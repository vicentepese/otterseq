import numpy as np
import subprocess
import json
import os
from os.path import isfile, join
from collections import defaultdict, OrderedDict
import matplotlib.pyplot as plt
import scipy 
from scipy.spatial.distance import euclidean 
import seaborn as sns
import pandas as pd


def binarizeFiles(settings):

    # Print 
    print("Binarizing files: \n")
    print(str(f) + ", " for f in os.listdir(settings['directory']['GWAS']))

    # Call bash script to binarize files in Data/GWAS
    subprocess.call(['bash', settings['sh_script']['binarize.sh']])

def get_SNP(settings, path):

    print("Parsing common SNPs")

    # Get .bim files of controls
    files = [file for file in os.listdir(path) if isfile(join(path, file)) and '.bim' in file]

    # Initialize
    totalSNPs = list()
    SNPs = list()

    # Geet SNPs
    for file in files:
        print("Parsing SNPs from " + file)
        with open(join(path, file), 'r') as inFile:
            for row in inFile:
                row = row.split()
                SNPs.append(row[1])
            totalSNPs.append(SNPs)
    print('Finished parsing')

    # Get intersection
    commonSNPs = set(totalSNPs[0])
    for i in range(1, len(totalSNPs)):
        commonSNPs = commonSNPs.intersection(set(totalSNPs[i]))
    commonSNPs = list(commonSNPs)

    # Write common SNPS
    with open(settings['file']['commonSNPs'], 'w') as outFile:
        for snp in commonSNPs:
            outFile.write(snp + '\n')

    # Print 
    print("Common SNPs stored in " + settings['file']['commonSNPs'])

def mergeFiles(settings):
    # Creates a merge file and uses it to merge files in GWAS_binaries

    # Get control files
    filePath = settings['directory']['GWAS_binaries']
    files = [f for f in os.listdir(filePath) if isfile(join(filePath,f)) and '.bim' in f]
    files = [f.split('.')[0] for f in files]
    files = [settings['directory']['GWAS_binaries'] + f for f in files]

    # Write mergelist
    with open(settings['file']['mergeList'],'w') as outFile:
        for f in files:
            outFile.write(f + '.bed ' + f + '.bim ' + f + '.fam' + '\n')

    # Merge files
    print('Merging files')
    subprocess.call(['bash', settings['sh_script']['mergeFiles.sh']])
    print('Files successfully merged') 

def QC(settings):
    # This function computes and IBD computation and creates 
    # A list of subjects to remove removeIDs.txt in Resources 

    # Run IBD computation
    print('Computing IBD')
    subprocess.call(['bash',settings['sh_script']['IBD.sh']])
    print("IBD successfully computed")

    # Get list of patients to be removed 
    IBD_IDs = list()
    with open(settings['file']['IBDGenome'], 'r') as inFile:
        next(inFile)
        for row in inFile:
            row = row.split('\t')
            row = [r for r in row if r is not '']
            PI_hat = row[8]
            if float(PI_hat) > settings['IBD_threshold']:
                IBD_IDs.append(row[0])
       
    # Get list of patients and cases (remove unknown pheno)
    patientList = list()
    with open(settings['plinkFiles']['GWAS'] + '.fam','r') as inFile:
        for row in inFile:
            row = row.split()
            pheno = int(row[5].split('\n')[0])
            if pheno != -9:
                patientList.append(row[0])
    
    # Get list of IBDS patient/cases (unknown pheno not included)
    high_IBD = [ID for ID in np.unique(IBD_IDs) if ID in patientList]
    
    # Check number of cases and controls
    pheno = pd.read_csv(settings['file']['pheno'], sep = ' ', header = None)
    pheno.columns = ['IID', 'FID', 'pheno']
    cases = pheno[pheno['pheno'] == 1]
    controls = pheno[pheno['pheno'] == 0]
    case_count = len([subj for subj in high_IBD if subj in cases['IID']])
    control_count = len([subj for subj in high_IBD if subj in controls['IID']])
    print('Number of cases to be excluded by IBD: ' + str(case_count))
    print('Number of controls to be excluded by IBD: ' + str(control_count))

    # Create exclusion file 
    with open(settings['file']['excludeID_IBD'], 'w') as outFile:
        for ID in high_IBD: 
            outFile.write(ID + ' ' + ID + ' ' + '\n')

    # Compute QC
    subprocess.call(['bash', settings['sh_script']['qc.sh']])

def computePCA(settings, type):

    # Compute PCA
    print("Compute PCA")
    subprocess.call('bash ' + settings['sh_script']['pca.sh'] + ' -t ' + type, shell=True)
    print("PCA computed. Ploting PCs")

    # Import PCA eigenvectors
    PCA = pd.read_csv(settings['file']['PCA_eigenvec'], header = None, delim_whitespace = True)
    PCA.columns = ['IID', 'FID'] + ['PC' + str(x) for x in range(1,21)]

    # Merge with pheno
    pheno = pd.read_csv(settings['file']['pheno'], sep = ' ', header = None)
    pheno.columns = ['IID', 'FID', 'pheno']
    pheno['IID']=pheno['IID'].astype(object)
    tst = pd.merge(PCA, pheno, on='IID')

    # Plot PCs
    Nsubjs = PCA.shape[0]
    sns.set_style()
    ax = sns.relplot(data=tst, x = 'PC1', y = 'PC2', hue = "pheno")
    ax.set(xlabel = "PC1", ylabel = "PC2", title = "PCA / #Subjects: " + str(Nsubjs))
    plt.show()

def patientMatching(settings):

    # Perform patient matching
    print("Computing patient matching as the Euclidean distance between each case and control")

    # Get PCs, cases and controls
    PCA = pd.read_csv(settings['file']['PCA_eigenvec'], header = None, delim_whitespace = True)
    PCA.columns = ['IID', 'FID'] + ['PC' + str(x) for x in range(1,21)]
    pheno = pheno = pd.read_csv(settings['file']['pheno'], sep = ' ', header = None)
    pheno.columns = ['IID', 'FID', 'pheno']
    cases = pheno[pheno['pheno'] == 1]
    controls = pheno[pheno['pheno'] == 0]

    # Get matching patients ratio
    ratio = settings['ControlCaseRatio']

    # For each case, compute euclidean distance
    # For each case, compute euclidean distance
    matched_controls = pd.DataFrame()
    for idx, (IID, FID) in enumerate(zip(cases['IID'], cases['FID'])):

        # Print
        print("\rMatching subject {0} of {1} \r".format(idx, cases.shape[0]))

        # Get the index and PCs of the case. Store.
        PC_case = PCA[PCA.IID.eq(IID) & PCA.FID.eq(FID)].drop(['IID', 'FID'], axis = 1)

        # Compute euclidean distance against all controls 
        patEuDist = defaultdict(list)
        for idx_ctrl, IID_ctrl, FID_ctrl in zip(controls.index, controls['IID'], controls['FID']):

            # Get controls PCs and compute euclidean distance
            PC_ctrl = PCA[PCA.IID.eq(IID_ctrl) & PCA.FID.eq(FID_ctrl)].drop(['IID', 'FID'], axis = 1) 
            patEuDist['IID_case'].append(IID); patEuDist['FID_case'].append(FID)
            patEuDist['IID_control'].append(IID_ctrl); patEuDist['FID_control'].append(FID_ctrl)
            patEuDist['dist'].append(euclidean(PC_case, PC_ctrl))

        # Sort DataFrame
        patEuDist = pd.DataFrame(patEuDist)
        patEuDist= patEuDist.sort_values(by = "dist",  ascending=True)

        # Get #ratio# of matched controls and append to DataFrame
        matched_controls = matched_controls.append(patEuDist[['IID_control','FID_control']][0:50])

    # Drop duplicated and append pheno
    matched_controls.drop_duplicates()
    matched_controls['pheno'] = [0]*matched_controls.shape[0]

    # Create patient list
    patList = matched_controls
    patList.columns = ['IID', 'FID', 'pheno']
    patList = patList.append(cases)
        
    # Write patient list to csv file 
    patList.to_csv('test.txt', sep = '\t', header = None)

def main():

    # Open settings
    with open('settings.json', 'r') as inFile:
        settings = json.load(inFile)

    # Add current directory to settings 
    settings["directory"]["main"] = os.getcwd()

    # Binarize GWAS files
    if (not os.listdir(settings['directory']['GWAS_binaries'])):
        binarizeFiles(settings)
    else:
        print("Files already binarized")

    # Get list of common SNPs across files 
    get_SNP(settings, path = settings['directory']['GWAS_binaries'])

    # Merge files based on common SNPs
    mergeFiles(settings)

    # Quality control (QC) + IBD filtering
    QC(settings)

    # Compute PCA 
    computePCA(settings, type = "batch")

    # Compute case-control matching 
    

    

if __name__ == "__main__":
    
    main()