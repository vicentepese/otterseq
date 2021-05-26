#!/usr/bin/env python
"""Performs GWAS: Merge, Quality Control, PCA, Case matching, Logistic regression.

- Merge: from multiple binary .ped files, merge files using only common SNPs
- Quality Control (QC): filter based on minimum allele frquency, and missing genotype 
    and phenotype. Removes duplicated variants and IIDs/FIDs.
- PCA: Principal Component Analysis - verifies that there is no batch biases.
- Case Matching: Matches cases with controls based on a predefined ratio, 
    using PCs and Euclidian distance. 
- Logistic Regression: Computes a logistic regression controling for PC
"""

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

__author__ = "Vicente Peris Sempere"
__credits__ = ["Vicente Peris Sempere"]
__license__ = "GPL"
__version__ = "3.0"
__maintainer__ = "Vicente Peris Sempere"
__email__ = "vipese@stanford.edu"
__status__ = "Finalized"


def binarizeFiles(settings):
    """ Convert .ped to .bed files through PLINK
    Input: 
        - settings: settings JSON file

    Output: 
        - Binarized files
    """

    # Print 
    print("Binarizing files: \n")
    print(str(f) + ", " for f in os.listdir(settings['directory']['GWAS']))

    # Call bash script to binarize files in Data/GWAS
    subprocess.call(['bash', settings['sh_script']['binarize.sh']])

def get_SNP(settings):
    """ Get common SNPs binarized files
    Input:
        - settings: settings JSON file 
    Output:
        - Common SNPs text file
    """

    print("Parsing common SNPs")

    # Get .bim files of controls
    path = settings['directory']['GWAS_binaries']
    files = [file for file in os.listdir(path) if isfile(join(path, file)) and '.bim' in file]

    # Initialize
    totalSNPs = list()

    # Geet SNPs
    for file in files:
        SNPs = list()
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
    """Creates merge file text and merges files through PLINK including 
    common SNPs only.
    Input: 
        - settings: settings JSON file

    Output:
        - Merged files 
    """

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
    """Quality Control: computes Inheritance By Descendance, and removes 
    subjetcs based on IBD filter, MAF, missingnes and duplicates throug PLINK.
    Inputs:
        - settings: settings JSON file
    Outputs:
        - QCed files through PLINK
    """

    # Run IBD computation
    print('Computing IBD')
    subprocess.call(['bash',settings['sh_script']['IBD.sh']])
    print("IBD successfully computed")

    # Get list of patients to be removed by IBD
    IBD_IDs = list()
    with open(settings['file']['IBDGenome'], 'r') as inFile:
        next(inFile)
        for row in inFile:
            row = row.split('\t')
            row = [r for r in row if r != '']
            PI_hat = float(row[8])
            if float(PI_hat) > settings['IBD_threshold']:
                IBD_IDs.append(row[0])
       
    # Get list of patients and cases 
    patientList = list()
    with open(settings['plinkFiles']['GWAS'] + settings['plinkFiles']['prefix'] +'.fam','r') as inFile:
        for row in inFile:
            row = row.split()
            pheno = int(row[5].split('\n')[0])
            patientList.append(row[0])
    
    # Get list of IBDS patient/cases 
    high_IBD = [ID for ID in np.unique(IBD_IDs) if ID in patientList]
    
    # Check number of cases and controls
    pheno = pd.read_csv(settings['file']['pheno'], sep = ' ', header = None)
    pheno.columns = ['FID', 'IID', 'pheno']
    cases = pheno[pheno.pheno.eq(2)]
    controls = pheno[pheno.pheno.eq(1)]
    case_count = cases[cases["IID"].isin(high_IBD)].shape[0]
    control_count = controls[controls["IID"].isin(high_IBD)].shape[0]
    print('Number of cases to be excluded by IBD: ' + str(case_count))
    print('Number of controls to be excluded by IBD: ' + str(control_count))

    # Create exclusion file 
    with open(settings['file']['excludeID_IBD'], 'w') as outFile:
        for ID in high_IBD: 
            outFile.write(ID + ' ' + ID + ' ' + '\n')

    # Compute QC
    subprocess.call(['bash', settings['sh_script']['qc.sh']])

def computePCA(settings):
    """Computes PCA removing HLA region 
    Inputs:
        - settings: settings JSON file 
    Outputs:
        - Principal Components (PCs)
    """

    # Compute PCA
    print("Compute PCA")
    subprocess.call('bash ' + settings['sh_script']['pca.sh'], shell=True)
    print("PCA computed. Ploting PCs")

def plotPCA(settings, type = "batch"):
    """ Plots PCA. If type is batch, uses pheno provided. If type is match, uses matched subjects
    Inputs: 
        - settings: settings JSON file 
        - type: "batch", "match". Specifies the type of plot to be computed (matched, or unmatched)
    Outputs:
        - PCA plot
    """

    # Import PCA 
    PCA = pd.read_csv(settings['file']['PCA_eigenvec'], header = None, delim_whitespace = True)
    PCA.columns = ['FID', 'IID'] + ['PC' + str(x) for x in range(1,21)]

    # Merge with pheno -- batch: study batch biases
    if type == "batch":

        # Import pheno 
        pheno = pd.read_csv(settings['file']['pheno'], sep = ' ', header = None)
        pheno.columns = ['FID', 'IID', 'pheno']
        pheno['FID']=pheno['FID'].astype(object)

    # match: plot PCA of matched cases /subjects
    elif type == "match":

        # Import matched pheno
        pheno = pd.read_csv(settings['file']['pheno_matched'], sep = ' ', header = None)
        pheno.columns = ['FID', 'IID', 'pheno']
        pheno['FID']=pheno['FID'].astype(object)

    # Merge with PCA
    PCA = pd.merge(PCA, pheno, on='FID')

    # Plot PCs
    Ncases = PCA[PCA.pheno.eq(2)].shape[0]
    Ncontrols = PCA[PCA.pheno.eq(1)].shape[0]
    sns.set_style()
    ax = sns.relplot(data=PCA, x = 'PC1', y = 'PC2', hue = "pheno")
    ax.set(xlabel = "PC1", ylabel = "PC2", title = "PCA -- Cases: " + str(Ncases) + " // Controls " + str(Ncontrols))
    plt.show()

def patientMatching(settings):
    """Matches cases with its closes controls based on provided ratio, and Euclidian Distance
    of PCs.
    Inputs:
        - settings: settings JSON file 
    Outputs:
        - File with list of cases and controls matched
    """

    # Perform patient matching
    print("Computing patient matching as the Euclidean distance between each case and control")

    # Get PCs, cases and controls
    PCA = pd.read_csv(settings['file']['PCA_eigenvec'], header = None, delim_whitespace = True)
    PCA.columns = ['FID', 'IID'] + ['PC' + str(x) for x in range(1,21)]
    pheno = pheno = pd.read_csv(settings['file']['pheno'], sep = ' ', header = None)
    pheno.columns = ['FID', 'IID', 'pheno']
    cases = pheno[pheno.pheno.eq(2)]
    controls = pheno[pheno.pheno.eq(1)]

    # Filter cases and controls present in the data 
    cases = cases[cases["FID"].isin(PCA["FID"]) & cases["IID"].isin(PCA["IID"])]
    controls = controls[controls["FID"].isin(PCA["FID"]) & controls["IID"].isin(PCA["IID"])]

    # Get matching patients ratio
    ratio = settings['ControlCaseRatio']

    # For each case, compute euclidean distance
    matched_controls = pd.DataFrame()
    for idx, (FID, IID) in enumerate(zip(cases['FID'], cases['IID'])):

        # Print
        print("Matching subject {0} of {1}".format(idx, cases.shape[0]-1), end='\r', flush=True)

        # Get the index and PCs of the case. Store.
        PC_case = PCA[PCA.FID.eq(FID) & PCA.IID.eq(IID)].drop(['FID', 'IID'], axis = 1)

        # Compute euclidean distance against all controls 
        patEuDist = defaultdict(list)
        for FID_ctrl, IID_ctrl in zip(controls['FID'], controls['IID']):

            # Get controls PCs and compute euclidean distance
            PC_ctrl = PCA[PCA.FID.eq(FID_ctrl) & PCA.IID.eq(IID_ctrl)].drop(['FID', 'IID'], axis = 1) 
            patEuDist['FID_case'].append(FID); patEuDist['IID_case'].append(IID)
            patEuDist['FID_control'].append(FID_ctrl); patEuDist['IID_control'].append(IID_ctrl)
            patEuDist['dist'].append(euclidean(PC_case[['PC1','PC2']], PC_ctrl[['PC1','PC2']]))

        # Sort DataFrame
        patEuDist = pd.DataFrame(patEuDist)
        patEuDist= patEuDist.sort_values(by = "dist",  ascending=True)

        # Get #ratio# of matched controls and append to DataFrame
        matched_controls = matched_controls.append(patEuDist[['FID_control','IID_control']][0:ratio])

    # Drop duplicated and append pheno
    matched_controls = matched_controls.drop_duplicates(subset = ["FID_control", "IID_control"])
    matched_controls['pheno'] = [1]*matched_controls.shape[0]

    # Create patient list
    patList = matched_controls
    patList.columns = ['FID', 'IID', 'pheno']
    patList = patList.append(cases)
        
    # Write patient list to csv file 
    patList.to_csv(settings['file']['pheno_matched'], sep = ' ', header = None, index = False)

    # Print 
    print ("Cases successfully matched. Matched controls: " +  str(patList[patList.pheno.eq(1)].shape[0]))

def logistic_regression(settings):
    """ Computes the logistic regression of QCed files and matched cases
    Inputs:
        - settings: settings JSON file 
    Outputs: 
        - PLINK file with logistic regression results
    """

    # Run logistic regression 
    print ("Computing logistic regression")
    subprocess.call("bash " + settings['sh_script']["logistic_regression.sh"], shell=True)
    print("Logistic regression successfully computed")


def main():
    """ Main function of the script. 
    """

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
    get_SNP(settings)

    # Merge files based on common SNPs
    mergeFiles(settings)

    # Quality control (QC) + IBD filtering
    QC(settings)

    # Compute PCA 
    computePCA(settings)

    # Plot PCA
    plotPCA(settings, type = "batch")

    # Compute case-control matching 
    patientMatching(settings)

    # Compute PCA 
    plotPCA(settings, type = "match")

    # Perform association analysis - logistic regression
    logistic_regression(settings)


if __name__ == "__main__":
    main()