import numpy as np
import subprocess
import json
import os
from os.path import isfile, join
from collections import defaultdict
import matplotlib.pyplot as plt
import scipy 
from scipy.spatial.distance import euclidean 

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
                row = row.split('\t')
                SNPs.append(row[1])
            totalSNPs.append(SNPs)
    print('Finished parsing')

    # Get intersection
    commonSNPs = set(totalSNPs[0])
    for i in range(1, len(totalSNPs)):
        commonSNPs = commonSNPs.intersection(set(totalSNPs[i]))
    commonSNPs = list(commonSNPs)

    # Write common SNPS
    with open(settings['file']['commonSNPs.txt'], 'w') as outFile:
        for snp in commonSNPs:
            outFile.write(snp + '\n')

    # Print 
    print("Common SNPs stored in " + settings['file']['commonSNPs.txt'])

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


def addPhenotype(settings, f, sep):
    # Opens filtLGI1 and adds the phenotype based on the list of 
    # controls and cases. If not in cases or controls, ID is added 
    # to list and written in exclude.txt
    # Output modLGI1.fam 

    # Read GWAS IDS
    controlsID = list() 
    with open(settings['file']['GWASIDsControls'], 'r') as inFile:
        for row in inFile:
            controlsID.append(row.split('"')[1])
    casesID = list()
    with open(settings['file']['GWASIDsCases'], 'r') as inFile:
        for row in inFile:
            casesID.append(row.split('"')[1])

    # Feed phenotype to fam file and create exclude list
    exclude = list()
    with open(settings['file']['modLGI1.fam'], 'w') as outFile:
        with open(f + '.fam', 'r') as inFile:
            for row in inFile:
                row = row.split(sep)
                subjectID = row[0].split('_')[3].split('.')[0]
                if subjectID in controlsID:
                    row[-1] = '1'
                    outFile.write('\t'.join(row) + '\n')
                elif subjectID in casesID:
                    row[-1] = '2'
                    outFile.write('\t'.join(row) + '\n')
                else:
                    exclude.append(row[0])
                    outFile.write('\t'.join(row))  

def IBDfilt(settings):
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
    patientList = list(); 
    case_count, control_count = 0, 0
    with open(settings['plinkFiles']['GWAS'] + '.fam','r') as inFile:
        for row in inFile:
            row = row.split('\t')
            pheno = int(row[5].split('\n')[0])
            if pheno != -9:
                patientList.append(row[0])
    
    # Get list of IBDS patient/cases (unknown pheno not included)
    high_IBD = [ID for ID in np.unique(IBD_IDs) if ID in patientList]
    
    # # Check number of cases and controls
    # case_count = 0 ; control_count = 0 
    # for ID in high_IBD:
    #     ID_pre = ID.split('_')[-1].split('.')[0]
    #     if ID_pre in cases:
    #         case_count +=1
    #     elif ID_pre in controls:
    #         control_count += 1
    # print('Number of cases to be excluded: ' + str(case_count))
    # print('Number of controls to be excluded: ' + str(control_count))

    # Create exclusion file 
    with open(settings['file']['excludeID.txt'], 'w') as outFile:
        for ID in high_IBD: 
            outFile.write(ID + ' ' + ID + ' ' + '\n')


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

    # Merge files
    mergeFiles(settings)

    # Filter IBD
    IBDfilt(settings)


    

if __name__ == "__main__":
    
    main()