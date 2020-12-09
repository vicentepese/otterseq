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

def main():

    # Open settings
    with open('settings.json', 'r') as inFile:
        settings = json.load(inFile)

    # Add current directory to settings 
    settings["directory"]["main"] = os.getcwd()

    # Binarize GWAS files
    if (os.listdir(settings['directory']['GWAS_binaries'])):
        binarizeFiles(settings)
    else:
        print("Files already binarized")

    # 


    

if __name__ == "__main__":
    
    main()