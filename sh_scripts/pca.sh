#!/bin/bash

###################################################################
#Script Name	: binarize.sh                                                                                            
#Description	: Computes Principal Component Analysis of .bed files                                                                          
#Args           : None                                                                                           
#Author       	: Vicente Peris Sempere                                                
#Email         	: vipese@stanford.edu                                        
###################################################################

# Read variables from settings
PHENOMISS=$(jq -r '.phenomiss' settings.json)
GENOMISS=$(jq -r '.genomiss' settings.json)
MAF=$(jq -r '.maf' settings.json)
GWASQC=$(jq -r '.plinkFiles.GWASQC' settings.json)
EXCLUDEIBD=$(jq -r '.file.excludeID_IBD' settings.json)
PCA=$(jq -r '.plinkFiles.PCA' settings.json)
PREFIX=$(jq -r '.plinkFiles.prefix' settings.json)


# Parse SNPs (rs_xxxxxxxx) belonging to HLA region (EXCLUDE VARS with ,?)
./bin/plink --bfile ${GWASQC}${PREFIX}_QC --allow-no-sex --chr 6 --from-bp 28477797 --to-bp 33448354 --make-bed --out temp > temp
awk '{print $2}' temp.bim > exclude.txt
rm -r temp*

# Compute PCA excluding by IBD, missingness, and HLA region
./bin/plink --bfile ${GWASQC}${PREFIX}_QC \
    --allow-no-sex \
    --exclude exclude.txt \
    --pca 20  --out $PCA
rm exclude.txt
