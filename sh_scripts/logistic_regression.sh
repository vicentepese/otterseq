#!/bin/bash

###################################################################
#Script Name	: binarize.sh                                                                                            
#Description	: Computes a logistic regression through PLINK                                                                       
#Args           : None                                                                                           
#Author       	: Vicente Peris Sempere                                                
#Email         	: vipese@stanford.edu                                        
###################################################################

# Read settings 
GWASFILE=$(jq -r '.plinkFiles.GWASQC' settings.json)
PHENOFILE=$(jq -r '.file.pheno_matched' settings.json)
PCA=$(jq -r '.file.PCA_eigenvec' settings.json)
OUTPUT=$(jq -r '.directory.GWAS_out' settings.json)
PREFIX=$(jq -r '.plinkFiles.prefix' settings.json)

# Create covar file 
HEADER=("FID" "IID")
for i in $(seq 1 20)
do 
    HEADER+=("PC"$i)
done
printf "%s " "${HEADER[@]}" > tst
printf "\n" >> tst
awk '{
    print $0
}' $PCA >> tst
awk -F " " '{
    gsub(/ +/, "\t"); print
}' tst > covartemp  

# Run logistic regression -- association analysis
./bin/plink --bfile ${GWASFILE}${PREFIX}_QC --pheno $PHENOFILE \
    --covar covartemp --covar-name PC1, PC2, PC3, PC4 \
    --logistic --allow-no-sex --out ${OUTPUT}${PREFIX}

# Remove temporary covariates
rm -r covartemp
rm tst

