#!/bin/bash

# Read variables from settings
PHENOMISS=$(jq -r '.phenomiss' settings.json)
GENOMISS=$(jq -r '.genomiss' settings.json)
MAF=$(jq -r '.maf' settings.json)
GWASQC=$(jq -r '.plinkFiles.GWASQC' settings.json)
EXCLUDEIBD=$(jq -r '.file.excludeID_IBD' settings.json)
PCA=$(jq -r '.plinkFiles.PCA' settings.json)
PREFIX=$(jq -r '.plinkFiles.prexi' settings.json)


# Parse SNPs (rs_xxxxxxxx) belonging to HLA region (EXCLUDE VARS with ,?)
plink --bfile ${GWASQC}${PREFIX} --allow-no-sex --chr 6 --from-bp 28477797 --to-bp 33448354 --make-bed --out temp > temp
rsExclude=$(awk '{ ORS=", "};{print $2}' temp.bim)
rm -r temp*

# Compute PCA excluding by IBD, missingness, and HLA region
plink --bfile ${GWASQC}${PREFIX} \
    --allow-no-sex \
    --exclude-snps $rsExclude --d ';' \
    --pca 20  --out $PCA


