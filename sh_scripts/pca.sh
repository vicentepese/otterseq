#!/bin/bash

# Read variables from settings
PHENOMISS=$(jq -r '.phenomiss' settings.json)
GENOMISS=$(jq -r '.genomiss' settings.json)
MAF=$(jq -r '.maf' settings.json)
GWASFILE=$(jq -r '.plinkFiles.GWAS' settings.json)
EXCLUDEIBD=$(jq -r '.file.excludeID_IBD' settings.json)
PCA=$(jq -r '.plinkFiles.PCA' settings.json)

# Parse SNPs (rs_xxxxxxxx) belonging to HLA region 
plink --bfile $GWASFILE --no-fid --no-sex --no-parents --chr 6 --from-bp 28477797 --to-bp 33448354 --make-bed --out temp > temp
rsExclude=$(awk '{ ORS=", "};{print $2}' temp.bim)
rm -r temp*

# Compute PCA excluding by IBD, missingness, and HLA region
plink --bfile $GWASFILE --remove $EXCLUDEIBD \
    --no-fid --no-sex --no-parents \
    --exclude-snps $rsExclude --d ';' --not-chr XY \
    --maf $MAF --geno $GENOMISS --mind $PHENOMISS \
    --pca 20  --out $PCA > $PCA.log

