#!bin/bash 

# Read from settings
    # Plink files
GWASDATA=$(jq -r '.plinkFiles.GWAS' settings.json)
GWASDATAQC=$(jq -r '.plinkFiles.GWASQC' settings.json)
    # Variables for QC
IBD_ID=$(jq -r '.file.excludeID_IBD' settings.json)
PHENOMISS=$(jq -r '.phenomiss' settings.json)
GENOMISS=$(jq -r '.genomiss' settings.json)
MAF=$(jq -r '.maf' settings.json)

# Perform Quality control
plink --bfile $GWASDATA --remove $IBD_ID \
    --no-fid --no-sex --no-parents --not-chr 25,26 \
    --maf $MAF --geno $GENOMISS --mind $PHENOMISS \
    --make-bed --out $GWASDATAQC >> $GWASDATAQC

