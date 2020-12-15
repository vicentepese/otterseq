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
DupSNPs=$(jq -r '.file.DupSNPs' settings.json)
TripSNPS=$(jq -r '.file.TripSNPs' settigs.json)

# Remove duplicate variants (based on position and allele codes)
plink --bfile $GWASDATA --list-duplicate-vars\
    --out temp
awk '{print $4}' temp.dupvar > $DupSNPs
rm -r temp*

# Perform Quality control
plink --bfile $GWASDATA --remove $IBD_ID --exclude $DupSNPs\
    --no-sex --no-parents --not-chr 25,26 \
    --maf $MAF --geno $GENOMISS --mind $PHENOMISS \
    --make-bed --out $GWASDATAQC >> $GWASDATAQC

# Remove triplicated variants / multiallelic variants
plink --bfile $GWASDATAQC --list-duplicate-vars\
    --out temp
awk '{print $4}' temp.dupvar > $TripSNPs
rm -r temp*

plink --bfile $GWASDATA --remove $IBD_ID --exclude $TripSNPs\
    --no-sex --no-parents --not-chr 25,26 \
    --maf $MAF --geno $GENOMISS --mind $PHENOMISS \
    --make-bed --out $GWASDATAQC >> $GWASDATAQC