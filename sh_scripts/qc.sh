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
TripSNPS=$(jq -r '.file.TripSNPs' settings.json)

# Parse duplicate variants (based on position and allele codes)
plink --bfile $GWASDATA --list-duplicate-vars suppress-first \
    --out temp > temp
awk '{print $4}' temp.dupvar > $DupSNPs
rm -r temp*

# Perform Quality control - Remove duplicated variants
plink --bfile $GWASDATA --remove $IBD_ID --exclude $DupSNPs\
    --no-sex --no-parents --not-chr 25,26 \
    --maf $MAF --geno $GENOMISS --mind $PHENOMISS \
    --make-bed --out gwastemp >> gwastemp

# Remove triplicated variants / multiallelic variants
plink --bfile gwastemp --list-duplicate-vars\
    --out temp > temp
awk '{print $4}' temp.dupvar > $TripSNPS
rm -r temp*

plink --bfile gwastemp --remove $IBD_ID --exclude $TripSNPS\
    --no-sex --no-parents --not-chr 25,26 \
    --maf $MAF --geno $GENOMISS --mind $PHENOMISS \
    --make-bed --out $GWASDATAQC >> $GWASDATAQC
rm -r gwastemp*