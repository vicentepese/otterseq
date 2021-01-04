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
DupIIDs=$(jq -r '.file.DupIIDs' settings.json)
TripSNPS=$(jq -r '.file.TripSNPs' settings.json)

# Parse duplicate variants (based on position and allele codes)
plink --bfile $GWASDATA --list-duplicate-vars suppress-first \
    --allow-no-sex --out temp > temp
awk '{print $4}' temp.dupvar > $DupSNPs
rm -r temp*

# Parse duplicated  IIDs
awk 'a[$1]++{$1=$1" "$1; print $1}' $GWASDATA.fam > $DupIIDs

# Perform Quality control - Remove duplicated variants
plink --bfile $GWASDATA --remove $IBD_ID --exclude $DupSNPs\
    --allow-no-sex \
    --maf $MAF --geno $GENOMISS --mind $PHENOMISS \
    --make-bed --out gwastemp > gwastemp

# Perform Quality Control - Remove duplicated IIDs
plink --bfile gwastemp --remove $DupIIDs \
    --allow-no-sex \
    --make-bed --out gwastempFilt > gwastempFilt
rm -r gwastemp.*

# Parse triplicated variants / multiallelic variants
plink --bfile gwastempFilt --list-duplicate-vars\
    --allow-no-sex --out temp > temp
awk '{print $4}' temp.dupvar > $TripSNPS
rm -r temp*

# Remove triplicated / multiallelic variants
plink --bfile gwastempFilt --exclude $TripSNPS\
    --allow-no-sex \
    --make-bed --out $GWASDATAQC >> $GWASDATAQC
rm -r gwastemp*
