#!bin/bash 

###################################################################
#Script Name	: binarize.sh                                                                                            
#Description	: Perform Quality Control of .bed files through PLINK                                                                      
#Args           : None                                                                                           
#Author       	: Vicente Peris Sempere                                                
#Email         	: vipese@stanford.edu                                        
###################################################################

# Read from settings
    # PLINK files
GWASDATA=$(jq -r '.plinkFiles.GWAS' settings.json)
GWASDATAQC=$(jq -r '.plinkFiles.GWASQC' settings.json)
PREFIX=$(jq -r '.plinkFiles.prefix' settings.json)
    # Variables for QC
IBD_ID=$(jq -r '.file.excludeID_IBD' settings.json)
PHENOMISS=$(jq -r '.phenomiss' settings.json)
GENOMISS=$(jq -r '.genomiss' settings.json)
MAF=$(jq -r '.maf' settings.json)
DupSNPs=$(jq -r '.file.DupSNPs' settings.json)
DupIIDs=$(jq -r '.file.DupIIDs' settings.json)
TripSNPS=$(jq -r '.file.TripSNPs' settings.json)

# Parse duplicate variants (based on position and allele codes)
./bin/plink --bfile ${GWASDATA}${PREFIX} --list-duplicate-vars suppress-first \
    --allow-no-sex --out temp > temp
awk '{print $4}' temp.dupvar > $DupSNPs
rm -r temp*

# Parse duplicated  FIDs
awk '{seen[$1,$2]++}' ${GWASDATA}${PREFIX}.fam > $DupIIDs

# Concat IBD remove and DupSNPS
cp $IBD_ID temp_remove
awk '{print $0}' $DupIIDs >> temp_remove

# Perform Quality control - Remove duplicated variants
./bin/plink --bfile ${GWASDATA}${PREFIX} --remove temp_remove --exclude $DupSNPs\
    --allow-no-sex \
    --maf $MAF --geno $GENOMISS --mind $PHENOMISS \
    --make-bed --out gwastempFilt > gwastempFilt
rm temp_remove

# Parse triplicated variants / multiallelic variants
./bin/plink --bfile gwastempFilt --list-duplicate-vars\
    --allow-no-sex --out temp > temp
awk '{print $4}' temp.dupvar > $TripSNPS
rm -r temp*

# Remove triplicated / multiallelic variants
./bin/plink --bfile gwastempFilt --exclude $TripSNPS\
    --allow-no-sex \
    --make-bed --out ${GWASDATAQC}${PREFIX}_QC >>${GWASDATAQC}${PREFIX}_QC
rm -r gwastemp*

# Clean .bim file to remove commas (PCA)
awk '{
    gsub(/\,/, ":", $2); print $0
    }' ${GWASDATAQC}${PREFIX}_QC.bim > temp.bim
mv temp.bim ${GWASDATAQC}${PREFIX}_QC.bim