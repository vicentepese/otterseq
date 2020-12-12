#!bin/bash 

# Read from settings
GWASDATA=$(jq -r '.plinkFiles.GWAS' settings.json)
GWASDATAQC=$(jq -r 'plinkFiles.GWAS' settings.json)

# Perform Quality control
plink 