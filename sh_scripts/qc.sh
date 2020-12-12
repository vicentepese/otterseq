#!bin/bash 

# Read from settings
GWASDATA=$(jq -r '.plinkFiles.GWAS' settings.json)
