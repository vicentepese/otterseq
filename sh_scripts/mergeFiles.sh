#!/bin/bash 

# Get directories, files and resources
GWASBIN=$(jq -r '.directory.GWAS_binaries' settings.json)
GWASMERGE=$(jq -r '.plinkFiles.GWAS' settings.json)
MERGELIST=$(jq -r '.file.mergeList.txt' settings.json)

# Get list of files 
BINFILES=$(ls $GWASBIN*.bim | head -1)
BINFILES=${BINFILES[1]}

# Merge 
plink --bfile $BINFILES --merge-list $MERGELIST --out $GWASMERGE

