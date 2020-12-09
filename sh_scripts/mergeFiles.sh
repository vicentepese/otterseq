#!/bin/bash 

# Get directories, files and resources
GWASBIN=$(jq -r '.directory.GWAS_binaries' settings.json)
GWASMERGE=$(jq -r '.plinkFiles.GWAS' settings.json)
MERGELIST=$(jq -r '.file.mergeList' settings.json)

# Get list of files 
BINFILES=$(find $GWASBIN -iname '*.bim' -type f -exec sh -c 'printf "%s\n" "${0%.*}"' {} ';'| head -1)

# Merge 
plink --bfile $BINFILES --merge-list $MERGELIST --out $GWASMERGE >> temp
rm temp

