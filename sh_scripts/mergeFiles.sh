#!/bin/bash 

# Get directories, files and resources
GWASBIN=$(jq -r '.directory.GWAS_binaries' settings.json)
GWASMERGE=$(jq -r '.plinkFiles.GWAS' settings.json)
MERGELIST=$(jq -r '.file.mergeList' settings.json)
COMMONSNPS=$(jq -r '.file.commonSNPs' settings.json)

# Get list of files 
BINFILES=$(find $GWASBIN -iname '*.bim' -type f -exec sh -c 'printf "%s\n" "${0%.*}"' {} ';'| head -1)

# Get number of files to merge 
NUMFILES=$(awk 'END{print NR}' $MERGELIST)

# If more than one dataset, merge and filter - else, filter
if [ $NUMFILES -lt 1 ];
then 
    plink --bfile $BINFILES \
    --merge-list $MERGELIST \
    --extract $COMMONSNPS \
    --allow-no-sex \
    --memory 4626791360 \
    --make-bed \
    --out $GWASMERGE > $GWASMERGE
else
    cp $BINFILES.bed $GWASMERGE.bed
    cp $BINFILES.fam $GWASMERGE.fam
    cp $BINFILES.bim $GWASMERGE.bim
fi

