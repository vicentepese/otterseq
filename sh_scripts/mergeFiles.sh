#!/bin/bash 

# Get directories, files and resources
GWASBIN=$(jq -r '.directory.GWAS_binaries' settings.json)
GWASMERGE=$(jq -r '.plinkFiles.GWAS' settings.json)
MERGELIST=$(jq -r '.file.mergeList' settings.json)
COMMONSNPS=$(jq -r '.file.commonSNPs' settings.json)

# Get list of files 
BINFILES=$(find $GWASBIN -iname '*.bim' -type f -exec sh -c 'printf "%s\n" "${0%.*}"' {} ';'| head -1)

# Merge and create temporal file
plink --bfile $BINFILES \
    --merge-list $MERGELIST \
    --extract $COMMONSNPS \
    --allow-no-sex \
    --memory 4604003776 \
    --make-bed \
    --out $GWASMERGE > $GWASMERGE

# Take only common SNPS from temporal file and create merged GWAS files
plink --bfile temp \
    --allow-no-sex  \
    --extract $COMMONSNPS \
    --memory 4604003776 \
    --make-bed --out $GWASMERGE >> $GWASMERGE




# Delete temporal files
rm *temp*