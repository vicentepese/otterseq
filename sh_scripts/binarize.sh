#!/bin/bash 

# Get directory and files
MAINDIR=$(pwd)
GWASDIR=$(jq -r '.directory.GWAS' settings.json)
GWASBINDIR=$(jq -r '.directory.GWAS_binaries' settings.json)

# For each file, binarize

GWASFILES=$(ls $GWASDIR)
for file in ${GWASFILES[@]} ; do
    if [[ $file == *'.ped'* ]] ; then
        echo "Converting $file to binary"
        IFS='.' read -a strarr <<< "$file"
        plink --file ${GWASDIR}${strarr[0]} --no-sex --no-pheno --no-fid --no-parents \
            --noweb --make-bed --out ${GWASBINDIR}${strarr[0]} >> ${GWASBINDIR}${strarr[0]}.log
    fi  
done 

# Print 
echo "Files binarized"


