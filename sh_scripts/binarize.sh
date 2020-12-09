#!/bin/bash 

# Get directory and files
MAINDIR=$(pwd)
GWASDIR=$(jq -r '.directory.GWAS' settings.json)
GWASBINDIR=$(jq -r '.directory.GWAS_binaries' settings.json)

# Check if different datasets
DBCOUNT=$(ls -l $GWASDIR | grep -c ^d)

# If subdatasets in GWAS folder, create directories in GWAS_binaries
if [[ $DBCOUNT -gt 1 ]] ; then 

    # Go to directory, get list of datasets
    cd $GWASDIR
    DBLIST=$(ls */ -d)
    cd $MAINDIR

    # Creat directories
    for DB in ${DBLIST[@]} ; do 
        [[ -d $GWASBINDIR/$DB ]] || mkdir $GWASBINDIR/$DB
    done 
fi 

# For each file, binarize
for file in ${GWASFILES[@]} ; do
    if [[ $file == *'.ped'* ]] ; then
        echo "Converting $file to binary"
        IFS='.' read -a strarr <<< "$file"
        plink --file ${GWASDIR}${strarr[0]} --no-sex --no-pheno --no-fid --no-parents \
            --noweb --make-bed --out ${GWASBINDIR}${strarr[0]} >> ${GWASBINDIR}${strarr[0]}.log
        mv ${GWASBINDIR}${strarr[0]}.log ${GWASBINDIR}Log/${strarr[0]}.log
    fi  
done 

# Print 
echo "Files binarized"


