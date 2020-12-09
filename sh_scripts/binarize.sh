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
if [[ $DBCOUNT -gt 1 ]] ; then 
    DBLIST=$(ls $GWASDIR* -d)
    for DB in ${DBLIST[@]} ; do 
        echo "Binarizing files in "$DB
        GWASFILES=$(ls $DB)
        for file in ${GWASFILES[@]} ; do
            if [[ $file == *'.ped'* ]] ; then
                echo "Converting $file to binary"
                IFS='.' read -a strarr <<< "$file"
                plink --file ${DB}/${strarr[0]} --no-sex --no-pheno --no-fid --no-parents \
                    --noweb --make-bed --out ${DB}/${strarr[0]} >> ${DB}/${strarr[0]}.log
            fi  
        done 
    done
else 
    GWASFILES=$(ls)
    for file in ${GWASFILES[@]} ; do
        if [[ $file == *'.ped'* ]] ; then
                echo "Converting $file to binary"
                IFS='.' read -a strarr <<< "$file"
                plink --file ${DB}/${strarr[0]} --no-sex --no-pheno --no-fid --no-parents \
                    --noweb --make-bed --out ${DB}/${strarr[0]} >> ${DB}/${strarr[0]}.log
        fi  
    done 
fi


# Print 
echo "Files binarized"


