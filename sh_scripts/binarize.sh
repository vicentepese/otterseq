#!/bin/bash 

###################################################################
#Script Name	: binarize.sh                                                                                            
#Description	: Converts .ped to .bed                                                                              
#Args           : None                                                                                           
#Author       	: Vicente Peris Sempere                                                
#Email         	: vipese@stanford.edu                                        
###################################################################

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
        ./bin/plink --file ${GWASDIR}${strarr[0]} --allow-no-sex \
            --noweb --make-bed --out ${GWASBINDIR}${strarr[0]} >> ${GWASBINDIR}${strarr[0]}.log
    fi  
done 

# Print 
echo "Files binarized"


