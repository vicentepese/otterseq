#!/bin/bash

# Read variables from settings
PHENOMISS=$(jq -r '.phenomiss' settings.json)
GENOMISS=$(jq -r '.genomiss' settings.json)
MAF=$(jq -r '.maf' settings.json)
GWASFILE=$(jq -r '.plinkFiles.GWASQC' settings.json)
EXCLUDEIBD=$(jq -r '.file.excludeID_IBD' settings.json)
PCA=$(jq -r '.plinkFiles.PCA' settings.json)

PROGNAME=$0

usage() {
  cat << EOF >&2
Usage: $PROGNAME [-p <path>]

One or more of the following flags is missing:

-t <type>: Type of PCA. Options:
        batch: Performs PCA removing CHR 6 to study any batches biases.
        match: Performs PCA after case-control matching.
EOF
  exit 1
}

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -h|--help)
    echo '-t <type>: Type of PCA. Options:'
    echo '        batch: Performs PCA removing CHR 6 to study any batches biases.'
    echo '        match: Performs PCA after case-control matching.'
    shift
    shift
    ;;
    -t|--type)
    if [ "$2" == "batch" ]
    then 
        TYPE=batch
    elif [ "$2" == "match" ]
    then
        TYPE=match
    else
        echo "Invalid option, please see --help for more information."
        exit 1
    fi
    shift # past argument
    shift # past value
    ;;
    *)    # unknown option
    usage # save it in an array for later 
esac
done

if [ "$TYPE" == batch ]
then 
    # Parse SNPs (rs_xxxxxxxx) belonging to HLA region (EXCLUDE VARS with ,?)
    plink --bfile $GWASFILE --allow-no-sex --chr 6 --from-bp 28477797 --to-bp 33448354 --make-bed --out temp > temp
    rsExclude=$(awk '{ ORS=", "};{print $2}' temp.bim)
    rm -r temp*

    # Compute PCA excluding by IBD, missingness, and HLA region
    plink --bfile $GWASFILE \
        --allow-no-sex \
        --exclude-snps $rsExclude --d ';' \
        --pca 20  --out $PCA
else

    # Parse SNPs (rs_xxxxxxxx) belonging to HLA region (EXCLUDE VARS with ,?)
    plink --bfile $GWASFILE --allow-no-sex --chr 6 --from-bp 28477797 --to-bp 33448354 --make-bed --out temp > temp
    rsExclude=$(awk '{ ORS=", "};{print $2}' temp.bim)
    rm -r temp*

    # Matched cases-controls 
    PATLIST=$(jq -r '.file.pheno_matched' settings.json)
    PCAMATCH=$(jq -r '.plinkFiles.PCA_matched' settings.json)

    # Compute PCA excluding by IBD, missingness, and HLA region
    plink --bfile $GWASFILE \
        --allow-no-sex \
        --pheno $PATLIST \
        --exclude-snps $rsExclude --d ';' \
        --pca 20 --out $PCAMATCH
fi

