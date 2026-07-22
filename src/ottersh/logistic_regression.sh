#!/bin/bash

# Parse named arguments
while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        --bfile)
        bfile="$2"
        shift # past argument
        shift # past value
        ;;
        --outpath)
        outpath="$2"
        shift # past argument
        shift # past value
        ;;
        --pheno)
        pheno="$2"
        shift # past argument
        shift # past value
        ;;
        --covar)
        covar="$2"
        shift # past argument
        shift # past value
        ;;
        *)    # unknown option
        echo "Unknown option $key"
        exit 1
        ;;
    esac
done


# Run logistic regression -- association analysis
bin/plink --bfile "${bfile}" \
      --pheno "${pheno}" \
      --covar "${covar}" \
      --logistic \
      --hide-covar \
      --allow-no-sex \
      --out "${outpath}"

