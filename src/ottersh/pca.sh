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
        --exclude-hla)
        exclude_hla="$2"
        shift # past argument
        shift # past value
        ;;
        --pcs)
        pcs="$2"
        shift # past argument
        shift # past value
        ;;
        --pheno)
        pheno="$2"
        shift # past argument
        shift # past value
        ;;
        --remove)
        remove="$2"
        shift # past argument
        shift # past value
        ;;
        *)    # unknown option
        echo "Unknown option $key"
        exit 1
        ;;
    esac
done

echo "${exclude_hla}"

if [ "${exclude_hla}" == "False" ]; then
    bin/plink --bfile "${bfile}" ${pheno:+--pheno "$pheno"} ${remove:+--remove "$remove"} --allow-no-sex --pca "${pcs}" --out "${outpath}"
else
    bin/plink --bfile "${bfile}" --allow-no-sex --chr 6 --from-bp 28477797 \
        --to-bp 33448354 --make-bed --out temp --silent
    awk '{print $2}' temp.bim > exclude.txt
    rm -r temp*
    bin/plink --bfile "${bfile}" \
        ${pheno:+--pheno "$pheno"} \
        ${remove:+--remove "$remove"} \
        --allow-no-sex \
        --exclude exclude.txt \
        --pca "${pcs}" --out "${outpath}"
    rm exclude.txt
fi
