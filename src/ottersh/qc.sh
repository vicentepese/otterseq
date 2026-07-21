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
        --maf)
        maf="$2"
        shift # past argument
        shift # past value
        ;;
        --geno-miss)
        geno_miss="$2"
        shift # past argument
        shift # past value
        ;;
        --indv-miss)
        indv_miss="$2"
        shift # past argument
        shift # past value
        ;;
        --rm-vars)
        rm_vars="$2"
        shift # past argument
        shift # past value
        ;;
        --rm-indv)
        rm_indv="$2"
        shift # past argument
        shift # past value
        ;;
        --pheno)
        pheno="$2"
        shift # past argument
        shift # past value
        ;;
        *)    # unknown option
        echo "Unknown option $key"
        exit 1
        ;;
    esac
done

# Extract basename from input file
basename=$(basename "$bfile")
outfile="${outpath}/${basename}"

# shellcheck disable=SC2046
bin/plink --bfile "${bfile}" \
      --pheno "${pheno}" \
      $([ -n "${geno_miss}" ] && echo "--geno ${geno_miss}") \
      $([ -n "${indv_miss}" ] && echo "--mind ${indv_miss}") \
      $([ -n "${maf}" ] && echo "--maf ${maf}") \
      $([ -n "${rm_indv}" ] && echo "--remove ${rm_indv}") \
      $([ -n "${rm_vars}" ] && echo "--exclude ${rm_vars}") \
      --make-bed \
      --out "${outfile}"