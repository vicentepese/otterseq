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
        --threshold)
        threshold="$2"
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

# Compute IBD
bin/plink2 --bfile "$bfile" --pheno "$pheno" ${remove:+--remove "$remove"} --make-king-table --out "$bfile"
bin/plink2 --bfile "$bfile" --king-cutoff-table "$bfile".kin0 "$threshold" --out "$bfile"
