#!/bin/bash

# Get directories
GWASDATA=$(jq -r .plinkFiles.GWAS settings.json)
IBDGEN=$(jq -r .file.IBDGenome settings.json)

king -b $GWASDATA.bed --ibd
mv king.seg $IBDGEN
rm -r *king*
