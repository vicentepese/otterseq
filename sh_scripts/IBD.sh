#!/bin/bash

# Get directories
GWASDATA=$(jq -r .plinkFiles.GWAS settings.json)
IBDGEN=$(jq -r .file.IBDGenome settings.json)
PREFIX=$(jq -r '.plinkFiles.prexi' settings.json)

king -b ${GWASDATA}${PREFIX}.bed --ibd 
mv king.seg $IBDGEN
rm -r *king*
