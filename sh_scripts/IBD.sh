#!/bin/bash

# Get directories
GWASDATA=$(jq -r .plinkFikes.GWAS)
IBDGEN=$(jq -r .files.IBDGenome)

king -b $GWASDATA.bed --ibd
mv king.seg $IBDGEN
rm -r *king*
