# GWAS_pipeline
GWAS pipeline 

# Requirements
Required Python version
> Python v3.6.8 

Required Python modules
> numpy v1.15.4 
> subprocess
> json 
> os
> collections
> matplotlib v3.0.2
> scipy v1.2.1
> king v2.2.4

# Steps / set settings

## Folders
_Data_: Master folder containing all data utilized in the project. 
_Data/GWAS_: Folder contaning non-binarized GWAS files. The pipeline will binarize the data. Files can be added alread binarized if decided.
_Data/GWAS\_binaries_: Folder containing binarized GWAS files. Files can be directly added into the folder, or in _Data/GWAS_ to be binarized. 

## Dataset modifications
_Data/GWAS/Stanford/Plates\_117\_118\_119\_PMRA.ped_: 
1. Removed line 194 and 289, COPY
