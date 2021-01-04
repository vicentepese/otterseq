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
2. Modified bim file from Oxford (removed
10	chr10-96521657:rs12248560-AC	115.27	96521657	A	C
10	chr10-96521657:rs12248560-TC	115.27	96521657	A	G 
for being biallelic)

## THINGS TO CHANGE 

### PREPROCESSING 
Non-standaradized, write about it 

### 1000 genome imputation 
(i)Non-standarized, write main code
    (i)+ make path standarized
(ii) Issue when converting from bgen to binary PED - IID is replaced by index --> fix code so that sample becomes .fam
(iii) Create SLURM log folder to save slurm, an shape it
(iv) Change name of SLURM logs
(v) Still cannot merge all binary PED files
(vi) Include QC there 

### GWAS PIPELINE
(i) QC does not remove multiallelic, maybe include
(ii) PCA having issues
(iii) Current version includes no pheno -- maybe allow option instead of list?
(iv) Remove content of folders in each iteration automatically?
