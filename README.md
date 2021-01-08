# GWAS_pipeline
GWAS pipeline 

# Requirements
Required Python version
> Python v3.6.8 

Required Python modules
> numpy v1.15.4 <br>
> subprocess <br>
> json <br>
> os <br>
> collections <br>
> matplotlib v3.0.2 <br>
> scipy v1.2.1 <br>
> king v2.2.4 <br>

# Steps / set settings

## Folders
_Data_: Master folder containing all data utilized in the project. 
_Data/GWAS_: Folder contaning non-binarized GWAS files. The pipeline will binarize the data. Files can be added alread binarized if decided.
_Data/GWAS\_binaries_: Folder containing binarized GWAS files. Files can be directly added into the folder, or in _Data/GWAS_ to be binarized. 

## Dataset modifications
1. Removed line 194 and 289, COPY
2. Modified bim file from Oxford (removed
10	chr10-96521657:rs12248560-AC	115.27	96521657	A	C
10	chr10-96521657:rs12248560-TC	115.27	96521657	A	G 
for being biallelic)

## THINGS TO CHANGE 

### PREPROCESSING 
Non-standaradized, write about it 

### 1000 genome imputation 
- [X] Non-standarized, write main code (V)
  - [X]   (i)+ make path standarized (V)
- [X] Issue when converting from bgen to binary PED - IID is replaced by index --> fix code so that sample becomes .fam
- [X] Create SLURM log folder to save slurm, an shape it
- [X] Change name of SLURM logs
- [X] Still cannot merge all binary PED files
- [ ] Include QC there 

### GWAS PIPELINE
- [X] QC does not remove multiallelic, maybe include
- [X] PCA having issues
- [X] Current version includes no pheno -- maybe allow option instead of list?
- [ ] Remove content of folders in each iteration automatically?
- [X] PCA cannot be performed if "," are in the SNPS ids 
