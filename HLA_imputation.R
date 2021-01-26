# Import libraries
library(jsonlite)
library(tidyverse)
library(readr)
library(data.table)
library(ggrepel)
library(viridis)
library(hrbrthemes)
library(HIBAG)
library(parallel)
library(ggplot2)
library(gridExtra)
library(rlist)

########## IMPORT ##########
setwd("~/Documents/GWAS_pipeline")

# Import settings
settings <- jsonlite::fromJSON('settings.json')

########## PREFIT MODELS ############

# Load pre-fit model and comvert to hlaMODEL
model.list <- get(load(settings$file$PMRA_HLA_model))
drb3 <- get(load(settings$file$PMRA_HLA_DRB3))
drb4 <- get(load(settings$file$PMRA_HLA_DRB4))
drb5 <- get(load(settings$file$PMRA_HLA_DRB5))
model.list[["DRB3"]] <- drb3
model.list[["DRB4"]] <- drb4
model.list[["DRB5"]] <- drb5
hla.id <- names(model.list)

# Import file
gname <- settings$plinkFiles$GWASQC
yourgeno <- hlaBED2Geno(bed.fn=paste(gname, ".bed", sep = ''), fam.fn=paste(gname, ".fam", sep='')
                        , bim.fn=paste(gname, ".bim", sep=''), assembly = 'hg19')
summary(yourgeno)

# Make cluster 
cl <- makeCluster(10)

# Make predictions
for (locus in hla.id){
  model.hla <- hlaModelFromObj(model.list[[locus]])
  summary(model.hla)
  pred.guess <- predict(model.hla, yourgeno, type="response+prob", nclassifiers=100, cl=cl)
  save(pred.guess, file = paste(settings$directory$HLA_Imputation, paste('HLA_', locus, sep = ''), '.RData', sep= ''))
}

################ MERGE #############

# Initialize loop# List of files to merge 
file.names <- list.files(settings$directory$HLA_Imputation) 
file.names <- file.names[which(grepl(x=file.names, pattern = '.RData'))]

# First pass
f <- file.names[1]
load(paste0(settings$directory$HLA_Imputation, f))
HLA.data <- pred.guess$value

# Get probs and alleles
probs.data <- HLA.data[,c(1,4)]
HLA.data <- HLA.data[,1:3]
colnames(HLA.data) <- c("sample.id", paste0(pred.guess$locus, '.1'), paste0(pred.guess$locus, '.2'))
colnames(probs.data) <- c("sample.id", paste0("prob.", pred.guess$locus))

# Merge rest of files
for (f in file.names[2:length(file.names)]){
  
  # Print for loop
  print(paste0("Current file: ", f))
  
  # Load file 
  load(paste0(settings$directory$HLA_Imputation, f))
  data <- pred.guess$value[,c(1:3)]
  probs <- pred.guess$value[,c(1,4)]
  colnames(data)[1:3] <- c("sample.id", paste0(pred.guess$locus, '.1'), paste0(pred.guess$locus, '.2'))
  colnames(probs) <- c("sample.id", paste0("prob.", pred.guess$locus))
  
  # Merge by sample.id
  HLA.data <- merge(HLA.data, data, by = "sample.id")
  
  # Merge probs by sample.id
  probs.data <- merge(probs.data, probs, by = "sample.id")
}

# Write 
write.table(HLA.data, file = paste0(settings$directory$HLA_Imputation, "HLA_calls.csv"), sep = ',', quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(probs.data, file = paste0(settings$directory$HLA_Imputation, "HLA_probs.csv"), sep = ',', quote = FALSE, row.names = FALSE, col.names = TRUE)

# Plot probabilities
pl <- vector('list', ncol(probs.data)-1)
idx <- 1
for (i in 2:ncol(probs.data)){
  
  pl[[i-1]] <- local({
    i <- i
    p1 <- ggplot(probs.data, aes(get(colnames(probs.data)[i]))) +
      geom_histogram() +
      xlab(colnames(probs.data)[i]) + xlim(0,1)
    print(p1)
  })
  
  
}
do.call(grid.arrange, pl)

############ COVARIATES #############

# Load PCA
pca <- read.table(file = settings$file$PCA_eigenvec, sep = " ", header = FALSE)
colnames(pca) <- c("FID", "IID", paste0(rep("PC", 20), as.character(1:20)))

# Load pheno
pheno <- read.table(file = settings$file$pheno_matched, sep = " ", header = FALSE)
colnames(pheno) <- c("FID","IID", "pheno")

# Create covars
covars <- merge(pca, pheno, by = c("FID", "IID"))

# Write covariates
write.table(covars, file = paste0(settings$directory$HLA_Imputation, "covars.txt"), quote=TRUE, sep = ",", row.names = FALSE, col.names = TRUE)
