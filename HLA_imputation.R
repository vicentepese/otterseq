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
library(plotly)

########## IMPORT ##########
setwd("~/Documents/GWAS_pipeline")

# Import settings
settings <- jsonlite::fromJSON('settings.json')

########## HLA IMPUTATION ############

# Impute HLA using individual models if there is an ethnicity file

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
gname <- paste0(settings$plinkFiles$GWASQC, settings$plinkFiles$prefix, "_QC")
yourgeno <- hlaBED2Geno(bed.fn=paste(gname, ".bed", sep = ''), fam.fn=paste(gname, ".fam", sep='')
                        , bim.fn=paste(gname, ".bim", sep=''), assembly = 'hg19')
summary(yourgeno)

# Make cluster 
cl <- makeCluster(10)

# Make predictions
for (locus in hla.id){
  model.hla <- hlaModelFromObj(model.list[[locus]])
  summary(model.hla)
  pred.guess <- predict(model.hla, yourgeno, type="response+prob", nclassifiers=100, cl=cl, match.type="Position")
  pred.guess$value$FID <- pred.guess$value$sample.id %>% lapply(function(x) strsplit(x,"-") %>% unlist() %>% .[1]) %>% unlist()
  pred.guess$value$IID <- pred.guess$value$sample.id %>% lapply(function(x) strsplit(x,"-") %>% unlist() %>% tail(n=1)) %>% unlist()
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
probs.data <- HLA.data[,c("FID","IID", "sample.id","prob")]
HLA.data <- HLA.data[,c("FID","IID", "sample.id","allele1", "allele2")]
colnames(HLA.data) <- c("FID","IID", "sample.id", paste0(pred.guess$locus, '.1'), paste0(pred.guess$locus, '.2'))
colnames(probs.data) <- c("FID","IID", "sample.id", paste0("prob.", pred.guess$locus))

# Merge rest of files
for (f in file.names[2:length(file.names)]){
  
  # Print for loop
  print(paste0("Current file: ", f))
  
  # Load data and change colnames
  load(paste0(settings$directory$HLA_Imputation, f))
  data <- pred.guess$value[,c("FID","IID", "sample.id","allele1", "allele2")]
  colnames(data) <- c("FID","IID", "sample.id", paste0(pred.guess$locus, '.1'), paste0(pred.guess$locus, '.2'))
  
  # Create probs and change colnames
  probs <- pred.guess$value[,c("FID","IID", "sample.id","prob")]
  colnames(probs) <- c("FID","IID", "sample.id", paste0("prob.", pred.guess$locus))
  
  # Merge by sample.id
  HLA.data <- merge(HLA.data, data, by =  c("FID", "IID", "sample.id"))
  
  # Merge probs by sample.id
  probs.data <- merge(probs.data, probs, by = c("FID", "IID", "sample.id"))
}

# Write 
write.table(HLA.data, file = paste0(settings$directory$HLA_Imputation, "HLA_calls.csv"), sep = ',', quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(probs.data, file = paste0(settings$directory$HLA_Imputation, "HLA_probs.csv"), sep = ',', quote = FALSE, row.names = FALSE, col.names = TRUE)

############## PLOT PROBABILITIES #############

## Plot probabilities distributions
pl <- vector('list', ncol(probs.data)-3)
idx <- 1
for (i in 4:ncol(probs.data)){
  
  pl[[i-3]] <- local({
    i <- i
    p1 <- ggplot(probs.data, aes(get(colnames(probs.data)[i]))) +
      geom_histogram() +
      xlab(colnames(probs.data)[i]) + xlim(0,1)
    print(p1)
  })
}
do.call(grid.arrange, pl)

##  Plot PCA probabilities
# Load PCA
pca <- read.table(file = settings$file$PCA_eigenvec, sep = " ", header = FALSE)
colnames(pca) <- c("FID", "IID", paste0(rep("PC", 20), as.character(1:20)))
pca$sample.id <- paste0(pca$FID, rep('-', nrow(pca)), pca$IID)

# Merge PCA with probs 
data_plot <- merge(pca, probs.data, by = c("FID", "IID","sample.id"))

# Plot PCAs with probs
pl <- vector('list', ncol(probs.data)-3)
idx <- 1
for (i in 4:ncol(probs.data)){

  pl[[i-3]] <- local({
    i <- i
    p1 <- ggplot(data_plot, aes(PC1, PC2, color = get(colnames(probs.data)[i]))) + 
      geom_point() +  
      scale_colour_gradientn(limits = c(0,1), colors =c("navyblue", "darkmagenta", "darkorange1"), oob = scales::squish) + 
      labs(color = colnames(probs.data)[i])
    print(p1)
  })
}
do.call(grid.arrange, pl)
g <- arrangeGrob(pl) #generates g
ggsave(file="pca_probs.pdf",g)

############ COVARIATES #############

# Load PCA
pca <- read.table(file = settings$file$PCA_eigenvec, sep = " ", header = FALSE)
colnames(pca) <- c("FID", "IID", paste0(rep("PC", 20), as.character(1:20)))

# Load pheno
pheno <- read.table(file = settings$file$pheno_matched, sep = " ", header = FALSE)
colnames(pheno) <- c("FID","IID", "pheno")

# Create covars
covars <- merge(pca, pheno, by = c("FID", "IID"))
covars$sample.id <- paste0(covars$FID, rep("-", nrow(covars)), covars$IID)

# Write covariates
write.table(covars, file = paste0(settings$directory$HLA_Imputation, "covars.txt"), quote=TRUE, sep = ",", row.names = FALSE, col.names = TRUE)
