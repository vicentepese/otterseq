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

######### IMPUTATION FUNCTION ############

imputeHLA = function(data, model.list, eth){
  
  # Get loci
  hla.id <- names(model.list)
  
  # Initialize output
  pred.guess.out <- data.frame(sample.id = data$sample.id)
  probs.out <- data.frame(sample.id = data$sample.id)
  
  # Impute for each locus
  for (locus in hla.id){
    
    # Get model
    model.hla <- hlaModelFromObj(model.list[[locus]])
    summary(model.hla)
    
    # Make predictions
    pred.guess <- predict(model.hla, data, type="response+prob", nclassifiers=100, cl=cl, match.type="Position")
    preds.all <- pred.guess$value[,c("sample.id", "allele1", "allele2")]
    colnames(preds.all) <- c("sample.id", paste0(locus, '.1'), paste0(locus, '.2'))
    
    # Get probabilities
    probs <- pred.guess$value[,c("sample.id","prob")]
    colnames(probs) <- c("sample.id", paste0("prob.", locus))
    
    # Merge 
    pred.guess.out <- pred.guess.out %>% merge(preds.all, by = "sample.id")
    probs.out <- probs.out %>% merge(probs, by = "sample.id")
    
  }
  
  # Add FID and IID
  pred.guess.out$FID <- pred.guess.out$sample.id %>% lapply(function(x) strsplit(x,"-") %>% unlist() %>% .[1]) %>% unlist()
  pred.guess.out$IID <-pred.guess.out$sample.id %>% lapply(function(x) strsplit(x,"-") %>% unlist() %>% tail(n=1)) %>% unlist()
  probs.out$FID <- probs.out$sample.id %>% lapply(function(x) strsplit(x,"-") %>% unlist() %>% .[1]) %>% unlist()
  probs.out$IID <-probs.out$sample.id %>% lapply(function(x) strsplit(x,"-") %>% unlist() %>% tail(n=1)) %>% unlist()
  
  # Write 
  save(pred.guess.out, file = paste(settings$directory$HLA_Imputation, paste('HLA_eth', eth, sep = '_'), '.RData', sep= ''))
  save(probs.out, file = paste(settings$directory$HLA_Imputation, paste('probs_eth', eth, sep = '_'), '.RData', sep= ''))
}

########## HLA IMPUTATION ############

# Import file
gname <- paste0(settings$plinkFiles$GWASQC, settings$plinkFiles$prefix, "_QC")
yourgeno <- hlaBED2Geno(bed.fn=paste(gname, ".bed", sep = ''), fam.fn=paste(gname, ".fam", sep='')
                        , bim.fn=paste(gname, ".bim", sep=''), assembly = 'hg19')
summary(yourgeno)

# Import ethnicity
eth.df <- read.table(file = settings$file$ethnicity, header = TRUE, sep = ",")
eth.df$sample.id <- paste0(eth.df$FID, rep("-", nrow(eth.df)), eth.df$IID)

# Load ethnicity models
EUR.list <- get(load(settings$file$PMRA_EUR))
AFR.list <- get(load(settings$file$PMRA_AFR))
AS.list <- get(load(settings$file$PMRA_AS))
AMR.list <- get(load(settings$file$PMRA_AMR))

# Make cluster 
cl <- makeCluster(10)

# Parse ethnicities
ethnicities <- unique(eth.df$Population)

# Make predictions for each ethnicity
for (eth in ethnicities){
  
  # Parse IDS
  eth_ids <- eth.df %>% filter(Population == eth) %>% .['sample.id'] %>% unlist()
  
  # Filter 
  data_eth.list <- yourgeno
  data_eth.list$sample.id <- data_eth.list$sample.id[which(data_eth.list$sample.id %in% eth_ids)]
  data_eth.list$genotype <- data_eth.list$genotype[,which(data_eth.list$sample.id %in% eth_ids)]
  
  # Impute for each locus 
  switch (eth,
    "EUR" = imputeHLA(data = data_eth.list, model.list = EUR.list, eth = "EUR"),
    "AFR" = imputeHLA(data = data_eth.list, model.list = AFR.list, eth = "AFR"),
    "SAS" = imputeHLA(data = data_eth.list, model.list = AS.list, eth = "SAS"),
    "EAS" = imputeHLA(data = data_eth.list, model.list = AS.list, eth = "EAS"),
    "AMR" = imputeHLA(data = data_eth.list, model.list = AMR.list, eth = "AMR")
  )
  
}


################ MERGE #############

# Initialize loop# List of files to merge 
file.names <- list.files(settings$directory$HLA_Imputation) 
file.names <- file.names[which(grepl(x=file.names, pattern = 'eth'))]
probs.file.names <- file.names[which(grepl(x = file.names, pattern = "probs"))]
HLA.file.names <- file.names[which(grepl(x = file.names, pattern = "HLA"))]

# Initialize HLA for loop
HLA.df <- data.frame()

# Merge HLA dataframes
for (HLA.file in HLA.file.names){
  
  # Print 
  print(paste0("Merging :", HLA.file))
  
  # Load file
  load(paste0(settings$directory$HLA_Imputation, HLA.file))
  
  # Bind by row
  HLA.df <- rbind(HLA.df, pred.guess.out)
  
}

# Initilize probabilities for loop
probs.df <- data.frame()
for (probs.file in probs.file.names){
  
  # Print 
  print(paste0("Merging :", probs.file))
  
  # Load file 
  load(paste0(settings$directory$HLA_Imputation, probs.file))
  
  # Bind by row 
  probs.df <- rbind(probs.df, probs.out)
  
}

# Write 
write.table(HLA.df, file = paste0(settings$directory$HLA_Imputation, "HLA_calls_eth.csv"),
            sep = ',', quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(probs.df, file = paste0(settings$directory$HLA_Imputation, "HLA_probs_eth.csv"),
            sep = ',', quote = FALSE, row.names = FALSE, col.names = TRUE)

############## PLOT PROBABILITIES #############

## Plot probabilities distributions
pl <- vector('list', ncol(probs.df)-3)
idx <- 1
for (i in 2:(ncol(probs.df)-2)){
  
  pl[[i-1]] <- local({
    i <- i
    p1 <- ggplot(probs.df, aes(get(colnames(probs.df)[i]))) +
      geom_histogram() +
      xlab(colnames(probs.df)[i]) + xlim(0,1)
    print(p1)
  })
}
do.call(grid.arrange, pl)

##  Plot PCA probabilities
# Load PCA
pca <- read.table(file = settings$file$PCA_eigenvec, sep = " ", header = FALSE)
colnames(pca) <- c("FID", "IID", paste0(rep("PC", 20), as.character(1:20)))
pca$sample.id <- paste0(pca$FID, rep('-', nrow(pca)), pca$IID)

# Merge PCA with probs and ethnicity
data_plot <- merge(pca, probs.df, by = c("sample.id"))
data_plot <-merge(data_plot, eth.df, by = "sample.id")

# Plot PCAs with probs
pl <- vector('list', ncol(probs.df)-3)
idx <- 1
for (i in 2:(ncol(probs.df)-2)){
  
  pl[[i-1]] <- local({
    i <- i
    p1 <- ggplot(data_plot, aes(PC1, PC2, color = get(colnames(probs.df)[i]))) + 
      geom_point(aes(shape = Population)) +  
      scale_colour_gradientn(limits = c(0,1), colors =c("navyblue", "darkmagenta", "darkorange1"), oob = scales::squish) + 
      labs(color = colnames(probs.df)[i])
    ggplotly(p1)
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
