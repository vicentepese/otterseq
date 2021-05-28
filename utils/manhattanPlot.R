# Import libraries
library(jsonlite)
library(tidyverse)
library(readr)
library(data.table)
library(ggrepel)
library(viridis)
library(hrbrthemes)

########## IMPORT ##########
# Import settings
settings <- jsonlite::fromJSON('settings.json')

# Import association file 
assoc.data <- read.table(paste0(settings$directory$GWAS_out,settings$plinkFiles$prefix ,".assoc.logistic"), header = TRUE, sep = '')
assoc.data <- assoc.data[which(assoc.data$TEST == 'ADD'),]

######### PLOT ###########

# Data table variables
pvals <- assoc.data$P
x <- c(1:length(pvals))
chr <- assoc.data$CHR
colors <- c()
flag <- 1
for (c in unique(chr)){
  l <- nrow(assoc.data[which(assoc.data$CHR == c),])
  if (flag == 1){
    colors <- c(colors, rep('#F8766D', l))
    flag = 2
  }
  else{
    colors <- c(colors, rep('#00BFC4', l))
    flag = 1
  }
}

# Create datatable
plot.data <- data.table('pval' =  -log10(pvals), 'pos' =  x, 'colors' = colors, 'chr' = chr)
axisdf <- plot.data %>% group_by(chr) %>% summarize(center=( max(pos) + min(pos)) / 2 )

# Plot
png("GWAS.png", width = 1920, height = 1080)
ggplot() + geom_point(data = plot.data, aes(pos, pval, color = as.factor(colors)), alpha = 1) +
  scale_x_continuous(label = axisdf$chr, breaks= axisdf$center ) + 
  labs(x='Chromosome',y='-log10(pval)', title='GWAS') +
  geom_hline(yintercept = -log10(5e-8)) + 
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
dev.off()

########## MANHATTAN PLOT TEMPLATE ############3

ggManh=function(metaout){
  df.temp= metaout%>% group_by(chr) %>% 
    summarise(chr.len= max(pos)) %>% 
    mutate(tot=cumsum(as.numeric(chr.len))-chr.len) %>% 
    select(-chr.len) %>% left_join(metaout, ., by=c("chr"="chr")) %>% arrange(chr, pos) %>% mutate(BPcum=pos+tot)
  axisdf <- df.temp %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  setDT(df.temp)
  a=ggplot(df.temp[pvals > 1e-20], aes(x=BPcum, y=-log10(pval))) +
    geom_point(aes(color=as.factor(chr)), alpha=0.3, size=1.4)+geom_hline(yintercept = -log10(5e-8))+
    scale_x_continuous(label = axisdf$chr, breaks= axisdf$center )+ labs(x = "Chromosome")+theme_bw(base_size = 22) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position="none",
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) 
  return(a)
}
