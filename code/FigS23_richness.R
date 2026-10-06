library(ggplot2)
library(randomForest)
library(pROC)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggpubr)

setwd("/Users/yrx/Downloads/FeatureSelection/feature_selection/code")

features <- read.csv(file = "../results/Real_feautures_filter.csv")

pdf1 = list()
index = 1

for (ffs in unique(features$ff)){
  dat1_sub <- read.csv(file = paste("../data/",ffs,sep = ""))
  ngroup <-  length(unique(features$frac[features$ff==ffs]))
  data_rich <- data.frame(species=colnames(dat1_sub),richness1=colSums(dat1_sub>0)/nrow(dat1_sub))
  if (ngroup == 10){
    pdf1[[index]] <- ggplot(data_rich, aes(richness1)) +
      geom_histogram(color="#0F7BA2FF",fill="#0F7BA2FF") +
      theme_bw()+xlab("Richness")+ylab("Frequency")+labs(title=ffs,size=2)+
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold.italic", size = 10),
        strip.background = element_blank(),
        axis.text.x = element_text(size = 8,color = 'black'),
        axis.text.y = element_text(size = 8,color = 'black'),
        #strip.text.x = element_blank(),
        axis.ticks = element_line(size = 0.2),
        legend.position = "none")
    
    index = index + 1
  }
}

p1 = do.call(ggarrange, c(pdf1, list(ncol = 5, nrow = 2),align="hv"))
print(p1)
