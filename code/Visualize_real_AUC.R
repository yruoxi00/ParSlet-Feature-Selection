library(ggplot2)
library(randomForest)
library(pROC)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggpubr)

setwd("/Users/yrx/Downloads/feature_selection/code")

features <- read.csv(file = "../results/Real_AUC_cross.csv")
df_features <- data.frame(AUCs=c(features$AUC_a,features$AUC_g),dataset=c(features$ff,features$ff),
                          frac=c(features$frac,features$frac),type=c(rep("Integrated",nrow(features)),rep("Original",nrow(features))))

g1 <- ggplot(data = df_features,aes(type,AUCs,color=type))+geom_boxplot()+
  facet_wrap(~frac+dataset,nrow=3,scale="free")+
  stat_compare_means(comparisons = list(c("Integrated", "Original")),paired = T,method.args = list(alternative = "greater"))+
  theme_bw()+xlab("")+ylab("AUC decrease")+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold.italic", size = 10),
    strip.background = element_blank(),
    #strip.text.x = element_blank(),
    axis.ticks = element_line(size = 0.2),
    legend.position = "none")
  
print(g1)
