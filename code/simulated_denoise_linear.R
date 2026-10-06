library(ggplot2)
library(randomForest)
library(pROC)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggpubr)

setwd("/Users/yrx/Downloads/feature_selection/code")

AUCs <- read.csv(file = "../results/Simulated0.01_AUC_cross.csv")
dat_AUC <- data.frame(AUC=c(AUCs$AUC_g,AUCs$AUC_a),frac=c(AUCs$frac,AUCs$frac),ff=c(AUCs$ff,AUCs$ff),
                      types=c(rep("Original",nrow(AUCs)),rep("Integrated",nrow(AUCs))))


pdf1<- ggplot(data = dat_AUC[dat_AUC$AUC>0,],aes(x=types,y=(AUC),color=factor(types),fill=factor(types)))+
  geom_boxplot(outlier.size = 2,outlier.shape = 21,lwd=0.75,alpha=0.2)+
  scale_fill_manual(values = c("#DD5129FF","#0F7BA2FF"))+scale_color_manual(values = c("#DD5129FF","#0F7BA2FF"))+
  facet_wrap(~ff+frac,nrow = 3,scales = "free")+
  theme_bw()+xlab("")+ylab("P-value difference")+stat_compare_means(method.args=list(alternative = "less"),label="p.format",label.y = 0.9)+
  labs(title = paste("Fraction of spiking-in: 0.05"))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 8),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text.x = element_text(size = 8,color = 'black',angle = 45,hjust = 1),
    axis.text.y = element_text(size = 8,color = 'black'),
    #strip.text.x = element_blank(),
    axis.ticks = element_line(size = 0.2),
    legend.position = "none")

pdf11<-ggplot(data = dat_AUC[dat_AUC$AUC<0,],aes(x=types,y=(AUC),color=factor(types),fill=factor(types)))+
  geom_boxplot(outlier.size = 2,outlier.shape = 21,lwd=0.75,alpha=0.2)+
  scale_fill_manual(values = c("#DD5129FF","#0F7BA2FF"))+scale_color_manual(values = c("#DD5129FF","#0F7BA2FF"))+
  facet_wrap(~ff+frac,nrow = 3,scales = "free")+
  theme_bw()+xlab("")+ylab("P-value difference")+stat_compare_means(method.args=list(alternative = "less"),label="p.format",label.y = -0.05)+
  labs(title = paste("Fraction of spiking-in: 0.05"))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 8),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text.x = element_text(size = 8,color = 'black',angle = 45,hjust = 1),
    axis.text.y = element_text(size = 8,color = 'black'),
    #strip.text.x = element_blank(),
    axis.ticks = element_line(size = 0.2),
    legend.position = "none")


AUCs <- read.csv(file = "../results/Simulated0.05_AUC_cross.csv")
dat_AUC <- data.frame(AUC=c(AUCs$AUC_a,AUCs$AUC_g),frac=c(AUCs$frac,AUCs$frac),ff=c(AUCs$ff,AUCs$ff),
                      types=c(rep("Original",nrow(AUCs)),rep("Integrated",nrow(AUCs))))

pdf2<- ggplot(data = dat_AUC[dat_AUC$AUC>0,],aes(x=types,y=(AUC),color=factor(types),fill=factor(types)))+
  geom_boxplot(outlier.size = 2,outlier.shape = 21,lwd=0.75,alpha=0.2)+
  scale_fill_manual(values = c("#DD5129FF","#0F7BA2FF"))+scale_color_manual(values = c("#DD5129FF","#0F7BA2FF"))+
  facet_wrap(~ff+frac,nrow = 3,scales = "free")+
  theme_bw()+xlab("")+ylab("P-value difference")+stat_compare_means(method.args=list(alternative = "less"),label="p.format",label.y = 0.9)+
  labs(title = paste("Fraction of spiking-in: 0.1"))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 8),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text.x = element_text(size = 8,color = 'black',angle = 45,hjust = 1),
    axis.text.y = element_text(size = 8,color = 'black'),
    #strip.text.x = element_blank(),
    axis.ticks = element_line(size = 0.2),
    legend.position = "none")

pdf21<- ggplot(data = dat_AUC[dat_AUC$AUC<0,],aes(x=types,y=(AUC),color=factor(types),fill=factor(types)))+
  geom_boxplot(outlier.size = 2,outlier.shape = 21,lwd=0.75,alpha=0.2)+
  scale_fill_manual(values = c("#DD5129FF","#0F7BA2FF"))+scale_color_manual(values = c("#DD5129FF","#0F7BA2FF"))+
  facet_wrap(~ff+frac,nrow = 3,scales = "free")+
  theme_bw()+xlab("")+ylab("P-value difference")+stat_compare_means(method.args=list(alternative = "less"),label="p.format",label.y = -0.05)+
  labs(title = paste("Fraction of spiking-in: 0.1"))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 8),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text.x = element_text(size = 8,color = 'black',angle = 45,hjust = 1),
    axis.text.y = element_text(size = 8,color = 'black'),
    #strip.text.x = element_blank(),
    axis.ticks = element_line(size = 0.2),
    legend.position = "none")

AUCs <- read.csv(file = "../results/Simulated0.2_AUC_cross.csv")
dat_AUC <- data.frame(AUC=c(AUCs$AUC_a,AUCs$AUC_g),frac=c(AUCs$frac,AUCs$frac),ff=c(AUCs$ff,AUCs$ff),
                      types=c(rep("Original",nrow(AUCs)),rep("Integrated",nrow(AUCs))))

pdf3<- ggplot(data = dat_AUC[dat_AUC$AUC>0,],aes(x=types,y=(AUC),color=factor(types),fill=factor(types)))+
  geom_boxplot(outlier.size = 2,outlier.shape = 21,lwd=0.75,alpha=0.2)+
  scale_fill_manual(values = c("#DD5129FF","#0F7BA2FF"))+scale_color_manual(values = c("#DD5129FF","#0F7BA2FF"))+
  facet_wrap(~ff+frac,nrow = 3,scales = "free")+
  theme_bw()+xlab("")+ylab("P-value difference")+stat_compare_means(method.args=list(alternative = "less"),label="p.format",label.y = 0.9)+
  labs(title = paste("Fraction of spiking-in: 0.2"))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 8),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text.x = element_text(size = 8,color = 'black',angle = 45,hjust = 1),
    axis.text.y = element_text(size = 8,color = 'black'),
    #strip.text.x = element_blank(),
    axis.ticks = element_line(size = 0.2),
    legend.position = "none")
  
pdf31<- ggplot(data = dat_AUC[dat_AUC$AUC<0,],aes(x=types,y=(AUC),color=factor(types),fill=factor(types)))+
  geom_boxplot(outlier.size = 2,outlier.shape = 21,lwd=0.75,alpha=0.2)+
  scale_fill_manual(values = c("#DD5129FF","#0F7BA2FF"))+scale_color_manual(values = c("#DD5129FF","#0F7BA2FF"))+
  facet_wrap(~ff+frac,nrow = 3,scales = "free")+
  theme_bw()+xlab("")+ylab("P-value difference")+stat_compare_means(method.args=list(alternative = "less"),label="p.format",label.y = -0.05)+
  labs(title = paste("Fraction of spiking-in: 0.2"))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 8),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text.x = element_text(size = 8,color = 'black',angle = 45,hjust = 1),
    axis.text.y = element_text(size = 8,color = 'black'),
    #strip.text.x = element_blank(),
    axis.ticks = element_line(size = 0.2),
    legend.position = "none")

p1 = ggarrange(pdf1,pdf2,pdf3,pdf11,pdf21,pdf31,ncol = 6,nrow=1)
print(p1)

ggsave(p1,file=paste("../figures/denoise_linear.pdf",sep = ""),width=8.5, height=7.5)

