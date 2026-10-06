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
pdf2 = list()
index = 1

for (ffs in unique(features$ff)){
  ngroup <-  length(unique(features$frac[features$ff==ffs]))
  print(ngroup)
  dat1_sub <- read.csv(file = paste("../data/",ffs,sep = ""))
  
  if (ngroup==10){
    overlap_matrix1 <- matrix(0, nrow = ngroup,ncol=ngroup)
    overlap_matrix2 <- matrix(0,nrow = ngroup,ncol=ngroup)
    overlap_matrix3 <- matrix(0, nrow = ngroup,ncol=ngroup)
    overlap_matrix4 <- matrix(0,nrow = ngroup,ncol=ngroup)
    
    for (rea in 1:10){
      for (i in 1:ngroup){
        for (j in 1:ngroup){
          features_sub <- features[features$rea==rea & features$ff==ffs,]
          i1 <- features_sub$frac == unique(features_sub$frac)[i]
          i2 <- features_sub$frac == unique(features_sub$frac)[j]
          
          g1 <-  features_sub$MeanDecreaseGini[i1]
          g2 <-  features_sub$MeanDecreaseGini[i2]
          
          a1 <-  features_sub$MeanDecreaseAccuracy[i1]
          a2 <-  features_sub$MeanDecreaseAccuracy[i2]
          
          sorted_indicesg1 <- order(g1,decreasing = T)
          sorted_indicesg2 <- order(g2,decreasing = T)
          
          sorted_indicesa1 <- order(a1,decreasing = T)
          sorted_indicesa2 <- order(a2,decreasing = T)
          
          set.seed(rea)
          train_table1 <- dat1_sub[sample(1:nrow(dat1_sub), size = unique(features_sub$frac)[i] * ncol(dat1_sub)),]
          train_table2 <- dat1_sub[sample(1:nrow(dat1_sub), size = unique(features_sub$frac)[j] * ncol(dat1_sub)),]
          
          data_rich <- data.frame(species=colnames(train_table1),richness1=colSums(train_table1>0)/nrow(train_table1),
                                  richness2=colSums(train_table2>0)/nrow(train_table2))
          
          
          g3 <- features_sub$MeanDecreaseGini[i1] *  as.numeric(data_rich$richness1[match(features_sub$Feature[i1],data_rich$species)])^2.5
          g4 <- features_sub$MeanDecreaseGini[i2] *  as.numeric(data_rich$richness2[match(features_sub$Feature[i2],data_rich$species)])^2.5
          
          a3 <- features_sub$MeanDecreaseAccuracy[i1] *  (data_rich$richness1[match(features_sub$Feature[i1],data_rich$species)])^2.5
          a4 <- features_sub$MeanDecreaseAccuracy[i2] *  (data_rich$richness2[match(features_sub$Feature[i2],data_rich$species)])^2.5
          
          sorted_indicesg3 <- order(g3,decreasing = T)
          sorted_indicesg4 <- order(g4,decreasing = T)
          
          sorted_indicesa3 <- order(a3,decreasing = T)
          sorted_indicesa4 <- order(a4,decreasing = T)
          
          top_n = floor(0.1*length(g1))
          
          overlap_matrix1[i,j] <- overlap_matrix1[i,j] + length(intersect(sorted_indicesg1[1:top_n],sorted_indicesg2[1:top_n]))/top_n
          overlap_matrix2[i,j] <- overlap_matrix2[i,j] + length(intersect(sorted_indicesa1[1:top_n],sorted_indicesa2[1:top_n]))/top_n
          overlap_matrix3[i,j] <- overlap_matrix3[i,j] + length(intersect(sorted_indicesg3[1:top_n],sorted_indicesg4[1:top_n]))/top_n
          overlap_matrix4[i,j] <- overlap_matrix4[i,j] + length(intersect(sorted_indicesa3[1:top_n],sorted_indicesa4[1:top_n]))/top_n
        }
      }
    }
    overlap_matrix1 <- overlap_matrix1 / 10
    overlap_matrix2 <- overlap_matrix2 / 10
    overlap_matrix3 <- overlap_matrix3 / 10
    overlap_matrix4 <- overlap_matrix4 / 10
    
    diag(overlap_matrix1) <- NA
    diag(overlap_matrix2) <- NA
    diag(overlap_matrix3) <- NA
    diag(overlap_matrix4) <- NA
    
    rownames(overlap_matrix1) <- unique(features_sub$frac)
    colnames(overlap_matrix1) <- unique(features_sub$frac)
    
    rownames(overlap_matrix2) <- unique(features_sub$frac)
    colnames(overlap_matrix2) <- unique(features_sub$frac)
    
    rownames(overlap_matrix3) <- unique(features_sub$frac)
    colnames(overlap_matrix3) <- unique(features_sub$frac)
    
    rownames(overlap_matrix4) <- unique(features_sub$frac)
    colnames(overlap_matrix4) <- unique(features_sub$frac)
    
    overlap_means1 <- apply(overlap_matrix1, 1, function(x) mean(x, na.rm = TRUE))
    overlap_means2 <- apply(overlap_matrix2, 1, function(x) mean(x, na.rm = TRUE))
    overlap_means3 <- apply(overlap_matrix3, 1, function(x) mean(x, na.rm = TRUE))
    overlap_means4 <- apply(overlap_matrix4, 1, function(x) mean(x, na.rm = TRUE))
    
    overlap_sd1 <- apply(overlap_matrix1, 1, function(x) sd(x, na.rm = TRUE))
    overlap_sd2 <- apply(overlap_matrix2, 1, function(x) sd(x, na.rm = TRUE))
    overlap_sd3 <- apply(overlap_matrix3, 1, function(x) sd(x, na.rm = TRUE))
    overlap_sd4 <- apply(overlap_matrix4, 1, function(x) sd(x, na.rm = TRUE))
    
    dat1 <- data.frame(sample_fraction = rep(rownames(overlap_matrix1),2),Overlap = c(overlap_means1,overlap_means2,overlap_means3,overlap_means4),
                       sds = c(overlap_sd1/sqrt(9),overlap_sd2/sqrt(9),overlap_sd3/sqrt(9),overlap_sd4/sqrt(9)),
                       types=c(rep("Gini_o",10),rep("MeanAccuracy_o",10),rep("Gini_i",10),rep("MeanAccuracy_i",10)))
    
    dat1$sample_fraction <- as.numeric(dat1$sample_fraction)
    
    # pdf1[[index]] <- ggplot(dat1, aes(sample_fraction, Overlap)) +
    #   geom_line(aes(group=types,color=types),linetype = "dotdash",size=0.4) +geom_point(aes(color=types,shape=types),size=2.5) +
    #   geom_errorbar(aes(ymin = Overlap - sds, ymax = Overlap + sds,color=types), width = 0.02) +
    #   scale_color_manual(values = c("#DD5129FF","#DD5129FF","#0F7BA2FF", "#0F7BA2FF"))+scale_shape_manual(values = c(16,21,16,21))+
    #   theme_bw()+xlab("")+ylab("")+
    #   theme(
    #     plot.title = element_text(hjust = 0.5, size = 8),
    #     panel.grid.major = element_blank(),
    #     panel.grid.minor = element_blank(),
    #     strip.text = element_text(face = "bold.italic", size = 10),
    #     strip.background = element_blank(),
    #     axis.text.x = element_text(size = 8,color = 'black'),
    #     axis.text.y = element_text(size = 8,color = 'black'),
    #     #strip.text.x = element_blank(),
    #     axis.ticks = element_line(size = 0.2),
    #     legend.position = "none")
    
    pdf1[[index]] <- ggplot(dat1, aes(sample_fraction, Overlap)) +
      geom_line(aes(group=types,color=types), linetype = "dotdash", size = 0.4) +
      geom_point(aes(color=types, shape=types), size = 2.5) +
      geom_errorbar(aes(ymin = Overlap - sds, ymax = Overlap + sds, color=types), width = 0.02) +
      scale_color_manual(values = c("#DD5129FF","#DD5129FF","#0F7BA2FF", "#0F7BA2FF")) +
      scale_shape_manual(values = c(16,21,16,21)) +
      ggtitle(gsub("_", " ", gsub("\\.csv", "", ffs))) +  # use ffs (filename) as the title
      theme_bw() +
      xlab("Sample size") + ylab("Overlap") +
      theme(
        plot.title = element_text(hjust = 0.5, size = 9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 8, color = 'black'),
        axis.text.y = element_text(size = 8, color = 'black'),
        axis.ticks = element_line(size = 0.2),
        legend.position = "none"
      )
    
    index = index + 1
  }
}

p1 = do.call(ggarrange, c(pdf1, list(ncol = 5, nrow = 2),align="hv"))
print(p1)

legend_plot <- ggplot(data.frame(x=1, y=1, types=factor(c("Gini_i", "Gini_o", "MeanAccuracy_i", "MeanAccuracy_o"),
                                                        levels=c("Gini_i", "Gini_o", "MeanAccuracy_i", "MeanAccuracy_o"))),
                      aes(x=x, y=y, color=types, shape=types)) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c("Gini_i"="#DD5129FF", "Gini_o"="#DD5129FF",
               "MeanAccuracy_i"="#0F7BA2FF", "MeanAccuracy_o"="#0F7BA2FF"),
    labels = c("Gini impurity (Integrated)", "Gini impurity (Original)",
               "Mean accuracy decrease (Integrated)", "Mean accuracy decrease (Original)")
  ) +
  scale_shape_manual(
    values = c("Gini_i"=16, "Gini_o"=21,
               "MeanAccuracy_i"=16, "MeanAccuracy_o"=21),
    labels = c("Gini impurity (Integrated)", "Gini impurity (Original)",
               "Mean accuracy decrease (Integrated)", "Mean accuracy decrease (Original)")
  ) +
  guides(
    color = guide_legend(ncol = 2, override.aes = list(size = 4)),
    shape = guide_legend(ncol = 2)
  ) +
  theme_void() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.box = "horizontal"
  )

legend <- get_legend(legend_plot)

p1_all <- ggarrange(
  do.call(ggarrange, c(pdf1, list(ncol = 5, 
                                  nrow = 2, 
                                  align = "hv",
                                  labels = letters[1:10],   
                                  label.x = 0,             
                                  label.y = 1))),
  legend,
  ncol = 1,
  heights = c(10, 1)
)
print(p1_all)

ggsave(p1_all,file=paste("../figs/Real_feature_power_filter2.png",sep = ""),width=12, height=5.5, dpi = 500)
