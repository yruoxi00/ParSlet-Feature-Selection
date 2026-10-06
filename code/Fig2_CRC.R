library(ggplot2)
library(randomForest)
library(pROC)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggpubr)
library(ggtext)
library(cowplot)

setwd("/Users/yrx/Downloads/FeatureSelection/feature_selection/code")

features <- read.csv(file = "../results/Real_feautures_filter.csv")

pdf1 = list()
pdf2 = list()
pdf3 = list()
pdf4 = list()

index = 1

for (ffs in unique(features$ff)[6:7]){
  ngroup <-  length(unique(features$frac[features$ff==ffs]))
  print(ngroup)
  
  taxa_table <- read.csv(file = paste("../data/",ffs,sep = ""))
  set.seed(1)
  train_indices1 <- sample(1:nrow(taxa_table), size = 0.1 * (ncol(taxa_table)-1))
  
  set.seed(8)
  train_indices2 <- sample(1:nrow(taxa_table), size = 0.33 * (ncol(taxa_table)-1))
  
  features_sub <- features[features$rea==1 & features$ff==ffs,]
  i1 <- features_sub$frac == unique(features_sub$frac)[1]
  i2 <- features_sub$frac == unique(features_sub$frac)[8]
  
  feature_sub1 <- features_sub[i1,]
  feature_sub2 <- features_sub[i2,]
  
  g1 <-  feature_sub1$MeanDecreaseGini
  g2 <-  feature_sub2$MeanDecreaseGini
  
  a1 <-  feature_sub1$MeanDecreaseAccuracy
  a2 <-  feature_sub2$MeanDecreaseAccuracy
  
  sorted_indicesg1 <- order(g1,decreasing = T)
  sorted_indicesg2 <- order(g2,decreasing = T)
  
  sorted_indicesa1 <- order(a1,decreasing = T)
  sorted_indicesa2 <- order(a2,decreasing = T)
  
  top_n = floor(0.1*length(g1))
  
  train_table1 <- taxa_table[train_indices1,]
  train_table2 <- taxa_table[train_indices2,]
  
  #data_rich1 <- data.frame(species = colnames(train_table1),richness = colSums(train_table1[train_table1$status=="H",]>0)-colSums(train_table1[train_table1$status=="CRC",]>0))
  #data_rich2 <- data.frame(species = colnames(train_table2),richness = colSums(train_table2[train_table2$status=="H",]>0)-colSums(train_table2[train_table2$status=="CRC",]>0))  
  
  data_rich1 <- data.frame(species = colnames(train_table1),richness = colSums(train_table1>0))
  data_rich2 <- data.frame(species = colnames(train_table2),richness = colSums(train_table2>0))  
  
  
  data_rich1$richness <- data_rich1$richness/nrow(train_table1)
  data_rich2$richness <- data_rich2$richness/nrow(train_table2)
  
  
  dat1 <- data.frame(features1 = feature_sub1$Feature[sorted_indicesg1[1:top_n]],rank=1:top_n,
                     richness1=data_rich1$richness[match(feature_sub1$Feature[sorted_indicesg1[1:top_n]],data_rich1$species)],
                     importance1 = feature_sub1$MeanDecreaseGini[sorted_indicesg1[1:top_n]],
                     
                     features2 = feature_sub2$Feature[sorted_indicesg2[1:top_n]],rank=1:top_n,
                     richness2=data_rich2$richness[match(feature_sub2$Feature[sorted_indicesg2[1:top_n]],data_rich2$species)],
                     importance2 = feature_sub2$MeanDecreaseGini[sorted_indicesg2[1:top_n]],
                     
                     features3 = feature_sub1$Feature[sorted_indicesa1[1:top_n]],rank=1:top_n,
                     richness3=data_rich1$richness[match(feature_sub1$Feature[sorted_indicesa1[1:top_n]],data_rich1$species)],
                     importance3 = feature_sub1$MeanDecreaseAccuracy[sorted_indicesa1[1:top_n]],
                     
                     features4 = feature_sub2$Feature[sorted_indicesa2[1:top_n]],rank=1:top_n,
                     richness4=data_rich2$richness[match(feature_sub2$Feature[sorted_indicesa2[1:top_n]],data_rich2$species)],
                     importance4 = feature_sub2$MeanDecreaseAccuracy[sorted_indicesa2[1:top_n]]
  )
  
  extract_genus <- function(x) {
    result <- "unclassified"
    
    genus <- stringr::str_match(x, "g__([^.; ]*)")[,2]
    if (!is.na(genus) && genus != "" && genus != "." && genus != "_") return(genus)
    
    family <- stringr::str_match(x, "f__([^.; ]*)")[,2]
    if (!is.na(family) && family != "" && family != "." && family != "_") return(family)
    
    order <- stringr::str_match(x, "o__([^.; ]*)")[,2]
    if (!is.na(order) && order != "" && order != "." && order != "_") return(order)
    
    class <- stringr::str_match(x, "c__([^.; ]*)")[,2]
    if (!is.na(class) && class != "" && class != "." && class != "_") return(class)
    
    phylum <- stringr::str_match(x, "p__([^.; ]*)")[,2]
    if (!is.na(phylum) && phylum != "" && phylum != "." && phylum != "_") return(phylum)
    
    kingdom <- stringr::str_match(x, "k__([^.; ]*)")[,2]
    if (!is.na(kingdom) && kingdom != "" && kingdom != "." && kingdom != "_") return(kingdom)
    
    return(result)
  }
  
  dat1$genus1 <- vapply(dat1$features1, extract_genus, "")
  dat1$genus2 <- vapply(dat1$features2, extract_genus, "")
  dat1$genus3 <- vapply(dat1$features3, extract_genus, "")
  dat1$genus4 <- vapply(dat1$features4, extract_genus, "")
  
  highlight_taxa <- c("Clostridium_XVIII","Ruminococcus", "Peptostreptococcus")
  dat1$genus1_lab <- ifelse(dat1$genus1 %in% highlight_taxa,
                            paste0("<span style='color:#33a02c'>", dat1$genus1, "</span>"),
                            dat1$genus1)
  dat1$genus2_lab <- ifelse(dat1$genus2 %in% highlight_taxa,
                            paste0("<span style='color:#33a02c'>", dat1$genus2, "</span>"),
                            dat1$genus2)
  dat1$genus3_lab <- ifelse(dat1$genus3 %in% highlight_taxa,
                            paste0("<span style='color:#33a02c'>", dat1$genus3, "</span>"),
                            dat1$genus3)
  dat1$genus4_lab <- ifelse(dat1$genus4 %in% highlight_taxa,
                            paste0("<span style='color:#33a02c'>", dat1$genus4, "</span>"),
                            dat1$genus4)
  
  # last_levels1 <- sapply(dat1$features1, function(taxonomy) {
  #   last_level <- unlist(strsplit(taxonomy, "__"))[length(unlist(strsplit(taxonomy, "__")))]
  # })
  # dat1$genus1 <- as.character(last_levels1) 
  # 
  # last_levels2 <- sapply(dat1$features2, function(taxonomy) {
  #   last_level <- unlist(strsplit(taxonomy, "__"))[length(unlist(strsplit(taxonomy, "__")))]
  # })
  # dat1$genus2 <- as.character(last_levels2) 
  # 
  # last_levels3 <- sapply(dat1$features3, function(taxonomy) {
  #   last_level <- unlist(strsplit(taxonomy, "__"))[length(unlist(strsplit(taxonomy, "__")))]
  # })
  # dat1$genus3 <- as.character(last_levels3) 
  # 
  # last_levels4 <- sapply(dat1$features4, function(taxonomy) {
  #   last_level <- unlist(strsplit(taxonomy, "__"))[length(unlist(strsplit(taxonomy, "__")))]
  # })
  # dat1$genus4 <- as.character(last_levels4) 
  
  
  pdf1[[index]] <- ggplot(dat1, aes(importance1, reorder(genus1_lab,importance1))) +
    geom_segment(aes(yend = genus1_lab, xend = 0), color = "gray70", size = 0.7) +
    geom_point(aes(size = abs(richness1)), color="#fb8072") +
    scale_size_continuous(range = c(0.5, 4)) +
    theme_bw() + xlab("Gini impurity") + ylab("") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_markdown(size=10),   
      axis.title.x = element_text(size = 8, face = "bold"),
      legend.position = "none"
    )
  
  pdf2[[index]] <- ggplot(dat1, aes(importance2, reorder(genus2_lab,importance2))) +
    geom_segment(aes(yend = genus2_lab, xend = 0), color = "gray70", size = 0.7) +
    geom_point(aes(size = abs(richness2)), color="#fb8072") +
    scale_size_continuous(range = c(0.5, 4)) +
    theme_bw() + xlab("Gini impurity") + ylab("") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_markdown(size=10),
      axis.title.x = element_text(size = 8, face = "bold"),
      legend.position = "none"
    )
  
  pdf3[[index]] <- ggplot(dat1, aes(importance3, reorder(genus3_lab,importance3))) +
    geom_segment(aes(yend = genus3_lab, xend = 0), color = "gray70", size = 0.7) +
    geom_point(aes(size = abs(richness3)), color="#fb8072") +
    scale_size_continuous(range = c(0.5, 4)) +
    theme_bw() + xlab("Mean decrease accuracy") + ylab("") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_markdown(size=10),
      axis.title.x = element_text(size = 8, face = "bold"),
      legend.position = "none"
    )
  
  pdf4[[index]] <- ggplot(dat1, aes(importance4, reorder(genus4_lab,importance4))) +
    geom_segment(aes(yend = genus4_lab, xend = 0), color = "gray70", size = 0.7) +
    geom_point(aes(size = abs(richness4)), color="#fb8072") +
    scale_size_continuous(range = c(0.5, 4)) +
    theme_bw() + xlab("Mean decrease accuracy") + ylab("") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_markdown(size=10),
      axis.title.x = element_text(size = 8, face = "bold"),
      legend.position = "none"
    )
  
  index = index + 1
}

plots1 <- c(pdf1, pdf2, pdf3, pdf4)
aligned1 <- align_plots(plotlist = plots1, align = "hv")

p1 <- plot_grid(plotlist = aligned1, ncol = 4, nrow = 2,
                   labels = c("a","b","c","d","e","f","g","h"))
print(p1)

header_row <- ggarrange(
  as_ggplot(text_grob("CRC_Baxter (10% training)", face = "bold", size = 10)),
  as_ggplot(text_grob("CRC_Zeller (10% training)", face = "bold", size = 10)),
  as_ggplot(text_grob("CRC_Baxter (33% training)", face = "bold", size = 10)),
  as_ggplot(text_grob("CRC_Zeller (33% training)", face = "bold", size = 10)),
  ncol = 4
)

final_plot <- ggarrange(
  header_row,
  p1,   
  ncol = 1, nrow = 2,
  heights = c(0.1, 1, 1)
)
print(final_plot)

#############################################################################

pdf1 = list()
pdf2 = list()
pdf3 = list()
pdf4 = list()

index = 1

for (ffs in unique(features$ff)[6:7]){
  ngroup <- length(unique(features$frac[features$ff == ffs]))
  print(ngroup)
  
  taxa_table <- read.csv(file = paste("../data/", ffs, sep = ""))
  set.seed(1)
  train_indices1 <- sample(1:nrow(taxa_table), size = 0.1 * (ncol(taxa_table) - 1))
  
  set.seed(8)
  train_indices2 <- sample(1:nrow(taxa_table), size = 0.33 * (ncol(taxa_table) - 1))
  
  features_sub <- features[features$rea == 1 & features$ff == ffs, ]
  i1 <- features_sub$frac == unique(features_sub$frac)[1]
  i2 <- features_sub$frac == unique(features_sub$frac)[8]
  
  train_table1 <- taxa_table[train_indices1, ]
  train_table2 <- taxa_table[train_indices2, ]
  
  data_rich1 <- data.frame(species = colnames(train_table1),
                           richness = colSums(train_table1 > 0) / nrow(train_table1))
  data_rich2 <- data.frame(species = colnames(train_table2),
                           richness = colSums(train_table2 > 0) / nrow(train_table2))
  
  feature_sub1 <- features_sub[i1, ]
  feature_sub2 <- features_sub[i2, ]
  
  feature_sub1$richness <- data_rich1$richness[match(feature_sub1$Feature, data_rich1$species)]
  feature_sub2$richness <- data_rich2$richness[match(feature_sub2$Feature, data_rich2$species)]
  
  g1 <- feature_sub1$MeanDecreaseGini * (feature_sub1$richness ^ 2.5)
  g2 <- feature_sub2$MeanDecreaseGini * (feature_sub2$richness ^ 2.5)
  a1 <- feature_sub1$MeanDecreaseAccuracy * (feature_sub1$richness ^ 2.5)
  a2 <- feature_sub2$MeanDecreaseAccuracy * (feature_sub2$richness ^ 2.5)
  
  sorted_indicesg1 <- order(g1, decreasing = TRUE)
  sorted_indicesg2 <- order(g2, decreasing = TRUE)
  sorted_indicesa1 <- order(a1, decreasing = TRUE)
  sorted_indicesa2 <- order(a2, decreasing = TRUE)
  
  top_n <- floor(0.1 * length(g1))
  
  dat1 <- data.frame(
    features1 = feature_sub1$Feature[sorted_indicesg1[1:top_n]], rank = 1:top_n,
    richness1 = feature_sub1$richness[sorted_indicesg1[1:top_n]],
    importance1 = g1[sorted_indicesg1[1:top_n]],
    
    features2 = feature_sub2$Feature[sorted_indicesg2[1:top_n]], rank = 1:top_n,
    richness2 = feature_sub2$richness[sorted_indicesg2[1:top_n]],
    importance2 = g2[sorted_indicesg2[1:top_n]],
    
    features3 = feature_sub1$Feature[sorted_indicesa1[1:top_n]], rank = 1:top_n,
    richness3 = feature_sub1$richness[sorted_indicesa1[1:top_n]],
    importance3 = a1[sorted_indicesa1[1:top_n]],
    
    features4 = feature_sub2$Feature[sorted_indicesa2[1:top_n]], rank = 1:top_n,
    richness4 = feature_sub2$richness[sorted_indicesa2[1:top_n]],
    importance4 = a2[sorted_indicesa2[1:top_n]]
  )
  
  last_levels1 <- sapply(dat1$features1, function(taxonomy) {
    last_level <- unlist(strsplit(taxonomy, "__"))[length(unlist(strsplit(taxonomy, "__")))]
  })
  dat1$genus1 <- as.character(last_levels1) 
  
  last_levels2 <- sapply(dat1$features2, function(taxonomy) {
    last_level <- unlist(strsplit(taxonomy, "__"))[length(unlist(strsplit(taxonomy, "__")))]
  })
  dat1$genus2 <- as.character(last_levels2) 
  
  last_levels3 <- sapply(dat1$features3, function(taxonomy) {
    last_level <- unlist(strsplit(taxonomy, "__"))[length(unlist(strsplit(taxonomy, "__")))]
  })
  dat1$genus3 <- as.character(last_levels3) 
  
  last_levels4 <- sapply(dat1$features4, function(taxonomy) {
    last_level <- unlist(strsplit(taxonomy, "__"))[length(unlist(strsplit(taxonomy, "__")))]
  })
  dat1$genus4 <- as.character(last_levels4) 
  
  
  pdf1[[index]] <-  ggplot(dat1, aes(importance1, reorder(genus1,importance1))) +
    geom_segment(aes(yend = genus1, xend = 0), color = "gray70", size = 0.7) +
    geom_point(aes(size = abs(richness1)),color="#fb8072")+
    scale_size_continuous(range = c(0.5, 4))+
    theme_bw()+xlab("Gini impurity")+ylab("")+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold.italic", size = 10),
      strip.background = element_blank(),
      #strip.text.x = element_blank(),
      axis.ticks = element_line(size = 0.2),
      legend.position = "none", 
      axis.title.x = element_text(size = 6, face = "bold"))
  
  
  pdf2[[index]] <- ggplot(dat1, aes(importance2, reorder(genus2,importance2))) +
    geom_segment(aes(yend = genus2, xend = 0), color = "gray70", size = 0.7) +
    geom_point(aes(size = abs(richness2)),color="#fb8072")+
    scale_size_continuous(range = c(0.5, 4))+
    theme_bw()+xlab("Gini impurity")+ylab("")+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold.italic", size = 10),
      strip.background = element_blank(),
      #strip.text.x = element_blank(),
      axis.ticks = element_line(size = 0.2),
      legend.position = "none", 
      axis.title.x = element_text(size = 6, face = "bold"))
  
  pdf3[[index]] <-  ggplot(dat1, aes(importance3, reorder(genus3,importance3))) +
    geom_segment(aes(yend = genus3, xend = 0), color = "gray70", size = 0.7) +
    geom_point(aes(size = abs(richness3)),color="#fb8072")+
    scale_size_continuous(range = c(0.5, 4))+
    theme_bw()+xlab("Mean decrease accuracy")+ylab("")+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold.italic", size = 10),
      strip.background = element_blank(),
      #strip.text.x = element_blank(),
      axis.ticks = element_line(size = 0.2),
      legend.position = "none", 
      axis.title.x = element_text(size = 6, face = "bold"))
  
  pdf4[[index]] <-  ggplot(dat1, aes(importance4, reorder(genus4,importance4))) +
    geom_segment(aes(yend = genus4, xend = 0), color = "gray70", size = 0.7) +
    geom_point(aes(size = abs(richness4)),color="#fb8072")+
    scale_size_continuous(range = c(0.5, 4))+
    theme_bw()+xlab("Mean decrease accuracy")+ylab("")+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold.italic", size = 10),
      strip.background = element_blank(),
      #strip.text.x = element_blank(),
      axis.ticks = element_line(size = 0.2),
      legend.position = "none", 
      axis.title.x = element_text(size = 6, face = "bold"))
  
  index = index + 1
}

plots2 <- c(pdf1, pdf2, pdf3, pdf4)
aligned2 <- align_plots(plotlist = plots2, align = "hv")

p2 <- plot_grid(plotlist = aligned2, ncol = 4, nrow = 2,
                   labels = c("a2","b2","c2","d2","e2","f2","g2","h2"))

header_row <- ggarrange(
  as_ggplot(text_grob("CRC_Baxter (10% training)", face = "bold", size = 10)),
  as_ggplot(text_grob("CRC_Zeller (10% training)", face = "bold", size = 10)),
  as_ggplot(text_grob("CRC_Baxter (33% training)", face = "bold", size = 10)),
  as_ggplot(text_grob("CRC_Zeller (33% training)", face = "bold", size = 10)),
  ncol = 4
)

final_plot <- ggarrange(
  header_row,
  p2,   
  ncol = 1, nrow = 2,
  heights = c(0.1, 1, 1)
)
print(final_plot)

ggsave(plot = final_plot, filename = "../figs/CRC_feature.pdf", width = 16, height = 6)
