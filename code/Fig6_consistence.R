library(ggplot2)
library(randomForest)
library(pROC)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggpubr)
library(unikn) 
library(gridExtra)

setwd("/Users/yrx/Downloads/FeatureSelection/feature_selection/code")

features <- read.csv(file = "../results/Curated_feautures.csv")
dat_features = NULL

for (threshold in seq(0.02,0.2,by=0.02)){
  for (phenotype in unique(features$phenotype)){
    dat1_sub <- features[features$phenotype==phenotype,]
    taxa_1 = c()
    taxa_2 = c()
    taxa_3 = c()
    for (study1 in unique(dat1_sub$study)){
      i1 <- dat1_sub$study == study1
      g1 <-  dat1_sub$MeanDecreaseGini[i1]
      sorted_indicesg1 <- order(g1,decreasing = T)
      g3 <- dat1_sub$MeanDecreaseGini[i1] *  (dat1_sub$richness[i1])^(2.5)
      sorted_indicesg3 <- order(g3,decreasing = T)
      g4 <- dat1_sub$MeanDecreaseGini[i1] *  (dat1_sub$richness[i1])
      sorted_indicesg4 <- order(g4,decreasing = T)
      
      
      top_n = floor(threshold*length(g1))
      taxa_1 = c(taxa_1,sorted_indicesg1[1:top_n])
      taxa_2 = c(taxa_2,sorted_indicesg3[1:top_n])
      taxa_3 = c(taxa_3,sorted_indicesg4[1:top_n])
    }
    freq_1 <- as.data.frame(table(taxa_1))
    freq_2 <- as.data.frame(table(taxa_2))
    freq_3 <- as.data.frame(table(taxa_3))
    count1 <- sum(freq_1$Freq==max(freq_1$Freq))/top_n/length(unique(dat1_sub$study))
    count2 <- sum(freq_2$Freq==max(freq_2$Freq))/top_n/length(unique(dat1_sub$study))
    count3 <- sum(freq_3$Freq==max(freq_3$Freq))/top_n/length(unique(dat1_sub$study))

    dat_features <- rbind(dat_features, data.frame(phenotype=c(phenotype,phenotype,phenotype),freq=c(count1,count2,count3),
                                                   threshold = c(threshold,threshold,threshold),
                                                   types = c("Original","Integrated (power)","Integrated (linear)")))
    print(paste("Phenotype:", phenotype, "Threshold:", threshold))
    cat("Richness values:", head(dat1_sub$richness[i1], 10), "\n")
    cat(" g1:", head(g1, 10), "\n")
    cat(" g3:", head(g3, 10), "\n")
    cat(" g4:", head(g4, 10), "\n\n")
    cat("  taxa_1 (Original):   ", head(taxa_1, 20), " ... length:", length(taxa_1), "\n")
    cat("  taxa_2 (Power):      ", head(taxa_2, 20), " ... length:", length(taxa_2), "\n")
    cat("  taxa_3 (Linear):     ", head(taxa_3, 20), " ... length:", length(taxa_3), "\n\n")
    print(paste("Original count:", count1))
    print(paste("Power count:", count2))
    print(paste("Linear count:", count3))
    print(top_n)
  }
}


p_values <- dat_features %>%
  filter(types %in% c("Original", "Integrated (power)")) %>%
  pivot_wider(
    id_cols = c(phenotype, threshold),
    names_from = types,
    values_from = freq
  ) %>%
  group_by(phenotype) %>%
  summarise(
    p_value = wilcox.test(`Integrated (power)`, Original, paired = TRUE)$p.value,
    higher_mean = mean(`Integrated (power)`) > mean(Original),  # quick check who is higher
    .groups = "drop"
  )


dat_features <- left_join(dat_features, p_values, by = "phenotype")

p_value_positions <- dat_features %>%
  group_by(phenotype) %>%
  summarise(
    x_pos = max(threshold) * 0.75,  
    y_pos = max(freq) * 1.05      
  )

p_values <- left_join(p_values, p_value_positions, by = "phenotype")

study_counts <- features %>%
  group_by(phenotype) %>%
  summarise(n_studies = n_distinct(study), .groups = "drop") %>%
  mutate(label = paste0(phenotype, " (", n_studies, ")"))

facet_labels <- setNames(study_counts$label, study_counts$phenotype)

gg1 <- ggplot(data = dat_features,aes(threshold,freq,group=types))+
  geom_line(aes(color=types),size=0.3)+
  geom_point(aes(color=types))+facet_wrap(~phenotype, scales = "free", nrow = 1,
                                                     labeller = labeller(phenotype = facet_labels))+
  scale_color_manual(values = c("#5495CFFF", "#DB4743FF","gray70"))+
  geom_text(data = p_values, aes(x = x_pos, y = y_pos, 
                                 label = paste0("p = ", round(p_value, 3))),
            inherit.aes = FALSE, color = "black", size = 3)+
  theme_bw()+ylab("TPR")+xlab("Fraction of selected features")+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 8),
    strip.background = element_blank(),
    axis.text.x = element_text(size = 8,color = 'black'),
    axis.text.y = element_text(size = 8,color = 'black'),
    #strip.text.x = element_blank(),
    axis.ticks = element_line(size = 0.2),
    legend.position = "none")



dat_features = NULL

for (threshold in seq(0.02,0.2,by=0.02)){
  for (phenotype in unique(features$phenotype)){
    dat1_sub <- features[features$phenotype==phenotype,]
    taxa_1 = c()
    taxa_2 = c()
    taxa_3 = c()
    for (study1 in unique(dat1_sub$study)){
      i1 <- dat1_sub$study == study1
      g1 <-  dat1_sub$MeanDecreaseGini[i1]
      sorted_indicesg1 <- order(g1,decreasing = T)
      g3 <- dat1_sub$MeanDecreaseGini[i1] *  (dat1_sub$richness[i1])^(2.5)
      sorted_indicesg3 <- order(g3,decreasing = T)
      g4 <- dat1_sub$MeanDecreaseGini[i1] *  (dat1_sub$richness[i1])
      sorted_indicesg4 <- order(g4,decreasing = T)
      
      
      top_n = floor(threshold*length(g1))
      taxa_1 = c(taxa_1,sorted_indicesg1[1:top_n])
      taxa_2 = c(taxa_2,sorted_indicesg3[1:top_n])
      taxa_3 = c(taxa_3,sorted_indicesg4[1:top_n])
    }
    freq_1 <- as.data.frame(table(taxa_1))
    freq_2 <- as.data.frame(table(taxa_2))
    freq_3 <- as.data.frame(table(taxa_3))
    count1 <- sum(freq_1$Freq==1)/top_n/length(unique(dat1_sub$study))
    count2 <- sum(freq_2$Freq==1)/top_n/length(unique(dat1_sub$study))
    count3 <- sum(freq_3$Freq==1)/top_n/length(unique(dat1_sub$study))
    
    dat_features <- rbind(dat_features, data.frame(phenotype=c(phenotype,phenotype,phenotype),freq=c(count1,count2,count3),
                                                   threshold = c(threshold,threshold,threshold),
                                                   types = c("Original","Integrated (power)","Integrated (linear)")))
  }
}


p_values <- dat_features %>%
  filter(types %in% c("Original", "Integrated (power)")) %>%
  pivot_wider(
    id_cols = c(phenotype, threshold),
    names_from = types,
    values_from = freq
  ) %>%
  group_by(phenotype) %>%
  summarise(
    p_value = wilcox.test(`Integrated (power)`, Original, paired = TRUE)$p.value,
    higher_mean = mean(`Integrated (power)`) > mean(Original),  # quick check who is higher
    .groups = "drop"
  )

dat_features <- left_join(dat_features, p_values, by = "phenotype")

p_value_positions <- dat_features %>%
  group_by(phenotype) %>%
  summarise(
    x_pos = max(threshold) * 0.75,  
    y_pos = max(freq) * 1.05       
  )

p_values <- left_join(p_values, p_value_positions, by = "phenotype")

gg2 <- ggplot(data = dat_features,aes(threshold,freq,group=types))+
  geom_line(aes(color=types),size=0.3)+
  geom_point(aes(color=types))+facet_wrap(~phenotype, scales = "free", nrow = 1,
                                          labeller = labeller(phenotype = facet_labels))+
  scale_color_manual(values = c("#5495CFFF", "#DB4743FF","gray70"))+
  geom_text(data = p_values, aes(x = x_pos, y = y_pos, 
                                 label = paste0("p = ", round(p_value, 3))),
            inherit.aes = FALSE, color = "black", size = 3)+
  theme_bw()+ylab("FPR")+xlab("Fraction of selected features")+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 8),
    strip.background = element_blank(),
    axis.text.x = element_text(size = 8,color = 'black'),
    axis.text.y = element_text(size = 8,color = 'black'),
    #strip.text.x = element_blank(),
    axis.ticks = element_line(size = 0.2),
    legend.position = "none")

gg_tmp <- ggplot(data = dat_features, aes(threshold, freq, group=types)) +
  geom_line(aes(color=types), size=0.3) +
  geom_point(aes(color=types)) +
  facet_wrap(~phenotype, scales = "free", nrow = 1) +
  scale_color_manual(
    values = c("#DB4743FF","#5495CFFF","gray70"),
    labels = c("Gini impurity (Integrated-power)",
               "Gini impurity (Integrated-linear)",
               "Gini impurity (original)")
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.justification = "center"
  ) +
  guides(color = guide_legend(
    override.aes = list(linetype = 0, shape = 16, size = 4) # no line, filled circle
  ))

legend <- cowplot::get_legend(gg_tmp)

gg1 <- gg1 + theme(legend.position="none")
gg2 <- gg2 + theme(legend.position="none")

p1 <- ggarrange(
  gg1,
  legend,
  gg2,
  ncol = 1,
  heights = c(1, 0.08, 1),
  labels = letters[1:2],   
  label.x = 0,                 
  label.y = 1,                 
  hjust = -0.5,                
  vjust = 1.5
)

print(p1)


ggsave(p1,file=paste("../figs/Curated_feature.png",sep = ""),width=10, height=5.5, dpi = 500,scale=0.9)
