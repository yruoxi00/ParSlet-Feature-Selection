library(ggplot2)
library(randomForest)
library(pROC)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggpubr)
library(paletteer)

setwd("/Users/yrx/Downloads/FeatureSelection/feature_selection/code")

features <- read.csv(file = "../results/Simulated0.01_features.csv")
ngroup <-  length(unique(features$frac))
pdf1 <- list()
index <- 1

for (ffs in c(0.05,0.1,0.2)){
  dat1 <- NULL
  dat1_sub <- read.csv(file = paste("../data/Stool0.01_f_",2,'_',ffs,'.csv',sep = ""),row.names = 1, header = T)
  for (rea in 1:10){
    for (i in 1:ngroup){
      features_sub <- features[features$rea==rea & features$ff==ffs,]
      i1 <- features_sub$frac == unique(features_sub$frac)[i]
      g1 <-  features_sub$MeanDecreaseGini[i1]
      a1 <-  features_sub$MeanDecreaseAccuracy[i1]
      sorted_indicesg1 <- order(g1,decreasing = T)
      sorted_indicesa1 <- order(a1,decreasing = T)
      
      set.seed(rea)
      train_table1 <- dat1_sub[sample(1:nrow(dat1_sub), size = unique(features_sub$frac)[i] * ncol(dat1_sub)),]
      
      data_rich <- data.frame(species=colnames(train_table1),richness1=colSums(train_table1>0)/nrow(train_table1))
      
      r1 <- as.numeric(data_rich$richness1[match(features_sub$Feature[i1],data_rich$species)])
      
      top_n = floor(0.1*length(g1))
      
      dat1 <- rbind(dat1, data.frame(features = sorted_indicesg1[1:top_n],richness=r1[sorted_indicesg1[1:top_n]]))
      
    }
  }
  
  new_df <- dat1 %>%
    group_by(features) %>%
    summarise(
      frequency = n(),           # Count the frequency of each feature
      mean_richness = mean(richness, na.rm = TRUE)  # Compute the mean richness
    )
  
  new_df$frequency <- new_df$frequency/10/ngroup
  
  fit <- nls(
    frequency ~ 1 / (1 + exp(-(b0 + b1 * mean_richness))),
    data = new_df,
    start = list(b0 = -2, b1 = 5),
    control = nls.control(maxiter = 200)
  )
  
  new_df$fitted_values <- predict(fit, newdata = new_df)

  pdf1[[index]] <- ggplot(new_df, aes(mean_richness, frequency)) +
    geom_point(color="#0F7BA2FF",size=2,alpha=0.6) +
    geom_line(aes(y = fitted_values), color = "#DD5129FF", size = 1) +
    theme_bw()+xlab("")+ylab("")+labs(title = paste("Fraction of spiking-in:", ffs))+
    theme(
      plot.title = element_text(hjust = 0.5, size = 8),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold.italic", size = 10),
      strip.background = element_blank(),
      axis.text.x = element_text(size = 8,color = 'black'),
      axis.text.y = element_text(size = 8,color = 'black'),
      #strip.text.x = element_blank(),
      axis.ticks = element_line(size = 0.2),
      legend.position = "none")
  
  index <- index + 1
}

pdf1[[1]] <- pdf1[[1]] + ylab("Frequency")

p1 = do.call(ggarrange, c(pdf1, list(labels = c("a", "b", "c"),ncol = 3, nrow = 1),align="hv"))

##############################
features <- read.csv(file = "../results/Simulated0.05_features.csv")
ngroup <-  length(unique(features$frac))
pdf1 <- list()
index <- 1

for (ffs in c(0.05,0.1,0.2)){
  dat1 <- NULL
  dat1_sub <- read.csv(file = paste("../data/Stool0.05_f_",2,'_',ffs,'.csv',sep = ""),row.names = 1, header = T)
  for (rea in 1:10){
    for (i in 1:ngroup){
      features_sub <- features[features$rea==rea & features$ff==ffs,]
      i1 <- features_sub$frac == unique(features_sub$frac)[i]
      g1 <-  features_sub$MeanDecreaseGini[i1]
      a1 <-  features_sub$MeanDecreaseAccuracy[i1]
      sorted_indicesg1 <- order(g1,decreasing = T)
      sorted_indicesa1 <- order(a1,decreasing = T)
      
      set.seed(rea)
      train_table1 <- dat1_sub[sample(1:nrow(dat1_sub), size = unique(features_sub$frac)[i] * ncol(dat1_sub)),]
      
      data_rich <- data.frame(species=colnames(train_table1),richness1=colSums(train_table1>0)/nrow(train_table1))
      
      r1 <- as.numeric(data_rich$richness1[match(features_sub$Feature[i1],data_rich$species)])
      
      top_n = floor(0.1*length(g1))
      
      dat1 <- rbind(dat1, data.frame(features = sorted_indicesg1[1:top_n],richness=r1[sorted_indicesg1[1:top_n]]))
      
    }
  }
  
  new_df <- dat1 %>%
    group_by(features) %>%
    summarise(
      frequency = n(),           # Count the frequency of each feature
      mean_richness = mean(richness, na.rm = TRUE)  # Compute the mean richness
    )
  
  new_df$frequency <- new_df$frequency/10/ngroup
  
  fit <- nls(
    frequency ~ 1 / (1 + exp(-(b0 + b1 * mean_richness))),
    data = new_df,
    start = list(b0 = -2, b1 = 5),
    control = nls.control(maxiter = 200)
  )
  
  new_df$fitted_values <- predict(fit, newdata = new_df)
  
  pdf1[[index]] <- ggplot(new_df, aes(mean_richness, frequency)) +
    geom_point(color="#0F7BA2FF",size=2,alpha=0.6) +
    geom_line(aes(y = fitted_values), color = "#DD5129FF", size = 1) +
    theme_bw()+xlab("")+ylab("")+
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
  
  index <- index + 1
}

pdf1[[1]] <- pdf1[[1]] + ylab("Frequency")

p2 = do.call(ggarrange, c(pdf1, list(labels = c("d", "e", "f"),ncol = 3, nrow = 1),align="hv"))


#############################################
features <- read.csv(file = "../results/Simulated0.2_features.csv")
ngroup <-  length(unique(features$frac))
pdf1 <- list()
index <- 1

for (ffs in c(0.05,0.1,0.2)){
  dat1 <- NULL
  dat1_sub <- read.csv(file = paste("../data/Stool0.2_f_",2,'_',ffs,'.csv',sep = ""),row.names = 1, header = T)
  for (rea in 1:10){
    for (i in 1:ngroup){
      features_sub <- features[features$rea==rea & features$ff==ffs,]
      i1 <- features_sub$frac == unique(features_sub$frac)[i]
      g1 <-  features_sub$MeanDecreaseGini[i1]
      a1 <-  features_sub$MeanDecreaseAccuracy[i1]
      sorted_indicesg1 <- order(g1,decreasing = T)
      sorted_indicesa1 <- order(a1,decreasing = T)
      
      set.seed(rea)
      train_table1 <- dat1_sub[sample(1:nrow(dat1_sub), size = unique(features_sub$frac)[i] * ncol(dat1_sub)),]
      
      data_rich <- data.frame(species=colnames(train_table1),richness1=colSums(train_table1>0)/nrow(train_table1))
      
      r1 <- as.numeric(data_rich$richness1[match(features_sub$Feature[i1],data_rich$species)])
      
      top_n = floor(0.1*length(g1))
      
      dat1 <- rbind(dat1, data.frame(features = sorted_indicesg1[1:top_n],richness=r1[sorted_indicesg1[1:top_n]]))
      
    }
  }
  
  new_df <- dat1 %>%
    group_by(features) %>%
    summarise(
      frequency = n(),           # Count the frequency of each feature
      mean_richness = mean(richness, na.rm = TRUE)  # Compute the mean richness
    )
  
  new_df$frequency <- new_df$frequency/10/ngroup
  
  fit <- nls(
    frequency ~ 1 / (1 + exp(-(b0 + b1 * mean_richness))),
    data = new_df,
    start = list(b0 = -2, b1 = 5),
    control = nls.control(maxiter = 200)
  )
  
  new_df$fitted_values <- predict(fit, newdata = new_df)
  
  pdf1[[index]] <- ggplot(new_df, aes(mean_richness, frequency)) +
    geom_point(color="#0F7BA2FF",size=2,alpha=0.6) +
    geom_line(aes(y = fitted_values), color = "#DD5129FF", size = 1) +
    theme_bw()+xlab("Prevalence")+ylab("")+
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
  
  index <- index + 1
}

pdf1[[1]] <- pdf1[[1]] + ylab("Frequency")

p3 = do.call(ggarrange, c(pdf1, list(labels = c("g", "h", "i"),ncol = 3, nrow = 1),align="hv"))

p4 = ggarrange(p1,p2,p3,ncol = 1,nrow=3,align = "hv")
print(p4)

row_label_1 <- text_grob("Effect Size: 0.01", rot = 90, size = 10, face = "plain")
row_label_2 <- text_grob("Effect Size: 0.05", rot = 90, size = 10, face = "plain")
row_label_3 <- text_grob("Effect Size: 0.2",  rot = 90, size = 10, face = "plain")

row1 <- ggarrange(row_label_1, p1, ncol = 2, widths = c(0.05, 1))
row2 <- ggarrange(row_label_2, p2, ncol = 2, widths = c(0.05, 1))
row3 <- ggarrange(row_label_3, p3, ncol = 2, widths = c(0.05, 1))

p4 <- ggarrange(row1, row2, row3, ncol = 1, nrow = 3, align = "v")
print(p4)

ggsave(p4,file=paste("../figs/test_correlation_logi.pdf",sep = ""))

# R2 <- function(y, yhat){
#   1 - sum((y - yhat)^2) / sum((y - mean(y))^2)
# }
# 
# fit_power <- lm(frequency ~ I(mean_richness^2.5), data = new_df)
# pred_power <- predict(fit_power, new_df)
# 
# R2_power <- R2(new_df$frequency, pred_power)
# R2_power
# 
# fit_logistic <- nls(
#   frequency ~ 1 / (1 + exp(-(b0 + b1 * mean_richness))),
#   data = new_df,
#   start = list(b0 = -2, b1 = 5)
# )
# pred_logistic <- predict(fit_logistic, new_df)
# 
# R2_logistic <- R2(new_df$frequency, pred_logistic)
# R2_logistic
# 
# deg <- 2
# fit_poly <- lm(
#   frequency ~ poly(mean_richness, deg, raw = TRUE),
#   data = new_df
# )
# pred_poly <- predict(fit_poly, newdata = new_df)
# 
# R2_poly <- 1 - sum((new_df$frequency - pred_poly)^2) / sum((new_df$frequency - mean(new_df$frequency))^2)
# R2_poly
# 
# fit_gam <- gam(frequency ~ s(mean_richness, k = 6),
#            data = new_df, method = "REML")
# pred_gam <- predict(fit_gam, newdata = new_df, type = "response")
# 
# R2_gam <- 1 - sum((new_df$frequency - pred_gam)^2) /
#   sum((new_df$frequency - mean(new_df$frequency))^2)
# 
# R2_gam
# 
# fit_spline <- glm(
#   frequency ~ ns(mean_richness, df = 4),
#   data = new_df,
#   family = binomial(link = "logit")
# )
# pred_spline <- predict(fit_spline, newdata = new_df, type = "response")
# 
# R2_spline <- 1 - sum((new_df$frequency - pred_spline)^2) /
#   sum((new_df$frequency - mean(new_df$frequency))^2)
# 
# R2_spline
# 
# r2_df <- data.frame(
#   Model = c("Power", "Polynomial", "GAM", "Logistic", "Spline"),
#   R2 = c(R2_power, R2_poly, R2_gam, R2_logistic, R2_spline)
# )
# 
# r2_df$Model <- factor(r2_df$Model,
#                       levels = c("Power", "Polynomial", "GAM", "Logistic", "Spline"))
# 
# p <- ggplot(r2_df, aes(x = Model, y = R2, group = 1)) +
#   geom_line(color = "gray60", linewidth = 0.8) +
#   geom_point(size = 3) +
#   geom_text(aes(label = round(R2, 3)), vjust = -1, size = 3) +
#   ylim(0, 1) +
#   theme_bw() +
#   labs(y = expression(R^2), x = "") +
#   theme(
#     axis.text.x = element_text(face = "bold"),
#     plot.title = element_text(hjust = 0.5)
#   )
# p
# 
# ggsave(p,file=paste("../figs/test_correlation_r2.pdf",sep = ""))
