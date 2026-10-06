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
  new_df <- new_df %>%
    mutate(transformed_richness = mean_richness^(3.0))
  
  fit <- lm(frequency ~ transformed_richness, data = new_df)
  new_df <- new_df %>%
    mutate(fitted_values = predict(fit, newdata = .))
  
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
  new_df <- new_df %>%
    mutate(transformed_richness = mean_richness^(3.0))
  
  fit <- lm(frequency ~ transformed_richness, data = new_df)
  new_df <- new_df %>%
    mutate(fitted_values = predict(fit, newdata = .))
  
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
  new_df <- new_df %>%
    mutate(transformed_richness = mean_richness^(3.0))
  
  fit <- lm(frequency ~ transformed_richness, data = new_df)
  new_df <- new_df %>%
    mutate(fitted_values = predict(fit, newdata = .))
  
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

alphas <- seq(1, 3, by = 0.5)

fit_stats <- lapply(alphas, function(alpha){
  fit <- lm(frequency ~ I(mean_richness^alpha), data = new_df)
  pred <- predict(fit, newdata = new_df)
  rmse <- sqrt(mean((new_df$frequency - pred)^2))
  r2   <- summary(fit)$r.squared
  aic  <- AIC(fit)
  data.frame(alpha = alpha, RMSE = rmse, R2 = r2, AIC = aic)
}) %>% bind_rows()

fit_stats %>% arrange(RMSE) %>% head(5)
best_alpha <- fit_stats$alpha[which.min(fit_stats$RMSE)]
best_alpha

p1 <- ggplot(fit_stats, aes(alpha, RMSE)) +
  geom_line() +
  geom_point(size = 2) +
  geom_vline(xintercept = 2.5, linetype = "dashed") +
  theme_bw() +
  labs(
    x = expression(alpha),
    y = "RMSE",
    title = expression(RMSE ~ "vs Exponent" ~ alpha)
  )

p2 <- ggplot(fit_stats, aes(alpha, R2)) +
  geom_line() +
  geom_point(size = 2) +
  geom_vline(xintercept = 2.5, linetype = "dashed") +
  theme_bw() +
  labs(
    x = expression(alpha),
    y = expression(R^2), 
    title = expression(R^2 ~ "vs Exponent" ~ alpha)
  )

p1 <- p1 +
  geom_vline(xintercept = 2.5, linetype = "dashed", linewidth = 0.4) +
  theme(plot.title = element_text(size = 10))

p2 <- p2 +
  geom_vline(xintercept = 2.5, linetype = "dashed", linewidth = 0.4) +
  theme(plot.title = element_text(size = 10))

p <- ggarrange(
  p1, p2,
  labels = c("a", "b"), 
  ncol = 2, nrow = 1,
  align = "hv"
)
p
