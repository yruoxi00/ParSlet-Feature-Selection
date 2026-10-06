library(ggplot2)
library(randomForest)

setwd("/Users/yrx/Downloads/FeatureSelection/feature_selection/code")

features = NULL
for (rea in 1:10){
  for (frac in seq(0.1, 1, length.out = 10)){
    for (ff in c(0.05,0.1,0.2)){
      set.seed(rea)
      dat1_sub <- read.csv(file = paste("../data/Stool0.01_f_",2,'_',ff,'.csv',sep = ""),row.names = 1, header = T)
      
      dat1_sub$status <- factor(dat1_sub$status)
      
      dat1_sub <- dat1_sub[sample(1:nrow(dat1_sub), size = frac * nrow(dat1_sub)),]
      train_data <- dat1_sub
      
      set.seed(rea)
      rf_model <- randomForest(
        status ~ .,
        data = train_data,
        importance = TRUE
      )
      
      importance_data <- as.data.frame(importance(rf_model))
      importance_data$Feature <- rownames(importance_data)
      importance_data$ff=ff
      importance_data$frac=frac
      importance_data$rank = rank(importance_data$MeanDecreaseGini)
      importance_data$rea = rea
      features = rbind(features, importance_data)
    }
  }
}
write.csv(features, file = "../results/Simulated0.01_features_all.csv", row.names = F)

features = NULL
for (rea in 1:10){
  for (frac in seq(0.1, 1, length.out = 10)){
    for (ff in c(0.05,0.1,0.2)){
      set.seed(rea)
      dat1_sub <- read.csv(file = paste("../data/Stool0.05_f_",2,'_',ff,'.csv',sep = ""),row.names = 1, header = T)
      
      dat1_sub$status <- factor(dat1_sub$status)
      
      dat1_sub <- dat1_sub[sample(1:nrow(dat1_sub), size = frac * nrow(dat1_sub)),]
      train_data <- dat1_sub 
      
      set.seed(rea)
      rf_model <- randomForest(
        status ~ .,
        data = train_data,
        importance = TRUE
      )
      
      importance_data <- as.data.frame(importance(rf_model))
      importance_data$Feature <- rownames(importance_data)
      importance_data$ff=ff
      importance_data$frac=frac
      importance_data$rank = rank(importance_data$MeanDecreaseGini)
      importance_data$rea = rea
      features = rbind(features, importance_data)
    }
  }
}
write.csv(features, file = "../results/Simulated0.05_features_all.csv", row.names = F)

features = NULL
for (rea in 1:10){
  for (frac in seq(0.1, 1, length.out = 10)){
    for (ff in c(0.05,0.1,0.2)){
      set.seed(rea)
      dat1_sub <- read.csv(file = paste("../data/Stool0.2_f_",2,'_',ff,'.csv',sep = ""),row.names = 1, header = T)
      
      dat1_sub$status <- factor(dat1_sub$status)
      
      dat1_sub <- dat1_sub[sample(1:nrow(dat1_sub), size = frac * nrow(dat1_sub)),]
      train_data <- dat1_sub 
      
      set.seed(rea)
      rf_model <- randomForest(
        status ~ .,
        data = train_data,
        importance = TRUE
      )
      
      importance_data <- as.data.frame(importance(rf_model))
      importance_data$Feature <- rownames(importance_data)
      importance_data$ff=ff
      importance_data$frac=frac
      importance_data$rank = rank(importance_data$MeanDecreaseGini)
      importance_data$rea = rea
      features = rbind(features, importance_data)
    }
  }
}
write.csv(features, file = "../results/Simulated0.2_features_all.csv", row.names = F)
