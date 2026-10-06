library(ggplot2)
library(randomForest)
library(pROC)

setwd("/Users/yrx/Downloads/FeatureSelection/feature_selection/code")
AUC = NULL
features = NULL

for (rea in 1:10){
  for (frac in seq(0.1, 1, length.out = 10)){
    for (ff in c(0.05,0.1,0.2)){
      set.seed(rea)
      dat1_sub <- read.csv(file = paste("../data/Stool0.01_f_",2,'_',ff,'.csv',sep = ""),row.names = 1, header = T)
      meta_data <- read.csv(file = paste('../data/Spike_Stool0.01__f_',2,'_',ff,'.csv',sep = ""),row.names = 1, header = T)
      
      dat1_sub$status <- factor(dat1_sub$status)
      
      dat1_sub <- dat1_sub[sample(1:nrow(dat1_sub), size = frac * ncol(dat1_sub)),]
      train_indices <- sample(1:nrow(dat1_sub), size = 0.8 * nrow(dat1_sub))
      train_data <- dat1_sub[train_indices, ]
      test_data <- dat1_sub[-train_indices, ]
      
      set.seed(rea)
      rf_model <- randomForest(
        status ~ .,
        data = train_data,
        importance = TRUE,
        probability = TRUE
      )
      predictions <- predict(rf_model, test_data, type = "prob")
      roc_curve <- roc(test_data$status, predictions[, 1])
      
      AUC = rbind(AUC,data.frame(AUC=as.numeric(auc(roc_curve)),ff=ff,frac=frac, rea=rea))
      
      importance_data <- as.data.frame(importance(rf_model))
      importance_data$Feature <- rownames(importance_data)
      importance_data$ff=ff
      importance_data$frac=frac
      importance_data$rank = rank(importance_data$MeanDecreaseGini)
      importance_data$rea = rea
      features = rbind(features,importance_data[,3:9])
    }
  }
}

write.csv(AUC,file = "Simulated0.01_AUC.csv",row.names = F)
write.csv(features,file = "Simulated0.01_features.csv",row.names = F)


AUC = NULL
features = NULL
for (rea in 1:10){
  for (frac in seq(0.1, 1, length.out = 10)){
    for (ff in c(0.05,0.1,0.2)){
      set.seed(rea)
      dat1_sub <- read.csv(file = paste("../data/Stool0.05_f_",2,'_',ff,'.csv',sep = ""),row.names = 1, header = T)
      meta_data <- read.csv(file = paste('../data/Spike_Stool0.05__f_',2,'_',ff,'.csv',sep = ""),row.names = 1, header = T)
      
      dat1_sub$status <- factor(dat1_sub$status)
      
      dat1_sub <- dat1_sub[sample(1:nrow(dat1_sub), size = frac * ncol(dat1_sub)),]
      train_indices <- sample(1:nrow(dat1_sub), size = 0.8 * nrow(dat1_sub))
      train_data <- dat1_sub[train_indices, ]
      test_data <- dat1_sub[-train_indices, ]
      
      set.seed(rea)
      rf_model <- randomForest(
        status ~ .,
        data = train_data,
        importance = TRUE,
        probability = TRUE
      )
      predictions <- predict(rf_model, test_data, type = "prob")
      roc_curve <- roc(test_data$status, predictions[, 1])
      
      AUC = rbind(AUC,data.frame(AUC=as.numeric(auc(roc_curve)),ff=ff,frac=frac, rea=rea))
      
      importance_data <- as.data.frame(importance(rf_model))
      importance_data$Feature <- rownames(importance_data)
      importance_data$ff=ff
      importance_data$frac=frac
      importance_data$rank = rank(importance_data$MeanDecreaseGini)
      importance_data$rea = rea
      features = rbind(features,importance_data[,3:9])
    }
  }
}

write.csv(AUC,file = "Simulated0.05_AUC.csv",row.names = F)
write.csv(features,file = "Simulated0.05_features.csv",row.names = F)



AUC = NULL
features = NULL
for (rea in 1:10){
  for (frac in seq(0.1, 1, length.out = 10)){
    for (ff in c(0.05,0.1,0.2)){
      set.seed(rea)
      dat1_sub <- read.csv(file = paste("../data/Stool0.2_f_",2,'_',ff,'.csv',sep = ""),row.names = 1, header = T)
      meta_data <- read.csv(file = paste('../data/Spike_Stool0.2__f_',2,'_',ff,'.csv',sep = ""),row.names = 1, header = T)
      
      dat1_sub$status <- factor(dat1_sub$status)
      
      dat1_sub <- dat1_sub[sample(1:nrow(dat1_sub), size = frac * ncol(dat1_sub)),]
      train_indices <- sample(1:nrow(dat1_sub), size = 0.8 * nrow(dat1_sub))
      train_data <- dat1_sub[train_indices, ]
      test_data <- dat1_sub[-train_indices, ]
      
      set.seed(rea)
      rf_model <- randomForest(
        status ~ .,
        data = train_data,
        importance = TRUE,
        probability = TRUE
      )
      predictions <- predict(rf_model, test_data, type = "prob")
      roc_curve <- roc(test_data$status, predictions[, 1])
      
      AUC = rbind(AUC,data.frame(AUC=as.numeric(auc(roc_curve)),ff=ff,frac=frac, rea=rea))
      
      importance_data <- as.data.frame(importance(rf_model))
      importance_data$Feature <- rownames(importance_data)
      importance_data$ff=ff
      importance_data$frac=frac
      importance_data$rank = rank(importance_data$MeanDecreaseGini)
      importance_data$rea = rea
      features = rbind(features,importance_data[,3:9])
    }
  }
}

write.csv(AUC,file = "Simulated0.2_AUC.csv",row.names = F)
write.csv(features,file = "Simulated0.2_features.csv",row.names = F)