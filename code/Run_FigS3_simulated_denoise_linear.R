library(ggplot2)
library(randomForest)
library(pROC)

setwd("/udd/spxuw/GMrepo/code")
AUC = NULL
for (ff in c(0.2)){
  for (frac in c(0.1,0.3,0.5)){
    for (rea in 1:50){
      set.seed(rea)
      dat1_sub <- read.csv(file = paste("../SparseDOSSA2/Stool0.01_f_",2,'_',ff,'.csv',sep = ""),row.names = 1, header = T)
      meta_data <- read.csv(file = paste('../SparseDOSSA2/Spike_Stool0.01__f_',2,'_',ff,'.csv',sep = ""),row.names = 1, header = T)
      
      dat1_sub$status <- factor(dat1_sub$status)
      train_indices1 <- sample(1:nrow(dat1_sub), size = frac * nrow(dat1_sub))
      test_indices <-  sample(setdiff(1:nrow(dat1_sub),c(train_indices1)), size = 200)
        
      train_data1 <- dat1_sub[train_indices1, ]
      test_data <- dat1_sub[test_indices, ]
      
      set.seed(rea)
      rf_model <- randomForest(
        status ~ .,
        data = train_data1,
        importance = TRUE,
        probability = TRUE
      )
      importance_data <- as.data.frame(importance(rf_model))
      importance_data$Feature <- rownames(importance_data)
      
      top_n = floor(0.1*nrow(importance_data))
      top_features1 <- importance_data$Feature[order(importance_data$MeanDecreaseGini,decreasing = T)[1:top_n]]
      
      pre <- data.frame(richness = colSums(train_data1>0)/nrow(train_data1),species=colnames(train_data1))
      importance_data$pre <- (pre$richness[match(importance_data$Feature,pre$species)])
      top_features2 <- importance_data$Feature[order(importance_data$MeanDecreaseGini*importance_data$pre,decreasing = T)[1:top_n]]
      
      train_data1[[ncol(train_data1)]] <- as.factor(train_data1[[ncol(train_data1)]])
      p_values1 <- sapply(top_features1, function(var) {
        test_result <- wilcox.test(train_data1[[var]] ~ train_data1$status)
        test_result$p.value
      })
      
      p_values11 <- sapply(top_features1, function(var) {
        test_result <- wilcox.test(test_data[[var]] ~ test_data$status)
        test_result$p.value
      })
      
      p_values2 <- sapply(top_features2, function(var) {
        test_result <- wilcox.test(train_data1[[var]] ~ train_data1$status)
        test_result$p.value
      })
      
      p_values21 <- sapply(top_features2, function(var) {
        test_result <- wilcox.test(test_data[[var]] ~ test_data$status)
        test_result$p.value
      })
      
      AUC = rbind(AUC,data.frame(AUC_a=p_values11-p_values1,AUC_g=p_values21-p_values2,ff=ff,frac=frac,rea=rea))
    }
  }
}
write.csv(AUC,file = "Simulated0.01_AUC_cross.csv",row.names = F)


AUC = NULL
for (ff in c(0.2)){
  for (frac in c(0.1,0.3,0.5)){
    for (rea in 1:50){
      set.seed(rea)
      dat1_sub <- read.csv(file = paste("../SparseDOSSA2/Stool0.05_f_",2,'_',ff,'.csv',sep = ""),row.names = 1, header = T)
      meta_data <- read.csv(file = paste('../SparseDOSSA2/Spike_Stool0.05__f_',2,'_',ff,'.csv',sep = ""),row.names = 1, header = T)
      
      dat1_sub$status <- factor(dat1_sub$status)
      train_indices1 <- sample(1:nrow(dat1_sub), size = frac * nrow(dat1_sub))
      test_indices <-  sample(setdiff(1:nrow(dat1_sub),c(train_indices1)), size = 200)
      
      train_data1 <- dat1_sub[train_indices1, ]
      test_data <- dat1_sub[test_indices, ]
      
      set.seed(rea)
      rf_model <- randomForest(
        status ~ .,
        data = train_data1,
        importance = TRUE,
        probability = TRUE
      )
      importance_data <- as.data.frame(importance(rf_model))
      importance_data$Feature <- rownames(importance_data)
      
      top_n = floor(0.1*nrow(importance_data))
      top_features1 <- importance_data$Feature[order(importance_data$MeanDecreaseGini,decreasing = T)[1:top_n]]
      
      pre <- data.frame(richness = colSums(train_data1>0)/nrow(train_data1),species=colnames(train_data1))
      importance_data$pre <- (pre$richness[match(importance_data$Feature,pre$species)])
      top_features2 <- importance_data$Feature[order(importance_data$MeanDecreaseGini*importance_data$pre,decreasing = T)[1:top_n]]
      
      train_data1[[ncol(train_data1)]] <- as.factor(train_data1[[ncol(train_data1)]])
      p_values1 <- sapply(top_features1, function(var) {
        test_result <- wilcox.test(train_data1[[var]] ~ train_data1$status)
        test_result$p.value
      })
      
      p_values11 <- sapply(top_features1, function(var) {
        test_result <- wilcox.test(test_data[[var]] ~ test_data$status)
        test_result$p.value
      })
      
      p_values2 <- sapply(top_features2, function(var) {
        test_result <- wilcox.test(train_data1[[var]] ~ train_data1$status)
        test_result$p.value
      })
      
      p_values21 <- sapply(top_features2, function(var) {
        test_result <- wilcox.test(test_data[[var]] ~ test_data$status)
        test_result$p.value
      })
      
      AUC = rbind(AUC,data.frame(AUC_a=p_values11-p_values1,AUC_g=p_values21-p_values2,ff=ff,frac=frac,rea=rea))
    }
  }
}
write.csv(AUC,file = "Simulated0.05_AUC_cross.csv",row.names = F)


AUC = NULL
for (ff in c(0.2)){
  for (frac in c(0.1,0.3,0.5)){
    for (rea in 1:50){
      set.seed(rea)
      dat1_sub <- read.csv(file = paste("../SparseDOSSA2/Stool0.2_f_",2,'_',ff,'.csv',sep = ""),row.names = 1, header = T)
      meta_data <- read.csv(file = paste('../SparseDOSSA2/Spike_Stool0.2__f_',2,'_',ff,'.csv',sep = ""),row.names = 1, header = T)
      
      dat1_sub$status <- factor(dat1_sub$status)
      train_indices1 <- sample(1:nrow(dat1_sub), size = frac * nrow(dat1_sub))
      test_indices <-  sample(setdiff(1:nrow(dat1_sub),c(train_indices1)), size = 200)
      
      train_data1 <- dat1_sub[train_indices1, ]
      test_data <- dat1_sub[test_indices, ]
      
      set.seed(rea)
      rf_model <- randomForest(
        status ~ .,
        data = train_data1,
        importance = TRUE,
        probability = TRUE
      )
      importance_data <- as.data.frame(importance(rf_model))
      importance_data$Feature <- rownames(importance_data)
      
      top_n = floor(0.1*nrow(importance_data))
      top_features1 <- importance_data$Feature[order(importance_data$MeanDecreaseGini,decreasing = T)[1:top_n]]
      
      pre <- data.frame(richness = colSums(train_data1>0)/nrow(train_data1),species=colnames(train_data1))
      importance_data$pre <- (pre$richness[match(importance_data$Feature,pre$species)])
      top_features2 <- importance_data$Feature[order(importance_data$MeanDecreaseGini*importance_data$pre,decreasing = T)[1:top_n]]
      
      train_data1[[ncol(train_data1)]] <- as.factor(train_data1[[ncol(train_data1)]])
      p_values1 <- sapply(top_features1, function(var) {
        test_result <- wilcox.test(train_data1[[var]] ~ train_data1$status)
        test_result$p.value
      })
      
      p_values11 <- sapply(top_features1, function(var) {
        test_result <- wilcox.test(test_data[[var]] ~ test_data$status)
        test_result$p.value
      })
      
      p_values2 <- sapply(top_features2, function(var) {
        test_result <- wilcox.test(train_data1[[var]] ~ train_data1$status)
        test_result$p.value
      })
      
      p_values21 <- sapply(top_features2, function(var) {
        test_result <- wilcox.test(test_data[[var]] ~ test_data$status)
        test_result$p.value
      })
      
      AUC = rbind(AUC,data.frame(AUC_a=p_values11-p_values1,AUC_g=p_values21-p_values2,ff=ff,frac=frac,rea=rea))
    }
  }
}
write.csv(AUC,file = "Simulated0.2_AUC_cross.csv",row.names = F)
