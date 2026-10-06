library(ggplot2)
library(randomForest)
library(pROC)

setwd("/udd/spxuw/GMrepo/code")
AUC = NULL
features = NULL
all_files <- list.files(paste("../processed_Hackathon/"))
all_files <- setdiff(all_files,c("Blueberry.csv","Exercise.csv","hiv_noguerajulian.csv","Ji_WTP_DS.csv","ob_ross.csv","ob_turnbaugh.csv",
                                 "ob_zupancic.csv","Office.csv","seston_plastic_mccormick.csv","wood_plastic_kesy.csv","asd_son.csv","t1d_alkanani.csv"))

for (ff in all_files){
  dat1_sub <- read.csv(file = paste("../processed_Hackathon/",ff,sep = ""),row.names = 1, header = T)
  upper_f <- as.numeric(substr((nrow(dat1_sub)-20)/(ncol(dat1_sub)), 1, 3))
  upper_f <- min(upper_f,1)
  for (frac in seq(0.1, upper_f, length.out = 10)){
    for (rea in 1:10){
      set.seed(rea)
      dat1_sub <- read.csv(file = paste("../processed_Hackathon/",ff,sep = ""),row.names = 1, header = T)
      
      dat1_sub$status <- factor(dat1_sub$status)
      train_indices <- sample(1:nrow(dat1_sub), size = frac * (ncol(dat1_sub)-1))
      test_indices <-  setdiff(1:nrow(dat1_sub),train_indices)
      
      train_data <- dat1_sub[train_indices, ]
      test_data <- dat1_sub[test_indices, ]
      
      set.seed(rea)
      rf_model <- randomForest(
        status ~ .,
        data = train_data,
        importance = TRUE,
        probability = TRUE
      )
      predictions <- predict(rf_model, test_data, type = "prob")
      roc_curve <- roc(test_data$status, predictions[, 1])
      
      AUC = rbind(AUC,data.frame(AUC=as.numeric(auc(roc_curve)),ff=ff,frac=frac, rea=rea,nfeatures=dim(dat1_sub)[1],nsample=dim(dat1_sub)[2]))
      
      importance_data <- as.data.frame(importance(rf_model))
      importance_data$Feature <- rownames(importance_data)

      importance_data = importance_data[,3:5]
      importance_data$ff = ff
      importance_data$frac = frac
      importance_data$rea = rea
      features = rbind(features,importance_data)
    }
  }
}

write.csv(features,file = "Real_feautures.csv",row.names = F)
write.csv(AUC,file = "Real_AUC.csv",row.names = F)

