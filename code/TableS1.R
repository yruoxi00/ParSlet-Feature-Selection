library(glmnet)
library(CORElearn) 
library(mRMRe)    
library(dplyr)
library(stringr)
library(tidyr)
library(PreLectR)

setwd("/Users/yrx/Downloads/FeatureSelection/feature_selection/code")

features_filter <- read.csv(file = "../results/Real_feautures_filter.csv")

extract_genus <- function(x) {
  result <- "unclassified"
  clean_x <- gsub("[._]$","", x)
  
  patterns <- c("g__", "f__", "o__", "c__", "p__", "k__")
  
  for (p in patterns) {
    regex_pattern <- paste0(p, "(.*?)(?=\\.[a-z]__|;|$)")
    match <- str_match(clean_x, regex_pattern)[,2]
    
    if (!is.na(match)) {
      match <- gsub("^[._]|[._]$", "", match)
      if (match != "" && match != ".") return(match)
    }
  }
  return(result)
}

target_taxa <- c("Ruminococcus", "Clostridium_XVIII", "Faecalibacterium", "Blautia")
results_list <- list()
run_id <- 1

for (ffs in unique(features_filter$ff)[6:7]) {
  taxa_table <- read.csv(file = paste0("../data/", ffs))
  
  taxa_table$status_num <- as.numeric(as.factor(taxa_table$status))
  
  total_features <- ncol(taxa_table) - 2
  top_n <- max(1, floor(0.10 * total_features)) 
  
  set.seed(1)
  train_idx1 <- sample(1:nrow(taxa_table), size = floor(0.1 * nrow(taxa_table)))
  set.seed(8)
  train_idx2 <- sample(1:nrow(taxa_table), size = floor(0.33 * nrow(taxa_table)))
  indices_list <- list("10%" = train_idx1, "33%" = train_idx2)
  
  for (label in names(indices_list)) {
    train_data_raw <- taxa_table[indices_list[[label]], ]
    
    feat_cols_idx <- which(!colnames(train_data_raw) %in% c("status", "status_num"))
    
    is_not_constant <- sapply(train_data_raw[, feat_cols_idx], function(x) length(unique(x)) > 1)
    valid_feats <- names(is_not_constant)[is_not_constant]
    
    train_data <- train_data_raw[, c(valid_feats, "status", "status_num")]
    
    current_total_features <- length(valid_feats)
    top_n <- max(1, floor(0.10 * current_total_features)) 
    
    x_train <- as.matrix(train_data[, valid_feats])
    y_train <- as.factor(train_data$status)
    
    # --- LASSO ---
    cv_lasso <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 1)
    lasso_coefs <- as.matrix(coef(cv_lasso, s = "lambda.min"))
    lasso_df <- data.frame(feat = rownames(lasso_coefs), coef = abs(lasso_coefs[,1])) %>%
      filter(feat != "(Intercept)" & coef > 0) %>%
      arrange(desc(coef))
    sel_lasso <- head(lasso_df$feat, top_n)
    
    # --- Elastic Net ---
    cv_enet <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 0.5)
    enet_coefs <- as.matrix(coef(cv_enet, s = "lambda.min"))
    enet_df <- data.frame(feat = rownames(enet_coefs), coef = abs(enet_coefs[,1])) %>%
      filter(feat != "(Intercept)" & coef > 0) %>%
      arrange(desc(coef))
    sel_enet <- head(enet_df$feat, top_n)
    
    # --- ReliefF ---
    relief_scores <- attrEval(status ~ ., data = train_data[, !colnames(train_data) == "status_num"], estimator = "ReliefFequalK")
    sel_relief <- names(sort(relief_scores, decreasing = TRUE)[1:top_n])
    
    # --- mRMRe ---
    mrmr_in <- train_data[, !colnames(train_data) == "status"]
    mrmr_data <- mRMR.data(data = mrmr_in)
    mrmr_res <- mRMR.classic(data = mrmr_data, target_indices = which(colnames(mrmr_in) == "status_num"), feature_count = top_n)
    sel_mrmr <- colnames(mrmr_in)[solutions(mrmr_res)[[1]]]
    
    methods <- list(LASSO = sel_lasso, Enet = sel_enet, Relief = sel_relief, mRMR = sel_mrmr)
    
    for (m_name in names(methods)) {
      selected_features <- methods[[m_name]]
      selected_genera <- unique(vapply(selected_features, extract_genus, ""))
      
      run_results <- data.frame(
        Dataset = ffs,
        TrainSize = label,
        Method = m_name,
        Taxon = target_taxa,
        Identified = target_taxa %in% selected_genera
      )
      results_list[[run_id]] <- run_results
      run_id <- run_id + 1
    }
  }
}

final_df <- do.call(rbind, results_list)
summary_table <- final_df %>%
  group_by(Method, Taxon) %>%
  summarize(Times_Identified = sum(Identified), .groups = 'drop') %>%
  pivot_wider(names_from = Method, values_from = Times_Identified)

print(summary_table)

output_dir <- "./prelect_output"
if(!dir.exists(output_dir)) dir.create(output_dir)

extract_genus <- function(x) {
  result <- "unclassified"
  clean_x <- gsub("[._]$","", x)
  patterns <- c("g__", "f__", "o__", "c__", "p__", "k__")
  for (p in patterns) {
    match <- str_match(clean_x, paste0(p, "([^.; ]*)"))[,2]
    if (!is.na(match) && match != "" && match != "." && match != "_") return(match)
  }
  return(result)
}

target_taxa <- c("Ruminococcus", "Clostridium_XVIII", "Faecalibacterium", "Blautia")
results_list <- list()
run_id <- 1

for (ffs in unique(features_filter$ff)[6:7]) {
  
  taxa_table <- read.csv(file = paste0("../data/", ffs))
  taxa_table$status_fac <- factor(taxa_table$status)
  
  set.seed(1)
  train_idx1 <- sample(1:nrow(taxa_table), size = floor(0.1 * nrow(taxa_table)))
  set.seed(8)
  train_idx2 <- sample(1:nrow(taxa_table), size = floor(0.33 * nrow(taxa_table)))
  indices_list <- list("10%" = train_idx1, "33%" = train_idx2)
  
  for (label in names(indices_list)) {
    train_data <- taxa_table[indices_list[[label]], ]
    
    feat_cols_all <- which(!colnames(train_data) %in% c("status", "status_fac"))
    is_constant <- sapply(train_data[, feat_cols_all], function(x) length(unique(x)) == 1)
    cols_to_keep <- names(is_constant)[!is_constant]
    
    top_n <- max(1, floor(0.10 * length(cols_to_keep)))
    
    X_raw <- t(as.matrix(train_data[, cols_to_keep]))
    X_scaled <- t(scale(t(X_raw)))
    diagnosis <- train_data$status_fac
    
    message(paste("Scanning Lambda for:", ffs, "| Split:", label))
    lrange <- AutoScanning(X_scaled, X_raw, diagnosis, task = "classification", step = 30)
    
    tuning_res <- LambdaTuning(X_scaled, X_raw, diagnosis, lrange, 
                               outpath = output_dir, 
                               spl_ratio = 0.7, 
                               task = "classification")
    
    lmbd_picking <- LambdaDecision(tuning_res$TuningResult, tuning_res$PvlDistSummary)
    
    prevalence <- GetPrevalence(X_raw)
    prelect_out <- PreLect(X_scaled, prevalence, diagnosis, 
                           lambda = lmbd_picking$opt_lmbd, 
                           task = "classification", 
                           max_iter = 100000)
    
    prelect_coefs <- unlist(prelect_out$coef)
    
    prelect_res_df <- data.frame(
      Feature = cols_to_keep, 
      Coef = as.numeric(prelect_coefs)
    ) %>%
      mutate(AbsCoef = abs(Coef)) %>%
      arrange(desc(AbsCoef))
    
    selected_features <- head(prelect_res_df$Feature, top_n)
    selected_genera <- unique(vapply(selected_features, extract_genus, ""))
    
    run_results <- data.frame(
      Dataset = ffs,
      TrainSize = label,
      Taxon = target_taxa,
      Identified = target_taxa %in% selected_genera
    )
    
    results_list[[run_id]] <- run_results
    run_id <- run_id + 1
  }
}

final_prelect_df <- do.call(rbind, results_list)

summary_prelect <- final_prelect_df %>%
  group_by(Taxon) %>%
  summarize(
    Times_Identified = sum(Identified),
    Total_Runs = n(),
    Consistency_Percentage = (sum(Identified) / n()) * 100
  )

print(as.data.frame(summary_prelect))