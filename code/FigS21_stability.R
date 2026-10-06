library(PreLectR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(cowplot)

setwd("/Users/yrx/Downloads/FeatureSelection/feature_selection/code")

effect_sizes <- c("0.01", "0.05", "0.2")
ffs_list     <- c(0.05, 0.1, 0.2)
n_rea        <- 10

drop_status <- function(dat){
  target_cols <- c("status", "label", "y", "class", "disease")
  idx <- which(tolower(colnames(dat)) %in% target_cols)
  if(length(idx) > 0) return(dat[, -idx, drop=FALSE])
  return(dat)
}

get_prelect_top <- function(dat_raw, top_n) {
  y <- factor(dat_raw[[ncol(dat_raw)]])
  dat_x <- drop_status(dat_raw)
  X_raw <- t(as.matrix(dat_x))
  storage.mode(X_raw) <- "numeric"
  X_scaled <- t(scale(t(log1p(X_raw))))
  X_scaled[is.na(X_scaled)] <- 0
  
  lrange <- AutoScanning(X_scaled, X_raw, y, task="classification", step=10)
  tuning_res <- LambdaTuningParallel(X_scaled, X_raw, y, lrange, n_cores=8, 
                                     outpath=tempdir(), spl_ratio=0.7, task="classification")
  pick <- LambdaDecision(tuning_res$TuningResult, tuning_res$PvlDistSummary)
  pre_out <- PreLect(X_scaled, GetPrevalence(X_raw), y, lambda=pick$opt_lmbd, task="classification")
  featprop <- FeatureProperty(X_raw, y, pre_out, task="classification")
  
  featprop <- featprop[order(abs(featprop$coef), decreasing = TRUE), ]
  return(as.character(featprop$FeatName[1:min(top_n, nrow(featprop))]))
}

stability_results <- data.frame()

for (ffs in ffs_list) {
  message("Calculating stability for Spike-in fraction: ", ffs)
  
  temp_jaccard_storage <- list()
  
  for (rea_id in 1:n_rea) {
    message("  - Processing Replicate: ", rea_id)
    
    method_features <- list(ParSlet = list(), PreLect = list())
    
    for (es in effect_sizes) {
      dat_path <- paste0("../data/Stool", es, "_f_", 2, "_", ffs, ".csv")
      spike_path <- paste0("../data/Spike_Stool", es, "__f_2_", ffs, ".csv")
      features_score_path <- paste0("../results/Simulated", es, "_features_all.csv")
      
      current_spike_data <- read.csv(spike_path)
      current_top_n <- nrow(current_spike_data)
      
      dat_full <- read.csv(dat_path, row.names = 1, header = TRUE)
      current_features_all <- read.csv(features_score_path)
      
      dat_x <- drop_status(dat_full)
      richness_full <- colSums(dat_x > 0) / nrow(dat_x)
      
      ps_feat <- current_features_all %>% 
        filter(ff == ffs, rea == rea_id, frac == max(frac)) %>%
        mutate(ps_score = MeanDecreaseGini * (richness_full[Feature] ^ 2.5)) %>%
        arrange(desc(ps_score)) %>% 
        head(current_top_n) %>% 
        pull(Feature)
      
      set.seed(rea_id)
      pl_feat <- get_prelect_top(dat_full, current_top_n)
      
      clean_name <- function(x) {
        x <- tolower(x)
        x <- gsub("[[:punct:]]+", ".", x)
        return(x)
      }
      
      method_features$ParSlet[[as.character(es)]] <- clean_name(ps_feat)
      method_features$PreLect[[as.character(es)]]  <- clean_name(pl_feat)
    }
    
    es_pairs <- list(c("0.01", "0.05"), c("0.05", "0.2"), c("0.01", "0.2"))
    for (pair in es_pairs) {
      es1 <- pair[1]; es2 <- pair[2]; comp_label <- paste0("ES ", es1, " vs ", es2)
      
      for (m in c("ParSlet", "PreLect")) {
        list1 <- method_features[[m]][[es1]]
        list2 <- method_features[[m]][[es2]]
        
        jaccard_val <- length(intersect(list1, list2)) / length(union(list1, list2))
        
        key <- paste(comp_label, m, sep = "_")
        temp_jaccard_storage[[key]] <- c(temp_jaccard_storage[[key]], jaccard_val)
      }
    }
  }
  
  for (key in names(temp_jaccard_storage)) {
    split_key <- strsplit(key, "_")[[1]]
    comp_name <- split_key[1]
    method_name <- split_key[2]
    
    stability_results <- rbind(stability_results, data.frame(
      SpikeIn = paste0("Spike: ", ffs),
      Comparison = comp_name,
      Method = method_name,
      Stability = mean(temp_jaccard_storage[[key]])
    ))
  }
}

wide_stability <- stability_results %>%
  pivot_wider(
    names_from = Method, 
    values_from = Stability
  )

p <- ggplot(wide_stability, aes(x = PreLect, y = ParSlet)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_point(aes(color = SpikeIn, shape = Comparison), size = 5, alpha = 0.8) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_bw() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Stability Comparison",
    x = "PreLect Stability",
    y = "ParSlet Stability",
    color = "Spike-in Ratio",
    shape = "ES Comparison"
  )
print(p)

ggsave("../figs/FigR4_stability.pdf", p, width = 7, height = 6)
