library(curatedMetagenomicData)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(randomForest)
library(pROC)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(PreLectR)

set.seed(1)

taxa_table <- curatedMetagenomicData(pattern = "relative_abundance", dryrun = FALSE)
merged_taxa_table <- mergeData(taxa_table)

X_all <- assay(merged_taxa_table, "relative_abundance")
meta  <- curatedMetagenomicData::sampleMetadata

target_studies <- c("GuptaA_2019", "HanniganGD_2017", "YachidaS_2019")

meta_sub <- meta %>%
  filter(disease %in% c("healthy", "CRC"), study_name %in% target_studies)
rownames(meta_sub) <- meta_sub$sample_id

keep_samples <- intersect(colnames(X_all), rownames(meta_sub))
X_all <- X_all[, keep_samples]
meta_sub <- meta_sub[keep_samples, ]

feat_names <- sub(".*s__", "", rownames(X_all))
feat_names <- make.unique(feat_names)
rownames(X_all) <- feat_names

X_scaled <- t(scale(t(X_all)))
X_scaled[is.na(X_scaled)] <- 0
diagnosis <- factor(meta_sub$disease, levels = c("healthy", "CRC"))
prevalence_vec <- rowSums(X_all > 0) / ncol(X_all)

lrange <- AutoScanning(X_scaled, X_all, diagnosis, task = "classification", step = 30)
out_dir_tmp <- tempdir()

tuning_res <- LambdaTuningParallel(
  X_scaled, X_all, diagnosis, lrange,
  outpath = out_dir_tmp,
  n_cores = 8,
  spl_ratio = 0.7,
  task = "classification"
)

lmbd_decision <- LambdaDecision(tuning_res$TuningResult, tuning_res$PvlDistSummary)

PreLect_out <- PreLect(X_scaled, prevalence_vec, diagnosis, 
                       lambda = lmbd_decision$opt_lmbd, task = "classification")
feat_prop <- FeatureProperty(X_all, diagnosis, PreLect_out)

selected_prop <- feat_prop %>% 
  filter(selected == "Selected") %>%
  mutate(
    richness = prevalence, 
    med_val = median(richness),
    prev_group = ifelse(richness >= med_val, "HighPrev", "LowPrev")
  )

make_X_aligned <- function(X_df, cols_train) {
  miss <- setdiff(cols_train, colnames(X_df))
  if (length(miss) > 0) { for (m in miss) X_df[[m]] <- 0 }
  X_df[, cols_train, drop = FALSE]
}

res_pl <- list()
idx <- 1
crc_studies <- target_studies 

for (st_train in crc_studies) {
  meta_tr <- meta_sub %>% filter(study_name == st_train)
  y_tr    <- factor(meta_tr$disease, levels = c("healthy", "CRC"))
  X_tr    <- as.data.frame(t(X_all[, meta_tr$sample_id, drop=FALSE]))
  
  high_feats_tr <- selected_prop %>% filter(prev_group == "HighPrev") %>% pull(FeatName) %>% intersect(colnames(X_tr))
  low_feats_tr  <- selected_prop %>% filter(prev_group == "LowPrev")  %>% pull(FeatName) %>% intersect(colnames(X_tr))
  
  if (length(high_feats_tr) < 2 || length(low_feats_tr) < 2) next
  
  rf_high <- randomForest(x = X_tr[, high_feats_tr], y = y_tr, ntree = 500)
  rf_low  <- randomForest(x = X_tr[, low_feats_tr],  y = y_tr, ntree = 500)
  
  for (st_test in setdiff(crc_studies, st_train)) {
    meta_te <- meta_sub %>% filter(study_name == st_test)
    if (nrow(meta_te) == 0) next
    
    X_te    <- as.data.frame(t(X_all[, meta_te$sample_id, drop=FALSE]))
    y_te_01 <- as.numeric(meta_te$disease == "CRC")
    
    prob_high <- predict(rf_high, make_X_aligned(X_te, high_feats_tr), type="prob")[, "CRC"]
    prob_low  <- predict(rf_low,  make_X_aligned(X_te, low_feats_tr),  type="prob")[, "CRC"]
    
    res_pl[[idx]] <- data.frame(
      train_study = st_train,
      test_study  = st_test,
      auc_high    = as.numeric(pROC::roc(y_te_01, prob_high, quiet=T)$auc),
      auc_low     = as.numeric(pROC::roc(y_te_01, prob_low,  quiet=T)$auc)
    )
    idx <- idx + 1
  }
}

auc_pl_long <- bind_rows(res_pl) %>%
  pivot_longer(cols = c(auc_high, auc_low), names_to = "group", values_to = "AUROC") %>%
  mutate(group = recode(group, auc_high = "HighPrev", auc_low = "LowPrev"))

auc_diffs <- res_pl %>%
  bind_rows() %>%
  mutate(diff = auc_high - auc_low)

obs_t <- t.test(auc_diffs$auc_high, auc_diffs$auc_low, paired = TRUE)$statistic

n_resample <- 999
resample_t <- numeric(n_resample)

for (i in 1:n_resample) {
  signs <- sample(c(-1, 1), nrow(auc_diffs), replace = TRUE)
  resample_t[i] <- t.test(auc_diffs$diff * signs)$statistic
}

p_val_resampled <- sum(abs(resample_t) >= abs(obs_t)) / n_resample

p_a <- ggplot(auc_pl_long, aes(x = group, y = AUROC, color = group, fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.25) +
  geom_jitter(width = 0.12, size = 1.6, alpha = 0.75) +
  scale_color_manual(values = c("HighPrev" = "#DD5129FF", "LowPrev" = "#0F7BA2FF")) +
  scale_fill_manual(values = c("HighPrev" = "#DD5129FF", "LowPrev" = "#0F7BA2FF")) +
  coord_cartesian(ylim = c(0.5, 1.1)) + # Increased limit to fit label
  theme_bw(base_size = 11) +
  labs(title = "PreLect Cross-study AUROC", x = NULL, y = "AUROC") +
  annotate("text", x = 1.5, y = 1.05, 
           label = paste0("Resampled p = ", format.pval(p_val_resampled, digits = 3)),
           size = 3.5, fontface = "italic")

study_prev <- meta_sub %>%
  group_by(study_name) %>%
  do({
    samples <- .$sample_id
    data.frame(
      study    = .$study_name[1],
      Feature  = rownames(X_all),
      richness = as.numeric(rowSums(X_all[, samples] > 0) / length(samples))
    )
  }) %>% ungroup()

plot_b_data <- study_prev %>%
  filter(study %in% target_studies) %>%
  inner_join(selected_prop %>% select(FeatName, prev_group), by = c("Feature" = "FeatName"))

get_resampled_p <- function(df, n = 999) {
  obs_diff <- mean(df$richness[df$prev_group == "HighPrev"]) - 
    mean(df$richness[df$prev_group == "LowPrev"])
  
  null_diffs <- replicate(n, {
    shuffled_labels <- sample(df$prev_group)
    mean(df$richness[shuffled_labels == "HighPrev"]) - 
      mean(df$richness[shuffled_labels == "LowPrev"])
  })
  
  return(sum(abs(null_diffs) >= abs(obs_diff)) / n)
}

facet_pvals <- plot_b_data %>%
  group_by(study) %>%
  do(data.frame(p_resamp = get_resampled_p(.))) %>%
  mutate(label = paste0("Resampled p = ", format.pval(p_resamp, digits = 2)))

p_b <- ggplot(plot_b_data, aes(x = prev_group, y = richness, color = prev_group, fill = prev_group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.25) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.5) +
  facet_wrap(~ study, scales = "free_y") +
  geom_text(data = facet_pvals, aes(x = 1.5, y = Inf, label = label), 
            vjust = 2, inherit.aes = FALSE, size = 3, fontface = "italic") +
  scale_color_manual(values = c("HighPrev" = "#DD5129FF", "LowPrev" = "#0F7BA2FF")) +
  scale_fill_manual(values = c("HighPrev" = "#DD5129FF", "LowPrev" = "#0F7BA2FF")) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", strip.text = element_text(face = "bold")) +
  labs(title = "Prevalence of PreLect Features", x = NULL, y = "Prevalence")

p <- plot_grid(p_a, p_b, labels = c("a", "b"), ncol = 2, rel_widths = c(1, 1.5))
print(p)

ggsave("../figs/r6b.pdf", plot = p, width = 12, height = 6)
