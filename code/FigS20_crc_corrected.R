library(curatedMetagenomicData)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(randomForest)
library(pROC)
library(ggplot2)
library(ggpubr)
library(cowplot)

setwd("/Users/yrx/Downloads/FeatureSelection/feature_selection/code")
set.seed(1)

gamma <- 2.5
top_prop <- 0.10
ntree_rf <- 500

corrected_t_test <- function(diffs, n_train, n_test) {
  n <- length(diffs)
  if (n < 2) return(NA)
  
  mean_d <- mean(diffs)
  var_d <- var(diffs)
  
  corrected_se <- sqrt((1/n + n_test/n_train) * var_d)
  
  t_stat <- mean_d / corrected_se
  p_val <- 2 * pt(-abs(t_stat), df = n - 1)
  return(p_val)
}

taxa_table <- curatedMetagenomicData(pattern = "relative_abundance", dryrun = FALSE)
merged_taxa_table <- mergeData(taxa_table)

X_all <- assay(merged_taxa_table, "relative_abundance")
meta  <- curatedMetagenomicData::sampleMetadata

meta_sub <- meta[meta$disease %in% c("healthy","IBD","CRC","T2D","RA","ACVD","adenoma","IGT"), ]

all_study <- c(
  "GuptaA_2019","HanniganGD_2017","HallAB_2017","HMP_2019_ibdmdb",
  "IjazUZ_2017","NielsenHB_2014","QinJ_2012","MetaCardis_2020_a",
  "KarlssonFH_2013","SankaranarayananK_2015","YachidaS_2019","ZellerG_2014","HMP_2019_t2d"
)

meta_sub <- meta_sub[meta_sub$study_name %in% all_study, ]
rownames(meta_sub) <- meta_sub$sample_id

keep_samples <- intersect(colnames(X_all), rownames(meta_sub))
meta_sub <- meta_sub[keep_samples, , drop=FALSE]
X_all    <- X_all[, keep_samples, drop=FALSE]

feat_names <- sub(".*s__", "", rownames(X_all))
feat_names <- make.unique(feat_names)
rownames(X_all) <- feat_names

features <- read.csv("../results/Curated_feautures.csv", stringsAsFactors = FALSE)

crc <- features %>%
  filter(phenotype == "CRC") %>%
  mutate(
    Feature = as.character(Feature),
    study   = as.character(study),
    score   = MeanDecreaseGini * (richness ^ gamma)
  ) %>%
  filter(Feature %in% rownames(X_all))

avail_tbl <- meta_sub %>%
  filter(disease %in% c("CRC","healthy")) %>%
  group_by(study_name, disease) %>%
  summarise(n = dplyr::n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = disease, values_from = n, values_fill = 0) %>%
  mutate(total = CRC + healthy)

crc_studies <- avail_tbl %>%
  filter(CRC > 0, healthy > 0) %>%
  pull(study_name) %>%
  intersect(unique(crc$study)) %>%
  sort()

crc_top <- crc %>%
  filter(study %in% crc_studies) %>%
  group_by(study) %>%
  arrange(desc(score), .by_group = TRUE) %>%
  mutate(n_top = max(1L, floor(top_prop * n()))) %>%
  filter(row_number() <= n_top) %>%
  mutate(
    med_prev   = median(richness, na.rm = TRUE),
    prev_group = ifelse(richness >= med_prev, "HighPrev", "LowPrev")
  ) %>%
  ungroup()

feat_lists <- crc_top %>%
  group_by(study, prev_group) %>%
  summarise(features = list(unique(Feature)), .groups="drop")

make_X_aligned <- function(X_df, cols_train) {
  miss <- setdiff(cols_train, colnames(X_df))
  if (length(miss) > 0) {
    for (m in miss) X_df[[m]] <- 0
  }
  X_df[, cols_train, drop = FALSE]
}

res <- list()
idx <- 1
crc_studies <- as.character(crc_studies)

for (st_train in crc_studies) {
  meta_tr <- meta_sub %>% filter(study_name == st_train, disease %in% c("CRC","healthy"))
  samp_tr <- meta_tr$sample_id
  y_tr <- factor(meta_tr$disease, levels = c("healthy","CRC"))
  X_tr_all <- as.data.frame(t(X_all[, samp_tr, drop=FALSE]))
  
  high_feats <- feat_lists %>% filter(study == st_train, prev_group == "HighPrev") %>% pull(features)
  low_feats  <- feat_lists %>% filter(study == st_train, prev_group == "LowPrev")  %>% pull(features)
  if (length(high_feats) == 0 || length(low_feats) == 0) next
  
  rf_high <- randomForest(x = X_tr_all[, high_feats[[1]], drop=FALSE], y = y_tr, ntree = ntree_rf)
  rf_low  <- randomForest(x = X_tr_all[, low_feats[[1]], drop=FALSE], y = y_tr, ntree = ntree_rf)
  
  for (st_test in setdiff(crc_studies, st_train)) {
    meta_te <- meta_sub %>% filter(study_name == st_test, disease %in% c("CRC","healthy"))
    y_te <- factor(meta_te$disease, levels = c("healthy","CRC"))
    if (nlevels(droplevels(y_te)) < 2) next
    
    X_te <- as.data.frame(t(X_all[, meta_te$sample_id, drop=FALSE]))
    prob_high <- predict(rf_high, make_X_aligned(X_te, high_feats[[1]]), type="prob")[, "CRC"]
    prob_low  <- predict(rf_low,  make_X_aligned(X_te, low_feats[[1]]),  type="prob")[, "CRC"]
    
    res[[idx]] <- data.frame(
      train_study = st_train, test_study = st_test,
      n_train = nrow(X_tr_all), n_test = nrow(X_te),
      auc_high = as.numeric(pROC::roc(y_te, prob_high, quiet=TRUE)$auc),
      auc_low  = as.numeric(pROC::roc(y_te, prob_low,  quiet=TRUE)$auc)
    )
    idx <- idx + 1
  }
}
auc_results <- bind_rows(res)

auc_long <- auc_results %>%
  pivot_longer(cols = c(auc_high, auc_low), names_to = "group", values_to = "AUROC") %>%
  mutate(group = recode(group, auc_high = "HighPrev", auc_low = "LowPrev"),
         pair  = paste0(train_study, " -> ", test_study))

auc_diffs <- auc_results %>% mutate(diff = auc_high - auc_low)
p_corrected_main <- corrected_t_test(auc_diffs$diff, mean(auc_results$n_train), mean(auc_results$n_test))
p_label_main <- paste0("Corrected Resampled t-test\np = ", formatC(p_corrected_main, format = "e", digits = 2))

p_auc <- ggplot(auc_long, aes(x = group, y = AUROC, color = group, fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.25) +
  geom_jitter(width = 0.12, size = 1.6, alpha = 0.75) +
  scale_color_manual(values = c("HighPrev" = "#DD5129FF", "LowPrev" = "#0F7BA2FF")) +
  scale_fill_manual(values = c("HighPrev" = "#DD5129FF", "LowPrev" = "#0F7BA2FF")) +
  coord_cartesian(ylim = c(0.45, 1.15)) +
  theme_bw(base_size = 11) +
  labs(title = "Cross-study AUROC", x = NULL, y = "AUROC") +
  annotate("text", x = 1.5, y = 1.1, label = p_label_main, size = 3.2, fontface = "italic") +
  theme(legend.position = "none", panel.grid.major.x = element_blank())

avg_ntr <- mean(auc_results$n_train)
avg_nte <- mean(auc_results$n_test)

pvals_prev <- crc_top %>%
  group_by(study) %>%
  summarise(
    p_val = {
      h <- richness[prev_group == "HighPrev"]
      l <- richness[prev_group == "LowPrev"]
      corrected_t_test(h - mean(l), avg_ntr, avg_nte)
    }, .groups = "drop"
  ) %>%
  mutate(label = paste0("corr. p = ", formatC(p_val, format = "e", digits = 2)))

label_pos <- crc_top %>% group_by(study) %>% summarise(y = max(richness, na.rm=T) * 1.2, .groups="drop")
pvals_plot <- left_join(pvals_prev, label_pos, by = "study")

p_prev <- ggplot(crc_top, aes(x = prev_group, y = richness, color = prev_group, fill = prev_group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.25) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.5) +
  facet_wrap(~ study, scales = "free_y") +
  geom_text(data = pvals_plot, aes(x = 1.5, y = y, label = label), inherit.aes = FALSE, size = 2.6, fontface = "italic") +
  scale_color_manual(values = c("HighPrev" = "#DD5129FF", "LowPrev" = "#0F7BA2FF")) +
  scale_fill_manual(values = c("HighPrev" = "#DD5129FF", "LowPrev" = "#0F7BA2FF")) +
  theme_bw(base_size = 12) +
  labs(title = "Prevalence Distribution of HighPrev vs LowPrev features (top 10%)", x = NULL, y = "Prevalence (richness)") +
  theme(legend.position = "none", strip.background = element_rect(fill = "grey95", color = NA), strip.text = element_text(face = "bold"))

final_plot <- plot_grid(p_auc, p_prev, labels = c("a", "b"), ncol = 2, rel_widths = c(1, 1.6))
print(final_plot)

ggsave("../figs/prev_corrected.pdf", plot = final_plot, width = 13, height = 6)

auroc_medians <- auc_long %>%
  group_by(group) %>%
  summarise(
    median_auroc = median(AUROC, na.rm = TRUE),
    mean_auroc   = mean(AUROC, na.rm = TRUE),
    n_observations = n()
  )

print(auroc_medians)
