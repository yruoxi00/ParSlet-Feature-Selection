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

taxa_table <- curatedMetagenomicData(pattern = "relative_abundance", dryrun = FALSE)
merged_taxa_table <- mergeData(taxa_table)

X_all <- assay(merged_taxa_table, "relative_abundance")
meta  <- curatedMetagenomicData::sampleMetadata

meta_sub <- meta[meta$disease %in% c("healthy","IBD","CRC","T2D","RA","ACVD","adenoma","IGT"), ]

all_study <- c(
  "GuptaA_2019","HanniganGD_2017",
  "HallAB_2017","HMP_2019_ibdmdb","IjazUZ_2017","NielsenHB_2014",
  "QinJ_2012","MetaCardis_2020_a","KarlssonFH_2013","SankaranarayananK_2015",
  "HanniganGD_2017","YachidaS_2019","ZellerG_2014",
  "KarlssonFH_2013","MetaCardis_2020_a","HMP_2019_t2d"
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

corrected_t_test <- function(diffs, n_train, n_test) {
  n <- length(diffs)
  mean_d <- mean(diffs)
  var_d <- var(diffs)
  
  corrected_se <- sqrt((1/n + n_test/n_train) * var_d)
  
  t_stat <- mean_d / corrected_se
  p_val <- 2 * pt(-abs(t_stat), df = n - 1)
  return(p_val)
}

crc <- features %>%
  filter(phenotype == "CRC") %>%
  mutate(
    Feature = as.character(Feature),
    study   = as.character(study),
    score   = MeanDecreaseGini * (richness ^ gamma)
  )

crc <- crc %>% filter(Feature %in% rownames(X_all))

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
  X_df <- X_df[, cols_train, drop = FALSE]
  X_df
}

res <- list()
idx <- 1

crc_studies <- as.character(crc_studies)

for (st_train in crc_studies) {
  st_train <- as.character(st_train)
  
  meta_tr <- meta_sub %>%
    filter(study_name == st_train, disease %in% c("CRC","healthy"))
  
  samp_tr <- meta_tr$sample_id
  y_tr <- factor(meta_tr$disease, levels = c("healthy","CRC"))
  X_tr_all <- as.data.frame(t(X_all[, samp_tr, drop=FALSE]))
  
  high_feats <- feat_lists %>% filter(study == st_train, prev_group == "HighPrev") %>% pull(features)
  low_feats  <- feat_lists %>% filter(study == st_train, prev_group == "LowPrev")  %>% pull(features)
  
  if (length(high_feats) == 0 || length(low_feats) == 0) next
  
  high_feats <- high_feats[[1]]
  low_feats  <- low_feats[[1]]
  
  high_feats_tr <- intersect(high_feats, colnames(X_tr_all))
  low_feats_tr  <- intersect(low_feats,  colnames(X_tr_all))
  
  if (length(high_feats_tr) < 2 || length(low_feats_tr) < 2) next
  
  rf_high <- randomForest(
    x = X_tr_all[, high_feats_tr, drop=FALSE],
    y = y_tr,
    ntree = ntree_rf,
    importance = FALSE
  )
  
  rf_low <- randomForest(
    x = X_tr_all[, low_feats_tr, drop=FALSE],
    y = y_tr,
    ntree = ntree_rf,
    importance = FALSE
  )
  
  for (st_test in setdiff(crc_studies, st_train)) {
    st_test <- as.character(st_test)
    
    meta_te <- meta_sub %>%
      filter(study_name == st_test, disease %in% c("CRC","healthy"))
    
    samp_te <- meta_te$sample_id
    y_te <- factor(meta_te$disease, levels = c("healthy","CRC"))
    if (nlevels(droplevels(y_te)) < 2) next
    
    X_te <- as.data.frame(t(X_all[, samp_te, drop=FALSE]))
    
    X_te_high <- make_X_aligned(X_te, high_feats_tr)
    X_te_low  <- make_X_aligned(X_te, low_feats_tr)
    
    prob_high <- predict(rf_high, X_te_high, type="prob")[, "CRC"]
    prob_low  <- predict(rf_low,  X_te_low,  type="prob")[, "CRC"]
    
    y01 <- as.numeric(y_te == "CRC")
    
    auc_high <- as.numeric(pROC::roc(response = y01, predictor = prob_high, quiet = TRUE)$auc)
    auc_low  <- as.numeric(pROC::roc(response = y01, predictor = prob_low,  quiet = TRUE)$auc)
    
    res[[idx]] <- data.frame(
      train_study = st_train,
      test_study  = st_test,
      n_train     = nrow(X_tr_all),
      n_test      = nrow(X_te),
      n_feat_high = length(high_feats_tr),
      n_feat_low  = length(low_feats_tr),
      auc_high    = auc_high,
      auc_low     = auc_low,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1
  }
}

auc_results <- bind_rows(res)

auc_long <- auc_results %>%
  filter(train_study != test_study) %>%
  pivot_longer(
    cols = c(auc_high, auc_low),
    names_to = "group",
    values_to = "AUROC"
  ) %>%
  mutate(
    group = recode(group, auc_high = "HighPrev", auc_low = "LowPrev"),
    pair  = paste0(train_study, " → ", test_study)
  )

auc_long$pair <- factor(auc_long$pair, levels = unique(auc_long$pair))

wilcox_test_res <- wilcox.test(auc_results$auc_high, auc_results$auc_low, paired = TRUE)
p_val_wilcox <- wilcox_test_res$p.value

p_label <- paste0("Paired Wilcoxon test\n p = ", formatC(p_val_wilcox, format = "e", digits = 2))

p <- ggplot(auc_long, aes(x = group, y = AUROC, color = group, fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.25) +
  geom_jitter(width = 0.12, size = 1.6, alpha = 0.75) +
  geom_line(aes(group = pair), color = "grey", alpha = 0.2, size = 0.3) + 
  scale_color_manual(values = c("HighPrev" = "#DD5129FF", "LowPrev" = "#0F7BA2FF")) +
  scale_fill_manual(values = c("HighPrev" = "#DD5129FF", "LowPrev" = "#0F7BA2FF")) +
  coord_cartesian(ylim = c(0.5, 1.1)) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Cross-study AUROC",
    x = NULL, y = "AUROC"
  ) +
  annotate("text", x = 1.5, y = 1.05, label = p_label, size = 3.5, fontface = "italic")

pvals_prev <- crc_top %>%
  group_by(study) %>%
  summarise(
    p_value = wilcox.test(
      richness[prev_group == "HighPrev"],
      richness[prev_group == "LowPrev"]
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(label = paste0("p = ", formatC(p_value, format = "e", digits = 2)))

label_pos <- crc_top %>%
  group_by(study) %>%
  summarise(y = max(richness, na.rm = TRUE) * 1.05, .groups = "drop")

pvals_plot <- left_join(pvals_prev, label_pos, by = "study")

p_prev <- ggplot(
  crc_top,
  aes(x = prev_group, y = richness, color = prev_group, fill = prev_group)
) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.25) +
  geom_jitter(width = 0.12, size = 1.8, alpha = 0.75) +
  facet_wrap(~ study, scales = "free_y") +
  geom_text(
    data = pvals_plot,
    aes(x = 1.5, y = y, label = label),
    inherit.aes = FALSE,
    size = 3.2
  ) +
  scale_color_manual(values = c("HighPrev" = "#DD5129FF",
                                "LowPrev"  = "#0F7BA2FF")) +
  scale_fill_manual(values = c("HighPrev" = "#DD5129FF",
                               "LowPrev"  = "#0F7BA2FF")) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Prevalence Distribution of HighPrev vs LowPrev features (top 10%)",
    x = NULL,
    y = "Prevalence (richness)"
  )

cp <- plot_grid(
  p, p_prev,
  labels = c("a", "b"),
  ncol = 2, 
  rel_widths = c(1, 1.5)
)
cp

ggsave("../figs/prev.pdf", plot = cp, width = 12, height = 6)
