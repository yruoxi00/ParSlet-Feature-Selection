library(curatedMetagenomicData)
library(SummarizedExperiment)
library(PreLectR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

setwd("/Users/yrx/Downloads/FeatureSelection/feature_selection/code")
set.seed(1)

thresholds <- c(0.05, 0.10)

features_all <- read.csv("../results/Curated_feautures.csv", stringsAsFactors = FALSE)
features_all <- as.data.frame(lapply(features_all, function(x) {
  if (is(x, "Rle")) return(as.vector(x))
  return(x)
}))

taxa_table <- curatedMetagenomicData(pattern = "relative_abundance", dryrun = FALSE)
merged_taxa_table <- mergeData(taxa_table)

X_all <- assay(merged_taxa_table, "relative_abundance")
meta  <- curatedMetagenomicData::sampleMetadata

all_study <- c(
  "GuptaA_2019","HanniganGD_2017","HallAB_2017","HMP_2019_ibdmdb",
  "IjazUZ_2017","NielsenHB_2014","QinJ_2012","MetaCardis_2020_a",
  "KarlssonFH_2013","SankaranarayananK_2015","YachidaS_2019","ZellerG_2014","HMP_2019_t2d"
)

meta_sub <- meta %>%
  filter(study_name %in% all_study, disease %in% c("CRC", "healthy"))

keep_samples <- intersect(colnames(X_all), meta_sub$sample_id)
X_all <- X_all[, keep_samples]
meta_sub <- meta_sub[match(keep_samples, meta_sub$sample_id), ]

feat_names <- sub(".*s__", "", rownames(X_all))
feat_names <- make.unique(feat_names)
rownames(X_all) <- feat_names

crc <- features_all %>%
  filter(phenotype == "CRC") %>%
  mutate(
    Feature = as.character(Feature),
    study   = as.character(study),
    score   = MeanDecreaseGini * (richness ^ 2.5)
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

overlap_results <- data.frame()

for (tp in thresholds) {
  parslet_lists <- list()
  prelect_lists <- list()
  
  for (st in crc_studies) {
    message("Processing study: ", st, " at Threshold: ", tp)
    
    idx <- which(meta_sub$study_name == st)
    sub_X <- as.data.frame(t(as.matrix(X_all[, idx])))
    sub_X$label <- factor(meta_sub$disease[idx])
    
    n_feat <- ncol(sub_X) - 1
    top_n <- max(1, floor(tp * n_feat))
    
    study_fs <- features_all %>% 
      filter(study == st, phenotype == "CRC") %>%
      mutate(ps_score = as.numeric(MeanDecreaseGini) * (as.numeric(richness) ^ 2.5)) %>%
      arrange(desc(ps_score)) %>%
      head(top_n)
    
    parslet_lists[[st]] <- as.character(study_fs$Feature)
    prelect_lists[[st]] <- get_prelect_top(sub_X, top_n)
  }
  
  for (i in 1:(length(crc_studies)-1)) {
    for (j in (i+1):length(crc_studies)) {
      s1 <- crc_studies[i]; s2 <- crc_studies[j]
      
      ps_ov <- length(intersect(parslet_lists[[s1]], parslet_lists[[s2]])) / top_n
      pl_ov <- length(intersect(prelect_lists[[s1]], prelect_lists[[s2]])) / top_n
      
      overlap_results <- rbind(overlap_results, data.frame(
        Comparison = paste(s1, "vs", s2),
        Threshold = factor(paste0(tp*100, "%"), levels=c("5%", "10%")),
        ParSlet_Overlap = ps_ov,
        PreLect_Overlap = pl_ov
      ))
    }
  }
}

p <- ggplot(overlap_results, aes(x = PreLect_Overlap, y = ParSlet_Overlap)) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_jitter(aes(color = Threshold, shape = Threshold), 
              size = 4, width = 0.005, height = 0.005, alpha = 0.8) +
  geom_text_repel(aes(label = Comparison), size = 2.5, color = "black", max.overlaps = 10) +
  scale_color_manual(values = c("5%" = "#E41A1C", "10%" = "#377EB8")) +
  theme_bw() + 
  labs(
    title = "Cross-study Pairwise Overlap Comparison",
    x = "PreLect Pairwise Overlap",
    y = "ParSlet Pairwise Overlap"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )
print(p)

ggsave("../figs/FigR4_CrossStudy_Stability.pdf", p, width=7, height=6)
