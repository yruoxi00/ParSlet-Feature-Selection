library(mRMRe)
library(ggplot2)
library(dplyr)
library(ggpubr)

setwd("/Users/yrx/Downloads/FeatureSelection/feature_selection/code")

mrmr_scores <- function(train_data, yname = "status", K = 200) {
  y <- factor(train_data[[yname]])
  
  X <- train_data[, setdiff(colnames(train_data), yname), drop = FALSE]
  X <- as.data.frame(lapply(X, function(z) log1p(as.numeric(z))))  # log1p
  
  y01 <- as.numeric(y == levels(y)[2])
  
  df <- cbind(y = y01, X)
  df <- as.data.frame(df)
  
  data_obj <- mRMRe::mRMR.data(data = df)
  
  p <- ncol(X)
  K <- min(K, p)
  
  res <- mRMRe::mRMR.classic(data = data_obj, target_indices = 1, feature_count = K)
  
  feat_idx <- mRMRe::solutions(res)[[1]]
  feat_idx_X <- feat_idx - 1
  feat_idx_X <- feat_idx_X[feat_idx_X >= 1 & feat_idx_X <= p]
  
  feat_names <- colnames(X)[feat_idx_X]
  
  ranks <- seq_along(feat_names)
  scores <- (length(feat_names) - ranks + 1)
  names(scores) <- feat_names
  
  scores
}

make_feat_df <- function(scores, ff, frac, rea, method = "mRMR") {
  scores <- scores[!is.na(scores)]
  data.frame(
    Feature = names(scores),
    score   = as.numeric(scores),
    rank    = rank(-as.numeric(scores), ties.method = "average"),
    ff = ff,
    frac = frac,
    rea = rea,
    method = method,
    stringsAsFactors = FALSE
  )
}

run_mrmr_pipeline_one_effect <- function(effect_size,
                                         f_fixed = 2,
                                         ffs = c(0.05, 0.1, 0.2),
                                         reas = 1:10,
                                         fracs = seq(0.1, 1, length.out = 10),
                                         K = 200,
                                         out_prefix = NULL) {
  features <- NULL
  
  for (rea in reas) {
    for (frac in fracs) {
      for (ff in ffs) {
        
        set.seed(rea)
        
        dat_path <- paste0("../data/Stool", effect_size, "_f_", f_fixed, "_", ff, ".csv")
        dat1_sub <- read.csv(dat_path, row.names = 1, header = TRUE)
        dat1_sub$status <- factor(dat1_sub$status)
        
        n_keep <- max(2, floor(frac * nrow(dat1_sub)))
        dat1_sub <- dat1_sub[sample(1:nrow(dat1_sub), size = n_keep), , drop = FALSE]
        
        n_train <- max(2, floor(0.8 * nrow(dat1_sub)))
        train_indices <- sample(1:nrow(dat1_sub), size = n_train)
        train_data <- dat1_sub[train_indices, , drop = FALSE]
        
        sc <- mrmr_scores(train_data, yname = "status", K = K)
        
        if (length(sc) > 0) {
          features <- rbind(features, make_feat_df(sc, ff = ff, frac = frac, rea = rea, method = "mRMR"))
        }
      }
    }
  }
  
  if (is.null(out_prefix)) out_prefix <- paste0("Simulated", effect_size)
  write.csv(features, file = paste0(out_prefix, "_features_mrmr.csv"), row.names = FALSE)
  invisible(features)
}

run_mrmr_pipeline_one_effect("0.01", out_prefix = "../data/Simulated0.01_mrmr", K = 200)
run_mrmr_pipeline_one_effect("0.05", out_prefix = "../data/Simulated0.05_mrmr", K = 200)
run_mrmr_pipeline_one_effect("0.2",  out_prefix = "../data/Simulated0.2_mrmr",  K = 200)

calc_prev_for_frac <- function(dat1_full, frac, rea_val, yname = "status") {
  stopifnot(frac > 0 && frac <= 1)
  
  set.seed(rea_val)
  
  feat_cols <- setdiff(colnames(dat1_full), yname)
  p <- length(feat_cols)
  
  n_keep_feat <- max(2, floor(frac * p))
  n_keep_feat <- min(n_keep_feat, p)
  
  keep_feats <- sample(feat_cols, size = n_keep_feat)
  
  mat <- dat1_full[, keep_feats, drop = FALSE]
  prev <- colSums(mat > 0, na.rm = TRUE) / nrow(mat)
  prev
}

calc_overlap_curve_mrmr_two_lines <- function(features_df,
                                              dat1_full,
                                              ffs,
                                              prev_power = 2.5,
                                              top_prop = 0.1,
                                              nrea = 10,
                                              yname = "status") {
  fracs <- sort(unique(features_df$frac))
  ngroup <- length(fracs)
  
  all_feats <- setdiff(colnames(dat1_full), yname)
  top_n <- max(1, floor(top_prop * length(all_feats)))
  
  overlap_o <- matrix(0, nrow = ngroup, ncol = ngroup)
  overlap_i <- matrix(0, nrow = ngroup, ncol = ngroup)
  
  for (rea_id in 1:nrea) {
    
    sub_rea <- features_df %>% filter(.data$rea == rea_id, .data$ff == ffs)
    
    for (ii in 1:ngroup) {
      for (jj in 1:ngroup) {
        
        fi <- fracs[ii]
        fj <- fracs[jj]
        
        si <- sub_rea %>% filter(.data$frac == fi) %>% select(Feature, score)
        sj <- sub_rea %>% filter(.data$frac == fj) %>% select(Feature, score)
        
        v1_o <- setNames(rep(0, length(all_feats)), all_feats)
        v2_o <- setNames(rep(0, length(all_feats)), all_feats)
        
        if (nrow(si) > 0) {
          hit <- intersect(si$Feature, all_feats)
          v1_o[hit] <- si$score[match(hit, si$Feature)]
        }
        if (nrow(sj) > 0) {
          hit <- intersect(sj$Feature, all_feats)
          v2_o[hit] <- sj$score[match(hit, sj$Feature)]
        }
        
        prev1 <- calc_prev_for_frac(dat1_full, frac = fi, rea_val = rea_id, yname = yname)
        prev2 <- calc_prev_for_frac(dat1_full, frac = fj, rea_val = rea_id, yname = yname)
        
        p1 <- setNames(rep(0, length(all_feats)), all_feats)
        p2 <- setNames(rep(0, length(all_feats)), all_feats)
        p1[names(prev1)] <- prev1
        p2[names(prev2)] <- prev2
        
        v1_i <- v1_o * (p1 ^ prev_power)
        v2_i <- v2_o * (p2 ^ prev_power)
        
        ord1_o <- order(v1_o, decreasing = TRUE)
        ord2_o <- order(v2_o, decreasing = TRUE)
        ord1_i <- order(v1_i, decreasing = TRUE)
        ord2_i <- order(v2_i, decreasing = TRUE)
        
        overlap_o[ii, jj] <- overlap_o[ii, jj] +
          length(intersect(ord1_o[1:top_n], ord2_o[1:top_n])) / top_n
        
        overlap_i[ii, jj] <- overlap_i[ii, jj] +
          length(intersect(ord1_i[1:top_n], ord2_i[1:top_n])) / top_n
      }
    }
  }
  
  overlap_o <- overlap_o / nrea
  overlap_i <- overlap_i / nrea
  diag(overlap_o) <- NA
  diag(overlap_i) <- NA
  
  means_o <- apply(overlap_o, 1, function(x) mean(x, na.rm = TRUE))
  means_i <- apply(overlap_i, 1, function(x) mean(x, na.rm = TRUE))
  sd_o    <- apply(overlap_o, 1, function(x) sd(x, na.rm = TRUE))
  sd_i    <- apply(overlap_i, 1, function(x) sd(x, na.rm = TRUE))
  
  data.frame(
    sample_fraction = rep(as.numeric(fracs), 2),
    Overlap = c(as.numeric(means_o), as.numeric(means_i)),
    sds     = c(as.numeric(sd_o),    as.numeric(sd_i)),
    types   = factor(rep(c("mRMR", "ParSlet"), each = length(fracs)),
                     levels = c("mRMR", "ParSlet")),
    stringsAsFactors = FALSE
  )
}

plot_one_effect_mrmr_two_lines <- function(effect_size,
                                           mrmr_csv,
                                           stool_prefix = "../data/Stool",
                                           prev_power = 2.5,
                                           panel_labels = c("a","b","c"),
                                           ylab_left = "Frequency") {
  
  features <- read.csv(mrmr_csv, stringsAsFactors = FALSE)
  features$Feature <- as.character(features$Feature)
  features$frac <- as.numeric(features$frac)
  features$ff   <- as.numeric(features$ff)
  features$rea  <- as.integer(features$rea)
  
  col_vals <- c(
    mRMR    = "#DD5129FF",
    ParSlet = "#0F7BA2FF"
  )
  
  shape_vals <- c(
    mRMR    = 16,
    ParSlet = 17
  )
  
  pdf_list <- list()
  idx <- 1
  
  for (ffs in c(0.05, 0.1, 0.2)) {
    
    dat1_sub <- read.csv(
      paste0(stool_prefix, effect_size, "_f_", 2, "_", ffs, ".csv"),
      row.names = 1, header = TRUE
    )
    
    dat1 <- calc_overlap_curve_mrmr_two_lines(
      features_df = features,
      dat1_full   = dat1_sub,
      ffs         = ffs,
      prev_power  = prev_power,
      top_prop    = 0.1,
      nrea        = 10,
      yname       = "status"
    )
    
    p_one <- ggplot(dat1, aes(sample_fraction, Overlap)) +
      geom_line(aes(group = types, color = types),
                linetype = "dotdash", linewidth = 0.4) +
      geom_point(aes(color = types, shape = types), size = 2.5) +
      geom_errorbar(aes(ymin = Overlap - sds,
                        ymax = Overlap + sds,
                        color = types),
                    width = 0.02, linewidth = 0.3) +
      scale_color_manual(values = col_vals) +
      scale_shape_manual(values = shape_vals) +
      theme_bw() +
      labs(title = paste("Fraction of spiking-in:", ffs), x = NULL, y = NULL) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.ticks  = element_line(linewidth = 0.2),
        legend.position = "none"
      )
    
    if (idx == 1) {
      p_one <- p_one + ylab(ylab_left)
    } else {
      p_one <- p_one + ylab(NULL)
    }
    
    pdf_list[[idx]] <- p_one
    idx <- idx + 1
  }
  
  ggarrange(
    plotlist = pdf_list,
    ncol = 3, nrow = 1,
    align = "hv",
    labels = panel_labels
  )
}

prev_power <- 2.5
type_i_name <- paste0("mRMR_i")
col_vals <- setNames(c("#DD5129FF", "#0F7BA2FF"), c("mRMR_o", type_i_name))
shape_vals <- setNames(c(16, 17), c("mRMR_o", type_i_name))

p1 <- plot_one_effect_mrmr_two_lines(
  "0.01", "../data/Simulated0.01_mrmr_features_mrmr.csv",
  stool_prefix = "../data/Stool",
  prev_power = prev_power,
  panel_labels = c("a","b","c"),
  ylab_left = "Frequency"
)

p2 <- plot_one_effect_mrmr_two_lines(
  "0.05", "../data/Simulated0.05_mrmr_features_mrmr.csv",
  stool_prefix = "../data/Stool",
  prev_power = prev_power,
  panel_labels = c("d","e","f"),
  ylab_left = "Frequency"
)

p3 <- plot_one_effect_mrmr_two_lines(
  "0.2", "../data/Simulated0.2_mrmr_features_mrmr.csv",
  stool_prefix = "../data/Stool",
  prev_power = prev_power,
  panel_labels = c("g","h","i"),
  ylab_left = "Frequency"
)

row1 <- annotate_figure(p1, left = text_grob("Effect size = 0.01", rot = 90, size = 10))
row2 <- annotate_figure(p2, left = text_grob("Effect size = 0.05", rot = 90, size = 10))
row3 <- annotate_figure(p3, left = text_grob("Effect size = 0.20", rot = 90, size = 10))

p_body <- ggarrange(row1, row2, row3, ncol = 1, nrow = 3, align = "hv")

col_vals <- c(
  mRMR    = "#DD5129FF",
  ParSlet = "#0F7BA2FF"
)

shape_vals <- c(
  mRMR    = 16,
  ParSlet = 17
)

legend_plot <- ggplot(
  data.frame(
    x = 1, y = 1,
    types = factor(c("mRMR", "ParSlet"), levels = c("mRMR", "ParSlet"))
  ),
  aes(x = x, y = y, color = types, shape = types)
) +
  geom_point(size = 3) +
  scale_color_manual(values = col_vals) +
  scale_shape_manual(values = shape_vals) +
  guides(
    color = guide_legend(ncol = 2, override.aes = list(size = 4)),
    shape = guide_legend(ncol = 2)
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  )

legend <- get_legend(legend_plot)

final_plot <- ggarrange(p_body, legend, ncol = 1, heights = c(1, 0.1))
final_plot

ggsave(final_plot, file = "../figs/Simulated_feature_mrmr.pdf", width = 7.2, height = 7.7)
