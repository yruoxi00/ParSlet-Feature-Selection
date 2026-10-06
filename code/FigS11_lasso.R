library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(glmnet)
library(cowplot)

setwd("/Users/yrx/Downloads/FeatureSelection/feature_selection/code")

lasso_scores <- function(train_data, yname = "status", nfolds = 5,
                         base = c("lambda.1se","lambda.min"),
                         lambda_scale = 0.7) {
  base <- match.arg(base)
  
  y <- factor(train_data[[yname]])
  X <- train_data[, setdiff(colnames(train_data), yname), drop = FALSE]
  Xmm <- model.matrix(~ . - 1, data = X)
  y01 <- as.numeric(y == levels(y)[2])
  
  cvfit <- cv.glmnet(
    x = Xmm, y = y01,
    family = "binomial",
    alpha = 1,
    nfolds = nfolds
  )
  
  lam0 <- if (base == "lambda.1se") cvfit$lambda.1se else cvfit$lambda.min
  lam  <- lambda_scale * lam0
  
  b <- as.matrix(coef(cvfit, s = lam))
  b <- b[rownames(b) != "(Intercept)", , drop = FALSE]
  
  sc <- abs(b[, 1])
  sc <- sc[sc > 0]
  sc
}

make_feat_df <- function(scores, ff, frac, rea, method = "LASSO") {
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

run_lasso_pipeline_one_effect <- function(effect_size,
                                          f_fixed = 2,
                                          ffs = c(0.05, 0.1, 0.2),
                                          reas = 1:10,
                                          fracs = seq(0.1, 1, length.out = 10),
                                          out_prefix = NULL) {
  features <- NULL
  
  for (rea in reas) {
    for (frac in fracs) {
      for (ff in ffs) {
        
        set.seed(rea)
        
        dat_path  <- paste0("../data/Stool", effect_size, "_f_", f_fixed, "_", ff, ".csv")
        meta_path <- paste0("../data/Spike_Stool", effect_size, "__f_", f_fixed, "_", ff, ".csv")
        
        dat1_sub  <- read.csv(dat_path, row.names = 1, header = TRUE)
        dat1_sub$status <- factor(dat1_sub$status)
        
        n_keep <- max(2, floor(frac * nrow(dat1_sub))) 
        dat1_sub <- dat1_sub[sample(1:nrow(dat1_sub), size = n_keep), , drop = FALSE]
        
        n_train <- max(2, floor(0.8 * nrow(dat1_sub)))
        train_indices <- sample(1:nrow(dat1_sub), size = n_train)
        train_data <- dat1_sub[train_indices, , drop = FALSE]
        
        sc <- lasso_scores(train_data, yname = "status", nfolds = 5,
                           base = "lambda.min", lambda_scale = 0.7)
        
        if (length(sc) > 0) {
          features <- rbind(features, make_feat_df(sc, ff = ff, frac = frac, rea = rea, method = "LASSO"))
        } else {
        }
      }
    }
  }
  
  if (is.null(out_prefix)) out_prefix <- paste0("Simulated", effect_size)
  write.csv(features, file = paste0(out_prefix, "_features_lasso.csv"), row.names = FALSE)
  
  invisible(features)
}

run_lasso_pipeline_one_effect(effect_size = "0.01", out_prefix = "../data/Simulated0.01_lasso")
run_lasso_pipeline_one_effect(effect_size = "0.05", out_prefix = "../data/Simulated0.05_lasso")
run_lasso_pipeline_one_effect(effect_size = "0.2",  out_prefix = "../data/Simulated0.2_lasso")

calc_prev_all <- function(dat1_full, yname="status") {
  feat_cols <- setdiff(colnames(dat1_full), yname)
  mat <- dat1_full[, feat_cols, drop=FALSE]
  colSums(mat > 0, na.rm=TRUE) / nrow(mat)
}

calc_overlap_curve_lasso_two_lines <- function(features_df,
                                               dat1_full,
                                               ffs,
                                               fracs_global,
                                               nrea = 10,
                                               top_prop = 0.1,
                                               power_alpha = 2.5,
                                               yname = "status") {
  all_feats <- setdiff(colnames(dat1_full), yname)
  p <- length(all_feats)
  top_n <- max(1, floor(top_prop * p))
  
  ngroup <- length(fracs_global)
  overlap_o <- matrix(0, nrow = ngroup, ncol = ngroup)
  overlap_i <- matrix(0, nrow = ngroup, ncol = ngroup)
  
  for (rea in 1:nrea) {
    
    sub_rea <- features_df %>%
      filter(.data$rea == rea, .data$ff == ffs)
    
    for (ii in 1:ngroup) {
      for (jj in 1:ngroup) {
        
        fi <- fracs_global[ii]
        fj <- fracs_global[jj]
        
        si <- sub_rea %>% filter(.data$frac == fi) %>% select(Feature, score)
        sj <- sub_rea %>% filter(.data$frac == fj) %>% select(Feature, score)
        
        v1_o <- setNames(rep(0, p), all_feats)
        v2_o <- setNames(rep(0, p), all_feats)
        
        if (nrow(si) > 0) {
          hit <- intersect(si$Feature, all_feats)
          v1_o[hit] <- si$score[match(hit, si$Feature)]
        }
        if (nrow(sj) > 0) {
          hit <- intersect(sj$Feature, all_feats)
          v2_o[hit] <- sj$score[match(hit, sj$Feature)]
        }
        
        prev_all <- calc_prev_all(dat1_full, yname=yname)
        
        p1 <- setNames(rep(0, p), all_feats)
        p2 <- setNames(rep(0, p), all_feats)
        p1[names(prev_all)] <- prev_all
        p2[names(prev_all)] <- prev_all
        
        v1_i <- v1_o * (p1 ^ power_alpha)
        v2_i <- v2_o * (p2 ^ power_alpha)
        
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
    sample_fraction = rep(as.numeric(fracs_global), 2),
    Overlap = c(as.numeric(means_o), as.numeric(means_i)),
    sds     = c(as.numeric(sd_o),    as.numeric(sd_i)),
    types   = factor(rep(c("LASSO", "ParSlet"), each = length(fracs_global)),
                     levels = c("LASSO", "ParSlet")),
    stringsAsFactors = FALSE
  )
}

plot_one_effect_lasso_like_rf <- function(effect_size,
                                          lasso_feat_csv,
                                          stool_prefix = "../data/Stool",
                                          f_fixed = 2,
                                          ffs_list = c(0.05, 0.1, 0.2),
                                          nrea = 10,
                                          top_prop = 0.1,
                                          power_alpha = 2.5,
                                          panel_labels = c("a","b","c")) {
  
  features <- read.csv(lasso_feat_csv, stringsAsFactors = FALSE)
  features$Feature <- as.character(features$Feature)
  features$frac <- as.numeric(features$frac)
  features$ff   <- as.numeric(features$ff)
  features$rea  <- as.integer(features$rea)
  
  fracs_global <- sort(unique(features$frac))
  all_panels <- list()
  
  col_vals <- c(
    LASSO   = "#DD5129FF",
    ParSlet = "#0F7BA2FF"
  )
  
  shape_vals <- c(
    LASSO   = 16,
    ParSlet = 17
  )
  
  idx <- 1
  for (ffs in ffs_list) {
    
    dat1_sub <- read.csv(
      file = paste0(stool_prefix, effect_size, "_f_", f_fixed, "_", ffs, ".csv"),
      row.names = 1, header = TRUE
    )
    
    dat1 <- calc_overlap_curve_lasso_two_lines(
      features_df   = features,
      dat1_full     = dat1_sub,
      ffs           = ffs,
      fracs_global  = fracs_global,
      nrea          = nrea,
      top_prop      = top_prop,
      power_alpha   = power_alpha,
      yname         = "status"
    )
    
    all_panels[[idx]] <- ggplot(dat1, aes(sample_fraction, Overlap)) +
      geom_line(aes(group = types, color = types), linetype = "dotdash", linewidth = 0.4) +
      geom_point(aes(color = types, shape = types), size = 2.5) +
      geom_errorbar(aes(ymin = Overlap - sds, ymax = Overlap + sds, color = types),
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
      all_panels[[idx]] <- all_panels[[idx]] + ylab("Frequency")
    } else {
      all_panels[[idx]] <- all_panels[[idx]] + ylab(NULL)
    }
    
    idx <- idx + 1
  }
  
  ggarrange(
    plotlist = all_panels,
    ncol = 3, nrow = 1, align = "hv",
    labels = panel_labels
  )
}

p1 <- plot_one_effect_lasso_like_rf(
  effect_size = "0.01",
  lasso_feat_csv = "../data/Simulated0.01_lasso_features_lasso.csv",
  power_alpha = 2.5,
  panel_labels = c("a","b","c")
)

p2 <- plot_one_effect_lasso_like_rf(
  effect_size = "0.05",
  lasso_feat_csv = "../data/Simulated0.05_lasso_features_lasso.csv",
  power_alpha = 2.5,
  panel_labels = c("d","e","f")
)

p3 <- plot_one_effect_lasso_like_rf(
  effect_size = "0.2",
  lasso_feat_csv = "../data/Simulated0.2_lasso_features_lasso.csv",
  power_alpha = 2.5,
  panel_labels = c("g","h","i")
)

row1 <- annotate_figure(p1, left = text_grob("Effect size = 0.01", rot = 90, size = 10))
row2 <- annotate_figure(p2, left = text_grob("Effect size = 0.05", rot = 90, size = 10))
row3 <- annotate_figure(p3, left = text_grob("Effect size = 0.20", rot = 90, size = 10))

p_body <- ggarrange(row1, row2, row3, ncol = 1, nrow = 3, align = "hv")

col_vals <- c(
  LASSO   = "#DD5129FF",
  ParSlet = "#0F7BA2FF"
)

shape_vals <- c(
  LASSO   = 16,
  ParSlet = 17
)

legend_plot <- ggplot(
  data.frame(x = 1, y = 1,
             types = factor(c("LASSO", "ParSlet"),
                            levels = c("LASSO", "ParSlet"))),
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

ggsave(final_plot, file = "../figs/Simulated_feature_lasso.pdf", width = 7.2, height = 7.7)
