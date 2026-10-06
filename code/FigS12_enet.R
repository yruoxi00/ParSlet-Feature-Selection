library(glmnet)
library(ggplot2)
library(dplyr)
library(ggpubr)

setwd("/Users/yrx/Downloads/FeatureSelection/feature_selection/code")

enet_scores <- function(train_data,
                        yname = "status",
                        alpha_en = 0.5,
                        nfolds = 5,
                        s_choice = c("lambda.1se", "lambda.min"),
                        lambda_scale = 0.5,
                        nlambda = 200,
                        lambda_min_ratio = 1e-4) {
  s_choice <- match.arg(s_choice)
  
  y <- factor(train_data[[yname]])
  X <- train_data[, setdiff(colnames(train_data), yname), drop = FALSE]
  Xmm <- model.matrix(~ . - 1, data = X)
  y01 <- as.numeric(y == levels(y)[2])
  
  cvfit <- cv.glmnet(
    x = Xmm, y = y01,
    family = "binomial",
    alpha = alpha_en,
    nfolds = nfolds,
    nlambda = nlambda,
    lambda.min.ratio = lambda_min_ratio
  )
  
  lam0 <- if (s_choice == "lambda.1se") cvfit$lambda.1se else cvfit$lambda.min
  
  lam <- lambda_scale * lam0
  
  b <- as.matrix(coef(cvfit, s = lam))
  b <- b[rownames(b) != "(Intercept)", , drop = FALSE]
  
  sc <- abs(b[, 1])
  sc <- sc[sc > 0]
  sc
}

make_feat_df <- function(scores, ff, frac, rea,
                         method = "ENET",
                         alpha_en = NA_real_) {
  scores <- scores[!is.na(scores)]
  data.frame(
    Feature = names(scores),
    score   = as.numeric(scores),
    rank    = rank(-as.numeric(scores), ties.method = "average"),
    ff = ff,
    frac = frac,
    rea = rea,
    method = method,
    alpha_en = alpha_en,
    stringsAsFactors = FALSE
  )
}

run_enet_pipeline_one_effect <- function(effect_size,
                                         f_fixed = 2,
                                         ffs = c(0.05, 0.1, 0.2),
                                         reas = 1:10,
                                         fracs = seq(0.1, 1, length.out = 10),
                                         alpha_grid = c(0.25, 0.5, 0.75),
                                         nfolds = 5,
                                         s_choice = c("lambda.1se", "lambda.min"),
                                         lambda_scale = 0.5,
                                         out_prefix = NULL) {
  
  s_choice <- match.arg(s_choice)
  features <- NULL
  
  for (alpha_en in alpha_grid) {
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
          
          sc <- enet_scores(
            train_data,
            yname = "status",
            alpha_en = alpha_en,
            nfolds = nfolds,
            s_choice = s_choice,
            lambda_scale = lambda_scale
          )
          
          if (length(sc) > 0) {
            features <- rbind(
              features,
              make_feat_df(sc, ff = ff, frac = frac, rea = rea,
                           method = "ENET", alpha_en = alpha_en)
            )
          }
        }
      }
    }
  }
  
  if (is.null(out_prefix)) out_prefix <- paste0("Simulated", effect_size, "_enet")
  
  write.csv(features, file = paste0(out_prefix, "_features_enet.csv"), row.names = FALSE)
  
  invisible(features)
}

run_enet_pipeline_one_effect(effect_size = "0.01", out_prefix = "../data/Simulated0.01_enet",
                             alpha_grid = c(0.25, 0.5, 0.75),
                             s_choice = "lambda.min", 
                             lambda_scale = 0.5)

run_enet_pipeline_one_effect(effect_size = "0.05", out_prefix = "../data/Simulated0.05_enet",
                             alpha_grid = c(0.25, 0.5, 0.75),
                             s_choice = "lambda.min", 
                             lambda_scale = 0.5)

run_enet_pipeline_one_effect(effect_size = "0.2",  out_prefix = "../data/Simulated0.2_enet",
                             alpha_grid = c(0.25, 0.5, 0.75),
                             s_choice = "lambda.min", 
                             lambda_scale = 0.5)

calc_prev_all <- function(dat1_full, yname="status") {
  feat_cols <- setdiff(colnames(dat1_full), yname)
  mat <- dat1_full[, feat_cols, drop = FALSE]
  colSums(mat > 0, na.rm = TRUE) / nrow(mat)
}

calc_overlap_curve_enet_two_lines <- function(features_df,
                                              dat1_full,
                                              ffs,
                                              alpha_keep = 0.5,
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
  
  for (rea_val in 1:nrea) {
    
    sub_rea <- features_df %>%
      filter(.data$rea == rea_val, .data$ff == ffs, .data$alpha_en == alpha_keep)
    
    for (ii in 1:ngroup) {
      for (jj in 1:ngroup) {
        
        fi <- fracs[ii]
        fj <- fracs[jj]
        
        si <- sub_rea %>% filter(.data$frac == fi) %>% select(Feature, score)
        sj <- sub_rea %>% filter(.data$frac == fj) %>% select(Feature, score)
        
        v1_o <- setNames(rep(0, length(all_feats)), all_feats)
        v2_o <- setNames(rep(0, length(all_feats)), all_feats)
        
        if (nrow(si) > 0) {
          keep <- intersect(si$Feature, all_feats)
          v1_o[keep] <- si$score[match(keep, si$Feature)]
        }
        if (nrow(sj) > 0) {
          keep <- intersect(sj$Feature, all_feats)
          v2_o[keep] <- sj$score[match(keep, sj$Feature)]
        }
        
        prev_all <- calc_prev_all(dat1_full, yname=yname)
        
        p1 <- setNames(prev_all, names(prev_all))
        p2 <- p1
        p1 <- setNames(rep(0, length(all_feats)), all_feats); p1[names(prev_all)] <- prev_all
        p2 <- p1
        
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
    types   = factor(rep(c("ENET", "ParSlet"), each = length(fracs)),
                     levels = c("ENET", "ParSlet")),
    stringsAsFactors = FALSE
  )
}

plot_one_effect_enet_two_lines <- function(effect_size,
                                           enet_csv,
                                           stool_prefix = "../data/Stool",
                                           alpha_keep = 0.5,
                                           prev_power = 2.5,
                                           panel_labels = c("a","b","c"),
                                           ylab_left = "Frequency") {
  
  features <- read.csv(enet_csv, stringsAsFactors = FALSE)
  if (!("alpha_en" %in% names(features))) stop("ENET feature file must contain column: alpha_en")
  features$alpha_en <- as.numeric(features$alpha_en)
  features$Feature  <- as.character(features$Feature)
  features$frac     <- as.numeric(features$frac)
  features$ff       <- as.numeric(features$ff)
  features$rea      <- as.integer(features$rea)
  
  col_vals <- c(
    ENET    = "#DD5129FF",
    ParSlet = "#0F7BA2FF"
  )
  
  shape_vals <- c(
    ENET    = 16,
    ParSlet = 17
  )
  
  pdf_list <- list()
  idx <- 1
  
  for (ffs in c(0.05, 0.1, 0.2)) {
    
    dat1_sub <- read.csv(
      file = paste0(stool_prefix, effect_size, "_f_", 2, "_", ffs, ".csv"),
      row.names = 1, header = TRUE
    )
    
    dat1 <- calc_overlap_curve_enet_two_lines(
      features_df = features,
      dat1_full   = dat1_sub,
      ffs         = ffs,
      alpha_keep  = alpha_keep,
      prev_power  = prev_power,
      top_prop    = 0.1,
      nrea        = 10
    )
    
    p_one <- ggplot(dat1, aes(sample_fraction, Overlap)) +
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
      p_one <- p_one + ylab(ylab_left)
    } else {
      p_one <- p_one + ylab(NULL)
    }
    
    pdf_list[[idx]] <- p_one
    idx <- idx + 1
  }
  
  ggarrange(
    plotlist = pdf_list, ncol = 3, nrow = 1, align = "hv",
    labels = panel_labels
  )
}

alpha_keep <- 0.5
prev_power <- 2.5
type_i_name <- paste0("ENET_i")
col_vals <- setNames(c("#DD5129FF", "#0F7BA2FF"), c("ENET_o", type_i_name))
shape_vals <- setNames(c(16, 17), c("ENET_o", type_i_name))

p1 <- plot_one_effect_enet_two_lines(
  "0.01", "../data/Simulated0.01_enet_features_enet.csv",
  stool_prefix = "../data/Stool",
  alpha_keep = alpha_keep, prev_power = prev_power,
  panel_labels = c("a","b","c"),
  ylab_left = "Frequency"
)

p2 <- plot_one_effect_enet_two_lines(
  "0.05", "../data/Simulated0.05_enet_features_enet.csv",
  stool_prefix = "../data/Stool",
  alpha_keep = alpha_keep, prev_power = prev_power,
  panel_labels = c("d","e","f"),
  ylab_left = "Frequency"
)

p3 <- plot_one_effect_enet_two_lines(
  "0.2", "../data/Simulated0.2_enet_features_enet.csv",
  stool_prefix = "../data/Stool",
  alpha_keep = alpha_keep, prev_power = prev_power,
  panel_labels = c("g","h","i"),
  ylab_left = "Frequency"
)

row1 <- annotate_figure(p1, left = text_grob("Effect size = 0.01", rot = 90, size = 10))
row2 <- annotate_figure(p2, left = text_grob("Effect size = 0.05", rot = 90, size = 10))
row3 <- annotate_figure(p3, left = text_grob("Effect size = 0.20", rot = 90, size = 10))
p_body <- ggarrange(row1, row2, row3, ncol = 1, nrow = 3, align = "hv")

col_vals <- c(
  ENET    = "#DD5129FF",
  ParSlet = "#0F7BA2FF"
)

shape_vals <- c(
  ENET    = 16,
  ParSlet = 17
)

legend_plot <- ggplot(
  data.frame(
    x = 1, y = 1,
    types = factor(c("ENET", "ParSlet"), levels = c("ENET", "ParSlet"))
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

ggsave(final_plot, file = "../figs/Simulated_feature_enet.pdf", width = 7.2, height = 7.7)
