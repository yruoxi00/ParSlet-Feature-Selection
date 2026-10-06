library(PreLectR)
library(dplyr)
library(ggplot2)
library(ggpubr)

setwd("/Users/yrx/Downloads/FeatureSelection/feature_selection/code")

features <- read.csv("../results/Simulated0.01_features.csv")

ffs_list <- c(0.05, 0.1, 0.2)
fracs    <- sort(unique(features$frac))
ngroup   <- length(fracs)

n_rea    <- 10
top_prop <- 0.10
gamma    <- 2.5

lambda_steps    <- 50
inner_spl_ratio <- 0.7
n_cores         <- 8
max_iter_fit    <- 10000

drop_status <- function(dat){
  last <- dat[[ncol(dat)]]
  if (!is.numeric(last) || tolower(colnames(dat)[ncol(dat)]) %in% c("status","label","y","class")) {
    dat[, -ncol(dat), drop=FALSE]
  } else {
    dat
  }
}

top_our_gini <- function(dat_sf, features_sub, frac_val, seed, gamma, top_n){
  set.seed(seed)
  dat_x <- drop_status(dat_sf)
  n     <- nrow(dat_x)
  n_sub <- max(2, floor(frac_val * n))
  train_table <- dat_x[sample(seq_len(n), size=n_sub), , drop=FALSE]
  
  richness <- colSums(train_table > 0) / nrow(train_table)
  names(richness) <- colnames(train_table)
  
  sub <- features_sub %>% filter(frac == frac_val)
  r <- as.numeric(richness[sub$Feature])
  r[is.na(r)] <- 0
  
  score <- sub$MeanDecreaseGini * (r ^ gamma)
  ord <- order(score, decreasing = TRUE)
  sub$Feature[ord][1:top_n]
}

top_prelect <- function(dat_sf, frac_val, seed,
                        lambda_steps, inner_spl_ratio, n_cores, max_iter_fit,
                        top_n){
  set.seed(seed)
  
  y <- factor(dat_sf[[ncol(dat_sf)]])
  if (length(levels(y)) != 2) stop("PreLect classification requires y with 2 levels.")
  
  dat_x <- drop_status(dat_sf)
  n     <- nrow(dat_x)
  n_sub <- max(2, floor(frac_val * n))
  idx   <- sample(seq_len(n), size=n_sub)
  
  X_sf  <- dat_x[idx, , drop=FALSE]
  y_sub <- y[idx]
  
  X_raw <- t(as.matrix(X_sf))
  storage.mode(X_raw) <- "numeric"
  X_scaled <- t(scale(t(log1p(X_raw))))
  X_scaled[is.na(X_scaled)] <- 0
  
  lrange <- AutoScanning(X_scaled, X_raw, y_sub, task="classification", step=lambda_steps)
  out_dir_tmp <- tempdir()
  
  tuning_res <- LambdaTuningParallel(
    X_scaled, X_raw, y_sub, lrange,
    n_cores = n_cores,
    outpath = out_dir_tmp,
    spl_ratio = inner_spl_ratio,
    task = "classification"
  )
  
  pick <- LambdaDecision(tuning_res$TuningResult, tuning_res$PvlDistSummary,
                         maxdepth=5, minbucket=3)
  opt_lmbd <- pick$opt_lmbd
  
  prev <- GetPrevalence(X_raw)
  pre_out <- PreLect(X_scaled, prev, y_sub,
                     lambda = opt_lmbd, task="classification",
                     max_iter = max_iter_fit)
  
  featprop <- FeatureProperty(X_raw, y_sub, pre_out, task="classification")
  
  score <- abs(featprop$coef)
  ord <- order(score, decreasing = TRUE)
  featprop$FeatName[ord][1:top_n]
}

pdf_list <- list()
index <- 1

for (ffs in ffs_list){
  
  dat1_sub <- read.csv(
    paste("../data/Stool0.01_f_", 2, "_", ffs, ".csv", sep=""),
    row.names=1, header=TRUE
  )
  
  overlap_ours <- matrix(0, nrow=ngroup, ncol=ngroup)
  overlap_pl   <- matrix(0, nrow=ngroup, ncol=ngroup)
  
  for (rea_id in 1:n_rea){
    
    features_sub <- features %>% filter(rea == rea_id, ff == ffs)
    
    any_frac <- fracs[1]
    n_feat <- sum(features_sub$frac == any_frac)
    top_n <- max(1, floor(top_prop * n_feat))
    
    ours_top_list <- vector("list", ngroup)
    pl_top_list   <- vector("list", ngroup)
    
    for (ii in seq_along(fracs)){
      fr <- fracs[ii]
      ours_top_list[[ii]] <- top_our_gini(dat1_sub, features_sub, fr, seed=rea_id, gamma=gamma, top_n=top_n)
      pl_top_list[[ii]] <- top_prelect(dat1_sub, fr, seed=rea_id,
                                       lambda_steps=lambda_steps,
                                       inner_spl_ratio=inner_spl_ratio,
                                       n_cores=n_cores,
                                       max_iter_fit=max_iter_fit,
                                       top_n=top_n)
    }
    
    for (i in seq_along(fracs)){
      for (j in seq_along(fracs)){
        overlap_ours[i,j] <- overlap_ours[i,j] + length(intersect(ours_top_list[[i]], ours_top_list[[j]])) / top_n
        overlap_pl[i,j]   <- overlap_pl[i,j]   + length(intersect(pl_top_list[[i]],   pl_top_list[[j]]))   / top_n
      }
    }
  }
  
  overlap_ours <- overlap_ours / n_rea
  overlap_pl   <- overlap_pl   / n_rea
  
  diag(overlap_ours) <- NA
  diag(overlap_pl)   <- NA
  
  rownames(overlap_ours) <- colnames(overlap_ours) <- fracs
  rownames(overlap_pl)   <- colnames(overlap_pl)   <- fracs
  
  mean_ours <- apply(overlap_ours, 1, function(x) mean(x, na.rm=TRUE))
  mean_pl   <- apply(overlap_pl,   1, function(x) mean(x, na.rm=TRUE))
  
  sd_ours <- apply(overlap_ours, 1, function(x) sd(x, na.rm=TRUE))
  sd_pl   <- apply(overlap_pl,   1, function(x) sd(x, na.rm=TRUE))
  
  dat_plot <- data.frame(
    sample_fraction = rep(as.numeric(fracs), 2),
    Overlap = c(mean_ours, mean_pl),
    se = c(sd_ours / sqrt(n_rea), sd_pl / sqrt(n_rea)),
    types = c(rep("Gini^2.5", ngroup),
              rep("PreLect", ngroup))
  )
  
  pdf_list[[index]] <- ggplot(dat_plot, aes(sample_fraction, Overlap)) +
    geom_line(aes(group=types, color=types), linetype="dotdash", size=0.4) +
    geom_point(aes(color=types, shape=types), size=2.5) +
    geom_errorbar(aes(ymin=Overlap-se, ymax=Overlap+se, color=types), width=0.02) +
    scale_color_manual(values = c("Gini^2.5" = "#0F7BA2FF", "PreLect" = "#DD5129FF")) + 
    scale_shape_manual(values = c("Gini^2.5" = 17, "PreLect" = 16)) +
    theme_bw() + xlab("") + ylab("") +
    labs(title=paste("Fraction of spiking-in:", ffs)) +
    theme(
      plot.title = element_text(hjust=0.5, size=8),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size=8, color='black'),
      axis.text.y = element_text(size=8, color='black'),
      axis.ticks  = element_line(size=0.2),
      legend.position = "none"
    )
  
  index <- index + 1
}

pdf_list[[1]] <- pdf_list[[1]] + ylab("Frequency")
p1 <- do.call(ggarrange, c(pdf_list, list(labels=c("a","b","c"), ncol=3, nrow=1), align="hv"))
p1

### 

features <- read.csv("../results/Simulated0.05_features.csv")

ffs_list <- c(0.05, 0.1, 0.2)
fracs    <- sort(unique(features$frac))
ngroup   <- length(fracs)

n_rea    <- 10
top_prop <- 0.10
gamma    <- 2.5

lambda_steps    <- 50
inner_spl_ratio <- 0.7
n_cores         <- 8
max_iter_fit    <- 10000

drop_status <- function(dat){
  last <- dat[[ncol(dat)]]
  if (!is.numeric(last) || tolower(colnames(dat)[ncol(dat)]) %in% c("status","label","y","class")) {
    dat[, -ncol(dat), drop=FALSE]
  } else {
    dat
  }
}

top_our_gini <- function(dat_sf, features_sub, frac_val, seed, gamma, top_n){
  set.seed(seed)
  dat_x <- drop_status(dat_sf)
  n     <- nrow(dat_x)
  n_sub <- max(2, floor(frac_val * n))
  train_table <- dat_x[sample(seq_len(n), size=n_sub), , drop=FALSE]
  
  richness <- colSums(train_table > 0) / nrow(train_table)
  names(richness) <- colnames(train_table)
  
  sub <- features_sub %>% filter(frac == frac_val)
  r <- as.numeric(richness[sub$Feature])
  r[is.na(r)] <- 0
  
  score <- sub$MeanDecreaseGini * (r ^ gamma)
  ord <- order(score, decreasing = TRUE)
  sub$Feature[ord][1:top_n]
}

top_prelect <- function(dat_sf, frac_val, seed,
                        lambda_steps, inner_spl_ratio, n_cores, max_iter_fit,
                        top_n){
  set.seed(seed)
  
  y <- factor(dat_sf[[ncol(dat_sf)]])
  if (length(levels(y)) != 2) stop("PreLect classification requires y with 2 levels.")
  
  dat_x <- drop_status(dat_sf)
  n     <- nrow(dat_x)
  n_sub <- max(2, floor(frac_val * n))
  idx   <- sample(seq_len(n), size=n_sub)
  
  X_sf  <- dat_x[idx, , drop=FALSE]
  y_sub <- y[idx]
  
  X_raw <- t(as.matrix(X_sf))
  storage.mode(X_raw) <- "numeric"
  X_scaled <- t(scale(t(log1p(X_raw))))
  X_scaled[is.na(X_scaled)] <- 0
  
  lrange <- AutoScanning(X_scaled, X_raw, y_sub, task="classification", step=lambda_steps)
  out_dir_tmp <- tempdir()
  
  tuning_res <- LambdaTuningParallel(
    X_scaled, X_raw, y_sub, lrange,
    n_cores = n_cores,
    outpath = out_dir_tmp,
    spl_ratio = inner_spl_ratio,
    task = "classification"
  )
  
  pick <- LambdaDecision(tuning_res$TuningResult, tuning_res$PvlDistSummary,
                         maxdepth=5, minbucket=3)
  opt_lmbd <- pick$opt_lmbd
  
  prev <- GetPrevalence(X_raw)
  pre_out <- PreLect(X_scaled, prev, y_sub,
                     lambda = opt_lmbd, task="classification",
                     max_iter = max_iter_fit)
  
  featprop <- FeatureProperty(X_raw, y_sub, pre_out, task="classification")
  
  score <- abs(featprop$coef)
  ord <- order(score, decreasing = TRUE)
  featprop$FeatName[ord][1:top_n]
}

pdf_list <- list()
index <- 1

for (ffs in ffs_list){
  
  dat1_sub <- read.csv(
    paste("../data/Stool0.05_f_", 2, "_", ffs, ".csv", sep=""),
    row.names=1, header=TRUE
  )
  
  overlap_ours <- matrix(0, nrow=ngroup, ncol=ngroup)
  overlap_pl   <- matrix(0, nrow=ngroup, ncol=ngroup)
  
  for (rea_id in 1:n_rea){
    
    features_sub <- features %>% filter(rea == rea_id, ff == ffs)
    
    any_frac <- fracs[1]
    n_feat <- sum(features_sub$frac == any_frac)
    top_n <- max(1, floor(top_prop * n_feat))
    
    ours_top_list <- vector("list", ngroup)
    pl_top_list   <- vector("list", ngroup)
    
    for (ii in seq_along(fracs)){
      fr <- fracs[ii]
      ours_top_list[[ii]] <- top_our_gini(dat1_sub, features_sub, fr, seed=rea_id, gamma=gamma, top_n=top_n)
      pl_top_list[[ii]] <- top_prelect(dat1_sub, fr, seed=rea_id,
                                       lambda_steps=lambda_steps,
                                       inner_spl_ratio=inner_spl_ratio,
                                       n_cores=n_cores,
                                       max_iter_fit=max_iter_fit,
                                       top_n=top_n)
    }
    
    for (i in seq_along(fracs)){
      for (j in seq_along(fracs)){
        overlap_ours[i,j] <- overlap_ours[i,j] + length(intersect(ours_top_list[[i]], ours_top_list[[j]])) / top_n
        overlap_pl[i,j]   <- overlap_pl[i,j]   + length(intersect(pl_top_list[[i]],   pl_top_list[[j]]))   / top_n
      }
    }
  }
  
  overlap_ours <- overlap_ours / n_rea
  overlap_pl   <- overlap_pl   / n_rea
  
  diag(overlap_ours) <- NA
  diag(overlap_pl)   <- NA
  
  rownames(overlap_ours) <- colnames(overlap_ours) <- fracs
  rownames(overlap_pl)   <- colnames(overlap_pl)   <- fracs
  
  mean_ours <- apply(overlap_ours, 1, function(x) mean(x, na.rm=TRUE))
  mean_pl   <- apply(overlap_pl,   1, function(x) mean(x, na.rm=TRUE))
  
  sd_ours <- apply(overlap_ours, 1, function(x) sd(x, na.rm=TRUE))
  sd_pl   <- apply(overlap_pl,   1, function(x) sd(x, na.rm=TRUE))
  
  dat_plot <- data.frame(
    sample_fraction = rep(as.numeric(fracs), 2),
    Overlap = c(mean_ours, mean_pl),
    se = c(sd_ours / sqrt(n_rea), sd_pl / sqrt(n_rea)),
    types = c(rep("Gini^2.5", ngroup),
              rep("PreLect", ngroup))
  )
  
  pdf_list[[index]] <- ggplot(dat_plot, aes(sample_fraction, Overlap)) +
    geom_line(aes(group=types, color=types), linetype="dotdash", size=0.4) +
    geom_point(aes(color=types, shape=types), size=2.5) +
    geom_errorbar(aes(ymin=Overlap-se, ymax=Overlap+se, color=types), width=0.02) +
    scale_color_manual(values = c("Gini^2.5" = "#0F7BA2FF", "PreLect" = "#DD5129FF")) + 
    scale_shape_manual(values = c("Gini^2.5" = 17, "PreLect" = 16)) +
    theme_bw() + xlab("") + ylab("") +
    labs(title=paste("Fraction of spiking-in:", ffs)) +
    theme(
      plot.title = element_text(hjust=0.5, size=8),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size=8, color='black'),
      axis.text.y = element_text(size=8, color='black'),
      axis.ticks  = element_line(size=0.2),
      legend.position = "none"
    )
  
  index <- index + 1
}

pdf_list[[1]] <- pdf_list[[1]] + ylab("Frequency")
p2 <- do.call(ggarrange, c(pdf_list, list(labels=c("d","e","f"), ncol=3, nrow=1), align="hv"))
p2

### 

features <- read.csv("../results/Simulated0.2_features.csv")

ffs_list <- c(0.05, 0.1, 0.2)
fracs    <- sort(unique(features$frac))
ngroup   <- length(fracs)

n_rea    <- 10
top_prop <- 0.10
gamma    <- 2.5

lambda_steps    <- 50
inner_spl_ratio <- 0.7
n_cores         <- 8
max_iter_fit    <- 10000

drop_status <- function(dat){
  last <- dat[[ncol(dat)]]
  if (!is.numeric(last) || tolower(colnames(dat)[ncol(dat)]) %in% c("status","label","y","class")) {
    dat[, -ncol(dat), drop=FALSE]
  } else {
    dat
  }
}

top_our_gini <- function(dat_sf, features_sub, frac_val, seed, gamma, top_n){
  set.seed(seed)
  dat_x <- drop_status(dat_sf)
  n     <- nrow(dat_x)
  n_sub <- max(2, floor(frac_val * n))
  train_table <- dat_x[sample(seq_len(n), size=n_sub), , drop=FALSE]
  
  richness <- colSums(train_table > 0) / nrow(train_table)
  names(richness) <- colnames(train_table)
  
  sub <- features_sub %>% filter(frac == frac_val)
  r <- as.numeric(richness[sub$Feature])
  r[is.na(r)] <- 0
  
  score <- sub$MeanDecreaseGini * (r ^ gamma)
  ord <- order(score, decreasing = TRUE)
  sub$Feature[ord][1:top_n]
}

top_prelect <- function(dat_sf, frac_val, seed,
                        lambda_steps, inner_spl_ratio, n_cores, max_iter_fit,
                        top_n){
  set.seed(seed)
  
  y <- factor(dat_sf[[ncol(dat_sf)]])
  if (length(levels(y)) != 2) stop("PreLect classification requires y with 2 levels.")
  
  dat_x <- drop_status(dat_sf)
  n     <- nrow(dat_x)
  n_sub <- max(2, floor(frac_val * n))
  idx   <- sample(seq_len(n), size=n_sub)
  
  X_sf  <- dat_x[idx, , drop=FALSE]
  y_sub <- y[idx]
  
  X_raw <- t(as.matrix(X_sf))
  storage.mode(X_raw) <- "numeric"
  X_scaled <- t(scale(t(log1p(X_raw))))
  X_scaled[is.na(X_scaled)] <- 0
  
  lrange <- AutoScanning(X_scaled, X_raw, y_sub, task="classification", step=lambda_steps)
  out_dir_tmp <- tempdir()
  
  tuning_res <- LambdaTuningParallel(
    X_scaled, X_raw, y_sub, lrange,
    n_cores = n_cores,
    outpath = out_dir_tmp,
    spl_ratio = inner_spl_ratio,
    task = "classification"
  )
  
  pick <- LambdaDecision(tuning_res$TuningResult, tuning_res$PvlDistSummary,
                         maxdepth=5, minbucket=3)
  opt_lmbd <- pick$opt_lmbd
  
  prev <- GetPrevalence(X_raw)
  pre_out <- PreLect(X_scaled, prev, y_sub,
                     lambda = opt_lmbd, task="classification",
                     max_iter = max_iter_fit)
  
  featprop <- FeatureProperty(X_raw, y_sub, pre_out, task="classification")
  
  score <- abs(featprop$coef)
  ord <- order(score, decreasing = TRUE)
  featprop$FeatName[ord][1:top_n]
}

pdf_list <- list()
index <- 1

for (ffs in ffs_list){
  
  dat1_sub <- read.csv(
    paste("../data/Stool0.2_f_", 2, "_", ffs, ".csv", sep=""),
    row.names=1, header=TRUE
  )
  
  overlap_ours <- matrix(0, nrow=ngroup, ncol=ngroup)
  overlap_pl   <- matrix(0, nrow=ngroup, ncol=ngroup)
  
  for (rea_id in 1:n_rea){
    
    features_sub <- features %>% filter(rea == rea_id, ff == ffs)
    
    any_frac <- fracs[1]
    n_feat <- sum(features_sub$frac == any_frac)
    top_n <- max(1, floor(top_prop * n_feat))
    
    ours_top_list <- vector("list", ngroup)
    pl_top_list   <- vector("list", ngroup)
    
    for (ii in seq_along(fracs)){
      fr <- fracs[ii]
      ours_top_list[[ii]] <- top_our_gini(dat1_sub, features_sub, fr, seed=rea_id, gamma=gamma, top_n=top_n)
      pl_top_list[[ii]] <- top_prelect(dat1_sub, fr, seed=rea_id,
                                       lambda_steps=lambda_steps,
                                       inner_spl_ratio=inner_spl_ratio,
                                       n_cores=n_cores,
                                       max_iter_fit=max_iter_fit,
                                       top_n=top_n)
    }
    
    for (i in seq_along(fracs)){
      for (j in seq_along(fracs)){
        overlap_ours[i,j] <- overlap_ours[i,j] + length(intersect(ours_top_list[[i]], ours_top_list[[j]])) / top_n
        overlap_pl[i,j]   <- overlap_pl[i,j]   + length(intersect(pl_top_list[[i]],   pl_top_list[[j]]))   / top_n
      }
    }
  }
  
  overlap_ours <- overlap_ours / n_rea
  overlap_pl   <- overlap_pl   / n_rea
  
  diag(overlap_ours) <- NA
  diag(overlap_pl)   <- NA
  
  rownames(overlap_ours) <- colnames(overlap_ours) <- fracs
  rownames(overlap_pl)   <- colnames(overlap_pl)   <- fracs
  
  mean_ours <- apply(overlap_ours, 1, function(x) mean(x, na.rm=TRUE))
  mean_pl   <- apply(overlap_pl,   1, function(x) mean(x, na.rm=TRUE))
  
  sd_ours <- apply(overlap_ours, 1, function(x) sd(x, na.rm=TRUE))
  sd_pl   <- apply(overlap_pl,   1, function(x) sd(x, na.rm=TRUE))
  
  dat_plot <- data.frame(
    sample_fraction = rep(as.numeric(fracs), 2),
    Overlap = c(mean_ours, mean_pl),
    se = c(sd_ours / sqrt(n_rea), sd_pl / sqrt(n_rea)),
    types = c(rep("Gini^2.5", ngroup),
              rep("PreLect", ngroup))
  )
  
  pdf_list[[index]] <- ggplot(dat_plot, aes(sample_fraction, Overlap)) +
    geom_line(aes(group=types, color=types), linetype="dotdash", size=0.4) +
    geom_point(aes(color=types, shape=types), size=2.5) +
    geom_errorbar(aes(ymin=Overlap-se, ymax=Overlap+se, color=types), width=0.02) +
    scale_color_manual(values = c("Gini^2.5" = "#0F7BA2FF", "PreLect" = "#DD5129FF")) + 
    scale_shape_manual(values = c("Gini^2.5" = 17, "PreLect" = 16)) +
    theme_bw() + xlab("") + ylab("") +
    labs(title=paste("Fraction of spiking-in:", ffs)) +
    theme(
      plot.title = element_text(hjust=0.5, size=8),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size=8, color='black'),
      axis.text.y = element_text(size=8, color='black'),
      axis.ticks  = element_line(size=0.2),
      legend.position = "none"
    )
  
  index <- index + 1
}

pdf_list[[1]] <- pdf_list[[1]] + ylab("Frequency")
p3 <- do.call(ggarrange, c(pdf_list, list(labels=c("g","h","i"), ncol=3, nrow=1), align="hv"))
p3

row1 <- annotate_figure(
  p1,
  left = text_grob("Effect size = 0.01", rot = 90, size = 10)
)

row2 <- annotate_figure(
  p2,
  left = text_grob("Effect size = 0.05", rot = 90, size = 10)
)

row3 <- annotate_figure(
  p3,
  left = text_grob("Effect size = 0.20", rot = 90, size = 10)
)
p4 <- ggarrange(row1, row2, row3, ncol = 1, nrow = 3, align = "hv")
print(p4)

legend_plot <- ggplot(data.frame(x=1, y=1, types=factor(c("PreLect", "Gini^2.5"),
                                                        levels=c("PreLect", "Gini^2.5"))),
                      aes(x=x, y=y, color=types, shape=types)) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c("PreLect" = "#DD5129FF", "Gini^2.5" = "#0F7BA2FF"),
    labels = c("PreLect", "ParSlet")
  ) +
  scale_shape_manual(
    values = c("PreLect"=16, "Gini^2.5"=17),
    labels = c("PreLect", "ParSlet")
  ) +
  guides(
    color = guide_legend(ncol = 2, override.aes = list(size = 4)),
    shape = guide_legend(ncol = 2)
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.box = "horizontal"
  )

legend <- get_legend(legend_plot)

final_plot <- ggarrange(p4, legend, ncol = 1, heights = c(1, 0.1))
print(final_plot)

ggsave(final_plot,file=paste("../figs/prelect_l50.pdf",sep = ""), width=7.2, height=7.7, dpi = 500,scale = 0.9)
