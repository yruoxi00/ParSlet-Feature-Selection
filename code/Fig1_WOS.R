library(ggplot2)
library(randomForest)
library(pROC)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggpubr)
library(unikn) 
library(treemap)
library(treemapify)
library(gridExtra)
library(cowplot)

setwd("/Users/yrx/Downloads/FeatureSelection/feature_selection/code")

df <- read.csv("../data/random_forest_microbiome_bib.csv")

df$Feature.Importance.Keywords <- as.character(df$Feature.Importance.Keywords)

df_keywords <- df %>%
  filter(Feature.Importance.Keywords != "") %>%
  separate_rows(Feature.Importance.Keywords, sep = ", ") 

keyword_mapping <- c(
  "biomarker" = "biomarker",
  "biomarkers" = "biomarker",
  "importance score" = "importance score",
  "feature importance" = "feature importance",
  "variable importance" = "variable importance"
)

df_keywords <- df_keywords %>%
  mutate(Feature.Importance.Keywords = recode(Feature.Importance.Keywords, !!!keyword_mapping))

df_keywords_summary <- df_keywords %>%
  group_by(Feature.Importance.Keywords) %>%
  summarise(Count = n())

years <- seq(min(df$Year, na.rm = TRUE), max(df$Year, na.rm = TRUE))

g1 <- ggplot(df, aes(x = factor(Year, levels = years))) +
  geom_bar(fill = "steelblue") +
  scale_x_discrete(drop = FALSE) +
  theme_bw() + xlab("Year") + ylab("Number of Papers") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold.italic", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    axis.ticks = element_line(size = 0.2),
    axis.title.x = element_text(size = 11, face = "bold"),
    axis.title.y = element_text(size = 11, face = "bold"), 
    legend.position = "none"
  )
print(g1)

df_keywords_summary$label2 <- paste0(df_keywords_summary$Feature.Importance.Keywords,
                                     "\n", df_keywords_summary$Count)
g2 <- ggplot(
  df_keywords_summary,
  aes(area = Count,
      fill = Feature.Importance.Keywords,
      label = label2)
) +
  geom_treemap() +
  geom_treemap_text(
    place = "centre",
    grow = FALSE,
    reflow = TRUE,
    colour = "black",
    fontface = "bold", 
    size = 10
  ) +
  scale_fill_brewer(palette = "Set3") +
  labs(y = "Feature importance 
indicators") +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    strip.text = element_text(face = "bold.italic", size = 10),
    axis.title.y = element_text(size = 11, face = "bold")
  )
print(g2)

g3 <- plot_grid(
  g1, g2,
  labels = c("a", "b"),
  label_size = 15,
  ncol = 1,
  align = "v",       
  axis = "lr",       
  rel_heights = c(1, 1) 
)
print(g3)

ggsave("../figs/WOS.png", g3, width = 6.5, height = 7, dpi = 500)
