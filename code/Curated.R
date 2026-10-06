library(curatedMetagenomicData)
library(dplyr)
library(DT)
library(randomForest)

taxa_table <- curatedMetagenomicData(pattern = "relative_abundance", dryrun = FALSE)

merged_taxa_table <- mergeData(taxa_table)

abundance_matrix <- assay(merged_taxa_table, "relative_abundance")
taxonomy_info <- rowData(merged_taxa_table)

# Create a taxonomic profile
taxonomic_profile <- data.frame(TaxonID = rownames(abundance_matrix), 
                                taxonomy_info, 
                                abundance_matrix)

sampleMetadata <- curatedMetagenomicData::sampleMetadata

sampleMetadata_sub <- sampleMetadata[sampleMetadata$disease%in%c("healthy","IBD", "CRC", "T2D", "RA","ACVD", "adenoma", "IGT"),]

freqs <- (table(sampleMetadata_sub$study_name,sampleMetadata_sub$disease))
freqs_df <- as.data.frame.matrix(freqs)
freqs_df <- freqs_df[rowSums(freqs_df>0)>=2,]


# CRC: c("FengQ_2015","GuptaA_2019","HanniganGD_2017")
# IBD: c("HallAB_2017","HMP_2019_ibdmdb","IjazUZ_2017","NielsenHB_2014")
# T2D: c("QinJ_2012","MetaCardis_2020_a","KarlssonFH_2013","SankaranarayananK_2015")
# adenoma: c("HanniganGD_2017","ThomasAM_2018a","YachidaS_2019","ZellerG_2014")
# IGT: c("KarlssonFH_2013","MetaCardis_2020_a","HMP_2019_t2d")

all_study <- c("FengQ_2015","GuptaA_2019","HanniganGD_2017","HallAB_2017","HMP_2019_ibdmdb","IjazUZ_2017","NielsenHB_2014",
               "QinJ_2012","MetaCardis_2020_a","KarlssonFH_2013","SankaranarayananK_2015",
               "HanniganGD_2017","ThomasAM_2018a","YachidaS_2019","ZellerG_2014",
               "KarlssonFH_2013","MetaCardis_2020_a","HMP_2019_t2d")


sampleMetadata_sub1 <- sampleMetadata_sub[sampleMetadata_sub$study_name%in%all_study,]
rownames(sampleMetadata_sub1) <- sampleMetadata_sub1$sample_id
overlaps <- intersect(colnames(taxonomic_profile),rownames(sampleMetadata_sub1))
sampleMetadata_sub1 <- sampleMetadata_sub1[overlaps,]
taxonomic_profile_sub <- taxonomic_profile[,overlaps] 
rownames(taxonomic_profile_sub) <- sub(".*s__", "", rownames(taxonomic_profile_sub))

features = NULL
for (study in unique(sampleMetadata_sub1$study_name)){
  for (phenotype in c("adenoma",  "CRC",  "IBD",  "IGT",  "T2D")){
    s1 <- which(sampleMetadata_sub1$study_name==study & sampleMetadata_sub1$disease == phenotype)
    s2 <- which(sampleMetadata_sub1$study_name==study & sampleMetadata_sub1$disease == "healthy")
    if ((length(s1)>0) & (length(s2)>0)){
      dat1_sub <- cbind(taxonomic_profile_sub[,s1],taxonomic_profile_sub[,s2])
      dat1_sub <- as.data.frame(t(dat1_sub))
      
      richness <- colSums(dat1_sub>0)
      dat1_sub$status <- c(rep(phenotype,length(s1)),rep("healthy",length(s2)))
      
      dat1_sub$status <- factor(dat1_sub$status)
      train_data <- dat1_sub
      rf_model <- randomForest(
        status ~ .,
        data = train_data,
        importance = TRUE,
        probability = TRUE
      )
      importance_data <- as.data.frame(importance(rf_model))
      importance_data$Feature <- rownames(importance_data)
      
      importance_data = importance_data[,3:5]
      importance_data$study = study
      importance_data$phenotype = phenotype
      importance_data$richness = richness
      features = rbind(features,importance_data)
    }
  }
}

write.csv(features,file = "Curated_feautures.csv",row.names = F)


library(dplyr)
library(tidyr)

phenos <- c("adenoma","CRC","IBD","IGT","T2D")
all_study_unique <- unique(all_study)

md <- sampleMetadata %>%
  filter(study_name %in% all_study_unique,
         disease %in% c("healthy", phenos))

counts_wide <- md %>%
  count(study_name, disease) %>%
  pivot_wider(names_from = disease, values_from = n, values_fill = 0)

baseline_by_study <- md %>%
  group_by(study_name) %>%
  summarise(
    PMID = paste(sort(unique(na.omit(PMID))), collapse="; "),
    country = paste(sort(unique(na.omit(country))), collapse="; "),
    body_site = paste(sort(unique(na.omit(body_site))), collapse="; "),
    sequencing_platform = paste(sort(unique(na.omit(sequencing_platform))), collapse="; "),
    N_total = n(),
    age_median = median(age, na.rm = TRUE),
    age_q1 = quantile(age, 0.25, na.rm = TRUE),
    age_q3 = quantile(age, 0.75, na.rm = TRUE),
    BMI_median = median(BMI, na.rm = TRUE),
    BMI_q1 = quantile(BMI, 0.25, na.rm = TRUE),
    BMI_q3 = quantile(BMI, 0.75, na.rm = TRUE),
    reads_median = median(number_reads, na.rm = TRUE),
    reads_q1 = quantile(number_reads, 0.25, na.rm = TRUE),
    reads_q3 = quantile(number_reads, 0.75, na.rm = TRUE),
    gender_levels = paste(sort(unique(na.omit(gender))), collapse="; "),
    .groups = "drop"
  )

baseline_table <- baseline_by_study %>%
  left_join(counts_wide, by = "study_name") %>%
  mutate(
    healthy = ifelse(is.na(healthy), 0, healthy),
    Total_cases = rowSums(across(all_of(phenos), ~ ifelse(is.na(.x), 0, .x))),
    Total_controls = healthy
  ) %>%
  arrange(desc(N_total))

baseline_table

case_control_table <- purrr::map_dfr(phenos, function(ph) {
  md %>%
    filter(disease %in% c("healthy", ph)) %>%
    group_by(study_name) %>%
    summarise(
      phenotype = ph,
      cases = sum(disease == ph, na.rm = TRUE),
      controls = sum(disease == "healthy", na.rm = TRUE),
      total = n(),
      .groups = "drop"
    ) %>%
    filter(cases > 0 & controls > 0)
}) %>%
  arrange(phenotype, desc(total))

case_control_table

write.csv(baseline_table, "MDAS.csv", row.names = FALSE)
write.csv(case_control_table, "MDAS_case_control.csv", row.names = FALSE)
