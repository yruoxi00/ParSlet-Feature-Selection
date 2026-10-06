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
