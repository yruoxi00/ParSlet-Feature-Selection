

# R 4.0.3

setwd("/Users/xuwenwang/Downloads/SparseDOSSA2")

library(SparseDOSSA2)

################################################################### fraction
for (Nsam in c(0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5)){
  for (frac in c(0.05,0.1,0.2)){
    sim <- SparseDOSSA2(template = "Stool", n_sample = round(Nsam*332), new_features = FALSE,spike_metadata="both",
                        perc_feature_spiked_metadata=frac,metadata_effect_size=0.01)
    relative_abudance = as.data.frame(t(sim$simulated_matrices$rel))
    spike_ined = sim$spike_metadata$feature_metadata_spike_df
    labels = as.numeric(sim$spike_metadata$metadata_matrix[,2])
    labels[labels==0] = "Health"
    labels[labels==1] = "Nonhealth"
    relative_abudance$status = labels
    write.table(relative_abudance,file = paste('Stool0.01_f_',Nsam,'_',frac,'.csv',sep = ""),row.names = T,col.names = T,sep = ",")
    write.table(spike_ined,file = paste('Spike_Stool0.01__f_',Nsam,'_',frac,'.csv',sep = ""),row.names = T,col.names = T,sep = ",")
  }
}


################################################################### fraction
for (Nsam in c(0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5)){
  for (frac in c(0.05,0.1,0.2)){
    sim <- SparseDOSSA2(template = "Stool", n_sample = round(Nsam*332), new_features = FALSE,spike_metadata="both",
                        perc_feature_spiked_metadata=frac,metadata_effect_size=0.05)
    relative_abudance = as.data.frame(t(sim$simulated_matrices$rel))
    spike_ined = sim$spike_metadata$feature_metadata_spike_df
    labels = as.numeric(sim$spike_metadata$metadata_matrix[,2])
    labels[labels==0] = "Health"
    labels[labels==1] = "Nonhealth"
    relative_abudance$status = labels
    write.table(relative_abudance,file = paste('Stool0.05_f_',Nsam,'_',frac,'.csv',sep = ""),row.names = T,col.names = T,sep = ",")
    write.table(spike_ined,file = paste('Spike_Stool0.05__f_',Nsam,'_',frac,'.csv',sep = ""),row.names = T,col.names = T,sep = ",")
  }
}


