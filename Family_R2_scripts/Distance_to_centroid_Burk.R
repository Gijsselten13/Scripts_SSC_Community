library("dplyr") #Version 1.1.2
library("phyloseq") #Version 1.44.0
library("vegan") #Version 2.6-4

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Distance to centroid per genus and SynCom - Burkholderiaceae =====

#KO profiles
all_iso_2 =read.table(paste(working_directory,"KO_genome/KO_SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Taxonomy
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")

colnames(all_iso_2)[grep("M.16",colnames(all_iso_2))] <- "M-16"
colnames(all_iso_2)[grep("M.6",colnames(all_iso_2))] <- "M-6"
colnames(all_iso_2)[grep("M.10",colnames(all_iso_2))] <- "M-10"
colnames(all_iso_2)[grep("M.11_2",colnames(all_iso_2))] <- "M-11_2"
colnames(all_iso_2)[grep("M.11",colnames(all_iso_2))] <- "M-11"
colnames(all_iso_2) <- gsub("X", "", colnames(all_iso_2))

hopla_2 <- data.frame()

SynComs <- c("AtSC", "HvSC", "LjSC")

for (syncom in SynComs){
  tax_sub <- tax_df_2[tax_df_2$SynCom == paste(syncom),]
  
  #Set the OTU, TAX and sample data for making phyloseq object
  OTU = otu_table(as.matrix(all_iso_2),taxa_are_rows = TRUE)
  samples = sample_data(tax_sub)
  
  phylo_sub = phyloseq(OTU,samples)
  
  phylo_sub_RA=microbiome::transform(x = phylo_sub, transform = "compositional" )
  
  #Bray Curtis distance matrix
  beta <- as.matrix(vegdist(t(phylo_sub_RA@otu_table@.Data), method = "bray", diag = T))
  row.names(beta) <- gsub("X", "", row.names(beta))
  
  #Make PCoA plot for Bray Curtis Distance matrix
  pcoa_tax = cmdscale(beta, k=3, eig=T)
  points_beta = as.data.frame(pcoa_tax$points)
  colnames(points_beta) = c("x", "y", "z") 
  eig = pcoa_tax$eig
  row.names(points_beta) <- gsub("X", "", row.names(points_beta))
  points_beta = merge(points_beta,tax_df_2, by = "row.names")
  row.names(points_beta) <- points_beta$Row.names
  points_beta <- points_beta %>% dplyr::select (-Row.names)
  
  # Define parameter for centroid calculation
  
  param="Genus"
  
  # Calculate centroids for each group
  centroids_tax <- points_beta %>%
    group_by(!!sym(param)) %>%
    dplyr::summarize(
      centroid_x = mean(x, na.rm = TRUE),
      centroid_y = mean(y, na.rm = TRUE),
      centroid_z = mean(z, na.rm = TRUE)
    )
  
  # Join centroids back to the original data
  data_with_centroids_tax <- left_join(points_beta, centroids_tax , by = param)
  # Calculate distance to centroid for each point
  data_with_centroids_tax  <- data_with_centroids_tax  %>% 
    rowwise() %>%
    mutate(distance_to_centroid = sqrt((x - centroid_x)^2 + 
                                         (y - centroid_y)^2 + 
                                         (z - centroid_z)^2))
  
  data_with_centroids_tax$distance_to_centroid
  data_with_centroids_tax_2 <- as.data.frame(data_with_centroids_tax)
  row.names(data_with_centroids_tax_2) <- row.names(points_beta)
  
  for (family in unique(data_with_centroids_tax_2$Genus)){
    data_with_centroids_tax_2_sub <- data_with_centroids_tax_2[data_with_centroids_tax_2$Genus == paste(family),]
    data_with_centroids_tax_2_sub <- data_with_centroids_tax_2_sub[!is.na(data_with_centroids_tax_2_sub$x),]
    average <- sum(data_with_centroids_tax_2_sub$distance_to_centroid)/length(data_with_centroids_tax_2_sub$distance_to_centroid)
    
    hopla <- data.frame(t(data.frame(c(paste(family), average, paste(syncom)))))
    hopla_2 <- rbind(hopla_2, hopla)
  }
}

tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")

fams_left <- c("Acidovorax","Cupriavidus","Pelomonas","Polaromonas","Rhizobacter","Variovorax")

row.names(hopla_2) <- NULL
colnames(hopla_2) <- c("Genus", "Distance_to_Centroid", "SynCom")

hopla_3 <- hopla_2[hopla_2$Genus %in% fams_left,]

write.table(hopla_3, paste(working_directory, "Family_R2/Dist_to_centroid_Burkholderiaceae_genus.txt", sep = ""), col.names =T, row.names =F, quote =F, sep = "\t")
