library("dplyr") #Version 1.1.2
library("phyloseq") #Version 1.44.0
library("vegan") #Version 2.6-4

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Script to generate R2_values_genus.txt - Necessary for S12 =====
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")

R2_values <- data.frame(matrix(NA, ncol = 3))
colnames(R2_values) <- c("run", "Tax_R2" ,"Func_R2")
R2_values_2 <- R2_values[-1,]

for (i in 1:1000){
  tax_table <- read.table(paste(working_directory,"simulation_R/run_",i,"/table_",i,".txt", sep = ""), sep = "\t", header = F)
  row.names(tax_table) <- tax_table$V1
  tax_table_2 <- tax_table %>% dplyr::select (-V1)
  vector <- c()
  for (j in 1:54){
    vector <- c(vector, paste("sample_", j, sep = ""))
  }
  colnames(tax_table_2) <- vector
  
  #Set the OTU, TAX and sample data for making phyloseq object
  OTU = otu_table(as.matrix(tax_table_2),taxa_are_rows = TRUE)
  #TAX = tax_table(tax_mat)
  TAX = tax_table(as.matrix(tax_df_2))
  
  phylo = phyloseq(OTU,TAX)
  
  phylo_RA=microbiome::transform(x = phylo, transform = "compositional" )
  
  #Agglomerate to phylum-level and rename
  phylo_RA <- phyloseq::tax_glom(phylo_RA, "Genus")
  phyloseq::taxa_names(phylo_RA) <- phyloseq::tax_table(phylo_RA)[, "Genus"]
  
  #Bray Curtis distance matrix
  beta_tax <- as.matrix(vegdist(t(phylo_RA@otu_table@.Data), method = "bray", diag = T))
  
  #Make PCoA plot for Bray Curtis Distance matrix
  pcoa = cmdscale(beta_tax, k=3, eig=T)
  points = as.data.frame(pcoa$points)
  colnames(points) = c("x", "y", "z") 
  eig = pcoa$eig
  
  set.seed(1)
  clusters <- kmeans(points, 9, iter.max = 10, nstart = 1)
  clusters_2 <- as.data.frame(clusters$cluster)
  colnames(clusters_2) <- "Cluster"
  
  set.seed(1)
  Tax_adonis <- adonis2(beta_tax ~ Cluster, data=clusters_2, method="bray", permutations=999)
  Tax_value <- Tax_adonis$R2[1]
  
  KO_table <- read.table(paste(working_directory,"simulation_R/run_",i,"/KO_table_",i,".txt", sep = ""), sep = "\t", header = T)
  row.names(KO_table) <- KO_table$function.
  KO_table_2 <- KO_table %>% dplyr::select (-function.)
  
  #Set the OTU, TAX and sample data for making phyloseq object
  OTU_KO = otu_table(as.matrix(KO_table_2),taxa_are_rows = TRUE)
  
  phylo_KO = phyloseq(OTU_KO)
  
  phylo_KO_RA=microbiome::transform(x = phylo_KO, transform = "compositional" )
  
  
  #Bray Curtis distance matrix
  beta_KO <- as.matrix(vegdist(t(phylo_KO_RA), method = "bray", diag = T))
  set.seed(1)
  KO_adonis <- adonis2(beta_KO ~ Cluster, data=clusters_2, method="bray", permutations=999)
  KO_value <- KO_adonis$R2[1]
  
  new <- c(paste("run_", i, sep = ""), Tax_value, KO_value)
  
  R2_values_2 <- rbind(R2_values_2, new)
}

colnames(R2_values_2) <- c("run", "Tax_R2" ,"Func_R2")
write.table(R2_values_2, paste(working_directory, "R2_values_genus.txt", sep= ""), quote = F, sep = "\t", col.names = T, row.names =T)
