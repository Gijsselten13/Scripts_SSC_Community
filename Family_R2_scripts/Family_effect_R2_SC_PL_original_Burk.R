library("phyloseq") #Version 1.44.0
library("vegan") #Version 2.6-4
library("dplyr") #Version 1.1.2

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Family effect - SynCom & Plant R2 - dom - Burkholderiaceae genera - generation of file =====
#otu table
KO_SSC =read.table(paste(working_directory,"KO_tables/Original/SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Metadata
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object
OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)

#Sample subsetting
cond="ES"

samples_df_sub <- subset(samples_df, samples_df$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

SynComs <- c("AtSC", "HvSC", "LjSC", "SSC")
Hosts <- c("At","Hv", "Lj")

for (syncom in SynComs){
  samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum == paste(syncom))
  
  samples_sub = sample_data(samples_df_sub_3)
  
  phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
  
  phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
  
  beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))
  
  Bray_curtis_df=beta_isolate_KO
  
  #Make PCoA plot for Bray Curtis Distance matrix
  pcoa = cmdscale(Bray_curtis_df, k=3, eig=T)
  points = as.data.frame(pcoa$points)
  colnames(points) = c("x", "y", "z") 
  eig = pcoa$eig
  
  points = merge(points,samples_df_sub_2, by = "row.names")
  rownames(points) <- points$Row.names
  points <- points %>% dplyr::select (-Row.names)
  
  points$Condition <- factor(points$Condition, levels = c("At","Hv", "Lj"))
  points$Nutrient <- factor(points$Nutrient, levels = c("low", "high"))
  points$Experiment  <- factor(points$Experiment, levels = c("R1", "R2"))
  
  metadata=points[,-c(1,2,3)]
  
  set.seed(1)
  SSC_bray_KO_adonis <- adonis2(beta_isolate_KO ~ Condition*Nutrient*Experiment, data=metadata, method="bray", permutations=999)
  
  colnames(points)[colnames(points) == "Condition"] <- "Plant"
  
  R2_Plant_KO <- SSC_bray_KO_adonis$R2[1]
  assign(paste("R2_Plant_KO_",syncom,sep =""),R2_Plant_KO)
}

for (host in Hosts){
  samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Condition == paste(host))
  
  samples_sub = sample_data(samples_df_sub_3)
  
  phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
  
  phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
  
  beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))
  
  Bray_curtis_df=beta_isolate_KO
  
  #Make PCoA plot for Bray Curtis Distance matrix
  pcoa = cmdscale(Bray_curtis_df, k=3, eig=T)
  points = as.data.frame(pcoa$points)
  colnames(points) = c("x", "y", "z") 
  eig = pcoa$eig
  
  points = merge(points,samples_df_sub_2, by = "row.names")
  rownames(points) <- points$Row.names
  points <- points %>% dplyr::select (-Row.names)
  
  points$Condition <- factor(points$Inoculum, levels = c("AtSC","HvSC", "LjSC", "SSC"))
  points$Nutrient <- factor(points$Nutrient, levels = c("low", "high"))
  points$Experiment  <- factor(points$Experiment, levels = c("R1", "R2"))
  
  metadata=points[,-c(1,2,3)]
  
  set.seed(1)
  SSC_bray_KO_adonis <- adonis2(beta_isolate_KO ~ Inoculum*Nutrient*Experiment, data=metadata, method="bray", permutations=999)
  
  colnames(points)[colnames(points) == "Inoculum"] <- "SynCom"
  
  R2_SynCom_KO <- SSC_bray_KO_adonis$R2[1]
  assign(paste("R2_SynCom_KO_",host,sep =""),R2_SynCom_KO)
}

#Family
fam_data <- data.frame()

tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")

fams_left <- c("Acidovorax","Cupriavidus","Pelomonas","Polaromonas","Rhizobacter","Variovorax")

SynComs <- c("AtSC","HvSC","LjSC","SSC")
Hosts <- c("At","Hv","Lj")

for (family in fams_left){
  #KOs
  genera_data <- read.table(paste(working_directory, "Family_R2/Genus_dom/",family,".tsv", sep =""), sep ="\t", header =T, row.names =1)
  
  #Samples TABLE
  samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t", row.names =1) #make the SampleID column into the row.names
  colnames(samples_df)[5]="Nutrient"
  samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
  samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")
  
  #Set the OTU, TAX and sample data for making phyloseq object
  OTU = otu_table(as.matrix(genera_data),taxa_are_rows = TRUE)
  
  #Sample subsetting
  
  cond="ES"
  samples_df_sub <- subset(samples_df, samples_df$Compartment == cond)
  samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")
  
  #Plant effect
  for (syncom in SynComs){
    samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum == paste(syncom))
    samples_sub = sample_data(samples_df_sub_3)
    
    phylo_sub = phyloseq(OTU, samples_sub)
    
    phylo_sub_RA=microbiome::transform(x = phylo_sub, transform = "compositional" )
    
    #Agglomerate to phylum-level and rename
    #Bray Curtis distance matrix
    beta_genus <- as.matrix(vegdist(t(phylo_sub_RA@otu_table@.Data), method = "bray", diag = T))
    mean_value_genus= mean(beta_genus)
    
    #Make PCoA plot for Bray Curtis Distance matrix
    pcoa = cmdscale(beta_genus, k=3, eig=T)
    points = as.data.frame(pcoa$points)
    colnames(points) = c("x", "y", "z") 
    eig = pcoa$eig
    
    points = merge(points,samples_df_sub_2, by = "row.names")
    rownames(points) <- points$Row.names
    points <- points %>% dplyr::select (-Row.names)
    
    points$Condition <- factor(points$Condition, levels = c("At","Hv", "Lj"))
    points$Nutrient <- factor(points$Nutrient, levels = c("low", "high"))
    points$Experiment  <- factor(points$Experiment, levels = c("R1", "R2"))
    
    metadata=points[,-c(1,2,3)]
    
    #  Run adonis PERMANOVA test
    set.seed(1)
    SSC_bray_adonis_fam <- adonis2(beta_genus ~ Condition*Nutrient*Experiment, data=metadata, method="bray", permutations=999)
    
    colnames(points)[colnames(points) == "Condition"] <- "Plant"
    
    if (paste(syncom) == "AtSC"){
      R2_Plant_KO <- R2_Plant_KO_AtSC
    } else if (paste(syncom) == "HvSC"){
      R2_Plant_KO <- R2_Plant_KO_HvSC
    } else if (paste(syncom) == "LjSC"){
      R2_Plant_KO <- R2_Plant_KO_LjSC
    } else {
      R2_Plant_KO <- R2_Plant_KO_SSC
    }
    
    R2_Plant_fam <- SSC_bray_adonis_fam$R2[1] - R2_Plant_KO
    
    fam_data_2 <- data.frame(t(data.frame(c(paste(family), paste(syncom),R2_Plant_fam,"KO"))))
    
    fam_data <- rbind(fam_data, fam_data_2)
  }
  
  for (host in Hosts){
    samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Condition == paste(host))
    samples_sub = sample_data(samples_df_sub_3)
    
    phylo_sub = phyloseq(OTU, samples_sub)
    
    phylo_sub_RA=microbiome::transform(x = phylo_sub, transform = "compositional" )
    
    #Agglomerate to phylum-level and rename
    #Bray Curtis distance matrix
    beta_genus <- as.matrix(vegdist(t(phylo_sub_RA@otu_table@.Data), method = "bray", diag = T))
    mean_value_genus= mean(beta_genus)
    
    #Make PCoA plot for Bray Curtis Distance matrix
    pcoa = cmdscale(beta_genus, k=3, eig=T)
    points = as.data.frame(pcoa$points)
    colnames(points) = c("x", "y", "z") 
    eig = pcoa$eig
    
    points = merge(points,samples_df_sub_2, by = "row.names")
    rownames(points) <- points$Row.names
    points <- points %>% dplyr::select (-Row.names)
    
    points$Inoculum <- factor(points$Inoculum, levels = c("AtSC","HvSC", "LjSC", "SSC"))
    points$Nutrient <- factor(points$Nutrient, levels = c("low", "high"))
    points$Experiment  <- factor(points$Experiment, levels = c("R1", "R2"))
    
    metadata=points[,-c(1,2,3)]
    
    #  Run adonis PERMANOVA test
    set.seed(1)
    SSC_bray_adonis_fam <- adonis2(beta_genus ~ Inoculum*Nutrient*Experiment, data=metadata, method="bray", permutations=999)
    
    colnames(points)[colnames(points) == "Inoculum"] <- "SynCom"
    
    if (paste(host) == "At"){
      R2_SynCom <- R2_SynCom_KO_At
    } else if (paste(host) == "Hv"){
      R2_SynCom <- R2_SynCom_KO_Hv
    } else {
      R2_SynCom <- R2_SynCom_KO_Lj
    }
    
    R2_SynCom_fam <- SSC_bray_adonis_fam$R2[1] - R2_SynCom
    fam_data_2 <- data.frame(t(data.frame(c(paste(family), paste(host),R2_SynCom_fam, "KO"))))
    
    fam_data <- rbind(fam_data, fam_data_2)
  }
}

row.names(fam_data) <- NULL
colnames(fam_data) <- c("Family", "Subset","R2_change","KO")

write.table(fam_data, paste(working_directory, "Family_R2/SSC_Gen_Burk_R2_effects_subs_with_dom.txt", sep = ""), sep ="\t", quote =F, col.names =T, row.names =T)
