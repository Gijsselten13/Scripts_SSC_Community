library("dplyr") #Version 1.1.2
library("phyloseq") #Version 1.44.0
library("vegan") #Version 2.6-4
library("reshape2") #Version 1.4.4
library("ggplot2") #Version 3.4.2
library("webr") #Version 0.1.6

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Figure 3ab & Table S4 - PieDonut plots =====
combined_df_syncom_4_dom <- data.frame(matrix(NA, ncol = 11))
colnames(combined_df_syncom_4_dom) <- c("Df", "SumofSqs", "R2", "F", "Pr(>F)", "Variable","Subset", "Test", "Rank", "Dominance", "Drop_out")
combined_df_syncom_4_dom <- combined_df_syncom_4_dom[-1,]

combined_df_syncom_5_dom <- combined_df_syncom_4_dom

#otu table
norm_SSC=read.table(paste(working_directory,"Isolate_tables/Original/SSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
#Taxonomy table
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")
#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
colnames(samples_df_2)[5]="Nutrient"
samples_df_2$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df_2$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

sapply(tax_df, function(x) length(unique(x)))

#Set the OTU, TAX and sample data for making phyloseq object
OTU = otu_table(as.matrix(norm_SSC),taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(tax_df_2))

#Sample subsetting
cond="ES"
samples_df_sub <- subset(samples_df_2, samples_df_2$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

samples_sub = sample_data(samples_df_sub_2)

phylo_sub = phyloseq(OTU,TAX, samples_sub)

phylo_sub_RA=microbiome::transform(x = phylo_sub, transform = "compositional" )

#Agglomerate to phylum-level and rename
ps_phylum <- phyloseq::tax_glom(phylo_sub_RA, "Phylum")
ps_class <- phyloseq::tax_glom(phylo_sub_RA, "Class")
ps_order <- phyloseq::tax_glom(phylo_sub_RA, "Order")
ps_family <- phyloseq::tax_glom(phylo_sub_RA, "Family")
ps_genus <- phyloseq::tax_glom(phylo_sub_RA, "Genus")

phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]
phyloseq::taxa_names(ps_class) <- phyloseq::tax_table(ps_class)[, "Class"]
phyloseq::taxa_names(ps_order) <- phyloseq::tax_table(ps_order)[, "Order"]
phyloseq::taxa_names(ps_family) <- phyloseq::tax_table(ps_family)[, "Family"]
phyloseq::taxa_names(ps_genus) <- phyloseq::tax_table(ps_genus)[, "Genus"]

#Bray Curtis distance matrix
beta_phylum <- as.matrix(vegdist(t(ps_phylum@otu_table@.Data), method = "bray", diag = T))
beta_class <- as.matrix(vegdist(t(ps_class@otu_table@.Data), method = "bray", diag = T))
beta_order <- as.matrix(vegdist(t(ps_order@otu_table@.Data), method = "bray", diag = T))
beta_family <- as.matrix(vegdist(t(ps_family@otu_table@.Data), method = "bray", diag = T))
beta_genus <- as.matrix(vegdist(t(ps_genus@otu_table@.Data), method = "bray", diag = T))
beta_isolate <- as.matrix(vegdist(t(phylo_sub_RA@otu_table@.Data), method = "bray", diag = T))

#otu table
KO_SSC=read.table(paste(working_directory,"KO_tables/Original/SSC.tsv", sep =""), header=TRUE,sep="\t", row.names = 1)

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object
OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)

#Sample subsetting
cond="ES"

samples_df_sub <- subset(samples_df_2, samples_df_2$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

samples_sub = sample_data(samples_df_sub_2)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)

phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )

beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

#Data without Dominances
#otu table
norm_SSC=read.table(paste(working_directory,"Isolate_tables/No_dominances/SSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
#Taxonomy table
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")
#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
colnames(samples_df_2)[5]="Nutrient"
samples_df_2$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df_2$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

sapply(tax_df, function(x) length(unique(x)))

#Set the OTU, TAX and sample data for making phyloseq object
OTU = otu_table(as.matrix(norm_SSC),taxa_are_rows = TRUE)
#TAX = tax_table(tax_mat)
TAX = tax_table(as.matrix(tax_df_2))

#Sample subsetting
cond="ES"
samples_df_sub <- subset(samples_df_2, samples_df_2$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

samples_sub = sample_data(samples_df_sub_2)

phylo_sub = phyloseq(OTU,TAX, samples_sub)

phylo_sub_RA=microbiome::transform(x = phylo_sub, transform = "compositional" )

#Agglomerate to phylum-level and rename
ps_phylum <- phyloseq::tax_glom(phylo_sub_RA, "Phylum")
ps_class <- phyloseq::tax_glom(phylo_sub_RA, "Class")
ps_order <- phyloseq::tax_glom(phylo_sub_RA, "Order")
ps_family <- phyloseq::tax_glom(phylo_sub_RA, "Family")
ps_genus <- phyloseq::tax_glom(phylo_sub_RA, "Genus")

phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]
phyloseq::taxa_names(ps_class) <- phyloseq::tax_table(ps_class)[, "Class"]
phyloseq::taxa_names(ps_order) <- phyloseq::tax_table(ps_order)[, "Order"]
phyloseq::taxa_names(ps_family) <- phyloseq::tax_table(ps_family)[, "Family"]
phyloseq::taxa_names(ps_genus) <- phyloseq::tax_table(ps_genus)[, "Genus"]

#Bray Curtis distance matrix
beta_phylum_dom <- as.matrix(vegdist(t(ps_phylum@otu_table@.Data), method = "bray", diag = T))
beta_class_dom <- as.matrix(vegdist(t(ps_class@otu_table@.Data), method = "bray", diag = T))
beta_order_dom <- as.matrix(vegdist(t(ps_order@otu_table@.Data), method = "bray", diag = T))
beta_family_dom <- as.matrix(vegdist(t(ps_family@otu_table@.Data), method = "bray", diag = T))
beta_genus_dom <- as.matrix(vegdist(t(ps_genus@otu_table@.Data), method = "bray", diag = T))
beta_isolate_dom <- as.matrix(vegdist(t(phylo_sub_RA@otu_table@.Data), method = "bray", diag = T))

#otu table
KO_SSC=read.table(paste(working_directory,"KO_tables/No_dominances/SSC.tsv", sep =""), header=TRUE,sep="\t", row.names = 1)

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object
OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)

#Sample subsetting
cond="ES"

samples_df_sub <- subset(samples_df_2, samples_df_2$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

samples_sub = sample_data(samples_df_sub_2)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)

phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )

beta_isolate_KO_dom <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

beta_distance_list=list(beta_isolate, beta_genus, beta_family, beta_order, beta_class, beta_phylum)

for (i in 1:length(beta_distance_list)) {
  Bray_curtis_df=beta_distance_list[[i]]
  
  #Subsetting Bray Curtis matrix based on the sample subsetting on the lines above
  Bray_match <- as.matrix(as.logical(match(row.names(Bray_curtis_df),row.names(samples_df_sub))))
  Bray_Curtis_SSC <- subset(Bray_curtis_df, Bray_match)
  Bray_Curtis_SSC_T <- t(Bray_Curtis_SSC)
  Bray_Curtis_SSC_T <- subset(Bray_Curtis_SSC_T,Bray_match)
  Bray_Curtis_SSC <- t(Bray_Curtis_SSC_T)
  
  #Make PCoA plot for Bray Curtis Distance matrix
  pcoa = cmdscale(Bray_curtis_df, k=3, eig=T)
  points = as.data.frame(pcoa$points)
  colnames(points) = c("x", "y", "z") 
  eig = pcoa$eig
  
  points = merge(points,samples_df_sub_2, by = "row.names")
  rownames(points) <- points$Row.names
  points <- points %>% dplyr::select (-Row.names)
  
  points$Condition <- factor(points$Condition, levels = c("At","Hv", "Lj"))
  points$Inoculum <- factor(points$Inoculum, levels = c("SSC","AtSC", "HvSC", "LjSC","NS"))
  points$Nutrient <- factor(points$Nutrient, levels = c("low", "high"))
  points$Experiment  <- factor(points$Experiment, levels = c("R1", "R2"))
  
  points$Compartment <- factor(points$Compartment, levels = c("NP", "RZ","ES","dom"))
  points$Growth_condition <- factor(points$Growth_condition, levels = c("GC1_Lj","GC2_At","GC3_Hv"))
  points$Plant_inoc=paste(points$Plant_species,"_",points$Inoculum)
  metadata=points[,-c(1,2,3)]
}

comb_plant_list=list()
comb_syncom_list=list()

Taxo_order <- c("Isolate", "Genus", "Family", "Order", "Phylum","KO_functions")

for (i in 1:6) {
  #Per SynCom
  metadata_AtSC <- metadata[metadata$Inoculum == "AtSC",]
  taxo_beta_AtSC <- beta_distance_list[[i]][row.names(beta_distance_list[[i]]) %in% row.names(metadata_AtSC),row.names(beta_distance_list[[i]]) %in% row.names(metadata_AtSC)]
  taxo_adonis_AtSC <- adonis2(taxo_beta_AtSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_AtSC, method="bray", permutations=999)
  
  lala=pairwiseAdonis::pairwise.adonis(x =taxo_beta_AtSC, factors = metadata_AtSC$Condition, perm = 999 )
  
  AtSC_taxo_df <- cbind(taxo_adonis_AtSC,Variable = rownames(taxo_adonis_AtSC), Subset = "AtSC", Test = "Taxonomy", Rank=Taxo_order[i])
  
  metadata_HvSC <- metadata[metadata$Inoculum == "HvSC",]
  taxo_beta_HvSC <- beta_distance_list[[i]][row.names(beta_distance_list[[i]]) %in% row.names(metadata_HvSC),row.names(beta_distance_list[[i]]) %in% row.names(metadata_HvSC)]
  taxo_adonis_HvSC <- adonis2(taxo_beta_HvSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_HvSC, method="bray", permutations=999)
  HvSC_taxo_df <- cbind(taxo_adonis_HvSC,Variable = rownames(taxo_adonis_HvSC), Subset = "HvSC", Test = "Taxonomy", Rank=Taxo_order[i])
  
  lala_HvSC=pairwiseAdonis::pairwise.adonis(x =taxo_beta_HvSC, factors = metadata_HvSC$Condition, perm = 999 )
  
  metadata_LjSC <- metadata[metadata$Inoculum == "LjSC",]
  taxo_beta_LjSC <- beta_distance_list[[i]][row.names(beta_distance_list[[i]]) %in% row.names(metadata_LjSC),row.names(beta_distance_list[[i]]) %in% row.names(metadata_LjSC)]
  taxo_adonis_LjSC <- adonis2(taxo_beta_LjSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_LjSC, method="bray", permutations=999)
  LjSC_taxo_df <- cbind(taxo_adonis_LjSC,Variable = rownames(taxo_adonis_LjSC), Subset = "LjSC", Test = "Taxonomy", Rank=Taxo_order[i])
  
  lala_LjSC=pairwiseAdonis::pairwise.adonis(x =taxo_beta_LjSC, factors = metadata_LjSC$Condition, perm = 999 )
  
  metadata_SSC <- metadata[metadata$Inoculum == "SSC",]
  taxo_beta_SSC <- beta_distance_list[[i]][row.names(beta_distance_list[[i]]) %in% row.names(metadata_SSC),row.names(beta_distance_list[[i]]) %in% row.names(metadata_SSC)]
  taxo_adonis_SSC <- adonis2(taxo_beta_SSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_SSC, method="bray", permutations=999)
  SSC_taxo_df <- cbind(taxo_adonis_SSC,Variable = rownames(taxo_adonis_SSC), Subset = "SSC", Test = "Taxonomy", Rank=Taxo_order[i])
  
  lala_SSC=pairwiseAdonis::pairwise.adonis(x =taxo_beta_SSC, factors = metadata_SSC$Condition, perm = 999 )
  
  #Per plant
  metadata_At <- metadata[metadata$Condition == "At",]
  taxo_beta_At <- beta_distance_list[[i]][row.names(beta_distance_list[[i]]) %in% row.names(metadata_At),row.names(beta_distance_list[[i]]) %in% row.names(metadata_At)]
  taxo_adonis_At <- adonis2(taxo_beta_At ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_At, method="bray", permutations=999)
  At_taxo_df <- cbind(taxo_adonis_At,Variable = rownames(taxo_adonis_At), Subset = "At", Test = "Taxonomy", Rank=Taxo_order[i])
  
  lala_At=pairwiseAdonis::pairwise.adonis(x =taxo_beta_At, factors = metadata_At$Inoculum, perm = 999 )
  
  metadata_Hv <- metadata[metadata$Condition == "Hv",]
  taxo_beta_Hv <- beta_distance_list[[i]][row.names(beta_distance_list[[i]]) %in% row.names(metadata_Hv),row.names(beta_distance_list[[i]]) %in% row.names(metadata_Hv)]
  taxo_adonis_Hv <- adonis2(taxo_beta_Hv ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_Hv, method="bray", permutations=999)
  Hv_taxo_df <- cbind(taxo_adonis_Hv,Variable = rownames(taxo_adonis_Hv), Subset = "Hv", Test = "Taxonomy", Rank=Taxo_order[i])
  
  lala_Hv=pairwiseAdonis::pairwise.adonis(x =taxo_beta_Hv, factors = metadata_Hv$Inoculum, perm = 999 )
  
  metadata_Lj <- metadata[metadata$Condition == "At",]
  taxo_beta_Lj <- beta_distance_list[[i]][row.names(beta_distance_list[[i]]) %in% row.names(metadata_Lj),row.names(beta_distance_list[[i]]) %in% row.names(metadata_Lj)]
  taxo_adonis_Lj <- adonis2(taxo_beta_Lj ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_Lj, method="bray", permutations=999)
  Lj_taxo_df <- cbind(taxo_adonis_Lj,Variable = rownames(taxo_adonis_Lj), Subset = "Lj", Test = "Taxonomy", Rank=Taxo_order[i])
  
  lala_Lj=pairwiseAdonis::pairwise.adonis(x =taxo_beta_Lj, factors = metadata_Lj$Inoculum, perm = 999 )
  
  Combined_df_Syncom=rbind(AtSC_taxo_df, HvSC_taxo_df,LjSC_taxo_df,SSC_taxo_df)
  Sub_comb_syncom <- Combined_df_Syncom[!(Combined_df_Syncom$Variable %in% c("Total", "Residual", "Condition:Nutrient", "Condition:Experiment")), ]
  comb_syncom_list[[i]]=Sub_comb_syncom
  
  Combined_df_Plant=rbind(At_taxo_df,Hv_taxo_df,Lj_taxo_df)
  Sub_comb_plant <- Combined_df_Plant[!(Combined_df_Plant$Variable %in% c("Total", "Residual", "Inoculum:Nutrient", "Inoculum:Experiment")), ]
  comb_plant_list[[i]]=Sub_comb_plant
  
}

#Separation by plants/ Syncom KO Functionality 
#Per SynCom
metadata_AtSC <- metadata[metadata$Inoculum == "AtSC",]
KO_beta_AtSC <- beta_isolate_KO[row.names(beta_isolate_KO) %in% row.names(metadata_AtSC),row.names(beta_isolate_KO) %in% row.names(metadata_AtSC)]
KO_adonis_AtSC <- adonis2(KO_beta_AtSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_AtSC, method="bray", permutations=999)
AtSC_KO_df <- cbind(KO_adonis_AtSC,Variable = rownames(KO_adonis_AtSC), Subset = "AtSC", Test = "Functions", Rank="Functions")

lala_KO_AtSC=pairwiseAdonis::pairwise.adonis(x =KO_beta_AtSC, factors = metadata_AtSC$Condition, perm = 999 )

metadata_HvSC <- metadata[metadata$Inoculum == "HvSC",]
KO_beta_HvSC <- beta_isolate_KO[row.names(beta_isolate_KO) %in% row.names(metadata_HvSC),row.names(beta_isolate_KO) %in% row.names(metadata_HvSC)]
KO_adonis_HvSC <- adonis2(KO_beta_HvSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_HvSC, method="bray", permutations=999)
HvSC_KO_df <- cbind(KO_adonis_HvSC,Variable = rownames(KO_adonis_HvSC), Subset = "HvSC", Test = "Functions", Rank="Functions")

lala_KO_HvSC=pairwiseAdonis::pairwise.adonis(x =KO_beta_HvSC, factors = metadata_HvSC$Condition, perm = 999 )

metadata_LjSC <- metadata[metadata$Inoculum == "LjSC",]
KO_beta_LjSC <- beta_isolate_KO[row.names(beta_isolate_KO) %in% row.names(metadata_LjSC),row.names(beta_isolate_KO) %in% row.names(metadata_LjSC)]
KO_adonis_LjSC <- adonis2(KO_beta_LjSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_LjSC, method="bray", permutations=999)
LjSC_KO_df <- cbind(KO_adonis_LjSC,Variable = rownames(KO_adonis_LjSC), Subset = "LjSC", Test = "Functions", Rank="Functions")

lala_KO_LjSC=pairwiseAdonis::pairwise.adonis(x =KO_beta_LjSC, factors = metadata_LjSC$Condition, perm = 999 )

metadata_SSC <- metadata[metadata$Inoculum == "SSC",]
KO_beta_SSC <- beta_isolate_KO[row.names(beta_isolate_KO) %in% row.names(metadata_SSC),row.names(beta_isolate_KO) %in% row.names(metadata_SSC)]
KO_adonis_SSC <- adonis2(KO_beta_SSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_SSC, method="bray", permutations=999)
SSC_KO_df <- cbind(KO_adonis_SSC,Variable = rownames(KO_adonis_SSC), Subset = "SSC", Test = "Functions", Rank="Functions")

lala_KO_SSC=pairwiseAdonis::pairwise.adonis(x =KO_beta_SSC, factors = metadata_SSC$Condition, perm = 999 )

#Per plant
metadata_At <- metadata[metadata$Condition == "At",]
KO_beta_At <- beta_isolate_KO[row.names(beta_isolate_KO) %in% row.names(metadata_At),row.names(beta_isolate_KO) %in% row.names(metadata_At)]
KO_adonis_At <- adonis2(KO_beta_At ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_At, method="bray", permutations=999)
At_KO_df <- cbind(KO_adonis_At,Variable = rownames(KO_adonis_At), Subset = "At", Test = "Functions", Rank="Functions")

metadata_Hv <- metadata[metadata$Condition == "Hv",]
KO_beta_Hv <- beta_isolate_KO[row.names(beta_isolate_KO) %in% row.names(metadata_Hv),row.names(beta_isolate_KO) %in% row.names(metadata_Hv)]
KO_adonis_Hv <- adonis2(KO_beta_Hv ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_Hv, method="bray", permutations=999)
Hv_KO_df <- cbind(KO_adonis_Hv,Variable = rownames(KO_adonis_Hv), Subset = "Hv", Test = "Functions", Rank="Functions")

metadata_Lj <- metadata[metadata$Condition == "Lj",]
KO_beta_Lj <- beta_isolate_KO[row.names(beta_isolate_KO) %in% row.names(metadata_Lj),row.names(beta_isolate_KO) %in% row.names(metadata_Lj)]
KO_adonis_Lj <- adonis2(KO_beta_Lj ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_Lj, method="bray", permutations=999)
Lj_KO_df <- cbind(KO_adonis_Lj,Variable = rownames(KO_adonis_Lj), Subset = "Lj", Test = "Functions", Rank="Functions")

#Combining Data
combined_df_plant <- do.call(rbind, comb_plant_list)
combined_df_plant <- rbind(combined_df_plant,At_KO_df,Hv_KO_df,Lj_KO_df )
combined_df_plant <- combined_df_plant[!(combined_df_plant$Variable %in% c("Total", "Residual", "Inoculum:Nutrient", "Inoculum:Experiment", "Experiment", "Nutrient")), ]
combined_df_plant$Rank=factor(combined_df_plant$Rank, levels = c("Isolate", "Genus", "Family", "Order", "Class", "Phylum","Functions"))
combined_df_plant$Test=factor(combined_df_plant$Test, levels = c("Taxonomy", "Functions"))

combined_df_syncom <- do.call(rbind, comb_syncom_list)
combined_df_syncom <- rbind(combined_df_syncom,AtSC_KO_df,HvSC_KO_df,LjSC_KO_df,SSC_KO_df)
combined_df_syncom <- combined_df_syncom[!(combined_df_syncom$Variable %in% c("Total", "Residual", "Condition:Nutrient", "Condition:Experiment","Experiment", "Nutrient")), ]
combined_df_syncom$Rank=factor(combined_df_syncom$Rank, levels = c("Isolate", "Genus", "Family", "Order", "Class", "Phylum","Functions"))
combined_df_syncom$Test=factor(combined_df_syncom$Test, levels = c("Taxonomy", "Functions"))

#Dominances drop out data
beta_distance_list_dom=list(beta_isolate_dom, beta_genus_dom, beta_family_dom, beta_order_dom, beta_class_dom, beta_phylum_dom)
comb_plant_list_dom=list()
comb_syncom_list_dom=list()

for (i in 1:6) {
  #Per SynCom
  metadata_AtSC <- metadata[metadata$Inoculum == "AtSC",]
  taxo_beta_AtSC <- beta_distance_list_dom[[i]][row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_AtSC),row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_AtSC)]
  taxo_adonis_AtSC <- adonis2(taxo_beta_AtSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_AtSC, method="bray", permutations=999)
  
  lala=pairwiseAdonis::pairwise.adonis(x =taxo_beta_AtSC, factors = metadata_AtSC$Condition, perm = 999 )
  
  AtSC_taxo_df <- cbind(taxo_adonis_AtSC,Variable = rownames(taxo_adonis_AtSC), Subset = "AtSC", Test = "Taxonomy", Rank=Taxo_order[i])
  
  metadata_HvSC <- metadata[metadata$Inoculum == "HvSC",]
  taxo_beta_HvSC <- beta_distance_list_dom[[i]][row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_HvSC),row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_HvSC)]
  taxo_adonis_HvSC <- adonis2(taxo_beta_HvSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_HvSC, method="bray", permutations=999)
  HvSC_taxo_df <- cbind(taxo_adonis_HvSC,Variable = rownames(taxo_adonis_HvSC), Subset = "HvSC", Test = "Taxonomy", Rank=Taxo_order[i])
  
  lala_HvSC=pairwiseAdonis::pairwise.adonis(x =taxo_beta_HvSC, factors = metadata_HvSC$Condition, perm = 999 )
  
  metadata_LjSC <- metadata[metadata$Inoculum == "LjSC",]
  taxo_beta_LjSC <- beta_distance_list_dom[[i]][row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_LjSC),row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_LjSC)]
  taxo_adonis_LjSC <- adonis2(taxo_beta_LjSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_LjSC, method="bray", permutations=999)
  LjSC_taxo_df <- cbind(taxo_adonis_LjSC,Variable = rownames(taxo_adonis_LjSC), Subset = "LjSC", Test = "Taxonomy", Rank=Taxo_order[i])
  
  lala_LjSC=pairwiseAdonis::pairwise.adonis(x =taxo_beta_LjSC, factors = metadata_LjSC$Condition, perm = 999 )
  
  metadata_SSC <- metadata[metadata$Inoculum == "SSC",]
  taxo_beta_SSC <- beta_distance_list_dom[[i]][row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_SSC),row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_SSC)]
  taxo_adonis_SSC <- adonis2(taxo_beta_SSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_SSC, method="bray", permutations=999)
  SSC_taxo_df <- cbind(taxo_adonis_SSC,Variable = rownames(taxo_adonis_SSC), Subset = "SSC", Test = "Taxonomy", Rank=Taxo_order[i])
  
  lala_SSC=pairwiseAdonis::pairwise.adonis(x =taxo_beta_SSC, factors = metadata_SSC$Condition, perm = 999 )
  
  #Per plant
  metadata_At <- metadata[metadata$Condition == "At",]
  taxo_beta_At <- beta_distance_list_dom[[i]][row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_At),row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_At)]
  taxo_adonis_At <- adonis2(taxo_beta_At ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_At, method="bray", permutations=999)
  At_taxo_df <- cbind(taxo_adonis_At,Variable = rownames(taxo_adonis_At), Subset = "At", Test = "Taxonomy", Rank=Taxo_order[i])
  
  lala_At=pairwiseAdonis::pairwise.adonis(x =taxo_beta_At, factors = metadata_At$Inoculum, perm = 999 )
  
  metadata_Hv <- metadata[metadata$Condition == "Hv",]
  taxo_beta_Hv <- beta_distance_list_dom[[i]][row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_Hv),row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_Hv)]
  taxo_adonis_Hv <- adonis2(taxo_beta_Hv ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_Hv, method="bray", permutations=999)
  Hv_taxo_df <- cbind(taxo_adonis_Hv,Variable = rownames(taxo_adonis_Hv), Subset = "Hv", Test = "Taxonomy", Rank=Taxo_order[i])
  
  lala_Hv=pairwiseAdonis::pairwise.adonis(x =taxo_beta_Hv, factors = metadata_Hv$Inoculum, perm = 999 )
  
  metadata_Lj <- metadata[metadata$Condition == "At",]
  taxo_beta_Lj <- beta_distance_list_dom[[i]][row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_Lj),row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_Lj)]
  taxo_adonis_Lj <- adonis2(taxo_beta_Lj ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_Lj, method="bray", permutations=999)
  Lj_taxo_df <- cbind(taxo_adonis_Lj,Variable = rownames(taxo_adonis_Lj), Subset = "Lj", Test = "Taxonomy", Rank=Taxo_order[i])
  
  lala_Lj=pairwiseAdonis::pairwise.adonis(x =taxo_beta_Lj, factors = metadata_Lj$Inoculum, perm = 999 )
  
  Combined_df_Syncom=rbind(AtSC_taxo_df, HvSC_taxo_df,LjSC_taxo_df,SSC_taxo_df)
  Sub_comb_syncom <- Combined_df_Syncom[!(Combined_df_Syncom$Variable %in% c("Total", "Residual", "Condition:Nutrient", "Condition:Experiment")), ]
  comb_syncom_list_dom[[i]]=Sub_comb_syncom
  
  Combined_df_Plant=rbind(At_taxo_df,Hv_taxo_df,Lj_taxo_df)
  Sub_comb_plant <- Combined_df_Plant[!(Combined_df_Plant$Variable %in% c("Total", "Residual", "Inoculum:Nutrient", "Inoculum:Experiment")), ]
  comb_plant_list_dom[[i]]=Sub_comb_plant
  
}

#Separation by plants/ Syncom KO Functionality 
#Per SynCom
metadata_AtSC <- metadata[metadata$Inoculum == "AtSC",]
KO_beta_AtSC <- beta_isolate_KO_dom[row.names(beta_isolate_KO_dom) %in% row.names(metadata_AtSC),row.names(beta_isolate_KO_dom) %in% row.names(metadata_AtSC)]
KO_adonis_AtSC <- adonis2(KO_beta_AtSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_AtSC, method="bray", permutations=999)
AtSC_KO_df <- cbind(KO_adonis_AtSC,Variable = rownames(KO_adonis_AtSC), Subset = "AtSC", Test = "Functions", Rank="Functions")

lala_KO_AtSC=pairwiseAdonis::pairwise.adonis(x =KO_beta_AtSC, factors = metadata_AtSC$Condition, perm = 999 )

metadata_HvSC <- metadata[metadata$Inoculum == "HvSC",]
KO_beta_HvSC <- beta_isolate_KO_dom[row.names(beta_isolate_KO_dom) %in% row.names(metadata_HvSC),row.names(beta_isolate_KO_dom) %in% row.names(metadata_HvSC)]
KO_adonis_HvSC <- adonis2(KO_beta_HvSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_HvSC, method="bray", permutations=999)
HvSC_KO_df <- cbind(KO_adonis_HvSC,Variable = rownames(KO_adonis_HvSC), Subset = "HvSC", Test = "Functions", Rank="Functions")

lala_KO_HvSC=pairwiseAdonis::pairwise.adonis(x =KO_beta_HvSC, factors = metadata_HvSC$Condition, perm = 999 )

metadata_LjSC <- metadata[metadata$Inoculum == "LjSC",]
KO_beta_LjSC <- beta_isolate_KO_dom[row.names(beta_isolate_KO_dom) %in% row.names(metadata_LjSC),row.names(beta_isolate_KO_dom) %in% row.names(metadata_LjSC)]
KO_adonis_LjSC <- adonis2(KO_beta_LjSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_LjSC, method="bray", permutations=999)
LjSC_KO_df <- cbind(KO_adonis_LjSC,Variable = rownames(KO_adonis_LjSC), Subset = "LjSC", Test = "Functions", Rank="Functions")

lala_KO_LjSC=pairwiseAdonis::pairwise.adonis(x =KO_beta_LjSC, factors = metadata_LjSC$Condition, perm = 999 )

metadata_SSC <- metadata[metadata$Inoculum == "SSC",]
KO_beta_SSC <- beta_isolate_KO_dom[row.names(beta_isolate_KO_dom) %in% row.names(metadata_SSC),row.names(beta_isolate_KO_dom) %in% row.names(metadata_SSC)]
KO_adonis_SSC <- adonis2(KO_beta_SSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_SSC, method="bray", permutations=999)
SSC_KO_df <- cbind(KO_adonis_SSC,Variable = rownames(KO_adonis_SSC), Subset = "SSC", Test = "Functions", Rank="Functions")

lala_KO_SSC=pairwiseAdonis::pairwise.adonis(x =KO_beta_SSC, factors = metadata_SSC$Condition, perm = 999 )

#Per plant
metadata_At <- metadata[metadata$Condition == "At",]
KO_beta_At <- beta_isolate_KO_dom[row.names(beta_isolate_KO_dom) %in% row.names(metadata_At),row.names(beta_isolate_KO_dom) %in% row.names(metadata_At)]
KO_adonis_At <- adonis2(KO_beta_At ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_At, method="bray", permutations=999)
At_KO_df <- cbind(KO_adonis_At,Variable = rownames(KO_adonis_At), Subset = "At", Test = "Functions", Rank="Functions")

metadata_Hv <- metadata[metadata$Condition == "Hv",]
KO_beta_Hv <- beta_isolate_KO_dom[row.names(beta_isolate_KO_dom) %in% row.names(metadata_Hv),row.names(beta_isolate_KO_dom) %in% row.names(metadata_Hv)]
KO_adonis_Hv <- adonis2(KO_beta_Hv ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_Hv, method="bray", permutations=999)
Hv_KO_df <- cbind(KO_adonis_Hv,Variable = rownames(KO_adonis_Hv), Subset = "Hv", Test = "Functions", Rank="Functions")

metadata_Lj <- metadata[metadata$Condition == "Lj",]
KO_beta_Lj <- beta_isolate_KO_dom[row.names(beta_isolate_KO_dom) %in% row.names(metadata_Lj),row.names(beta_isolate_KO_dom) %in% row.names(metadata_Lj)]
KO_adonis_Lj <- adonis2(KO_beta_Lj ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_Lj, method="bray", permutations=999)
Lj_KO_df <- cbind(KO_adonis_Lj,Variable = rownames(KO_adonis_Lj), Subset = "Lj", Test = "Functions", Rank="Functions")

#Combining Data
combined_df_plant_dom <- do.call(rbind, comb_plant_list_dom)
combined_df_plant_dom <- rbind(combined_df_plant_dom,At_KO_df,Hv_KO_df,Lj_KO_df )
combined_df_plant_dom <- combined_df_plant_dom[!(combined_df_plant_dom$Variable %in% c("Total", "Residual", "Inoculum:Nutrient", "Inoculum:Experiment", "Experiment", "Nutrient")), ]
combined_df_plant_dom$Rank=factor(combined_df_plant_dom$Rank, levels = c("Isolate", "Genus", "Family", "Order", "Class", "Phylum","Functions"))
combined_df_plant_dom$Test=factor(combined_df_plant_dom$Test, levels = c("Taxonomy", "Functions"))

combined_df_syncom_dom <- do.call(rbind, comb_syncom_list_dom)
combined_df_syncom_dom <- rbind(combined_df_syncom_dom,AtSC_KO_df,HvSC_KO_df,LjSC_KO_df,SSC_KO_df)
combined_df_syncom_dom <- combined_df_syncom_dom[!(combined_df_syncom_dom$Variable %in% c("Total", "Residual", "Condition:Nutrient", "Condition:Experiment","Experiment", "Nutrient")), ]
combined_df_syncom_dom$Rank=factor(combined_df_syncom_dom$Rank, levels = c("Isolate", "Genus", "Family", "Order", "Class", "Phylum","Functions"))
combined_df_syncom_dom$Test=factor(combined_df_syncom_dom$Test, levels = c("Taxonomy", "Functions"))

beta_distance_list_dom=list(beta_isolate_dom, beta_genus_dom, beta_family_dom, beta_order_dom, beta_class_dom)
comb_plant_list_dom=list()
comb_syncom_list_dom=list()

beta_distance_list=list(beta_isolate, beta_genus, beta_family, beta_order, beta_class)
Plant_list=c("None","Arabidopsis", "Barley", "Lotus")
comb_plant_list=list()
comb_syncom_list=list()
line_plot_list=list()

metadata$Condition <- as.character(metadata$Condition)
metadata$Condition[metadata$Condition == "At"] <- "Arabidopsis"
metadata$Condition[metadata$Condition == "Hv"] <- "Barley"
metadata$Condition[metadata$Condition == "Lj"] <- "Lotus"

#original table
for (j in 1:4) {
  for (i in 1:5) {
    #Per SynCom
    metadata_AtSC <- metadata[metadata$Inoculum == "AtSC",]
    metadata_AtSC <- metadata_AtSC[metadata_AtSC$Condition != Plant_list[j],]
    taxo_beta_AtSC <- beta_distance_list[[i]][row.names(beta_distance_list[[i]]) %in% row.names(metadata_AtSC),row.names(beta_distance_list[[i]]) %in% row.names(metadata_AtSC)]
    taxo_adonis_AtSC <- adonis2(taxo_beta_AtSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_AtSC, method="bray", permutations=999)
    AtSC_taxo_df <- cbind(taxo_adonis_AtSC,Variable = rownames(taxo_adonis_AtSC), Subset = "AtSC", Test = "Taxonomy", Rank=Taxo_order[i])
    
    metadata_HvSC <- metadata[metadata$Inoculum == "HvSC",]
    metadata_HvSC <- metadata_HvSC[metadata_HvSC$Condition != Plant_list[j],]
    taxo_beta_HvSC <- beta_distance_list[[i]][row.names(beta_distance_list[[i]]) %in% row.names(metadata_HvSC),row.names(beta_distance_list[[i]]) %in% row.names(metadata_HvSC)]
    taxo_adonis_HvSC <- adonis2(taxo_beta_HvSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_HvSC, method="bray", permutations=999)
    HvSC_taxo_df <- cbind(taxo_adonis_HvSC,Variable = rownames(taxo_adonis_HvSC), Subset = "HvSC", Test = "Taxonomy", Rank=Taxo_order[i])
    
    metadata_LjSC <- metadata[metadata$Inoculum == "LjSC",]
    metadata_LjSC <- metadata_LjSC[metadata_LjSC$Condition != Plant_list[j],]
    taxo_beta_LjSC <- beta_distance_list[[i]][row.names(beta_distance_list[[i]]) %in% row.names(metadata_LjSC),row.names(beta_distance_list[[i]]) %in% row.names(metadata_LjSC)]
    taxo_adonis_LjSC <- adonis2(taxo_beta_LjSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_LjSC, method="bray", permutations=999)
    LjSC_taxo_df <- cbind(taxo_adonis_LjSC,Variable = rownames(taxo_adonis_LjSC), Subset = "LjSC", Test = "Taxonomy", Rank=Taxo_order[i])
    
    metadata_SSC <- metadata[metadata$Inoculum == "SSC",]
    metadata_SSC <- metadata_SSC[metadata_SSC$Condition != Plant_list[j],]
    taxo_beta_SSC <- beta_distance_list[[i]][row.names(beta_distance_list[[i]]) %in% row.names(metadata_SSC),row.names(beta_distance_list[[i]]) %in% row.names(metadata_SSC)]
    taxo_adonis_SSC <- adonis2(taxo_beta_SSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_SSC, method="bray", permutations=999)
    SSC_taxo_df <- cbind(taxo_adonis_SSC,Variable = rownames(taxo_adonis_SSC), Subset = "SSC", Test = "Taxonomy", Rank=Taxo_order[i])
    
    Combined_df_Syncom=rbind(AtSC_taxo_df, HvSC_taxo_df,LjSC_taxo_df,SSC_taxo_df)
    Sub_comb_syncom <- Combined_df_Syncom[!(Combined_df_Syncom$Variable %in% c("Total", "Residual", "Condition:Nutrient", "Condition:Experiment")), ]
    comb_syncom_list[[i]]=Sub_comb_syncom
  }
  
  for (i in 1:5) {
    #Per SynCom
    metadata_AtSC <- metadata[metadata$Inoculum == "AtSC",]
    metadata_AtSC <- metadata_AtSC[metadata_AtSC$Condition != Plant_list[j],]
    taxo_beta_AtSC <- beta_distance_list_dom[[i]][row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_AtSC),row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_AtSC)]
    taxo_adonis_AtSC <- adonis2(taxo_beta_AtSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_AtSC, method="bray", permutations=999)
    AtSC_taxo_df <- cbind(taxo_adonis_AtSC,Variable = rownames(taxo_adonis_AtSC), Subset = "AtSC", Test = "Taxonomy", Rank=Taxo_order[i])
    
    metadata_HvSC <- metadata[metadata$Inoculum == "HvSC",]
    metadata_HvSC <- metadata_HvSC[metadata_HvSC$Condition != Plant_list[j],]
    taxo_beta_HvSC <- beta_distance_list_dom[[i]][row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_HvSC),row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_HvSC)]
    taxo_adonis_HvSC <- adonis2(taxo_beta_HvSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_HvSC, method="bray", permutations=999)
    HvSC_taxo_df <- cbind(taxo_adonis_HvSC,Variable = rownames(taxo_adonis_HvSC), Subset = "HvSC", Test = "Taxonomy", Rank=Taxo_order[i])
    
    metadata_LjSC <- metadata[metadata$Inoculum == "LjSC",]
    metadata_LjSC <- metadata_LjSC[metadata_LjSC$Condition != Plant_list[j],]
    taxo_beta_LjSC <- beta_distance_list_dom[[i]][row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_LjSC),row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_LjSC)]
    taxo_adonis_LjSC <- adonis2(taxo_beta_LjSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_LjSC, method="bray", permutations=999)
    LjSC_taxo_df <- cbind(taxo_adonis_LjSC,Variable = rownames(taxo_adonis_LjSC), Subset = "LjSC", Test = "Taxonomy", Rank=Taxo_order[i])
    
    metadata_SSC <- metadata[metadata$Inoculum == "SSC",]
    metadata_SSC <- metadata_SSC[metadata_SSC$Condition != Plant_list[j],]
    taxo_beta_SSC <- beta_distance_list_dom[[i]][row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_SSC),row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_SSC)]
    taxo_adonis_SSC <- adonis2(taxo_beta_SSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_SSC, method="bray", permutations=999)
    SSC_taxo_df <- cbind(taxo_adonis_SSC,Variable = rownames(taxo_adonis_SSC), Subset = "SSC", Test = "Taxonomy", Rank=Taxo_order[i])
    
    Combined_df_Syncom=rbind(AtSC_taxo_df, HvSC_taxo_df,LjSC_taxo_df,SSC_taxo_df)
    Sub_comb_syncom <- Combined_df_Syncom[!(Combined_df_Syncom$Variable %in% c("Total", "Residual", "Condition:Nutrient", "Condition:Experiment")), ]
    comb_syncom_list_dom[[i]]=Sub_comb_syncom
    
  }
  
  #Per SynCom
  metadata_AtSC <- metadata[metadata$Inoculum == "AtSC",]
  metadata_AtSC <- metadata_AtSC[metadata_AtSC$Condition != Plant_list[j],]
  KO_beta_AtSC <- beta_isolate_KO[row.names(beta_isolate_KO) %in% row.names(metadata_AtSC),row.names(beta_isolate_KO) %in% row.names(metadata_AtSC)]
  KO_adonis_AtSC <- adonis2(KO_beta_AtSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_AtSC, method="bray", permutations=999)
  AtSC_KO_df <- cbind(KO_adonis_AtSC,Variable = rownames(KO_adonis_AtSC), Subset = "AtSC", Test = "Functions", Rank="Functions")
  
  metadata_HvSC <- metadata[metadata$Inoculum == "HvSC",]
  metadata_HvSC <- metadata_HvSC[metadata_HvSC$Condition != Plant_list[j],]
  KO_beta_HvSC <- beta_isolate_KO[row.names(beta_isolate_KO) %in% row.names(metadata_HvSC),row.names(beta_isolate_KO) %in% row.names(metadata_HvSC)]
  KO_adonis_HvSC <- adonis2(KO_beta_HvSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_HvSC, method="bray", permutations=999)
  HvSC_KO_df <- cbind(KO_adonis_HvSC,Variable = rownames(KO_adonis_HvSC), Subset = "HvSC", Test = "Functions", Rank="Functions")
  
  metadata_LjSC <- metadata[metadata$Inoculum == "LjSC",]
  metadata_LjSC <- metadata_LjSC[metadata_LjSC$Condition != Plant_list[j],]
  KO_beta_LjSC <- beta_isolate_KO[row.names(beta_isolate_KO) %in% row.names(metadata_LjSC),row.names(beta_isolate_KO) %in% row.names(metadata_LjSC)]
  KO_adonis_LjSC <- adonis2(KO_beta_LjSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_LjSC, method="bray", permutations=999)
  LjSC_KO_df <- cbind(KO_adonis_LjSC,Variable = rownames(KO_adonis_LjSC), Subset = "LjSC", Test = "Functions", Rank="Functions")
  
  metadata_SSC <- metadata[metadata$Inoculum == "SSC",]
  metadata_SSC <- metadata_SSC[metadata_SSC$Condition != Plant_list[j],]
  KO_beta_SSC <- beta_isolate_KO[row.names(beta_isolate_KO) %in% row.names(metadata_SSC),row.names(beta_isolate_KO) %in% row.names(metadata_SSC)]
  KO_adonis_SSC <- adonis2(KO_beta_SSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_SSC, method="bray", permutations=999)
  SSC_KO_df <- cbind(KO_adonis_SSC,Variable = rownames(KO_adonis_SSC), Subset = "SSC", Test = "Functions", Rank="Functions")
  
  #Per SynCom
  metadata_AtSC <- metadata[metadata$Inoculum == "AtSC",]
  metadata_AtSC <- metadata_AtSC[metadata_AtSC$Condition != Plant_list[j],]
  KO_beta_AtSC <- beta_isolate_KO_dom[row.names(beta_isolate_KO_dom) %in% row.names(metadata_AtSC),row.names(beta_isolate_KO_dom) %in% row.names(metadata_AtSC)]
  KO_adonis_AtSC <- adonis2(KO_beta_AtSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_AtSC, method="bray", permutations=999)
  AtSC_KO_df_dom <- cbind(KO_adonis_AtSC,Variable = rownames(KO_adonis_AtSC), Subset = "AtSC", Test = "Functions", Rank="Functions")
  
  metadata_HvSC <- metadata[metadata$Inoculum == "HvSC",]
  metadata_HvSC <- metadata_HvSC[metadata_HvSC$Condition != Plant_list[j],]
  KO_beta_HvSC <- beta_isolate_KO_dom[row.names(beta_isolate_KO_dom) %in% row.names(metadata_HvSC),row.names(beta_isolate_KO_dom) %in% row.names(metadata_HvSC)]
  KO_adonis_HvSC <- adonis2(KO_beta_HvSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_HvSC, method="bray", permutations=999)
  HvSC_KO_df_dom <- cbind(KO_adonis_HvSC,Variable = rownames(KO_adonis_HvSC), Subset = "HvSC", Test = "Functions", Rank="Functions")
  
  metadata_LjSC <- metadata[metadata$Inoculum == "LjSC",]
  metadata_LjSC <- metadata_LjSC[metadata_LjSC$Condition != Plant_list[j],]
  KO_beta_LjSC <- beta_isolate_KO_dom[row.names(beta_isolate_KO_dom) %in% row.names(metadata_LjSC),row.names(beta_isolate_KO_dom) %in% row.names(metadata_LjSC)]
  KO_adonis_LjSC <- adonis2(KO_beta_LjSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_LjSC, method="bray", permutations=999)
  LjSC_KO_df_dom <- cbind(KO_adonis_LjSC,Variable = rownames(KO_adonis_LjSC), Subset = "LjSC", Test = "Functions", Rank="Functions")
  
  metadata_SSC <- metadata[metadata$Inoculum == "SSC",]
  metadata_SSC <- metadata_SSC[metadata_SSC$Condition != Plant_list[j],]
  KO_beta_SSC <- beta_isolate_KO_dom[row.names(beta_isolate_KO_dom) %in% row.names(metadata_SSC),row.names(beta_isolate_KO_dom) %in% row.names(metadata_SSC)]
  KO_adonis_SSC <- adonis2(KO_beta_SSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_SSC, method="bray", permutations=999)
  SSC_KO_df_dom <- cbind(KO_adonis_SSC,Variable = rownames(KO_adonis_SSC), Subset = "SSC", Test = "Functions", Rank="Functions")
  
  combined_df_syncom <- do.call(rbind, comb_syncom_list)
  combined_df_syncom <- rbind(combined_df_syncom,AtSC_KO_df,HvSC_KO_df,LjSC_KO_df,SSC_KO_df)
  
  combined_df_syncom_dom <- do.call(rbind, comb_syncom_list_dom)
  combined_df_syncom_dom <- rbind(combined_df_syncom_dom,AtSC_KO_df_dom,HvSC_KO_df_dom,LjSC_KO_df_dom,SSC_KO_df_dom)
  
  combined_df_syncom <- combined_df_syncom[!(combined_df_syncom$Variable %in% c("Total", "Residual", "Condition:Nutrient", "Condition:Experiment","Experiment", "Nutrient")), ]
  combined_df_syncom$Rank=factor(combined_df_syncom$Rank, levels = c("Isolate", "Genus", "Family", "Order", "Class","Functions"))
  combined_df_syncom$Test=factor(combined_df_syncom$Test, levels = c("Taxonomy", "Functions"))
  combined_df_syncom$Dominance <- "Yes"
  
  combined_df_syncom_dom <- combined_df_syncom_dom[!(combined_df_syncom_dom$Variable %in% c("Total", "Residual", "Condition:Nutrient", "Condition:Experiment","Experiment", "Nutrient")), ]
  combined_df_syncom_dom$Rank=factor(combined_df_syncom_dom$Rank, levels = c("Isolate", "Genus", "Family", "Order", "Class","Functions"))
  combined_df_syncom_dom$Test=factor(combined_df_syncom_dom$Test, levels = c("Taxonomy", "Functions"))
  combined_df_syncom_dom$Dominance <- "No"
  
  combined_df_syncom <- rbind(combined_df_syncom,combined_df_syncom_dom)
  
  combined_df_syncom_2=subset.data.frame(x = combined_df_syncom, subset = combined_df_syncom$Rank!="Functions")
  combined_df_syncom_3=subset.data.frame(x = combined_df_syncom, subset = combined_df_syncom$Rank=="Functions")
  
  #For heatmap
  combined_df_syncom_3$Drop_out <- Plant_list[j]
  row.names(combined_df_syncom_3) <- NULL
  combined_df_syncom_4_dom <- rbind(combined_df_syncom_4_dom, combined_df_syncom_3)
  
  #For Table S4
  combined_df_syncom$Drop_out <- Plant_list[j]
  row.names(combined_df_syncom) <- NULL
  combined_df_syncom_5_dom <- rbind(combined_df_syncom_5_dom, combined_df_syncom)
  
}

write.table(combined_df_syncom_5_dom, paste(results.dir, "Table_S4_R2_no_dominators.tsv"), col.names =T, row.names =F, sep = "\t", quote =F)

#PieDonut plot
source(paste(working_directory, "PieDonutCustom_SSC_FL.R", sep = ""))

# SynCom colors
SynCom_colors <- data.frame(c("AtSC", "HvSC", "LjSC", "SSC"),c("#A3A500","#00B0F6","#00BF7D","#F8766D"))
colnames(SynCom_colors) <- c("SynCom", "Colour")

SynComs <- c("AtSC", "LjSC", "HvSC", "SSC")

hop_comp <- data.frame(matrix(NA, ncol =4))
colnames(hop_comp) <- c("SynCom", "Arabidopsis", "Barley", "Lotus")
hop_comp <- hop_comp[-1,]

for (syncom in SynComs){
  combined_df_syncom_4_dom_sub <- unique(combined_df_syncom_4_dom[combined_df_syncom_4_dom$Subset == paste(syncom),])
  combined_df_syncom_4_dom_sub_2 <- combined_df_syncom_4_dom_sub[combined_df_syncom_4_dom_sub$Dominance == "Yes",]
  
  At <- abs(combined_df_syncom_4_dom_sub_2$R2[combined_df_syncom_4_dom_sub_2$Drop_out == "None"] - combined_df_syncom_4_dom_sub_2$R2[combined_df_syncom_4_dom_sub_2$Drop_out == "Arabidopsis"])
  Hv <- abs(combined_df_syncom_4_dom_sub_2$R2[combined_df_syncom_4_dom_sub_2$Drop_out == "None"] - combined_df_syncom_4_dom_sub_2$R2[combined_df_syncom_4_dom_sub_2$Drop_out == "Barley"])
  Lj <- abs(combined_df_syncom_4_dom_sub_2$R2[combined_df_syncom_4_dom_sub_2$Drop_out == "None"] - combined_df_syncom_4_dom_sub_2$R2[combined_df_syncom_4_dom_sub_2$Drop_out == "Lotus"])
  
  hop_sub <- t(data.frame(c(paste(syncom), At, Hv, Lj)))
  row.names(hop_sub) <- NULL
  hop_comp <- rbind(hop_comp, hop_sub)
}

colnames(hop_comp) <- c("SynCom", "Arabidopsis", "Barley", "Lotus")

hop_comp$Arabidopsis <- as.numeric(hop_comp$Arabidopsis)
hop_comp$Barley <- as.numeric(hop_comp$Barley)
hop_comp$Lotus <- as.numeric(hop_comp$Lotus)

hop_comp_2 <- melt(hop_comp)

hop_comp_2$Rel <- round(hop_comp_2$value/sum(hop_comp_2$value)*10000,0)

hop_comp_2$Combination <- paste(hop_comp_2$SynCom, hop_comp_2$variable,sep = "_")

pie_data_2 <- data.frame()

for (combi in hop_comp_2$Combination){
  new <- hop_comp_2$Rel[hop_comp_2$Combination == paste(combi)]
  syncom <- hop_comp_2$SynCom[hop_comp_2$Combination == paste(combi)]
  plant <- hop_comp_2$variable[hop_comp_2$Combination == paste(combi)]
  
  for (i in 1:new){
    pie_data <- data.frame(paste(syncom), paste(plant))
    pie_data_2 <- rbind(pie_data_2, pie_data)
  }
}

colnames(pie_data_2) <- c("Inoculum", "Host")

print(PieDonutCustom_SSC(pie_data_2,aes(pies=Inoculum,donuts=Host),showRatioThreshold = 0.001, r0 = getOption("PieDonut.r0", 0.3), r1 = getOption("PieDonut.r1", 0.7), r2 = getOption("PieDonut.r2", 1.1), color = "gray50"))

pdf(paste(results.dir,"Figure_3a_PieDonut_with_dom.pdf", sep=""), width=6, height=6)
print(PieDonutCustom_SSC(pie_data_2,aes(pies=SynCom,donuts=Host),showRatioThreshold = 0.001, r0 = getOption("PieDonut.r0", 0.3), r1 = getOption("PieDonut.r1", 0.7), r2 = getOption("PieDonut.r2", 1.1), color = "gray50"))
dev.off()

#Figure 3b - piedonut without dominances
SynComs <- c("AtSC", "LjSC", "HvSC", "SSC")

hop_comp_B <- data.frame(matrix(NA, ncol =4))
colnames(hop_comp_B) <- c("SynCom", "Arabidopsis", "Barley", "Lotus")
hop_comp_B <- hop_comp_B[-1,]

for (syncom in SynComs){
  combined_df_syncom_4_dom_sub <- combined_df_syncom_4_dom[combined_df_syncom_4_dom$Subset == paste(syncom),]
  combined_df_syncom_4_dom_sub_2 <- combined_df_syncom_4_dom_sub[combined_df_syncom_4_dom_sub$Dominance == "No",]
  
  At <- abs(combined_df_syncom_4_dom_sub_2$R2[combined_df_syncom_4_dom_sub_2$Drop_out == "None"] - combined_df_syncom_4_dom_sub_2$R2[combined_df_syncom_4_dom_sub_2$Drop_out == "Arabidopsis"])
  Hv <- abs(combined_df_syncom_4_dom_sub_2$R2[combined_df_syncom_4_dom_sub_2$Drop_out == "None"] - combined_df_syncom_4_dom_sub_2$R2[combined_df_syncom_4_dom_sub_2$Drop_out == "Barley"])
  Lj <- abs(combined_df_syncom_4_dom_sub_2$R2[combined_df_syncom_4_dom_sub_2$Drop_out == "None"] - combined_df_syncom_4_dom_sub_2$R2[combined_df_syncom_4_dom_sub_2$Drop_out == "Lotus"])
  
  hop_sub_B <- t(data.frame(c(paste(syncom), At, Hv, Lj)))
  row.names(hop_sub_B) <- NULL
  hop_comp_B <- rbind(hop_comp_B, hop_sub_B)
}

colnames(hop_comp_B) <- c("SynCom", "Arabidopsis", "Barley", "Lotus")

hop_comp_B$Arabidopsis <- as.numeric(hop_comp_B$Arabidopsis)
hop_comp_B$Barley <- as.numeric(hop_comp_B$Barley)
hop_comp_B$Lotus <- as.numeric(hop_comp_B$Lotus)

hop_comp_2_B <- melt(hop_comp_B)

hop_comp_2_B$Rel <- round(hop_comp_2_B$value/sum(hop_comp_2_B$value)*10000,0)

hop_comp_2_B$Combination <- paste(hop_comp_2_B$SynCom, hop_comp_2_B$variable,sep = "_")

pie_data_2_B <- data.frame()

for (combi in hop_comp_2_B$Combination){
  new <- hop_comp_2_B$Rel[hop_comp_2_B$Combination == paste(combi)]
  syncom <- hop_comp_2_B$SynCom[hop_comp_2_B$Combination == paste(combi)]
  plant <- hop_comp_2$variable[hop_comp_2_B$Combination == paste(combi)]
  
  for (i in 1:new){
    pie_data_B <- data.frame(paste(syncom), paste(plant))
    pie_data_2_B <- rbind(pie_data_2_B, pie_data_B)
  }
}

colnames(pie_data_2_B) <- c("Inoculum", "Host")

print(PieDonutCustom_SSC(pie_data_2_B,aes(pies=Inoculum,donuts=Host),showRatioThreshold = 0.001, r0 = getOption("PieDonut.r0", 0.3), r1 = getOption("PieDonut.r1", 0.7), r2 = getOption("PieDonut.r2", 1.1), color = "gray50"))

pdf(paste(results.dir,"Figure_3b_PieDonut_no_dom.pdf", sep=""), width=6, height=6)
print(PieDonutCustom_SSC(pie_data_2_B,aes(pies=Inoculum,donuts=Host),showRatioThreshold = 0.001, r0 = getOption("PieDonut.r0", 0.3), r1 = getOption("PieDonut.r1", 0.7), r2 = getOption("PieDonut.r2", 1.1), color = "gray50"))
dev.off()
