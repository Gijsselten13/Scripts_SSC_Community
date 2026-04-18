library("dplyr") #Version 1.1.2
library("phyloseq") #Version 1.44.0
library("vegan") #Version 2.6-4
library("ggplot2") #Version 3.4.2
library("ggpubr") #Version 0.6.0
library("ggrepel") #Version 0.9.3

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Figure S12 and Table S3 - R2 Simulations =====
R2_values_2 <- read.table(paste(working_directory,"R2_values_genus.txt", sep = ""), sep = "\t", header = T)

#Add our own R2 value to the simulated data
#otu table
norm_SSC =read.table(paste(working_directory,"Isolate_tables/Original/SSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
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
ps_genus <- phyloseq::tax_glom(phylo_sub_RA, "Genus")
phyloseq::taxa_names(ps_genus) <- phyloseq::tax_table(ps_genus)[, "Genus"]

#Bray Curtis distance matrix
beta_genus <- as.matrix(vegdist(t(ps_genus@otu_table@.Data), method = "bray", diag = T))

#Make PCoA plot for Bray Curtis Distance matrix
pcoa = cmdscale(beta_genus, k=3, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z") 
eig = pcoa$eig

points = merge(points,samples_df_sub_2, by = "row.names")
rownames(points) <- points$Row.names
points <- points %>% dplyr::select (-Row.names)
metadata=points[,-c(1,2,3)]

set.seed(1)
SSC_bray_adonis <- adonis2(beta_genus ~ Inoculum, data=metadata, method="bray", permutations=999)
SSC_Tax_SC_value <- SSC_bray_adonis$R2[1]

SSC_bray_adonis <- adonis2(beta_genus ~ Condition, data=metadata, method="bray", permutations=999)
SSC_Tax_PL_value <- SSC_bray_adonis$R2[1]

#KO table
KO_SSC =read.table(paste(working_directory,"KO_tables/Original/SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

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

#Make PCoA plot for Bray Curtis Distance matrix
pcoa = cmdscale(beta_isolate_KO, k=3, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z") 
eig = pcoa$eig

points = merge(points,samples_df_sub_2, by = "row.names")
rownames(points) <- points$Row.names
points <- points %>% dplyr::select (-Row.names)

metadata=points[,-c(1,2,3)]

set.seed(1)
SSC_bray_KO_adonis <- adonis2(beta_isolate_KO ~ Inoculum, data=metadata, method="bray", permutations=999)
SSC_KO_SC_value <- SSC_bray_KO_adonis$R2[1]

SSC_bray_KO_adonis <- adonis2(beta_isolate_KO ~ Condition, data=metadata, method="bray", permutations=999)
SSC_KO_PL_value <- SSC_bray_KO_adonis$R2[1]

SSC_values_SC <- c("Inoculum", SSC_Tax_SC_value,SSC_KO_SC_value)
SSC_values_PL <- c("Host", SSC_Tax_PL_value,SSC_KO_PL_value)

#otu table - no nodulators
norm_SSC =read.table(paste(working_directory,"Isolate_tables/No_nodulators/SSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
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
ps_genus <- phyloseq::tax_glom(phylo_sub_RA, "Genus")
phyloseq::taxa_names(ps_genus) <- phyloseq::tax_table(ps_genus)[, "Genus"]

#Bray Curtis distance matrix
beta_genus <- as.matrix(vegdist(t(ps_genus@otu_table@.Data), method = "bray", diag = T))

#Make PCoA plot for Bray Curtis Distance matrix
pcoa = cmdscale(beta_genus, k=3, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z") 
eig = pcoa$eig

points = merge(points,samples_df_sub_2, by = "row.names")
rownames(points) <- points$Row.names
points <- points %>% dplyr::select (-Row.names)
metadata=points[,-c(1,2,3)]

set.seed(1)
SSC_bray_adonis <- adonis2(beta_genus ~ Inoculum, data=metadata, method="bray", permutations=999)
SSC_Tax_SC_value <- SSC_bray_adonis$R2[1]

SSC_bray_adonis <- adonis2(beta_genus ~ Condition, data=metadata, method="bray", permutations=999)
SSC_Tax_PL_value <- SSC_bray_adonis$R2[1]

#KO table
KO_SSC =read.table(paste(working_directory,"KO_tables/No_nodulators/SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

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

#Make PCoA plot for Bray Curtis Distance matrix
pcoa = cmdscale(beta_isolate_KO, k=3, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z") 
eig = pcoa$eig

points = merge(points,samples_df_sub_2, by = "row.names")
rownames(points) <- points$Row.names
points <- points %>% dplyr::select (-Row.names)

metadata=points[,-c(1,2,3)]

set.seed(1)
SSC_bray_KO_adonis <- adonis2(beta_isolate_KO ~ Inoculum, data=metadata, method="bray", permutations=999)
SSC_KO_SC_value <- SSC_bray_KO_adonis$R2[1]

SSC_bray_KO_adonis <- adonis2(beta_isolate_KO ~ Condition, data=metadata, method="bray", permutations=999)
SSC_KO_PL_value <- SSC_bray_KO_adonis$R2[1]

SSC_values_SC_nod <- c("Inoculum - no nodulators", SSC_Tax_SC_value,SSC_KO_SC_value)
SSC_values_PL_nod <- c("Host - no nodulators", SSC_Tax_PL_value,SSC_KO_PL_value)

#otu table - no Rhizobacter
norm_SSC =read.table(paste(working_directory,"Isolate_tables/No_rhizobacter/SSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
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
ps_genus <- phyloseq::tax_glom(phylo_sub_RA, "Genus")
phyloseq::taxa_names(ps_genus) <- phyloseq::tax_table(ps_genus)[, "Genus"]

#Bray Curtis distance matrix
beta_genus <- as.matrix(vegdist(t(ps_genus@otu_table@.Data), method = "bray", diag = T))

#Make PCoA plot for Bray Curtis Distance matrix
pcoa = cmdscale(beta_genus, k=3, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z") 
eig = pcoa$eig

points = merge(points,samples_df_sub_2, by = "row.names")
rownames(points) <- points$Row.names
points <- points %>% dplyr::select (-Row.names)
metadata=points[,-c(1,2,3)]

set.seed(1)
SSC_bray_adonis <- adonis2(beta_genus ~ Inoculum, data=metadata, method="bray", permutations=999)
SSC_Tax_SC_value <- SSC_bray_adonis$R2[1]

SSC_bray_adonis <- adonis2(beta_genus ~ Condition, data=metadata, method="bray", permutations=999)
SSC_Tax_PL_value <- SSC_bray_adonis$R2[1]

#KO table
KO_SSC =read.table(paste(working_directory,"KO_tables/No_rhizobacter/SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

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

#Make PCoA plot for Bray Curtis Distance matrix
pcoa = cmdscale(beta_isolate_KO, k=3, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z") 
eig = pcoa$eig

points = merge(points,samples_df_sub_2, by = "row.names")
rownames(points) <- points$Row.names
points <- points %>% dplyr::select (-Row.names)

metadata=points[,-c(1,2,3)]

set.seed(1)
SSC_bray_KO_adonis <- adonis2(beta_isolate_KO ~ Inoculum, data=metadata, method="bray", permutations=999)
SSC_KO_SC_value <- SSC_bray_KO_adonis$R2[1]

SSC_bray_KO_adonis <- adonis2(beta_isolate_KO ~ Condition, data=metadata, method="bray", permutations=999)
SSC_KO_PL_value <- SSC_bray_KO_adonis$R2[1]

SSC_values_SC_rhizo <- c("Inoculum - no ", SSC_Tax_SC_value,SSC_KO_SC_value)
SSC_values_PL_rhizo <- c("Host - no ", SSC_Tax_PL_value,SSC_KO_PL_value)

#otu table - no Rhizobacter & Nodulators
norm_SSC =read.table(paste(working_directory,"Isolate_tables/No_dominances/SSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
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
ps_genus <- phyloseq::tax_glom(phylo_sub_RA, "Genus")
phyloseq::taxa_names(ps_genus) <- phyloseq::tax_table(ps_genus)[, "Genus"]

#Bray Curtis distance matrix
beta_genus <- as.matrix(vegdist(t(ps_genus@otu_table@.Data), method = "bray", diag = T))

#Make PCoA plot for Bray Curtis Distance matrix
pcoa = cmdscale(beta_genus, k=3, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z") 
eig = pcoa$eig

points = merge(points,samples_df_sub_2, by = "row.names")
rownames(points) <- points$Row.names
points <- points %>% dplyr::select (-Row.names)
metadata=points[,-c(1,2,3)]

set.seed(1)
SSC_bray_adonis <- adonis2(beta_genus ~ Inoculum, data=metadata, method="bray", permutations=999)
SSC_Tax_SC_value <- SSC_bray_adonis$R2[1]

SSC_bray_adonis <- adonis2(beta_genus ~ Condition, data=metadata, method="bray", permutations=999)
SSC_Tax_PL_value <- SSC_bray_adonis$R2[1]

#KO table
KO_SSC =read.table(paste(working_directory,"KO_tables/No_dominances/SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

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

#Make PCoA plot for Bray Curtis Distance matrix
pcoa = cmdscale(beta_isolate_KO, k=3, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z") 
eig = pcoa$eig

points = merge(points,samples_df_sub_2, by = "row.names")
rownames(points) <- points$Row.names
points <- points %>% dplyr::select (-Row.names)

metadata=points[,-c(1,2,3)]

set.seed(1)
SSC_bray_KO_adonis <- adonis2(beta_isolate_KO ~ Inoculum, data=metadata, method="bray", permutations=999)
SSC_KO_SC_value <- SSC_bray_KO_adonis$R2[1]

SSC_bray_KO_adonis <- adonis2(beta_isolate_KO ~ Condition, data=metadata, method="bray", permutations=999)
SSC_KO_PL_value <- SSC_bray_KO_adonis$R2[1]

SSC_values_SC_rhizo_nod <- c("Inoculum - no dominators", SSC_Tax_SC_value,SSC_KO_SC_value)
SSC_values_PL_rhizo_nod <- c("Host - no dominators", SSC_Tax_PL_value,SSC_KO_PL_value)

R2_values_2 <- rbind(R2_values_2, SSC_values_SC,SSC_values_PL,SSC_values_SC_nod,SSC_values_PL_nod, SSC_values_SC_rhizo, SSC_values_PL_rhizo, SSC_values_SC_rhizo_nod,SSC_values_PL_rhizo_nod)
R2_values_2$Tax_R2 <- as.numeric(R2_values_2$Tax_R2)
R2_values_2$Func_R2 <- as.numeric(R2_values_2$Func_R2)
R2_values_2$Simulation <- "Simulation"
R2_values_2$Simulation[R2_values_2$run == "Inoculum"] <- "SSC_SC"
R2_values_2$Simulation[R2_values_2$run == "Inoculum - no nodulators"] <- "SSC_no_nod_SC"
R2_values_2$Simulation[R2_values_2$run == "Host"] <- "SSC_PL"
R2_values_2$Simulation[R2_values_2$run == "Host - no nodulators"] <- "SSC_no_nod_PL"
R2_values_2$Simulation[R2_values_2$run == "Inoculum - no "] <- "SSC_SC_no_rhizo"
R2_values_2$Simulation[R2_values_2$run == "Inoculum - no dominators"] <- "SSC_no_rhizo_and_nod_SC"
R2_values_2$Simulation[R2_values_2$run == "Host - no "] <- "SSC_PL_no_rhizo"
R2_values_2$Simulation[R2_values_2$run == "Host - no dominators"] <- "SSC_no_rhizo_and_nod_PL"

R2_values_2$Simulation_2 <- "Simulation"
R2_values_2$Simulation_2[R2_values_2$run == "Inoculum"] <- "SSC_SC"
R2_values_2$Simulation_2[R2_values_2$run == "Inoculum - no nodulators"] <- "SSC_SC"
R2_values_2$Simulation_2[R2_values_2$run == "Host"] <- "SSC_PL"
R2_values_2$Simulation_2[R2_values_2$run == "Host - no nodulators"] <- "SSC_PL"
R2_values_2$Simulation_2[R2_values_2$run == "Inoculum - no "] <- "SSC_SC"
R2_values_2$Simulation_2[R2_values_2$run == "Inoculum - no dominators"] <- "SSC_SC"
R2_values_2$Simulation_2[R2_values_2$run == "Host - no "] <- "SSC_PL"
R2_values_2$Simulation_2[R2_values_2$run == "Host - no dominators"] <- "SSC_PL"

R2_values_2$Simulation_2 <- factor(R2_values_2$Simulation_2, levels = c("Simulation", "SSC_SC", "SSC_PL"))

checkies <- c(3,6,9)

R2_values_3 <- R2_values_2[R2_values_2$kmeans %in% checkies, ]
R2_values_4 <- unique(R2_values_2[R2_values_2$Simulation != "Simulation", ])

R2_values_5 <- rbind(R2_values_3, R2_values_4)

data_stat <- data.frame()

groups <- R2_values_5$run[3001:3008]

for (group in groups){
  model <- lm(Tax_R2 ~ Func_R2, data = R2_values_5)
  R2_values_5$residuals <- residuals(model)
  std_error <- summary(model)$sigma
  specific_point <- R2_values_5[R2_values_5$run == paste(group), ]
  predicted_value <- predict(model, newdata = specific_point)
  residual <- specific_point$Tax_R2 - predicted_value
  standardized_residual <- residual / std_error
  p_value <- 2 * (1 - pnorm(abs(standardized_residual)) ) # Two-tailed test
  hop <- data.frame(t(data.frame(c(paste(group), p_value))))
  row.names(hop) <- NULL
  colnames(hop) <- c("Group","pvalue")
  data_stat <- rbind(data_stat,hop )
}

#Statistical deviations of trend
write.table(data_stat, paste(working_directory, "Table_S3_R2_simulation_stats.txt", sep =""), sep= "\t", quote =F, row.names =F, col.names=T)

R2_correlations <- ggscatter(R2_values_5, x="Tax_R2", y="Func_R2", color = "Simulation_2", conf.int = F,alpha = 0.7,
                             palette = c(Simulation_2 = "black", SSC_SC = "#F8766D",SSC_PL = "#80b006")) +
  stat_cor(label.x = 0.05, label.y = 0.4) +
  geom_smooth(method=lm, level = 0.99999999999999) +
  geom_text_repel(aes(label=ifelse(R2_values_5$Simulation_2 == "SSC_SC",yes = as.character(R2_values_5$run), '')),size=4,max.overlaps = Inf) +
  geom_text_repel(aes(label=ifelse(R2_values_5$Simulation_2 == "SSC_PL",yes = as.character(R2_values_5$run), '')),size=4, max.overlaps = Inf) +
  ggtitle("R2 Correlation between taxonomy and functionality") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  ylab(expression(paste("Functionality R"^2, sep =""))) + 
  xlab(expression(paste("Taxonomy R"^2, sep = ""))) + 
  theme(axis.text.x = element_text(size = 14), axis.title = element_text(size = 18), axis.text.y = element_text(size=14), plot.title = element_text(size=24)) +
  theme(legend.position = "none") 
R2_correlations

pdf(paste(results.dir,"Figure_S12_R2_correlation_plot_simulation_2.pdf", sep=""), width=10, height=7)
print(R2_correlations)
dev.off()
