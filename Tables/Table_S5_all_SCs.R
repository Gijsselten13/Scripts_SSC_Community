library("dplyr") #Version 1.1.2
library("phyloseq") #Version 1.44.0
library("vegan") #Version 2.6-4
library("reshape2") #Version 1.4.4

#Set working directory and results folder here
working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Table S5 - R2 values at different taxonomic levels - all SCs =====
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

mean_value_phylum= mean(beta_phylum)
mean_value_class= mean(beta_class)
mean_value_order= mean(beta_order)
mean_value_family= mean(beta_family)
mean_value_genus= mean(beta_genus)
mean_value_isolate= mean(beta_isolate)

Taxo_order=c("Isolate","Genus", "Family", "Order", "Class")
beta_distance_list=list(beta_isolate, beta_genus, beta_family, beta_order, beta_phylum)
Adonis_list=list()
Plot_list=list()
Plot_list_2=list()

for (i in 1:length(beta_distance_list)) {
  
  Bray_curtis_df=beta_distance_list[[i]]
  
  #Make PCoA plot for Bray Curtis Distance matrix
  pcoa = cmdscale(Bray_curtis_df, k=3, eig=T)
  points = as.data.frame(pcoa$points)
  colnames(points) = c("x", "y", "z") 
  eig = pcoa$eig
  
  points = merge(points,samples_df_sub_2, by = "row.names")
  rownames(points) <- points$Row.names
  points <- points %>% dplyr::select (-Row.names)
  
  points$Condition <- factor(points$Condition, levels = c("At","Hv", "Lj"))
  points$Inoculum <- factor(points$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC","NS"))
  points$Nutrient <- factor(points$Nutrient, levels = c("low", "high"))
  points$Experiment  <- factor(points$Experiment, levels = c("R1", "R2"))
  
  metadata=points[,-c(1,2,3)]
  
  #  Run adonis PERMANOVA test
  set.seed(1)
  SSC_bray_adonis <- adonis2(beta_distance_list[[i]] ~ Inoculum*Condition*Nutrient*Experiment, data=metadata, method="bray", permutations=999)
  
  Adonis_list[[i]]=SSC_bray_adonis
  
}

#otu table
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
points$Inoculum <- factor(points$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC","NS"))
points$Nutrient <- factor(points$Nutrient, levels = c("low", "high"))
points$Experiment  <- factor(points$Experiment, levels = c("R1", "R2"))

metadata=points[,-c(1,2,3)]

set.seed(1)
SSC_bray_KO_adonis <- adonis2(beta_isolate_KO ~ Inoculum*Condition*Nutrient*Experiment, data=metadata, method="bray", permutations=999)

#You need to run the script of figure 5 before running this one
# Adonis results from KO functions
Adonis_list[[length(beta_distance_list)+1]]=SSC_bray_KO_adonis

# Make a dataframe out of all these results with R2 values
R2_values <- matrix(nrow = length(Adonis_list), ncol = 5)
for (i in 1:length(Adonis_list)) {
  R2_values[i, ] <- Adonis_list[[i]]$R2[1:5]
}

# Row and column names
Taxo_order <- c("Isolate", "Genus", "Family", "Order", "Phylum","KO_functions")
col_names <- c("Syncom", "Plant", "Nutrient", "Experiment", paste("Syncom:", "Plant", sep = ""))

# Create the dataframe with R2 values
df <- data.frame(Taxo_order, R2_values)
colnames(df)[-1] <- col_names
df

# Melt the dataframe
melted_df <- melt(df, id.vars = "Taxo_order")
melted_df <- melted_df %>%
  mutate(Test = ifelse(Taxo_order == "KO_functions", "Function", "Taxonomy"))

colnames(melted_df) <- c("Taxo_order", "Variable", "R2", "Test")
melted_df$Variable <- as.character(melted_df$Variable)

melted_df$Variable[melted_df$Variable == "Syncom"] <- "Inoculum"
melted_df$Variable[melted_df$Variable == "Plant"] <- "Host"
melted_df$Variable[melted_df$Variable == "Syncom:Plant"] <- "Inoculum:Host"

melted_df$Taxo_order=factor(melted_df$Taxo_order, levels = c("Isolate", "Genus", "Family", "Order", "Class", "Phylum","KO_functions"))
melted_df$Variable=factor(melted_df$Variable, levels = c("Inoculum", "Host","Inoculum:Host", "Nutrient", "Experiment"))
melted_df$Test=factor(melted_df$Test, levels = c("Taxonomy", "Function"))

# Create extra dfs for the line plot

# Taxonomy
melted_df_2=subset.data.frame(x = melted_df, subset = melted_df$Taxo_order!="KO_functions")
# Functions for dashed line
melted_df_3=subset.data.frame(x = melted_df, subset = melted_df$Taxo_order=="KO_functions")

melted_df_2$Dominance <- "Yes"
melted_df_3$Dominance <- "Yes"

#ADDING THE DATA WITHOUT THE NODULATORS & RHIZOBACTER
#otu table
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
# row.names(samples_df_sub_2) <- sub('HL_orig', 'HL',row.names(samples_df_sub_2))

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

mean_value_phylum_dom = mean(beta_phylum_dom)
mean_value_class_dom = mean(beta_class_dom)
mean_value_order_dom = mean(beta_order_dom)
mean_value_family_dom = mean(beta_family_dom)
mean_value_genus_dom = mean(beta_genus_dom)
mean_value_isolate_dom = mean(beta_isolate_dom)

Taxo_order=c("Isolate","Genus", "Family", "Order", "Class")
beta_distance_list_dom=list(beta_isolate_dom, beta_genus_dom, beta_family_dom, beta_order_dom, beta_phylum_dom)
Adonis_list_dom=list()
Plot_list_dom=list()
Plot_list_2_dom=list()

for (i in 1:length(beta_distance_list_dom)) {
  
  Bray_curtis_df=beta_distance_list_dom[[i]]
  
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
  
  metadata=points[,-c(1,2,3)]
  
  #  Run adonis PERMANOVA test
  set.seed(1)
  SSC_bray_adonis <- adonis2(beta_distance_list_dom[[i]] ~ Inoculum*Condition*Nutrient*Experiment, data=metadata, method="bray", permutations=999)
  
  Adonis_list_dom[[i]]=SSC_bray_adonis
}

#otu table
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

beta_isolate_KO_dom <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

Bray_curtis_df=beta_isolate_KO_dom

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

metadata=points[,-c(1,2,3)]

set.seed(1)
SSC_bray_KO_adonis <- adonis2(beta_isolate_KO_dom ~ Inoculum*Condition*Nutrient*Experiment, data=metadata, method="bray", permutations=999)

# Adonis results from KO functions
Adonis_list_dom[[length(beta_distance_list_dom)+1]]=SSC_bray_KO_adonis

# Make a dataframe out of all these results with R2 values
R2_values_dom <- matrix(nrow = length(Adonis_list_dom), ncol = 5)
for (i in 1:length(Adonis_list_dom)) {
  R2_values_dom[i, ] <- Adonis_list_dom[[i]]$R2[1:5]
}

# Row and column names
Taxo_order <- c("Isolate", "Genus", "Family", "Order", "Phylum","KO_functions")
col_names <- c("Syncom", "Plant", "Nutrient", "Experiment", paste("Syncom:", "Plant", sep = ""))

# Create the dataframe with R2 values
df <- data.frame(Taxo_order, R2_values_dom)
colnames(df)[-1] <- col_names
df

# Melt the dataframe
melted_df <- melt(df, id.vars = "Taxo_order")
melted_df <- melted_df %>%
  mutate(Test = ifelse(Taxo_order == "KO_functions", "Function", "Taxonomy"))

colnames(melted_df) <- c("Taxo_order", "Variable", "R2", "Test")
melted_df$Variable <- as.character(melted_df$Variable)

melted_df$Variable[melted_df$Variable == "Syncom"] <- "Inoculum"
melted_df$Variable[melted_df$Variable == "Plant"] <- "Host"
melted_df$Variable[melted_df$Variable == "Syncom:Plant"] <- "Inoculum:Host"

melted_df$Taxo_order=factor(melted_df$Taxo_order, levels = c("Isolate", "Genus", "Family", "Order", "Class", "Phylum","KO_functions"))
melted_df$Variable=factor(melted_df$Variable, levels = c("Inoculum", "Host","Inoculum:Host", "Nutrient", "Experiment"))
melted_df$Test=factor(melted_df$Test, levels = c("Taxonomy", "Function"))

write.table(melted_df, paste(results.dir,"Table_S5_R2_values_all_SCs.txt", sep =""), sep = "\t", quote =F,col.names =T, row.names =F)
