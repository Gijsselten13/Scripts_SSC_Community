library("phyloseq") #Version 1.44.0
library("dplyr") #Version 1.1.2
library("vegan") #Version 2.6-4
library("ggplot2") #Version 3.4.2

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Figure S5 - Rarefaction curves =====

# OTU TABLE 
norm_SSC =read.table(paste(working_directory,"Isolate_tables/Original/SSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
# Taxonomy TABLE 
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
# Samples TABLE 
samples_df = read.table(paste(working_directory,"SSC_R2_metadata.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
colnames(samples_df)[6]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Set the OTU, TAX and sample data for making phyloseq object
OTU = otu_table(as.matrix(norm_SSC),taxa_are_rows = TRUE)
# TAX = tax_table(tax_mat)
TAX = tax_table(as.matrix(tax_df_2))

cond="ES"

#Subset for root (endosphere) samples
samples_df_sub <- subset(samples_df, samples_df$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

samples_sub = sample_data(samples_df_sub_2)

phylo_sub = phyloseq(OTU,TAX, samples_sub)

racur_data=rarecurve(x =floor(t(phylo_sub@otu_table@.Data)), step = 250, label = F, tidy = T)
racur_data=merge.data.frame(x = racur_data, y = samples_df_sub, by.x = "Site", by.y = 0)
racur_data$Inoculum <- factor(racur_data$Inoculum, levels = c("SSC","AtSC", "LjSC", "HvSC"))
racur_data$Condition <- gsub("At", "Arabidopsis", racur_data$Condition)
racur_data$Condition <- gsub("Hv", "Barley", racur_data$Condition)
racur_data$Condition <- gsub("Lj", "Lotus", racur_data$Condition)

racu=ggplot(data = racur_data, aes(x = Sample, y = Species, color = Experiment))
racu_all <- racu+geom_point( size = 0.05) + xlim(0, 50000)+ 
  xlab("No. of pseudoaligned reads") + ylab("No of isolates") +
  facet_grid(facets = c("Condition","Inoculum")) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 25)) +
  theme(
    panel.grid.major = element_line(color = "grey80"),  
    panel.grid.minor = element_line(color = "grey90")   
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))
racu_all

pdf(paste(results.dir,"Figure_S5_rarefaction_curves.pdf", sep=""), width=8, height=6)
print(racu_all)
dev.off()
