library("dplyr") #Version 1.1.2
library("phyloseq") #Version 1.44.0
library("reshape2") #Version 1.4.4
library("tidyr") #Version 1.3.0
library("ggplot2") #Version 3.4.2

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Figure S9 - Relative abundances dominators & nodulators in nodules =====

#Figure S9 Dominance by rhizobacter and mesorhizobium
#OTU TABLE
norm_SSC=read.table(paste(working_directory,"Isolate_tables/Original/SSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
round_SSC=floor(x = norm_SSC)

#Taxonomy TABLE
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
colnames(samples_df)[6]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparation
#Set the OTU, TAX and sample data for making phyloseq object
OTU = otu_table(as.matrix(round_SSC),taxa_are_rows = TRUE)
# TAX = tax_table(tax_mat)
TAX = tax_table(as.matrix(tax_df_2))

#Sample subsetting
samples_df_sub <- subset(samples_df, samples_df$Compartment != "RZ")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "AM")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "NOD")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "Input")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Condition != "NP")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Inoculum != "AtSC")

samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

samples <- sample_data(samples_df_sub_2)

phylo_sub = phyloseq(OTU,TAX, samples)

subsetted_table <- otu_table(phylo_sub)
subsetted_table_long <- melt(subsetted_table)

Hank_the_normalizer <- function(df,group,amount){
  df_2 <- df %>% dplyr::group_by_at(group) %>% dplyr::summarise(total=sum(.data[[amount]]))
  df_3 <- df_2$total
  names(df_3) <- df_2[[group]]
  df$total <- df_3[as.character(df[[group]])]
  df$Rel <- df[[amount]] / df$total
  return(df)
}

subsetted_table_long_2 <- Hank_the_normalizer(subsetted_table_long,"Var2","value")
subsetted_table_long_2$value[subsetted_table_long_2$Rel < 0.0005] <- 0
subsetted_table_long_3 <- subsetted_table_long_2[1:3]
data_wide <- spread(subsetted_table_long_3, Var2, value)
row.names(data_wide) <- data_wide$Var1
data_wide_2 <- data_wide %>% dplyr::select (-Var1)

OTU = otu_table(as.matrix(data_wide_2),taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(tax_df_2))
samples <- sample_data(samples_df_sub_2)
phylo_sub = phyloseq(OTU,TAX, samples)

# Transform to relative abundance
phylo_sub_RA <- microbiome::transform(phylo_sub, "compositional")

# Melt phyloseq object into a dataframe
data_melt <- psmelt(phylo_sub_RA)

# Create a new OTU column where non-top5 are grouped
data_melt$OTU_grouped <- ifelse(
  data_melt$OTU %in% c("P2_G4", "LjNodule214", "P1_H10", "P2_A12", "P2_D6"), 
  data_melt$OTU, 
  ifelse(!is.na(data_melt$genus) & data_melt$genus == "Rhizobacter", "Other Rhizobacter", 
         ifelse(!is.na(data_melt$genus) & data_melt$genus == "Mesorhizobium", "Other Mesorhizobium", "Other strains"))
)

# Ensure OTU_grouped is a factor with correct legend order
otu_order <- c("P2_G4","Other Rhizobacter", "Other strains", "LjNodule214", "P1_H10", "P2_A12", "P2_D6","Other Mesorhizobium")
data_melt$OTU_grouped <- factor(data_melt$OTU_grouped, levels = otu_order)

# Define colors
otu_colors <- c(
  "P2_G4" = "#C85AC8",  # Purple
  "LjNodule214" = "#66E1D0",  # Blue
  "P1_H10" = "#00C1C8",       # Orange
  "P2_A12" = "#00AA95",       # Green
  "P2_D6" = "#00C18C",        # Red
  "Other Rhizobacter" = "#F096F0",  # Light purple
  "Other Mesorhizobium" = "#00773E",  # Greenish
  "Other strains" = "lightgray"  # Generic category
)

# Stacked barplot with grouped OTUs
bar <- ggplot(data_melt, aes(fill=OTU_grouped, y=Abundance, x=sample_id)) + 
  geom_bar(position="stack", stat="identity", colour = "darkgray", linewidth = 0.01) +  
  scale_fill_manual(values = otu_colors) +  
  ggtitle("Nodulators and Rhizobacter Root colonization profiles") + 
  theme_classic() +
  labs(y = "Relative abundance", fill = "Isolates") +  
  theme(
    plot.title = element_text(hjust = 0.5, size = 20), 
    axis.text.x = element_blank(),
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 12)
  )

# Add facet wrap by inoculum if needed
bar2 <- bar + facet_wrap(~Inoculum+Condition+Nutrient, scales = "free_x", nrow = 1) +
  theme(strip.text.x = element_text(size = 12))
bar2

pdf(paste(results.dir,"Figure_S9_dominator_taxonomic_profile.pdf", sep=""), width=15, height=4.5)
print(bar2)
dev.off()
