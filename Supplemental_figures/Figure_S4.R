library("dplyr") #Version 1.1.2
library("stringr") #Version 1.5.0
library("rstatix") #Version 0.7.2
library("multcompView") #Version 0.1-9
library("plyr") #Version 1.8.8
library("ggplot2") #Version 3.4.2
library("ggpubr") #Version 0.6.0
library("phyloseq") #Version 1.44.0
library("reshape2") #Version 1.4.4

working_directory <- ""
results.dir <- paste(working_directory,"results/", sep = "")

###Figure S4 - Nodule Numbers =====
importdat_nod <- read.table(file = paste(working_directory,"SSC_nodules.txt", sep = ""), sep = "\t", header=TRUE, dec = ".")

# Manipulate dataframe so it can be processed by ggplot2 package
nod_good_table=reshape2::melt(importdat_nod, measure.vars = c("Pink","White","Total"), value.name = "number", id.vars = c("condition","replicate"), variable.name = "type")
nod_good_table=mutate(nod_good_table,inoculum=word(condition,start = 2, sep=fixed("_")))
nod_good_table=mutate(nod_good_table,nutrient=word(condition,start = -1, sep=fixed("_")))

nod_good_table_R1 <- nod_good_table[grepl(pattern = "R1_",nod_good_table$condition),]
new <- str_split(nod_good_table_R1$condition, pattern = "_")
nod_good_table_R1$inoculum <- data.table::transpose(new)[[3]]

nod_good_table_R2 <- nod_good_table[!grepl(pattern = "R1_",nod_good_table$condition),]
nod_good_table_2 <- rbind(nod_good_table_R2,nod_good_table_R1)

nod_good_table_2$inoc_nut=paste0(nod_good_table_2$inoculum,"_",nod_good_table_2$nutrient)

# Remove NA values because of different number of samples being tested
nod_good_table_2=na.omit(nod_good_table_2)
nod_good_table_2$inoc_nut <- factor(nod_good_table_2$inoc_nut, levels = c("AtSC_Low", "HvSC_Low", "LjSC_Low", "SSC_Low","NS_Low", "AtSC_High","HvSC_High", "LjSC_High", "SSC_High","NS_High"))
nod_good_table_2$inoculum <- factor(nod_good_table_2$inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC","NS"))

# Visualize the data using box plots. Plot weight by groups.
nod_good_table_1=subset(x = nod_good_table_2, nod_good_table_2$type==c("Total"))

#  Kruskal wallis test followed by Dunn post hoc test
res.kruskal <- nod_good_table_1 %>% kruskal_test(number ~ inoc_nut)

# Dunn posthoc Pairwise comparisons
pwc <- nod_good_table_1 %>% 
  dunn_test(number ~ inoc_nut, p.adjust.method = "bonferroni") 

# Generating letters for kruskal+dunn pairwise comparisons
tukey_values= data.frame()
fit=aov(data=nod_good_table_1,number ~ inoc_nut)
anova(fit)
res=TukeyHSD(fit)
res[[1]][,4]=pwc$p.adj
Tukey.levels <- res[[1]][,4]
Labels_pairwise <- multcompLetters(Tukey.levels)['Letters']
inoc_nut <- names(Labels_pairwise[['Letters']])

boxplot.df <- ddply(nod_good_table_1, .(inoc_nut), function (x) max(fivenum(x$number)+0.04*(max(x$number))))

# Create a data frame out of the factor levels and Tukey's homogenous group letters
plot.levels <- data.frame(inoc_nut, labels = Labels_pairwise[['Letters']],
                          stringsAsFactors = FALSE)

# Merge it with the labels
labels.df <- merge(plot.levels, boxplot.df, by= "inoc_nut", sort = FALSE)

nod_good_table_1$inoc_nut <- factor(nod_good_table_1$inoc_nut, levels = c("AtSC_Low", "HvSC_Low", "LjSC_Low", "SSC_Low","NS_Low", "AtSC_High","HvSC_High", "LjSC_High", "SSC_High","NS_High"))

nod_total <- ggplot(nod_good_table_1, aes(x = inoc_nut, y = number)) +
  geom_boxplot(notch = FALSE, size=1, aes(color = inoculum),outlier.shape = NA) +
  theme_pubr() + ylab(label = "Number of nodule")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
        axis.title.x = element_blank(),
        title = element_text(hjust = 0.5, size = 12),
        plot.title = element_text(hjust = 0.5, size = 12), 
        legend.position = "none")+
  ylim(0,18) +
  scale_color_manual(values = c("#A3A500","#00B0F6","#00BF7D","#F8766D","white"))+
  geom_jitter(position=position_jitter(0.2))+
  geom_text(data = labels.df,size=5, aes(x = inoc_nut, y = max(V1), label = labels))

nod_total

pdf(paste(results.dir,"Figure_S4b_Nodules.pdf", sep=""), width=8, height=6)
print(nod_total)
dev.off()

#Figure S4c - nodulators - Taxonomic profile - nodules
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
cond="NOD"

samples_df_sub <- subset(samples_df, samples_df$Compartment != "RZ")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "AM")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "ES")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Condition != "NP")

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
samples_df_sub <- subset(samples_df, samples_df$Compartment != "RZ")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "AM")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "ES")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Condition != "NP")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Condition != "Input")

samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")
samples <- sample_data(samples_df_sub_2)
phylo_sub = phyloseq(OTU,TAX, samples)

# Transform to relative abundance
phylo_sub_RA <- microbiome::transform(phylo_sub, "compositional")

# Select the top 5 most abundant OTUs
top5 <- names(sort(taxa_sums(phylo_sub_RA), decreasing=TRUE))[1:4]

# Melt phyloseq object into a dataframe
data_melt <- psmelt(phylo_sub_RA)

# Create a new OTU column where non-top5 are grouped as "Other strains"
data_melt$OTU_grouped <- ifelse(
  data_melt$OTU %in% c("LjNodule214", "P1_H10", "P2_A12", "P2_D6"), 
  data_melt$OTU, 
  ifelse(!is.na(data_melt$genus) & data_melt$genus == "Mesorhizobium", "Other Mesorhizobium", "Other strains"))

# Manually set factor levels for ordered legend
otu_order <- c("LjNodule214", "P1_H10", "P2_A12", "P2_D6","Other Mesorhizobium", "Other strains")
data_melt$OTU_grouped <- factor(data_melt$OTU_grouped, levels = otu_order)

# Custom color palette (forcing "Other" to dark gray)
otu_colors <- c(
  "LjNodule214" = "#66E1D0",  # Blue
  "P1_H10" = "#00C1C8",       # Orange
  "P2_A12" = "#00AA95",       # Green
  "P2_D6" = "#00C18C",        # Red
  "Other Mesorhizobium" = "#00773E",  # Assign a unique color for Other Mesorhizobium
  "Other strains" = "darkgray"        # Forced dark gray
)

# Improved stacked barplot
bar <- ggplot(data_melt, aes(fill=OTU_grouped, y=Abundance, x=Biorep)) + 
  geom_bar(position="stack", stat="identity", colour = "darkgray", linewidth = 0.01) +  
  scale_fill_manual(values = otu_colors) +  # Apply custom colors
  ggtitle("Lotus nodule colonization profiles") + 
  theme_classic() +
  labs(y = "Relative abundance", fill = "OTU") +  
  labs(fill = "Isolate") +
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
bar2 <- bar + facet_wrap(~Inoculum+Nutrient, scales = "free_x", nrow = 1) +
  theme(strip.text.x = element_text(size = 12))

bar2

pdf(paste(results.dir,"Figure_S4c_nodule_taxonomic_profile.pdf", sep=""), width=7, height=4.5)
print(bar2)
dev.off()
