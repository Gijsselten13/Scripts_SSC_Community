library("dplyr") #Version 1.1.2
library("phyloseq") #Version 1.44.0
library("vegan") #Version 2.6-4
library("ggplot2") #Version 3.4.2
library("ggpubr") #Version 0.6.0

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Figure S11 & 2i - Beta Diversity - distance to centroid =====
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
samples_df_2$Plant_inoc_compartment_nutrient=paste(samples_df$Condition, samples_df$Inoculum, samples_df$Compartment, samples_df$Nutrient, sep ="_")
samples_df_2$Plant_inoc_compartment=paste(samples_df$Condition, samples_df$Inoculum, samples_df$Compartment, sep ="_")

sapply(tax_df, function(x) length(unique(x)))
#  Class  number is not very different from Phylum number, skip class rank that gives weird results because nclass=6 and nphylum=5 

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

mean_value_genus= mean(beta_genus)

Bray_curtis_df=beta_genus

#Make PCoA plot for Bray Curtis Distance matrix
pcoa_tax = cmdscale(Bray_curtis_df, k=3, eig=T)
points_beta = as.data.frame(pcoa_tax$points)
colnames(points_beta) = c("x", "y", "z") 
eig = pcoa_tax$eig

points_beta = merge(points_beta,samples_df_sub_2, by = "row.names")
rownames(points_beta) <- points_beta$Row.names
points_beta <- points_beta %>% dplyr::select (-Row.names)

points_beta$Condition <- factor(points_beta$Condition, levels = c("At","Hv", "Lj"))
points_beta$Inoculum <- factor(points_beta$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC","NS"))
points_beta$Nutrient <- factor(points_beta$Nutrient, levels = c("low", "high"))
points_beta$Experiment  <- factor(points_beta$Experiment, levels = c("R1", "R2"))

# KO_otu table
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

Bray_curtis_KO_df=beta_isolate_KO

#Make PCoA plot for Bray Curtis Distance matrix
pcoa_KO = cmdscale(Bray_curtis_KO_df, k=3, eig=T)
points_beta_KO = as.data.frame(pcoa_KO$points)
colnames(points_beta_KO) = c("x", "y", "z") 
eig_KO = pcoa_KO$eig

points_beta_KO = merge(points_beta_KO,samples_df_sub_2, by = "row.names")
rownames(points_beta_KO) <- points_beta_KO$Row.names
points_beta_KO <- points_beta_KO %>% dplyr::select (-Row.names)

points_beta_KO$Condition <- factor(points_beta_KO$Condition, levels = c("At","Hv", "Lj"))
points_beta_KO$Inoculum <- factor(points_beta_KO$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC","NS"))
points_beta_KO$Nutrient <- factor(points_beta_KO$Nutrient, levels = c("low", "high"))
points_beta_KO$Experiment  <- factor(points_beta_KO$Experiment, levels = c("R1", "R2"))

# Subset tables if needed, by experiment
points_beta_tax_sub=subset.data.frame(x = points_beta, subset = points_beta$Experiment=="R2")
points_beta_KO_sub=subset.data.frame(x = points_beta_KO, subset = points_beta$Experiment=="R2")

# Define parameter for centroid calculation
param="Plant_inoc_compartment_nutrient"

# Calculate centroids for each group
centroids_tax <- points_beta_tax_sub %>%
  group_by(!!sym(param)) %>%
  dplyr::summarize(
    centroid_x = mean(x, na.rm = TRUE),
    centroid_y = mean(y, na.rm = TRUE),
    centroid_z = mean(z, na.rm = TRUE)
  )

# Join centroids back to the original data
data_with_centroids_tax <- left_join(points_beta_tax_sub, centroids_tax , by = param)

# Calculate distance to centroid for each point
data_with_centroids_tax  <- data_with_centroids_tax  %>% 
  rowwise() %>%
  mutate(distance_to_centroid = sqrt((x - centroid_x)^2 + 
                                       (y - centroid_y)^2 + 
                                       (z - centroid_z)^2))

custom_labels_plant <- c(At = "A.thaliana", Hv = "H. vulgare", Lj = "L. japonicus")

# Calculate centroids for each group
centroids_KO <- points_beta_KO_sub %>%
  group_by(!!sym(param)) %>%
  dplyr::summarize(
    centroid_x = mean(x, na.rm = TRUE),
    centroid_y = mean(y, na.rm = TRUE),
    centroid_z = mean(z, na.rm = TRUE)
  )

# Join centroids back to the original data
data_with_centroids_KO <- left_join(points_beta_KO_sub, centroids_KO , by = param)

# Calculate distance to centroid for each point
data_with_centroids_KO  <- data_with_centroids_KO  %>% 
  rowwise() %>%
  mutate(distance_to_centroid = sqrt((x - centroid_x)^2 + 
                                       (y - centroid_y)^2 + 
                                       (z - centroid_z)^2))

custom_labels_plant <- c(At = "A.thaliana", Hv = "H. vulgare", Lj = "L. japonicus")

#  Combine both datasets to compare functional vs taxonomic distance from centroids
data_with_centroids_KO$Type  <- "Functional" 
data_with_centroids_tax$Type  <- "Taxonomic" 

centroid_merge=rbind(data_with_centroids_KO,data_with_centroids_tax)

dist_cen <- ggplot(centroid_merge, aes(x = param, y = distance_to_centroid, color=Type)) +
  theme_classic()+
  scale_color_manual(values = c("gray70","black"))+
  scale_shape_manual(values = c(0,3))+
  scale_fill_manual(values = c("white","gray70"))+
  geom_boxplot(outlier.shape = NA) + # Hide outliers since jitter will show all points
  facet_wrap(as.formula(paste(".~", param)), scales = "free_x", ncol = nlevels(factor(centroid_merge[[param]]))/3) +
  geom_jitter()+
  theme( axis.text.x=element_blank(), 
         axis.title.x=element_blank(), 
         title=element_text(hjust=0.5, size=15), 
         axis.ticks.x=element_blank(),
         strip.background=element_rect(colour="gray50", size=0.3), # Change 'size' for thickness
         axis.text=element_text(color="gray50"),
         axis.line = element_line(color="gray50", size=0.3)) +
  labs(title = "Distance to Centroid Functional vs Taxonomic",
       x = param,
       y = "Distance to Centroid") +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), label = "p.signif", vjust = 0.7)
dist_cen

pdf(paste(results.dir,"Figure_S11_distance_to_centroid.pdf", sep=""), width=13, height=6)
print(dist_cen)
dev.off()

#Figure 2i - Distance to centroid stats 
centroid_merge$Type <- factor(centroid_merge$Type, levels = c("Functional","Taxonomic"))

centroid_merge$Type <- factor(centroid_merge$Type, levels = c("Taxonomic","Functional"))

centro_stats_hor <- ggpaired(centroid_merge, x = "Type", y = "distance_to_centroid",
                             color = "Type", line.color = "transparent", line.size = 0.1, linetype = "solid")+
  labs(x = "", y = "Distance to centroid") +
  scale_color_manual(values = c("black","gray70"))+
  stat_compare_means(paired = TRUE, size =5, hjust=-0.5) +
  theme(legend.position = "none") +
  theme(axis.title = element_text(size = 15, hjust = 0.5), axis.text.x = element_text(size = 15, angle = 45, vjust=1, hjust=1), title = element_text(size =15)) +
  ggtitle("")
centro_stats_hor

pdf(paste(results.dir,"Figure_2i_stats_distance_to_centroid_hor.pdf", sep=""), width=3, height=4)
print(centro_stats_hor)
dev.off()
