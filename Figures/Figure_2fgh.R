library("dplyr") #Version 1.1.2
library("phyloseq") #Version 1.44.0
library("vegan") #Version 2.6-4
library("ggplot2") #Version 3.4.2
library("ggpubr") #Version 0.6.0
library("grid") #Version 4.4.1

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = ""

###Figure 2fgh - Beta Diversity plots =====
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

samples_df_2$Condition <- gsub("At", "Arabidopsis", samples_df_2$Condition)
samples_df_2$Condition <- gsub("Hv", "Barley", samples_df_2$Condition)
samples_df_2$Condition <- gsub("Lj", "Lotus", samples_df_2$Condition)
samples_df_2$Nutrient <- gsub("high", "High", samples_df_2$Nutrient)
samples_df_2$Nutrient <- gsub("low", "Low", samples_df_2$Nutrient)

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

Bray_curtis_df=beta_genus

#Make PCoA plot for Bray Curtis Distance matrix
pcoa_tax = cmdscale(Bray_curtis_df, k=3, eig=T)
points_beta = as.data.frame(pcoa_tax$points)
colnames(points_beta) = c("x", "y", "z") 
eig = pcoa_tax$eig

points_beta = merge(points_beta,samples_df_sub_2, by = "row.names")
rownames(points_beta) <- points_beta$Row.names
points_beta <- points_beta %>% dplyr::select (-Row.names)

points_beta$Condition <- factor(points_beta$Condition, levels = c("Arabidopsis","Barley", "Lotus"))
points_beta$Inoculum <- factor(points_beta$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC","NS"))
points_beta$Nutrient <- factor(points_beta$Nutrient, levels = c("Low", "High"))
points_beta$Experiment  <- factor(points_beta$Experiment, levels = c("R1", "R2"))

#  Run adonis PERMANOVA test
metadata=points_beta[,-c(1,2,3)]
set.seed(1)
SSC_bray_adonis <- adonis2(beta_distance_list[[2]] ~ Inoculum*Condition*Nutrient*Experiment, data=metadata, method="bray", permutations=999)

#Bray Curtis - SynCom
if (SSC_bray_adonis$`Pr(>F)`[1] <= 0.001) {
  sig_syncom <- "***"
} else if (SSC_bray_adonis$`Pr(>F)`[1] <= 0.01) {
  sig_syncom <- "**"
} else if (SSC_bray_adonis$`Pr(>F)`[1] <= 0.05) {
  sig_syncom <- "*"
} else {
  sig_syncom <- "ns"
}

#Bray Curtis - Plant
if (SSC_bray_adonis$`Pr(>F)`[2] <= 0.001) {
  sig_plant <- "***"
} else if (SSC_bray_adonis$`Pr(>F)`[2] <= 0.01) {
  sig_plant <- "**"
} else if (SSC_bray_adonis$`Pr(>F)`[2] <= 0.05) {
  sig_plant <- "*"
} else {
  sig_plant <- "ns"
}

#Bray Curtis - Nutrient
if (SSC_bray_adonis$`Pr(>F)`[3] <= 0.001) {
  sig_nut <- "***"
} else if (SSC_bray_adonis$`Pr(>F)`[3] <= 0.01) {
  sig_nut <- "**"
} else if (SSC_bray_adonis$`Pr(>F)`[3] <= 0.05) {
  sig_nut <- "*"
} else {
  sig_nut <- "ns"
}

#Bray Curtis - Experiment
if (SSC_bray_adonis$`Pr(>F)`[4] <= 0.001) {
  sig_exp <- "***"
} else if (SSC_bray_adonis$`Pr(>F)`[4] <= 0.01) {
  sig_exp <- "**"
} else if (SSC_bray_adonis$`Pr(>F)`[4] <= 0.05) {
  sig_exp <- "*"
} else {
  sig_exp <- "ns"
}

#Bray Curtis - SynCom:Plant Interaction
if (SSC_bray_adonis$`Pr(>F)`[5] <= 0.001) {
  sig_inter <- "***"
} else if (SSC_bray_adonis$`Pr(>F)`[5] <= 0.01) {
  sig_inter <- "**"
} else if (SSC_bray_adonis$`Pr(>F)`[5] <= 0.05) {
  sig_inter <- "*"
} else {
  sig_inter <- "ns"
}

colnames(points_beta)[colnames(points_beta) == "Condition"] <- "Host"

pbray_tax_1 <- ggplot(points_beta, aes(x=x, y=y, fill=Inoculum, color=Nutrient, shape=Host, stroke=0.3))+
  scale_shape_manual(values = c(23,24,21))+
  scale_color_manual(values = c("gray70","black"))+
  scale_fill_manual(values = c("#A3A500","#00B0F6","#00BF7D","#F8766D","White"))+
  theme(panel.background=element_blank(),panel.grid=element_blank(),axis.line.x=element_line(size=.5, color="black"),axis.line.y=element_line(size=.5, color="black"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black", size=7),legend.position="right",legend.background=element_blank(),legend.key=element_blank(),legend.text= element_text(size=10),text=element_text(family="sans", size=10))+
  geom_point(alpha=1, size=4)+
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""))+
  theme(axis.title = element_text(size = 14, hjust = 1), axis.text = element_text(size = 12, hjust = 1))+
  guides(fill = guide_legend(override.aes=list(shape=21, size=4)))+
  guides(color = guide_legend(override.aes=list(shape=21, size = 4))) +
  guides(shape = guide_legend(override.aes=list(size = 4))) +
  stat_ellipse(aes(group = Inoculum ), size=0.1)+
  ggtitle(paste0("Taxonomic (",Taxo_order[2], ") composition"))
pbray_tax_1

text <- c(paste0("PCoA Bray Curtis - ",Taxo_order[2], " composition"),
          paste("Inoculum - R2 = ",round(SSC_bray_adonis$R2[1],3),sig_syncom),
          paste("Host - R2 = ",round(SSC_bray_adonis$R2[2],3),sig_plant),
          paste("Nutrient - R2 = ",round(SSC_bray_adonis$R2[3],3),sig_nut),
          paste("Experiment - R2 = ",round(SSC_bray_adonis$R2[4],3),sig_exp),
          paste("Inoculum:Host - R2 = ",round(SSC_bray_adonis$R2[5],3),sig_inter))

tl <- textGrob(paste(strwrap(text, 40), collapse="\n"), hjust=0, x=0)
pbray_tax_1_2 <- ggarrange(pbray_tax_1,tl, ncol = 2, nrow =1, widths = c(2,1))

# KO_otu table
KO_SSC =read.table(paste(working_directory,"KO_tables/Original/SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Taxonomy table

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object
OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)

#Sample subsetting
cond="ES"

samples_df_sub <- subset(samples_df_2, samples_df_2$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")
samples_df_sub_2$Condition[samples_df_sub_2$Condition == "At"] <- "Arabidopsis"
samples_df_sub_2$Condition[samples_df_sub_2$Condition == "Hv"] <- "Barley"
samples_df_sub_2$Condition[samples_df_sub_2$Condition == "Lj"] <- "Lotus"
samples_df_sub_2$Nutrient[samples_df_sub_2$Nutrient == "high"] <- "High"
samples_df_sub_2$Nutrient[samples_df_sub_2$Nutrient == "low"] <- "Low"
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

points_beta_KO$Condition <- factor(points_beta_KO$Condition, levels = c("Arabidopsis","Barley", "Lotus"))
points_beta_KO$Inoculum <- factor(points_beta_KO$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC","NS"))
points_beta_KO$Nutrient <- factor(points_beta_KO$Nutrient, levels = c("Low", "High"))
points_beta_KO$Experiment  <- factor(points_beta_KO$Experiment, levels = c("R1", "R2"))

# PERMANOVA Adonis test
set.seed(1)
metadata=points_beta_KO[,-c(1,2,3)]
SSC_bray_KO_adonis <- adonis2(beta_isolate_KO ~ Inoculum*Condition*Nutrient*Experiment, data=metadata, method="bray", permutations=999)

#Bray Curtis - SynCom
if (SSC_bray_KO_adonis$`Pr(>F)`[1] <= 0.001) {
  sig_syncom <- "***"
} else if (SSC_bray_KO_adonis$`Pr(>F)`[1] <= 0.01) {
  sig_syncom <- "**"
} else if (SSC_bray_KO_adonis$`Pr(>F)`[1] <= 0.05) {
  sig_syncom <- "*"
} else {
  sig_syncom <- "ns"
}

#Bray Curtis - Plant
if (SSC_bray_KO_adonis$`Pr(>F)`[2] <= 0.001) {
  sig_plant <- "***"
} else if (SSC_bray_KO_adonis$`Pr(>F)`[2] <= 0.01) {
  sig_plant <- "**"
} else if (SSC_bray_KO_adonis$`Pr(>F)`[2] <= 0.05) {
  sig_plant <- "*"
} else {
  sig_plant <- "ns"
}

#Bray Curtis - Nutrient
if (SSC_bray_KO_adonis$`Pr(>F)`[3] <= 0.001) {
  sig_nut <- "***"
} else if (SSC_bray_KO_adonis$`Pr(>F)`[3] <= 0.01) {
  sig_nut <- "**"
} else if (SSC_bray_KO_adonis$`Pr(>F)`[3] <= 0.05) {
  sig_nut <- "*"
} else {
  sig_nut <- "ns"
}

#Bray Curtis - Experiment
if (SSC_bray_KO_adonis$`Pr(>F)`[4] <= 0.001) {
  sig_exp <- "***"
} else if (SSC_bray_KO_adonis$`Pr(>F)`[4] <= 0.01) {
  sig_exp <- "**"
} else if (SSC_bray_KO_adonis$`Pr(>F)`[4] <= 0.05) {
  sig_exp <- "*"
} else {
  sig_exp <- "ns"
}

#Bray Curtis - SynCom:Plant Interaction
if (SSC_bray_KO_adonis$`Pr(>F)`[5] <= 0.001) {
  sig_inter <- "***"
} else if (SSC_bray_KO_adonis$`Pr(>F)`[5] <= 0.01) {
  sig_inter <- "**"
} else if (SSC_bray_KO_adonis$`Pr(>F)`[5] <= 0.05) {
  sig_inter <- "*"
} else {
  sig_inter <- "ns"
}

colnames(points_beta_KO)[colnames(points_beta_KO) == "Condition"] <- "Host"

pbray_KO_1 <- ggplot(points_beta_KO, aes(x=x, y=y, fill=Inoculum, color=Nutrient, shape=Host,size=12,stroke=0.3))+
  scale_shape_manual(values = c(23,24,21))+
  scale_color_manual(values = c("gray70","black"))+
  scale_fill_manual(values = c("#A3A500","#00B0F6","#00BF7D","#F8766D","White"))+
  theme(panel.background=element_blank(),panel.grid=element_blank(),axis.line.x=element_line(size=.5, color="black"),axis.line.y=element_line(size=.5, color="black"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black", size=7),legend.position="none",legend.background=element_blank(),legend.key=element_blank(),legend.text= element_text(size=10),text=element_text(family="sans", size=10))+
  geom_point(alpha=1, size=4,)+
  labs(x=paste("PCoA 1 (", format(100 * eig_KO[1] / sum(eig_KO), digits=4), "%)", sep=""),y=paste("PCoA 2 (", format(100 * eig_KO[2] / sum(eig_KO), digits=4), "%)", sep=""))+
  theme(axis.title = element_text(size = 14, hjust = 1), axis.text = element_text(size = 12, hjust = 1))+
  guides(fill = guide_legend(override.aes=list(shape=21, size=4)))+
  guides(color = guide_legend(override.aes=list(shape=21, size = 4))) +
  guides(shape = guide_legend(override.aes=list(size = 4))) +
  stat_ellipse(aes(group = Inoculum ), size=0.1)+
  ggtitle("Functional (KO) composition")
pbray_KO_1

#Plot text next to it
text_KO <- c(paste0("PCoA Bray Curtis - KO composition"), 
             paste("Inoculum - R2 = ",round(SSC_bray_KO_adonis$R2[1],3),sig_syncom),
             paste("Host - R2 = ",round(SSC_bray_KO_adonis$R2[2],3),sig_plant),
             paste("Nutrient - R2 = ",round(SSC_bray_KO_adonis$R2[3],3),sig_nut),
             paste("Experiment - R2 = ",round(SSC_bray_KO_adonis$R2[4],3),sig_exp),
             paste("Inoculum:Host - R2 = ",round(SSC_bray_KO_adonis$R2[5],3),sig_inter))

tl <- textGrob(paste(strwrap(text_KO, 40), collapse="\n"), hjust=0, x=0)
plot_KO <- ggarrange(pbray_KO_1,tl, ncol = 2, nrow =1, widths = c(2,1))

genus_plot <- pbray_tax_1
genus_plot<- genus_plot+theme(legend.position = "none")

# Define legend text size and key size as variables
text_size <- 14
key_size <- 1.2

# Create a theme object with your legend modifications
legend_theme <- theme(legend.text = element_text(size = text_size),  # Apply text size
                      legend.key.size = unit(key_size, "cm"), # Apply key size
                      strip.background=element_rect(colour="gray50", size=0.3), # Change 'size' for thickness
                      axis.text=element_text(color="gray50"),
                      axis.line = element_line(color="gray50", size=0.3),
                      axis.line.x = element_line(color="gray50", size=0.3), 
                      axis.line.y = element_line(color="gray50", size=0.3),
                      axis.ticks.x = element_line(color="gray50", size=0.3),
                      axis.ticks.y =element_line(color="gray50", size=0.3)
)

# Apply the legend_theme to each ggplot object
genus_plot <- genus_plot + legend_theme
pbray_KO_1 <- pbray_KO_1 + legend_theme

#Adonis test on taxonomy tables
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
  
  points$Condition <- factor(points$Condition, levels = c("Arabidopsis","Barley", "Lotus"))
  points$Inoculum <- factor(points$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC","NS"))
  points$Nutrient <- factor(points$Nutrient, levels = c("Low", "High"))
  points$Experiment  <- factor(points$Experiment, levels = c("R1", "R2"))
  
  metadata=points[,-c(1,2,3)]
  
  #  Run adonis PERMANOVA test
  set.seed(1)
  SSC_bray_adonis <- adonis2(beta_distance_list[[i]] ~ Inoculum*Condition*Nutrient*Experiment, data=metadata, method="bray", permutations=999)
  
  Adonis_list[[i]]=SSC_bray_adonis
  
}

# Add values of the adonis tests in dataframes to make a table:
Taxo_R2=as.data.frame(Adonis_list[[2]])
Func_R2=as.data.frame(SSC_bray_KO_adonis)

R2_table=merge.data.frame(Taxo_R2,Func_R2, by = 0)
R2_table=R2_table[,-c(2,3,5,6,7,8,10,11)]

# Subset specific rows
R2_table <- R2_table[c(1,4,5,6,11), ]  # Select rows 2 to 4 and rows 11 to 13
R2_table <- R2_table[order(-R2_table$R2.x), ]

# Round numeric values to 2 digits after the decimal point
R2_table$R2.x <- round(R2_table$R2.x, digits = 2)
R2_table$R2.y <- round(R2_table$R2.y, digits = 2)

# Rename columns
colnames(R2_table)=c("Taxonomy R²","Parameter","Function R²")

# Replace "Condition" with "Plant"
R2_table$Parameter <- gsub("Condition", "Host", R2_table$Parameter)

# my_table <- tableGrob(R2_table)
ggtexttable(R2_table, rows = NULL, theme = ttheme("classic"))

beta_figures <- ggarrange(genus_plot,"", pbray_KO_1+scale_y_continuous(position="right"), nrow =1, ncol =3, common.legend = F, widths = c(1,0.3,1))
beta_figures

pdf(paste(results.dir,"Figure_2fg_Beta_div_plots.pdf", sep=""), width=width_in, height=height_in)
print(beta_figures)
dev.off()

pdf(paste(results.dir,"Figure_2h_table_R2.pdf", sep=""), width=width_in, height=height_in)
ggtexttable(R2_table, rows = NULL, theme = ttheme("classic"))
dev.off()
