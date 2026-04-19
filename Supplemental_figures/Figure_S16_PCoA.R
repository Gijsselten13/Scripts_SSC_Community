library("phyloseq") #Version 1.44.0
library("vegan") #Version 2.6-4
library("ggplot2") #Version 3.4.2
library("ggpubr") #Version 0.6.0

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Figure S16 - PCoAs to compare sPLS-DA results =====
#Rhizobiaceae - with dom - HvSC, LjSC, SSC
#otu table
KO_SSC_only=read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Rhizobiaceae_KO_with_dom.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t", row.names =1) #make the SampleID column into the row.names
colnames(samples_df)[5]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object
#Sample subsetting
samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

#Subset for all SynComs but AtSC
samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum != "AtSC")

OTU_KO = otu_table(as.matrix(KO_SSC_only),taxa_are_rows = TRUE)
samples_sub = sample_data(samples_df_sub_3)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

bray_2 <- as.matrix(beta_isolate_KO)

str(samples_df_sub_3)
str(bray_2)

#Bind metadata with distance matrix
pcoa = cmdscale(bray_2, k=10, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z", "a", "b", "c", "d", "e", "f", "g") 
eig = pcoa$eig
points_2 <- points[order(row.names(points)), ]
samples_df_sub_6 <- samples_df_sub_3[row.names(samples_df_sub_3) %in% row.names(points),]
samples_df_sub_7 <- samples_df_sub_6[order(row.names(samples_df_sub_6)), ]
points_3 <- cbind(points_2,samples_df_sub_7)
colnames(points_3) <- c("PCoA Axis 1", "PCoA Axis 2", "PCoA Axis 3","PCoA Axis 4", "b", "c", "d", "e", "f", "g",colnames(samples_df_sub_7))

groups <- c("PCoA Axis 1", "PCoA Axis 2", "PCoA Axis 3","PCoA Axis 4")

for (group in groups) {
  for (group_2 in groups) {
    if (group != group_2) {  # Avoid comparing the same group with itself
      points_sub <- points_3[, colnames(points_3) %in% c(group, group_2, "Condition")]
      if (group == "PCoA Axis 1"){
        axis_x <- paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep="")
      } else if (group == "PCoA Axis 2"){
        axis_x <- paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")
      } else if (group == "PCoA Axis 3"){
        axis_x <- paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep="")
      } else {
        axis_x <- paste("PCoA 4 (", format(100 * eig[4] / sum(eig), digits=4), "%)", sep="")
      }
      
      if (group_2 == "PCoA Axis 1"){
        axis_y <- paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep="")
      } else if (group_2 == "PCoA Axis 2"){
        axis_y <- paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")
      } else if (group_2 == "PCoA Axis 3"){
        axis_y <- paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep="")
      } else {
        axis_y <- paste("PCoA 4 (", format(100 * eig[4] / sum(eig), digits=4), "%)", sep="")
      }
      
      # PCoA plot with vectors
      pcoa_plot <- ggplot(points_sub, aes_string(x = paste("`", group, "`", sep = ""), 
                                                 y = paste("`", group_2, "`", sep = ""), 
                                                 color = "Condition")) +
        geom_point(size = 3) +
        theme_classic() +
        scale_colour_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +  # Center title
        labs(title = "",
             x = axis_x, 
             y = axis_y)
      
      # Store the plot dynamically
      assign(paste(group, group_2, sep = "_"), pcoa_plot)
    }
  }
}


one <- ggarrange(`PCoA Axis 1_PCoA Axis 2`, `PCoA Axis 1_PCoA Axis 3`, `PCoA Axis 1_PCoA Axis 4`, nrow =3, ncol =1, common.legend = T)
two <- ggarrange(NULL, `PCoA Axis 2_PCoA Axis 3`, `PCoA Axis 2_PCoA Axis 4`, nrow =3, ncol =1,common.legend = T)
three <- ggarrange(NULL, NULL, `PCoA Axis 3_PCoA Axis 4`, nrow =3, ncol =1,common.legend = T)

all <- ggarrange(one, two, three, ncol = 3, nrow =1,common.legend = T)
all

pdf(paste(results.dir,"Figure_S16_Rhizobiaceae_PC4.pdf", sep=""), width=10, height=10)
print(all)
dev.off()

#Xanthomonadaceae - SynCom - no dom 

#otu table
KO_SSC_only=read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Xanthomonadaceae_KO_no_dom.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t", row.names =1) #make the SampleID column into the row.names
colnames(samples_df)[5]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object

#Sample subsetting
samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

OTU_KO = otu_table(as.matrix(KO_SSC_only),taxa_are_rows = TRUE)
samples_sub = sample_data(samples_df_sub_2)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

bray_2 <- as.matrix(beta_isolate_KO)

str(samples_df_sub_2)
str(bray_2)

#Bind metadata with distance matrix
pcoa = cmdscale(bray_2, k=10, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z", "a", "b", "c", "d", "e", "f", "g") 
eig = pcoa$eig
points_2 <- points[order(row.names(points)), ]
samples_df_sub_6 <- samples_df_sub_2[row.names(samples_df_sub_2) %in% row.names(points),]
samples_df_sub_7 <- samples_df_sub_6[order(row.names(samples_df_sub_6)), ]
points_3 <- cbind(points_2,samples_df_sub_7)
colnames(points_3) <- c("PCoA Axis 1", "PCoA Axis 2", "PCoA Axis 3","PCoA Axis 4", "b", "c", "d", "e", "f", "g",colnames(samples_df_sub_7))

groups <- c("PCoA Axis 1", "PCoA Axis 2", "PCoA Axis 3","PCoA Axis 4")

for (group in groups) {
  for (group_2 in groups) {
    if (group != group_2) {  # Avoid comparing the same group with itself
      points_sub <- points_3[, colnames(points_3) %in% c(group, group_2, "Inoculum")]
      
      if (group == "PCoA Axis 1"){
        axis_x <- paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep="")
      } else if (group == "PCoA Axis 2"){
        axis_x <- paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")
      } else if (group == "PCoA Axis 3"){
        axis_x <- paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep="")
      } else {
        axis_x <- paste("PCoA 4 (", format(100 * eig[4] / sum(eig), digits=4), "%)", sep="")
      }
      
      if (group_2 == "PCoA Axis 1"){
        axis_y <- paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep="")
      } else if (group_2 == "PCoA Axis 2"){
        axis_y <- paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")
      } else if (group_2 == "PCoA Axis 3"){
        axis_y <- paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep="")
      } else {
        axis_y <- paste("PCoA 4 (", format(100 * eig[4] / sum(eig), digits=4), "%)", sep="")
      }
      
      # PCoA plot with vectors
      pcoa_plot <- ggplot(points_sub, aes_string(x = paste("`", group, "`", sep = ""), 
                                                 y = paste("`", group_2, "`", sep = ""), 
                                                 color = "Inoculum")) +
        geom_point(size = 3) +
        theme_classic() +
        scale_colour_manual(values = c("#A3A500","#00B0F6","#00BF7D","#F8766D")) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +  # Center title
        labs(title = "",
             x = axis_x, 
             y = axis_y)
      
      # Store the plot dynamically
      assign(paste(group, group_2, sep = "_"), pcoa_plot)
    }
  }
}

one <- ggarrange(`PCoA Axis 1_PCoA Axis 2`, `PCoA Axis 1_PCoA Axis 3`, `PCoA Axis 1_PCoA Axis 4`, nrow =3, ncol =1, common.legend = T)
two <- ggarrange(NULL, `PCoA Axis 2_PCoA Axis 3`, `PCoA Axis 2_PCoA Axis 4`, nrow =3, ncol =1,common.legend = T)
three <- ggarrange(NULL, NULL, `PCoA Axis 3_PCoA Axis 4`, nrow =3, ncol =1,common.legend = T)

all <- ggarrange(one, two, three, ncol = 3, nrow =1,common.legend = T)
all

pdf(paste(results.dir,"Figure_S16_Xanthomonadaceae_PC4.pdf", sep=""), width=10, height=10)
print(all)
dev.off()

#Burkholderiaceae - SynCom - no dom
#otu table
KO_SSC_only=read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Burkholderiaceae_KO_no_dom.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t", row.names =1) #make the SampleID column into the row.names
colnames(samples_df)[5]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object

#Sample subsetting
samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

#At subset - Enterobacteriaceae - dom vs no-dom 
OTU_KO = otu_table(as.matrix(KO_SSC_only),taxa_are_rows = TRUE)
samples_sub = sample_data(samples_df_sub_2)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

bray_2 <- as.matrix(beta_isolate_KO)

str(samples_df_sub_2)
str(bray_2)

#Bind metadata with distance matrix
pcoa = cmdscale(bray_2, k=10, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z", "a", "b", "c", "d", "e", "f", "g") 
eig = pcoa$eig
points_2 <- points[order(row.names(points)), ]
samples_df_sub_6 <- samples_df_sub_2[row.names(samples_df_sub_2) %in% row.names(points),]
samples_df_sub_7 <- samples_df_sub_6[order(row.names(samples_df_sub_6)), ]
points_3 <- cbind(points_2,samples_df_sub_7)
colnames(points_3) <- c("PCoA Axis 1", "PCoA Axis 2", "PCoA Axis 3","a", "b", "c", "d", "e", "f", "g",colnames(samples_df_sub_7))

groups <- c("PCoA Axis 1", "PCoA Axis 2", "PCoA Axis 3")

for (group in groups) {
  for (group_2 in groups) {
    if (group != group_2) {  # Avoid comparing the same group with itself
      points_sub <- points_3[, colnames(points_3) %in% c(group, group_2, "Inoculum")]
      
      if (group == "PCoA Axis 1"){
        axis_x <- paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep="")
      } else if (group == "PCoA Axis 2"){
        axis_x <- paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")
      } else if (group == "PCoA Axis 3"){
        axis_x <- paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep="")
      } 
      
      if (group_2 == "PCoA Axis 1"){
        axis_y <- paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep="")
      } else if (group_2 == "PCoA Axis 2"){
        axis_y <- paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")
      } else if (group_2 == "PCoA Axis 3"){
        axis_y <- paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep="")
      } 
      
      # PCoA plot with vectors
      pcoa_plot <- ggplot(points_sub, aes_string(x = paste("`", group, "`", sep = ""), 
                                                 y = paste("`", group_2, "`", sep = ""), 
                                                 color = "Inoculum")) +
        geom_point(size = 3) +
        theme_classic() +
        scale_colour_manual(values = c("#A3A500","#00B0F6","#00BF7D","#F8766D")) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +  # Center title
        labs(title = "",
             x = axis_x, 
             y = axis_y)
      
      # Store the plot dynamically
      assign(paste(group, group_2, sep = "_"), pcoa_plot)
    }
  }
}

one <- ggarrange(`PCoA Axis 1_PCoA Axis 2`, `PCoA Axis 1_PCoA Axis 3`, nrow =2, ncol =1, common.legend = T)
two <- ggarrange(NULL, `PCoA Axis 2_PCoA Axis 3`, nrow =2, ncol =1,common.legend = T)

all <- ggarrange(one, two,ncol = 2, nrow =1,common.legend = T)
all

pdf(paste(results.dir,"Figure_S16_Burkholderiaceae_PC3.pdf", sep=""), width=8, height=8)
print(all)
dev.off()

#Xanthomondaceae - Plant - LjSC - no dom
#otu table
KO_SSC_only=read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Xanthomonadaceae_KO_no_dom.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t", row.names =1) #make the SampleID column into the row.names
colnames(samples_df)[5]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object

#Sample subsetting
samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

#Subset for LjSC
samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum == "LjSC")

OTU_KO = otu_table(as.matrix(KO_SSC_only),taxa_are_rows = TRUE)
samples_sub = sample_data(samples_df_sub_3)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

bray_2 <- as.matrix(beta_isolate_KO)

str(samples_df_sub_3)
str(bray_2)

#Bind metadata with distance matrix
pcoa = cmdscale(bray_2, k=10, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z", "a", "b", "c", "d", "e", "f", "g") 
eig = pcoa$eig
points_2 <- points[order(row.names(points)), ]
samples_df_sub_6 <- samples_df_sub_3[row.names(samples_df_sub_3) %in% row.names(points),]
samples_df_sub_7 <- samples_df_sub_6[order(row.names(samples_df_sub_6)), ]
points_3 <- cbind(points_2,samples_df_sub_7)
colnames(points_3) <- c("PCoA Axis 1", "PCoA Axis 2", "z","a", "b", "c", "d", "e", "f", "g",colnames(samples_df_sub_7))
groups <- c("PCoA Axis 1", "PCoA Axis 2")

# Avoid comparing the same group with itself
points_sub <- points_3[, colnames(points_3) %in% groups]
points_sub$Condition <- samples_df_sub_7$Condition[match(row.names(points_sub), row.names(samples_df_sub_7))]

axis_x <- paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep="")
axis_y <- paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")

pcoa_plot <- ggplot(points_sub, aes_string(x = "`PCoA Axis 1`", 
                                           y = "`PCoA Axis 2`", 
                                           color = "Condition")) +
  geom_point(size = 3) +
  theme_classic() +
  scale_colour_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +  # Center title
  labs(title = "",
       x = axis_x, 
       y = axis_y)
pcoa_plot


pdf(paste(results.dir,"Figure_S16_Xanthomonadaceae_PC2_plant.pdf", sep=""), width=6, height=5)
print(pcoa_plot)
dev.off()

#Pseudomonadaceae - Plant - LjSC - no dom 

#otu table
KO_SSC_only=read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Pseudomonadaceae_KO_no_dom.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t", row.names =1) #make the SampleID column into the row.names
colnames(samples_df)[5]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object

#Sample subsetting
samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

#Subset for LjSC
samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum == "LjSC")

OTU_KO = otu_table(as.matrix(KO_SSC_only),taxa_are_rows = TRUE)
samples_sub = sample_data(samples_df_sub_3)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

bray_2 <- as.matrix(beta_isolate_KO)

str(samples_df_sub_3)
str(bray_2)

#Bind metadata with distance matrix
pcoa = cmdscale(bray_2, k=10, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z", "a", "b", "c", "d", "e", "f", "g") 
eig = pcoa$eig
points_2 <- points[order(row.names(points)), ]
samples_df_sub_6 <- samples_df_sub_3[row.names(samples_df_sub_3) %in% row.names(points),]
samples_df_sub_7 <- samples_df_sub_6[order(row.names(samples_df_sub_6)), ]
points_3 <- cbind(points_2,samples_df_sub_7)
colnames(points_3) <- c("PCoA Axis 1", "PCoA Axis 2", "z","a", "b", "c", "d", "e", "f", "g",colnames(samples_df_sub_7))
groups <- c("PCoA Axis 1", "PCoA Axis 2")

# Avoid comparing the same group with itself
points_sub <- points_3[, colnames(points_3) %in% groups]
points_sub$Condition <- samples_df_sub_7$Condition[match(row.names(points_sub), row.names(samples_df_sub_7))]

axis_x <- paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep="")
axis_y <- paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")

pcoa_plot <- ggplot(points_sub, aes_string(x = "`PCoA Axis 1`", 
                                           y = "`PCoA Axis 2`", 
                                           color = "Condition")) +
  geom_point(size = 3) +
  theme_classic() +
  scale_colour_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +  # Center title
  labs(title = "",
       x = axis_x, 
       y = axis_y)
pcoa_plot

pdf(paste(results.dir,"Figure_S16_Pseudomonadaceae_PC2_plant.pdf", sep=""), width=6, height=5)
print(pcoa_plot)
dev.off()

#Caulobacteraceae - Plant - AtSC - no dom 

#otu table
KO_SSC_only=read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Caulobacteraceae.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t", row.names =1) #make the SampleID column into the row.names
colnames(samples_df)[5]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object

#Sample subsetting
samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

#Subset for AtSC
samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum == "AtSC")

OTU_KO = otu_table(as.matrix(KO_SSC_only),taxa_are_rows = TRUE)
samples_sub = sample_data(samples_df_sub_3)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

bray_2 <- as.matrix(beta_isolate_KO)

str(samples_df_sub_3)
str(bray_2)

#Bind metadata with distance matrix
pcoa = cmdscale(bray_2, k=10, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z", "a", "b", "c", "d", "e", "f", "g") 
eig = pcoa$eig
points_2 <- points[order(row.names(points)), ]
samples_df_sub_6 <- samples_df_sub_3[row.names(samples_df_sub_3) %in% row.names(points),]
samples_df_sub_7 <- samples_df_sub_6[order(row.names(samples_df_sub_6)), ]
points_3 <- cbind(points_2,samples_df_sub_7)
colnames(points_3) <- c("PCoA Axis 1", "PCoA Axis 2", "z","a", "b", "c", "d", "e", "f", "g",colnames(samples_df_sub_7))
groups <- c("PCoA Axis 1", "PCoA Axis 2")

# Avoid comparing the same group with itself
points_sub <- points_3[, colnames(points_3) %in% groups]
points_sub$Condition <- samples_df_sub_7$Condition[match(row.names(points_sub), row.names(samples_df_sub_7))]

axis_x <- paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep="")
axis_y <- paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")

pcoa_plot <- ggplot(points_sub, aes_string(x = "`PCoA Axis 1`", 
                                           y = "`PCoA Axis 2`", 
                                           color = "Condition")) +
  geom_point(size = 3) +
  theme_classic() +
  scale_colour_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +  # Center title
  labs(title = "",
       x = axis_x, 
       y = axis_y)
pcoa_plot

pdf(paste(results.dir,"Figure_S16_Caulobacteraceae_PC2_plant.pdf", sep=""), width=6, height=5)
print(pcoa_plot)
dev.off()

#Burkholderiaceae - Plant - AtSC - no dom
#otu table
KO_SSC_only=read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Burkholderiaceae.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t", row.names =1) #make the SampleID column into the row.names
colnames(samples_df)[5]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object

#Sample subsetting
samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

#subset for AtSC
samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum == "AtSC")

OTU_KO = otu_table(as.matrix(KO_SSC_only),taxa_are_rows = TRUE)
samples_sub = sample_data(samples_df_sub_3)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

bray_2 <- as.matrix(beta_isolate_KO)

str(samples_df_sub_3)
str(bray_2)

#Bind metadata with distance matrix
pcoa = cmdscale(bray_2, k=10, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z", "a", "b", "c", "d", "e", "f", "g") 
eig = pcoa$eig
points_2 <- points[order(row.names(points)), ]
samples_df_sub_6 <- samples_df_sub_3[row.names(samples_df_sub_3) %in% row.names(points),]
samples_df_sub_7 <- samples_df_sub_6[order(row.names(samples_df_sub_6)), ]
points_3 <- cbind(points_2,samples_df_sub_7)
colnames(points_3) <- c("PCoA Axis 1", "PCoA Axis 2", "z","a", "b", "c", "d", "e", "f", "g",colnames(samples_df_sub_7))
groups <- c("PCoA Axis 1", "PCoA Axis 2")

# Avoid comparing the same group with itself
points_sub <- points_3[, colnames(points_3) %in% groups]
points_sub$Condition <- samples_df_sub_7$Condition[match(row.names(points_sub), row.names(samples_df_sub_7))]

axis_x <- paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep="")
axis_y <- paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")

pcoa_plot <- ggplot(points_sub, aes_string(x = "`PCoA Axis 1`", 
                                           y = "`PCoA Axis 2`", 
                                           color = "Condition")) +
  geom_point(size = 3) +
  theme_classic() +
  scale_colour_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +  # Center title
  labs(title = "",
       x = axis_x, 
       y = axis_y)
pcoa_plot

pdf(paste(results.dir,"Figure_S16_Burkholderiaceae_PC2_plant.pdf", sep=""), width=6, height=5)
print(pcoa_plot)
dev.off()
