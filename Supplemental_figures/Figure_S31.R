library("ggplot2") #Version 3.4.2
library("dplyr") #Version 1.1.2

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Figure S31 - Niche replacement score ======
#With nodulators - plant-specific pathways
pathway_selection <- c("ABC transporters", "Quorum sensing", "Two-component system", "Purine metabolism", "Oxidative phosphorylation", "Pentose and glucuronate interconversions", "Porphyrin metabolism", "Galactose metabolism", "Cell cycle - Caulobacter", "Exopolysaccharide biosynthesis")

Families <- c("Burkholderiaceae", "Caulobacteraceae", "Pseudomonadaceae", "other", "all")
#Make figure
data_fam <- read.table(paste(working_directory, "LjSC_Family_drop_out_experiment/Figure_5d_validation_in_vivo_Fam_drop_out_with_nod.txt", sep = ""), header =T, sep = "\t")

data_fam$Family[data_fam$Family == "LDT1"] <- "Burkholderiaceae"
data_fam$Family[data_fam$Family == "LDT2"] <- "Caulobacteraceae"
data_fam$Family[data_fam$Family == "LDT3"] <- "Pseudomonadaceae"
data_fam$Family[data_fam$Family == "LDT4"] <- "Rhizobiaceae"
data_fam$Family[data_fam$Family == "LDT5"] <- "other"
data_fam$Family[data_fam$Family == "LDT6"] <- "No_drop_out"
data_fam$Family[data_fam$Family == "all"] <- "No_drop_out"

data_fam$Log <- log2(data_fam$Present_fold_change)

data_fam_3 <- data_fam[data_fam$Pathway %in% pathway_selection,]
data_fam_3$Family <- factor(data_fam_3$Family, levels = c("No_drop_out", "Burkholderiaceae", "Caulobacteraceae", "Pseudomonadaceae", "Rhizobiaceae", "other"))
groups_pathways <- unique(data_fam_3$Pathway)
groups_pathways_2 <- groups_pathways[order(groups_pathways,decreasing = T)]
data_fam_3$Pathway <- factor(data_fam_3$Pathway, levels = groups_pathways_2)

#Comparison - Rhizobiaceae removed as there was contamination in the experiment
Families <- c("Burkholderiaceae", "Caulobacteraceae", "Pseudomonadaceae", "other")
Pathway <- unique(data_fam_3$Pathway)

new_data <- data.frame()

for (fam in Families){
  for (path in Pathway){
    data_fam_sub <- data_fam_3[data_fam_3$Pathway == paste(path) & data_fam_3$Family == paste(fam),]
    data_sub <- data_fam_3[data_fam_3$Pathway == paste(path) & data_fam_3$Family == "No_drop_out",]
    
    value <- data_fam_sub$Present_fold_change/data_sub$Present_fold_change
    
    new_data_2 <- data.frame(t(data.frame(c(paste(path), paste(fam), value))))
    new_data <- rbind(new_data, new_data_2)
  }
}


row.names(new_data) <- NULL
colnames(new_data) <- c("Pathway", "Family", "Niche_replacement_score")

new_data$Family <- factor(new_data$Family, levels = c("Burkholderiaceae", "Caulobacteraceae", "Pseudomonadaceae", "other"))
new_data$Pathway <- factor(new_data$Pathway, levels = groups_pathways_2)
new_data$Log <- log2(as.numeric(new_data$Niche_replacement_score))

plot_2 <- new_data %>% ggplot() + geom_tile(aes(y=Pathway, x=Family, fill = as.numeric(Log))) +
  geom_text(aes(y=Pathway, x=Family, label = round(as.numeric(Log),2))) + 
  scale_fill_gradient2(low = "red", mid = "white", high = "green", midpoint = 0) +
  labs(y="Pathway", x="Dropped-out family", fill="", title="Niche replacement score") +
  theme_bw(base_size = 14) %+replace% theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  theme(plot.title = element_text(hjust =0.5))
plot_2

pdf(paste(results.dir,"Figure_S31a_Niche_replacement_with_dom.pdf", sep=""), width=7.5, height=6)
print(plot_2)
dev.off()

#Without nodulators - core pathways
pathway_selection <- c("Vitamin B6 metabolism","Transcriptional regulator","Secretion","Quorum sensing",
                       "Phenylalanine, tyrosine and tryptophan biosynthesis","Pantothenate and CoA biosynthesis",
                       "Oxidative phosphorylation","Membrane protein","Glycerophospholipid metabolism",
                       "Folate biosynthesis","Flagellar assembly","Exopolysaccharide biosynthesis","Cysteine and methionine metabolism",
                       "Arginine biosynthesis", "Methane metabolism", "Aminobenzoate degradation", "Ascorbate and aldarate metabolism")

Families <- c("Burkholderiaceae", "Caulobacteraceae", "Pseudomonadaceae", "other", "all")

data_fam <- read.table(paste(working_directory, "LjSC_Family_drop_out_experiment/Figure_5f_validation_in_vivo_Fam_drop_out_no_nod.txt", sep = ""), header =T, sep = "\t")

data_fam$Family[data_fam$Family == "LDT1"] <- "Burkholderiaceae"
data_fam$Family[data_fam$Family == "LDT2"] <- "Caulobacteraceae"
data_fam$Family[data_fam$Family == "LDT3"] <- "Pseudomonadaceae"
data_fam$Family[data_fam$Family == "LDT4"] <- "Rhizobiaceae"
data_fam$Family[data_fam$Family == "LDT5"] <- "other"
data_fam$Family[data_fam$Family == "LDT6"] <- "No_drop_out"

colnames(data_fam) <- c("RA", "Present_fold_change", "Absent_fold_change", "Family", "Pathway", "Percentage")

data_fam$Log <- log2(data_fam$Present_fold_change)

data_fam_3 <- data_fam[data_fam$Pathway %in% pathway_selection,]
data_fam_3$Family <- factor(data_fam_3$Family, levels = c("No_drop_out", "Burkholderiaceae", "Caulobacteraceae", "Pseudomonadaceae", "other"))
groups_pathways <- unique(data_fam_3$Pathway)
groups_pathways_2 <- groups_pathways[order(groups_pathways,decreasing = T)]
data_fam_3$Pathway <- factor(data_fam_3$Pathway, levels = groups_pathways_2)

#Comparison - Rhizobiaceae not taken along as there was contamination in those samples
Families <- c("Burkholderiaceae", "Caulobacteraceae", "Pseudomonadaceae", "other")
Pathway <- unique(data_fam_3$Pathway)

new_data <- data.frame()

for (fam in Families){
  for (path in Pathway){
    data_fam_sub <- na.omit(data_fam_3[data_fam_3$Pathway == paste(path) & data_fam_3$Family == paste(fam),])
    data_sub <- na.omit(data_fam_3[data_fam_3$Pathway == paste(path) & data_fam_3$Family == "No_drop_out",])
    
    value <- data_fam_sub$Present_fold_change/data_sub$Present_fold_change
    
    new_data_2 <- data.frame(t(data.frame(c(paste(path), paste(fam), value))))
    new_data <- rbind(new_data, new_data_2)
  }
}

row.names(new_data) <- NULL
colnames(new_data) <- c("Pathway", "Family", "Niche_replacement_score")

new_data$Family <- factor(new_data$Family, levels = c("Burkholderiaceae", "Caulobacteraceae", "Pseudomonadaceae", "other"))
new_data$Pathway <- factor(new_data$Pathway, levels = groups_pathways_2)
new_data$Log <- log2(as.numeric(new_data$Niche_replacement_score))

plot_2 <- new_data %>% ggplot() + geom_tile(aes(y=Pathway, x=Family, fill = as.numeric(Log))) +
  geom_text(aes(y=Pathway, x=Family, label = round(as.numeric(Log),2))) + 
  scale_fill_gradient2(low = "red", mid = "white", high = "green", midpoint = 0) +
  labs(y="Pathway", x="Dropped-out family", fill="", title="Niche replacement score") +
  theme_bw(base_size = 14) %+replace% theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  theme(plot.title = element_text(hjust =0.5))
plot_2

pdf(paste(results.dir,"Figure_S31b_Niche_replacement_no_dom.pdf", sep=""), width=9, height=8)
print(plot_2)
dev.off()
