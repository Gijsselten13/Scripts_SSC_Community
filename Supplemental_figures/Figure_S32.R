library("multcompView") #Version 0.1-9
library("ggplot2") #Version 3.4.2
library("scales") #Version 1.2.1
library("ggpubr") #Version 0.6.0

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Figure S32 - Bacterial/plant ratio =====
#With nodulators
plant_reads <- read.table(paste(working_directory, "LjSC_Family_drop_out_experiment/plant_reads.tsv", sep = ""), header =T, sep = "\t")
plant_reads_2 <- plant_reads[grep("ROOT", plant_reads$Sample),]

norm_SSC =read.table(paste(working_directory,"LjSC_Family_drop_out_experiment/Data_with_synthetic_input.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
norm_SSC_3 <- data.frame(colSums(norm_SSC))
colnames(norm_SSC_3) <- "No_of_bacterial_reads"
norm_SSC_3$Sample <- row.names(norm_SSC_3)
norm_SSC_4 <- norm_SSC_3[grep("ROOT", row.names(norm_SSC_3)),]
norm_SSC_4$No_of_plant_reads <- plant_reads_2$Host_mapped_Infr[match(norm_SSC_4$Sample, plant_reads_2$Sample)]
norm_SSC_4$Bacterial_read_per_plant <- norm_SSC_4$No_of_bacterial_reads/norm_SSC_4$No_of_plant_reads

norm_SSC_4_2 <- norm_SSC_4[!grepl("LDT7", norm_SSC_4$Sample),]
norm_SSC_5 <- norm_SSC_4_2[!grepl("LDT4", norm_SSC_4_2$Sample),]

norm_SSC_5$SynCom <- sapply(strsplit(as.vector(row.names(norm_SSC_5)), "_"),`[`, 1)

#Rhizobiaceae was not taken along as it was found to have a contamination
norm_SSC_5$SynCom <- gsub("LDT1", "Burkholderiaceae drop out", norm_SSC_5$SynCom)
norm_SSC_5$SynCom <- gsub("LDT2", "Caulobacteraceae drop out", norm_SSC_5$SynCom)
norm_SSC_5$SynCom <- gsub("LDT3", "Pseudomonadaceae drop out", norm_SSC_5$SynCom)
#norm_SSC_5$SynCom <- gsub("LDT4", "Rhizobiaceae drop out", norm_SSC_5$SynCom)
norm_SSC_5$SynCom <- gsub("LDT5", "All other families drop out", norm_SSC_5$SynCom)
norm_SSC_5$SynCom <- gsub("LDT6", "Full LjSC", norm_SSC_5$SynCom)

# Perform ANOVA
fitAnova <- aov(Bacterial_read_per_plant ~ SynCom, data = norm_SSC_5)

# Perform Tukey's post-hoc test
Tukey <- TukeyHSD(fitAnova)

# Get letters for significance
letters_anova <- multcompLetters4(fitAnova, Tukey)$SynCom$Letters

# Combine results into a data frame for plotting
ltlbl_combined <- data.frame(
  SynCom = names(letters_anova),
  Letters = letters_anova
)

# Ensure the order of the factor levels is correct
ltlbl_combined$SynCom <- factor(ltlbl_combined$SynCom, levels = unique(norm_SSC_5$SynCom))

# Apply the manually defined order to the factor
ltlbl_combined <- ltlbl_combined[order(ltlbl_combined$SynCom), ]

# Merge the letters with the data frame
norm_SSC_5 <- merge(norm_SSC_5, ltlbl_combined, by = "SynCom")

#Add number of 266 KOs
KO_table = read.table(paste(working_directory, "KO_genome/KO_LjSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

taxonomy <- read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",row.names =1)
taxonomy_LjSC <- taxonomy[taxonomy$SynCom == "LjSC",]

families <- c("Pseudomonadaceae", "Rhizobiaceae", "Burkholderiaceae", "Caulobacteraceae", "other")

input_table <- read.table(paste(working_directory,"DESeq2/Sig_KO_all_no_nod_rhizo.txt", sep = ""), header=T, sep="\t")
input_table_2 <- table(input_table$KO)
input_table_3 <- names(input_table_2)[input_table_2 == 12]

#No dominator dataset
new_fam_data <- data.frame()

for (group in families){
  if(group == "other"){
    new_tax <- row.names(taxonomy_LjSC)[!taxonomy_LjSC$family %in% families]
  } else {
    new_tax <- row.names(taxonomy_LjSC)[taxonomy_LjSC$family == paste(group)]
  }
  
  KO_table_sub <- KO_table[,colnames(KO_table) %in% new_tax]
  KO_table_sub_2 <- rowSums(KO_table_sub)
  total_KOs <- length(names(KO_table_sub_2)[KO_table_sub_2 != 0])
  Family_KOs <- names(KO_table_sub_2)[KO_table_sub_2 != 0]
  
  total_266_Fam <- length(Family_KOs[Family_KOs %in% input_table_3])
  
  all_data <- data.frame(t(data.frame(c(paste(group), total_266_Fam))))
  new_fam_data <- rbind(new_fam_data, all_data)
}

row.names(new_fam_data) <- NULL
colnames(new_fam_data) <- c("Family", "Total no of 266 KOs")

norm_SSC_5$Group <- NA
norm_SSC_5$Group[norm_SSC_5$SynCom == "All other families drop out"] <- new_fam_data$`Total no of 266 KOs`[new_fam_data$Family == "other"]
norm_SSC_5$Group[norm_SSC_5$SynCom == "Pseudomonadaceae drop out"] <- new_fam_data$`Total no of 266 KOs`[new_fam_data$Family == "Pseudomonadaceae"]
norm_SSC_5$Group[norm_SSC_5$SynCom == "Caulobacteraceae drop out"] <- new_fam_data$`Total no of 266 KOs`[new_fam_data$Family == "Caulobacteraceae"]
norm_SSC_5$Group[norm_SSC_5$SynCom == "Burkholderiaceae drop out"] <- new_fam_data$`Total no of 266 KOs`[new_fam_data$Family == "Burkholderiaceae"]
norm_SSC_5$Group[norm_SSC_5$SynCom == "Full LjSC"] <- 266

norm_SSC_5$Group_2 <- paste(norm_SSC_5$SynCom, " (no of KO: " , norm_SSC_5$Group, ")",sep = "")

norm_SSC_5$Group_2 <- factor(norm_SSC_5$Group_2, levels = c("Burkholderiaceae drop out (no of KO: 141)", "Caulobacteraceae drop out (no of KO: 130)", "Pseudomonadaceae drop out (no of KO: 84)", "All other families drop out (no of KO: 242)", "Full LjSC (no of KO: 266)"))
# Plot the results
plot_col <- ggplot(norm_SSC_5, aes(x = Group_2, y = Bacterial_read_per_plant, colour = Group_2)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  ylab("Bacterial read/plant read") +
  geom_jitter(shape = 16, position = position_jitter(0.2), aes(colour = Group_2), show.legend = TRUE) +
  theme(
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    strip.text.x = element_text(size = 10)
  ) +
  guides(colour = FALSE) +
  scale_y_continuous(labels = scales::scientific,  ) +
  stat_summary(aes(label = Letters, y = max(Bacterial_read_per_plant)*1.03), fun = max, geom = "text") +
  ggtitle("Bacterial colonization with nodulators")
plot_col

plot_col_2 <- ggarrange(NULL, plot_col, ncol =2, widths=c(1,5))

#Without nodulators

norm_SSC =read.table(paste(working_directory,"LjSC_Family_drop_out_experiment/Data_with_synthetic_input.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
norm_SSC_2 <- norm_SSC[!row.names(norm_SSC) %in% c("LjNodule214", "LjRoot234", "LjRoot228"),]
norm_SSC_3 <- data.frame(colSums(norm_SSC_2))
colnames(norm_SSC_3) <- "No_of_bacterial_reads"
norm_SSC_3$Sample <- row.names(norm_SSC_3)
norm_SSC_4 <- norm_SSC_3[grep("ROOT", row.names(norm_SSC_3)),]
norm_SSC_4$No_of_plant_reads <- plant_reads_2$Host_mapped_Infr[match(norm_SSC_4$Sample, plant_reads_2$Sample)]
norm_SSC_4$Bacterial_read_per_plant <- norm_SSC_4$No_of_bacterial_reads/norm_SSC_4$No_of_plant_reads

norm_SSC_4_2 <- norm_SSC_4[!grepl("LDT7", norm_SSC_4$Sample),]
norm_SSC_5 <- norm_SSC_4_2[!grepl("LDT4", norm_SSC_4_2$Sample),]

norm_SSC_5$SynCom <- sapply(strsplit(as.vector(row.names(norm_SSC_5)), "_"),`[`, 1)

#Rhizobiaceae was not taken along as it was found to have a contamination
norm_SSC_5$SynCom <- gsub("LDT1", "Burkholderiaceae drop out", norm_SSC_5$SynCom)
norm_SSC_5$SynCom <- gsub("LDT2", "Caulobacteraceae drop out", norm_SSC_5$SynCom)
norm_SSC_5$SynCom <- gsub("LDT3", "Pseudomonadaceae drop out", norm_SSC_5$SynCom)
#norm_SSC_5$SynCom <- gsub("LDT4", "Rhizobiaceae drop out", norm_SSC_5$SynCom)
norm_SSC_5$SynCom <- gsub("LDT5", "All other families drop out", norm_SSC_5$SynCom)
norm_SSC_5$SynCom <- gsub("LDT6", "Full LjSC", norm_SSC_5$SynCom)

# Perform ANOVA
fitAnova <- aov(Bacterial_read_per_plant ~ SynCom, data = norm_SSC_5)

# Perform Tukey's post-hoc test
Tukey <- TukeyHSD(fitAnova)

# Get letters for significance
letters_anova <- multcompLetters4(fitAnova, Tukey)$SynCom$Letters

# Combine results into a data frame for plotting
ltlbl_combined <- data.frame(
  SynCom = names(letters_anova),
  Letters = letters_anova
)

# Ensure the order of the factor levels is correct
ltlbl_combined$SynCom <- factor(ltlbl_combined$SynCom, levels = unique(norm_SSC_5$SynCom))

# Apply the manually defined order to the factor
ltlbl_combined <- ltlbl_combined[order(ltlbl_combined$SynCom), ]

# Merge the letters with the data frame
norm_SSC_5 <- merge(norm_SSC_5, ltlbl_combined, by = "SynCom")

#Add number of 266 KOs
KO_table = read.table(paste(working_directory, "KO_genome/KO_LjSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

taxonomy <- read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",row.names =1)
taxonomy_LjSC <- taxonomy[taxonomy$SynCom == "LjSC",]

families <- c("Pseudomonadaceae", "Rhizobiaceae", "Burkholderiaceae", "Caulobacteraceae", "other")

input_table <- read.table(paste(working_directory,"DESeq2/Sig_KO_all_no_nod_rhizo.txt", sep = ""), header=T, sep="\t")
input_table_2 <- table(input_table$KO)
input_table_3 <- names(input_table_2)[input_table_2 == 12]

#No dominator dataset
new_fam_data <- data.frame()

for (group in families){
  if(group == "other"){
    new_tax <- row.names(taxonomy_LjSC)[!taxonomy_LjSC$family %in% families]
  } else {
    new_tax <- row.names(taxonomy_LjSC)[taxonomy_LjSC$family == paste(group)]
  }
  
  KO_table_sub <- KO_table[,colnames(KO_table) %in% new_tax]
  KO_table_sub_2 <- rowSums(KO_table_sub)
  total_KOs <- length(names(KO_table_sub_2)[KO_table_sub_2 != 0])
  Family_KOs <- names(KO_table_sub_2)[KO_table_sub_2 != 0]
  
  total_266_Fam <- length(Family_KOs[Family_KOs %in% input_table_3])
  
  all_data <- data.frame(t(data.frame(c(paste(group), total_266_Fam))))
  new_fam_data <- rbind(new_fam_data, all_data)
}

row.names(new_fam_data) <- NULL
colnames(new_fam_data) <- c("Family", "Total no of 266 KOs")

norm_SSC_5$Group <- NA
norm_SSC_5$Group[norm_SSC_5$SynCom == "All other families drop out"] <- new_fam_data$`Total no of 266 KOs`[new_fam_data$Family == "other"]
norm_SSC_5$Group[norm_SSC_5$SynCom == "Pseudomonadaceae drop out"] <- new_fam_data$`Total no of 266 KOs`[new_fam_data$Family == "Pseudomonadaceae"]
norm_SSC_5$Group[norm_SSC_5$SynCom == "Caulobacteraceae drop out"] <- new_fam_data$`Total no of 266 KOs`[new_fam_data$Family == "Caulobacteraceae"]
norm_SSC_5$Group[norm_SSC_5$SynCom == "Burkholderiaceae drop out"] <- new_fam_data$`Total no of 266 KOs`[new_fam_data$Family == "Burkholderiaceae"]
norm_SSC_5$Group[norm_SSC_5$SynCom == "Full LjSC"] <- 266

norm_SSC_5$Group_2 <- paste(norm_SSC_5$SynCom, " (no of KO: " , norm_SSC_5$Group, ")",sep = "")

norm_SSC_5$Group_2 <- factor(norm_SSC_5$Group_2, levels = c("Burkholderiaceae drop out (no of KO: 141)", "Caulobacteraceae drop out (no of KO: 130)", "Pseudomonadaceae drop out (no of KO: 84)", "All other families drop out (no of KO: 242)", "Full LjSC (no of KO: 266)"))
# Plot the results
plot_col_3 <- ggplot(norm_SSC_5, aes(x = Group_2, y = Bacterial_read_per_plant, colour = Group_2)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  ylab("Bacterial read/plant read") +
  geom_jitter(shape = 16, position = position_jitter(0.2), aes(colour = Group_2), show.legend = TRUE) +
  theme(
    axis.text.x = element_text(size = 12, hjust = 1, angle = 25),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    strip.text.x = element_text(size = 10)
  ) +
  guides(colour = FALSE) +
  scale_y_continuous(labels = scales::scientific,  ) +
  stat_summary(aes(label = Letters, y = max(Bacterial_read_per_plant)*1.03), fun = max, geom = "text") +
  ggtitle("Bacterial colonization without nodulators")
plot_col_3

plot_col_4 <- ggarrange(NULL, plot_col_3, ncol =2, widths=c(1,5))

plot_col_5 <- ggarrange(plot_col_2, plot_col_4, nrow =2, ncol =1, heights=c(2,3))

pdf(paste(results.dir,"Figure_S32_Bac_col_vs_plant.pdf", sep=""), width=10, height=10)
print(plot_col_5)
dev.off()
