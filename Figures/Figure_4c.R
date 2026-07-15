library("multcompView") #Version 0.1-9
library("ggplot2") #Version 3.4.2
library("grid") #Version 4.4.1
library("ggtern") #Version 3.4.1

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Figure 4c boxplots file preparation =====

# First Load sPLS-DA results
empty_vector_all_2_rhizo <- read.table(paste(working_directory, "sPLS-DA/output/PLSDA_Rhizobiaceae_with_dom.tsv", sep=""), sep = "\t",  header = T)
empty_vector_all_2_rhizo$Group <- "Rhizobiaceae_dominators_HvSC_LjSC_SSC"
empty_vector_all_2_Burk_SC <- read.table(paste(working_directory, "sPLS-DA/output/PLSDA_Burkholderiaceae_no_dom_SynCom.tsv", sep=""), sep = "\t",  header = T)
empty_vector_all_2_Burk_SC$Group <- "Burkholderiaceae_no_dominators"
empty_vector_all_2_Xant_SC <- read.table(paste(working_directory, "sPLS-DA/output/PLSDA_Xanthomonadaceae_no_dom_SynCom.tsv", sep=""), sep = "\t",  header = T)
empty_vector_all_2_Xant_SC$Group <- "Xanthomonadaceae_no_dominators"
empty_vector_all_2_Xant_PL <- read.table(paste(working_directory, "sPLS-DA/output/PLSDA_Xanthomonadaceae_no_dom_Plant.tsv", sep=""), sep = "\t",  header = T)
empty_vector_all_2_Xant_PL$Group <- "Xanthomonadaceae_no_dominators_LjSC"
empty_vector_all_2_Pseud <- read.table(paste(working_directory, "sPLS-DA/output/PLSDA_Pseudomonadaceae_no_dom_Plant.tsv", sep=""), sep = "\t",  header = T)
empty_vector_all_2_Pseud$Group <- "Pseudomonadacaea_no_dominators_LjSC"
empty_vector_all_2_Caulo <- read.table(paste(working_directory, "sPLS-DA/output/PLSDA_Caulobacteraceae_no_dom_Plant.tsv", sep=""), sep = "\t",  header = T)
empty_vector_all_2_Caulo$Group <- "Caulobacteraceae_no_dominators_AtSC"
empty_vector_all_2_Burk_PL <- read.table(paste(working_directory, "sPLS-DA/output/PLSDA_Burkholderiaceae_no_dom_Plant.tsv", sep=""), sep = "\t",  header = T)
empty_vector_all_2_Burk_PL$Group <- "Burkholderiaceae_no_dominators_AtSC"

empty_vector_all_2_all <- rbind(empty_vector_all_2_rhizo,empty_vector_all_2_Burk_SC,empty_vector_all_2_Xant_SC,empty_vector_all_2_Xant_PL,empty_vector_all_2_Pseud,empty_vector_all_2_Caulo,empty_vector_all_2_Burk_PL)
empty_vector_all_2_all$No_of_isolates <- NA
KO_SSC =read.table(paste(working_directory,"KO_genome/KO_SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

for (KO in unique(empty_vector_all_2_all$KO)){
  KO_SSC_sub <- KO_SSC[row.names(KO_SSC) == paste(KO),]
  KO_SSC_sub[KO_SSC_sub > 0] <- 1
  group <- sum(KO_SSC_sub)
  empty_vector_all_2_all$No_of_isolates[empty_vector_all_2_all$KO == paste(KO)] <- group
}

KOS <- table(empty_vector_all_2_all$KO)
KO_groups <- names(KOS)[KOS > 1]
empty_vector_all_2_all$Contribution <- (empty_vector_all_2_all$contrib * empty_vector_all_2_all$Variance) * 100
empty_vector_all_2_all <- empty_vector_all_2_all[empty_vector_all_2_all$Contribution >= 10,]

empty_vector_all_2_all$KO_2 <-  empty_vector_all_2_all$KO
empty_vector_all_2_all_sub <- empty_vector_all_2_all[!empty_vector_all_2_all$KO %in% KO_groups,]

for (KO in KO_groups){
  empty_vector_all_2_all_subber <- empty_vector_all_2_all[empty_vector_all_2_all$KO == paste(KO),]
  if(length(unique(empty_vector_all_2_all_subber$Data)) == 1 & length(unique(empty_vector_all_2_all_subber$Group)) == 1){
    empty_vector_all_2_all_subber_2 <- empty_vector_all_2_all_subber[1,]
    empty_vector_all_2_all_subber_2$contrib <- sum(empty_vector_all_2_all_subber$contrib)/length(empty_vector_all_2_all_subber$contrib)
    empty_vector_all_2_all_subber_2$Variance <- sum(empty_vector_all_2_all_subber$Variance)/length(empty_vector_all_2_all_subber$Variance)
    empty_vector_all_2_all_subber_2$Contribution <- sum(empty_vector_all_2_all_subber$Contribution)/length(empty_vector_all_2_all_subber$Contribution)
    empty_vector_all_2_all_sub <- rbind(empty_vector_all_2_all_sub,empty_vector_all_2_all_subber_2)
  } else if (length(unique(empty_vector_all_2_all_subber$Data)) == 1 & length(unique(empty_vector_all_2_all_subber$Group)) != 1){
    empty_vector_all_2_all_subber$KO <- paste(empty_vector_all_2_all_subber$KO, empty_vector_all_2_all_subber$Group, sep = "_")
    empty_vector_all_2_all_sub <- rbind(empty_vector_all_2_all_sub,empty_vector_all_2_all_subber)
  } else {
    empty_vector_all_2_all_subber$KO <- paste(empty_vector_all_2_all_subber$KO, empty_vector_all_2_all_subber$Data, sep = "_")
    empty_vector_all_2_all_sub <- rbind(empty_vector_all_2_all_sub,empty_vector_all_2_all_subber)
  }
}

###Figure 4c boxplots =====

#Rhizobiaceae with dominators
# Load the KO abundance tables (Rhizobiaceae-inclusive and exclusive datasets)
KO_SSC_only <- read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Rhizobiaceae_KO_with_dom.tsv", sep = ""),
                          header = TRUE, sep = "\t", row.names = 1)
KO_SSC_ex_only <- read.table(paste(working_directory, "sPLS-DA/isolate_drop_out_data/No_rhizobiaceae_KO_with_dom.tsv", sep = ""),
                             header = TRUE, sep = "\t", row.names = 1)

# Load metadata and filter for relevant samples
samples_df <- read.table(paste(working_directory, "SSC_R2_metadata_no_HL.tsv", sep = ""), header = TRUE, sep = "\t", row.names = 1)
colnames(samples_df)[5] <- "Nutrient"

samples_df_sub <- subset(samples_df, Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, Inoculum != "NS")
samples_df_sub_3 <- subset(samples_df_sub_2, Inoculum != "AtSC")

# Normalize KO tables by row sums (relative abundance)
KO_SSC_only_2 <- t(t(KO_SSC_only) / rowSums(t(KO_SSC_only)))
KO_SSC_ex_only_2 <- t(t(KO_SSC_ex_only) / rowSums(t(KO_SSC_ex_only)))

# Extract KOs associated with Rhizobiaceae and dominant families
KOs <- empty_vector_all_2_all_sub[empty_vector_all_2_all_sub$Group == "Rhizobiaceae_dominators_HvSC_LjSC_SSC", ]
KOs_Lj <- KOs$KO[KOs$Data == "Lj"]
KOs_Hv <- KOs$KO[KOs$Data == "Hv"]

# Subset for barley KOs
KO_SSC_only_2_Hv <- data.frame(KO_SSC_only_2[row.names(KO_SSC_only_2) %in% KOs_Hv, ])
colnames(KO_SSC_only_2_Hv) <- "KO"
KO_SSC_ex_only_2_Hv <- data.frame(KO_SSC_ex_only_2[row.names(KO_SSC_ex_only_2) %in% KOs_Hv, ])
colnames(KO_SSC_ex_only_2_Hv) <- "KO"

KO_SSC_only_2_Hv$Dataset <- "Rhizobiaceae"
KO_SSC_ex_only_2_Hv$Dataset <- "other families"

# Subset for Lotus KOs and sum across rows
KO_SSC_only_2_Lj <- KO_SSC_only_2[row.names(KO_SSC_only_2) %in% KOs_Lj, ]
KO_SSC_ex_only_2_Lj <- KO_SSC_ex_only_2[row.names(KO_SSC_ex_only_2) %in% KOs_Lj, ]

KO_SSC_only_2_Lj_2 <- data.frame(colSums(KO_SSC_only_2_Lj))
KO_SSC_ex_only_2_Lj_2 <- data.frame(colSums(KO_SSC_ex_only_2_Lj))

colnames(KO_SSC_only_2_Lj_2) <- "KO"
colnames(KO_SSC_ex_only_2_Lj_2) <- "KO"

KO_SSC_only_2_Lj_2$Dataset <- "Rhizobiaceae"
KO_SSC_ex_only_2_Lj_2$Dataset <- "Other families"

# Assign contributing plant information
KO_SSC_only_2_Hv$Contributing_plant <- "Barley KOs"
KO_SSC_ex_only_2_Hv$Contributing_plant <- "Barley KOs"
KO_SSC_only_2_Lj_2$Contributing_plant <- "Lotus KOs"
KO_SSC_ex_only_2_Lj_2$Contributing_plant <- "Lotus KOs"

# Assign sample names
KO_SSC_only_2_Hv$Sample <- row.names(KO_SSC_only_2_Hv)
KO_SSC_ex_only_2_Hv$Sample <- row.names(KO_SSC_ex_only_2_Hv)
KO_SSC_only_2_Lj_2$Sample <- row.names(KO_SSC_only_2_Lj_2)
KO_SSC_ex_only_2_Lj_2$Sample <- row.names(KO_SSC_ex_only_2_Lj_2)

# Function to subset and add plant information
subset_and_add_plant <- function(df) {
  df <- df[df$Sample %in% row.names(samples_df_sub_3), ]
  df$Plant <- samples_df_sub_3$Condition[match(df$Sample, row.names(samples_df_sub_3))]
  df$Plant <- as.factor(df$Plant)
  return(df)
}

# Apply the function to each data frame
KO_SSC_only_2_Hv <- subset_and_add_plant(KO_SSC_only_2_Hv)
KO_SSC_ex_only_2_Hv <- subset_and_add_plant(KO_SSC_ex_only_2_Hv)
KO_SSC_only_2_Lj_2 <- subset_and_add_plant(KO_SSC_only_2_Lj_2)
KO_SSC_ex_only_2_Lj_2 <- subset_and_add_plant(KO_SSC_ex_only_2_Lj_2)

# Function to perform ANOVA and Tukey test, and plot the results
anova_and_plot <- function(df, title) {
  # Perform ANOVA
  fitAnova <- aov(KO ~ Plant, data = df)
  
  # Perform Tukey's post-hoc test
  Tukey <- TukeyHSD(fitAnova)
  
  # Get letters for significance
  letters_anova <- multcompLetters4(fitAnova, Tukey)$Plant$Letters
  
  # Combine results into a data frame for plotting
  ltlbl_combined <- data.frame(
    Plant = names(letters_anova),
    Letters = letters_anova
  )
  
  # Ensure the order of the factor levels is correct
  ltlbl_combined$Plant <- factor(ltlbl_combined$Plant, levels = unique(df$Plant))
  
  # Apply the manually defined order to the factor
  ltlbl_combined <- ltlbl_combined[order(ltlbl_combined$Plant), ]
  
  # Merge the letters with the data frame
  df <- merge(df, ltlbl_combined, by = "Plant")
  
  # Plot the results
  plot <- ggplot(df, aes(x = Plant, y = KO, colour = Plant)) +
    geom_boxplot(outlier.shape = NA) +
    theme_classic() +
    ylab("Cumulative RA of KOs") +
    geom_jitter(shape = 16, position = position_jitter(0.2), aes(colour = Plant), show.legend = TRUE) +
    theme(
      axis.text.x = element_blank(),
      axis.title.y = element_text(size = 14),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14),
      plot.title = element_text(size = 14),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      strip.text.x = element_text(size = 10)
    ) +
    scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
    guides(colour = FALSE) +
    facet_wrap(~Contributing_plant + Dataset, scales = "free") +
    scale_y_continuous(labels = scales::scientific,  ) +
    stat_summary(aes(label = Letters, y = max(KO)*1.03), fun = max, geom = "text")
  
  return(plot)
}

# Generate plots for each data frame
plot_Hv_Rhiz <- anova_and_plot(KO_SSC_only_2_Hv, "Barley KOs - Rhizobiaceae")
plot_Hv_other <- anova_and_plot(KO_SSC_ex_only_2_Hv, "Barley KOs - Other Families")
plot_Lj_Rhiz <- anova_and_plot(KO_SSC_only_2_Lj_2, "Lotus KOs - Rhizobiaceae")
plot_Lj_other <- anova_and_plot(KO_SSC_ex_only_2_Lj_2, "Lotus KOs - Other Families")

# Arrange plots using ggarrange
rhizo_plot=grid.arrange(plot_Hv_Rhiz, plot_Hv_other, plot_Lj_Rhiz, plot_Lj_other, nrow = 2)

### Panel C-2 Xanthomonadaceae without dominators

# Load the KO abundance tables (Xanthomonadaceae-inclusive and exclusive datasets)
KO_SSC_only=read.table(paste(working_directory,"sPLS-DA/isolate_subset_data/Xanthomonadaceae_KO_no_dom.tsv",sep =""), header=TRUE,sep="\t", row.names = 1)
KO_SSC_ex_only=read.table(paste(working_directory,"sPLS-DA/isolate_drop_out_data/No_xanthomonadaceae_KO_no_dom.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

# Load metadata and filter for relevant samples
samples_df <- read.table(paste(working_directory, "SSC_R2_metadata_no_HL.tsv", sep = ""), header = TRUE, sep = "\t", row.names = 1)
colnames(samples_df)[5]="Nutrient"
samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

# Normalize KO tables by row sums (relative abundance)
KO_SSC_only_2 <- t(t(KO_SSC_only)/rowSums(t(KO_SSC_only)))
KO_SSC_ex_only_2 <- t(t(KO_SSC_ex_only)/rowSums(t(KO_SSC_ex_only)))

# Extract KOs associated with Xanthomonadaceae and without dominators families
KOs <- empty_vector_all_2_all_sub[empty_vector_all_2_all_sub$Group == "Xanthomonadaceae_no_dominators", ]
KOs_LjSC <- KOs$KO[KOs$Data == "LjSC"]
KOs_HvSC <- KOs$KO[KOs$Data == "HvSC"]

#Subset for HvSC related KOs
KO_SSC_only_2_HvSC <- data.frame(KO_SSC_only_2[row.names(KO_SSC_only_2) %in% KOs_HvSC,])
colnames(KO_SSC_only_2_HvSC) <- "KO"
KO_SSC_ex_only_2_HvSC <- data.frame(KO_SSC_ex_only_2[row.names(KO_SSC_ex_only_2) %in% KOs_HvSC,])
colnames(KO_SSC_ex_only_2_HvSC) <- "KO"

KO_SSC_ex_only_2_HvSC$Dataset <- "other families"
KO_SSC_only_2_HvSC$Dataset <- "Xanthomonadaceae"

#Subset for LjSC related KOs
KO_SSC_only_2_LjSC <- KO_SSC_only_2[row.names(KO_SSC_only_2) %in% KOs_LjSC,]
KO_SSC_ex_only_2_LjSC <- KO_SSC_ex_only_2[row.names(KO_SSC_ex_only_2) %in% KOs_LjSC,]

# sum across rows
KO_SSC_only_2_LjSC_2 <- data.frame(colSums(KO_SSC_only_2_LjSC))
KO_SSC_ex_only_2_LjSC_2 <- data.frame(colSums(KO_SSC_ex_only_2_LjSC))

colnames(KO_SSC_only_2_LjSC_2) <- "KO"
colnames(KO_SSC_ex_only_2_LjSC_2) <- "KO"

KO_SSC_only_2_LjSC_2$Dataset <- "Xanthomonadaceae"
KO_SSC_ex_only_2_LjSC_2$Dataset <- "other families"

# Assign contributing inoculum information
KO_SSC_only_2_HvSC$Contributing_plant <- "HvSC KOs"
KO_SSC_ex_only_2_HvSC$Contributing_plant <- "HvSC KOs"
KO_SSC_only_2_LjSC_2$Contributing_plant <- "LjSC KOs"
KO_SSC_ex_only_2_LjSC_2$Contributing_plant <- "LjSC KOs"

# Assign sample names
KO_SSC_only_2_HvSC$Sample <- row.names(KO_SSC_only_2_HvSC)
KO_SSC_ex_only_2_HvSC$Sample <- row.names(KO_SSC_ex_only_2_HvSC)
KO_SSC_only_2_LjSC_2$Sample <- row.names(KO_SSC_only_2_LjSC_2)
KO_SSC_ex_only_2_LjSC_2$Sample <- row.names(KO_SSC_ex_only_2_LjSC_2)

# Function to subset and add plant information
subset_and_add_inoc <- function(df) {
  df <- df[df$Sample %in% row.names(samples_df_sub_2), ]
  df$Plant <- samples_df_sub_2$Inoculum[match(df$Sample, row.names(samples_df_sub_2))]
  df$Plant <- as.factor(df$Plant)
  return(df)
}

# Apply the function to each data frame
KO_SSC_only_2_HvSC <- subset_and_add_inoc(KO_SSC_only_2_HvSC)
KO_SSC_only_2_LjSC_2 <- subset_and_add_inoc(KO_SSC_only_2_LjSC_2)
KO_SSC_ex_only_2_HvSC <- subset_and_add_inoc(KO_SSC_ex_only_2_HvSC)
KO_SSC_ex_only_2_LjSC_2 <- subset_and_add_inoc(KO_SSC_ex_only_2_LjSC_2)

# Function to perform ANOVA and Tukey test, and plot the results
anova_and_plot <- function(df, title) {
  # Perform ANOVA
  fitAnova <- aov(KO ~ Plant, data = df)
  
  # Perform Tukey's post-hoc test
  Tukey <- TukeyHSD(fitAnova)
  
  # Get letters for significance
  letters_anova <- multcompLetters4(fitAnova, Tukey)$Plant$Letters
  
  # Combine results into a data frame for plotting
  ltlbl_combined <- data.frame(
    Plant = names(letters_anova),
    Letters = letters_anova
  )
  
  # Ensure the order of the factor levels is correct
  ltlbl_combined$Plant <- factor(ltlbl_combined$Plant, levels = unique(df$Plant))
  
  # Apply the manually defined order to the factor
  ltlbl_combined <- ltlbl_combined[order(ltlbl_combined$Plant), ]
  
  # Merge the letters with the data frame
  df <- merge(df, ltlbl_combined, by = "Plant")
  
  # Plot the results
  plot <- ggplot(df, aes(x = Plant, y = KO, colour = Plant)) +
    geom_boxplot(outlier.shape = NA) +
    theme_classic() +
    ylab("Cumulative RA of KOs") +
    geom_jitter(shape = 16, position = position_jitter(0.2), aes(colour = Plant), show.legend = TRUE) +
    theme(
      axis.text.x = element_blank(),
      axis.title.y = element_text(size = 14),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14),
      plot.title = element_text(size = 14),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      strip.text.x = element_text(size = 10)
    ) +
    scale_color_manual(values = c("#A3A500","#00B0F6","#00BF7D","#F8766D" )) +
    guides(colour = FALSE) +
    facet_wrap(~Contributing_plant + Dataset, scales = "free") +
    scale_y_continuous(labels = scales::scientific,  ) +
    stat_summary(aes(label = Letters, y = max(KO)*1.03), fun = max, geom = "text")
  
  return(plot)
}

# Generate plots for each data frame
plot_HvSC <- anova_and_plot(KO_SSC_only_2_HvSC, "HvSC KOs - Xanthomonadaceae")
plot_HvSC_other <- anova_and_plot(KO_SSC_ex_only_2_HvSC, "HvSC KOs - Other Families")
plot_LjSC <- anova_and_plot(KO_SSC_only_2_LjSC_2, "Lotus KOs - Rhizobiaceae")
plot_LjSC_other <- anova_and_plot(KO_SSC_ex_only_2_LjSC_2, "Lotus KOs - Other Families")

# Arrange plots using ggarrange
grid.arrange(plot_HvSC, plot_HvSC_other, plot_LjSC, plot_LjSC_other, nrow = 2)

### Panel C-3 Pseudomonadaceae LjSC without dominators (not follwing expected trend)

KO_SSC_only=read.table(paste(working_directory,"sPLS-DA/isolate_subset_data/Pseudomonadaceae_KO_no_dom.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
KO_SSC_ex_only=read.table(paste(working_directory,"sPLS-DA/isolate_drop_out_data/No_pseudomonadaceae_KO_no_dom.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Samples TABLE
samples_df <- read.table(paste(working_directory, "SSC_R2_metadata_no_HL.tsv", sep = ""), header = TRUE, sep = "\t", row.names = 1)
colnames(samples_df)[5]="Nutrient"
samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")
samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum == "LjSC")

KO_SSC_only_2 <- t(t(KO_SSC_only)/rowSums(t(KO_SSC_only)))
KO_SSC_ex_only_2 <- t(t(KO_SSC_ex_only)/rowSums(t(KO_SSC_ex_only)))

KOs <- empty_vector_all_2_all_sub[empty_vector_all_2_all_sub$Group == "Pseudomonadacaea_no_dominators_LjSC", ]
KOs_Lj <- KOs$KO[KOs$Data == "Lj"]
KOs_At <- KOs$KO[KOs$Data == "At"]

#Arabidopsis 
KO_SSC_only_2_At <- data.frame(KO_SSC_only_2[row.names(KO_SSC_only_2) %in% KOs_At,])
colnames(KO_SSC_only_2_At) <- "KO"
KO_SSC_ex_only_2_At <- data.frame(KO_SSC_ex_only_2[row.names(KO_SSC_ex_only_2) %in% KOs_At,])
colnames(KO_SSC_ex_only_2_At) <- "KO"

KO_SSC_ex_only_2_At$Dataset <- "Other families LjSC"
KO_SSC_only_2_At$Dataset <- "Pseudomonadaceae LjSC"

#Lotus
KO_SSC_only_2_Lj <- KO_SSC_only_2[row.names(KO_SSC_only_2) %in% KOs_Lj,]
KO_SSC_ex_only_2_Lj <- KO_SSC_ex_only_2[row.names(KO_SSC_ex_only_2) %in% KOs_Lj,]

KO_SSC_only_2_Lj_2 <- data.frame(colSums(KO_SSC_only_2_Lj))
KO_SSC_ex_only_2_Lj_2 <- data.frame(colSums(KO_SSC_ex_only_2_Lj))

colnames(KO_SSC_only_2_Lj_2) <- "KO"
colnames(KO_SSC_ex_only_2_Lj_2) <- "KO"

KO_SSC_only_2_Lj_2$Dataset <- "Pseudomonadaceae LjSC"
KO_SSC_ex_only_2_Lj_2$Dataset <- "Other families LjSC"

KO_SSC_only_2_At$Contributing_plant <- "Arabidopsis KOs"
KO_SSC_ex_only_2_At$Contributing_plant <- "Arabidopsis KOs"
KO_SSC_only_2_Lj_2$Contributing_plant <- "Lotus KOs"
KO_SSC_ex_only_2_Lj_2$Contributing_plant <- "Lotus KOs"

KO_SSC_only_2_At$Sample <- row.names(KO_SSC_only_2_At)
KO_SSC_ex_only_2_At$Sample <- row.names(KO_SSC_ex_only_2_At)
KO_SSC_only_2_Lj_2$Sample <- row.names(KO_SSC_only_2_Lj_2)
KO_SSC_ex_only_2_Lj_2$Sample <- row.names(KO_SSC_ex_only_2_Lj_2)

# Function to subset and add plant information
subset_and_add_plant <- function(df) {
  df <- df[df$Sample %in% row.names(samples_df_sub_3), ]
  df$Plant <- samples_df_sub_3$Condition[match(df$Sample, row.names(samples_df_sub_3))]
  df$Plant <- as.factor(df$Plant)
  return(df)
}

# Apply the function to each data frame
KO_SSC_only_2_At <- subset_and_add_plant(KO_SSC_only_2_At)
KO_SSC_ex_only_2_At <- subset_and_add_plant(KO_SSC_ex_only_2_At)
KO_SSC_only_2_Lj_2 <- subset_and_add_plant(KO_SSC_only_2_Lj_2)
KO_SSC_ex_only_2_Lj_2 <- subset_and_add_plant(KO_SSC_ex_only_2_Lj_2)

# Function to perform ANOVA and Tukey test, and plot the results
anova_and_plot <- function(df, title) {
  # Perform ANOVA
  fitAnova <- aov(KO ~ Plant, data = df)
  
  # Perform Tukey's post-hoc test
  Tukey <- TukeyHSD(fitAnova)
  
  # Get letters for significance
  letters_anova <- multcompLetters4(fitAnova, Tukey)$Plant$Letters
  
  # Combine results into a data frame for plotting
  ltlbl_combined <- data.frame(
    Plant = names(letters_anova),
    Letters = letters_anova
  )
  
  # Ensure the order of the factor levels is correct
  ltlbl_combined$Plant <- factor(ltlbl_combined$Plant, levels = unique(df$Plant))
  
  # Apply the manually defined order to the factor
  ltlbl_combined <- ltlbl_combined[order(ltlbl_combined$Plant), ]
  
  # Merge the letters with the data frame
  df <- merge(df, ltlbl_combined, by = "Plant")
  
  # Plot the results
  plot <- ggplot(df, aes(x = Plant, y = KO, colour = Plant)) +
    geom_boxplot(outlier.shape = NA) +
    theme_classic() +
    ylab("Cumulative RA of KOs") +
    geom_jitter(shape = 16, position = position_jitter(0.2), aes(colour = Plant), show.legend = TRUE) +
    theme(
      axis.text.x = element_blank(),
      axis.title.y = element_text(size = 14),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14),
      plot.title = element_text(size = 14),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      strip.text.x = element_text(size = 10)
    ) +
    scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
    guides(colour = FALSE) +
    facet_wrap(~Contributing_plant + Dataset, scales = "free") +
    scale_y_continuous(labels = scales::scientific,  ) +
    stat_summary(aes(label = Letters, y = max(KO)*1.03), fun = max, geom = "text")
  
  return(plot)
}

# Generate plots for each data frame
plot_At_Pseudo <- anova_and_plot(KO_SSC_only_2_At, "Arabidopsis KOs - Pseudomonadaceae")
plot_At_other <- anova_and_plot(KO_SSC_ex_only_2_At, "Arabidopsis KOs - Other Families")
plot_Lj_Pseudo <- anova_and_plot(KO_SSC_only_2_Lj_2, "Lotus KOs - Pseudomonadaceae")
plot_Lj_2_other <- anova_and_plot(KO_SSC_ex_only_2_Lj_2, "Lotus KOs - Other Families")

# Arrange plots using ggarrange
grid.arrange(plot_At_Pseudo, plot_At_other, plot_Lj_Pseudo, plot_Lj_2_other, nrow = 2)

### Combine the 3 panels
combined_plot_2 <- 
  (plot_Hv_Rhiz + plot_Hv_other + plot_Lj_Rhiz + plot_Lj_other) / 
  (plot_HvSC + plot_HvSC_other + plot_LjSC + plot_LjSC_other) / 
  (plot_At_Pseudo + plot_At_other + plot_Lj_Pseudo + plot_Lj_2_other) +
  plot_layout(ncol = 3)

print(combined_plot_2)

#Save plot
pdf(paste(results.dir, "Figure_4c_boxplots.pdf", sep = ""), width = 21, height = 8)
print(combined_plot_2)
dev.off()
