library("multcompView") #Version 0.1-9
library("ggplot2") #Version 3.4.2
library("ggtern") #Version 3.4.1

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Figure S16 boxplots file preparation =====

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

###Figure S16 - boxplots ======

#Burkholderiaceae AtSC no dominators - Plant

# Load the KO abundance tables (Burkholderiaceae-inclusive and exclusive datasets)
KO_SSC_only <- read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Burkholderiaceae_KO_no_dom.tsv", sep = ""),
                          header = TRUE, sep = "\t", row.names = 1)
KO_SSC_ex_only <- read.table(paste(working_directory, "sPLS-DA/isolate_drop_out_data/No_burkholderiaceae_KO_no_dom.tsv", sep = ""),
                             header = TRUE, sep = "\t", row.names = 1)

# Load metadata and filter for relevant samples
samples_df <- read.table(paste(working_directory, "SSC_R2_metadata_no_HL.tsv", sep = ""), header = TRUE, sep = "\t", row.names = 1)
colnames(samples_df)[5] <- "Nutrient"

samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")
samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum == "AtSC")

# Normalize KO tables by row sums (relative abundance)
KO_SSC_only_2 <- t(t(KO_SSC_only) / rowSums(t(KO_SSC_only)))
KO_SSC_ex_only_2 <- t(t(KO_SSC_ex_only) / rowSums(t(KO_SSC_ex_only)))

# Extract KOs associated with Burkholderiaceae in AtSC
KOs <- empty_vector_all_2_all_sub[empty_vector_all_2_all_sub$Group == "Burkholderiaceae_no_dominators_AtSC", ]
KOs_Lj <- KOs$KO[KOs$Data == "Lj"]
KOs_Hv <- KOs$KO[KOs$Data == "Hv"]
KOs_At <- KOs$KO[KOs$Data == "At"]

# Subset for Arabidopsis KOs
KO_SSC_only_2_At <- data.frame(KO_SSC_only_2[row.names(KO_SSC_only_2) %in% KOs_At, ])
colnames(KO_SSC_only_2_At) <- "KO"
KO_SSC_ex_only_2_At <- data.frame(KO_SSC_ex_only_2[row.names(KO_SSC_ex_only_2) %in% KOs_At, ])
colnames(KO_SSC_ex_only_2_At) <- "KO"

KO_SSC_only_2_At$Dataset <- "Burkholderiaceae AtSC"
KO_SSC_ex_only_2_At$Dataset <- "Other families AtSC" 

# Subset for Barley KOs and sum across rows
KO_SSC_only_2_Hv <- KO_SSC_only_2[row.names(KO_SSC_only_2) %in% KOs_Hv, ]
KO_SSC_ex_only_2_Hv <- KO_SSC_ex_only_2[row.names(KO_SSC_ex_only_2) %in% KOs_Hv, ]

KO_SSC_only_2_Hv_2 <- data.frame(colSums(KO_SSC_only_2_Hv))
KO_SSC_ex_only_2_Hv_2 <- data.frame(colSums(KO_SSC_ex_only_2_Hv))

colnames(KO_SSC_only_2_Hv_2) <- "KO"
colnames(KO_SSC_ex_only_2_Hv_2) <- "KO"

KO_SSC_only_2_Hv_2$Dataset <- "Burkholderiaceae AtSC"
KO_SSC_ex_only_2_Hv_2$Dataset <- "Other families AtSC"

# Subset for Lotus KOs and sum across rows
KO_SSC_only_2_Lj <- KO_SSC_only_2[row.names(KO_SSC_only_2) %in% KOs_Lj, ]
KO_SSC_ex_only_2_Lj <- KO_SSC_ex_only_2[row.names(KO_SSC_ex_only_2) %in% KOs_Lj, ]

KO_SSC_only_2_Lj_2 <- data.frame(colSums(KO_SSC_only_2_Lj))
KO_SSC_ex_only_2_Lj_2 <- data.frame(colSums(KO_SSC_ex_only_2_Lj))

colnames(KO_SSC_only_2_Lj_2) <- "KO"
colnames(KO_SSC_ex_only_2_Lj_2) <- "KO"

KO_SSC_only_2_Lj_2$Dataset <- "Burkholderiaceae AtSC"
KO_SSC_ex_only_2_Lj_2$Dataset <- "Other families AtSC"

# Assign contributing plant information
KO_SSC_only_2_At$Contributing_plant <- "Arabidopsis KOs"
KO_SSC_ex_only_2_At$Contributing_plant <- "Arabidopsis KOs"
KO_SSC_only_2_Hv_2$Contributing_plant <- "Barley KOs"
KO_SSC_ex_only_2_Hv_2$Contributing_plant <- "Barley KOs"
KO_SSC_only_2_Lj_2$Contributing_plant <- "Lotus KOs"
KO_SSC_ex_only_2_Lj_2$Contributing_plant <- "Lotus KOs"

# Assign sample names
KO_SSC_only_2_At$Sample <- row.names(KO_SSC_only_2_At)
KO_SSC_ex_only_2_At$Sample <- row.names(KO_SSC_ex_only_2_At)
KO_SSC_only_2_Hv_2$Sample <- row.names(KO_SSC_only_2_Hv_2)
KO_SSC_ex_only_2_Hv_2$Sample <- row.names(KO_SSC_ex_only_2_Hv_2)
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
KO_SSC_only_2_Hv_2 <- subset_and_add_plant(KO_SSC_only_2_Hv_2)
KO_SSC_ex_only_2_Hv_2 <- subset_and_add_plant(KO_SSC_ex_only_2_Hv_2)
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
plot_At_Burk <- anova_and_plot(KO_SSC_only_2_At, "Arabidopsis KOs - Burkholderiaceae AtSC")
plot_At_other <- anova_and_plot(KO_SSC_ex_only_2_At, "Arabidopsis KOs - Other families AtSC")
plot_Hv_Burk <- anova_and_plot(KO_SSC_only_2_Hv_2, "Barley KOs - Burkholderiaceae AtSC")
plot_Hv_other <- anova_and_plot(KO_SSC_ex_only_2_Hv_2, "Barley KOs - Other families AtSC")
plot_Lj_Burk <- anova_and_plot(KO_SSC_only_2_Lj_2, "Lotus KOs - Burkholderiaceae AtSC")
plot_Lj_other <- anova_and_plot(KO_SSC_ex_only_2_Lj_2, "Lotus KOs - Other families AtSC")

# Arrange plots using ggarrange
Burk_plot=grid.arrange(plot_At_Burk,plot_At_other,plot_Hv_Burk, plot_Hv_other, plot_Lj_Burk, plot_Lj_other, nrow = 3)

#Caulobacteraceae AtSC without dominators - Plant
# Load the KO abundance tables (Caulobacteraceae-inclusive and exclusive datasets)
KO_SSC_only <- read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Caulobacteraceae_KO_no_dom.tsv", sep = ""),
                          header = TRUE, sep = "\t", row.names = 1)
KO_SSC_ex_only <- read.table(paste(working_directory, "sPLS-DA/isolate_drop_out_data/No_caulobacteraceae_KO_no_dom.tsv", sep = ""),
                             header = TRUE, sep = "\t", row.names = 1)

# Load metadata and filter for relevant samples
samples_df <- read.table(paste(working_directory, "SSC_R2_metadata_no_HL.tsv", sep = ""), header = TRUE, sep = "\t", row.names = 1)
colnames(samples_df)[5] <- "Nutrient"

samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")
samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum == "AtSC")

# Normalize KO tables by row sums (relative abundance)
KO_SSC_only_2 <- t(t(KO_SSC_only) / rowSums(t(KO_SSC_only)))
KO_SSC_ex_only_2 <- t(t(KO_SSC_ex_only) / rowSums(t(KO_SSC_ex_only)))

# Extract KOs associated with Caulobacteraceae in AtSC
KOs <- empty_vector_all_2_all_sub[empty_vector_all_2_all_sub$Group == "Caulobacteraceae_no_dominators_AtSC", ]
KOs_Hv <- KOs$KO[KOs$Data == "Hv"]
KOs_At <- KOs$KO[KOs$Data == "At"]

# Subset for Arabidopsis KOs and sum across rows
KO_SSC_only_2_At <- KO_SSC_only_2[row.names(KO_SSC_only_2) %in% KOs_At, ]
KO_SSC_ex_only_2_At <- KO_SSC_ex_only_2[row.names(KO_SSC_ex_only_2) %in% KOs_At, ]

KO_SSC_only_2_At_2 <- data.frame(colSums(KO_SSC_only_2_At))
KO_SSC_ex_only_2_At_2 <- data.frame(colSums(KO_SSC_ex_only_2_At))

colnames(KO_SSC_only_2_At_2) <- "KO"
colnames(KO_SSC_ex_only_2_At_2) <- "KO"

KO_SSC_only_2_At_2$Dataset <- "Caulobacteraceae AtSC"
KO_SSC_ex_only_2_At_2$Dataset <- "Other families AtSC"

# Subset for Barley KOs and sum across rows
KO_SSC_only_2_Hv <- KO_SSC_only_2[row.names(KO_SSC_only_2) %in% KOs_Hv, ]
KO_SSC_ex_only_2_Hv <- KO_SSC_ex_only_2[row.names(KO_SSC_ex_only_2) %in% KOs_Hv, ]

KO_SSC_only_2_Hv_2 <- data.frame(colSums(KO_SSC_only_2_Hv))
KO_SSC_ex_only_2_Hv_2 <- data.frame(colSums(KO_SSC_ex_only_2_Hv))

colnames(KO_SSC_only_2_Hv_2) <- "KO"
colnames(KO_SSC_ex_only_2_Hv_2) <- "KO"

KO_SSC_only_2_Hv_2$Dataset <- "Caulobacteraceae AtSC"
KO_SSC_ex_only_2_Hv_2$Dataset <- "Other families AtSC"

# Assign contributing plant information
KO_SSC_only_2_At_2$Contributing_plant <- "Arabidopsis KOs"
KO_SSC_ex_only_2_At_2$Contributing_plant <- "Arabidopsis KOs"
KO_SSC_only_2_Hv_2$Contributing_plant <- "Barley KOs"
KO_SSC_ex_only_2_Hv_2$Contributing_plant <- "Barley KOs"

# Assign sample names
KO_SSC_only_2_At_2$Sample <- row.names(KO_SSC_only_2_At_2)
KO_SSC_ex_only_2_At_2$Sample <- row.names(KO_SSC_ex_only_2_At_2)
KO_SSC_only_2_Hv_2$Sample <- row.names(KO_SSC_only_2_Hv_2)
KO_SSC_ex_only_2_Hv_2$Sample <- row.names(KO_SSC_ex_only_2_Hv_2)

# Function to subset and add plant information
subset_and_add_plant <- function(df) {
  df <- df[df$Sample %in% row.names(samples_df_sub_3), ]
  df$Plant <- samples_df_sub_3$Condition[match(df$Sample, row.names(samples_df_sub_3))]
  df$Plant <- as.factor(df$Plant)
  return(df)
}

# Apply the function to each data frame
KO_SSC_only_2_At_2 <- subset_and_add_plant(KO_SSC_only_2_At_2)
KO_SSC_ex_only_2_At_2 <- subset_and_add_plant(KO_SSC_ex_only_2_At_2)
KO_SSC_only_2_Hv_2 <- subset_and_add_plant(KO_SSC_only_2_Hv_2)
KO_SSC_ex_only_2_Hv_2 <- subset_and_add_plant(KO_SSC_ex_only_2_Hv_2)

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
plot_At_Caulo <- anova_and_plot(KO_SSC_only_2_At_2, "Arabidopsis KOs - Caulobacteraceae AtSC")
plot_At_other_caulo <- anova_and_plot(KO_SSC_ex_only_2_At_2, "Arabidopsis KOs - Other families AtSC")
plot_Hv_Caulo <- anova_and_plot(KO_SSC_only_2_Hv_2, "Barley KOs - Caulobacteraceae AtSC")
plot_Hv_other_caulo <- anova_and_plot(KO_SSC_ex_only_2_Hv_2, "Barley KOs - Other families AtSC")

# Arrange plots using ggarrange
Caulo_plot=grid.arrange(plot_At_Caulo,plot_At_other_caulo,plot_Hv_Caulo, plot_Hv_other_caulo,nrow = 2)

#Xanthomonadaceae LjSC without dominators - Plant
# Load the KO abundance tables (Xanthomonadaceae-inclusive and exclusive datasets)
KO_SSC_only <- read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Xanthomonadaceae_KO_no_dom.tsv", sep = ""),
                          header = TRUE, sep = "\t", row.names = 1)
KO_SSC_ex_only <- read.table(paste(working_directory, "sPLS-DA/isolate_drop_out_data/No_xanthomonadaceae_KO_no_dom.tsv", sep = ""),
                             header = TRUE, sep = "\t", row.names = 1)

# Load metadata and filter for relevant samples
samples_df <- read.table(paste(working_directory, "SSC_R2_metadata_no_HL.tsv", sep = ""), header = TRUE, sep = "\t", row.names = 1)
colnames(samples_df)[5] <- "Nutrient"

samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")
samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum == "LjSC")

# Normalize KO tables by row sums (relative abundance)
KO_SSC_only_2 <- t(t(KO_SSC_only) / rowSums(t(KO_SSC_only)))
KO_SSC_ex_only_2 <- t(t(KO_SSC_ex_only) / rowSums(t(KO_SSC_ex_only)))

# Extract KOs associated with Xanthomonadaceae and dominant families
KOs <- empty_vector_all_2_all_sub[empty_vector_all_2_all_sub$Group == "Xanthomonadaceae_no_dominators_LjSC", ]
KOs_Hv <- KOs$KO[KOs$Data == "Hv"]
KOs_At <- KOs$KO[KOs$Data == "At"]

# Subset for Arabidopsis KOs
KO_SSC_only_2_At <- data.frame(KO_SSC_only_2[row.names(KO_SSC_only_2) %in% KOs_At, ])
colnames(KO_SSC_only_2_At) <- "KO"
KO_SSC_ex_only_2_At <- data.frame(KO_SSC_ex_only_2[row.names(KO_SSC_ex_only_2) %in% KOs_At, ])
colnames(KO_SSC_ex_only_2_At) <- "KO"

KO_SSC_only_2_At$Dataset <- "Xanthomonadaceae LjSC"
KO_SSC_ex_only_2_At$Dataset <- "other families LjSC"

# Subset for Barley KOs and sum across rows
KO_SSC_only_2_Hv <- KO_SSC_only_2[row.names(KO_SSC_only_2) %in% KOs_Hv, ]
KO_SSC_ex_only_2_Hv <- KO_SSC_ex_only_2[row.names(KO_SSC_ex_only_2) %in% KOs_Hv, ]

KO_SSC_only_2_Hv_2 <- data.frame(colSums(KO_SSC_only_2_Hv))
KO_SSC_ex_only_2_Hv_2 <- data.frame(colSums(KO_SSC_ex_only_2_Hv))

colnames(KO_SSC_only_2_Hv_2) <- "KO"
colnames(KO_SSC_ex_only_2_Hv_2) <- "KO"

KO_SSC_only_2_Hv_2$Dataset <- "Xanthomonadaceae LjSC"
KO_SSC_ex_only_2_Hv_2$Dataset <- "Other families LjSC"

# Assign contributing plant information
KO_SSC_only_2_At$Contributing_plant <- "Arabidopsis KOs"
KO_SSC_ex_only_2_At$Contributing_plant <- "Arabidopsis KOs"
KO_SSC_only_2_Hv_2$Contributing_plant <- "Barley KOs"
KO_SSC_ex_only_2_Hv_2$Contributing_plant <- "Barley KOs"

# Assign sample names
KO_SSC_only_2_At$Sample <- row.names(KO_SSC_only_2_At)
KO_SSC_ex_only_2_At$Sample <- row.names(KO_SSC_ex_only_2_At)
KO_SSC_only_2_Hv_2$Sample <- row.names(KO_SSC_only_2_Hv_2)
KO_SSC_ex_only_2_Hv_2$Sample <- row.names(KO_SSC_ex_only_2_Hv_2)

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
KO_SSC_only_2_Hv_2 <- subset_and_add_plant(KO_SSC_only_2_Hv_2)
KO_SSC_ex_only_2_Hv_2 <- subset_and_add_plant(KO_SSC_ex_only_2_Hv_2)

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
plot_At_Xanth <- anova_and_plot(KO_SSC_only_2_At, "Arabidopsis KOs - Xanthomonadaceae LjSC")
plot_At_other_xanth <- anova_and_plot(KO_SSC_ex_only_2_At, "Arabidopsis KOs - Other families LjSC")
plot_Hv_Xanth <- anova_and_plot(KO_SSC_only_2_Hv_2, "Barley KOs - Xanthomonadaceae LjSC")
plot_Hv_other_xanth <- anova_and_plot(KO_SSC_ex_only_2_Hv_2, "Barley KOs - Other families LjSC")

# Arrange plots using ggarrange
Xanth_plot=grid.arrange(plot_At_Xanth,plot_At_other_xanth,plot_Hv_Xanth, plot_Hv_other_xanth, nrow = 2)

#Burkholderiaceae no dominators - SynComs
KO_SSC_only <- read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Burkholderiaceae_KO_no_dom.tsv", sep = ""),
                          header = TRUE, sep = "\t", row.names = 1)
KO_SSC_ex_only <- read.table(paste(working_directory, "sPLS-DA/isolate_drop_out_data/No_burkholderiaceae_KO_no_dom.tsv", sep = ""),
                             header = TRUE, sep = "\t", row.names = 1)

# Load metadata and filter for relevant samples
samples_df <- read.table(paste(working_directory, "SSC_R2_metadata_no_HL.tsv", sep = ""), header = TRUE, sep = "\t", row.names = 1)
colnames(samples_df)[5] <- "Nutrient"

samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

# Normalize KO tables by row sums (relative abundance)
KO_SSC_only_2 <- t(t(KO_SSC_only) / rowSums(t(KO_SSC_only)))
KO_SSC_ex_only_2 <- t(t(KO_SSC_ex_only) / rowSums(t(KO_SSC_ex_only)))

# Extract KOs associated with Burkholderiaceae and dominant families
KOs <- empty_vector_all_2_all_sub[empty_vector_all_2_all_sub$Group == "Burkholderiaceae_no_dominators", ]
KOs_AtSC <- KOs$KO[KOs$Data == "AtSC"]
KOs_LjSC <- KOs$KO[KOs$Data == "LjSC"]

# Subset for AtSC KOs
KO_SSC_only_2_AtSC <- data.frame(KO_SSC_only_2[row.names(KO_SSC_only_2) %in% KOs_AtSC, ])
colnames(KO_SSC_only_2_AtSC ) <- "KO"
KO_SSC_ex_only_2_AtSC <- data.frame(KO_SSC_ex_only_2[row.names(KO_SSC_ex_only_2) %in% KOs_AtSC, ])
colnames(KO_SSC_ex_only_2_AtSC) <- "KO"

KO_SSC_only_2_AtSC$Dataset <- "Burkholderiaceae"
KO_SSC_ex_only_2_AtSC$Dataset <- "Other families" 

# Subset for LjSC KOs
KO_SSC_only_2_LjSC <- data.frame(KO_SSC_only_2[row.names(KO_SSC_only_2) %in% KOs_LjSC, ])
colnames(KO_SSC_only_2_LjSC ) <- "KO"
KO_SSC_ex_only_2_LjSC <- data.frame(KO_SSC_ex_only_2[row.names(KO_SSC_ex_only_2) %in% KOs_LjSC, ])
colnames(KO_SSC_ex_only_2_LjSC) <- "KO"

KO_SSC_only_2_LjSC$Dataset <- "Burkholderiaceae"
KO_SSC_ex_only_2_LjSC$Dataset <- "Other families" 

# Assign contributing plant information
KO_SSC_only_2_AtSC$Contributing_SynCom <- "AtSC KOs"
KO_SSC_ex_only_2_AtSC$Contributing_SynCom <- "AtSC KOs"
KO_SSC_only_2_LjSC$Contributing_SynCom <- "LjSC KOs"
KO_SSC_ex_only_2_LjSC$Contributing_SynCom <- "LjSC KOs"

# Assign sample names
KO_SSC_only_2_AtSC$Sample <- row.names(KO_SSC_only_2_AtSC)
KO_SSC_ex_only_2_AtSC$Sample <- row.names(KO_SSC_ex_only_2_AtSC)
KO_SSC_only_2_LjSC$Sample <- row.names(KO_SSC_only_2_LjSC)
KO_SSC_ex_only_2_LjSC$Sample <- row.names(KO_SSC_ex_only_2_LjSC)

# Function to subset and add plant information
subset_and_add_syncom <- function(df) {
  df <- df[df$Sample %in% row.names(samples_df_sub_2), ]
  df$SynCom <- samples_df_sub_2$Inoculum[match(df$Sample, row.names(samples_df_sub_2))]
  df$SynCom <- as.factor(df$SynCom)
  return(df)
}

# Apply the function to each data frame
KO_SSC_only_2_AtSC <- subset_and_add_syncom(KO_SSC_only_2_AtSC)
KO_SSC_ex_only_2_AtSC <- subset_and_add_syncom(KO_SSC_ex_only_2_AtSC)
KO_SSC_only_2_LjSC <- subset_and_add_syncom(KO_SSC_only_2_LjSC)
KO_SSC_ex_only_2_LjSC <- subset_and_add_syncom(KO_SSC_ex_only_2_LjSC)

# Function to perform ANOVA and Tukey test, and plot the results
anova_and_plot <- function(df, title) {
  # Perform ANOVA
  fitAnova <- aov(KO ~ SynCom, data = df)
  
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
  ltlbl_combined$SynCom <- factor(ltlbl_combined$SynCom, levels = unique(df$SynCom))
  
  # Apply the manually defined order to the factor
  ltlbl_combined <- ltlbl_combined[order(ltlbl_combined$SynCom), ]
  
  # Merge the letters with the data frame
  df <- merge(df, ltlbl_combined, by = "SynCom")
  
  # Plot the results
  plot <- ggplot(df, aes(x = SynCom, y = KO, colour = SynCom)) +
    geom_boxplot(outlier.shape = NA) +
    theme_classic() +
    ylab("Cumulative RA of KOs") +
    geom_jitter(shape = 16, position = position_jitter(0.2), aes(colour = SynCom), show.legend = TRUE) +
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
    facet_wrap(~Contributing_SynCom + Dataset, scales = "free") +
    scale_y_continuous(labels = scales::scientific,  ) +
    stat_summary(aes(label = Letters, y = max(KO)*1.03), fun = max, geom = "text")
  
  return(plot)
}

# Generate plots for each data frame
plot_AtSC_Burk <- anova_and_plot(KO_SSC_only_2_AtSC, "AtSC KOs - Burkholderiaceae")
plot_AtSC_other <- anova_and_plot(KO_SSC_ex_only_2_AtSC, "AtSC KOs - Other families")
plot_LjSC_Burk <- anova_and_plot(KO_SSC_only_2_LjSC, "LjSC KOs - Burkholderiaceae")
plot_LjSC_other <- anova_and_plot(KO_SSC_ex_only_2_LjSC, "LjSC KOs - Other families")

# Arrange plots using ggarrange
Burk_plot_2=grid.arrange(plot_AtSC_Burk,plot_AtSC_other,plot_LjSC_Burk, plot_LjSC_other, nrow = 2)

### Combine the the 3 with 2 host-specific or SynCom-specific KOs
combined_plot_2 <- 
  (plot_At_Caulo + plot_At_other_caulo + plot_Hv_Caulo +  plot_Hv_other_caulo) / 
  (plot_At_Xanth + plot_At_other_xanth + plot_Hv_Xanth + plot_Hv_other_xanth) /
  (plot_AtSC_Burk + plot_AtSC_other + plot_LjSC_Burk +  plot_LjSC_other) /
  plot_layout(ncol = 3)

print(combined_plot_2)

##Save plot
pdf(paste0(results.dir, "Figure_S16_boxplots_1.pdf"), width = 21, height = 8)
print(combined_plot_2)
dev.off()

#Just the Burkholderiaceae AtSC without dominators Plant plots
combined_plot_3 <- 
  (plot_At_Burk + plot_At_other) /
  (plot_Hv_Burk + plot_Hv_other) /
  (plot_Lj_Burk + plot_Lj_other) /
  plot_layout(ncol = 3)

print(combined_plot_3)

## Save plot
pdf(paste0(results.dir, "Figure_S16_boxplots_2.pdf"), width = 21, height = 4)
print(combined_plot_3)
dev.off()
