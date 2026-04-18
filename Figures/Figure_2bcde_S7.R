library("dplyr") #Version 1.1.2
library("phyloseq") #Version 1.44.0
library("reshape2") #Version 1.4.4
library("tidyr") #Version 1.3.0
library("ggplot2") #Version 3.4.2

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Figure 2bcde & S7 - Alpha Diversity plots =====
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
cond="ES"

samples_df_sub <- subset(samples_df, samples_df$Compartment != "RZ")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "AM")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "NOD")
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
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "NOD")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Condition != "NP")

samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")
samples <- sample_data(samples_df_sub_2)
phylo_sub = phyloseq(OTU,TAX, samples)

#Observed strains
method="Observed"

obs_df=estimate_richness(physeq = phylo_sub, measures = method)
obs_df=merge(x = obs_df, y = samples_df_sub_2, by = "row.names" )

obs_df=estimate_richness(physeq = phylo_sub, measures = method)
obs_df=merge(x = obs_df, y = samples_df_sub_2, by = "row.names" )
obs_df$inoculum_experiment=paste(obs_df$Condition,obs_df$Inoculum, sep ="_")

obs_df$Inoculum <- factor(obs_df$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC"))
obs_df$Condition <- factor(obs_df$Condition, levels = c("Input","At","Hv", "Lj"))
obs_df$Nutrient <- factor(obs_df$Nutrient, levels = c("high","low", "Input"))

colnames(obs_df)[colnames(obs_df) == "Condition"] <- "Host"

plants <- unique(obs_df$Host)
custom_labels_plant <- c(Input="Input", At = "Arabidopsis", Hv = "Barley", Lj = "Lotus")

# Initialize an empty list for storing results
anova_results <- list()
# Loop through each plant type to perform ANOVA, Tukey's test, and get letters
for(plant in plants) {
  # Subset data for the current plant
  subset_data <- obs_df[obs_df$Host == plant, ]
  
  # Perform ANOVA
  fitAnova <- aov(Observed ~ Inoculum, data=subset_data)
  
  # Perform Tukey's post-hoc test
  Tukey <- TukeyHSD(fitAnova)
  
  # Get letters
  letters_anova <- multcompView::multcompLetters4(fitAnova, Tukey)$Inoculum$Letters
  
  # Store results
  anova_results[[plant]] <- letters_anova
}
# Combine results into a data frame for plotting
ltlbl_combined <- do.call(rbind, lapply(names(anova_results), function(plant) {
  data.frame(Host = plant, Inoculum = names(anova_results[[plant]]), Letters = anova_results[[plant]])
}))
# Order factors based on your original setup
ltlbl_combined$Inoculum <- factor(ltlbl_combined$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC"))
ltlbl_combined$Host <- factor(ltlbl_combined$Host, levels = c("Input", "At", "Hv", "Lj"))

ltlbl_combined <- ltlbl_combined[order(ltlbl_combined$Host, ltlbl_combined$Inoculum), ]

alpha_tax_plot_obs_1 <- ggplot(data = obs_df, aes(x=Inoculum, y = Observed, color = Inoculum)) +
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+scale_color_manual(values = c("#A3A500","#00B0F6","#00BF7D","#F8766D", "black","gray70","black"))+
  scale_shape_manual(values = c(0,3))+
  geom_jitter(position=position_jitter(0.2), size =0.5, aes(color =Nutrient, shape=Experiment))+
  facet_wrap(~Host, scales="free_x", nrow=1, labeller=as_labeller(custom_labels_plant)) +
  theme(
    axis.text.x=element_blank(), 
    axis.title.x=element_blank(), 
    title=element_text(hjust=0.5, size=15), 
    axis.ticks.x=element_blank(),
    strip.background=element_rect(colour="gray50", size=0.3), # Change 'size' for thickness
    axis.text=element_text(color="gray50"),
    axis.line = element_line(color="gray50", size=0.3)
  ) +
  ylab(paste(method,"isolates"))+
  stat_summary(geom = 'text', label = ltlbl_combined$Letters, fun.y = max, aes(y = max(Observed)*1.05), show.legend=FALSE)

alpha_tax_plot_obs_1

#Plant subset
inocs <- unique(obs_df$Inoculum)

# Initialize an empty list for storing results
anova_results <- list()
# Loop through each plant type to perform ANOVA, Tukey's test, and get letters
for(inoc in inocs) {
  # Subset data for the current inoc
  subset_data <- obs_df[obs_df$Inoculum == inoc, ]
  
  # Perform ANOVA
  fitAnova <- aov(Observed ~ Host, data=subset_data)
  
  # Perform Tukey's post-hoc test
  Tukey <- TukeyHSD(fitAnova)
  
  # Get letters
  letters_anova <- multcompView::multcompLetters4(fitAnova, Tukey)$Host$Letters
  
  # Store results
  anova_results[[inoc]] <- letters_anova
}
# Combine results into a data frame for plotting
ltlbl_combined_2 <- do.call(rbind, lapply(names(anova_results), function(inoc) {
  data.frame(Inoculum = inoc, Host = names(anova_results[[inoc]]) , Letters = anova_results[[inoc]])
}))


# Order factors based on your original setup
ltlbl_combined_2$Inoculum <- factor(ltlbl_combined_2$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC"))
ltlbl_combined_2$Host <- factor(ltlbl_combined_2$Host, levels = c("Input", "At", "Hv", "Lj"))

#Alpha diversity plot, facet by plant/Input
alpha_tax_plot_obs_2 <- ggplot(data = obs_df, aes(x=Host, y = Observed, color = Host)) +
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+scale_color_manual(values = c("#1b9e77","#d95f02", "#e7298a","#7570b3","black", "gray70"))+
  scale_shape_manual(values = c(0,3))+
  geom_jitter(position=position_jitter(0.2), size =0.5, aes(color =Nutrient, shape=Experiment))+
  facet_wrap(~Inoculum, scales="free_x", nrow=1) +
  theme(
    axis.text.x=element_blank(), 
    axis.title.x=element_blank(), 
    title=element_text(hjust=0.5, size=15), 
    axis.ticks.x=element_blank(),
    strip.background=element_rect(colour="gray50", size=0.3), # Change 'size' for thickness
    axis.text=element_text(color="gray50"),
    axis.line = element_line(color="gray50", size=0.3)
  ) +
  ylab(paste(method,"isolates"))+
  stat_summary(geom = 'text', label = ltlbl_combined_2$Letters, fun.y = max, aes(y = max(Observed)*1.05), show.legend=FALSE)
alpha_tax_plot_obs_2

#Shannon diversity analysis
method="Shannon"

shannon_df=estimate_richness(physeq = phylo_sub, measures = method)
shannon_df=merge(x = shannon_df, y = samples_df_sub_2, by = "row.names" )

shannon_df=estimate_richness(physeq = phylo_sub, measures = method)
shannon_df=merge(x = shannon_df, y = samples_df_sub_2, by = "row.names" )
colnames(shannon_df)[colnames(shannon_df) == "Condition"] <- "Host"
shannon_df$inoculum_experiment=paste(shannon_df$Host,shannon_df$Inoculum, sep ="_")

shannon_df$Inoculum <- factor(shannon_df$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC"))
shannon_df$Host <- factor(shannon_df$Host, levels = c("Input","At","Hv", "Lj"))
shannon_df$Nutrient <- factor(shannon_df$Nutrient, levels = c("high","low", "Input"))

plants <- unique(shannon_df$Host)

# Initialize an empty list for storing results
anova_results <- list()
# Loop through each plant type to perform ANOVA, Tukey's test, and get letters
for(plant in plants) {
  # Subset data for the current plant
  subset_data <- shannon_df[shannon_df$Host == plant, ]
  
  # Perform ANOVA
  fitAnova <- aov(Shannon ~ Inoculum, data=subset_data)
  
  # Perform Tukey's post-hoc test
  Tukey <- TukeyHSD(fitAnova)
  
  # Get letters
  letters_anova <- multcompView::multcompLetters4(fitAnova, Tukey)$Inoculum$Letters
  
  # Store results
  anova_results[[plant]] <- letters_anova
}
# Combine results into a data frame for plotting
ltlbl_combined <- do.call(rbind, lapply(names(anova_results), function(plant) {
  data.frame(Host = plant, Inoculum = names(anova_results[[plant]]), Letters = anova_results[[plant]])
}))
# Order factors based on your original setup
ltlbl_combined$Inoculum <- factor(ltlbl_combined$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC"))
ltlbl_combined$Host <- factor(ltlbl_combined$Host, levels = c("Input", "At", "Hv", "Lj"))

ltlbl_combined <- ltlbl_combined[order(ltlbl_combined$Host, ltlbl_combined$Inoculum), ]

#Alpha diversity plot, facet by plant/Input
alpha_tax_plot_shan_1 <- ggplot(data = shannon_df, aes(x=Inoculum, y = Shannon, color = Inoculum)) +
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+scale_color_manual(values = c("#A3A500","#00B0F6","#00BF7D","#F8766D", "black","gray70","black"))+
  scale_shape_manual(values = c(0,3))+
  geom_jitter(position=position_jitter(0.2), size =0.5, aes(color =Nutrient, shape=Experiment))+
  facet_wrap(~Host, scales="free_x", nrow=1, labeller=as_labeller(custom_labels_plant)) +
  theme(
    axis.text.x=element_blank(), 
    axis.title.x=element_blank(), 
    title=element_text(hjust=0.5, size=15), 
    axis.ticks.x=element_blank(),
    strip.background=element_rect(colour="gray50", size=0.3), # Change 'size' for thickness
    axis.text=element_text(color="gray50"),
    axis.line = element_line(color="gray50", size=0.3)
  ) +
  ylab(paste(method,"isolates"))+
  stat_summary(geom = 'text', label = ltlbl_combined$Letters, fun.y = max, aes(y = max(Shannon)*1.05), show.legend=FALSE)

alpha_tax_plot_shan_1

#Alpha diversity plot, facet by Syncom

#Plant subset
inocs <- unique(shannon_df$Inoculum)

# Initialize an empty list for storing results
anova_results <- list()
# Loop through each plant type to perform ANOVA, Tukey's test, and get letters
for(inoc in inocs) {
  # Subset data for the current inoc
  subset_data <- shannon_df[shannon_df$Inoculum == inoc, ]
  
  # Perform ANOVA
  fitAnova <- aov(Shannon ~ Host, data=subset_data)
  
  # Perform Tukey's post-hoc test
  Tukey <- TukeyHSD(fitAnova)
  
  # Get letters
  letters_anova <- multcompView::multcompLetters4(fitAnova, Tukey)$Host$Letters
  
  # Store results
  anova_results[[inoc]] <- letters_anova
}
# Combine results into a data frame for plotting
ltlbl_combined_2 <- do.call(rbind, lapply(names(anova_results), function(inoc) {
  data.frame(Inoculum = inoc, Host = names(anova_results[[inoc]]) , Letters = anova_results[[inoc]])
}))

# Order factors based on your original setup
ltlbl_combined_2$Inoculum <- factor(ltlbl_combined_2$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC"))
ltlbl_combined_2$Host <- factor(ltlbl_combined_2$Host, levels = c("Input", "At", "Hv", "Lj"))

alpha_tax_plot_shan_2 <- ggplot(data = shannon_df, aes(x=Host, y = Shannon, color = Host)) +
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+scale_color_manual(values = c("#1b9e77","#d95f02", "#e7298a","#7570b3", "black", "gray70"))+
  scale_shape_manual(values = c(0,3))+
  geom_jitter(position=position_jitter(0.2), size =0.5, aes(color =Nutrient, shape=Experiment))+
  facet_wrap(~Inoculum, scales="free_x", nrow=1) +
  theme(
    axis.text.x=element_blank(), 
    axis.title.x=element_blank(), 
    title=element_text(hjust=0.5, size=15), 
    axis.ticks.x=element_blank(),
    strip.background=element_rect(colour="gray50", size=0.3), # Change 'size' for thickness
    axis.text=element_text(color="gray50"),
    axis.line = element_line(color="gray50", size=0.3)
  ) +
  ylab(paste(method,"isolates"))+
  stat_summary(geom = 'text', label = ltlbl_combined_2$Letters, fun.y = max, aes(y = max(Shannon)*1.05), show.legend=FALSE)
alpha_tax_plot_shan_2

#KO OTU TABLE
norm_SSC_KO =read.table(paste(working_directory,"KO_tables/Original/SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
round_SSC_KO=floor(x = norm_SSC_KO)

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
colnames(samples_df)[6]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparation
#Set the OTU, TAX and sample data for making phyloseq object
OTU = otu_table(as.matrix(round_SSC_KO),taxa_are_rows = TRUE)

#Sample subsetting
cond="ES"

samples_df_sub <- subset(samples_df, samples_df$Compartment != "RZ")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "AM")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "NOD")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Condition != "NP")

samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")
samples <- sample_data(samples_df_sub_2)

round_SSC_KO_2 <- round_SSC_KO[,colnames(round_SSC_KO) %in% row.names(samples_df_sub_2)]
round_SSC_2 <- round_SSC[,colnames(round_SSC) %in% row.names(samples_df_sub_2)]

average_mult_tax_to_KO <- sum(colSums(round_SSC_KO_2)/colSums(round_SSC_2))/length(colSums(round_SSC_KO_2))

phylo_KO_sub = phyloseq(OTU,samples)

subsetted_table <- otu_table(phylo_KO_sub)
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
subsetted_table_long_2$value[subsetted_table_long_2$Rel < 0.0005/average_mult_tax_to_KO] <- 0
subsetted_table_long_3 <- subsetted_table_long_2[1:3]
data_wide <- spread(subsetted_table_long_3, Var2, value)
row.names(data_wide) <- data_wide$Var1
data_wide_2 <- data_wide %>% dplyr::select (-Var1)

OTU = otu_table(as.matrix(data_wide_2),taxa_are_rows = TRUE)
samples_df_sub <- subset(samples_df, samples_df$Compartment != "RZ")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "AM")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "NOD")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Condition != "NP")

samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")
samples <- sample_data(samples_df_sub_2)
phylo_KO_sub = phyloseq(OTU, samples)

#Observed KOs
method="Observed"

obs_df_KO=estimate_richness(physeq = phylo_KO_sub, measures = method)
obs_df_KO=merge(x = obs_df_KO, y = samples_df_sub_2, by = "row.names" )

obs_df_KO=estimate_richness(physeq = phylo_KO_sub, measures = method)
obs_df_KO=merge(x = obs_df_KO, y = samples_df_sub_2, by = "row.names" )
obs_df_KO$inoculum_experiment=paste(obs_df_KO$Condition,obs_df_KO$Inoculum, sep ="_")

obs_df_KO$Inoculum <- factor(obs_df_KO$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC"))
obs_df_KO$Condition <- factor(obs_df_KO$Condition, levels = c("Input","At","Hv", "Lj"))
obs_df_KO$Nutrient <- factor(obs_df_KO$Nutrient, levels = c("high","low", "Input"))

colnames(obs_df_KO)[colnames(obs_df_KO) == "Condition"] <- "Host"

plants <- unique(obs_df_KO$Host)

# Initialize an empty list for storing results
anova_results <- list()
# Loop through each plant type to perform ANOVA, Tukey's test, and get letters
for(plant in plants) {
  # Subset data for the current plant
  subset_data <- obs_df_KO[obs_df_KO$Host == plant, ]
  
  # Perform ANOVA
  fitAnova <- aov(Observed ~ Inoculum, data=subset_data)
  
  # Perform Tukey's post-hoc test
  Tukey <- TukeyHSD(fitAnova)
  
  # Get letters
  letters_anova <- multcompView::multcompLetters4(fitAnova, Tukey)$Inoculum$Letters
  
  # Store results
  anova_results[[plant]] <- letters_anova
}
# Combine results into a data frame for plotting
ltlbl_combined <- do.call(rbind, lapply(names(anova_results), function(plant) {
  data.frame(Host = plant, Inoculum = names(anova_results[[plant]]), Letters = anova_results[[plant]])
}))
# Order factors based on your original setup
ltlbl_combined$Inoculum <- factor(ltlbl_combined$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC"))
ltlbl_combined$Host <- factor(ltlbl_combined$Host, levels = c("Input", "At", "Hv", "Lj"))

ltlbl_combined <- ltlbl_combined[order(ltlbl_combined$Host, ltlbl_combined$Inoculum), ]

#Alpha diversity plot, facet by plant/Input
alpha_plot_KO_obs_1 <- ggplot(data = obs_df_KO, aes(x=Inoculum, y = Observed, color = Inoculum)) +
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+scale_color_manual(values = c("#A3A500","#00B0F6","#00BF7D","#F8766D", "black","gray70","black"))+
  scale_shape_manual(values = c(0,3))+
  geom_jitter(position=position_jitter(0.2), size =0.5, aes(color =Nutrient, shape=Experiment))+
  facet_wrap(~Host, scales="free_x", nrow=1, labeller=as_labeller(custom_labels_plant)) +
  theme(
    axis.text.x=element_blank(), 
    axis.title.x=element_blank(), 
    title=element_text(hjust=0.5, size=15), 
    axis.ticks.x=element_blank(),
    strip.background=element_rect(colour="gray50", size=0.3), # Change 'size' for thickness
    axis.text=element_text(color="gray50"),
    axis.line = element_line(color="gray50", size=0.3)
  ) +
  ylab(paste(method,"KOs"))+
  stat_summary(geom = 'text', label = ltlbl_combined$Letters, fun.y = max, aes(y = max(Observed)*1.05), show.legend=FALSE)

alpha_plot_KO_obs_1

#Alpha diversity plot, facet by Syncom
#Plant subset
inocs <- unique(obs_df_KO$Inoculum)

# Initialize an empty list for storing results
anova_results <- list()
# Loop through each plant type to perform ANOVA, Tukey's test, and get letters
for(inoc in inocs) {
  # Subset data for the current inoc
  subset_data <- obs_df_KO[obs_df_KO$Inoculum == inoc, ]
  
  # Perform ANOVA
  fitAnova <- aov(Observed ~ Host, data=subset_data)
  
  # Perform Tukey's post-hoc test
  Tukey <- TukeyHSD(fitAnova)
  
  # Get letters
  letters_anova <- multcompView::multcompLetters4(fitAnova, Tukey)$Host$Letters
  
  # Store results
  anova_results[[inoc]] <- letters_anova
}
# Combine results into a data frame for plotting
ltlbl_combined_2 <- do.call(rbind, lapply(names(anova_results), function(inoc) {
  data.frame(Inoculum = inoc, Host = names(anova_results[[inoc]]) , Letters = anova_results[[inoc]])
}))

# Order factors based on your original setup
ltlbl_combined_2$Inoculum <- factor(ltlbl_combined_2$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC"))
ltlbl_combined_2$Host <- factor(ltlbl_combined_2$Host, levels = c("Input", "At", "Hv", "Lj"))

alpha_plot_KO_obs_2 <- ggplot(data = obs_df_KO, aes(x=Host, y = Observed, color = Host)) +
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+scale_color_manual(values = c("#1b9e77", "#d95f02", "#e7298a","#7570b3","black", "gray70"))+
  scale_shape_manual(values = c(0,3))+
  geom_jitter(position=position_jitter(0.2), size =0.5, aes(color =Nutrient, shape=Experiment))+
  facet_wrap(~Inoculum, scales="free_x", nrow=1) +
  theme(
    axis.text.x=element_blank(), 
    axis.title.x=element_blank(), 
    title=element_text(hjust=0.5, size=15), 
    axis.ticks.x=element_blank(),
    strip.background=element_rect(colour="gray50", size=0.3), # Change 'size' for thickness
    axis.text=element_text(color="gray50"),
    axis.line = element_line(color="gray50", size=0.3)
  ) +
  ylab(paste(method,"KOs"))+
  stat_summary(geom = 'text', label = ltlbl_combined_2$Letters, fun.y = max, aes(y = max(Observed)*1.05), show.legend=FALSE)
alpha_plot_KO_obs_2

#Shannon diversity KOs
method="Shannon"

shannon_df_KO=estimate_richness(physeq = phylo_KO_sub, measures = method)
shannon_df_KO=merge(x = shannon_df_KO, y = samples_df_sub_2, by = "row.names" )
shannon_df_KO$inoculum_experiment=paste(shannon_df_KO$Condition,shannon_df_KO$Inoculum, sep ="_")

shannon_df_KO$Inoculum <- factor(shannon_df_KO$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC"))
shannon_df_KO$Condition <- factor(shannon_df_KO$Condition, levels = c("Input","At","Hv", "Lj"))
shannon_df_KO$Nutrient <- factor(shannon_df_KO$Nutrient, levels = c("high","low", "Input"))

colnames(shannon_df_KO)[colnames(shannon_df_KO) == "Condition"] <- "Host"

plants <- unique(shannon_df_KO$Host)

# Initialize an empty list for storing results
anova_results <- list()
# Loop through each plant type to perform ANOVA, Tukey's test, and get letters
for(plant in plants) {
  # Subset data for the current plant
  subset_data <- shannon_df_KO[shannon_df_KO$Host == plant, ]
  
  # Perform ANOVA
  fitAnova <- aov(Shannon ~ Inoculum, data=subset_data)
  
  # Perform Tukey's post-hoc test
  Tukey <- TukeyHSD(fitAnova)
  
  # Get letters
  letters_anova <- multcompView::multcompLetters4(fitAnova, Tukey)$Inoculum$Letters
  
  # Store results
  anova_results[[plant]] <- letters_anova
}
# Combine results into a data frame for plotting
ltlbl_combined <- do.call(rbind, lapply(names(anova_results), function(plant) {
  data.frame(Host = plant, Inoculum = names(anova_results[[plant]]), Letters = anova_results[[plant]])
}))
# Order factors based on your original setup
ltlbl_combined$Inoculum <- factor(ltlbl_combined$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC"))
ltlbl_combined$Host <- factor(ltlbl_combined$Host, levels = c("Input", "At", "Hv", "Lj"))

ltlbl_combined <- ltlbl_combined[order(ltlbl_combined$Host, ltlbl_combined$Inoculum), ]

#Alpha diversity plot, facet by plant/Input
alpha_plot_KO_shan_1 <- ggplot(data = shannon_df_KO, aes(x=Inoculum, y = Shannon, color = Inoculum)) +
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+scale_color_manual(values = c("#A3A500","#00B0F6","#00BF7D","#F8766D", "black","gray70","black"))+
  scale_shape_manual(values = c(0,3))+
  geom_jitter(position=position_jitter(0.2), size =0.5, aes(color =Nutrient, shape=Experiment))+
  facet_wrap(~Host, scales="free_x", nrow=1, labeller=as_labeller(custom_labels_plant)) +
  theme(
    axis.text.x=element_blank(), 
    axis.title.x=element_blank(), 
    title=element_text(hjust=0.5, size=15), 
    axis.ticks.x=element_blank(),
    strip.background=element_rect(colour="gray50", size=0.3), # Change 'size' for thickness
    axis.text=element_text(color="gray50"),
    axis.line = element_line(color="gray50", size=0.3)
  ) +
  ylab(paste(method,"KOs"))+
  stat_summary(geom = 'text', label = ltlbl_combined$Letters, fun.y = max, aes(y = max(Shannon)*1.05), show.legend=FALSE)

alpha_plot_KO_shan_1

#Alpha diversity plot, facet by Syncom
#Plant subset
inocs <- unique(shannon_df_KO$Inoculum)

# Initialize an empty list for storing results
anova_results <- list()
# Loop through each plant type to perform ANOVA, Tukey's test, and get letters
for(inoc in inocs) {
  # Subset data for the current inoc
  subset_data <- shannon_df_KO[shannon_df_KO$Inoculum == inoc, ]
  
  # Perform ANOVA
  fitAnova <- aov(Shannon ~ Host, data=subset_data)
  
  # Perform Tukey's post-hoc test
  Tukey <- TukeyHSD(fitAnova)
  
  # Get letters
  letters_anova <- multcompView::multcompLetters4(fitAnova, Tukey)$Host$Letters
  
  # Store results
  anova_results[[inoc]] <- letters_anova
}
# Combine results into a data frame for plotting
ltlbl_combined_2 <- do.call(rbind, lapply(names(anova_results), function(inoc) {
  data.frame(Inoculum = inoc, Host = names(anova_results[[inoc]]) , Letters = anova_results[[inoc]])
}))


# Order factors based on your original setup
ltlbl_combined_2$Inoculum <- factor(ltlbl_combined_2$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC"))
ltlbl_combined_2$Host <- factor(ltlbl_combined_2$Host, levels = c("Input", "At", "Hv", "Lj"))

alpha_plot_KO_shan_2 <- ggplot(data = shannon_df_KO, aes(x=Host, y = Shannon, color = Host)) +
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+scale_color_manual(values = c("#1b9e77", "#d95f02", "#e7298a","#7570b3","black", "gray70"))+
  scale_shape_manual(values = c(0,3))+
  geom_jitter(position=position_jitter(0.2), size =0.5, aes(color =Nutrient, shape=Experiment))+
  facet_wrap(~Inoculum, scales="free_x", nrow=1) +
  theme(
    axis.text.x=element_blank(), 
    axis.title.x=element_blank(), 
    title=element_text(hjust=0.5, size=15), 
    axis.ticks.x=element_blank(),
    strip.background=element_rect(colour="gray50", size=0.3), # Change 'size' for thickness
    axis.text=element_text(color="gray50"),
    axis.line = element_line(color="gray50", size=0.3)
  ) +
  ylab(paste(method,"KOs"))+
  stat_summary(geom = 'text', label = ltlbl_combined_2$Letters, fun.y = max, aes(y = max(Shannon)*1.05), show.legend=FALSE)
alpha_plot_KO_shan_2

#Printing all alpha diversity plots
pdf(paste(results.dir,"alpha_diversity_2b.pdf", sep=""), width=7, height=3)
print(alpha_tax_plot_obs_2)
dev.off()

pdf(paste(results.dir,"alpha_diversity_2c.pdf", sep=""), width=7, height=3)
print(alpha_tax_plot_obs_1)
dev.off()

pdf(paste(results.dir,"alpha_diversity_2d.pdf", sep=""), width=7, height=3)
print(alpha_plot_KO_obs_2)
dev.off()

pdf(paste(results.dir,"alpha_diversity_2e.pdf", sep=""), width=7, height=3)
print(alpha_plot_KO_obs_1)
dev.off()

pdf(paste(results.dir,"alpha_diversity_S7a.pdf", sep=""), width=7, height=3)
print(alpha_tax_plot_shan_2)
dev.off()

pdf(paste(results.dir,"alpha_diversity_S7b.pdf", sep=""), width=7, height=3)
print(alpha_tax_plot_shan_1)
dev.off()

pdf(paste(results.dir,"alpha_diversity_S7c.pdf", sep=""), width=7, height=3)
print(alpha_plot_KO_shan_2)
dev.off()

pdf(paste(results.dir,"alpha_diversity_S7d.pdf", sep=""), width=7, height=3)
print(alpha_plot_KO_shan_1)
dev.off()
