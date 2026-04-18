library("dplyr") #Version 1.1.2
library("phyloseq") #Version 1.44.0
library("reshape2") #Version 1.4.4
library("tidyr") #Version 1.3.0
library("ggplot2") #Version 3.4.2
library("ggpubr") #Version 0.6.0

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Figure S10 - Alpha diversity - nutrient condition subset =====
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
samples_df_sub_2$Condition[samples_df_sub_2$Condition == "At"] <- "Arabidopsis"
samples_df_sub_2$Condition[samples_df_sub_2$Condition == "Hv"] <- "Barley"
samples_df_sub_2$Condition[samples_df_sub_2$Condition == "Lj"] <- "Lotus"

samples <- sample_data(samples_df_sub_2)
phylo_sub = phyloseq(OTU,TAX, samples)

#Observed isolates
method="Observed"

obs_df=estimate_richness(physeq = phylo_sub, measures = method)
obs_df=merge(x = obs_df, y = samples_df_sub_2, by = "row.names" )

obs_df=estimate_richness(physeq = phylo_sub, measures = method)
obs_df=merge(x = obs_df, y = samples_df_sub_2, by = "row.names" )
obs_df$inoculum_experiment=paste(obs_df$Condition,obs_df$Inoculum, sep =" ")

obs_df$Inoculum <- factor(obs_df$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC"))
obs_df$Condition <- factor(obs_df$Condition, levels = c("Input","Arabidopsis","Barley", "Lotus"))
obs_df$Nutrient <- factor(obs_df$Nutrient, levels = c("high","low", "Input"))

colnames(obs_df)[colnames(obs_df) == "Condition"] <- "Host"

obs_df_2 <- obs_df[obs_df$Host != "Input",]

obs_nut_sub <- ggplot(obs_df_2, aes(x = Nutrient, y = Observed, color=Nutrient)) +
  theme_classic()+
  scale_color_manual(values = c("black", "gray70"))+
  scale_shape_manual(values = c(0,3))+
  scale_fill_manual(values = c("white","gray70"))+
  geom_jitter()+
  geom_boxplot(outlier.shape = NA) + # Hide outliers since jitter will show all points
  facet_wrap(~inoculum_experiment, scales="free_x", nrow=2) +
  theme( axis.text.x=element_blank(), 
         axis.title.x=element_blank(), 
         title=element_text(hjust=0.5, size=15), 
         axis.ticks.x=element_blank(),
         strip.background=element_rect(colour="gray50", size=0.3), # Change 'size' for thickness
         axis.text=element_text(color="gray50"),
         axis.line = element_line(color="gray50", size=0.3)) +
  labs(title = "",
       x = "",
       y = "Observed isolates") +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), label = "p.signif", vjust = 0.7, label.x = 1.5)
obs_nut_sub

#Shannon diversity isolates
method="Shannon"

shannon_df=estimate_richness(physeq = phylo_sub, measures = method)
shannon_df=merge(x = shannon_df, y = samples_df_sub_2, by = "row.names" )

shannon_df=estimate_richness(physeq = phylo_sub, measures = method)
shannon_df=merge(x = shannon_df, y = samples_df_sub_2, by = "row.names" )
shannon_df$inoculum_experiment=paste(shannon_df$Condition,shannon_df$Inoculum, sep =" ")

shannon_df$Inoculum <- factor(shannon_df$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC"))
shannon_df$Condition <- factor(shannon_df$Condition, levels = c("Input","Arabidopsis","Barley", "Lotus"))
shannon_df$Nutrient <- factor(shannon_df$Nutrient, levels = c("high","low", "Input"))

colnames(shannon_df)[colnames(shannon_df) == "Condition"] <- "Host"

shannon_df_2 <- shannon_df[shannon_df$Host != "Input",]

shan_nut_sub <- ggplot(shannon_df_2, aes(x = Nutrient, y = Shannon, color=Nutrient)) +
  theme_classic()+
  scale_color_manual(values = c("black", "gray70"))+
  scale_shape_manual(values = c(0,3))+
  scale_fill_manual(values = c("white","gray70"))+
  geom_jitter()+
  geom_boxplot(outlier.shape = NA) + # Hide outliers since jitter will show all points
  facet_wrap(~inoculum_experiment, scales="free_x", nrow=2) +
  theme( axis.text.x=element_blank(), 
         axis.title.x=element_blank(), 
         title=element_text(hjust=0.5, size=15), 
         axis.ticks.x=element_blank(),
         strip.background=element_rect(colour="gray50", size=0.3), # Change 'size' for thickness
         axis.text=element_text(color="gray50"),
         axis.line = element_line(color="gray50", size=0.3)) +
  labs(title = "",
       x = "",
       y = "Shannon isolates") +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), label = "p.signif", vjust = 0.7, label.x = 1.5)
shan_nut_sub

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
samples_df_sub_2$Condition[samples_df_sub_2$Condition == "At"] <- "Arabidopsis"
samples_df_sub_2$Condition[samples_df_sub_2$Condition == "Hv"] <- "Barley"
samples_df_sub_2$Condition[samples_df_sub_2$Condition == "Lj"] <- "Lotus"

samples <- sample_data(samples_df_sub_2)
phylo_KO_sub = phyloseq(OTU, samples)

#Observed KOs
method="Observed"

obs_df_KO=estimate_richness(physeq = phylo_KO_sub, measures = method)
obs_df_KO=merge(x = obs_df_KO, y = samples_df_sub_2, by = "row.names" )

obs_df_KO=estimate_richness(physeq = phylo_KO_sub, measures = method)
obs_df_KO=merge(x = obs_df_KO, y = samples_df_sub_2, by = "row.names" )
obs_df_KO$inoculum_experiment=paste(obs_df_KO$Condition,obs_df_KO$Inoculum, sep =" ")

obs_df_KO$Inoculum <- factor(obs_df_KO$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC"))
obs_df_KO$Condition <- factor(obs_df_KO$Condition, levels = c("Input","Arabidopsis","Barley", "Lotus"))
obs_df_KO$Nutrient <- factor(obs_df_KO$Nutrient, levels = c("high","low", "Input"))

colnames(obs_df_KO)[colnames(obs_df_KO) == "Condition"] <- "Host"

obs_df_KO_2 <- obs_df_KO[obs_df_KO$Host != "Input",]

obs_nut_sub_KO <- ggplot(obs_df_KO_2, aes(x = Nutrient, y = Observed, color=Nutrient)) +
  theme_classic()+
  scale_color_manual(values = c("black", "gray70"))+
  scale_shape_manual(values = c(0,3))+
  scale_fill_manual(values = c("white","gray70"))+
  geom_jitter()+
  geom_boxplot(outlier.shape = NA) + # Hide outliers since jitter will show all points
  facet_wrap(~inoculum_experiment, scales="free_x", nrow=2) +
  theme( axis.text.x=element_blank(), 
         axis.title.x=element_blank(), 
         title=element_text(hjust=0.5, size=15), 
         axis.ticks.x=element_blank(),
         strip.background=element_rect(colour="gray50", size=0.3), # Change 'size' for thickness
         axis.text=element_text(color="gray50"),
         axis.line = element_line(color="gray50", size=0.3)) +
  labs(title = "",
       x = "",
       y = "Observed KOs") +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), label = "p.signif", vjust = 0.7, label.x = 1.5)
obs_nut_sub_KO

#Shannon diversity KOs
method="Shannon"

shannon_KO_df=estimate_richness(physeq = phylo_sub, measures = method)
shannon_KO_df=merge(x = shannon_KO_df, y = samples_df_sub_2, by = "row.names" )

shannon_KO_df=estimate_richness(physeq = phylo_sub, measures = method)
shannon_KO_df=merge(x = shannon_KO_df, y = samples_df_sub_2, by = "row.names" )
shannon_KO_df$inoculum_experiment=paste(shannon_KO_df$Condition,shannon_KO_df$Inoculum, sep =" ")

shannon_KO_df$Inoculum <- factor(shannon_KO_df$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC"))
shannon_KO_df$Condition <- factor(shannon_KO_df$Condition, levels = c("Input","Arabidopsis","Barley", "Lotus"))
shannon_KO_df$Nutrient <- factor(shannon_KO_df$Nutrient, levels = c("high","low", "Input"))

colnames(shannon_KO_df)[colnames(shannon_KO_df) == "Condition"] <- "Host"

shannon_KO_df_2 <- shannon_KO_df[shannon_KO_df$Host != "Input",]

shan_nut_sub_KO <- ggplot(shannon_KO_df_2, aes(x = Nutrient, y = Shannon, color=Nutrient)) +
  theme_classic()+
  scale_color_manual(values = c("black", "gray70"))+
  scale_shape_manual(values = c(0,3))+
  scale_fill_manual(values = c("white","gray70"))+
  geom_jitter()+
  geom_boxplot(outlier.shape = NA) + # Hide outliers since jitter will show all points
  facet_wrap(~inoculum_experiment, scales="free_x", nrow=2) +
  theme( axis.text.x=element_blank(), 
         axis.title.x=element_blank(), 
         title=element_text(hjust=0.5, size=15), 
         axis.ticks.x=element_blank(),
         strip.background=element_rect(colour="gray50", size=0.3), # Change 'size' for thickness
         axis.text=element_text(color="gray50"),
         axis.line = element_line(color="gray50", size=0.3)) +
  labs(title = "",
       x = "",
       y = "Shannon KOs") +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), label = "p.signif", vjust = 0.7, label.x = 1.5)
shan_nut_sub_KO

alpha_shan <- ggarrange(obs_nut_sub, obs_nut_sub_KO, shan_nut_sub,shan_nut_sub_KO,ncol = 1, nrow = 4, labels = c("a", "b", "c", "d"), common.legend = T)

pdf(paste(results.dir,"Figure_S10_alpha_diversity_nutrients.pdf", sep=""), width=8, height=18)
print(alpha_shan)
dev.off()
