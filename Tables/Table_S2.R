library("dplyr") #Version 1.1.2
library("phyloseq") #Version 1.44.0
library("reshape2") #Version 1.4.4
library("tidyr") #Version 1.3.0

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Table S2 - abundance categories R2 simulation =====
norm_SSC=read.table(paste(working_directory,"Isolate_tables/Original/SSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
round_SSC=floor(x = norm_SSC)

#Taxonomy TABLE
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
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

data_wide_3 <- colnames(data_wide_2)[grep("Input",colnames(data_wide_2))]
data_wide_4 <- data_wide_2[,!colnames(data_wide_2) %in% data_wide_3]

#Count no of isolates
data_wide_4[data_wide_4 > 0] <- 1
Av_no_of_isolates <- round(sum(colSums(data_wide_4))/length(colnames(data_wide_4)))
Av_no_of_isolates

#Getting average per group
new <- data.frame(colSums(data_wide_4))
colnames(new) <- "number"
new$Plant <- samples_df$Condition[match(row.names(new), row.names(samples_df))]
new$SynCom <- samples_df$Inoculum[match(row.names(new), row.names(samples_df))]
new$Sample_group <- paste(new$Plant, new$SynCom, sep = "_")

hop_2 <- data.frame()

for (group in unique(new$Plant)){
  new_sub <- new[new$Plant == paste(group),]
  value <- sum(new_sub$number)/length(new_sub$number)
  hop <- t(data.frame(c(paste(group), value)))
  hop_2 <- rbind(hop_2, hop)
}

#Distribution
data_wide_5 <- as.numeric(unlist(data_wide_2[,!colnames(data_wide_2) %in% data_wide_3]))
data_wide_6 <- data_wide_5[data_wide_5 != 0]
value <- length(data_wide_6[data_wide_6 > 500])
data_wide_7 <- data_wide_6[data_wide_6 < 500]

hist_plot <- hist(data_wide_7, breaks = 20)
new_list <- as.data.frame(hist_plot$breaks)
new_list_2 <- as.data.frame(new_list[-1,])
new_list_2$counts <- hist_plot$counts

new_list_2 <- rbind(new_list_2, c(1000,value))
val <- sum(new_list_2$counts)
new_list_2$Rel <- new_list_2$counts/val
new_val <- 100/Av_no_of_isolates
new_list_2$new_counts <- (new_list_2$Rel/new_val)*100

new_list_3 <- as.data.frame(new_list_2[,1])
new_list_3$counts <- round(new_list_2$new_counts)
colnames(new_list_3) <- c("Category", "Counts")

write.table(new_list_3, paste(working_directory, "Table_S2_abundance_categories.txt", sep = ""), sep = "\t", quote = F)
