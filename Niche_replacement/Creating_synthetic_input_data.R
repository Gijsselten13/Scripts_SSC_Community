library("phyloseq") #Version 1.44.0
library("vegan") #Version 2.6-4

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Creating Data_with_synthetic_input.tsv =====

#Since only one input samples was sequenced for each SynCom, we will subdivide these into four samples by random subsampling at a lower sequencing depth
otu_table <- read.table(paste(working_directory,"LjSC_Family_drop_out_experiment/original_data_norm.tsv", sep = ""), sep = "\t", row.names =1, header =T)

Input_data <- otu_table[,grep("INPUT", colnames(otu_table))]
OTU = otu_table(as.matrix(Input_data),taxa_are_rows = TRUE)
phylo_sub = phyloseq(OTU)

racur_data=rarecurve(x =floor(t(phylo_sub@.Data)), step = 250, label = F, tidy = T)
racur_data$Site <- gsub("LDT1_INPUT", "Burkholderiaceae drop out", racur_data$Site)
racur_data$Site <- gsub("LDT2_INPUT", "Caulobacteraceae drop out", racur_data$Site)
racur_data$Site <- gsub("LDT3_INPUT", "Pseudomonadaceae drop out", racur_data$Site)
racur_data$Site <- gsub("LDT4_INPUT", "Rhizobiaceae drop out", racur_data$Site)
racur_data$Site <- gsub("LDT5_INPUT", "All other families drop out", racur_data$Site)
racur_data$Site <- gsub("LDT6_INPUT", "Full LjSC", racur_data$Site)

values <- c()

for (group in unique(racur_data$Site)){
  racur_sub <- racur_data[racur_data$Site == paste(group),]
  number <- tail(racur_sub$Species, n=1) * 0.95
  depth <- racur_sub$Sample[racur_sub$Species > number][1]
  
  values <- c(values, depth)
}

#Selected rarefaction depth
rarefaction_depth <- sum(values)/length(values)

#creating subsets
otu_table_no_input <- otu_table[,!grepl("INPUT", colnames(otu_table))]
Input_data <- otu_table[,grep("INPUT", colnames(otu_table))]

OTU = otu_table(as.matrix(Input_data),taxa_are_rows = TRUE)
phylo_sub = phyloseq(OTU)
phylo_sub_2 <- rarefy_even_depth(phylo_sub, sample.size = rarefaction_depth, rngseed = 1)
phylo_sub_3 <- rarefy_even_depth(phylo_sub, sample.size = rarefaction_depth, rngseed = 2)
phylo_sub_4 <- rarefy_even_depth(phylo_sub, sample.size = rarefaction_depth, rngseed = 3)
phylo_sub_5 <- rarefy_even_depth(phylo_sub, sample.size = rarefaction_depth, rngseed = 4)

otu_2 <- phylo_sub_2@.Data
otu_3 <- phylo_sub_3@.Data
otu_4 <- phylo_sub_4@.Data
otu_5 <- phylo_sub_5@.Data

colnames(otu_2) <- gsub("INPUT","INPUT_B1", colnames(otu_2))
colnames(otu_3) <- gsub("INPUT","INPUT_B2", colnames(otu_3))
colnames(otu_4) <- gsub("INPUT","INPUT_B3", colnames(otu_4))
colnames(otu_5) <- gsub("INPUT","INPUT_B4", colnames(otu_5))

reset_names_2 <- row.names(otu_table_no_input)[!row.names(otu_table_no_input) %in% row.names(otu_2) ]
for (group in reset_names_2){
  new <- data.frame(t(data.frame(c(0,0,0,0,0,0))))
  row.names(new) <- paste(group)
  colnames(new) <- colnames(otu_2)
  otu_2 <- rbind(otu_2, new)
}

reset_names_3 <- row.names(otu_table_no_input)[!row.names(otu_table_no_input) %in% row.names(otu_3) ]
for (group in reset_names_3){
  new <- data.frame(t(data.frame(c(0,0,0,0,0,0))))
  row.names(new) <- paste(group)
  colnames(new) <- colnames(otu_3)
  otu_3 <- rbind(otu_3, new)
}

reset_names_4 <- row.names(otu_table_no_input)[!row.names(otu_table_no_input) %in% row.names(otu_4) ]
for (group in reset_names_4){
  new <- data.frame(t(data.frame(c(0,0,0,0,0,0))))
  row.names(new) <- paste(group)
  colnames(new) <- colnames(otu_4)
  otu_4 <- rbind(otu_4, new)
}

reset_names_5 <- row.names(otu_table_no_input)[!row.names(otu_table_no_input) %in% row.names(otu_5) ]
for (group in reset_names_5){
  new <- data.frame(t(data.frame(c(0,0,0,0,0,0))))
  row.names(new) <- paste(group)
  colnames(new) <- colnames(otu_5)
  otu_5 <- rbind(otu_5, new)
}

otu_2_sub <- otu_2[row.names(otu_table_no_input),]
otu_3_sub <- otu_3[row.names(otu_table_no_input),]
otu_4_sub <- otu_4[row.names(otu_table_no_input),]
otu_5_sub <- otu_5[row.names(otu_table_no_input),]

new_table <- cbind(otu_table_no_input,otu_2, otu_3, otu_4, otu_5)

write.table(new_table,paste(working_directory, "LjSC_Family_drop_out_experiment/Data_with_synthetic_input.tsv", sep =""), sep = "\t", quote =F, col.names = T, row.names =T)
