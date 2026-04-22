library("reshape2") #Version 1.4.4
library("tidyr") #Version 1.3.0

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Table S8 - KO overview =====
KO_table = read.table(paste(working_directory, "KO_genome/KO_LjSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
sigtab_col_all_2 <- read.table(paste(working_directory, "LjSC_Family_drop_out_experiment/DESeq2_Root_vs_input_Fam_drop_no_nod.txt", sep = ""), header =T,  sep = "\t")
sigtab_col_all_2_with_nod <- read.table(paste(working_directory, "LjSC_Family_drop_out_experiment/DESeq2_Root_vs_input_Fam_drop.txt", sep = ""), header =T,  sep = "\t")

taxonomy <- read.table("/media/wolf1/InRoot/Gijs_Selten/SSC_R2_11_2022/SSC_Community_2024/SSC_taxonomy_GTDB.tsv", header =T, row.names =1)
taxonomy_LjSC <- taxonomy[taxonomy$SynCom == "LjSC",]

families <- c("Pseudomonadaceae", "Rhizobiaceae", "Burkholderiaceae", "Caulobacteraceae", "other")
Nodulators <- c("LjRoot228", "LjNodule214", "LjRoot234")

input_table <- read.table(paste(working_directory,"DESeq2/Sig_KO_all_no_nod_rhizo.txt", sep = ""), header=T, sep="\t")
input_table_2 <- table(input_table$KO)
input_table_3 <- names(input_table_2)[input_table_2 == 12]

#No nodulators dataset
new_fam_data <- data.frame()

for (group in families){
  if(group == "other"){
    new_tax <- row.names(taxonomy_LjSC)[!taxonomy_LjSC$family %in% families]
    new_tax <- new_tax[!new_tax %in% Nodulators]
    sigtab_col_all_2_sub <- sigtab_col_all_2[sigtab_col_all_2$Subset == "All other families drop out",]
  } else {
    new_tax <- row.names(taxonomy_LjSC)[taxonomy_LjSC$family == paste(group)]
    new_tax <- new_tax[!new_tax %in% Nodulators]
    sigtab_col_all_2_sub <- sigtab_col_all_2[sigtab_col_all_2$Subset == paste(group, "drop out", sep = " "),]
  }
  
  KO_table_sub <- KO_table[,colnames(KO_table) %in% new_tax]
  KO_table_sub_2 <- rowSums(KO_table_sub)
  total_KOs <- length(names(KO_table_sub_2)[KO_table_sub_2 != 0])
  Family_KOs <- names(KO_table_sub_2)[KO_table_sub_2 != 0]
  
  sigtab_col_all_2_full <- sigtab_col_all_2[sigtab_col_all_2_with_nod$Subset == "Full LjSC",]
  sigtab_col_all_2_full_2 <- sigtab_col_all_2_full[sigtab_col_all_2_full$padj < 0.05,]
  sigtab_col_all_2_full_3 <- sigtab_col_all_2_full_2[sigtab_col_all_2_full_2$log2FoldChange > 0,]
  contrib_full <- length(sigtab_col_all_2_full_3$KO[sigtab_col_all_2_full_3$KO %in% Family_KOs])
  
  sigtab_col_all_2_2 <- sigtab_col_all_2_sub[sigtab_col_all_2_sub$padj < 0.05,]
  sigtab_col_all_2_3 <- sigtab_col_all_2_2[sigtab_col_all_2_2$log2FoldChange > 0,]
  lost_KOs <- length(sigtab_col_all_2_full_3$KO[!sigtab_col_all_2_full_3$KO %in% sigtab_col_all_2_3$KO])
  
  total_266_Fam <- length(Family_KOs[Family_KOs %in% input_table_3])
  sigtab_col_all_2_full_3_sig <- sigtab_col_all_2_full_3$KO[!sigtab_col_all_2_full_3$KO %in% sigtab_col_all_2_3$KO]
  lost_266 <- length(sigtab_col_all_2_full_3_sig[sigtab_col_all_2_full_3_sig %in% input_table_3])
  
  all_data <- data.frame(t(data.frame(c(paste(group), total_KOs, contrib_full, lost_KOs, total_266_Fam, lost_266 ))))
  new_fam_data <- rbind(new_fam_data, all_data)
}

row.names(new_fam_data) <- NULL
colnames(new_fam_data) <- c("Family", "Total no of KOs", "Total sig KOs", "Lost KOs in drop out", "Total no of 266 KOs", "Lost 266 KOs in drop out")

new_fam_data$`Total no of KOs` <- as.numeric(new_fam_data$`Total no of KOs`)
new_fam_data$`Total sig KOs` <- as.numeric(new_fam_data$`Total sig KOs`)
new_fam_data$`Lost KOs in drop out` <- as.numeric(new_fam_data$`Lost KOs in drop out`)
new_fam_data$`Total no of 266 KOs` <- as.numeric(new_fam_data$`Total no of 266 KOs`)
new_fam_data$`Lost 266 KOs in drop out` <- as.numeric(new_fam_data$`Lost 266 KOs in drop out`)

new_fam_data_3 <- melt(new_fam_data)

for (group in unique(new_fam_data_3$Family)){
  new_fam_data_3_sub <- new_fam_data_3[new_fam_data_3$Family == paste(group),]
  prop <- as.numeric(new_fam_data_3_sub$value[new_fam_data_3_sub$variable == "Lost KOs in drop out"])/as.numeric(new_fam_data_3_sub$value[new_fam_data_3_sub$variable == "Total sig KOs"])
  have_266 <- as.numeric(new_fam_data_3_sub$value[new_fam_data_3_sub$variable == "Lost 266 KOs in drop out"])
  not_have_266 <- as.numeric(new_fam_data_3_sub$value[new_fam_data_3_sub$variable == "Total no of 266 KOs"])
  not_have_266_2 <- (not_have_266 - have_266)
  
  binom_out <- binom.test(have_266, not_have_266, prop)
  pval <- binom_out$p.value
  
  new_data <- data.frame(t(data.frame(c(paste(group), "Percentage lost KOs vs sig", paste(prop)))))
  colnames(new_data) <- colnames(new_fam_data_3)
  row.names(new_data) <- NULL
  
  new_fam_data_3 <- rbind(new_fam_data_3,new_data)
  
  new_data <- data.frame(t(data.frame(c(paste(group), "Percentage lost KOs vs sig - 266", have_266/not_have_266))))
  colnames(new_data) <- colnames(new_fam_data_3)
  row.names(new_data) <- NULL
  
  new_fam_data_3 <- rbind(new_fam_data_3,new_data)
  
  new_data <- data.frame(t(data.frame(c(paste(group), "Binomial p-value", paste(round(pval,5))))))
  colnames(new_data) <- colnames(new_fam_data_3)
  row.names(new_data) <- NULL
  
  new_fam_data_3 <- rbind(new_fam_data_3,new_data)
}

data_wide_4 <- spread(new_fam_data_3, variable, value)
#Rhizobiaceae were eventually excluded as there was contamination in the drop out
data_wide_4_2 <- data_wide_4[data_wide_4$Family != "Rhizobiaceae",]
data_wide_4_2$Dataset <- "No nodulator"
#With dominator dataset
new_fam_data <- data.frame()

for (group in families){
  if(group == "other"){
    new_tax <- row.names(taxonomy_LjSC)[!taxonomy_LjSC$family %in% families]
    sigtab_col_all_2_with_nod_sub <- sigtab_col_all_2_with_nod[sigtab_col_all_2_with_nod$Subset == "All other families drop out",]
  } else {
    new_tax <- row.names(taxonomy_LjSC)[taxonomy_LjSC$family == paste(group)]
    sigtab_col_all_2_with_nod_sub <- sigtab_col_all_2_with_nod[sigtab_col_all_2_with_nod$Subset == paste(group, "drop out", sep = " "),]
  }
  
  KO_table_sub <- KO_table[,colnames(KO_table) %in% new_tax]
  KO_table_sub_2 <- rowSums(KO_table_sub)
  total_KOs <- length(names(KO_table_sub_2)[KO_table_sub_2 != 0])
  Family_KOs <- names(KO_table_sub_2)[KO_table_sub_2 != 0]
  
  sigtab_col_all_2_with_nod_full <- sigtab_col_all_2_with_nod[sigtab_col_all_2_with_nod$Subset == "Full LjSC",]
  sigtab_col_all_2_with_nod_full_2 <- sigtab_col_all_2_with_nod_full[sigtab_col_all_2_with_nod_full$padj < 0.05,]
  sigtab_col_all_2_with_nod_full_3 <- sigtab_col_all_2_with_nod_full_2[sigtab_col_all_2_with_nod_full_2$log2FoldChange > 0,]
  contrib_full <- length(sigtab_col_all_2_with_nod_full_3$KO[sigtab_col_all_2_with_nod_full_3$KO %in% Family_KOs])
  
  sigtab_col_all_2_with_nod_2 <- sigtab_col_all_2_with_nod_sub[sigtab_col_all_2_with_nod_sub$padj < 0.05,]
  sigtab_col_all_2_with_nod_3 <- sigtab_col_all_2_with_nod_2[sigtab_col_all_2_with_nod_2$log2FoldChange > 0,]
  lost_KOs <- length(sigtab_col_all_2_with_nod_full_3$KO[!sigtab_col_all_2_with_nod_full_3$KO %in% sigtab_col_all_2_with_nod_3$KO])
  
  total_266_Fam <- length(Family_KOs[Family_KOs %in% input_table_3])
  sigtab_col_all_2_with_nod_full_3_sig <- sigtab_col_all_2_with_nod_full_3$KO[!sigtab_col_all_2_with_nod_full_3$KO %in% sigtab_col_all_2_with_nod_3$KO]
  lost_266 <- length(sigtab_col_all_2_with_nod_full_3_sig[sigtab_col_all_2_with_nod_full_3_sig %in% input_table_3])
  
  all_data <- data.frame(t(data.frame(c(paste(group), total_KOs, contrib_full, lost_KOs, total_266_Fam, lost_266 ))))
  new_fam_data <- rbind(new_fam_data, all_data)
}

row.names(new_fam_data) <- NULL
colnames(new_fam_data) <- c("Family", "Total no of KOs", "Total sig KOs", "Lost KOs in drop out", "Total no of 266 KOs", "Lost 266 KOs in drop out")

new_fam_data$`Total no of KOs` <- as.numeric(new_fam_data$`Total no of KOs`)
new_fam_data$`Total sig KOs` <- as.numeric(new_fam_data$`Total sig KOs`)
new_fam_data$`Lost KOs in drop out` <- as.numeric(new_fam_data$`Lost KOs in drop out`)
new_fam_data$`Total no of 266 KOs` <- as.numeric(new_fam_data$`Total no of 266 KOs`)
new_fam_data$`Lost 266 KOs in drop out` <- as.numeric(new_fam_data$`Lost 266 KOs in drop out`)

new_fam_data_3 <- melt(new_fam_data)

for (group in unique(new_fam_data_3$Family)){
  new_fam_data_3_sub <- new_fam_data_3[new_fam_data_3$Family == paste(group),]
  prop <- as.numeric(new_fam_data_3_sub$value[new_fam_data_3_sub$variable == "Lost KOs in drop out"])/as.numeric(new_fam_data_3_sub$value[new_fam_data_3_sub$variable == "Total sig KOs"])
  have_266 <- as.numeric(new_fam_data_3_sub$value[new_fam_data_3_sub$variable == "Lost 266 KOs in drop out"])
  not_have_266 <- as.numeric(new_fam_data_3_sub$value[new_fam_data_3_sub$variable == "Total no of 266 KOs"])
  not_have_266_2 <- (not_have_266 - have_266)
  
  binom_out <- binom.test(have_266, not_have_266, prop)
  pval <- binom_out$p.value
  
  new_data <- data.frame(t(data.frame(c(paste(group), "Percentage lost KOs vs sig", paste(prop)))))
  colnames(new_data) <- colnames(new_fam_data_3)
  row.names(new_data) <- NULL
  
  new_fam_data_3 <- rbind(new_fam_data_3,new_data)
  
  new_data <- data.frame(t(data.frame(c(paste(group), "Percentage lost KOs vs sig - 266", have_266/not_have_266))))
  colnames(new_data) <- colnames(new_fam_data_3)
  row.names(new_data) <- NULL
  
  new_fam_data_3 <- rbind(new_fam_data_3,new_data)
  
  new_data <- data.frame(t(data.frame(c(paste(group), "Binomial p-value", paste(round(pval,5))))))
  colnames(new_data) <- colnames(new_fam_data_3)
  row.names(new_data) <- NULL
  
  new_fam_data_3 <- rbind(new_fam_data_3,new_data)
}

data_wide_4_with_nod <- spread(new_fam_data_3, variable, value)
data_wide_4_2_with_nod <- data_wide_4_with_nod[data_wide_4_with_nod$Family != "Rhizobiaceae",]
data_wide_4_2_with_nod$Dataset <- "Nodulator"

data_together <- rbind(data_wide_4_2, data_wide_4_2_with_nod)

write.table(data_together, paste(results.dir, "Table_S8_KO_overview_Drop_out_exp.txt", sep = ""), sep = "\t", quote =F, col.names =T, row.names =F)
