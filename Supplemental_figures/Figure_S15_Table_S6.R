library("dplyr") #Version 1.1.2
library("reshape2") #Version 1.4.4
library("mclust") #Version 6.1.1
library("ggplot2") #Version 3.4.2

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Figure S15 & Table S6 - panel B Line diagram =====
Families <- c("Burkholderiaceae","Rhizobiaceae","Xanthomonadaceae","Pseudomonadaceae","Caulobacteraceae","Beijerinckiaceae")

tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
KO_iso <- read.table(paste(working_directory,"KO_genome/KO_SSC.tsv", sep = ""), row.names =1, header =T)
tax_df <- tax_df[tax_df$SynCom == "AtSC",]

colnames(KO_iso) <- gsub("X","", colnames(KO_iso))
colnames(KO_iso) <- gsub("M.11_2","M-11_2", colnames(KO_iso))
colnames(KO_iso) <- gsub("M.11","M-11", colnames(KO_iso))
colnames(KO_iso) <- gsub("M.16","M-16", colnames(KO_iso))
colnames(KO_iso) <- gsub("M.10","M-10", colnames(KO_iso))
colnames(KO_iso) <- gsub("M.6","M-6", colnames(KO_iso))

#Filter out samples with high contamination
plant_reads <- read.table(paste(working_directory, "AtSC_7_SynCom_experiment/plant_reads.tsv", sep = ""), header =T)  
SynCom_reads <- read.table(paste(working_directory, "AtSC_7_SynCom_experiment/mapped_reads.tsv", sep = ""), header =T)

plant_reads$SynCom_reads_2 <- SynCom_reads$Full_200[match(plant_reads$sample_id, SynCom_reads$Sample)]
plant_reads$SynCom_reads <- SynCom_reads$Map_200[match(plant_reads$sample_id, SynCom_reads$Sample)]

plant_reads$non_plant_pairs <- plant_reads$total_reads-plant_reads$plant_reads

plant_reads$other_isolate_reads <- plant_reads$SynCom_reads_2 - plant_reads$SynCom_reads
plant_reads$Contaminant_reads <- plant_reads$non_plant_pairs - plant_reads$SynCom_reads_2

Hank_the_normalizer <- function(df,group,amount){
  df_2 <- df %>% dplyr::group_by_at(group) %>% dplyr::summarise(total=sum(.data[[amount]]))
  df_3 <- df_2$total
  names(df_3) <- df_2[[group]]
  df$total <- df_3[as.character(df[[group]])]
  df$Rel <- df[[amount]] / df$total
  return(df)
}

plant_reads_2 <- plant_reads[,c(1,3,5,7,8)]

plant_reads_melt <- melt(plant_reads_2)
plant_reads_melt_2 <- Hank_the_normalizer(plant_reads_melt,"sample_id","value")

removal_sam <- c()

for (sample in unique(plant_reads_melt_2$sample_id)){
  plant_reads_melt_2_sub <- plant_reads_melt_2[plant_reads_melt_2$sample_id == paste(sample),]
  value <- plant_reads_melt_2_sub$Rel[plant_reads_melt_2_sub$variable == "SynCom_reads"]/sum(plant_reads_melt_2_sub$Rel[plant_reads_melt_2_sub$variable != "plant_reads"])
  if (value <= 0.6){
    removal_sam <- c(removal_sam, paste(sample))
  }
}

#Microbiome data
table_bac <- read.table(paste(working_directory,"AtSC_7_SynCom_experiment/isolate_norm.txt", sep = ""), row.names =1, header =T)
table_bac_1 <- t(t(table_bac) / rowSums(t(table_bac)))
table_bac_2 <- table_bac_1[,!colnames(table_bac_1) %in% removal_sam] 

#metadata
metadata <- read.table(paste(working_directory,"AtSC_7_SynCom_experiment/metadata.txt", sep= ""), row.names =1, header =T)
metadata_2 <- metadata[grep("Syncom",metadata$Syncom),]
SynComs <- unique(metadata_2$Syncom)

#isolates
groups <- read.table(paste(working_directory,"AtSC_7_SynCom_experiment/SynCom_isolates.txt", sep = ""), row.names =1, header =F)

new_datafr <- data.frame()

for (syncom in SynComs){
  #Observed
  metadata_sub <- metadata_2[metadata_2$Syncom == paste(syncom),]
  metadata_sub_2 <- row.names(metadata_sub)[metadata_sub$Compartment != "inoculum"]
  
  isolates_in_SynCom <- as.vector(unlist(as.vector(groups[row.names(groups) == paste(syncom),])))
  tax_df_sub <- tax_df[tax_df$isolate %in% isolates_in_SynCom,]
  
  table_bac_sub <- data.frame(table_bac_2[,colnames(table_bac_2) %in% metadata_sub_2])
  
  for (family in Families){
    tax_df_sub_2 <- tax_df_sub$isolate[tax_df_sub$family == paste(family)]
    
    #Expected
    KO_iso_2 <- KO_iso[,colnames(KO_iso) %in% tax_df_sub_2]
    KO_iso_2[KO_iso_2 > 0] <- 1
    
    KO_iso_3 <- data.frame(colSums(KO_iso_2))
    colnames(KO_iso_3) <- "No_of_KOs"
    KO_iso_4 <- KO_iso_3 %>% dplyr::arrange(desc(No_of_KOs))
    
    m <- Mclust(KO_iso_4$No_of_KOs)     
    KO_iso_4$group <- m$classification
    
    #Observed
    table_bac_sub_2 <- table_bac_sub[row.names(table_bac_sub) %in% tax_df_sub_2, ]
    out <- data.frame(rowSums(table_bac_sub_2)/length(colnames(table_bac_sub_2)))
    colnames(out) <- "RA"
    
    for (cluster in unique(KO_iso_4$group)){
      KO_iso_4_sub <- row.names(KO_iso_4)[KO_iso_4$group == paste(cluster)]
      out_sub <- sum(out[row.names(out) %in% KO_iso_4_sub,])
      KO_iso_sub <- KO_iso_4[KO_iso_4$group == paste(cluster),]
      value_KO <- sum(KO_iso_sub$No_of_KOs)/length(KO_iso_sub$No_of_KOs)
      new_datafr <- rbind(new_datafr, data.frame(t(data.frame(c(paste(syncom),paste(family),paste(cluster),value_KO, out_sub)))))
    }
  }
}

row.names(new_datafr) <- NULL
colnames(new_datafr) <- c("SynCom","Family","Cluster","No_of_KO","Cum_RA")

new_datafr$SynFam <- paste(new_datafr$SynCom, new_datafr$Family, sep = "_")

new_datafr_2 <- data.frame()

for (group in unique(new_datafr$SynFam)){
  new_datafr_sub <- new_datafr[new_datafr$SynFam == paste(group),]
  new_datafr_sub_2 <- new_datafr_sub[order(as.numeric(new_datafr_sub$Cum_RA), decreasing = F), ]
  new_datafr_sub_2$order_RA <- 1:length(new_datafr_sub_2$SynFam)
  new_datafr_2 <- rbind(new_datafr_2,new_datafr_sub_2 )
}

row.names(new_datafr_2) <- NULL
colnames(new_datafr_2) <- c("SynCom","Family","Exp_order","No_of_KO","Cum_RA","SynFam", "Obs_order")

stat_data <- data.frame()

for (group in unique(new_datafr_2$Family)){
  new_datafr_2_sub <- new_datafr_2[new_datafr_2$Family == paste(group),]
  check_table <- table(new_datafr_2_sub$SynFam)
  check_table_2 <- names(check_table)[check_table == 1]
  new_datafr_2_sub <- new_datafr_2_sub[!new_datafr_2_sub$SynFam %in% check_table_2,]
  
  if (length(new_datafr_2_sub$Family) > 2){
    stat_out <- cor.test(as.numeric(new_datafr_2_sub$Exp_order), as.numeric(new_datafr_2_sub$Obs_order), method = "kendall")
    pval <- stat_out$p.value
    stat_data <- rbind(stat_data, data.frame(t(data.frame(c(paste(group), pval)))))
  } else {
    pval <- NA
    SC <- unique(new_datafr_2$SynCom[new_datafr_2$SynFam == paste(group)])
    FAM <- unique(new_datafr_2$Family[new_datafr_2$SynFam == paste(group)])
    stat_data <- rbind(stat_data, data.frame(t(data.frame(c(paste(SC), paste(FAM), pval)))))
  }
}

row.names(stat_data) <- NULL
colnames(stat_data) <- c("Family","Kendall correlation p-value")

new_datafr_2$SynCom <- gsub("Syncom", "SynCom", new_datafr_2$SynCom)
new_datafr_2$SynCom[new_datafr_2$SynCom == "SynCom11"] <- "SynCom1"

plot <- ggplot(new_datafr_2, aes(x = Exp_order, y = Obs_order, color = as.factor(Family), group = 1)) +
  geom_point(size = 3) +
  #geom_line() +   
  theme_classic() +
  labs(x = "Functional diversity rank",
       y = "Relative abundance rank",
       color = "Family") +
  ggtitle("Functional diversity vs Relative abundance") +
  theme(axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 14),
        axis.text.y = element_text(size=10),
        legend.title = element_text(size=16),
        legend.text = element_text(size=12,  face = rep("italic")),
        plot.title = element_text(size=18)) +
  facet_grid(SynCom ~ Family ,scales='free')
plot <- plot +  geom_smooth(method = "lm", se = FALSE) + theme(strip.text.x = element_text(face=rep("italic")))
plot

pdf(paste(results.dir,"Figure_S15_Line_plot_Func_vs_RA.pdf", sep=""), width=12, height=9)
print(plot)
dev.off()

stat_data$`SynCom number` <- c("4/7 SynComs", "5/5 SynComs", "3/7 SynComs", "7/7 SynComs", "7/7 SynComs", "5/5 SynComs")

write.table(new_datafr_2,paste(results.dir,"Data_func_vs_RA_output.txt", sep = ""), col.names =T, row.names =F, quote =F, sep = "\t")
write.table(stat_data,paste(results.dir,"Table_S6_func_RA_stats.txt", sep = ""), col.names =T, row.names =F, quote =F, sep = "\t")
