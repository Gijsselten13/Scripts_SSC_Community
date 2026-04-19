library("reshape2") #Version 1.4.4
library("dplyr") #Version 1.1.2
library("ggplot2") #Version 3.4.2

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Figure S15 - Panel A Heatmap =====
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

new_data_next <- data.frame()

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
    
    KO_iso_2_sum <- rowSums(KO_iso_2)
    KO_iso_2_sum_val <- length(names(KO_iso_2_sum)[KO_iso_2_sum != 0])
    
    KO_iso_3 <- data.frame(colSums(KO_iso_2))
    colnames(KO_iso_3) <- "No_of_KOs"
    
    #Observed
    table_bac_sub_2 <- table_bac_sub[row.names(table_bac_sub) %in% tax_df_sub_2, ]
    out <- data.frame(rowSums(table_bac_sub_2)/length(colnames(table_bac_sub_2)))
    colnames(out) <- "RA"
    
    KO_iso_3$RA <- out$RA[match(row.names(KO_iso_3), row.names(out))]
    KO_iso_3 <- na.omit(KO_iso_3)
    KO_iso_3$RA_prop <- KO_iso_3$RA/sum(KO_iso_3$RA)
    KO_iso_3$KO_prop <- KO_iso_3$No_of_KOs/KO_iso_2_sum_val
    
    KO_iso_3$SynCom <- paste(syncom)
    KO_iso_3$Family <- paste(family)
    
    
    new_data_next <- rbind(new_data_next,KO_iso_3)
  }
}

table_cor_2 <- data.frame()

for (family in unique(new_data_next$Family)){
  new_data_next_2 <- new_data_next[new_data_next$Family == paste(family),]
  for (syncom in unique(new_data_next$SynCom)){
    new_data_next_3 <- new_data_next_2[new_data_next_2$SynCom == paste(syncom),]
    new_data_next_3 <- new_data_next_3[new_data_next_3$RA_prop != 0,]
    
    if (length(row.names(new_data_next_3)) < 6){
      pval <- NA
      cor <- NA
    } else {
      val <- lm(formula = KO_prop ~ log(RA_prop), data = new_data_next_3)
      
      sum <- summary(val)
      cor <- sum$adj.r.squared
      pval_sum <- sum$coefficients
      pval <- pval_sum[8]
      
    }
    
    table_cor <- data.frame(t(data.frame(c(paste(syncom), cor, pval, paste(family)))))
    table_cor_2 <- rbind(table_cor_2, table_cor)
  }
}

row.names(table_cor_2) <- NULL
colnames(table_cor_2) <- c("SynCom", "R2", "Pvalue", "Family")

table_cor_2$R2 <- as.numeric(table_cor_2$R2)
table_cor_2$Pvalue <- as.numeric(table_cor_2$Pvalue)

table_cor_2$Sig <- ""
table_cor_2$Sig[table_cor_2$Pvalue < 0.05] <- "*"

table_cor_2$SynCom <- gsub("Syncom", "SynCom", table_cor_2$SynCom)

Plot_fam <- ggplot(table_cor_2, aes(SynCom, Family)) +
  geom_tile(height=0.98, mapping = aes(fill = R2)) +
  geom_text(aes(label = Sig), size =4) +
  scale_fill_gradient2(low = "#D55e00", mid = "white", high = "#56b4e9", midpoint =0, na.value = "lightgrey")+
  theme_classic() +
  labs(x ="Inoculum", y = "Family", fill = "R2") +
  theme(panel.background=element_blank(),panel.grid=element_blank(),axis.line.x=element_line(size=.5, colour="black"),axis.line.y=element_line(size=.5, colour="black"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black", size=7),legend.position="right",legend.text= element_text(size=10),text=element_text(family="sans", size=10))+
  theme(axis.text.x = element_text(size = 14, angle = 25,hjust=1),axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18), axis.text.y = element_text(size=14, face = rep("italic")), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=18)) +
  ggtitle("Root colonization versus functional diversity") +
  theme(plot.title = element_text(hjust = 0.5))
Plot_fam

pdf(paste(results.dir,"Figure_S15_AtSC_7_200_member_SCs_heatmap.pdf", sep=""), width=8, height=6)
print(Plot_fam)
dev.off()
