library("dplyr") #Version 1.1.2

working_directory <- ""
results_sa.dir <- paste(working_directory, "Shiny_app/", sep = "")

###Generation of Abundances_full.txt and Family_piedonuts_full.tsv - Gene/KO, Isolate, and family abundances of KOs (on basis of S20-S24) - Original dataset =====
hop_4 <- read.table(paste(working_directory, "Shiny_app/ternary_KOs_av_med_all.txt", sep = ""), header =T, sep ="\t", row.names = 1)

tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)

samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
samples_df_2$Condition[samples_df_2$Condition == "At"] <- "Arabidopsis"
samples_df_2$Condition[samples_df_2$Condition == "Lj"] <- "Lotus"
samples_df_2$Condition[samples_df_2$Condition == "Hv"] <- "Barley"

samples_df_3 <- samples_df_2[samples_df_2$Compartment == "ES",]
samples_df_4 <- samples_df_3[samples_df_3$Inoculum != "NS",]

SynComs <- c("AtSC","HvSC","LjSC", "SSC")

KO_table = read.table(paste(working_directory, "KO_genome/KO_SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
colnames(KO_table) <- gsub("X", "", colnames(KO_table))
colnames(KO_table) <- gsub("M.1", "M-1", colnames(KO_table))
colnames(KO_table) <- gsub("M.6", "M-6", colnames(KO_table))

together_2 <- data.frame()
fams <- data.frame()

for (gene in row.names(KO_table)){
  hop_4_sub <- hop_4[hop_4$KO == paste(gene),]
  print(gene)
  
  for (syncom in SynComs) {
    norm_KO = read.table(paste(working_directory, "KO_tables/Original/", syncom, ".tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
    norm_iso = read.table(paste(working_directory, "Isolate_tables/Original/", syncom, "_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
    
    if (syncom == "SSC"){
      KO_table_2 <- KO_table
    } else {
      KO_table_2 <- KO_table[,colnames(KO_table) %in% row.names(tax_df_2)[tax_df_2$SynCom == paste(syncom)]]
    }
    
    KO_table_3 <- KO_table_2[row.names(KO_table_2) == paste(gene),]
    KO_table_4 <- colSums(KO_table_3)
    KO_table_5 <- na.omit(names(KO_table_4)[KO_table_4 != 0])
    
    samples_df_5 <- samples_df_4[samples_df_4$Inoculum == paste(syncom),]
    samples_df_8 <- samples_df_5[!grepl("HL_orig", row.names(samples_df_5)),]
    samples_df_7 <- samples_df_5[!grepl("HL_orig", row.names(samples_df_5)),]
    
    norm_KO_2 <- norm_KO[,colnames(norm_KO) %in% row.names(samples_df_7)]
    
    norm_KO_3 <- t(t(norm_KO_2)/rowSums(t(norm_KO_2)))
    norm_KO_4 <- norm_KO_3[row.names(norm_KO_3) == paste(gene),]
    
    if (length(row.names(KO_table_3)) != 0){
      norm_iso_5 <- norm_iso[,colnames(norm_iso) %in% row.names(samples_df_7)]
      norm_iso_6 <- t(t(norm_iso_5)/rowSums(t(norm_iso_5)))
      norm_iso_7 <- norm_iso_6[row.names(norm_iso_6) %in% KO_table_5,]
      families <- unique(tax_df_2$family[row.names(tax_df_2) %in% KO_table_5])
      
      norm_iso_sub <- data.frame()
      if(length(families) > 0){
        for (fam in families){
          if (length(KO_table_5) == 1){
            norm_iso_8 <- norm_iso_7
          } else {
            if (is.matrix(norm_iso_7 == TRUE)){
              norm_iso_8 <- norm_iso_7[row.names(norm_iso_7) %in% KO_table_5[KO_table_5 %in% row.names(tax_df_2)[tax_df_2$family == paste(fam)]],]
            } else {
              norm_iso_8 <- norm_iso_7
            }
          }
          
          
          if (is.matrix(norm_iso_8) == FALSE){
            norm_iso_9 <- norm_iso_8
          } else {
            norm_iso_9 <- colSums(norm_iso_8)
          }
          
          fam_value <- sum(norm_iso_9)/length(norm_iso_9)
          fam_data <- t(data.frame(c(paste(fam), fam_value, paste(syncom))))
          norm_iso_sub <- rbind(norm_iso_sub, fam_data)
        } 
        
        norm_iso_sub$V4 <- paste(gene)
        fams <- rbind(fams, norm_iso_sub)
      }
      
      
      if (length(row.names(KO_table_3)) != 0){
        norm_iso_10 <- norm_iso[,colnames(norm_iso) %in% row.names(samples_df_8)]
        norm_iso_11 <- t(t(norm_iso_10)/rowSums(t(norm_iso_10)))
        norm_iso_12 <- norm_iso_11[row.names(norm_iso_11) %in% KO_table_5,]
      } 
      
      Groups <- c("Arabidopsis", "Barley", "Lotus")
      
      norm_KO_2 <- norm_KO[,colnames(norm_KO) %in% row.names(samples_df_8)]
      
      norm_KO_3 <- t(t(norm_KO_2)/rowSums(t(norm_KO_2)))
      norm_KO_4 <- norm_KO_3[row.names(norm_KO_3) == paste(gene),]
      
      for (group in Groups){
        
        samples_df_9 <- samples_df_8[samples_df_8$Condition == paste(group),]
        
        norm_KO_5 <- norm_KO_4[names(norm_KO_4) %in% row.names(samples_df_9)]
        norm_KO_6 <- norm_KO_5
        
        if (syncom == "SSC"){
          if (length(KO_table_5) > 1){
            At_iso <- KO_table_5[KO_table_5 %in% row.names(tax_df_2)[tax_df_2$SynCom == "AtSC"]] 
            Lj_iso <- KO_table_5[KO_table_5 %in% row.names(tax_df_2)[tax_df_2$SynCom == "LjSC"]] 
            Hv_iso <- KO_table_5[KO_table_5 %in% row.names(tax_df_2)[tax_df_2$SynCom == "HvSC"]] 
            
            if (is.matrix(norm_iso_12) == TRUE){
              norm_iso_12_At <- norm_iso_12[row.names(norm_iso_12) %in% At_iso,]
              norm_iso_12_Lj <- norm_iso_12[row.names(norm_iso_12) %in% Lj_iso,]
              norm_iso_12_Hv <- norm_iso_12[row.names(norm_iso_12) %in% Hv_iso,]
            } else {
              
              if (length(At_iso) > 0){
                norm_iso_12_At <- norm_iso_12
              } else {
                norm_iso_12_At <- NULL
              }
              if (length(Lj_iso) > 0){
                norm_iso_12_Lj <- norm_iso_12
              } else {
                norm_iso_12_Lj <- NULL
              }
              if (length(Hv_iso) > 0){
                norm_iso_12_Hv <- norm_iso_12
              } else {
                norm_iso_12_Hv <- NULL
              }
            }
            
            if (length(At_iso) > 1 & is.matrix(norm_iso_12_At) == TRUE){
              norm_iso_13_At <- colSums(norm_iso_12_At)
            } else {
              norm_iso_13_At <- sum(norm_iso_12_At)
            }
            
            if (length(Hv_iso) > 1 & is.matrix(norm_iso_12_Hv) == TRUE){
              norm_iso_13_Hv <- colSums(norm_iso_12_Hv)
            } else {
              norm_iso_13_Hv <- sum(norm_iso_12_Hv)
            }
            
            if (length(Lj_iso) > 1 & is.matrix(norm_iso_12_Lj) == TRUE){
              norm_iso_13_Lj <- colSums(norm_iso_12_Lj)
            } else {
              norm_iso_13_Lj <- sum(norm_iso_12_Lj)
            }
            
            norm_iso_14_At <- norm_iso_13_At[names(norm_iso_13_At) %in% row.names(samples_df_9)]
            norm_iso_14_Lj <- norm_iso_13_Lj[names(norm_iso_13_Lj) %in% row.names(samples_df_9)]
            norm_iso_14_Hv <- norm_iso_13_Hv[names(norm_iso_13_Hv) %in% row.names(samples_df_9)]
            
            value_iso_At <- sum(norm_iso_14_At)/length(names(norm_iso_14_At))
            value_iso_Lj <- sum(norm_iso_14_Lj)/length(names(norm_iso_14_Lj))
            value_iso_Hv <- sum(norm_iso_14_Hv)/length(names(norm_iso_14_Hv))
            
          } else if (length(KO_table_5) == 1) {
            plant_sel <- tax_df_2$SynCom[row.names(tax_df_2) == paste(KO_table_5)]
            norm_iso_14 <- norm_iso_12[names(norm_iso_12) %in% row.names(samples_df_9)]
            
            if (plant_sel == "AtSC"){
              value_iso_At <- sum(norm_iso_14)/length(names(norm_iso_14))
              value_iso_Lj <- 0
              value_iso_Hv <- 0
            } else if (plant_sel == "HvSC"){
              value_iso_At <- 0
              value_iso_Lj <- 0
              value_iso_Hv <- sum(norm_iso_14)/length(names(norm_iso_14))
            } else if (plant_sel == "LjSC"){
              value_iso_At <- 0
              value_iso_Lj <- sum(norm_iso_14)/length(names(norm_iso_14))
              value_iso_Hv <- 0
            }
          } else {
            value_iso_Hv <- 0
            value_iso_At <- 0
            value_iso_Lj <- 0
          }
        } else {
          if (length(KO_table_5) > 1){
            if(is.matrix(norm_iso_12) == TRUE){
              norm_iso_13 <- colSums(norm_iso_12)
              norm_iso_14 <- norm_iso_13[names(norm_iso_13) %in% row.names(samples_df_9)]
              value_iso <- sum(norm_iso_14)/length(names(norm_iso_14))
            } else {
              norm_iso_13 <- norm_iso_12[names(norm_iso_12) %in% row.names(samples_df_9)]
              norm_iso_14 <- sum(norm_iso_13)
              value_iso <- norm_iso_14
            }
          } else if (length(KO_table_5) == 1) {
            norm_iso_14 <- norm_iso_12[names(norm_iso_12) %in% row.names(samples_df_9)]
            value_iso <- sum(norm_iso_14)/length(names(norm_iso_14))
          } else {
            value_iso <- 0
          }
        }
        
        norm_KO_7 <- norm_KO_6[names(norm_KO_6) %in% row.names(samples_df_9)]
        
        value <- sum(norm_KO_7)/length(names(norm_KO_7))
        
        if (syncom == "SSC"){
          together_3 <- t(data.frame(c(paste(gene), paste(syncom), paste(group), as.numeric(value), "AtSC", as.numeric(value_iso_At))))
          together_4 <- t(data.frame(c(paste(gene), paste(syncom), paste(group), as.numeric(value), "HvSC", as.numeric(value_iso_Hv))))
          together_5 <- t(data.frame(c(paste(gene), paste(syncom), paste(group), as.numeric(value), "LjSC", as.numeric(value_iso_Lj))))
          together <- rbind(together_3, together_4, together_5)
        } else {
          together <- t(data.frame(c(paste(gene), paste(syncom), paste(group), as.numeric(value), paste(syncom), as.numeric(value_iso))))
        }
        
        row.names(together) <- NULL
        
        together_2 <- rbind(together_2, together)
      }
    }
  }
}

row.names(together_2) <- NULL
colnames(together_2) <- c("Gene", "SynCom", "Plant", "RA_KO", "Origin", "RA_Iso")

together_2$RA_KO[together_2$RA_KO == NaN] <- 0
together_2$RA_Iso[together_2$RA_Iso == NaN] <- 0

row.names(fams) <- NULL
colnames(fams) <- c("Family", "RA", "SynCom", "Gene")

write.table(together_2, paste(results_sa.dir, "Abundances_full.tsv", col.names =T, row.names =F, sep = "\t", quote =F))
write.table(fams, paste(results_sa.dir, "Family_piedonuts_full.tsv", col.names =T, row.names =F, sep = "\t", quote =F))
