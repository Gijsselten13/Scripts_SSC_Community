working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Creating Figure_5d_validation_in_vivo_Fam_drop_out_with_nod.txt =====
pathway_selection <- c("ABC transporters", "Quorum sensing", "Two-component system", "Purine metabolism", "Oxidative phosphorylation", "Pentose and glucuronate interconversions", "Porphyrin metabolism", "Galactose metabolism", "Cell cycle - Caulobacter", "Exopolysaccharide biosynthesis")

#Insilico drop-out - recalculating pathways
Families <- c("Burkholderiaceae", "Caulobacteraceae", "Pseudomonadaceae", "other", "all")

#Invivo drop-out - recalculating pathways
top <- read.table(paste(working_directory, "Annotations/pathway_top.txt", sep = ""), header=F, sep="\t")
KO_to_pathway <- read.table(paste(working_directory, "Annotations/KO_to_pathway.txt", sep = ""), header=T, sep="\t")
KO_to_pathway$V3 <- top$V2[match(KO_to_pathway$V2, top$V1)]

KO_to_pathway_2 <- read.table(paste(working_directory, "Annotations/KO_to_pathway_unannotated_2.txt", sep = ""), header=F, sep="\t")
colnames(KO_to_pathway_2) <- c("KO","new_category")

for (KO in KO_to_pathway_2$KO){
  KO_to_pathway$V3[KO_to_pathway$V1 == paste(KO)] <- KO_to_pathway_2$new_category[KO_to_pathway_2$KO == paste(KO)]
}

KO_table = read.table(paste(working_directory, "KO_genome/KO_LjSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
KO_table_2 <- t(KO_table)

#Significant KOs from SSC DEseq2 analysis
input_table <- read.table(paste(working_directory, "DESeq2/Sig_KO_all.txt", sep = ""), header=T, sep="\t")
input_table_2 <- input_table$KO[input_table$Plant == "Lotus" & input_table$SynCom == "LjSC"]

SynComs <- c("LDT1","LDT2","LDT3","LDT4","LDT5","LDT6")
hop_4 <- data.frame()

for (path in pathway_selection){
  KOs <- as.vector(na.omit(KO_to_pathway$V1[KO_to_pathway$V3 == paste(path)]))
  KOs_2 <- KOs[KOs %in% input_table_2]
  if(length(KOs_2) > 0){
    hop_2 <- data.frame()
    if(length(KOs_2) > 0 ){
      for (KO in KOs_2){
        for (syncom in SynComs){
          norm_SSC =read.table(paste(working_directory, "LjSC_Family_drop_out_experiment/Data_with_synthetic_input.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
          norm_SSC_2 <- norm_SSC[,grep("ROOT", colnames(norm_SSC))]
          norm_SSC_3 <- norm_SSC_2[,grep(paste(syncom), colnames(norm_SSC_2))]
          norm_SSC_input <- norm_SSC[,grep("INPUT", colnames(norm_SSC))]
          
          KO_table_3 <- data.frame(KO_table_2[,colnames(KO_table_2) == paste(KO)])
          colnames(KO_table_3) <- "KO"
          KO_table_4 <- KO_table_3[row.names(KO_table_3) %in% row.names(norm_SSC_3),]
          names(KO_table_4) <- row.names(KO_table_3)[row.names(KO_table_3) %in% row.names(norm_SSC_3)]
          
          norm_SSC_5 <- t(t(norm_SSC_3)/rowSums(t(norm_SSC_3)))
          
          No <- names(KO_table_4)[KO_table_4 == 0]
          Yes <- names(KO_table_4)[KO_table_4 != 0]
          
          norm_SSC_yes <- norm_SSC_5[row.names(norm_SSC_5) %in% Yes,]
          norm_SSC_no <- norm_SSC_5[row.names(norm_SSC_5) %in% No,]
          
          if (length(Yes) > 1){
            Yes_sum <- sum(colSums(norm_SSC_yes))/length(colSums(norm_SSC_yes))
          } else if (length(Yes) == 1){
            Yes_sum <- sum(norm_SSC_yes)/length(norm_SSC_yes)
          } else {
            Yes_sum <- 0
          }
          
          if (length(No) > 1){
            No_sum <- sum(colSums(norm_SSC_no))/length(colSums(norm_SSC_no))
          } else if (length(No) == 1){
            No_sum <- sum(norm_SSC_no)/length(norm_SSC_no)
          } else {
            No_sum <- 0
          }
          
          norm_SSC_input_2 <- t(t(norm_SSC_input)/rowSums(t(norm_SSC_input)))
          norm_SSC_input_3 <- norm_SSC_input_2[,grep(paste(syncom), colnames(norm_SSC_input_2))]
          norm_SSC_input_yes <- norm_SSC_input_3[row.names(norm_SSC_input_3) %in% Yes,]
          norm_SSC_input_no <- norm_SSC_input_3[row.names(norm_SSC_input_3) %in% No,]
          
          if (length(Yes) > 1){
            Yes_sum_input <- sum(colSums(norm_SSC_input_yes))/length(colSums(norm_SSC_input_yes))
          } else if (length(Yes) == 1){
            Yes_sum_input <- sum(norm_SSC_input_yes)/length(norm_SSC_input_yes)
          } else {
            Yes_sum_input <- 0
          }
          
          if (length(No) > 1){
            No_sum_input <- sum(colSums(norm_SSC_input_no))/length(colSums(norm_SSC_input_no))
          } else if (length(No) == 1){
            No_sum_input <- sum(norm_SSC_input_no)/length(norm_SSC_input_no)
          } else {
            No_sum_input <- 0
          }
          
          hop <- t(data.frame(c(paste(KO), Yes_sum, No_sum, Yes_sum_input, No_sum_input, length(Yes), length(No),paste(syncom))))
          hop_2 <- rbind(hop_2, hop)
        }
      }
      
      hop_2$V2 <- as.numeric(hop_2$V2)
      hop_2$V3 <- as.numeric(hop_2$V3)
      hop_2$V4 <- as.numeric(hop_2$V4)
      hop_2$V5 <- as.numeric(hop_2$V5)
      hop_2$V6 <- as.numeric(hop_2$V6)
      hop_2$V7 <- as.numeric(hop_2$V7)
      
      int_value <- min(hop_2$V4[hop_2$V4 != 0])
      int_value_2 <- min(hop_2$V5[hop_2$V5 != 0])
      
      hop_2$V4[hop_2$V4 == 0] <- int_value
      hop_2$V5[hop_2$V5 == 0] <- int_value_2
      
      hop_2$V9 <- (hop_2$V2/hop_2$V4)
      hop_2$V10 <- (hop_2$V3/hop_2$V5)
      
      vec <- hop_2$V6/(hop_2$V6 + hop_2$V7)
      vec[is.nan(vec)] <- 0
      
      value <- sum(vec)/length(hop_2$V6)
      
      for (syncom in SynComs){
        hop_2_6 <- hop_2[hop_2$V8 == paste(syncom),]
        value_path_yes <- sum(as.numeric(hop_2_6$V9))/length(hop_2_6$V9)
        value_path_no <- sum(as.numeric(hop_2_6$V10))/length(hop_2_6$V10)
        value_RA <-  sum(as.numeric(hop_2_6$V2))/length(hop_2_6$V2)
        hop_3 <- t(data.frame(c(value_RA, value_path_yes, value_path_no, paste(syncom), paste(path), value)))
        hop_4 <- rbind(hop_4, hop_3)
      }
    }
  }
}

row.names(hop_4) <- NULL
colnames(hop_4) <- c("RA", "Present_fold_change", "Absent_fold_change", "Family", "Pathway", "Percentage")

write.table(hop_4, paste(working_directory, "LjSC_Family_drop_out_experiment/Figure_5d_validation_in_vivo_Fam_drop_out_with_nod.txt", sep = ""), quote = F, col.names = T, row.names = F, sep = "\t")
