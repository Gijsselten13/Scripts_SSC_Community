working_directory <- ""
results_sa.dir <- paste(working_directory, "Shiny_app/", sep = "")

###Generation of boxplots_full.txt - FC of isolates with pathway - Original dataset ======
top <- read.table(paste(working_directory,"Annotations/pathway_top.txt", sep = ""), header=F, sep="\t")
KO_to_pathway <- read.table(paste(working_directory,"Annotations/KO_to_pathway.txt", sep = ""), header=T, sep="\t")
KO_to_pathway$V3 <- top$V2[match(KO_to_pathway$V2, top$V1)]

KO_to_pathway_2 <- read.table(paste(working_directory,"Annotations/KO_to_pathway_unannotated_2.txt", sep = ""), header=F, sep="\t")
colnames(KO_to_pathway_2) <- c("KO","new_category")

for (KO in KO_to_pathway_2$KO){
  KO_to_pathway$V3[KO_to_pathway$V1 == paste(KO)] <- KO_to_pathway_2$new_category[KO_to_pathway_2$KO == paste(KO)]
}

new_table <- read.table(paste(working_directory, "Shiny_app/KOs.txt", sep = ""), header=T, sep="\t",quote ="")
KO_to_pathway$V4 <- new_table$Gene[match(KO_to_pathway$V1, new_table$KO)]
KO_to_pathway_3 <- KO_to_pathway[,c(1,4,2,3)]

KO_to_pathway_3$V5 <- NA

for (KO in top$V2){
  KO_to_pathway_3$V5[KO_to_pathway_3$V3 == paste(KO)] <- top$V3[top$V2 == paste(KO)]
}

KO_to_pathway_3$V1[is.na(KO_to_pathway_3$V1)] <- "Unknown"
KO_to_pathway_3$V2[is.na(KO_to_pathway_3$V2)] <- "Unknown"
KO_to_pathway_3$V3[is.na(KO_to_pathway_3$V3)] <- "Unknown"
KO_to_pathway_3$V4[is.na(KO_to_pathway_3$V4)] <- "Unknown"
KO_to_pathway_3$V5[is.na(KO_to_pathway_3$V5)] <- "Unknown"

colnames(KO_to_pathway_3) <- c("KO", "Description", "pathway", "Pathway_Description", "Category")

write.table(KO_to_pathway_3,paste(results_sa.dir, "Gene_descriptions.txt", sep = ""), col.names =T, row.names = F, sep = "\t", quote =F)

KO_table = read.table(paste(working_directory,"KO_genome/KO_SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
colnames(KO_table) <- gsub("X", "", colnames(KO_table))

pathways <- unique(KO_to_pathway_3$Pathway_Description)
pathways_2 <- pathways[pathways != "Unknown"]

Categories <- na.omit(unique(top$V3[top$V2 %in% pathways]))

groups <- Categories
SynComs <- c("AtSC","LjSC", "HvSC", "SSC")
Plant <- c("At", "Hv", "Lj")
hop_4 <- data.frame()

for (cat in groups){
  paths <- top$V2[top$V3 == paste(cat)]
  for (path in paths){
    KOs <- na.omit(KO_to_pathway$V1[KO_to_pathway$V3 == paste(path)])
    
    hop_2 <- data.frame()
    if(length(KOs) > 0 ){
      for (KO in KOs){
        for (syncom in SynComs){
          norm_SSC =read.table(paste(working_directory, "Isolate_tables/Original/", syncom,"_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
          norm_SSC_2 <- norm_SSC[,grep("ES", colnames(norm_SSC))]
          norm_SSC_3 <- norm_SSC_2[,grep(paste(syncom), colnames(norm_SSC_2))]
          norm_SSC_input <- norm_SSC[,grep("Input", colnames(norm_SSC))]
          
          KO_table_2 <- KO_table[row.names(KO_table) == paste(KO),]
          KO_table_3 <- KO_table_2[,colnames(KO_table_2) %in% row.names(norm_SSC_3)]
          
          for (plant in Plant){
            norm_SSC_4 <- norm_SSC_3[,grep(paste(plant, "_", sep = ""), colnames(norm_SSC_3))]
            norm_SSC_5 <- t(t(norm_SSC_4)/rowSums(t(norm_SSC_4)))
            
            No <- names(KO_table_3)[KO_table_3 == 0]
            Yes <- names(KO_table_3)[KO_table_3 != 0]
            
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
            
            norm_SSC_input_yes <- norm_SSC_input_2[row.names(norm_SSC_input_2) %in% Yes,]
            norm_SSC_input_no <- norm_SSC_input_2[row.names(norm_SSC_input_2) %in% No,]
            
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
            
            
            
            hop <- t(data.frame(c(paste(KO), Yes_sum, No_sum, Yes_sum_input, No_sum_input, length(Yes), length(No), paste(plant), paste(syncom))))
            hop_2 <- rbind(hop_2, hop)
          }
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
      
      hop_2$V10 <- (hop_2$V2/hop_2$V4)
      hop_2$V11 <- (hop_2$V3/hop_2$V5)
      
      value <- sum(hop_2$V6/(hop_2$V6 + hop_2$V7))/length(hop_2$V6)
      
      for (plant in Plant){
        for (syncom in SynComs){
          hop_2_5 <- hop_2[hop_2$V8 == paste(plant),]
          hop_2_6 <- hop_2_5[hop_2_5$V9 == paste(syncom),]
          value_path_yes <- sum(as.numeric(hop_2_6$V10))/length(hop_2_6$V10)
          value_path_no <- sum(as.numeric(hop_2_6$V11))/length(hop_2_6$V11)
          value_RA <-  sum(as.numeric(hop_2_6$V2))/length(hop_2_6$V2)
          hop_3 <- t(data.frame(c(paste(plant),value_RA, value_path_yes, value_path_no, paste(syncom), paste(path), paste(cat), value)))
          hop_4 <- rbind(hop_4, hop_3)
        }
      }
    }
  }
}

row.names(hop_4) <- NULL

write.table(hop_4, paste(results_sa.dir,"boxplots_full.txt", sep = ""), quote = F, col.names = T, row.names = T, sep = "\t")
