working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Table S6 - DESeq2 overlap across SynComs per plant ======
#With dominators
input_table <- read.table(paste(working_directory, "/DESeq2/Sig_KO_all.txt", sep = ""), header=T, sep="\t")

plants <- c("Arabidopsis", "Barley", "Lotus")
SynComs <- c("AtSC", "HvSC", "LjSC", "SSC")

KO_table_2 <- data.frame()

for (plant in plants){
  for (syncom in SynComs){
    
    if (syncom == "AtSC"){
      SynComs2 <- c("HvSC","LjSC","SSC")
    } else if (syncom == "HvSC") {
      SynComs2 <- c("AtSC","LjSC","SSC")
    } else if (syncom == "LjSC") {
      SynComs2 <- c("AtSC","HvSC","SSC")
    } else {
      SynComs2 <- c("AtSC","HvSC","LjSC")
    }
    
    input_table_2 <- input_table[input_table$Plant == paste(plant),]
    input_table_3 <- input_table_2[input_table_2$SynCom == paste(syncom),]
    
    value1 <- length(input_table_3$KO)
    
    for (SC2 in SynComs2){
      norm_KO = read.table(paste(working_directory,"KO_tables/Original/", SC2, ".tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
      value2 <- length(row.names(norm_KO)[row.names(norm_KO) %in% input_table_3$KO])
      input_table_4 <- input_table_2[input_table_2$SynCom == paste(SC2),]
      value3 <- length(input_table_3$KO[input_table_3$KO %in% input_table_4$KO])
      
      new_data <- data.frame(t(data.frame(c(paste(plant), paste(syncom), value1, paste(SC2), value2, value3))))
      
      KO_table_2 <- rbind(KO_table_2, new_data)
    }
  }
}

row.names(KO_table_2) <- NULL
colnames(KO_table_2) <- c("Plant", "SynCom", "No_of_sig_KOs", "SynCom_comparison", "Overlap_present", "Overlap_significant")

write.table(KO_table_2, paste(results.dir, "Table_S6_KO_overlap_with_dom.txt", sep = ""), sep = "\t", row.names = T, col.names =T, quote =F)

#Deseq2 file without dominators
input_table <- read.table(paste(working_directory, "/DESeq2/Sig_KO_all_no_nod_rhizo.txt", sep = ""), header=T, sep="\t")

plants <- c("Arabidopsis", "Barley", "Lotus")
SynComs <- c("AtSC", "HvSC", "LjSC", "SSC")

KO_table_2 <- data.frame()

for (plant in plants){
  for (syncom in SynComs){
    
    if (syncom == "AtSC"){
      SynComs2 <- c("HvSC","LjSC","SSC")
    } else if (syncom == "HvSC") {
      SynComs2 <- c("AtSC","LjSC","SSC")
    } else if (syncom == "LjSC") {
      SynComs2 <- c("AtSC","HvSC","SSC")
    } else {
      SynComs2 <- c("AtSC","HvSC","LjSC")
    }
    
    input_table_2 <- input_table[input_table$plant == paste(plant),]
    input_table_3 <- input_table_2[input_table_2$SynCom == paste(syncom),]
    
    value1 <- length(input_table_3$KO)
    
    for (SC2 in SynComs2){
      norm_KO = read.table(paste(working_directory,"KO_tables/Original/", SC2, ".tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
      value2 <- length(row.names(norm_KO)[row.names(norm_KO) %in% input_table_3$KO])
      input_table_4 <- input_table_2[input_table_2$SynCom == paste(SC2),]
      value3 <- length(input_table_3$KO[input_table_3$KO %in% input_table_4$KO])
      
      new_data <- data.frame(t(data.frame(c(paste(plant), paste(syncom), value1, paste(SC2), value2, value3))))
      
      KO_table_2 <- rbind(KO_table_2, new_data)
    }
  }
}

row.names(KO_table_2) <- NULL
colnames(KO_table_2) <- c("Plant", "SynCom", "No_of_sig_KOs", "SynCom_comparison", "Overlap_present", "Overlap_significant")

write.table(KO_table_2, paste(results.dir, "Table_S6_KO_overlap_without_dom.txt", sep = ""), sep = "\t", row.names = T, col.names =T, quote =F)
