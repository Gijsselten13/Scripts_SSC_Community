library("dplyr") #Version 1.1.2

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Script to generate pangenome_order.tsv for SSC - Necessary for 1d =====
SynComs <- c("AtSC", "HvSC", "LjSC", "SSC")
table_4 <- data.frame()

for (inoculum in SynComs) {
  table <- read.table(paste(working_directory, "KO_genome/KO_",inoculum,".tsv", sep = ""), sep= "\t", header =T) 
  row.names(table) <- table$sequence
  table_2 <- table %>% dplyr::select (-sequence)
  
  table_2[is.na(table_2)] <- 0
  
  vector <- colnames(table_2)
  new_list <- list()
  
  new_order <- c()
  
  vector_3 <- vector
  i = 1
  removing_KOs <- c()
  
  while(i <=length(vector)){
    vector_2 <- vector_3
    new_list <- list()
    for (isolate in vector_2) {
      sub <- table_2[colnames(table_2) == paste(isolate)]
      sub_2 <- row.names(sub)[sub>0]
      new_list[[paste(isolate)]] <- sub_2
    }
    
    new <- as.data.frame(names(new_list))
    colnames(new) <- "isolate"
    
    no_table <- as.data.frame(unique(unlist(new_list)))
    colnames(no_table) <- "KO_term"
    no_table_2 <- no_table[!no_table$KO_term %in% removing_KOs,]
    new_list_2 <- unlist(new_list)
    new_list_3 <- new_list_2[!new_list_2 %in% removing_KOs]
    
    new_list_4 <- as.data.frame(table(new_list_3))
    
    isolates <- as.vector(new)
    
    for (isolate in isolates$isolate) {
      hop <- new_list[paste(isolate)]
      new_vector <- as.vector(unlist(hop))
      new_vector_2 <- new_vector[!new_vector %in% removing_KOs]
      other_new <- as.data.frame(new_vector_2)
      colnames(other_new) <- "KO_term"
      for (another in new_vector_2){
        value <- new_list_4$Freq[new_list_4$new_list_3 == paste(another)]
        other_new$number[other_new$KO_term == paste(another)] <- value
      }
      new_value <- length(other_new$KO_term[other_new$number == 1])
      new$number[new$isolate == paste(isolate)] <- new_value
    }
    new_2 <- new[order(new$number, decreasing = T),]
    new_order <- c(new_order,new_2$isolate[1])
    
    isolate_3 <- new_2$isolate[1]
    sub_3 <- table_2[colnames(table_2) == paste(isolate_3)]
    sub_4 <- row.names(sub_3)[sub_3>0]
    removing_KOs <- c(removing_KOs, sub_4)
    removing_KOs <- unique(removing_KOs)
    vector_3 <- vector_2[!vector_2 == new_2$isolate[1]]
    i <- i +1
    
  } 
  
  SynCom_vec <- vector()
  other_vec <- vector()
  
  for (column in new_order) {
    table_sub <- table_2[,colnames(table_2) == paste(column)]
    table_sub_2 <- as.data.frame(table_sub)
    row.names(table_sub_2) <- row.names(table_2)
    table_sub_3 <- row.names(table_sub_2)[table_sub_2[1] > 0]
    SynCom_vec <- c(SynCom_vec,table_sub_3)
    new_value <- length(unique(SynCom_vec))
    other_vec <- c(other_vec,new_value)
  }
  
  table_3 <- t(as.data.frame(other_vec))
  colnames(table_3) <- 1:length(other_vec)
  row.names(table_3) <- paste(inoculum)
  
  late_vector <- rep(NA,1000 -length(table_3))
  table_3 <- c(table_3,late_vector)
  table_4 <- rbind(table_4,table_3)
}

row.names(table_4) <- SynComs
colnames(table_4) <- 1:1000
table_4$SynCom <- row.names(table_4)
table_5 <- melt(table_4)
table_6 <- na.omit(table_5)

write.table(table_6, paste(working_directory, "KO_intravariability/pangenome_order.tsv", sep=""), sep = "\t", quote = F, row.names = F, col.names = T)
