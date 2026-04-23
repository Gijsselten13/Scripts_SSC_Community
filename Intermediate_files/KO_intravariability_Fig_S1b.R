library("dplyr") #Version 1.1.2
library("ape") #Version 5.7-1

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Script to generate KO_intravarability.tsv for 1000 one-strain-per-family SynComs - Necessary for S1b ======
KO_to_gene = read.table(paste(working_directory,"KO_to_gene.txt", sep = ""), header=TRUE,sep="\t")

table_all <- read.table(paste(working_directory, "KO_intravariability/simulations/random_family_sets.csv", sep = ""), header =T, row.names =1)
SSC_data <-   read.table(paste(working_directory, "KO_genome/KO_SSC.tsv",sep = ""), sep= "\t", header =T, row.names =1) 
colnames(SSC_data)[grep("M.16",colnames(SSC_data))] <- "M-16"
colnames(SSC_data)[grep("M.6",colnames(SSC_data))] <- "M-6"
colnames(SSC_data)[grep("M.10",colnames(SSC_data))] <- "M-10"
colnames(SSC_data)[grep("M.11_2",colnames(SSC_data))] <- "M-11_2"
colnames(SSC_data)[grep("^M.11$",colnames(SSC_data))] <- "M-11"
colnames(SSC_data) <- gsub("X", "", colnames(SSC_data))

new_KO_list_2 <- data.frame()

for (syncom in 1:1000){
  table_all_sub <- as.vector(unlist(as.vector(table_all[syncom,])))
  SSC_data_sub <- SSC_data[,colnames(SSC_data) %in% table_all_sub]
  
  KO_to_gene_sub <- KO_to_gene[KO_to_gene$isolate %in% table_all_sub,]
  
  KO_list <- row.names(SSC_data_sub)[rowSums(SSC_data_sub) != 0]
  
  for (KO_term in KO_list){
    print(syncom)
    print(KO_term)
    #Phylogenetic tree
    subset <- KO_to_gene_sub[ which(KO_to_gene_sub$kegg == paste0(KO_term)),]
    subset_2 <- subset[!duplicated(subset), ]
    subset_2$name <- paste(subset_2$isolate, subset_2$gene, sep = "_")
    subset_3 <- subset_2[,4:5]
    
    if (length(na.omit(subset_3$sequence)) > 1){
      subset_4 <- subset_3[!(is.na(subset_3$sequence) | subset_3$sequence==""), ]
      
      y <- strsplit(subset_4[,1],"")
      names(y) <- subset_4[,2]
      
      df.fasta <- ape::as.DNAbin(y)
      dist_matrix <- kmer::kdistance(df.fasta)
      
      # Produce dendrogram
      hclust = hclust(dist_matrix)
      
      my_tree <- as.phylo(hclust) 
      
      value_tree <- sum(my_tree$edge.length)/length(names(y))
      
    } else if (length(na.omit(subset_3$sequence)) > 0 ) { 
      value_tree <- 1
    } else {
      value_tree <- 0
    }
    
    vector <- c(paste(syncom),paste(KO_term), value_tree)
    new_KO_list_2 <- rbind(new_KO_list_2, vector)
  }
}

colnames(new_KO_list_2) <- c("SynCom", "KO" ,"value")
KO_out <- new_KO_list_2

new_KO_lists <- data.frame()

for (syncom in 1:1000){
  table_all_sub <- as.vector(unlist(as.vector(table_all[syncom,])))
  SSC_data_sub <- SSC_data[,colnames(SSC_data) %in% table_all_sub]
  
  KO_to_gene_sub <- KO_to_gene[KO_to_gene$isolate %in% table_all_sub,]
  
  KO_list <- row.names(SSC_data_sub)[rowSums(SSC_data_sub) != 0]
  
  KO_out_out <- KO_out[KO_out$SynCom == paste(syncom),]
  
  for (KO_term in KO_list){
    print(syncom)
    print(KO_term)
    subset <- KO_to_gene_sub[ which(KO_to_gene_sub$kegg == paste0(KO_term)),]
    subset_2 <- subset[!duplicated(subset), ]
    length_value <- length(subset_2$kegg)
    
    KO_out_out_2 <- KO_out_out[KO_out_out$KO == paste(KO_term),]
    
    total_branch_length_value <- as.numeric(KO_out_out_2$value) * length_value
    
    vector <- c(paste(syncom),paste(KO_term), total_branch_length_value, length_value)
    new_KO_lists <- rbind(new_KO_lists, vector)
  }
}

colnames(new_KO_lists) <- c("SynCom", "KO" ,"total_branch_length","no_of_genes")
write.table(new_KO_lists, paste(working_directory, "KO_intravariability/simulations/KO_intravariability.txt", sep = ""),col.names = T, row.names =F, quote =F, sep = "\t")
