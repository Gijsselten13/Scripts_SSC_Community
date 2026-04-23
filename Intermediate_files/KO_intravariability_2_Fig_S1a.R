library("dplyr") #Version 1.1.2
library("ape") #Version 5.7-1

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Script to generate KO_intravariability_2.tsv for SSC KO intravariability - Necessary for S1a =====
KO_to_gene = read.table(paste(working_directory,"KO_to_gene.txt", sep = ""), header=TRUE,sep="\t")
KO_to_gene$SynCom <- NA
list_AtSC =read.table(paste(working_directory,"Isolate_tables/Original/AtSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
list_AtSC <- row.names(list_AtSC)
list_LjSC =read.table(paste(working_directory,"Isolate_tables/Original/LjSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
list_LjSC <- row.names(list_LjSC)
list_HvSC =read.table(paste(working_directory,"Isolate_tables/Original/HvSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
list_HvSC <- row.names(list_HvSC)

KO_to_gene$SynCom[KO_to_gene$isolate %in% list_AtSC] <- "AtSC"
KO_to_gene$SynCom[KO_to_gene$isolate %in% list_LjSC] <- "LjSC"
KO_to_gene$SynCom[KO_to_gene$isolate %in% list_HvSC] <- "HvSC"
KO_list <- unique(KO_to_gene$kegg)

new_KO_list_3 <- data.frame(matrix(NA, ncol = 9))
colnames(new_KO_list_3) <- c("KO","AtSC_total_branch_length","LjSC_total_branch_length","HvSC_total_branch_length","SSC_total_branch_length", "No_of_genes_AtSC", "No_of_genes_LjSC","No_of_genes_HvSC", "No_of_genes_SSC")
new_KO_list_4 <- new_KO_list_3[-1,]

for (KO_term in KO_list){
  print(KO_term)
  #Phylogenetic tree
  subset <- KO_to_gene[ which(KO_to_gene$kegg == paste0(KO_term)),]
  subset_2 <- subset[!duplicated(subset), ]
  subset_2$name <- paste(subset_2$isolate, subset_2$gene, sep = "_")
  
  subset_AtSC <- subset_2[subset_2$SynCom == "AtSC",]
  subset_LjSC <- subset_2[subset_2$SynCom == "LjSC",]
  subset_HvSC <- subset_2[subset_2$SynCom == "HvSC",]
  
  subset_3 <- subset_AtSC[,4:6]
  subset_3 <- subset_3 %>% dplyr::select (-SynCom)
  
  if (length(na.omit(subset_3$sequence)) > 1){
    subset_4 <- subset_3[!(is.na(subset_3$sequence) | subset_3$sequence==""), ]
    
    y <- strsplit(subset_4[,1],"")
    names(y) <- subset_4[,2]
    
    df.fasta <- ape::as.DNAbin(y)
    dist_matrix <- kmer::kdistance(df.fasta)
    
    # Produce dendrogram
    hclust = hclust(dist_matrix)
    
    my_tree <- as.phylo(hclust) 
    len_AtSC <- sum(my_tree$edge.length)
    tree_AtSC <- length(names(y))
    
  } else if (length(na.omit(subset_3$sequence)) > 0 ) { 
    len_AtSC <- 1
    tree_AtSC <- 1
  } else {
    len_AtSC <- 0
    tree_AtSC <- 0
  }
  
  subset_3 <- subset_LjSC[,4:6]
  subset_3 <- subset_3 %>% dplyr::select (-SynCom)
  
  if (length(na.omit(subset_3$sequence)) > 1){
    subset_4 <- subset_3[!(is.na(subset_3$sequence) | subset_3$sequence==""), ]
    
    y <- strsplit(subset_4[,1],"")
    names(y) <- subset_4[,2]
    
    df.fasta <- ape::as.DNAbin(y)
    dist_matrix <- kmer::kdistance(df.fasta)
    
    # Produce dendrogram
    hclust = hclust(dist_matrix)
    
    my_tree <- as.phylo(hclust) 
    len_LjSC <- sum(my_tree$edge.length)
    tree_LjSC <- length(names(y))
    
  } else if (length(na.omit(subset_3$sequence)) > 0) { 
    len_LjSC <- 1
    tree_LjSC <- 1
  } else {
    len_LjSC <- 0
    tree_LjSC <- 0
  }
  
  subset_3 <- subset_HvSC[,4:6]
  subset_3 <- subset_3 %>% dplyr::select (-SynCom)
  
  if (length(na.omit(subset_3$sequence)) > 1){
    subset_4 <- subset_3[!(is.na(subset_3$sequence) | subset_3$sequence==""), ]
    
    y <- strsplit(subset_4[,1],"")
    names(y) <- subset_4[,2]
    
    df.fasta <- ape::as.DNAbin(y)
    dist_matrix <- kmer::kdistance(df.fasta)
    
    # Produce dendrogram
    hclust = hclust(dist_matrix)
    
    my_tree <- as.phylo(hclust) 
    len_HvSC <- sum(my_tree$edge.length)
    tree_HvSC <- length(names(y))
    
  } else if (length(na.omit(subset_3$sequence)) > 0 ){ 
    len_HvSC <- 1
    tree_HvSC <- 1
  } else {
    len_HvSC <- 0
    tree_HvSC <- 0
  }
  
  subset_3 <- subset_2[,4:6]
  subset_3 <- subset_3 %>% dplyr::select (-SynCom)
  
  if (length(na.omit(subset_3$sequence)) > 2){
    subset_4 <- subset_3[!(is.na(subset_3$sequence) | subset_3$sequence==""), ]
    
    y <- strsplit(subset_4[,1],"")
    names(y) <- subset_4[,2]
    
    df.fasta <- ape::as.DNAbin(y)
    dist_matrix <- kmer::kdistance(df.fasta)
    
    # Produce dendrogram
    hclust = hclust(dist_matrix)
    
    my_tree <- as.phylo(hclust) 
    len_SSC <- sum(my_tree$edge.length)
    tree_SSC <- length(names(y))
  } else if (length(na.omit(subset_3$sequence)) > 0 ) { 
    len_SSC <- 1
    tree_SSC <- 1
  } else {
    len_SSC <- 0
    tree_SSC <- 0
  }

  vector_2 <- c(paste(KO_term, sep = ""), len_AtSC,len_LjSC,len_HvSC,len_SSC,tree_AtSC, tree_LjSC, tree_HvSC, tree_SSC)
  new_KO_list_4 <- rbind(new_KO_list_4, vector_2)
}

colnames(new_KO_list_4) <- c("KO","AtSC_total_branch_length","LjSC_total_branch_length","HvSC_total_branch_length","SSC_total_branch_length", "No_of_genes_AtSC", "No_of_genes_LjSC","No_of_genes_HvSC", "No_of_genes_SSC")
write.table(new_KO_list_4, paste(working_directory,"KO_intravariability/KO_intrafunctionality_2.tsv", sep = ""), sep = "\t", quote = F)
