working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Script to generate top70_isolates_no_dom.txt =====
SynComs <- c("AtSC", "HvSC", "LjSC", "SSC")
Plant <- c("At", "Hv", "Lj")

samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names

list_70 <- list()

uno <- 1

for (syncom in SynComs){
  norm_table = read.table(paste(working_directory, "Isolate_tables/No_dominances/",syncom, "_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
  
  #Filter for root samples
  samples_df_sub <- samples_df[samples_df$Inoculum == paste(syncom),]
  samples_df_sub_2 <- samples_df_sub[samples_df_sub$Compartment == "ES",]
  
  for (plant in Plant){
    samples_df_sub_3 <- samples_df_sub_2[samples_df_sub_2$Condition == paste(plant),]
    norm_table_2 <- norm_table[, colnames(norm_table) %in% samples_df_sub_3$sample_id]
    norm_table_3 <- t(t(norm_table_2)/rowSums(t(norm_table_2)))
    
    hop <- data.frame(rowSums(norm_table_3)/length(colnames(norm_table_3)))
    colnames(hop) <- "Rel"
    
    hop_2 <- data.frame(hop[order(hop$Rel, decreasing=T),])
    row.names(hop_2) <- row.names(hop)[order(hop$Rel, decreasing=T)]
    
    colnames(hop_2) <- "Rel"
    hop_2$Cum <- NA
    
    for (i in 1:length(hop_2$Rel)){
      hop_sub <- hop_2[1:i,]
      
      if (length(hop_sub$Rel) == 1){
        hop_2$Cum[i] <- hop_sub$Rel
      } else {
        hop_2$Cum[i] <- sum(hop_sub$Rel)
      }
    }
    
    group_70 <- row.names(hop_2)[1:(length(row.names(hop_2)[hop_2$Cum < 0.7]) +1)]
    
    if (length(group_70) != 0){
      list_70[[uno]] <- group_70
      names(list_70[[uno]]) <- paste(plant,syncom, sep ="_")
    }
    
    uno <- uno + 1
  }
}

# Function to pad entries to the maximum length
pad_to_max_length <- function(x, max_len) {
  length(x) <- max_len
  return(x)
}

max_length_70 <- max(sapply(list_70, length))
df_70 <- do.call(rbind, lapply(list_70, pad_to_max_length, max_len = max_length_70))
df_70 <- as.data.frame(df_70)

write.table(df_70, paste(working_directory, "top70_isolates_no_dom.txt", sep = ""), sep = "\t", quote =F, row.names =F, col.names =F)
