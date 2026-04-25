working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Script to generate top90_isolates_no_dom.txt - Necessary for Figure S30 =====
SynComs <- c("AtSC", "HvSC", "LjSC", "SSC")
Plant <- c("At", "Hv", "Lj")

samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names

list_90 <- list()

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
    
    group_90 <- row.names(hop_2)[1:(length(row.names(hop_2)[hop_2$Cum < 0.9]) +1)]
    
    if (length(group_90) != 0){
      list_90[[uno]] <- group_90
      names(list_90[[uno]]) <- paste(plant,syncom, sep ="_")
    }
    
    uno <- uno + 1
  }
}

# Function to pad entries to the maximum length
pad_to_max_length <- function(x, max_len) {
  length(x) <- max_len
  return(x)
}

max_length_90 <- max(sapply(list_90, length))
df_90 <- do.call(rbind, lapply(list_90, pad_to_max_length, max_len = max_length_90))
df_90 <- as.data.frame(df_90)

write.table(df_90, paste(working_directory, "top90_isolates_no_dom.txt", sep = ""), sep = "\t", quote =F, row.names =F, col.names =F)
