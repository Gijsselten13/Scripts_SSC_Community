library("dplyr") #Version 1.1.2

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Frequency isolates per Family across SynComs =====
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")

top <- read.table(paste(working_directory, "top70_isolates_no_dom.txt", sep = ""), sep = "\t", header =F)
top_2 <- unlist(as.vector(top))
top_3 <- top_2[!is.na(top_2)]
top_4 <- unique(top_3)

fams_left <- unique(tax_df_2$Family[row.names(tax_df_2) %in% top_4])

SynComs <- c("AtSC","HvSC", "LjSC")

all_together <- data.frame()

for(syncom in SynComs){
  tax_df_2_SC <- tax_df_2[tax_df_2$SynCom == paste(syncom),]
  
  new <- data.frame(table(tax_df_2_SC$Family))
  new$SynCom <- paste(syncom)
  
  all_together <- rbind(all_together, new)
}

all_together_2 <- all_together[all_together$Var1 %in% fams_left,]

family_order <- c("Chitinophagaceae", "Microbacteriaceae","Micrococcaceae","Xanthobacteraceae","Sphingobacteriaceae","Rhodanobacteraceae","Sphingomonadaceae","Flavobacteriaceae","Devosiaceae","Beijerinckiaceae","Enterobacteriaceae","Caulobacteraceae","Pseudomonadaceae","Xanthomonadaceae","Burkholderiaceae","Rhizobiaceae")

all_together_2$SynCom <- factor(all_together_2$SynCom, levels = c("LjSC","HvSC","AtSC"))
all_together_2$Var1 <- factor(all_together_2$Var1, levels = family_order)

write.table(all_together_2, paste(working_directory, "Family_R2/No_of_isolates_per_fam.txt", sep = ""), col.names =T, row.names =F, quote =F, sep = "\t")
