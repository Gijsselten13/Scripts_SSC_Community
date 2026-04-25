library("dplyr") #Version 1.1.2

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Frequency isolates per Genus across SynComs - Burkholderiaceae genera =====
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")

fams_left <- c("Acidovorax","Cupriavidus","Pelomonas","Polaromonas","Rhizobacter","Variovorax")

SynComs <- c("AtSC","HvSC", "LjSC")

all_together <- data.frame()

for(syncom in SynComs){
  tax_df_2_SC <- tax_df_2[tax_df_2$SynCom == paste(syncom),]
  
  new <- data.frame(table(tax_df_2_SC$Genus))
  new$SynCom <- paste(syncom)
  
  all_together <- rbind(all_together, new)
}

all_together_2 <- all_together[all_together$Var1 %in% fams_left,]

family_order <- c("Pelomonas","Cupriavidus","Polaromonas","Variovorax","Rhizobacter","Acidovorax")

all_together_2$SynCom <- factor(all_together_2$SynCom, levels = c("LjSC","HvSC","AtSC"))
all_together_2$Var1 <- factor(all_together_2$Var1, levels = family_order)

write.table(all_together_2, paste(working_directory, "Family_R2/No_of_isolates_per_Burkholderiaceae_genus.txt", sep = ""), col.names =T, row.names =F, quote =F, sep = "\t")
