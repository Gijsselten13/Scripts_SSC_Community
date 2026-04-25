library("ggplot2") #Version 3.4.2
library("dplyr") #Version 1.1.2
library("ggpp") #Version 0.5.5

working_directory <- ""
working_directory_SA = paste(working_directory, "Shiny_app/", sep = "")
dir.create(paste(working_directory_SA, "results", sep = ""))
results.dir_SA <- paste(working_directory_SA,"results/", sep = "")

###Similar to figure S19-S23 - Gene abundance bar chart =====
table = read.table(paste(working_directory_SA,"Abundances_full.tsv", sep = ""), header =T, sep = "\t")
genes = table$Gene

for (gene in genes){
  table_sub <- table [table$Gene == paste(gene),]
  table_sub$RA_KO <- as.numeric(table_sub$RA_KO)
  
  table_sub$New_column <- paste(table_sub$SynCom, table_sub$Origin, sep = "-")
  table_sub$New_column <- gsub("HvSC-HvSC", "HvSC", table_sub$New_column)
  table_sub$New_column <- gsub("AtSC-AtSC", "AtSC", table_sub$New_column)
  table_sub$New_column <- gsub("LjSC-LjSC", "LjSC", table_sub$New_column)
  table_sub$New_column <- factor(table_sub$New_column, levels = c("AtSC", "HvSC", "LjSC", "SSC-AtSC", "SSC-HvSC", "SSC-LjSC"))
  
  g1 <- ggplot(table_sub %>% filter(SynCom != "SSC"),
               aes(x= Plant, weight=RA_KO, fill=New_column)) +
    theme_classic() +
    geom_bar(position = "dodge", width=0.5, just = 0.5) + 
    labs(x="Plant") +
    scale_fill_manual(values = c("#A3A500","#00B0F6","#00BF7D")) +
    ggtitle(paste(gene)) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    ylab("Relative abundance - Gene") + 
    xlab("Plant") +
    labs(fill = "SynCom") +
    theme(axis.text.x = element_text(size = 14, hjust =0.2), axis.title.y = element_text(size = 18),axis.title.x = element_blank(), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=24))
  g1
  
  g2 <- g1 + geom_bar(data=table_sub %>% filter(SynCom == "SSC"),
                      aes(x=Plant, fill=New_column),
                      position=position_stacknudge(x = 0.335), width=0.17) +
    scale_fill_manual(values = c("#A3A500","#00B0F6","#00BF7D","#fcddd9", "#fabab3","#F8766D"))
  
  pdf(paste(results.dir_SA,"RA_", gene,".pdf", sep=""), width=8, height=4)
  print(g2)
  dev.off()
}
