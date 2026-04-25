library("ggplot2") #Version 3.4.2
library("ggtern") #Version 3.4.1

working_directory <- ""
working_directory_SA = paste(working_directory, "Shiny_app/", sep = "")
dir.create(paste(working_directory_SA, "results", sep = ""))
results.dir_SA <- paste(working_directory_SA,"results/", sep = "")

###Similar to figure S19-S23 - Ternary plots - per gene and SynCom =====
table <- read.table(paste(working_directory_SA,"ternary_KOs_av_med_all.txt", sep = ""), header =T, sep ="\t", row.names = 1)
SynComs <- c("AtSC", "LjSC", "HvSC", "SSC")

for (gene in unique(table$KO)){
  
  table_2 <- table[table$KO == paste(gene),]
  
  table_gene_2 <-data.frame()
  for (syncom in SynComs){
    table_sub_2 <- table_2[table_2$SynCom == paste(syncom),]
    At_val <- sum(table_sub_2$Arabidopsis)/length(table_sub_2$Arabidopsis)
    Hv_val <- sum(table_sub_2$Barley)/length(table_sub_2$Barley)
    Lj_val <- sum(table_sub_2$Lotus)/length(table_sub_2$Lotus)
    Prop_val <- sum(table_sub_2$Proportion_of_strains)/length(table_sub_2$Proportion_of_strains)
    
    table_gene <- t(data.frame(c(paste(gene), At_val, Hv_val, Lj_val, Prop_val, paste(syncom))))
    table_gene_2 <- rbind(table_gene_2, table_gene)
  }
  
  row.names(table_gene_2) <- NULL
  colnames(table_gene_2) <- c("gene", "Arabidopsis", "Barley", "Lotus", "Proportion_of_strains", "SynCom")
  
  table_gene_2$Arabidopsis <- as.numeric(table_gene_2$Arabidopsis)
  table_gene_2$Barley <- as.numeric(table_gene_2$Barley)
  table_gene_2$Lotus <- as.numeric(table_gene_2$Lotus)
  table_gene_2$Proportion_of_strains <- as.numeric(table_gene_2$Proportion_of_strains)
  nv = 0.005
  pn = position_nudge_tern(y=nv,x=-nv/2,z=-nv/2)
  
  ternary <- ggtern(data=table_gene_2,aes(x=Arabidopsis,y=Barley, z=Lotus, color = SynCom)) +
    geom_point(size = 6) +
    theme_bw()+
    scale_color_manual(values =c("#A3A500","#00B0F6","#00BF7D","#F8766D") ) +
    ggtitle(paste(gene)) + 
    theme(plot.title = element_text(hjust = 0.5, size = 20)) + 
    labs(color = "Gene") +
    theme(text = element_text(size=18)) + theme(legend.text=element_text(size=16)) +
    theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
  ternary 
  
  pdf(paste(results.dir_SA,"ternary_", gene,".pdf", sep=""), width=8, height=6)
  print(ternary)
  dev.off()
}
