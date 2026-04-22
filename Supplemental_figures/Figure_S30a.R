library("reshape2") #Version 1.4.4
library("ggplot2") #Version 3.4.2
library("ggpubr") #Version 0.6.0

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Figure S30a - Overview of 266 across families =====
df_90 <- read.table(paste(working_directory, "top90_isolates_no_dom.txt", sep = ""), sep = "\t",header =F)
df_90_2 <- unique(na.omit(unlist(as.vector(df_90))))

tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE, row.names =1)
top <- tax_df[row.names(tax_df) %in% df_90_2,]
non_top <- tax_df[!row.names(tax_df) %in% df_90_2,]

input_table <- read.table(paste(working_directory, "DESeq2/Sig_KO_all_no_nod_rhizo.txt", sep = ""), header=T, sep="\t")
input_table_2 <- table(input_table$KO)
input_table_3 <- names(input_table_2)[input_table_2 == 12]

table <- read.table(paste(working_directory, "KO_genome/KO_SSC.tsv", sep = ""), sep= "\t", header =T,row.names =1) 
colnames(table)[grep("M.16",colnames(table))] <- "M-16"
colnames(table)[grep("M.6",colnames(table))] <- "M-6"
colnames(table)[grep("M.10",colnames(table))] <- "M-10"
colnames(table)[grep("M.11_2",colnames(table))] <- "M-11_2"
colnames(table)[grep("^M.11$",colnames(table))] <- "M-11"
colnames(table) <- gsub("X", "", colnames(table))
table[table > 0] <- 1

table_2 <- table[row.names(table) %in% input_table_3,]

Tax_Fam_top <- table(top$family)
Tax_Fam_top_2 <- names(Tax_Fam_top)[Tax_Fam_top > 2]

top_iso <- top[top$family %in% Tax_Fam_top_2,]
No_top_iso <- non_top[non_top$family %in% Tax_Fam_top_2,]

table_3 <- data.frame()

for (family in Tax_Fam_top_2 ){
  top_iso_sub <- row.names(top_iso)[top_iso$family == paste(family)]
  table_sub <- table_2[,colnames(table_2) %in% top_iso_sub]
  table_sub_2 <- rowSums(table_sub)/length(colnames(table_sub))
  
  value_100 <- length(names(table_sub_2)[table_sub_2 >= 0.9])
  value_90 <- length(names(table_sub_2)[table_sub_2 >= 0.6 & table_sub_2 < 0.9])
  value_60 <- length(names(table_sub_2)[table_sub_2 >= 0.3 & table_sub_2 < 0.6])
  value_30 <- length(names(table_sub_2)[table_sub_2 > 0 & table_sub_2 < 0.3])
  value_0 <- length(names(table_sub_2)[table_sub_2 == 0])
  
  table_top <- data.frame(t(data.frame(c(paste(family), length(top_iso_sub), "Top", value_100, value_90, value_60, value_30, value_0))))
  
  No_top_iso_sub <- row.names(No_top_iso)[No_top_iso$family == paste(family)]
  table_sub <- table_2[,colnames(table_2) %in% No_top_iso_sub]
  table_sub_2 <- rowSums(table_sub)/length(colnames(table_sub))
  
  value_100 <- length(names(table_sub_2)[table_sub_2 >= 0.9])
  value_90 <- length(names(table_sub_2)[table_sub_2 >= 0.6 & table_sub_2 < 0.9])
  value_60 <- length(names(table_sub_2)[table_sub_2 >= 0.3 & table_sub_2 < 0.6])
  value_30 <- length(names(table_sub_2)[table_sub_2 > 0 & table_sub_2 < 0.3])
  value_0 <- length(names(table_sub_2)[table_sub_2 == 0])
  
  table_non_top <- data.frame(t(data.frame(c(paste(family),length(No_top_iso_sub), "No top", value_100, value_90, value_60, value_30, value_0))))
  
  table_3 <- rbind(table_3,table_top,table_non_top)
}

row.names(table_3) <- paste(table_3$X1, " (n = ",table_3$X2,") - ", table_3$X3, sep ="")
table_4 <- table_3[,c(4,5,6,7,8)]
colnames(table_4) <- c("100%-90%", "90%-60%","60%-30%","30%-1%", "0%")
table_4$`100%-90%` <- as.numeric(table_4$`100%-90%`)
table_4$`90%-60%` <- as.numeric(table_4$`90%-60%`)
table_4$`60%-30%` <- as.numeric(table_4$`60%-30%`)
table_4$`30%-1%` <- as.numeric(table_4$`30%-1%`)
table_4$`0%` <- as.numeric(table_4$`0%`)

table_4$Group <- row.names(table_4)
table_5 <- melt(table_4)

table_5$variable <- factor(table_5$variable, levels = c("0%","30%-1%","60%-30%","90%-60%","100%-90%"))

table_5$Group <- factor(table_5$Group, levels =c("Rhizobiaceae (n = 54) - Top", "Rhizobiaceae (n = 61) - No top", 
                                                 "Xanthobacteraceae (n = 12) - Top" ,"Xanthobacteraceae (n = 3) - No top",
                                                 "Beijerinckiaceae (n = 20) - Top", "Beijerinckiaceae (n = 26) - No top",
                                                 "Caulobacteraceae (n = 26) - Top", "Caulobacteraceae (n = 14) - No top",
                                                 "Cellulomonadaceae (n = 3) - Top","Cellulomonadaceae (n = 11) - No top", 
                                                 "Micrococcaceae (n = 8) - Top", "Micrococcaceae (n = 74) - No top",
                                                 "Sphingomonadaceae (n = 11) - Top", "Sphingomonadaceae (n = 3) - No top",
                                                 "Pseudomonadaceae (n = 42) - Top", "Pseudomonadaceae (n = 41) - No top",
                                                 "Microbacteriaceae (n = 6) - Top",  "Microbacteriaceae (n = 23) - No top",
                                                 "Sphingobacteriaceae (n = 11) - Top", "Sphingobacteriaceae (n = 21) - No top",
                                                 "Burkholderiaceae (n = 138) - Top","Burkholderiaceae (n = 129) - No top", 
                                                 "Rhodanobacteraceae (n = 6) - Top","Rhodanobacteraceae (n = 9) - No top",
                                                 "Xanthomonadaceae (n = 24) - Top", "Xanthomonadaceae (n = 40) - No top"))

g1 <- ggplot(table_5, aes(x= Group, weight=value, fill=variable)) +
  theme_classic() +
  geom_bar(position = "stack") + 
  scale_fill_manual(values = c("black","#d3d3d3","#fa8072","#ffd700","#66bd63")) +
  ggtitle("266 KO distribution in 90% CRA top vs non-top") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  ylab("No of KOs") + 
  xlab("Family - top vs non-top") +
  labs(fill = "Proportion of KO in Family") +
  theme(axis.text.x = element_text(size = 14, hjust =1, angle = 45), axis.title.y = element_text(size = 18),axis.title.x = element_blank(), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=24))
g1 

g2 <- ggarrange(NA, g1, ncol = 2, widths = c(1,6))

pdf(paste(results.dir, "Figure_S30a_overview.pdf", sep=""), width=20, height=7)
print(g2)
dev.off()
