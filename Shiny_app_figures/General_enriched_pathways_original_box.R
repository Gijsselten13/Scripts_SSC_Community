library("ggplot2") #Version 3.4.2
library("ggpubr") #Version 0.6.0
library("ggrepel") #Version 0.9.3
library("scales") #Version 1.2.1

working_directory <- ""
working_directory_SA = paste(working_directory, "Shiny_app/", sep = "")
dir.create(paste(working_directory_SA, "results", sep = ""))
results.dir_SA <- paste(working_directory_SA,"results/", sep = "")

###Similar to figure 5f - General enriched pathways - Box plots strains with vs strains without =====
table <- read.table(paste(working_directory_SA,"boxplots_full.txt", sep = ""),header = T, sep = "\t", row.names =1)

colnames(table) <- c("Plant", "RA", "Have", "No_Have", "SynCom", "pathway", "Category_2", "No_of_strains")
table$Category <- paste(table$pathway, " (n = ",format(round(as.numeric(table$No_of_strains),2), nsmall = 2), ")", sep = "")

table_2 <- table[,colnames(table) != "No_Have"]
table_3 <- table[,colnames(table) != "Have"]

colnames(table_2) <- c("Plant", "RA", "value","SynCom", "pathway", "Category_2","No_of_strains", "Category")
colnames(table_3) <- c("Plant", "RA", "value","SynCom", "pathway", "Category_2","No_of_strains", "Category")

table_2$Present <-"Present"
table_3$Present <-"Absent"

table_4 <- rbind(table_2, table_3)

table_4$Present <- factor(table_4$Present, levels = c("Present", "Absent"))
table_5 <- na.omit(table_4)

list_of_cats <- unique(table_5$Category)
new_2 <- data.frame()

for (cat in list_of_cats){
  table_sub <- table_5[table_5$Category == paste(cat),]
  table_sub$value <- as.numeric(table_sub$value)
  
  stat <- compare_means(value~Present, table_sub, method = "t.test")
  
  Yes_RA_change <- sum(table_sub$value[table_sub$Present == "Present"])/length(table_sub$value[table_sub$Present == "Present"])
  No_RA_change <- sum(table_sub$value[table_sub$Present == "Absent"])/length(table_sub$value[table_sub$Present == "Absent"])
  
  Yes_RA_change <- sum(table_sub$value[table_sub$Present == "Present"])/length(table_sub$value[table_sub$Present == "Present"])
  No_RA_change <- sum(table_sub$value[table_sub$Present == "Absent"])/length(table_sub$value[table_sub$Present == "Absent"])
  
  Category_2 <- unique(table_sub$Category_2)
  description <- unique(table_sub$pathway)
  No_of_strains <- unique(table_sub$No_of_strains)
  RA <- sum(table_sub$RA)/length(table_sub$RA)
  
  new <- t(data.frame(c(Category_2, description, round(as.numeric(No_of_strains),2), RA,Yes_RA_change, No_RA_change, stat$p.format)))
  new_2 <- rbind(new_2, new)
}

row.names(new_2) <- NULL
colnames(new_2) <- c("Category", "Pathway", "No_of_strains", "RA","RA change", "RA change w/o pathway", "Adjusted_p_value")
new_2$Adjusted_p_value <- as.numeric(new_2$Adjusted_p_value)

new_3 <- new_2[order(new_2$Adjusted_p_value, decreasing =F),]
new_3$`RA change` <- round(as.numeric(new_3$`RA change`), 2)
new_3$`RA change w/o pathway` <- round(as.numeric(new_3$`RA change w/o pathway`), 2)

new_4 <- new_3[new_3$`RA change` > new_3$`RA change w/o pathway`,]
new_5 <- new_4[new_4$Adjusted_p_value < 0.05,]

#Here pathways of interest can be put, e.g: 
pathways_of_interest <- c("Amino sugar and nucleotide sugar metabolism","Inositol phosphate metabolism", "Biofilm formation - Escherichia coli","Biosynthesis of various antibiotics","Quorum sensing", "Flagellar assembly", "Sulfur metabolism", "Cysteine and methionine metabolism", "Starch and sucrose metabolism", "ABC transporters", "Transcriptional regulator", "Biosynthesis of enediyne antibiotics")

table_6 <- table_5[table_5$pathway %in% pathways_of_interest,]
table_7 <- table_6[table_6$Present == "Present",] 
table_8 <- table_6[table_6$Present == "Absent",] 
table_7$Present <- "Present"
table_8$Present <- "Absent"

table_9 <- rbind(table_7, table_8)

hex <- hue_pal()(length(unique(table_9$pathway))) 

swr = function(string, nwrap=23) {
  paste(strwrap(string, width=nwrap), collapse="\n")
}
swr = Vectorize(swr)

# Create line breaks in Year
table_9$Category = swr(table_9$Category)

plot_top_sig <- ggplot(table_9, aes(x=Present , y=value, color = pathway)) + 
  geom_boxplot(outlier.shape = NA)+
  ggtitle("Pathways") + 
  geom_jitter(aes(shape=Plant,size = 3), cex =3)+
  theme(legend.position="top")+ 
  scale_color_manual(values = hex) +
  guides(color = FALSE) +
  guides(size = FALSE) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(size =12), axis.title = element_text(size = 14), axis.text.y = element_text(size=12), legend.title = element_text(size=16), legend.text = element_text(size=14), plot.title = element_text(size=20)) +  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  ylab("Fold change vs Input") + 
  xlab("") +
  facet_wrap(~Category,scales = "free_y") +
  theme(strip.text = element_text(size=12))
plot_top_sig_2 <- plot_top_sig + stat_compare_means(comparisons = list(c("Absent", "Present")), method = "t.test", vjust = 1.5) +
  guides(shape = guide_legend(override.aes = list(size = 5)))
plot_top_sig_2

pdf(paste(results.dir_SA,"Box_plots.pdf", sep=""), width=16, height=8)
print(plot_top_sig_2)
dev.off()
