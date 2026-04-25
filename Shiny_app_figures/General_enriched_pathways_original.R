library("ggplot2") #Version 3.4.2
library("ggpubr") #Version 0.6.0
library("ggrepel") #Version 0.9.3

working_directory <- ""
working_directory_SA = paste(working_directory, "Shiny_app/", sep = "")
dir.create(paste(working_directory_SA, "results", sep = ""))
results.dir_SA <- paste(working_directory_SA,"results/", sep = "")

###Similar to figure S27b - General enriched pathways =====

table <- read.table(paste(working_directory_SA, "boxplots_full.txt", sep = ""), header = T, sep = "\t", row.names =1)

colnames(table) <- c("Plant","RA", "Have", "No_Have", "SynCom", "pathway", "Category_2", "No_of_strains")
table$Category <- paste(table$pathway, " (n = ",format(round(as.numeric(table$No_of_strains),2), nsmall = 2), ")", sep = "")

table$Log <- log10(table$Have)

table_2 <- table[,colnames(table) != "No_Have"]
table_3 <- table[,colnames(table) != "Have"]

colnames(table_2) <- c("Plant","RA", "value", "SynCom", "pathway", "Category_2","No_of_strains","Category")
colnames(table_3) <- c("Plant", "RA", "value", "SynCom", "pathway", "Category_2","No_of_strains","Category")

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
  
  RA <- sum(table_sub$RA[table_sub$Present == "Present"])/length(table_sub$RA[table_sub$Present == "Present"])
  
  
  new <- t(data.frame(c(Category_2, description, round(as.numeric(No_of_strains),2),RA, Yes_RA_change, No_RA_change, stat$p.format)))
  new_2 <- rbind(new_2, new)
}

row.names(new_2) <- NULL
colnames(new_2) <- c("Category", "Pathway", "No_of_strains", "RA", "RA_change", "RA change w/o pathway", "Adjusted_p_value")
new_2$Adjusted_p_value <- as.numeric(new_2$Adjusted_p_value)

new_3 <- new_2[order(new_2$Adjusted_p_value, decreasing =F),]
new_3$RA_change <- round(as.numeric(new_3$RA_change), 2)
new_3$`RA change w/o pathway` <- round(as.numeric(new_3$`RA change w/o pathway`), 2)

new_4 <- new_3[new_3$RA_change > new_3$`RA change w/o pathway`,]
new_5 <- new_4[new_4$Adjusted_p_value < 0.05,]
new_5$No_of_strains <- as.numeric(new_5$No_of_strains)
new_5$RA <- as.numeric(new_5$RA)

Encoding(new_5$Pathway) <- 'latin1'

new_5$Square <- sqrt(((new_5$RA/max(new_5$RA)) * (new_5$RA/max(new_5$RA))) + ((new_5$RA_change/max(new_5$RA_change)) * (new_5$RA_change/max(new_5$RA_change))))
pathway_selection <- new_5$Pathway[order(new_5$Square, decreasing = T)][1:9]
pathway_selection_2 <- new_5$Pathway[order(new_5$RA_change, decreasing = T)][3:9]
pathway_selection <- c(pathway_selection, pathway_selection_2)

plot_cor <- ggscatter(new_5, x = "RA", y = "RA_change", color = "Category", size  = "No_of_strains") + 
  ggtitle("General plant-selected pathways") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  ylab("Fold change vs Input") + 
  xlab("Relative abundance") +
  theme(legend.position = "right") +
  guides(color=guide_legend(ncol =1)) +
  geom_text_repel(aes(label=ifelse(new_5$RA > 0.1 & new_5$RA_change > 2 & !new_5$Pathway %in% pathway_selection,as.character(new_5$Pathway), '')),size=5,max.overlaps = 5) +
  geom_text_repel(aes(label=ifelse(new_5$RA_change > 6.5 & !new_5$Pathway %in% pathway_selection,as.character(new_5$Pathway), '')),size=5,max.overlaps = 3) +
  geom_text_repel(aes(label=ifelse(new_5$Pathway %in% pathway_selection,as.character(new_5$Pathway), '')),size=5,max.overlaps = Inf,fontface = "bold") +
  labs(color = "Category", size = "Proportion of strains") +
  theme(axis.text.x = element_text(size = 14), axis.title = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=24))
plot_cor

pdf(paste(results.dir_SA,"CRAC_vs_RA_KO.pdf", sep=""), width=16, height=8)
print(plot_cor_KO)
dev.off()
