library("ggplot2") #Version 3.4.2
library("ggpubr") #Version 0.6.0
library("see") #Version 0.11.0

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Figure S18 - DESeq2 bar plots - No nodulators or no Rhizobacter =====
#Data no nodulators
input_table <- read.table(paste(working_directory,"DESeq2/Sig_KO_all_no_nod.txt", sep = ""), header=T, sep="\t")

#Make everything in one table
input_table$combination <- paste(input_table$Plant, input_table$SynCom, sep = "_")
input_table_1 <- input_table[input_table$log2FoldChange >0,]
input_table_2 <- input_table_1[input_table_1$padj < 0.05,]

vector_plant <- c("Arabidopsis", "Barley", "Lotus")
vector_sc <- c("AtSC", "LjSC", "HvSC", "SSC")

new_4 <- data.frame(matrix(NA, ncol = 2))
colnames(new_4) <- c("No_of_SynComs", "Frequency")
new_4 <- new_4[-1,]

for (plant in vector_plant){
  input_table_sub <- input_table_2[input_table_2$plant == paste(plant),]
  new_vector <- unique(input_table_sub$KO)
  
  new_2 <- data.frame(matrix(NA, ncol = 2))
  colnames(new_2) <- c("KO", "No_of_sig_occurences")
  new_2 <- new_2[-1,]
  
  for (ko in new_vector){
    value <- length(input_table_sub$KO[input_table_sub$KO == paste(ko)])
    new <- t(as.data.frame(c(paste(ko), value)))
    new_2 <- rbind(new_2, new)
  }
  colnames(new_2) <- c("KO", "No_of_sig_occurences")
  new_3 <- data.frame(table(new_2$No_of_sig_occurences))
  colnames(new_3) <- c("No_of_SynComs", "Frequency")
  new_3$Plant <- paste(plant)
  
  new_4 <- rbind(new_4,new_3)
}

legend_theme <- theme(legend.text = element_text(size = 16),  # Apply text size
                      legend.key.size = unit(0.5, "cm"), # Apply key size
                      strip.background=element_rect(colour="gray50", size=0.3), # Change 'size' for thickness
                      axis.text=element_text(color="gray50"),
                      axis.line = element_line(color="gray50", size=0.3),
                      axis.line.x = element_line(color="gray50", size=0.3), 
                      axis.line.y = element_line(color="gray50", size=0.3),
                      axis.ticks.x = element_line(color="gray50", size=0.3),
                      axis.ticks.y =element_line(color="gray50", size=0.3),
                      axis.text.x = element_text(size = 16), 
                      axis.title = element_text(size = 16), 
                      axis.text.y = element_text(size = 16), 
                      legend.title = element_text(size = 16), 
                      legend.position = "top")

legend_theme_2 <- theme(legend.text = element_text(size = 16),  # Apply text size
                        legend.key.size = unit(0.5, "cm"), # Apply key size
                        strip.background=element_rect(colour="gray50", size=0.3), # Change 'size' for thickness
                        axis.text=element_text(color="gray50"),
                        axis.line = element_line(color="gray50", size=0.3),
                        axis.line.x = element_line(color="gray50", size=0.3), 
                        axis.line.y = element_line(color="gray50", size=0.3),
                        axis.ticks.x = element_line(color="gray50", size=0.3),
                        axis.ticks.y =element_line(color="gray50", size=0.3),
                        axis.text.x = element_text(size = 16), 
                        axis.title = element_text(size = 16), 
                        axis.text.y = element_text(size = 16), 
                        legend.title = element_text(size = 16), 
                        legend.position = "right"
)

# Function to replace underscores with intersection symbol
replace_underscore_with_intersection <- function(labels) {
  sapply(labels, function(label) gsub("_", " ∩ ", label, fixed = TRUE))
}

bar_no_nod <- ggplot(new_4, aes(fill=No_of_SynComs, y=Frequency, x=Plant)) + 
  geom_bar(stat="identity") +  ggtitle("Functional overlap - excl nodulators") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_classic()+
  scale_fill_manual(values = c("grey95", "grey75", "grey55","grey40")) +  # Manual scale for grays
  geom_text(aes(label = Frequency), position = position_stack(vjust = 0.8), size = 5) +  # Add text labels
  labs(x ="", y = "Number of KOs", fill = "No of Inocula") + legend_theme+
  guides(fill = guide_legend(nrow=1, title="Number of Inocula"))
bar_no_nod

new_6 <- data.frame(matrix(NA, ncol = 3))
colnames(new_6) <- c("Plant", "Frequency", "SynCom")
new_6 <- new_6[-1,]

for (syncom in vector_sc){
  input_table_sub <- input_table_2[input_table_2$SynCom == paste(syncom),]
  new_vector <- unique(input_table_sub$KO)
  
  new_7 <- data.frame(matrix(NA, ncol = 2))
  colnames(new_7) <- c("KO", "Plant")
  new_7 <- new_7[-1,]
  
  for (ko in new_vector){
    value <- input_table_sub$plant[input_table_sub$KO == paste(ko)]
    if (length(value) > 1){
      value_2 <- paste(value, collapse = "_")
      value <- value_2
    }
    new_8 <- t(as.data.frame(c(paste(ko), value)))
    new_7 <- rbind(new_7, new_8)
  }
  colnames(new_7) <- c("KO", "Plant")
  new_9 <- data.frame(table(new_7$Plant))
  colnames(new_9) <- c("Plant", "Frequency")
  new_9$SynCom <- paste(syncom)
  
  new_6 <- rbind(new_6,new_9)
}

new_6$Plant <- factor(new_6$Plant, levels = c("Arabidopsis", "Barley", "Lotus", "Barley_Arabidopsis", "Lotus_Arabidopsis", "Lotus_Barley", "Lotus_Barley_Arabidopsis"))

bar_2_no_nod <- ggplot(new_6, aes(fill=Plant, y=Frequency, x=SynCom)) + 
  geom_bar(stat="identity") +  ggtitle("KO specificity - excl nodulators") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_classic() +  scale_fill_okabeito(labels = replace_underscore_with_intersection)+
  geom_text(aes(label = Frequency), position = position_stack(vjust = .4), size = 5) +  # Add text labels
  labs(x ="", y = "Number of KOs", fill = "Significant in Hosts") + legend_theme_2+
  guides(fill = guide_legend(ncol = 1, title="Significant KOs in:"))  # Change nrow to desired number of rows
bar_2_no_nod

#Figure S18c and d - Data No Rhizobacter
input_table <- read.table(paste(working_directory,"DESeq2/Sig_KO_all_no_rhizobacter.txt", sep = ""), header=T, sep="\t")

#Make everything in one table
input_table$combination <- paste(input_table$Plant, input_table$SynCom, sep = "_")
input_table_1 <- input_table[input_table$log2FoldChange >0,]
input_table_2 <- input_table_1[input_table_1$padj < 0.05,]

vector_plant <- c("Arabidopsis", "Barley", "Lotus")
vector_sc <- c("AtSC", "LjSC", "HvSC", "SSC")

new_4 <- data.frame(matrix(NA, ncol = 2))
colnames(new_4) <- c("No_of_SynComs", "Frequency")
new_4 <- new_4[-1,]

for (plant in vector_plant){
  input_table_sub <- input_table_2[input_table_2$plant == paste(plant),]
  new_vector <- unique(input_table_sub$KO)
  
  new_2 <- data.frame(matrix(NA, ncol = 2))
  colnames(new_2) <- c("KO", "No_of_sig_occurences")
  new_2 <- new_2[-1,]
  
  for (ko in new_vector){
    value <- length(input_table_sub$KO[input_table_sub$KO == paste(ko)])
    new <- t(as.data.frame(c(paste(ko), value)))
    new_2 <- rbind(new_2, new)
  }
  colnames(new_2) <- c("KO", "No_of_sig_occurences")
  new_3 <- data.frame(table(new_2$No_of_sig_occurences))
  colnames(new_3) <- c("No_of_SynComs", "Frequency")
  new_3$Plant <- paste(plant)
  
  new_4 <- rbind(new_4,new_3)
}

bar_3_no_rhizo <- ggplot(new_4, aes(fill=No_of_SynComs, y=Frequency, x=Plant)) + 
  geom_bar(stat="identity") +  ggtitle("Functional overlap - excl Rhizobacter") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_classic()+
  scale_fill_manual(values = c("grey95", "grey75", "grey55","grey40")) +  # Manual scale for grays
  geom_text(aes(label = Frequency), position = position_stack(vjust = 0.8), size = 5) +  # Add text labels
  labs(x ="", y = "Number of KOs", fill = "No of Inocula") + legend_theme+
  guides(fill = guide_legend(nrow=1, title="Number of Inocula"))
bar_3_no_rhizo

new_6 <- data.frame(matrix(NA, ncol = 3))
colnames(new_6) <- c("Plant", "Frequency", "SynCom")
new_6 <- new_6[-1,]

for (syncom in vector_sc){
  input_table_sub <- input_table_2[input_table_2$SynCom == paste(syncom),]
  new_vector <- unique(input_table_sub$KO)
  
  new_7 <- data.frame(matrix(NA, ncol = 2))
  colnames(new_7) <- c("KO", "Plant")
  new_7 <- new_7[-1,]
  
  for (ko in new_vector){
    value <- input_table_sub$plant[input_table_sub$KO == paste(ko)]
    if (length(value) > 1){
      value_2 <- paste(value, collapse = "_")
      value <- value_2
    }
    new_8 <- t(as.data.frame(c(paste(ko), value)))
    new_7 <- rbind(new_7, new_8)
  }
  colnames(new_7) <- c("KO", "Plant")
  new_9 <- data.frame(table(new_7$Plant))
  colnames(new_9) <- c("Plant", "Frequency")
  new_9$SynCom <- paste(syncom)
  
  new_6 <- rbind(new_6,new_9)
}

new_6$Plant <- factor(new_6$Plant, levels = c("Arabidopsis", "Barley", "Lotus", "Barley_Arabidopsis", "Lotus_Arabidopsis", "Lotus_Barley", "Lotus_Barley_Arabidopsis"))

bar_4_no_rhizo <- ggplot(new_6, aes(fill=Plant, y=Frequency, x=SynCom)) + 
  geom_bar(stat="identity") +  ggtitle("KO specificity - excl Rhizobacter") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_classic() +  scale_fill_okabeito(labels = replace_underscore_with_intersection)+
  geom_text(aes(label = Frequency), position = position_stack(vjust = .4), size = 5) +  # Add text labels
  labs(x ="", y = "Number of KOs", fill = "Significant in Hosts") + legend_theme_2+
  guides(fill = guide_legend(ncol = 1, title="Significant KOs in:"))  # Change nrow to desired number of rows
bar_4_no_rhizo

bar_plot_sup_fig_nod <- ggarrange(bar_no_nod,bar_2_no_nod,ncol = 2, nrow = 1, labels = c("a", "b"),font.label = list(size = 20), widths = c(2,3))
bar_plot_sup_fig_rhizo <- ggarrange(bar_3_no_rhizo,bar_4_no_rhizo,ncol = 2, nrow = 1, labels = c("c", "d"),font.label = list(size = 20), widths = c(2,3))
bar_plot_sup_fig_S18 <- ggarrange(bar_plot_sup_fig_nod,bar_plot_sup_fig_rhizo,ncol = 1, nrow = 2)

pdf(paste(results.dir,"Figure_S18_DESeq2_bars_nod_rhizo.pdf", sep=""), width=12, height=12)
print(bar_plot_sup_fig_S18) 
dev.off()
