library("ggplot2") #Version 3.4.2
library("ggpubr") #Version 0.6.0

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Figure 5a - DESeq2 bar plots - SynCom overlap =====

#First the DESeq2 results of the original dataset
input_table <- read.table(paste(working_directory,"DESeq2/Sig_KO_all.txt", sep = ""), header=T, sep="\t")

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
  input_table_sub <- input_table_2[input_table_2$Plant == paste(plant),]
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

bar <- ggplot(new_4, aes(fill=No_of_SynComs, y=Frequency, x=Plant)) + 
  geom_bar(stat="identity") +  ggtitle("Functional overlap in SynComs") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_classic() +
  scale_fill_manual(values = c("grey95", "grey75", "grey55","grey40")) +  # Manual scale for grays
  geom_text(aes(label = Frequency), position = position_stack(vjust = 0.8), size = 5) +  # Add text labels
  labs(x ="", y = "Number of KOs", fill = "No of Inocula") + legend_theme+
  guides(fill = guide_legend(nrow=1, title="Number of Inocula"))
bar

#Dataset without dominances
input_table <- read.table(paste(working_directory,"DESeq2/Sig_KO_all_no_nod_rhizo.txt", sep = ""), header=T, sep="\t")

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

bar_3_no_dom <- ggplot(new_4, aes(fill=No_of_SynComs, y=Frequency, x=Plant)) + 
  geom_bar(stat="identity") +  ggtitle("Functional overlap - excl dominators") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_classic()+
  scale_fill_manual(values = c("grey95", "grey75", "grey55","grey40")) +  # Manual scale for grays
  geom_text(aes(label = Frequency), position = position_stack(vjust = 0.8), size = 5) +  # Add text labels
  labs(x ="", y = "Number of KOs", fill = "No of Inocula") + legend_theme+
  guides(fill = guide_legend(nrow=1, title="Number of Inocula"))
bar_3_no_dom

bar_plot_main_fig <- ggarrange(bar,bar_3_no_dom,ncol = 2, nrow = 1, labels = c("a", "b"),font.label = list(size = 20))

pdf(paste(results.dir,"Figure_5a_Deseq2_bars.pdf", sep=""), width=12, height=6)
print(bar_plot_main_fig) 
dev.off()
