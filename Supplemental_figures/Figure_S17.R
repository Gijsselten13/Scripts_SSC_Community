library("ggplot2") #Version 3.4.2
library("ggpubr") #Version 0.6.0
library("see") #Version 0.11.0

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Figure S17 - DESeq2 bar plots - Plant comparison =====
#Figure S17A - data with dominators
#Loading in tables
input_table <- read.table(paste(working_directory,"DESeq2/Sig_KO_all.txt", sep = ""), header=T, sep="\t")

#Make everything in one table
input_table$combination <- paste(input_table$Plant, input_table$SynCom, sep = "_")
input_table_1 <- input_table[input_table$log2FoldChange >0,]
input_table_2 <- input_table_1[input_table_1$padj < 0.05,]

vector_plant <- c("Arabidopsis", "Barley", "Lotus")
vector_sc <- c("AtSC", "LjSC", "HvSC", "SSC")

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
    value <- input_table_sub$Plant[input_table_sub$KO == paste(ko)]
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

bar_2 <- ggplot(new_6, aes(fill=Plant, y=Frequency, x=SynCom)) + 
  geom_bar(stat="identity") +  ggtitle("KO specificity in each SynCom") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_classic() +  scale_fill_okabeito(labels = replace_underscore_with_intersection)+
  geom_text(aes(label = Frequency), position = position_stack(vjust = .4), size = 5) +  # Add text labels
  labs(x ="", y = "Number of KOs", fill = "Significant in Hosts") + legend_theme_2+
  guides(fill = guide_legend(ncol = 1, title="Significant KOs in:"))  # Change nrow to desired number of rows
bar_2

#Data No Dominances - figure S17b
input_table <- read.table(paste(working_directory,"DESeq2/Sig_KO_all_no_nod_rhizo.txt", sep = ""), header=T, sep="\t")

#Make everything in one table
input_table$combination <- paste(input_table$Plant, input_table$SynCom, sep = "_")
input_table_1 <- input_table[input_table$log2FoldChange >0,]
input_table_2 <- input_table_1[input_table_1$padj < 0.05,]

vector_plant <- c("Arabidopsis", "Barley", "Lotus")
vector_sc <- c("AtSC", "LjSC", "HvSC", "SSC")

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

bar_4_no_dom <- ggplot(new_6, aes(fill=Plant, y=Frequency, x=SynCom)) + 
  geom_bar(stat="identity") +  ggtitle("KO specificity - excl dominators") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_classic() +  scale_fill_okabeito(labels = replace_underscore_with_intersection)+
  geom_text(aes(label = Frequency), position = position_stack(vjust = .4), size = 5) +  # Add text labels
  labs(x ="", y = "Number of KOs", fill = "Significant in Hosts") + legend_theme_2+
  guides(fill = guide_legend(ncol = 1, title="Significant KOs in:"))  # Change nrow to desired number of rows
bar_4_no_dom

bar_plot_sup_fig <- ggarrange(bar_2,bar_4_no_dom,ncol = 1, nrow = 2, labels = c("a", "b"),font.label = list(size = 20))

pdf(paste(results.dir,"Figure_S17_DESeq2_Bars_plant_comp.pdf", sep=""), width=10, height=12)
print(bar_plot_sup_fig) 
dev.off()
