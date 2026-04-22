library("ggplot2") #Version 3.4.2

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Figure 5d - bar plot & ABC transporter diversity =====
input_table <- read.table(paste(working_directory, "/DESeq2/Sig_KO_all.txt", sep = ""), header=T, sep="\t")

top <- read.table(paste(working_directory,"/Annotations/pathway_top.txt", sep = ""), header=F, sep="\t")
KO_to_pathway <- read.table(paste(working_directory,"/Annotations/KO_to_pathway.txt", sep = ""), header=T, sep="\t")
KO_to_pathway$V3 <- top$V2[match(KO_to_pathway$V2, top$V1)]

input_table$pathway <- KO_to_pathway$V3[match(input_table$KO, KO_to_pathway$V1)]
input_table_2 <- input_table[!is.na(input_table$pathway),]

input_table_3 <- input_table_2[input_table_2$padj < 0.05,]

Plant <- c("Arabidopsis","Barley",  "Lotus")
Plant_count <- data.frame()

for (plant in Plant){
  input_table_sub <- input_table_3[input_table_3$Plant == paste(plant),]
  data_count <- data.frame(table(input_table_sub$KO))
  data_count_KOs <- data_count$Var1[data_count$Freq == 4]
  input_table_sub_2 <- input_table_sub[input_table_sub$KO %in% data_count_KOs,]
  data_count_together <- data.frame(table(input_table_sub_2$pathway))
  data_count_together$Plant <- paste(plant)
  Plant_count <- rbind(Plant_count, data_count_together)
}

#Take top 10 per data
Plant_count_3 <- data.frame()

for (pathway in unique(Plant_count$Var1)){
  Plant_count_2 <- Plant_count[Plant_count$Var1 == paste(pathway),]
  count_value <- sum(Plant_count_2$Count)
  count_data <- data.frame(t(data.frame(c(paste(pathway), count_value))))
  Plant_count_3 <- rbind(Plant_count_3, count_data)
}

row.names(Plant_count_3) <- NULL
Plant_count_3$X2 <- as.numeric(Plant_count_3$X2)/4
Plant_count_4 <- Plant_count_3[order(Plant_count_3$X2, decreasing=TRUE),]
Plant_count_5 <- Plant_count_4$X1[1:10]

#plotting
colnames(Plant_count) <- c("Pathway", "Count", "Plant")

Plant_count_2 <- Plant_count[Plant_count$Pathway %in% Plant_count_5,]
Plant_count_2$Count <- as.numeric(Plant_count_2$Count)/4

for (pathway in unique(Plant_count_2$Pathway)){
  Plant_count_sub <- Plant_count_2[Plant_count_2$Pathway == paste(pathway),]
  if (length(Plant_count_sub != 3)){
    plant_missing <- Plant[!Plant %in% Plant_count_sub$Plant]
    for (plant_hop in plant_missing){
      new_data <- data.frame(t(data.frame(c(paste(pathway),0,paste(plant_hop)))))
      colnames(new_data) <- colnames(Plant_count_2)
      Plant_count_2 <- rbind(Plant_count_2, new_data)
    }
  }
}
row.names(Plant_count_2) <- NULL

Plant_count_2$Count <- as.numeric(Plant_count_2$Count)

# Calculate total counts for each pathway
total_counts <- aggregate(Count ~ Pathway, data = Plant_count_2, sum)
# Order the pathways by total counts in descending order
ordered_pathways <- total_counts$Pathway[order(total_counts$Count, decreasing = TRUE)]

# Adjust the factor levels of the Pathway variable before wrapping text
Plant_count_2$Pathway <- factor(Plant_count_2$Pathway, levels = ordered_pathways)

# Define a function to wrap text to a specified width
wrap_text <- function(s, width = 20) {
  # Use strwrap to wrap text and then paste it back together with \n to create multi-line text
  sapply(s, function(x) paste(strwrap(x, width = width), collapse = "\n"))
}

# Preprocess the labels to insert line breaks after adjusting the factor levels
Plant_count_2$Pathway <- wrap_text(as.character(Plant_count_2$Pathway))

# Ensure the factor levels are maintained after wrapping text
Plant_count_2$Pathway <- factor(Plant_count_2$Pathway, levels = wrap_text(as.character(ordered_pathways)))

# Define the custom labels for the legend
custom_labels_plant <- c(Arabidopsis = "italic('A. thaliana')", Barley = "italic('H. vulgare')", Lotus = "italic('L. japonicus')")

# Define the theme
legend_theme <- theme(
  legend.key.size = unit(0.2, "in"),      # Adjust legend key size for better fit
  legend.text = element_text(size = 10),   # Ensure the legend text is readable
  strip.background = element_rect(colour = "gray50", fill = "transparent", size = 1), # Change 'size' for thickness
  axis.text = element_text(color = "gray50"),
  axis.line = element_line(color = "gray50", size = 0.3),
  axis.line.x = element_line(color = "gray50", size = 0.3), 
  axis.line.y = element_line(color = "gray50", size = 0.3),
  axis.ticks.x = element_line(color = "gray50", size = 0.3),
  axis.ticks.y = element_line(color = "gray50", size = 0.3),
  axis.text.x = element_text(size = 10, angle = 45, hjust = 1), 
  axis.title = element_text(size = 12), 
  axis.text.y = element_text(size = 10), 
  legend.title = element_text(size = 10),
  strip.text.x = element_text(size = 10),
  strip.text.y = element_text(size = 10),
  axis.title.x = element_blank(), 
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  legend.position = c(0.9, 0.6),          # Position the legend inside the plot
  legend.box = "vertical"              # Arrange legends horizontally
)

# Plotting
g2 <- ggplot(Plant_count_2, aes(x = Pathway, weight = Count, fill = Plant)) +
  theme_classic() +
  geom_bar(position = "dodge", width = 0.5, just = 0.5) + 
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3"), 
                    labels = parse(text = custom_labels_plant)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ylab("Number of significant KOs") + 
  xlab("Pathway") +
  labs(fill = "Plant species") +
  legend_theme

# Display the plot
g2

# Define the dimensions in centimeters
width_cm <- 21
height_cm <- 29.7/2

# Convert dimensions to inches (1 inch = 2.54 cm)
width_in <- width_cm / 2.54
height_in <- height_cm / 2.54

pdf(paste(results.dir, "Figure_5d_ABC_bar_chart.pdf", sep=""), width=width_in, height=height_in)
print(g2)
dev.off()

#Heatmap in figure 5d 
Plants <- c("At", "Hv", "Lj")
SynComs <- c("AtSC", "HvSC", "LjSC", "SSC")

KO_gene <- read.table(paste(working_directory,"/Annotations/Ann_ABC_2.txt", sep = ""), header=T, sep="\t")

new_2 <- data.frame()

for (syncom in SynComs){
  table <- read.table(paste(working_directory,"/KO_genome/KO_",syncom,".tsv", sep = ""), sep= "\t", header =T, row.names =1) 
  
  table_2 <- table[row.names(table) %in% KO_gene$KO,]
  table_3 <- rowSums(table_2)
  value <- sum(table_3)/length(colnames(table_2))
  print(sum(table_3))
  new <- data.frame(t(data.frame(c(paste(syncom), value))))
  new_2 <- rbind(new_2,new)
}

row.names(new_2) <- NULL
colnames(new_2) <- c("SynCom", "No_of_ABC_KO")

data_together <- data.frame()

for (plant in Plants){
  for (syncom in SynComs){
    norm_SSC=read.table(paste(working_directory, "/Isolate_tables/Original/", syncom,"_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
    
    norm_SSC_2 <- norm_SSC[, grep(paste(syncom),colnames(norm_SSC))]
    norm_SSC_3 <- norm_SSC_2[, grep(paste(plant, "_", sep = ""),colnames(norm_SSC_2))]
    norm_SSC_4 <- norm_SSC_3[, grep("ES_",colnames(norm_SSC_3))]
    norm_SSC_5 <- norm_SSC_4[, !grepl("HL",colnames(norm_SSC_4))]
    
    norm_SSC_6 <- t(t(norm_SSC_5)/rowSums(t(norm_SSC_5)))
    norm_SSC_7 <- data.frame(rowSums(norm_SSC_6)/length(colnames(norm_SSC_6)))
    colnames(norm_SSC_7) <- "Abun"
    norm_SSC_8 <- norm_SSC_7[order(norm_SSC_7$Abun, decreasing =T),]
    norm_SSC_9 <- row.names(norm_SSC_7)[order(norm_SSC_7$Abun, decreasing =T)]
    
    isolates <- c()
    for (j in 1:length(norm_SSC_8)){
      value <- sum(norm_SSC_8[1:j])
      isolates <- c(isolates, value)
    }
    isolates_2 <- length(isolates[isolates < 0.75])
    
    norm_SSC_10 <- norm_SSC_9[1:isolates_2]
    
    table <- read.table(paste(working_directory, "KO_genome/KO_",syncom,".tsv", sep = ""), sep= "\t", header =T, row.names =1) 
    table_2 <- table[row.names(table) %in% KO_gene$KO,]
    colnames(table_2) <- gsub("X", "", colnames(table_2))
    table_3 <- table_2[,colnames(table_2) %in% norm_SSC_10]
    
    table_4 <- rowSums(table_3)
    value <- sum(table_3)/isolates_2
    print(sum(table_3))
    
    value_all <- new_2$No_of_ABC_KO[new_2$SynCom == paste(syncom)]
    
    prop <- value/as.numeric(value_all)
    
    data_together <- rbind(data_together, data.frame(t(data.frame(c(paste(plant), paste(syncom), prop)))))
  }
}

row.names(data_together) <- NULL
colnames(data_together) <- c("Plant", "SynCom", "Proportion")
data_together$Proportion <- as.numeric(data_together$Proportion)

data_together$Plant

plot <- ggplot(data_together, aes(SynCom, Plant)) +
  geom_tile(aes(fill = Proportion)) +
  geom_text(aes(label = round(as.numeric(Proportion),2)), size =6) +
  scale_fill_gradient2(low = "red", mid = "white", high = "green", midpoint =1, na.value = "lightgrey")+
  theme_classic() +
  labs(x ="", y = "", fill = "Enrichment") +
  theme(panel.background=element_blank(),panel.grid=element_blank(),axis.line.x=element_line(size=.5, colour="black"),axis.line.y=element_line(size=.5, colour="black"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black", size=7),legend.position="right",legend.text= element_text(size=10),text=element_text(family="sans", size=10))+
  theme(axis.text.x = element_text(size = 14, angle = 25,hjust=1),axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=18)) +
  ggtitle("ABC transporter enrichment") +
  theme(plot.title = element_text(hjust = 0.5))
plot

pdf(paste(results.dir, "Figure_5d_ABC_heatmap.pdf", sep=""), width=8, height= 4)
print(plot)
dev.off()
