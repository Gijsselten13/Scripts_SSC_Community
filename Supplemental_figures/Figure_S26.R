library("poolr") #Version 1.1-1
library("ggplot2") #Version 3.4.2
library("ggpubr") #Version 0.6.0

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Figure S26 - Metabolite dot plot =====
input_table <- read.table(paste(working_directory, "/DESeq2/Sig_KO_all.txt", sep = ""), header=T, sep="\t")

top <- read.table(paste(working_directory,"/Annotations/pathway_top.txt", sep = ""), header=F, sep="\t")
KO_to_pathway <- read.table(paste(working_directory,"/Annotations/KO_to_pathway.txt", sep = ""), header=T, sep="\t")
KO_to_pathway$V3 <- top$V2[match(KO_to_pathway$V2, top$V1)]

input_table$pathway <- KO_to_pathway$V3[match(input_table$KO, KO_to_pathway$V1)]
input_table_2 <- input_table[!is.na(input_table$pathway),]

input_table_3 <- input_table_2[input_table_2$padj < 0.05,]

Plant <- c("Arabidopsis","Barley",  "Lotus")
Iso <- read.table(paste(working_directory,"/Shiny_app/Abundances_full.tsv", sep = ""), header=T, sep ="\t")

for (plant in Plant){
  input_table_sub <- input_table_3[input_table_3$Plant == paste(plant),]
  data_count <- data.frame(table(input_table_sub$KO))
  data_count_KOs <- data_count$Var1[data_count$Freq == 4]
  
  if (plant == "Arabidopsis"){
    Arabidopsis <- data_count_KOs
  } else if (plant == "Barley"){
    Barley <- data_count_KOs
  } else {
    Lotus <- data_count_KOs
  }
}

input_table <- read.table(paste(working_directory, "/DESeq2/Sig_KO_all.txt", sep = ""), header=T, sep="\t")
KO_gene <- read.table(paste(working_directory,"/Annotations/Ann_ABC_5_FL.txt", sep = ""), header=T, sep="\t")

selection <- table(KO_gene$Category)

selection_2 <- names(selection)[selection > 10]

KO_gene_2 <- KO_gene[KO_gene$Category %in% selection_2,]

At_KO <- Arabidopsis[Arabidopsis %in% KO_gene_2$KO]
Hv_KO <- Barley[Barley %in% KO_gene_2$KO]
Lj_KO <- Lotus[Lotus %in% KO_gene_2$KO]

At_KO_2 <- unique(KO_gene_2$Substrate[KO_gene_2$KO %in% At_KO])
Hv_KO_2 <- unique(KO_gene_2$Substrate[KO_gene_2$KO %in% Hv_KO])
Lj_KO_2 <- unique(KO_gene_2$Substrate[KO_gene_2$KO %in% Lj_KO])

SynComs <- c("AtSC", "HvSC", "LjSC", "SSC")
Plant <- c("Arabidopsis", "Barley", "Lotus")

new_data_2 <- data.frame()

for (plant in Plant){
  if (plant == "Arabidopsis"){
    genes <- At_KO_2
  } else if (plant == "Barley"){
    genes <- Hv_KO_2
  } else if (plant == "Lotus"){
    genes <- Lj_KO_2
  }
  
  for (gene in unique(genes)){
    KOs <- KO_gene$KO[KO_gene$Substrate == paste(gene)]
    input_table_2 <- input_table[input_table$Plant == paste(plant),]
    input_table_3 <- input_table_2[input_table_2$KO %in% KOs,]
    
    Iso_2 <- Iso[Iso$Gene %in% KOs,]
    Iso_3 <- Iso_2[Iso_2$Plant == paste(plant),]
    
    data_sub <- data.frame()
    
    for (syncom in SynComs){
      if (syncom != "SSC"){
        Iso_sub <- Iso_3[Iso_3$SynCom == paste(syncom),]
        data_sub <- rbind(data_sub, Iso_sub)
      } else {
        Iso_sub <- Iso_3[Iso_3$SynCom == paste(syncom),]
        for (KO in unique(Iso_sub$Gene)){
          Iso_sub_2 <- Iso_sub[Iso_sub$Gene == paste(KO),]
          inbetween <- data.frame(t(data.frame(c(unique(Iso_sub_2$Gene),paste(syncom), paste(plant), sum(as.numeric(Iso_sub_2$RA_KO))/length(Iso_sub_2$RA_KO), "SSC", sum(as.numeric(Iso_sub_2$RA_Iso))/length(Iso_sub_2$RA_Iso)))))
          row.names(inbetween) <- NULL
          colnames(inbetween) <- colnames(Iso_sub)
          data_sub <- rbind(data_sub,inbetween )
        }
      }
    }
    
    if (length(input_table_3$KO) != 0){
      value_com <- stouffer(input_table_3$padj)
      value <- value_com$p
      value_RA <- sum(as.numeric(data_sub$RA_Iso))/length(data_sub$RA_Iso)
      
      if (value == 0){
        value <- "-Inf"
      }
    } else {
      value <- 0
      value_RA <- 0
    }
    length_value <- length(unique(input_table_3$KO))/length(KOs)
    
    new_data <- data.frame(t(data.frame(c(paste(gene), paste(plant), value, value_RA,length_value))))
    
    new_data_2 <- rbind(new_data_2, new_data)
  }
}

row.names(new_data_2) <- NULL
colnames(new_data_2) <- c("Substrate", "Plant", "padj","RA_Iso","no_of_genes")

new_data_sub_sub <- new_data_2[new_data_2$padj != "-Inf",]
new_data_sub_sub_2 <- new_data_sub_sub[new_data_sub_sub$padj != 0,]

new_data_2$padj[new_data_2$padj == "-Inf"] <- min(as.numeric(new_data_sub_sub_2$padj))
new_data_2$neglogp <- -log10(as.numeric(new_data_2$padj))
new_data_2$neglogp[new_data_2$neglogp == "Inf"] <- NA
new_data_2$no_of_genes[is.na(new_data_2$neglogp)] <- NA

new_data_2$no_of_genes <- as.numeric(new_data_2$no_of_genes)

order_2 <- data.frame()
order_3 <- data.frame()
order_4 <- data.frame()

for (order in unique(new_data_2$Substrate)){
  new_data_sub <- new_data_2[new_data_2$Substrate == paste(order),]
  if (length(new_data_sub$padj) == 3){
    order_4 <- rbind(order_4, data.frame(t(data.frame(c(paste(order), new_data_sub$padj, new_data_sub$neglogp)))))
  } else if (length(new_data_sub$padj) == 2){
    order_3 <- rbind(order_3,data.frame(t(data.frame(c(paste(order), sum(as.numeric(new_data_sub$padj))/2, sum(as.numeric(new_data_sub$neglogp))/2)))))
  } else {
    order_2 <- rbind(order_2,data.frame(t(data.frame(c(paste(order), sum(as.numeric(new_data_sub$padj))/3, sum(as.numeric(new_data_sub$neglogp))/3)))))
  }
}

order_2_sub <- order_2$X1[order(as.numeric(order_2$X2), decreasing =F)]
order_3_sub <- order_3$X1[order(as.numeric(order_3$X2), decreasing =F)]
order_4_sub <- order_4$X1[order(as.numeric(order_4$X2), decreasing =F)]

final_order <- c(order_4_sub,order_3_sub, order_2_sub)

new_data_2$Category <- KO_gene$Category[match(new_data_2$Substrate, KO_gene$Substrate)]

new_data_3 <- new_data_2[new_data_2$Category %in% selection_2,]

new_data_sub <- new_data_3[,c(1,6)]
new_data_sub_2 <- unique(new_data_sub)
new_data_sub_3 <- table(new_data_sub_2$Category)

new_data_3$Substrate <- factor(new_data_3$Substrate, levels = final_order)

# Curation of KO substrate and categorization

#These are removed as they are double in the dataset (Manganese as in Iron(II)/Manganese and Histidine as in Neutral amino acid/Histidine )
missing <- c("Histidine", "Manganese" )
new_data_3_sub <- new_data_3[!new_data_3$Substrate %in% missing,]

new_data_3_sub$Category_better <- KO_gene$Category_better[match(new_data_3_sub$Substrate, KO_gene$Substrate)]
new_data_3_sub$Substrate_2 <- KO_gene$Substrate_2[match(new_data_3_sub$Substrate, KO_gene$Substrate)]
new_data_3_sub$Column_category <- KO_gene$Column_category[match(new_data_3_sub$Substrate, KO_gene$Substrate)]

legend_theme <- theme(legend.key.size = unit(0.3, "in"),      # Adjust legend key size for better fit
                      legend.text = element_text(size = 8),   # Ensure the legend text is readable
                      strip.background=element_rect(colour="gray50",fill = "transparent", size=1), # Change 'size' for thickness
                      axis.text=element_text(color="gray50"),
                      axis.line = element_line(color="gray50", size=0.3),
                      axis.line.x = element_line(color="gray50", size=0.3), 
                      axis.line.y = element_line(color="gray50", size=0.3),
                      axis.ticks.x = element_line(color="gray50", size=0.3),
                      axis.ticks.y =element_line(color="gray50", size=0.3),
                      axis.text.x = element_text(size = 10, angle = 45, hjust = 1), 
                      axis.title = element_text(size = 10), 
                      axis.text.y = element_text(size = 10), 
                      legend.title = element_text(size = 10),
                      strip.text.x = element_text(size = 10),
                      strip.text.y = element_text(size = 9),
                      panel.border = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      legend.position = "right",             # Move legend to the bottom
                      legend.box = "vertical"              # Arrange legends horizontally
)
plot_list=list()
for (category in na.omit(unique(new_data_3_sub$Column_category))){
  
  new_data_4 <- na.omit(new_data_3_sub[new_data_3_sub$Column_category == paste(category),])
  new_data_4$RA_Iso <- as.numeric(new_data_4$RA_Iso)
  
  dot_plot <- ggplot(new_data_4, aes(y=Substrate_2, x=Plant, color=neglogp, size = RA_Iso)) + 
    geom_point() + 
    labs(x = "",y = "", size = "RA Isolates") +
    scale_color_gradient2(midpoint=mean(as.numeric(na.omit(new_data_3$neglogp))), low="blue",high="red", mid = "purple", space ="Lab", limits = c(0,max(as.numeric(na.omit(new_data_3$neglogp)))))+
    scale_size_continuous(range = c(1, 4))+
    facet_grid(Category_better ~ . ,scales='free', space = 'free') +
    legend_theme
  dot_plot
  
  plot_list[[category]]=dot_plot
  
}

dotplot_final= ggarrange(plot_list[[2]], plot_list[[1]], nrow = 1, common.legend = T, legend = "right", widths = c(1,0.95))
dotplot_final

pdf(paste(results.dir,"Figure_S26_metabolite_dotplot.pdf", sep=""), width=10, height= 8)
print(dotplot_final)
dev.off()
