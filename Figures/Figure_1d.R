library("ggplot2") #Version 3.4.2
library("ggvenn") #Version 0.1.10
library("patchwork") #Version 1.2.0

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Figure 1d - KOs in SSC =====

table_6 <- read.table(paste(working_directory, "KO_intravariability/pangenome_order.tsv", sep=""), sep = "\t", header = T)
SynComs <- c("AtSC", "HvSC", "LjSC")

for (inoculum in SynComs) {
  table <- read.table(paste(working_directory,"KO_genome/KO_",inoculum,".tsv", sep = ""), sep= "\t", header =T) 
  if (inoculum == "AtSC"){
    genes_AtSC <- table$sequence
  } else if (inoculum == "LjSC"){
    genes_LjSC <- table$sequence
  } else {
    genes_HvSC <- table$sequence
  }
}
x <- list(
  AtSC = sample(genes_AtSC), 
  LjSC = sample(genes_LjSC), 
  HvSC = sample(genes_HvSC)
)

hop <- data.frame(unique(c(genes_AtSC, genes_LjSC, genes_HvSC)))
colnames(hop) <- "KO"
hop$AtSC <- FALSE
hop$HvSC <- FALSE
hop$LjSC <- FALSE
hop$AtSC[hop$KO %in% genes_AtSC] <- TRUE
hop$HvSC[hop$KO %in% genes_HvSC] <- TRUE
hop$LjSC[hop$KO %in% genes_LjSC] <- TRUE

Venn <- ggplot(hop) +
  geom_venn(aes(A = AtSC, B = HvSC, C =LjSC ), fill_color = c("#A3A500", "#00B0F6", "#00BF7D"), show_percentage = FALSE, set_name_size = 6,text_size = 4) + 
  theme_void() + ggtitle("KO overlap") +  theme(plot.title = element_text(hjust = 0.5))

#plot
plot <- ggplot(table_6,aes(x = variable, y = value, group = SynCom)) +
  geom_line(data = table_6, 
            mapping = aes(x = variable, y = value, color = SynCom),size=1.2) + 
  ylab("Total number of KOs") + 
  xlab("No. of isolates") +
  scale_color_manual(values = c("#A3A500","#00B0F6","#00BF7D","#F8766D"))+
  theme(axis.text.x = element_text(size = 14), axis.title = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=24)) +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) +
  ggtitle("Pangenome openness curve") +  theme(plot.title = element_text(hjust = 0.5))
plot

plot_2 <- plot + inset_element(Venn, left = 0.4, bottom = 0.05, right = 1, top = 0.75)
plot_2

pdf(paste(results.dir,"Figure_1d_SSC_curve.pdf", sep=""), width=8, height=6)
print(plot_2)
dev.off()
