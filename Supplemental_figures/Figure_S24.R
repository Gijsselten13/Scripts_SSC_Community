library("ggplot2") #Version 3.4.2

###Figure S24 - Cologne AtCC vs LjCC =====
AtCC_Cologne <- read.table(paste(working_directory, "KO_genome/KO_AtCC_Cologne.tsv", sep = ""), header =T, row.names =1)
AtCC_Cologne[AtCC_Cologne > 0] <- 1
AtCC_Cologne_summary <- data.frame(colSums(AtCC_Cologne))
colnames(AtCC_Cologne_summary) <- "No_of_KOs"
AtCC_Cologne_summary$Collection <- "AtCC"

LjSC <- read.table(paste(working_directory, "KO_genome/KO_LjSC.tsv", sep = ""), header =T, row.names =1)
LjSC[LjSC > 0] <- 1
LjSC_summary <- data.frame(colSums(LjSC))
colnames(LjSC_summary) <- "No_of_KOs"
LjSC_summary$Collection <- "LjCC"

combined <- rbind(AtCC_Cologne_summary,LjSC_summary)

plot_LjCC_vs_AtSC <- ggplot()+
  geom_density(data = combined, aes(x = No_of_KOs, fill = Collection),
               alpha = 0.5, size = 0.2)+
  scale_fill_manual(values = c("purple","#00BF7D"))+
  theme_classic() +
  xlab("No of KOs") +
  ylab("Density") +
  ggtitle("Functional diversity LjCC vs AtCC Cologne soil") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), axis.text =element_text(size = 16),axis.title =element_text(size = 18) )

plot_LjCC_vs_AtSC

ks.test(AtCC_Cologne_summary$No_of_KOs, LjSC_summary$No_of_KOs)

pdf(paste(results.dir,"Figure_S24_AtCC_vs_LjCC_Cologne.pdf", sep=""), width=10, height=6)
print(plot_LjCC_vs_AtSC)
dev.off()
