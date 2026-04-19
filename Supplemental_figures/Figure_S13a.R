library("dplyr") #Version 1.1.2
library("ggplot2") #Version 3.4.2
library("ggpubr") #Version 0.6.0

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Figure S13a - Intrafamily functional variation =====
SSC_KO_profiles <- read.table(paste(working_directory,"KO_genome/KO_SSC.tsv", sep = ""), sep= "\t", header =T, row.names =1) 
colnames(SSC_KO_profiles) <- gsub("X","", colnames(SSC_KO_profiles))

colnames(SSC_KO_profiles)[grep("M.16",colnames(SSC_KO_profiles))] <- "M-16"
colnames(SSC_KO_profiles)[grep("M.6",colnames(SSC_KO_profiles))] <- "M-6"
colnames(SSC_KO_profiles)[grep("M.10",colnames(SSC_KO_profiles))] <- "M-10"
colnames(SSC_KO_profiles)[grep("M.11_2",colnames(SSC_KO_profiles))] <- "M-11_2"
colnames(SSC_KO_profiles)[571]  <- "M-11"

tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)

tax_df_3 <- table(tax_df_2$family)
Families <- names(tax_df_3)[tax_df_3 > 10]

together_2 <- data.frame()
order_data_2 <- data.frame()

for (family in Families){
  isolates <- row.names(tax_df_2)[tax_df_2$family == paste(family)]
  SSC_KO_profiles_sub <- SSC_KO_profiles[,colnames(SSC_KO_profiles) %in% isolates]
  SSC_KO_profiles_sub[SSC_KO_profiles_sub > 0] <- 1
  new_data <- data.frame(rowSums(SSC_KO_profiles_sub)/length(colnames(SSC_KO_profiles_sub))*100)
  colnames(new_data) <- "Proportion"
  
  # Define KO occurrence categories
  Core <- data.frame(group = "Core (present in 100% of genomes)",
                     No_of_KOs = sum(new_data$Proportion == 100))
  Soft_Core <- data.frame(group = "Soft core (present in 90–99.9%)",
                          No_of_KOs = sum(new_data$Proportion < 100 & new_data$Proportion >= 90))
  Shell <- data.frame(group = "Shell (present in 50–89.9%)",
                      No_of_KOs = sum(new_data$Proportion < 90 & new_data$Proportion >= 50))
  Cloud <- data.frame(group = "Cloud (present in <50%)",
                      No_of_KOs = sum(new_data$Proportion < 50 & new_data$Proportion > 0))
  
  together <- rbind(Core,Soft_Core,Shell,Cloud)
  
  order_data <- data.frame(t(data.frame(c(paste(family),sum(together$No_of_KOs)))))
  row.names(order_data) <- NULL
  colnames(order_data) <- c("Family","No_of_KOs")
  
  order_data_2 <- rbind(order_data_2,order_data)
  
  together$Family <- paste(family)
  
  together_2 <- rbind(together_2, together)
}

order_data_3 <- order_data_2$Family[order(order_data_2$No_of_KOs, decreasing =T)]

together_2$Family <- factor(together_2$Family, levels = order_data_3)
together_2$group <- factor(together_2$group, levels = c("Cloud (present in <50%)","Shell (present in 50–89.9%)","Soft core (present in 90–99.9%)","Core (present in 100% of genomes)"))

g2 <- ggplot(together_2, aes(x = Family, weight = No_of_KOs, fill = group)) +
  theme_classic() +
  geom_bar(position = "stack", width = 0.5, just = 0.5) + 
  scale_fill_manual(values =  c("#eaebfe","#b0b8ce","#505a74","#022954")) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ylab("Number of KOs") + 
  xlab("Family") +
  ylim(0,6000) +
  labs(fill = "Distribution across family") +
  ggtitle("Family pangenome") +
  theme( axis.text.x=element_text(size = 12, angle = 25, hjust =1,face = "italic",), 
         axis.title.x=element_blank(), 
         title=element_text(hjust=0.5, size=15), 
         strip.background=element_rect(colour="gray50", size=0.3), # Change 'size' for thickness
         axis.text.y=element_text(color="gray50", size =12),
         axis.line = element_line(color="gray50", size=0.3))

# Display the plot
g3 <- ggarrange(NULL, g2, ncol=2, widths=c(1,14))

pdf(paste(results.dir,"Figure_S13a_Family_pangenome.pdf", sep=""), width=14, height=8)
print(g3)
dev.off()
