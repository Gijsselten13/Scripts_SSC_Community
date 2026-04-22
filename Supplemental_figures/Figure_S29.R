library("ggplot2") #Version 3.4.2
library("ggvenn") #Version 0.1.10
library("patchwork") #Version 1.2.0

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Figure S29 - KOs in SSC and Levy et al. (2018) data =====
input_table <- read.table(paste(working_directory, "/DESeq2/Sig_KO_all_no_nod_rhizo.txt", sep = ""), header=T, sep="\t")

input_table_2 <- table(input_table$KO)
input_table_3 <- names(input_table_2)[input_table_2 == 12]

SynComs <- c("AtSC", "LjSC", "HvSC")

disttrib_plot_df <- data.frame()

for (inoculum in SynComs){
  table <- read.table(paste(working_directory,"KO_genome/KO_",inoculum,".tsv", sep = ""), sep= "\t", header =T, row.names =1) 
  colnames(table) <- gsub("X", "", colnames(table))
  table[table > 0] <- 1
  table_2 <- table[row.names(table) %in% input_table_3,]
  
  for (isolate in colnames(table)){
    table_3 <- table_2[,colnames(table_2) == paste(isolate)]
    value_isolate <- sum(table_3)
    
    hop <- data.frame(t(data.frame(c(paste(isolate), paste(inoculum), value_isolate))))
    
    disttrib_plot_df <- rbind(disttrib_plot_df, hop)
  }
}

row.names(disttrib_plot_df) <- NULL
colnames(disttrib_plot_df) <- c("Isolate", "SynCom", "No_of_KOs")
disttrib_plot_df$No_of_KOs <- as.numeric(disttrib_plot_df$No_of_KOs)

plot_SSC <- ggplot()+
  geom_density(data = disttrib_plot_df, aes(x = No_of_KOs, fill = SynCom),
               alpha = 0.5, size = 0.2)+
  scale_fill_manual(values = c("#A3A500","#00B0F6","#00BF7D"))+
  theme_classic() +
  xlab("No of KOs") +
  ylab("Density") +
  xlim(0,266) +
  ggtitle("Distribution 266 KOs") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), axis.text =element_text(size = 16),axis.title =element_text(size = 18) )

plot_SSC

pdf(paste(results.dir, "Figure_S29a_Dist_266.pdf", sep=""), width=7, height=5)
print(plot_SSC)
dev.off()

#Figure S29b - Levy data - 266 KO dataset 
input_table <- read.table(paste(working_directory, "DESeq2/Sig_KO_all_no_nod_rhizo.txt", sep = ""), header=T, sep="\t")

input_table_2 <- table(input_table$KO)
input_table_3 <- names(input_table_2)[input_table_2 == 12]

all_SSC_KOs <- read.table(paste(working_directory,"KO_tables/Original/SSC.tsv", sep = ""), header=T, sep="\t")
all_SSC_KOs_2 <- all_SSC_KOs$function.

Levy_KO <- read.table(paste(working_directory, "Levy/Levy_genomes_ko.tsv", sep = ""), sep = "\t", header =T, row.names=1)

metadata <- read.table(paste(working_directory,"Levy/metadata.txt", sep = ""), header =T, sep = "\t", row.names=1)
metadata$Classification[metadata$Classification == "soil"] <- "Soil"

hop_2 <- data.frame()

colnames(Levy_KO) <- gsub("X", "", colnames(Levy_KO))
Levy_KO[Levy_KO > 0] <- 1
Levy_KO_2 <- Levy_KO[row.names(Levy_KO) %in% input_table_3,]

for (isolate in colnames(Levy_KO_2)){
  Levy_KO_3 <- Levy_KO_2[,colnames(Levy_KO_2) == paste(isolate)]
  value_isolate <- sum(Levy_KO_3)
  status <- metadata$Classification[row.names(metadata) == paste(isolate)]
  
  hop <- data.frame(t(data.frame(c(paste(isolate), paste(status), value_isolate))))
  
  hop_2 <- rbind(hop_2, hop)
}

row.names(hop_2) <- NULL
colnames(hop_2) <- c("Isolate", "Status", "No_of_KOs")
hop_2$No_of_KOs <- as.numeric(hop_2$No_of_KOs)
hop_2$Status <- factor(hop_2$Status, levels = c("NPA", "Soil", "PA"))

plot <- ggplot()+
  geom_density(data = hop_2, aes(x = No_of_KOs, fill = Status),
               alpha = 0.5, size = 0.2)+
  theme_classic() +
  xlab("No of KOs") +
  ylab("Density") +
  xlim(0,266) +
  ggtitle("Distribution 266 KOs in Levy et al. (2018)") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), axis.text =element_text(size = 16),axis.title =element_text(size = 18) )
plot

hop_venn <- as.vector(read.table(paste(working_directory,"Levy/Levy_list_PA_sig_KOs.txt", sep = ""), sep = "\t", header = F))
hop_venn_1 <- hop_venn$V1[hop_venn$V1 %in% all_SSC_KOs_2]
input_table_4 <- input_table_3[input_table_3 %in% row.names(Levy_KO)]

x <- list(
  Levy = sample(hop_venn_1), 
  SSC = sample(input_table_4)
)

hop_venn_2 <- data.frame(unique(c(x$Levy, x$SSC)))
colnames(hop_venn_2) <- "KO"
hop_venn_2$Levy <- FALSE
hop_venn_2$SSC <- FALSE

hop_venn_2$Levy[hop_venn_2$KO %in% hop_venn_1] <- TRUE
hop_venn_2$SSC[hop_venn_2$KO %in% input_table_4] <- TRUE

colnames(hop_venn_2) <- c("OG", "Levy et al. (2018)", "SSC")

Venn <- ggplot(hop_venn_2) +
  geom_venn(aes(A = `Levy et al. (2018)`, B = SSC), fill_color = c("gray50", "gray80"), show_percentage = FALSE, set_name_size = 6,text_size = 6) + 
  theme_void() +  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size = 10))

Venn
plot_2 <- plot + inset_element(Venn, left = 0.2, bottom = 0.25, right = 0.95, top = 0.95)
plot_2

pdf(paste(results.dir, "Figure_S26b_Levy_dist_266.pdf", sep=""), width=7, height=5)
print(plot_2)
dev.off()

#Enrichment p-value
Levy_table <- read.table(paste(working_directory, "Levy/Levy_genomes_ko.tsv", sep = ""), header =T, row.names =1)
Levy_table_2 <- rowSums(Levy_table)
no_of_KOs <- length(names(Levy_table_2)[Levy_table_2 !=0])
unique_Levy_PA_KOs <- 3566

KO_table =read.table(paste(working_directory,"KO_genome/KO_SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
value_all <- length(row.names(KO_table)[row.names(KO_table) %in% names(Levy_table_2)[Levy_table_2 !=0]])
Expectancy <- unique_Levy_PA_KOs/value_all

KOs_present_in_Levy_from_266 <- 228
overlap <- 177
together <- c(overlap, KOs_present_in_Levy_from_266 - overlap)

binom_out_RA <- binom.test(together,KOs_present_in_Levy_from_266, Expectancy)

Enrichment <- binom_out_RA$p.value
Enrichment
