#Libraries - R version 4.4
library("ggplot2") #Version 3.4.2
library("grid") #Version 4.4.1
library("dplyr") #Version 1.1.2
library("ggvenn") #Version 0.1.10
library("patchwork") #Version 1.2.0
library("cowplot") #Version 1.1.3
library("ape") #Version 5.7-1
library("reshape2") #Version 1.4.4
library("stringr") #Version 1.5.0
library("ggpubr") #Version 0.6.0
library("rstatix") #Version 0.7.2
library("multcompView") #Version 0.1-9
library("plyr") #Version 1.8.8
library("phyloseq") #Version 1.44.0
library("vegan") #Version 2.6-4
library("ggrepel") #Version 0.9.3
library("forcats") #Version 1.0.0
library("moonBook") #Version 0.3.3
library("webr") #Version 0.1.6
library("tidyheatmaps") #Version 0.1.0
library("ggtern") #Version 3.4.1
library("scales") #Version 1.2.1
library("poolr") #Version 1.1-1
library("KEGGREST") #Version 1.40.0
library("ggpp") #Version 0.5.5
library("DESeq2") #Version 1.40.0
library("mixOmics") #Version 6.30.0
library("see") #Version 0.11.0
library("multcomp") #Version 1.4-28
library("mclust") #Version 6.1.1
library("tidyr") #Version 1.3.0
library("nloptr") #Version 2.1.1
library("ANCOMBC") #Version 2.6.0

#Set working directory and results folder here
working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")
dir.create(paste(results.dir, "/Small_plots", sep = ""))
results.dir_2 <- paste(results.dir,"Small_plots/", sep = "")

############################# Figures #############################################
###Figures 1b and 1c - AtCC and HvCC Coverage of natural communities =====
AtCC <- read.table(paste(working_directory, "NatCom/AtCC_NatCom_table.tsv", sep = ""), sep = "\t", header =T)
AtCC_kaiju_tax <- read.table(paste(working_directory,"NatCom/SSC_taxonomy_kaiju.tsv", sep =""), header=T,sep="\t",quote="\"", fill = FALSE)

data_wide <- spread(AtCC, sample, percent)
row.names(data_wide) <- data_wide$taxon_name
data_wide_2 <- data_wide %>% dplyr::select(-taxon_name)

data_wide_3 <- data.frame(rowSums(data_wide_2)/length(colnames(data_wide_2)))
colnames(data_wide_3) <- "Percent"
data_wide_3$Taxon <- row.names(data_wide_3)

data_wide_4 <- data_wide_3[data_wide_3$Taxon != "unclassified",]
data_wide_5 <- data_wide_4[data_wide_4$Taxon != "cannot be assigned to a (non-viral) genus",]

#Filter out low abundant genera
data_wide_5_sub <- data_wide_5[data_wide_5$Percent > 0.05,]
data_wide_5 <- data_wide_5_sub

data_wide_5$Rel <- data_wide_5$Percent/sum(data_wide_5$Percent)

data_wide_6 <- data_wide_5[order(data_wide_5$Rel, decreasing =TRUE),]
data_wide_7 <- data_wide_6[1:70,]

AtCC_kaiju_tax_2 <- AtCC_kaiju_tax[AtCC_kaiju_tax$SynCom == "AtSC",]
AtCC_Genera <- unique(AtCC_kaiju_tax_2$genus)

data_wide_7$Retrieved <- "Not recovered"
data_wide_7$Retrieved[data_wide_7$Taxon %in% AtCC_Genera] <- "Recovered"

data_wide_7$Taxon <- factor(data_wide_7$Taxon, levels = data_wide_7$Taxon)

#Calculating recovery rate
data_wide_6$Retrieved <- "Not recovered"
data_wide_6$Retrieved[data_wide_6$Taxon %in% AtCC_Genera] <- "Recovered"

recovery_rate <- sum(data_wide_6$Rel[data_wide_6$Retrieved == "Recovered"])
recovery_rate_2 <- paste("recovery rate: ",round(recovery_rate*100, 2), "%", sep ="")

data_wide_7$Cum <- 0

for (i in 1:length(row.names(data_wide_7))){
  new <- data_wide_7[1:i,]
  new_2 <- new[new$Retrieved != "Not recovered",]
  data_wide_7$Cum[i] <- sum(new_2$Rel)
}

data_wide_7$Cum <- as.numeric(data_wide_7$Cum)
data_wide_8=rbind(data.frame(Percent =0, Taxon="None", Cum=0,
                             Rel = 0, Retrieved = "Recovered"),
                  data_wide_7)

data_wide_8$order <- 1:length(data_wide_8$Percent)

# Get max values for scaling
max_bar <- max(data_wide_8$Rel)
max_line <- max(data_wide_8$Cum)

# Define scaling factor
scale_factor <- max_bar / max_line

g_dual_AtCC <- ggplot(data_wide_8, aes(x = order)) +
  # Barplot (left y-axis)
  geom_bar(aes(y = Rel, fill = Retrieved), stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("lightgrey", "black")) +
  
  # Line plot (scaled to bar height, right y-axis)
  geom_line(aes(y = Cum * scale_factor), size = 1, color = "#A3A500") +
  
  # Left y-axis
  scale_y_continuous(
    name = "Relative abundance in NatCom",
    # Right y-axis
    sec.axis = sec_axis(~ . / scale_factor, name = "Cumulative relative abundance")
  ) +
  annotate("text", x = 50, y=0.02, label = paste(recovery_rate_2))+
  theme_classic() +
  labs(fill = "") +
  theme(
    axis.title.y.left = element_text(size = 14),
    axis.title.y.right = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 18)
  ) +
  ggtitle("AtCC recovery")

g_dual_AtCC

pdf(paste(results.dir, "Figure_1b_AtCC_cov_combined_merged.pdf", sep = ""), width = 8, height = 5)
print(g_dual_AtCC)
dev.off()

#Now HvCC - real clustering
HvCC <- read.table(paste(working_directory, "NatCom/HvCC_NatCom_table.tsv", sep =""), sep = "\t", header =T, row.names =1)
HvCC_meta <- read.table(paste(working_directory, "NatCom/HvCC_metadata.txt", sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)

#Take only Barley - Golden Promise samples
HvCC_meta_GP <- HvCC_meta$SampleID[HvCC_meta$Genotype == "GP"]

HvCC_2 <- HvCC[, colnames(HvCC) %in% HvCC_meta_GP]

HvCC_3 <- data.frame(rowSums(HvCC_2)/length(colnames(HvCC_2)))
colnames(HvCC_3) <- "Counts"
HvCC_3$Relative_abundance <- HvCC_3$Counts/sum(HvCC_3$Counts)

HvCC_3$Retrieved <- "Not recovered"
HvCC_3$Retrieved[1:116] <- "Recovered"

HvCC_4 <- HvCC_3[HvCC_3$Counts != 0,]
HvCC_5 <- HvCC_4[order(HvCC_4$Relative_abundance, decreasing =T),]

recovery_rate <- sum(HvCC_5$Relative_abundance[HvCC_5$Retrieved == "Recovered"])
recovery_rate_2 <- paste("recovery rate: ",round(recovery_rate*100, 2), "%", sep ="")

HvCC_5$isolates <- row.names(HvCC_5)
HvCC_6 <- HvCC_5[1:70,]

HvCC_6$isolates <- factor(HvCC_6$isolates, levels = HvCC_6$isolates)

HvCC_6$Cum <- 0

for (i in 1:length(row.names(HvCC_6))){
  new <- HvCC_6[1:i,]
  new_2 <- new[new$Retrieved != "Not recovered",]
  HvCC_6$Cum[i] <- sum(new_2$Rel)
}

HvCC_6$Cum <- as.numeric(HvCC_6$Cum)
HvCC_7=rbind(data.frame(Counts =0, Cum=0,
                        Relative_abundance = 0, Retrieved = "Recovered", isolates="Null"),
             HvCC_6)


HvCC_7$order <- 1:length(HvCC_7$Relative_abundance)

# Get max values for scaling
max_bar <- max(HvCC_7$Relative_abundance)
max_line <- max(HvCC_7$Cum)

# Define scaling factor
scale_factor <- max_bar / max_line

g_dual_HvCC <- ggplot(HvCC_7, aes(x = order)) +
  # Barplot (left y-axis)
  geom_bar(aes(y = Relative_abundance, fill = Retrieved), stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("lightgrey", "black")) +
  
  # Line plot (scaled to bar height, right y-axis)
  geom_line(aes(y = Cum * scale_factor), size = 1.1, color = "#00B0F6") +
  
  # Left y-axis
  scale_y_continuous(
    name = "Relative abundance in NatCom",
    # Right y-axis
    sec.axis = sec_axis(~ . / scale_factor, name = "Cumulative relative abundance")
  ) +
  annotate("text", x = 50, y=0.02, label = paste(recovery_rate_2))+
  theme_classic() +
  labs(fill = "") +
  theme(
    axis.title.y.left = element_text(size = 14),
    axis.title.y.right = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 18)
  ) +
  ggtitle("HvCC recovery")

g_dual_HvCC

pdf(paste(results.dir, "Figure_1c_HvCC_cov_combined_merged.pdf", sep = ""), width = 8, height = 5)
print(g_dual_HvCC)
dev.off()

grid.arrange(g_dual_AtCC, g_dual_HvCC,  ncol = 2, nrow = 1)


pdf(paste(results.dir, "Figure_1bc_All_recovery_rates_merged.pdf", sep = ""), width = 16, height = 5)
grid.arrange(g_dual_AtCC, g_dual_HvCC,  ncol = 2, nrow = 1)
dev.off()
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

###Figure 1e - KOs in simulated one-strain-per-family-SynCom =====
table_6 <- read.table(paste(working_directory, "KO_intravariability/simulations/pangenome_order.tsv", sep = ""), sep = "\t", header = T)
SynComs <- 1:25

table_SC <- data.frame()

for (syncom in SynComs){
  table_6_sub <- na.omit(table_6[table_6$variable == paste(syncom),])
  average <- round(sum(table_6_sub$value)/length(table_6_sub$value),0)
  
  low_95 <-  table_6_sub$value[order(table_6_sub$value)][50]
  high_95 <-  table_6_sub$value[order(table_6_sub$value)][950]
  
  hop <- data.frame(t(data.frame(c("Average", paste(syncom),average))))
  hop_2 <- data.frame(t(data.frame(c("Low confidence interval", paste(syncom),low_95))))
  hop_3 <- data.frame(t(data.frame(c("High confidence interval", paste(syncom),high_95))))
  
  row.names(hop) <- NULL
  row.names(hop_2) <- NULL
  row.names(hop_3) <- NULL
  
  colnames(hop) <- colnames(table_6)
  colnames(hop_2) <- colnames(table_6)
  colnames(hop_3) <- colnames(table_6)
  
  table_SC <- rbind(table_SC,hop, hop_2, hop_3 )
}

table_SC$variable <- as.numeric(table_SC$variable)
table_SC$value <- as.numeric(table_SC$value)

table_SC$SynCom <- factor(table_SC$SynCom, levels = c("High confidence interval", "Average", "Low confidence interval"))

#plot
plot_1e <- ggplot(table_SC,aes(x = variable, y = value, group = SynCom)) +
  geom_line(data = table_SC, 
            mapping = aes(x = variable, y = value, color = SynCom),size=1.2) + 
  ylab("Total number of KOs") + 
  xlab("No. of isolates") +
  scale_y_continuous(limits = c(2000,8000))+
  scale_color_manual(values = c("grey", "black","grey"))+
  theme(axis.text.x = element_text(size = 14), axis.title = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=24)) +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) +
  ggtitle("One isolate per family simulation") +  theme(plot.title = element_text(hjust = 0.5, size = 18))
plot_1e

pdf(paste(results.dir,"Figure_1E_SynCom_simulation.pdf", sep=""), width=8, height=6)
print(plot_1e)
dev.off()

###Figure 1e - Pie chart for KO occurrence in 1000 syncom simulations =====
overview <- read.table(file.path(working_directory, "KO_intravariability/simulations/Sim_SynComs_per_KO.txt"), 
                       header = TRUE, sep = "\t")

# Define KO occurrence categories
Core <- data.frame(group = "Core (present in 100% of simulations)",
                   No_of_KOs = sum(overview$Perc == 100))
Soft_Core <- data.frame(group = "Soft core (present in 90–99.9%)",
                        No_of_KOs = sum(overview$Perc < 100 & overview$Perc >= 90))
Shell <- data.frame(group = "Shell (present in 50–89.9%)",
                    No_of_KOs = sum(overview$Perc < 90 & overview$Perc >= 50))
Cloud <- data.frame(group = "Cloud (present in <50%)",
                    No_of_KOs = sum(overview$Perc < 50))

# Combine and format
all <- rbind(Core, Soft_Core, Shell, Cloud)
all$Perc <- round(all$No_of_KOs / sum(all$No_of_KOs) * 100, 1)
all$Perc_label <- paste0(all$Perc, "%")

# Define order and color
all$group <- factor(all$group, levels = c("Core (present in 100% of simulations)",
                                          "Soft core (present in 90–99.9%)",
                                          "Shell (present in 50–89.9%)",
                                          "Cloud (present in <50%)"))

gray_palette <- c("Core (present in 100% of simulations)" = "#555555",
                  "Soft core (present in 90–99.9%)" = "#777777",
                  "Shell (present in 50–89.9%)" = "#AAAAAA",
                  "Cloud (present in <50%)" = "#DDDDDD")

# Plot

pie <- ggplot(all, aes(x = "", y = No_of_KOs, fill = group)) +
  geom_col(color = "black") +
  geom_label_repel(aes(label = Perc_label), 
                   color = c("white", "white", "white","black"),
                   position = position_stack(vjust = 0.5), 
                   show.legend = FALSE, size = 5) +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = gray_palette) +
  coord_polar(theta = "y") +
  ggtitle("KO occurrence frequency across 1000 random \nSynCom simulations with taxonomic diversity") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 16))
pie

# Save plot
pdf(file.path(results.dir, "Figure_1E_pie_plot.pdf"), width = 6, height = 4)
print(pie)
dev.off()

###Figure S1a & S1b - KO_intravariability SSC and simulated SynComs =====
KO_out <- read.table(paste(working_directory,"KO_intravariability/KO_intrafunctionality_2.tsv", sep = ""), sep = "\t", header =T)
row.names(KO_out) <- KO_out$KO

SynComs <- c("AtSC", "LjSC","HvSC", "SSC")

KO_out_sub <- data.frame(matrix(NA, ncol = 4))
colnames(KO_out_sub) <- c("Total_branch_length", "No_of_genes" ,"SynCom", "KO")
KO_out_sub_2 <- KO_out_sub[-1,]


for (syncom in SynComs){
  KO_out_sub <- KO_out[,grep(paste(syncom),colnames(KO_out))]
  KO_out_sub$SynCom <- paste(syncom)
  colnames(KO_out_sub) <- c("Total_branch_length", "No_of_genes" ,"SynCom")
  KO_out_sub$KO <- row.names(KO_out_sub)
  
  KO_out_sub_2 <- rbind(KO_out_sub_2, KO_out_sub)
}


KO_out_sub_2$log_genes <- log10(KO_out_sub_2$No_of_genes)

KO_out_sub_2$log_branch <- KO_out_sub_2$Total_branch_length
KO_out_sub_3 <- KO_out_sub_2[KO_out_sub_2$log_branch > 1,]
KO_out_sub_4 <- KO_out_sub_2[KO_out_sub_2$log_branch <= 1,]
KO_out_sub_3$log_branch <- log10(KO_out_sub_3$Total_branch_length)

KO_out_sub_5 <- rbind(KO_out_sub_3,KO_out_sub_4)

KO_out_sub_6 <- KO_out_sub_5[is.finite(KO_out_sub_5$log_genes),]

KO_out_sub_7 <- KO_out_sub_6[KO_out_sub_6$log_genes != 0,]
KO_out_sub_8 <- KO_out_sub_7[KO_out_sub_7$log_branch != 0,]
KO_out_sub_8$log_genes <- as.numeric(KO_out_sub_8$log_genes)
KO_out_sub_8$log_branch <- as.numeric(KO_out_sub_8$log_branch)

# Main plot
pmain <- ggplot(KO_out_sub_8, aes(x = log_branch, y = log_genes, color = SynCom, alpha =0.3))+
  scale_color_manual(values = c("#A3A500","#00B0F6","#00BF7D","#F8766D")) +
  geom_point() +
  theme_classic()+
  ylab("Log10(number of genes)") +
  xlab("Log10(total gene diversity)") +
  ggtitle("KO intravariability") +
  guides(alpha = FALSE) +
  theme(axis.text.x = element_text(size = 14), axis.title = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=24)) +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(hjust = 0.5))

# Marginal densities along x axis
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = KO_out_sub_8, aes(x = log_branch, fill = SynCom),
               alpha = 0.5, size = 0.2)+
  scale_fill_manual(values = c("#A3A500","#00B0F6","#00BF7D","#F8766D"))
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = KO_out_sub_8, aes(x = log_genes, fill = SynCom),
               alpha = 0.5, size = 0.2)+
  scale_fill_manual(values = c("#A3A500","#00B0F6","#00BF7D","#F8766D")) +
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")

ggdraw(p2)

pdf(paste(results.dir,"Figure_S1A_KO_intraspecificity.pdf", sep=""), width=8, height=6)
print(ggdraw(p2))
dev.off()

#lower quality but also lower size
ggsave(paste(results.dir,"Figure_S1A_KO_intraspecificity.pdf", sep=""), plot = ggdraw(p2), width = 8, height = 6, dpi = 300)

#S1B - KO intravariability of 1000 one-strain-per-family SynComs
KO_out <- read.table(paste(working_directory,"KO_intravariability/simulations/KO_intravariability.txt", sep = ""), sep = "\t", header =T)

KO_out_sub <- data.frame()

for (syncom in 1:1000){
  KO_out_2 <- KO_out[KO_out$SynCom == paste(syncom),]
  value <- sum(KO_out_2$Total_branch_length)
  hop <- data.frame(t(data.frame(c(paste(syncom), value))))
  row.names(hop) <- NULL
  colnames(hop) <- c("SynCom","value")
  
  KO_out_sub <- rbind(KO_out_sub, hop)
}

KO_out_sub_2 <- KO_out_sub[order(as.numeric(KO_out_sub$value), decreasing =T),]
TopSynComs <- KO_out_sub_2$SynCom[1:50]

KO_out$Conf <- "Normal"
KO_out$Conf[KO_out$SynCom %in% TopSynComs] <- "Top 50 functionally diverse"

KO_out$log_genes <- log10(KO_out$No_of_genes)

KO_out$log_branch <- KO_out$Total_branch_length
KO_out_sub_3 <- KO_out[KO_out$log_branch > 1,]
KO_out_sub_4 <- KO_out[KO_out$log_branch <= 1,]
KO_out_sub_3$log_branch <- log10(KO_out_sub_3$Total_branch_length)

KO_out_sub_5 <- rbind(KO_out_sub_3,KO_out_sub_4)

KO_out_sub_6 <- KO_out_sub_5[is.finite(KO_out_sub_5$log_genes),]

KO_out_sub_7 <- KO_out_sub_6[KO_out_sub_6$log_genes != 0,]
KO_out_sub_8 <- KO_out_sub_7[KO_out_sub_7$log_branch != 0,]
KO_out_sub_8$log_genes <- as.numeric(KO_out_sub_8$log_genes)
KO_out_sub_8$log_branch <- as.numeric(KO_out_sub_8$log_branch)
KO_out_sub_9 <- KO_out_sub_8[order(KO_out_sub_8$Conf),]

# Main plot
pmain <- ggplot(KO_out_sub_9, aes(x = log_branch, y = log_genes, color = Conf, alpha =0.3))+
  scale_color_manual(values = c("lightgray","darkgrey")) +
  geom_point() +
  theme_classic()+
  ylab("Log10(number of genes)") +
  xlab("Log10(total gene diversity)") +
  ggtitle("KO intravariability") +
  guides(alpha = FALSE) +
  theme(axis.text.x = element_text(size = 14), axis.title = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=24)) +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(hjust = 0.5))

# Marginal densities along x axis
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = KO_out_sub_9, aes(x = log_branch, fill = Conf),
               alpha = 0.5, size = 0.2)+
  scale_fill_manual(values = c("lightgray","darkgrey"))
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = KO_out_sub_9, aes(x = log_genes, fill = Conf),
               alpha = 0.5, size = 0.2)+
  scale_fill_manual(values = c("lightgray","darkgrey")) +
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")

ggdraw(p2)

pdf(paste(results.dir,"Figure_S1B_KO_intraspecificity_simulation.pdf", sep=""), width=8, height=6)
print(ggdraw(p2))
dev.off()

#Smaller size plot
ggsave(paste(results.dir,"Figure_S1B_KO_intraspecificity_simulation.pdf", sep=""), plot = ggdraw(p2), width = 8, height = 6, dpi = 300)

###Figure S3 - Shoot weights =====
importdat <- read.table(paste(working_directory,"SSC_R3_shoot_weights.txt", sep =""), sep = "\t", header=TRUE, dec = ".")

# Manipulate dataframe so it can be processed by ggplot2 package
good_table=t(importdat)
good_table=melt(good_table)
colnames(good_table)=c("local_condition", "replicate", "mass")

good_table=mutate(good_table,plant=word(local_condition,start = 1, sep=fixed("_")))
good_table=mutate(good_table,inoculum=word(local_condition,start = 2, sep=fixed("_")))
good_table=mutate(good_table,nutrient=word(local_condition,start = 4, sep=fixed("_")))
good_table=mutate(good_table,condition=word(local_condition,start = 1, sep=fixed("_"), end = 4))

# Remove NA values because of different number of samples being tested
good_table=na.omit(good_table)
good_table$mass=as.numeric(good_table$mass)

#Reformat results from first experiment
good_table$condition <-gsub("R1_", "",good_table$condition)
good_table_R1 <- good_table[grepl(pattern = "R1_",good_table$local_condition),]
good_table_R1$condition <- paste(good_table_R1$condition, "_Low", sep ="")
good_table_R1$plant <- good_table_R1$inoculum
good_table_R1$nutrient <- "Low"
new <- str_split(good_table_R1$local_condition, pattern = "_")
good_table_R1$inoculum <- data.table::transpose(new)[[3]]

good_table_R2 <- good_table[!grepl(pattern = "R1_",good_table$local_condition),]

good_table_2 <- good_table_R2

good_table_2$inoc <- factor(good_table_2$inoc, levels = c("SSC","AtSC","LjSC","HvSC","NS"))
colnames(good_table_2)[colnames(good_table_2) == "inoculum"] <- "SynCom"

PLANTS <- c("At", "Hv", "Lj")
FULL_NAMES <- c("Arabidopsis", "Barley", "Lotus")
NUTRIENT=c("Low", "High")
Plot_list=list()

max_scale <- data.frame(c("Arabidopsis", "Barley", "Lotus"), c(max(good_table$mass[good_table$plant == "At"])+0.1*max(good_table$mass[good_table$plant == "At"]),
                                                               max(good_table$mass[good_table$plant == "Hv"])+0.1*max(good_table$mass[good_table$plant == "Hv"]),
                                                               max(good_table$mass[good_table$plant == "Lj"])+0.1*max(good_table$mass[good_table$plant == "Lj"])))
colnames(max_scale) <- c("Host","Max")

for (i in 0:2) {
  for (j in 1:2){
    
    good_table_3=subset (x = good_table_2, subset = plant==PLANTS[i+1])
    good_table_3$SynCom <- factor(good_table_3$SynCom, levels = c("SSC","AtSC", "LjSC", "HvSC","NS"))
    good_table_3$Experiment <- "R2"
    good_table_3$Experiment[grepl("R1_", good_table_3$local_condition)] <- "R1"
    
    max_value <- max_scale$Max[max_scale$Host == FULL_NAMES[i+1]]
    good_table_3=subset (x = good_table_3, subset = nutrient==NUTRIENT[j])
    
    nutrient.labs <- paste(NUTRIENT[j],"nutrient")
    names(nutrient.labs) <-NUTRIENT[j]
    
    colnames(good_table_3)[colnames(good_table_3) == "SynCom"] <- "Inoculum"
    good_table_3$Inoculum <- factor(good_table_3$Inoculum, levels = c("AtSC", "HvSC", "LjSC", "SSC", "NS"))
    
    # Visualization the data using box plots. Plot weight by groups.
    
    bxpot <- ggboxplot(
      data = good_table_3,
      x = "Inoculum", 
      y = "mass",
      combine = FALSE,
      merge = FALSE,
      legend="left",
      fill =  "Inoculum",
      color = "black",
      palette = c("#A3A500","#00B0F6","#00BF7D","#F8766D","white"),
      title = ifelse(test = j==1,yes =  paste("R2 - ", FULL_NAMES[i+1], sep = ""),no =  ""),
      font.title=c(14,"italic"),
      xlab = "Inoculum",
      ylab = "Shoot mass (g)",
      bxp.errorbar = FALSE,
      bxp.errorbar.width = 0.4,
      facet.by = NULL,
      scales="free",
      panel.labs = NULL,
      short.panel.labs = TRUE,
      linetype = "solid",
      size = NULL,
      width = 0.8,
      notch = FALSE,
      outlier.shape = NA,
      select = NULL,
      remove = NULL,
      order = NULL,
      add = "none",
      add.params = list(),
      error.plot = "pointrange",
      label = NULL,
      font.label = list(size = 11, color = "black"),
      label.select = NULL,
      repel = FALSE,
      label.rectangle = FALSE,
      ggtheme = theme_pubr()
    ) +theme(axis.text.x = element_blank(),axis.title.x = element_blank(), title = element_text(hjust = 0.5, size = 12),plot.title = element_text(hjust = 0.5, size=12))+
      facet_wrap(facets = "nutrient", labeller = labeller(nutrient = nutrient.labs)
      ) + theme_classic() +
      coord_cartesian(ylim = c(0, max_value))      
    
    
    #  Kruskal wallis test followed by Dunn post hoc test
    # Test computation 
    res.kruskal <- good_table_3 %>% kruskal_test(mass ~ Inoculum)
    res.kruskal
    
    # Dunn posthoc Pairwise comparisons
    pwc <- good_table_3 %>% 
      # group_by(growth) %>%
      dunn_test(mass ~ Inoculum, p.adjust.method = "BH")
    
    # Generating letters for kruskal+dunn pairwise comparisons
    
    # Make fake ANOVA TUKEYHSD test only to replace the non parametric p-values
    tukey_values= data.frame()
    fit=aov(data=good_table_3,mass ~ Inoculum)
    anova(fit)
    res=TukeyHSD(fit)
    
    # Replace the adjusted pvalues from Kruskal+Dunn
    res[[1]][,4]=pwc$p.adj
    Tukey.levels <- res[[1]][,4]
    Labels_pairwise <- multcompLetters(Tukey.levels)['Letters']
    Inoculum <- names(Labels_pairwise[['Letters']])
    
    boxplot.df <- ddply(good_table_3, .(Inoculum), function (x) fivenum(max_scale$Max[max_scale$Host == FULL_NAMES[i+1]]*0.95))
    
    # Create a data frame out of the factor levels and Tukey's homogenous group letters
    plot.levels <- data.frame(Inoculum, labels = Labels_pairwise[['Letters']],
                              stringsAsFactors = FALSE)
    
    # Merge it with the labels
    labels.df <- merge(plot.levels, boxplot.df, by = "Inoculum" , sort = FALSE)
    
    p1=bxpot+
      geom_text(data = labels.df, aes(x = Inoculum, y = V1, label = labels))+
      geom_jitter(position=position_jitter(0.2)) +
      scale_shape_manual(values = c(15))
    
    Plot_list[[2*i+j]]=p1+rremove("x.ticks")
    
    
  }
}

#R1 data
Plot_list_R1=list()

for (i in 0:2) {
  
  good_table_3=subset (x = good_table_R1, subset = plant==PLANTS[i+1])
  colnames(good_table_3)[colnames(good_table_3) == "inoculum"] <- "SynCom"
  good_table_3$SynCom <- factor(good_table_3$SynCom, levels = c("SSC","AtSC", "LjSC", "HvSC","NS"))
  good_table_3$Experiment <- "R1"
  
  nutrient.labs <- paste(NUTRIENT[1],"nutrient")
  names(nutrient.labs) <-NUTRIENT[1]
  max <- max_scale$Max[max_scale$Host == FULL_NAMES[i+1]]
  
  colnames(good_table_3)[colnames(good_table_3) == "SynCom"] <- "Inoculum"
  good_table_3$Inoculum <- factor(good_table_3$Inoculum, levels = c("AtSC", "HvSC", "LjSC", "SSC", "NS"))
  
  # Visualization the data using box plots. Plot weight by groups.
  
  bxpot <- ggboxplot(
    data = good_table_3,
    x = "Inoculum", 
    y = "mass",
    combine = FALSE,
    merge = FALSE,
    legend="left",
    fill =  "Inoculum",
    color = "black",
    palette = c("#A3A500","#00B0F6","#00BF7D","#F8766D","white"),
    title = paste("R1 - ",FULL_NAMES[i+1],sep = ""),
    font.title=c(14,"italic"),
    xlab = "Inoculum",
    ylab = "Shoot mass (g)",
    bxp.errorbar = FALSE,
    bxp.errorbar.width = 0.4,
    facet.by = NULL,
    scales="free",
    panel.labs = NULL,
    short.panel.labs = TRUE,
    linetype = "solid",
    size = NULL,
    width = 0.8,
    notch = FALSE,
    outlier.shape = NA,
    select = NULL,
    remove = NULL,
    order = NULL,
    add = "none",
    add.params = list(),
    error.plot = "pointrange",
    label = NULL,
    font.label = list(size = 11, color = "black"),
    label.select = NULL,
    repel = FALSE,
    label.rectangle = FALSE,
    ggtheme = theme_pubr()
  ) +theme(axis.text.x = element_blank(),axis.title.x = element_blank(), title = element_text(hjust = 0.5, size = 12),plot.title = element_text(hjust = 0.5, size=12))+
    facet_wrap(facets = "nutrient", labeller = labeller(nutrient = nutrient.labs)
    ) + theme_classic() +
    scale_y_continuous(limits = c(0, max, na.rm = TRUE))
  
  #  Kruskal wallis test followed by Dunn post hoc test
  # Test computation 
  res.kruskal <- good_table_3 %>% kruskal_test(mass ~ Inoculum)
  res.kruskal
  
  # Dunn posthoc Pairwise comparisons
  pwc <- good_table_3 %>% 
    # group_by(growth) %>%
    dunn_test(mass ~ Inoculum, p.adjust.method = "BH")
  
  # Generating letters for kruskal+dunn pairwise comparisons
  
  # Make fake ANOVA TUKEYHSD test only to replace the non parametric p-values
  tukey_values= data.frame()
  fit=aov(data=good_table_3,mass ~ Inoculum)
  anova(fit)
  res=TukeyHSD(fit)
  
  # Replace the adjusted pvalues from Kruskal+Dunn
  res[[1]][,4]=pwc$p.adj
  Tukey.levels <- res[[1]][,4]
  Labels_pairwise <- multcompLetters(Tukey.levels)['Letters']
  Inoculum <- names(Labels_pairwise[['Letters']])
  
  boxplot.df <- ddply(good_table_3, .(Inoculum), function (x) fivenum(max_scale$Max[max_scale$Host == FULL_NAMES[i+1]]*0.95))
  
  # Create a data frame out of the factor levels and Tukey's homogenous group letters
  plot.levels <- data.frame(Inoculum, labels = Labels_pairwise[['Letters']],
                            stringsAsFactors = FALSE)
  
  # Merge it with the labels
  labels.df <- merge(plot.levels, boxplot.df, by = "Inoculum" , sort = FALSE)
  
  
  p1=bxpot+
    geom_text(data = labels.df, aes(x = Inoculum, y = V1, label = labels))+
    geom_jitter(position=position_jitter(0.2)) +
    scale_shape_manual(values = c(15))
  
  
  Plot_list_R1[[i+1]]=p1+rremove("x.ticks")
  
}

plot_shoots <- ggarrange(Plot_list[[1]],Plot_list[[3]]+rremove("ylab"),Plot_list[[5]]+rremove("ylab"),Plot_list[[2]],Plot_list[[4]]+rremove("ylab"),Plot_list[[6]]+rremove("ylab"), 
                         Plot_list_R1[[1]],Plot_list_R1[[2]]+rremove("ylab"),Plot_list_R1[[3]]+rremove("ylab"), ncol = 3,nrow = 3, common.legend = T, legend = "right" )

pdf(paste(results.dir,"Figure_S3_Shoot_weights.pdf", sep=""), width=15, height=15)
print(plot_shoots)
dev.off()
###Figure S4 - Nodule Numbers =====
importdat_nod <- read.table(file = paste(working_directory,"SSC_nodules.txt", sep = ""), sep = "\t", header=TRUE, dec = ".")

# Manipulate dataframe so it can be processed by ggplot2 package
nod_good_table=reshape2::melt(importdat_nod, measure.vars = c("Pink","White","Total"), value.name = "number", id.vars = c("condition","replicate"), variable.name = "type")
nod_good_table=mutate(nod_good_table,inoculum=word(condition,start = 2, sep=fixed("_")))
nod_good_table=mutate(nod_good_table,nutrient=word(condition,start = -1, sep=fixed("_")))

nod_good_table_R1 <- nod_good_table[grepl(pattern = "R1_",nod_good_table$condition),]
new <- str_split(nod_good_table_R1$condition, pattern = "_")
nod_good_table_R1$inoculum <- data.table::transpose(new)[[3]]

nod_good_table_R2 <- nod_good_table[!grepl(pattern = "R1_",nod_good_table$condition),]
nod_good_table_2 <- rbind(nod_good_table_R2,nod_good_table_R1)

nod_good_table_2$inoc_nut=paste0(nod_good_table_2$inoculum,"_",nod_good_table_2$nutrient)

# Remove NA values because of different number of samples being tested
nod_good_table_2=na.omit(nod_good_table_2)
nod_good_table_2$inoc_nut <- factor(nod_good_table_2$inoc_nut, levels = c("AtSC_Low", "HvSC_Low", "LjSC_Low", "SSC_Low","NS_Low", "AtSC_High","HvSC_High", "LjSC_High", "SSC_High","NS_High"))
nod_good_table_2$inoculum <- factor(nod_good_table_2$inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC","NS"))

# Visualize the data using box plots. Plot weight by groups.
nod_good_table_1=subset(x = nod_good_table_2, nod_good_table_2$type==c("Total"))

#  Kruskal wallis test followed by Dunn post hoc test
res.kruskal <- nod_good_table_1 %>% kruskal_test(number ~ inoc_nut)

# Dunn posthoc Pairwise comparisons
pwc <- nod_good_table_1 %>% 
  dunn_test(number ~ inoc_nut, p.adjust.method = "bonferroni") 

# Generating letters for kruskal+dunn pairwise comparisons
tukey_values= data.frame()
fit=aov(data=nod_good_table_1,number ~ inoc_nut)
anova(fit)
res=TukeyHSD(fit)
res[[1]][,4]=pwc$p.adj
Tukey.levels <- res[[1]][,4]
Labels_pairwise <- multcompLetters(Tukey.levels)['Letters']
inoc_nut <- names(Labels_pairwise[['Letters']])

boxplot.df <- ddply(nod_good_table_1, .(inoc_nut), function (x) max(fivenum(x$number)+0.04*(max(x$number))))

# Create a data frame out of the factor levels and Tukey's homogenous group letters
plot.levels <- data.frame(inoc_nut, labels = Labels_pairwise[['Letters']],
                          stringsAsFactors = FALSE)

# Merge it with the labels
labels.df <- merge(plot.levels, boxplot.df, by= "inoc_nut", sort = FALSE)

nod_good_table_1$inoc_nut <- factor(nod_good_table_1$inoc_nut, levels = c("AtSC_Low", "HvSC_Low", "LjSC_Low", "SSC_Low","NS_Low", "AtSC_High","HvSC_High", "LjSC_High", "SSC_High","NS_High"))

nod_total <- ggplot(nod_good_table_1, aes(x = inoc_nut, y = number)) +
  geom_boxplot(notch = FALSE, size=1, aes(color = inoculum),outlier.shape = NA) +
  theme_pubr() + ylab(label = "Number of nodule")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
        axis.title.x = element_blank(),
        title = element_text(hjust = 0.5, size = 12),
        plot.title = element_text(hjust = 0.5, size = 12), 
        legend.position = "none")+
  ylim(0,18) +
  scale_color_manual(values = c("#A3A500","#00B0F6","#00BF7D","#F8766D","white"))+
  geom_jitter(position=position_jitter(0.2))+
  geom_text(data = labels.df,size=5, aes(x = inoc_nut, y = max(V1), label = labels))

nod_total

pdf(paste(results.dir,"Figure_S4b_Nodules.pdf", sep=""), width=8, height=6)
print(nod_total)
dev.off()

#Figure S4c - nodulators - Taxonomic profile - nodules
#OTU TABLE
norm_SSC=read.table(paste(working_directory,"Isolate_tables/Original/SSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
round_SSC=floor(x = norm_SSC)

#Taxonomy TABLE
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
colnames(samples_df)[6]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparation
#Set the OTU, TAX and sample data for making phyloseq object
OTU = otu_table(as.matrix(round_SSC),taxa_are_rows = TRUE)
# TAX = tax_table(tax_mat)
TAX = tax_table(as.matrix(tax_df_2))

#Sample subsetting
cond="NOD"

samples_df_sub <- subset(samples_df, samples_df$Compartment != "RZ")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "AM")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "ES")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Condition != "NP")

samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")
samples <- sample_data(samples_df_sub_2)

phylo_sub = phyloseq(OTU,TAX, samples)

subsetted_table <- otu_table(phylo_sub)
subsetted_table_long <- melt(subsetted_table)

Hank_the_normalizer <- function(df,group,amount){
  df_2 <- df %>% dplyr::group_by_at(group) %>% dplyr::summarise(total=sum(.data[[amount]]))
  df_3 <- df_2$total
  names(df_3) <- df_2[[group]]
  df$total <- df_3[as.character(df[[group]])]
  df$Rel <- df[[amount]] / df$total
  return(df)
}

subsetted_table_long_2 <- Hank_the_normalizer(subsetted_table_long,"Var2","value")
subsetted_table_long_2$value[subsetted_table_long_2$Rel < 0.0005] <- 0
subsetted_table_long_3 <- subsetted_table_long_2[1:3]
data_wide <- spread(subsetted_table_long_3, Var2, value)
row.names(data_wide) <- data_wide$Var1
data_wide_2 <- data_wide %>% dplyr::select (-Var1)

OTU = otu_table(as.matrix(data_wide_2),taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(tax_df_2))
samples_df_sub <- subset(samples_df, samples_df$Compartment != "RZ")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "AM")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "ES")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Condition != "NP")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Condition != "Input")

samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")
samples <- sample_data(samples_df_sub_2)
phylo_sub = phyloseq(OTU,TAX, samples)

# Transform to relative abundance
phylo_sub_RA <- microbiome::transform(phylo_sub, "compositional")

# Select the top 5 most abundant OTUs
top5 <- names(sort(taxa_sums(phylo_sub_RA), decreasing=TRUE))[1:4]

# Melt phyloseq object into a dataframe
data_melt <- psmelt(phylo_sub_RA)

# Create a new OTU column where non-top5 are grouped as "Other strains"
data_melt$OTU_grouped <- ifelse(
  data_melt$OTU %in% c("LjNodule214", "P1_H10", "P2_A12", "P2_D6"), 
  data_melt$OTU, 
  ifelse(!is.na(data_melt$genus) & data_melt$genus == "Mesorhizobium", "Other Mesorhizobium", "Other strains"))

# Manually set factor levels for ordered legend
otu_order <- c("LjNodule214", "P1_H10", "P2_A12", "P2_D6","Other Mesorhizobium", "Other strains")
data_melt$OTU_grouped <- factor(data_melt$OTU_grouped, levels = otu_order)

# Custom color palette (forcing "Other" to dark gray)
otu_colors <- c(
  "LjNodule214" = "#66E1D0",  # Blue
  "P1_H10" = "#00C1C8",       # Orange
  "P2_A12" = "#00AA95",       # Green
  "P2_D6" = "#00C18C",        # Red
  "Other Mesorhizobium" = "#00773E",  # Assign a unique color for Other Mesorhizobium
  "Other strains" = "darkgray"        # Forced dark gray
)

# Improved stacked barplot
bar <- ggplot(data_melt, aes(fill=OTU_grouped, y=Abundance, x=Biorep)) + 
  geom_bar(position="stack", stat="identity", colour = "darkgray", linewidth = 0.01) +  
  scale_fill_manual(values = otu_colors) +  # Apply custom colors
  ggtitle("Lotus nodule colonization profiles") + 
  theme_classic() +
  labs(y = "Relative abundance", fill = "OTU") +  
  labs(fill = "Isolate") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20), 
    axis.text.x = element_blank(),
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 12)
  )

# Add facet wrap by inoculum if needed
bar2 <- bar + facet_wrap(~Inoculum+Nutrient, scales = "free_x", nrow = 1) +
  theme(strip.text.x = element_text(size = 12))

bar2

pdf(paste(results.dir,"Figure_S4c_nodule_taxonomic_profile.pdf", sep=""), width=7, height=4.5)
print(bar2)
dev.off()


###Figure S5 - Rarefaction curves =====

# OTU TABLE 
norm_SSC =read.table(paste(working_directory,"Isolate_tables/Original/SSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
# Taxonomy TABLE 
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
# Samples TABLE 
samples_df = read.table(paste(working_directory,"SSC_R2_metadata.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
colnames(samples_df)[6]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Set the OTU, TAX and sample data for making phyloseq object
OTU = otu_table(as.matrix(norm_SSC),taxa_are_rows = TRUE)
# TAX = tax_table(tax_mat)
TAX = tax_table(as.matrix(tax_df_2))

cond="ES"

#Subset for root (endosphere) samples
samples_df_sub <- subset(samples_df, samples_df$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

samples_sub = sample_data(samples_df_sub_2)

phylo_sub = phyloseq(OTU,TAX, samples_sub)

racur_data=rarecurve(x =floor(t(phylo_sub@otu_table@.Data)), step = 250, label = F, tidy = T)
racur_data=merge.data.frame(x = racur_data, y = samples_df_sub, by.x = "Site", by.y = 0)
racur_data$Inoculum <- factor(racur_data$Inoculum, levels = c("SSC","AtSC", "LjSC", "HvSC"))
racur_data$Condition <- gsub("At", "Arabidopsis", racur_data$Condition)
racur_data$Condition <- gsub("Hv", "Barley", racur_data$Condition)
racur_data$Condition <- gsub("Lj", "Lotus", racur_data$Condition)

racu=ggplot(data = racur_data, aes(x = Sample, y = Species, color = Experiment))
racu_all <- racu+geom_point( size = 0.05) + xlim(0, 50000)+ 
  xlab("No. of pseudoaligned reads") + ylab("No of isolates") +
  facet_grid(facets = c("Condition","Inoculum")) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 25)) +
  theme(
    panel.grid.major = element_line(color = "grey80"),  
    panel.grid.minor = element_line(color = "grey90")   
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))
racu_all

pdf(paste(results.dir,"Figure_S5_rarefaction_curves.pdf", sep=""), width=8, height=6)
print(racu_all)
dev.off()

###Figure S6 - Taxonomic profile - SynCom color =====
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
samples_df = read.table(paste(working_directory,"SSC_R2_metadata.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
colnames(samples_df)[6]="Nutrient"

#Fuse and rename
samples_df$Sample=paste(samples_df$Inoculum,samples_df$Condition,samples_df$Compartment, samples_df$Nutrient, samples_df$Experiment, sep ="_")
samples_df$Plant_Inoculum_compartment=paste(samples_df$Condition, samples_df$Inoculum, samples_df$Compartment, sep ="_")
samples_df$Sample[samples_df$Sample == "SSC_Input_Input_Input_R1"] <- "SSC_Input_R1"
samples_df$Sample[samples_df$Sample == "AtSC_Input_Input_Input_R1"] <- "AtSC_Input_R1"
samples_df$Sample[samples_df$Sample == "HvSC_Input_Input_Input_R1"] <- "HvSC_Input_R1"
samples_df$Sample[samples_df$Sample == "LjSC_Input_Input_Input_R1"] <- "LjSC_Input_R1"
samples_df$Sample[samples_df$Sample == "SSC_Input_Input_Input_R2"] <- "SSC_Input_R2"
samples_df$Sample[samples_df$Sample == "AtSC_Input_Input_Input_R2"] <- "AtSC_Input_R2"
samples_df$Sample[samples_df$Sample == "HvSC_Input_Input_Input_R2"] <- "HvSC_Input_R2"
samples_df$Sample[samples_df$Sample == "LjSC_Input_Input_Input_R2"] <- "LjSC_Input_R2"

# TAX = tax_table(tax_mat)
TAX = tax_table(as.matrix(tax_df_2))
samples = sample_data(samples_df)

SynComs <- c("AtSC", "LjSC", "HvSC", "SSC")

barplot_df_4 <- data.frame()

for (inoculum in SynComs){
  #Set the OTU, TAX and sample data for making phyloseq object
  norm_otu=read.table(paste(working_directory,"Isolate_tables/Original/SSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
  OTU = otu_table(as.matrix(norm_otu),taxa_are_rows = TRUE)
  
  samples_df_2 <-subset(samples_df, samples_df$Inoculum == paste(inoculum))
  samples_df_inp <- subset(samples_df_2, samples_df_2$Compartment == "Input")
  samples_df_nod <- subset(samples_df_2, samples_df_2$Compartment == "NOD")
  samples_df_es <- subset(samples_df_2, samples_df_2$Compartment == "ES")
  samples_df_sub <- rbind(samples_df_inp,samples_df_nod, samples_df_es)
  
  df <- as.data.frame(lapply(sample_data(samples_df_sub),function (y) if(class(y)!="factor" ) as.factor(y) else y),stringsAsFactors=T)
  row.names(df) <- row.names(samples_df_sub)
  
  samples_sub = sample_data(df)
  phylo_inoculum = phyloseq(OTU,TAX, samples_sub)
  
  #order taxa by abundance
  top100 <- names(sort(taxa_sums(phylo_inoculum), decreasing=TRUE))[1:1000]
  phylo_fraction <- prune_taxa(taxa = top100, x = phylo_inoculum)
  
  phylo_fraction_2 <- merge_samples(phylo_fraction, group = "Sample")
  Lalala=t(as.data.frame(phylo_fraction_2@otu_table))
  
  lololo=merge.data.frame(x = Lalala, tax_df, by = 0 )
  rownames(lololo)=lololo$Row.names
  lololo=lololo[,-1]
  
  barplot_df=melt(data = lololo)
  Hank_the_normalizer <- function(df,group,amount){
    df_2 <- df %>% dplyr::group_by_at(group) %>% dplyr::summarise(total=sum(.data[[amount]]))
    df_3 <- df_2$total
    names(df_3) <- df_2[[group]]
    df$total <- df_3[as.character(df[[group]])]
    df$Rel <- df[[amount]] / df$total
    return(df)
  }
  barplot_df_2 <- Hank_the_normalizer(barplot_df,"variable","value")
  barplot_df_3 <- cbind(barplot_df_2,str_split_fixed(barplot_df_2$variable,'_',5))
  colnames(barplot_df_3)[colnames(barplot_df_3) == "1"] <- "SynCom_2"
  colnames(barplot_df_3)[colnames(barplot_df_3) == "2"] <- "Plant"
  colnames(barplot_df_3)[colnames(barplot_df_3) == "3"] <- "Compartment"
  colnames(barplot_df_3)[colnames(barplot_df_3) == "4"] <- "Nutrient"
  colnames(barplot_df_3)[colnames(barplot_df_3) == "5"] <- "Experiment"
  
  barplot_df_4 <- rbind(barplot_df_4,barplot_df_3)
}

barplot_df_5 <- barplot_df_4[barplot_df_4$Nutrient == "low",]
barplot_df_5_sub <- barplot_df_4[barplot_df_4$Plant == "Input",]
barplot_df_5 <- rbind(barplot_df_5, barplot_df_5_sub)
barplot_df_5$variable <- factor(barplot_df_5$variable, levels = c("AtSC_Input_R1", "AtSC_Input_R2", "AtSC_At_ES_low_R1", "AtSC_At_RZ_low_R2", "AtSC_Hv_ES_low_R1", "AtSC_Hv_ES_low_R2", "AtSC_Lj_ES_low_R1", "AtSC_Lj_ES_low_R2", "HvSC_Input_R1", "HvSC_Input_R2", "HvSC_At_ES_low_R1", "HvSC_At_ES_low_R2", "HvSC_Hv_ES_low_R1", "HvSC_Hv_ES_low_R2", "HvSC_Lj_ES_low_R1", "HvSC_Lj_ES_low_R2", "HvSC_Lj_NOD_low_R2", "LjSC_Input_R1", "LjSC_Input_R2", "LjSC_At_ES_low_R1", "LjSC_At_ES_low_R2", "LjSC_Hv_ES_low_R1", "LjSC_Hv_ES_low_R2", "LjSC_Lj_ES_low_R1", "LjSC_Lj_ES_low_R2", "LjSC_Lj_NOD_low_R2", "SSC_Input_R1", "SSC_Input_R2", "SSC_At_ES_low_R1", "SSC_At_ES_low_R2", "SSC_Hv_ES_low_R1", "SSC_Hv_ES_low_R2", "SSC_Lj_ES_low_R1", "SSC_Lj_ES_low_R2", "SSC_Lj_NOD_low_R2"))

barplot_df_6 <- barplot_df_4[barplot_df_4$Nutrient == "high",]
barplot_df_6_sub <- barplot_df_4[barplot_df_4$Plant == "Input",]
barplot_df_6 <- rbind(barplot_df_6, barplot_df_6_sub)
barplot_df_6$variable <- factor(barplot_df_6$variable, levels = c("AtSC_Input_R1", "AtSC_Input_R2", "AtSC_At_ES_high_R2", "AtSC_Hv_ES_high_R2", "AtSC_Lj_ES_high_R2", "HvSC_Input_R1", "HvSC_Input_R2", "HvSC_At_ES_high_R2", "HvSC_Hv_ES_high_R2", "HvSC_Lj_ES_high_R2", "HvSC_Lj_NOD_high_R2", "LjSC_Input_R1", "LjSC_Input_R2", "LjSC_At_ES_high_R2", "LjSC_Hv_ES_high_R2", "LjSC_Lj_ES_high_R2", "LjSC_Lj_NOD_high_R2", "SSC_Input_R1", "SSC_Input_R2", "SSC_At_ES_high_R2", "SSC_Hv_ES_high_R2", "SSC_Lj_ES_high_R2", "SSC_Lj_NOD_high_R2"))
barplot_df_6$Plant <- factor(barplot_df_6$Plant, levels = c("Input", "At", "Hv", "Lj"))

bar <- ggplot(barplot_df_5, aes(fill=SynCom, y=Rel, x=variable)) + 
  geom_bar(position="stack", stat="identity") +  ggtitle("Taxonomic profile") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_classic() +
  labs(x ="Sample", y = "Relative abundance - low nutrient samples", fill = "Inoculum") +
  theme(panel.background=element_blank(),panel.grid=element_blank(),axis.line.x=element_line(size=.5, colour="black"),axis.line.y=element_line(size=.5, colour="black"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black", size=7),legend.position="right",legend.background=element_blank(),legend.key=element_blank(),legend.text= element_text(size=10),text=element_text(family="sans", size=10))+
  theme(axis.text.x = element_text(size = 14, angle = 25,hjust=1),axis.title.x = element_blank(), axis.title.y = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=18))
bar_2 <- bar+facet_wrap(~SynCom_2, scales = "free", nrow = 1 ) +theme(strip.text.x = element_text(size = 18))
bar_2

bar_3 <- ggplot(barplot_df_6, aes(fill=SynCom, y=Rel, x=variable)) + 
  geom_bar(position="stack", stat="identity") +  ggtitle("Taxonomic profile") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_classic() +
  labs(x ="Sample", y = "Relative abundance - high nutrient samples", fill = "Inoculum") +
  theme(panel.background=element_blank(),panel.grid=element_blank(),axis.line.x=element_line(size=.5, colour="black"),axis.line.y=element_line(size=.5, colour="black"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black", size=7),legend.position="right",legend.background=element_blank(),legend.key=element_blank(),legend.text= element_text(size=10),text=element_text(family="sans", size=10))+
  theme(axis.text.x = element_text(size = 14, angle = 25,hjust=1), axis.title.x = element_blank(), axis.title.y = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=18))
bar_4 <- bar_3+facet_wrap(~SynCom_2, scales = "free", nrow = 1 ) + theme(strip.text.x = element_text(size = 18))
bar_4

bar_5 <- ggarrange(bar_2, bar_4, ncol =1, nrow =2, common.legend = T, legend = "right")

bar_5

pdf(paste(results.dir,"Figure_S6_Tax_profile_cont.pdf", sep=""), width=35, height=15)
print(bar_5)
dev.off()



###Figure 2bcde & S7 - Alpha Diversity plots =====
#OTU TABLE
norm_SSC=read.table(paste(working_directory,"Isolate_tables/Original/SSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
round_SSC=floor(x = norm_SSC)

#Taxonomy TABLE
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
colnames(samples_df)[6]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparation
#Set the OTU, TAX and sample data for making phyloseq object
OTU = otu_table(as.matrix(round_SSC),taxa_are_rows = TRUE)
# TAX = tax_table(tax_mat)
TAX = tax_table(as.matrix(tax_df_2))

#Sample subsetting
cond="ES"

samples_df_sub <- subset(samples_df, samples_df$Compartment != "RZ")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "AM")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "NOD")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Condition != "NP")

samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")
samples <- sample_data(samples_df_sub_2)

phylo_sub = phyloseq(OTU,TAX, samples)

subsetted_table <- otu_table(phylo_sub)
subsetted_table_long <- melt(subsetted_table)

Hank_the_normalizer <- function(df,group,amount){
  df_2 <- df %>% dplyr::group_by_at(group) %>% dplyr::summarise(total=sum(.data[[amount]]))
  df_3 <- df_2$total
  names(df_3) <- df_2[[group]]
  df$total <- df_3[as.character(df[[group]])]
  df$Rel <- df[[amount]] / df$total
  return(df)
}

subsetted_table_long_2 <- Hank_the_normalizer(subsetted_table_long,"Var2","value")
subsetted_table_long_2$value[subsetted_table_long_2$Rel < 0.0005] <- 0
subsetted_table_long_3 <- subsetted_table_long_2[1:3]
data_wide <- spread(subsetted_table_long_3, Var2, value)
row.names(data_wide) <- data_wide$Var1
data_wide_2 <- data_wide %>% dplyr::select (-Var1)

OTU = otu_table(as.matrix(data_wide_2),taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(tax_df_2))
samples_df_sub <- subset(samples_df, samples_df$Compartment != "RZ")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "AM")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "NOD")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Condition != "NP")

samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")
samples <- sample_data(samples_df_sub_2)
phylo_sub = phyloseq(OTU,TAX, samples)

#Observed strains
method="Observed"

obs_df=estimate_richness(physeq = phylo_sub, measures = method)
obs_df=merge(x = obs_df, y = samples_df_sub_2, by = "row.names" )

obs_df=estimate_richness(physeq = phylo_sub, measures = method)
obs_df=merge(x = obs_df, y = samples_df_sub_2, by = "row.names" )
obs_df$inoculum_experiment=paste(obs_df$Condition,obs_df$Inoculum, sep ="_")

obs_df$Inoculum <- factor(obs_df$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC"))
obs_df$Condition <- factor(obs_df$Condition, levels = c("Input","At","Hv", "Lj"))
obs_df$Nutrient <- factor(obs_df$Nutrient, levels = c("high","low", "Input"))

colnames(obs_df)[colnames(obs_df) == "Condition"] <- "Host"

plants <- unique(obs_df$Host)
custom_labels_plant <- c(Input="Input", At = "Arabidopsis", Hv = "Barley", Lj = "Lotus")

# Initialize an empty list for storing results
anova_results <- list()
# Loop through each plant type to perform ANOVA, Tukey's test, and get letters
for(plant in plants) {
  # Subset data for the current plant
  subset_data <- obs_df[obs_df$Host == plant, ]
  
  # Perform ANOVA
  fitAnova <- aov(Observed ~ Inoculum, data=subset_data)
  
  # Perform Tukey's post-hoc test
  Tukey <- TukeyHSD(fitAnova)
  
  # Get letters
  letters_anova <- multcompView::multcompLetters4(fitAnova, Tukey)$Inoculum$Letters
  
  # Store results
  anova_results[[plant]] <- letters_anova
}
# Combine results into a data frame for plotting
ltlbl_combined <- do.call(rbind, lapply(names(anova_results), function(plant) {
  data.frame(Host = plant, Inoculum = names(anova_results[[plant]]), Letters = anova_results[[plant]])
}))
# Order factors based on your original setup
ltlbl_combined$Inoculum <- factor(ltlbl_combined$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC"))
ltlbl_combined$Host <- factor(ltlbl_combined$Host, levels = c("Input", "At", "Hv", "Lj"))

ltlbl_combined <- ltlbl_combined[order(ltlbl_combined$Host, ltlbl_combined$Inoculum), ]

alpha_tax_plot_obs_1 <- ggplot(data = obs_df, aes(x=Inoculum, y = Observed, color = Inoculum)) +
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+scale_color_manual(values = c("#A3A500","#00B0F6","#00BF7D","#F8766D", "black","gray70","black"))+
  scale_shape_manual(values = c(0,3))+
  geom_jitter(position=position_jitter(0.2), size =0.5, aes(color =Nutrient, shape=Experiment))+
  facet_wrap(~Host, scales="free_x", nrow=1, labeller=as_labeller(custom_labels_plant)) +
  theme(
    axis.text.x=element_blank(), 
    axis.title.x=element_blank(), 
    title=element_text(hjust=0.5, size=15), 
    axis.ticks.x=element_blank(),
    strip.background=element_rect(colour="gray50", size=0.3), # Change 'size' for thickness
    axis.text=element_text(color="gray50"),
    axis.line = element_line(color="gray50", size=0.3)
  ) +
  ylab(paste(method,"isolates"))+
  stat_summary(geom = 'text', label = ltlbl_combined$Letters, fun.y = max, aes(y = max(Observed)*1.05), show.legend=FALSE)

alpha_tax_plot_obs_1

#Plant subset
inocs <- unique(obs_df$Inoculum)

# Initialize an empty list for storing results
anova_results <- list()
# Loop through each plant type to perform ANOVA, Tukey's test, and get letters
for(inoc in inocs) {
  # Subset data for the current inoc
  subset_data <- obs_df[obs_df$Inoculum == inoc, ]
  
  # Perform ANOVA
  fitAnova <- aov(Observed ~ Host, data=subset_data)
  
  # Perform Tukey's post-hoc test
  Tukey <- TukeyHSD(fitAnova)
  
  # Get letters
  letters_anova <- multcompView::multcompLetters4(fitAnova, Tukey)$Host$Letters
  
  # Store results
  anova_results[[inoc]] <- letters_anova
}
# Combine results into a data frame for plotting
ltlbl_combined_2 <- do.call(rbind, lapply(names(anova_results), function(inoc) {
  data.frame(Inoculum = inoc, Host = names(anova_results[[inoc]]) , Letters = anova_results[[inoc]])
}))


# Order factors based on your original setup
ltlbl_combined_2$Inoculum <- factor(ltlbl_combined_2$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC"))
ltlbl_combined_2$Host <- factor(ltlbl_combined_2$Host, levels = c("Input", "At", "Hv", "Lj"))

#Alpha diversity plot, facet by plant/Input
alpha_tax_plot_obs_2 <- ggplot(data = obs_df, aes(x=Host, y = Observed, color = Host)) +
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+scale_color_manual(values = c("#1b9e77","#d95f02", "#e7298a","#7570b3","black", "gray70"))+
  scale_shape_manual(values = c(0,3))+
  geom_jitter(position=position_jitter(0.2), size =0.5, aes(color =Nutrient, shape=Experiment))+
  facet_wrap(~Inoculum, scales="free_x", nrow=1) +
  theme(
    axis.text.x=element_blank(), 
    axis.title.x=element_blank(), 
    title=element_text(hjust=0.5, size=15), 
    axis.ticks.x=element_blank(),
    strip.background=element_rect(colour="gray50", size=0.3), # Change 'size' for thickness
    axis.text=element_text(color="gray50"),
    axis.line = element_line(color="gray50", size=0.3)
  ) +
  ylab(paste(method,"isolates"))+
  stat_summary(geom = 'text', label = ltlbl_combined_2$Letters, fun.y = max, aes(y = max(Observed)*1.05), show.legend=FALSE)
alpha_tax_plot_obs_2

#Shannon diversity analysis
method="Shannon"

shannon_df=estimate_richness(physeq = phylo_sub, measures = method)
shannon_df=merge(x = shannon_df, y = samples_df_sub_2, by = "row.names" )

shannon_df=estimate_richness(physeq = phylo_sub, measures = method)
shannon_df=merge(x = shannon_df, y = samples_df_sub_2, by = "row.names" )
colnames(shannon_df)[colnames(shannon_df) == "Condition"] <- "Host"
shannon_df$inoculum_experiment=paste(shannon_df$Host,shannon_df$Inoculum, sep ="_")

shannon_df$Inoculum <- factor(shannon_df$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC"))
shannon_df$Host <- factor(shannon_df$Host, levels = c("Input","At","Hv", "Lj"))
shannon_df$Nutrient <- factor(shannon_df$Nutrient, levels = c("high","low", "Input"))

plants <- unique(shannon_df$Host)

# Initialize an empty list for storing results
anova_results <- list()
# Loop through each plant type to perform ANOVA, Tukey's test, and get letters
for(plant in plants) {
  # Subset data for the current plant
  subset_data <- shannon_df[shannon_df$Host == plant, ]
  
  # Perform ANOVA
  fitAnova <- aov(Shannon ~ Inoculum, data=subset_data)
  
  # Perform Tukey's post-hoc test
  Tukey <- TukeyHSD(fitAnova)
  
  # Get letters
  letters_anova <- multcompView::multcompLetters4(fitAnova, Tukey)$Inoculum$Letters
  
  # Store results
  anova_results[[plant]] <- letters_anova
}
# Combine results into a data frame for plotting
ltlbl_combined <- do.call(rbind, lapply(names(anova_results), function(plant) {
  data.frame(Host = plant, Inoculum = names(anova_results[[plant]]), Letters = anova_results[[plant]])
}))
# Order factors based on your original setup
ltlbl_combined$Inoculum <- factor(ltlbl_combined$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC"))
ltlbl_combined$Host <- factor(ltlbl_combined$Host, levels = c("Input", "At", "Hv", "Lj"))

ltlbl_combined <- ltlbl_combined[order(ltlbl_combined$Host, ltlbl_combined$Inoculum), ]

#Alpha diversity plot, facet by plant/Input
alpha_tax_plot_shan_1 <- ggplot(data = shannon_df, aes(x=Inoculum, y = Shannon, color = Inoculum)) +
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+scale_color_manual(values = c("#A3A500","#00B0F6","#00BF7D","#F8766D", "black","gray70","black"))+
  scale_shape_manual(values = c(0,3))+
  geom_jitter(position=position_jitter(0.2), size =0.5, aes(color =Nutrient, shape=Experiment))+
  facet_wrap(~Host, scales="free_x", nrow=1, labeller=as_labeller(custom_labels_plant)) +
  theme(
    axis.text.x=element_blank(), 
    axis.title.x=element_blank(), 
    title=element_text(hjust=0.5, size=15), 
    axis.ticks.x=element_blank(),
    strip.background=element_rect(colour="gray50", size=0.3), # Change 'size' for thickness
    axis.text=element_text(color="gray50"),
    axis.line = element_line(color="gray50", size=0.3)
  ) +
  ylab(paste(method,"isolates"))+
  stat_summary(geom = 'text', label = ltlbl_combined$Letters, fun.y = max, aes(y = max(Shannon)*1.05), show.legend=FALSE)

alpha_tax_plot_shan_1

#Alpha diversity plot, facet by Syncom

#Plant subset
inocs <- unique(shannon_df$Inoculum)

# Initialize an empty list for storing results
anova_results <- list()
# Loop through each plant type to perform ANOVA, Tukey's test, and get letters
for(inoc in inocs) {
  # Subset data for the current inoc
  subset_data <- shannon_df[shannon_df$Inoculum == inoc, ]
  
  # Perform ANOVA
  fitAnova <- aov(Shannon ~ Host, data=subset_data)
  
  # Perform Tukey's post-hoc test
  Tukey <- TukeyHSD(fitAnova)
  
  # Get letters
  letters_anova <- multcompView::multcompLetters4(fitAnova, Tukey)$Host$Letters
  
  # Store results
  anova_results[[inoc]] <- letters_anova
}
# Combine results into a data frame for plotting
ltlbl_combined_2 <- do.call(rbind, lapply(names(anova_results), function(inoc) {
  data.frame(Inoculum = inoc, Host = names(anova_results[[inoc]]) , Letters = anova_results[[inoc]])
}))

# Order factors based on your original setup
ltlbl_combined_2$Inoculum <- factor(ltlbl_combined_2$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC"))
ltlbl_combined_2$Host <- factor(ltlbl_combined_2$Host, levels = c("Input", "At", "Hv", "Lj"))

alpha_tax_plot_shan_2 <- ggplot(data = shannon_df, aes(x=Host, y = Shannon, color = Host)) +
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+scale_color_manual(values = c("#1b9e77","#d95f02", "#e7298a","#7570b3", "black", "gray70"))+
  scale_shape_manual(values = c(0,3))+
  geom_jitter(position=position_jitter(0.2), size =0.5, aes(color =Nutrient, shape=Experiment))+
  facet_wrap(~Inoculum, scales="free_x", nrow=1) +
  theme(
    axis.text.x=element_blank(), 
    axis.title.x=element_blank(), 
    title=element_text(hjust=0.5, size=15), 
    axis.ticks.x=element_blank(),
    strip.background=element_rect(colour="gray50", size=0.3), # Change 'size' for thickness
    axis.text=element_text(color="gray50"),
    axis.line = element_line(color="gray50", size=0.3)
  ) +
  ylab(paste(method,"isolates"))+
  stat_summary(geom = 'text', label = ltlbl_combined_2$Letters, fun.y = max, aes(y = max(Shannon)*1.05), show.legend=FALSE)
alpha_tax_plot_shan_2

#KO OTU TABLE
norm_SSC_KO =read.table(paste(working_directory,"KO_tables/Original/SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
round_SSC_KO=floor(x = norm_SSC_KO)

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
colnames(samples_df)[6]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparation
#Set the OTU, TAX and sample data for making phyloseq object
OTU = otu_table(as.matrix(round_SSC_KO),taxa_are_rows = TRUE)

#Sample subsetting
cond="ES"

samples_df_sub <- subset(samples_df, samples_df$Compartment != "RZ")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "AM")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "NOD")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Condition != "NP")

samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")
samples <- sample_data(samples_df_sub_2)

round_SSC_KO_2 <- round_SSC_KO[,colnames(round_SSC_KO) %in% row.names(samples_df_sub_2)]
round_SSC_2 <- round_SSC[,colnames(round_SSC) %in% row.names(samples_df_sub_2)]

average_mult_tax_to_KO <- sum(colSums(round_SSC_KO_2)/colSums(round_SSC_2))/length(colSums(round_SSC_KO_2))

phylo_KO_sub = phyloseq(OTU,samples)

subsetted_table <- otu_table(phylo_KO_sub)
subsetted_table_long <- melt(subsetted_table)

Hank_the_normalizer <- function(df,group,amount){
  df_2 <- df %>% dplyr::group_by_at(group) %>% dplyr::summarise(total=sum(.data[[amount]]))
  df_3 <- df_2$total
  names(df_3) <- df_2[[group]]
  df$total <- df_3[as.character(df[[group]])]
  df$Rel <- df[[amount]] / df$total
  return(df)
}

subsetted_table_long_2 <- Hank_the_normalizer(subsetted_table_long,"Var2","value")
subsetted_table_long_2$value[subsetted_table_long_2$Rel < 0.0005/average_mult_tax_to_KO] <- 0
subsetted_table_long_3 <- subsetted_table_long_2[1:3]
data_wide <- spread(subsetted_table_long_3, Var2, value)
row.names(data_wide) <- data_wide$Var1
data_wide_2 <- data_wide %>% dplyr::select (-Var1)

OTU = otu_table(as.matrix(data_wide_2),taxa_are_rows = TRUE)
samples_df_sub <- subset(samples_df, samples_df$Compartment != "RZ")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "AM")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "NOD")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Condition != "NP")

samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")
samples <- sample_data(samples_df_sub_2)
phylo_KO_sub = phyloseq(OTU, samples)

#Observed KOs
method="Observed"

obs_df_KO=estimate_richness(physeq = phylo_KO_sub, measures = method)
obs_df_KO=merge(x = obs_df_KO, y = samples_df_sub_2, by = "row.names" )

obs_df_KO=estimate_richness(physeq = phylo_KO_sub, measures = method)
obs_df_KO=merge(x = obs_df_KO, y = samples_df_sub_2, by = "row.names" )
obs_df_KO$inoculum_experiment=paste(obs_df_KO$Condition,obs_df_KO$Inoculum, sep ="_")

obs_df_KO$Inoculum <- factor(obs_df_KO$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC"))
obs_df_KO$Condition <- factor(obs_df_KO$Condition, levels = c("Input","At","Hv", "Lj"))
obs_df_KO$Nutrient <- factor(obs_df_KO$Nutrient, levels = c("high","low", "Input"))

colnames(obs_df_KO)[colnames(obs_df_KO) == "Condition"] <- "Host"

plants <- unique(obs_df_KO$Host)

# Initialize an empty list for storing results
anova_results <- list()
# Loop through each plant type to perform ANOVA, Tukey's test, and get letters
for(plant in plants) {
  # Subset data for the current plant
  subset_data <- obs_df_KO[obs_df_KO$Host == plant, ]
  
  # Perform ANOVA
  fitAnova <- aov(Observed ~ Inoculum, data=subset_data)
  
  # Perform Tukey's post-hoc test
  Tukey <- TukeyHSD(fitAnova)
  
  # Get letters
  letters_anova <- multcompView::multcompLetters4(fitAnova, Tukey)$Inoculum$Letters
  
  # Store results
  anova_results[[plant]] <- letters_anova
}
# Combine results into a data frame for plotting
ltlbl_combined <- do.call(rbind, lapply(names(anova_results), function(plant) {
  data.frame(Host = plant, Inoculum = names(anova_results[[plant]]), Letters = anova_results[[plant]])
}))
# Order factors based on your original setup
ltlbl_combined$Inoculum <- factor(ltlbl_combined$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC"))
ltlbl_combined$Host <- factor(ltlbl_combined$Host, levels = c("Input", "At", "Hv", "Lj"))

ltlbl_combined <- ltlbl_combined[order(ltlbl_combined$Host, ltlbl_combined$Inoculum), ]

#Alpha diversity plot, facet by plant/Input
alpha_plot_KO_obs_1 <- ggplot(data = obs_df_KO, aes(x=Inoculum, y = Observed, color = Inoculum)) +
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+scale_color_manual(values = c("#A3A500","#00B0F6","#00BF7D","#F8766D", "black","gray70","black"))+
  scale_shape_manual(values = c(0,3))+
  geom_jitter(position=position_jitter(0.2), size =0.5, aes(color =Nutrient, shape=Experiment))+
  facet_wrap(~Host, scales="free_x", nrow=1, labeller=as_labeller(custom_labels_plant)) +
  theme(
    axis.text.x=element_blank(), 
    axis.title.x=element_blank(), 
    title=element_text(hjust=0.5, size=15), 
    axis.ticks.x=element_blank(),
    strip.background=element_rect(colour="gray50", size=0.3), # Change 'size' for thickness
    axis.text=element_text(color="gray50"),
    axis.line = element_line(color="gray50", size=0.3)
  ) +
  ylab(paste(method,"KOs"))+
  stat_summary(geom = 'text', label = ltlbl_combined$Letters, fun.y = max, aes(y = max(Observed)*1.05), show.legend=FALSE)

alpha_plot_KO_obs_1

#Alpha diversity plot, facet by Syncom
#Plant subset
inocs <- unique(obs_df_KO$Inoculum)

# Initialize an empty list for storing results
anova_results <- list()
# Loop through each plant type to perform ANOVA, Tukey's test, and get letters
for(inoc in inocs) {
  # Subset data for the current inoc
  subset_data <- obs_df_KO[obs_df_KO$Inoculum == inoc, ]
  
  # Perform ANOVA
  fitAnova <- aov(Observed ~ Host, data=subset_data)
  
  # Perform Tukey's post-hoc test
  Tukey <- TukeyHSD(fitAnova)
  
  # Get letters
  letters_anova <- multcompView::multcompLetters4(fitAnova, Tukey)$Host$Letters
  
  # Store results
  anova_results[[inoc]] <- letters_anova
}
# Combine results into a data frame for plotting
ltlbl_combined_2 <- do.call(rbind, lapply(names(anova_results), function(inoc) {
  data.frame(Inoculum = inoc, Host = names(anova_results[[inoc]]) , Letters = anova_results[[inoc]])
}))

# Order factors based on your original setup
ltlbl_combined_2$Inoculum <- factor(ltlbl_combined_2$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC"))
ltlbl_combined_2$Host <- factor(ltlbl_combined_2$Host, levels = c("Input", "At", "Hv", "Lj"))

alpha_plot_KO_obs_2 <- ggplot(data = obs_df_KO, aes(x=Host, y = Observed, color = Host)) +
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+scale_color_manual(values = c("#1b9e77", "#d95f02", "#e7298a","#7570b3","black", "gray70"))+
  scale_shape_manual(values = c(0,3))+
  geom_jitter(position=position_jitter(0.2), size =0.5, aes(color =Nutrient, shape=Experiment))+
  facet_wrap(~Inoculum, scales="free_x", nrow=1) +
  theme(
    axis.text.x=element_blank(), 
    axis.title.x=element_blank(), 
    title=element_text(hjust=0.5, size=15), 
    axis.ticks.x=element_blank(),
    strip.background=element_rect(colour="gray50", size=0.3), # Change 'size' for thickness
    axis.text=element_text(color="gray50"),
    axis.line = element_line(color="gray50", size=0.3)
  ) +
  ylab(paste(method,"KOs"))+
  stat_summary(geom = 'text', label = ltlbl_combined_2$Letters, fun.y = max, aes(y = max(Observed)*1.05), show.legend=FALSE)
alpha_plot_KO_obs_2

#Shannon diversity KOs
method="Shannon"

shannon_df_KO=estimate_richness(physeq = phylo_KO_sub, measures = method)
shannon_df_KO=merge(x = shannon_df_KO, y = samples_df_sub_2, by = "row.names" )
shannon_df_KO$inoculum_experiment=paste(shannon_df_KO$Condition,shannon_df_KO$Inoculum, sep ="_")

shannon_df_KO$Inoculum <- factor(shannon_df_KO$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC"))
shannon_df_KO$Condition <- factor(shannon_df_KO$Condition, levels = c("Input","At","Hv", "Lj"))
shannon_df_KO$Nutrient <- factor(shannon_df_KO$Nutrient, levels = c("high","low", "Input"))

colnames(shannon_df_KO)[colnames(shannon_df_KO) == "Condition"] <- "Host"

plants <- unique(shannon_df_KO$Host)

# Initialize an empty list for storing results
anova_results <- list()
# Loop through each plant type to perform ANOVA, Tukey's test, and get letters
for(plant in plants) {
  # Subset data for the current plant
  subset_data <- shannon_df_KO[shannon_df_KO$Host == plant, ]
  
  # Perform ANOVA
  fitAnova <- aov(Shannon ~ Inoculum, data=subset_data)
  
  # Perform Tukey's post-hoc test
  Tukey <- TukeyHSD(fitAnova)
  
  # Get letters
  letters_anova <- multcompView::multcompLetters4(fitAnova, Tukey)$Inoculum$Letters
  
  # Store results
  anova_results[[plant]] <- letters_anova
}
# Combine results into a data frame for plotting
ltlbl_combined <- do.call(rbind, lapply(names(anova_results), function(plant) {
  data.frame(Host = plant, Inoculum = names(anova_results[[plant]]), Letters = anova_results[[plant]])
}))
# Order factors based on your original setup
ltlbl_combined$Inoculum <- factor(ltlbl_combined$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC"))
ltlbl_combined$Host <- factor(ltlbl_combined$Host, levels = c("Input", "At", "Hv", "Lj"))

ltlbl_combined <- ltlbl_combined[order(ltlbl_combined$Host, ltlbl_combined$Inoculum), ]

#Alpha diversity plot, facet by plant/Input
alpha_plot_KO_shan_1 <- ggplot(data = shannon_df_KO, aes(x=Inoculum, y = Shannon, color = Inoculum)) +
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+scale_color_manual(values = c("#A3A500","#00B0F6","#00BF7D","#F8766D", "black","gray70","black"))+
  scale_shape_manual(values = c(0,3))+
  geom_jitter(position=position_jitter(0.2), size =0.5, aes(color =Nutrient, shape=Experiment))+
  facet_wrap(~Host, scales="free_x", nrow=1, labeller=as_labeller(custom_labels_plant)) +
  theme(
    axis.text.x=element_blank(), 
    axis.title.x=element_blank(), 
    title=element_text(hjust=0.5, size=15), 
    axis.ticks.x=element_blank(),
    strip.background=element_rect(colour="gray50", size=0.3), # Change 'size' for thickness
    axis.text=element_text(color="gray50"),
    axis.line = element_line(color="gray50", size=0.3)
  ) +
  ylab(paste(method,"KOs"))+
  stat_summary(geom = 'text', label = ltlbl_combined$Letters, fun.y = max, aes(y = max(Shannon)*1.05), show.legend=FALSE)

alpha_plot_KO_shan_1

#Alpha diversity plot, facet by Syncom
#Plant subset
inocs <- unique(shannon_df_KO$Inoculum)

# Initialize an empty list for storing results
anova_results <- list()
# Loop through each plant type to perform ANOVA, Tukey's test, and get letters
for(inoc in inocs) {
  # Subset data for the current inoc
  subset_data <- shannon_df_KO[shannon_df_KO$Inoculum == inoc, ]
  
  # Perform ANOVA
  fitAnova <- aov(Shannon ~ Host, data=subset_data)
  
  # Perform Tukey's post-hoc test
  Tukey <- TukeyHSD(fitAnova)
  
  # Get letters
  letters_anova <- multcompView::multcompLetters4(fitAnova, Tukey)$Host$Letters
  
  # Store results
  anova_results[[inoc]] <- letters_anova
}
# Combine results into a data frame for plotting
ltlbl_combined_2 <- do.call(rbind, lapply(names(anova_results), function(inoc) {
  data.frame(Inoculum = inoc, Host = names(anova_results[[inoc]]) , Letters = anova_results[[inoc]])
}))


# Order factors based on your original setup
ltlbl_combined_2$Inoculum <- factor(ltlbl_combined_2$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC"))
ltlbl_combined_2$Host <- factor(ltlbl_combined_2$Host, levels = c("Input", "At", "Hv", "Lj"))

alpha_plot_KO_shan_2 <- ggplot(data = shannon_df_KO, aes(x=Host, y = Shannon, color = Host)) +
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+scale_color_manual(values = c("#1b9e77", "#d95f02", "#e7298a","#7570b3","black", "gray70"))+
  scale_shape_manual(values = c(0,3))+
  geom_jitter(position=position_jitter(0.2), size =0.5, aes(color =Nutrient, shape=Experiment))+
  facet_wrap(~Inoculum, scales="free_x", nrow=1) +
  theme(
    axis.text.x=element_blank(), 
    axis.title.x=element_blank(), 
    title=element_text(hjust=0.5, size=15), 
    axis.ticks.x=element_blank(),
    strip.background=element_rect(colour="gray50", size=0.3), # Change 'size' for thickness
    axis.text=element_text(color="gray50"),
    axis.line = element_line(color="gray50", size=0.3)
  ) +
  ylab(paste(method,"KOs"))+
  stat_summary(geom = 'text', label = ltlbl_combined_2$Letters, fun.y = max, aes(y = max(Shannon)*1.05), show.legend=FALSE)
alpha_plot_KO_shan_2

#Printing all alpha diversity plots
pdf(paste(results.dir,"alpha_diversity_2b.pdf", sep=""), width=7, height=3)
print(alpha_tax_plot_obs_2)
dev.off()

pdf(paste(results.dir,"alpha_diversity_2c.pdf", sep=""), width=7, height=3)
print(alpha_tax_plot_obs_1)
dev.off()

pdf(paste(results.dir,"alpha_diversity_2d.pdf", sep=""), width=7, height=3)
print(alpha_plot_KO_obs_2)
dev.off()

pdf(paste(results.dir,"alpha_diversity_2e.pdf", sep=""), width=7, height=3)
print(alpha_plot_KO_obs_1)
dev.off()

pdf(paste(results.dir,"alpha_diversity_S7a.pdf", sep=""), width=7, height=3)
print(alpha_tax_plot_shan_2)
dev.off()

pdf(paste(results.dir,"alpha_diversity_S7b.pdf", sep=""), width=7, height=3)
print(alpha_tax_plot_shan_1)
dev.off()

pdf(paste(results.dir,"alpha_diversity_S7c.pdf", sep=""), width=7, height=3)
print(alpha_plot_KO_shan_2)
dev.off()

pdf(paste(results.dir,"alpha_diversity_S7d.pdf", sep=""), width=7, height=3)
print(alpha_plot_KO_shan_1)
dev.off()

###Figure S8 - Taxonomic profile - Genus color =====
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
samples_df = read.table(paste(working_directory,"SSC_R2_metadata.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
colnames(samples_df)[6]="Nutrient"

#Fuse and rename
samples_df$Sample=paste(samples_df$Inoculum,samples_df$Condition,samples_df$Compartment, samples_df$Nutrient, samples_df$Experiment, sep ="_")
samples_df$Plant_Inoculum_compartment=paste(samples_df$Condition, samples_df$Inoculum, samples_df$Compartment, sep ="_")
samples_df$Sample[samples_df$Sample == "SSC_Input_Input_Input_R1"] <- "SSC_Input_R1"
samples_df$Sample[samples_df$Sample == "AtSC_Input_Input_Input_R1"] <- "AtSC_Input_R1"
samples_df$Sample[samples_df$Sample == "HvSC_Input_Input_Input_R1"] <- "HvSC_Input_R1"
samples_df$Sample[samples_df$Sample == "LjSC_Input_Input_Input_R1"] <- "LjSC_Input_R1"
samples_df$Sample[samples_df$Sample == "SSC_Input_Input_Input_R2"] <- "SSC_Input_R2"
samples_df$Sample[samples_df$Sample == "AtSC_Input_Input_Input_R2"] <- "AtSC_Input_R2"
samples_df$Sample[samples_df$Sample == "HvSC_Input_Input_Input_R2"] <- "HvSC_Input_R2"
samples_df$Sample[samples_df$Sample == "LjSC_Input_Input_Input_R2"] <- "LjSC_Input_R2"

# TAX = tax_table(tax_mat)
TAX = tax_table(as.matrix(tax_df_2))
samples = sample_data(samples_df)

SynComs <- c("AtSC", "LjSC", "HvSC", "SSC")

barplot_df_4 <- data.frame()

for (inoculum in SynComs){
  #Set the OTU, TAX and sample data for making phyloseq object
  norm_otu=read.table(paste(working_directory,"Isolate_tables/Original/SSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
  OTU = otu_table(as.matrix(norm_otu),taxa_are_rows = TRUE)
  
  samples_df_2 <-subset(samples_df, samples_df$Inoculum == paste(inoculum))
  samples_df_inp <- subset(samples_df_2, samples_df_2$Compartment == "Input")
  samples_df_nod <- subset(samples_df_2, samples_df_2$Compartment == "NOD")
  samples_df_es <- subset(samples_df_2, samples_df_2$Compartment == "ES")
  samples_df_sub <- rbind(samples_df_inp,samples_df_nod, samples_df_es)
  
  df <- as.data.frame(lapply(sample_data(samples_df_sub),function (y) if(class(y)!="factor" ) as.factor(y) else y),stringsAsFactors=T)
  row.names(df) <- row.names(samples_df_sub)
  
  samples_sub = sample_data(df)
  phylo_inoculum = phyloseq(OTU,TAX, samples_sub)
  
  #order taxa by abundance
  top100 <- names(sort(taxa_sums(phylo_inoculum), decreasing=TRUE))[1:1000]
  phylo_fraction <- prune_taxa(taxa = top100, x = phylo_inoculum)
  
  phylo_fraction_2 <- merge_samples(phylo_fraction, group = "Sample")
  Lalala=t(as.data.frame(phylo_fraction_2@otu_table))
  
  lololo=merge.data.frame(x = Lalala, tax_df, by = 0 )
  rownames(lololo)=lololo$Row.names
  lololo=lololo[,-1]
  
  barplot_df=melt(data = lololo)
  Hank_the_normalizer <- function(df,group,amount){
    df_2 <- df %>% dplyr::group_by_at(group) %>% dplyr::summarise(total=sum(.data[[amount]]))
    df_3 <- df_2$total
    names(df_3) <- df_2[[group]]
    df$total <- df_3[as.character(df[[group]])]
    df$Rel <- df[[amount]] / df$total
    return(df)
  }
  barplot_df_2 <- Hank_the_normalizer(barplot_df,"variable","value")
  barplot_df_3 <- cbind(barplot_df_2,str_split_fixed(barplot_df_2$variable,'_',5))
  colnames(barplot_df_3)[colnames(barplot_df_3) == "1"] <- "SynCom_2"
  colnames(barplot_df_3)[colnames(barplot_df_3) == "2"] <- "Plant"
  colnames(barplot_df_3)[colnames(barplot_df_3) == "3"] <- "Compartment"
  colnames(barplot_df_3)[colnames(barplot_df_3) == "4"] <- "Nutrient"
  colnames(barplot_df_3)[colnames(barplot_df_3) == "5"] <- "Experiment"
  
  barplot_df_4 <- rbind(barplot_df_4,barplot_df_3)
}

barplot_df_5 <- barplot_df_4[barplot_df_4$Nutrient == "low",]
barplot_df_5_sub <- barplot_df_4[barplot_df_4$Plant == "Input",]
barplot_df_5 <- rbind(barplot_df_5, barplot_df_5_sub)
barplot_df_5$variable <- factor(barplot_df_5$variable, levels = c("AtSC_Input_R1", "AtSC_Input_R2", "AtSC_At_ES_low_R1", "AtSC_At_ES_low_R2", "AtSC_Hv_ES_low_R1", "AtSC_Hv_ES_low_R2", "AtSC_Lj_ES_low_R1", "AtSC_Lj_ES_low_R2", "HvSC_Input_R1", "HvSC_Input_R2", "HvSC_At_ES_low_R1", "HvSC_At_ES_low_R2", "HvSC_Hv_ES_low_R1", "HvSC_Hv_ES_low_R2", "HvSC_Lj_ES_low_R1", "HvSC_Lj_ES_low_R2", "HvSC_Lj_NOD_low_R2", "LjSC_Input_R1", "LjSC_Input_R2", "LjSC_At_ES_low_R1", "LjSC_At_ES_low_R2", "LjSC_Hv_ES_low_R1", "LjSC_Hv_ES_low_R2", "LjSC_Lj_ES_low_R1", "LjSC_Lj_ES_low_R2", "LjSC_Lj_NOD_low_R2", "SSC_Input_R1", "SSC_Input_R2", "SSC_At_ES_low_R1", "SSC_At_ES_low_R2", "SSC_Hv_ES_low_R1", "SSC_Hv_ES_low_R2", "SSC_Lj_ES_low_R1", "SSC_Lj_ES_low_R2", "SSC_Lj_NOD_low_R2"))

barplot_df_6 <- barplot_df_4[barplot_df_4$Nutrient == "high",]
barplot_df_6_sub <- barplot_df_4[barplot_df_4$Plant == "Input",]
barplot_df_6 <- rbind(barplot_df_6, barplot_df_6_sub)
barplot_df_6$variable <- factor(barplot_df_6$variable, levels = c("AtSC_Input_R1", "AtSC_Input_R2", "AtSC_At_ES_high_R2", "AtSC_Hv_ES_high_R2", "AtSC_Lj_ES_high_R2", "HvSC_Input_R1", "HvSC_Input_R2", "HvSC_At_ES_high_R2", "HvSC_Hv_ES_high_R2", "HvSC_Lj_ES_high_R2", "HvSC_Lj_NOD_high_R2", "LjSC_Input_R1", "LjSC_Input_R2", "LjSC_At_ES_high_R2", "LjSC_Hv_ES_high_R2", "LjSC_Lj_ES_high_R2", "LjSC_Lj_NOD_high_R2", "SSC_Input_R1", "SSC_Input_R2", "SSC_At_ES_high_R2", "SSC_Hv_ES_high_R2", "SSC_Lj_ES_high_R2", "SSC_Lj_NOD_high_R2"))
barplot_df_6$Plant <- factor(barplot_df_6$Plant, levels = c("Input", "At", "Hv", "Lj"))

#plot
bar <- ggplot(barplot_df_5, aes(fill=genus, y=Rel, x=variable)) + 
  geom_bar(position="stack", stat="identity") +  ggtitle("Taxonomic profile") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_classic() +
  labs(x ="Sample", y = "Relative abundance - low nutrient samples", fill = "Genus") +
  theme(panel.background=element_blank(),panel.grid=element_blank(),axis.line.x=element_line(size=.5, colour="black"),axis.line.y=element_line(size=.5, colour="black"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black", size=7),legend.position="right",legend.background=element_blank(),legend.key=element_blank(),legend.text= element_text(size=10),text=element_text(family="sans", size=10))+
  theme(axis.text.x = element_text(size = 14, angle = 25,hjust=1),axis.title.x = element_blank(), axis.title.y = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14, face = rep("italic")), plot.title = element_text(size=18))
bar_2 <- bar+facet_wrap(~SynCom_2, scales = "free", nrow = 1 ) +theme(strip.text.x = element_text(size = 18))
bar_2

bar_3 <- ggplot(barplot_df_6, aes(fill=genus, y=Rel, x=variable)) + 
  geom_bar(position="stack", stat="identity") +  ggtitle("Taxonomic profile") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_classic() +
  labs(x ="Sample", y = "Relative abundance - high nutrient samples", fill = "Genus") +
  theme(panel.background=element_blank(),panel.grid=element_blank(),axis.line.x=element_line(size=.5, colour="black"),axis.line.y=element_line(size=.5, colour="black"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black", size=7),legend.position="right",legend.background=element_blank(),legend.key=element_blank(),legend.text= element_text(size=10),text=element_text(family="sans", size=10))+
  theme(axis.text.x = element_text(size = 14, angle = 25,hjust=1), axis.title.x = element_blank(), axis.title.y = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14,face = rep("italic")), plot.title = element_text(size=18))
bar_4 <- bar_3+facet_wrap(~SynCom_2, scales = "free", nrow = 1 ) + theme(strip.text.x = element_text(size = 18))
bar_4

bar_5 <- ggarrange(bar_2, bar_4, ncol =1, nrow =2, common.legend = T, legend = "right")

bar_5

pdf(paste(results.dir,"Figure_S8_Tax_profile.pdf", sep=""), width=35, height=15)
print(bar_5)
dev.off()


###Figure S9 - Relative abundances dominators & nodulators in nodules =====

#Figure S9 Dominance by rhizobacter and mesorhizobium
#OTU TABLE
norm_SSC=read.table(paste(working_directory,"Isolate_tables/Original/SSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
round_SSC=floor(x = norm_SSC)

#Taxonomy TABLE
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
colnames(samples_df)[6]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparation
#Set the OTU, TAX and sample data for making phyloseq object
OTU = otu_table(as.matrix(round_SSC),taxa_are_rows = TRUE)
# TAX = tax_table(tax_mat)
TAX = tax_table(as.matrix(tax_df_2))

#Sample subsetting
samples_df_sub <- subset(samples_df, samples_df$Compartment != "RZ")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "AM")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "NOD")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "Input")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Condition != "NP")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Inoculum != "AtSC")

samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

samples <- sample_data(samples_df_sub_2)

phylo_sub = phyloseq(OTU,TAX, samples)

subsetted_table <- otu_table(phylo_sub)
subsetted_table_long <- melt(subsetted_table)

Hank_the_normalizer <- function(df,group,amount){
  df_2 <- df %>% dplyr::group_by_at(group) %>% dplyr::summarise(total=sum(.data[[amount]]))
  df_3 <- df_2$total
  names(df_3) <- df_2[[group]]
  df$total <- df_3[as.character(df[[group]])]
  df$Rel <- df[[amount]] / df$total
  return(df)
}

subsetted_table_long_2 <- Hank_the_normalizer(subsetted_table_long,"Var2","value")
subsetted_table_long_2$value[subsetted_table_long_2$Rel < 0.0005] <- 0
subsetted_table_long_3 <- subsetted_table_long_2[1:3]
data_wide <- spread(subsetted_table_long_3, Var2, value)
row.names(data_wide) <- data_wide$Var1
data_wide_2 <- data_wide %>% dplyr::select (-Var1)

OTU = otu_table(as.matrix(data_wide_2),taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(tax_df_2))
samples <- sample_data(samples_df_sub_2)
phylo_sub = phyloseq(OTU,TAX, samples)

# Transform to relative abundance
phylo_sub_RA <- microbiome::transform(phylo_sub, "compositional")

# Melt phyloseq object into a dataframe
data_melt <- psmelt(phylo_sub_RA)

# Create a new OTU column where non-top5 are grouped
data_melt$OTU_grouped <- ifelse(
  data_melt$OTU %in% c("P2_G4", "LjNodule214", "P1_H10", "P2_A12", "P2_D6"), 
  data_melt$OTU, 
  ifelse(!is.na(data_melt$genus) & data_melt$genus == "Rhizobacter", "Other Rhizobacter", 
         ifelse(!is.na(data_melt$genus) & data_melt$genus == "Mesorhizobium", "Other Mesorhizobium", "Other strains"))
)

# Ensure OTU_grouped is a factor with correct legend order
otu_order <- c("P2_G4","Other Rhizobacter", "Other strains", "LjNodule214", "P1_H10", "P2_A12", "P2_D6","Other Mesorhizobium")
data_melt$OTU_grouped <- factor(data_melt$OTU_grouped, levels = otu_order)

# Define colors
otu_colors <- c(
  "P2_G4" = "#C85AC8",  # Purple
  "LjNodule214" = "#66E1D0",  # Blue
  "P1_H10" = "#00C1C8",       # Orange
  "P2_A12" = "#00AA95",       # Green
  "P2_D6" = "#00C18C",        # Red
  "Other Rhizobacter" = "#F096F0",  # Light purple
  "Other Mesorhizobium" = "#00773E",  # Greenish
  "Other strains" = "lightgray"  # Generic category
)

# Stacked barplot with grouped OTUs
bar <- ggplot(data_melt, aes(fill=OTU_grouped, y=Abundance, x=sample_id)) + 
  geom_bar(position="stack", stat="identity", colour = "darkgray", linewidth = 0.01) +  
  scale_fill_manual(values = otu_colors) +  
  ggtitle("Nodulators and Rhizobacter Root colonization profiles") + 
  theme_classic() +
  labs(y = "Relative abundance", fill = "Isolates") +  
  theme(
    plot.title = element_text(hjust = 0.5, size = 20), 
    axis.text.x = element_blank(),
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 12)
  )

# Add facet wrap by inoculum if needed
bar2 <- bar + facet_wrap(~Inoculum+Condition+Nutrient, scales = "free_x", nrow = 1) +
  theme(strip.text.x = element_text(size = 12))
bar2

pdf(paste(results.dir,"Figure_S9_dominator_taxonomic_profile.pdf", sep=""), width=15, height=4.5)
print(bar2)
dev.off()



###Figure S10 - Alpha diversity - nutrient condition subset =====
#OTU TABLE
norm_SSC=read.table(paste(working_directory,"Isolate_tables/Original/SSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
round_SSC=floor(x = norm_SSC)

#Taxonomy TABLE
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
colnames(samples_df)[6]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparation
#Set the OTU, TAX and sample data for making phyloseq object
OTU = otu_table(as.matrix(round_SSC),taxa_are_rows = TRUE)
# TAX = tax_table(tax_mat)
TAX = tax_table(as.matrix(tax_df_2))

#Sample subsetting
cond="ES"

samples_df_sub <- subset(samples_df, samples_df$Compartment != "RZ")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "AM")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "NOD")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Condition != "NP")

samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")
samples <- sample_data(samples_df_sub_2)

phylo_sub = phyloseq(OTU,TAX, samples)

subsetted_table <- otu_table(phylo_sub)
subsetted_table_long <- melt(subsetted_table)

Hank_the_normalizer <- function(df,group,amount){
  df_2 <- df %>% dplyr::group_by_at(group) %>% dplyr::summarise(total=sum(.data[[amount]]))
  df_3 <- df_2$total
  names(df_3) <- df_2[[group]]
  df$total <- df_3[as.character(df[[group]])]
  df$Rel <- df[[amount]] / df$total
  return(df)
}

subsetted_table_long_2 <- Hank_the_normalizer(subsetted_table_long,"Var2","value")
subsetted_table_long_2$value[subsetted_table_long_2$Rel < 0.0005] <- 0
subsetted_table_long_3 <- subsetted_table_long_2[1:3]
data_wide <- spread(subsetted_table_long_3, Var2, value)
row.names(data_wide) <- data_wide$Var1
data_wide_2 <- data_wide %>% dplyr::select (-Var1)

OTU = otu_table(as.matrix(data_wide_2),taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(tax_df_2))
samples_df_sub <- subset(samples_df, samples_df$Compartment != "RZ")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "AM")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "NOD")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Condition != "NP")

samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")
samples_df_sub_2$Condition[samples_df_sub_2$Condition == "At"] <- "Arabidopsis"
samples_df_sub_2$Condition[samples_df_sub_2$Condition == "Hv"] <- "Barley"
samples_df_sub_2$Condition[samples_df_sub_2$Condition == "Lj"] <- "Lotus"

samples <- sample_data(samples_df_sub_2)
phylo_sub = phyloseq(OTU,TAX, samples)

#Observed isolates
method="Observed"

obs_df=estimate_richness(physeq = phylo_sub, measures = method)
obs_df=merge(x = obs_df, y = samples_df_sub_2, by = "row.names" )

obs_df=estimate_richness(physeq = phylo_sub, measures = method)
obs_df=merge(x = obs_df, y = samples_df_sub_2, by = "row.names" )
obs_df$inoculum_experiment=paste(obs_df$Condition,obs_df$Inoculum, sep =" ")

obs_df$Inoculum <- factor(obs_df$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC"))
obs_df$Condition <- factor(obs_df$Condition, levels = c("Input","Arabidopsis","Barley", "Lotus"))
obs_df$Nutrient <- factor(obs_df$Nutrient, levels = c("high","low", "Input"))

colnames(obs_df)[colnames(obs_df) == "Condition"] <- "Host"

obs_df_2 <- obs_df[obs_df$Host != "Input",]

obs_nut_sub <- ggplot(obs_df_2, aes(x = Nutrient, y = Observed, color=Nutrient)) +
  theme_classic()+
  scale_color_manual(values = c("black", "gray70"))+
  scale_shape_manual(values = c(0,3))+
  scale_fill_manual(values = c("white","gray70"))+
  geom_jitter()+
  geom_boxplot(outlier.shape = NA) + # Hide outliers since jitter will show all points
  facet_wrap(~inoculum_experiment, scales="free_x", nrow=2) +
  theme( axis.text.x=element_blank(), 
         axis.title.x=element_blank(), 
         title=element_text(hjust=0.5, size=15), 
         axis.ticks.x=element_blank(),
         strip.background=element_rect(colour="gray50", size=0.3), # Change 'size' for thickness
         axis.text=element_text(color="gray50"),
         axis.line = element_line(color="gray50", size=0.3)) +
  labs(title = "",
       x = "",
       y = "Observed isolates") +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), label = "p.signif", vjust = 0.7, label.x = 1.5)
obs_nut_sub

#Shannon diversity isolates
method="Shannon"

shannon_df=estimate_richness(physeq = phylo_sub, measures = method)
shannon_df=merge(x = shannon_df, y = samples_df_sub_2, by = "row.names" )

shannon_df=estimate_richness(physeq = phylo_sub, measures = method)
shannon_df=merge(x = shannon_df, y = samples_df_sub_2, by = "row.names" )
shannon_df$inoculum_experiment=paste(shannon_df$Condition,shannon_df$Inoculum, sep =" ")

shannon_df$Inoculum <- factor(shannon_df$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC"))
shannon_df$Condition <- factor(shannon_df$Condition, levels = c("Input","Arabidopsis","Barley", "Lotus"))
shannon_df$Nutrient <- factor(shannon_df$Nutrient, levels = c("high","low", "Input"))

colnames(shannon_df)[colnames(shannon_df) == "Condition"] <- "Host"

shannon_df_2 <- shannon_df[shannon_df$Host != "Input",]

shan_nut_sub <- ggplot(shannon_df_2, aes(x = Nutrient, y = Shannon, color=Nutrient)) +
  theme_classic()+
  scale_color_manual(values = c("black", "gray70"))+
  scale_shape_manual(values = c(0,3))+
  scale_fill_manual(values = c("white","gray70"))+
  geom_jitter()+
  geom_boxplot(outlier.shape = NA) + # Hide outliers since jitter will show all points
  facet_wrap(~inoculum_experiment, scales="free_x", nrow=2) +
  theme( axis.text.x=element_blank(), 
         axis.title.x=element_blank(), 
         title=element_text(hjust=0.5, size=15), 
         axis.ticks.x=element_blank(),
         strip.background=element_rect(colour="gray50", size=0.3), # Change 'size' for thickness
         axis.text=element_text(color="gray50"),
         axis.line = element_line(color="gray50", size=0.3)) +
  labs(title = "",
       x = "",
       y = "Shannon isolates") +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), label = "p.signif", vjust = 0.7, label.x = 1.5)
shan_nut_sub

#KO OTU TABLE
norm_SSC_KO =read.table(paste(working_directory,"KO_tables/Original/SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
round_SSC_KO=floor(x = norm_SSC_KO)

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
colnames(samples_df)[6]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparation
#Set the OTU, TAX and sample data for making phyloseq object
OTU = otu_table(as.matrix(round_SSC_KO),taxa_are_rows = TRUE)

#Sample subsetting
cond="ES"

samples_df_sub <- subset(samples_df, samples_df$Compartment != "RZ")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "AM")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "NOD")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Condition != "NP")

samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")
samples <- sample_data(samples_df_sub_2)

round_SSC_KO_2 <- round_SSC_KO[,colnames(round_SSC_KO) %in% row.names(samples_df_sub_2)]
round_SSC_2 <- round_SSC[,colnames(round_SSC) %in% row.names(samples_df_sub_2)]

average_mult_tax_to_KO <- sum(colSums(round_SSC_KO_2)/colSums(round_SSC_2))/length(colSums(round_SSC_KO_2))

phylo_KO_sub = phyloseq(OTU,samples)

subsetted_table <- otu_table(phylo_KO_sub)
subsetted_table_long <- melt(subsetted_table)

Hank_the_normalizer <- function(df,group,amount){
  df_2 <- df %>% dplyr::group_by_at(group) %>% dplyr::summarise(total=sum(.data[[amount]]))
  df_3 <- df_2$total
  names(df_3) <- df_2[[group]]
  df$total <- df_3[as.character(df[[group]])]
  df$Rel <- df[[amount]] / df$total
  return(df)
}

subsetted_table_long_2 <- Hank_the_normalizer(subsetted_table_long,"Var2","value")
subsetted_table_long_2$value[subsetted_table_long_2$Rel < 0.0005/average_mult_tax_to_KO] <- 0
subsetted_table_long_3 <- subsetted_table_long_2[1:3]
data_wide <- spread(subsetted_table_long_3, Var2, value)
row.names(data_wide) <- data_wide$Var1
data_wide_2 <- data_wide %>% dplyr::select (-Var1)

OTU = otu_table(as.matrix(data_wide_2),taxa_are_rows = TRUE)
samples_df_sub <- subset(samples_df, samples_df$Compartment != "RZ")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "AM")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "NOD")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Condition != "NP")

samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")
samples_df_sub_2$Condition[samples_df_sub_2$Condition == "At"] <- "Arabidopsis"
samples_df_sub_2$Condition[samples_df_sub_2$Condition == "Hv"] <- "Barley"
samples_df_sub_2$Condition[samples_df_sub_2$Condition == "Lj"] <- "Lotus"

samples <- sample_data(samples_df_sub_2)
phylo_KO_sub = phyloseq(OTU, samples)

#Observed KOs
method="Observed"

obs_df_KO=estimate_richness(physeq = phylo_KO_sub, measures = method)
obs_df_KO=merge(x = obs_df_KO, y = samples_df_sub_2, by = "row.names" )

obs_df_KO=estimate_richness(physeq = phylo_KO_sub, measures = method)
obs_df_KO=merge(x = obs_df_KO, y = samples_df_sub_2, by = "row.names" )
obs_df_KO$inoculum_experiment=paste(obs_df_KO$Condition,obs_df_KO$Inoculum, sep =" ")

obs_df_KO$Inoculum <- factor(obs_df_KO$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC"))
obs_df_KO$Condition <- factor(obs_df_KO$Condition, levels = c("Input","Arabidopsis","Barley", "Lotus"))
obs_df_KO$Nutrient <- factor(obs_df_KO$Nutrient, levels = c("high","low", "Input"))

colnames(obs_df_KO)[colnames(obs_df_KO) == "Condition"] <- "Host"

obs_df_KO_2 <- obs_df_KO[obs_df_KO$Host != "Input",]

obs_nut_sub_KO <- ggplot(obs_df_KO_2, aes(x = Nutrient, y = Observed, color=Nutrient)) +
  theme_classic()+
  scale_color_manual(values = c("black", "gray70"))+
  scale_shape_manual(values = c(0,3))+
  scale_fill_manual(values = c("white","gray70"))+
  geom_jitter()+
  geom_boxplot(outlier.shape = NA) + # Hide outliers since jitter will show all points
  facet_wrap(~inoculum_experiment, scales="free_x", nrow=2) +
  theme( axis.text.x=element_blank(), 
         axis.title.x=element_blank(), 
         title=element_text(hjust=0.5, size=15), 
         axis.ticks.x=element_blank(),
         strip.background=element_rect(colour="gray50", size=0.3), # Change 'size' for thickness
         axis.text=element_text(color="gray50"),
         axis.line = element_line(color="gray50", size=0.3)) +
  labs(title = "",
       x = "",
       y = "Observed KOs") +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), label = "p.signif", vjust = 0.7, label.x = 1.5)
obs_nut_sub_KO

#Shannon diversity KOs
method="Shannon"

shannon_KO_df=estimate_richness(physeq = phylo_sub, measures = method)
shannon_KO_df=merge(x = shannon_KO_df, y = samples_df_sub_2, by = "row.names" )

shannon_KO_df=estimate_richness(physeq = phylo_sub, measures = method)
shannon_KO_df=merge(x = shannon_KO_df, y = samples_df_sub_2, by = "row.names" )
shannon_KO_df$inoculum_experiment=paste(shannon_KO_df$Condition,shannon_KO_df$Inoculum, sep =" ")

shannon_KO_df$Inoculum <- factor(shannon_KO_df$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC"))
shannon_KO_df$Condition <- factor(shannon_KO_df$Condition, levels = c("Input","Arabidopsis","Barley", "Lotus"))
shannon_KO_df$Nutrient <- factor(shannon_KO_df$Nutrient, levels = c("high","low", "Input"))

colnames(shannon_KO_df)[colnames(shannon_KO_df) == "Condition"] <- "Host"

shannon_KO_df_2 <- shannon_KO_df[shannon_KO_df$Host != "Input",]

shan_nut_sub_KO <- ggplot(shannon_KO_df_2, aes(x = Nutrient, y = Shannon, color=Nutrient)) +
  theme_classic()+
  scale_color_manual(values = c("black", "gray70"))+
  scale_shape_manual(values = c(0,3))+
  scale_fill_manual(values = c("white","gray70"))+
  geom_jitter()+
  geom_boxplot(outlier.shape = NA) + # Hide outliers since jitter will show all points
  facet_wrap(~inoculum_experiment, scales="free_x", nrow=2) +
  theme( axis.text.x=element_blank(), 
         axis.title.x=element_blank(), 
         title=element_text(hjust=0.5, size=15), 
         axis.ticks.x=element_blank(),
         strip.background=element_rect(colour="gray50", size=0.3), # Change 'size' for thickness
         axis.text=element_text(color="gray50"),
         axis.line = element_line(color="gray50", size=0.3)) +
  labs(title = "",
       x = "",
       y = "Shannon KOs") +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), label = "p.signif", vjust = 0.7, label.x = 1.5)
shan_nut_sub_KO

alpha_shan <- ggarrange(obs_nut_sub, obs_nut_sub_KO, shan_nut_sub,shan_nut_sub_KO,ncol = 1, nrow = 4, labels = c("a", "b", "c", "d"), common.legend = T)

pdf(paste(results.dir,"Figure_S10_alpha_diversity_nutrients.pdf", sep=""), width=8, height=18)
print(alpha_shan)
dev.off()

###Figure 2fgh - Beta Diversity plots =====
#otu table
norm_SSC =read.table(paste(working_directory,"Isolate_tables/Original/SSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
#Taxonomy table
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")
#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
colnames(samples_df_2)[5]="Nutrient"
samples_df_2$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df_2$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

samples_df_2$Condition <- gsub("At", "Arabidopsis", samples_df_2$Condition)
samples_df_2$Condition <- gsub("Hv", "Barley", samples_df_2$Condition)
samples_df_2$Condition <- gsub("Lj", "Lotus", samples_df_2$Condition)
samples_df_2$Nutrient <- gsub("high", "High", samples_df_2$Nutrient)
samples_df_2$Nutrient <- gsub("low", "Low", samples_df_2$Nutrient)

sapply(tax_df, function(x) length(unique(x)))
#  Class  number is not very different from Phylum number, skip class rank that gives weird results because nclass=6 and nphylum=5 

#Set the OTU, TAX and sample data for making phyloseq object
OTU = otu_table(as.matrix(norm_SSC),taxa_are_rows = TRUE)
#TAX = tax_table(tax_mat)
TAX = tax_table(as.matrix(tax_df_2))

#Sample subsetting
cond="ES"
samples_df_sub <- subset(samples_df_2, samples_df_2$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")
samples_sub = sample_data(samples_df_sub_2)

phylo_sub = phyloseq(OTU,TAX, samples_sub)

phylo_sub_RA=microbiome::transform(x = phylo_sub, transform = "compositional" )

#Agglomerate to phylum-level and rename
ps_phylum <- phyloseq::tax_glom(phylo_sub_RA, "Phylum")
ps_class <- phyloseq::tax_glom(phylo_sub_RA, "Class")
ps_order <- phyloseq::tax_glom(phylo_sub_RA, "Order")
ps_family <- phyloseq::tax_glom(phylo_sub_RA, "Family")
ps_genus <- phyloseq::tax_glom(phylo_sub_RA, "Genus")

phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]
phyloseq::taxa_names(ps_class) <- phyloseq::tax_table(ps_class)[, "Class"]
phyloseq::taxa_names(ps_order) <- phyloseq::tax_table(ps_order)[, "Order"]
phyloseq::taxa_names(ps_family) <- phyloseq::tax_table(ps_family)[, "Family"]
phyloseq::taxa_names(ps_genus) <- phyloseq::tax_table(ps_genus)[, "Genus"]

#Bray Curtis distance matrix
beta_phylum <- as.matrix(vegdist(t(ps_phylum@otu_table@.Data), method = "bray", diag = T))
beta_class <- as.matrix(vegdist(t(ps_class@otu_table@.Data), method = "bray", diag = T))
beta_order <- as.matrix(vegdist(t(ps_order@otu_table@.Data), method = "bray", diag = T))
beta_family <- as.matrix(vegdist(t(ps_family@otu_table@.Data), method = "bray", diag = T))
beta_genus <- as.matrix(vegdist(t(ps_genus@otu_table@.Data), method = "bray", diag = T))
beta_isolate <- as.matrix(vegdist(t(phylo_sub_RA@otu_table@.Data), method = "bray", diag = T))

mean_value_phylum= mean(beta_phylum)
mean_value_class= mean(beta_class)
mean_value_order= mean(beta_order)
mean_value_family= mean(beta_family)
mean_value_genus= mean(beta_genus)
mean_value_isolate= mean(beta_isolate)

Taxo_order=c("Isolate","Genus", "Family", "Order", "Class")
beta_distance_list=list(beta_isolate, beta_genus, beta_family, beta_order, beta_phylum)

Bray_curtis_df=beta_genus

#Make PCoA plot for Bray Curtis Distance matrix
pcoa_tax = cmdscale(Bray_curtis_df, k=3, eig=T)
points_beta = as.data.frame(pcoa_tax$points)
colnames(points_beta) = c("x", "y", "z") 
eig = pcoa_tax$eig

points_beta = merge(points_beta,samples_df_sub_2, by = "row.names")
rownames(points_beta) <- points_beta$Row.names
points_beta <- points_beta %>% dplyr::select (-Row.names)

points_beta$Condition <- factor(points_beta$Condition, levels = c("Arabidopsis","Barley", "Lotus"))
points_beta$Inoculum <- factor(points_beta$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC","NS"))
points_beta$Nutrient <- factor(points_beta$Nutrient, levels = c("Low", "High"))
points_beta$Experiment  <- factor(points_beta$Experiment, levels = c("R1", "R2"))

#  Run adonis PERMANOVA test
metadata=points_beta[,-c(1,2,3)]
set.seed(1)
SSC_bray_adonis <- adonis2(beta_distance_list[[2]] ~ Inoculum*Condition*Nutrient*Experiment, data=metadata, method="bray", permutations=999)

#Bray Curtis - SynCom
if (SSC_bray_adonis$`Pr(>F)`[1] <= 0.001) {
  sig_syncom <- "***"
} else if (SSC_bray_adonis$`Pr(>F)`[1] <= 0.01) {
  sig_syncom <- "**"
} else if (SSC_bray_adonis$`Pr(>F)`[1] <= 0.05) {
  sig_syncom <- "*"
} else {
  sig_syncom <- "ns"
}

#Bray Curtis - Plant
if (SSC_bray_adonis$`Pr(>F)`[2] <= 0.001) {
  sig_plant <- "***"
} else if (SSC_bray_adonis$`Pr(>F)`[2] <= 0.01) {
  sig_plant <- "**"
} else if (SSC_bray_adonis$`Pr(>F)`[2] <= 0.05) {
  sig_plant <- "*"
} else {
  sig_plant <- "ns"
}

#Bray Curtis - Nutrient
if (SSC_bray_adonis$`Pr(>F)`[3] <= 0.001) {
  sig_nut <- "***"
} else if (SSC_bray_adonis$`Pr(>F)`[3] <= 0.01) {
  sig_nut <- "**"
} else if (SSC_bray_adonis$`Pr(>F)`[3] <= 0.05) {
  sig_nut <- "*"
} else {
  sig_nut <- "ns"
}

#Bray Curtis - Experiment
if (SSC_bray_adonis$`Pr(>F)`[4] <= 0.001) {
  sig_exp <- "***"
} else if (SSC_bray_adonis$`Pr(>F)`[4] <= 0.01) {
  sig_exp <- "**"
} else if (SSC_bray_adonis$`Pr(>F)`[4] <= 0.05) {
  sig_exp <- "*"
} else {
  sig_exp <- "ns"
}

#Bray Curtis - SynCom:Plant Interaction
if (SSC_bray_adonis$`Pr(>F)`[5] <= 0.001) {
  sig_inter <- "***"
} else if (SSC_bray_adonis$`Pr(>F)`[5] <= 0.01) {
  sig_inter <- "**"
} else if (SSC_bray_adonis$`Pr(>F)`[5] <= 0.05) {
  sig_inter <- "*"
} else {
  sig_inter <- "ns"
}

colnames(points_beta)[colnames(points_beta) == "Condition"] <- "Host"

pbray_tax_1 <- ggplot(points_beta, aes(x=x, y=y, fill=Inoculum, color=Nutrient, shape=Host, stroke=0.3))+
  scale_shape_manual(values = c(23,24,21))+
  scale_color_manual(values = c("gray70","black"))+
  scale_fill_manual(values = c("#A3A500","#00B0F6","#00BF7D","#F8766D","White"))+
  theme(panel.background=element_blank(),panel.grid=element_blank(),axis.line.x=element_line(size=.5, color="black"),axis.line.y=element_line(size=.5, color="black"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black", size=7),legend.position="right",legend.background=element_blank(),legend.key=element_blank(),legend.text= element_text(size=10),text=element_text(family="sans", size=10))+
  geom_point(alpha=1, size=4)+
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""))+
  theme(axis.title = element_text(size = 14, hjust = 1), axis.text = element_text(size = 12, hjust = 1))+
  guides(fill = guide_legend(override.aes=list(shape=21, size=4)))+
  guides(color = guide_legend(override.aes=list(shape=21, size = 4))) +
  guides(shape = guide_legend(override.aes=list(size = 4))) +
  stat_ellipse(aes(group = Inoculum ), size=0.1)+
  ggtitle(paste0("Taxonomic (",Taxo_order[2], ") composition"))
pbray_tax_1

text <- c(paste0("PCoA Bray Curtis - ",Taxo_order[2], " composition"),
          paste("Inoculum - R2 = ",round(SSC_bray_adonis$R2[1],3),sig_syncom),
          paste("Host - R2 = ",round(SSC_bray_adonis$R2[2],3),sig_plant),
          paste("Nutrient - R2 = ",round(SSC_bray_adonis$R2[3],3),sig_nut),
          paste("Experiment - R2 = ",round(SSC_bray_adonis$R2[4],3),sig_exp),
          paste("Inoculum:Host - R2 = ",round(SSC_bray_adonis$R2[5],3),sig_inter))

tl <- textGrob(paste(strwrap(text, 40), collapse="\n"), hjust=0, x=0)
pbray_tax_1_2 <- ggarrange(pbray_tax_1,tl, ncol = 2, nrow =1, widths = c(2,1))

# KO_otu table
KO_SSC =read.table(paste(working_directory,"KO_tables/Original/SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Taxonomy table

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object
OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)

#Sample subsetting
cond="ES"

samples_df_sub <- subset(samples_df_2, samples_df_2$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")
samples_df_sub_2$Condition[samples_df_sub_2$Condition == "At"] <- "Arabidopsis"
samples_df_sub_2$Condition[samples_df_sub_2$Condition == "Hv"] <- "Barley"
samples_df_sub_2$Condition[samples_df_sub_2$Condition == "Lj"] <- "Lotus"
samples_df_sub_2$Nutrient[samples_df_sub_2$Nutrient == "high"] <- "High"
samples_df_sub_2$Nutrient[samples_df_sub_2$Nutrient == "low"] <- "Low"
samples_sub = sample_data(samples_df_sub_2)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)

phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )

beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

Bray_curtis_KO_df=beta_isolate_KO

#Make PCoA plot for Bray Curtis Distance matrix
pcoa_KO = cmdscale(Bray_curtis_KO_df, k=3, eig=T)
points_beta_KO = as.data.frame(pcoa_KO$points)
colnames(points_beta_KO) = c("x", "y", "z") 
eig_KO = pcoa_KO$eig

points_beta_KO = merge(points_beta_KO,samples_df_sub_2, by = "row.names")
rownames(points_beta_KO) <- points_beta_KO$Row.names
points_beta_KO <- points_beta_KO %>% dplyr::select (-Row.names)

points_beta_KO$Condition <- factor(points_beta_KO$Condition, levels = c("Arabidopsis","Barley", "Lotus"))
points_beta_KO$Inoculum <- factor(points_beta_KO$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC","NS"))
points_beta_KO$Nutrient <- factor(points_beta_KO$Nutrient, levels = c("Low", "High"))
points_beta_KO$Experiment  <- factor(points_beta_KO$Experiment, levels = c("R1", "R2"))

# PERMANOVA Adonis test
set.seed(1)
metadata=points_beta_KO[,-c(1,2,3)]
SSC_bray_KO_adonis <- adonis2(beta_isolate_KO ~ Inoculum*Condition*Nutrient*Experiment, data=metadata, method="bray", permutations=999)

#Bray Curtis - SynCom
if (SSC_bray_KO_adonis$`Pr(>F)`[1] <= 0.001) {
  sig_syncom <- "***"
} else if (SSC_bray_KO_adonis$`Pr(>F)`[1] <= 0.01) {
  sig_syncom <- "**"
} else if (SSC_bray_KO_adonis$`Pr(>F)`[1] <= 0.05) {
  sig_syncom <- "*"
} else {
  sig_syncom <- "ns"
}

#Bray Curtis - Plant
if (SSC_bray_KO_adonis$`Pr(>F)`[2] <= 0.001) {
  sig_plant <- "***"
} else if (SSC_bray_KO_adonis$`Pr(>F)`[2] <= 0.01) {
  sig_plant <- "**"
} else if (SSC_bray_KO_adonis$`Pr(>F)`[2] <= 0.05) {
  sig_plant <- "*"
} else {
  sig_plant <- "ns"
}

#Bray Curtis - Nutrient
if (SSC_bray_KO_adonis$`Pr(>F)`[3] <= 0.001) {
  sig_nut <- "***"
} else if (SSC_bray_KO_adonis$`Pr(>F)`[3] <= 0.01) {
  sig_nut <- "**"
} else if (SSC_bray_KO_adonis$`Pr(>F)`[3] <= 0.05) {
  sig_nut <- "*"
} else {
  sig_nut <- "ns"
}

#Bray Curtis - Experiment
if (SSC_bray_KO_adonis$`Pr(>F)`[4] <= 0.001) {
  sig_exp <- "***"
} else if (SSC_bray_KO_adonis$`Pr(>F)`[4] <= 0.01) {
  sig_exp <- "**"
} else if (SSC_bray_KO_adonis$`Pr(>F)`[4] <= 0.05) {
  sig_exp <- "*"
} else {
  sig_exp <- "ns"
}

#Bray Curtis - SynCom:Plant Interaction
if (SSC_bray_KO_adonis$`Pr(>F)`[5] <= 0.001) {
  sig_inter <- "***"
} else if (SSC_bray_KO_adonis$`Pr(>F)`[5] <= 0.01) {
  sig_inter <- "**"
} else if (SSC_bray_KO_adonis$`Pr(>F)`[5] <= 0.05) {
  sig_inter <- "*"
} else {
  sig_inter <- "ns"
}

colnames(points_beta_KO)[colnames(points_beta_KO) == "Condition"] <- "Host"

pbray_KO_1 <- ggplot(points_beta_KO, aes(x=x, y=y, fill=Inoculum, color=Nutrient, shape=Host,size=12,stroke=0.3))+
  scale_shape_manual(values = c(23,24,21))+
  scale_color_manual(values = c("gray70","black"))+
  scale_fill_manual(values = c("#A3A500","#00B0F6","#00BF7D","#F8766D","White"))+
  theme(panel.background=element_blank(),panel.grid=element_blank(),axis.line.x=element_line(size=.5, color="black"),axis.line.y=element_line(size=.5, color="black"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black", size=7),legend.position="none",legend.background=element_blank(),legend.key=element_blank(),legend.text= element_text(size=10),text=element_text(family="sans", size=10))+
  geom_point(alpha=1, size=4,)+
  labs(x=paste("PCoA 1 (", format(100 * eig_KO[1] / sum(eig_KO), digits=4), "%)", sep=""),y=paste("PCoA 2 (", format(100 * eig_KO[2] / sum(eig_KO), digits=4), "%)", sep=""))+
  theme(axis.title = element_text(size = 14, hjust = 1), axis.text = element_text(size = 12, hjust = 1))+
  guides(fill = guide_legend(override.aes=list(shape=21, size=4)))+
  guides(color = guide_legend(override.aes=list(shape=21, size = 4))) +
  guides(shape = guide_legend(override.aes=list(size = 4))) +
  stat_ellipse(aes(group = Inoculum ), size=0.1)+
  ggtitle("Functional (KO) composition")
pbray_KO_1

#Plot text next to it
text_KO <- c(paste0("PCoA Bray Curtis - KO composition"), 
             paste("Inoculum - R2 = ",round(SSC_bray_KO_adonis$R2[1],3),sig_syncom),
             paste("Host - R2 = ",round(SSC_bray_KO_adonis$R2[2],3),sig_plant),
             paste("Nutrient - R2 = ",round(SSC_bray_KO_adonis$R2[3],3),sig_nut),
             paste("Experiment - R2 = ",round(SSC_bray_KO_adonis$R2[4],3),sig_exp),
             paste("Inoculum:Host - R2 = ",round(SSC_bray_KO_adonis$R2[5],3),sig_inter))

tl <- textGrob(paste(strwrap(text_KO, 40), collapse="\n"), hjust=0, x=0)
plot_KO <- ggarrange(pbray_KO_1,tl, ncol = 2, nrow =1, widths = c(2,1))

genus_plot <- pbray_tax_1
genus_plot<- genus_plot+theme(legend.position = "none")

# Define legend text size and key size as variables
text_size <- 14
key_size <- 1.2

# Create a theme object with your legend modifications
legend_theme <- theme(legend.text = element_text(size = text_size),  # Apply text size
                      legend.key.size = unit(key_size, "cm"), # Apply key size
                      strip.background=element_rect(colour="gray50", size=0.3), # Change 'size' for thickness
                      axis.text=element_text(color="gray50"),
                      axis.line = element_line(color="gray50", size=0.3),
                      axis.line.x = element_line(color="gray50", size=0.3), 
                      axis.line.y = element_line(color="gray50", size=0.3),
                      axis.ticks.x = element_line(color="gray50", size=0.3),
                      axis.ticks.y =element_line(color="gray50", size=0.3)
)

# Apply the legend_theme to each ggplot object
genus_plot <- genus_plot + legend_theme
pbray_KO_1 <- pbray_KO_1 + legend_theme

#Adonis test on taxonomy tables
Adonis_list=list()
Plot_list=list()
Plot_list_2=list()

for (i in 1:length(beta_distance_list)) {
  
  Bray_curtis_df=beta_distance_list[[i]]
  
  #Make PCoA plot for Bray Curtis Distance matrix
  pcoa = cmdscale(Bray_curtis_df, k=3, eig=T)
  points = as.data.frame(pcoa$points)
  colnames(points) = c("x", "y", "z") 
  eig = pcoa$eig
  
  points = merge(points,samples_df_sub_2, by = "row.names")
  rownames(points) <- points$Row.names
  points <- points %>% dplyr::select (-Row.names)
  
  points$Condition <- factor(points$Condition, levels = c("Arabidopsis","Barley", "Lotus"))
  points$Inoculum <- factor(points$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC","NS"))
  points$Nutrient <- factor(points$Nutrient, levels = c("Low", "High"))
  points$Experiment  <- factor(points$Experiment, levels = c("R1", "R2"))
  
  metadata=points[,-c(1,2,3)]
  
  #  Run adonis PERMANOVA test
  set.seed(1)
  SSC_bray_adonis <- adonis2(beta_distance_list[[i]] ~ Inoculum*Condition*Nutrient*Experiment, data=metadata, method="bray", permutations=999)
  
  Adonis_list[[i]]=SSC_bray_adonis
  
}

# Add values of the adonis tests in dataframes to make a table:
Taxo_R2=as.data.frame(Adonis_list[[2]])
Func_R2=as.data.frame(SSC_bray_KO_adonis)

R2_table=merge.data.frame(Taxo_R2,Func_R2, by = 0)
R2_table=R2_table[,-c(2,3,5,6,7,8,10,11)]

# Subset specific rows
R2_table <- R2_table[c(1,4,5,6,11), ]  # Select rows 2 to 4 and rows 11 to 13
R2_table <- R2_table[order(-R2_table$R2.x), ]

# Round numeric values to 2 digits after the decimal point
R2_table$R2.x <- round(R2_table$R2.x, digits = 2)
R2_table$R2.y <- round(R2_table$R2.y, digits = 2)

# Rename columns
colnames(R2_table)=c("Taxonomy R²","Parameter","Function R²")

# Replace "Condition" with "Plant"
R2_table$Parameter <- gsub("Condition", "Host", R2_table$Parameter)

# my_table <- tableGrob(R2_table)
ggtexttable(R2_table, rows = NULL, theme = ttheme("classic"))

beta_figures <- ggarrange(genus_plot,"", pbray_KO_1+scale_y_continuous(position="right"), nrow =1, ncol =3, common.legend = F, widths = c(1,0.3,1))
beta_figures

pdf(paste(results.dir,"Figure_2fg_Beta_div_plots.pdf", sep=""), width=width_in, height=height_in)
print(beta_figures)
dev.off()

pdf(paste(results.dir,"Figure_2h_table_R2.pdf", sep=""), width=width_in, height=height_in)
ggtexttable(R2_table, rows = NULL, theme = ttheme("classic"))
dev.off()

###Figure S11 & 2i - Beta Diversity - distance to centroid =====
#otu table
norm_SSC =read.table(paste(working_directory,"Isolate_tables/Original/SSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
#Taxonomy table
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")
#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
colnames(samples_df_2)[5]="Nutrient"
samples_df_2$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df_2$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")
samples_df_2$Plant_inoc_compartment_nutrient=paste(samples_df$Condition, samples_df$Inoculum, samples_df$Compartment, samples_df$Nutrient, sep ="_")
samples_df_2$Plant_inoc_compartment=paste(samples_df$Condition, samples_df$Inoculum, samples_df$Compartment, sep ="_")

sapply(tax_df, function(x) length(unique(x)))
#  Class  number is not very different from Phylum number, skip class rank that gives weird results because nclass=6 and nphylum=5 

#Set the OTU, TAX and sample data for making phyloseq object
OTU = otu_table(as.matrix(norm_SSC),taxa_are_rows = TRUE)
#TAX = tax_table(tax_mat)
TAX = tax_table(as.matrix(tax_df_2))

#Sample subsetting

cond="ES"
samples_df_sub <- subset(samples_df_2, samples_df_2$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

samples_sub = sample_data(samples_df_sub_2)

phylo_sub = phyloseq(OTU,TAX, samples_sub)

phylo_sub_RA=microbiome::transform(x = phylo_sub, transform = "compositional" )

#Agglomerate to phylum-level and rename
ps_genus <- phyloseq::tax_glom(phylo_sub_RA, "Genus")

phyloseq::taxa_names(ps_genus) <- phyloseq::tax_table(ps_genus)[, "Genus"]

#Bray Curtis distance matrix
beta_genus <- as.matrix(vegdist(t(ps_genus@otu_table@.Data), method = "bray", diag = T))

mean_value_genus= mean(beta_genus)

Bray_curtis_df=beta_genus

#Make PCoA plot for Bray Curtis Distance matrix
pcoa_tax = cmdscale(Bray_curtis_df, k=3, eig=T)
points_beta = as.data.frame(pcoa_tax$points)
colnames(points_beta) = c("x", "y", "z") 
eig = pcoa_tax$eig

points_beta = merge(points_beta,samples_df_sub_2, by = "row.names")
rownames(points_beta) <- points_beta$Row.names
points_beta <- points_beta %>% dplyr::select (-Row.names)

points_beta$Condition <- factor(points_beta$Condition, levels = c("At","Hv", "Lj"))
points_beta$Inoculum <- factor(points_beta$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC","NS"))
points_beta$Nutrient <- factor(points_beta$Nutrient, levels = c("low", "high"))
points_beta$Experiment  <- factor(points_beta$Experiment, levels = c("R1", "R2"))

# KO_otu table
KO_SSC =read.table(paste(working_directory,"KO_tables/Original/SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object
OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)

#Sample subsetting
cond="ES"

samples_df_sub <- subset(samples_df_2, samples_df_2$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

samples_sub = sample_data(samples_df_sub_2)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)

phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )

beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

Bray_curtis_KO_df=beta_isolate_KO

#Make PCoA plot for Bray Curtis Distance matrix
pcoa_KO = cmdscale(Bray_curtis_KO_df, k=3, eig=T)
points_beta_KO = as.data.frame(pcoa_KO$points)
colnames(points_beta_KO) = c("x", "y", "z") 
eig_KO = pcoa_KO$eig

points_beta_KO = merge(points_beta_KO,samples_df_sub_2, by = "row.names")
rownames(points_beta_KO) <- points_beta_KO$Row.names
points_beta_KO <- points_beta_KO %>% dplyr::select (-Row.names)

points_beta_KO$Condition <- factor(points_beta_KO$Condition, levels = c("At","Hv", "Lj"))
points_beta_KO$Inoculum <- factor(points_beta_KO$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC","NS"))
points_beta_KO$Nutrient <- factor(points_beta_KO$Nutrient, levels = c("low", "high"))
points_beta_KO$Experiment  <- factor(points_beta_KO$Experiment, levels = c("R1", "R2"))

# Subset tables if needed, by experiment
points_beta_tax_sub=subset.data.frame(x = points_beta, subset = points_beta$Experiment=="R2")
points_beta_KO_sub=subset.data.frame(x = points_beta_KO, subset = points_beta$Experiment=="R2")

# Define parameter for centroid calculation
param="Plant_inoc_compartment_nutrient"

# Calculate centroids for each group
centroids_tax <- points_beta_tax_sub %>%
  group_by(!!sym(param)) %>%
  dplyr::summarize(
    centroid_x = mean(x, na.rm = TRUE),
    centroid_y = mean(y, na.rm = TRUE),
    centroid_z = mean(z, na.rm = TRUE)
  )

# Join centroids back to the original data
data_with_centroids_tax <- left_join(points_beta_tax_sub, centroids_tax , by = param)

# Calculate distance to centroid for each point
data_with_centroids_tax  <- data_with_centroids_tax  %>% 
  rowwise() %>%
  mutate(distance_to_centroid = sqrt((x - centroid_x)^2 + 
                                       (y - centroid_y)^2 + 
                                       (z - centroid_z)^2))

custom_labels_plant <- c(At = "A.thaliana", Hv = "H. vulgare", Lj = "L. japonicus")

# Calculate centroids for each group
centroids_KO <- points_beta_KO_sub %>%
  group_by(!!sym(param)) %>%
  dplyr::summarize(
    centroid_x = mean(x, na.rm = TRUE),
    centroid_y = mean(y, na.rm = TRUE),
    centroid_z = mean(z, na.rm = TRUE)
  )

# Join centroids back to the original data
data_with_centroids_KO <- left_join(points_beta_KO_sub, centroids_KO , by = param)

# Calculate distance to centroid for each point
data_with_centroids_KO  <- data_with_centroids_KO  %>% 
  rowwise() %>%
  mutate(distance_to_centroid = sqrt((x - centroid_x)^2 + 
                                       (y - centroid_y)^2 + 
                                       (z - centroid_z)^2))

custom_labels_plant <- c(At = "A.thaliana", Hv = "H. vulgare", Lj = "L. japonicus")

#  Combine both datasets to compare functional vs taxonomic distance from centroids
data_with_centroids_KO$Type  <- "Functional" 
data_with_centroids_tax$Type  <- "Taxonomic" 

centroid_merge=rbind(data_with_centroids_KO,data_with_centroids_tax)

dist_cen <- ggplot(centroid_merge, aes(x = param, y = distance_to_centroid, color=Type)) +
  theme_classic()+
  scale_color_manual(values = c("gray70","black"))+
  scale_shape_manual(values = c(0,3))+
  scale_fill_manual(values = c("white","gray70"))+
  geom_boxplot(outlier.shape = NA) + # Hide outliers since jitter will show all points
  facet_wrap(as.formula(paste(".~", param)), scales = "free_x", ncol = nlevels(factor(centroid_merge[[param]]))/3) +
  geom_jitter()+
  theme( axis.text.x=element_blank(), 
         axis.title.x=element_blank(), 
         title=element_text(hjust=0.5, size=15), 
         axis.ticks.x=element_blank(),
         strip.background=element_rect(colour="gray50", size=0.3), # Change 'size' for thickness
         axis.text=element_text(color="gray50"),
         axis.line = element_line(color="gray50", size=0.3)) +
  labs(title = "Distance to Centroid Functional vs Taxonomic",
       x = param,
       y = "Distance to Centroid") +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), label = "p.signif", vjust = 0.7)
dist_cen

pdf(paste(results.dir,"Figure_S11_distance_to_centroid.pdf", sep=""), width=13, height=6)
print(dist_cen)
dev.off()

#Figure 2i - Distance to centroid stats 
centroid_merge$Type <- factor(centroid_merge$Type, levels = c("Functional","Taxonomic"))

centroid_merge$Type <- factor(centroid_merge$Type, levels = c("Taxonomic","Functional"))

centro_stats_hor <- ggpaired(centroid_merge, x = "Type", y = "distance_to_centroid",
                             color = "Type", line.color = "transparent", line.size = 0.1, linetype = "solid")+
  labs(x = "", y = "Distance to centroid") +
  scale_color_manual(values = c("black","gray70"))+
  stat_compare_means(paired = TRUE, size =5, hjust=-0.5) +
  theme(legend.position = "none") +
  theme(axis.title = element_text(size = 15, hjust = 0.5), axis.text.x = element_text(size = 15, angle = 45, vjust=1, hjust=1), title = element_text(size =15)) +
  ggtitle("")
centro_stats_hor

pdf(paste(results.dir,"Figure_2i_stats_distance_to_centroid_hor.pdf", sep=""), width=3, height=4)
print(centro_stats_hor)
dev.off()


###Table S2 - abundance categories R2 simulation =====
norm_SSC=read.table(paste(working_directory,"Isolate_tables/Original/SSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
round_SSC=floor(x = norm_SSC)

#Taxonomy TABLE
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
colnames(samples_df)[6]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparation
#Set the OTU, TAX and sample data for making phyloseq object
OTU = otu_table(as.matrix(round_SSC),taxa_are_rows = TRUE)
# TAX = tax_table(tax_mat)
TAX = tax_table(as.matrix(tax_df_2))

#Sample subsetting
cond="ES"

samples_df_sub <- subset(samples_df, samples_df$Compartment != "RZ")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "AM")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Compartment != "NOD")
samples_df_sub <- subset(samples_df_sub, samples_df_sub$Condition != "NP")

samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")
samples <- sample_data(samples_df_sub_2)

phylo_sub = phyloseq(OTU,TAX, samples)

subsetted_table <- otu_table(phylo_sub)
subsetted_table_long <- melt(subsetted_table)

Hank_the_normalizer <- function(df,group,amount){
  df_2 <- df %>% dplyr::group_by_at(group) %>% dplyr::summarise(total=sum(.data[[amount]]))
  df_3 <- df_2$total
  names(df_3) <- df_2[[group]]
  df$total <- df_3[as.character(df[[group]])]
  df$Rel <- df[[amount]] / df$total
  return(df)
}

subsetted_table_long_2 <- Hank_the_normalizer(subsetted_table_long,"Var2","value")
subsetted_table_long_2$value[subsetted_table_long_2$Rel < 0.0005] <- 0
subsetted_table_long_3 <- subsetted_table_long_2[1:3]
data_wide <- spread(subsetted_table_long_3, Var2, value)
row.names(data_wide) <- data_wide$Var1
data_wide_2 <- data_wide %>% dplyr::select (-Var1)

data_wide_3 <- colnames(data_wide_2)[grep("Input",colnames(data_wide_2))]
data_wide_4 <- data_wide_2[,!colnames(data_wide_2) %in% data_wide_3]

#Count no of isolates
data_wide_4[data_wide_4 > 0] <- 1
Av_no_of_isolates <- round(sum(colSums(data_wide_4))/length(colnames(data_wide_4)))
Av_no_of_isolates

#Getting average per group
new <- data.frame(colSums(data_wide_4))
colnames(new) <- "number"
new$Plant <- samples_df$Condition[match(row.names(new), row.names(samples_df))]
new$SynCom <- samples_df$Inoculum[match(row.names(new), row.names(samples_df))]
new$Sample_group <- paste(new$Plant, new$SynCom, sep = "_")

hop_2 <- data.frame()

for (group in unique(new$Plant)){
  new_sub <- new[new$Plant == paste(group),]
  value <- sum(new_sub$number)/length(new_sub$number)
  hop <- t(data.frame(c(paste(group), value)))
  hop_2 <- rbind(hop_2, hop)
}

#Distribution
data_wide_5 <- as.numeric(unlist(data_wide_2[,!colnames(data_wide_2) %in% data_wide_3]))
data_wide_6 <- data_wide_5[data_wide_5 != 0]
value <- length(data_wide_6[data_wide_6 > 500])
data_wide_7 <- data_wide_6[data_wide_6 < 500]

hist_plot <- hist(data_wide_7, breaks = 20)
new_list <- as.data.frame(hist_plot$breaks)
new_list_2 <- as.data.frame(new_list[-1,])
new_list_2$counts <- hist_plot$counts

new_list_2 <- rbind(new_list_2, c(1000,value))
val <- sum(new_list_2$counts)
new_list_2$Rel <- new_list_2$counts/val
new_val <- 100/Av_no_of_isolates
new_list_2$new_counts <- (new_list_2$Rel/new_val)*100

new_list_3 <- as.data.frame(new_list_2[,1])
new_list_3$counts <- round(new_list_2$new_counts)
colnames(new_list_3) <- c("Category", "Counts")

write.table(new_list_3, paste(working_directory, "Table_S2_abundance_categories.txt", sep = ""), sep = "\t", quote = F)

###Figure S12 and Table S3 - R2 Simulations =====
R2_values_2 <- read.table(paste(working_directory,"R2_values_genus.txt", sep = ""), sep = "\t", header = T)

#Add our own R2 value to the simulated data
#otu table
norm_SSC =read.table(paste(working_directory,"Isolate_tables/Original/SSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
#Taxonomy table
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")
#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
colnames(samples_df_2)[5]="Nutrient"
samples_df_2$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df_2$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

sapply(tax_df, function(x) length(unique(x)))

#Set the OTU, TAX and sample data for making phyloseq object
OTU = otu_table(as.matrix(norm_SSC),taxa_are_rows = TRUE)
#TAX = tax_table(tax_mat)
TAX = tax_table(as.matrix(tax_df_2))

#Sample subsetting
cond="ES"
samples_df_sub <- subset(samples_df_2, samples_df_2$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

samples_sub = sample_data(samples_df_sub_2)

phylo_sub = phyloseq(OTU,TAX, samples_sub)

phylo_sub_RA=microbiome::transform(x = phylo_sub, transform = "compositional" )

#Agglomerate to phylum-level and rename
ps_genus <- phyloseq::tax_glom(phylo_sub_RA, "Genus")
phyloseq::taxa_names(ps_genus) <- phyloseq::tax_table(ps_genus)[, "Genus"]

#Bray Curtis distance matrix
beta_genus <- as.matrix(vegdist(t(ps_genus@otu_table@.Data), method = "bray", diag = T))

#Make PCoA plot for Bray Curtis Distance matrix
pcoa = cmdscale(beta_genus, k=3, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z") 
eig = pcoa$eig

points = merge(points,samples_df_sub_2, by = "row.names")
rownames(points) <- points$Row.names
points <- points %>% dplyr::select (-Row.names)
metadata=points[,-c(1,2,3)]

set.seed(1)
SSC_bray_adonis <- adonis2(beta_genus ~ Inoculum, data=metadata, method="bray", permutations=999)
SSC_Tax_SC_value <- SSC_bray_adonis$R2[1]

SSC_bray_adonis <- adonis2(beta_genus ~ Condition, data=metadata, method="bray", permutations=999)
SSC_Tax_PL_value <- SSC_bray_adonis$R2[1]

#KO table
KO_SSC =read.table(paste(working_directory,"KO_tables/Original/SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object
OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)

#Sample subsetting
cond="ES"

samples_df_sub <- subset(samples_df_2, samples_df_2$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

samples_sub = sample_data(samples_df_sub_2)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)

phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )

beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

#Make PCoA plot for Bray Curtis Distance matrix
pcoa = cmdscale(beta_isolate_KO, k=3, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z") 
eig = pcoa$eig

points = merge(points,samples_df_sub_2, by = "row.names")
rownames(points) <- points$Row.names
points <- points %>% dplyr::select (-Row.names)

metadata=points[,-c(1,2,3)]

set.seed(1)
SSC_bray_KO_adonis <- adonis2(beta_isolate_KO ~ Inoculum, data=metadata, method="bray", permutations=999)
SSC_KO_SC_value <- SSC_bray_KO_adonis$R2[1]

SSC_bray_KO_adonis <- adonis2(beta_isolate_KO ~ Condition, data=metadata, method="bray", permutations=999)
SSC_KO_PL_value <- SSC_bray_KO_adonis$R2[1]

SSC_values_SC <- c("Inoculum", SSC_Tax_SC_value,SSC_KO_SC_value)
SSC_values_PL <- c("Host", SSC_Tax_PL_value,SSC_KO_PL_value)

#otu table - no nodulators
norm_SSC =read.table(paste(working_directory,"Isolate_tables/No_nodulators/SSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
#Taxonomy table
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")
#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
colnames(samples_df_2)[5]="Nutrient"
samples_df_2$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df_2$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

sapply(tax_df, function(x) length(unique(x)))

#Set the OTU, TAX and sample data for making phyloseq object
OTU = otu_table(as.matrix(norm_SSC),taxa_are_rows = TRUE)
#TAX = tax_table(tax_mat)
TAX = tax_table(as.matrix(tax_df_2))

#Sample subsetting
cond="ES"
samples_df_sub <- subset(samples_df_2, samples_df_2$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

samples_sub = sample_data(samples_df_sub_2)

phylo_sub = phyloseq(OTU,TAX, samples_sub)

phylo_sub_RA=microbiome::transform(x = phylo_sub, transform = "compositional" )

#Agglomerate to phylum-level and rename
ps_genus <- phyloseq::tax_glom(phylo_sub_RA, "Genus")
phyloseq::taxa_names(ps_genus) <- phyloseq::tax_table(ps_genus)[, "Genus"]

#Bray Curtis distance matrix
beta_genus <- as.matrix(vegdist(t(ps_genus@otu_table@.Data), method = "bray", diag = T))

#Make PCoA plot for Bray Curtis Distance matrix
pcoa = cmdscale(beta_genus, k=3, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z") 
eig = pcoa$eig

points = merge(points,samples_df_sub_2, by = "row.names")
rownames(points) <- points$Row.names
points <- points %>% dplyr::select (-Row.names)
metadata=points[,-c(1,2,3)]

set.seed(1)
SSC_bray_adonis <- adonis2(beta_genus ~ Inoculum, data=metadata, method="bray", permutations=999)
SSC_Tax_SC_value <- SSC_bray_adonis$R2[1]

SSC_bray_adonis <- adonis2(beta_genus ~ Condition, data=metadata, method="bray", permutations=999)
SSC_Tax_PL_value <- SSC_bray_adonis$R2[1]

#KO table
KO_SSC =read.table(paste(working_directory,"KO_tables/No_nodulators/SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object
OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)

#Sample subsetting
cond="ES"

samples_df_sub <- subset(samples_df_2, samples_df_2$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

samples_sub = sample_data(samples_df_sub_2)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)

phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )

beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

#Make PCoA plot for Bray Curtis Distance matrix
pcoa = cmdscale(beta_isolate_KO, k=3, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z") 
eig = pcoa$eig

points = merge(points,samples_df_sub_2, by = "row.names")
rownames(points) <- points$Row.names
points <- points %>% dplyr::select (-Row.names)

metadata=points[,-c(1,2,3)]

set.seed(1)
SSC_bray_KO_adonis <- adonis2(beta_isolate_KO ~ Inoculum, data=metadata, method="bray", permutations=999)
SSC_KO_SC_value <- SSC_bray_KO_adonis$R2[1]

SSC_bray_KO_adonis <- adonis2(beta_isolate_KO ~ Condition, data=metadata, method="bray", permutations=999)
SSC_KO_PL_value <- SSC_bray_KO_adonis$R2[1]

SSC_values_SC_nod <- c("Inoculum - no nodulators", SSC_Tax_SC_value,SSC_KO_SC_value)
SSC_values_PL_nod <- c("Host - no nodulators", SSC_Tax_PL_value,SSC_KO_PL_value)

#otu table - no Rhizobacter
norm_SSC =read.table(paste(working_directory,"Isolate_tables/No_rhizobacter/SSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
#Taxonomy table
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")
#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
colnames(samples_df_2)[5]="Nutrient"
samples_df_2$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df_2$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

sapply(tax_df, function(x) length(unique(x)))

#Set the OTU, TAX and sample data for making phyloseq object
OTU = otu_table(as.matrix(norm_SSC),taxa_are_rows = TRUE)
#TAX = tax_table(tax_mat)
TAX = tax_table(as.matrix(tax_df_2))

#Sample subsetting
cond="ES"
samples_df_sub <- subset(samples_df_2, samples_df_2$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

samples_sub = sample_data(samples_df_sub_2)

phylo_sub = phyloseq(OTU,TAX, samples_sub)

phylo_sub_RA=microbiome::transform(x = phylo_sub, transform = "compositional" )

#Agglomerate to phylum-level and rename
ps_genus <- phyloseq::tax_glom(phylo_sub_RA, "Genus")
phyloseq::taxa_names(ps_genus) <- phyloseq::tax_table(ps_genus)[, "Genus"]

#Bray Curtis distance matrix
beta_genus <- as.matrix(vegdist(t(ps_genus@otu_table@.Data), method = "bray", diag = T))

#Make PCoA plot for Bray Curtis Distance matrix
pcoa = cmdscale(beta_genus, k=3, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z") 
eig = pcoa$eig

points = merge(points,samples_df_sub_2, by = "row.names")
rownames(points) <- points$Row.names
points <- points %>% dplyr::select (-Row.names)
metadata=points[,-c(1,2,3)]

set.seed(1)
SSC_bray_adonis <- adonis2(beta_genus ~ Inoculum, data=metadata, method="bray", permutations=999)
SSC_Tax_SC_value <- SSC_bray_adonis$R2[1]

SSC_bray_adonis <- adonis2(beta_genus ~ Condition, data=metadata, method="bray", permutations=999)
SSC_Tax_PL_value <- SSC_bray_adonis$R2[1]

#KO table
KO_SSC =read.table(paste(working_directory,"KO_tables/No_rhizobacter/SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object
OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)

#Sample subsetting
cond="ES"

samples_df_sub <- subset(samples_df_2, samples_df_2$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

samples_sub = sample_data(samples_df_sub_2)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)

phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )

beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

#Make PCoA plot for Bray Curtis Distance matrix
pcoa = cmdscale(beta_isolate_KO, k=3, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z") 
eig = pcoa$eig

points = merge(points,samples_df_sub_2, by = "row.names")
rownames(points) <- points$Row.names
points <- points %>% dplyr::select (-Row.names)

metadata=points[,-c(1,2,3)]

set.seed(1)
SSC_bray_KO_adonis <- adonis2(beta_isolate_KO ~ Inoculum, data=metadata, method="bray", permutations=999)
SSC_KO_SC_value <- SSC_bray_KO_adonis$R2[1]

SSC_bray_KO_adonis <- adonis2(beta_isolate_KO ~ Condition, data=metadata, method="bray", permutations=999)
SSC_KO_PL_value <- SSC_bray_KO_adonis$R2[1]

SSC_values_SC_rhizo <- c("Inoculum - no ", SSC_Tax_SC_value,SSC_KO_SC_value)
SSC_values_PL_rhizo <- c("Host - no ", SSC_Tax_PL_value,SSC_KO_PL_value)

#otu table - no Rhizobacter & Nodulators
norm_SSC =read.table(paste(working_directory,"Isolate_tables/No_dominances/SSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
#Taxonomy table
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")
#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
colnames(samples_df_2)[5]="Nutrient"
samples_df_2$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df_2$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

sapply(tax_df, function(x) length(unique(x)))

#Set the OTU, TAX and sample data for making phyloseq object
OTU = otu_table(as.matrix(norm_SSC),taxa_are_rows = TRUE)
#TAX = tax_table(tax_mat)
TAX = tax_table(as.matrix(tax_df_2))

#Sample subsetting

cond="ES"
samples_df_sub <- subset(samples_df_2, samples_df_2$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

samples_sub = sample_data(samples_df_sub_2)

phylo_sub = phyloseq(OTU,TAX, samples_sub)

phylo_sub_RA=microbiome::transform(x = phylo_sub, transform = "compositional" )

#Agglomerate to phylum-level and rename
ps_genus <- phyloseq::tax_glom(phylo_sub_RA, "Genus")
phyloseq::taxa_names(ps_genus) <- phyloseq::tax_table(ps_genus)[, "Genus"]

#Bray Curtis distance matrix
beta_genus <- as.matrix(vegdist(t(ps_genus@otu_table@.Data), method = "bray", diag = T))

#Make PCoA plot for Bray Curtis Distance matrix
pcoa = cmdscale(beta_genus, k=3, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z") 
eig = pcoa$eig

points = merge(points,samples_df_sub_2, by = "row.names")
rownames(points) <- points$Row.names
points <- points %>% dplyr::select (-Row.names)
metadata=points[,-c(1,2,3)]

set.seed(1)
SSC_bray_adonis <- adonis2(beta_genus ~ Inoculum, data=metadata, method="bray", permutations=999)
SSC_Tax_SC_value <- SSC_bray_adonis$R2[1]

SSC_bray_adonis <- adonis2(beta_genus ~ Condition, data=metadata, method="bray", permutations=999)
SSC_Tax_PL_value <- SSC_bray_adonis$R2[1]

#KO table
KO_SSC =read.table(paste(working_directory,"KO_tables/No_dominances/SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object
OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)

#Sample subsetting
cond="ES"

samples_df_sub <- subset(samples_df_2, samples_df_2$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

samples_sub = sample_data(samples_df_sub_2)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)

phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )

beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

#Make PCoA plot for Bray Curtis Distance matrix
pcoa = cmdscale(beta_isolate_KO, k=3, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z") 
eig = pcoa$eig

points = merge(points,samples_df_sub_2, by = "row.names")
rownames(points) <- points$Row.names
points <- points %>% dplyr::select (-Row.names)

metadata=points[,-c(1,2,3)]

set.seed(1)
SSC_bray_KO_adonis <- adonis2(beta_isolate_KO ~ Inoculum, data=metadata, method="bray", permutations=999)
SSC_KO_SC_value <- SSC_bray_KO_adonis$R2[1]

SSC_bray_KO_adonis <- adonis2(beta_isolate_KO ~ Condition, data=metadata, method="bray", permutations=999)
SSC_KO_PL_value <- SSC_bray_KO_adonis$R2[1]

SSC_values_SC_rhizo_nod <- c("Inoculum - no dominators", SSC_Tax_SC_value,SSC_KO_SC_value)
SSC_values_PL_rhizo_nod <- c("Host - no dominators", SSC_Tax_PL_value,SSC_KO_PL_value)

R2_values_2 <- rbind(R2_values_2, SSC_values_SC,SSC_values_PL,SSC_values_SC_nod,SSC_values_PL_nod, SSC_values_SC_rhizo, SSC_values_PL_rhizo, SSC_values_SC_rhizo_nod,SSC_values_PL_rhizo_nod)
R2_values_2$Tax_R2 <- as.numeric(R2_values_2$Tax_R2)
R2_values_2$Func_R2 <- as.numeric(R2_values_2$Func_R2)
R2_values_2$Simulation <- "Simulation"
R2_values_2$Simulation[R2_values_2$run == "Inoculum"] <- "SSC_SC"
R2_values_2$Simulation[R2_values_2$run == "Inoculum - no nodulators"] <- "SSC_no_nod_SC"
R2_values_2$Simulation[R2_values_2$run == "Host"] <- "SSC_PL"
R2_values_2$Simulation[R2_values_2$run == "Host - no nodulators"] <- "SSC_no_nod_PL"
R2_values_2$Simulation[R2_values_2$run == "Inoculum - no "] <- "SSC_SC_no_rhizo"
R2_values_2$Simulation[R2_values_2$run == "Inoculum - no dominators"] <- "SSC_no_rhizo_and_nod_SC"
R2_values_2$Simulation[R2_values_2$run == "Host - no "] <- "SSC_PL_no_rhizo"
R2_values_2$Simulation[R2_values_2$run == "Host - no dominators"] <- "SSC_no_rhizo_and_nod_PL"

R2_values_2$Simulation_2 <- "Simulation"
R2_values_2$Simulation_2[R2_values_2$run == "Inoculum"] <- "SSC_SC"
R2_values_2$Simulation_2[R2_values_2$run == "Inoculum - no nodulators"] <- "SSC_SC"
R2_values_2$Simulation_2[R2_values_2$run == "Host"] <- "SSC_PL"
R2_values_2$Simulation_2[R2_values_2$run == "Host - no nodulators"] <- "SSC_PL"
R2_values_2$Simulation_2[R2_values_2$run == "Inoculum - no "] <- "SSC_SC"
R2_values_2$Simulation_2[R2_values_2$run == "Inoculum - no dominators"] <- "SSC_SC"
R2_values_2$Simulation_2[R2_values_2$run == "Host - no "] <- "SSC_PL"
R2_values_2$Simulation_2[R2_values_2$run == "Host - no dominators"] <- "SSC_PL"

R2_values_2$Simulation_2 <- factor(R2_values_2$Simulation_2, levels = c("Simulation", "SSC_SC", "SSC_PL"))

checkies <- c(3,6,9)

R2_values_3 <- R2_values_2[R2_values_2$kmeans %in% checkies, ]
R2_values_4 <- unique(R2_values_2[R2_values_2$Simulation != "Simulation", ])

R2_values_5 <- rbind(R2_values_3, R2_values_4)

data_stat <- data.frame()

groups <- R2_values_5$run[3001:3008]

for (group in groups){
  model <- lm(Tax_R2 ~ Func_R2, data = R2_values_5)
  R2_values_5$residuals <- residuals(model)
  std_error <- summary(model)$sigma
  specific_point <- R2_values_5[R2_values_5$run == paste(group), ]
  predicted_value <- predict(model, newdata = specific_point)
  residual <- specific_point$Tax_R2 - predicted_value
  standardized_residual <- residual / std_error
  p_value <- 2 * (1 - pnorm(abs(standardized_residual)) ) # Two-tailed test
  hop <- data.frame(t(data.frame(c(paste(group), p_value))))
  row.names(hop) <- NULL
  colnames(hop) <- c("Group","pvalue")
  data_stat <- rbind(data_stat,hop )
}

#Statistical deviations of trend
write.table(data_stat, paste(working_directory, "Table_S3_R2_simulation_stats.txt", sep =""), sep= "\t", quote =F, row.names =F, col.names=T)

R2_correlations <- ggscatter(R2_values_5, x="Tax_R2", y="Func_R2", color = "Simulation_2", conf.int = F,alpha = 0.7,
                             palette = c(Simulation_2 = "black", SSC_SC = "#F8766D",SSC_PL = "#80b006")) +
  stat_cor(label.x = 0.05, label.y = 0.4) +
  geom_smooth(method=lm, level = 0.99999999999999) +
  geom_text_repel(aes(label=ifelse(R2_values_5$Simulation_2 == "SSC_SC",yes = as.character(R2_values_5$run), '')),size=4,max.overlaps = Inf) +
  geom_text_repel(aes(label=ifelse(R2_values_5$Simulation_2 == "SSC_PL",yes = as.character(R2_values_5$run), '')),size=4, max.overlaps = Inf) +
  ggtitle("R2 Correlation between taxonomy and functionality") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  ylab(expression(paste("Functionality R"^2, sep =""))) + 
  xlab(expression(paste("Taxonomy R"^2, sep = ""))) + 
  theme(axis.text.x = element_text(size = 14), axis.title = element_text(size = 18), axis.text.y = element_text(size=14), plot.title = element_text(size=24)) +
  theme(legend.position = "none") 
R2_correlations

pdf(paste(results.dir,"Figure_S12_R2_correlation_plot_simulation_2.pdf", sep=""), width=10, height=7)
print(R2_correlations)
dev.off()

###Figure 3ab & Table S5 - PieDonut plots =====
combined_df_syncom_4_dom <- data.frame(matrix(NA, ncol = 11))
colnames(combined_df_syncom_4_dom) <- c("Df", "SumofSqs", "R2", "F", "Pr(>F)", "Variable","Subset", "Test", "Rank", "Dominance", "Drop_out")
combined_df_syncom_4_dom <- combined_df_syncom_4_dom[-1,]

combined_df_syncom_5_dom <- combined_df_syncom_4_dom

#otu table
norm_SSC=read.table(paste(working_directory,"Isolate_tables/Original/SSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
#Taxonomy table
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")
#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
colnames(samples_df_2)[5]="Nutrient"
samples_df_2$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df_2$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

sapply(tax_df, function(x) length(unique(x)))

#Set the OTU, TAX and sample data for making phyloseq object
OTU = otu_table(as.matrix(norm_SSC),taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(tax_df_2))

#Sample subsetting
cond="ES"
samples_df_sub <- subset(samples_df_2, samples_df_2$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

samples_sub = sample_data(samples_df_sub_2)

phylo_sub = phyloseq(OTU,TAX, samples_sub)

phylo_sub_RA=microbiome::transform(x = phylo_sub, transform = "compositional" )

#Agglomerate to phylum-level and rename
ps_phylum <- phyloseq::tax_glom(phylo_sub_RA, "Phylum")
ps_class <- phyloseq::tax_glom(phylo_sub_RA, "Class")
ps_order <- phyloseq::tax_glom(phylo_sub_RA, "Order")
ps_family <- phyloseq::tax_glom(phylo_sub_RA, "Family")
ps_genus <- phyloseq::tax_glom(phylo_sub_RA, "Genus")

phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]
phyloseq::taxa_names(ps_class) <- phyloseq::tax_table(ps_class)[, "Class"]
phyloseq::taxa_names(ps_order) <- phyloseq::tax_table(ps_order)[, "Order"]
phyloseq::taxa_names(ps_family) <- phyloseq::tax_table(ps_family)[, "Family"]
phyloseq::taxa_names(ps_genus) <- phyloseq::tax_table(ps_genus)[, "Genus"]

#Bray Curtis distance matrix
beta_phylum <- as.matrix(vegdist(t(ps_phylum@otu_table@.Data), method = "bray", diag = T))
beta_class <- as.matrix(vegdist(t(ps_class@otu_table@.Data), method = "bray", diag = T))
beta_order <- as.matrix(vegdist(t(ps_order@otu_table@.Data), method = "bray", diag = T))
beta_family <- as.matrix(vegdist(t(ps_family@otu_table@.Data), method = "bray", diag = T))
beta_genus <- as.matrix(vegdist(t(ps_genus@otu_table@.Data), method = "bray", diag = T))
beta_isolate <- as.matrix(vegdist(t(phylo_sub_RA@otu_table@.Data), method = "bray", diag = T))

#otu table
KO_SSC=read.table(paste(working_directory,"KO_tables/Original/SSC.tsv", sep =""), header=TRUE,sep="\t", row.names = 1)

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object
OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)

#Sample subsetting
cond="ES"

samples_df_sub <- subset(samples_df_2, samples_df_2$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

samples_sub = sample_data(samples_df_sub_2)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)

phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )

beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

#Data without Dominances
#otu table
norm_SSC=read.table(paste(working_directory,"Isolate_tables/No_dominances/SSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
#Taxonomy table
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")
#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
colnames(samples_df_2)[5]="Nutrient"
samples_df_2$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df_2$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

sapply(tax_df, function(x) length(unique(x)))

#Set the OTU, TAX and sample data for making phyloseq object
OTU = otu_table(as.matrix(norm_SSC),taxa_are_rows = TRUE)
#TAX = tax_table(tax_mat)
TAX = tax_table(as.matrix(tax_df_2))

#Sample subsetting
cond="ES"
samples_df_sub <- subset(samples_df_2, samples_df_2$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

samples_sub = sample_data(samples_df_sub_2)

phylo_sub = phyloseq(OTU,TAX, samples_sub)

phylo_sub_RA=microbiome::transform(x = phylo_sub, transform = "compositional" )

#Agglomerate to phylum-level and rename
ps_phylum <- phyloseq::tax_glom(phylo_sub_RA, "Phylum")
ps_class <- phyloseq::tax_glom(phylo_sub_RA, "Class")
ps_order <- phyloseq::tax_glom(phylo_sub_RA, "Order")
ps_family <- phyloseq::tax_glom(phylo_sub_RA, "Family")
ps_genus <- phyloseq::tax_glom(phylo_sub_RA, "Genus")

phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]
phyloseq::taxa_names(ps_class) <- phyloseq::tax_table(ps_class)[, "Class"]
phyloseq::taxa_names(ps_order) <- phyloseq::tax_table(ps_order)[, "Order"]
phyloseq::taxa_names(ps_family) <- phyloseq::tax_table(ps_family)[, "Family"]
phyloseq::taxa_names(ps_genus) <- phyloseq::tax_table(ps_genus)[, "Genus"]

#Bray Curtis distance matrix
beta_phylum_dom <- as.matrix(vegdist(t(ps_phylum@otu_table@.Data), method = "bray", diag = T))
beta_class_dom <- as.matrix(vegdist(t(ps_class@otu_table@.Data), method = "bray", diag = T))
beta_order_dom <- as.matrix(vegdist(t(ps_order@otu_table@.Data), method = "bray", diag = T))
beta_family_dom <- as.matrix(vegdist(t(ps_family@otu_table@.Data), method = "bray", diag = T))
beta_genus_dom <- as.matrix(vegdist(t(ps_genus@otu_table@.Data), method = "bray", diag = T))
beta_isolate_dom <- as.matrix(vegdist(t(phylo_sub_RA@otu_table@.Data), method = "bray", diag = T))

#otu table
KO_SSC=read.table(paste(working_directory,"KO_tables/No_dominances/SSC.tsv", sep =""), header=TRUE,sep="\t", row.names = 1)

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object
OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)

#Sample subsetting
cond="ES"

samples_df_sub <- subset(samples_df_2, samples_df_2$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

samples_sub = sample_data(samples_df_sub_2)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)

phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )

beta_isolate_KO_dom <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

beta_distance_list=list(beta_isolate, beta_genus, beta_family, beta_order, beta_class, beta_phylum)

for (i in 1:length(beta_distance_list)) {
  Bray_curtis_df=beta_distance_list[[i]]
  
  #Subsetting Bray Curtis matrix based on the sample subsetting on the lines above
  Bray_match <- as.matrix(as.logical(match(row.names(Bray_curtis_df),row.names(samples_df_sub))))
  Bray_Curtis_SSC <- subset(Bray_curtis_df, Bray_match)
  Bray_Curtis_SSC_T <- t(Bray_Curtis_SSC)
  Bray_Curtis_SSC_T <- subset(Bray_Curtis_SSC_T,Bray_match)
  Bray_Curtis_SSC <- t(Bray_Curtis_SSC_T)
  
  #Make PCoA plot for Bray Curtis Distance matrix
  pcoa = cmdscale(Bray_curtis_df, k=3, eig=T)
  points = as.data.frame(pcoa$points)
  colnames(points) = c("x", "y", "z") 
  eig = pcoa$eig
  
  points = merge(points,samples_df_sub_2, by = "row.names")
  rownames(points) <- points$Row.names
  points <- points %>% dplyr::select (-Row.names)
  
  points$Condition <- factor(points$Condition, levels = c("At","Hv", "Lj"))
  points$Inoculum <- factor(points$Inoculum, levels = c("SSC","AtSC", "HvSC", "LjSC","NS"))
  points$Nutrient <- factor(points$Nutrient, levels = c("low", "high"))
  points$Experiment  <- factor(points$Experiment, levels = c("R1", "R2"))
  
  points$Compartment <- factor(points$Compartment, levels = c("NP", "RZ","ES","dom"))
  points$Growth_condition <- factor(points$Growth_condition, levels = c("GC1_Lj","GC2_At","GC3_Hv"))
  points$Plant_inoc=paste(points$Plant_species,"_",points$Inoculum)
  metadata=points[,-c(1,2,3)]
}

comb_plant_list=list()
comb_syncom_list=list()

Taxo_order <- c("Isolate", "Genus", "Family", "Order", "Phylum","KO_functions")

for (i in 1:6) {
  #Per SynCom
  metadata_AtSC <- metadata[metadata$Inoculum == "AtSC",]
  taxo_beta_AtSC <- beta_distance_list[[i]][row.names(beta_distance_list[[i]]) %in% row.names(metadata_AtSC),row.names(beta_distance_list[[i]]) %in% row.names(metadata_AtSC)]
  taxo_adonis_AtSC <- adonis2(taxo_beta_AtSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_AtSC, method="bray", permutations=999)
  
  lala=pairwiseAdonis::pairwise.adonis(x =taxo_beta_AtSC, factors = metadata_AtSC$Condition, perm = 999 )
  
  AtSC_taxo_df <- cbind(taxo_adonis_AtSC,Variable = rownames(taxo_adonis_AtSC), Subset = "AtSC", Test = "Taxonomy", Rank=Taxo_order[i])
  
  metadata_HvSC <- metadata[metadata$Inoculum == "HvSC",]
  taxo_beta_HvSC <- beta_distance_list[[i]][row.names(beta_distance_list[[i]]) %in% row.names(metadata_HvSC),row.names(beta_distance_list[[i]]) %in% row.names(metadata_HvSC)]
  taxo_adonis_HvSC <- adonis2(taxo_beta_HvSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_HvSC, method="bray", permutations=999)
  HvSC_taxo_df <- cbind(taxo_adonis_HvSC,Variable = rownames(taxo_adonis_HvSC), Subset = "HvSC", Test = "Taxonomy", Rank=Taxo_order[i])
  
  lala_HvSC=pairwiseAdonis::pairwise.adonis(x =taxo_beta_HvSC, factors = metadata_HvSC$Condition, perm = 999 )
  
  metadata_LjSC <- metadata[metadata$Inoculum == "LjSC",]
  taxo_beta_LjSC <- beta_distance_list[[i]][row.names(beta_distance_list[[i]]) %in% row.names(metadata_LjSC),row.names(beta_distance_list[[i]]) %in% row.names(metadata_LjSC)]
  taxo_adonis_LjSC <- adonis2(taxo_beta_LjSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_LjSC, method="bray", permutations=999)
  LjSC_taxo_df <- cbind(taxo_adonis_LjSC,Variable = rownames(taxo_adonis_LjSC), Subset = "LjSC", Test = "Taxonomy", Rank=Taxo_order[i])
  
  lala_LjSC=pairwiseAdonis::pairwise.adonis(x =taxo_beta_LjSC, factors = metadata_LjSC$Condition, perm = 999 )
  
  metadata_SSC <- metadata[metadata$Inoculum == "SSC",]
  taxo_beta_SSC <- beta_distance_list[[i]][row.names(beta_distance_list[[i]]) %in% row.names(metadata_SSC),row.names(beta_distance_list[[i]]) %in% row.names(metadata_SSC)]
  taxo_adonis_SSC <- adonis2(taxo_beta_SSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_SSC, method="bray", permutations=999)
  SSC_taxo_df <- cbind(taxo_adonis_SSC,Variable = rownames(taxo_adonis_SSC), Subset = "SSC", Test = "Taxonomy", Rank=Taxo_order[i])
  
  lala_SSC=pairwiseAdonis::pairwise.adonis(x =taxo_beta_SSC, factors = metadata_SSC$Condition, perm = 999 )
  
  #Per plant
  metadata_At <- metadata[metadata$Condition == "At",]
  taxo_beta_At <- beta_distance_list[[i]][row.names(beta_distance_list[[i]]) %in% row.names(metadata_At),row.names(beta_distance_list[[i]]) %in% row.names(metadata_At)]
  taxo_adonis_At <- adonis2(taxo_beta_At ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_At, method="bray", permutations=999)
  At_taxo_df <- cbind(taxo_adonis_At,Variable = rownames(taxo_adonis_At), Subset = "At", Test = "Taxonomy", Rank=Taxo_order[i])
  
  lala_At=pairwiseAdonis::pairwise.adonis(x =taxo_beta_At, factors = metadata_At$Inoculum, perm = 999 )
  
  metadata_Hv <- metadata[metadata$Condition == "Hv",]
  taxo_beta_Hv <- beta_distance_list[[i]][row.names(beta_distance_list[[i]]) %in% row.names(metadata_Hv),row.names(beta_distance_list[[i]]) %in% row.names(metadata_Hv)]
  taxo_adonis_Hv <- adonis2(taxo_beta_Hv ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_Hv, method="bray", permutations=999)
  Hv_taxo_df <- cbind(taxo_adonis_Hv,Variable = rownames(taxo_adonis_Hv), Subset = "Hv", Test = "Taxonomy", Rank=Taxo_order[i])
  
  lala_Hv=pairwiseAdonis::pairwise.adonis(x =taxo_beta_Hv, factors = metadata_Hv$Inoculum, perm = 999 )
  
  metadata_Lj <- metadata[metadata$Condition == "At",]
  taxo_beta_Lj <- beta_distance_list[[i]][row.names(beta_distance_list[[i]]) %in% row.names(metadata_Lj),row.names(beta_distance_list[[i]]) %in% row.names(metadata_Lj)]
  taxo_adonis_Lj <- adonis2(taxo_beta_Lj ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_Lj, method="bray", permutations=999)
  Lj_taxo_df <- cbind(taxo_adonis_Lj,Variable = rownames(taxo_adonis_Lj), Subset = "Lj", Test = "Taxonomy", Rank=Taxo_order[i])
  
  lala_Lj=pairwiseAdonis::pairwise.adonis(x =taxo_beta_Lj, factors = metadata_Lj$Inoculum, perm = 999 )
  
  Combined_df_Syncom=rbind(AtSC_taxo_df, HvSC_taxo_df,LjSC_taxo_df,SSC_taxo_df)
  Sub_comb_syncom <- Combined_df_Syncom[!(Combined_df_Syncom$Variable %in% c("Total", "Residual", "Condition:Nutrient", "Condition:Experiment")), ]
  comb_syncom_list[[i]]=Sub_comb_syncom
  
  Combined_df_Plant=rbind(At_taxo_df,Hv_taxo_df,Lj_taxo_df)
  Sub_comb_plant <- Combined_df_Plant[!(Combined_df_Plant$Variable %in% c("Total", "Residual", "Inoculum:Nutrient", "Inoculum:Experiment")), ]
  comb_plant_list[[i]]=Sub_comb_plant
  
}

#Separation by plants/ Syncom KO Functionality 
#Per SynCom
metadata_AtSC <- metadata[metadata$Inoculum == "AtSC",]
KO_beta_AtSC <- beta_isolate_KO[row.names(beta_isolate_KO) %in% row.names(metadata_AtSC),row.names(beta_isolate_KO) %in% row.names(metadata_AtSC)]
KO_adonis_AtSC <- adonis2(KO_beta_AtSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_AtSC, method="bray", permutations=999)
AtSC_KO_df <- cbind(KO_adonis_AtSC,Variable = rownames(KO_adonis_AtSC), Subset = "AtSC", Test = "Functions", Rank="Functions")

lala_KO_AtSC=pairwiseAdonis::pairwise.adonis(x =KO_beta_AtSC, factors = metadata_AtSC$Condition, perm = 999 )

metadata_HvSC <- metadata[metadata$Inoculum == "HvSC",]
KO_beta_HvSC <- beta_isolate_KO[row.names(beta_isolate_KO) %in% row.names(metadata_HvSC),row.names(beta_isolate_KO) %in% row.names(metadata_HvSC)]
KO_adonis_HvSC <- adonis2(KO_beta_HvSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_HvSC, method="bray", permutations=999)
HvSC_KO_df <- cbind(KO_adonis_HvSC,Variable = rownames(KO_adonis_HvSC), Subset = "HvSC", Test = "Functions", Rank="Functions")

lala_KO_HvSC=pairwiseAdonis::pairwise.adonis(x =KO_beta_HvSC, factors = metadata_HvSC$Condition, perm = 999 )

metadata_LjSC <- metadata[metadata$Inoculum == "LjSC",]
KO_beta_LjSC <- beta_isolate_KO[row.names(beta_isolate_KO) %in% row.names(metadata_LjSC),row.names(beta_isolate_KO) %in% row.names(metadata_LjSC)]
KO_adonis_LjSC <- adonis2(KO_beta_LjSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_LjSC, method="bray", permutations=999)
LjSC_KO_df <- cbind(KO_adonis_LjSC,Variable = rownames(KO_adonis_LjSC), Subset = "LjSC", Test = "Functions", Rank="Functions")

lala_KO_LjSC=pairwiseAdonis::pairwise.adonis(x =KO_beta_LjSC, factors = metadata_LjSC$Condition, perm = 999 )

metadata_SSC <- metadata[metadata$Inoculum == "SSC",]
KO_beta_SSC <- beta_isolate_KO[row.names(beta_isolate_KO) %in% row.names(metadata_SSC),row.names(beta_isolate_KO) %in% row.names(metadata_SSC)]
KO_adonis_SSC <- adonis2(KO_beta_SSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_SSC, method="bray", permutations=999)
SSC_KO_df <- cbind(KO_adonis_SSC,Variable = rownames(KO_adonis_SSC), Subset = "SSC", Test = "Functions", Rank="Functions")

lala_KO_SSC=pairwiseAdonis::pairwise.adonis(x =KO_beta_SSC, factors = metadata_SSC$Condition, perm = 999 )

#Per plant
metadata_At <- metadata[metadata$Condition == "At",]
KO_beta_At <- beta_isolate_KO[row.names(beta_isolate_KO) %in% row.names(metadata_At),row.names(beta_isolate_KO) %in% row.names(metadata_At)]
KO_adonis_At <- adonis2(KO_beta_At ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_At, method="bray", permutations=999)
At_KO_df <- cbind(KO_adonis_At,Variable = rownames(KO_adonis_At), Subset = "At", Test = "Functions", Rank="Functions")

metadata_Hv <- metadata[metadata$Condition == "Hv",]
KO_beta_Hv <- beta_isolate_KO[row.names(beta_isolate_KO) %in% row.names(metadata_Hv),row.names(beta_isolate_KO) %in% row.names(metadata_Hv)]
KO_adonis_Hv <- adonis2(KO_beta_Hv ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_Hv, method="bray", permutations=999)
Hv_KO_df <- cbind(KO_adonis_Hv,Variable = rownames(KO_adonis_Hv), Subset = "Hv", Test = "Functions", Rank="Functions")

metadata_Lj <- metadata[metadata$Condition == "Lj",]
KO_beta_Lj <- beta_isolate_KO[row.names(beta_isolate_KO) %in% row.names(metadata_Lj),row.names(beta_isolate_KO) %in% row.names(metadata_Lj)]
KO_adonis_Lj <- adonis2(KO_beta_Lj ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_Lj, method="bray", permutations=999)
Lj_KO_df <- cbind(KO_adonis_Lj,Variable = rownames(KO_adonis_Lj), Subset = "Lj", Test = "Functions", Rank="Functions")

#Combining Data
combined_df_plant <- do.call(rbind, comb_plant_list)
combined_df_plant <- rbind(combined_df_plant,At_KO_df,Hv_KO_df,Lj_KO_df )
combined_df_plant <- combined_df_plant[!(combined_df_plant$Variable %in% c("Total", "Residual", "Inoculum:Nutrient", "Inoculum:Experiment", "Experiment", "Nutrient")), ]
combined_df_plant$Rank=factor(combined_df_plant$Rank, levels = c("Isolate", "Genus", "Family", "Order", "Class", "Phylum","Functions"))
combined_df_plant$Test=factor(combined_df_plant$Test, levels = c("Taxonomy", "Functions"))

combined_df_syncom <- do.call(rbind, comb_syncom_list)
combined_df_syncom <- rbind(combined_df_syncom,AtSC_KO_df,HvSC_KO_df,LjSC_KO_df,SSC_KO_df)
combined_df_syncom <- combined_df_syncom[!(combined_df_syncom$Variable %in% c("Total", "Residual", "Condition:Nutrient", "Condition:Experiment","Experiment", "Nutrient")), ]
combined_df_syncom$Rank=factor(combined_df_syncom$Rank, levels = c("Isolate", "Genus", "Family", "Order", "Class", "Phylum","Functions"))
combined_df_syncom$Test=factor(combined_df_syncom$Test, levels = c("Taxonomy", "Functions"))

#Dominances drop out data
beta_distance_list_dom=list(beta_isolate_dom, beta_genus_dom, beta_family_dom, beta_order_dom, beta_class_dom, beta_phylum_dom)
comb_plant_list_dom=list()
comb_syncom_list_dom=list()

for (i in 1:6) {
  #Per SynCom
  metadata_AtSC <- metadata[metadata$Inoculum == "AtSC",]
  taxo_beta_AtSC <- beta_distance_list_dom[[i]][row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_AtSC),row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_AtSC)]
  taxo_adonis_AtSC <- adonis2(taxo_beta_AtSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_AtSC, method="bray", permutations=999)
  
  lala=pairwiseAdonis::pairwise.adonis(x =taxo_beta_AtSC, factors = metadata_AtSC$Condition, perm = 999 )
  
  AtSC_taxo_df <- cbind(taxo_adonis_AtSC,Variable = rownames(taxo_adonis_AtSC), Subset = "AtSC", Test = "Taxonomy", Rank=Taxo_order[i])
  
  metadata_HvSC <- metadata[metadata$Inoculum == "HvSC",]
  taxo_beta_HvSC <- beta_distance_list_dom[[i]][row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_HvSC),row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_HvSC)]
  taxo_adonis_HvSC <- adonis2(taxo_beta_HvSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_HvSC, method="bray", permutations=999)
  HvSC_taxo_df <- cbind(taxo_adonis_HvSC,Variable = rownames(taxo_adonis_HvSC), Subset = "HvSC", Test = "Taxonomy", Rank=Taxo_order[i])
  
  lala_HvSC=pairwiseAdonis::pairwise.adonis(x =taxo_beta_HvSC, factors = metadata_HvSC$Condition, perm = 999 )
  
  metadata_LjSC <- metadata[metadata$Inoculum == "LjSC",]
  taxo_beta_LjSC <- beta_distance_list_dom[[i]][row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_LjSC),row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_LjSC)]
  taxo_adonis_LjSC <- adonis2(taxo_beta_LjSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_LjSC, method="bray", permutations=999)
  LjSC_taxo_df <- cbind(taxo_adonis_LjSC,Variable = rownames(taxo_adonis_LjSC), Subset = "LjSC", Test = "Taxonomy", Rank=Taxo_order[i])
  
  lala_LjSC=pairwiseAdonis::pairwise.adonis(x =taxo_beta_LjSC, factors = metadata_LjSC$Condition, perm = 999 )
  
  metadata_SSC <- metadata[metadata$Inoculum == "SSC",]
  taxo_beta_SSC <- beta_distance_list_dom[[i]][row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_SSC),row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_SSC)]
  taxo_adonis_SSC <- adonis2(taxo_beta_SSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_SSC, method="bray", permutations=999)
  SSC_taxo_df <- cbind(taxo_adonis_SSC,Variable = rownames(taxo_adonis_SSC), Subset = "SSC", Test = "Taxonomy", Rank=Taxo_order[i])
  
  lala_SSC=pairwiseAdonis::pairwise.adonis(x =taxo_beta_SSC, factors = metadata_SSC$Condition, perm = 999 )
  
  #Per plant
  metadata_At <- metadata[metadata$Condition == "At",]
  taxo_beta_At <- beta_distance_list_dom[[i]][row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_At),row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_At)]
  taxo_adonis_At <- adonis2(taxo_beta_At ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_At, method="bray", permutations=999)
  At_taxo_df <- cbind(taxo_adonis_At,Variable = rownames(taxo_adonis_At), Subset = "At", Test = "Taxonomy", Rank=Taxo_order[i])
  
  lala_At=pairwiseAdonis::pairwise.adonis(x =taxo_beta_At, factors = metadata_At$Inoculum, perm = 999 )
  
  metadata_Hv <- metadata[metadata$Condition == "Hv",]
  taxo_beta_Hv <- beta_distance_list_dom[[i]][row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_Hv),row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_Hv)]
  taxo_adonis_Hv <- adonis2(taxo_beta_Hv ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_Hv, method="bray", permutations=999)
  Hv_taxo_df <- cbind(taxo_adonis_Hv,Variable = rownames(taxo_adonis_Hv), Subset = "Hv", Test = "Taxonomy", Rank=Taxo_order[i])
  
  lala_Hv=pairwiseAdonis::pairwise.adonis(x =taxo_beta_Hv, factors = metadata_Hv$Inoculum, perm = 999 )
  
  metadata_Lj <- metadata[metadata$Condition == "At",]
  taxo_beta_Lj <- beta_distance_list_dom[[i]][row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_Lj),row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_Lj)]
  taxo_adonis_Lj <- adonis2(taxo_beta_Lj ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_Lj, method="bray", permutations=999)
  Lj_taxo_df <- cbind(taxo_adonis_Lj,Variable = rownames(taxo_adonis_Lj), Subset = "Lj", Test = "Taxonomy", Rank=Taxo_order[i])
  
  lala_Lj=pairwiseAdonis::pairwise.adonis(x =taxo_beta_Lj, factors = metadata_Lj$Inoculum, perm = 999 )
  
  Combined_df_Syncom=rbind(AtSC_taxo_df, HvSC_taxo_df,LjSC_taxo_df,SSC_taxo_df)
  Sub_comb_syncom <- Combined_df_Syncom[!(Combined_df_Syncom$Variable %in% c("Total", "Residual", "Condition:Nutrient", "Condition:Experiment")), ]
  comb_syncom_list_dom[[i]]=Sub_comb_syncom
  
  Combined_df_Plant=rbind(At_taxo_df,Hv_taxo_df,Lj_taxo_df)
  Sub_comb_plant <- Combined_df_Plant[!(Combined_df_Plant$Variable %in% c("Total", "Residual", "Inoculum:Nutrient", "Inoculum:Experiment")), ]
  comb_plant_list_dom[[i]]=Sub_comb_plant
  
}

#Separation by plants/ Syncom KO Functionality 
#Per SynCom
metadata_AtSC <- metadata[metadata$Inoculum == "AtSC",]
KO_beta_AtSC <- beta_isolate_KO_dom[row.names(beta_isolate_KO_dom) %in% row.names(metadata_AtSC),row.names(beta_isolate_KO_dom) %in% row.names(metadata_AtSC)]
KO_adonis_AtSC <- adonis2(KO_beta_AtSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_AtSC, method="bray", permutations=999)
AtSC_KO_df <- cbind(KO_adonis_AtSC,Variable = rownames(KO_adonis_AtSC), Subset = "AtSC", Test = "Functions", Rank="Functions")

lala_KO_AtSC=pairwiseAdonis::pairwise.adonis(x =KO_beta_AtSC, factors = metadata_AtSC$Condition, perm = 999 )

metadata_HvSC <- metadata[metadata$Inoculum == "HvSC",]
KO_beta_HvSC <- beta_isolate_KO_dom[row.names(beta_isolate_KO_dom) %in% row.names(metadata_HvSC),row.names(beta_isolate_KO_dom) %in% row.names(metadata_HvSC)]
KO_adonis_HvSC <- adonis2(KO_beta_HvSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_HvSC, method="bray", permutations=999)
HvSC_KO_df <- cbind(KO_adonis_HvSC,Variable = rownames(KO_adonis_HvSC), Subset = "HvSC", Test = "Functions", Rank="Functions")

lala_KO_HvSC=pairwiseAdonis::pairwise.adonis(x =KO_beta_HvSC, factors = metadata_HvSC$Condition, perm = 999 )

metadata_LjSC <- metadata[metadata$Inoculum == "LjSC",]
KO_beta_LjSC <- beta_isolate_KO_dom[row.names(beta_isolate_KO_dom) %in% row.names(metadata_LjSC),row.names(beta_isolate_KO_dom) %in% row.names(metadata_LjSC)]
KO_adonis_LjSC <- adonis2(KO_beta_LjSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_LjSC, method="bray", permutations=999)
LjSC_KO_df <- cbind(KO_adonis_LjSC,Variable = rownames(KO_adonis_LjSC), Subset = "LjSC", Test = "Functions", Rank="Functions")

lala_KO_LjSC=pairwiseAdonis::pairwise.adonis(x =KO_beta_LjSC, factors = metadata_LjSC$Condition, perm = 999 )

metadata_SSC <- metadata[metadata$Inoculum == "SSC",]
KO_beta_SSC <- beta_isolate_KO_dom[row.names(beta_isolate_KO_dom) %in% row.names(metadata_SSC),row.names(beta_isolate_KO_dom) %in% row.names(metadata_SSC)]
KO_adonis_SSC <- adonis2(KO_beta_SSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_SSC, method="bray", permutations=999)
SSC_KO_df <- cbind(KO_adonis_SSC,Variable = rownames(KO_adonis_SSC), Subset = "SSC", Test = "Functions", Rank="Functions")

lala_KO_SSC=pairwiseAdonis::pairwise.adonis(x =KO_beta_SSC, factors = metadata_SSC$Condition, perm = 999 )

#Per plant
metadata_At <- metadata[metadata$Condition == "At",]
KO_beta_At <- beta_isolate_KO_dom[row.names(beta_isolate_KO_dom) %in% row.names(metadata_At),row.names(beta_isolate_KO_dom) %in% row.names(metadata_At)]
KO_adonis_At <- adonis2(KO_beta_At ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_At, method="bray", permutations=999)
At_KO_df <- cbind(KO_adonis_At,Variable = rownames(KO_adonis_At), Subset = "At", Test = "Functions", Rank="Functions")

metadata_Hv <- metadata[metadata$Condition == "Hv",]
KO_beta_Hv <- beta_isolate_KO_dom[row.names(beta_isolate_KO_dom) %in% row.names(metadata_Hv),row.names(beta_isolate_KO_dom) %in% row.names(metadata_Hv)]
KO_adonis_Hv <- adonis2(KO_beta_Hv ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_Hv, method="bray", permutations=999)
Hv_KO_df <- cbind(KO_adonis_Hv,Variable = rownames(KO_adonis_Hv), Subset = "Hv", Test = "Functions", Rank="Functions")

metadata_Lj <- metadata[metadata$Condition == "Lj",]
KO_beta_Lj <- beta_isolate_KO_dom[row.names(beta_isolate_KO_dom) %in% row.names(metadata_Lj),row.names(beta_isolate_KO_dom) %in% row.names(metadata_Lj)]
KO_adonis_Lj <- adonis2(KO_beta_Lj ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_Lj, method="bray", permutations=999)
Lj_KO_df <- cbind(KO_adonis_Lj,Variable = rownames(KO_adonis_Lj), Subset = "Lj", Test = "Functions", Rank="Functions")

#Combining Data
combined_df_plant_dom <- do.call(rbind, comb_plant_list_dom)
combined_df_plant_dom <- rbind(combined_df_plant_dom,At_KO_df,Hv_KO_df,Lj_KO_df )
combined_df_plant_dom <- combined_df_plant_dom[!(combined_df_plant_dom$Variable %in% c("Total", "Residual", "Inoculum:Nutrient", "Inoculum:Experiment", "Experiment", "Nutrient")), ]
combined_df_plant_dom$Rank=factor(combined_df_plant_dom$Rank, levels = c("Isolate", "Genus", "Family", "Order", "Class", "Phylum","Functions"))
combined_df_plant_dom$Test=factor(combined_df_plant_dom$Test, levels = c("Taxonomy", "Functions"))

combined_df_syncom_dom <- do.call(rbind, comb_syncom_list_dom)
combined_df_syncom_dom <- rbind(combined_df_syncom_dom,AtSC_KO_df,HvSC_KO_df,LjSC_KO_df,SSC_KO_df)
combined_df_syncom_dom <- combined_df_syncom_dom[!(combined_df_syncom_dom$Variable %in% c("Total", "Residual", "Condition:Nutrient", "Condition:Experiment","Experiment", "Nutrient")), ]
combined_df_syncom_dom$Rank=factor(combined_df_syncom_dom$Rank, levels = c("Isolate", "Genus", "Family", "Order", "Class", "Phylum","Functions"))
combined_df_syncom_dom$Test=factor(combined_df_syncom_dom$Test, levels = c("Taxonomy", "Functions"))

beta_distance_list_dom=list(beta_isolate_dom, beta_genus_dom, beta_family_dom, beta_order_dom, beta_class_dom)
comb_plant_list_dom=list()
comb_syncom_list_dom=list()

beta_distance_list=list(beta_isolate, beta_genus, beta_family, beta_order, beta_class)
Plant_list=c("None","Arabidopsis", "Barley", "Lotus")
comb_plant_list=list()
comb_syncom_list=list()
line_plot_list=list()

metadata$Condition <- as.character(metadata$Condition)
metadata$Condition[metadata$Condition == "At"] <- "Arabidopsis"
metadata$Condition[metadata$Condition == "Hv"] <- "Barley"
metadata$Condition[metadata$Condition == "Lj"] <- "Lotus"

#original table
for (j in 1:4) {
  for (i in 1:5) {
    #Per SynCom
    metadata_AtSC <- metadata[metadata$Inoculum == "AtSC",]
    metadata_AtSC <- metadata_AtSC[metadata_AtSC$Condition != Plant_list[j],]
    taxo_beta_AtSC <- beta_distance_list[[i]][row.names(beta_distance_list[[i]]) %in% row.names(metadata_AtSC),row.names(beta_distance_list[[i]]) %in% row.names(metadata_AtSC)]
    taxo_adonis_AtSC <- adonis2(taxo_beta_AtSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_AtSC, method="bray", permutations=999)
    AtSC_taxo_df <- cbind(taxo_adonis_AtSC,Variable = rownames(taxo_adonis_AtSC), Subset = "AtSC", Test = "Taxonomy", Rank=Taxo_order[i])
    
    metadata_HvSC <- metadata[metadata$Inoculum == "HvSC",]
    metadata_HvSC <- metadata_HvSC[metadata_HvSC$Condition != Plant_list[j],]
    taxo_beta_HvSC <- beta_distance_list[[i]][row.names(beta_distance_list[[i]]) %in% row.names(metadata_HvSC),row.names(beta_distance_list[[i]]) %in% row.names(metadata_HvSC)]
    taxo_adonis_HvSC <- adonis2(taxo_beta_HvSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_HvSC, method="bray", permutations=999)
    HvSC_taxo_df <- cbind(taxo_adonis_HvSC,Variable = rownames(taxo_adonis_HvSC), Subset = "HvSC", Test = "Taxonomy", Rank=Taxo_order[i])
    
    metadata_LjSC <- metadata[metadata$Inoculum == "LjSC",]
    metadata_LjSC <- metadata_LjSC[metadata_LjSC$Condition != Plant_list[j],]
    taxo_beta_LjSC <- beta_distance_list[[i]][row.names(beta_distance_list[[i]]) %in% row.names(metadata_LjSC),row.names(beta_distance_list[[i]]) %in% row.names(metadata_LjSC)]
    taxo_adonis_LjSC <- adonis2(taxo_beta_LjSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_LjSC, method="bray", permutations=999)
    LjSC_taxo_df <- cbind(taxo_adonis_LjSC,Variable = rownames(taxo_adonis_LjSC), Subset = "LjSC", Test = "Taxonomy", Rank=Taxo_order[i])
    
    metadata_SSC <- metadata[metadata$Inoculum == "SSC",]
    metadata_SSC <- metadata_SSC[metadata_SSC$Condition != Plant_list[j],]
    taxo_beta_SSC <- beta_distance_list[[i]][row.names(beta_distance_list[[i]]) %in% row.names(metadata_SSC),row.names(beta_distance_list[[i]]) %in% row.names(metadata_SSC)]
    taxo_adonis_SSC <- adonis2(taxo_beta_SSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_SSC, method="bray", permutations=999)
    SSC_taxo_df <- cbind(taxo_adonis_SSC,Variable = rownames(taxo_adonis_SSC), Subset = "SSC", Test = "Taxonomy", Rank=Taxo_order[i])
    
    Combined_df_Syncom=rbind(AtSC_taxo_df, HvSC_taxo_df,LjSC_taxo_df,SSC_taxo_df)
    Sub_comb_syncom <- Combined_df_Syncom[!(Combined_df_Syncom$Variable %in% c("Total", "Residual", "Condition:Nutrient", "Condition:Experiment")), ]
    comb_syncom_list[[i]]=Sub_comb_syncom
  }
  
  for (i in 1:5) {
    #Per SynCom
    metadata_AtSC <- metadata[metadata$Inoculum == "AtSC",]
    metadata_AtSC <- metadata_AtSC[metadata_AtSC$Condition != Plant_list[j],]
    taxo_beta_AtSC <- beta_distance_list_dom[[i]][row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_AtSC),row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_AtSC)]
    taxo_adonis_AtSC <- adonis2(taxo_beta_AtSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_AtSC, method="bray", permutations=999)
    AtSC_taxo_df <- cbind(taxo_adonis_AtSC,Variable = rownames(taxo_adonis_AtSC), Subset = "AtSC", Test = "Taxonomy", Rank=Taxo_order[i])
    
    metadata_HvSC <- metadata[metadata$Inoculum == "HvSC",]
    metadata_HvSC <- metadata_HvSC[metadata_HvSC$Condition != Plant_list[j],]
    taxo_beta_HvSC <- beta_distance_list_dom[[i]][row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_HvSC),row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_HvSC)]
    taxo_adonis_HvSC <- adonis2(taxo_beta_HvSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_HvSC, method="bray", permutations=999)
    HvSC_taxo_df <- cbind(taxo_adonis_HvSC,Variable = rownames(taxo_adonis_HvSC), Subset = "HvSC", Test = "Taxonomy", Rank=Taxo_order[i])
    
    metadata_LjSC <- metadata[metadata$Inoculum == "LjSC",]
    metadata_LjSC <- metadata_LjSC[metadata_LjSC$Condition != Plant_list[j],]
    taxo_beta_LjSC <- beta_distance_list_dom[[i]][row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_LjSC),row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_LjSC)]
    taxo_adonis_LjSC <- adonis2(taxo_beta_LjSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_LjSC, method="bray", permutations=999)
    LjSC_taxo_df <- cbind(taxo_adonis_LjSC,Variable = rownames(taxo_adonis_LjSC), Subset = "LjSC", Test = "Taxonomy", Rank=Taxo_order[i])
    
    metadata_SSC <- metadata[metadata$Inoculum == "SSC",]
    metadata_SSC <- metadata_SSC[metadata_SSC$Condition != Plant_list[j],]
    taxo_beta_SSC <- beta_distance_list_dom[[i]][row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_SSC),row.names(beta_distance_list_dom[[i]]) %in% row.names(metadata_SSC)]
    taxo_adonis_SSC <- adonis2(taxo_beta_SSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_SSC, method="bray", permutations=999)
    SSC_taxo_df <- cbind(taxo_adonis_SSC,Variable = rownames(taxo_adonis_SSC), Subset = "SSC", Test = "Taxonomy", Rank=Taxo_order[i])
    
    Combined_df_Syncom=rbind(AtSC_taxo_df, HvSC_taxo_df,LjSC_taxo_df,SSC_taxo_df)
    Sub_comb_syncom <- Combined_df_Syncom[!(Combined_df_Syncom$Variable %in% c("Total", "Residual", "Condition:Nutrient", "Condition:Experiment")), ]
    comb_syncom_list_dom[[i]]=Sub_comb_syncom
    
  }
  
  #Per SynCom
  metadata_AtSC <- metadata[metadata$Inoculum == "AtSC",]
  metadata_AtSC <- metadata_AtSC[metadata_AtSC$Condition != Plant_list[j],]
  KO_beta_AtSC <- beta_isolate_KO[row.names(beta_isolate_KO) %in% row.names(metadata_AtSC),row.names(beta_isolate_KO) %in% row.names(metadata_AtSC)]
  KO_adonis_AtSC <- adonis2(KO_beta_AtSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_AtSC, method="bray", permutations=999)
  AtSC_KO_df <- cbind(KO_adonis_AtSC,Variable = rownames(KO_adonis_AtSC), Subset = "AtSC", Test = "Functions", Rank="Functions")
  
  metadata_HvSC <- metadata[metadata$Inoculum == "HvSC",]
  metadata_HvSC <- metadata_HvSC[metadata_HvSC$Condition != Plant_list[j],]
  KO_beta_HvSC <- beta_isolate_KO[row.names(beta_isolate_KO) %in% row.names(metadata_HvSC),row.names(beta_isolate_KO) %in% row.names(metadata_HvSC)]
  KO_adonis_HvSC <- adonis2(KO_beta_HvSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_HvSC, method="bray", permutations=999)
  HvSC_KO_df <- cbind(KO_adonis_HvSC,Variable = rownames(KO_adonis_HvSC), Subset = "HvSC", Test = "Functions", Rank="Functions")
  
  metadata_LjSC <- metadata[metadata$Inoculum == "LjSC",]
  metadata_LjSC <- metadata_LjSC[metadata_LjSC$Condition != Plant_list[j],]
  KO_beta_LjSC <- beta_isolate_KO[row.names(beta_isolate_KO) %in% row.names(metadata_LjSC),row.names(beta_isolate_KO) %in% row.names(metadata_LjSC)]
  KO_adonis_LjSC <- adonis2(KO_beta_LjSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_LjSC, method="bray", permutations=999)
  LjSC_KO_df <- cbind(KO_adonis_LjSC,Variable = rownames(KO_adonis_LjSC), Subset = "LjSC", Test = "Functions", Rank="Functions")
  
  metadata_SSC <- metadata[metadata$Inoculum == "SSC",]
  metadata_SSC <- metadata_SSC[metadata_SSC$Condition != Plant_list[j],]
  KO_beta_SSC <- beta_isolate_KO[row.names(beta_isolate_KO) %in% row.names(metadata_SSC),row.names(beta_isolate_KO) %in% row.names(metadata_SSC)]
  KO_adonis_SSC <- adonis2(KO_beta_SSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_SSC, method="bray", permutations=999)
  SSC_KO_df <- cbind(KO_adonis_SSC,Variable = rownames(KO_adonis_SSC), Subset = "SSC", Test = "Functions", Rank="Functions")
  
  #Per SynCom
  metadata_AtSC <- metadata[metadata$Inoculum == "AtSC",]
  metadata_AtSC <- metadata_AtSC[metadata_AtSC$Condition != Plant_list[j],]
  KO_beta_AtSC <- beta_isolate_KO_dom[row.names(beta_isolate_KO_dom) %in% row.names(metadata_AtSC),row.names(beta_isolate_KO_dom) %in% row.names(metadata_AtSC)]
  KO_adonis_AtSC <- adonis2(KO_beta_AtSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_AtSC, method="bray", permutations=999)
  AtSC_KO_df_dom <- cbind(KO_adonis_AtSC,Variable = rownames(KO_adonis_AtSC), Subset = "AtSC", Test = "Functions", Rank="Functions")
  
  metadata_HvSC <- metadata[metadata$Inoculum == "HvSC",]
  metadata_HvSC <- metadata_HvSC[metadata_HvSC$Condition != Plant_list[j],]
  KO_beta_HvSC <- beta_isolate_KO_dom[row.names(beta_isolate_KO_dom) %in% row.names(metadata_HvSC),row.names(beta_isolate_KO_dom) %in% row.names(metadata_HvSC)]
  KO_adonis_HvSC <- adonis2(KO_beta_HvSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_HvSC, method="bray", permutations=999)
  HvSC_KO_df_dom <- cbind(KO_adonis_HvSC,Variable = rownames(KO_adonis_HvSC), Subset = "HvSC", Test = "Functions", Rank="Functions")
  
  metadata_LjSC <- metadata[metadata$Inoculum == "LjSC",]
  metadata_LjSC <- metadata_LjSC[metadata_LjSC$Condition != Plant_list[j],]
  KO_beta_LjSC <- beta_isolate_KO_dom[row.names(beta_isolate_KO_dom) %in% row.names(metadata_LjSC),row.names(beta_isolate_KO_dom) %in% row.names(metadata_LjSC)]
  KO_adonis_LjSC <- adonis2(KO_beta_LjSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_LjSC, method="bray", permutations=999)
  LjSC_KO_df_dom <- cbind(KO_adonis_LjSC,Variable = rownames(KO_adonis_LjSC), Subset = "LjSC", Test = "Functions", Rank="Functions")
  
  metadata_SSC <- metadata[metadata$Inoculum == "SSC",]
  metadata_SSC <- metadata_SSC[metadata_SSC$Condition != Plant_list[j],]
  KO_beta_SSC <- beta_isolate_KO_dom[row.names(beta_isolate_KO_dom) %in% row.names(metadata_SSC),row.names(beta_isolate_KO_dom) %in% row.names(metadata_SSC)]
  KO_adonis_SSC <- adonis2(KO_beta_SSC ~ Inoculum*Condition*Nutrient*Experiment, data=metadata_SSC, method="bray", permutations=999)
  SSC_KO_df_dom <- cbind(KO_adonis_SSC,Variable = rownames(KO_adonis_SSC), Subset = "SSC", Test = "Functions", Rank="Functions")
  
  combined_df_syncom <- do.call(rbind, comb_syncom_list)
  combined_df_syncom <- rbind(combined_df_syncom,AtSC_KO_df,HvSC_KO_df,LjSC_KO_df,SSC_KO_df)
  
  combined_df_syncom_dom <- do.call(rbind, comb_syncom_list_dom)
  combined_df_syncom_dom <- rbind(combined_df_syncom_dom,AtSC_KO_df_dom,HvSC_KO_df_dom,LjSC_KO_df_dom,SSC_KO_df_dom)
  
  combined_df_syncom <- combined_df_syncom[!(combined_df_syncom$Variable %in% c("Total", "Residual", "Condition:Nutrient", "Condition:Experiment","Experiment", "Nutrient")), ]
  combined_df_syncom$Rank=factor(combined_df_syncom$Rank, levels = c("Isolate", "Genus", "Family", "Order", "Class","Functions"))
  combined_df_syncom$Test=factor(combined_df_syncom$Test, levels = c("Taxonomy", "Functions"))
  combined_df_syncom$Dominance <- "Yes"
  
  combined_df_syncom_dom <- combined_df_syncom_dom[!(combined_df_syncom_dom$Variable %in% c("Total", "Residual", "Condition:Nutrient", "Condition:Experiment","Experiment", "Nutrient")), ]
  combined_df_syncom_dom$Rank=factor(combined_df_syncom_dom$Rank, levels = c("Isolate", "Genus", "Family", "Order", "Class","Functions"))
  combined_df_syncom_dom$Test=factor(combined_df_syncom_dom$Test, levels = c("Taxonomy", "Functions"))
  combined_df_syncom_dom$Dominance <- "No"
  
  combined_df_syncom <- rbind(combined_df_syncom,combined_df_syncom_dom)
  
  combined_df_syncom_2=subset.data.frame(x = combined_df_syncom, subset = combined_df_syncom$Rank!="Functions")
  combined_df_syncom_3=subset.data.frame(x = combined_df_syncom, subset = combined_df_syncom$Rank=="Functions")
  
  #For heatmap
  combined_df_syncom_3$Drop_out <- Plant_list[j]
  row.names(combined_df_syncom_3) <- NULL
  combined_df_syncom_4_dom <- rbind(combined_df_syncom_4_dom, combined_df_syncom_3)
  
  #For Table S4
  combined_df_syncom$Drop_out <- Plant_list[j]
  row.names(combined_df_syncom) <- NULL
  combined_df_syncom_5_dom <- rbind(combined_df_syncom_5_dom, combined_df_syncom)
  
}

write.table(combined_df_syncom_5_dom, paste(results.dir, "Table_S5_R2_no_dominators.tsv"), col.names =T, row.names =F, sep = "\t", quote =F)

#PieDonut plot
source(paste(working_directory, "PieDonutCustom_SSC_FL.R", sep = ""))

# SynCom colors
SynCom_colors <- data.frame(c("AtSC", "HvSC", "LjSC", "SSC"),c("#A3A500","#00B0F6","#00BF7D","#F8766D"))
colnames(SynCom_colors) <- c("SynCom", "Colour")

SynComs <- c("AtSC", "LjSC", "HvSC", "SSC")

hop_comp <- data.frame(matrix(NA, ncol =4))
colnames(hop_comp) <- c("SynCom", "Arabidopsis", "Barley", "Lotus")
hop_comp <- hop_comp[-1,]

for (syncom in SynComs){
  combined_df_syncom_4_dom_sub <- unique(combined_df_syncom_4_dom[combined_df_syncom_4_dom$Subset == paste(syncom),])
  combined_df_syncom_4_dom_sub_2 <- combined_df_syncom_4_dom_sub[combined_df_syncom_4_dom_sub$Dominance == "Yes",]
  
  At <- abs(combined_df_syncom_4_dom_sub_2$R2[combined_df_syncom_4_dom_sub_2$Drop_out == "None"] - combined_df_syncom_4_dom_sub_2$R2[combined_df_syncom_4_dom_sub_2$Drop_out == "Arabidopsis"])
  Hv <- abs(combined_df_syncom_4_dom_sub_2$R2[combined_df_syncom_4_dom_sub_2$Drop_out == "None"] - combined_df_syncom_4_dom_sub_2$R2[combined_df_syncom_4_dom_sub_2$Drop_out == "Barley"])
  Lj <- abs(combined_df_syncom_4_dom_sub_2$R2[combined_df_syncom_4_dom_sub_2$Drop_out == "None"] - combined_df_syncom_4_dom_sub_2$R2[combined_df_syncom_4_dom_sub_2$Drop_out == "Lotus"])
  
  hop_sub <- t(data.frame(c(paste(syncom), At, Hv, Lj)))
  row.names(hop_sub) <- NULL
  hop_comp <- rbind(hop_comp, hop_sub)
}

colnames(hop_comp) <- c("SynCom", "Arabidopsis", "Barley", "Lotus")

hop_comp$Arabidopsis <- as.numeric(hop_comp$Arabidopsis)
hop_comp$Barley <- as.numeric(hop_comp$Barley)
hop_comp$Lotus <- as.numeric(hop_comp$Lotus)

hop_comp_2 <- melt(hop_comp)

hop_comp_2$Rel <- round(hop_comp_2$value/sum(hop_comp_2$value)*10000,0)

hop_comp_2$Combination <- paste(hop_comp_2$SynCom, hop_comp_2$variable,sep = "_")

pie_data_2 <- data.frame()

for (combi in hop_comp_2$Combination){
  new <- hop_comp_2$Rel[hop_comp_2$Combination == paste(combi)]
  syncom <- hop_comp_2$SynCom[hop_comp_2$Combination == paste(combi)]
  plant <- hop_comp_2$variable[hop_comp_2$Combination == paste(combi)]
  
  for (i in 1:new){
    pie_data <- data.frame(paste(syncom), paste(plant))
    pie_data_2 <- rbind(pie_data_2, pie_data)
  }
}

colnames(pie_data_2) <- c("Inoculum", "Host")

print(PieDonutCustom_SSC(pie_data_2,aes(pies=Inoculum,donuts=Host),showRatioThreshold = 0.001, r0 = getOption("PieDonut.r0", 0.3), r1 = getOption("PieDonut.r1", 0.7), r2 = getOption("PieDonut.r2", 1.1), color = "gray50"))

pdf(paste(results.dir,"Figure_3a_PieDonut_with_dom.pdf", sep=""), width=6, height=6)
print(PieDonutCustom_SSC(pie_data_2,aes(pies=SynCom,donuts=Host),showRatioThreshold = 0.001, r0 = getOption("PieDonut.r0", 0.3), r1 = getOption("PieDonut.r1", 0.7), r2 = getOption("PieDonut.r2", 1.1), color = "gray50"))
dev.off()

#Figure 3b - piedonut without dominances
SynComs <- c("AtSC", "LjSC", "HvSC", "SSC")

hop_comp_B <- data.frame(matrix(NA, ncol =4))
colnames(hop_comp_B) <- c("SynCom", "Arabidopsis", "Barley", "Lotus")
hop_comp_B <- hop_comp_B[-1,]

for (syncom in SynComs){
  combined_df_syncom_4_dom_sub <- combined_df_syncom_4_dom[combined_df_syncom_4_dom$Subset == paste(syncom),]
  combined_df_syncom_4_dom_sub_2 <- combined_df_syncom_4_dom_sub[combined_df_syncom_4_dom_sub$Dominance == "No",]
  
  At <- abs(combined_df_syncom_4_dom_sub_2$R2[combined_df_syncom_4_dom_sub_2$Drop_out == "None"] - combined_df_syncom_4_dom_sub_2$R2[combined_df_syncom_4_dom_sub_2$Drop_out == "Arabidopsis"])
  Hv <- abs(combined_df_syncom_4_dom_sub_2$R2[combined_df_syncom_4_dom_sub_2$Drop_out == "None"] - combined_df_syncom_4_dom_sub_2$R2[combined_df_syncom_4_dom_sub_2$Drop_out == "Barley"])
  Lj <- abs(combined_df_syncom_4_dom_sub_2$R2[combined_df_syncom_4_dom_sub_2$Drop_out == "None"] - combined_df_syncom_4_dom_sub_2$R2[combined_df_syncom_4_dom_sub_2$Drop_out == "Lotus"])
  
  hop_sub_B <- t(data.frame(c(paste(syncom), At, Hv, Lj)))
  row.names(hop_sub_B) <- NULL
  hop_comp_B <- rbind(hop_comp_B, hop_sub_B)
}

colnames(hop_comp_B) <- c("SynCom", "Arabidopsis", "Barley", "Lotus")

hop_comp_B$Arabidopsis <- as.numeric(hop_comp_B$Arabidopsis)
hop_comp_B$Barley <- as.numeric(hop_comp_B$Barley)
hop_comp_B$Lotus <- as.numeric(hop_comp_B$Lotus)

hop_comp_2_B <- melt(hop_comp_B)

hop_comp_2_B$Rel <- round(hop_comp_2_B$value/sum(hop_comp_2_B$value)*10000,0)

hop_comp_2_B$Combination <- paste(hop_comp_2_B$SynCom, hop_comp_2_B$variable,sep = "_")

pie_data_2_B <- data.frame()

for (combi in hop_comp_2_B$Combination){
  new <- hop_comp_2_B$Rel[hop_comp_2_B$Combination == paste(combi)]
  syncom <- hop_comp_2_B$SynCom[hop_comp_2_B$Combination == paste(combi)]
  plant <- hop_comp_2$variable[hop_comp_2_B$Combination == paste(combi)]
  
  for (i in 1:new){
    pie_data_B <- data.frame(paste(syncom), paste(plant))
    pie_data_2_B <- rbind(pie_data_2_B, pie_data_B)
  }
}

colnames(pie_data_2_B) <- c("Inoculum", "Host")

print(PieDonutCustom_SSC(pie_data_2_B,aes(pies=Inoculum,donuts=Host),showRatioThreshold = 0.001, r0 = getOption("PieDonut.r0", 0.3), r1 = getOption("PieDonut.r1", 0.7), r2 = getOption("PieDonut.r2", 1.1), color = "gray50"))

pdf(paste(results.dir,"Figure_3b_PieDonut_no_dom.pdf", sep=""), width=6, height=6)
print(PieDonutCustom_SSC(pie_data_2_B,aes(pies=Inoculum,donuts=Host),showRatioThreshold = 0.001, r0 = getOption("PieDonut.r0", 0.3), r1 = getOption("PieDonut.r1", 0.7), r2 = getOption("PieDonut.r2", 1.1), color = "gray50"))
dev.off()


###Figure 3c - R2 Line diagram - with and without Dominances =====
#otu table
norm_SSC =read.table(paste(working_directory,"Isolate_tables/Original/SSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
#Taxonomy table
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")
#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
colnames(samples_df_2)[5]="Nutrient"
samples_df_2$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df_2$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

sapply(tax_df, function(x) length(unique(x)))

#Set the OTU, TAX and sample data for making phyloseq object
OTU = otu_table(as.matrix(norm_SSC),taxa_are_rows = TRUE)
#TAX = tax_table(tax_mat)
TAX = tax_table(as.matrix(tax_df_2))

#Sample subsetting

cond="ES"
samples_df_sub <- subset(samples_df_2, samples_df_2$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

samples_sub = sample_data(samples_df_sub_2)

phylo_sub = phyloseq(OTU,TAX, samples_sub)

phylo_sub_RA=microbiome::transform(x = phylo_sub, transform = "compositional" )

#Agglomerate to phylum-level and rename
ps_phylum <- phyloseq::tax_glom(phylo_sub_RA, "Phylum")
ps_class <- phyloseq::tax_glom(phylo_sub_RA, "Class")
ps_order <- phyloseq::tax_glom(phylo_sub_RA, "Order")
ps_family <- phyloseq::tax_glom(phylo_sub_RA, "Family")
ps_genus <- phyloseq::tax_glom(phylo_sub_RA, "Genus")

phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]
phyloseq::taxa_names(ps_class) <- phyloseq::tax_table(ps_class)[, "Class"]
phyloseq::taxa_names(ps_order) <- phyloseq::tax_table(ps_order)[, "Order"]
phyloseq::taxa_names(ps_family) <- phyloseq::tax_table(ps_family)[, "Family"]
phyloseq::taxa_names(ps_genus) <- phyloseq::tax_table(ps_genus)[, "Genus"]

#Bray Curtis distance matrix
beta_phylum <- as.matrix(vegdist(t(ps_phylum@otu_table@.Data), method = "bray", diag = T))
beta_class <- as.matrix(vegdist(t(ps_class@otu_table@.Data), method = "bray", diag = T))
beta_order <- as.matrix(vegdist(t(ps_order@otu_table@.Data), method = "bray", diag = T))
beta_family <- as.matrix(vegdist(t(ps_family@otu_table@.Data), method = "bray", diag = T))
beta_genus <- as.matrix(vegdist(t(ps_genus@otu_table@.Data), method = "bray", diag = T))
beta_isolate <- as.matrix(vegdist(t(phylo_sub_RA@otu_table@.Data), method = "bray", diag = T))

mean_value_phylum= mean(beta_phylum)
mean_value_class= mean(beta_class)
mean_value_order= mean(beta_order)
mean_value_family= mean(beta_family)
mean_value_genus= mean(beta_genus)
mean_value_isolate= mean(beta_isolate)

Taxo_order=c("Isolate","Genus", "Family", "Order", "Class")
beta_distance_list=list(beta_isolate, beta_genus, beta_family, beta_order, beta_phylum)
Adonis_list=list()
Plot_list=list()
Plot_list_2=list()

for (i in 1:length(beta_distance_list)) {
  
  Bray_curtis_df=beta_distance_list[[i]]
  
  #Make PCoA plot for Bray Curtis Distance matrix
  pcoa = cmdscale(Bray_curtis_df, k=3, eig=T)
  points = as.data.frame(pcoa$points)
  colnames(points) = c("x", "y", "z") 
  eig = pcoa$eig
  
  points = merge(points,samples_df_sub_2, by = "row.names")
  rownames(points) <- points$Row.names
  points <- points %>% dplyr::select (-Row.names)
  
  points$Condition <- factor(points$Condition, levels = c("At","Hv", "Lj"))
  points$Inoculum <- factor(points$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC","NS"))
  points$Nutrient <- factor(points$Nutrient, levels = c("low", "high"))
  points$Experiment  <- factor(points$Experiment, levels = c("R1", "R2"))
  
  metadata=points[,-c(1,2,3)]
  
  #  Run adonis PERMANOVA test
  set.seed(1)
  SSC_bray_adonis <- adonis2(beta_distance_list[[i]] ~ Inoculum*Condition, data=metadata, method="bray", permutations=999)
  
  Adonis_list[[i]]=SSC_bray_adonis
}

#otu table
KO_SSC =read.table(paste(working_directory,"KO_tables/Original/SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object
OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)

#Sample subsetting
cond="ES"

samples_df_sub <- subset(samples_df_2, samples_df_2$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

samples_sub = sample_data(samples_df_sub_2)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)

phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )

beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

Bray_curtis_df=beta_isolate_KO

#Make PCoA plot for Bray Curtis Distance matrix
pcoa = cmdscale(Bray_curtis_df, k=3, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z") 
eig = pcoa$eig

points = merge(points,samples_df_sub_2, by = "row.names")
rownames(points) <- points$Row.names
points <- points %>% dplyr::select (-Row.names)

points$Condition <- factor(points$Condition, levels = c("At","Hv", "Lj"))
points$Inoculum <- factor(points$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC","NS"))
points$Nutrient <- factor(points$Nutrient, levels = c("low", "high"))
points$Experiment  <- factor(points$Experiment, levels = c("R1", "R2"))

metadata=points[,-c(1,2,3)]

set.seed(1)
SSC_bray_KO_adonis <- adonis2(beta_isolate_KO ~ Inoculum*Condition, data=metadata, method="bray", permutations=999)

#You need to run the script of figure 5 before running this one
# Adonis results from KO functions
Adonis_list[[length(beta_distance_list)+1]]=SSC_bray_KO_adonis

# Make a dataframe out of all these results with R2 values
R2_values <- matrix(nrow = length(Adonis_list), ncol = 2)
for (i in 1:length(Adonis_list)) {
  R2_values[i, ] <- Adonis_list[[i]]$R2[1:2]
}

# Row and column names
Taxo_order <- c("Isolate", "Genus", "Family", "Order", "Phylum","KO_functions")
col_names <- c("Syncom", "Plant")

# Create the dataframe with R2 values
df <- data.frame(Taxo_order, R2_values)
colnames(df)[-1] <- col_names
df

# Melt the dataframe
melted_df <- melt(df, id.vars = "Taxo_order")
melted_df <- melted_df %>%
  mutate(Test = ifelse(Taxo_order == "KO_functions", "Function", "Taxonomy"))

colnames(melted_df) <- c("Taxo_order", "Variable", "R2", "Test")
melted_df$Variable <- as.character(melted_df$Variable)

melted_df$Variable[melted_df$Variable == "Syncom"] <- "Inoculum"
melted_df$Variable[melted_df$Variable == "Plant"] <- "Host"

melted_df$Taxo_order=factor(melted_df$Taxo_order, levels = c("Isolate", "Genus", "Family", "Order", "Class", "Phylum","KO_functions"))
melted_df$Variable=factor(melted_df$Variable, levels = c("Inoculum", "Host"))
melted_df$Test=factor(melted_df$Test, levels = c("Taxonomy", "Function"))

# Create extra dfs for the line plot
# Taxonomy
melted_df_2=subset.data.frame(x = melted_df, subset = melted_df$Taxo_order!="KO_functions")
# Functions for dashed line
melted_df_3=subset.data.frame(x = melted_df, subset = melted_df$Taxo_order=="KO_functions")

melted_df_2$Dominance <- "Yes"
melted_df_3$Dominance <- "Yes"

#Data without dominances
#otu table
norm_SSC =read.table(paste(working_directory,"Isolate_tables/No_dominances/SSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
#Taxonomy table
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")
#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
colnames(samples_df_2)[5]="Nutrient"
samples_df_2$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df_2$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

sapply(tax_df, function(x) length(unique(x)))

#Set the OTU, TAX and sample data for making phyloseq object
OTU = otu_table(as.matrix(norm_SSC),taxa_are_rows = TRUE)
#TAX = tax_table(tax_mat)
TAX = tax_table(as.matrix(tax_df_2))

#Sample subsetting

cond="ES"
samples_df_sub <- subset(samples_df_2, samples_df_2$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

samples_sub = sample_data(samples_df_sub_2)

phylo_sub = phyloseq(OTU,TAX, samples_sub)

phylo_sub_RA=microbiome::transform(x = phylo_sub, transform = "compositional" )

#Agglomerate to phylum-level and rename
ps_phylum <- phyloseq::tax_glom(phylo_sub_RA, "Phylum")
ps_class <- phyloseq::tax_glom(phylo_sub_RA, "Class")
ps_order <- phyloseq::tax_glom(phylo_sub_RA, "Order")
ps_family <- phyloseq::tax_glom(phylo_sub_RA, "Family")
ps_genus <- phyloseq::tax_glom(phylo_sub_RA, "Genus")

phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]
phyloseq::taxa_names(ps_class) <- phyloseq::tax_table(ps_class)[, "Class"]
phyloseq::taxa_names(ps_order) <- phyloseq::tax_table(ps_order)[, "Order"]
phyloseq::taxa_names(ps_family) <- phyloseq::tax_table(ps_family)[, "Family"]
phyloseq::taxa_names(ps_genus) <- phyloseq::tax_table(ps_genus)[, "Genus"]

#Bray Curtis distance matrix
beta_phylum_dom <- as.matrix(vegdist(t(ps_phylum@otu_table@.Data), method = "bray", diag = T))
beta_class_dom <- as.matrix(vegdist(t(ps_class@otu_table@.Data), method = "bray", diag = T))
beta_order_dom <- as.matrix(vegdist(t(ps_order@otu_table@.Data), method = "bray", diag = T))
beta_family_dom <- as.matrix(vegdist(t(ps_family@otu_table@.Data), method = "bray", diag = T))
beta_genus_dom <- as.matrix(vegdist(t(ps_genus@otu_table@.Data), method = "bray", diag = T))
beta_isolate_dom <- as.matrix(vegdist(t(phylo_sub_RA@otu_table@.Data), method = "bray", diag = T))

mean_value_phylum_dom = mean(beta_phylum_dom)
mean_value_class_dom = mean(beta_class_dom)
mean_value_order_dom = mean(beta_order_dom)
mean_value_family_dom = mean(beta_family_dom)
mean_value_genus_dom = mean(beta_genus_dom)
mean_value_isolate_dom = mean(beta_isolate_dom)

Taxo_order=c("Isolate","Genus", "Family", "Order", "Class")
beta_distance_list_dom=list(beta_isolate_dom, beta_genus_dom, beta_family_dom, beta_order_dom, beta_phylum_dom)
Adonis_list_dom=list()
Plot_list_dom=list()
Plot_list_2_dom=list()

for (i in 1:length(beta_distance_list_dom)) {
  
  Bray_curtis_df=beta_distance_list_dom[[i]]
  
  #Make PCoA plot for Bray Curtis Distance matrix
  pcoa = cmdscale(Bray_curtis_df, k=3, eig=T)
  points = as.data.frame(pcoa$points)
  colnames(points) = c("x", "y", "z") 
  eig = pcoa$eig
  
  points = merge(points,samples_df_sub_2, by = "row.names")
  rownames(points) <- points$Row.names
  points <- points %>% dplyr::select (-Row.names)
  
  points$Condition <- factor(points$Condition, levels = c("At","Hv", "Lj"))
  points$Inoculum <- factor(points$Inoculum, levels = c("SSC","AtSC", "HvSC", "LjSC","NS"))
  points$Nutrient <- factor(points$Nutrient, levels = c("low", "high"))
  points$Experiment  <- factor(points$Experiment, levels = c("R1", "R2"))
  
  metadata=points[,-c(1,2,3)]
  
  #  Run adonis PERMANOVA test
  set.seed(1)
  SSC_bray_adonis <- adonis2(beta_distance_list_dom[[i]] ~ Inoculum*Condition, data=metadata, method="bray", permutations=999)
  
  Adonis_list_dom[[i]]=SSC_bray_adonis
}

#otu table
KO_SSC =read.table(paste(working_directory,"KO_tables/No_dominances/SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object
OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)

#Sample subsetting
cond="ES"

samples_df_sub <- subset(samples_df_2, samples_df_2$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

samples_sub = sample_data(samples_df_sub_2)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)

phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )

beta_isolate_KO_dom <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

Bray_curtis_df=beta_isolate_KO_dom

#Make PCoA plot for Bray Curtis Distance matrix
pcoa = cmdscale(Bray_curtis_df, k=3, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z") 
eig = pcoa$eig

points = merge(points,samples_df_sub_2, by = "row.names")
rownames(points) <- points$Row.names
points <- points %>% dplyr::select (-Row.names)

points$Condition <- factor(points$Condition, levels = c("At","Hv", "Lj"))
points$Inoculum <- factor(points$Inoculum, levels = c("SSC","AtSC", "HvSC", "LjSC","NS"))
points$Nutrient <- factor(points$Nutrient, levels = c("low", "high"))
points$Experiment  <- factor(points$Experiment, levels = c("R1", "R2"))

metadata=points[,-c(1,2,3)]

set.seed(1)
SSC_bray_KO_adonis <- adonis2(beta_isolate_KO_dom ~ Inoculum*Condition, data=metadata, method="bray", permutations=999)

# Adonis results from KO functions
Adonis_list_dom[[length(beta_distance_list_dom)+1]]=SSC_bray_KO_adonis

# Make a dataframe out of all these results with R2 values
R2_values_dom <- matrix(nrow = length(Adonis_list_dom), ncol = 2)
for (i in 1:length(Adonis_list_dom)) {
  R2_values_dom[i, ] <- Adonis_list_dom[[i]]$R2[1:2]
}

# Row and column names
Taxo_order <- c("Isolate", "Genus", "Family", "Order", "Phylum","KO_functions")
col_names <- c("Syncom", "Plant")

# Create the dataframe with R2 values
df <- data.frame(Taxo_order, R2_values_dom)
colnames(df)[-1] <- col_names
df

# Melt the dataframe
melted_df <- melt(df, id.vars = "Taxo_order")
melted_df <- melted_df %>%
  mutate(Test = ifelse(Taxo_order == "KO_functions", "Function", "Taxonomy"))

colnames(melted_df) <- c("Taxo_order", "Variable", "R2", "Test")
melted_df$Variable <- as.character(melted_df$Variable)

melted_df$Variable[melted_df$Variable == "Syncom"] <- "Inoculum"
melted_df$Variable[melted_df$Variable == "Plant"] <- "Host"

melted_df$Taxo_order=factor(melted_df$Taxo_order, levels = c("Isolate", "Genus", "Family", "Order", "Class", "Phylum","KO_functions"))
melted_df$Variable=factor(melted_df$Variable, levels = c("Inoculum", "Host"))
melted_df$Test=factor(melted_df$Test, levels = c("Taxonomy", "Function"))

# Create extra dfs for the line plot
# Taxonomy
melted_df_dom_2=subset.data.frame(x = melted_df, subset = melted_df$Taxo_order!="KO_functions")
# Functions for dashed line
melted_df_dom_3=subset.data.frame(x = melted_df, subset = melted_df$Taxo_order=="KO_functions")

melted_df_dom_2$Dominance <- "No"
melted_df_dom_3$Dominance <- "No"

melted_df_2 <- rbind(melted_df_2, melted_df_dom_2)
melted_df_3 <- rbind(melted_df_3, melted_df_dom_3)

colnames(melted_df_2) <- c("Taxo_order", "Variable", "R2", "Test", "Dominance")
colnames(melted_df_3) <- c("Taxo_order", "Variable", "R2", "Test", "Dominance")

melted_df_2_2 <- melted_df_2[melted_df_2$Variable != "Experiment",]
melted_df_3_2 <- melted_df_3[melted_df_3$Variable != "Experiment",]

line_plot1=ggplot(melted_df_2_2, aes(x = Taxo_order, y = R2, color = Dominance)) +
  theme_bw() +
  geom_point(size = 5) +
  geom_line(aes(group = Dominance), size = 2) + # Ensure grouping is correct
  ylim(0,0.6)+
  scale_color_manual(values = c("gray70","black"))+
  labs(x = "Taxonomic order", y = "R2 values") +
  geom_line(aes(group = Taxo_order), color = "gray30", alpha = 0.5, linetype = "dashed") + # Distance between dominance and non dominance conditions
  geom_hline(data = melted_df_3_2, aes(yintercept = R2, linetype = "R2 KO_functions", color=Dominance), size = 1) + # Functional R2 horizontal dashed line
  scale_linetype_manual(values = c("dashed", "dashed", "dashed", "dashed", "dashed", "dashed")) +
  facet_wrap(~ Variable, nrow=1)+theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "right")+
  theme(
    strip.background=element_rect(colour="gray50", fill = "transparent", size=0.3), # Change 'size' for thickness
    axis.text=element_text(color="gray50"),
    axis.line = element_line(color="gray50", size=0.3)
  ) +
  
  guides(linetype = guide_legend(title = "Dashed line",nrow = 2 ))

line_plot1

pdf(paste(results.dir,"Figure_3c_R2_line_plot.pdf", sep=""), width=12, height=8)
print(line_plot1)
dev.off()

###Table S4 - R2 values at different taxonomic levels (drop out R2 values are associated with 3c) =====
#otu table
norm_SSC =read.table(paste(working_directory,"Isolate_tables/Original/SSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
#Taxonomy table
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")
#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
colnames(samples_df_2)[5]="Nutrient"
samples_df_2$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df_2$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

sapply(tax_df, function(x) length(unique(x)))

#Set the OTU, TAX and sample data for making phyloseq object
OTU = otu_table(as.matrix(norm_SSC),taxa_are_rows = TRUE)
#TAX = tax_table(tax_mat)
TAX = tax_table(as.matrix(tax_df_2))

#Sample subsetting

cond="ES"
samples_df_sub <- subset(samples_df_2, samples_df_2$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

samples_sub = sample_data(samples_df_sub_2)

phylo_sub = phyloseq(OTU,TAX, samples_sub)

phylo_sub_RA=microbiome::transform(x = phylo_sub, transform = "compositional" )

#Agglomerate to phylum-level and rename
ps_phylum <- phyloseq::tax_glom(phylo_sub_RA, "Phylum")
ps_class <- phyloseq::tax_glom(phylo_sub_RA, "Class")
ps_order <- phyloseq::tax_glom(phylo_sub_RA, "Order")
ps_family <- phyloseq::tax_glom(phylo_sub_RA, "Family")
ps_genus <- phyloseq::tax_glom(phylo_sub_RA, "Genus")

phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]
phyloseq::taxa_names(ps_class) <- phyloseq::tax_table(ps_class)[, "Class"]
phyloseq::taxa_names(ps_order) <- phyloseq::tax_table(ps_order)[, "Order"]
phyloseq::taxa_names(ps_family) <- phyloseq::tax_table(ps_family)[, "Family"]
phyloseq::taxa_names(ps_genus) <- phyloseq::tax_table(ps_genus)[, "Genus"]

#Bray Curtis distance matrix
beta_phylum <- as.matrix(vegdist(t(ps_phylum@otu_table@.Data), method = "bray", diag = T))
beta_class <- as.matrix(vegdist(t(ps_class@otu_table@.Data), method = "bray", diag = T))
beta_order <- as.matrix(vegdist(t(ps_order@otu_table@.Data), method = "bray", diag = T))
beta_family <- as.matrix(vegdist(t(ps_family@otu_table@.Data), method = "bray", diag = T))
beta_genus <- as.matrix(vegdist(t(ps_genus@otu_table@.Data), method = "bray", diag = T))
beta_isolate <- as.matrix(vegdist(t(phylo_sub_RA@otu_table@.Data), method = "bray", diag = T))

mean_value_phylum= mean(beta_phylum)
mean_value_class= mean(beta_class)
mean_value_order= mean(beta_order)
mean_value_family= mean(beta_family)
mean_value_genus= mean(beta_genus)
mean_value_isolate= mean(beta_isolate)

Taxo_order=c("Isolate","Genus", "Family", "Order", "Class")
beta_distance_list=list(beta_isolate, beta_genus, beta_family, beta_order, beta_phylum)
Adonis_list=list()
Plot_list=list()
Plot_list_2=list()

for (i in 1:length(beta_distance_list)) {
  
  Bray_curtis_df=beta_distance_list[[i]]
  
  #Make PCoA plot for Bray Curtis Distance matrix
  pcoa = cmdscale(Bray_curtis_df, k=3, eig=T)
  points = as.data.frame(pcoa$points)
  colnames(points) = c("x", "y", "z") 
  eig = pcoa$eig
  
  points = merge(points,samples_df_sub_2, by = "row.names")
  rownames(points) <- points$Row.names
  points <- points %>% dplyr::select (-Row.names)
  
  points$Condition <- factor(points$Condition, levels = c("At","Hv", "Lj"))
  points$Inoculum <- factor(points$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC","NS"))
  points$Nutrient <- factor(points$Nutrient, levels = c("low", "high"))
  points$Experiment  <- factor(points$Experiment, levels = c("R1", "R2"))
  
  metadata=points[,-c(1,2,3)]
  
  #  Run adonis PERMANOVA test
  set.seed(1)
  SSC_bray_adonis <- adonis2(beta_distance_list[[i]] ~ Inoculum*Condition*Nutrient*Experiment, data=metadata, method="bray", permutations=999)
  
  Adonis_list[[i]]=SSC_bray_adonis
  
}

#otu table
KO_SSC =read.table(paste(working_directory,"KO_tables/Original/SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object
OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)

#Sample subsetting
cond="ES"

samples_df_sub <- subset(samples_df_2, samples_df_2$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

samples_sub = sample_data(samples_df_sub_2)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)

phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )

beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

Bray_curtis_df=beta_isolate_KO

#Make PCoA plot for Bray Curtis Distance matrix
pcoa = cmdscale(Bray_curtis_df, k=3, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z") 
eig = pcoa$eig

points = merge(points,samples_df_sub_2, by = "row.names")
rownames(points) <- points$Row.names
points <- points %>% dplyr::select (-Row.names)

points$Condition <- factor(points$Condition, levels = c("At","Hv", "Lj"))
points$Inoculum <- factor(points$Inoculum, levels = c("AtSC", "HvSC", "LjSC","SSC","NS"))
points$Nutrient <- factor(points$Nutrient, levels = c("low", "high"))
points$Experiment  <- factor(points$Experiment, levels = c("R1", "R2"))

metadata=points[,-c(1,2,3)]

set.seed(1)
SSC_bray_KO_adonis <- adonis2(beta_isolate_KO ~ Inoculum*Condition*Nutrient*Experiment, data=metadata, method="bray", permutations=999)

#You need to run the script of figure 5 before running this one
# Adonis results from KO functions
Adonis_list[[length(beta_distance_list)+1]]=SSC_bray_KO_adonis

# Make a dataframe out of all these results with R2 values
R2_values <- matrix(nrow = length(Adonis_list), ncol = 5)
for (i in 1:length(Adonis_list)) {
  R2_values[i, ] <- Adonis_list[[i]]$R2[1:5]
}

# Row and column names
Taxo_order <- c("Isolate", "Genus", "Family", "Order", "Phylum","KO_functions")
col_names <- c("Syncom", "Plant", "Nutrient", "Experiment", paste("Syncom:", "Plant", sep = ""))

# Create the dataframe with R2 values
df <- data.frame(Taxo_order, R2_values)
colnames(df)[-1] <- col_names
df

# Melt the dataframe
melted_df <- melt(df, id.vars = "Taxo_order")
melted_df <- melted_df %>%
  mutate(Test = ifelse(Taxo_order == "KO_functions", "Function", "Taxonomy"))

colnames(melted_df) <- c("Taxo_order", "Variable", "R2", "Test")
melted_df$Variable <- as.character(melted_df$Variable)

melted_df$Variable[melted_df$Variable == "Syncom"] <- "Inoculum"
melted_df$Variable[melted_df$Variable == "Plant"] <- "Host"
melted_df$Variable[melted_df$Variable == "Syncom:Plant"] <- "Inoculum:Host"

melted_df$Taxo_order=factor(melted_df$Taxo_order, levels = c("Isolate", "Genus", "Family", "Order", "Class", "Phylum","KO_functions"))
melted_df$Variable=factor(melted_df$Variable, levels = c("Inoculum", "Host","Inoculum:Host", "Nutrient", "Experiment"))
melted_df$Test=factor(melted_df$Test, levels = c("Taxonomy", "Function"))

# Create extra dfs for the line plot

# Taxonomy
melted_df_2=subset.data.frame(x = melted_df, subset = melted_df$Taxo_order!="KO_functions")
# Functions for dashed line
melted_df_3=subset.data.frame(x = melted_df, subset = melted_df$Taxo_order=="KO_functions")

melted_df_2$Dominance <- "Yes"
melted_df_3$Dominance <- "Yes"

#ADDING THE DATA WITHOUT THE NODULATORS & RHIZOBACTER
#otu table
norm_SSC =read.table(paste(working_directory,"Isolate_tables/No_dominances/SSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
#Taxonomy table
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")
#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
colnames(samples_df_2)[5]="Nutrient"
samples_df_2$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df_2$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

sapply(tax_df, function(x) length(unique(x)))

#Set the OTU, TAX and sample data for making phyloseq object
OTU = otu_table(as.matrix(norm_SSC),taxa_are_rows = TRUE)
#TAX = tax_table(tax_mat)
TAX = tax_table(as.matrix(tax_df_2))

#Sample subsetting

cond="ES"
samples_df_sub <- subset(samples_df_2, samples_df_2$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")
# row.names(samples_df_sub_2) <- sub('HL_orig', 'HL',row.names(samples_df_sub_2))

samples_sub = sample_data(samples_df_sub_2)

phylo_sub = phyloseq(OTU,TAX, samples_sub)

phylo_sub_RA=microbiome::transform(x = phylo_sub, transform = "compositional" )

#Agglomerate to phylum-level and rename
ps_phylum <- phyloseq::tax_glom(phylo_sub_RA, "Phylum")
ps_class <- phyloseq::tax_glom(phylo_sub_RA, "Class")
ps_order <- phyloseq::tax_glom(phylo_sub_RA, "Order")
ps_family <- phyloseq::tax_glom(phylo_sub_RA, "Family")
ps_genus <- phyloseq::tax_glom(phylo_sub_RA, "Genus")

phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]
phyloseq::taxa_names(ps_class) <- phyloseq::tax_table(ps_class)[, "Class"]
phyloseq::taxa_names(ps_order) <- phyloseq::tax_table(ps_order)[, "Order"]
phyloseq::taxa_names(ps_family) <- phyloseq::tax_table(ps_family)[, "Family"]
phyloseq::taxa_names(ps_genus) <- phyloseq::tax_table(ps_genus)[, "Genus"]

#Bray Curtis distance matrix
beta_phylum_dom <- as.matrix(vegdist(t(ps_phylum@otu_table@.Data), method = "bray", diag = T))
beta_class_dom <- as.matrix(vegdist(t(ps_class@otu_table@.Data), method = "bray", diag = T))
beta_order_dom <- as.matrix(vegdist(t(ps_order@otu_table@.Data), method = "bray", diag = T))
beta_family_dom <- as.matrix(vegdist(t(ps_family@otu_table@.Data), method = "bray", diag = T))
beta_genus_dom <- as.matrix(vegdist(t(ps_genus@otu_table@.Data), method = "bray", diag = T))
beta_isolate_dom <- as.matrix(vegdist(t(phylo_sub_RA@otu_table@.Data), method = "bray", diag = T))

mean_value_phylum_dom = mean(beta_phylum_dom)
mean_value_class_dom = mean(beta_class_dom)
mean_value_order_dom = mean(beta_order_dom)
mean_value_family_dom = mean(beta_family_dom)
mean_value_genus_dom = mean(beta_genus_dom)
mean_value_isolate_dom = mean(beta_isolate_dom)

Taxo_order=c("Isolate","Genus", "Family", "Order", "Class")
beta_distance_list_dom=list(beta_isolate_dom, beta_genus_dom, beta_family_dom, beta_order_dom, beta_phylum_dom)
Adonis_list_dom=list()
Plot_list_dom=list()
Plot_list_2_dom=list()

for (i in 1:length(beta_distance_list_dom)) {
  
  Bray_curtis_df=beta_distance_list_dom[[i]]
  
  #Make PCoA plot for Bray Curtis Distance matrix
  pcoa = cmdscale(Bray_curtis_df, k=3, eig=T)
  points = as.data.frame(pcoa$points)
  colnames(points) = c("x", "y", "z") 
  eig = pcoa$eig
  
  points = merge(points,samples_df_sub_2, by = "row.names")
  rownames(points) <- points$Row.names
  points <- points %>% dplyr::select (-Row.names)
  
  points$Condition <- factor(points$Condition, levels = c("At","Hv", "Lj"))
  points$Inoculum <- factor(points$Inoculum, levels = c("SSC","AtSC", "HvSC", "LjSC","NS"))
  points$Nutrient <- factor(points$Nutrient, levels = c("low", "high"))
  points$Experiment  <- factor(points$Experiment, levels = c("R1", "R2"))
  
  metadata=points[,-c(1,2,3)]
  
  #  Run adonis PERMANOVA test
  set.seed(1)
  SSC_bray_adonis <- adonis2(beta_distance_list_dom[[i]] ~ Inoculum*Condition*Nutrient*Experiment, data=metadata, method="bray", permutations=999)
  
  Adonis_list_dom[[i]]=SSC_bray_adonis
}

#otu table
KO_SSC =read.table(paste(working_directory,"KO_tables/No_dominances/SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object
OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)

#Sample subsetting
cond="ES"

samples_df_sub <- subset(samples_df_2, samples_df_2$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

samples_sub = sample_data(samples_df_sub_2)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)

phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )

beta_isolate_KO_dom <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

Bray_curtis_df=beta_isolate_KO_dom

#Make PCoA plot for Bray Curtis Distance matrix
pcoa = cmdscale(Bray_curtis_df, k=3, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z") 
eig = pcoa$eig

points = merge(points,samples_df_sub_2, by = "row.names")
rownames(points) <- points$Row.names
points <- points %>% dplyr::select (-Row.names)

points$Condition <- factor(points$Condition, levels = c("At","Hv", "Lj"))
points$Inoculum <- factor(points$Inoculum, levels = c("SSC","AtSC", "HvSC", "LjSC","NS"))
points$Nutrient <- factor(points$Nutrient, levels = c("low", "high"))
points$Experiment  <- factor(points$Experiment, levels = c("R1", "R2"))

metadata=points[,-c(1,2,3)]

set.seed(1)
SSC_bray_KO_adonis <- adonis2(beta_isolate_KO_dom ~ Inoculum*Condition*Nutrient*Experiment, data=metadata, method="bray", permutations=999)

# Adonis results from KO functions
Adonis_list_dom[[length(beta_distance_list_dom)+1]]=SSC_bray_KO_adonis

# Make a dataframe out of all these results with R2 values
R2_values_dom <- matrix(nrow = length(Adonis_list_dom), ncol = 5)
for (i in 1:length(Adonis_list_dom)) {
  R2_values_dom[i, ] <- Adonis_list_dom[[i]]$R2[1:5]
}

# Row and column names
Taxo_order <- c("Isolate", "Genus", "Family", "Order", "Phylum","KO_functions")
col_names <- c("Syncom", "Plant", "Nutrient", "Experiment", paste("Syncom:", "Plant", sep = ""))

# Create the dataframe with R2 values
df <- data.frame(Taxo_order, R2_values_dom)
colnames(df)[-1] <- col_names
df

# Melt the dataframe
melted_df <- melt(df, id.vars = "Taxo_order")
melted_df <- melted_df %>%
  mutate(Test = ifelse(Taxo_order == "KO_functions", "Function", "Taxonomy"))

colnames(melted_df) <- c("Taxo_order", "Variable", "R2", "Test")
melted_df$Variable <- as.character(melted_df$Variable)

melted_df$Variable[melted_df$Variable == "Syncom"] <- "Inoculum"
melted_df$Variable[melted_df$Variable == "Plant"] <- "Host"
melted_df$Variable[melted_df$Variable == "Syncom:Plant"] <- "Inoculum:Host"

melted_df$Taxo_order=factor(melted_df$Taxo_order, levels = c("Isolate", "Genus", "Family", "Order", "Class", "Phylum","KO_functions"))
melted_df$Variable=factor(melted_df$Variable, levels = c("Inoculum", "Host","Inoculum:Host", "Nutrient", "Experiment"))
melted_df$Test=factor(melted_df$Test, levels = c("Taxonomy", "Function"))

write.table(melted_df, paste(results.dir,"Table_S4_R2_values_all_SCs.txt", sep =""), sep = "\t", quote =F,col.names =T, row.names =F)


###Figure 3d - Intrafamily diversity heatmap =====
SynComs <- c("AtSC", "LjSC", "HvSC", "SSC")

fam_5 <- data.frame(matrix(NA, ncol =9))
colnames(fam_5) <- c("Isolate", "Rel", "Rel_prop", "KO", "KO_prop", "Family", "Rel_prop_Z", "KO_prop_Z", "SynCom")
fam_6 <- fam_5[-1,]

for (syncom in SynComs){
  norm_table =read.table(paste(working_directory,"Isolate_tables/Original/", syncom, "_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
  #KO table
  KO_table =read.table(paste(working_directory,"KO_genome/KO_", syncom, ".tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
  colnames(KO_table) <- gsub("X", "", colnames(KO_table))
  
  if (syncom == "AtSC"){
    colnames(KO_table)[grep("M.16",colnames(KO_table))] <- "M-16"
    colnames(KO_table)[grep("M.6",colnames(KO_table))] <- "M-6"
    colnames(KO_table)[grep("M.10",colnames(KO_table))] <- "M-10"
    colnames(KO_table)[grep("M.11_2",colnames(KO_table))] <- "M-11_2"
    colnames(KO_table)[394] <- "M-11"
  }
  
  #taxonomy
  tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
  rownames(tax_df) <- tax_df$isolate
  tax_df_2 <- tax_df %>% dplyr::select (-isolate)
  colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")
  
  tax_df_3 <- table(tax_df_2$Family)
  tax_df_4 <- names(tax_df_3)[tax_df_3 > 10]
  
  #Samples TABLE
  samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
  rownames(samples_df) <- samples_df$sample_id
  samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
  
  #Subset for the right SynCom
  samples_df_3 <- subset(samples_df_2, samples_df_2$Compartment == "ES")
  samples_df_4 <- subset(samples_df_3, samples_df_3$Inoculum == paste(syncom))
  
  #Subset microbiome table for the right SynCom
  norm_table_2 <- norm_table[,colnames(norm_table) %in% row.names(samples_df_4)]
  
  #Subset taxonomy accordingly 
  if (syncom == "SSC"){
    tax_df_3 <- tax_df_2[tax_df_2$Family %in% tax_df_4,]
  } else {
    tax_df_3 <- tax_df_2[tax_df_2$Family %in% tax_df_4,]
    tax_df_3 <- tax_df_3[tax_df_3$SynCom == paste(syncom),]
  }
  
  #Set the OTU, TAX and sample data for making phyloseq object
  OTU = otu_table(as.matrix(norm_table_2),taxa_are_rows = TRUE)
  #TAX = tax_table(tax_mat)
  TAX = tax_table(as.matrix(tax_df_3))
  samples_sub = sample_data(samples_df_4)
  
  phylo = phyloseq(OTU,TAX, samples_sub)
  
  phylo_RA=microbiome::transform(x = phylo, transform = "compositional" )
  ps_family <- phyloseq::tax_glom(phylo, "Family")
  phylo_RA_fam=microbiome::transform(x = ps_family, transform = "compositional" )
  
  isolate_tab <- phylo_RA@otu_table
  OTU1 = as(otu_table(phylo_RA_fam), "matrix")
  TAX1 = as.data.frame(as(tax_table(phylo_RA_fam), "matrix"))
  
  row.names(OTU1) <- TAX1$Family
  Families <- unique(tax_df_3$Family)
  
  fam_3 <- data.frame(matrix(NA, ncol =7))
  colnames(fam_3) <- c("Isolate", "Rel", "Rel_prop", "KO", "KO_prop", "Family","Family_KO")
  fam_4 <- fam_3[-1,]
  
  for (family in Families) {
    isolate_set <- row.names(tax_df_3)[tax_df_3$Family == paste(family)]
    
    isolate_set_2 <- isolate_set[isolate_set %in% row.names(isolate_tab)]
    
    fam <- data.frame(matrix(NA, ncol = 3))
    colnames(fam) <- c("Isolate", "Rel", "KO")
    fam_2 <- fam[-1,]
    
    KO_table_2 <- KO_table[,colnames(KO_table) %in% isolate_set_2]
    veg_dist <- as.matrix(vegdist(t(KO_table_2)), method = "bray", diag = T)
    veg_dist_2 <- 1-veg_dist
    
    for (isolate in isolate_set_2){
      isolate_tab_2 <-isolate_tab[row.names(isolate_tab) == paste(isolate),]
      isolate_value <- rowSums(isolate_tab_2)/length(isolate_tab_2)
      names(isolate_value) <- NULL
      
      if (length(isolate_set_2) > 1){
        KO_table_3 <- KO_table_2[, colnames(KO_table_2) == paste(isolate)]
        KO_table_4 <- KO_table_3[KO_table_3 != 0]
        KO_value <- length(KO_table_4)
      } else {
        KO_table_4 <- KO_table_2[KO_table_2 != 0]
        KO_value <- length(KO_table_4)
      }
      
      new <- t(data.frame(c(paste(isolate), isolate_value, KO_value)))
      fam_2 <- rbind(fam_2, new)
    }
    
    fam_tab_2 <- OTU1[row.names(OTU1) == paste(family),]
    fam_value <- sum(fam_tab_2)/length(fam_tab_2)
    fam_2$V2 <- as.numeric(fam_2$V2)
    fam_2$V4 <- fam_2$V2/fam_value
    
    KO_table_fam <- KO_table[,colnames(KO_table) %in% isolate_set]
    if (length(isolate_set) >1){
      KO_table_fam_2 <- rowSums(KO_table_fam)
      KO_table_fam_3 <- KO_table_fam_2[KO_table_fam_2 != 0]
      fam_KO <- length(KO_table_fam_3)
    } else {
      fam_KO <- sum(KO_table_fam)
    }
    
    fam_2$V3 <- as.numeric(fam_2$V3)
    fam_2$V5 <- fam_2$V3/fam_KO
    fam_2$V6 <- paste(family)
    fam_2$V7 <- fam_KO
    fam_2$V8 <- scale(fam_2$V4)
    fam_2$V9 <- scale(fam_2$V5)
    
    colnames(fam_2) <- c("Isolate", "Rel", "KO", "Rel_Prop", "KO_prop", "Family", "Family_KO", "Rel_Prop_Z","KO_prop_Z")
    fam_4 <- rbind(fam_4, fam_2)
  }
  row.names(fam_4) <- NULL
  fam_4$SynCom <- paste(syncom)
  fam_6 <- rbind(fam_6, fam_4)
}

fam_6$Rel <- as.numeric(fam_6$Rel)

families <- unique(fam_6$Family)

table_cor_2 <- data.frame()

for (family in families){
  fam_7 <- fam_6[fam_6$Family == paste(family),]
  average <- sum(fam_7$KO_prop)/length(fam_7$KO_prop)
  
  SynCom_colors <- data.frame(c("AtSC", "HvSC", "LjSC"),c("#A3A500","#00B0F6","#00BF7D"))
  colnames(SynCom_colors) <- c("SynCom", "color")            
  
  fam_8 <- fam_7[order(fam_7$SynCom),]
  fam_colors <- SynCom_colors$color[SynCom_colors$SynCom %in% fam_8$SynCom]
  
  fam_8$Rel_log <- log2(fam_8$Rel_Prop)
  fam_8$Rel_log[fam_8$Rel_log == "-Inf"] <- 0
  
  for (syncom in SynComs){
    fam_sub <- fam_8[fam_8$SynCom == paste(syncom),]
    
    if (length(fam_sub$Isolate) < 6){
      pval <- NA
      cor <- NA
    } else {
      fam_sub_2 <- fam_sub[fam_sub$Rel_Prop != 0, ]
      val <- lm(formula = KO_prop ~ log(Rel_Prop), data = fam_sub_2)
      sum <- summary(val)
      cor <- sum$adj.r.squared
      pval_sum <- sum$coefficients
      pval <- pval_sum[8]
    }
    
    table_cor <- data.frame(t(data.frame(c(paste(syncom), cor, pval, paste(family)))))
    table_cor_2 <- rbind(table_cor_2, table_cor)
  }
}

row.names(table_cor_2) <- NULL
colnames(table_cor_2) <- c("SynCom", "R2", "Pvalue", "Family")

table_cor_2$R2 <- as.numeric(table_cor_2$R2)
table_cor_2$Pvalue <- as.numeric(table_cor_2$Pvalue)

table_cor_2$Sig <- ""
table_cor_2$Sig[table_cor_2$Pvalue < 0.05] <- "*"

new_fam <- data.frame()
#to get the order
for (family in families){
  table_cor_sub <- table_cor_2[table_cor_2$Family == paste(family),]
  table_cor_sub_2 <- table_cor_sub[!is.na(table_cor_sub$R2),]
  sum_val <- sum(table_cor_sub_2$R2)
  
  new_fam_2 <- data.frame(t(data.frame(c(paste(family), sum_val))))
  new_fam <- rbind(new_fam, new_fam_2)
}

new_fam_2 <- new_fam$X1[order(new_fam$X2, decreasing = F)]

table_cor_2$Family <- factor(table_cor_2$Family, levels = new_fam_2)

#Get the order based on family abundance in the data. 
#taxonomy
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv", sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")

tax_df_3 <- table(tax_df_2$Family)
tax_df_4 <- names(tax_df_3)[tax_df_3 > 10]

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep = ""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)

#Subset for the right SynCom
samples_df_3 <- subset(samples_df_2, samples_df_2$Compartment == "ES")
samples_df_4 <- subset(samples_df_3, samples_df_3$Condition != "NS")

#Subset microbiome table for the right SynCom
norm_table_2 <- norm_table[,colnames(norm_table) %in% row.names(samples_df_4)]

#Subset taxonomy accordingly 
tax_df_3 <- tax_df_2[tax_df_2$Family %in% tax_df_4,]

#Set the OTU, TAX and sample data for making phyloseq object
OTU = otu_table(as.matrix(norm_table_2),taxa_are_rows = TRUE)
#TAX = tax_table(tax_mat)
TAX = tax_table(as.matrix(tax_df_3))
samples_sub = sample_data(samples_df_4)

phylo = phyloseq(OTU,TAX, samples_sub)

phylo_RA=microbiome::transform(x = phylo, transform = "compositional" )
isolate_tab <- phylo_RA@otu_table

fam_order_2 <- data.frame()

for (family in families){
  tax_df_isolates <- row.names(tax_df_3)[tax_df_3$Family == paste(family)]
  isolate_tab_2 <- isolate_tab[row.names(isolate_tab) %in% tax_df_isolates,]
  average <- sum(rowSums(isolate_tab_2)/length(colnames(isolate_tab_2)))
  
  fam_order <- data.frame(t(data.frame(c(paste(family), average))))
  
  fam_order_2 <- rbind(fam_order_2, fam_order)
}

fam_order_3 <- fam_order_2$X1[order(fam_order_2$X2, decreasing =F)]
table_cor_2$Family <- factor(table_cor_2$Family, levels = fam_order_3)

fam_order_2$X3 <- "RA"
fam_order_2$X4 <- ""
fam_order_2$X5 <- "-"

fam_order_4 <- fam_order_2[,c(3,2,4,1,5)]
colnames(fam_order_4) <- colnames(table_cor_2)
table_cor_2 <- rbind(table_cor_2, fam_order_4)
row.names(table_cor_2) <- NULL
table_cor_2$R2 <- as.numeric(table_cor_2$R2)

table_cor_2$Sig[is.na(table_cor_2$Pvalue)] <- "-"
table_cor_2$Sig[table_cor_2$Pvalue < 0.05] <- "*"

table_cor_3 <- table_cor_2[table_cor_2$SynCom != "RA",]

table_cor_3$SynCom <- factor(table_cor_3$SynCom, levels = c("AtSC", "HvSC", "LjSC", "SSC"))

# Count number of isolates which belong to every family
Family_count <- data.frame(table(tax_df_2$Family))
Family_count$V2 <- paste(Family_count$Var1, " (n = ", Family_count$Freq, ")", sep ="")
# Create a named vector for mapping V1 to V2
mapping <- setNames(Family_count$V2, Family_count$Var1)
# Use the X_subset to reorder and subset the mapping to generate Y
fam_order_final <- mapping[fam_order_3]

table_cor_3$Family_2 <- Family_count$V2[match(table_cor_3$Family, Family_count$Var1)]
table_cor_3$Family_2 <- factor(table_cor_3$Family_2, levels = fam_order_final)

Plot_fam <- ggplot(table_cor_3, aes(SynCom, Family_2)) +
  geom_tile(height=0.98, mapping = aes(fill = R2)) +
  geom_text(aes(label = Sig), size =4) +
  scale_fill_gradient2(low = "#D55e00", mid = "white", high = "#56b4e9", midpoint =0, na.value = "lightgrey")+
  theme_classic() +
  labs(x ="Inoculum", y = "Family", fill = "R2") +
  theme(panel.background=element_blank(),panel.grid=element_blank(),axis.line.x=element_line(size=.5, colour="black"),axis.line.y=element_line(size=.5, colour="black"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black", size=7),legend.position="right",legend.text= element_text(size=10),text=element_text(family="sans", size=10))+
  theme(axis.text.x = element_text(size = 14, angle = 25,hjust=1),axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18), axis.text.y = element_text(size=14, face = rep("italic")), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=18)) +
  ggtitle("Root colonization versus functional diversity") +
  theme(plot.title = element_text(hjust = 0.5))
Plot_fam

pdf(paste(results.dir, "Figure_3d_family_heatmap.pdf", sep=""), width=6, height=6)
print(Plot_fam)
dev.off()

table_cor_4 <- table_cor_2[table_cor_2$SynCom == "RA",]
table_cor_4$Family <- factor(table_cor_4$Family, levels = fam_order_3)

plot_bar_fam <- ggplot(table_cor_4, aes(x=Family, y=R2)) + 
  geom_bar(stat = "identity", width = 0.98) +
  coord_flip() +
  theme(panel.background=element_blank(),panel.grid=element_blank(),axis.line.x=element_line(size=.5, colour="black"),axis.line.y=element_line(size=.5, colour="black"),axis.ticks=element_line(color="black")) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 14,angle = 25,hjust=1)) +
  ylab("Relative Abundance")
plot_bar_fam

pdf(paste(results.dir, "Figure_3d_family_bar.pdf", sep=""), width=1.5, height=6)
print(plot_bar_fam)
dev.off()

#Figure 3d - Genera heatmap
SynComs <- c("AtSC", "LjSC", "HvSC", "SSC")

genus_5 <- data.frame(matrix(NA, ncol =7))
colnames(genus_5) <- c("Isolate", "Rel", "Rel_prop", "KO", "KO_prop", "Genus", "SynCom")
genus_6 <- genus_5[-1,]

for (syncom in SynComs){
  norm_table =read.table(paste(working_directory, "Isolate_tables/Original/",syncom,"_norm.tsv",sep = ""), header=TRUE,sep="\t", row.names = 1)
  #KO table
  KO_table = read.table(paste(working_directory, "KO_genome/KO_",syncom, ".tsv", sep = ""), header = TRUE, sep = "\t", row.names =1)
  colnames(KO_table) <- gsub("X", "", colnames(KO_table))
  
  if (syncom == "AtSC"){
    colnames(KO_table)[grep("M.16",colnames(KO_table))] <- "M-16"
    colnames(KO_table)[grep("M.6",colnames(KO_table))] <- "M-6"
    colnames(KO_table)[grep("M.10",colnames(KO_table))] <- "M-10"
    colnames(KO_table)[grep("M.11_2",colnames(KO_table))] <- "M-11_2"
    colnames(KO_table)[394] <- "M-11"
  }
  
  if (syncom == "SSC"){
    colnames(KO_table)[grep("M.16",colnames(KO_table))] <- "M-16"
    colnames(KO_table)[grep("M.6",colnames(KO_table))] <- "M-6"
    colnames(KO_table)[grep("M.10",colnames(KO_table))] <- "M-10"
    colnames(KO_table)[grep("M.11_2",colnames(KO_table))] <- "M-11_2"
    colnames(KO_table)[571] <- "M-11"
  }
  
  #taxonomy
  tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv", sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
  rownames(tax_df) <- tax_df$isolate
  tax_df_2 <- tax_df %>% dplyr::select (-isolate)
  colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")
  
  tax_df_3 <- table(tax_df_2$Genus)
  tax_df_4 <- na.omit(names(tax_df_3)[tax_df_3 > 5])
  
  #Samples TABLE
  samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep = ""), header=TRUE,sep="\t") #make the SampleID column into the row.names
  rownames(samples_df) <- samples_df$sample_id
  samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
  
  #Subset for the right SynCom
  samples_df_3 <- subset(samples_df_2, samples_df_2$Compartment == "ES")
  samples_df_4 <- subset(samples_df_3, samples_df_3$Inoculum == paste(syncom))
  
  #Subset microbiome table for the right SynCom
  norm_table_2 <- norm_table[,colnames(norm_table) %in% row.names(samples_df_4)]
  
  #Subset taxonomy accordingly 
  if (syncom == "SSC"){
    tax_df_3 <- tax_df_2[tax_df_2$Genus %in% tax_df_4,]
  } else {
    tax_df_3 <- tax_df_2[tax_df_2$Genus %in% tax_df_4,]
    tax_df_3 <- tax_df_3[tax_df_3$SynCom == paste(syncom),]
  }
  
  #Set the OTU, TAX and sample data for making phyloseq object
  OTU = otu_table(as.matrix(norm_table_2),taxa_are_rows = TRUE)
  #TAX = tax_table(tax_mat)
  TAX = tax_table(as.matrix(tax_df_3))
  samples_sub = sample_data(samples_df_4)
  
  phylo = phyloseq(OTU,TAX, samples_sub)
  
  phylo_RA=microbiome::transform(x = phylo, transform = "compositional" )
  ps_genus <- phyloseq::tax_glom(phylo, "Genus")
  phylo_RA_gen=microbiome::transform(x = ps_genus, transform = "compositional" )
  
  isolate_tab <- phylo_RA@otu_table
  OTU1 = as(otu_table(phylo_RA_gen), "matrix")
  TAX1 = as.data.frame(as(tax_table(phylo_RA_gen), "matrix"))
  
  row.names(OTU1) <- TAX1$Genus
  Genera <- na.omit(unique(tax_df_3$Genus))
  
  genus_3 <- data.frame(matrix(NA, ncol =7))
  colnames(genus_3) <- c("Isolate", "Rel", "Rel_prop", "KO", "KO_prop", "Genus","Family_KO")
  genus_4 <- genus_3[-1,]
  
  for (genera in Genera) {
    isolate_set <- row.names(tax_df_3)[tax_df_3$Genus == paste(genera)]
    
    isolate_set_2 <- isolate_set[isolate_set %in% row.names(isolate_tab)]
    
    genus <- data.frame(matrix(NA, ncol = 3))
    colnames(genus) <- c("Isolate", "Rel", "KO")
    genus_2 <- genus[-1,]
    
    KO_table_2 <- KO_table[,colnames(KO_table) %in% isolate_set_2]
    veg_dist <- as.matrix(vegdist(t(KO_table_2)), method = "bray", diag = T)
    veg_dist_2 <- 1-veg_dist
    
    for (isolate in isolate_set_2){
      isolate_tab_2 <-isolate_tab[row.names(isolate_tab) == paste(isolate),]
      isolate_value <- rowSums(isolate_tab_2)/length(isolate_tab_2)
      names(isolate_value) <- NULL
      
      if (length(na.omit(isolate_set_2)) > 1){
        KO_table_3 <- KO_table_2[, colnames(KO_table_2) == paste(isolate)]
        KO_table_4 <- KO_table_3[KO_table_3 != 0]
        KO_value <- length(KO_table_4)
      } else {
        KO_table_4 <- KO_table_2[KO_table_2 != 0]
        KO_value <- length(KO_table_4)
      }
      
      new <- t(data.frame(c(paste(isolate), isolate_value, KO_value)))
      genus_2 <- rbind(genus_2, new)
    }
    
    genus_tab_2 <- OTU1[row.names(OTU1) == paste(genera),]
    genus_value <- sum(genus_tab_2)/length(genus_tab_2)
    genus_2$V2 <- as.numeric(genus_2$V2)
    genus_2$V4 <- genus_2$V2/genus_value
    
    KO_table_genus <- KO_table[,colnames(KO_table) %in% isolate_set]
    if (length(na.omit(isolate_set)) >1){
      KO_table_genus_2 <- rowSums(KO_table_genus)
      KO_table_genus_3 <- KO_table_genus_2[KO_table_genus_2 != 0]
      genus_KO <- length(KO_table_genus_3)
    } else {
      genus_KO <- sum(KO_table_genus)
    }
    
    if (length(genus_2$V1) != 0){
      genus_2$V3 <- as.numeric(genus_2$V3)
      genus_2$V5 <- genus_2$V3/genus_KO
      genus_2$V6 <- paste(genera)
      genus_2$V7 <- genus_KO
      colnames(genus_2) <- c("Isolate", "Rel", "KO", "Rel_Prop", "KO_prop", "Genus", "Genus_KO")
      genus_4 <- rbind(genus_4, genus_2)
    }
  }
  
  row.names(genus_4) <- NULL
  genus_4$SynCom <- paste(syncom)
  genus_6 <- rbind(genus_6, genus_4)
}

genus_6$Rel <- as.numeric(genus_6$Rel)

Genera <- unique(genus_6$Genus)

table_cor_2 <- data.frame()

for (genera in Genera){
  genus_7 <- genus_6[genus_6$Genus == paste(genera),]
  average <- sum(genus_7$KO_prop)/length(genus_7$KO_prop)
  
  SynCom_colors <- data.frame(c("AtSC", "HvSC", "LjSC"),c("#A3A500","#00B0F6","#00BF7D"))
  colnames(SynCom_colors) <- c("SynCom", "color")            
  
  genus_8 <- genus_7[order(genus_7$SynCom),]
  genus_colors <- SynCom_colors$color[SynCom_colors$SynCom %in% genus_8$SynCom]
  
  genus_8$Rel_log <- log2(genus_8$Rel_Prop)
  genus_8$Rel_log[genus_8$Rel_log == "-Inf"] <- 0
  
  for (syncom in SynComs){
    genus_sub <- genus_8[genus_8$SynCom == paste(syncom),]
    
    if (length(genus_sub$Isolate) < 6){
      pval <- NA
      cor <- NA
    } else {
      genus_sub_2 <- genus_sub[genus_sub$Rel_Prop != 0, ]
      val <- lm(formula = KO_prop ~ log(Rel_Prop), data = genus_sub_2)
      sum <- summary(val)
      cor <- sum$adj.r.squared
      pval_sum <- sum$coefficients
      pval <- pval_sum[8]
    }
    
    table_cor <- data.frame(t(data.frame(c(paste(syncom), cor, pval, paste(genera)))))
    table_cor_2 <- rbind(table_cor_2, table_cor)
  }
}

row.names(table_cor_2) <- NULL
colnames(table_cor_2) <- c("SynCom", "R2", "Pvalue", "Genus")

table_cor_2$R2 <- as.numeric(table_cor_2$R2)
table_cor_2$Pvalue <- as.numeric(table_cor_2$Pvalue)

table_cor_2$Sig <- ""
table_cor_2$Sig[table_cor_2$Pvalue < 0.05] <- "*"

new_genus <- data.frame()
#to get the order
for (genera in Genera){
  table_cor_sub <- table_cor_2[table_cor_2$Genus == paste(genera),]
  table_cor_sub_2 <- table_cor_sub[!is.na(table_cor_sub$R2),]
  sum_val <- sum(table_cor_sub_2$R2)
  
  new_genus_2 <- data.frame(t(data.frame(c(paste(genera), sum_val))))
  new_genus <- rbind(new_genus, new_genus_2)
}

new_genus_2 <- new_genus$X1[order(new_genus$X2, decreasing = F)]

table_cor_2$Genus <- factor(table_cor_2$Genus, levels = new_genus_2)

#Get the order based on family abundance in the data. 
#taxonomy
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv", sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")

tax_df_3 <- table(tax_df_2$Genus)
tax_df_4 <- names(tax_df_3)[tax_df_3 > 5]

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep = ""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)

#Subset for the right SynCom
samples_df_3 <- subset(samples_df_2, samples_df_2$Compartment == "ES")
samples_df_4 <- subset(samples_df_3, samples_df_3$Condition != "NS")

#Subset microbiome table for the right SynCom
norm_table_2 <- norm_table[,colnames(norm_table) %in% row.names(samples_df_4)]

#Subset taxonomy accordingly 
tax_df_3 <- tax_df_2[tax_df_2$Genus %in% tax_df_4,]

#Set the OTU, TAX and sample data for making phyloseq object
OTU = otu_table(as.matrix(norm_table_2),taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(tax_df_3))
samples_sub = sample_data(samples_df_4)

phylo = phyloseq(OTU,TAX, samples_sub)

phylo_RA=microbiome::transform(x = phylo, transform = "compositional" )
isolate_tab <- phylo_RA@otu_table

genus_order_2 <- data.frame()

for (genera in Genera){
  tax_df_isolates <- row.names(tax_df_3)[tax_df_3$Genus == paste(genera)]
  isolate_tab_2 <- isolate_tab[row.names(isolate_tab) %in% tax_df_isolates,]
  average <- sum(rowSums(isolate_tab_2)/length(colnames(isolate_tab_2)))
  
  genus_order <- data.frame(t(data.frame(c(paste(genera), average))))
  
  genus_order_2 <- rbind(genus_order_2, genus_order)
}

genus_order_3 <- genus_order_2$X1[order(genus_order_2$X2, decreasing =F)]

genus_order_2$X3 <- "RA"
genus_order_2$X4 <- ""
genus_order_2$X5 <- ""

genus_order_4 <- genus_order_2[,c(3,2,4,1,5)]
colnames(genus_order_4) <- colnames(table_cor_2)
table_cor_2 <- rbind(table_cor_2, genus_order_4)
row.names(table_cor_2) <- NULL
table_cor_2$R2 <- as.numeric(table_cor_2$R2)

table_cor_2$Sig[is.na(table_cor_2$Pvalue)] <- "-"
table_cor_2$Sig[table_cor_2$Pvalue < 0.05] <- "*"

table_cor_3 <- table_cor_2[table_cor_2$SynCom != "RA",]

table_cor_3$Genus <- factor(table_cor_3$Genus, levels = genus_order_3)
table_cor_3$SynCom <- factor(table_cor_3$SynCom, levels = c("AtSC", "HvSC", "LjSC", "SSC"))

table_cor_4 <- table_cor_3[table_cor_3$Genus %in% Genera,]
table_cor_4$Genus <- tax_df_2$Genus[match(table_cor_4$Genus, tax_df_2$Genus)]

Genera_order <- c("Acidovorax", "Rhizobacter", "Variovorax", "Polaromonas", "Cupriavidus", "Pelomonas")

# Count number of isolates which belong to every family
Genus_count <- data.frame(table(tax_df_2$Genus))
Genus_count$V2 <- paste(Genus_count$Var1, " (n = ", Genus_count$Freq, ")", sep ="")
# Create a named vector for mapping V1 to V2
mapping <- setNames(Genus_count$V2, Genus_count$Var1)
# Use the X_subset to reorder and subset the mapping to generate Y
gen_order_final <- mapping[Genera_order]

table_cor_4$Genus_2 <- Genus_count$V2[match(table_cor_3$Genus, Genus_count$Var1)]
table_cor_4 <- table_cor_4[table_cor_4$Genus_2 %in% gen_order_final,]
table_cor_4$Genus_2 <- factor(table_cor_4$Genus_2, levels = rev(gen_order_final))

Plot_genus <- ggplot(table_cor_4, aes(SynCom, Genus_2)) +
  geom_tile(aes(fill = R2)) +
  geom_text(aes(label = Sig), size =4) +
  scale_fill_gradient2(low = "#D55e00", mid = "white", high = "#56b4e9", midpoint =0, na.value = "lightgrey")+
  theme_classic() +
  labs(x ="SynCom", y = "Genus", fill = "R2") +
  theme(panel.background=element_blank(),panel.grid=element_blank(),axis.line.x=element_line(size=.5, colour="black"),axis.line.y=element_line(size=.5, colour="black"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black", size=7),legend.position="right",legend.text= element_text(size=10),text=element_text(family="sans", size=10))+
  theme(axis.text.x = element_text(size = 14, angle = 25,hjust=1),axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=18)) +
  theme(plot.title = element_text(hjust = 0.5))
Plot_genus

pdf(paste(results.dir, "Figure_3d_genus_heatmap.pdf", sep=""), width=6, height=3)
print(Plot_genus)
dev.off()

table_cor_5 <- table_cor_2[table_cor_2$SynCom == "RA",]
table_cor_5$Genus <- factor(table_cor_5$Genus, levels = genus_order_3)
table_cor_6 <- table_cor_5[table_cor_5$Genus %in% Genera_order,]
table_cor_6 <- table_cor_6[order(table_cor_6$R2, decreasing = TRUE), ]

plot_bar_genus <- ggplot(table_cor_6, aes(x=Genus, y=R2)) + 
  geom_bar(stat = "identity", width = 0.98) +
  coord_flip() +
  ylim(0,0.4) +
  theme(panel.background=element_blank(),panel.grid=element_blank(),axis.line.x=element_line(size=.5, colour="black"),axis.line.y=element_line(size=.5, colour="black"),axis.ticks=element_line(color="black")) +
  theme(axis.text.y = element_blank(),axis.title.x = element_text(size = 18), axis.title.y = element_blank()) +
  ylab("Relative Abundance")
plot_bar_genus

pdf(paste(results.dir, "Figure_3d_genus_bar.pdf", sep=""), width=1.5, height=2.117)
print(plot_bar_genus)
dev.off()

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

###Figure S13b - Intrafamily diversity =====
SynComs <- c("AtSC", "LjSC", "HvSC")

fam_5 <- data.frame(matrix(NA, ncol =7))
colnames(fam_5) <- c("Isolate", "Rel", "Rel_prop", "KO", "KO_prop", "Family", "SynCom")
fam_6 <- fam_5[-1,]

for (syncom in SynComs){
  norm_table =read.table(paste(working_directory,"Isolate_tables/Original/", syncom, "_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
  #KO table
  KO_table =read.table(paste(working_directory,"KO_genome/KO_", syncom, ".tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
  colnames(KO_table) <- gsub("X", "", colnames(KO_table))
  
  if (syncom == "AtSC"){
    colnames(KO_table)[grep("M.16",colnames(KO_table))] <- "M-16"
    colnames(KO_table)[grep("M.6",colnames(KO_table))] <- "M-6"
    colnames(KO_table)[grep("M.10",colnames(KO_table))] <- "M-10"
    colnames(KO_table)[grep("M.11_2",colnames(KO_table))] <- "M-11_2"
    colnames(KO_table)[394] <- "M-11"
  }
  
  #taxonomy
  tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
  rownames(tax_df) <- tax_df$isolate
  tax_df_2 <- tax_df %>% dplyr::select (-isolate)
  colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")
  
  tax_df_3 <- table(tax_df_2$Family)
  tax_df_4 <- names(tax_df_3)[tax_df_3 > 10]
  
  #Samples TABLE
  samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
  rownames(samples_df) <- samples_df$sample_id
  samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
  
  #Subset for the right SynCom
  samples_df_3 <- subset(samples_df_2, samples_df_2$Compartment == "ES")
  samples_df_4 <- subset(samples_df_3, samples_df_3$Inoculum == paste(syncom))
  
  #Subset microbiome table for the right SynCom
  norm_table_2 <- norm_table[,colnames(norm_table) %in% row.names(samples_df_4)]
  
  #Subset taxonomy accordingly 
  if (syncom == "SSC"){
    tax_df_3 <- tax_df_2[tax_df_2$Family %in% tax_df_4,]
  } else {
    tax_df_3 <- tax_df_2[tax_df_2$Family %in% tax_df_4,]
    tax_df_3 <- tax_df_3[tax_df_3$SynCom == paste(syncom),]
  }
  
  #Set the OTU, TAX and sample data for making phyloseq object
  OTU = otu_table(as.matrix(norm_table_2),taxa_are_rows = TRUE)
  #TAX = tax_table(tax_mat)
  TAX = tax_table(as.matrix(tax_df_3))
  samples_sub = sample_data(samples_df_4)
  
  phylo = phyloseq(OTU,TAX, samples_sub)
  
  phylo_RA=microbiome::transform(x = phylo, transform = "compositional" )
  ps_family <- phyloseq::tax_glom(phylo, "Family")
  phylo_RA_fam=microbiome::transform(x = ps_family, transform = "compositional" )
  
  isolate_tab <- phylo_RA@otu_table
  OTU1 = as(otu_table(phylo_RA_fam), "matrix")
  TAX1 = as.data.frame(as(tax_table(phylo_RA_fam), "matrix"))
  
  row.names(OTU1) <- TAX1$Family
  Families <- unique(tax_df_3$Family)
  
  fam_3 <- data.frame(matrix(NA, ncol =7))
  colnames(fam_3) <- c("Isolate", "Rel", "Rel_prop", "KO", "KO_prop", "Family","Family_KO")
  fam_4 <- fam_3[-1,]
  
  for (family in Families) {
    isolate_set <- row.names(tax_df_3)[tax_df_3$Family == paste(family)]
    
    isolate_set_2 <- isolate_set[isolate_set %in% row.names(isolate_tab)]
    
    fam <- data.frame(matrix(NA, ncol = 3))
    colnames(fam) <- c("Isolate", "Rel", "KO")
    fam_2 <- fam[-1,]
    
    KO_table_2 <- KO_table[,colnames(KO_table) %in% isolate_set_2]
    veg_dist <- as.matrix(vegdist(t(KO_table_2)), method = "bray", diag = T)
    veg_dist_2 <- 1-veg_dist
    
    for (isolate in isolate_set_2){
      isolate_tab_2 <-isolate_tab[row.names(isolate_tab) == paste(isolate),]
      isolate_value <- rowSums(isolate_tab_2)/length(isolate_tab_2)
      names(isolate_value) <- NULL
      
      if (length(isolate_set_2) > 1){
        KO_table_3 <- KO_table_2[, colnames(KO_table_2) == paste(isolate)]
        KO_table_4 <- KO_table_3[KO_table_3 != 0]
        KO_value <- length(KO_table_4)
      } else {
        KO_table_4 <- KO_table_2[KO_table_2 != 0]
        KO_value <- length(KO_table_4)
      }
      
      new <- t(data.frame(c(paste(isolate), isolate_value, KO_value)))
      fam_2 <- rbind(fam_2, new)
    }
    
    fam_tab_2 <- OTU1[row.names(OTU1) == paste(family),]
    fam_value <- sum(fam_tab_2)/length(fam_tab_2)
    fam_2$V2 <- as.numeric(fam_2$V2)
    fam_2$V4 <- fam_2$V2/fam_value
    
    KO_table_fam <- KO_table[,colnames(KO_table) %in% isolate_set]
    if (length(isolate_set) >1){
      KO_table_fam_2 <- rowSums(KO_table_fam)
      KO_table_fam_3 <- KO_table_fam_2[KO_table_fam_2 != 0]
      fam_KO <- length(KO_table_fam_3)
    } else {
      fam_KO <- sum(KO_table_fam)
    }
    
    fam_2$V3 <- as.numeric(fam_2$V3)
    fam_2$V5 <- fam_2$V3/fam_KO
    fam_2$V6 <- paste(family)
    fam_2$V7 <- fam_KO
    
    colnames(fam_2) <- c("Isolate", "Rel", "KO", "Rel_Prop", "KO_prop", "Family", "Family_KO")
    fam_4 <- rbind(fam_4, fam_2)
  }
  row.names(fam_4) <- NULL
  fam_4$SynCom <- paste(syncom)
  fam_6 <- rbind(fam_6, fam_4)
}

fam_6$Rel <- as.numeric(fam_6$Rel)

families <- unique(fam_6$Family)

plot_list <- list()
i <- 1

fam_6

for (family in families){
  fam_7 <- fam_6[fam_6$Family == paste(family),]
  average <- sum(fam_7$KO_prop)/length(fam_7$KO_prop)
  
  SynCom_colors <- data.frame(c("AtSC", "HvSC", "LjSC"),c("#A3A500","#00B0F6","#00BF7D"))
  colnames(SynCom_colors) <- c("SynCom", "color")            
  
  fam_8 <- fam_7[order(fam_7$SynCom),]
  fam_colors <- SynCom_colors$color[SynCom_colors$SynCom %in% fam_8$SynCom]
  
  plot <- ggplot(fam_8, aes(x = Rel_Prop, y = KO_prop, 
                            color = as.factor(SynCom))) +
    geom_point(size = 3) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) + 
    labs(x = "Relative abundance proportion",y = paste("KO Proportion (n = ", unique(fam_8$Family_KO), ")", sep = ""), color = "Inoculum") +
    ggtitle(paste(family, " - (n = ", length(fam_8$Isolate), ")", sep = "")) +
    ylim(0.4,1) +
    scale_color_manual(values = fam_colors)+
    geom_smooth(method="nls", se=FALSE, formula=y~a*log(x)+k,
                method.args=list(start=c(a=1, k=1))) +
    scale_size_continuous(limits = c(0,0.37)) +
    theme(axis.text.x = element_text(size = 10), axis.title = element_text(size = 14), axis.text.y = element_text(size=10), legend.title = element_text(size=16), legend.text = element_text(size=12), plot.title = element_text(size=18)) +
    guides(shape = guide_legend(override.aes = list(size = 5)))
  plot_list[[i]] <- plot
  i <- i + 1
}

all_plots <- ggarrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],plot_list[[6]],plot_list[[7]],plot_list[[8]],plot_list[[9]],plot_list[[10]],plot_list[[11]],plot_list[[12]],plot_list[[13]],plot_list[[14]],plot_list[[15]],plot_list[[16]],plot_list[[17]], nrow =5, ncol = 4, common.legend = T)

pdf(paste(results.dir,"Figure_S13b_Family_plot_loglinear.pdf", sep=""), width=18, height=20)
print(all_plots)
dev.off()

###Figure S14 - Intragenus diversity =====
SynComs <- c("AtSC", "LjSC", "HvSC")

genus_5 <- data.frame(matrix(NA, ncol =7))
colnames(genus_5) <- c("Isolate", "Rel", "Rel_prop", "KO", "KO_prop", "Genus", "SynCom")
genus_6 <- genus_5[-1,]

for (syncom in SynComs){
  norm_table =read.table(paste(working_directory,"Isolate_tables/Original/", syncom, "_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
  #KO table
  KO_table =read.table(paste(working_directory,"KO_genome/KO_", syncom, ".tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
  colnames(KO_table) <- gsub("X", "", colnames(KO_table))
  
  if (syncom == "AtSC"){
    colnames(KO_table)[grep("M.16",colnames(KO_table))] <- "M-16"
    colnames(KO_table)[grep("M.6",colnames(KO_table))] <- "M-6"
    colnames(KO_table)[grep("M.10",colnames(KO_table))] <- "M-10"
    colnames(KO_table)[grep("M.11_2",colnames(KO_table))] <- "M-11_2"
    colnames(KO_table)[394] <- "M-11"
  }
  
  #taxonomy
  tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
  rownames(tax_df) <- tax_df$isolate
  tax_df_2 <- tax_df %>% dplyr::select (-isolate)
  colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")
  tax_df_2_burk <- tax_df_2[tax_df_2$Family == "Burkholderiaceae",]
  
  tax_df_3 <- table(tax_df_2_burk$Genus)
  tax_df_4 <- names(tax_df_3)[tax_df_3 > 10]
  
  #Samples TABLE
  samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
  rownames(samples_df) <- samples_df$sample_id
  samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
  
  #Subset for the right SynCom
  samples_df_3 <- subset(samples_df_2, samples_df_2$Compartment == "ES")
  samples_df_4 <- subset(samples_df_3, samples_df_3$Inoculum == paste(syncom))
  
  #Subset microbiome table for the right SynCom
  norm_table_2 <- norm_table[,colnames(norm_table) %in% row.names(samples_df_4)]
  
  #Subset taxonomy accordingly 
  if (syncom == "SSC"){
    tax_df_3 <- tax_df_2[tax_df_2$Genus %in% tax_df_4,]
  } else {
    tax_df_3 <- tax_df_2[tax_df_2$Genus %in% tax_df_4,]
    tax_df_3 <- tax_df_3[tax_df_3$SynCom == paste(syncom),]
  }
  
  #Set the OTU, TAX and sample data for making phyloseq object
  OTU = otu_table(as.matrix(norm_table_2),taxa_are_rows = TRUE)
  #TAX = tax_table(tax_mat)
  TAX = tax_table(as.matrix(tax_df_3))
  samples_sub = sample_data(samples_df_4)
  
  phylo = phyloseq(OTU,TAX, samples_sub)
  
  phylo_RA=microbiome::transform(x = phylo, transform = "compositional" )
  ps_genus <- phyloseq::tax_glom(phylo, "Genus")
  phylo_RA_genus=microbiome::transform(x = ps_genus, transform = "compositional" )
  
  isolate_tab <- phylo_RA@otu_table
  OTU1 = as(otu_table(phylo_RA_genus), "matrix")
  TAX1 = as.data.frame(as(tax_table(phylo_RA_genus), "matrix"))
  
  row.names(OTU1) <- TAX1$Genus
  Genera <- unique(tax_df_3$Genus)
  
  genus_3 <- data.frame(matrix(NA, ncol =7))
  colnames(genus_3) <- c("Isolate", "Rel", "Rel_prop", "KO", "KO_prop", "Genus","Genus_KO")
  genus_4 <- genus_3[-1,]
  
  for (genus in Genera) {
    isolate_set <- row.names(tax_df_3)[tax_df_3$Genus == paste(genus)]
    
    isolate_set_2 <- isolate_set[isolate_set %in% row.names(isolate_tab)]
    
    genus_1 <- data.frame(matrix(NA, ncol = 3))
    colnames(genus_1) <- c("Isolate", "Rel", "KO")
    genus_2 <- genus_1[-1,]
    
    KO_table_2 <- KO_table[,colnames(KO_table) %in% isolate_set_2]
    veg_dist <- as.matrix(vegdist(t(KO_table_2)), method = "bray", diag = T)
    veg_dist_2 <- 1-veg_dist
    
    for (isolate in isolate_set_2){
      isolate_tab_2 <-isolate_tab[row.names(isolate_tab) == paste(isolate),]
      isolate_value <- rowSums(isolate_tab_2)/length(isolate_tab_2)
      names(isolate_value) <- NULL
      
      if (length(isolate_set_2) > 1){
        KO_table_3 <- KO_table_2[, colnames(KO_table_2) == paste(isolate)]
        KO_table_4 <- KO_table_3[KO_table_3 != 0]
        KO_value <- length(KO_table_4)
      } else {
        KO_table_4 <- KO_table_2[KO_table_2 != 0]
        KO_value <- length(KO_table_4)
      }
      
      new <- t(data.frame(c(paste(isolate), isolate_value, KO_value)))
      genus_2 <- rbind(genus_2, new)
    }
    
    genus_tab_2 <- OTU1[row.names(OTU1) == paste(genus),]
    genus_value <- sum(genus_tab_2)/length(genus_tab_2)
    genus_2$V2 <- as.numeric(genus_2$V2)
    genus_2$V4 <- genus_2$V2/genus_value
    
    KO_table_genus <- KO_table[,colnames(KO_table) %in% isolate_set]
    if (length(isolate_set) >1){
      KO_table_genus_2 <- rowSums(KO_table_genus)
      KO_table_genus_3 <- KO_table_genus_2[KO_table_genus_2 != 0]
      genus_KO <- length(KO_table_genus_3)
    } else {
      genus_KO <- sum(KO_table_genus)
    }
    
    genus_2$V3 <- as.numeric(genus_2$V3)
    genus_2$V5 <- genus_2$V3/genus_KO
    genus_2$V6 <- paste(genus)
    genus_2$V7 <- genus_KO
    
    colnames(genus_2) <- c("Isolate", "Rel", "KO", "Rel_Prop", "KO_prop", "Genus", "Genus_KO")
    genus_4 <- rbind(genus_4, genus_2)
  }
  row.names(genus_4) <- NULL
  genus_4$SynCom <- paste(syncom)
  genus_6 <- rbind(genus_6, genus_4)
}

genus_6$Rel <- as.numeric(genus_6$Rel)

genera <- unique(genus_6$Genus)

plot_list <- list()
i <- 1

for (genus in genera){
  genus_7 <- genus_6[genus_6$Genus == paste(genus),]
  average <- sum(genus_7$KO_prop)/length(genus_7$KO_prop)
  
  SynCom_colors <- data.frame(c("AtSC", "HvSC", "LjSC"),c("#A3A500","#00B0F6","#00BF7D"))
  colnames(SynCom_colors) <- c("SynCom", "color")            
  
  genus_8 <- genus_7[order(genus_7$SynCom),]
  genus_colors <- SynCom_colors$color[SynCom_colors$SynCom %in% genus_8$SynCom]
  
  plot <- ggplot(genus_8, aes(x = Rel_Prop, y = KO_prop, 
                              color = as.factor(SynCom))) +
    geom_point(size = 3) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) + 
    labs(x = "Relative abundance proportion",y = paste("KO Proportion (n = ", unique(genus_8$Genus_KO), ")", sep = ""), color = "Inoculum") +
    ggtitle(paste(genus, " - (n = ", length(genus_8$Isolate), ")", sep = "")) +
    ylim(0.4,1) +
    scale_color_manual(values = genus_colors)+
    geom_smooth(method="nls", se=FALSE, formula=y~a*log(x)+k,
                method.args=list(start=c(a=1, k=1))) +
    scale_size_continuous(limits = c(0,0.37)) +
    theme(axis.text.x = element_text(size = 10), axis.title = element_text(size = 14), axis.text.y = element_text(size=10), legend.title = element_text(size=16), legend.text = element_text(size=12), plot.title = element_text(size=14)) +
    guides(shape = guide_legend(override.aes = list(size = 5)))
  plot_list[[i]] <- plot
  i <- i + 1
}

all_plots <- ggarrange(plot_list[[2]],plot_list[[1]],plot_list[[3]],plot_list[[4]],plot_list[[5]],plot_list[[6]], nrow =2, ncol = 3, common.legend = T)

pdf(paste(results.dir,"Figure_S14_Genus_plot_loglinear.pdf", sep=""), width=12, height=8)
print(all_plots)
dev.off()

###Table S5 - Statistical randomization test =====
SynComs <- c("AtSC", "LjSC", "HvSC", "SSC")

fam_5 <- data.frame(matrix(NA, ncol =9))
colnames(fam_5) <- c("Isolate", "Rel", "Rel_prop", "KO", "KO_prop", "Family", "Rel_prop_Z", "KO_prop_Z", "SynCom")
fam_6 <- fam_5[-1,]

for (syncom in SynComs){
  norm_table =read.table(paste(working_directory,"Isolate_tables/Original/", syncom, "_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
  #KO table
  KO_table =read.table(paste(working_directory,"KO_genome/KO_", syncom, ".tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
  colnames(KO_table) <- gsub("X", "", colnames(KO_table))
  
  if (syncom == "AtSC"){
    colnames(KO_table)[grep("M.16",colnames(KO_table))] <- "M-16"
    colnames(KO_table)[grep("M.6",colnames(KO_table))] <- "M-6"
    colnames(KO_table)[grep("M.10",colnames(KO_table))] <- "M-10"
    colnames(KO_table)[grep("M.11_2",colnames(KO_table))] <- "M-11_2"
    colnames(KO_table)[394] <- "M-11"
  }
  
  #taxonomy
  tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
  rownames(tax_df) <- tax_df$isolate
  tax_df_2 <- tax_df %>% dplyr::select (-isolate)
  colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")
  
  tax_df_3 <- table(tax_df_2$Family)
  tax_df_4 <- names(tax_df_3)[tax_df_3 > 10]
  
  #Samples TABLE
  samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
  rownames(samples_df) <- samples_df$sample_id
  samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
  
  #Subset for the right SynCom
  samples_df_3 <- subset(samples_df_2, samples_df_2$Compartment == "ES")
  samples_df_4 <- subset(samples_df_3, samples_df_3$Inoculum == paste(syncom))
  
  #Subset microbiome table for the right SynCom
  norm_table_2 <- norm_table[,colnames(norm_table) %in% row.names(samples_df_4)]
  
  #Subset taxonomy accordingly 
  if (syncom == "SSC"){
    tax_df_3 <- tax_df_2[tax_df_2$Family %in% tax_df_4,]
  } else {
    tax_df_3 <- tax_df_2[tax_df_2$Family %in% tax_df_4,]
    tax_df_3 <- tax_df_3[tax_df_3$SynCom == paste(syncom),]
  }
  
  #Set the OTU, TAX and sample data for making phyloseq object
  OTU = otu_table(as.matrix(norm_table_2),taxa_are_rows = TRUE)
  #TAX = tax_table(tax_mat)
  TAX = tax_table(as.matrix(tax_df_3))
  samples_sub = sample_data(samples_df_4)
  
  phylo = phyloseq(OTU,TAX, samples_sub)
  
  phylo_RA=microbiome::transform(x = phylo, transform = "compositional" )
  ps_family <- phyloseq::tax_glom(phylo, "Family")
  phylo_RA_fam=microbiome::transform(x = ps_family, transform = "compositional" )
  
  isolate_tab <- phylo_RA@otu_table
  OTU1 = as(otu_table(phylo_RA_fam), "matrix")
  TAX1 = as.data.frame(as(tax_table(phylo_RA_fam), "matrix"))
  
  row.names(OTU1) <- TAX1$Family
  Families <- unique(tax_df_3$Family)
  
  fam_3 <- data.frame(matrix(NA, ncol =7))
  colnames(fam_3) <- c("Isolate", "Rel", "Rel_prop", "KO", "KO_prop", "Family","Family_KO")
  fam_4 <- fam_3[-1,]
  
  for (family in Families) {
    isolate_set <- row.names(tax_df_3)[tax_df_3$Family == paste(family)]
    
    isolate_set_2 <- isolate_set[isolate_set %in% row.names(isolate_tab)]
    
    fam <- data.frame(matrix(NA, ncol = 3))
    colnames(fam) <- c("Isolate", "Rel", "KO")
    fam_2 <- fam[-1,]
    
    KO_table_2 <- KO_table[,colnames(KO_table) %in% isolate_set_2]
    veg_dist <- as.matrix(vegdist(t(KO_table_2)), method = "bray", diag = T)
    veg_dist_2 <- 1-veg_dist
    
    for (isolate in isolate_set_2){
      isolate_tab_2 <-isolate_tab[row.names(isolate_tab) == paste(isolate),]
      isolate_value <- rowSums(isolate_tab_2)/length(isolate_tab_2)
      names(isolate_value) <- NULL
      
      if (length(isolate_set_2) > 1){
        KO_table_3 <- KO_table_2[, colnames(KO_table_2) == paste(isolate)]
        KO_table_4 <- KO_table_3[KO_table_3 != 0]
        KO_value <- length(KO_table_4)
      } else {
        KO_table_4 <- KO_table_2[KO_table_2 != 0]
        KO_value <- length(KO_table_4)
      }
      
      new <- t(data.frame(c(paste(isolate), isolate_value, KO_value)))
      fam_2 <- rbind(fam_2, new)
    }
    
    fam_tab_2 <- OTU1[row.names(OTU1) == paste(family),]
    fam_value <- sum(fam_tab_2)/length(fam_tab_2)
    fam_2$V2 <- as.numeric(fam_2$V2)
    fam_2$V4 <- fam_2$V2/fam_value
    
    KO_table_fam <- KO_table[,colnames(KO_table) %in% isolate_set]
    if (length(isolate_set) >1){
      KO_table_fam_2 <- rowSums(KO_table_fam)
      KO_table_fam_3 <- KO_table_fam_2[KO_table_fam_2 != 0]
      fam_KO <- length(KO_table_fam_3)
    } else {
      fam_KO <- sum(KO_table_fam)
    }
    
    fam_2$V3 <- as.numeric(fam_2$V3)
    fam_2$V5 <- fam_2$V3/fam_KO
    fam_2$V6 <- paste(family)
    fam_2$V7 <- fam_KO
    fam_2$V8 <- scale(fam_2$V4)
    fam_2$V9 <- scale(fam_2$V5)
    
    colnames(fam_2) <- c("Isolate", "Rel", "KO", "Rel_Prop", "KO_prop", "Family", "Family_KO", "Rel_Prop_Z","KO_prop_Z")
    fam_4 <- rbind(fam_4, fam_2)
  }
  row.names(fam_4) <- NULL
  fam_4$SynCom <- paste(syncom)
  fam_6 <- rbind(fam_6, fam_4)
}

fam_6$Rel <- as.numeric(fam_6$Rel)

families <- unique(fam_6$Family)

table_cor_2 <- data.frame()

for (family in families){
  fam_7 <- fam_6[fam_6$Family == paste(family),]
  average <- sum(fam_7$KO_prop)/length(fam_7$KO_prop)
  
  SynCom_colors <- data.frame(c("AtSC", "HvSC", "LjSC"),c("#A3A500","#00B0F6","#00BF7D"))
  colnames(SynCom_colors) <- c("SynCom", "color")            
  
  fam_8 <- fam_7[order(fam_7$SynCom),]
  fam_colors <- SynCom_colors$color[SynCom_colors$SynCom %in% fam_8$SynCom]
  
  fam_8$Rel_log <- log2(fam_8$Rel_Prop)
  fam_8$Rel_log[fam_8$Rel_log == "-Inf"] <- 0
  
  for (syncom in SynComs){
    fam_sub <- fam_8[fam_8$SynCom == paste(syncom),]
    if (length(fam_sub$Isolate) < 6){
      pval <- NA
      cor <- NA
    } else {
      fam_sub_2 <- fam_sub[fam_sub$Rel_Prop != 0, ]
      
      for (i in 1:1000){
        fam_sub_3 <- sample(fam_sub_2$KO_prop)
        fam_sub_2$KO_prop <- fam_sub_3
        val <- lm(formula = KO_prop ~ log(Rel_Prop), data = fam_sub_2)
        sum <- summary(val)
        cor <- sum$adj.r.squared
        pval_sum <- sum$coefficients
        pval <- pval_sum[8]
        table_cor <- data.frame(t(data.frame(c(paste(syncom), cor, pval, paste(family), paste(i)))))
        colnames(table_cor) <- c("SynCom", "R2", "Pvalue", "Family", "Simulation")
        
        table_cor_2 <- rbind(table_cor_2, table_cor)
      }
    }
  }
}

#Real data 
for (family in families){
  fam_7 <- fam_6[fam_6$Family == paste(family),]
  average <- sum(fam_7$KO_prop)/length(fam_7$KO_prop)
  
  SynCom_colors <- data.frame(c("AtSC", "HvSC", "LjSC"),c("#A3A500","#00B0F6","#00BF7D"))
  colnames(SynCom_colors) <- c("SynCom", "color")            
  
  fam_8 <- fam_7[order(fam_7$SynCom),]
  fam_colors <- SynCom_colors$color[SynCom_colors$SynCom %in% fam_8$SynCom]
  
  fam_8$Rel_log <- log2(fam_8$Rel_Prop)
  fam_8$Rel_log[fam_8$Rel_log == "-Inf"] <- 0
  
  for (syncom in SynComs){
    fam_sub <- fam_8[fam_8$SynCom == paste(syncom),]
    
    if (length(fam_sub$Isolate) < 6){
      pval <- NA
      cor <- NA
    } else {
      fam_sub_2 <- fam_sub[fam_sub$Rel_Prop != 0, ]
      val <- lm(formula = KO_prop ~ log(Rel_Prop), data = fam_sub_2)
      sum <- summary(val)
      cor <- sum$adj.r.squared
      pval_sum <- sum$coefficients
      pval <- pval_sum[8]
    }
    
    table_cor <- data.frame(t(data.frame(c(paste(syncom), cor, pval, paste(family), "Real"))))
    colnames(table_cor) <- c("SynCom", "R2", "Pvalue", "Family", "Simulation")
    table_cor_2 <- rbind(table_cor_2, table_cor)
  }
}

row.names(table_cor_2) <- NULL
colnames(table_cor_2) <- c("SynCom", "R2", "Pvalue", "Family", "Simulation")

table_cor_2$R2 <- as.numeric(table_cor_2$R2)
table_cor_2$Pvalue <- as.numeric(table_cor_2$Pvalue)

table_cor_2$Sig <- ""
table_cor_2$Sig[table_cor_2$Pvalue < 0.05] <- "*"
new_data <- data.frame()

for (family in families){
  for (syncom in SynComs){
    table_cor_2_sub <- table_cor_2[table_cor_2$SynCom == paste(syncom) & table_cor_2$Family == paste(family) & table_cor_2$Simulation != "Real",]
    table_cor_2_sub$Pvalue
    threshold <- quantile(as.numeric(table_cor_2_sub$Pvalue), probs = 0.025)
    real <- table_cor_2$Pvalue[table_cor_2$SynCom == paste(syncom) & table_cor_2$Family == paste(family) & table_cor_2$Simulation == "Real"]
    
    table_cor <- data.frame(t(data.frame(c(paste(syncom), paste(family), threshold, real))))
    new_data <- rbind(new_data, table_cor)
  }
}

row.names(new_data) <- NULL
colnames(new_data) <- c("SynCom","Family","pvalue_conf_interval_0.025_simulations","real_pvalue")

write.table(new_data, paste(results.dir, "Table_S5_Statistical_randomization_test.tsv",sep =""), row.names =F, col.names =T, quote =F, sep = "\t")


###Figure 3e - AtSC 7x200 validation of functional diversity vs root colonization =====
Families <- c("Burkholderiaceae","Rhizobiaceae","Xanthomonadaceae","Pseudomonadaceae","Caulobacteraceae","Beijerinckiaceae")

tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
KO_iso <- read.table(paste(working_directory,"KO_genome/KO_SSC.tsv", sep = ""), row.names =1, header =T)
tax_df <- tax_df[tax_df$SynCom == "AtSC",]

colnames(KO_iso) <- gsub("X","", colnames(KO_iso))
colnames(KO_iso) <- gsub("M.11_2","M-11_2", colnames(KO_iso))
colnames(KO_iso) <- gsub("M.11","M-11", colnames(KO_iso))
colnames(KO_iso) <- gsub("M.16","M-16", colnames(KO_iso))
colnames(KO_iso) <- gsub("M.10","M-10", colnames(KO_iso))
colnames(KO_iso) <- gsub("M.6","M-6", colnames(KO_iso))

#Filter out samples with high contamination
plant_reads <- read.table(paste(working_directory, "AtSC_7_SynCom_experiment/plant_reads.tsv", sep = ""), header =T)  
SynCom_reads <- read.table(paste(working_directory, "AtSC_7_SynCom_experiment/mapped_reads.tsv", sep = ""), header =T)

plant_reads$SynCom_reads_2 <- SynCom_reads$Full_200[match(plant_reads$sample_id, SynCom_reads$Sample)]
plant_reads$SynCom_reads <- SynCom_reads$Map_200[match(plant_reads$sample_id, SynCom_reads$Sample)]

plant_reads$non_plant_pairs <- plant_reads$total_reads-plant_reads$plant_reads

plant_reads$other_isolate_reads <- plant_reads$SynCom_reads_2 - plant_reads$SynCom_reads
plant_reads$Contaminant_reads <- plant_reads$non_plant_pairs - plant_reads$SynCom_reads_2

Hank_the_normalizer <- function(df,group,amount){
  df_2 <- df %>% dplyr::group_by_at(group) %>% dplyr::summarise(total=sum(.data[[amount]]))
  df_3 <- df_2$total
  names(df_3) <- df_2[[group]]
  df$total <- df_3[as.character(df[[group]])]
  df$Rel <- df[[amount]] / df$total
  return(df)
}

plant_reads_2 <- plant_reads[,c(1,3,5,7,8)]

plant_reads_melt <- melt(plant_reads_2)
plant_reads_melt_2 <- Hank_the_normalizer(plant_reads_melt,"sample_id","value")

removal_sam <- c()

for (sample in unique(plant_reads_melt_2$sample_id)){
  plant_reads_melt_2_sub <- plant_reads_melt_2[plant_reads_melt_2$sample_id == paste(sample),]
  value <- plant_reads_melt_2_sub$Rel[plant_reads_melt_2_sub$variable == "SynCom_reads"]/sum(plant_reads_melt_2_sub$Rel[plant_reads_melt_2_sub$variable != "plant_reads"])
  if (value <= 0.6){
    removal_sam <- c(removal_sam, paste(sample))
  }
}

#Microbiome data
table_bac <- read.table(paste(working_directory,"AtSC_7_SynCom_experiment/isolate_norm.txt", sep = ""), row.names =1, header =T)
table_bac_1 <- t(t(table_bac) / rowSums(t(table_bac)))
table_bac_2 <- table_bac_1[,!colnames(table_bac_1) %in% removal_sam] 

#metadata
metadata <- read.table(paste(working_directory,"AtSC_7_SynCom_experiment/metadata.txt", sep= ""), row.names =1, header =T)
metadata_2 <- metadata[grep("Syncom",metadata$Syncom),]
SynComs <- unique(metadata_2$Syncom)

#isolates
groups <- read.table(paste(working_directory,"AtSC_7_SynCom_experiment/SynCom_isolates.txt", sep = ""), row.names =1, header =F)

new_datafr <- data.frame()

for (syncom in SynComs){
  #Observed
  metadata_sub <- metadata_2[metadata_2$Syncom == paste(syncom),]
  metadata_sub_2 <- row.names(metadata_sub)[metadata_sub$Compartment != "inoculum"]
  
  isolates_in_SynCom <- as.vector(unlist(as.vector(groups[row.names(groups) == paste(syncom),])))
  tax_df_sub <- tax_df[tax_df$isolate %in% isolates_in_SynCom,]
  
  table_bac_sub <- data.frame(table_bac_2[,colnames(table_bac_2) %in% metadata_sub_2])
  
  for (family in Families){
    tax_df_sub_2 <- tax_df_sub$isolate[tax_df_sub$family == paste(family)]
    
    #Expected
    KO_iso_2 <- KO_iso[,colnames(KO_iso) %in% tax_df_sub_2]
    KO_iso_2[KO_iso_2 > 0] <- 1
    
    KO_iso_3 <- data.frame(colSums(KO_iso_2))
    colnames(KO_iso_3) <- "No_of_KOs"
    KO_iso_4 <- KO_iso_3 %>% dplyr::arrange(desc(No_of_KOs))
    
    m <- Mclust(KO_iso_4$No_of_KOs)     
    KO_iso_4$group <- m$classification
    
    #Observed
    table_bac_sub_2 <- table_bac_sub[row.names(table_bac_sub) %in% tax_df_sub_2, ]
    out <- data.frame(rowSums(table_bac_sub_2)/length(colnames(table_bac_sub_2)))
    colnames(out) <- "RA"
    
    for (cluster in unique(KO_iso_4$group)){
      KO_iso_4_sub <- row.names(KO_iso_4)[KO_iso_4$group == paste(cluster)]
      out_sub <- sum(out[row.names(out) %in% KO_iso_4_sub,])
      KO_iso_sub <- KO_iso_4[KO_iso_4$group == paste(cluster),]
      value_KO <- sum(KO_iso_sub$No_of_KOs)/length(KO_iso_sub$No_of_KOs)
      new_datafr <- rbind(new_datafr, data.frame(t(data.frame(c(paste(syncom),paste(family),paste(cluster),value_KO, out_sub)))))
    }
  }
}

row.names(new_datafr) <- NULL
colnames(new_datafr) <- c("SynCom","Family","Cluster","No_of_KO","Cum_RA")

new_datafr$SynFam <- paste(new_datafr$SynCom, new_datafr$Family, sep = "_")

new_datafr_2 <- data.frame()

for (group in unique(new_datafr$SynFam)){
  new_datafr_sub <- new_datafr[new_datafr$SynFam == paste(group),]
  new_datafr_sub_2 <- new_datafr_sub[order(as.numeric(new_datafr_sub$Cum_RA), decreasing = F), ]
  new_datafr_sub_2$order_RA <- 1:length(new_datafr_sub_2$SynFam)
  new_datafr_2 <- rbind(new_datafr_2,new_datafr_sub_2 )
}

row.names(new_datafr_2) <- NULL
colnames(new_datafr_2) <- c("SynCom","Family","Exp_order","No_of_KO","Cum_RA","SynFam", "Obs_order")

stat_data <- data.frame()

for (group in unique(new_datafr_2$Family)){
  new_datafr_2_sub <- new_datafr_2[new_datafr_2$Family == paste(group),]
  check_table <- table(new_datafr_2_sub$SynFam)
  check_table_2 <- names(check_table)[check_table == 1]
  new_datafr_2_sub <- new_datafr_2_sub[!new_datafr_2_sub$SynFam %in% check_table_2,]
  
  if (length(new_datafr_2_sub$Family) > 2){
    stat_out <- cor.test(as.numeric(new_datafr_2_sub$Exp_order), as.numeric(new_datafr_2_sub$Obs_order), method = "kendall")
    pval <- stat_out$p.value
    stat_data <- rbind(stat_data, data.frame(t(data.frame(c(paste(group), pval)))))
  } else {
    pval <- NA
    SC <- unique(new_datafr_2$SynCom[new_datafr_2$SynFam == paste(group)])
    FAM <- unique(new_datafr_2$Family[new_datafr_2$SynFam == paste(group)])
    stat_data <- rbind(stat_data, data.frame(t(data.frame(c(paste(SC), paste(FAM), pval)))))
  }
}

row.names(stat_data) <- NULL
colnames(stat_data) <- c("Family","Kendall correlation p-value")

stat_data$Pos_cor <- c(4/7,5/5,3/7,7/7,7/7,5/5)
stat_data$neglogp <- -1*log10(as.numeric(stat_data$`Kendall correlation p-value`))

new_datafr$Cum_RA <- as.numeric(new_datafr$Cum_RA)
new_datafr_new <- data.frame()

for (group in unique(new_datafr$SynFam)){
  new_datafr_sub <- new_datafr[new_datafr$SynFam == paste(group),]
  
  if (length(new_datafr_sub$SynFam > 1)){
    cum_value <- sum(new_datafr_sub$Cum_RA)
  } else {
    cum_value <- new_datafr_sub$Cum_RA
  }
  fam_val <- unique(new_datafr_sub$Family)
  set <- data.frame(t(data.frame(c(paste(group), cum_value,fam_val))))
  new_datafr_new <- rbind(new_datafr_new, set)
}

row.names(new_datafr_new) <- NULL
colnames(new_datafr_new) <- c("SynFam", "Cum_RA", "Family")
new_datafr_new$Cum_RA <- as.numeric(new_datafr_new$Cum_RA)

order <- c("Pseudomonadaceae", "Beijerinckiaceae", "Xanthomonadaceae", "Burkholderiaceae", "Rhizobiaceae", "Caulobacteraceae")
new_datafr_new$Family <- factor(new_datafr_new$Family, levels = order)

plot_col <- ggplot(new_datafr_new, aes(x = Family, y = Cum_RA),color="#515455") +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  ylab("Cumulative RA") +
  geom_jitter(shape = 16,color="#515455", position = position_jitter(0.2), show.legend = TRUE) +
  theme(panel.background=element_blank(),panel.grid=element_blank(),axis.line.x=element_line(size=.5, colour="black"),axis.line.y=element_line(size=.5, colour="black"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black", size=7),legend.position="right",legend.text= element_text(size=10),text=element_text(family="sans", size=10))+
  theme(axis.text.x = element_text(size = 14),axis.title.x = element_text(size = 18), axis.title.y = element_blank(), axis.text.y = element_blank(), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=18, hjust = 0.5)) +
  coord_flip() + ggtitle("Cumulative RA Family isolates - 7x200 AtSC") 
plot_col

stat_data$Family <- factor(stat_data$Family, levels = order)

dot_plot <- ggplot(stat_data, aes(y=Family, x=neglogp)) + geom_point(color="#515455",size =5 ) + theme_classic() + ggtitle("Kendall rank correlation") + theme(plot.title = element_text(hjust = 0.5)) +  
  labs(x ="-log10(p-value)", y = "Family") +
  geom_vline(xintercept=-1*log10(0.05), linetype = "dashed" ) +
  theme(panel.background=element_blank(),panel.grid=element_blank(),axis.line.x=element_line(size=.5, colour="black"),axis.line.y=element_line(size=.5, colour="black"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black", size=7),legend.position="right",legend.text= element_text(size=10),text=element_text(family="sans", size=10))+
  theme(axis.text.x = element_text(size = 14),axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18), axis.text.y = element_text(size=14, face = rep("italic")), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=18))
print(dot_plot)

tog_plot <- ggarrange(dot_plot,plot_col, ncol =2)

pdf(paste(results.dir,"Figure_3e_AtSC_7x200_validation.pdf", sep=""), width=12, height=6)
print(tog_plot)
dev.off()

###Figure S15 - Panel A Heatmap =====
Families <- c("Burkholderiaceae","Rhizobiaceae","Xanthomonadaceae","Pseudomonadaceae","Caulobacteraceae","Beijerinckiaceae")

tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
KO_iso <- read.table(paste(working_directory,"KO_genome/KO_SSC.tsv", sep = ""), row.names =1, header =T)
tax_df <- tax_df[tax_df$SynCom == "AtSC",]

colnames(KO_iso) <- gsub("X","", colnames(KO_iso))
colnames(KO_iso) <- gsub("M.11_2","M-11_2", colnames(KO_iso))
colnames(KO_iso) <- gsub("M.11","M-11", colnames(KO_iso))
colnames(KO_iso) <- gsub("M.16","M-16", colnames(KO_iso))
colnames(KO_iso) <- gsub("M.10","M-10", colnames(KO_iso))
colnames(KO_iso) <- gsub("M.6","M-6", colnames(KO_iso))

#Filter out samples with high contamination
plant_reads <- read.table(paste(working_directory, "AtSC_7_SynCom_experiment/plant_reads.tsv", sep = ""), header =T)  
SynCom_reads <- read.table(paste(working_directory, "AtSC_7_SynCom_experiment/mapped_reads.tsv", sep = ""), header =T)

plant_reads$SynCom_reads_2 <- SynCom_reads$Full_200[match(plant_reads$sample_id, SynCom_reads$Sample)]
plant_reads$SynCom_reads <- SynCom_reads$Map_200[match(plant_reads$sample_id, SynCom_reads$Sample)]

plant_reads$non_plant_pairs <- plant_reads$total_reads-plant_reads$plant_reads

plant_reads$other_isolate_reads <- plant_reads$SynCom_reads_2 - plant_reads$SynCom_reads
plant_reads$Contaminant_reads <- plant_reads$non_plant_pairs - plant_reads$SynCom_reads_2

Hank_the_normalizer <- function(df,group,amount){
  df_2 <- df %>% dplyr::group_by_at(group) %>% dplyr::summarise(total=sum(.data[[amount]]))
  df_3 <- df_2$total
  names(df_3) <- df_2[[group]]
  df$total <- df_3[as.character(df[[group]])]
  df$Rel <- df[[amount]] / df$total
  return(df)
}

plant_reads_2 <- plant_reads[,c(1,3,5,7,8)]

plant_reads_melt <- melt(plant_reads_2)
plant_reads_melt_2 <- Hank_the_normalizer(plant_reads_melt,"sample_id","value")

removal_sam <- c()

for (sample in unique(plant_reads_melt_2$sample_id)){
  plant_reads_melt_2_sub <- plant_reads_melt_2[plant_reads_melt_2$sample_id == paste(sample),]
  value <- plant_reads_melt_2_sub$Rel[plant_reads_melt_2_sub$variable == "SynCom_reads"]/sum(plant_reads_melt_2_sub$Rel[plant_reads_melt_2_sub$variable != "plant_reads"])
  if (value <= 0.6){
    removal_sam <- c(removal_sam, paste(sample))
  }
}

#Microbiome data
table_bac <- read.table(paste(working_directory,"AtSC_7_SynCom_experiment/isolate_norm.txt", sep = ""), row.names =1, header =T)
table_bac_1 <- t(t(table_bac) / rowSums(t(table_bac)))
table_bac_2 <- table_bac_1[,!colnames(table_bac_1) %in% removal_sam] 

#metadata
metadata <- read.table(paste(working_directory,"AtSC_7_SynCom_experiment/metadata.txt", sep= ""), row.names =1, header =T)
metadata_2 <- metadata[grep("Syncom",metadata$Syncom),]
SynComs <- unique(metadata_2$Syncom)

#isolates
groups <- read.table(paste(working_directory,"AtSC_7_SynCom_experiment/SynCom_isolates.txt", sep = ""), row.names =1, header =F)

new_data_next <- data.frame()

for (syncom in SynComs){
  #Observed
  metadata_sub <- metadata_2[metadata_2$Syncom == paste(syncom),]
  metadata_sub_2 <- row.names(metadata_sub)[metadata_sub$Compartment != "inoculum"]
  
  isolates_in_SynCom <- as.vector(unlist(as.vector(groups[row.names(groups) == paste(syncom),])))
  tax_df_sub <- tax_df[tax_df$isolate %in% isolates_in_SynCom,]
  
  table_bac_sub <- data.frame(table_bac_2[,colnames(table_bac_2) %in% metadata_sub_2])
  
  for (family in Families){
    tax_df_sub_2 <- tax_df_sub$isolate[tax_df_sub$family == paste(family)]
    
    #Expected
    KO_iso_2 <- KO_iso[,colnames(KO_iso) %in% tax_df_sub_2]
    KO_iso_2[KO_iso_2 > 0] <- 1
    
    KO_iso_2_sum <- rowSums(KO_iso_2)
    KO_iso_2_sum_val <- length(names(KO_iso_2_sum)[KO_iso_2_sum != 0])
    
    KO_iso_3 <- data.frame(colSums(KO_iso_2))
    colnames(KO_iso_3) <- "No_of_KOs"
    
    #Observed
    table_bac_sub_2 <- table_bac_sub[row.names(table_bac_sub) %in% tax_df_sub_2, ]
    out <- data.frame(rowSums(table_bac_sub_2)/length(colnames(table_bac_sub_2)))
    colnames(out) <- "RA"
    
    KO_iso_3$RA <- out$RA[match(row.names(KO_iso_3), row.names(out))]
    KO_iso_3 <- na.omit(KO_iso_3)
    KO_iso_3$RA_prop <- KO_iso_3$RA/sum(KO_iso_3$RA)
    KO_iso_3$KO_prop <- KO_iso_3$No_of_KOs/KO_iso_2_sum_val
    
    KO_iso_3$SynCom <- paste(syncom)
    KO_iso_3$Family <- paste(family)
    
    
    new_data_next <- rbind(new_data_next,KO_iso_3)
  }
}

table_cor_2 <- data.frame()

for (family in unique(new_data_next$Family)){
  new_data_next_2 <- new_data_next[new_data_next$Family == paste(family),]
  for (syncom in unique(new_data_next$SynCom)){
    new_data_next_3 <- new_data_next_2[new_data_next_2$SynCom == paste(syncom),]
    new_data_next_3 <- new_data_next_3[new_data_next_3$RA_prop != 0,]
    
    if (length(row.names(new_data_next_3)) < 6){
      pval <- NA
      cor <- NA
    } else {
      val <- lm(formula = KO_prop ~ log(RA_prop), data = new_data_next_3)
      
      sum <- summary(val)
      cor <- sum$adj.r.squared
      pval_sum <- sum$coefficients
      pval <- pval_sum[8]
      
    }
    
    table_cor <- data.frame(t(data.frame(c(paste(syncom), cor, pval, paste(family)))))
    table_cor_2 <- rbind(table_cor_2, table_cor)
  }
}

row.names(table_cor_2) <- NULL
colnames(table_cor_2) <- c("SynCom", "R2", "Pvalue", "Family")

table_cor_2$R2 <- as.numeric(table_cor_2$R2)
table_cor_2$Pvalue <- as.numeric(table_cor_2$Pvalue)

table_cor_2$Sig <- ""
table_cor_2$Sig[table_cor_2$Pvalue < 0.05] <- "*"

table_cor_2$SynCom <- gsub("Syncom", "SynCom", table_cor_2$SynCom)

Plot_fam <- ggplot(table_cor_2, aes(SynCom, Family)) +
  geom_tile(height=0.98, mapping = aes(fill = R2)) +
  geom_text(aes(label = Sig), size =4) +
  scale_fill_gradient2(low = "#D55e00", mid = "white", high = "#56b4e9", midpoint =0, na.value = "lightgrey")+
  theme_classic() +
  labs(x ="Inoculum", y = "Family", fill = "R2") +
  theme(panel.background=element_blank(),panel.grid=element_blank(),axis.line.x=element_line(size=.5, colour="black"),axis.line.y=element_line(size=.5, colour="black"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black", size=7),legend.position="right",legend.text= element_text(size=10),text=element_text(family="sans", size=10))+
  theme(axis.text.x = element_text(size = 14, angle = 25,hjust=1),axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18), axis.text.y = element_text(size=14, face = rep("italic")), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=18)) +
  ggtitle("Root colonization versus functional diversity") +
  theme(plot.title = element_text(hjust = 0.5))
Plot_fam

pdf(paste(results.dir,"Figure_S15_AtSC_7_200_member_SCs_heatmap.pdf", sep=""), width=8, height=6)
print(Plot_fam)
dev.off()

###Figure S15 & Table S6 - panel B Line diagram =====
Families <- c("Burkholderiaceae","Rhizobiaceae","Xanthomonadaceae","Pseudomonadaceae","Caulobacteraceae","Beijerinckiaceae")

tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
KO_iso <- read.table(paste(working_directory,"KO_genome/KO_SSC.tsv", sep = ""), row.names =1, header =T)
tax_df <- tax_df[tax_df$SynCom == "AtSC",]

colnames(KO_iso) <- gsub("X","", colnames(KO_iso))
colnames(KO_iso) <- gsub("M.11_2","M-11_2", colnames(KO_iso))
colnames(KO_iso) <- gsub("M.11","M-11", colnames(KO_iso))
colnames(KO_iso) <- gsub("M.16","M-16", colnames(KO_iso))
colnames(KO_iso) <- gsub("M.10","M-10", colnames(KO_iso))
colnames(KO_iso) <- gsub("M.6","M-6", colnames(KO_iso))

#Filter out samples with high contamination
plant_reads <- read.table(paste(working_directory, "AtSC_7_SynCom_experiment/plant_reads.tsv", sep = ""), header =T)  
SynCom_reads <- read.table(paste(working_directory, "AtSC_7_SynCom_experiment/mapped_reads.tsv", sep = ""), header =T)

plant_reads$SynCom_reads_2 <- SynCom_reads$Full_200[match(plant_reads$sample_id, SynCom_reads$Sample)]
plant_reads$SynCom_reads <- SynCom_reads$Map_200[match(plant_reads$sample_id, SynCom_reads$Sample)]

plant_reads$non_plant_pairs <- plant_reads$total_reads-plant_reads$plant_reads

plant_reads$other_isolate_reads <- plant_reads$SynCom_reads_2 - plant_reads$SynCom_reads
plant_reads$Contaminant_reads <- plant_reads$non_plant_pairs - plant_reads$SynCom_reads_2

Hank_the_normalizer <- function(df,group,amount){
  df_2 <- df %>% dplyr::group_by_at(group) %>% dplyr::summarise(total=sum(.data[[amount]]))
  df_3 <- df_2$total
  names(df_3) <- df_2[[group]]
  df$total <- df_3[as.character(df[[group]])]
  df$Rel <- df[[amount]] / df$total
  return(df)
}

plant_reads_2 <- plant_reads[,c(1,3,5,7,8)]

plant_reads_melt <- melt(plant_reads_2)
plant_reads_melt_2 <- Hank_the_normalizer(plant_reads_melt,"sample_id","value")

removal_sam <- c()

for (sample in unique(plant_reads_melt_2$sample_id)){
  plant_reads_melt_2_sub <- plant_reads_melt_2[plant_reads_melt_2$sample_id == paste(sample),]
  value <- plant_reads_melt_2_sub$Rel[plant_reads_melt_2_sub$variable == "SynCom_reads"]/sum(plant_reads_melt_2_sub$Rel[plant_reads_melt_2_sub$variable != "plant_reads"])
  if (value <= 0.6){
    removal_sam <- c(removal_sam, paste(sample))
  }
}

#Microbiome data
table_bac <- read.table(paste(working_directory,"AtSC_7_SynCom_experiment/isolate_norm.txt", sep = ""), row.names =1, header =T)
table_bac_1 <- t(t(table_bac) / rowSums(t(table_bac)))
table_bac_2 <- table_bac_1[,!colnames(table_bac_1) %in% removal_sam] 

#metadata
metadata <- read.table(paste(working_directory,"AtSC_7_SynCom_experiment/metadata.txt", sep= ""), row.names =1, header =T)
metadata_2 <- metadata[grep("Syncom",metadata$Syncom),]
SynComs <- unique(metadata_2$Syncom)

#isolates
groups <- read.table(paste(working_directory,"AtSC_7_SynCom_experiment/SynCom_isolates.txt", sep = ""), row.names =1, header =F)

new_datafr <- data.frame()

for (syncom in SynComs){
  #Observed
  metadata_sub <- metadata_2[metadata_2$Syncom == paste(syncom),]
  metadata_sub_2 <- row.names(metadata_sub)[metadata_sub$Compartment != "inoculum"]
  
  isolates_in_SynCom <- as.vector(unlist(as.vector(groups[row.names(groups) == paste(syncom),])))
  tax_df_sub <- tax_df[tax_df$isolate %in% isolates_in_SynCom,]
  
  table_bac_sub <- data.frame(table_bac_2[,colnames(table_bac_2) %in% metadata_sub_2])
  
  for (family in Families){
    tax_df_sub_2 <- tax_df_sub$isolate[tax_df_sub$family == paste(family)]
    
    #Expected
    KO_iso_2 <- KO_iso[,colnames(KO_iso) %in% tax_df_sub_2]
    KO_iso_2[KO_iso_2 > 0] <- 1
    
    KO_iso_3 <- data.frame(colSums(KO_iso_2))
    colnames(KO_iso_3) <- "No_of_KOs"
    KO_iso_4 <- KO_iso_3 %>% dplyr::arrange(desc(No_of_KOs))
    
    m <- Mclust(KO_iso_4$No_of_KOs)     
    KO_iso_4$group <- m$classification
    
    #Observed
    table_bac_sub_2 <- table_bac_sub[row.names(table_bac_sub) %in% tax_df_sub_2, ]
    out <- data.frame(rowSums(table_bac_sub_2)/length(colnames(table_bac_sub_2)))
    colnames(out) <- "RA"
    
    for (cluster in unique(KO_iso_4$group)){
      KO_iso_4_sub <- row.names(KO_iso_4)[KO_iso_4$group == paste(cluster)]
      out_sub <- sum(out[row.names(out) %in% KO_iso_4_sub,])
      KO_iso_sub <- KO_iso_4[KO_iso_4$group == paste(cluster),]
      value_KO <- sum(KO_iso_sub$No_of_KOs)/length(KO_iso_sub$No_of_KOs)
      new_datafr <- rbind(new_datafr, data.frame(t(data.frame(c(paste(syncom),paste(family),paste(cluster),value_KO, out_sub)))))
    }
  }
}

row.names(new_datafr) <- NULL
colnames(new_datafr) <- c("SynCom","Family","Cluster","No_of_KO","Cum_RA")

new_datafr$SynFam <- paste(new_datafr$SynCom, new_datafr$Family, sep = "_")

new_datafr_2 <- data.frame()

for (group in unique(new_datafr$SynFam)){
  new_datafr_sub <- new_datafr[new_datafr$SynFam == paste(group),]
  new_datafr_sub_2 <- new_datafr_sub[order(as.numeric(new_datafr_sub$Cum_RA), decreasing = F), ]
  new_datafr_sub_2$order_RA <- 1:length(new_datafr_sub_2$SynFam)
  new_datafr_2 <- rbind(new_datafr_2,new_datafr_sub_2 )
}

row.names(new_datafr_2) <- NULL
colnames(new_datafr_2) <- c("SynCom","Family","Exp_order","No_of_KO","Cum_RA","SynFam", "Obs_order")

stat_data <- data.frame()

for (group in unique(new_datafr_2$Family)){
  new_datafr_2_sub <- new_datafr_2[new_datafr_2$Family == paste(group),]
  check_table <- table(new_datafr_2_sub$SynFam)
  check_table_2 <- names(check_table)[check_table == 1]
  new_datafr_2_sub <- new_datafr_2_sub[!new_datafr_2_sub$SynFam %in% check_table_2,]
  
  if (length(new_datafr_2_sub$Family) > 2){
    stat_out <- cor.test(as.numeric(new_datafr_2_sub$Exp_order), as.numeric(new_datafr_2_sub$Obs_order), method = "kendall")
    pval <- stat_out$p.value
    stat_data <- rbind(stat_data, data.frame(t(data.frame(c(paste(group), pval)))))
  } else {
    pval <- NA
    SC <- unique(new_datafr_2$SynCom[new_datafr_2$SynFam == paste(group)])
    FAM <- unique(new_datafr_2$Family[new_datafr_2$SynFam == paste(group)])
    stat_data <- rbind(stat_data, data.frame(t(data.frame(c(paste(SC), paste(FAM), pval)))))
  }
}

row.names(stat_data) <- NULL
colnames(stat_data) <- c("Family","Kendall correlation p-value")

new_datafr_2$SynCom <- gsub("Syncom", "SynCom", new_datafr_2$SynCom)
new_datafr_2$SynCom[new_datafr_2$SynCom == "SynCom11"] <- "SynCom1"

plot <- ggplot(new_datafr_2, aes(x = Exp_order, y = Obs_order, color = as.factor(Family), group = 1)) +
  geom_point(size = 3) +
  #geom_line() +   
  theme_classic() +
  labs(x = "Functional diversity rank",
       y = "Relative abundance rank",
       color = "Family") +
  ggtitle("Functional diversity vs Relative abundance") +
  theme(axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 14),
        axis.text.y = element_text(size=10),
        legend.title = element_text(size=16),
        legend.text = element_text(size=12,  face = rep("italic")),
        plot.title = element_text(size=18)) +
  facet_grid(SynCom ~ Family ,scales='free')
plot <- plot +  geom_smooth(method = "lm", se = FALSE) + theme(strip.text.x = element_text(face=rep("italic")))
plot

pdf(paste(results.dir,"Figure_S15_Line_plot_Func_vs_RA.pdf", sep=""), width=12, height=9)
print(plot)
dev.off()

stat_data$`SynCom number` <- c("4/7 SynComs", "5/5 SynComs", "3/7 SynComs", "7/7 SynComs", "7/7 SynComs", "5/5 SynComs")

write.table(new_datafr_2,paste(results.dir,"Data_func_vs_RA_output.txt", sep = ""), col.names =T, row.names =F, quote =F, sep = "\t")
write.table(stat_data,paste(results.dir,"Table_S6_func_RA_stats.txt", sep = ""), col.names =T, row.names =F, quote =F, sep = "\t")

###Figure 4a - Family KO R2 effects =====
##Panel A subsetted no Burkho genera, and <0.05
#Contribution of the different bacterial families to the observed community-level split by inoculum
# Define hosts
Hosts <- c("At", "Hv", "Lj")
PL_colors <- c("Lj" = "#7570b3", "Hv" = "#d95f02",  "At"= "#1b9e77")

# Define SynComs
SynComs <- c("AtSC", "HvSC", "LjSC", "SSC")
SC_colors <- c("#F8766D", "#00BF7D", "#00B0F6", "#A3A500")


# Read data - With Dominators
fam_data <- read.table(paste0(working_directory, "Family_R2/SSC_Fam_R2_effects_with_dom.txt"),
                       sep="\t", header=TRUE, row.names=1)

combined_syncom_with_dom <- rbind(fam_data) %>%
  mutate(dataset = "Dominators")

# Read data - No Dominators
fam_data5 <- read.table(paste0(working_directory, "Family_R2/SSC_Fam_R2_effects_no_dom.txt"),
                        sep="\t", header=TRUE, row.names=1)


combined_syncom_no_dom <- rbind(fam_data5) %>%
  mutate(dataset = "No_Dominators")

# Merge datasets
pyramid_data_fam <- bind_rows(combined_syncom_no_dom, combined_syncom_with_dom)

family_order <- c("Chitinophagaceae", "Microbacteriaceae", "Micrococcaceae",
                  "Xanthobacteraceae", "Sphingobacteriaceae", "Rhodanobacteraceae",
                  "Sphingomonadaceae", "Flavobacteriaceae", "Devosiaceae",
                  "Beijerinckiaceae", "Enterobacteriaceae", "Caulobacteraceae",
                  "Pseudomonadaceae", "Xanthomonadaceae", "Burkholderiaceae",
                  "Rhizobiaceae")

# Filter for KO only and subset by hosts
pyramid_data_fam_host <- pyramid_data_fam %>% filter(Subset %in% Hosts)

# Factor levels
pyramid_data_fam_host$Subset <- factor(pyramid_data_fam_host$Subset, levels = c("Lj", "Hv", "At"))

pyramid_data_fam_host$Family <- factor(pyramid_data_fam_host$Family, levels = family_order)
pyramid_data_fam_host$R2_change <- abs(pyramid_data_fam_host$R2_change)
pyramid_data_fam_host$dataset <- factor(pyramid_data_fam_host$dataset, levels = c("Dominators", "No_Dominators"))

# Compute sum of R2_change per Family and dataset
sum_per_group_host <- aggregate(R2_change ~ Family + dataset, data = pyramid_data_fam_host, sum, na.rm = TRUE)

#Contribution of the different bacterial families to the observed community-level split by plant host

# Filter for KO only and subset by Inoculum
pyramid_data_fam_inoc <- pyramid_data_fam %>% filter(Subset %in% SynComs)

# Factor levels
pyramid_data_fam_inoc$Subset <- factor(pyramid_data_fam_inoc$Subset, levels = c("SSC", "LjSC", "HvSC", "AtSC"))
pyramid_data_fam_inoc$Family <- factor(pyramid_data_fam_inoc$Family, levels = family_order)
pyramid_data_fam_inoc$R2_change <- abs(pyramid_data_fam_inoc$R2_change)
pyramid_data_fam_inoc$dataset <- factor(pyramid_data_fam_inoc$dataset, levels = c("Dominators", "No_Dominators"))

# Compute sum of R2_change per Family and dataset
sum_per_group_inoc <- aggregate(R2_change ~ Family + dataset, data = pyramid_data_fam_inoc, sum, na.rm = TRUE)

# Get families with R2_change > 0.05 in either dataset
families_keep <- union(
  unique(sum_per_group_host$Family[sum_per_group_host$R2_change > 0.05]),
  unique(sum_per_group_inoc$Family[sum_per_group_inoc$R2_change > 0.05])
)

# Apply the same family filter to both datasets
pyramid_data_fam_host <- pyramid_data_fam_host %>%
  filter(Family %in% families_keep) %>%
  droplevels()

pyramid_data_fam_inoc <- pyramid_data_fam_inoc %>%
  filter(Family %in% families_keep) %>%
  droplevels()

# Define bar width
barwidth <- 0.3

# Calculate positions for Dominators and Non-Dominators
Dominators <- filter(pyramid_data_fam_host, dataset == "Dominators") %>%
  group_by(Family) %>%
  arrange(-Subset) 

No_Dominators <- filter(pyramid_data_fam_host, dataset == "No_Dominators") %>%
  group_by(Family) %>%
  arrange(-Subset) 

# Plot with the calculated positions
plot8 <- ggplot() +
  geom_bar(data = Dominators,
           mapping = aes(x = as.numeric(Family) + barwidth/2, y = R2_change, fill = Subset, alpha = "Dominators"),
           stat = "identity",
           position = 'stack',
           color = "black",
           size = 0.2,
           width = barwidth) +
  geom_bar(data = No_Dominators,
           mapping = aes(x = as.numeric(Family) - barwidth/2, y = R2_change, fill = Subset, color = Subset, alpha = "No_Dominators"),
           stat = "identity",
           position = 'stack',
           size = 0.2,
           width = barwidth) +
  ggtitle("Family R² Effects With vs. Without Dominators") +
  theme(plot.title = element_text(hjust = 0.5, size = 10)) +
  theme_classic() +
  labs(x = "Family", y = "Effect on Inoculum R²", fill = "Plant Host", alpha="Dataset") +
  scale_fill_manual(values = PL_colors) +
  scale_color_manual(values = PL_colors) +
  scale_alpha_manual(values = c("Dominators" = 1,"No_Dominators" = 0.3)) +
  scale_y_continuous(limits = c(0.32, 0), 
                     expand = expansion(mult = c(-0.05, 0)), 
                     trans = "reverse",
                     breaks = seq(0.3, 0, by = -0.1),  # Ensure 0.3 is included in breaks
                     labels = scales::label_number(accuracy = 0.01))+  # Format labels with two decimal places
  scale_x_continuous(breaks = 1:length(levels(pyramid_data_fam_host$Family)), labels = levels(pyramid_data_fam_host$Family), position = "top") + # Adjust the expand parameter here
  coord_flip() + # Flip the coordinates
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black", size = 7),
        legend.position = "right",
        legend.background = element_blank(),
        legend.key = element_blank(),
        text = element_text(family = "sans", size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(angle = 0, vjust = 0.5), # Adjust y-axis title position
        axis.text.y = element_text(face = "italic", size = 10, angle = 0, hjust = 0)) # Adjust y-axis text position

# Display the plot
print(plot8)

# Generate Dominators and Non-Dominators datasets
Dominators <- filter(pyramid_data_fam_inoc, dataset == "Dominators") %>%
  group_by(Family) %>%
  arrange(-Subset)

No_Dominators <- filter(pyramid_data_fam_inoc, dataset == "No_Dominators") %>%
  group_by(Family) %>%
  arrange(-Subset)

# Plot native
plot9 <- ggplot() +
  geom_bar(data = Dominators,
           mapping = aes(x = as.numeric(Family) + barwidth/2, y = R2_change, fill = Subset, alpha = "Dominators"),
           stat = "identity",
           position = 'stack',
           color = "black",
           size = 0.2,
           width = barwidth) +
  geom_bar(data = No_Dominators,
           mapping = aes(x = as.numeric(Family) - barwidth/2, y = R2_change, fill = Subset, alpha = "No_Dominators", color=Subset),
           stat = "identity",
           position = 'stack',
           size = 0.2,
           width = barwidth) +
  ggtitle("Family R² Effects With vs. Without Dominators") +
  theme(plot.title = element_text(hjust = 0.5, size = 10)) +
  theme_classic() +
  labs(x = "Family", y = "Effect on SynCom R²", fill = "Inoculum", alpha="Dataset") +
  scale_fill_manual(values = SC_colors) +
  scale_color_manual(values = SC_colors) +
  scale_alpha_manual(values = c("Dominators" = 1, "No_Dominators" = 0.3)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black", size = 7),
        legend.position = "right",
        legend.background = element_blank(),
        legend.key = element_blank(),
        text = element_text(family = "sans", size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = "italic", size = 10)) +
  scale_x_continuous(breaks = 1:length(levels(pyramid_data_fam_inoc$Family)), labels = levels(pyramid_data_fam_inoc$Family))+
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01), expand = expansion(mult = c(0, 0.05)), breaks = seq(0,0.3, 0.1)) +
  coord_flip()

# Display the plot
plot9

# Plot with Rhizobiaceae break
plot9_lim <- ggplot() +
  geom_bar(data = Dominators,
           mapping = aes(x = as.numeric(Family) + barwidth/2, y = R2_change, fill = Subset, alpha = "Dominators"),
           stat = "identity",
           position = 'stack',
           color = "black",
           size = 0.2,
           width = barwidth) +
  geom_bar(data = No_Dominators,
           mapping = aes(x = as.numeric(Family) - barwidth/2, y = R2_change, fill = Subset, alpha = "No_Dominators", color=Subset),
           stat = "identity",
           position = 'stack',
           size = 0.2,
           width = barwidth) +
  ggtitle("Family R² Effects With vs. Without Dominators") +
  theme(plot.title = element_text(hjust = 0.5, size = 10)) +
  theme_classic() +
  labs(x = "Family", y = "Effect on Host R²", fill = "Inoculum", alpha="Dataset") +
  scale_fill_manual(values = SC_colors) +
  scale_color_manual(values = SC_colors) +
  scale_alpha_manual(values = c("Dominators" = 1, "No_Dominators" = 0.3)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black", size = 7),
        legend.position = "right",
        legend.background = element_blank(),
        legend.key = element_blank(),
        text = element_text(family = "sans", size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = "italic", size = 10)) +
  scale_x_continuous(breaks = 1:length(levels(pyramid_data_fam_inoc$Family)), labels = levels(pyramid_data_fam_inoc$Family))+
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01), expand = expansion(mult = c(0, 0.05)), breaks = seq(0,0.3, 0.1), limits = c(0,0.32) ) +
  coord_flip()

plot9_lim




# Combine syncom and plant R² family effect (NATIVE)
pyramid_dom=ggarrange(print(plot8), print(plot9), 
                      ncol = 2, nrow = 1,    
                      common.legend = F,  
                      legend = "none")      
pyramid_dom

# Combine syncom and plant R² family effect with limits
pyramid_dom_lim=ggarrange(print(plot8), print(plot9_lim), 
                          ncol = 2, nrow = 1,    
                          common.legend = F,  
                          legend = "none")      
pyramid_dom_lim

# Extract legends from both plots
legend8 <- get_legend(plot8 + theme(legend.position = "right"))
legend9 <- get_legend(plot9 + theme(legend.position = "right"))

# Arrange legends in a single plot
legend_plot <- ggarrange(legend8, legend9, ncol = 2, nrow = 1)

# Display the legend plot
legend_plot

# # Save plot
pdf(paste0(results.dir, "Figure_4a_Pyramid_Plot_main_lim.pdf"), width = 10, height = 8)
print(pyramid_dom_lim)
dev.off()

pdf(paste0(results.dir, "Figure_4a_Pyramid_Plot_legends.pdf"), width = 10, height = 8)
print(legend_plot)
dev.off()


###Figure 4b - sPLS-DA dot plot =====
empty_vector_all_2_all_sub <- read.table(paste(working_directory, "sPLS-DA/data_table_sPLSDA_main.txt", sep =""), header =T, sep = "\t")

dot_plot <- ggplot(empty_vector_all_2_all_sub, aes(y=Family, x=Contribution, size=No_of_isolates)) + 
  geom_point(aes(fill = factor(Data)),shape=21, stroke =1) + 
  ggtitle("KO contribution") + 
  theme_classic() +
  scale_fill_manual(values = c("#1b9e77","#A3A500","#d95f02","#00B0F6","#7570b3","#00BF7D","#F8766D" ))+
  theme(plot.title = element_text(hjust = 0.5)) + theme(text = element_text(size=12)) + theme(legend.text=element_text(size=12)) + 
  labs(x = "Contribution",y = "Family", size = "No of isolates", fill = "Group") +
  geom_text_repel(aes(label=as.character(empty_vector_all_2_all_sub$KO_2)),size=4, max.overlaps = 8) +
  theme(panel.border = element_blank(),panel.grid.major = element_line(color = "gray70", size = 0.5),panel.grid.minor = element_line(color = "gray90", size = 0.5),panel.background = element_blank(),axis.line = element_line(colour = "black")) +
  facet_wrap(~Facet_group,scales ="free_y",nrow =4) +
  theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 12))
dot_plot

dot_plot_rotated <- ggplot(empty_vector_all_2_all_sub, aes(x=Family, y=Contribution, size=No_of_isolates)) + 
  geom_point(aes(fill = factor(Data)), shape=21, stroke=1) + 
  ggtitle("KO contribution") + 
  theme_classic() +
  scale_fill_manual(values = c("#1b9e77","#A3A500","#d95f02","#00B0F6","#7570b3","#00BF7D","#F8766D")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size=14)) + 
  theme(legend.text = element_text(size=14)) +
  labs(y = "Contribution", x = "Family", size = "No of isolates", fill = "Group") +
  geom_text_repel(aes(label = as.character(KO_2)), size = 6, max.overlaps = 4, angle = 0, direction = "both") +
  theme(panel.border = element_blank(),
        panel.grid.major = element_line(color = "gray70", size = 0.5),
        panel.grid.minor = element_line(color = "gray90", size = 0.5),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  facet_wrap(~Dataset+Subset, scales = "free_x", ncol = 4) +  # horizontal facets now
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14))
dot_plot_rotated
dot_plot_rotated+theme(legend.position = "none")

#Save plot
pdf(paste(results.dir, "Figure_4b_dotplot.pdf", sep =""), width = 21, height = 11, )
print(dot_plot_rotated+theme(legend.position = "none"))
dev.off()

# Extract legends 
legend10 <- get_legend(dot_plot_rotated )

#Save plot
pdf(paste(results.dir, "Figure_4b_dotplot legend.pdf", sep = ""), width = 10, height = 8)
print(legend10)
dev.off()


###Figure 4c and S16 boxplots file preparation =====

# First Load sPLS-DA results
empty_vector_all_2_rhizo <- read.table(paste(working_directory, "sPLS-DA/output/PLSDA_Rhizobiaceae_with_dom.tsv", sep=""), sep = "\t",  header = T)
empty_vector_all_2_rhizo$Group <- "Rhizobiaceae_dominators_HvSC_LjSC_SSC"
empty_vector_all_2_Burk_SC <- read.table(paste(working_directory, "sPLS-DA/output/PLSDA_Burkholderiaceae_no_dom_SynCom.tsv", sep=""), sep = "\t",  header = T)
empty_vector_all_2_Burk_SC$Group <- "Burkholderiaceae_no_dominators"
empty_vector_all_2_Xant_SC <- read.table(paste(working_directory, "sPLS-DA/output/PLSDA_Xanthomonadaceae_no_dom_SynCom.tsv", sep=""), sep = "\t",  header = T)
empty_vector_all_2_Xant_SC$Group <- "Xanthomonadaceae_no_dominators"
empty_vector_all_2_Xant_PL <- read.table(paste(working_directory, "sPLS-DA/output/PLSDA_Xanthomonadaceae_no_dom_Plant.tsv", sep=""), sep = "\t",  header = T)
empty_vector_all_2_Xant_PL$Group <- "Xanthomonadaceae_no_dominators_LjSC"
empty_vector_all_2_Pseud <- read.table(paste(working_directory, "sPLS-DA/output/PLSDA_Pseudomonadaceae_no_dom_Plant.tsv", sep=""), sep = "\t",  header = T)
empty_vector_all_2_Pseud$Group <- "Pseudomonadacaea_no_dominators_LjSC"
empty_vector_all_2_Caulo <- read.table(paste(working_directory, "sPLS-DA/output/PLSDA_Caulobacteraceae_no_dom_Plant.tsv", sep=""), sep = "\t",  header = T)
empty_vector_all_2_Caulo$Group <- "Caulobacteraceae_no_dominators_AtSC"
empty_vector_all_2_Burk_PL <- read.table(paste(working_directory, "sPLS-DA/output/PLSDA_Burkholderiaceae_no_dom_Plant.tsv", sep=""), sep = "\t",  header = T)
empty_vector_all_2_Burk_PL$Group <- "Burkholderiaceae_no_dominators_AtSC"

empty_vector_all_2_all <- rbind(empty_vector_all_2_rhizo,empty_vector_all_2_Burk_SC,empty_vector_all_2_Xant_SC,empty_vector_all_2_Xant_PL,empty_vector_all_2_Pseud,empty_vector_all_2_Caulo,empty_vector_all_2_Burk_PL)
empty_vector_all_2_all$No_of_isolates <- NA
KO_SSC =read.table(paste(working_directory,"KO_genome/KO_SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

for (KO in unique(empty_vector_all_2_all$KO)){
  KO_SSC_sub <- KO_SSC[row.names(KO_SSC) == paste(KO),]
  KO_SSC_sub[KO_SSC_sub > 0] <- 1
  group <- sum(KO_SSC_sub)
  empty_vector_all_2_all$No_of_isolates[empty_vector_all_2_all$KO == paste(KO)] <- group
}

KOS <- table(empty_vector_all_2_all$KO)
KO_groups <- names(KOS)[KOS > 1]
empty_vector_all_2_all$Contribution <- (empty_vector_all_2_all$contrib * empty_vector_all_2_all$Variance) * 100
empty_vector_all_2_all <- empty_vector_all_2_all[empty_vector_all_2_all$Contribution >= 10,]

empty_vector_all_2_all$KO_2 <-  empty_vector_all_2_all$KO
empty_vector_all_2_all_sub <- empty_vector_all_2_all[!empty_vector_all_2_all$KO %in% KO_groups,]

for (KO in KO_groups){
  empty_vector_all_2_all_subber <- empty_vector_all_2_all[empty_vector_all_2_all$KO == paste(KO),]
  if(length(unique(empty_vector_all_2_all_subber$Data)) == 1 & length(unique(empty_vector_all_2_all_subber$Group)) == 1){
    empty_vector_all_2_all_subber_2 <- empty_vector_all_2_all_subber[1,]
    empty_vector_all_2_all_subber_2$contrib <- sum(empty_vector_all_2_all_subber$contrib)/length(empty_vector_all_2_all_subber$contrib)
    empty_vector_all_2_all_subber_2$Variance <- sum(empty_vector_all_2_all_subber$Variance)/length(empty_vector_all_2_all_subber$Variance)
    empty_vector_all_2_all_subber_2$Contribution <- sum(empty_vector_all_2_all_subber$Contribution)/length(empty_vector_all_2_all_subber$Contribution)
    empty_vector_all_2_all_sub <- rbind(empty_vector_all_2_all_sub,empty_vector_all_2_all_subber_2)
  } else if (length(unique(empty_vector_all_2_all_subber$Data)) == 1 & length(unique(empty_vector_all_2_all_subber$Group)) != 1){
    empty_vector_all_2_all_subber$KO <- paste(empty_vector_all_2_all_subber$KO, empty_vector_all_2_all_subber$Group, sep = "_")
    empty_vector_all_2_all_sub <- rbind(empty_vector_all_2_all_sub,empty_vector_all_2_all_subber)
  } else {
    empty_vector_all_2_all_subber$KO <- paste(empty_vector_all_2_all_subber$KO, empty_vector_all_2_all_subber$Data, sep = "_")
    empty_vector_all_2_all_sub <- rbind(empty_vector_all_2_all_sub,empty_vector_all_2_all_subber)
  }
}

###Figure 4c boxplots =====

#Rhizobiaceae with dominators
# Load the KO abundance tables (Rhizobiaceae-inclusive and exclusive datasets)
KO_SSC_only <- read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Rhizobiaceae_KO_with_dom.tsv", sep = ""),
                          header = TRUE, sep = "\t", row.names = 1)
KO_SSC_ex_only <- read.table(paste(working_directory, "sPLS-DA/isolate_drop_out_data/No_rhizobiaceae_KO_with_dom.tsv", sep = ""),
                             header = TRUE, sep = "\t", row.names = 1)

# Load metadata and filter for relevant samples
samples_df <- read.table(paste(working_directory, "SSC_R2_metadata_no_HL.tsv", sep = ""), header = TRUE, sep = "\t", row.names = 1)
colnames(samples_df)[5] <- "Nutrient"

samples_df_sub <- subset(samples_df, Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, Inoculum != "NS")
samples_df_sub_3 <- subset(samples_df_sub_2, Inoculum != "AtSC")

# Normalize KO tables by row sums (relative abundance)
KO_SSC_only_2 <- t(t(KO_SSC_only) / rowSums(t(KO_SSC_only)))
KO_SSC_ex_only_2 <- t(t(KO_SSC_ex_only) / rowSums(t(KO_SSC_ex_only)))

# Extract KOs associated with Rhizobiaceae and dominant families
KOs <- empty_vector_all_2_all_sub[empty_vector_all_2_all_sub$Group == "Rhizobiaceae_dominators_HvSC_LjSC_SSC", ]
KOs_Lj <- KOs$KO[KOs$Data == "Lj"]
KOs_Hv <- KOs$KO[KOs$Data == "Hv"]

# Subset for barley KOs
KO_SSC_only_2_Hv <- data.frame(KO_SSC_only_2[row.names(KO_SSC_only_2) %in% KOs_Hv, ])
colnames(KO_SSC_only_2_Hv) <- "KO"
KO_SSC_ex_only_2_Hv <- data.frame(KO_SSC_ex_only_2[row.names(KO_SSC_ex_only_2) %in% KOs_Hv, ])
colnames(KO_SSC_ex_only_2_Hv) <- "KO"

KO_SSC_only_2_Hv$Dataset <- "Rhizobiaceae"
KO_SSC_ex_only_2_Hv$Dataset <- "other families"

# Subset for Lotus KOs and sum across rows
KO_SSC_only_2_Lj <- KO_SSC_only_2[row.names(KO_SSC_only_2) %in% KOs_Lj, ]
KO_SSC_ex_only_2_Lj <- KO_SSC_ex_only_2[row.names(KO_SSC_ex_only_2) %in% KOs_Lj, ]

KO_SSC_only_2_Lj_2 <- data.frame(colSums(KO_SSC_only_2_Lj))
KO_SSC_ex_only_2_Lj_2 <- data.frame(colSums(KO_SSC_ex_only_2_Lj))

colnames(KO_SSC_only_2_Lj_2) <- "KO"
colnames(KO_SSC_ex_only_2_Lj_2) <- "KO"

KO_SSC_only_2_Lj_2$Dataset <- "Rhizobiaceae"
KO_SSC_ex_only_2_Lj_2$Dataset <- "Other families"

# Assign contributing plant information
KO_SSC_only_2_Hv$Contributing_plant <- "Barley KOs"
KO_SSC_ex_only_2_Hv$Contributing_plant <- "Barley KOs"
KO_SSC_only_2_Lj_2$Contributing_plant <- "Lotus KOs"
KO_SSC_ex_only_2_Lj_2$Contributing_plant <- "Lotus KOs"

# Assign sample names
KO_SSC_only_2_Hv$Sample <- row.names(KO_SSC_only_2_Hv)
KO_SSC_ex_only_2_Hv$Sample <- row.names(KO_SSC_ex_only_2_Hv)
KO_SSC_only_2_Lj_2$Sample <- row.names(KO_SSC_only_2_Lj_2)
KO_SSC_ex_only_2_Lj_2$Sample <- row.names(KO_SSC_ex_only_2_Lj_2)

# Function to subset and add plant information
subset_and_add_plant <- function(df) {
  df <- df[df$Sample %in% row.names(samples_df_sub_3), ]
  df$Plant <- samples_df_sub_3$Condition[match(df$Sample, row.names(samples_df_sub_3))]
  df$Plant <- as.factor(df$Plant)
  return(df)
}

# Apply the function to each data frame
KO_SSC_only_2_Hv <- subset_and_add_plant(KO_SSC_only_2_Hv)
KO_SSC_ex_only_2_Hv <- subset_and_add_plant(KO_SSC_ex_only_2_Hv)
KO_SSC_only_2_Lj_2 <- subset_and_add_plant(KO_SSC_only_2_Lj_2)
KO_SSC_ex_only_2_Lj_2 <- subset_and_add_plant(KO_SSC_ex_only_2_Lj_2)

# Function to perform ANOVA and Tukey test, and plot the results
anova_and_plot <- function(df, title) {
  # Perform ANOVA
  fitAnova <- aov(KO ~ Plant, data = df)
  
  # Perform Tukey's post-hoc test
  Tukey <- TukeyHSD(fitAnova)
  
  # Get letters for significance
  letters_anova <- multcompLetters4(fitAnova, Tukey)$Plant$Letters
  
  # Combine results into a data frame for plotting
  ltlbl_combined <- data.frame(
    Plant = names(letters_anova),
    Letters = letters_anova
  )
  
  # Ensure the order of the factor levels is correct
  ltlbl_combined$Plant <- factor(ltlbl_combined$Plant, levels = unique(df$Plant))
  
  # Apply the manually defined order to the factor
  ltlbl_combined <- ltlbl_combined[order(ltlbl_combined$Plant), ]
  
  # Merge the letters with the data frame
  df <- merge(df, ltlbl_combined, by = "Plant")
  
  # Plot the results
  plot <- ggplot(df, aes(x = Plant, y = KO, colour = Plant)) +
    geom_boxplot(outlier.shape = NA) +
    theme_classic() +
    ylab("Cumulative RA of KOs") +
    geom_jitter(shape = 16, position = position_jitter(0.2), aes(colour = Plant), show.legend = TRUE) +
    theme(
      axis.text.x = element_blank(),
      axis.title.y = element_text(size = 14),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14),
      plot.title = element_text(size = 14),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      strip.text.x = element_text(size = 10)
    ) +
    scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
    guides(colour = FALSE) +
    facet_wrap(~Contributing_plant + Dataset, scales = "free") +
    scale_y_continuous(labels = scales::scientific,  ) +
    stat_summary(aes(label = Letters, y = max(KO)*1.03), fun = max, geom = "text")
  
  return(plot)
}

# Generate plots for each data frame
plot_Hv_Rhiz <- anova_and_plot(KO_SSC_only_2_Hv, "Barley KOs - Rhizobiaceae")
plot_Hv_other <- anova_and_plot(KO_SSC_ex_only_2_Hv, "Barley KOs - Other Families")
plot_Lj_Rhiz <- anova_and_plot(KO_SSC_only_2_Lj_2, "Lotus KOs - Rhizobiaceae")
plot_Lj_other <- anova_and_plot(KO_SSC_ex_only_2_Lj_2, "Lotus KOs - Other Families")

# Arrange plots using ggarrange
rhizo_plot=grid.arrange(plot_Hv_Rhiz, plot_Hv_other, plot_Lj_Rhiz, plot_Lj_other, nrow = 2)

### Panel C-2 Xanthomonadaceae without dominators

# Load the KO abundance tables (Xanthomonadaceae-inclusive and exclusive datasets)
KO_SSC_only=read.table(paste(working_directory,"sPLS-DA/isolate_subset_data/Xanthomonadaceae_KO_no_dom.tsv",sep =""), header=TRUE,sep="\t", row.names = 1)
KO_SSC_ex_only=read.table(paste(working_directory,"sPLS-DA/isolate_drop_out_data/No_xanthomonadaceae_KO_no_dom.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

# Load metadata and filter for relevant samples
samples_df <- read.table(paste(working_directory, "SSC_R2_metadata_no_HL.tsv", sep = ""), header = TRUE, sep = "\t", row.names = 1)
colnames(samples_df)[5]="Nutrient"
samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

# Normalize KO tables by row sums (relative abundance)
KO_SSC_only_2 <- t(t(KO_SSC_only)/rowSums(t(KO_SSC_only)))
KO_SSC_ex_only_2 <- t(t(KO_SSC_ex_only)/rowSums(t(KO_SSC_ex_only)))

# Extract KOs associated with Xanthomonadaceae and without dominators families
KOs <- empty_vector_all_2_all_sub[empty_vector_all_2_all_sub$Group == "Xanthomonadaceae_no_dominators", ]
KOs_LjSC <- KOs$KO[KOs$Data == "LjSC"]
KOs_HvSC <- KOs$KO[KOs$Data == "HvSC"]

#Subset for HvSC related KOs
KO_SSC_only_2_HvSC <- data.frame(KO_SSC_only_2[row.names(KO_SSC_only_2) %in% KOs_HvSC,])
colnames(KO_SSC_only_2_HvSC) <- "KO"
KO_SSC_ex_only_2_HvSC <- data.frame(KO_SSC_ex_only_2[row.names(KO_SSC_ex_only_2) %in% KOs_HvSC,])
colnames(KO_SSC_ex_only_2_HvSC) <- "KO"

KO_SSC_ex_only_2_HvSC$Dataset <- "other families"
KO_SSC_only_2_HvSC$Dataset <- "Xanthomonadaceae"

#Subset for LjSC related KOs
KO_SSC_only_2_LjSC <- KO_SSC_only_2[row.names(KO_SSC_only_2) %in% KOs_LjSC,]
KO_SSC_ex_only_2_LjSC <- KO_SSC_ex_only_2[row.names(KO_SSC_ex_only_2) %in% KOs_LjSC,]

# sum across rows
KO_SSC_only_2_LjSC_2 <- data.frame(colSums(KO_SSC_only_2_LjSC))
KO_SSC_ex_only_2_LjSC_2 <- data.frame(colSums(KO_SSC_ex_only_2_LjSC))

colnames(KO_SSC_only_2_LjSC_2) <- "KO"
colnames(KO_SSC_ex_only_2_LjSC_2) <- "KO"

KO_SSC_only_2_LjSC_2$Dataset <- "Xanthomonadaceae"
KO_SSC_ex_only_2_LjSC_2$Dataset <- "other families"

# Assign contributing inoculum information
KO_SSC_only_2_HvSC$Contributing_plant <- "HvSC KOs"
KO_SSC_ex_only_2_HvSC$Contributing_plant <- "HvSC KOs"
KO_SSC_only_2_LjSC_2$Contributing_plant <- "LjSC KOs"
KO_SSC_ex_only_2_LjSC_2$Contributing_plant <- "LjSC KOs"

# Assign sample names
KO_SSC_only_2_HvSC$Sample <- row.names(KO_SSC_only_2_HvSC)
KO_SSC_ex_only_2_HvSC$Sample <- row.names(KO_SSC_ex_only_2_HvSC)
KO_SSC_only_2_LjSC_2$Sample <- row.names(KO_SSC_only_2_LjSC_2)
KO_SSC_ex_only_2_LjSC_2$Sample <- row.names(KO_SSC_ex_only_2_LjSC_2)

# Function to subset and add plant information
subset_and_add_inoc <- function(df) {
  df <- df[df$Sample %in% row.names(samples_df_sub_2), ]
  df$Plant <- samples_df_sub_2$Inoculum[match(df$Sample, row.names(samples_df_sub_2))]
  df$Plant <- as.factor(df$Plant)
  return(df)
}

# Apply the function to each data frame
KO_SSC_only_2_HvSC <- subset_and_add_inoc(KO_SSC_only_2_HvSC)
KO_SSC_only_2_LjSC_2 <- subset_and_add_inoc(KO_SSC_only_2_LjSC_2)
KO_SSC_ex_only_2_HvSC <- subset_and_add_inoc(KO_SSC_ex_only_2_HvSC)
KO_SSC_ex_only_2_LjSC_2 <- subset_and_add_inoc(KO_SSC_ex_only_2_LjSC_2)

# Function to perform ANOVA and Tukey test, and plot the results
anova_and_plot <- function(df, title) {
  # Perform ANOVA
  fitAnova <- aov(KO ~ Plant, data = df)
  
  # Perform Tukey's post-hoc test
  Tukey <- TukeyHSD(fitAnova)
  
  # Get letters for significance
  letters_anova <- multcompLetters4(fitAnova, Tukey)$Plant$Letters
  
  # Combine results into a data frame for plotting
  ltlbl_combined <- data.frame(
    Plant = names(letters_anova),
    Letters = letters_anova
  )
  
  # Ensure the order of the factor levels is correct
  ltlbl_combined$Plant <- factor(ltlbl_combined$Plant, levels = unique(df$Plant))
  
  # Apply the manually defined order to the factor
  ltlbl_combined <- ltlbl_combined[order(ltlbl_combined$Plant), ]
  
  # Merge the letters with the data frame
  df <- merge(df, ltlbl_combined, by = "Plant")
  
  # Plot the results
  plot <- ggplot(df, aes(x = Plant, y = KO, colour = Plant)) +
    geom_boxplot(outlier.shape = NA) +
    theme_classic() +
    ylab("Cumulative RA of KOs") +
    geom_jitter(shape = 16, position = position_jitter(0.2), aes(colour = Plant), show.legend = TRUE) +
    theme(
      axis.text.x = element_blank(),
      axis.title.y = element_text(size = 14),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14),
      plot.title = element_text(size = 14),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      strip.text.x = element_text(size = 10)
    ) +
    scale_color_manual(values = c("#A3A500","#00B0F6","#00BF7D","#F8766D" )) +
    guides(colour = FALSE) +
    facet_wrap(~Contributing_plant + Dataset, scales = "free") +
    scale_y_continuous(labels = scales::scientific,  ) +
    stat_summary(aes(label = Letters, y = max(KO)*1.03), fun = max, geom = "text")
  
  return(plot)
}

# Generate plots for each data frame
plot_HvSC <- anova_and_plot(KO_SSC_only_2_HvSC, "HvSC KOs - Xanthomonadaceae")
plot_HvSC_other <- anova_and_plot(KO_SSC_ex_only_2_HvSC, "HvSC KOs - Other Families")
plot_LjSC <- anova_and_plot(KO_SSC_only_2_LjSC_2, "Lotus KOs - Rhizobiaceae")
plot_LjSC_other <- anova_and_plot(KO_SSC_ex_only_2_LjSC_2, "Lotus KOs - Other Families")

# Arrange plots using ggarrange
grid.arrange(plot_HvSC, plot_HvSC_other, plot_LjSC, plot_LjSC_other, nrow = 2)

### Panel C-3 Pseudomonadaceae LjSC without dominators (not follwing expected trend)

KO_SSC_only=read.table(paste(working_directory,"sPLS-DA/isolate_subset_data/Pseudomonadaceae_KO_no_dom.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
KO_SSC_ex_only=read.table(paste(working_directory,"sPLS-DA/isolate_drop_out_data/No_pseudomonadaceae_KO_no_dom.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Samples TABLE
samples_df <- read.table(paste(working_directory, "SSC_R2_metadata_no_HL.tsv", sep = ""), header = TRUE, sep = "\t", row.names = 1)
colnames(samples_df)[5]="Nutrient"
samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")
samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum == "LjSC")

KO_SSC_only_2 <- t(t(KO_SSC_only)/rowSums(t(KO_SSC_only)))
KO_SSC_ex_only_2 <- t(t(KO_SSC_ex_only)/rowSums(t(KO_SSC_ex_only)))

KOs <- empty_vector_all_2_all_sub[empty_vector_all_2_all_sub$Group == "Pseudomonadacaea no dominators LjSC", ]
KOs_Lj <- KOs$KO[KOs$Data == "Lj"]
KOs_At <- KOs$KO[KOs$Data == "At"]

#Arabidopsis 
KO_SSC_only_2_At <- data.frame(KO_SSC_only_2[row.names(KO_SSC_only_2) %in% KOs_At,])
colnames(KO_SSC_only_2_At) <- "KO"
KO_SSC_ex_only_2_At <- data.frame(KO_SSC_ex_only_2[row.names(KO_SSC_ex_only_2) %in% KOs_At,])
colnames(KO_SSC_ex_only_2_At) <- "KO"

KO_SSC_ex_only_2_At$Dataset <- "Other families LjSC"
KO_SSC_only_2_At$Dataset <- "Pseudomonadaceae LjSC"

#Lotus
KO_SSC_only_2_Lj <- KO_SSC_only_2[row.names(KO_SSC_only_2) %in% KOs_Lj,]
KO_SSC_ex_only_2_Lj <- KO_SSC_ex_only_2[row.names(KO_SSC_ex_only_2) %in% KOs_Lj,]

KO_SSC_only_2_Lj_2 <- data.frame(colSums(KO_SSC_only_2_Lj))
KO_SSC_ex_only_2_Lj_2 <- data.frame(colSums(KO_SSC_ex_only_2_Lj))

colnames(KO_SSC_only_2_Lj_2) <- "KO"
colnames(KO_SSC_ex_only_2_Lj_2) <- "KO"

KO_SSC_only_2_Lj_2$Dataset <- "Pseudomonadaceae LjSC"
KO_SSC_ex_only_2_Lj_2$Dataset <- "Other families LjSC"

KO_SSC_only_2_At$Contributing_plant <- "Arabidopsis KOs"
KO_SSC_ex_only_2_At$Contributing_plant <- "Arabidopsis KOs"
KO_SSC_only_2_Lj_2$Contributing_plant <- "Lotus KOs"
KO_SSC_ex_only_2_Lj_2$Contributing_plant <- "Lotus KOs"

KO_SSC_only_2_At$Sample <- row.names(KO_SSC_only_2_At)
KO_SSC_ex_only_2_At$Sample <- row.names(KO_SSC_ex_only_2_At)
KO_SSC_only_2_Lj_2$Sample <- row.names(KO_SSC_only_2_Lj_2)
KO_SSC_ex_only_2_Lj_2$Sample <- row.names(KO_SSC_ex_only_2_Lj_2)

# Function to subset and add plant information
subset_and_add_plant <- function(df) {
  df <- df[df$Sample %in% row.names(samples_df_sub_3), ]
  df$Plant <- samples_df_sub_3$Condition[match(df$Sample, row.names(samples_df_sub_3))]
  df$Plant <- as.factor(df$Plant)
  return(df)
}

# Apply the function to each data frame
KO_SSC_only_2_At <- subset_and_add_plant(KO_SSC_only_2_At)
KO_SSC_ex_only_2_At <- subset_and_add_plant(KO_SSC_ex_only_2_At)
KO_SSC_only_2_Lj_2 <- subset_and_add_plant(KO_SSC_only_2_Lj_2)
KO_SSC_ex_only_2_Lj_2 <- subset_and_add_plant(KO_SSC_ex_only_2_Lj_2)

# Function to perform ANOVA and Tukey test, and plot the results
anova_and_plot <- function(df, title) {
  # Perform ANOVA
  fitAnova <- aov(KO ~ Plant, data = df)
  
  # Perform Tukey's post-hoc test
  Tukey <- TukeyHSD(fitAnova)
  
  # Get letters for significance
  letters_anova <- multcompLetters4(fitAnova, Tukey)$Plant$Letters
  
  # Combine results into a data frame for plotting
  ltlbl_combined <- data.frame(
    Plant = names(letters_anova),
    Letters = letters_anova
  )
  
  # Ensure the order of the factor levels is correct
  ltlbl_combined$Plant <- factor(ltlbl_combined$Plant, levels = unique(df$Plant))
  
  # Apply the manually defined order to the factor
  ltlbl_combined <- ltlbl_combined[order(ltlbl_combined$Plant), ]
  
  # Merge the letters with the data frame
  df <- merge(df, ltlbl_combined, by = "Plant")
  
  # Plot the results
  plot <- ggplot(df, aes(x = Plant, y = KO, colour = Plant)) +
    geom_boxplot(outlier.shape = NA) +
    theme_classic() +
    ylab("Cumulative RA of KOs") +
    geom_jitter(shape = 16, position = position_jitter(0.2), aes(colour = Plant), show.legend = TRUE) +
    theme(
      axis.text.x = element_blank(),
      axis.title.y = element_text(size = 14),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14),
      plot.title = element_text(size = 14),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      strip.text.x = element_text(size = 10)
    ) +
    scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
    guides(colour = FALSE) +
    facet_wrap(~Contributing_plant + Dataset, scales = "free") +
    scale_y_continuous(labels = scales::scientific,  ) +
    stat_summary(aes(label = Letters, y = max(KO)*1.03), fun = max, geom = "text")
  
  return(plot)
}

# Generate plots for each data frame
plot_At_Pseudo <- anova_and_plot(KO_SSC_only_2_At, "Barley KOs - Pseudomonadaceae")
plot_At_other <- anova_and_plot(KO_SSC_ex_only_2_At, "Barley KOs - Other Families")
plot_Lj_Pseudo <- anova_and_plot(KO_SSC_only_2_Lj_2, "Lotus KOs - Pseudomonadaceae")
plot_Lj_2_other <- anova_and_plot(KO_SSC_ex_only_2_Lj_2, "Lotus KOs - Other Families")

# Arrange plots using ggarrange
grid.arrange(plot_At_Pseudo, plot_At_other, plot_Lj_Pseudo, plot_Lj_2_other, nrow = 2)

### Combine the 3 panels
combined_plot_2 <- 
  (plot_Hv_Rhiz + plot_Hv_other + plot_Lj_Rhiz + plot_Lj_other) / 
  (plot_HvSC + plot_HvSC_other + plot_LjSC + plot_LjSC_other) / 
  (plot_At_Pseudo + plot_At_other + plot_Lj_Pseudo + plot_Lj_2_other) +
  plot_layout(ncol = 3)

print(combined_plot_2)

#Save plot
pdf(paste(results.dir, "Figure_4c_boxplots.pdf", sep = ""), width = 21, height = 8)
print(combined_plot_2)
dev.off()

###Figure S16a - Family KO R2 effects - full =====
# Define hosts
Hosts <- c("At", "Hv", "Lj")
PL_colors <- c("Lj" = "#7570b3", "Hv" = "#d95f02",  "At"= "#1b9e77")

# Read data - With Dominators
fam_data <- read.table(paste0(working_directory, "Family_R2/SSC_Fam_R2_effects_with_dom.txt"),
                       sep="\t", header=TRUE, row.names=1)

fam_data2 <- read.table(paste0(working_directory, "Family_R2/SSC_Gen_Burk_R2_effects_with_dom.txt"),
                        sep="\t", header=TRUE, row.names=1)

combined_syncom_with_dom <- rbind(fam_data, fam_data2) %>%
  mutate(dataset = "Dominators")

# Read data - No Dominators
fam_data5 <- read.table(paste0(working_directory, "Family_R2/SSC_Fam_R2_effects_no_dom.txt"),
                        sep="\t", header=TRUE, row.names=1)

fam_data6 <- read.table(paste0(working_directory, "Family_R2/SSC_Gen_Burk_R2_effects_no_dom.txt"),
                        sep="\t", header=TRUE, row.names=1)

combined_syncom_no_dom <- rbind(fam_data5, fam_data6) %>%
  mutate(dataset = "No_Dominators")

# Merge datasets
pyramid_data_fam <- bind_rows(combined_syncom_no_dom, combined_syncom_with_dom)

# Filter for KO only
pyramid_data_fam <- pyramid_data_fam %>% filter(KO == "KO", Subset %in% Hosts)

# Factor levels
pyramid_data_fam$Subset <- factor(pyramid_data_fam$Subset, levels = c("Lj", "Hv", "At"))
family_order <- c("Pelomonas", "Cupriavidus", "Polaromonas",
                  "Variovorax", "Rhizobacter", "Acidovorax", "Chitinophagaceae", "Microbacteriaceae", "Micrococcaceae",
                  "Xanthobacteraceae", "Sphingobacteriaceae", "Rhodanobacteraceae",
                  "Sphingomonadaceae", "Flavobacteriaceae", "Devosiaceae",
                  "Beijerinckiaceae", "Enterobacteriaceae", "Caulobacteraceae",
                  "Pseudomonadaceae", "Xanthomonadaceae", "Burkholderiaceae",
                  "Rhizobiaceae")

pyramid_data_fam$Family <- factor(pyramid_data_fam$Family, levels = family_order)
pyramid_data_fam$R2_change <- abs(pyramid_data_fam$R2_change)
pyramid_data_fam$dataset <- factor(pyramid_data_fam$dataset, levels = c("Dominators", "No_Dominators"))
# Define bar width
barwidth <- 0.4

# Calculate positions for Dominators and Non-Dominators
Dominators <- filter(pyramid_data_fam, dataset == "Dominators") %>%
  group_by(Family) %>%
  arrange(-Subset) 

No_Dominators <- filter(pyramid_data_fam, dataset == "No_Dominators") %>%
  group_by(Family) %>%
  arrange(-Subset) 

# Plot with the calculated positions
plot8 <- ggplot() +
  geom_bar(data = Dominators,
           mapping = aes(x = as.numeric(Family) + barwidth/2, y = R2_change, fill = Subset, alpha = "Dominators"),
           stat = "identity",
           position = 'stack',
           color = "black",
           size = 0.2,
           width = barwidth) +
  geom_bar(data = No_Dominators,
           mapping = aes(x = as.numeric(Family) - barwidth/2, y = R2_change, fill = Subset, color = Subset, alpha = "No_Dominators"),
           stat = "identity",
           position = 'stack',
           size = 0.2,
           width = barwidth) +
  ggtitle("Family R² Effects With vs. Without Dominators") +
  theme(plot.title = element_text(hjust = 0.5, size = 10)) +
  theme_classic() +
  labs(x = "Family", y = "Effect on SynCom R²", fill = "Plant Host", alpha="Dataset") +
  scale_fill_manual(values = PL_colors) +
  scale_color_manual(values = PL_colors) +
  scale_alpha_manual(values = c("Dominators" = 1,"No_Dominators" = 0.3)) +
  scale_y_continuous(limits = c(0.32, 0), 
                     expand = expansion(mult = c(-0.05, 0)), 
                     trans = "reverse",
                     breaks = seq(0.3, 0, by = -0.1),  # Ensure 0.3 is included in breaks
                     labels = scales::label_number(accuracy = 0.01))+  # Format labels with two decimal places
  scale_x_continuous(breaks = 1:length(levels(pyramid_data_fam$Family)), labels = levels(pyramid_data_fam$Family), position = "top") + # Adjust the expand parameter here
  coord_flip() + # Flip the coordinates
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black", size = 7),
        legend.position = "right",
        legend.background = element_blank(),
        legend.key = element_blank(),
        text = element_text(family = "sans", size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(angle = 0, vjust = 0.5), # Adjust y-axis title position
        axis.text.y = element_text(face = "italic", size = 10, angle = 0, hjust = 0)) # Adjust y-axis text position

# Display the plot
print(plot8)


# Define SynComs
SynComs <- c("AtSC", "HvSC", "LjSC", "SSC")
SC_colors <- c("#F8766D", "#00BF7D", "#00B0F6", "#A3A500")

# Read data - With Dominators
fam_data3 <- read.table(paste0(working_directory, "Family_R2/SSC_Fam_R2_effects_with_dom.txt"),
                        sep="\t", header=TRUE, row.names=1)

fam_data4 <- read.table(paste0(working_directory, "Family_R2/SSC_Gen_Burk_R2_effects_with_dom.txt"),
                        sep="\t", header=TRUE, row.names=1)

combined_syncom_with_dom <- rbind(fam_data3, fam_data4) %>%
  mutate(dataset = "Dominators")

# Read data - No Dominators
fam_data7 <- read.table(paste0(working_directory, "Family_R2/SSC_Fam_R2_effects_no_dom.txt"),
                        sep="\t", header=TRUE, row.names=1)

fam_data8 <- read.table(paste0(working_directory, "Family_R2/SSC_Gen_Burk_R2_effects_no_dom.txt"),
                        sep="\t", header=TRUE, row.names=1)

combined_syncom_no_dom <- rbind(fam_data7, fam_data8) %>%
  mutate(dataset = "No_Dominators")

# Merge datasets
pyramid_data_fam <- bind_rows(combined_syncom_with_dom, combined_syncom_no_dom)

# Filter for KO only
pyramid_data_fam <- pyramid_data_fam %>% filter(KO == "KO", Subset %in% SynComs)

# Factor levels
pyramid_data_fam$Subset <- factor(pyramid_data_fam$Subset, levels = c("SSC", "LjSC", "HvSC", "AtSC"))
family_order <- c("Pelomonas", "Cupriavidus", "Polaromonas",
                  "Variovorax", "Rhizobacter", "Acidovorax", "Chitinophagaceae", "Microbacteriaceae", "Micrococcaceae",
                  "Xanthobacteraceae", "Sphingobacteriaceae", "Rhodanobacteraceae",
                  "Sphingomonadaceae", "Flavobacteriaceae", "Devosiaceae",
                  "Beijerinckiaceae", "Enterobacteriaceae", "Caulobacteraceae",
                  "Pseudomonadaceae", "Xanthomonadaceae", "Burkholderiaceae", "Rhizobiaceae")

pyramid_data_fam$Family <- factor(pyramid_data_fam$Family, levels = family_order)
pyramid_data_fam$R2_change <- abs(pyramid_data_fam$R2_change)
pyramid_data_fam$dataset <- factor(pyramid_data_fam$dataset, levels = c("Dominators", "No_Dominators"))

# Define bar width
barwidth <- 0.4

# Generate Dominators and Non-Dominators datasets
Dominators <- filter(pyramid_data_fam, dataset == "Dominators") %>%
  group_by(Family) %>%
  arrange(-Subset)

No_Dominators <- filter(pyramid_data_fam, dataset == "No_Dominators") %>%
  group_by(Family) %>%
  arrange(-Subset)

# Plot native
plot9 <- ggplot() +
  geom_bar(data = Dominators,
           mapping = aes(x = as.numeric(Family) + barwidth/2, y = R2_change, fill = Subset, alpha = "Dominators"),
           stat = "identity",
           position = 'stack',
           color = "black",
           size = 0.2,
           width = barwidth) +
  geom_bar(data = No_Dominators,
           mapping = aes(x = as.numeric(Family) - barwidth/2, y = R2_change, fill = Subset, alpha = "No_Dominators", color=Subset),
           stat = "identity",
           position = 'stack',
           size = 0.2,
           width = barwidth) +
  ggtitle("Family R² Effects With vs. Without Dominators") +
  theme(plot.title = element_text(hjust = 0.5, size = 10)) +
  theme_classic() +
  labs(x = "Family", y = "Effect on SynCom R²", fill = "Inoculum", alpha="Dataset") +
  scale_fill_manual(values = SC_colors) +
  scale_color_manual(values = SC_colors) +
  scale_alpha_manual(values = c("Dominators" = 1, "No_Dominators" = 0.3)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black", size = 7),
        legend.position = "right",
        legend.background = element_blank(),
        legend.key = element_blank(),
        text = element_text(family = "sans", size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = "italic", size = 10)) +
  scale_x_continuous(breaks = 1:length(levels(pyramid_data_fam$Family)), labels = levels(pyramid_data_fam$Family))+
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01), expand = expansion(mult = c(0, 0.05)), breaks = seq(0,0.3, 0.1)) +
  coord_flip()

# Display the plot
plot9

# Plot with Rhizobiaceae break
plot9_lim <- ggplot() +
  geom_bar(data = Dominators,
           mapping = aes(x = as.numeric(Family) + barwidth/2, y = R2_change, fill = Subset, alpha = "Dominators"),
           stat = "identity",
           position = 'stack',
           color = "black",
           size = 0.2,
           width = barwidth) +
  geom_bar(data = No_Dominators,
           mapping = aes(x = as.numeric(Family) - barwidth/2, y = R2_change, fill = Subset, alpha = "No_Dominators", color=Subset),
           stat = "identity",
           position = 'stack',
           size = 0.2,
           width = barwidth) +
  ggtitle("Family R² Effects With vs. Without Dominators") +
  theme(plot.title = element_text(hjust = 0.5, size = 10)) +
  theme_classic() +
  labs(x = "Family", y = "Effect on SynCom R²", fill = "Inoculum", alpha="Dataset") +
  scale_fill_manual(values = SC_colors) +
  scale_color_manual(values = SC_colors) +
  scale_alpha_manual(values = c("Dominators" = 1, "No_Dominators" = 0.3)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black", size = 7),
        legend.position = "right",
        legend.background = element_blank(),
        legend.key = element_blank(),
        text = element_text(family = "sans", size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = "italic", size = 10)) +
  scale_x_continuous(breaks = 1:length(levels(pyramid_data_fam$Family)), labels = levels(pyramid_data_fam$Family))+
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01), expand = expansion(mult = c(0, 0.05)), breaks = seq(0,0.3, 0.1), limits = c(0,0.32) ) +
  coord_flip()

plot9_lim

# Combine syncom and plant R² family effect (NATIVE)
pyramid_dom=ggarrange(print(plot8), print(plot9), 
                      ncol = 2, nrow = 1,    
                      common.legend = F,  
                      legend = "none")      
pyramid_dom

# Combine syncom and plant R² family effect with limits
pyramid_dom_lim=ggarrange(print(plot8), print(plot9_lim), 
                          ncol = 2, nrow = 1,    
                          common.legend = F,  
                          legend = "none")      
pyramid_dom_lim

# Extract legends from both plots
legend8 <- get_legend(plot8 + theme(legend.position = "right"))
legend9 <- get_legend(plot9 + theme(legend.position = "right"))

# Arrange legends in a single plot
legend_plot <- ggarrange(legend8, legend9, ncol = 2, nrow = 1)

# Display the legend plot
legend_plot


# # Save plot
pdf(paste0(results.dir, "Figure_S16a_Pyramid_Plot_main_lim.pdf"), width = 10, height = 8)
print(pyramid_dom_lim)
dev.off()

pdf(paste0(results.dir, "Figure_S16a_Pyramid_Plot_legends.pdf"), width = 10, height = 8)
print(legend_plot)
dev.off()


###Figure S16 - boxplots ======

#Burkholderiaceae AtSC no dominators - Plant

# Load the KO abundance tables (Burkholderiaceae-inclusive and exclusive datasets)
KO_SSC_only <- read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Burkholderiaceae_KO_no_dom.tsv", sep = ""),
                          header = TRUE, sep = "\t", row.names = 1)
KO_SSC_ex_only <- read.table(paste(working_directory, "sPLS-DA/isolate_drop_out_data/No_burkholderiaceae_KO_no_dom.tsv", sep = ""),
                             header = TRUE, sep = "\t", row.names = 1)

# Load metadata and filter for relevant samples
samples_df <- read.table(paste(working_directory, "SSC_R2_metadata_no_HL.tsv", sep = ""), header = TRUE, sep = "\t", row.names = 1)
colnames(samples_df)[5] <- "Nutrient"

samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")
samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum == "AtSC")

# Normalize KO tables by row sums (relative abundance)
KO_SSC_only_2 <- t(t(KO_SSC_only) / rowSums(t(KO_SSC_only)))
KO_SSC_ex_only_2 <- t(t(KO_SSC_ex_only) / rowSums(t(KO_SSC_ex_only)))

# Extract KOs associated with Burkholderiaceae in AtSC
KOs <- empty_vector_all_2_all_sub[empty_vector_all_2_all_sub$Group == "Burkholderiaceae_no_dominators_AtSC", ]
KOs_Lj <- KOs$KO[KOs$Data == "Lj"]
KOs_Hv <- KOs$KO[KOs$Data == "Hv"]
KOs_At <- KOs$KO[KOs$Data == "At"]

# Subset for Arabidopsis KOs
KO_SSC_only_2_At <- data.frame(KO_SSC_only_2[row.names(KO_SSC_only_2) %in% KOs_At, ])
colnames(KO_SSC_only_2_At) <- "KO"
KO_SSC_ex_only_2_At <- data.frame(KO_SSC_ex_only_2[row.names(KO_SSC_ex_only_2) %in% KOs_At, ])
colnames(KO_SSC_ex_only_2_At) <- "KO"

KO_SSC_only_2_At$Dataset <- "Burkholderiaceae AtSC"
KO_SSC_ex_only_2_At$Dataset <- "Other families AtSC" 

# Subset for Barley KOs and sum across rows
KO_SSC_only_2_Hv <- KO_SSC_only_2[row.names(KO_SSC_only_2) %in% KOs_Hv, ]
KO_SSC_ex_only_2_Hv <- KO_SSC_ex_only_2[row.names(KO_SSC_ex_only_2) %in% KOs_Hv, ]

KO_SSC_only_2_Hv_2 <- data.frame(colSums(KO_SSC_only_2_Hv))
KO_SSC_ex_only_2_Hv_2 <- data.frame(colSums(KO_SSC_ex_only_2_Hv))

colnames(KO_SSC_only_2_Hv_2) <- "KO"
colnames(KO_SSC_ex_only_2_Hv_2) <- "KO"

KO_SSC_only_2_Hv_2$Dataset <- "Burkholderiaceae AtSC"
KO_SSC_ex_only_2_Hv_2$Dataset <- "Other families AtSC"

# Subset for Lotus KOs and sum across rows
KO_SSC_only_2_Lj <- KO_SSC_only_2[row.names(KO_SSC_only_2) %in% KOs_Lj, ]
KO_SSC_ex_only_2_Lj <- KO_SSC_ex_only_2[row.names(KO_SSC_ex_only_2) %in% KOs_Lj, ]

KO_SSC_only_2_Lj_2 <- data.frame(colSums(KO_SSC_only_2_Lj))
KO_SSC_ex_only_2_Lj_2 <- data.frame(colSums(KO_SSC_ex_only_2_Lj))

colnames(KO_SSC_only_2_Lj_2) <- "KO"
colnames(KO_SSC_ex_only_2_Lj_2) <- "KO"

KO_SSC_only_2_Lj_2$Dataset <- "Burkholderiaceae AtSC"
KO_SSC_ex_only_2_Lj_2$Dataset <- "Other families AtSC"

# Assign contributing plant information
KO_SSC_only_2_At$Contributing_plant <- "Arabidopsis KOs"
KO_SSC_ex_only_2_At$Contributing_plant <- "Arabidopsis KOs"
KO_SSC_only_2_Hv_2$Contributing_plant <- "Barley KOs"
KO_SSC_ex_only_2_Hv_2$Contributing_plant <- "Barley KOs"
KO_SSC_only_2_Lj_2$Contributing_plant <- "Lotus KOs"
KO_SSC_ex_only_2_Lj_2$Contributing_plant <- "Lotus KOs"

# Assign sample names
KO_SSC_only_2_At$Sample <- row.names(KO_SSC_only_2_At)
KO_SSC_ex_only_2_At$Sample <- row.names(KO_SSC_ex_only_2_At)
KO_SSC_only_2_Hv_2$Sample <- row.names(KO_SSC_only_2_Hv_2)
KO_SSC_ex_only_2_Hv_2$Sample <- row.names(KO_SSC_ex_only_2_Hv_2)
KO_SSC_only_2_Lj_2$Sample <- row.names(KO_SSC_only_2_Lj_2)
KO_SSC_ex_only_2_Lj_2$Sample <- row.names(KO_SSC_ex_only_2_Lj_2)

# Function to subset and add plant information
subset_and_add_plant <- function(df) {
  df <- df[df$Sample %in% row.names(samples_df_sub_3), ]
  df$Plant <- samples_df_sub_3$Condition[match(df$Sample, row.names(samples_df_sub_3))]
  df$Plant <- as.factor(df$Plant)
  return(df)
}

# Apply the function to each data frame
KO_SSC_only_2_At <- subset_and_add_plant(KO_SSC_only_2_At)
KO_SSC_ex_only_2_At <- subset_and_add_plant(KO_SSC_ex_only_2_At)
KO_SSC_only_2_Hv_2 <- subset_and_add_plant(KO_SSC_only_2_Hv_2)
KO_SSC_ex_only_2_Hv_2 <- subset_and_add_plant(KO_SSC_ex_only_2_Hv_2)
KO_SSC_only_2_Lj_2 <- subset_and_add_plant(KO_SSC_only_2_Lj_2)
KO_SSC_ex_only_2_Lj_2 <- subset_and_add_plant(KO_SSC_ex_only_2_Lj_2)

# Function to perform ANOVA and Tukey test, and plot the results
anova_and_plot <- function(df, title) {
  # Perform ANOVA
  fitAnova <- aov(KO ~ Plant, data = df)
  
  # Perform Tukey's post-hoc test
  Tukey <- TukeyHSD(fitAnova)
  
  # Get letters for significance
  letters_anova <- multcompLetters4(fitAnova, Tukey)$Plant$Letters
  
  # Combine results into a data frame for plotting
  ltlbl_combined <- data.frame(
    Plant = names(letters_anova),
    Letters = letters_anova
  )
  
  # Ensure the order of the factor levels is correct
  ltlbl_combined$Plant <- factor(ltlbl_combined$Plant, levels = unique(df$Plant))
  
  # Apply the manually defined order to the factor
  ltlbl_combined <- ltlbl_combined[order(ltlbl_combined$Plant), ]
  
  # Merge the letters with the data frame
  df <- merge(df, ltlbl_combined, by = "Plant")
  
  # Plot the results
  plot <- ggplot(df, aes(x = Plant, y = KO, colour = Plant)) +
    geom_boxplot(outlier.shape = NA) +
    theme_classic() +
    ylab("Cumulative RA of KOs") +
    geom_jitter(shape = 16, position = position_jitter(0.2), aes(colour = Plant), show.legend = TRUE) +
    theme(
      axis.text.x = element_blank(),
      axis.title.y = element_text(size = 14),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14),
      plot.title = element_text(size = 14),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      strip.text.x = element_text(size = 10)
    ) +
    scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
    guides(colour = FALSE) +
    facet_wrap(~Contributing_plant + Dataset, scales = "free") +
    scale_y_continuous(labels = scales::scientific,  ) +
    stat_summary(aes(label = Letters, y = max(KO)*1.03), fun = max, geom = "text")
  
  return(plot)
}

# Generate plots for each data frame
plot_At_Burk <- anova_and_plot(KO_SSC_only_2_At, "Arabidopsis KOs - Burkholderiaceae AtSC")
plot_At_other <- anova_and_plot(KO_SSC_ex_only_2_At, "Arabidopsis KOs - Other families AtSC")
plot_Hv_Burk <- anova_and_plot(KO_SSC_only_2_Hv_2, "Barley KOs - Burkholderiaceae AtSC")
plot_Hv_other <- anova_and_plot(KO_SSC_ex_only_2_Hv_2, "Barley KOs - Other families AtSC")
plot_Lj_Burk <- anova_and_plot(KO_SSC_only_2_Lj_2, "Lotus KOs - Burkholderiaceae AtSC")
plot_Lj_other <- anova_and_plot(KO_SSC_ex_only_2_Lj_2, "Lotus KOs - Other families AtSC")

# Arrange plots using ggarrange
Burk_plot=grid.arrange(plot_At_Burk,plot_At_other,plot_Hv_Burk, plot_Hv_other, plot_Lj_Burk, plot_Lj_other, nrow = 3)

#Caulobacteraceae AtSC without dominators - Plant
# Load the KO abundance tables (Caulobacteraceae-inclusive and exclusive datasets)
KO_SSC_only <- read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Caulobacteraceae_KO_no_dom.tsv", sep = ""),
                          header = TRUE, sep = "\t", row.names = 1)
KO_SSC_ex_only <- read.table(paste(working_directory, "sPLS-DA/isolate_drop_out_data/No_caulobacteraceae_KO_no_dom.tsv", sep = ""),
                             header = TRUE, sep = "\t", row.names = 1)

# Load metadata and filter for relevant samples
samples_df <- read.table(paste(working_directory, "SSC_R2_metadata_no_HL.tsv", sep = ""), header = TRUE, sep = "\t", row.names = 1)
colnames(samples_df)[5] <- "Nutrient"

samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")
samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum == "AtSC")

# Normalize KO tables by row sums (relative abundance)
KO_SSC_only_2 <- t(t(KO_SSC_only) / rowSums(t(KO_SSC_only)))
KO_SSC_ex_only_2 <- t(t(KO_SSC_ex_only) / rowSums(t(KO_SSC_ex_only)))

# Extract KOs associated with Caulobacteraceae in AtSC
KOs <- empty_vector_all_2_all_sub[empty_vector_all_2_all_sub$Group == "Caulobacteraceae_no_dominators_AtSC", ]
KOs_Hv <- KOs$KO[KOs$Data == "Hv"]
KOs_At <- KOs$KO[KOs$Data == "At"]

# Subset for Arabidopsis KOs and sum across rows
KO_SSC_only_2_At <- KO_SSC_only_2[row.names(KO_SSC_only_2) %in% KOs_At, ]
KO_SSC_ex_only_2_At <- KO_SSC_ex_only_2[row.names(KO_SSC_ex_only_2) %in% KOs_At, ]

KO_SSC_only_2_At_2 <- data.frame(colSums(KO_SSC_only_2_At))
KO_SSC_ex_only_2_At_2 <- data.frame(colSums(KO_SSC_ex_only_2_At))

colnames(KO_SSC_only_2_At_2) <- "KO"
colnames(KO_SSC_ex_only_2_At_2) <- "KO"

KO_SSC_only_2_At_2$Dataset <- "Caulobacteraceae AtSC"
KO_SSC_ex_only_2_At_2$Dataset <- "Other families AtSC"

# Subset for Barley KOs and sum across rows
KO_SSC_only_2_Hv <- KO_SSC_only_2[row.names(KO_SSC_only_2) %in% KOs_Hv, ]
KO_SSC_ex_only_2_Hv <- KO_SSC_ex_only_2[row.names(KO_SSC_ex_only_2) %in% KOs_Hv, ]

KO_SSC_only_2_Hv_2 <- data.frame(colSums(KO_SSC_only_2_Hv))
KO_SSC_ex_only_2_Hv_2 <- data.frame(colSums(KO_SSC_ex_only_2_Hv))

colnames(KO_SSC_only_2_Hv_2) <- "KO"
colnames(KO_SSC_ex_only_2_Hv_2) <- "KO"

KO_SSC_only_2_Hv_2$Dataset <- "Caulobacteraceae AtSC"
KO_SSC_ex_only_2_Hv_2$Dataset <- "Other families AtSC"

# Assign contributing plant information
KO_SSC_only_2_At_2$Contributing_plant <- "Arabidopsis KOs"
KO_SSC_ex_only_2_At_2$Contributing_plant <- "Arabidopsis KOs"
KO_SSC_only_2_Hv_2$Contributing_plant <- "Barley KOs"
KO_SSC_ex_only_2_Hv_2$Contributing_plant <- "Barley KOs"

# Assign sample names
KO_SSC_only_2_At_2$Sample <- row.names(KO_SSC_only_2_At_2)
KO_SSC_ex_only_2_At_2$Sample <- row.names(KO_SSC_ex_only_2_At_2)
KO_SSC_only_2_Hv_2$Sample <- row.names(KO_SSC_only_2_Hv_2)
KO_SSC_ex_only_2_Hv_2$Sample <- row.names(KO_SSC_ex_only_2_Hv_2)

# Function to subset and add plant information
subset_and_add_plant <- function(df) {
  df <- df[df$Sample %in% row.names(samples_df_sub_3), ]
  df$Plant <- samples_df_sub_3$Condition[match(df$Sample, row.names(samples_df_sub_3))]
  df$Plant <- as.factor(df$Plant)
  return(df)
}

# Apply the function to each data frame
KO_SSC_only_2_At_2 <- subset_and_add_plant(KO_SSC_only_2_At_2)
KO_SSC_ex_only_2_At_2 <- subset_and_add_plant(KO_SSC_ex_only_2_At_2)
KO_SSC_only_2_Hv_2 <- subset_and_add_plant(KO_SSC_only_2_Hv_2)
KO_SSC_ex_only_2_Hv_2 <- subset_and_add_plant(KO_SSC_ex_only_2_Hv_2)

# Function to perform ANOVA and Tukey test, and plot the results
anova_and_plot <- function(df, title) {
  # Perform ANOVA
  fitAnova <- aov(KO ~ Plant, data = df)
  
  # Perform Tukey's post-hoc test
  Tukey <- TukeyHSD(fitAnova)
  
  # Get letters for significance
  letters_anova <- multcompLetters4(fitAnova, Tukey)$Plant$Letters
  
  # Combine results into a data frame for plotting
  ltlbl_combined <- data.frame(
    Plant = names(letters_anova),
    Letters = letters_anova
  )
  
  # Ensure the order of the factor levels is correct
  ltlbl_combined$Plant <- factor(ltlbl_combined$Plant, levels = unique(df$Plant))
  
  # Apply the manually defined order to the factor
  ltlbl_combined <- ltlbl_combined[order(ltlbl_combined$Plant), ]
  
  # Merge the letters with the data frame
  df <- merge(df, ltlbl_combined, by = "Plant")
  
  # Plot the results
  plot <- ggplot(df, aes(x = Plant, y = KO, colour = Plant)) +
    geom_boxplot(outlier.shape = NA) +
    theme_classic() +
    ylab("Cumulative RA of KOs") +
    geom_jitter(shape = 16, position = position_jitter(0.2), aes(colour = Plant), show.legend = TRUE) +
    theme(
      axis.text.x = element_blank(),
      axis.title.y = element_text(size = 14),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14),
      plot.title = element_text(size = 14),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      strip.text.x = element_text(size = 10)
    ) +
    scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
    guides(colour = FALSE) +
    facet_wrap(~Contributing_plant + Dataset, scales = "free") +
    scale_y_continuous(labels = scales::scientific,  ) +
    stat_summary(aes(label = Letters, y = max(KO)*1.03), fun = max, geom = "text")
  
  return(plot)
}

# Generate plots for each data frame
plot_At_Caulo <- anova_and_plot(KO_SSC_only_2_At_2, "Arabidopsis KOs - Caulobacteraceae AtSC")
plot_At_other_caulo <- anova_and_plot(KO_SSC_ex_only_2_At_2, "Arabidopsis KOs - Other families AtSC")
plot_Hv_Caulo <- anova_and_plot(KO_SSC_only_2_Hv_2, "Barley KOs - Caulobacteraceae AtSC")
plot_Hv_other_caulo <- anova_and_plot(KO_SSC_ex_only_2_Hv_2, "Barley KOs - Other families AtSC")

# Arrange plots using ggarrange
Caulo_plot=grid.arrange(plot_At_Caulo,plot_At_other_caulo,plot_Hv_Caulo, plot_Hv_other_caulo,nrow = 2)

#Xanthomonadaceae LjSC without dominators - Plant
# Load the KO abundance tables (Xanthomonadaceae-inclusive and exclusive datasets)
KO_SSC_only <- read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Xanthomonadaceae_KO_no_dom.tsv", sep = ""),
                          header = TRUE, sep = "\t", row.names = 1)
KO_SSC_ex_only <- read.table(paste(working_directory, "sPLS-DA/isolate_drop_out_data/No_xanthomonadaceae_KO_no_dom.tsv", sep = ""),
                             header = TRUE, sep = "\t", row.names = 1)

# Load metadata and filter for relevant samples
samples_df <- read.table(paste(working_directory, "SSC_R2_metadata_no_HL.tsv", sep = ""), header = TRUE, sep = "\t", row.names = 1)
colnames(samples_df)[5] <- "Nutrient"

samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")
samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum == "LjSC")

# Normalize KO tables by row sums (relative abundance)
KO_SSC_only_2 <- t(t(KO_SSC_only) / rowSums(t(KO_SSC_only)))
KO_SSC_ex_only_2 <- t(t(KO_SSC_ex_only) / rowSums(t(KO_SSC_ex_only)))

# Extract KOs associated with Xanthomonadaceae and dominant families
KOs <- empty_vector_all_2_all_sub[empty_vector_all_2_all_sub$Group == "Xanthomonadaceae_no_dominators_LjSC", ]
KOs_Hv <- KOs$KO[KOs$Data == "Hv"]
KOs_At <- KOs$KO[KOs$Data == "At"]

# Subset for Arabidopsis KOs
KO_SSC_only_2_At <- data.frame(KO_SSC_only_2[row.names(KO_SSC_only_2) %in% KOs_At, ])
colnames(KO_SSC_only_2_At) <- "KO"
KO_SSC_ex_only_2_At <- data.frame(KO_SSC_ex_only_2[row.names(KO_SSC_ex_only_2) %in% KOs_At, ])
colnames(KO_SSC_ex_only_2_At) <- "KO"

KO_SSC_only_2_At$Dataset <- "Xanthomonadaceae LjSC"
KO_SSC_ex_only_2_At$Dataset <- "other families LjSC"

# Subset for Barley KOs and sum across rows
KO_SSC_only_2_Hv <- KO_SSC_only_2[row.names(KO_SSC_only_2) %in% KOs_Hv, ]
KO_SSC_ex_only_2_Hv <- KO_SSC_ex_only_2[row.names(KO_SSC_ex_only_2) %in% KOs_Hv, ]

KO_SSC_only_2_Hv_2 <- data.frame(colSums(KO_SSC_only_2_Hv))
KO_SSC_ex_only_2_Hv_2 <- data.frame(colSums(KO_SSC_ex_only_2_Hv))

colnames(KO_SSC_only_2_Hv_2) <- "KO"
colnames(KO_SSC_ex_only_2_Hv_2) <- "KO"

KO_SSC_only_2_Hv_2$Dataset <- "Xanthomonadaceae LjSC"
KO_SSC_ex_only_2_Hv_2$Dataset <- "Other families LjSC"

# Assign contributing plant information
KO_SSC_only_2_At$Contributing_plant <- "Arabidopsis KOs"
KO_SSC_ex_only_2_At$Contributing_plant <- "Arabidopsis KOs"
KO_SSC_only_2_Hv_2$Contributing_plant <- "Barley KOs"
KO_SSC_ex_only_2_Hv_2$Contributing_plant <- "Barley KOs"

# Assign sample names
KO_SSC_only_2_At$Sample <- row.names(KO_SSC_only_2_At)
KO_SSC_ex_only_2_At$Sample <- row.names(KO_SSC_ex_only_2_At)
KO_SSC_only_2_Hv_2$Sample <- row.names(KO_SSC_only_2_Hv_2)
KO_SSC_ex_only_2_Hv_2$Sample <- row.names(KO_SSC_ex_only_2_Hv_2)

# Function to subset and add plant information
subset_and_add_plant <- function(df) {
  df <- df[df$Sample %in% row.names(samples_df_sub_3), ]
  df$Plant <- samples_df_sub_3$Condition[match(df$Sample, row.names(samples_df_sub_3))]
  df$Plant <- as.factor(df$Plant)
  return(df)
}

# Apply the function to each data frame
KO_SSC_only_2_At <- subset_and_add_plant(KO_SSC_only_2_At)
KO_SSC_ex_only_2_At <- subset_and_add_plant(KO_SSC_ex_only_2_At)
KO_SSC_only_2_Hv_2 <- subset_and_add_plant(KO_SSC_only_2_Hv_2)
KO_SSC_ex_only_2_Hv_2 <- subset_and_add_plant(KO_SSC_ex_only_2_Hv_2)

# Function to perform ANOVA and Tukey test, and plot the results
anova_and_plot <- function(df, title) {
  # Perform ANOVA
  fitAnova <- aov(KO ~ Plant, data = df)
  
  # Perform Tukey's post-hoc test
  Tukey <- TukeyHSD(fitAnova)
  
  # Get letters for significance
  letters_anova <- multcompLetters4(fitAnova, Tukey)$Plant$Letters
  
  # Combine results into a data frame for plotting
  ltlbl_combined <- data.frame(
    Plant = names(letters_anova),
    Letters = letters_anova
  )
  
  # Ensure the order of the factor levels is correct
  ltlbl_combined$Plant <- factor(ltlbl_combined$Plant, levels = unique(df$Plant))
  
  # Apply the manually defined order to the factor
  ltlbl_combined <- ltlbl_combined[order(ltlbl_combined$Plant), ]
  
  # Merge the letters with the data frame
  df <- merge(df, ltlbl_combined, by = "Plant")
  
  # Plot the results
  plot <- ggplot(df, aes(x = Plant, y = KO, colour = Plant)) +
    geom_boxplot(outlier.shape = NA) +
    theme_classic() +
    ylab("Cumulative RA of KOs") +
    geom_jitter(shape = 16, position = position_jitter(0.2), aes(colour = Plant), show.legend = TRUE) +
    theme(
      axis.text.x = element_blank(),
      axis.title.y = element_text(size = 14),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14),
      plot.title = element_text(size = 14),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      strip.text.x = element_text(size = 10)
    ) +
    scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
    guides(colour = FALSE) +
    facet_wrap(~Contributing_plant + Dataset, scales = "free") +
    scale_y_continuous(labels = scales::scientific,  ) +
    stat_summary(aes(label = Letters, y = max(KO)*1.03), fun = max, geom = "text")
  
  return(plot)
}

# Generate plots for each data frame
plot_At_Xanth <- anova_and_plot(KO_SSC_only_2_At, "Arabidopsis KOs - Xanthomonadaceae LjSC")
plot_At_other_xanth <- anova_and_plot(KO_SSC_ex_only_2_At, "Arabidopsis KOs - Other families LjSC")
plot_Hv_Xanth <- anova_and_plot(KO_SSC_only_2_Hv_2, "Barley KOs - Xanthomonadaceae LjSC")
plot_Hv_other_xanth <- anova_and_plot(KO_SSC_ex_only_2_Hv_2, "Barley KOs - Other families LjSC")

# Arrange plots using ggarrange
Xanth_plot=grid.arrange(plot_At_Xanth,plot_At_other_xanth,plot_Hv_Xanth, plot_Hv_other_xanth, nrow = 2)

#Burkholderiaceae no dominators - SynComs
KO_SSC_only <- read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Burkholderiaceae_KO_no_dom.tsv", sep = ""),
                          header = TRUE, sep = "\t", row.names = 1)
KO_SSC_ex_only <- read.table(paste(working_directory, "sPLS-DA/isolate_drop_out_data/No_burkholderiaceae_KO_no_dom.tsv", sep = ""),
                             header = TRUE, sep = "\t", row.names = 1)

# Load metadata and filter for relevant samples
samples_df <- read.table(paste(working_directory, "SSC_R2_metadata_no_HL.tsv", sep = ""), header = TRUE, sep = "\t", row.names = 1)
colnames(samples_df)[5] <- "Nutrient"

samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

# Normalize KO tables by row sums (relative abundance)
KO_SSC_only_2 <- t(t(KO_SSC_only) / rowSums(t(KO_SSC_only)))
KO_SSC_ex_only_2 <- t(t(KO_SSC_ex_only) / rowSums(t(KO_SSC_ex_only)))

# Extract KOs associated with Burkholderiaceae and dominant families
KOs <- empty_vector_all_2_all_sub[empty_vector_all_2_all_sub$Group == "Burkholderiaceae_no_dominators", ]
KOs_AtSC <- KOs$KO[KOs$Data == "AtSC"]
KOs_LjSC <- KOs$KO[KOs$Data == "LjSC"]

# Subset for AtSC KOs
KO_SSC_only_2_AtSC <- data.frame(KO_SSC_only_2[row.names(KO_SSC_only_2) %in% KOs_AtSC, ])
colnames(KO_SSC_only_2_AtSC ) <- "KO"
KO_SSC_ex_only_2_AtSC <- data.frame(KO_SSC_ex_only_2[row.names(KO_SSC_ex_only_2) %in% KOs_AtSC, ])
colnames(KO_SSC_ex_only_2_AtSC) <- "KO"

KO_SSC_only_2_AtSC$Dataset <- "Burkholderiaceae"
KO_SSC_ex_only_2_AtSC$Dataset <- "Other families" 

# Subset for LjSC KOs
KO_SSC_only_2_LjSC <- data.frame(KO_SSC_only_2[row.names(KO_SSC_only_2) %in% KOs_LjSC, ])
colnames(KO_SSC_only_2_LjSC ) <- "KO"
KO_SSC_ex_only_2_LjSC <- data.frame(KO_SSC_ex_only_2[row.names(KO_SSC_ex_only_2) %in% KOs_LjSC, ])
colnames(KO_SSC_ex_only_2_LjSC) <- "KO"

KO_SSC_only_2_LjSC$Dataset <- "Burkholderiaceae"
KO_SSC_ex_only_2_LjSC$Dataset <- "Other families" 

# Assign contributing plant information
KO_SSC_only_2_AtSC$Contributing_SynCom <- "AtSC KOs"
KO_SSC_ex_only_2_AtSC$Contributing_SynCom <- "AtSC KOs"
KO_SSC_only_2_LjSC$Contributing_SynCom <- "LjSC KOs"
KO_SSC_ex_only_2_LjSC$Contributing_SynCom <- "LjSC KOs"

# Assign sample names
KO_SSC_only_2_AtSC$Sample <- row.names(KO_SSC_only_2_AtSC)
KO_SSC_ex_only_2_AtSC$Sample <- row.names(KO_SSC_ex_only_2_AtSC)
KO_SSC_only_2_LjSC$Sample <- row.names(KO_SSC_only_2_LjSC)
KO_SSC_ex_only_2_LjSC$Sample <- row.names(KO_SSC_ex_only_2_LjSC)

# Function to subset and add plant information
subset_and_add_syncom <- function(df) {
  df <- df[df$Sample %in% row.names(samples_df_sub_2), ]
  df$SynCom <- samples_df_sub_2$Inoculum[match(df$Sample, row.names(samples_df_sub_2))]
  df$SynCom <- as.factor(df$SynCom)
  return(df)
}

# Apply the function to each data frame
KO_SSC_only_2_AtSC <- subset_and_add_syncom(KO_SSC_only_2_AtSC)
KO_SSC_ex_only_2_AtSC <- subset_and_add_syncom(KO_SSC_ex_only_2_AtSC)
KO_SSC_only_2_LjSC <- subset_and_add_syncom(KO_SSC_only_2_LjSC)
KO_SSC_ex_only_2_LjSC <- subset_and_add_syncom(KO_SSC_ex_only_2_LjSC)

# Function to perform ANOVA and Tukey test, and plot the results
anova_and_plot <- function(df, title) {
  # Perform ANOVA
  fitAnova <- aov(KO ~ SynCom, data = df)
  
  # Perform Tukey's post-hoc test
  Tukey <- TukeyHSD(fitAnova)
  
  # Get letters for significance
  letters_anova <- multcompLetters4(fitAnova, Tukey)$SynCom$Letters
  
  # Combine results into a data frame for plotting
  ltlbl_combined <- data.frame(
    SynCom = names(letters_anova),
    Letters = letters_anova
  )
  
  # Ensure the order of the factor levels is correct
  ltlbl_combined$SynCom <- factor(ltlbl_combined$SynCom, levels = unique(df$SynCom))
  
  # Apply the manually defined order to the factor
  ltlbl_combined <- ltlbl_combined[order(ltlbl_combined$SynCom), ]
  
  # Merge the letters with the data frame
  df <- merge(df, ltlbl_combined, by = "SynCom")
  
  # Plot the results
  plot <- ggplot(df, aes(x = SynCom, y = KO, colour = SynCom)) +
    geom_boxplot(outlier.shape = NA) +
    theme_classic() +
    ylab("Cumulative RA of KOs") +
    geom_jitter(shape = 16, position = position_jitter(0.2), aes(colour = SynCom), show.legend = TRUE) +
    theme(
      axis.text.x = element_blank(),
      axis.title.y = element_text(size = 14),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14),
      plot.title = element_text(size = 14),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      strip.text.x = element_text(size = 10)
    ) +
    scale_color_manual(values = c("#A3A500","#00B0F6","#00BF7D","#F8766D" )) +
    guides(colour = FALSE) +
    facet_wrap(~Contributing_SynCom + Dataset, scales = "free") +
    scale_y_continuous(labels = scales::scientific,  ) +
    stat_summary(aes(label = Letters, y = max(KO)*1.03), fun = max, geom = "text")
  
  return(plot)
}

# Generate plots for each data frame
plot_AtSC_Burk <- anova_and_plot(KO_SSC_only_2_AtSC, "AtSC KOs - Burkholderiaceae")
plot_AtSC_other <- anova_and_plot(KO_SSC_ex_only_2_AtSC, "AtSC KOs - Other families")
plot_LjSC_Burk <- anova_and_plot(KO_SSC_only_2_LjSC, "LjSC KOs - Burkholderiaceae")
plot_LjSC_other <- anova_and_plot(KO_SSC_ex_only_2_LjSC, "LjSC KOs - Other families")

# Arrange plots using ggarrange
Burk_plot_2=grid.arrange(plot_AtSC_Burk,plot_AtSC_other,plot_LjSC_Burk, plot_LjSC_other, nrow = 2)

### Combine the the 3 with 2 host-specific or SynCom-specific KOs
combined_plot_2 <- 
  (plot_At_Caulo + plot_At_other_caulo + plot_Hv_Caulo +  plot_Hv_other_caulo) / 
  (plot_At_Xanth + plot_At_other_xanth + plot_Hv_Xanth + plot_Hv_other_xanth) /
  (plot_AtSC_Burk + plot_AtSC_other + plot_LjSC_Burk +  plot_LjSC_other) /
  plot_layout(ncol = 3)

print(combined_plot_2)

##Save plot
pdf(paste0(results.dir, "Figure_S16_boxplots_1.pdf"), width = 21, height = 8)
print(combined_plot_2)
dev.off()

#Just the Burkholderiaceae AtSC without dominators Plant plots
combined_plot_3 <- 
  (plot_At_Burk + plot_At_other) /
  (plot_Hv_Burk + plot_Hv_other) /
  (plot_Lj_Burk + plot_Lj_other) /
  plot_layout(ncol = 3)

print(combined_plot_3)

## Save plot
pdf(paste0(results.dir, "Figure_S16_boxplots_2.pdf"), width = 21, height = 4)
print(combined_plot_3)
dev.off()

###Figure S16 - PCoAs to compare sPLS-DA results =====
#Rhizobiaceae - with dom - HvSC, LjSC, SSC
#otu table
KO_SSC_only=read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Rhizobiaceae_KO_with_dom.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t", row.names =1) #make the SampleID column into the row.names
colnames(samples_df)[5]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object
#Sample subsetting
samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

#Subset for all SynComs but AtSC
samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum != "AtSC")

OTU_KO = otu_table(as.matrix(KO_SSC_only),taxa_are_rows = TRUE)
samples_sub = sample_data(samples_df_sub_3)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

bray_2 <- as.matrix(beta_isolate_KO)

str(samples_df_sub_3)
str(bray_2)

#Bind metadata with distance matrix
pcoa = cmdscale(bray_2, k=10, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z", "a", "b", "c", "d", "e", "f", "g") 
eig = pcoa$eig
points_2 <- points[order(row.names(points)), ]
samples_df_sub_6 <- samples_df_sub_3[row.names(samples_df_sub_3) %in% row.names(points),]
samples_df_sub_7 <- samples_df_sub_6[order(row.names(samples_df_sub_6)), ]
points_3 <- cbind(points_2,samples_df_sub_7)
colnames(points_3) <- c("PCoA Axis 1", "PCoA Axis 2", "PCoA Axis 3","PCoA Axis 4", "b", "c", "d", "e", "f", "g",colnames(samples_df_sub_7))

groups <- c("PCoA Axis 1", "PCoA Axis 2", "PCoA Axis 3","PCoA Axis 4")

for (group in groups) {
  for (group_2 in groups) {
    if (group != group_2) {  # Avoid comparing the same group with itself
      points_sub <- points_3[, colnames(points_3) %in% c(group, group_2, "Condition")]
      if (group == "PCoA Axis 1"){
        axis_x <- paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep="")
      } else if (group == "PCoA Axis 2"){
        axis_x <- paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")
      } else if (group == "PCoA Axis 3"){
        axis_x <- paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep="")
      } else {
        axis_x <- paste("PCoA 4 (", format(100 * eig[4] / sum(eig), digits=4), "%)", sep="")
      }
      
      if (group_2 == "PCoA Axis 1"){
        axis_y <- paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep="")
      } else if (group_2 == "PCoA Axis 2"){
        axis_y <- paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")
      } else if (group_2 == "PCoA Axis 3"){
        axis_y <- paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep="")
      } else {
        axis_y <- paste("PCoA 4 (", format(100 * eig[4] / sum(eig), digits=4), "%)", sep="")
      }
      
      # PCoA plot with vectors
      pcoa_plot <- ggplot(points_sub, aes_string(x = paste("`", group, "`", sep = ""), 
                                                 y = paste("`", group_2, "`", sep = ""), 
                                                 color = "Condition")) +
        geom_point(size = 3) +
        theme_classic() +
        scale_colour_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +  # Center title
        labs(title = "",
             x = axis_x, 
             y = axis_y)
      
      # Store the plot dynamically
      assign(paste(group, group_2, sep = "_"), pcoa_plot)
    }
  }
}


one <- ggarrange(`PCoA Axis 1_PCoA Axis 2`, `PCoA Axis 1_PCoA Axis 3`, `PCoA Axis 1_PCoA Axis 4`, nrow =3, ncol =1, common.legend = T)
two <- ggarrange(NULL, `PCoA Axis 2_PCoA Axis 3`, `PCoA Axis 2_PCoA Axis 4`, nrow =3, ncol =1,common.legend = T)
three <- ggarrange(NULL, NULL, `PCoA Axis 3_PCoA Axis 4`, nrow =3, ncol =1,common.legend = T)

all <- ggarrange(one, two, three, ncol = 3, nrow =1,common.legend = T)
all

pdf(paste(results.dir,"Figure_S16_Rhizobiaceae_PC4.pdf", sep=""), width=10, height=10)
print(all)
dev.off()

#Xanthomonadaceae - SynCom - no dom 

#otu table
KO_SSC_only=read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Xanthomonadaceae_KO_no_dom.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t", row.names =1) #make the SampleID column into the row.names
colnames(samples_df)[5]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object

#Sample subsetting
samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

OTU_KO = otu_table(as.matrix(KO_SSC_only),taxa_are_rows = TRUE)
samples_sub = sample_data(samples_df_sub_2)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

bray_2 <- as.matrix(beta_isolate_KO)

str(samples_df_sub_2)
str(bray_2)

#Bind metadata with distance matrix
pcoa = cmdscale(bray_2, k=10, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z", "a", "b", "c", "d", "e", "f", "g") 
eig = pcoa$eig
points_2 <- points[order(row.names(points)), ]
samples_df_sub_6 <- samples_df_sub_2[row.names(samples_df_sub_2) %in% row.names(points),]
samples_df_sub_7 <- samples_df_sub_6[order(row.names(samples_df_sub_6)), ]
points_3 <- cbind(points_2,samples_df_sub_7)
colnames(points_3) <- c("PCoA Axis 1", "PCoA Axis 2", "PCoA Axis 3","PCoA Axis 4", "b", "c", "d", "e", "f", "g",colnames(samples_df_sub_7))

groups <- c("PCoA Axis 1", "PCoA Axis 2", "PCoA Axis 3","PCoA Axis 4")

for (group in groups) {
  for (group_2 in groups) {
    if (group != group_2) {  # Avoid comparing the same group with itself
      points_sub <- points_3[, colnames(points_3) %in% c(group, group_2, "Inoculum")]
      
      if (group == "PCoA Axis 1"){
        axis_x <- paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep="")
      } else if (group == "PCoA Axis 2"){
        axis_x <- paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")
      } else if (group == "PCoA Axis 3"){
        axis_x <- paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep="")
      } else {
        axis_x <- paste("PCoA 4 (", format(100 * eig[4] / sum(eig), digits=4), "%)", sep="")
      }
      
      if (group_2 == "PCoA Axis 1"){
        axis_y <- paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep="")
      } else if (group_2 == "PCoA Axis 2"){
        axis_y <- paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")
      } else if (group_2 == "PCoA Axis 3"){
        axis_y <- paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep="")
      } else {
        axis_y <- paste("PCoA 4 (", format(100 * eig[4] / sum(eig), digits=4), "%)", sep="")
      }
      
      # PCoA plot with vectors
      pcoa_plot <- ggplot(points_sub, aes_string(x = paste("`", group, "`", sep = ""), 
                                                 y = paste("`", group_2, "`", sep = ""), 
                                                 color = "Inoculum")) +
        geom_point(size = 3) +
        theme_classic() +
        scale_colour_manual(values = c("#A3A500","#00B0F6","#00BF7D","#F8766D")) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +  # Center title
        labs(title = "",
             x = axis_x, 
             y = axis_y)
      
      # Store the plot dynamically
      assign(paste(group, group_2, sep = "_"), pcoa_plot)
    }
  }
}

one <- ggarrange(`PCoA Axis 1_PCoA Axis 2`, `PCoA Axis 1_PCoA Axis 3`, `PCoA Axis 1_PCoA Axis 4`, nrow =3, ncol =1, common.legend = T)
two <- ggarrange(NULL, `PCoA Axis 2_PCoA Axis 3`, `PCoA Axis 2_PCoA Axis 4`, nrow =3, ncol =1,common.legend = T)
three <- ggarrange(NULL, NULL, `PCoA Axis 3_PCoA Axis 4`, nrow =3, ncol =1,common.legend = T)

all <- ggarrange(one, two, three, ncol = 3, nrow =1,common.legend = T)
all

pdf(paste(results.dir,"Figure_S16_Xanthomonadaceae_PC4.pdf", sep=""), width=10, height=10)
print(all)
dev.off()

#Burkholderiaceae - SynCom - no dom
#otu table
KO_SSC_only=read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Burkholderiaceae_KO_no_dom.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t", row.names =1) #make the SampleID column into the row.names
colnames(samples_df)[5]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object

#Sample subsetting
samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

#At subset - Enterobacteriaceae - dom vs no-dom 
OTU_KO = otu_table(as.matrix(KO_SSC_only),taxa_are_rows = TRUE)
samples_sub = sample_data(samples_df_sub_2)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

bray_2 <- as.matrix(beta_isolate_KO)

str(samples_df_sub_2)
str(bray_2)

#Bind metadata with distance matrix
pcoa = cmdscale(bray_2, k=10, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z", "a", "b", "c", "d", "e", "f", "g") 
eig = pcoa$eig
points_2 <- points[order(row.names(points)), ]
samples_df_sub_6 <- samples_df_sub_2[row.names(samples_df_sub_2) %in% row.names(points),]
samples_df_sub_7 <- samples_df_sub_6[order(row.names(samples_df_sub_6)), ]
points_3 <- cbind(points_2,samples_df_sub_7)
colnames(points_3) <- c("PCoA Axis 1", "PCoA Axis 2", "PCoA Axis 3","a", "b", "c", "d", "e", "f", "g",colnames(samples_df_sub_7))

groups <- c("PCoA Axis 1", "PCoA Axis 2", "PCoA Axis 3")

for (group in groups) {
  for (group_2 in groups) {
    if (group != group_2) {  # Avoid comparing the same group with itself
      points_sub <- points_3[, colnames(points_3) %in% c(group, group_2, "Inoculum")]
      
      if (group == "PCoA Axis 1"){
        axis_x <- paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep="")
      } else if (group == "PCoA Axis 2"){
        axis_x <- paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")
      } else if (group == "PCoA Axis 3"){
        axis_x <- paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep="")
      } 
      
      if (group_2 == "PCoA Axis 1"){
        axis_y <- paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep="")
      } else if (group_2 == "PCoA Axis 2"){
        axis_y <- paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")
      } else if (group_2 == "PCoA Axis 3"){
        axis_y <- paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep="")
      } 
      
      # PCoA plot with vectors
      pcoa_plot <- ggplot(points_sub, aes_string(x = paste("`", group, "`", sep = ""), 
                                                 y = paste("`", group_2, "`", sep = ""), 
                                                 color = "Inoculum")) +
        geom_point(size = 3) +
        theme_classic() +
        scale_colour_manual(values = c("#A3A500","#00B0F6","#00BF7D","#F8766D")) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +  # Center title
        labs(title = "",
             x = axis_x, 
             y = axis_y)
      
      # Store the plot dynamically
      assign(paste(group, group_2, sep = "_"), pcoa_plot)
    }
  }
}

one <- ggarrange(`PCoA Axis 1_PCoA Axis 2`, `PCoA Axis 1_PCoA Axis 3`, nrow =2, ncol =1, common.legend = T)
two <- ggarrange(NULL, `PCoA Axis 2_PCoA Axis 3`, nrow =2, ncol =1,common.legend = T)

all <- ggarrange(one, two,ncol = 2, nrow =1,common.legend = T)
all

pdf(paste(results.dir,"Figure_S16_Burkholderiaceae_PC3.pdf", sep=""), width=8, height=8)
print(all)
dev.off()

#Xanthomondaceae - Plant - LjSC - no dom
#otu table
KO_SSC_only=read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Xanthomonadaceae_KO_no_dom.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t", row.names =1) #make the SampleID column into the row.names
colnames(samples_df)[5]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object

#Sample subsetting
samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

#Subset for LjSC
samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum == "LjSC")

OTU_KO = otu_table(as.matrix(KO_SSC_only),taxa_are_rows = TRUE)
samples_sub = sample_data(samples_df_sub_3)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

bray_2 <- as.matrix(beta_isolate_KO)

str(samples_df_sub_3)
str(bray_2)

#Bind metadata with distance matrix
pcoa = cmdscale(bray_2, k=10, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z", "a", "b", "c", "d", "e", "f", "g") 
eig = pcoa$eig
points_2 <- points[order(row.names(points)), ]
samples_df_sub_6 <- samples_df_sub_3[row.names(samples_df_sub_3) %in% row.names(points),]
samples_df_sub_7 <- samples_df_sub_6[order(row.names(samples_df_sub_6)), ]
points_3 <- cbind(points_2,samples_df_sub_7)
colnames(points_3) <- c("PCoA Axis 1", "PCoA Axis 2", "z","a", "b", "c", "d", "e", "f", "g",colnames(samples_df_sub_7))
groups <- c("PCoA Axis 1", "PCoA Axis 2")

# Avoid comparing the same group with itself
points_sub <- points_3[, colnames(points_3) %in% groups]
points_sub$Condition <- samples_df_sub_7$Condition[match(row.names(points_sub), row.names(samples_df_sub_7))]

axis_x <- paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep="")
axis_y <- paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")

pcoa_plot <- ggplot(points_sub, aes_string(x = "`PCoA Axis 1`", 
                                           y = "`PCoA Axis 2`", 
                                           color = "Condition")) +
  geom_point(size = 3) +
  theme_classic() +
  scale_colour_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +  # Center title
  labs(title = "",
       x = axis_x, 
       y = axis_y)
pcoa_plot


pdf(paste(results.dir,"Figure_S16_Xanthomonadaceae_PC2_plant.pdf", sep=""), width=6, height=5)
print(pcoa_plot)
dev.off()

#Pseudomonadaceae - Plant - LjSC - no dom 

#otu table
KO_SSC_only=read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Pseudomonadaceae_KO_no_dom.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t", row.names =1) #make the SampleID column into the row.names
colnames(samples_df)[5]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object

#Sample subsetting
samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

#Subset for LjSC
samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum == "LjSC")

OTU_KO = otu_table(as.matrix(KO_SSC_only),taxa_are_rows = TRUE)
samples_sub = sample_data(samples_df_sub_3)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

bray_2 <- as.matrix(beta_isolate_KO)

str(samples_df_sub_3)
str(bray_2)

#Bind metadata with distance matrix
pcoa = cmdscale(bray_2, k=10, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z", "a", "b", "c", "d", "e", "f", "g") 
eig = pcoa$eig
points_2 <- points[order(row.names(points)), ]
samples_df_sub_6 <- samples_df_sub_3[row.names(samples_df_sub_3) %in% row.names(points),]
samples_df_sub_7 <- samples_df_sub_6[order(row.names(samples_df_sub_6)), ]
points_3 <- cbind(points_2,samples_df_sub_7)
colnames(points_3) <- c("PCoA Axis 1", "PCoA Axis 2", "z","a", "b", "c", "d", "e", "f", "g",colnames(samples_df_sub_7))
groups <- c("PCoA Axis 1", "PCoA Axis 2")

# Avoid comparing the same group with itself
points_sub <- points_3[, colnames(points_3) %in% groups]
points_sub$Condition <- samples_df_sub_7$Condition[match(row.names(points_sub), row.names(samples_df_sub_7))]

axis_x <- paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep="")
axis_y <- paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")

pcoa_plot <- ggplot(points_sub, aes_string(x = "`PCoA Axis 1`", 
                                           y = "`PCoA Axis 2`", 
                                           color = "Condition")) +
  geom_point(size = 3) +
  theme_classic() +
  scale_colour_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +  # Center title
  labs(title = "",
       x = axis_x, 
       y = axis_y)
pcoa_plot

pdf(paste(results.dir,"Figure_S16_Pseudomonadaceae_PC2_plant.pdf", sep=""), width=6, height=5)
print(pcoa_plot)
dev.off()

#Caulobacteraceae - Plant - AtSC - no dom 

#otu table
KO_SSC_only=read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Caulobacteraceae.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t", row.names =1) #make the SampleID column into the row.names
colnames(samples_df)[5]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object

#Sample subsetting
samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

#Subset for AtSC
samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum == "AtSC")

OTU_KO = otu_table(as.matrix(KO_SSC_only),taxa_are_rows = TRUE)
samples_sub = sample_data(samples_df_sub_3)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

bray_2 <- as.matrix(beta_isolate_KO)

str(samples_df_sub_3)
str(bray_2)

#Bind metadata with distance matrix
pcoa = cmdscale(bray_2, k=10, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z", "a", "b", "c", "d", "e", "f", "g") 
eig = pcoa$eig
points_2 <- points[order(row.names(points)), ]
samples_df_sub_6 <- samples_df_sub_3[row.names(samples_df_sub_3) %in% row.names(points),]
samples_df_sub_7 <- samples_df_sub_6[order(row.names(samples_df_sub_6)), ]
points_3 <- cbind(points_2,samples_df_sub_7)
colnames(points_3) <- c("PCoA Axis 1", "PCoA Axis 2", "z","a", "b", "c", "d", "e", "f", "g",colnames(samples_df_sub_7))
groups <- c("PCoA Axis 1", "PCoA Axis 2")

# Avoid comparing the same group with itself
points_sub <- points_3[, colnames(points_3) %in% groups]
points_sub$Condition <- samples_df_sub_7$Condition[match(row.names(points_sub), row.names(samples_df_sub_7))]

axis_x <- paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep="")
axis_y <- paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")

pcoa_plot <- ggplot(points_sub, aes_string(x = "`PCoA Axis 1`", 
                                           y = "`PCoA Axis 2`", 
                                           color = "Condition")) +
  geom_point(size = 3) +
  theme_classic() +
  scale_colour_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +  # Center title
  labs(title = "",
       x = axis_x, 
       y = axis_y)
pcoa_plot

pdf(paste(results.dir,"Figure_S16_Caulobacteraceae_PC2_plant.pdf", sep=""), width=6, height=5)
print(pcoa_plot)
dev.off()

#Burkholderiaceae - Plant - AtSC - no dom
#otu table
KO_SSC_only=read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Burkholderiaceae.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t", row.names =1) #make the SampleID column into the row.names
colnames(samples_df)[5]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object

#Sample subsetting
samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

#subset for AtSC
samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum == "AtSC")

OTU_KO = otu_table(as.matrix(KO_SSC_only),taxa_are_rows = TRUE)
samples_sub = sample_data(samples_df_sub_3)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

bray_2 <- as.matrix(beta_isolate_KO)

str(samples_df_sub_3)
str(bray_2)

#Bind metadata with distance matrix
pcoa = cmdscale(bray_2, k=10, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z", "a", "b", "c", "d", "e", "f", "g") 
eig = pcoa$eig
points_2 <- points[order(row.names(points)), ]
samples_df_sub_6 <- samples_df_sub_3[row.names(samples_df_sub_3) %in% row.names(points),]
samples_df_sub_7 <- samples_df_sub_6[order(row.names(samples_df_sub_6)), ]
points_3 <- cbind(points_2,samples_df_sub_7)
colnames(points_3) <- c("PCoA Axis 1", "PCoA Axis 2", "z","a", "b", "c", "d", "e", "f", "g",colnames(samples_df_sub_7))
groups <- c("PCoA Axis 1", "PCoA Axis 2")

# Avoid comparing the same group with itself
points_sub <- points_3[, colnames(points_3) %in% groups]
points_sub$Condition <- samples_df_sub_7$Condition[match(row.names(points_sub), row.names(samples_df_sub_7))]

axis_x <- paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep="")
axis_y <- paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")

pcoa_plot <- ggplot(points_sub, aes_string(x = "`PCoA Axis 1`", 
                                           y = "`PCoA Axis 2`", 
                                           color = "Condition")) +
  geom_point(size = 3) +
  theme_classic() +
  scale_colour_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +  # Center title
  labs(title = "",
       x = axis_x, 
       y = axis_y)
pcoa_plot

pdf(paste(results.dir,"Figure_S16_Burkholderiaceae_PC2_plant.pdf", sep=""), width=6, height=5)
print(pcoa_plot)
dev.off()


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

###Table S7 - DESeq2 overlap across SynComs per plant ======
#With dominators
input_table <- read.table(paste(working_directory, "/DESeq2/Sig_KO_all.txt", sep = ""), header=T, sep="\t")

plants <- c("Arabidopsis", "Barley", "Lotus")
SynComs <- c("AtSC", "HvSC", "LjSC", "SSC")

KO_table_2 <- data.frame()

for (plant in plants){
  for (syncom in SynComs){
    
    if (syncom == "AtSC"){
      SynComs2 <- c("HvSC","LjSC","SSC")
    } else if (syncom == "HvSC") {
      SynComs2 <- c("AtSC","LjSC","SSC")
    } else if (syncom == "LjSC") {
      SynComs2 <- c("AtSC","HvSC","SSC")
    } else {
      SynComs2 <- c("AtSC","HvSC","LjSC")
    }
    
    input_table_2 <- input_table[input_table$Plant == paste(plant),]
    input_table_3 <- input_table_2[input_table_2$SynCom == paste(syncom),]
    
    value1 <- length(input_table_3$KO)
    
    for (SC2 in SynComs2){
      norm_KO = read.table(paste(working_directory,"KO_tables/Original/", SC2, ".tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
      value2 <- length(row.names(norm_KO)[row.names(norm_KO) %in% input_table_3$KO])
      input_table_4 <- input_table_2[input_table_2$SynCom == paste(SC2),]
      value3 <- length(input_table_3$KO[input_table_3$KO %in% input_table_4$KO])
      
      new_data <- data.frame(t(data.frame(c(paste(plant), paste(syncom), value1, paste(SC2), value2, value3))))
      
      KO_table_2 <- rbind(KO_table_2, new_data)
    }
  }
}

row.names(KO_table_2) <- NULL
colnames(KO_table_2) <- c("Plant", "SynCom", "No_of_sig_KOs", "SynCom_comparison", "Overlap_present", "Overlap_significant")

write.table(KO_table_2, paste(results.dir, "Table_S7_KO_overlap_with_dom.txt", sep = ""), sep = "\t", row.names = T, col.names =T, quote =F)

#Deseq2 file without dominators
input_table <- read.table(paste(working_directory, "/DESeq2/Sig_KO_all_no_nod_rhizo.txt", sep = ""), header=T, sep="\t")

plants <- c("Arabidopsis", "Barley", "Lotus")
SynComs <- c("AtSC", "HvSC", "LjSC", "SSC")

KO_table_2 <- data.frame()

for (plant in plants){
  for (syncom in SynComs){
    
    if (syncom == "AtSC"){
      SynComs2 <- c("HvSC","LjSC","SSC")
    } else if (syncom == "HvSC") {
      SynComs2 <- c("AtSC","LjSC","SSC")
    } else if (syncom == "LjSC") {
      SynComs2 <- c("AtSC","HvSC","SSC")
    } else {
      SynComs2 <- c("AtSC","HvSC","LjSC")
    }
    
    input_table_2 <- input_table[input_table$plant == paste(plant),]
    input_table_3 <- input_table_2[input_table_2$SynCom == paste(syncom),]
    
    value1 <- length(input_table_3$KO)
    
    for (SC2 in SynComs2){
      norm_KO = read.table(paste(working_directory,"KO_tables/Original/", SC2, ".tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
      value2 <- length(row.names(norm_KO)[row.names(norm_KO) %in% input_table_3$KO])
      input_table_4 <- input_table_2[input_table_2$SynCom == paste(SC2),]
      value3 <- length(input_table_3$KO[input_table_3$KO %in% input_table_4$KO])
      
      new_data <- data.frame(t(data.frame(c(paste(plant), paste(syncom), value1, paste(SC2), value2, value3))))
      
      KO_table_2 <- rbind(KO_table_2, new_data)
    }
  }
}

row.names(KO_table_2) <- NULL
colnames(KO_table_2) <- c("Plant", "SynCom", "No_of_sig_KOs", "SynCom_comparison", "Overlap_present", "Overlap_significant")

write.table(KO_table_2, paste(results.dir, "Table_S7_KO_overlap_without_dom.txt", sep = ""), sep = "\t", row.names = T, col.names =T, quote =F)

###Figure 5b - Ternary plot - lenient selection =====
hop_4 <- read.table(paste(working_directory, "Functionality/852/ternary_KOs_av_med.txt", sep = ""), header =T, row.names =1, sep ="\t")

top <- read.table(paste(working_directory,"Annotations/pathway_top.txt", sep = ""), header=F, sep="\t")
KO_to_pathway <- read.table(paste(working_directory,"Annotations/KO_to_pathway.txt", sep = ""), header=T, sep="\t")
KO_to_pathway$V3 <- top$V2[match(KO_to_pathway$V2, top$V1)]

KO_to_pathway_2 <- read.table(paste(working_directory,"Annotations/KO_to_pathway_unannotated_2.txt", sep = ""), header=F, sep="\t")
colnames(KO_to_pathway_2) <- c("KO","new_category")

for (KO in KO_to_pathway_2$KO){
  KO_to_pathway$V3[KO_to_pathway$V1 == paste(KO)] <- KO_to_pathway_2$new_category[KO_to_pathway_2$KO == paste(KO)]
}

KO_to_pathway$V4 <- top$V3[match(KO_to_pathway$V3, top$V2)]

hop_4$pathway <- KO_to_pathway$V3[match(hop_4$KO, KO_to_pathway$V1)]
hop_4$category <- KO_to_pathway$V4[match(hop_4$KO, KO_to_pathway$V1)]

hop_4$pathway[is.na(hop_4$pathway)] <- "Unknown"
hop_4$category[is.na(hop_4$category)] <- "Unknown"

hop_4$Arabidopsis <- as.numeric(hop_4$Arabidopsis)
hop_4$Barley <- as.numeric(hop_4$Barley)
hop_4$Lotus <- as.numeric(hop_4$Lotus)
hop_4$Proportion_of_strains <- as.numeric(hop_4$Proportion_of_strains)

#Ternary plot
nv = 0.005
pn = position_nudge_tern(y=nv,x=-nv/2,z=-nv/2)

input_table <- read.table(paste(working_directory, "DESeq2/Sig_KO_all.txt", sep = ""), header=T, sep="\t")
input_table_2 <- table(input_table$KO)
input_table_3 <- names(input_table_2)[input_table_2 > 6]

hop_5 <- hop_4[hop_4$KO %in% input_table_3,]

ternary_dens <- ggplot(data=hop_5,aes(x=Arabidopsis,y=Barley, z=Lotus, size = Proportion_of_strains)) +
  geom_point(color = 'grey', alpha = 0.01) +
  coord_tern() +
  stat_density_tern(geom = 'polygon',n= 200, aes(fill = ..level.., alpha = ..level..), bins =750) +
  scale_fill_gradient(low = "blue",high = "red") +
  theme_bw()+
  guides(size = FALSE) +
  guides(fill=FALSE, alpha=FALSE) +
  ggtitle("") + theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  theme(text = element_text(size=25), axis.title.x = element_blank(), axis.text.x = element_blank() ) +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
ternary_dens

pdf(paste(results.dir,"Figure_5b_ternary_density.pdf", sep=""), width=12, height=12)
print(ternary_dens)
dev.off()

###Figure 5c - Ternary plot - host-specificity =====
hop_4 <- read.table(paste(working_directory, "Functionality/ternary_KOs_av_med_no_nod.txt", sep = ""), header =T, row.names =1, sep ="\t")

top <- read.table(paste(working_directory,"Annotations/pathway_top.txt", sep = ""), header=F, sep="\t")
KO_to_pathway <- read.table(paste(working_directory,"Annotations/KO_to_pathway.txt", sep = ""), header=T, sep="\t")
KO_to_pathway$V3 <- top$V2[match(KO_to_pathway$V2, top$V1)]

KO_to_pathway_2 <- read.table(paste(working_directory,"Annotations/KO_to_pathway_unannotated_2.txt", sep = ""), header=F, sep="\t")
colnames(KO_to_pathway_2) <- c("KO","new_category")

for (KO in KO_to_pathway_2$KO){
  KO_to_pathway$V3[KO_to_pathway$V1 == paste(KO)] <- KO_to_pathway_2$new_category[KO_to_pathway_2$KO == paste(KO)]
}

KO_to_pathway$V4 <- top$V3[match(KO_to_pathway$V3, top$V2)]

hop_4$pathway <- KO_to_pathway$V3[match(row.names(hop_4), KO_to_pathway$V1)]
hop_4$category <- KO_to_pathway$V4[match(row.names(hop_4), KO_to_pathway$V1)]

hop_4$pathway[is.na(hop_4$pathway)] <- "Unknown"
hop_4$category[is.na(hop_4$category)] <- "Unknown"

hop_4$Arabidopsis <- as.numeric(hop_4$Arabidopsis)
hop_4$Barley <- as.numeric(hop_4$Barley)
hop_4$Lotus <- as.numeric(hop_4$Lotus)
hop_4$Proportion_of_strains <- as.numeric(hop_4$Proportion_of_strains)

#Ternary plot
nv = 0.005
pn = position_nudge_tern(y=nv,x=-nv/2,z=-nv/2)

hex <- hue_pal()(length(unique(hop_4$category))) 

hop_4$Text <- "No"
hop_4$Plant <- "No"
hop_4$Text[hop_4$Arabidopsis > 3 & hop_4$Lotus < 3 & hop_4$Barley < 3] <- "Yes"
hop_4$Plant[hop_4$Arabidopsis > 3 & hop_4$Lotus < 3 & hop_4$Barley < 3] <- "Arabidopsis"

hop_4$Text[hop_4$Barley > 3 & hop_4$Lotus < 3 & hop_4$Arabidopsis < 3] <- "Yes"
hop_4$Plant[hop_4$Barley > 3 & hop_4$Lotus < 3 & hop_4$Arabidopsis < 3] <- "Barley"

hop_4$Text[hop_4$Lotus > 3 & hop_4$Arabidopsis < 3 & hop_4$Barley < 3] <- "Yes"
hop_4$Plant[hop_4$Lotus > 3 & hop_4$Arabidopsis < 3 & hop_4$Barley < 3] <- "Lotus"

hop_6 <- hop_4[hop_4$Text == "Yes",]

hop_Lj <- hop_6$KO[hop_6$Plant == "Lotus"]

#Rerun the same script but now with the original dataset
hop_4 <- read.table(paste(working_directory, "Functionality/852/ternary_KOs_av_med.txt", sep = ""), header =T, row.names =1, sep ="\t")

top <- read.table(paste(working_directory,"Annotations/pathway_top.txt", sep = ""), header=F, sep="\t")
KO_to_pathway <- read.table(paste(working_directory,"Annotations/KO_to_pathway.txt", sep = ""), header=T, sep="\t")
KO_to_pathway$V3 <- top$V2[match(KO_to_pathway$V2, top$V1)]

KO_to_pathway_2 <- read.table(paste(working_directory,"Annotations/KO_to_pathway_unannotated_2.txt", sep = ""), header=F, sep="\t")
colnames(KO_to_pathway_2) <- c("KO","new_category")

for (KO in KO_to_pathway_2$KO){
  KO_to_pathway$V3[KO_to_pathway$V1 == paste(KO)] <- KO_to_pathway_2$new_category[KO_to_pathway_2$KO == paste(KO)]
}

KO_to_pathway$V4 <- top$V3[match(KO_to_pathway$V3, top$V2)]

hop_4$pathway <- KO_to_pathway$V3[match(hop_4$KO, KO_to_pathway$V1)]
hop_4$category <- KO_to_pathway$V4[match(hop_4$KO, KO_to_pathway$V1)]

hop_4$pathway[is.na(hop_4$pathway)] <- "Unknown"
hop_4$category[is.na(hop_4$category)] <- "Unknown"

hop_4$Arabidopsis <- as.numeric(hop_4$Arabidopsis)
hop_4$Barley <- as.numeric(hop_4$Barley)
hop_4$Lotus <- as.numeric(hop_4$Lotus)
hop_4$Proportion_of_strains <- as.numeric(hop_4$Proportion_of_strains)

#Ternary plot
nv = 0.005
pn = position_nudge_tern(y=nv,x=-nv/2,z=-nv/2)

hex <- hue_pal()(length(unique(hop_4$category))) 

hop_4$Text <- "No"
hop_4$Plant <- "No"
hop_4$Text[hop_4$Arabidopsis > 3 & hop_4$Lotus < 3 & hop_4$Barley < 3] <- "Yes"
hop_4$Plant[hop_4$Arabidopsis > 3 & hop_4$Lotus < 3 & hop_4$Barley < 3] <- "Arabidopsis"

hop_4$Text[hop_4$Barley > 3 & hop_4$Lotus < 3 & hop_4$Arabidopsis < 3] <- "Yes"
hop_4$Plant[hop_4$Barley > 3 & hop_4$Lotus < 3 & hop_4$Arabidopsis < 3] <- "Barley"

hop_4$Text[hop_4$Lotus > 3 & hop_4$Arabidopsis < 3 & hop_4$Barley < 3] <- "Yes"
hop_4$Plant[hop_4$Lotus > 3 & hop_4$Arabidopsis < 3 & hop_4$Barley < 3] <- "Lotus"

hop_6 <- hop_4[hop_4$Text == "Yes",]
hop_6$Nodulator_associated <- "No"
hop_6$Nodulator_associated[hop_6$Plant == "Lotus" & !hop_6$KO %in% hop_Lj] <- "Yes"

ternary_ps <- ggplot(data=hop_6,aes(x=Arabidopsis,y=Barley, z=Lotus, size = Proportion_of_strains, color = Plant, shape = Nodulator_associated)) +
  geom_point() +
  coord_tern() +
  theme_bw()+
  labs(x = "", y = "", z = "") + # Removing axis names
  scale_color_manual(values = c("#1b9e77","#d95f02", "#7570b3")) +
  scale_shape_manual(values = c(16, 21)) +
  guides(size = FALSE) +
  guides(color = FALSE) +
  ggtitle("") + theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  labs(x = "Arabidopsis", y = "Barley", z = "Lotus",color = "Category", size = "Proportion of strains with KO", shape = "Nodulator KOs") +
  theme(text = element_text(size=25)) + theme(legend.text=element_text(size=25))+
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
ternary_ps

pdf(paste(results.dir,"Figure_5c_ternary_plant_spec.pdf", sep=""), width=12, height=12)
print(ternary_ps)
dev.off()

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

###Figure S19,S20,S21,S22,S23 - Creation of files necessary to create the small plots =====
hop_4 <- read.table(paste(working_directory,"Shiny_app/ternary_KOs_av_med_all.txt", sep = ""), header =T, sep ="\t", row.names = 1)

tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv", sep = ""),header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)

Genes_with_cassettes <- read.table(paste(working_directory,"Functionality/Genes_with_cassettes_2.tsv", sep = ""), sep = "\t", header =T)

#List of genes
gene_selection <- unique(Genes_with_cassettes$Gene)

samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep = ""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
samples_df_2$Condition[samples_df_2$Condition == "At"] <- "Arabidopsis"
samples_df_2$Condition[samples_df_2$Condition == "Lj"] <- "Lotus"
samples_df_2$Condition[samples_df_2$Condition == "Hv"] <- "Barley"

samples_df_3 <- samples_df_2[samples_df_2$Compartment == "ES",]
samples_df_4 <- samples_df_3[samples_df_3$Inoculum != "NS",]

SynComs <- c("AtSC","HvSC","LjSC", "SSC")

KO_table <- read.table(paste(working_directory,"KO_genome/KO_SSC.tsv", sep = ""), header=T, sep = "\t", row.names =1)
colnames(KO_table) <- gsub("X", "", colnames(KO_table))

together_2 <- data.frame()
fams <- data.frame()

for (gene in gene_selection){
  Genes_with_cassettes_2 <- Genes_with_cassettes$KO[Genes_with_cassettes$Gene == paste(gene)]
  
  if(gene == "bch"){
    Genes_with_cassettes_2 <- Genes_with_cassettes_2[!Genes_with_cassettes_2 %in% c("K13604","K13605", "K13601", "K13603", "K13602", "K04034", "K04396")]
  }
  
  hop_4_sub <- hop_4[hop_4$KO %in% Genes_with_cassettes_2,]
  
  for (syncom in SynComs) {
    norm_KO = read.table(paste(working_directory,"KO_tables/Original/", syncom, ".tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
    norm_iso = read.table(paste(working_directory,"Isolate_tables/Original/", syncom,"_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
    
    if (syncom == "SSC"){
      KO_table_2 <- KO_table
    } else {
      KO_table_2 <- KO_table[,colnames(KO_table) %in% row.names(tax_df_2)[tax_df_2$SynCom == paste(syncom)]]
    }
    
    KO_table_3 <- KO_table_2[row.names(KO_table_2) %in% Genes_with_cassettes_2,]
    KO_table_4 <- colSums(KO_table_3)
    KO_table_5 <- na.omit(names(KO_table_4)[KO_table_4 != 0])
    
    samples_df_5 <- samples_df_4[samples_df_4$Inoculum == paste(syncom),]
    samples_df_8 <- samples_df_5[!grepl("HL_orig", row.names(samples_df_5)),]
    
    #For family
    Genes_with_cassettes_plant <- unique(Genes_with_cassettes$Plant[Genes_with_cassettes$Gene == paste(gene)])
    samples_df_6 <- samples_df_5[samples_df_5$Condition == paste(Genes_with_cassettes_plant),]
    samples_df_7 <- samples_df_6[!grepl("HL_orig", row.names(samples_df_6)),]
    
    norm_KO_2 <- norm_KO[,colnames(norm_KO) %in% row.names(samples_df_7)]
    
    norm_KO_3 <- t(t(norm_KO_2)/rowSums(t(norm_KO_2)))
    norm_KO_4 <- norm_KO_3[row.names(norm_KO_3) %in% Genes_with_cassettes_2,]
    
    if (length(row.names(KO_table_3)) != 0){
      norm_iso_5 <- norm_iso[,colnames(norm_iso) %in% row.names(samples_df_7)]
      norm_iso_6 <- t(t(norm_iso_5)/rowSums(t(norm_iso_5)))
      norm_iso_7 <- norm_iso_6[row.names(norm_iso_6) %in% KO_table_5,]
      families <- unique(tax_df_2$family[row.names(tax_df_2) %in% KO_table_5])
      
      norm_iso_sub <- data.frame()
      if(length(families) > 0){
        for (fam in families){
          if (length(KO_table_5) == 1){
            norm_iso_8 <- norm_iso_7
          } else {
            norm_iso_8 <- norm_iso_7[row.names(norm_iso_7) %in% KO_table_5[KO_table_5 %in% row.names(tax_df_2)[tax_df_2$family == paste(fam)]],]
          }
          
          if (length(KO_table_5[KO_table_5 %in% row.names(tax_df_2)[tax_df_2$family == paste(fam)]]) == 1){
            norm_iso_9 <- norm_iso_8
          } else {
            norm_iso_9 <- colSums(norm_iso_8)
          }
          
          fam_value <- sum(norm_iso_9)/length(norm_iso_9)
          fam_data <- t(data.frame(c(paste(fam), fam_value, paste(syncom))))
          norm_iso_sub <- rbind(norm_iso_sub, fam_data)
        } 
        
        norm_iso_sub$V4 <- paste(gene)
        fams <- rbind(fams, norm_iso_sub)
      }
      
      if (length(row.names(KO_table_3)) != 0){
        norm_iso_10 <- norm_iso[,colnames(norm_iso) %in% row.names(samples_df_8)]
        norm_iso_11 <- t(t(norm_iso_10)/rowSums(t(norm_iso_10)))
        norm_iso_12 <- norm_iso_11[row.names(norm_iso_11) %in% KO_table_5,]
      } 
      
      Groups <- c("Arabidopsis", "Barley", "Lotus")
      
      norm_KO_2 <- norm_KO[,colnames(norm_KO) %in% row.names(samples_df_8)]
      
      norm_KO_3 <- t(t(norm_KO_2)/rowSums(t(norm_KO_2)))
      norm_KO_4 <- norm_KO_3[row.names(norm_KO_3) %in% Genes_with_cassettes_2,]
      
      for (group in Groups){
        
        samples_df_9 <- samples_df_8[samples_df_8$Condition == paste(group),]
        
        if (gene == "bch" & syncom == "LjSC"){
          norm_KO_5 <- norm_KO_4[names(norm_KO_4) %in% row.names(samples_df_9)]
        } else {
          norm_KO_5 <- norm_KO_4[,colnames(norm_KO_4) %in% row.names(samples_df_9)]
        }
        if (length(Genes_with_cassettes_2) == 1 | gene == "bch" & syncom == "LjSC"){
          norm_KO_5 <- norm_KO_4[names(norm_KO_4) %in% row.names(samples_df_9)]
          norm_KO_6 <- norm_KO_5
        } else {
          norm_KO_5 <- norm_KO_4[,colnames(norm_KO_4) %in% row.names(samples_df_9)]
          norm_KO_6 <- colSums(norm_KO_5)
        }
        
        if (syncom == "SSC"){
          if (length(KO_table_5) > 1){
            At_iso <- KO_table_5[KO_table_5 %in% row.names(tax_df_2)[tax_df_2$SynCom == "AtSC"]] 
            Lj_iso <- KO_table_5[KO_table_5 %in% row.names(tax_df_2)[tax_df_2$SynCom == "LjSC"]] 
            Hv_iso <- KO_table_5[KO_table_5 %in% row.names(tax_df_2)[tax_df_2$SynCom == "HvSC"]] 
            
            norm_iso_12_At <- norm_iso_12[row.names(norm_iso_12) %in% At_iso,]
            norm_iso_12_Lj <- norm_iso_12[row.names(norm_iso_12) %in% Lj_iso,]
            norm_iso_12_Hv <- norm_iso_12[row.names(norm_iso_12) %in% Hv_iso,]
            
            norm_iso_13_At <- colSums(norm_iso_12_At)
            norm_iso_13_Lj <- colSums(norm_iso_12_Lj)
            norm_iso_13_Hv <- colSums(norm_iso_12_Hv)
            
            norm_iso_14_At <- norm_iso_13_At[names(norm_iso_13_At) %in% row.names(samples_df_9)]
            norm_iso_14_Lj <- norm_iso_13_Lj[names(norm_iso_13_Lj) %in% row.names(samples_df_9)]
            norm_iso_14_Hv <- norm_iso_13_Hv[names(norm_iso_13_Hv) %in% row.names(samples_df_9)]
            
            value_iso_At <- sum(norm_iso_14_At)/length(names(norm_iso_14_At))
            value_iso_Lj <- sum(norm_iso_14_Lj)/length(names(norm_iso_14_Lj))
            value_iso_Hv <- sum(norm_iso_14_Hv)/length(names(norm_iso_14_Hv))
            
          } else if (length(KO_table_5) == 1) {
            plant_sel <- tax_df_2$SynCom[row.names(tax_df_2) == paste(KO_table_5)]
            norm_iso_14 <- norm_iso_12[names(norm_iso_12) %in% row.names(samples_df_9)]
            
            if (plant_sel == "AtSC"){
              value_iso_At <- sum(norm_iso_14)/length(names(norm_iso_14))
              value_iso_Lj <- 0
              value_iso_Hv <- 0
            } else if (plant_sel == "HvSC"){
              value_iso_At <- 0
              value_iso_Lj <- 0
              value_iso_Hv <- sum(norm_iso_14)/length(names(norm_iso_14))
            } else if (plant_sel == "LjSC"){
              value_iso_At <- 0
              value_iso_Lj <- sum(norm_iso_14)/length(names(norm_iso_14))
              value_iso_Hv <- 0
            }
          } else {
            value_iso_Hv <- 0
            value_iso_At <- 0
            value_iso_Lj <- 0
          }
        } else {
          if (length(KO_table_5) > 1){
            norm_iso_13 <- colSums(norm_iso_12)
            norm_iso_14 <- norm_iso_13[names(norm_iso_13) %in% row.names(samples_df_9)]
            value_iso <- sum(norm_iso_14)/length(names(norm_iso_14))
          } else if (length(KO_table_5) == 1) {
            norm_iso_14 <- norm_iso_12[names(norm_iso_12) %in% row.names(samples_df_9)]
            value_iso <- sum(norm_iso_14)/length(names(norm_iso_14))
          } else {
            value_iso <- 0
          }
        }
        
        norm_KO_7 <- norm_KO_6[names(norm_KO_6) %in% row.names(samples_df_9)]
        
        value <- sum(norm_KO_7)/length(names(norm_KO_7))
        
        if (syncom == "SSC"){
          together_3 <- t(data.frame(c(paste(gene), paste(syncom), paste(group), as.numeric(value), "AtSC", as.numeric(value_iso_At))))
          together_4 <- t(data.frame(c(paste(gene), paste(syncom), paste(group), as.numeric(value), "HvSC", as.numeric(value_iso_Hv))))
          together_5 <- t(data.frame(c(paste(gene), paste(syncom), paste(group), as.numeric(value), "LjSC", as.numeric(value_iso_Lj))))
          together <- rbind(together_3, together_4, together_5)
        } else {
          together <- t(data.frame(c(paste(gene), paste(syncom), paste(group), as.numeric(value), paste(syncom), as.numeric(value_iso))))
        }
        
        row.names(together) <- NULL
        
        together_2 <- rbind(together_2, together)
      }
    }
  }
}

row.names(together_2) <- NULL
colnames(together_2) <- c("Gene", "Inoculum", "Plant", "RA_KO", "Origin", "RA_Iso")

together_2$RA_KO[together_2$RA_KO == NaN] <- 0
together_2$RA_Iso[together_2$RA_Iso == NaN] <- 0

row.names(fams) <- NULL
colnames(fams) <- c("Family", "RA", "Inoculum", "Gene")

#Bar plots - Isolates
for (gene in gene_selection){
  together_2_sub <- together_2[together_2$Gene == paste(gene),]
  together_2_sub$RA_KO <- as.numeric(together_2_sub$RA_KO)
  
  together_2_sub$RA_Iso <- as.numeric(together_2_sub$RA_Iso)
  
  together_2_sub$New_column <- paste(together_2_sub$Inoculum, together_2_sub$Origin, sep = "-")
  together_2_sub$New_column <- gsub("HvSC-HvSC", "HvSC", together_2_sub$New_column)
  together_2_sub$New_column <- gsub("AtSC-AtSC", "AtSC", together_2_sub$New_column)
  together_2_sub$New_column <- gsub("LjSC-LjSC", "LjSC", together_2_sub$New_column)
  together_2_sub$New_column <- factor(together_2_sub$New_column, levels = c("AtSC", "HvSC", "LjSC", "SSC-AtSC", "SSC-HvSC", "SSC-LjSC"))
  
  g1 <- ggplot(together_2_sub %>% filter(Inoculum != "SSC"),
               aes(x= Plant, weight=RA_Iso, fill=New_column)) +
    theme_classic() +
    geom_bar(position = "dodge", width=0.5, just = 0.5) + 
    labs(x="Plant") +
    ylim(0,1)+
    scale_fill_manual(values = c("#A3A500","#00B0F6","#00BF7D")) +
    ggtitle(paste(gene)) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    ylab("Relative abundance - Isolate") + 
    xlab("Plant") +
    labs(fill = "Inoculum") +
    theme(axis.text.x = element_text(size = 14, hjust =0.2), axis.title.y = element_text(size = 18),axis.title.x = element_blank(), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=24))
  g1
  
  g2 <- g1 + geom_bar(data=together_2_sub %>% filter(Inoculum == "SSC"),
                      aes(x=Plant, fill=New_column),
                      position=position_stacknudge(x = 0.335), width=0.17) +
    scale_fill_manual(values = c("#A3A500","#00B0F6","#00BF7D","#fcddd9", "#fabab3","#F8766D"))
  
  pdf(paste(results.dir_2,"RA_", gene,".pdf", sep=""), width=8, height=4)
  print(g2)
  dev.off()
}

#Small ternaries - isolate fold change
hop_4 <- read.table(paste(working_directory,"Shiny_app/ternary_KOs_av_med_all.txt", sep = ""), header =T, sep ="\t", row.names = 1)
Genes_with_cassettes <- read.table(paste(working_directory,"Functionality/Genes_with_cassettes_2.tsv", sep = ""), sep = "\t", header =T)

KO_table <- read.table(paste(working_directory,"KO_genome/KO_SSC.tsv", sep = ""), header=T, sep = "\t", row.names =1)
colnames(KO_table) <- gsub("X", "", colnames(KO_table))

tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv", sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)

#List of genes
for (gene in Genes_with_cassettes$Gene){
  Genes_with_cassettes_KOs <- Genes_with_cassettes$KO[Genes_with_cassettes$Gene == paste(gene)]
  
  if(gene == "bch"){
    Genes_with_cassettes_KOs <- Genes_with_cassettes_KOs[Genes_with_cassettes_KOs != "K13604"]
  }
  
  hop_2 <- hop_4[hop_4$KO %in% Genes_with_cassettes_KOs,]
  
  table_gene_2 <-data.frame()
  for (syncom in SynComs){
    hop_sub_2 <- hop_2[hop_2$SynCom == paste(syncom),]
    At_val <- sum(hop_sub_2$Arabidopsis)/length(hop_sub_2$Arabidopsis)
    Hv_val <- sum(hop_sub_2$Barley)/length(hop_sub_2$Barley)
    Lj_val <- sum(hop_sub_2$Lotus)/length(hop_sub_2$Lotus)
    Prop_val <- sum(hop_sub_2$Proportion_of_strains)/length(hop_sub_2$Proportion_of_strains)
    
    table_gene <- t(data.frame(c(paste(gene), At_val, Hv_val, Lj_val, Prop_val, paste(syncom))))
    table_gene_2 <- rbind(table_gene_2, table_gene)
  }
  
  row.names(table_gene_2) <- NULL
  colnames(table_gene_2) <- c("gene", "Arabidopsis", "Barley", "Lotus", "Proportion_of_strains", "SynCom")
  
  table_gene_2$Arabidopsis <- as.numeric(table_gene_2$Arabidopsis)
  table_gene_2$Barley <- as.numeric(table_gene_2$Barley)
  table_gene_2$Lotus <- as.numeric(table_gene_2$Lotus)
  table_gene_2$Proportion_of_strains <- as.numeric(table_gene_2$Proportion_of_strains)
  nv = 0.005
  pn = position_nudge_tern(y=nv,x=-nv/2,z=-nv/2)
  
  ternary <- ggtern(data=table_gene_2,aes(x=Arabidopsis,y=Barley, z=Lotus, color = SynCom)) +
    geom_point(size = 6) +
    theme_bw()+
    scale_color_manual(values =c("#A3A500","#00B0F6","#00BF7D","#F8766D") ) +
    ggtitle(paste(gene)) + 
    theme(plot.title = element_text(hjust = 0.5, size = 20)) + 
    labs(color = "Gene") +
    theme(text = element_text(size=18)) + theme(legend.text=element_text(size=16)) +
    theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
  ternary 
  
  pdf(paste(results.dir_2,"ternary_", gene,".pdf", sep=""), width=8, height=6)
  print(ternary)
  dev.off()
}

#Family PieDonut plot
Fam_colors <- data.frame(unique(tax_df_2$family))
colnames(Fam_colors) <- "Family"
hex <- hue_pal()(length(Fam_colors$Family)) 
Fam_colors$Colors <- hex

source(paste(working_directory, "PieDonutCustom_fams_GS.R", sep = ""))

for (gene in unique(Genes_with_cassettes$Gene)){
  fams_2 <- fams[fams$Gene == paste(gene),]
  fams_3 <- fams_2 %>% dplyr::select (-Gene)
  fams_4 <- fams_3[c(3,1,2)]
  fams_4$RA <- as.numeric(fams_4$RA)
  
  fams_sub_2 <- fams_4
  fams_sub_2$Color <- Fam_colors$Colors[match(fams_sub_2$Family, Fam_colors$Family)]
  
  fams_2$Rel_RA_2 <- round(as.numeric(fams_2$RA)/sum(as.numeric(fams_2$RA))*10000,0)
  fams_2$Combination <- paste(fams_2$Inoculum, fams_2$Family,sep = "_")
  
  pie_data_2 <- data.frame()
  
  for (combi in fams_2$Combination){
    new <- fams_2$Rel_RA_2[fams_2$Combination == paste(combi)]
    syncom <- fams_2$Inoculum[fams_2$Combination == paste(combi)]
    family <- fams_2$Family[fams_2$Combination == paste(combi)]
    
    for (i in 1:new){
      pie_data <- data.frame(paste(syncom), paste(family))
      pie_data_2 <- rbind(pie_data_2, pie_data)
    }
  }
  
  colnames(pie_data_2) <- c("Inoculum", "Family")
  
  fams_sub_2$Combination <- paste(fams_sub_2$Inoculum, fams_sub_2$Family, sep = "_")
  fams_sub_3 <- fams_sub_2[order(fams_sub_2$Combination),]
  colors <- fams_sub_3$Color
  
  SynCom_colors <- data.frame(c("AtSC", "HvSC", "LjSC", "SSC"),c("#A3A500","#00B0F6","#00BF7D","#F8766D"))
  colnames(SynCom_colors) <- c("Inoculum", "Colour")
  SynCom_colors_2 <- SynCom_colors$Colour[SynCom_colors$Inoculum %in% unique(fams_2$Inoculum)]
  
  print(PieDonutCustom_fams(pie_data_2,aes(pies=Inoculum,donuts=Family),showRatioThreshold = 0.02))
  
  pdf(paste(results.dir_2,"pie_", gene,".pdf", sep=""), width=8, height=8)
  print(PieDonutCustom_fams(pie_data_2,aes(pies=Inoculum,donuts=Family),showRatioThreshold = 0.02))
  dev.off()
}

#Recreating the file for individual bar plots
hop_4 <- read.table(paste(working_directory,"Shiny_app/ternary_KOs_av_med_all.txt", sep = ""), header =T, sep ="\t", row.names = 1)

tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv", sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)

Genes_with_cassettes <- read.table(paste(working_directory,"Functionality/Genes_with_cassettes_2.tsv", sep = ""), sep = "\t", header =T)
#List of genes
gene_selection <- unique(Genes_with_cassettes$Gene)

samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv",sep = ""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
samples_df_2$Condition[samples_df_2$Condition == "At"] <- "Arabidopsis"
samples_df_2$Condition[samples_df_2$Condition == "Lj"] <- "Lotus"
samples_df_2$Condition[samples_df_2$Condition == "Hv"] <- "Barley"

samples_df_3 <- samples_df_2[samples_df_2$Compartment == "ES",]
samples_df_4 <- samples_df_3[samples_df_3$Inoculum != "NS",]

SynComs <- c("AtSC","HvSC","LjSC", "SSC")

KO_table <- read.table(paste(working_directory,"KO_genome/KO_SSC.tsv", sep = ""), header=T, sep = "\t", row.names =1)
colnames(KO_table) <- gsub("X", "", colnames(KO_table))

together_2 <- data.frame()
KO_table_sub_6 <- data.frame()

for (gene in gene_selection){
  Genes_with_cassettes_2 <- Genes_with_cassettes$KO[Genes_with_cassettes$Gene == paste(gene)]
  
  for (syncom in SynComs) {
    norm_KO = read.table(paste(working_directory,"KO_tables/Original/", syncom, ".tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
    norm_iso = read.table(paste(working_directory,"Isolate_tables/Original/", syncom, "_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
    
    if (syncom == "SSC"){
      KO_table_2 <- KO_table
    } else {
      KO_table_2 <- KO_table[,colnames(KO_table) %in% row.names(tax_df_2)[tax_df_2$SynCom == paste(syncom)]]
    }
    
    for (KO in Genes_with_cassettes_2) {
      KO_table_3 <- KO_table_2[row.names(KO_table_2) == paste(KO),]
      KO_table_5 <- na.omit(names(KO_table_3)[KO_table_3 != 0])
      
      samples_df_5 <- samples_df_4[samples_df_4$Inoculum == paste(syncom),]
      samples_df_8 <- samples_df_5[!grepl("HL_orig", row.names(samples_df_5)),]
      
      if (length(row.names(KO_table_3)) != 0){
        norm_iso_10 <- norm_iso[,colnames(norm_iso) %in% row.names(samples_df_8)]
        norm_iso_11 <- t(t(norm_iso_10)/rowSums(t(norm_iso_10)))
        norm_iso_12 <- norm_iso_11[row.names(norm_iso_11) %in% KO_table_5,]
      } 
      
      Groups <- c("Arabidopsis", "Barley", "Lotus")
      
      norm_KO_2 <- norm_KO[,colnames(norm_KO) %in% row.names(samples_df_8)]
      
      norm_KO_3 <- t(t(norm_KO_2)/rowSums(t(norm_KO_2)))
      norm_KO_4 <- norm_KO_3[row.names(norm_KO_3) == paste(KO),]
      
      for (group in Groups){
        
        samples_df_9 <- samples_df_8[samples_df_8$Condition == paste(group),]
        
        norm_KO_5 <- norm_KO_4[names(norm_KO_4) %in% row.names(samples_df_9)]
        
        if (length(KO_table_5) > 1){
          norm_iso_13 <- colSums(norm_iso_12)
          norm_iso_14 <- norm_iso_13[names(norm_iso_13) %in% row.names(samples_df_9)]
          value_iso <- sum(norm_iso_14)/length(names(norm_iso_14))
        } else if (length(KO_table_5) == 1) {
          norm_iso_14 <- norm_iso_12[names(norm_iso_12) %in% row.names(samples_df_9)]
          value_iso <- sum(norm_iso_14)/length(names(norm_iso_14))
        } else {
          value_iso <- 0
        }
        
        norm_KO_7 <- norm_KO_5[names(norm_KO_5) %in% row.names(samples_df_9)]
        
        if(length(row.names(KO_table_3)) > 0) {
          value <- sum(norm_KO_7)/length(names(norm_KO_7))
        } else {
          value <- 0
        }
        
        together <- t(data.frame(c(paste(gene), paste(KO),paste(syncom), paste(group), as.numeric(value), as.numeric(value_iso))))
        row.names(together) <- NULL
        
        together_2 <- rbind(together_2, together)
      }
    }
  }
  KO_table_sub <- KO_table[row.names(KO_table) %in% Genes_with_cassettes_2,]
  KO_table_sub[KO_table_sub > 0] <- 1
  KO_table_sub_2 <- colSums(KO_table_sub)
  KO_table_sub_3 <- KO_table_sub_2[order(KO_table_sub_2, decreasing =T)]
  max_value <- max(KO_table_sub_3)
  KO_table_sub_4 <- names(KO_table_sub_3)[KO_table_sub_3 == max_value]
  KO_table_sub_5 <- data.frame(KO_table_sub_4)
  colnames(KO_table_sub_5) <- "isolate"
  KO_table_sub_5$Gene <- paste(gene)
  
  KO_table_sub_6 <- rbind(KO_table_sub_6, KO_table_sub_5)
}

row.names(together_2) <- NULL
colnames(together_2) <- c("Gene", "KO", "SynCom", "Plant", "RA_KO", "RA_Iso")

together_2$RA_KO[together_2$RA_KO == NaN] <- 0
together_2$RA_Iso[together_2$RA_Iso == NaN] <- 0

together_2$Gene_cassette <- Genes_with_cassettes$Gene_cassette[match(together_2$KO, Genes_with_cassettes$KO)]
gene_sel_cassette <- unique(Genes_with_cassettes$Gene_cassette)

table <- read.table(paste(working_directory,"Functionality/Gene_viz_2.txt", sep = ""), header =T, sep = "\t")
together_2$Significance <- table$Significant[match(together_2$Gene_cassette,table$Gene_cassette )]

#Individual gene bar plots

plot_list <- list()

i <- 1

for (gene in gene_sel_cassette){
  together_2_sub <- together_2[together_2$Gene_cassette == paste(gene),]
  together_2_sub$RA_KO <- as.numeric(together_2_sub$RA_KO)
  
  together_2_sub <- na.omit(together_2_sub)
  
  if (length(together_2_sub$KO) > 0){
    if(unique(together_2_sub$Significance) != "Yes") {
      bar_plot_KO <- ggplot(together_2_sub, aes(fill=SynCom, y=RA_KO, x=Plant)) + 
        theme_classic() +
        geom_bar(position="stack", stat="identity") +  ggtitle(paste(gene)) + 
        theme(plot.title = element_text(hjust = 0.5)) + 
        ylab("Relative abundance - Gene") + 
        xlab("Host") +
        labs(fill = "Inoculum") +
        scale_fill_manual(values = c("#A3A500","#00B0F6","#00BF7D","#F8766D")) +
        theme(axis.text.x =  element_blank(), legend.position = "none",axis.title.y = element_blank(),axis.title.x = element_blank(), legend.title =  element_blank(), legend.text =  element_blank(), plot.title = element_text(size=24))
      bar_plot_KO
    } else {
      bar_plot_KO <- ggplot(together_2_sub, aes(fill=SynCom, y=RA_KO, x=Plant)) + 
        theme_classic() +
        geom_bar(position="stack", stat="identity") +  ggtitle(paste(gene)) + 
        theme(plot.title = element_text(hjust = 0.5)) + 
        ylab("Relative abundance - Gene") + 
        xlab("Host") +
        labs(fill = "Inoculum") +
        scale_fill_manual(values = c("#A3A500","#00B0F6","#00BF7D","#F8766D")) +
        theme(axis.text.x =  element_blank(), legend.position = "none",axis.title.y = element_blank(),axis.title.x = element_blank(), legend.title =  element_blank(), legend.text =  element_blank(), plot.title = element_text(size=24, colour = "red"))
      bar_plot_KO
    }
    
    plot_list[[i]] <- bar_plot_KO  
    i <- i + 1
    
    pdf(paste(results.dir_2,"RA_", gene,".pdf", sep=""), width=2, height=3)
    print(bar_plot_KO)
    dev.off()
  }
}

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

###Figure S25 - ABC transporter diversity violin plot =====

#getting list of isolates per plant-SynCom combination
list_of_isolates <- list()

Plants <- c("At", "Hv", "Lj")
SynComs <- c("AtSC","HvSC", "LjSC", "SSC")

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
    
    list_of_isolates[[paste(plant, syncom, sep = "_")]] <- norm_SSC_10
    
  }
}

#Getting data.frame of number of ABC transporter KOs per isolate
data_together <- data.frame()

table <- read.table(paste(working_directory, "/KO_genome/KO_SSC.tsv",sep = ""), sep= "\t", header =T, row.names =1) 
KO_gene <- read.table(paste(working_directory,"Annotations/Ann_ABC_2.txt", sep = ""),header=T, sep="\t")
table_2 <- table[row.names(table) %in% KO_gene$KO,]
colnames(table_2) <- gsub("X", "", colnames(table_2))
table_4 <- data.frame(colSums(table_2))
colnames(table_4) <- "Num_of_KOs"    

tax_df = read.table(paste(working_directory,"/SSC_taxonomy_GTDB.tsv", sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)

table_4$SynCom <- tax_df_2$SynCom[match(row.names(table_4) , row.names(tax_df_2))]
table_4$SynCom[is.na(table_4$SynCom)] <- "AtSC"

SynComs <- c("AtSC","HvSC", "LjSC", "SSC")
Plants <- c("Arabidopsis", "Barley", "Lotus")
Dominators <- c("LjNodule214", "P1_H10", "P2_A12", "P2_D6", "P2_G4")

table_sub_3 <- data.frame()

for (syncom in SynComs){
  if (syncom != "SSC"){
    table_sub <- table_4[table_4$SynCom == paste(syncom),]
  } else {
    table_sub <- table_4
  }
  
  table_sub$isolate <- row.names(table_sub)
  row.names(table_sub) <- NULL
  table_sub_2 <- table_sub[,c(3,1,2)]
  table_sub_2$Top <- "No"
  table_sub_2$Dominator <- "No"
  
  for (plant in Plants){
    if (plant == "Arabidopsis"){
      plant_2 <- "At"
    } else if (plant == "Barley"){
      plant_2 <- "Hv"
    } else {
      plant_2 <- "Lj"
    }
    
    selection <- paste(plant_2, syncom, sep = "_")
    
    list_of_isolates_2 <- list_of_isolates[[paste(selection)]]
    
    table_sub_2_5 <- table_sub_2[table_sub_2$isolate %in% list_of_isolates_2,]
    table_sub_2_6 <- table_sub_2[!table_sub_2$isolate %in% list_of_isolates_2,]
    
    table_sub_2_5$Top <- "Yes"
    table_sub_2_5$Dominator[table_sub_2_5$isolate %in% Dominators] <- "Yes"
    table_sub_2_6$Dominator[table_sub_2_6$isolate %in% Dominators] <- "Yes"
    table_sub_2_7 <- rbind(table_sub_2_5,table_sub_2_6)
    table_sub_2_7$Plant <- paste(plant)
    table_sub_2_7$SynCom_2 <- paste(syncom)
    table_sub_3 <- rbind(table_sub_3,table_sub_2_7)
  }
}

#Getting SynCom average
SynComs <- c("AtSC", "HvSC", "LjSC", "SSC")
KO_gene <- read.table(paste(working_directory,"Annotations/Ann_ABC_2.txt", sep = ""),header=T, sep="\t")

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
colnames(new_2) <- c("SynCom_2", "average")
new_2$average <- as.numeric(new_2$average)

#Numbers tables
Table_75 <- data.frame()

for (syncom in SynComs){
  for (plant in Plants){
    table_sub_sub <- table_sub_3[table_sub_3$SynCom_2 == paste(syncom),]
    table_sub_sub_2 <- table_sub_sub[table_sub_sub$Plant == paste(plant),]
    table_sub_sub_3 <- table(table_sub_sub_2$Top)
    table_75_2 <- data.frame(t(data.frame(c(paste(syncom), paste(plant), table_sub_sub_3[names(table_sub_sub_3) == "Yes"],table_sub_sub_3[names(table_sub_sub_3) == "No"]))))
    Table_75 <- rbind(Table_75,table_75_2)
  }
}

row.names(Table_75) <- NULL
colnames(Table_75) <- c("SynCom_2", "Plant", "Yes", "No")

Table_av <- data.frame()

for (syncom in SynComs){
  for (plant in Plants){
    table_sub_sub <- table_sub_3[table_sub_3$SynCom_2 == paste(syncom),]
    table_sub_sub_2 <- table_sub_sub[table_sub_sub$Plant == paste(plant),]
    table_sub_sub_3 <- table_sub_sub_2[table_sub_sub_2$Top == "Yes",]
    table_sub_sub_4 <- table_sub_sub_2[table_sub_sub_2$Top == "No",]
    
    av_value <- new_2$average[new_2$SynCom == paste(syncom)]
    above_val_yes <- length(table_sub_sub_3$isolate[table_sub_sub_3$Num_of_KOs > av_value])
    above_val_no <- length(table_sub_sub_4$isolate[table_sub_sub_4$Num_of_KOs > av_value])
    
    all_yes <- length(table_sub_sub_3$isolate)
    all_no <- length(table_sub_sub_4$isolate)
    
    Table_av_2 <- data.frame(t(data.frame(c(paste(syncom), paste(plant), paste(round(100*above_val_yes/all_yes,0), "%", sep = "") ,paste(round(100*above_val_no/all_no,0), "%", sep = "")))))
    Table_av <- rbind(Table_av,Table_av_2)
  }
}

row.names(Table_av) <- NULL
colnames(Table_av) <- c("SynCom_2", "Plant", "perc_yes", "perc_no")

Table_average <- data.frame(SynCom_2=c("AtSC", "HvSC", "LjSC", "SSC"), Plant = c("Arabidopsis", "Arabidopsis", "Arabidopsis", "Arabidopsis"),label = c("Above average: ", "", "Above average: ",""))
Table_average_2 <- data.frame(SynCom_2=c("AtSC", "HvSC", "LjSC", "SSC"), Plant = c("Arabidopsis", "Arabidopsis", "Arabidopsis", "Arabidopsis"),label = c("No of strains: ", "", "No of strains: ",""))

plot_violin <- ggplot(table_sub_3, aes(x=Plant, y=Num_of_KOs)) + 
  geom_violin(scale="width") +
  geom_jitter(data =table_sub_3[table_sub_3$Top == "No" & table_sub_3$Dominator == "No",], shape=16, position=position_jitter(0.2), colour = "black", size = 1.5, alpha =0.2,show.legend = F) +
  geom_jitter(data =table_sub_3[table_sub_3$Top == "Yes" & table_sub_3$Dominator == "Yes" & table_sub_3$isolate == "P2_G4" & table_sub_3$Plant != "Lotus",], shape=17, position=position_jitter(0.2), colour = "red", size = 3.5, alpha =1,show.legend = F) +
  geom_jitter(data =table_sub_3[table_sub_3$Top == "Yes" & table_sub_3$Dominator == "Yes" & table_sub_3$Plant == "Lotus",], shape=17, position=position_jitter(0.2), colour = "red", size = 3.5, alpha =1,show.legend = F) +
  geom_jitter(data =table_sub_3[table_sub_3$Top == "Yes" & table_sub_3$Dominator == "No",], shape=16, position=position_jitter(0.2), aes(colour = Plant), size = 2.5, alpha =1,show.legend = T) +
  ylab("No of ABC transporter KOs") + 
  xlab("Plant") +
  scale_colour_manual(values = c("#1b9e77","#d95f02","#7570b3")) +
  theme(axis.text.x = element_text(size = 14,angle = 25, hjust =1), axis.title.y = element_text(size = 18), axis.title.x = element_blank(),axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=24)) +
  ggtitle("ABC transporter distribution") + 
  geom_text_repel(aes(label=ifelse(Dominator == "Yes" & Plant == "Lotus", yes = as.character(isolate), '')),size=4,max.overlaps = Inf) +
  geom_text_repel(aes(label=ifelse(Dominator == "Yes" & SynCom == "HvSC" & Plant == "Arabidopsis" & isolate == "P2_G4", yes = as.character(isolate), '')),size=4,max.overlaps = Inf) +
  geom_text_repel(aes(label=ifelse(Dominator == "Yes" & SynCom == "HvSC" & Plant == "Barley" & isolate == "P2_G4", yes = as.character(isolate), '')),size=4,max.overlaps = Inf) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) +
  theme(strip.text.x = element_text(size = 14))+
  guides(colour=guide_legend(title="Top 75% root colonizers")) +
  geom_hline(data = new_2, aes(yintercept = average)) +
  facet_wrap(~SynCom_2,scales = "free_x") +
  geom_text(data=Table_75,aes(x = Plant, y = -2, label=No),color="darkgrey", hjust = 1.2)+
  geom_text(data=Table_75,aes(x = Plant, y = -2, label=Yes,color= Plant), hjust = -0.25)+
  geom_text(data=Table_av,aes(x = Plant, y = 525, label=perc_no),color="darkgrey", hjust = 1.2)+
  geom_text(data=Table_av,aes(x = Plant, y = 525, label=perc_yes,color= Plant), hjust = -0.25) +
  geom_text(data=Table_average,label=Table_average$label,y = 525, hjust = 1.35) +
  geom_text(data=Table_average_2,label=Table_average_2$label,y = -2, hjust = 1.35) +
  theme_classic()

plot_violin

#This extra plot is produced to have the text enlarged
pdf(paste(results.dir, "Figure_S25_violin_plot_ABC_text.pdf", sep=""), width=22, height= 12)
print(plot_violin)
dev.off()

pdf(paste(results.dir, "Figure_S25_violin_plot_ABC.pdf", sep=""), width=12, height= 7)
print(plot_violin)
dev.off()

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

###Figure 5e - Ternary plot - strict selection =====
hop_4 <- read.table(paste(working_directory, "Functionality/266/ternary_KOs_av_med_no_dom_266.txt", sep = ""), header =T, row.names =1, sep ="\t")

top <- read.table(paste(working_directory,"Annotations/pathway_top.txt", sep = ""), header=F, sep="\t")
KO_to_pathway <- read.table(paste(working_directory,"Annotations/KO_to_pathway.txt", sep = ""), header=T, sep="\t")
KO_to_pathway$V3 <- top$V2[match(KO_to_pathway$V2, top$V1)]

KO_to_pathway_2 <- read.table(paste(working_directory,"Annotations/KO_to_pathway_unannotated_2.txt", sep = ""), header=F, sep="\t")
colnames(KO_to_pathway_2) <- c("KO","new_category")

for (KO in KO_to_pathway_2$KO){
  KO_to_pathway$V3[KO_to_pathway$V1 == paste(KO)] <- KO_to_pathway_2$new_category[KO_to_pathway_2$KO == paste(KO)]
}

KO_to_pathway$V4 <- top$V3[match(KO_to_pathway$V3, top$V2)]

hop_4$pathway <- KO_to_pathway$V3[match(hop_4$KO, KO_to_pathway$V1)]
hop_4$category <- KO_to_pathway$V4[match(hop_4$KO, KO_to_pathway$V1)]

hop_4$pathway[is.na(hop_4$pathway)] <- "Unknown"
hop_4$category[is.na(hop_4$category)] <- "Unknown"

hop_4$Arabidopsis <- as.numeric(hop_4$Arabidopsis)
hop_4$Barley <- as.numeric(hop_4$Barley)
hop_4$Lotus <- as.numeric(hop_4$Lotus)
hop_4$Proportion_of_strains <- as.numeric(hop_4$Proportion_of_strains)

#Ternary plot
nv = 0.005
pn = position_nudge_tern(y=nv,x=-nv/2,z=-nv/2)

hex <- hue_pal()(length(unique(hop_4$category))) 

hop_4$Text <- "No"
hop_4$Text[hop_4$Arabidopsis > hop_4$Lotus & hop_4$Arabidopsis > hop_4$Barley] <- "Yes"
hop_4$Text[hop_4$Barley > hop_4$Lotus & hop_4$Barley > hop_4$Arabidopsis] <- "Yes"
hop_4$Text[hop_4$Lotus > 2 * hop_4$Barley & hop_4$Lotus > 2 *hop_4$Arabidopsis] <- "Yes"

input_table <- read.table(paste(working_directory, "DESeq2/Sig_KO_all_no_nod_rhizo.txt", sep = ""), header=T, sep="\t")
input_table_2 <- table(input_table$KO)
input_table_3 <- names(input_table_2)[input_table_2 == 12]

ternary_dens <- ggplot(data=hop_4,aes(x=Arabidopsis,y=Barley, z=Lotus, size = Proportion_of_strains)) +
  geom_point(color = 'grey', alpha = 0.01) +
  coord_tern() +
  stat_density_tern(geom = 'polygon',n= 200, aes(fill = ..level.., alpha = ..level..), bins =750) +
  scale_fill_gradient(low = "blue",high = "red") +
  theme_bw()+
  guides(size = FALSE) +
  guides(fill=FALSE, alpha=FALSE) +
  scale_color_manual(values = hex) +
  ggtitle("") + theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  theme(text = element_text(size=25), axis.title.x = element_blank(), axis.text.x = element_blank() ) +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
ternary_dens

pdf(paste(results.dir,"Figure_5e_ternary_density_266.pdf", sep=""), width=12, height=12)
print(ternary_dens)
dev.off()

###Figure S27b - Fold change curve =====
hop_4 <- read.table(paste(working_directory,"Functionality/266/boxplots.txt", sep = ""), header = T, sep = "\t", row.names =1)

colnames(hop_4) <- c("Plant","RA", "Have", "No_Have", "SynCom", "pathway", "Category_2", "No_of_strains")
hop_4$Category <- paste(hop_4$pathway, " (n% = ",format(round(as.numeric(hop_4$No_of_strains*100),1), nsmall = 1), ")", sep = "")

hop_4$Log <- log10(hop_4$Have)

hop_6 <- hop_4[,colnames(hop_4) != "No_Have"]
hop_7 <- hop_4[,colnames(hop_4) != "Have"]

colnames(hop_6) <- c("Plant","RA", "FC", "SynCom", "pathway", "Category_2","No_of_strains","Category","unused")
colnames(hop_7) <- c("Plant", "RA", "FC", "SynCom", "pathway", "Category_2","No_of_strains","Category","unused")

hop_6$Present <-"Present"
hop_7$Present <-"Absent"

hop_8 <- rbind(hop_6, hop_7)

hop_8$Present <- factor(hop_8$Present, levels = c("Present", "Absent"))
hop_9 <- na.omit(hop_8)

list_of_cats <- unique(hop_9$pathway)
new_scat_2 <- data.frame()

# Transform fold changes into LFC
hop_9 <- hop_9 %>%
  mutate(value = ifelse(FC == 0, 0, log2(FC)))

for (cat in list_of_cats){
  hop_sub <- hop_9[hop_9$pathway == paste(cat),]
  hop_sub$value <- as.numeric(hop_sub$value)
  
  stat <- compare_means(value~Present, hop_sub, method = "wilcox.test")
  
  # Calculate average Log2foldchange plant vs input of absent and present subpopulations
  Yes_LFC <- sum(hop_sub$value[hop_sub$Present == "Present"])/length(hop_sub$value[hop_sub$Present == "Present"])
  No_LFC <- sum(hop_sub$value[hop_sub$Present == "Absent"])/length(hop_sub$value[hop_sub$Present == "Absent"])
  
  Category_2 <- unique(hop_sub$Category_2)
  description <- unique(hop_sub$pathway)
  No_of_strains <- unique(hop_sub$No_of_strains)
  
  Average_RA <- sum(hop_sub$RA[hop_sub$Present == "Present"])/length(hop_sub$RA[hop_sub$Present == "Present"])
  
  new_scat <- t(data.frame(c(Category_2, description, round(as.numeric(No_of_strains),1),Average_RA, Yes_LFC, No_LFC, stat$p.format)))
  new_scat_2 <- rbind(new_scat_2, new_scat)
  
}

row.names(new_scat_2) <- NULL
colnames(new_scat_2) <- c("Category", "Pathway", "No_of_strains", "Average_RA", "LFC_vs_input_pathway_PRESENT_subpop", "LFC_vs_input_pathway_ABSENT_subpop", "Raw p-values")
new_scat_2$`Raw p-values` <-  as.numeric(new_scat_2$`Raw p-values`)

# Applying multiple adjustment methods for comparison
new_scat_2$padj_bonf = p.adjust(new_scat_2$`Raw p-values`, method = "bonferroni")  

#  Table manipulation, sorting by lowest pvals, rounding LFCs
new_scat_3 <- new_scat_2[order(new_scat_2$padj_bonf, decreasing =F),]
new_scat_3$LFC_vs_input_pathway_PRESENT_subpop <- round(as.numeric(new_scat_3$LFC_vs_input_pathway_PRESENT_subpop), 2)
new_scat_3$`LFC_vs_input_pathway_ABSENT_subpop` <- round(as.numeric(new_scat_3$`LFC_vs_input_pathway_ABSENT_subpop`), 2)

new_scat_4 <- new_scat_3[new_scat_3$LFC_vs_input_pathway_PRESENT_subpop > new_scat_3$`LFC_vs_input_pathway_ABSENT_subpop`,]
new_scat_5 <- new_scat_4[new_scat_4$padj_bonf < 0.05,]
new_scat_5$No_of_strains <- as.numeric(new_scat_5$No_of_strains)
new_scat_5$Average_RA <- as.numeric(new_scat_5$Average_RA)

Encoding(new_scat_5$Pathway) <- 'latin1'

# Create a "Square" metrics that select the best tradeoff between RA and fold change
new_scat_5$Square <- sqrt(((new_scat_5$Average_RA/max(new_scat_5$Average_RA)) * (new_scat_5$Average_RA/max(new_scat_5$Average_RA))) + ((new_scat_5$LFC_vs_input_pathway_PRESENT_subpop/max(new_scat_5$LFC_vs_input_pathway_PRESENT_subpop)) * (new_scat_5$LFC_vs_input_pathway_PRESENT_subpop/max(new_scat_5$LFC_vs_input_pathway_PRESENT_subpop))))

#  Create my other metric to best select pathways of interest
new_scat_5$diff=new_scat_5$LFC_vs_input_pathway_PRESENT_subpop-new_scat_5$LFC_vs_input_pathway_ABSENT_subpop
new_scat_5$diff_RA=new_scat_5$diff*new_scat_5$No_of_strains

# Define a legend theme, with custom text size and key size
text_size <- 12
key_size <- 0

legend_theme <- theme(legend.text = element_text(size = text_size),  # Apply text size
                      legend.key.size = unit(key_size, "in"), # Apply key size
                      strip.background=element_rect(colour="gray50", size=0.3), # Change 'size' for thickness
                      axis.text=element_text(color="gray50"),
                      axis.line = element_line(color="gray50", size=0.3),
                      axis.line.x = element_line(color="gray50", size=0.3), 
                      axis.line.y = element_line(color="gray50", size=0.3),
                      axis.ticks.x = element_line(color="gray50", size=0.3),
                      axis.ticks.y =element_line(color="gray50", size=0.3),
                      axis.text.x = element_text(size = 14), 
                      axis.title = element_text(size = 18), 
                      axis.text.y = element_text(size = 14), 
                      legend.title = element_text(size = 18), 
)

# Get 'Category' finite number of unique values
categories <- unique(new_scat_5$Category)
# Generate a color palette
color_palette <- scales::hue_pal()(length(categories))

# Create a named vector where names are categories and values are colors
named_colors <- setNames(color_palette, categories)

# Add a new column to new_scat_5 that maps each category to its color
new_scat_5$Color <- named_colors[new_scat_5$Category]

# Snippet for data point selection using exponential decay function
# Exponential decay function
A <- 3  
B <- 0.93
exp_decay <- function(x) {
  A * exp(-(log(A) / B) * x)
}

# Calculate decay values
new_scat_5$decay_value <- exp_decay(new_scat_5$Average_RA)

# Identify points above the exponential decay line
new_scat_5$above_decay <- new_scat_5$diff > new_scat_5$decay_value

# Compute distance from decay line for filtering
new_scat_5$distance = new_scat_5$diff - new_scat_5$decay_value
new_scat_5 <- new_scat_5[order(-new_scat_5$diff),]

# Select top points
top_points <- new_scat_5[new_scat_5$above_decay, ]
top_points <- top_points[order(-top_points$diff), ]

# Define non-informative pathways to remove
non_informative_pathways <- c("Biosynthesis", "Thermogenesis", "RNA polymerase", "RNA protein")  # Add more as needed

# Filter these out from the top_points
top_points <- top_points[!top_points$Pathway %in% non_informative_pathways, ]

# First, make sure 'alphabet_capital' is defined and has enough letters
alphabet_capital <- LETTERS  # Using predefined R variable LETTERS for capital alphabets
# Number of letters you want to assign
num_letters_to_assign <- 30  # Adjust this number based on your specific need
top_points$dist_category <- NA

if ("dist_category" %in% names(top_points) && "diff" %in% names(top_points) && "decay_value" %in% names(top_points)) {
  # Get indices where the condition is TRUE
  indices_to_change <- which(top_points$diff > top_points$decay_value)
  
  # Only take as many as are available
  num_to_assign <- min(num_letters_to_assign, length(indices_to_change))
  
  # Subset the indices safely
  selected_indices <- indices_to_change[1:num_to_assign]
  
  # Assign letters only to those positions
  top_points$dist_category[selected_indices] <- as.character(alphabet_capital[1:num_to_assign])
  
} else {
  cat("One or more specified columns do not exist in the dataframe top_points")
}

bottom_points <- new_scat_5[!new_scat_5$above_decay, ]
bottom_points$dist_category <- NA

top_points <- rbind(top_points, bottom_points)

plot_cor <- ggscatter(top_points, x = "Average_RA", y = "diff", color = "Category", size = "No_of_strains") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))+
  ylab("Log2Foldchange difference Present vs Absent") + 
  xlab("Relative abundance") +
  theme(legend.position = "right") +
  guides(color = guide_legend(ncol = 1)) +   
  labs(color = "Category", size = "Proportion of strains") +
  legend_theme +
  geom_text_repel(aes(label = as.character(dist_category), color= as.factor(Category)), show.legend = FALSE, size = 5, max.overlaps = Inf) 
plot_cor

pdf(paste(results.dir,"Figure_S27b_General_plot_exp_decay.pdf", sep=""), width=15, height=7)
print(plot_cor)
dev.off()


###Figure 5f & S28 - Fold change box plot =====
hop_4 <- read.table(paste(working_directory,"Functionality/266/boxplots.txt", sep = ""), header = T, sep = "\t", row.names =1)

colnames(hop_4) <- c("Plant","RA", "Have", "No_Have", "SynCom", "pathway", "Category_2", "No_of_strains")
hop_4$Category <- paste(hop_4$pathway, " (n% = ",format(round(as.numeric(hop_4$No_of_strains*100),1), nsmall = 1), ")", sep = "")

hop_4$Log <- log10(hop_4$Have)

hop_6 <- hop_4[,colnames(hop_4) != "No_Have"]
hop_7 <- hop_4[,colnames(hop_4) != "Have"]

colnames(hop_6) <- c("Plant","RA", "FC", "SynCom", "pathway", "Category_2","No_of_strains","Category","unused")
colnames(hop_7) <- c("Plant", "RA", "FC", "SynCom", "pathway", "Category_2","No_of_strains","Category","unused")

hop_6$Present <-"Present"
hop_7$Present <-"Absent"

hop_8 <- rbind(hop_6, hop_7)

hop_8$Present <- factor(hop_8$Present, levels = c("Present", "Absent"))
hop_9 <- na.omit(hop_8)

list_of_cats <- unique(hop_9$pathway)
new_scat_2 <- data.frame()

# Transform fold changes into LFC
hop_9 <- hop_9 %>%
  mutate(value = ifelse(FC == 0, 0, log2(FC)))


for (cat in list_of_cats){
  
  # cat="DNA replication"
  hop_sub <- hop_9[hop_9$pathway == paste(cat),]
  hop_sub$value <- as.numeric(hop_sub$value)
  
  
  stat <- compare_means(value~Present, hop_sub, method = "wilcox.test")
  
  # Calculate average Log2foldchange plant vs input of absent and present subpopulations
  Yes_LFC <- sum(hop_sub$value[hop_sub$Present == "Present"])/length(hop_sub$value[hop_sub$Present == "Present"])
  No_LFC <- sum(hop_sub$value[hop_sub$Present == "Absent"])/length(hop_sub$value[hop_sub$Present == "Absent"])
  
  Category_2 <- unique(hop_sub$Category_2)
  description <- unique(hop_sub$pathway)
  No_of_strains <- unique(hop_sub$No_of_strains)
  
  Average_RA <- sum(hop_sub$RA[hop_sub$Present == "Present"])/length(hop_sub$RA[hop_sub$Present == "Present"])
  
  
  new_scat <- t(data.frame(c(Category_2, description, round(as.numeric(No_of_strains),1),Average_RA, Yes_LFC, No_LFC, stat$p.format)))
  new_scat_2 <- rbind(new_scat_2, new_scat)
  
}


row.names(new_scat_2) <- NULL
colnames(new_scat_2) <- c("Category", "Pathway", "No_of_strains", "Average_RA", "LFC_vs_input_pathway_PRESENT_subpop", "LFC_vs_input_pathway_ABSENT_subpop", "Raw p-values")
new_scat_2$`Raw p-values` <-  as.numeric(new_scat_2$`Raw p-values`)

# Applying multiple adjustment methods for comparison
new_scat_2$padj_bonf = p.adjust(new_scat_2$`Raw p-values`, method = "bonferroni")  
new_scat_2$padj_BH = p.adjust(new_scat_2$`Raw p-values`, method = "BH")  

#  Table manipulation, sorting by lowest pvals, rounding LFCs
new_scat_3 <- new_scat_2[order(new_scat_2$padj_BH, decreasing =F),]
new_scat_3$LFC_vs_input_pathway_PRESENT_subpop <- round(as.numeric(new_scat_3$LFC_vs_input_pathway_PRESENT_subpop), 2)
new_scat_3$`LFC_vs_input_pathway_ABSENT_subpop` <- round(as.numeric(new_scat_3$`LFC_vs_input_pathway_ABSENT_subpop`), 2)

new_scat_4 <- new_scat_3[new_scat_3$LFC_vs_input_pathway_PRESENT_subpop > new_scat_3$`LFC_vs_input_pathway_ABSENT_subpop`,]
new_scat_5 <- new_scat_4[new_scat_4$padj_bonf < 0.05,]
new_scat_5$No_of_strains <- as.numeric(new_scat_5$No_of_strains)
new_scat_5$Average_RA <- as.numeric(new_scat_5$Average_RA)

Encoding(new_scat_5$Pathway) <- 'latin1'

# Create a "Square" metrics that select the best tradeoff between RA and fold change
new_scat_5$Square <- sqrt(((new_scat_5$Average_RA/max(new_scat_5$Average_RA)) * (new_scat_5$Average_RA/max(new_scat_5$Average_RA))) + ((new_scat_5$LFC_vs_input_pathway_PRESENT_subpop/max(new_scat_5$LFC_vs_input_pathway_PRESENT_subpop)) * (new_scat_5$LFC_vs_input_pathway_PRESENT_subpop/max(new_scat_5$LFC_vs_input_pathway_PRESENT_subpop))))

#  Create my other metric to best select pathways of interest
new_scat_5$diff=new_scat_5$LFC_vs_input_pathway_PRESENT_subpop-new_scat_5$LFC_vs_input_pathway_ABSENT_subpop
new_scat_5$diff_RA=new_scat_5$diff*new_scat_5$No_of_strains

# Define a legend theme, with custom text size and key size

text_size <- 12
key_size <- 0

legend_theme <- theme(legend.text = element_text(size = text_size),  # Apply text size
                      legend.key.size = unit(key_size, "in"), # Apply key size
                      strip.background=element_rect(colour="gray50", size=0.3), # Change 'size' for thickness
                      axis.text=element_text(color="gray50"),
                      axis.line = element_line(color="gray50", size=0.3),
                      axis.line.x = element_line(color="gray50", size=0.3), 
                      axis.line.y = element_line(color="gray50", size=0.3),
                      axis.ticks.x = element_line(color="gray50", size=0.3),
                      axis.ticks.y =element_line(color="gray50", size=0.3),
                      axis.text.x = element_text(size = 14), 
                      axis.title = element_text(size = 18), 
                      axis.text.y = element_text(size = 14), 
                      legend.title = element_text(size = 18), 
)


# Get 'Category' finite number of unique values
categories <- unique(new_scat_5$Category)
# Generate a color palette
color_palette <- scales::hue_pal()(length(categories))

# Create a named vector where names are categories and values are colors
named_colors <- setNames(color_palette, categories)

# Add a new column to new_scat_5 that maps each category to its color
new_scat_5$Color <- named_colors[new_scat_5$Category]

# Snippet for data point selection using exponential decay function
# Exponential decay function
A <- 2.8  
B <- 0.95
exp_decay <- function(x) {
  A * exp(-(log(A) / B) * x)
}

# Calculate decay values
new_scat_5$decay_value <- exp_decay(new_scat_5$Average_RA)

# Identify points above the exponential decay line
new_scat_5$above_decay <- new_scat_5$diff > new_scat_5$decay_value

# Compute distance from decay line for filtering
new_scat_5$distance = new_scat_5$diff - new_scat_5$decay_value
new_scat_5 <- new_scat_5[order(-new_scat_5$diff),]

# Select top points
top_points <- new_scat_5[new_scat_5$above_decay, ]
top_points <- top_points[order(-top_points$diff), ]
# top_20_points <- head(top_points, 30)

# Define non-informative pathways to remove
non_informative_pathways <- c("Biosynthesis", "Thermogenesis", "RNA polymerase", "RNA protein")  # Add more as needed

# Filter these out from the top_points
top_points <- top_points[!top_points$Pathway %in% non_informative_pathways, ]

# First, make sure 'alphabet_capital' is defined and has enough letters
alphabet_capital <- LETTERS  # Using predefined R variable LETTERS for capital alphabets
# Number of letters you want to assign
num_letters_to_assign <- 30  # Adjust this number based on your specific need
top_points$dist_category <- NA

if ("dist_category" %in% names(top_points) && "diff" %in% names(top_points) && "decay_value" %in% names(top_points)) {
  # Assign letters to 'dist_category' only where the condition is true
  top_points$dist_category[top_points$diff > top_points$decay_value][1:num_letters_to_assign] <- as.character(alphabet_capital[1:num_letters_to_assign])
} else {
  cat("One or more specified columns do not exist in the dataframe top_points")
}

new_scat_6 <- new_scat_5 %>%
  left_join(top_points %>% dplyr::select(Pathway, dist_category), by = "Pathway", suffix = c("", "_new"))

# Coalesce the categories to ensure that non-updated categories stay as they are
new_scat_6$dist_category <- coalesce(new_scat_6$dist_category_new, new_scat_6$dist_category)

# Clean up the temporary column
new_scat_6$dist_category_new <- NULL

#266 KO dataset - Figure S23f
hop_10=hop_9
# swr function 1 separate pathway name vs proportion
swr = function(string) {
  pos = regexpr("\\(", string) # Find the position of the first '('
  if (pos > 0) { # Check if '(' is found
    # Insert a newline before the first '('
    string = paste0(substr(string, 1, pos - 1), "\n", substr(string, pos, nchar(string)))
  }
  return(string)
}
swr = Vectorize(swr)

# Create line breaks
hop_10$Category = swr(hop_10$Category)

select_Y <- subset(new_scat_6, !is.na(dist_category))

# Category selection decreasing foldchange diff
pathway_selection_Y = na.omit(merge.data.frame(x = select_Y , y = hop_10, by.x = "Pathway", by.y = "pathway" , all.y = T))
pathway_selection_Y$padj_bonf <- sprintf("%.1e", pathway_selection_Y$padj_bonf)
pathway_selection_Y$dist_category_full= paste0(pathway_selection_Y$dist_category,": ",pathway_selection_Y$Category.y," (p=",pathway_selection_Y$padj_bonf,")")

# Create a theme object with your legend modifications
text_size=12
key_size=1
legend_theme <- theme(legend.text = element_text(size = text_size),  # Apply text size
                      legend.key.size = unit(key_size, "cm"), # Apply key size
                      strip.background=element_rect(colour="gray50", size=0.3), # Change 'size' for thickness
                      axis.text=element_text(color="black", size =text_size ),
                      axis.line = element_line(color="gray50", size=0.3),
                      axis.line.x = element_line(color="gray50", size=0.3), 
                      axis.line.y = element_line(color="gray50", size=0.3),
                      axis.ticks.x = element_line(color="gray50", size=0.3),
                      axis.ticks.y =element_line(color="gray50", size=0.3)
)

plot_top_sig_RA_hor <- ggplot(pathway_selection_Y, aes(x=fct_rev(dist_category_full) , y=value, color = Present)) +
  geom_boxplot(outlier.shape = NA)+ theme_classic()+
  scale_shape_manual(values = c(5,2,1))+
  geom_jitter(aes(shape=Plant,color= SynCom ,size = 1), alpha=1, cex =1,position = position_jitter(width = 0.2))+
  theme(legend.position="left")+
  scale_fill_manual(values = named_colors) +
  scale_color_manual(values = c("black","gray70","#A3A500","#00B0F6","#00BF7D","#F8766D"))+
  ylab("Log2FoldChange vs Input") + 
  xlab("") +
  legend_theme+
  coord_flip()
plot_top_sig_RA_hor

width_cm <- 21
height_cm <- 29.7

# Convert dimensions to inches (1 inch = 2.54 cm)
width_in <- width_cm / 2.54
height_in <- height_cm / 2.54

pdf(paste(results.dir,"Figure_S28_General_plot_266.pdf", sep=""), width=width_in, height=height_in)
ggarrange(plot_top_sig_RA_hor+theme(legend.position = "none"))
dev.off()

#With legend
pdf(paste(results.dir,"Figure_S28_General_plot_266_legend.pdf", sep=""), width=width_in, height=height_in)
print(plot_top_sig_RA_hor)
dev.off()

pathway_selection_Y_main_sel <- unique(pathway_selection_Y$dist_category_full)[order(unique(pathway_selection_Y$dist_category_full))][1:10]

pathway_selection_Y_main <- pathway_selection_Y[pathway_selection_Y$dist_category_full %in% pathway_selection_Y_main_sel,]

plot_top_sig_RA_hor_main <- ggplot(pathway_selection_Y_main, aes(x=fct_rev(dist_category_full) , y=value, color = Present)) +
  geom_boxplot(outlier.shape = NA)+ theme_classic()+
  scale_shape_manual(values = c(5,2,1))+
  geom_jitter(aes(shape=Plant,color= SynCom ,size = 1), alpha=1, cex =1,position = position_jitter(width = 0.2))+
  theme(legend.position="left")+
  scale_fill_manual(values = named_colors) +
  scale_color_manual(values = c("black","gray70","#A3A500","#00B0F6","#00BF7D","#F8766D"))+
  ylab("Log2FoldChange vs Input") + 
  xlab("") +
  legend_theme+
  coord_flip()
plot_top_sig_RA_hor_main

width_cm <- 10
height_cm <- 15

# Convert dimensions to inches (1 inch = 2.54 cm)
width_in <- width_cm / 2.54
height_in <- height_cm / 2.54

pdf(paste(results.dir,"Figure_5f_General_plot_266.pdf", sep=""), width=width_in, height=height_in)
ggarrange(plot_top_sig_RA_hor_main+theme(legend.position = "none"))
dev.off()

#With legend
pdf(paste(results.dir,"Figure_5f_General_plot_266_legend.pdf", sep=""), width=width_in, height=height_in)
print(plot_top_sig_RA_hor_main)
dev.off()

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

###Figure S30a - Overview of 266 across families =====
df_90 <- read.table(paste(working_directory, "top90_isolates_no_dom.txt", sep = ""), sep = "\t",header =F)
df_90_2 <- unique(na.omit(unlist(as.vector(df_90))))

tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE, row.names =1)
top <- tax_df[row.names(tax_df) %in% df_90_2,]
non_top <- tax_df[!row.names(tax_df) %in% df_90_2,]

input_table <- read.table(paste(working_directory, "DESeq2/Sig_KO_all_no_nod_rhizo.txt", sep = ""), header=T, sep="\t")
input_table_2 <- table(input_table$KO)
input_table_3 <- names(input_table_2)[input_table_2 == 12]

table <- read.table(paste(working_directory, "KO_genome/KO_SSC.tsv", sep = ""), sep= "\t", header =T,row.names =1) 
colnames(table)[grep("M.16",colnames(table))] <- "M-16"
colnames(table)[grep("M.6",colnames(table))] <- "M-6"
colnames(table)[grep("M.10",colnames(table))] <- "M-10"
colnames(table)[grep("M.11_2",colnames(table))] <- "M-11_2"
colnames(table)[grep("^M.11$",colnames(table))] <- "M-11"
colnames(table) <- gsub("X", "", colnames(table))
table[table > 0] <- 1

table_2 <- table[row.names(table) %in% input_table_3,]

Tax_Fam_top <- table(top$family)
Tax_Fam_top_2 <- names(Tax_Fam_top)[Tax_Fam_top > 2]

top_iso <- top[top$family %in% Tax_Fam_top_2,]
No_top_iso <- non_top[non_top$family %in% Tax_Fam_top_2,]

table_3 <- data.frame()

for (family in Tax_Fam_top_2 ){
  top_iso_sub <- row.names(top_iso)[top_iso$family == paste(family)]
  table_sub <- table_2[,colnames(table_2) %in% top_iso_sub]
  table_sub_2 <- rowSums(table_sub)/length(colnames(table_sub))
  
  value_100 <- length(names(table_sub_2)[table_sub_2 >= 0.9])
  value_90 <- length(names(table_sub_2)[table_sub_2 >= 0.6 & table_sub_2 < 0.9])
  value_60 <- length(names(table_sub_2)[table_sub_2 >= 0.3 & table_sub_2 < 0.6])
  value_30 <- length(names(table_sub_2)[table_sub_2 > 0 & table_sub_2 < 0.3])
  value_0 <- length(names(table_sub_2)[table_sub_2 == 0])
  
  table_top <- data.frame(t(data.frame(c(paste(family), length(top_iso_sub), "Top", value_100, value_90, value_60, value_30, value_0))))
  
  No_top_iso_sub <- row.names(No_top_iso)[No_top_iso$family == paste(family)]
  table_sub <- table_2[,colnames(table_2) %in% No_top_iso_sub]
  table_sub_2 <- rowSums(table_sub)/length(colnames(table_sub))
  
  value_100 <- length(names(table_sub_2)[table_sub_2 >= 0.9])
  value_90 <- length(names(table_sub_2)[table_sub_2 >= 0.6 & table_sub_2 < 0.9])
  value_60 <- length(names(table_sub_2)[table_sub_2 >= 0.3 & table_sub_2 < 0.6])
  value_30 <- length(names(table_sub_2)[table_sub_2 > 0 & table_sub_2 < 0.3])
  value_0 <- length(names(table_sub_2)[table_sub_2 == 0])
  
  table_non_top <- data.frame(t(data.frame(c(paste(family),length(No_top_iso_sub), "No top", value_100, value_90, value_60, value_30, value_0))))
  
  table_3 <- rbind(table_3,table_top,table_non_top)
}

row.names(table_3) <- paste(table_3$X1, " (n = ",table_3$X2,") - ", table_3$X3, sep ="")
table_4 <- table_3[,c(4,5,6,7,8)]
colnames(table_4) <- c("100%-90%", "90%-60%","60%-30%","30%-1%", "0%")
table_4$`100%-90%` <- as.numeric(table_4$`100%-90%`)
table_4$`90%-60%` <- as.numeric(table_4$`90%-60%`)
table_4$`60%-30%` <- as.numeric(table_4$`60%-30%`)
table_4$`30%-1%` <- as.numeric(table_4$`30%-1%`)
table_4$`0%` <- as.numeric(table_4$`0%`)

table_4$Group <- row.names(table_4)
table_5 <- melt(table_4)

table_5$variable <- factor(table_5$variable, levels = c("0%","30%-1%","60%-30%","90%-60%","100%-90%"))

table_5$Group <- factor(table_5$Group, levels =c("Rhizobiaceae (n = 54) - Top", "Rhizobiaceae (n = 61) - No top", 
                                                 "Xanthobacteraceae (n = 12) - Top" ,"Xanthobacteraceae (n = 3) - No top",
                                                 "Beijerinckiaceae (n = 20) - Top", "Beijerinckiaceae (n = 26) - No top",
                                                 "Caulobacteraceae (n = 26) - Top", "Caulobacteraceae (n = 14) - No top",
                                                 "Cellulomonadaceae (n = 3) - Top","Cellulomonadaceae (n = 11) - No top", 
                                                 "Micrococcaceae (n = 8) - Top", "Micrococcaceae (n = 74) - No top",
                                                 "Sphingomonadaceae (n = 11) - Top", "Sphingomonadaceae (n = 3) - No top",
                                                 "Pseudomonadaceae (n = 42) - Top", "Pseudomonadaceae (n = 41) - No top",
                                                 "Microbacteriaceae (n = 6) - Top",  "Microbacteriaceae (n = 23) - No top",
                                                 "Sphingobacteriaceae (n = 11) - Top", "Sphingobacteriaceae (n = 21) - No top",
                                                 "Burkholderiaceae (n = 138) - Top","Burkholderiaceae (n = 129) - No top", 
                                                 "Rhodanobacteraceae (n = 6) - Top","Rhodanobacteraceae (n = 9) - No top",
                                                 "Xanthomonadaceae (n = 24) - Top", "Xanthomonadaceae (n = 40) - No top"))

g1 <- ggplot(table_5, aes(x= Group, weight=value, fill=variable)) +
  theme_classic() +
  geom_bar(position = "stack") + 
  scale_fill_manual(values = c("black","#d3d3d3","#fa8072","#ffd700","#66bd63")) +
  ggtitle("266 KO distribution in 90% CRA top vs non-top") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  ylab("No of KOs") + 
  xlab("Family - top vs non-top") +
  labs(fill = "Proportion of KO in Family") +
  theme(axis.text.x = element_text(size = 14, hjust =1, angle = 45), axis.title.y = element_text(size = 18),axis.title.x = element_blank(), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=24))
g1 

g2 <- ggarrange(NA, g1, ncol = 2, widths = c(1,6))

pdf(paste(results.dir, "Figure_S30a_overview.pdf", sep=""), width=20, height=7)
print(g2)
dev.off()


###Figure S31 - Niche replacement score ======
#With nodulators - plant-specific pathways
pathway_selection <- c("ABC transporters", "Quorum sensing", "Two-component system", "Purine metabolism", "Oxidative phosphorylation", "Pentose and glucuronate interconversions", "Porphyrin metabolism", "Galactose metabolism", "Cell cycle - Caulobacter", "Exopolysaccharide biosynthesis")

Families <- c("Burkholderiaceae", "Caulobacteraceae", "Pseudomonadaceae", "other", "all")
#Make figure
data_fam <- read.table(paste(working_directory, "LjSC_Family_drop_out_experiment/Figure_5d_validation_in_vivo_Fam_drop_out_with_nod.txt", sep = ""), header =T, sep = "\t")

data_fam$Family[data_fam$Family == "LDT1"] <- "Burkholderiaceae"
data_fam$Family[data_fam$Family == "LDT2"] <- "Caulobacteraceae"
data_fam$Family[data_fam$Family == "LDT3"] <- "Pseudomonadaceae"
data_fam$Family[data_fam$Family == "LDT4"] <- "Rhizobiaceae"
data_fam$Family[data_fam$Family == "LDT5"] <- "other"
data_fam$Family[data_fam$Family == "LDT6"] <- "No_drop_out"
data_fam$Family[data_fam$Family == "all"] <- "No_drop_out"

data_fam$Log <- log2(data_fam$Present_fold_change)

data_fam_3 <- data_fam[data_fam$Pathway %in% pathway_selection,]
data_fam_3$Family <- factor(data_fam_3$Family, levels = c("No_drop_out", "Burkholderiaceae", "Caulobacteraceae", "Pseudomonadaceae", "Rhizobiaceae", "other"))
groups_pathways <- unique(data_fam_3$Pathway)
groups_pathways_2 <- groups_pathways[order(groups_pathways,decreasing = T)]
data_fam_3$Pathway <- factor(data_fam_3$Pathway, levels = groups_pathways_2)

#Comparison - Rhizobiaceae removed as there was contamination in the experiment
Families <- c("Burkholderiaceae", "Caulobacteraceae", "Pseudomonadaceae", "other")
Pathway <- unique(data_fam_3$Pathway)

new_data <- data.frame()

for (fam in Families){
  for (path in Pathway){
    data_fam_sub <- data_fam_3[data_fam_3$Pathway == paste(path) & data_fam_3$Family == paste(fam),]
    data_sub <- data_fam_3[data_fam_3$Pathway == paste(path) & data_fam_3$Family == "No_drop_out",]
    
    value <- data_fam_sub$Present_fold_change/data_sub$Present_fold_change
    
    new_data_2 <- data.frame(t(data.frame(c(paste(path), paste(fam), value))))
    new_data <- rbind(new_data, new_data_2)
  }
}


row.names(new_data) <- NULL
colnames(new_data) <- c("Pathway", "Family", "Niche_replacement_score")

new_data$Family <- factor(new_data$Family, levels = c("Burkholderiaceae", "Caulobacteraceae", "Pseudomonadaceae", "other"))
new_data$Pathway <- factor(new_data$Pathway, levels = groups_pathways_2)
new_data$Log <- log2(as.numeric(new_data$Niche_replacement_score))

plot_2 <- new_data %>% ggplot() + geom_tile(aes(y=Pathway, x=Family, fill = as.numeric(Log))) +
  geom_text(aes(y=Pathway, x=Family, label = round(as.numeric(Log),2))) + 
  scale_fill_gradient2(low = "red", mid = "white", high = "green", midpoint = 0) +
  labs(y="Pathway", x="Dropped-out family", fill="", title="Niche replacement score") +
  theme_bw(base_size = 14) %+replace% theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  theme(plot.title = element_text(hjust =0.5))
plot_2

pdf(paste(results.dir,"Figure_S31a_Niche_replacement_with_dom.pdf", sep=""), width=7.5, height=6)
print(plot_2)
dev.off()

#Without nodulators - core pathways
pathway_selection <- c("Vitamin B6 metabolism","Transcriptional regulator","Secretion","Quorum sensing",
                       "Phenylalanine, tyrosine and tryptophan biosynthesis","Pantothenate and CoA biosynthesis",
                       "Oxidative phosphorylation","Membrane protein","Glycerophospholipid metabolism",
                       "Folate biosynthesis","Flagellar assembly","Exopolysaccharide biosynthesis","Cysteine and methionine metabolism",
                       "Arginine biosynthesis", "Methane metabolism", "Aminobenzoate degradation", "Ascorbate and aldarate metabolism")

Families <- c("Burkholderiaceae", "Caulobacteraceae", "Pseudomonadaceae", "other", "all")

data_fam <- read.table(paste(working_directory, "LjSC_Family_drop_out_experiment/Figure_5f_validation_in_vivo_Fam_drop_out_no_nod.txt", sep = ""), header =T, sep = "\t")

data_fam$Family[data_fam$Family == "LDT1"] <- "Burkholderiaceae"
data_fam$Family[data_fam$Family == "LDT2"] <- "Caulobacteraceae"
data_fam$Family[data_fam$Family == "LDT3"] <- "Pseudomonadaceae"
data_fam$Family[data_fam$Family == "LDT4"] <- "Rhizobiaceae"
data_fam$Family[data_fam$Family == "LDT5"] <- "other"
data_fam$Family[data_fam$Family == "LDT6"] <- "No_drop_out"

colnames(data_fam) <- c("RA", "Present_fold_change", "Absent_fold_change", "Family", "Pathway", "Percentage")

data_fam$Log <- log2(data_fam$Present_fold_change)

data_fam_3 <- data_fam[data_fam$Pathway %in% pathway_selection,]
data_fam_3$Family <- factor(data_fam_3$Family, levels = c("No_drop_out", "Burkholderiaceae", "Caulobacteraceae", "Pseudomonadaceae", "other"))
groups_pathways <- unique(data_fam_3$Pathway)
groups_pathways_2 <- groups_pathways[order(groups_pathways,decreasing = T)]
data_fam_3$Pathway <- factor(data_fam_3$Pathway, levels = groups_pathways_2)

#Comparison - Rhizobiaceae not taken along as there was contamination in those samples
Families <- c("Burkholderiaceae", "Caulobacteraceae", "Pseudomonadaceae", "other")
Pathway <- unique(data_fam_3$Pathway)

new_data <- data.frame()

for (fam in Families){
  for (path in Pathway){
    data_fam_sub <- na.omit(data_fam_3[data_fam_3$Pathway == paste(path) & data_fam_3$Family == paste(fam),])
    data_sub <- na.omit(data_fam_3[data_fam_3$Pathway == paste(path) & data_fam_3$Family == "No_drop_out",])
    
    value <- data_fam_sub$Present_fold_change/data_sub$Present_fold_change
    
    new_data_2 <- data.frame(t(data.frame(c(paste(path), paste(fam), value))))
    new_data <- rbind(new_data, new_data_2)
  }
}

row.names(new_data) <- NULL
colnames(new_data) <- c("Pathway", "Family", "Niche_replacement_score")

new_data$Family <- factor(new_data$Family, levels = c("Burkholderiaceae", "Caulobacteraceae", "Pseudomonadaceae", "other"))
new_data$Pathway <- factor(new_data$Pathway, levels = groups_pathways_2)
new_data$Log <- log2(as.numeric(new_data$Niche_replacement_score))

plot_2 <- new_data %>% ggplot() + geom_tile(aes(y=Pathway, x=Family, fill = as.numeric(Log))) +
  geom_text(aes(y=Pathway, x=Family, label = round(as.numeric(Log),2))) + 
  scale_fill_gradient2(low = "red", mid = "white", high = "green", midpoint = 0) +
  labs(y="Pathway", x="Dropped-out family", fill="", title="Niche replacement score") +
  theme_bw(base_size = 14) %+replace% theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  theme(plot.title = element_text(hjust =0.5))
plot_2

pdf(paste(results.dir,"Figure_S31b_Niche_replacement_no_dom.pdf", sep=""), width=9, height=8)
print(plot_2)
dev.off()


###Figure S32 - Bacterial/plant ratio =====
#With nodulators
plant_reads <- read.table(paste(working_directory, "LjSC_Family_drop_out_experiment/plant_reads.tsv", sep = ""), header =T, sep = "\t")
plant_reads_2 <- plant_reads[grep("ROOT", plant_reads$Sample),]

norm_SSC =read.table(paste(working_directory,"LjSC_Family_drop_out_experiment/Data_with_synthetic_input.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
norm_SSC_3 <- data.frame(colSums(norm_SSC))
colnames(norm_SSC_3) <- "No_of_bacterial_reads"
norm_SSC_3$Sample <- row.names(norm_SSC_3)
norm_SSC_4 <- norm_SSC_3[grep("ROOT", row.names(norm_SSC_3)),]
norm_SSC_4$No_of_plant_reads <- plant_reads_2$Host_mapped_Infr[match(norm_SSC_4$Sample, plant_reads_2$Sample)]
norm_SSC_4$Bacterial_read_per_plant <- norm_SSC_4$No_of_bacterial_reads/norm_SSC_4$No_of_plant_reads

norm_SSC_4_2 <- norm_SSC_4[!grepl("LDT7", norm_SSC_4$Sample),]
norm_SSC_5 <- norm_SSC_4_2[!grepl("LDT4", norm_SSC_4_2$Sample),]

norm_SSC_5$SynCom <- sapply(strsplit(as.vector(row.names(norm_SSC_5)), "_"),`[`, 1)

#Rhizobiaceae was not taken along as it was found to have a contamination
norm_SSC_5$SynCom <- gsub("LDT1", "Burkholderiaceae drop out", norm_SSC_5$SynCom)
norm_SSC_5$SynCom <- gsub("LDT2", "Caulobacteraceae drop out", norm_SSC_5$SynCom)
norm_SSC_5$SynCom <- gsub("LDT3", "Pseudomonadaceae drop out", norm_SSC_5$SynCom)
#norm_SSC_5$SynCom <- gsub("LDT4", "Rhizobiaceae drop out", norm_SSC_5$SynCom)
norm_SSC_5$SynCom <- gsub("LDT5", "All other families drop out", norm_SSC_5$SynCom)
norm_SSC_5$SynCom <- gsub("LDT6", "Full LjSC", norm_SSC_5$SynCom)

# Perform ANOVA
fitAnova <- aov(Bacterial_read_per_plant ~ SynCom, data = norm_SSC_5)

# Perform Tukey's post-hoc test
Tukey <- TukeyHSD(fitAnova)

# Get letters for significance
letters_anova <- multcompLetters4(fitAnova, Tukey)$SynCom$Letters

# Combine results into a data frame for plotting
ltlbl_combined <- data.frame(
  SynCom = names(letters_anova),
  Letters = letters_anova
)

# Ensure the order of the factor levels is correct
ltlbl_combined$SynCom <- factor(ltlbl_combined$SynCom, levels = unique(norm_SSC_5$SynCom))

# Apply the manually defined order to the factor
ltlbl_combined <- ltlbl_combined[order(ltlbl_combined$SynCom), ]

# Merge the letters with the data frame
norm_SSC_5 <- merge(norm_SSC_5, ltlbl_combined, by = "SynCom")

#Add number of 266 KOs
KO_table = read.table(paste(working_directory, "KO_genome/KO_LjSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

taxonomy <- read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",row.names =1)
taxonomy_LjSC <- taxonomy[taxonomy$SynCom == "LjSC",]

families <- c("Pseudomonadaceae", "Rhizobiaceae", "Burkholderiaceae", "Caulobacteraceae", "other")

input_table <- read.table(paste(working_directory,"DESeq2/Sig_KO_all_no_nod_rhizo.txt", sep = ""), header=T, sep="\t")
input_table_2 <- table(input_table$KO)
input_table_3 <- names(input_table_2)[input_table_2 == 12]

#No dominator dataset
new_fam_data <- data.frame()

for (group in families){
  if(group == "other"){
    new_tax <- row.names(taxonomy_LjSC)[!taxonomy_LjSC$family %in% families]
  } else {
    new_tax <- row.names(taxonomy_LjSC)[taxonomy_LjSC$family == paste(group)]
  }
  
  KO_table_sub <- KO_table[,colnames(KO_table) %in% new_tax]
  KO_table_sub_2 <- rowSums(KO_table_sub)
  total_KOs <- length(names(KO_table_sub_2)[KO_table_sub_2 != 0])
  Family_KOs <- names(KO_table_sub_2)[KO_table_sub_2 != 0]
  
  total_266_Fam <- length(Family_KOs[Family_KOs %in% input_table_3])
  
  all_data <- data.frame(t(data.frame(c(paste(group), total_266_Fam))))
  new_fam_data <- rbind(new_fam_data, all_data)
}

row.names(new_fam_data) <- NULL
colnames(new_fam_data) <- c("Family", "Total no of 266 KOs")

norm_SSC_5$Group <- NA
norm_SSC_5$Group[norm_SSC_5$SynCom == "All other families drop out"] <- new_fam_data$`Total no of 266 KOs`[new_fam_data$Family == "other"]
norm_SSC_5$Group[norm_SSC_5$SynCom == "Pseudomonadaceae drop out"] <- new_fam_data$`Total no of 266 KOs`[new_fam_data$Family == "Pseudomonadaceae"]
norm_SSC_5$Group[norm_SSC_5$SynCom == "Caulobacteraceae drop out"] <- new_fam_data$`Total no of 266 KOs`[new_fam_data$Family == "Caulobacteraceae"]
norm_SSC_5$Group[norm_SSC_5$SynCom == "Burkholderiaceae drop out"] <- new_fam_data$`Total no of 266 KOs`[new_fam_data$Family == "Burkholderiaceae"]
norm_SSC_5$Group[norm_SSC_5$SynCom == "Full LjSC"] <- 266

norm_SSC_5$Group_2 <- paste(norm_SSC_5$SynCom, " (no of KO: " , norm_SSC_5$Group, ")",sep = "")

norm_SSC_5$Group_2 <- factor(norm_SSC_5$Group_2, levels = c("Burkholderiaceae drop out (no of KO: 141)", "Caulobacteraceae drop out (no of KO: 130)", "Pseudomonadaceae drop out (no of KO: 84)", "All other families drop out (no of KO: 242)", "Full LjSC (no of KO: 266)"))
# Plot the results
plot_col <- ggplot(norm_SSC_5, aes(x = Group_2, y = Bacterial_read_per_plant, colour = Group_2)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  ylab("Bacterial read/plant read") +
  geom_jitter(shape = 16, position = position_jitter(0.2), aes(colour = Group_2), show.legend = TRUE) +
  theme(
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    strip.text.x = element_text(size = 10)
  ) +
  guides(colour = FALSE) +
  scale_y_continuous(labels = scales::scientific,  ) +
  stat_summary(aes(label = Letters, y = max(Bacterial_read_per_plant)*1.03), fun = max, geom = "text") +
  ggtitle("Bacterial colonization with nodulators")
plot_col

plot_col_2 <- ggarrange(NULL, plot_col, ncol =2, widths=c(1,5))

#Without nodulators

norm_SSC =read.table(paste(working_directory,"LjSC_Family_drop_out_experiment/Data_with_synthetic_input.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
norm_SSC_2 <- norm_SSC[!row.names(norm_SSC) %in% c("LjNodule214", "LjRoot234", "LjRoot228"),]
norm_SSC_3 <- data.frame(colSums(norm_SSC_2))
colnames(norm_SSC_3) <- "No_of_bacterial_reads"
norm_SSC_3$Sample <- row.names(norm_SSC_3)
norm_SSC_4 <- norm_SSC_3[grep("ROOT", row.names(norm_SSC_3)),]
norm_SSC_4$No_of_plant_reads <- plant_reads_2$Host_mapped_Infr[match(norm_SSC_4$Sample, plant_reads_2$Sample)]
norm_SSC_4$Bacterial_read_per_plant <- norm_SSC_4$No_of_bacterial_reads/norm_SSC_4$No_of_plant_reads

norm_SSC_4_2 <- norm_SSC_4[!grepl("LDT7", norm_SSC_4$Sample),]
norm_SSC_5 <- norm_SSC_4_2[!grepl("LDT4", norm_SSC_4_2$Sample),]

norm_SSC_5$SynCom <- sapply(strsplit(as.vector(row.names(norm_SSC_5)), "_"),`[`, 1)

#Rhizobiaceae was not taken along as it was found to have a contamination
norm_SSC_5$SynCom <- gsub("LDT1", "Burkholderiaceae drop out", norm_SSC_5$SynCom)
norm_SSC_5$SynCom <- gsub("LDT2", "Caulobacteraceae drop out", norm_SSC_5$SynCom)
norm_SSC_5$SynCom <- gsub("LDT3", "Pseudomonadaceae drop out", norm_SSC_5$SynCom)
#norm_SSC_5$SynCom <- gsub("LDT4", "Rhizobiaceae drop out", norm_SSC_5$SynCom)
norm_SSC_5$SynCom <- gsub("LDT5", "All other families drop out", norm_SSC_5$SynCom)
norm_SSC_5$SynCom <- gsub("LDT6", "Full LjSC", norm_SSC_5$SynCom)

# Perform ANOVA
fitAnova <- aov(Bacterial_read_per_plant ~ SynCom, data = norm_SSC_5)

# Perform Tukey's post-hoc test
Tukey <- TukeyHSD(fitAnova)

# Get letters for significance
letters_anova <- multcompLetters4(fitAnova, Tukey)$SynCom$Letters

# Combine results into a data frame for plotting
ltlbl_combined <- data.frame(
  SynCom = names(letters_anova),
  Letters = letters_anova
)

# Ensure the order of the factor levels is correct
ltlbl_combined$SynCom <- factor(ltlbl_combined$SynCom, levels = unique(norm_SSC_5$SynCom))

# Apply the manually defined order to the factor
ltlbl_combined <- ltlbl_combined[order(ltlbl_combined$SynCom), ]

# Merge the letters with the data frame
norm_SSC_5 <- merge(norm_SSC_5, ltlbl_combined, by = "SynCom")

#Add number of 266 KOs
KO_table = read.table(paste(working_directory, "KO_genome/KO_LjSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

taxonomy <- read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",row.names =1)
taxonomy_LjSC <- taxonomy[taxonomy$SynCom == "LjSC",]

families <- c("Pseudomonadaceae", "Rhizobiaceae", "Burkholderiaceae", "Caulobacteraceae", "other")

input_table <- read.table(paste(working_directory,"DESeq2/Sig_KO_all_no_nod_rhizo.txt", sep = ""), header=T, sep="\t")
input_table_2 <- table(input_table$KO)
input_table_3 <- names(input_table_2)[input_table_2 == 12]

#No dominator dataset
new_fam_data <- data.frame()

for (group in families){
  if(group == "other"){
    new_tax <- row.names(taxonomy_LjSC)[!taxonomy_LjSC$family %in% families]
  } else {
    new_tax <- row.names(taxonomy_LjSC)[taxonomy_LjSC$family == paste(group)]
  }
  
  KO_table_sub <- KO_table[,colnames(KO_table) %in% new_tax]
  KO_table_sub_2 <- rowSums(KO_table_sub)
  total_KOs <- length(names(KO_table_sub_2)[KO_table_sub_2 != 0])
  Family_KOs <- names(KO_table_sub_2)[KO_table_sub_2 != 0]
  
  total_266_Fam <- length(Family_KOs[Family_KOs %in% input_table_3])
  
  all_data <- data.frame(t(data.frame(c(paste(group), total_266_Fam))))
  new_fam_data <- rbind(new_fam_data, all_data)
}

row.names(new_fam_data) <- NULL
colnames(new_fam_data) <- c("Family", "Total no of 266 KOs")

norm_SSC_5$Group <- NA
norm_SSC_5$Group[norm_SSC_5$SynCom == "All other families drop out"] <- new_fam_data$`Total no of 266 KOs`[new_fam_data$Family == "other"]
norm_SSC_5$Group[norm_SSC_5$SynCom == "Pseudomonadaceae drop out"] <- new_fam_data$`Total no of 266 KOs`[new_fam_data$Family == "Pseudomonadaceae"]
norm_SSC_5$Group[norm_SSC_5$SynCom == "Caulobacteraceae drop out"] <- new_fam_data$`Total no of 266 KOs`[new_fam_data$Family == "Caulobacteraceae"]
norm_SSC_5$Group[norm_SSC_5$SynCom == "Burkholderiaceae drop out"] <- new_fam_data$`Total no of 266 KOs`[new_fam_data$Family == "Burkholderiaceae"]
norm_SSC_5$Group[norm_SSC_5$SynCom == "Full LjSC"] <- 266

norm_SSC_5$Group_2 <- paste(norm_SSC_5$SynCom, " (no of KO: " , norm_SSC_5$Group, ")",sep = "")

norm_SSC_5$Group_2 <- factor(norm_SSC_5$Group_2, levels = c("Burkholderiaceae drop out (no of KO: 141)", "Caulobacteraceae drop out (no of KO: 130)", "Pseudomonadaceae drop out (no of KO: 84)", "All other families drop out (no of KO: 242)", "Full LjSC (no of KO: 266)"))
# Plot the results
plot_col_3 <- ggplot(norm_SSC_5, aes(x = Group_2, y = Bacterial_read_per_plant, colour = Group_2)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  ylab("Bacterial read/plant read") +
  geom_jitter(shape = 16, position = position_jitter(0.2), aes(colour = Group_2), show.legend = TRUE) +
  theme(
    axis.text.x = element_text(size = 12, hjust = 1, angle = 25),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    strip.text.x = element_text(size = 10)
  ) +
  guides(colour = FALSE) +
  scale_y_continuous(labels = scales::scientific,  ) +
  stat_summary(aes(label = Letters, y = max(Bacterial_read_per_plant)*1.03), fun = max, geom = "text") +
  ggtitle("Bacterial colonization without nodulators")
plot_col_3

plot_col_4 <- ggarrange(NULL, plot_col_3, ncol =2, widths=c(1,5))

plot_col_5 <- ggarrange(plot_col_2, plot_col_4, nrow =2, ncol =1, heights=c(2,3))

pdf(paste(results.dir,"Figure_S32_Bac_col_vs_plant.pdf", sep=""), width=10, height=10)
print(plot_col_5)
dev.off()

###Table S8 - KO overview =====
KO_table = read.table(paste(working_directory, "KO_genome/KO_LjSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
sigtab_col_all_2 <- read.table(paste(working_directory, "LjSC_Family_drop_out_experiment/DESeq2_Root_vs_input_Fam_drop_no_nod.txt", sep = ""), header =T,  sep = "\t")
sigtab_col_all_2_with_nod <- read.table(paste(working_directory, "LjSC_Family_drop_out_experiment/DESeq2_Root_vs_input_Fam_drop.txt", sep = ""), header =T,  sep = "\t")

taxonomy <- read.table("/media/wolf1/InRoot/Gijs_Selten/SSC_R2_11_2022/SSC_Community_2024/SSC_taxonomy_GTDB.tsv", header =T, row.names =1)
taxonomy_LjSC <- taxonomy[taxonomy$SynCom == "LjSC",]

families <- c("Pseudomonadaceae", "Rhizobiaceae", "Burkholderiaceae", "Caulobacteraceae", "other")
Nodulators <- c("LjRoot228", "LjNodule214", "LjRoot234")

input_table <- read.table(paste(working_directory,"DESeq2/Sig_KO_all_no_nod_rhizo.txt", sep = ""), header=T, sep="\t")
input_table_2 <- table(input_table$KO)
input_table_3 <- names(input_table_2)[input_table_2 == 12]

#No nodulators dataset
new_fam_data <- data.frame()

for (group in families){
  if(group == "other"){
    new_tax <- row.names(taxonomy_LjSC)[!taxonomy_LjSC$family %in% families]
    new_tax <- new_tax[!new_tax %in% Nodulators]
    sigtab_col_all_2_sub <- sigtab_col_all_2[sigtab_col_all_2$Subset == "All other families drop out",]
  } else {
    new_tax <- row.names(taxonomy_LjSC)[taxonomy_LjSC$family == paste(group)]
    new_tax <- new_tax[!new_tax %in% Nodulators]
    sigtab_col_all_2_sub <- sigtab_col_all_2[sigtab_col_all_2$Subset == paste(group, "drop out", sep = " "),]
  }
  
  KO_table_sub <- KO_table[,colnames(KO_table) %in% new_tax]
  KO_table_sub_2 <- rowSums(KO_table_sub)
  total_KOs <- length(names(KO_table_sub_2)[KO_table_sub_2 != 0])
  Family_KOs <- names(KO_table_sub_2)[KO_table_sub_2 != 0]
  
  sigtab_col_all_2_full <- sigtab_col_all_2[sigtab_col_all_2_with_nod$Subset == "Full LjSC",]
  sigtab_col_all_2_full_2 <- sigtab_col_all_2_full[sigtab_col_all_2_full$padj < 0.05,]
  sigtab_col_all_2_full_3 <- sigtab_col_all_2_full_2[sigtab_col_all_2_full_2$log2FoldChange > 0,]
  contrib_full <- length(sigtab_col_all_2_full_3$KO[sigtab_col_all_2_full_3$KO %in% Family_KOs])
  
  sigtab_col_all_2_2 <- sigtab_col_all_2_sub[sigtab_col_all_2_sub$padj < 0.05,]
  sigtab_col_all_2_3 <- sigtab_col_all_2_2[sigtab_col_all_2_2$log2FoldChange > 0,]
  lost_KOs <- length(sigtab_col_all_2_full_3$KO[!sigtab_col_all_2_full_3$KO %in% sigtab_col_all_2_3$KO])
  
  total_266_Fam <- length(Family_KOs[Family_KOs %in% input_table_3])
  sigtab_col_all_2_full_3_sig <- sigtab_col_all_2_full_3$KO[!sigtab_col_all_2_full_3$KO %in% sigtab_col_all_2_3$KO]
  lost_266 <- length(sigtab_col_all_2_full_3_sig[sigtab_col_all_2_full_3_sig %in% input_table_3])
  
  all_data <- data.frame(t(data.frame(c(paste(group), total_KOs, contrib_full, lost_KOs, total_266_Fam, lost_266 ))))
  new_fam_data <- rbind(new_fam_data, all_data)
}

row.names(new_fam_data) <- NULL
colnames(new_fam_data) <- c("Family", "Total no of KOs", "Total sig KOs", "Lost KOs in drop out", "Total no of 266 KOs", "Lost 266 KOs in drop out")

new_fam_data$`Total no of KOs` <- as.numeric(new_fam_data$`Total no of KOs`)
new_fam_data$`Total sig KOs` <- as.numeric(new_fam_data$`Total sig KOs`)
new_fam_data$`Lost KOs in drop out` <- as.numeric(new_fam_data$`Lost KOs in drop out`)
new_fam_data$`Total no of 266 KOs` <- as.numeric(new_fam_data$`Total no of 266 KOs`)
new_fam_data$`Lost 266 KOs in drop out` <- as.numeric(new_fam_data$`Lost 266 KOs in drop out`)

new_fam_data_3 <- melt(new_fam_data)

for (group in unique(new_fam_data_3$Family)){
  new_fam_data_3_sub <- new_fam_data_3[new_fam_data_3$Family == paste(group),]
  prop <- as.numeric(new_fam_data_3_sub$value[new_fam_data_3_sub$variable == "Lost KOs in drop out"])/as.numeric(new_fam_data_3_sub$value[new_fam_data_3_sub$variable == "Total sig KOs"])
  have_266 <- as.numeric(new_fam_data_3_sub$value[new_fam_data_3_sub$variable == "Lost 266 KOs in drop out"])
  not_have_266 <- as.numeric(new_fam_data_3_sub$value[new_fam_data_3_sub$variable == "Total no of 266 KOs"])
  not_have_266_2 <- (not_have_266 - have_266)
  
  binom_out <- binom.test(have_266, not_have_266, prop)
  pval <- binom_out$p.value
  
  new_data <- data.frame(t(data.frame(c(paste(group), "Percentage lost KOs vs sig", paste(prop)))))
  colnames(new_data) <- colnames(new_fam_data_3)
  row.names(new_data) <- NULL
  
  new_fam_data_3 <- rbind(new_fam_data_3,new_data)
  
  new_data <- data.frame(t(data.frame(c(paste(group), "Percentage lost KOs vs sig - 266", have_266/not_have_266))))
  colnames(new_data) <- colnames(new_fam_data_3)
  row.names(new_data) <- NULL
  
  new_fam_data_3 <- rbind(new_fam_data_3,new_data)
  
  new_data <- data.frame(t(data.frame(c(paste(group), "Binomial p-value", paste(round(pval,5))))))
  colnames(new_data) <- colnames(new_fam_data_3)
  row.names(new_data) <- NULL
  
  new_fam_data_3 <- rbind(new_fam_data_3,new_data)
}

data_wide_4 <- spread(new_fam_data_3, variable, value)
#Rhizobiaceae were eventually excluded as there was contamination in the drop out
data_wide_4_2 <- data_wide_4[data_wide_4$Family != "Rhizobiaceae",]
data_wide_4_2$Dataset <- "No nodulator"
#With dominator dataset
new_fam_data <- data.frame()

for (group in families){
  if(group == "other"){
    new_tax <- row.names(taxonomy_LjSC)[!taxonomy_LjSC$family %in% families]
    sigtab_col_all_2_with_nod_sub <- sigtab_col_all_2_with_nod[sigtab_col_all_2_with_nod$Subset == "All other families drop out",]
  } else {
    new_tax <- row.names(taxonomy_LjSC)[taxonomy_LjSC$family == paste(group)]
    sigtab_col_all_2_with_nod_sub <- sigtab_col_all_2_with_nod[sigtab_col_all_2_with_nod$Subset == paste(group, "drop out", sep = " "),]
  }
  
  KO_table_sub <- KO_table[,colnames(KO_table) %in% new_tax]
  KO_table_sub_2 <- rowSums(KO_table_sub)
  total_KOs <- length(names(KO_table_sub_2)[KO_table_sub_2 != 0])
  Family_KOs <- names(KO_table_sub_2)[KO_table_sub_2 != 0]
  
  sigtab_col_all_2_with_nod_full <- sigtab_col_all_2_with_nod[sigtab_col_all_2_with_nod$Subset == "Full LjSC",]
  sigtab_col_all_2_with_nod_full_2 <- sigtab_col_all_2_with_nod_full[sigtab_col_all_2_with_nod_full$padj < 0.05,]
  sigtab_col_all_2_with_nod_full_3 <- sigtab_col_all_2_with_nod_full_2[sigtab_col_all_2_with_nod_full_2$log2FoldChange > 0,]
  contrib_full <- length(sigtab_col_all_2_with_nod_full_3$KO[sigtab_col_all_2_with_nod_full_3$KO %in% Family_KOs])
  
  sigtab_col_all_2_with_nod_2 <- sigtab_col_all_2_with_nod_sub[sigtab_col_all_2_with_nod_sub$padj < 0.05,]
  sigtab_col_all_2_with_nod_3 <- sigtab_col_all_2_with_nod_2[sigtab_col_all_2_with_nod_2$log2FoldChange > 0,]
  lost_KOs <- length(sigtab_col_all_2_with_nod_full_3$KO[!sigtab_col_all_2_with_nod_full_3$KO %in% sigtab_col_all_2_with_nod_3$KO])
  
  total_266_Fam <- length(Family_KOs[Family_KOs %in% input_table_3])
  sigtab_col_all_2_with_nod_full_3_sig <- sigtab_col_all_2_with_nod_full_3$KO[!sigtab_col_all_2_with_nod_full_3$KO %in% sigtab_col_all_2_with_nod_3$KO]
  lost_266 <- length(sigtab_col_all_2_with_nod_full_3_sig[sigtab_col_all_2_with_nod_full_3_sig %in% input_table_3])
  
  all_data <- data.frame(t(data.frame(c(paste(group), total_KOs, contrib_full, lost_KOs, total_266_Fam, lost_266 ))))
  new_fam_data <- rbind(new_fam_data, all_data)
}

row.names(new_fam_data) <- NULL
colnames(new_fam_data) <- c("Family", "Total no of KOs", "Total sig KOs", "Lost KOs in drop out", "Total no of 266 KOs", "Lost 266 KOs in drop out")

new_fam_data$`Total no of KOs` <- as.numeric(new_fam_data$`Total no of KOs`)
new_fam_data$`Total sig KOs` <- as.numeric(new_fam_data$`Total sig KOs`)
new_fam_data$`Lost KOs in drop out` <- as.numeric(new_fam_data$`Lost KOs in drop out`)
new_fam_data$`Total no of 266 KOs` <- as.numeric(new_fam_data$`Total no of 266 KOs`)
new_fam_data$`Lost 266 KOs in drop out` <- as.numeric(new_fam_data$`Lost 266 KOs in drop out`)

new_fam_data_3 <- melt(new_fam_data)

for (group in unique(new_fam_data_3$Family)){
  new_fam_data_3_sub <- new_fam_data_3[new_fam_data_3$Family == paste(group),]
  prop <- as.numeric(new_fam_data_3_sub$value[new_fam_data_3_sub$variable == "Lost KOs in drop out"])/as.numeric(new_fam_data_3_sub$value[new_fam_data_3_sub$variable == "Total sig KOs"])
  have_266 <- as.numeric(new_fam_data_3_sub$value[new_fam_data_3_sub$variable == "Lost 266 KOs in drop out"])
  not_have_266 <- as.numeric(new_fam_data_3_sub$value[new_fam_data_3_sub$variable == "Total no of 266 KOs"])
  not_have_266_2 <- (not_have_266 - have_266)
  
  binom_out <- binom.test(have_266, not_have_266, prop)
  pval <- binom_out$p.value
  
  new_data <- data.frame(t(data.frame(c(paste(group), "Percentage lost KOs vs sig", paste(prop)))))
  colnames(new_data) <- colnames(new_fam_data_3)
  row.names(new_data) <- NULL
  
  new_fam_data_3 <- rbind(new_fam_data_3,new_data)
  
  new_data <- data.frame(t(data.frame(c(paste(group), "Percentage lost KOs vs sig - 266", have_266/not_have_266))))
  colnames(new_data) <- colnames(new_fam_data_3)
  row.names(new_data) <- NULL
  
  new_fam_data_3 <- rbind(new_fam_data_3,new_data)
  
  new_data <- data.frame(t(data.frame(c(paste(group), "Binomial p-value", paste(round(pval,5))))))
  colnames(new_data) <- colnames(new_fam_data_3)
  row.names(new_data) <- NULL
  
  new_fam_data_3 <- rbind(new_fam_data_3,new_data)
}

data_wide_4_with_nod <- spread(new_fam_data_3, variable, value)
data_wide_4_2_with_nod <- data_wide_4_with_nod[data_wide_4_with_nod$Family != "Rhizobiaceae",]
data_wide_4_2_with_nod$Dataset <- "Nodulator"

data_together <- rbind(data_wide_4_2, data_wide_4_2_with_nod)

write.table(data_together, paste(results.dir, "Table_S8_KO_overview_Drop_out_exp.txt", sep = ""), sep = "\t", quote =F, col.names =T, row.names =F)


###Table S9 - KOs in Soil  =====
#KO overlap in natural rhizospheres vs natural soils

groups <- c("cucumber_rhizosphere","cucumber_soil","juanjo_arabidopsis_soil", "juanjo_arabidopsis_rhizosphere","stringlis_arabidopsis_rhizosphere","stringlis_arabidopsis_soil","wheat_ofek_rhizosphere","wheat_ofek_soil")

compile_dat <- data.frame()

for (group in groups){
  date_set <- read.table(paste(working_directory, "Natural_soil/KO_lists/KO_list_", group, ".tsv", sep =""), sep = " ", header =T)
  date_set$V3 <- paste(group)
  compile_dat <- rbind(compile_dat, date_set)
}

Soil_ofek <- unique(c(compile_dat$KEGG_ko[compile_dat$V3 == "cucumber_soil"], compile_dat$KEGG_ko[compile_dat$V3 == "wheat_ofek_soil"]))
Soil_jj_string <- unique(c(compile_dat$KEGG_ko[compile_dat$V3 == "juanjo_arabidopsis_soil"], compile_dat$KEGG_ko[compile_dat$V3 == "stringlis_arabidopsis_soil"]))
Rhizosphere_ofek <- unique(c(compile_dat$KEGG_ko[compile_dat$V3 == "cucumber_rhizosphere"], compile_dat$KEGG_ko[compile_dat$V3 == "wheat_ofek_rhizosphere"]))
Rhizosphere_jj_string <- unique(c(compile_dat$KEGG_ko[compile_dat$V3 == "juanjo_arabidopsis_rhizosphere"], compile_dat$KEGG_ko[compile_dat$V3 == "stringlis_arabidopsis_rhizosphere"]))

Soil_overlap <- table(c(Soil_ofek, Soil_jj_string))
Expectancy <- length(names(Soil_overlap)[Soil_overlap == 2])/length(names(Soil_overlap))

Root_overlap <- table(c(Rhizosphere_ofek, Rhizosphere_jj_string))
Root_overlap_2 <- length(names(Root_overlap)[Root_overlap == 2])

together <- c(Root_overlap_2, length(names(Root_overlap)) - Root_overlap_2)

binom_out_RA <- binom.test(together,length(names(Root_overlap_2)), Expectancy)

Enrichment <- binom_out_RA$p.value

dataset <- data.frame(c(length(Soil_ofek), length(Soil_jj_string)), 
                      c(length(names(Soil_overlap)[Soil_overlap == 2]),length(names(Soil_overlap)[Soil_overlap == 2])),
                      c(Expectancy, Expectancy),
                      c(length(Rhizosphere_ofek), length(Rhizosphere_jj_string)),
                      c(Root_overlap_2,Root_overlap_2),
                      c(Root_overlap_2/length(names(Root_overlap))),
                      c(Enrichment,Enrichment))
row.names(dataset) <- c("Maon soil", "Reijerscamp soil")
colnames(dataset) <- c("Soil KOs", "Soil KO overlap", "Soil KO overlap (%)", "Rhizosphere KOs", "Rhizosphere KO overlap", "Rhizosphere KO overlap (%)", "Enrichment root vs soil")

write.table(dataset, paste(results.dir,"Table_S9_Soil_Rhizosphere_KO_overlap.tsv", sep = ""), sep = "\t", quote =F, row.names =T, col.names =T)

#Table S9 - Stats of 266 and 852 KOs in natural rhizospheres
Jose_266 <- read.table(paste(working_directory, "Natural_soil/266_sig_in_natural_rhizospheres.tsv", sep = ""), header =T)
Jose_852 <- read.table(paste(working_directory, "Natural_soil/852_sig_in_natural_rhizospheres.tsv", sep = ""), header =T)

#266 KOs - overlap
input_table <- read.table(paste(working_directory,"DESeq2/Sig_KO_all_no_nod_rhizo.txt", sep = ""), header=T, sep="\t")
input_table_2 <- table(input_table$KO)
input_table_3 <- names(input_table_2)[input_table_2 == 12]

Jose_both <- rbind(Jose_266, Jose_852)
Jose_both_2 <- unique(Jose_both)
Jose_266_2 <- Jose_both_2[Jose_both_2$Orthogroup_Id %in% input_table_3,]

JJ_overlap_sig_266 <- length(na.omit(Jose_266_2$Orthogroup_Id[Jose_266_2$Sanchez == "YES"]))
Jose_overlap_sig_266 <- length(na.omit(Jose_266_2$Orthogroup_Id[Jose_266_2$Lopez == "YES"]))

#852 KO overlap
JJ_overlap_sig_852 <- length(na.omit(Jose_852$Orthogroup_Id[Jose_852$Sanchez == "YES"]))
Jose_overlap_sig_852 <- length(na.omit(Jose_852$Orthogroup_Id[Jose_852$Lopez == "YES"]))

#Enrichment tests
stat_test <- function(vector1, vector2){
  # Calculate the sizes of the two vectors
  size_vector1 <- length(vector1)
  size_vector2 <- length(vector2)
  # Calculate the overlap
  overlap <- length(intersect(vector1, vector2))  # Number of shared identifiers
  # Hypergeometric test
  p_value <- phyper(
    q = overlap - 1,            # Overlap - 1 for "probability of more extreme overlaps"
    m = size_vector1,           # Size of the first group
    n = universe_size - size_vector1,  # Size of the remaining universe
    k = size_vector2,           # Size of the second group
    lower.tail = FALSE          # Use the upper tail for significance of overlap
  )
  #p_value
  fisher_test <- fisher.test(matrix(c(overlap,
                                      size_vector1 - overlap,
                                      size_vector2 - overlap,
                                      universe_size - size_vector1 - size_vector2 + overlap),
                                    nrow = 2))
  #fisher_test$p.value
  
  # Output
  #cat("Overlap:", overlap, "\n")
  #cat("P-value:", p_value, "\n")
  return(data.frame(HG_pvalue = p_value, Fisher_pvalue = fisher_test$p.value, 
                    length_overlap = overlap))
}

#Load KO lists
ST <- read.table(paste(working_directory,"KO_genome/KO_SSC.tsv", sep = ""),header = TRUE)
stringent_KO <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/266KO.tsv", sep =""))
general_KO <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/852KO.tsv", sep =""))

### Compare to Jose KOs - 266
Jose_KO <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/Lopez_et_al_2023_output.tsv", sep= ""))
Jose_KO_R <- filter(Jose_KO, niche_association == "Rhizosphere")
Jose_KO_R_2 <- unique(dplyr::select(Jose_KO_R, Orthogroup_Id, molecular_function))
Jose_KO_R_2$Lopez <- "YES"

vector1_Jose_266 <- unique(Jose_KO_R_2$Orthogroup_Id)
vector2 <- unique(stringent_KO$KO)
KO_universe <- union(Jose_KO$Orthogroup_Id, ST$sequence)
universe_size <- length(KO_universe)
universe_size_Jose_266 <- length(KO_universe)
Jose_266_stat <- stat_test(vector1_Jose_266, vector2)

### Compare to Jose KOs - 852
vector1_Jose_852 <- unique(Jose_KO_R_2$Orthogroup_Id)
vector2 <- unique(general_KO$KO)
KO_universe <- union(Jose_KO$Orthogroup_Id, ST$sequence)
universe_size <- length(KO_universe)
universe_size_Jose_852 <- length(KO_universe)
Jose_852_stat <- stat_test(vector1_Jose_852, vector2)

### Compare to JJ KOs - 266
res_ara <- readRDS(paste(working_directory,"Natural_soil/ANCOMBC_output/ANCOMBC2_ara.rds", sep = ""))
res_ara <- res_ara%>%dplyr::select(taxon, lfc_TreatmentSoil, p_TreatmentSoil)
colnames(res_ara) <- c("KO", "LFC", "padj")
res_ara$LFC <- -1*res_ara$LFC
res_ara_JJ <- readRDS(paste(working_directory,"Natural_soil/ANCOMBC_output/ANCOMBC2_ara_JJ.rds", sep = ""))
res_ara_JJ <- res_ara_JJ%>%dplyr::select(taxon, lfc_TreatmentSoil, p_TreatmentSoil)
colnames(res_ara_JJ) <- c("KO", "LFC", "padj")
res_ara_JJ$LFC <- -1*res_ara_JJ$LFC

res_ara$padj <- p.adjust(res_ara$padj, method = "fdr")
res_ara_JJ$padj <- p.adjust(res_ara_JJ$padj, method = "fdr")

ancombc2_araC <- res_ara
ancombc2_araJ <- res_ara_JJ

#Since the ANCOMBC analyses on ancombc2_araC led to no significant results, it was later discarded, only ancombc2_araJ was taken
ancombc2_araC_R <- filter(ancombc2_araC, padj<0.05&LFC>0)
ancombc2_araJ_R <- filter(ancombc2_araJ, padj<0.05&LFC>0)

vector1_JJ_266 <- unique(ancombc2_araJ_R$KO)
vector2 <- unique(stringent_KO$KO)
KO_universe <- union(ancombc2_araJ_R$KO, ST$sequence)
universe_size <- length(KO_universe)
universe_size_JJ_266 <- length(KO_universe)
JJ_266_stat <- stat_test(vector1_JJ_266, vector2)

### Compare to JJ KOs - 852
vector1_JJ_852 <- unique(ancombc2_araJ_R$KO)
vector2 <- unique(general_KO$KO)
KO_universe <- union(ancombc2_araJ_R$KO, ST$sequence)
universe <- length(KO_universe)
universe_size_JJ_852 <- length(KO_universe)
JJ_852_stat <- stat_test(vector1_JJ_852, vector2)

dataset <- data.frame(c("Lopez et al. (2023)", "Lopez et al. (2023)", "Sanchez et al. unpublished", "Sanchez et al. unpublished"),
                      c("266", "852", "266", "852"),
                      c(length(vector1_Jose_266), length(vector1_Jose_852), length(vector1_JJ_266), length(vector1_JJ_852)), 
                      c(universe_size_Jose_266,universe_size_Jose_852, universe_size_JJ_266,universe_size_JJ_852),
                      c(Jose_overlap_sig_266,Jose_overlap_sig_852,JJ_overlap_sig_266, JJ_overlap_sig_852),
                      c(Jose_266_stat$HG_pvalue,Jose_852_stat$HG_pvalue,JJ_266_stat$HG_pvalue,JJ_852_stat$HG_pvalue),
                      c(Jose_266_stat$Fisher_pvalue,Jose_852_stat$Fisher_pvalue,JJ_266_stat$Fisher_pvalue,JJ_852_stat$Fisher_pvalue))

colnames(dataset) <- c("Data", "KO set", "Number of Rhizosphere KOs", "Union of KO set and Rhizosphere KOs","Overlap KOs", "Hypergeometric test p-value", "Fisher test p-value")

write.table(dataset, paste(results.dir,"Table_S9_266_852_KOs_in_Natural_soil.tsv", sep = ""), sep = "\t", quote =F, row.names =T, col.names =T)

############################# Scripts to generate files ###########################
###Script to generate pangenome_order.tsv for SSC - Necessary for 1d =====
SynComs <- c("AtSC", "HvSC", "LjSC", "SSC")
table_4 <- data.frame()

for (inoculum in SynComs) {
  table <- read.table(paste(working_directory, "KO_genome/KO_",inoculum,".tsv", sep = ""), sep= "\t", header =T) 
  row.names(table) <- table$sequence
  table_2 <- table %>% dplyr::select (-sequence)
  
  table_2[is.na(table_2)] <- 0
  
  vector <- colnames(table_2)
  new_list <- list()
  
  new_order <- c()
  
  vector_3 <- vector
  i = 1
  removing_KOs <- c()
  
  while(i <=length(vector)){
    vector_2 <- vector_3
    new_list <- list()
    for (isolate in vector_2) {
      sub <- table_2[colnames(table_2) == paste(isolate)]
      sub_2 <- row.names(sub)[sub>0]
      new_list[[paste(isolate)]] <- sub_2
    }
    
    new <- as.data.frame(names(new_list))
    colnames(new) <- "isolate"
    
    no_table <- as.data.frame(unique(unlist(new_list)))
    colnames(no_table) <- "KO_term"
    no_table_2 <- no_table[!no_table$KO_term %in% removing_KOs,]
    new_list_2 <- unlist(new_list)
    new_list_3 <- new_list_2[!new_list_2 %in% removing_KOs]
    
    new_list_4 <- as.data.frame(table(new_list_3))
    
    isolates <- as.vector(new)
    
    for (isolate in isolates$isolate) {
      hop <- new_list[paste(isolate)]
      new_vector <- as.vector(unlist(hop))
      new_vector_2 <- new_vector[!new_vector %in% removing_KOs]
      other_new <- as.data.frame(new_vector_2)
      colnames(other_new) <- "KO_term"
      for (another in new_vector_2){
        value <- new_list_4$Freq[new_list_4$new_list_3 == paste(another)]
        other_new$number[other_new$KO_term == paste(another)] <- value
      }
      new_value <- length(other_new$KO_term[other_new$number == 1])
      new$number[new$isolate == paste(isolate)] <- new_value
    }
    new_2 <- new[order(new$number, decreasing = T),]
    new_order <- c(new_order,new_2$isolate[1])
    
    isolate_3 <- new_2$isolate[1]
    sub_3 <- table_2[colnames(table_2) == paste(isolate_3)]
    sub_4 <- row.names(sub_3)[sub_3>0]
    removing_KOs <- c(removing_KOs, sub_4)
    removing_KOs <- unique(removing_KOs)
    vector_3 <- vector_2[!vector_2 == new_2$isolate[1]]
    i <- i +1
    
  } 
  
  SynCom_vec <- vector()
  other_vec <- vector()
  
  for (column in new_order) {
    table_sub <- table_2[,colnames(table_2) == paste(column)]
    table_sub_2 <- as.data.frame(table_sub)
    row.names(table_sub_2) <- row.names(table_2)
    table_sub_3 <- row.names(table_sub_2)[table_sub_2[1] > 0]
    SynCom_vec <- c(SynCom_vec,table_sub_3)
    new_value <- length(unique(SynCom_vec))
    other_vec <- c(other_vec,new_value)
  }
  
  table_3 <- t(as.data.frame(other_vec))
  colnames(table_3) <- 1:length(other_vec)
  row.names(table_3) <- paste(inoculum)
  
  late_vector <- rep(NA,1000 -length(table_3))
  table_3 <- c(table_3,late_vector)
  table_4 <- rbind(table_4,table_3)
}

row.names(table_4) <- SynComs
colnames(table_4) <- 1:1000
table_4$SynCom <- row.names(table_4)
table_5 <- melt(table_4)
table_6 <- na.omit(table_5)

write.table(table_6, paste(working_directory, "KO_intravariability/pangenome_order.tsv", sep=""), sep = "\t", quote = F, row.names = F, col.names = T)

###Script to generate pangenome_order.tsv for simulations - Necessary for 1e ====
table_4 <- data.frame()

table_all <- read.table(paste(working_directory, "KO_intravariability/simulations/random_family_sets.csv", sep = ""), header =T, row.names =1)
SSC_data <- read.table(paste(working_directory, "KO_genome/KO_SSC.tsv",sep = ""), sep= "\t", header =T, row.names =1) 

colnames(SSC_data)[grep("M.16",colnames(SSC_data))] <- "M-16"
colnames(SSC_data)[grep("M.6",colnames(SSC_data))] <- "M-6"
colnames(SSC_data)[grep("M.10",colnames(SSC_data))] <- "M-10"
colnames(SSC_data)[grep("M.11_2",colnames(SSC_data))] <- "M-11_2"
colnames(SSC_data)[grep("^M.11$",colnames(SSC_data))] <- "M-11"
colnames(SSC_data) <- gsub("X", "", colnames(SSC_data))

for (simulation in 1:1000) {
  table_all_sub <- unlist(as.vector(table_all[simulation,]))
  
  SSC_data_2 <- SSC_data[,colnames(SSC_data) %in% table_all_sub]
  
  vector <- colnames(SSC_data_2)
  new_list <- list()
  
  new_order <- c()
  
  vector_3 <- vector
  i = 1
  removing_KOs <- c()
  
  while(i <=length(vector)){
    vector_2 <- vector_3
    new_list <- list()
    for (isolate in vector_2) {
      sub <- SSC_data_2[colnames(SSC_data_2) == paste(isolate)]
      sub_2 <- row.names(sub)[sub>0]
      new_list[[paste(isolate)]] <- sub_2
    }
    
    new <- as.data.frame(names(new_list))
    colnames(new) <- "isolate"
    
    no_table <- as.data.frame(unique(unlist(new_list)))
    colnames(no_table) <- "KO_term"
    no_table_2 <- no_table[!no_table$KO_term %in% removing_KOs,]
    new_list_2 <- unlist(new_list)
    new_list_3 <- new_list_2[!new_list_2 %in% removing_KOs]
    
    new_list_4 <- as.data.frame(table(new_list_3))
    
    isolates <- as.vector(new)
    
    for (isolate in isolates$isolate) {
      hop <- new_list[paste(isolate)]
      new_vector <- as.vector(unlist(hop))
      new_vector_2 <- new_vector[!new_vector %in% removing_KOs]
      other_new <- as.data.frame(new_vector_2)
      colnames(other_new) <- "KO_term"
      for (another in new_vector_2){
        value <- new_list_4$Freq[new_list_4$new_list_3 == paste(another)]
        other_new$number[other_new$KO_term == paste(another)] <- value
      }
      new_value <- length(other_new$KO_term[other_new$number == 1])
      new$number[new$isolate == paste(isolate)] <- new_value
    }
    new_2 <- new[order(new$number, decreasing = T),]
    new_order <- c(new_order,new_2$isolate[1])
    
    isolate_3 <- new_2$isolate[1]
    sub_3 <- SSC_data_2[colnames(SSC_data_2) == paste(isolate_3)]
    sub_4 <- row.names(sub_3)[sub_3>0]
    removing_KOs <- c(removing_KOs, sub_4)
    removing_KOs <- unique(removing_KOs)
    vector_3 <- vector_2[!vector_2 == new_2$isolate[1]]
    i <- i +1
    
  } 
  
  SynCom_vec <- vector()
  other_vec <- vector()
  
  for (column in new_order) {
    table_sub <- SSC_data_2[,colnames(SSC_data_2) == paste(column)]
    table_sub_2 <- as.data.frame(table_sub)
    row.names(table_sub_2) <- row.names(SSC_data_2)
    table_sub_3 <- row.names(table_sub_2)[table_sub_2[1] > 0]
    SynCom_vec <- c(SynCom_vec,table_sub_3)
    new_value <- length(unique(SynCom_vec))
    other_vec <- c(other_vec,new_value)
  }
  
  table_3 <- t(as.data.frame(other_vec))
  colnames(table_3) <- 1:length(other_vec)
  row.names(table_3) <- paste(simulation)
  
  table_4 <- rbind(table_4,table_3)
}

row.names(table_4) <- NULL
colnames(table_4) <- 1:25
table_4$SynCom <- row.names(table_4)
table_5 <- melt(table_4)
table_6 <- na.omit(table_5)

write.table(table_6, paste(working_directory,"KO_intravariability/simulations/pangenome_order.tsv", sep=""), sep = "\t", quote = F, row.names = F, col.names = T)


###Script to generate KO_intravariability_2.tsv for SSC KO intravariability - Necessary for S1a =====
KO_to_gene = read.table(paste(working_directory,"KO_to_gene.txt", sep = ""), header=TRUE,sep="\t")
KO_to_gene$SynCom <- NA
list_AtSC =read.table(paste(working_directory,"Isolate_tables/Original/AtSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
list_AtSC <- row.names(list_AtSC)
list_LjSC =read.table(paste(working_directory,"Isolate_tables/Original/LjSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
list_LjSC <- row.names(list_LjSC)
list_HvSC =read.table(paste(working_directory,"Isolate_tables/Original/HvSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
list_HvSC <- row.names(list_HvSC)

KO_to_gene$SynCom[KO_to_gene$isolate %in% list_AtSC] <- "AtSC"
KO_to_gene$SynCom[KO_to_gene$isolate %in% list_LjSC] <- "LjSC"
KO_to_gene$SynCom[KO_to_gene$isolate %in% list_HvSC] <- "HvSC"
KO_list <- unique(KO_to_gene$kegg)

new_KO_list_3 <- data.frame(matrix(NA, ncol = 9))
colnames(new_KO_list_3) <- c("KO","AtSC_total_branch_length","LjSC_total_branch_length","HvSC_total_branch_length","SSC_total_branch_length", "No_of_genes_AtSC", "No_of_genes_LjSC","No_of_genes_HvSC", "No_of_genes_SSC")
new_KO_list_4 <- new_KO_list_3[-1,]

for (KO_term in KO_list){
  print(KO_term)
  #Phylogenetic tree
  subset <- KO_to_gene[ which(KO_to_gene$kegg == paste0(KO_term)),]
  subset_2 <- subset[!duplicated(subset), ]
  subset_2$name <- paste(subset_2$isolate, subset_2$gene, sep = "_")
  
  subset_AtSC <- subset_2[subset_2$SynCom == "AtSC",]
  subset_LjSC <- subset_2[subset_2$SynCom == "LjSC",]
  subset_HvSC <- subset_2[subset_2$SynCom == "HvSC",]
  
  subset_3 <- subset_AtSC[,4:6]
  subset_3 <- subset_3 %>% dplyr::select (-SynCom)
  
  if (length(na.omit(subset_3$sequence)) > 1){
    subset_4 <- subset_3[!(is.na(subset_3$sequence) | subset_3$sequence==""), ]
    
    y <- strsplit(subset_4[,1],"")
    names(y) <- subset_4[,2]
    
    df.fasta <- ape::as.DNAbin(y)
    dist_matrix <- kmer::kdistance(df.fasta)
    
    # Produce dendrogram
    hclust = hclust(dist_matrix)
    
    my_tree <- as.phylo(hclust) 
    len_AtSC <- sum(my_tree$edge.length)
    tree_AtSC <- length(names(y))
    
  } else if (length(na.omit(subset_3$sequence)) > 0 ) { 
    len_AtSC <- 1
    tree_AtSC <- 1
  } else {
    len_AtSC <- 0
    tree_AtSC <- 0
  }
  
  subset_3 <- subset_LjSC[,4:6]
  subset_3 <- subset_3 %>% dplyr::select (-SynCom)
  
  if (length(na.omit(subset_3$sequence)) > 1){
    subset_4 <- subset_3[!(is.na(subset_3$sequence) | subset_3$sequence==""), ]
    
    y <- strsplit(subset_4[,1],"")
    names(y) <- subset_4[,2]
    
    df.fasta <- ape::as.DNAbin(y)
    dist_matrix <- kmer::kdistance(df.fasta)
    
    # Produce dendrogram
    hclust = hclust(dist_matrix)
    
    my_tree <- as.phylo(hclust) 
    len_LjSC <- sum(my_tree$edge.length)
    tree_LjSC <- length(names(y))
    
  } else if (length(na.omit(subset_3$sequence)) > 0) { 
    len_LjSC <- 1
    tree_LjSC <- 1
  } else {
    len_LjSC <- 0
    tree_LjSC <- 0
  }
  
  subset_3 <- subset_HvSC[,4:6]
  subset_3 <- subset_3 %>% dplyr::select (-SynCom)
  
  if (length(na.omit(subset_3$sequence)) > 1){
    subset_4 <- subset_3[!(is.na(subset_3$sequence) | subset_3$sequence==""), ]
    
    y <- strsplit(subset_4[,1],"")
    names(y) <- subset_4[,2]
    
    df.fasta <- ape::as.DNAbin(y)
    dist_matrix <- kmer::kdistance(df.fasta)
    
    # Produce dendrogram
    hclust = hclust(dist_matrix)
    
    my_tree <- as.phylo(hclust) 
    len_HvSC <- sum(my_tree$edge.length)
    tree_HvSC <- length(names(y))
    
  } else if (length(na.omit(subset_3$sequence)) > 0 ){ 
    len_HvSC <- 1
    tree_HvSC <- 1
  } else {
    len_HvSC <- 0
    tree_HvSC <- 0
  }
  
  subset_3 <- subset_2[,4:6]
  subset_3 <- subset_3 %>% dplyr::select (-SynCom)
  
  if (length(na.omit(subset_3$sequence)) > 2){
    subset_4 <- subset_3[!(is.na(subset_3$sequence) | subset_3$sequence==""), ]
    
    y <- strsplit(subset_4[,1],"")
    names(y) <- subset_4[,2]
    
    df.fasta <- ape::as.DNAbin(y)
    dist_matrix <- kmer::kdistance(df.fasta)
    
    # Produce dendrogram
    hclust = hclust(dist_matrix)
    
    my_tree <- as.phylo(hclust) 
    len_SSC <- sum(my_tree$edge.length)
    tree_SSC <- length(names(y))
  } else if (length(na.omit(subset_3$sequence)) > 0 ) { 
    len_SSC <- 1
    tree_SSC <- 1
  } else {
    len_SSC <- 0
    tree_SSC <- 0
  }

  vector_2 <- c(paste(KO_term, sep = ""), len_AtSC,len_LjSC,len_HvSC,len_SSC,tree_AtSC, tree_LjSC, tree_HvSC, tree_SSC)
  new_KO_list_4 <- rbind(new_KO_list_4, vector_2)
}

colnames(new_KO_list_4) <- c("KO","AtSC_total_branch_length","LjSC_total_branch_length","HvSC_total_branch_length","SSC_total_branch_length", "No_of_genes_AtSC", "No_of_genes_LjSC","No_of_genes_HvSC", "No_of_genes_SSC")
write.table(new_KO_list_4, paste(working_directory,"KO_intravariability/KO_intrafunctionality_2.tsv", sep = ""), sep = "\t", quote = F)

###Script to generate KO_intravariability_2_simulations.tsv - Necessary for S1B ====!!!!!!!!!!=
###Script to generate KO_intravarability.tsv for 1000 one-strain-per-family SynComs - Necessary for S1b ======
KO_to_gene = read.table(paste(working_directory,"KO_to_gene.txt", sep = ""), header=TRUE,sep="\t")

table_all <- read.table(paste(working_directory, "KO_intravariability/simulations/random_family_sets.csv", sep = ""), header =T, row.names =1)
SSC_data <-   read.table(paste(working_directory, "KO_genome/KO_SSC.tsv",sep = ""), sep= "\t", header =T, row.names =1) 
colnames(SSC_data)[grep("M.16",colnames(SSC_data))] <- "M-16"
colnames(SSC_data)[grep("M.6",colnames(SSC_data))] <- "M-6"
colnames(SSC_data)[grep("M.10",colnames(SSC_data))] <- "M-10"
colnames(SSC_data)[grep("M.11_2",colnames(SSC_data))] <- "M-11_2"
colnames(SSC_data)[grep("^M.11$",colnames(SSC_data))] <- "M-11"
colnames(SSC_data) <- gsub("X", "", colnames(SSC_data))

new_KO_list_2 <- data.frame()

for (syncom in 1:1000){
  table_all_sub <- as.vector(unlist(as.vector(table_all[syncom,])))
  SSC_data_sub <- SSC_data[,colnames(SSC_data) %in% table_all_sub]
  
  KO_to_gene_sub <- KO_to_gene[KO_to_gene$isolate %in% table_all_sub,]
  
  KO_list <- row.names(SSC_data_sub)[rowSums(SSC_data_sub) != 0]
  
  for (KO_term in KO_list){
    print(syncom)
    print(KO_term)
    #Phylogenetic tree
    subset <- KO_to_gene_sub[ which(KO_to_gene_sub$kegg == paste0(KO_term)),]
    subset_2 <- subset[!duplicated(subset), ]
    subset_2$name <- paste(subset_2$isolate, subset_2$gene, sep = "_")
    subset_3 <- subset_2[,4:5]
    
    if (length(na.omit(subset_3$sequence)) > 1){
      subset_4 <- subset_3[!(is.na(subset_3$sequence) | subset_3$sequence==""), ]
      
      y <- strsplit(subset_4[,1],"")
      names(y) <- subset_4[,2]
      
      df.fasta <- ape::as.DNAbin(y)
      dist_matrix <- kmer::kdistance(df.fasta)
      
      # Produce dendrogram
      hclust = hclust(dist_matrix)
      
      my_tree <- as.phylo(hclust) 
      
      value_tree <- sum(my_tree$edge.length)/length(names(y))
      
    } else if (length(na.omit(subset_3$sequence)) > 0 ) { 
      value_tree <- 1
    } else {
      value_tree <- 0
    }
    
    vector <- c(paste(syncom),paste(KO_term), value_tree)
    new_KO_list_2 <- rbind(new_KO_list_2, vector)
  }
}

colnames(new_KO_list_2) <- c("SynCom", "KO" ,"value")
KO_out <- new_KO_list_2

new_KO_lists <- data.frame()

for (syncom in 1:1000){
  table_all_sub <- as.vector(unlist(as.vector(table_all[syncom,])))
  SSC_data_sub <- SSC_data[,colnames(SSC_data) %in% table_all_sub]
  
  KO_to_gene_sub <- KO_to_gene[KO_to_gene$isolate %in% table_all_sub,]
  
  KO_list <- row.names(SSC_data_sub)[rowSums(SSC_data_sub) != 0]
  
  KO_out_out <- KO_out[KO_out$SynCom == paste(syncom),]
  
  for (KO_term in KO_list){
    print(syncom)
    print(KO_term)
    subset <- KO_to_gene_sub[ which(KO_to_gene_sub$kegg == paste0(KO_term)),]
    subset_2 <- subset[!duplicated(subset), ]
    length_value <- length(subset_2$kegg)
    
    KO_out_out_2 <- KO_out_out[KO_out_out$KO == paste(KO_term),]
    
    total_branch_length_value <- as.numeric(KO_out_out_2$value) * length_value
    
    vector <- c(paste(syncom),paste(KO_term), total_branch_length_value, length_value)
    new_KO_lists <- rbind(new_KO_lists, vector)
  }
}

colnames(new_KO_lists) <- c("SynCom", "KO" ,"total_branch_length","no_of_genes")
write.table(new_KO_lists, paste(working_directory, "KO_intravariability/simulations/KO_intravariability.txt", sep = ""),col.names = T, row.names =F, quote =F, sep = "\t")

###Script to generate R2_values_genus.txt - Necessary for S12 =====
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")

R2_values <- data.frame(matrix(NA, ncol = 3))
colnames(R2_values) <- c("run", "Tax_R2" ,"Func_R2")
R2_values_2 <- R2_values[-1,]

for (i in 1:1000){
  tax_table <- read.table(paste(working_directory,"simulation_R/run_",i,"/table_",i,".txt", sep = ""), sep = "\t", header = F)
  row.names(tax_table) <- tax_table$V1
  tax_table_2 <- tax_table %>% dplyr::select (-V1)
  vector <- c()
  for (j in 1:54){
    vector <- c(vector, paste("sample_", j, sep = ""))
  }
  colnames(tax_table_2) <- vector
  
  #Set the OTU, TAX and sample data for making phyloseq object
  OTU = otu_table(as.matrix(tax_table_2),taxa_are_rows = TRUE)
  #TAX = tax_table(tax_mat)
  TAX = tax_table(as.matrix(tax_df_2))
  
  phylo = phyloseq(OTU,TAX)
  
  phylo_RA=microbiome::transform(x = phylo, transform = "compositional" )
  
  #Agglomerate to phylum-level and rename
  phylo_RA <- phyloseq::tax_glom(phylo_RA, "Genus")
  phyloseq::taxa_names(phylo_RA) <- phyloseq::tax_table(phylo_RA)[, "Genus"]
  
  #Bray Curtis distance matrix
  beta_tax <- as.matrix(vegdist(t(phylo_RA@otu_table@.Data), method = "bray", diag = T))
  
  #Make PCoA plot for Bray Curtis Distance matrix
  pcoa = cmdscale(beta_tax, k=3, eig=T)
  points = as.data.frame(pcoa$points)
  colnames(points) = c("x", "y", "z") 
  eig = pcoa$eig
  
  set.seed(1)
  clusters <- kmeans(points, 9, iter.max = 10, nstart = 1)
  clusters_2 <- as.data.frame(clusters$cluster)
  colnames(clusters_2) <- "Cluster"
  
  set.seed(1)
  Tax_adonis <- adonis2(beta_tax ~ Cluster, data=clusters_2, method="bray", permutations=999)
  Tax_value <- Tax_adonis$R2[1]
  
  KO_table <- read.table(paste(working_directory,"simulation_R/run_",i,"/KO_table_",i,".txt", sep = ""), sep = "\t", header = T)
  row.names(KO_table) <- KO_table$function.
  KO_table_2 <- KO_table %>% dplyr::select (-function.)
  
  #Set the OTU, TAX and sample data for making phyloseq object
  OTU_KO = otu_table(as.matrix(KO_table_2),taxa_are_rows = TRUE)
  
  phylo_KO = phyloseq(OTU_KO)
  
  phylo_KO_RA=microbiome::transform(x = phylo_KO, transform = "compositional" )
  
  
  #Bray Curtis distance matrix
  beta_KO <- as.matrix(vegdist(t(phylo_KO_RA), method = "bray", diag = T))
  set.seed(1)
  KO_adonis <- adonis2(beta_KO ~ Cluster, data=clusters_2, method="bray", permutations=999)
  KO_value <- KO_adonis$R2[1]
  
  new <- c(paste("run_", i, sep = ""), Tax_value, KO_value)
  
  R2_values_2 <- rbind(R2_values_2, new)
}

colnames(R2_values_2) <- c("run", "Tax_R2" ,"Func_R2")
write.table(R2_values_2, paste(working_directory, "R2_values_genus.txt", sep= ""), quote = F, sep = "\t", col.names = T, row.names =T)

###Script to generate ternary_KOs_av_med.txt with dom and without dom - Necessary for 5b, c, and e - lenient and strict thresholds =====

KO_table = read.table(paste(working_directory, "KO_genome/KO_SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
colnames(KO_table) <- gsub("X", "", colnames(KO_table))

tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)

rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)

input_table <- read.table(paste(working_directory, "DESeq2/Sig_KO_all.txt", sep = ""), header=T, sep="\t")
input_table_2 <- table(input_table$KO)
input_table_3 <- names(input_table_2)[input_table_2 > 6]

input_table_At <- input_table[input_table$Plant == "Arabidopsis",]
input_table_At_2 <- table(input_table_At$KO)
input_table_At_3 <- names(input_table_At_2)[input_table_At_2 > 2]

input_table_Hv <- input_table[input_table$Plant == "Barley",]
input_table_Hv_2 <- table(input_table_Hv$KO)
input_table_Hv_3 <- names(input_table_Hv_2)[input_table_Hv_2 > 2]

input_table_Lj <- input_table[input_table$Plant == "Lotus",]
input_table_Lj_2 <- table(input_table_Lj$KO)
input_table_Lj_3 <- names(input_table_Lj_2)[input_table_Lj_2 > 2]

input_table_3 <- unique(c(input_table_3,input_table_At_3, input_table_Hv_3, input_table_Lj_3))

SynComs <- c("AtSC","LjSC", "HvSC", "SSC")
hop_2 <- data.frame()

for (KO in input_table_3){
  for (syncom in SynComs){
    syncom_table <- read.table(paste(working_directory,"Isolate_tables/Original/", syncom,"_norm.tsv", sep= ""), sep = "\t", header =T, row.names =1)
    syncom_table_2 <- syncom_table[,grep("ES", colnames(syncom_table))]
    syncom_table_3 <- syncom_table_2[,grep(paste(syncom), colnames(syncom_table_2))]
    syncom_table_4 <- syncom_table_3[,!grepl("HL", colnames(syncom_table_3))]
    
    KO_table_sub <- KO_table[row.names(KO_table) == paste(KO),]
    
    if (syncom == "SSC"){
      KO_table_sub_2 <- KO_table_sub
    } else {
      KO_table_sub_2 <- KO_table_sub[,colnames(KO_table_sub) %in% row.names(tax_df_2)[tax_df_2$SynCom == paste(syncom)]]
    }
    
    KO_table_sub_yes <- names(KO_table_sub_2)[KO_table_sub_2 > 0]
    KO_table_sub_no <- names(KO_table_sub_2)[KO_table_sub_2 == 0]
    
    syncom_table_5 <- t(t(syncom_table_4)/rowSums(t(syncom_table_4)))
    
    syncom_table_At <- syncom_table_5[,grep("At_", colnames(syncom_table_5))]
    syncom_table_Hv <- syncom_table_5[,grep("Hv_", colnames(syncom_table_5))]
    syncom_table_Lj <- syncom_table_5[,grep("Lj_", colnames(syncom_table_5))]
    
    #Averages
    syncom_table_At_2 <- rowSums(syncom_table_At)/length(colnames(syncom_table_At))
    syncom_table_Hv_2 <- rowSums(syncom_table_Hv)/length(colnames(syncom_table_Hv))
    syncom_table_Lj_2 <- rowSums(syncom_table_Lj)/length(colnames(syncom_table_Lj))
    
    At_RA <- sum(syncom_table_At_2[names(syncom_table_At_2) %in% KO_table_sub_yes])
    Hv_RA <- sum(syncom_table_Hv_2[names(syncom_table_Hv_2) %in% KO_table_sub_yes])
    Lj_RA <- sum(syncom_table_Lj_2[names(syncom_table_Lj_2) %in% KO_table_sub_yes])
    
    syncom_table_inp <- syncom_table[,grep("Input", colnames(syncom_table))]
    syncom_table_inp_2 <- t(t(syncom_table_inp)/rowSums(t(syncom_table_inp)))
    
    syncom_table_inp_3  <- rowSums(syncom_table_inp_2)/length(colnames(syncom_table_inp_2))
    
    Input_RA <- sum(syncom_table_inp_3[names(syncom_table_inp_3) %in% KO_table_sub_yes])
    
    if(Input_RA == 0){
      At_val <- At_RA
      Hv_val <- Hv_RA
      Lj_val <- Lj_RA
    } else {
      At_val <- At_RA/Input_RA
      Hv_val <- Hv_RA/Input_RA
      Lj_val <- Lj_RA/Input_RA
    }
    
    len_val <- length(KO_table_sub_yes)/length(KO_table_sub_2)
    
    hop <- t(data.frame(c(paste(KO), At_val, Hv_val, Lj_val, len_val, paste(syncom), length(colnames(syncom_table_At)), length(colnames(syncom_table_Hv)),length(colnames(syncom_table_Lj)))))
    
    hop_2 <- rbind(hop_2, hop)
  }
}

row.names(hop_2) <- NULL
colnames(hop_2) <- c("KO", "At_val", "Hv_val", "Lj_val", "No_of_strains", "SynCom", "No_of_samples_At", "No_of_samples_Hv", "No_of_samples_Lj")

hop_4 <- data.frame()

for (KO in input_table_3){
  hop_sub <- hop_2[hop_2$KO == paste(KO),]
  
  new_3 <- data.frame()
  for (syncom in unique(hop_sub$SynCom)){
    hop_sub_2 <- hop_sub[hop_sub$SynCom == paste(syncom),]
    At_val <- as.numeric(hop_sub_2$At_val)
    Hv_val <- as.numeric(hop_sub_2$Hv_val)
    Lj_val <- as.numeric(hop_sub_2$Lj_val)
    new_2 <- data.frame(At_val,Hv_val, Lj_val,hop_sub_2$No_of_samples_At,hop_sub_2$No_of_samples_Hv, hop_sub_2$No_of_samples_Lj, paste(syncom))
    new_3 <- rbind(new_3,new_2)
  }
  
  syncom_table_inp_3 = apply(syncom_table_inp_2, 1, median, na.rm=TRUE)
  
  #Medians
  At_val <- median(as.numeric(new_3$At_val))
  Hv_val <- median(as.numeric(new_3$Hv_val))
  Lj_val <- median(as.numeric(new_3$Lj_val))
  
  No_of_strains <- sum(as.numeric(hop_sub$No_of_strains))/length(hop_sub$No_of_strains)
  
  hop_3 <- t(data.frame(c(paste(KO), At_val,Hv_val, Lj_val, No_of_strains)))
  
  hop_4 <- rbind(hop_4, hop_3)
}

row.names(hop_4) <- NULL
colnames(hop_4) <- c("KO", "Arabidopsis", "Barley", "Lotus", "Proportion_of_strains")

write.table(hop_4, paste(working_directory,"Functionality/852/ternary_KOs_av_med.txt", sep = ""), quote =F, col.names =T, row.names =T, sep ="\t")

#Generating file for ternary plots
KO_table = read.table(paste(working_directory, "KO_genome/KO_SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
colnames(KO_table) <- gsub("X", "", colnames(KO_table))

tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)

input_table <- read.table(paste(working_directory, "DESeq2/Sig_KO_all_no_nod_rhizo.txt", sep = ""), header=T, sep="\t")
input_table_2 <- table(input_table$KO)
input_table_3 <- names(input_table_2)[input_table_2 == 12]

input_table_At <- input_table[input_table$Plant == "Arabidopsis",]
input_table_At_2 <- table(input_table_At$KO)
input_table_At_3 <- names(input_table_At_2)[input_table_At_2 > 2]

input_table_Hv <- input_table[input_table$Plant == "Barley",]
input_table_Hv_2 <- table(input_table_Hv$KO)
input_table_Hv_3 <- names(input_table_Hv_2)[input_table_Hv_2 > 2]

input_table_Lj <- input_table[input_table$Plant == "Lotus",]
input_table_Lj_2 <- table(input_table_Lj$KO)
input_table_Lj_3 <- names(input_table_Lj_2)[input_table_Lj_2 > 2]

input_table_3 <- unique(c(input_table_3,input_table_At_3, input_table_Hv_3, input_table_Lj_3))

SynComs <- c("AtSC","LjSC", "HvSC", "SSC")
hop_2 <- data.frame()

for (KO in input_table_3){
  for (syncom in SynComs){
    syncom_table <- read.table(paste(working_directory, "Isolate_tables/No_dominances/", syncom,"_norm.tsv", sep= ""), sep = "\t", header =T, row.names =1)
    syncom_table_2 <- syncom_table[,grep("ES", colnames(syncom_table))]
    syncom_table_3 <- syncom_table_2[,grep(paste(syncom), colnames(syncom_table_2))]
    syncom_table_4 <- syncom_table_3[,!grepl("HL", colnames(syncom_table_3))]
    
    KO_table_sub <- KO_table[row.names(KO_table) == paste(KO),]
    
    if (syncom == "SSC"){
      KO_table_sub_2 <- KO_table_sub
    } else {
      KO_table_sub_2 <- KO_table_sub[,colnames(KO_table_sub) %in% row.names(tax_df_2)[tax_df_2$SynCom == paste(syncom)]]
    }
    
    KO_table_sub_yes <- names(KO_table_sub_2)[KO_table_sub_2 > 0]
    KO_table_sub_no <- names(KO_table_sub_2)[KO_table_sub_2 == 0]
    
    syncom_table_5 <- t(t(syncom_table_4)/rowSums(t(syncom_table_4)))
    
    syncom_table_At <- syncom_table_5[,grep("At_", colnames(syncom_table_5))]
    syncom_table_Hv <- syncom_table_5[,grep("Hv_", colnames(syncom_table_5))]
    syncom_table_Lj <- syncom_table_5[,grep("Lj_", colnames(syncom_table_5))]
    
    #Averages
    syncom_table_At_2 <- rowSums(syncom_table_At)/length(colnames(syncom_table_At))
    syncom_table_Hv_2 <- rowSums(syncom_table_Hv)/length(colnames(syncom_table_Hv))
    syncom_table_Lj_2 <- rowSums(syncom_table_Lj)/length(colnames(syncom_table_Lj))
    
    At_RA <- sum(syncom_table_At_2[names(syncom_table_At_2) %in% KO_table_sub_yes])
    Hv_RA <- sum(syncom_table_Hv_2[names(syncom_table_Hv_2) %in% KO_table_sub_yes])
    Lj_RA <- sum(syncom_table_Lj_2[names(syncom_table_Lj_2) %in% KO_table_sub_yes])
    
    syncom_table_inp <- syncom_table[,grep("Input", colnames(syncom_table))]
    syncom_table_inp_2 <- t(t(syncom_table_inp)/rowSums(t(syncom_table_inp)))
    
    syncom_table_inp_3  <- rowSums(syncom_table_inp_2)/length(colnames(syncom_table_inp_2))
    
    Input_RA <- sum(syncom_table_inp_3[names(syncom_table_inp_3) %in% KO_table_sub_yes])
    
    if(Input_RA == 0){
      At_val <- At_RA
      Hv_val <- Hv_RA
      Lj_val <- Lj_RA
    } else {
      At_val <- At_RA/Input_RA
      Hv_val <- Hv_RA/Input_RA
      Lj_val <- Lj_RA/Input_RA
    }
    
    len_val <- length(KO_table_sub_yes)/length(KO_table_sub_2)
    
    hop <- t(data.frame(c(paste(KO), At_val, Hv_val, Lj_val, len_val, paste(syncom), length(colnames(syncom_table_At)), length(colnames(syncom_table_Hv)),length(colnames(syncom_table_Lj)))))
    
    hop_2 <- rbind(hop_2, hop)
  }
}

row.names(hop_2) <- NULL
colnames(hop_2) <- c("KO", "At_val", "Hv_val", "Lj_val", "No_of_strains", "SynCom", "No_of_samples_At", "No_of_samples_Hv", "No_of_samples_Lj")

hop_4 <- data.frame()

for (KO in input_table_3){
  hop_sub <- hop_2[hop_2$KO == paste(KO),]
  
  new_3 <- data.frame()
  for (syncom in unique(hop_sub$SynCom)){
    hop_sub_2 <- hop_sub[hop_sub$SynCom == paste(syncom),]
    At_val <- as.numeric(hop_sub_2$At_val)
    Hv_val <- as.numeric(hop_sub_2$Hv_val)
    Lj_val <- as.numeric(hop_sub_2$Lj_val)
    new_2 <- data.frame(At_val,Hv_val, Lj_val,hop_sub_2$No_of_samples_At,hop_sub_2$No_of_samples_Hv, hop_sub_2$No_of_samples_Lj, paste(syncom))
    new_3 <- rbind(new_3,new_2)
  }
  
  syncom_table_inp_3 = apply(syncom_table_inp_2, 1, median, na.rm=TRUE)
  
  
  #Medians
  At_val <- median(as.numeric(new_3$At_val))
  Hv_val <- median(as.numeric(new_3$Hv_val))
  Lj_val <- median(as.numeric(new_3$Lj_val))
  
  No_of_strains <- sum(as.numeric(hop_sub$No_of_strains))/length(hop_sub$No_of_strains)
  
  hop_3 <- t(data.frame(c(paste(KO), At_val,Hv_val, Lj_val, No_of_strains)))
  
  hop_4 <- rbind(hop_4, hop_3)
}

row.names(hop_4) <- NULL
colnames(hop_4) <- c("KO", "Arabidopsis", "Barley", "Lotus", "Proportion_of_strains")

write.table(hop_4, paste(working_directory, "Functionality/266/ternary_KOs_av_med_no_dom_266.txt", sep = ""), quote =F, col.names =T, row.names =T, sep ="\t")


###Script to generate boxplots.txt - - Necessary for 5f and S25 =====
top <- read.table(paste(working_directory, "Annotations/pathway_top.txt", sep = ""), header=F, sep="\t")
KO_to_pathway <- read.table(paste(working_directory, "Annotations/KO_to_pathway.txt", sep = ""), header=T, sep="\t")
KO_to_pathway$V3 <- top$V2[match(KO_to_pathway$V2, top$V1)]

KO_to_pathway_2 <- read.table(paste(working_directory, "Annotations/KO_to_pathway_unannotated_2.txt", sep = ""), header=F, sep="\t")
colnames(KO_to_pathway_2) <- c("KO","new_category")

for (KO in KO_to_pathway_2$KO){
  KO_to_pathway$V3[KO_to_pathway$V1 == paste(KO)] <- KO_to_pathway_2$new_category[KO_to_pathway_2$KO == paste(KO)]
}

KO_table = read.table(paste(working_directory, "KO_genome/KO_SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
colnames(KO_table) <- gsub("X", "", colnames(KO_table))

input_table <- read.table(paste(working_directory, "DESeq2/Sig_KO_all_no_nod_rhizo.txt", sep = ""), header=T, sep="\t")
input_table_2 <- table(input_table$KO)
input_table_3 <- names(input_table_2)[input_table_2 == 12]

pathways <- unique(KO_to_pathway$V3[KO_to_pathway$V1 %in% input_table_3])
pathways_2 <- pathways[pathways != "Unknown"]

Categories <- na.omit(unique(top$V3[top$V2 %in% pathways]))

groups <- Categories
SynComs <- c("AtSC","LjSC", "HvSC", "SSC")
Plant <- c("At", "Hv", "Lj")
hop_4 <- data.frame()

for (cat in groups){
  paths <- top$V2[top$V3 == paste(cat)]
  for (path in paths){
    KOs <- as.vector(na.omit(KO_to_pathway$V1[KO_to_pathway$V3 == paste(path)]))
    KOs_2 <- KOs[KOs %in% input_table_3]
    
    hop_2 <- data.frame()
    if(length(KOs_2) > 0 ){
      for (KO in KOs_2){
        input_table_sub <- input_table[input_table$KO == paste(KO),]
        SynComs <- names(table(input_table_sub$SynCom))[table(input_table_sub$SynCom) > 0]
        for (syncom in SynComs){
          norm_SSC =read.table(paste(working_directory,"Isolate_tables/No_dominances/",syncom,"_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
          norm_SSC_2 <- norm_SSC[,grep("ES", colnames(norm_SSC))]
          norm_SSC_3 <- norm_SSC_2[,grep(paste(syncom), colnames(norm_SSC_2))]
          norm_SSC_input <- norm_SSC[,grep("Input", colnames(norm_SSC))]
          
          KO_table_2 <- KO_table[row.names(KO_table) == paste(KO),]
          KO_table_3 <- KO_table_2[,colnames(KO_table_2) %in% row.names(norm_SSC_3)]
          
          for (plant in Plant){
            norm_SSC_4 <- norm_SSC_3[,grep(paste(plant, "_", sep = ""), colnames(norm_SSC_3))]
            norm_SSC_5 <- t(t(norm_SSC_4)/rowSums(t(norm_SSC_4)))
            
            No <- names(KO_table_3)[KO_table_3 == 0]
            Yes <- names(KO_table_3)[KO_table_3 != 0]
            
            norm_SSC_yes <- norm_SSC_5[row.names(norm_SSC_5) %in% Yes,]
            norm_SSC_no <- norm_SSC_5[row.names(norm_SSC_5) %in% No,]
            
            if (length(Yes) > 1){
              Yes_sum <- sum(colSums(norm_SSC_yes))/length(colSums(norm_SSC_yes))
            } else if (length(Yes) == 1){
              Yes_sum <- sum(norm_SSC_yes)/length(norm_SSC_yes)
            } else {
              Yes_sum <- 0
            }
            
            
            if (length(No) > 1){
              No_sum <- sum(colSums(norm_SSC_no))/length(colSums(norm_SSC_no))
            } else if (length(No) == 1){
              No_sum <- sum(norm_SSC_no)/length(norm_SSC_no)
            } else {
              No_sum <- 0
            }
            
            norm_SSC_input_2 <- t(t(norm_SSC_input)/rowSums(t(norm_SSC_input)))
            
            norm_SSC_input_yes <- norm_SSC_input_2[row.names(norm_SSC_input_2) %in% Yes,]
            norm_SSC_input_no <- norm_SSC_input_2[row.names(norm_SSC_input_2) %in% No,]
            
            if (length(Yes) > 1){
              Yes_sum_input <- sum(colSums(norm_SSC_input_yes))/length(colSums(norm_SSC_input_yes))
            } else if (length(Yes) == 1){
              Yes_sum_input <- sum(norm_SSC_input_yes)/length(norm_SSC_input_yes)
            } else {
              Yes_sum_input <- 0
            }
            
            if (length(No) > 1){
              No_sum_input <- sum(colSums(norm_SSC_input_no))/length(colSums(norm_SSC_input_no))
            } else if (length(No) == 1){
              No_sum_input <- sum(norm_SSC_input_no)/length(norm_SSC_input_no)
            } else {
              No_sum_input <- 0
            }
            
            
            
            hop <- t(data.frame(c(paste(KO), Yes_sum, No_sum, Yes_sum_input, No_sum_input, length(Yes), length(No), paste(plant), paste(syncom))))
            hop_2 <- rbind(hop_2, hop)
          }
        }
      }
      hop_2$V2 <- as.numeric(hop_2$V2)
      hop_2$V3 <- as.numeric(hop_2$V3)
      hop_2$V4 <- as.numeric(hop_2$V4)
      hop_2$V5 <- as.numeric(hop_2$V5)
      hop_2$V6 <- as.numeric(hop_2$V6)
      hop_2$V7 <- as.numeric(hop_2$V7)
      
      int_value <- min(hop_2$V4[hop_2$V4 != 0])
      int_value_2 <- min(hop_2$V5[hop_2$V5 != 0])
      
      hop_2$V4[hop_2$V4 == 0] <- int_value
      hop_2$V5[hop_2$V5 == 0] <- int_value_2
      
      hop_2$V10 <- (hop_2$V2/hop_2$V4)
      hop_2$V11 <- (hop_2$V3/hop_2$V5)
      
      value <- sum(hop_2$V6/(hop_2$V6 + hop_2$V7))/length(hop_2$V6)
      
      for (plant in Plant){
        for (syncom in SynComs){
          hop_2_5 <- hop_2[hop_2$V8 == paste(plant),]
          hop_2_6 <- hop_2_5[hop_2_5$V9 == paste(syncom),]
          value_path_yes <- sum(as.numeric(hop_2_6$V10))/length(hop_2_6$V10)
          value_path_no <- sum(as.numeric(hop_2_6$V11))/length(hop_2_6$V11)
          value_RA <-  sum(as.numeric(hop_2_6$V2))/length(hop_2_6$V2)
          hop_3 <- t(data.frame(c(paste(plant),value_RA, value_path_yes, value_path_no, paste(syncom), paste(path), paste(cat), value)))
          hop_4 <- rbind(hop_4, hop_3)
        }
      }
    }
  }
}

row.names(hop_4) <- NULL

write.table(hop_4, paste(working_directory,"Functionality/266/boxplots.txt", sep = ""), quote = F, col.names = T, row.names = T, sep = "\t")


###Script to generate top70_isolates_no_dom.txt =====
SynComs <- c("AtSC", "HvSC", "LjSC", "SSC")
Plant <- c("At", "Hv", "Lj")

samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names

list_70 <- list()

uno <- 1

for (syncom in SynComs){
  norm_table = read.table(paste(working_directory, "Isolate_tables/No_dominances/",syncom, "_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
  
  #Filter for root samples
  samples_df_sub <- samples_df[samples_df$Inoculum == paste(syncom),]
  samples_df_sub_2 <- samples_df_sub[samples_df_sub$Compartment == "ES",]
  
  for (plant in Plant){
    samples_df_sub_3 <- samples_df_sub_2[samples_df_sub_2$Condition == paste(plant),]
    norm_table_2 <- norm_table[, colnames(norm_table) %in% samples_df_sub_3$sample_id]
    norm_table_3 <- t(t(norm_table_2)/rowSums(t(norm_table_2)))
    
    hop <- data.frame(rowSums(norm_table_3)/length(colnames(norm_table_3)))
    colnames(hop) <- "Rel"
    
    hop_2 <- data.frame(hop[order(hop$Rel, decreasing=T),])
    row.names(hop_2) <- row.names(hop)[order(hop$Rel, decreasing=T)]
    
    colnames(hop_2) <- "Rel"
    hop_2$Cum <- NA
    
    for (i in 1:length(hop_2$Rel)){
      hop_sub <- hop_2[1:i,]
      
      if (length(hop_sub$Rel) == 1){
        hop_2$Cum[i] <- hop_sub$Rel
      } else {
        hop_2$Cum[i] <- sum(hop_sub$Rel)
      }
    }
    
    group_70 <- row.names(hop_2)[1:(length(row.names(hop_2)[hop_2$Cum < 0.7]) +1)]
    
    if (length(group_70) != 0){
      list_70[[uno]] <- group_70
      names(list_70[[uno]]) <- paste(plant,syncom, sep ="_")
    }
    
    uno <- uno + 1
  }
}

# Function to pad entries to the maximum length
pad_to_max_length <- function(x, max_len) {
  length(x) <- max_len
  return(x)
}

max_length_70 <- max(sapply(list_70, length))
df_70 <- do.call(rbind, lapply(list_70, pad_to_max_length, max_len = max_length_70))
df_70 <- as.data.frame(df_70)

write.table(df_70, paste(working_directory, "top70_isolates_no_dom.txt", sep = ""), sep = "\t", quote =F, row.names =F, col.names =F)
###Script to generate top90_isolates_no_dom.txt - Necessary for Figure S30 =====
SynComs <- c("AtSC", "HvSC", "LjSC", "SSC")
Plant <- c("At", "Hv", "Lj")

samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names

list_90 <- list()

uno <- 1

for (syncom in SynComs){
  norm_table = read.table(paste(working_directory, "Isolate_tables/No_dominances/",syncom, "_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
  
  #Filter for root samples
  samples_df_sub <- samples_df[samples_df$Inoculum == paste(syncom),]
  samples_df_sub_2 <- samples_df_sub[samples_df_sub$Compartment == "ES",]
  
  for (plant in Plant){
    samples_df_sub_3 <- samples_df_sub_2[samples_df_sub_2$Condition == paste(plant),]
    norm_table_2 <- norm_table[, colnames(norm_table) %in% samples_df_sub_3$sample_id]
    norm_table_3 <- t(t(norm_table_2)/rowSums(t(norm_table_2)))
    
    hop <- data.frame(rowSums(norm_table_3)/length(colnames(norm_table_3)))
    colnames(hop) <- "Rel"
    
    hop_2 <- data.frame(hop[order(hop$Rel, decreasing=T),])
    row.names(hop_2) <- row.names(hop)[order(hop$Rel, decreasing=T)]
    
    colnames(hop_2) <- "Rel"
    hop_2$Cum <- NA
    
    for (i in 1:length(hop_2$Rel)){
      hop_sub <- hop_2[1:i,]
      
      if (length(hop_sub$Rel) == 1){
        hop_2$Cum[i] <- hop_sub$Rel
      } else {
        hop_2$Cum[i] <- sum(hop_sub$Rel)
      }
    }
    
    group_90 <- row.names(hop_2)[1:(length(row.names(hop_2)[hop_2$Cum < 0.9]) +1)]
    
    if (length(group_90) != 0){
      list_90[[uno]] <- group_90
      names(list_90[[uno]]) <- paste(plant,syncom, sep ="_")
    }
    
    uno <- uno + 1
  }
}

# Function to pad entries to the maximum length
pad_to_max_length <- function(x, max_len) {
  length(x) <- max_len
  return(x)
}

max_length_90 <- max(sapply(list_90, length))
df_90 <- do.call(rbind, lapply(list_90, pad_to_max_length, max_len = max_length_90))
df_90 <- as.data.frame(df_90)

write.table(df_90, paste(working_directory, "top90_isolates_no_dom.txt", sep = ""), sep = "\t", quote =F, row.names =F, col.names =F)

############################# DESeq2 scripts #######################################
###DESeq2 on original dataset=====
SynCom <- c("AtSC", "LjSC", "HvSC", "SSC")

sigtab_col_all_2 <- data.frame(matrix(NA, ncol = 8))
colnames(sigtab_col_all_2) <- c("baseMean", "log2FoldChange","lfcSE", "pvalue", "padj", "plant", "KO","SynCom")
sigtab_col_all_2 <- sigtab_col_all_2[-1,]

for (inoculum in SynCom) {
  #Insert OTU table
  otu_mat = read.table(paste(working_directory,"KO_tables/Original/", inoculum,".tsv", sep = ""), header=T, sep="\t")
  #Insert metadata file
  samples_df = read.table(paste(working_directory,"SSC_R2_metadata.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
  #make the SampleID column into the row.names
  row.names(samples_df) <- samples_df$sample_id
  samples_df <- samples_df %>% dplyr::select (-sample_id)
  row.names(otu_mat) <- otu_mat$function.
  otu_mat <- otu_mat %>% dplyr::select (-function.)
  
  ##DESeq2 Analysis
  otu_mat <- as.matrix(otu_mat)
  
  OTU = otu_table(otu_mat,taxa_are_rows = TRUE)
  samples = sample_data(samples_df)
  samples$Condition <- relevel(factor(samples$Condition), ref = "Input")
  
  #import as phyloseq object
  phylo <- phyloseq(OTU,samples)
  
  phylo_1 = subset_samples(phylo, Compartment != "RZ")
  phylo_2 = subset_samples(phylo_1, Compartment != "BS")
  phylo_3 = subset_samples(phylo_2, Inoculum == paste(inoculum))
  
  #make Deseq2 object, where the column 'genotype' will be used as comparison
  DEseq2_meta_col = phyloseq_to_deseq2(phylo_3, ~Condition)
  
  #work around for error: every gene contains at least one zero, cannot compute log geometric means
  # calculate geometric means prior to estimate size factors
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geoMeans = apply(counts(DEseq2_meta_col), 1, gm_mean)
  
  DEseq2_meta_col = estimateSizeFactors(DEseq2_meta_col, geoMeans = geoMeans)
  
  DEseq2_meta_col = DESeq(DEseq2_meta_col, fitType="local")
  
  #check model fitting of dispersion
  dispersion_plot_col <- plotDispEsts(DEseq2_meta_col)
  
  #investigate DESeq2 results
  resultsNames(DEseq2_meta_col)
  
  #Here is the input for each comparison between mutant vs WT
  res_col_Hv <- lfcShrink(DEseq2_meta_col, coef="Condition_Hv_vs_Input")
  res_col_At <- lfcShrink(DEseq2_meta_col, coef="Condition_At_vs_Input")
  res_col_Lj <- lfcShrink(DEseq2_meta_col, coef="Condition_Lj_vs_Input")
  
  alpha = 0.05
  sigtab_col_Hv = res_col_Hv[which(res_col_Hv$padj < alpha), ]
  sigtab_col_At = res_col_At[which(res_col_At$padj < alpha), ]
  sigtab_col_Lj = res_col_Lj[which(res_col_Lj$padj < alpha), ]
  
  sigtab_col_Hv <- sigtab_col_Hv[sigtab_col_Hv$log2FoldChange >= 0, ]
  sigtab_col_At <- sigtab_col_At[sigtab_col_At$log2FoldChange >= 0, ]
  sigtab_col_Lj <- sigtab_col_Lj[sigtab_col_Lj$log2FoldChange >= 0, ]
  
  sigtab_col_Hv$plant <- "Barley"
  sigtab_col_At$plant <- "Arabidopsis"
  sigtab_col_Lj$plant <- "Lotus"
  
  sigtab_col_Lj$KO <- row.names(sigtab_col_Lj)
  sigtab_col_At$KO <- row.names(sigtab_col_At)
  sigtab_col_Hv$KO <- row.names(sigtab_col_Hv)
  
  sigtab_col_all <- rbind(sigtab_col_Lj, sigtab_col_Hv, sigtab_col_At)
  sigtab_col_all$SynCom <- paste(inoculum)
  sigtab_col_all_2 <- rbind(sigtab_col_all_2,as.data.frame(sigtab_col_all))
}

write.table(sigtab_col_all_2 , file = paste(working_directory, "DESeq2/Sig_KO_all.txt", sep = ""), quote=F, sep="\t", col.names=T, row.names=F)

###DESeq2 on the dataset without dominators =====
SynCom <- c("AtSC", "LjSC", "HvSC", "SSC")

sigtab_col_all_2 <- data.frame(matrix(NA, ncol = 8))
colnames(sigtab_col_all_2) <- c("baseMean", "log2FoldChange","lfcSE", "pvalue", "padj", "plant", "KO","SynCom")
sigtab_col_all_2 <- sigtab_col_all_2[-1,]

for (inoculum in SynCom) {
  #Insert OTU table
  otu_mat = read.table(paste(working_directory,"KO_tables/No_dominances/", inoculum,".tsv", sep = ""), header=T, sep="\t")
  #Insert metadata file
  samples_df = read.table(paste(working_directory,"SSC_R2_metadata.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
  #make the SampleID column into the row.names
  row.names(samples_df) <- samples_df$sample_id
  samples_df <- samples_df %>% dplyr::select (-sample_id)
  row.names(otu_mat) <- otu_mat$function.
  otu_mat <- otu_mat %>% dplyr::select (-function.)
  
  ##DESeq2 Analysis
  otu_mat <- as.matrix(otu_mat)
  
  OTU = otu_table(otu_mat,taxa_are_rows = TRUE)
  samples = sample_data(samples_df)
  samples$Condition <- relevel(factor(samples$Condition), ref = "Input")
  
  #import as phyloseq object
  phylo <- phyloseq(OTU,samples)
  
  phylo_1 = subset_samples(phylo, Compartment != "RZ")
  phylo_2 = subset_samples(phylo_1, Compartment != "BS")
  phylo_3 = subset_samples(phylo_2, Inoculum == paste(inoculum))
  
  #make Deseq2 object, where the column 'genotype' will be used as comparison
  DEseq2_meta_col = phyloseq_to_deseq2(phylo_3, ~Condition)
  
  #work around for error: every gene contains at least one zero, cannot compute log geometric means
  # calculate geometric means prior to estimate size factors
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geoMeans = apply(counts(DEseq2_meta_col), 1, gm_mean)
  
  DEseq2_meta_col = estimateSizeFactors(DEseq2_meta_col, geoMeans = geoMeans)
  
  DEseq2_meta_col = DESeq(DEseq2_meta_col, fitType="local")
  
  #check model fitting of dispersion
  dispersion_plot_col <- plotDispEsts(DEseq2_meta_col)
  
  #investigate DESeq2 results
  resultsNames(DEseq2_meta_col)
  
  #Here is the input for each comparison between mutant vs WT
  res_col_Hv <- lfcShrink(DEseq2_meta_col, coef="Condition_Hv_vs_Input")
  res_col_At <- lfcShrink(DEseq2_meta_col, coef="Condition_At_vs_Input")
  res_col_Lj <- lfcShrink(DEseq2_meta_col, coef="Condition_Lj_vs_Input")
  
  alpha = 0.05
  sigtab_col_Hv = res_col_Hv[which(res_col_Hv$padj < alpha), ]
  sigtab_col_At = res_col_At[which(res_col_At$padj < alpha), ]
  sigtab_col_Lj = res_col_Lj[which(res_col_Lj$padj < alpha), ]
  
  sigtab_col_Hv <- sigtab_col_Hv[sigtab_col_Hv$log2FoldChange >= 0, ]
  sigtab_col_At <- sigtab_col_At[sigtab_col_At$log2FoldChange >= 0, ]
  sigtab_col_Lj <- sigtab_col_Lj[sigtab_col_Lj$log2FoldChange >= 0, ]
  
  sigtab_col_Hv$plant <- "Barley"
  sigtab_col_At$plant <- "Arabidopsis"
  sigtab_col_Lj$plant <- "Lotus"
  
  sigtab_col_Lj$KO <- row.names(sigtab_col_Lj)
  sigtab_col_At$KO <- row.names(sigtab_col_At)
  sigtab_col_Hv$KO <- row.names(sigtab_col_Hv)
  
  sigtab_col_all <- rbind(sigtab_col_Lj, sigtab_col_Hv, sigtab_col_At)
  sigtab_col_all$SynCom <- paste(inoculum)
  sigtab_col_all_2 <- rbind(sigtab_col_all_2,as.data.frame(sigtab_col_all))
}

write.table(sigtab_col_all_2 , file = paste(working_directory, "DESeq2/Sig_KO_all_no_nod_rhizo.txt", sep = ""), quote=F, sep="\t", col.names=T, row.names=F)


###DESeq2 on the dataset without nodulators =====
SynCom <- c("AtSC", "LjSC", "HvSC", "SSC")

sigtab_col_all_2 <- data.frame(matrix(NA, ncol = 8))
colnames(sigtab_col_all_2) <- c("baseMean", "log2FoldChange","lfcSE", "pvalue", "padj", "plant", "KO","SynCom")
sigtab_col_all_2 <- sigtab_col_all_2[-1,]

for (inoculum in SynCom) {
  #Insert OTU table
  otu_mat = read.table(paste(working_directory,"KO_tables/No_nodulators/", inoculum,".tsv", sep = ""), header=T, sep="\t")
  #Insert metadata file
  samples_df = read.table(paste(working_directory,"SSC_R2_metadata.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
  #make the SampleID column into the row.names
  row.names(samples_df) <- samples_df$sample_id
  samples_df <- samples_df %>% dplyr::select (-sample_id)
  row.names(otu_mat) <- otu_mat$function.
  otu_mat <- otu_mat %>% dplyr::select (-function.)
  
  ##DESeq2 Analysis
  otu_mat <- as.matrix(otu_mat)
  
  OTU = otu_table(otu_mat,taxa_are_rows = TRUE)
  samples = sample_data(samples_df)
  samples$Condition <- relevel(factor(samples$Condition), ref = "Input")
  
  #import as phyloseq object
  phylo <- phyloseq(OTU,samples)
  
  phylo_1 = subset_samples(phylo, Compartment != "RZ")
  phylo_2 = subset_samples(phylo_1, Compartment != "BS")
  phylo_3 = subset_samples(phylo_2, Inoculum == paste(inoculum))
  
  #make Deseq2 object, where the column 'genotype' will be used as comparison
  DEseq2_meta_col = phyloseq_to_deseq2(phylo_3, ~Condition)
  
  #work around for error: every gene contains at least one zero, cannot compute log geometric means
  # calculate geometric means prior to estimate size factors
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geoMeans = apply(counts(DEseq2_meta_col), 1, gm_mean)
  
  DEseq2_meta_col = estimateSizeFactors(DEseq2_meta_col, geoMeans = geoMeans)
  
  DEseq2_meta_col = DESeq(DEseq2_meta_col, fitType="local")
  
  #check model fitting of dispersion
  dispersion_plot_col <- plotDispEsts(DEseq2_meta_col)
  
  #investigate DESeq2 results
  resultsNames(DEseq2_meta_col)
  
  #Here is the input for each comparison between mutant vs WT
  res_col_Hv <- lfcShrink(DEseq2_meta_col, coef="Condition_Hv_vs_Input")
  res_col_At <- lfcShrink(DEseq2_meta_col, coef="Condition_At_vs_Input")
  res_col_Lj <- lfcShrink(DEseq2_meta_col, coef="Condition_Lj_vs_Input")
  
  alpha = 0.05
  sigtab_col_Hv = res_col_Hv[which(res_col_Hv$padj < alpha), ]
  sigtab_col_At = res_col_At[which(res_col_At$padj < alpha), ]
  sigtab_col_Lj = res_col_Lj[which(res_col_Lj$padj < alpha), ]
  
  sigtab_col_Hv <- sigtab_col_Hv[sigtab_col_Hv$log2FoldChange >= 0, ]
  sigtab_col_At <- sigtab_col_At[sigtab_col_At$log2FoldChange >= 0, ]
  sigtab_col_Lj <- sigtab_col_Lj[sigtab_col_Lj$log2FoldChange >= 0, ]
  
  sigtab_col_Hv$plant <- "Barley"
  sigtab_col_At$plant <- "Arabidopsis"
  sigtab_col_Lj$plant <- "Lotus"
  
  sigtab_col_Lj$KO <- row.names(sigtab_col_Lj)
  sigtab_col_At$KO <- row.names(sigtab_col_At)
  sigtab_col_Hv$KO <- row.names(sigtab_col_Hv)
  
  sigtab_col_all <- rbind(sigtab_col_Lj, sigtab_col_Hv, sigtab_col_At)
  sigtab_col_all$SynCom <- paste(inoculum)
  sigtab_col_all_2 <- rbind(sigtab_col_all_2,as.data.frame(sigtab_col_all))
}

write.table(sigtab_col_all_2 , file = paste(working_directory, "DESeq2/Sig_KO_all_no_nod.txt", sep = ""), quote=F, sep="\t", col.names=T, row.names=F)

###DESeq2 on the dataset without Rhizobacter =====
SynCom <- c("AtSC", "LjSC", "HvSC", "SSC")

sigtab_col_all_2 <- data.frame(matrix(NA, ncol = 8))
colnames(sigtab_col_all_2) <- c("baseMean", "log2FoldChange","lfcSE", "pvalue", "padj", "plant", "KO","SynCom")
sigtab_col_all_2 <- sigtab_col_all_2[-1,]

for (inoculum in SynCom) {
  #Insert OTU table
  otu_mat = read.table(paste(working_directory,"KO_tables/No_rhizobacter/", inoculum,".tsv", sep = ""), header=T, sep="\t")
  #Insert metadata file
  samples_df = read.table(paste(working_directory,"SSC_R2_metadata.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
  #make the SampleID column into the row.names
  row.names(samples_df) <- samples_df$sample_id
  samples_df <- samples_df %>% dplyr::select (-sample_id)
  row.names(otu_mat) <- otu_mat$function.
  otu_mat <- otu_mat %>% dplyr::select (-function.)
  
  ##DESeq2 Analysis
  otu_mat <- as.matrix(otu_mat)
  
  OTU = otu_table(otu_mat,taxa_are_rows = TRUE)
  samples = sample_data(samples_df)
  samples$Condition <- relevel(factor(samples$Condition), ref = "Input")
  
  #import as phyloseq object
  phylo <- phyloseq(OTU,samples)
  
  phylo_1 = subset_samples(phylo, Compartment != "RZ")
  phylo_2 = subset_samples(phylo_1, Compartment != "BS")
  phylo_3 = subset_samples(phylo_2, Inoculum == paste(inoculum))
  
  #make Deseq2 object, where the column 'genotype' will be used as comparison
  DEseq2_meta_col = phyloseq_to_deseq2(phylo_3, ~Condition)
  
  #work around for error: every gene contains at least one zero, cannot compute log geometric means
  # calculate geometric means prior to estimate size factors
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geoMeans = apply(counts(DEseq2_meta_col), 1, gm_mean)
  
  DEseq2_meta_col = estimateSizeFactors(DEseq2_meta_col, geoMeans = geoMeans)
  
  DEseq2_meta_col = DESeq(DEseq2_meta_col, fitType="local")
  
  #check model fitting of dispersion
  dispersion_plot_col <- plotDispEsts(DEseq2_meta_col)
  
  #investigate DESeq2 results
  resultsNames(DEseq2_meta_col)
  
  #Here is the input for each comparison between mutant vs WT
  res_col_Hv <- lfcShrink(DEseq2_meta_col, coef="Condition_Hv_vs_Input")
  res_col_At <- lfcShrink(DEseq2_meta_col, coef="Condition_At_vs_Input")
  res_col_Lj <- lfcShrink(DEseq2_meta_col, coef="Condition_Lj_vs_Input")
  
  alpha = 0.05
  sigtab_col_Hv = res_col_Hv[which(res_col_Hv$padj < alpha), ]
  sigtab_col_At = res_col_At[which(res_col_At$padj < alpha), ]
  sigtab_col_Lj = res_col_Lj[which(res_col_Lj$padj < alpha), ]
  
  sigtab_col_Hv <- sigtab_col_Hv[sigtab_col_Hv$log2FoldChange >= 0, ]
  sigtab_col_At <- sigtab_col_At[sigtab_col_At$log2FoldChange >= 0, ]
  sigtab_col_Lj <- sigtab_col_Lj[sigtab_col_Lj$log2FoldChange >= 0, ]
  
  sigtab_col_Hv$plant <- "Barley"
  sigtab_col_At$plant <- "Arabidopsis"
  sigtab_col_Lj$plant <- "Lotus"
  
  sigtab_col_Lj$KO <- row.names(sigtab_col_Lj)
  sigtab_col_At$KO <- row.names(sigtab_col_At)
  sigtab_col_Hv$KO <- row.names(sigtab_col_Hv)
  
  sigtab_col_all <- rbind(sigtab_col_Lj, sigtab_col_Hv, sigtab_col_At)
  sigtab_col_all$SynCom <- paste(inoculum)
  sigtab_col_all_2 <- rbind(sigtab_col_all_2,as.data.frame(sigtab_col_all))
}

write.table(sigtab_col_all_2 , file = paste(working_directory, "DESeq2/Sig_KO_all_no_rhizobacter.txt", sep = ""), quote=F, sep="\t", col.names=T, row.names=F)

###DESeq2 on the LjSC Family drop out experiment =====
KO_table <- read.table(paste(working_directory, "LjSC_Family_drop_out_experiment/KO_LjSC_Family_drop_out.tsv", sep = ""), header =T, row.names =1)
KO_table_2 <- KO_table[,!grepl("NOD", colnames(KO_table))]
KO_table_3 <- KO_table_2[, !grepl("LDT7", colnames(KO_table_2))]

#Metadata creation
metadata <- data.frame(colnames(KO_table_3))
colnames(metadata) <- "Sample"
metadata$Compartment <- "Root"
metadata$Compartment[grep("INPUT", metadata$Sample)] <- "Input"

metadata$Subset <- "Full LjSC"
metadata$Subset[grep("LDT1", metadata$Sample)] <- "Burkholderiaceae drop out"
metadata$Subset[grep("LDT2", metadata$Sample)] <- "Caulobacteraceae drop out"
metadata$Subset[grep("LDT3", metadata$Sample)] <- "Pseudomonadaceae drop out"
metadata$Subset[grep("LDT4", metadata$Sample)] <- "Rhizobiaceae drop out"
metadata$Subset[grep("LDT5", metadata$Sample)] <- "All other families drop out"

#DESeq2
sigtab_col_all_2 <- data.frame(matrix(NA, ncol = 7))
colnames(sigtab_col_all_2) <- c("baseMean", "log2FoldChange","lfcSE", "pvalue", "padj", "KO","Subset")
sigtab_col_all_2 <- sigtab_col_all_2[-1,]

for (inoculum in unique(metadata$Subset)) {
  KO_table_2_sub <- KO_table_3[, colnames(KO_table_3) %in% metadata$Sample[metadata$Subset == paste(inoculum)]]
  
  otu_mat <- as.matrix(KO_table_2_sub)
  
  OTU = otu_table(otu_mat,taxa_are_rows = TRUE)
  metadata_2 <- metadata[metadata$Subset == paste(inoculum),]
  row.names(metadata_2) <- metadata_2$Sample
  samples = sample_data(metadata_2)
  samples$Compartment <- relevel(factor(samples$Compartment), ref = "Input")
  
  #import as phyloseq object
  phylo <- phyloseq(OTU,samples)
  
  #make Deseq2 object, where the column 'genotype' will be used as comparison
  DEseq2_meta_col = phyloseq_to_deseq2(phylo, ~Compartment)
  
  #work around for error: every gene contains at least one zero, cannot compute log geometric means
  # calculate geometric means prior to estimate size factors
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geoMeans = apply(counts(DEseq2_meta_col), 1, gm_mean)
  
  DEseq2_meta_col = estimateSizeFactors(DEseq2_meta_col, geoMeans = geoMeans)
  
  DEseq2_meta_col = DESeq(DEseq2_meta_col, fitType="local")
  
  #check model fitting of dispersion
  dispersion_plot_col <- plotDispEsts(DEseq2_meta_col)
  
  #investigate DESeq2 results
  resultsNames(DEseq2_meta_col)
  
  #Here is the input for each comparison between mutant vs WT
  res_col <- lfcShrink(DEseq2_meta_col, coef="Compartment_Root_vs_Input")

  res_col$KO <- row.names(res_col)
  
  res_col$Subset <- paste(inoculum)
  sigtab_col_all_2 <- rbind(sigtab_col_all_2,as.data.frame(res_col))
}

write.table(sigtab_col_all_2 , file = paste(working_directory, "LjSC_Family_drop_out_experiment/DESeq2_Root_vs_input_Fam_drop.txt", sep = ""), quote=F, sep="\t", col.names=T, row.names=F)

###DEseq2 on the LjSC Family drop out experiment - without nodulators =====
KO_table <- read.table(paste(working_directory, "LjSC_Family_drop_out_experiment/KO_LjSC_Family_drop_out_no_nodulators.tsv", sep = ""), header =T, row.names =1)
KO_table_2 <- KO_table[,!grepl("NOD", colnames(KO_table))]
KO_table_3 <- KO_table_2[, !grepl("LDT7", colnames(KO_table_2))]

#Metadata creation
metadata <- data.frame(colnames(KO_table_3))
colnames(metadata) <- "Sample"
metadata$Compartment <- "Root"
metadata$Compartment[grep("INPUT", metadata$Sample)] <- "Input"

metadata$Subset <- "Full LjSC"
metadata$Subset[grep("LDT1", metadata$Sample)] <- "Burkholderiaceae drop out"
metadata$Subset[grep("LDT2", metadata$Sample)] <- "Caulobacteraceae drop out"
metadata$Subset[grep("LDT3", metadata$Sample)] <- "Pseudomonadaceae drop out"
metadata$Subset[grep("LDT4", metadata$Sample)] <- "Rhizobiaceae drop out"
metadata$Subset[grep("LDT5", metadata$Sample)] <- "All other families drop out"

#DESeq2
sigtab_col_all_2 <- data.frame(matrix(NA, ncol = 7))
colnames(sigtab_col_all_2) <- c("baseMean", "log2FoldChange","lfcSE", "pvalue", "padj", "KO","Subset")
sigtab_col_all_2 <- sigtab_col_all_2[-1,]

for (inoculum in unique(metadata$Subset)) {
  KO_table_2_sub <- KO_table_3[, colnames(KO_table_3) %in% metadata$Sample[metadata$Subset == paste(inoculum)]]
  
  otu_mat <- as.matrix(KO_table_2_sub)
  
  OTU = otu_table(otu_mat,taxa_are_rows = TRUE)
  metadata_2 <- metadata[metadata$Subset == paste(inoculum),]
  row.names(metadata_2) <- metadata_2$Sample
  samples = sample_data(metadata_2)
  samples$Compartment <- relevel(factor(samples$Compartment), ref = "Input")
  
  #import as phyloseq object
  phylo <- phyloseq(OTU,samples)
  
  #make Deseq2 object, where the column 'genotype' will be used as comparison
  DEseq2_meta_col = phyloseq_to_deseq2(phylo, ~Compartment)
  
  #work around for error: every gene contains at least one zero, cannot compute log geometric means
  # calculate geometric means prior to estimate size factors
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geoMeans = apply(counts(DEseq2_meta_col), 1, gm_mean)
  
  DEseq2_meta_col = estimateSizeFactors(DEseq2_meta_col, geoMeans = geoMeans)
  
  DEseq2_meta_col = DESeq(DEseq2_meta_col, fitType="local")
  
  #check model fitting of dispersion
  dispersion_plot_col <- plotDispEsts(DEseq2_meta_col)
  
  #investigate DESeq2 results
  resultsNames(DEseq2_meta_col)
  
  #Here is the input for each comparison between mutant vs WT
  res_col <- lfcShrink(DEseq2_meta_col, coef="Compartment_Root_vs_Input")
  
  res_col$KO <- row.names(res_col)
  
  res_col$Subset <- paste(inoculum)
  sigtab_col_all_2 <- rbind(sigtab_col_all_2,as.data.frame(res_col))
}

write.table(sigtab_col_all_2 , file = paste(working_directory, "LjSC_Family_drop_out_experiment/DESeq2_Root_vs_input_Fam_drop_no_nod.txt", sep = ""), quote=F, sep="\t", col.names=T, row.names=F)

###################### Family R2 scripts - Figure 4a ##############################
###Frequency isolates per Family across SynComs =====
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")

top <- read.table(paste(working_directory, "top70_isolates_no_dom.txt", sep = ""), sep = "\t", header =F)
top_2 <- unlist(as.vector(top))
top_3 <- top_2[!is.na(top_2)]
top_4 <- unique(top_3)

fams_left <- unique(tax_df_2$Family[row.names(tax_df_2) %in% top_4])

SynComs <- c("AtSC","HvSC", "LjSC")

all_together <- data.frame()

for(syncom in SynComs){
  tax_df_2_SC <- tax_df_2[tax_df_2$SynCom == paste(syncom),]
  
  new <- data.frame(table(tax_df_2_SC$Family))
  new$SynCom <- paste(syncom)
  
  all_together <- rbind(all_together, new)
}

all_together_2 <- all_together[all_together$Var1 %in% fams_left,]

family_order <- c("Chitinophagaceae", "Microbacteriaceae","Micrococcaceae","Xanthobacteraceae","Sphingobacteriaceae","Rhodanobacteraceae","Sphingomonadaceae","Flavobacteriaceae","Devosiaceae","Beijerinckiaceae","Enterobacteriaceae","Caulobacteraceae","Pseudomonadaceae","Xanthomonadaceae","Burkholderiaceae","Rhizobiaceae")

all_together_2$SynCom <- factor(all_together_2$SynCom, levels = c("LjSC","HvSC","AtSC"))
all_together_2$Var1 <- factor(all_together_2$Var1, levels = family_order)

write.table(all_together_2, paste(working_directory, "Family_R2/No_of_isolates_per_fam.txt", sep = ""), col.names =T, row.names =F, quote =F, sep = "\t")

###Frequency isolates per Genus across SynComs - Burkholderiaceae genera =====
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")

fams_left <- c("Acidovorax","Cupriavidus","Pelomonas","Polaromonas","Rhizobacter","Variovorax")

SynComs <- c("AtSC","HvSC", "LjSC")

all_together <- data.frame()

for(syncom in SynComs){
  tax_df_2_SC <- tax_df_2[tax_df_2$SynCom == paste(syncom),]
  
  new <- data.frame(table(tax_df_2_SC$Genus))
  new$SynCom <- paste(syncom)
  
  all_together <- rbind(all_together, new)
}

all_together_2 <- all_together[all_together$Var1 %in% fams_left,]

family_order <- c("Pelomonas","Cupriavidus","Polaromonas","Variovorax","Rhizobacter","Acidovorax")

all_together_2$SynCom <- factor(all_together_2$SynCom, levels = c("LjSC","HvSC","AtSC"))
all_together_2$Var1 <- factor(all_together_2$Var1, levels = family_order)

write.table(all_together_2, paste(working_directory, "Family_R2/No_of_isolates_per_Burkholderiaceae_genus.txt", sep = ""), col.names =T, row.names =F, quote =F, sep = "\t")

###Distance to centroid per family and SynCom =====

#KO profiles
all_iso_2 =read.table(paste(working_directory,"KO_genome/KO_SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Taxonomy
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")

colnames(all_iso_2)[grep("M.16",colnames(all_iso_2))] <- "M-16"
colnames(all_iso_2)[grep("M.6",colnames(all_iso_2))] <- "M-6"
colnames(all_iso_2)[grep("M.10",colnames(all_iso_2))] <- "M-10"
colnames(all_iso_2)[grep("M.11_2",colnames(all_iso_2))] <- "M-11_2"
colnames(all_iso_2)[grep("M.11",colnames(all_iso_2))] <- "M-11"
colnames(all_iso_2) <- gsub("X", "", colnames(all_iso_2))

hopla_2 <- data.frame()

SynComs <- c("AtSC", "HvSC", "LjSC")

for (syncom in SynComs){
  tax_sub <- tax_df_2[tax_df_2$SynCom == paste(syncom),]
  
  #Set the OTU, TAX and sample data for making phyloseq object
  OTU = otu_table(as.matrix(all_iso_2),taxa_are_rows = TRUE)
  samples = sample_data(tax_sub)
  
  phylo_sub = phyloseq(OTU,samples)
  
  phylo_sub_RA=microbiome::transform(x = phylo_sub, transform = "compositional" )
  
  #Bray Curtis distance matrix
  beta <- as.matrix(vegdist(t(phylo_sub_RA@otu_table@.Data), method = "bray", diag = T))
  row.names(beta) <- gsub("X", "", row.names(beta))
  
  #Make PCoA plot for Bray Curtis Distance matrix
  pcoa_tax = cmdscale(beta, k=3, eig=T)
  points_beta = as.data.frame(pcoa_tax$points)
  colnames(points_beta) = c("x", "y", "z") 
  eig = pcoa_tax$eig
  row.names(points_beta) <- gsub("X", "", row.names(points_beta))
  points_beta = merge(points_beta,tax_df_2, by = "row.names")
  row.names(points_beta) <- points_beta$Row.names
  points_beta <- points_beta %>% dplyr::select (-Row.names)
  
  # Define parameter for centroid calculation
  
  param="Family"
  
  # Calculate centroids for each group
  centroids_tax <- points_beta %>%
    group_by(!!sym(param)) %>%
    dplyr::summarize(
      centroid_x = mean(x, na.rm = TRUE),
      centroid_y = mean(y, na.rm = TRUE),
      centroid_z = mean(z, na.rm = TRUE)
    )
  
  # Join centroids back to the original data
  data_with_centroids_tax <- left_join(points_beta, centroids_tax , by = param)
  # Calculate distance to centroid for each point
  data_with_centroids_tax  <- data_with_centroids_tax  %>% 
    rowwise() %>%
    mutate(distance_to_centroid = sqrt((x - centroid_x)^2 + 
                                         (y - centroid_y)^2 + 
                                         (z - centroid_z)^2))
  
  data_with_centroids_tax$distance_to_centroid
  data_with_centroids_tax_2 <- as.data.frame(data_with_centroids_tax)
  row.names(data_with_centroids_tax_2) <- row.names(points_beta)
  
  for (family in unique(data_with_centroids_tax_2$Family)){
    data_with_centroids_tax_2_sub <- data_with_centroids_tax_2[data_with_centroids_tax_2$Family == paste(family),]
    average <- sum(data_with_centroids_tax_2_sub$distance_to_centroid)/length(data_with_centroids_tax_2_sub$distance_to_centroid)
    
    hopla <- data.frame(t(data.frame(c(paste(family), average, paste(syncom)))))
    hopla_2 <- rbind(hopla_2, hopla)
  }
}

#70% CRA filter
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")

top <- read.table(paste(working_directory, "top70_isolates_no_dom.txt", sep = ""), sep = "\t", header =F)
top_2 <- unlist(as.vector(top))
top_3 <- top_2[!is.na(top_2)]
top_4 <- unique(top_3)

fams_left <- unique(tax_df_2$Family[row.names(tax_df_2) %in% top_4])

row.names(hopla_2) <- NULL
colnames(hopla_2) <- c("Family", "Distance_to_Centroid", "SynCom")

hopla_3 <- hopla_2[hopla_2$Family %in% fams_left,]

write.table(hopla_3, paste(working_directory, "Family_R2/Dist_to_centroid_Fam.txt", sep = ""), col.names =T, row.names =F, quote =F, sep = "\t")

###Distance to centroid per genus and SynCom - Burkholderiaceae =====

#KO profiles
all_iso_2 =read.table(paste(working_directory,"KO_genome/KO_SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Taxonomy
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")

colnames(all_iso_2)[grep("M.16",colnames(all_iso_2))] <- "M-16"
colnames(all_iso_2)[grep("M.6",colnames(all_iso_2))] <- "M-6"
colnames(all_iso_2)[grep("M.10",colnames(all_iso_2))] <- "M-10"
colnames(all_iso_2)[grep("M.11_2",colnames(all_iso_2))] <- "M-11_2"
colnames(all_iso_2)[grep("M.11",colnames(all_iso_2))] <- "M-11"
colnames(all_iso_2) <- gsub("X", "", colnames(all_iso_2))

hopla_2 <- data.frame()

SynComs <- c("AtSC", "HvSC", "LjSC")

for (syncom in SynComs){
  tax_sub <- tax_df_2[tax_df_2$SynCom == paste(syncom),]
  
  #Set the OTU, TAX and sample data for making phyloseq object
  OTU = otu_table(as.matrix(all_iso_2),taxa_are_rows = TRUE)
  samples = sample_data(tax_sub)
  
  phylo_sub = phyloseq(OTU,samples)
  
  phylo_sub_RA=microbiome::transform(x = phylo_sub, transform = "compositional" )
  
  #Bray Curtis distance matrix
  beta <- as.matrix(vegdist(t(phylo_sub_RA@otu_table@.Data), method = "bray", diag = T))
  row.names(beta) <- gsub("X", "", row.names(beta))
  
  #Make PCoA plot for Bray Curtis Distance matrix
  pcoa_tax = cmdscale(beta, k=3, eig=T)
  points_beta = as.data.frame(pcoa_tax$points)
  colnames(points_beta) = c("x", "y", "z") 
  eig = pcoa_tax$eig
  row.names(points_beta) <- gsub("X", "", row.names(points_beta))
  points_beta = merge(points_beta,tax_df_2, by = "row.names")
  row.names(points_beta) <- points_beta$Row.names
  points_beta <- points_beta %>% dplyr::select (-Row.names)
  
  # Define parameter for centroid calculation
  
  param="Genus"
  
  # Calculate centroids for each group
  centroids_tax <- points_beta %>%
    group_by(!!sym(param)) %>%
    dplyr::summarize(
      centroid_x = mean(x, na.rm = TRUE),
      centroid_y = mean(y, na.rm = TRUE),
      centroid_z = mean(z, na.rm = TRUE)
    )
  
  # Join centroids back to the original data
  data_with_centroids_tax <- left_join(points_beta, centroids_tax , by = param)
  # Calculate distance to centroid for each point
  data_with_centroids_tax  <- data_with_centroids_tax  %>% 
    rowwise() %>%
    mutate(distance_to_centroid = sqrt((x - centroid_x)^2 + 
                                         (y - centroid_y)^2 + 
                                         (z - centroid_z)^2))
  
  data_with_centroids_tax$distance_to_centroid
  data_with_centroids_tax_2 <- as.data.frame(data_with_centroids_tax)
  row.names(data_with_centroids_tax_2) <- row.names(points_beta)
  
  for (family in unique(data_with_centroids_tax_2$Genus)){
    data_with_centroids_tax_2_sub <- data_with_centroids_tax_2[data_with_centroids_tax_2$Genus == paste(family),]
    data_with_centroids_tax_2_sub <- data_with_centroids_tax_2_sub[!is.na(data_with_centroids_tax_2_sub$x),]
    average <- sum(data_with_centroids_tax_2_sub$distance_to_centroid)/length(data_with_centroids_tax_2_sub$distance_to_centroid)
    
    hopla <- data.frame(t(data.frame(c(paste(family), average, paste(syncom)))))
    hopla_2 <- rbind(hopla_2, hopla)
  }
}

tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")

fams_left <- c("Acidovorax","Cupriavidus","Pelomonas","Polaromonas","Rhizobacter","Variovorax")

row.names(hopla_2) <- NULL
colnames(hopla_2) <- c("Genus", "Distance_to_Centroid", "SynCom")

hopla_3 <- hopla_2[hopla_2$Genus %in% fams_left,]

write.table(hopla_3, paste(working_directory, "Family_R2/Dist_to_centroid_Burkholderiaceae_genus.txt", sep = ""), col.names =T, row.names =F, quote =F, sep = "\t")

###Family effect - SynCom & Plant R2 - dom - generation of file =====
#otu table
KO_SSC =read.table(paste(working_directory,"KO_tables/Original/SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Metadata
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object
OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)

#Sample subsetting
cond="ES"

samples_df_sub <- subset(samples_df, samples_df$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

SynComs <- c("AtSC", "HvSC", "LjSC", "SSC")
Hosts <- c("At","Hv", "Lj")

for (syncom in SynComs){
  samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum == paste(syncom))
  
  samples_sub = sample_data(samples_df_sub_3)
  
  phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
  
  phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
  
  beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))
  
  Bray_curtis_df=beta_isolate_KO
  
  #Make PCoA plot for Bray Curtis Distance matrix
  pcoa = cmdscale(Bray_curtis_df, k=3, eig=T)
  points = as.data.frame(pcoa$points)
  colnames(points) = c("x", "y", "z") 
  eig = pcoa$eig
  
  points = merge(points,samples_df_sub_2, by = "row.names")
  rownames(points) <- points$Row.names
  points <- points %>% dplyr::select (-Row.names)
  
  points$Condition <- factor(points$Condition, levels = c("At","Hv", "Lj"))
  points$Nutrient <- factor(points$Nutrient, levels = c("low", "high"))
  points$Experiment  <- factor(points$Experiment, levels = c("R1", "R2"))
  
  metadata=points[,-c(1,2,3)]
  
  set.seed(1)
  SSC_bray_KO_adonis <- adonis2(beta_isolate_KO ~ Condition*Nutrient*Experiment, data=metadata, method="bray", permutations=999)
  
  colnames(points)[colnames(points) == "Condition"] <- "Plant"
  
  R2_Plant_KO <- SSC_bray_KO_adonis$R2[1]
  assign(paste("R2_Plant_KO_",syncom,sep =""),R2_Plant_KO)
}

for (host in Hosts){
  samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Condition == paste(host))
  
  samples_sub = sample_data(samples_df_sub_3)
  
  phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
  
  phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
  
  beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))
  
  Bray_curtis_df=beta_isolate_KO
  
  #Make PCoA plot for Bray Curtis Distance matrix
  pcoa = cmdscale(Bray_curtis_df, k=3, eig=T)
  points = as.data.frame(pcoa$points)
  colnames(points) = c("x", "y", "z") 
  eig = pcoa$eig
  
  points = merge(points,samples_df_sub_2, by = "row.names")
  rownames(points) <- points$Row.names
  points <- points %>% dplyr::select (-Row.names)
  
  points$Condition <- factor(points$Inoculum, levels = c("AtSC","HvSC", "LjSC", "SSC"))
  points$Nutrient <- factor(points$Nutrient, levels = c("low", "high"))
  points$Experiment  <- factor(points$Experiment, levels = c("R1", "R2"))
  
  metadata=points[,-c(1,2,3)]
  
  set.seed(1)
  SSC_bray_KO_adonis <- adonis2(beta_isolate_KO ~ Inoculum*Nutrient*Experiment, data=metadata, method="bray", permutations=999)
  
  colnames(points)[colnames(points) == "Inoculum"] <- "SynCom"
  
  R2_SynCom_KO <- SSC_bray_KO_adonis$R2[1]
  assign(paste("R2_SynCom_KO_",host,sep =""),R2_SynCom_KO)
}

#Family
list_SSC <- read.table(paste(working_directory, "Family_R2/SSC_list.txt", sep = ""), sep = "\t", header =F)

list_SSC_or <- list_SSC[list_SSC$V1 == "SSC",]
list_SSC_rm <- list_SSC[list_SSC$V1 != "SSC",]

fam_data <- data.frame()
list_SSC_or$V2

tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")

top <- read.table(paste(working_directory, "top70_isolates_no_dom.txt", sep = ""), sep = "\t", header =F)
top_2 <- unlist(as.vector(top))
top_3 <- top_2[!is.na(top_2)]
top_4 <- unique(top_3)

fams_left <- unique(tax_df_2$Family[row.names(tax_df_2) %in% top_4])

SynComs <- c("AtSC","HvSC","LjSC","SSC")
Hosts <- c("At","Hv","Lj")

for (family in fams_left){
  #KOs
  genera_data <- read.table(paste(working_directory, "Family_R2/Family_dom/",family,".tsv", sep =""), sep ="\t", header =T, row.names =1)
  
  #Samples TABLE
  samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t", row.names =1) #make the SampleID column into the row.names
  colnames(samples_df)[5]="Nutrient"
  samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
  samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")
  
  #Set the OTU, TAX and sample data for making phyloseq object
  OTU = otu_table(as.matrix(genera_data),taxa_are_rows = TRUE)
  
  #Sample subsetting
  
  cond="ES"
  samples_df_sub <- subset(samples_df, samples_df$Compartment == cond)
  samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")
  
  #Plant effect
  for (syncom in SynComs){
    samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum == paste(syncom))
    samples_sub = sample_data(samples_df_sub_3)
    
    phylo_sub = phyloseq(OTU, samples_sub)
    
    phylo_sub_RA=microbiome::transform(x = phylo_sub, transform = "compositional" )
    
    #Agglomerate to phylum-level and rename
    #Bray Curtis distance matrix
    beta_genus <- as.matrix(vegdist(t(phylo_sub_RA@otu_table@.Data), method = "bray", diag = T))
    mean_value_genus= mean(beta_genus)
    
    #Make PCoA plot for Bray Curtis Distance matrix
    pcoa = cmdscale(beta_genus, k=3, eig=T)
    points = as.data.frame(pcoa$points)
    colnames(points) = c("x", "y", "z") 
    eig = pcoa$eig
    
    points = merge(points,samples_df_sub_2, by = "row.names")
    rownames(points) <- points$Row.names
    points <- points %>% dplyr::select (-Row.names)
    
    points$Condition <- factor(points$Condition, levels = c("At","Hv", "Lj"))
    points$Nutrient <- factor(points$Nutrient, levels = c("low", "high"))
    points$Experiment  <- factor(points$Experiment, levels = c("R1", "R2"))
    
    metadata=points[,-c(1,2,3)]
    
    #  Run adonis PERMANOVA test
    set.seed(1)
    SSC_bray_adonis_fam <- adonis2(beta_genus ~ Condition*Nutrient*Experiment, data=metadata, method="bray", permutations=999)
    
    colnames(points)[colnames(points) == "Condition"] <- "Plant"
    
    if (paste(syncom) == "AtSC"){
      R2_Plant_KO <- R2_Plant_KO_AtSC
    } else if (paste(syncom) == "HvSC"){
      R2_Plant_KO <- R2_Plant_KO_HvSC
    } else if (paste(syncom) == "LjSC"){
      R2_Plant_KO <- R2_Plant_KO_LjSC
    } else {
      R2_Plant_KO <- R2_Plant_KO_SSC
    }
    
    R2_Plant_fam <- SSC_bray_adonis_fam$R2[1] - R2_Plant_KO
    
    fam_data_2 <- data.frame(t(data.frame(c(paste(family), paste(syncom),R2_Plant_fam,"KO"))))
    
    fam_data <- rbind(fam_data, fam_data_2)
  }
  
  for (host in Hosts){
    samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Condition == paste(host))
    samples_sub = sample_data(samples_df_sub_3)
    
    phylo_sub = phyloseq(OTU, samples_sub)
    
    phylo_sub_RA=microbiome::transform(x = phylo_sub, transform = "compositional" )
    
    #Agglomerate to phylum-level and rename
    #Bray Curtis distance matrix
    beta_genus <- as.matrix(vegdist(t(phylo_sub_RA@otu_table@.Data), method = "bray", diag = T))
    mean_value_genus= mean(beta_genus)
    
    #Make PCoA plot for Bray Curtis Distance matrix
    pcoa = cmdscale(beta_genus, k=3, eig=T)
    points = as.data.frame(pcoa$points)
    colnames(points) = c("x", "y", "z") 
    eig = pcoa$eig
    
    points = merge(points,samples_df_sub_2, by = "row.names")
    rownames(points) <- points$Row.names
    points <- points %>% dplyr::select (-Row.names)
    
    points$Inoculum <- factor(points$Inoculum, levels = c("AtSC","HvSC", "LjSC", "SSC"))
    points$Nutrient <- factor(points$Nutrient, levels = c("low", "high"))
    points$Experiment  <- factor(points$Experiment, levels = c("R1", "R2"))
    
    metadata=points[,-c(1,2,3)]
    
    #  Run adonis PERMANOVA test
    set.seed(1)
    SSC_bray_adonis_fam <- adonis2(beta_genus ~ Inoculum*Nutrient*Experiment, data=metadata, method="bray", permutations=999)
    
    colnames(points)[colnames(points) == "Inoculum"] <- "SynCom"
    
    if (paste(host) == "At"){
      R2_SynCom <- R2_SynCom_KO_At
    } else if (paste(host) == "Hv"){
      R2_SynCom <- R2_SynCom_KO_Hv
    } else {
      R2_SynCom <- R2_SynCom_KO_Lj
    }
    
    R2_SynCom_fam <- SSC_bray_adonis_fam$R2[1] - R2_SynCom
    fam_data_2 <- data.frame(t(data.frame(c(paste(family), paste(host),R2_SynCom_fam, "KO"))))
    
    fam_data <- rbind(fam_data, fam_data_2)
  }
}

row.names(fam_data) <- NULL
colnames(fam_data) <- c("Family", "Subset","R2_change","KO")

write.table(fam_data, paste(working_directory, "Family_R2/SSC_Fam_R2_effects_subs_with_dom.txt", sep = ""), sep ="\t", quote =F, col.names =T, row.names =T)

###Family effect - SynCom & Plant R2 - no dom - generation of file =====
#otu table
KO_SSC =read.table(paste(working_directory,"KO_tables/Original/SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Metadata
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object
OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)

#Sample subsetting
cond="ES"

samples_df_sub <- subset(samples_df, samples_df$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

SynComs <- c("AtSC", "HvSC", "LjSC", "SSC")
Hosts <- c("At","Hv", "Lj")

for (syncom in SynComs){
  samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum == paste(syncom))
  
  samples_sub = sample_data(samples_df_sub_3)
  
  phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
  
  phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
  
  beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))
  
  Bray_curtis_df=beta_isolate_KO
  
  #Make PCoA plot for Bray Curtis Distance matrix
  pcoa = cmdscale(Bray_curtis_df, k=3, eig=T)
  points = as.data.frame(pcoa$points)
  colnames(points) = c("x", "y", "z") 
  eig = pcoa$eig
  
  points = merge(points,samples_df_sub_2, by = "row.names")
  rownames(points) <- points$Row.names
  points <- points %>% dplyr::select (-Row.names)
  
  points$Condition <- factor(points$Condition, levels = c("At","Hv", "Lj"))
  points$Nutrient <- factor(points$Nutrient, levels = c("low", "high"))
  points$Experiment  <- factor(points$Experiment, levels = c("R1", "R2"))
  
  metadata=points[,-c(1,2,3)]
  
  set.seed(1)
  SSC_bray_KO_adonis <- adonis2(beta_isolate_KO ~ Condition*Nutrient*Experiment, data=metadata, method="bray", permutations=999)
  
  colnames(points)[colnames(points) == "Condition"] <- "Plant"
  
  R2_Plant_KO <- SSC_bray_KO_adonis$R2[1]
  assign(paste("R2_Plant_KO_",syncom,sep =""),R2_Plant_KO)
}

for (host in Hosts){
  samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Condition == paste(host))
  
  samples_sub = sample_data(samples_df_sub_3)
  
  phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
  
  phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
  
  beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))
  
  Bray_curtis_df=beta_isolate_KO
  
  #Make PCoA plot for Bray Curtis Distance matrix
  pcoa = cmdscale(Bray_curtis_df, k=3, eig=T)
  points = as.data.frame(pcoa$points)
  colnames(points) = c("x", "y", "z") 
  eig = pcoa$eig
  
  points = merge(points,samples_df_sub_2, by = "row.names")
  rownames(points) <- points$Row.names
  points <- points %>% dplyr::select (-Row.names)
  
  points$Condition <- factor(points$Inoculum, levels = c("AtSC","HvSC", "LjSC", "SSC"))
  points$Nutrient <- factor(points$Nutrient, levels = c("low", "high"))
  points$Experiment  <- factor(points$Experiment, levels = c("R1", "R2"))
  
  metadata=points[,-c(1,2,3)]
  
  set.seed(1)
  SSC_bray_KO_adonis <- adonis2(beta_isolate_KO ~ Inoculum*Nutrient*Experiment, data=metadata, method="bray", permutations=999)
  
  colnames(points)[colnames(points) == "Inoculum"] <- "SynCom"
  
  R2_SynCom_KO <- SSC_bray_KO_adonis$R2[1]
  assign(paste("R2_SynCom_KO_",host,sep =""),R2_SynCom_KO)
}

#Family
list_SSC <- read.table(paste(working_directory, "Family_R2/SSC_list.txt", sep = ""), sep = "\t", header =F)

list_SSC_or <- list_SSC[list_SSC$V1 == "SSC",]
list_SSC_rm <- list_SSC[list_SSC$V1 != "SSC",]

fam_data <- data.frame()
list_SSC_or$V2

tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")

top <- read.table(paste(working_directory, "top70_isolates_no_dom.txt", sep = ""), sep = "\t", header =F)
top_2 <- unlist(as.vector(top))
top_3 <- top_2[!is.na(top_2)]
top_4 <- unique(top_3)

fams_left <- unique(tax_df_2$Family[row.names(tax_df_2) %in% top_4])

SynComs <- c("AtSC","HvSC","LjSC","SSC")
Hosts <- c("At","Hv","Lj")

for (family in fams_left){
  #KOs
  genera_data <- read.table(paste(working_directory, "Family_R2/Family_no_dom/",family,".tsv", sep =""), sep ="\t", header =T, row.names =1)
  
  #Samples TABLE
  samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t", row.names =1) #make the SampleID column into the row.names
  colnames(samples_df)[5]="Nutrient"
  samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
  samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")
  
  #Set the OTU, TAX and sample data for making phyloseq object
  OTU = otu_table(as.matrix(genera_data),taxa_are_rows = TRUE)
  
  #Sample subsetting
  
  cond="ES"
  samples_df_sub <- subset(samples_df, samples_df$Compartment == cond)
  samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")
  
  #Plant effect
  for (syncom in SynComs){
    samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum == paste(syncom))
    samples_sub = sample_data(samples_df_sub_3)
    
    phylo_sub = phyloseq(OTU, samples_sub)
    
    phylo_sub_RA=microbiome::transform(x = phylo_sub, transform = "compositional" )
    
    #Agglomerate to phylum-level and rename
    #Bray Curtis distance matrix
    beta_genus <- as.matrix(vegdist(t(phylo_sub_RA@otu_table@.Data), method = "bray", diag = T))
    mean_value_genus= mean(beta_genus)
    
    #Make PCoA plot for Bray Curtis Distance matrix
    pcoa = cmdscale(beta_genus, k=3, eig=T)
    points = as.data.frame(pcoa$points)
    colnames(points) = c("x", "y", "z") 
    eig = pcoa$eig
    
    points = merge(points,samples_df_sub_2, by = "row.names")
    rownames(points) <- points$Row.names
    points <- points %>% dplyr::select (-Row.names)
    
    points$Condition <- factor(points$Condition, levels = c("At","Hv", "Lj"))
    points$Nutrient <- factor(points$Nutrient, levels = c("low", "high"))
    points$Experiment  <- factor(points$Experiment, levels = c("R1", "R2"))
    
    metadata=points[,-c(1,2,3)]
    
    #  Run adonis PERMANOVA test
    set.seed(1)
    SSC_bray_adonis_fam <- adonis2(beta_genus ~ Condition*Nutrient*Experiment, data=metadata, method="bray", permutations=999)
    
    colnames(points)[colnames(points) == "Condition"] <- "Plant"
    
    if (paste(syncom) == "AtSC"){
      R2_Plant_KO <- R2_Plant_KO_AtSC
    } else if (paste(syncom) == "HvSC"){
      R2_Plant_KO <- R2_Plant_KO_HvSC
    } else if (paste(syncom) == "LjSC"){
      R2_Plant_KO <- R2_Plant_KO_LjSC
    } else {
      R2_Plant_KO <- R2_Plant_KO_SSC
    }
    
    R2_Plant_fam <- SSC_bray_adonis_fam$R2[1] - R2_Plant_KO
    
    fam_data_2 <- data.frame(t(data.frame(c(paste(family), paste(syncom),R2_Plant_fam,"KO"))))
    
    fam_data <- rbind(fam_data, fam_data_2)
  }
  
  for (host in Hosts){
    samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Condition == paste(host))
    samples_sub = sample_data(samples_df_sub_3)
    
    phylo_sub = phyloseq(OTU, samples_sub)
    
    phylo_sub_RA=microbiome::transform(x = phylo_sub, transform = "compositional" )
    
    #Agglomerate to phylum-level and rename
    #Bray Curtis distance matrix
    beta_genus <- as.matrix(vegdist(t(phylo_sub_RA@otu_table@.Data), method = "bray", diag = T))
    mean_value_genus= mean(beta_genus)
    
    #Make PCoA plot for Bray Curtis Distance matrix
    pcoa = cmdscale(beta_genus, k=3, eig=T)
    points = as.data.frame(pcoa$points)
    colnames(points) = c("x", "y", "z") 
    eig = pcoa$eig
    
    points = merge(points,samples_df_sub_2, by = "row.names")
    rownames(points) <- points$Row.names
    points <- points %>% dplyr::select (-Row.names)
    
    points$Inoculum <- factor(points$Inoculum, levels = c("AtSC","HvSC", "LjSC", "SSC"))
    points$Nutrient <- factor(points$Nutrient, levels = c("low", "high"))
    points$Experiment  <- factor(points$Experiment, levels = c("R1", "R2"))
    
    metadata=points[,-c(1,2,3)]
    
    #  Run adonis PERMANOVA test
    set.seed(1)
    SSC_bray_adonis_fam <- adonis2(beta_genus ~ Inoculum*Nutrient*Experiment, data=metadata, method="bray", permutations=999)
    
    colnames(points)[colnames(points) == "Inoculum"] <- "SynCom"
    
    if (paste(host) == "At"){
      R2_SynCom <- R2_SynCom_KO_At
    } else if (paste(host) == "Hv"){
      R2_SynCom <- R2_SynCom_KO_Hv
    } else {
      R2_SynCom <- R2_SynCom_KO_Lj
    }
    
    R2_SynCom_fam <- SSC_bray_adonis_fam$R2[1] - R2_SynCom
    fam_data_2 <- data.frame(t(data.frame(c(paste(family), paste(host),R2_SynCom_fam, "KO"))))
    
    fam_data <- rbind(fam_data, fam_data_2)
  }
}

row.names(fam_data) <- NULL
colnames(fam_data) <- c("Family", "Subset","R2_change","KO")

write.table(fam_data, paste(working_directory, "Family_R2/SSC_Fam_R2_effects_subs.txt", sep = ""), sep ="\t", quote =F, col.names =T, row.names =T)

###Family effect - SynCom & Plant R2 - dom - Burkholderiaceae genera - generation of file =====
#otu table
KO_SSC =read.table(paste(working_directory,"KO_tables/Original/SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Metadata
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object
OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)

#Sample subsetting
cond="ES"

samples_df_sub <- subset(samples_df, samples_df$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

SynComs <- c("AtSC", "HvSC", "LjSC", "SSC")
Hosts <- c("At","Hv", "Lj")

for (syncom in SynComs){
  samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum == paste(syncom))
  
  samples_sub = sample_data(samples_df_sub_3)
  
  phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
  
  phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
  
  beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))
  
  Bray_curtis_df=beta_isolate_KO
  
  #Make PCoA plot for Bray Curtis Distance matrix
  pcoa = cmdscale(Bray_curtis_df, k=3, eig=T)
  points = as.data.frame(pcoa$points)
  colnames(points) = c("x", "y", "z") 
  eig = pcoa$eig
  
  points = merge(points,samples_df_sub_2, by = "row.names")
  rownames(points) <- points$Row.names
  points <- points %>% dplyr::select (-Row.names)
  
  points$Condition <- factor(points$Condition, levels = c("At","Hv", "Lj"))
  points$Nutrient <- factor(points$Nutrient, levels = c("low", "high"))
  points$Experiment  <- factor(points$Experiment, levels = c("R1", "R2"))
  
  metadata=points[,-c(1,2,3)]
  
  set.seed(1)
  SSC_bray_KO_adonis <- adonis2(beta_isolate_KO ~ Condition*Nutrient*Experiment, data=metadata, method="bray", permutations=999)
  
  colnames(points)[colnames(points) == "Condition"] <- "Plant"
  
  R2_Plant_KO <- SSC_bray_KO_adonis$R2[1]
  assign(paste("R2_Plant_KO_",syncom,sep =""),R2_Plant_KO)
}

for (host in Hosts){
  samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Condition == paste(host))
  
  samples_sub = sample_data(samples_df_sub_3)
  
  phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
  
  phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
  
  beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))
  
  Bray_curtis_df=beta_isolate_KO
  
  #Make PCoA plot for Bray Curtis Distance matrix
  pcoa = cmdscale(Bray_curtis_df, k=3, eig=T)
  points = as.data.frame(pcoa$points)
  colnames(points) = c("x", "y", "z") 
  eig = pcoa$eig
  
  points = merge(points,samples_df_sub_2, by = "row.names")
  rownames(points) <- points$Row.names
  points <- points %>% dplyr::select (-Row.names)
  
  points$Condition <- factor(points$Inoculum, levels = c("AtSC","HvSC", "LjSC", "SSC"))
  points$Nutrient <- factor(points$Nutrient, levels = c("low", "high"))
  points$Experiment  <- factor(points$Experiment, levels = c("R1", "R2"))
  
  metadata=points[,-c(1,2,3)]
  
  set.seed(1)
  SSC_bray_KO_adonis <- adonis2(beta_isolate_KO ~ Inoculum*Nutrient*Experiment, data=metadata, method="bray", permutations=999)
  
  colnames(points)[colnames(points) == "Inoculum"] <- "SynCom"
  
  R2_SynCom_KO <- SSC_bray_KO_adonis$R2[1]
  assign(paste("R2_SynCom_KO_",host,sep =""),R2_SynCom_KO)
}

#Family
fam_data <- data.frame()

tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")

fams_left <- c("Acidovorax","Cupriavidus","Pelomonas","Polaromonas","Rhizobacter","Variovorax")

SynComs <- c("AtSC","HvSC","LjSC","SSC")
Hosts <- c("At","Hv","Lj")

for (family in fams_left){
  #KOs
  genera_data <- read.table(paste(working_directory, "Family_R2/Genus_dom/",family,".tsv", sep =""), sep ="\t", header =T, row.names =1)
  
  #Samples TABLE
  samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t", row.names =1) #make the SampleID column into the row.names
  colnames(samples_df)[5]="Nutrient"
  samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
  samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")
  
  #Set the OTU, TAX and sample data for making phyloseq object
  OTU = otu_table(as.matrix(genera_data),taxa_are_rows = TRUE)
  
  #Sample subsetting
  
  cond="ES"
  samples_df_sub <- subset(samples_df, samples_df$Compartment == cond)
  samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")
  
  #Plant effect
  for (syncom in SynComs){
    samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum == paste(syncom))
    samples_sub = sample_data(samples_df_sub_3)
    
    phylo_sub = phyloseq(OTU, samples_sub)
    
    phylo_sub_RA=microbiome::transform(x = phylo_sub, transform = "compositional" )
    
    #Agglomerate to phylum-level and rename
    #Bray Curtis distance matrix
    beta_genus <- as.matrix(vegdist(t(phylo_sub_RA@otu_table@.Data), method = "bray", diag = T))
    mean_value_genus= mean(beta_genus)
    
    #Make PCoA plot for Bray Curtis Distance matrix
    pcoa = cmdscale(beta_genus, k=3, eig=T)
    points = as.data.frame(pcoa$points)
    colnames(points) = c("x", "y", "z") 
    eig = pcoa$eig
    
    points = merge(points,samples_df_sub_2, by = "row.names")
    rownames(points) <- points$Row.names
    points <- points %>% dplyr::select (-Row.names)
    
    points$Condition <- factor(points$Condition, levels = c("At","Hv", "Lj"))
    points$Nutrient <- factor(points$Nutrient, levels = c("low", "high"))
    points$Experiment  <- factor(points$Experiment, levels = c("R1", "R2"))
    
    metadata=points[,-c(1,2,3)]
    
    #  Run adonis PERMANOVA test
    set.seed(1)
    SSC_bray_adonis_fam <- adonis2(beta_genus ~ Condition*Nutrient*Experiment, data=metadata, method="bray", permutations=999)
    
    colnames(points)[colnames(points) == "Condition"] <- "Plant"
    
    if (paste(syncom) == "AtSC"){
      R2_Plant_KO <- R2_Plant_KO_AtSC
    } else if (paste(syncom) == "HvSC"){
      R2_Plant_KO <- R2_Plant_KO_HvSC
    } else if (paste(syncom) == "LjSC"){
      R2_Plant_KO <- R2_Plant_KO_LjSC
    } else {
      R2_Plant_KO <- R2_Plant_KO_SSC
    }
    
    R2_Plant_fam <- SSC_bray_adonis_fam$R2[1] - R2_Plant_KO
    
    fam_data_2 <- data.frame(t(data.frame(c(paste(family), paste(syncom),R2_Plant_fam,"KO"))))
    
    fam_data <- rbind(fam_data, fam_data_2)
  }
  
  for (host in Hosts){
    samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Condition == paste(host))
    samples_sub = sample_data(samples_df_sub_3)
    
    phylo_sub = phyloseq(OTU, samples_sub)
    
    phylo_sub_RA=microbiome::transform(x = phylo_sub, transform = "compositional" )
    
    #Agglomerate to phylum-level and rename
    #Bray Curtis distance matrix
    beta_genus <- as.matrix(vegdist(t(phylo_sub_RA@otu_table@.Data), method = "bray", diag = T))
    mean_value_genus= mean(beta_genus)
    
    #Make PCoA plot for Bray Curtis Distance matrix
    pcoa = cmdscale(beta_genus, k=3, eig=T)
    points = as.data.frame(pcoa$points)
    colnames(points) = c("x", "y", "z") 
    eig = pcoa$eig
    
    points = merge(points,samples_df_sub_2, by = "row.names")
    rownames(points) <- points$Row.names
    points <- points %>% dplyr::select (-Row.names)
    
    points$Inoculum <- factor(points$Inoculum, levels = c("AtSC","HvSC", "LjSC", "SSC"))
    points$Nutrient <- factor(points$Nutrient, levels = c("low", "high"))
    points$Experiment  <- factor(points$Experiment, levels = c("R1", "R2"))
    
    metadata=points[,-c(1,2,3)]
    
    #  Run adonis PERMANOVA test
    set.seed(1)
    SSC_bray_adonis_fam <- adonis2(beta_genus ~ Inoculum*Nutrient*Experiment, data=metadata, method="bray", permutations=999)
    
    colnames(points)[colnames(points) == "Inoculum"] <- "SynCom"
    
    if (paste(host) == "At"){
      R2_SynCom <- R2_SynCom_KO_At
    } else if (paste(host) == "Hv"){
      R2_SynCom <- R2_SynCom_KO_Hv
    } else {
      R2_SynCom <- R2_SynCom_KO_Lj
    }
    
    R2_SynCom_fam <- SSC_bray_adonis_fam$R2[1] - R2_SynCom
    fam_data_2 <- data.frame(t(data.frame(c(paste(family), paste(host),R2_SynCom_fam, "KO"))))
    
    fam_data <- rbind(fam_data, fam_data_2)
  }
}

row.names(fam_data) <- NULL
colnames(fam_data) <- c("Family", "Subset","R2_change","KO")

write.table(fam_data, paste(working_directory, "Family_R2/SSC_Gen_Burk_R2_effects_subs_with_dom.txt", sep = ""), sep ="\t", quote =F, col.names =T, row.names =T)

###Family effect - SynCom & Plant R2 - no dom - Burkholderiaceae genera - generation of file =====
#otu table
KO_SSC =read.table(paste(working_directory,"KO_tables/Original/SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Metadata
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object
OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)

#Sample subsetting
cond="ES"

samples_df_sub <- subset(samples_df, samples_df$Compartment == cond)
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

SynComs <- c("AtSC", "HvSC", "LjSC", "SSC")
Hosts <- c("At","Hv", "Lj")

for (syncom in SynComs){
  samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum == paste(syncom))
  
  samples_sub = sample_data(samples_df_sub_3)
  
  phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
  
  phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
  
  beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))
  
  Bray_curtis_df=beta_isolate_KO
  
  #Make PCoA plot for Bray Curtis Distance matrix
  pcoa = cmdscale(Bray_curtis_df, k=3, eig=T)
  points = as.data.frame(pcoa$points)
  colnames(points) = c("x", "y", "z") 
  eig = pcoa$eig
  
  points = merge(points,samples_df_sub_2, by = "row.names")
  rownames(points) <- points$Row.names
  points <- points %>% dplyr::select (-Row.names)
  
  points$Condition <- factor(points$Condition, levels = c("At","Hv", "Lj"))
  points$Nutrient <- factor(points$Nutrient, levels = c("low", "high"))
  points$Experiment  <- factor(points$Experiment, levels = c("R1", "R2"))
  
  metadata=points[,-c(1,2,3)]
  
  set.seed(1)
  SSC_bray_KO_adonis <- adonis2(beta_isolate_KO ~ Condition*Nutrient*Experiment, data=metadata, method="bray", permutations=999)
  
  colnames(points)[colnames(points) == "Condition"] <- "Plant"
  
  R2_Plant_KO <- SSC_bray_KO_adonis$R2[1]
  assign(paste("R2_Plant_KO_",syncom,sep =""),R2_Plant_KO)
}

for (host in Hosts){
  samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Condition == paste(host))
  
  samples_sub = sample_data(samples_df_sub_3)
  
  phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
  
  phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
  
  beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))
  
  Bray_curtis_df=beta_isolate_KO
  
  #Make PCoA plot for Bray Curtis Distance matrix
  pcoa = cmdscale(Bray_curtis_df, k=3, eig=T)
  points = as.data.frame(pcoa$points)
  colnames(points) = c("x", "y", "z") 
  eig = pcoa$eig
  
  points = merge(points,samples_df_sub_2, by = "row.names")
  rownames(points) <- points$Row.names
  points <- points %>% dplyr::select (-Row.names)
  
  points$Condition <- factor(points$Inoculum, levels = c("AtSC","HvSC", "LjSC", "SSC"))
  points$Nutrient <- factor(points$Nutrient, levels = c("low", "high"))
  points$Experiment  <- factor(points$Experiment, levels = c("R1", "R2"))
  
  metadata=points[,-c(1,2,3)]
  
  set.seed(1)
  SSC_bray_KO_adonis <- adonis2(beta_isolate_KO ~ Inoculum*Nutrient*Experiment, data=metadata, method="bray", permutations=999)
  
  colnames(points)[colnames(points) == "Inoculum"] <- "SynCom"
  
  R2_SynCom_KO <- SSC_bray_KO_adonis$R2[1]
  assign(paste("R2_SynCom_KO_",host,sep =""),R2_SynCom_KO)
}

#Family
fam_data <- data.frame()

tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")

fams_left <- c("Acidovorax","Cupriavidus","Pelomonas","Polaromonas","Rhizobacter","Variovorax")

SynComs <- c("AtSC","HvSC","LjSC","SSC")
Hosts <- c("At","Hv","Lj")

for (family in fams_left){
  #KOs
  genera_data <- read.table(paste(working_directory, "Family_R2/Genus_no_dom/",family,".tsv", sep =""), sep ="\t", header =T, row.names =1)
  
  #Samples TABLE
  samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t", row.names =1) #make the SampleID column into the row.names
  colnames(samples_df)[5]="Nutrient"
  samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
  samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")
  
  #Set the OTU, TAX and sample data for making phyloseq object
  OTU = otu_table(as.matrix(genera_data),taxa_are_rows = TRUE)
  
  #Sample subsetting
  
  cond="ES"
  samples_df_sub <- subset(samples_df, samples_df$Compartment == cond)
  samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")
  
  #Plant effect
  for (syncom in SynComs){
    samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum == paste(syncom))
    samples_sub = sample_data(samples_df_sub_3)
    
    phylo_sub = phyloseq(OTU, samples_sub)
    
    phylo_sub_RA=microbiome::transform(x = phylo_sub, transform = "compositional" )
    
    #Agglomerate to phylum-level and rename
    #Bray Curtis distance matrix
    beta_genus <- as.matrix(vegdist(t(phylo_sub_RA@otu_table@.Data), method = "bray", diag = T))
    mean_value_genus= mean(beta_genus)
    
    #Make PCoA plot for Bray Curtis Distance matrix
    pcoa = cmdscale(beta_genus, k=3, eig=T)
    points = as.data.frame(pcoa$points)
    colnames(points) = c("x", "y", "z") 
    eig = pcoa$eig
    
    points = merge(points,samples_df_sub_2, by = "row.names")
    rownames(points) <- points$Row.names
    points <- points %>% dplyr::select (-Row.names)
    
    points$Condition <- factor(points$Condition, levels = c("At","Hv", "Lj"))
    points$Nutrient <- factor(points$Nutrient, levels = c("low", "high"))
    points$Experiment  <- factor(points$Experiment, levels = c("R1", "R2"))
    
    metadata=points[,-c(1,2,3)]
    
    #  Run adonis PERMANOVA test
    set.seed(1)
    SSC_bray_adonis_fam <- adonis2(beta_genus ~ Condition*Nutrient*Experiment, data=metadata, method="bray", permutations=999)
    
    colnames(points)[colnames(points) == "Condition"] <- "Plant"
    
    if (paste(syncom) == "AtSC"){
      R2_Plant_KO <- R2_Plant_KO_AtSC
    } else if (paste(syncom) == "HvSC"){
      R2_Plant_KO <- R2_Plant_KO_HvSC
    } else if (paste(syncom) == "LjSC"){
      R2_Plant_KO <- R2_Plant_KO_LjSC
    } else {
      R2_Plant_KO <- R2_Plant_KO_SSC
    }
    
    R2_Plant_fam <- SSC_bray_adonis_fam$R2[1] - R2_Plant_KO
    
    fam_data_2 <- data.frame(t(data.frame(c(paste(family), paste(syncom),R2_Plant_fam,"KO"))))
    
    fam_data <- rbind(fam_data, fam_data_2)
  }
  
  for (host in Hosts){
    samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Condition == paste(host))
    samples_sub = sample_data(samples_df_sub_3)
    
    phylo_sub = phyloseq(OTU, samples_sub)
    
    phylo_sub_RA=microbiome::transform(x = phylo_sub, transform = "compositional" )
    
    #Agglomerate to phylum-level and rename
    #Bray Curtis distance matrix
    beta_genus <- as.matrix(vegdist(t(phylo_sub_RA@otu_table@.Data), method = "bray", diag = T))
    mean_value_genus= mean(beta_genus)
    
    #Make PCoA plot for Bray Curtis Distance matrix
    pcoa = cmdscale(beta_genus, k=3, eig=T)
    points = as.data.frame(pcoa$points)
    colnames(points) = c("x", "y", "z") 
    eig = pcoa$eig
    
    points = merge(points,samples_df_sub_2, by = "row.names")
    rownames(points) <- points$Row.names
    points <- points %>% dplyr::select (-Row.names)
    
    points$Inoculum <- factor(points$Inoculum, levels = c("AtSC","HvSC", "LjSC", "SSC"))
    points$Nutrient <- factor(points$Nutrient, levels = c("low", "high"))
    points$Experiment  <- factor(points$Experiment, levels = c("R1", "R2"))
    
    metadata=points[,-c(1,2,3)]
    
    #  Run adonis PERMANOVA test
    set.seed(1)
    SSC_bray_adonis_fam <- adonis2(beta_genus ~ Inoculum*Nutrient*Experiment, data=metadata, method="bray", permutations=999)
    
    colnames(points)[colnames(points) == "Inoculum"] <- "SynCom"
    
    if (paste(host) == "At"){
      R2_SynCom <- R2_SynCom_KO_At
    } else if (paste(host) == "Hv"){
      R2_SynCom <- R2_SynCom_KO_Hv
    } else {
      R2_SynCom <- R2_SynCom_KO_Lj
    }
    
    R2_SynCom_fam <- SSC_bray_adonis_fam$R2[1] - R2_SynCom
    fam_data_2 <- data.frame(t(data.frame(c(paste(family), paste(host),R2_SynCom_fam, "KO"))))
    
    fam_data <- rbind(fam_data, fam_data_2)
  }
}

row.names(fam_data) <- NULL
colnames(fam_data) <- c("Family", "Subset","R2_change","KO")

write.table(fam_data, paste(working_directory, "Family_R2/SSC_Gen_Burk_R2_effects_subs.txt", sep = ""), sep ="\t", quote =F, col.names =T, row.names =T)

############################# sPLS-DA scripts ######################################
###sPLS-DA - Rhizobiaceae with dominators and excl AtSC - Plants =====
KO_SSC=read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Rhizobiaceae_KO_with_dom.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t", row.names =1) #make the SampleID column into the row.names
colnames(samples_df)[5]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object

#Sample subsetting
samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

#At subset - Enterobacteriaceae - dom vs no-dom 
samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum != "AtSC")

OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)
samples_sub = sample_data(samples_df_sub_3)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

bray_2 <- as.matrix(beta_isolate_KO)

str(samples_df_sub_3)
str(bray_2)

#Bind metadata with distance matrix
pcoa = cmdscale(bray_2, k=10, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z", "a", "b", "c", "d", "e", "f", "g") 
eig = pcoa$eig
points_2 <- points[order(row.names(points)), ]
samples_df_sub_6 <- samples_df_sub_3[row.names(samples_df_sub_3) %in% row.names(points),]
samples_df_sub_7 <- samples_df_sub_6[order(row.names(samples_df_sub_6)), ]
points_3 <- cbind(points_2,samples_df_sub_7)
colnames(points_3) <- c("x", "y", "z","a", "b", "c", "d", "e", "f", "g",colnames(samples_df_sub_7))

# Principal Coordinates Analysis (PCoA)
pcoa_result <- cmdscale(bray_2, eig = TRUE, k = 2)  # k = number of dimensions

# Extract coordinates for plotting
pcoa_coords <- as.data.frame(pcoa_result$points)

# Adding metadata environmental data
pcoa_coords$Inoculum <- samples_df_sub_7$Inoculum[match(row.names(pcoa_coords), row.names(samples_df_sub_7))]

KO_SSC_3 <- KO_SSC[,colnames(KO_SSC) %in% row.names(pcoa_coords)]
KO_SSC_4 <- t(t(KO_SSC_3)/rowSums(t(KO_SSC_3)))

data_plot <- data.frame(matrix(NA, ncol = 5))
colnames(data_plot) <- c("KO", "contrib", "component", "Cluster", "Study")
data_plot_2 <- data_plot[-1,]

#mixOmics
KO_SSC_5 <- t(KO_SSC_4)
KO_SSC_6 <- KO_SSC_5[match(row.names(points_3),row.names(KO_SSC_5)),]
KO_SSC_7 <- KO_SSC_6[, colSums(KO_SSC_6 != 0) > 0]

#remove columns with 0's
final.plsda <- plsda(KO_SSC_7,points_3$Condition, ncomp = 10)

set.seed(30) # For reproducibility with this handbook, remove otherwise
perf.plsda <- perf(final.plsda, validation = 'Mfold', folds = 3, 
                   progressBar = FALSE,  # Set to TRUE to track progress
                   nrepeat = 10)     

list.keepX <- c(1:10,  seq(20, 100, 10))
list.keepX

tune.splsda_data <- tune.splsda(KO_SSC_7,points_3$Condition, ncomp = 4, validation = 'Mfold', 
                                folds = 5, dist = 'max.dist', 
                                test.keepX = list.keepX, nrepeat = 10)
ncomp <- tune.splsda_data$choice.ncomp$ncomp 
select.keepX <- tune.splsda_data$choice.keepX[1:ncomp]  

splsda.data <- splsda(KO_SSC_7,points_3$Condition, ncomp = ncomp, keepX = select.keepX) 
var.name.short <- colnames(KO_SSC_7)

empty_vector <- vector()

empty_vector_all <- data.frame(matrix(NA, ncol = 4))
colnames(empty_vector_all) <- c("KO", "contrib", "component", "Data")
empty_vector_all_2 <- empty_vector_all[-1,]

empty_vector_contrib <- data_frame()

for (comp in 1:ncomp) {
  list <- as.data.frame(splsda.data$loadings$X)
  nonredun <- row.names(list)[rowSums(list)!=0]
  nonredun_2 <- list[row.names(list) %in% nonredun,]
  nonredun_3 <- as.data.frame(cbind(nonredun, nonredun_2))
  comp_2 <- comp+1
  if (comp == 1) {
    PC1 <- nonredun_3[comp_2]
    PC1$KO <- nonredun_3$nonredun
  } else {
    PC1 <- nonredun_3[comp_2]
    PC1$KO <- row.names(PC1)
  }
  colnames(PC1) <- c(paste("comp",comp,sep=""),"KO")
  PC1_2 <- PC1[order(PC1[,paste("comp",comp,sep="")]),]
  PC1_3 <- PC1_2[PC1_2[,-1] !=0,]
  PC1_3 <- PC1_3[PC1_3[1] != 0,]
  PC1_3$extra <- abs(as.numeric(unlist(PC1_3[1])))
  
  empty_vector <- c(empty_vector,PC1_3$KO)
  
  table <- data.frame(PC1_3$KO)
  inter <- as.data.frame(do.call(cbind, PC1_3[3]))
  
  table$contrib <- unlist(inter[1])
  colnames(table) <- c("KO", "contrib")
  table$component <- paste("PC", comp, sep="")
  
  plot <- plotLoadings(splsda.data, comp = comp, method = 'mean', contrib = 'max', 
                       name.var = var.name.short)
  
  new_data <- data.frame(plot$X$GroupContrib)
  row.names(new_data) <- row.names(plot$X)
  empty_vector_contrib <- rbind(empty_vector_contrib,new_data)
  
  table$Data <- new_data$plot.X.GroupContrib[match(table$KO,row.names(new_data))]
  empty_vector_all_2 <- rbind(empty_vector_all_2, table)
  
}

#otu table
KO_SSC=read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Rhizobiaceae_KO_with_dom.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t", row.names =1) #make the SampleID column into the row.names
colnames(samples_df)[5]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object

#Sample subsetting
samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

#At subset - Enterobacteriaceae - dom vs no-dom 
samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum == "AtSC")

OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)
samples_sub = sample_data(samples_df_sub_3)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

bray_2 <- as.matrix(beta_isolate_KO)

str(samples_df_sub_3)
str(bray_2)

#Bind metadata with distance matrix
pcoa = cmdscale(bray_2, k=10, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z", "a", "b", "c", "d", "e", "f", "g") 
eig = pcoa$eig

empty_vector_all_2$Variance <- NA
empty_vector_all_2$Variance[empty_vector_all_2$component == "PC1"] <- eig[1]/sum(eig)
empty_vector_all_2$Variance[empty_vector_all_2$component == "PC2"] <- eig[2]/sum(eig)
empty_vector_all_2$Variance[empty_vector_all_2$component == "PC3"] <- eig[3]/sum(eig)
empty_vector_all_2$Variance[empty_vector_all_2$component == "PC4"] <- eig[4]/sum(eig)
empty_vector_all_2$Group <- "Rhizobiaceae_Plant"

write.table(empty_vector_all_2, paste(working_directory, "sPLS-DA/output/PLSDA_Rhizobiaceae_with_dom.tsv", sep=""), sep = "\t", quote = F, row.names = F, col.names = T)

###sPLS-DA - Burkholderiaceae without dominators - SynComs =====
#otu table
KO_SSC=read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Burkholderiaceae_KO_no_dom.tsv", sep =""), header=TRUE,sep="\t", row.names = 1)

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t", row.names =1) #make the SampleID column into the row.names
colnames(samples_df)[5]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object

#Sample subsetting
samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

#At subset - Enterobacteriaceae - dom vs no-dom 

OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)
samples_sub = sample_data(samples_df_sub_2)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

bray_2 <- as.matrix(beta_isolate_KO)

str(samples_df_sub_4)
str(bray_2)

#Bind metadata with distance matrix
pcoa = cmdscale(bray_2, k=10, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z", "a", "b", "c", "d", "e", "f", "g") 
eig = pcoa$eig
points_2 <- points[order(row.names(points)), ]
samples_df_sub_3 <- samples_df_sub_2[row.names(samples_df_sub_2) %in% row.names(points),]
samples_df_sub_4 <- samples_df_sub_3[order(row.names(samples_df_sub_3)), ]
points_3 <- cbind(points_2,samples_df_sub_4)
colnames(points_3) <- c("x", "y", "z","a", "b", "c", "d", "e", "f", "g",colnames(samples_df_sub_4))


# Principal Coordinates Analysis (PCoA)
pcoa_result <- cmdscale(bray_2, eig = TRUE, k = 2)  # k = number of dimensions

# Extract coordinates for plotting
pcoa_coords <- as.data.frame(pcoa_result$points)

# Adding metadata environmental data
pcoa_coords$Inoculum <- samples_df_sub_7$Inoculum[match(row.names(pcoa_coords), row.names(samples_df_sub_7))]

KO_SSC_3 <- KO_SSC[,colnames(KO_SSC) %in% row.names(pcoa_coords)]
KO_SSC_4 <- t(t(KO_SSC_3)/rowSums(t(KO_SSC_3)))

data_plot <- data.frame(matrix(NA, ncol = 5))
colnames(data_plot) <- c("KO", "contrib", "component", "Cluster", "Study")
data_plot_2 <- data_plot[-1,]

#mixOmics
KO_SSC_5 <- t(KO_SSC_4)
KO_SSC_6 <- KO_SSC_5[match(row.names(points_3),row.names(KO_SSC_5)),]
KO_SSC_7 <- KO_SSC_6[, colSums(KO_SSC_6 != 0) > 0]

#remove columns with 0's
final.plsda <- plsda(KO_SSC_7,points_3$Inoculum, ncomp = 10)

set.seed(30) # For reproducibility with this handbook, remove otherwise
perf.plsda <- perf(final.plsda, validation = 'Mfold', folds = 3, 
                   progressBar = FALSE,  # Set to TRUE to track progress
                   nrepeat = 10)     

list.keepX <- c(1:10,  seq(20, 100, 10))
list.keepX

tune.splsda_data <- tune.splsda(KO_SSC_7,points_3$Inoculum, ncomp = 4, validation = 'Mfold', 
                                folds = 5, dist = 'max.dist', 
                                test.keepX = list.keepX, nrepeat = 10)
ncomp <- tune.splsda_data$choice.ncomp$ncomp 
select.keepX <- tune.splsda_data$choice.keepX[1:ncomp]  

splsda.data <- splsda(KO_SSC_7,points_3$Inoculum, ncomp = ncomp, keepX = select.keepX) 
var.name.short <- colnames(KO_SSC_7)

empty_vector <- vector()

empty_vector_all <- data.frame(matrix(NA, ncol = 4))
colnames(empty_vector_all) <- c("KO", "contrib", "component", "Data")
empty_vector_all_2 <- empty_vector_all[-1,]

empty_vector_contrib <- data_frame()

for (comp in 1:ncomp) {
  list <- as.data.frame(splsda.data$loadings$X)
  nonredun <- row.names(list)[rowSums(list)!=0]
  nonredun_2 <- list[row.names(list) %in% nonredun,]
  nonredun_3 <- as.data.frame(cbind(nonredun, nonredun_2))
  comp_2 <- comp+1
  if (comp == 1) {
    PC1 <- nonredun_3[comp_2]
    PC1$KO <- nonredun_3$nonredun
  } else {
    PC1 <- nonredun_3[comp_2]
    PC1$KO <- row.names(PC1)
  }
  colnames(PC1) <- c(paste("comp",comp,sep=""),"KO")
  PC1_2 <- PC1[order(PC1[,paste("comp",comp,sep="")]),]
  PC1_3 <- PC1_2[PC1_2[,-1] !=0,]
  PC1_3 <- PC1_3[PC1_3[1] != 0,]
  PC1_3$extra <- abs(as.numeric(unlist(PC1_3[1])))
  
  empty_vector <- c(empty_vector,PC1_3$KO)
  
  table <- data.frame(PC1_3$KO)
  inter <- as.data.frame(do.call(cbind, PC1_3[3]))
  
  table$contrib <- unlist(inter[1])
  colnames(table) <- c("KO", "contrib")
  table$component <- paste("PC", comp, sep="")
  
  plot <- plotLoadings(splsda.data, comp = comp, method = 'mean', contrib = 'max', 
                       name.var = var.name.short)
  
  new_data <- data.frame(plot$X$GroupContrib)
  row.names(new_data) <- row.names(plot$X)
  empty_vector_contrib <- rbind(empty_vector_contrib,new_data)
  
  table$Data <- new_data$plot.X.GroupContrib[match(table$KO,row.names(new_data))]
  empty_vector_all_2 <- rbind(empty_vector_all_2, table)
  
}

#otu table
KO_SSC=read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Burkholderiaceae_KO_no_dom.tsv", sep =""), header=TRUE,sep="\t", row.names = 1)

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t", row.names =1) #make the SampleID column into the row.names
colnames(samples_df)[5]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object

#Sample subsetting
samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

#At subset - Enterobacteriaceae - dom vs no-dom 
OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)
samples_sub = sample_data(samples_df_sub_2)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

bray_2 <- as.matrix(beta_isolate_KO)

str(samples_df_sub_2)
str(bray_2)

#Bind metadata with distance matrix
pcoa = cmdscale(bray_2, k=10, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z", "a", "b", "c", "d", "e", "f", "g") 
eig = pcoa$eig

empty_vector_all_2$Variance <- NA
empty_vector_all_2$Variance[empty_vector_all_2$component == "PC1"] <- eig[1]/sum(eig)
empty_vector_all_2$Variance[empty_vector_all_2$component == "PC2"] <- eig[2]/sum(eig)
empty_vector_all_2$Variance[empty_vector_all_2$component == "PC3"] <- eig[3]/sum(eig)
empty_vector_all_2$Group <- "Burkholderiaceae_SynCom"

write.table(empty_vector_all_2, paste(working_directory, "sPLS-DA/output/PLSDA_Burkholderiaceae_no_dom_SynCom.tsv", sep=""), sep = "\t", quote = F, row.names = F, col.names = T)

###sPLS-DA - Xanthomonadaceae without dominators - SynComs =====
#otu table
KO_SSC=read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Xanthomonadaceae_KO_no_dom.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t", row.names =1) #make the SampleID column into the row.names
colnames(samples_df)[5]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object

#Sample subsetting
samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

#At subset - Enterobacteriaceae - dom vs no-dom 

OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)
samples_sub = sample_data(samples_df_sub_2)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

bray_2 <- as.matrix(beta_isolate_KO)

str(samples_df_sub_4)
str(bray_2)

#Bind metadata with distance matrix
pcoa = cmdscale(bray_2, k=10, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z", "a", "b", "c", "d", "e", "f", "g") 
eig = pcoa$eig
points_2 <- points[order(row.names(points)), ]
samples_df_sub_3 <- samples_df_sub_2[row.names(samples_df_sub_2) %in% row.names(points),]
samples_df_sub_4 <- samples_df_sub_3[order(row.names(samples_df_sub_3)), ]
points_3 <- cbind(points_2,samples_df_sub_4)
colnames(points_3) <- c("x", "y", "z","a", "b", "c", "d", "e", "f", "g",colnames(samples_df_sub_4))

# Principal Coordinates Analysis (PCoA)
pcoa_result <- cmdscale(bray_2, eig = TRUE, k = 2)  # k = number of dimensions

# Extract coordinates for plotting
pcoa_coords <- as.data.frame(pcoa_result$points)

# Adding metadata environmental data
pcoa_coords$Inoculum <- samples_df_sub_7$Inoculum[match(row.names(pcoa_coords), row.names(samples_df_sub_7))]

KO_SSC_3 <- KO_SSC[,colnames(KO_SSC) %in% row.names(pcoa_coords)]
KO_SSC_4 <- t(t(KO_SSC_3)/rowSums(t(KO_SSC_3)))

data_plot <- data.frame(matrix(NA, ncol = 5))
colnames(data_plot) <- c("KO", "contrib", "component", "Cluster", "Study")
data_plot_2 <- data_plot[-1,]

#mixOmics
KO_SSC_5 <- t(KO_SSC_4)
KO_SSC_6 <- KO_SSC_5[match(row.names(points_3),row.names(KO_SSC_5)),]
KO_SSC_7 <- KO_SSC_6[, colSums(KO_SSC_6 != 0) > 0]

#remove columns with 0's
final.plsda <- plsda(KO_SSC_7,points_3$Inoculum, ncomp = 10)

set.seed(30) # For reproducibility with this handbook, remove otherwise
perf.plsda <- perf(final.plsda, validation = 'Mfold', folds = 3, 
                   progressBar = FALSE,  # Set to TRUE to track progress
                   nrepeat = 10)     

list.keepX <- c(1:10,  seq(20, 100, 10))
list.keepX

tune.splsda_data <- tune.splsda(KO_SSC_7,points_3$Inoculum, ncomp = 4, validation = 'Mfold', 
                                folds = 5, dist = 'max.dist', 
                                test.keepX = list.keepX, nrepeat = 10)
ncomp <- tune.splsda_data$choice.ncomp$ncomp 
select.keepX <- tune.splsda_data$choice.keepX[1:ncomp]  

splsda.data <- splsda(KO_SSC_7,points_3$Inoculum, ncomp = ncomp, keepX = select.keepX) 
var.name.short <- colnames(KO_SSC_7)

empty_vector <- vector()

empty_vector_all <- data.frame(matrix(NA, ncol = 4))
colnames(empty_vector_all) <- c("KO", "contrib", "component", "Data")
empty_vector_all_2 <- empty_vector_all[-1,]

empty_vector_contrib <- data_frame()

for (comp in 1:ncomp) {
  list <- as.data.frame(splsda.data$loadings$X)
  nonredun <- row.names(list)[rowSums(list)!=0]
  nonredun_2 <- list[row.names(list) %in% nonredun,]
  nonredun_3 <- as.data.frame(cbind(nonredun, nonredun_2))
  comp_2 <- comp+1
  if (comp == 1) {
    PC1 <- nonredun_3[comp_2]
    PC1$KO <- nonredun_3$nonredun
  } else {
    PC1 <- nonredun_3[comp_2]
    PC1$KO <- row.names(PC1)
  }
  colnames(PC1) <- c(paste("comp",comp,sep=""),"KO")
  PC1_2 <- PC1[order(PC1[,paste("comp",comp,sep="")]),]
  PC1_3 <- PC1_2[PC1_2[,-1] !=0,]
  PC1_3 <- PC1_3[PC1_3[1] != 0,]
  PC1_3$extra <- abs(as.numeric(unlist(PC1_3[1])))
  
  empty_vector <- c(empty_vector,PC1_3$KO)
  
  table <- data.frame(PC1_3$KO)
  inter <- as.data.frame(do.call(cbind, PC1_3[3]))
  
  table$contrib <- unlist(inter[1])
  colnames(table) <- c("KO", "contrib")
  table$component <- paste("PC", comp, sep="")
  
  plot <- plotLoadings(splsda.data, comp = comp, method = 'mean', contrib = 'max', 
                       name.var = var.name.short)
  
  new_data <- data.frame(plot$X$GroupContrib)
  row.names(new_data) <- row.names(plot$X)
  empty_vector_contrib <- rbind(empty_vector_contrib,new_data)
  
  table$Data <- new_data$plot.X.GroupContrib[match(table$KO,row.names(new_data))]
  empty_vector_all_2 <- rbind(empty_vector_all_2, table)
}

#otu table
KO_SSC=read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Xanthomonadaceae_KO_no_dom.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t", row.names =1) #make the SampleID column into the row.names
colnames(samples_df)[5]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object

#Sample subsetting
samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

#At subset - Enterobacteriaceae - dom vs no-dom 
OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)
samples_sub = sample_data(samples_df_sub_2)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

bray_2 <- as.matrix(beta_isolate_KO)

str(samples_df_sub_2)
str(bray_2)

#Bind metadata with distance matrix
pcoa = cmdscale(bray_2, k=10, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z", "a", "b", "c", "d", "e", "f", "g") 
eig = pcoa$eig

empty_vector_all_2$Variance <- NA
empty_vector_all_2$Variance[empty_vector_all_2$component == "PC1"] <- eig[1]/sum(eig)
empty_vector_all_2$Variance[empty_vector_all_2$component == "PC2"] <- eig[2]/sum(eig)
empty_vector_all_2$Variance[empty_vector_all_2$component == "PC3"] <- eig[3]/sum(eig)
empty_vector_all_2$Variance[empty_vector_all_2$component == "PC4"] <- eig[4]/sum(eig)
empty_vector_all_2$Group <- "Xanthomonadaceae_SynCom"

write.table(empty_vector_all_2, paste(working_directory, "sPLS-DA/output/PLSDA_Xanthomonadaceae_no_dom_SynCom.tsv", sep=""), sep = "\t", quote = F, row.names = F, col.names = T)

###sPLS-DA - Xanthomonadaceae without dominators and only LjSC - Plants =====
#otu table
KO_SSC=read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Xanthomonadaceae_KO_no_dom.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t", row.names =1) #make the SampleID column into the row.names
colnames(samples_df)[5]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object

#Sample subsetting
samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

#At subset - Enterobacteriaceae - dom vs no-dom 
samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum == "LjSC")

OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)
samples_sub = sample_data(samples_df_sub_3)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

bray_2 <- as.matrix(beta_isolate_KO)

str(samples_df_sub_3)
str(bray_2)

#Bind metadata with distance matrix
pcoa = cmdscale(bray_2, k=10, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z", "a", "b", "c", "d", "e", "f", "g") 
eig = pcoa$eig
points_2 <- points[order(row.names(points)), ]
samples_df_sub_6 <- samples_df_sub_3[row.names(samples_df_sub_3) %in% row.names(points),]
samples_df_sub_7 <- samples_df_sub_6[order(row.names(samples_df_sub_6)), ]
points_3 <- cbind(points_2,samples_df_sub_7)
colnames(points_3) <- c("x", "y", "z","a", "b", "c", "d", "e", "f", "g",colnames(samples_df_sub_7))

# Principal Coordinates Analysis (PCoA)
pcoa_result <- cmdscale(bray_2, eig = TRUE, k = 2)  # k = number of dimensions

# Extract coordinates for plotting
pcoa_coords <- as.data.frame(pcoa_result$points)

# Adding metadata environmental data
pcoa_coords$Condition <- samples_df_sub_7$Condition[match(row.names(pcoa_coords), row.names(samples_df_sub_7))]

KO_SSC_3 <- KO_SSC[,colnames(KO_SSC) %in% row.names(pcoa_coords)]
KO_SSC_4 <- t(t(KO_SSC_3)/rowSums(t(KO_SSC_3)))

data_plot <- data.frame(matrix(NA, ncol = 5))
colnames(data_plot) <- c("KO", "contrib", "component", "Cluster", "Study")
data_plot_2 <- data_plot[-1,]

#mixOmics
KO_SSC_5 <- t(KO_SSC_4)
KO_SSC_7 <- na.omit(KO_SSC_5[match(row.names(points_3),row.names(KO_SSC_5)),])
points_4 <- points_3[row.names(points_3) %in% row.names(KO_SSC_7),]
points_3 <- points_4

#remove columns with 0's
final.plsda <- plsda(KO_SSC_7,points_3$Condition, ncomp = 10)

set.seed(30) # For reproducibility with this handbook, remove otherwise
perf.plsda <- perf(final.plsda, validation = 'Mfold', folds = 3, 
                   progressBar = FALSE,  # Set to TRUE to track progress
                   nrepeat = 10)     

list.keepX <- c(1:10,  seq(20, 100, 10))
list.keepX

tune.splsda_data <- tune.splsda(KO_SSC_7,points_3$Condition, ncomp = 4, validation = 'Mfold', 
                                folds = 5, dist = 'max.dist', 
                                test.keepX = list.keepX, nrepeat = 10)
ncomp <- tune.splsda_data$choice.ncomp$ncomp 
select.keepX <- tune.splsda_data$choice.keepX[1:ncomp]  

splsda.data <- splsda(KO_SSC_7,points_3$Condition, ncomp = ncomp, keepX = select.keepX) 
var.name.short <- colnames(KO_SSC_7)

empty_vector <- vector()

empty_vector_all <- data.frame(matrix(NA, ncol = 4))
colnames(empty_vector_all) <- c("KO", "contrib", "component", "Data")
empty_vector_all_2 <- empty_vector_all[-1,]

empty_vector_contrib <- data_frame()

for (comp in 1:ncomp) {
  list <- as.data.frame(splsda.data$loadings$X)
  nonredun <- row.names(list)[rowSums(list)!=0]
  nonredun_2 <- list[row.names(list) %in% nonredun,]
  nonredun_3 <- as.data.frame(cbind(nonredun, nonredun_2))
  comp_2 <- comp+1
  if (comp == 1) {
    PC1 <- nonredun_3[comp_2]
    PC1$KO <- nonredun_3$nonredun
  } else {
    PC1 <- nonredun_3[comp_2]
    PC1$KO <- row.names(PC1)
  }
  colnames(PC1) <- c(paste("comp",comp,sep=""),"KO")
  PC1_2 <- PC1[order(PC1[,paste("comp",comp,sep="")]),]
  PC1_3 <- PC1_2[PC1_2[,-1] !=0,]
  PC1_3 <- PC1_3[PC1_3[1] != 0,]
  PC1_3$extra <- abs(as.numeric(unlist(PC1_3[1])))
  
  empty_vector <- c(empty_vector,PC1_3$KO)
  
  table <- data.frame(PC1_3$KO)
  inter <- as.data.frame(do.call(cbind, PC1_3[3]))
  
  table$contrib <- unlist(inter[1])
  colnames(table) <- c("KO", "contrib")
  table$component <- paste("PC", comp, sep="")
  
  plot <- plotLoadings(splsda.data, comp = comp, method = 'mean', contrib = 'max', 
                       name.var = var.name.short)
  
  new_data <- data.frame(plot$X$GroupContrib)
  row.names(new_data) <- row.names(plot$X)
  empty_vector_contrib <- rbind(empty_vector_contrib,new_data)
  
  table$Data <- new_data$plot.X.GroupContrib[match(table$KO,row.names(new_data))]
  empty_vector_all_2 <- rbind(empty_vector_all_2, table)
}

#otu table
KO_SSC=read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Xanthomonadaceae_KO_no_dom.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t", row.names =1) #make the SampleID column into the row.names
colnames(samples_df)[5]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object

#Sample subsetting
samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

#At subset - Enterobacteriaceae - dom vs no-dom 
samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum == "LjSC")

OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)
samples_sub = sample_data(samples_df_sub_3)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

bray_2 <- as.matrix(beta_isolate_KO)

str(samples_df_sub_3)
str(bray_2)

#Bind metadata with distance matrix
pcoa = cmdscale(bray_2, k=10, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z", "a", "b", "c", "d", "e", "f", "g") 
eig = pcoa$eig

empty_vector_all_2$Variance <- NA
empty_vector_all_2$Variance[empty_vector_all_2$component == "PC1"] <- eig[1]/sum(eig)
empty_vector_all_2$Variance[empty_vector_all_2$component == "PC2"] <- eig[2]/sum(eig)
empty_vector_all_2$Group <- "Xanthomonadaceae_Plant"

write.table(empty_vector_all_2, paste(working_directory, "sPLS-DA/output/PLSDA_Xanthomonadaceae_no_dom_Plant.tsv", sep=""), sep = "\t", quote = F, row.names = F, col.names = T)

###sPLS-DA - Pseudomonadaceae without dominators and only LjSC - Plants =====
#otu table
KO_SSC=read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Pseudomonadaceae_KO_no_dom.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t", row.names =1) #make the SampleID column into the row.names
colnames(samples_df)[5]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object

#Sample subsetting
samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

#At subset - Enterobacteriaceae - dom vs no-dom 
samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum == "LjSC")

OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)
samples_sub = sample_data(samples_df_sub_3)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

bray_2 <- as.matrix(beta_isolate_KO)

str(samples_df_sub_3)
str(bray_2)

#Bind metadata with distance matrix
pcoa = cmdscale(bray_2, k=10, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z", "a", "b", "c", "d", "e", "f", "g") 
eig = pcoa$eig
points_2 <- points[order(row.names(points)), ]
samples_df_sub_6 <- samples_df_sub_3[row.names(samples_df_sub_3) %in% row.names(points),]
samples_df_sub_7 <- samples_df_sub_6[order(row.names(samples_df_sub_6)), ]
points_3 <- cbind(points_2,samples_df_sub_7)
colnames(points_3) <- c("x", "y", "z","a", "b", "c", "d", "e", "f", "g",colnames(samples_df_sub_7))

# Principal Coordinates Analysis (PCoA)
pcoa_result <- cmdscale(bray_2, eig = TRUE, k = 2)  # k = number of dimensions

# Extract coordinates for plotting
pcoa_coords <- as.data.frame(pcoa_result$points)

# Adding metadata environmental data
pcoa_coords$Condition <- samples_df_sub_7$Condition[match(row.names(pcoa_coords), row.names(samples_df_sub_7))]

KO_SSC_3 <- KO_SSC[,colnames(KO_SSC) %in% row.names(pcoa_coords)]
KO_SSC_4 <- t(t(KO_SSC_3)/rowSums(t(KO_SSC_3)))

data_plot <- data.frame(matrix(NA, ncol = 5))
colnames(data_plot) <- c("KO", "contrib", "component", "Cluster", "Study")
data_plot_2 <- data_plot[-1,]

#mixOmics
KO_SSC_5 <- t(KO_SSC_4)
KO_SSC_7 <- na.omit(KO_SSC_5[match(row.names(points_3),row.names(KO_SSC_5)),])
points_4 <- points_3[row.names(points_3) %in% row.names(KO_SSC_7),]
points_3 <- points_4

#remove columns with 0's
final.plsda <- plsda(KO_SSC_7,points_3$Condition, ncomp = 10)

set.seed(30) # For reproducibility with this handbook, remove otherwise
perf.plsda <- perf(final.plsda, validation = 'Mfold', folds = 3, 
                   progressBar = FALSE,  # Set to TRUE to track progress
                   nrepeat = 10)     

list.keepX <- c(1:10,  seq(20, 100, 10))
list.keepX

tune.splsda_data <- tune.splsda(KO_SSC_7,points_3$Condition, ncomp = 4, validation = 'Mfold', 
                                folds = 5, dist = 'max.dist', 
                                test.keepX = list.keepX, nrepeat = 10)
ncomp <- tune.splsda_data$choice.ncomp$ncomp 
select.keepX <- tune.splsda_data$choice.keepX[1:ncomp]  

splsda.data <- splsda(KO_SSC_7,points_3$Condition, ncomp = ncomp, keepX = select.keepX) 
var.name.short <- colnames(KO_SSC_7)

empty_vector <- vector()

empty_vector_all <- data.frame(matrix(NA, ncol = 4))
colnames(empty_vector_all) <- c("KO", "contrib", "component", "Data")
empty_vector_all_2 <- empty_vector_all[-1,]

empty_vector_contrib <- data_frame()

for (comp in 1:ncomp) {
  list <- as.data.frame(splsda.data$loadings$X)
  nonredun <- row.names(list)[rowSums(list)!=0]
  nonredun_2 <- list[row.names(list) %in% nonredun,]
  nonredun_3 <- as.data.frame(cbind(nonredun, nonredun_2))
  comp_2 <- comp+1
  if (comp == 1) {
    PC1 <- nonredun_3[comp_2]
    PC1$KO <- nonredun_3$nonredun
  } else {
    PC1 <- nonredun_3[comp_2]
    PC1$KO <- row.names(PC1)
  }
  colnames(PC1) <- c(paste("comp",comp,sep=""),"KO")
  PC1_2 <- PC1[order(PC1[,paste("comp",comp,sep="")]),]
  PC1_3 <- PC1_2[PC1_2[,-1] !=0,]
  PC1_3 <- PC1_3[PC1_3[1] != 0,]
  PC1_3$extra <- abs(as.numeric(unlist(PC1_3[1])))
  
  empty_vector <- c(empty_vector,PC1_3$KO)
  
  table <- data.frame(PC1_3$KO)
  inter <- as.data.frame(do.call(cbind, PC1_3[3]))
  
  table$contrib <- unlist(inter[1])
  colnames(table) <- c("KO", "contrib")
  table$component <- paste("PC", comp, sep="")
  
  plot <- plotLoadings(splsda.data, comp = comp, method = 'mean', contrib = 'max', 
                       name.var = var.name.short)
  
  new_data <- data.frame(plot$X$GroupContrib)
  row.names(new_data) <- row.names(plot$X)
  empty_vector_contrib <- rbind(empty_vector_contrib,new_data)
  
  table$Data <- new_data$plot.X.GroupContrib[match(table$KO,row.names(new_data))]
  empty_vector_all_2 <- rbind(empty_vector_all_2, table)
  
}

#otu table
KO_SSC=read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Pseudomonadaceae_KO_no_dom.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t", row.names =1) #make the SampleID column into the row.names
colnames(samples_df)[5]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object

#Sample subsetting
samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

#At subset - Enterobacteriaceae - dom vs no-dom 
samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum == "LjSC")

OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)
samples_sub = sample_data(samples_df_sub_3)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

bray_2 <- as.matrix(beta_isolate_KO)

str(samples_df_sub_3)
str(bray_2)

#Bind metadata with distance matrix
pcoa = cmdscale(bray_2, k=10, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z", "a", "b", "c", "d", "e", "f", "g") 
eig = pcoa$eig

empty_vector_all_2$Variance <- NA
empty_vector_all_2$Variance[empty_vector_all_2$component == "PC1"] <- eig[1]/sum(eig)
empty_vector_all_2$Variance[empty_vector_all_2$component == "PC2"] <- eig[2]/sum(eig)
empty_vector_all_2$Group <- "Pseudomonadaceae_Plant"

write.table(empty_vector_all_2, paste(working_directory, "sPLS-DA/output/PLSDA_Pseudomonadaceae_no_dom_Plant.tsv", sep=""), sep = "\t", quote = F, row.names = F, col.names = T)

###sPLS-DA - Caulobacteraceae without dominators and only AtSC - Plants =====

#otu table
KO_SSC=read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Caulobacteraceae_KO_no_dom.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t", row.names =1) #make the SampleID column into the row.names
colnames(samples_df)[5]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object

#Sample subsetting
samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

#At subset - Enterobacteriaceae - dom vs no-dom 
samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum == "AtSC")

OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)
samples_sub = sample_data(samples_df_sub_3)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

bray_2 <- as.matrix(beta_isolate_KO)

str(samples_df_sub_3)
str(bray_2)

#Bind metadata with distance matrix
pcoa = cmdscale(bray_2, k=10, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z", "a", "b", "c", "d", "e", "f", "g") 
eig = pcoa$eig
points_2 <- points[order(row.names(points)), ]
samples_df_sub_6 <- samples_df_sub_3[row.names(samples_df_sub_3) %in% row.names(points),]
samples_df_sub_7 <- samples_df_sub_6[order(row.names(samples_df_sub_6)), ]
points_3 <- cbind(points_2,samples_df_sub_7)
colnames(points_3) <- c("x", "y", "z","a", "b", "c", "d", "e", "f", "g",colnames(samples_df_sub_7))


# Principal Coordinates Analysis (PCoA)
pcoa_result <- cmdscale(bray_2, eig = TRUE, k = 2)  # k = number of dimensions

# Extract coordinates for plotting
pcoa_coords <- as.data.frame(pcoa_result$points)

# Adding metadata environmental data
pcoa_coords$Condition <- samples_df_sub_7$Condition[match(row.names(pcoa_coords), row.names(samples_df_sub_7))]

KO_SSC_3 <- KO_SSC[,colnames(KO_SSC) %in% row.names(pcoa_coords)]
KO_SSC_4 <- t(t(KO_SSC_3)/rowSums(t(KO_SSC_3)))

data_plot <- data.frame(matrix(NA, ncol = 5))
colnames(data_plot) <- c("KO", "contrib", "component", "Cluster", "Study")
data_plot_2 <- data_plot[-1,]

#mixOmics
KO_SSC_5 <- t(KO_SSC_4)
KO_SSC_7 <- na.omit(KO_SSC_5[match(row.names(points_3),row.names(KO_SSC_5)),])
points_4 <- points_3[row.names(points_3) %in% row.names(KO_SSC_7),]
points_3 <- points_4

#remove columns with 0's
final.plsda <- plsda(KO_SSC_7,points_3$Condition, ncomp = 10)

set.seed(30) # For reproducibility with this handbook, remove otherwise
perf.plsda <- perf(final.plsda, validation = 'Mfold', folds = 3, 
                   progressBar = FALSE,  # Set to TRUE to track progress
                   nrepeat = 10)     

list.keepX <- c(1:10,  seq(20, 100, 10))
list.keepX

tune.splsda_data <- tune.splsda(KO_SSC_7,points_3$Condition, ncomp = 4, validation = 'Mfold', 
                                folds = 5, dist = 'max.dist', 
                                test.keepX = list.keepX, nrepeat = 10)
ncomp <- tune.splsda_data$choice.ncomp$ncomp 
select.keepX <- tune.splsda_data$choice.keepX[1:ncomp]  

splsda.data <- splsda(KO_SSC_7,points_3$Condition, ncomp = ncomp, keepX = select.keepX) 
var.name.short <- colnames(KO_SSC_7)

empty_vector <- vector()

empty_vector_all <- data.frame(matrix(NA, ncol = 4))
colnames(empty_vector_all) <- c("KO", "contrib", "component", "Data")
empty_vector_all_2 <- empty_vector_all[-1,]

empty_vector_contrib <- data_frame()

for (comp in 1:ncomp) {
  list <- as.data.frame(splsda.data$loadings$X)
  nonredun <- row.names(list)[rowSums(list)!=0]
  nonredun_2 <- list[row.names(list) %in% nonredun,]
  nonredun_3 <- as.data.frame(cbind(nonredun, nonredun_2))
  comp_2 <- comp+1
  if (comp == 1) {
    PC1 <- nonredun_3[comp_2]
    PC1$KO <- nonredun_3$nonredun
  } else {
    PC1 <- nonredun_3[comp_2]
    PC1$KO <- row.names(PC1)
  }
  colnames(PC1) <- c(paste("comp",comp,sep=""),"KO")
  PC1_2 <- PC1[order(PC1[,paste("comp",comp,sep="")]),]
  PC1_3 <- PC1_2[PC1_2[,-1] !=0,]
  PC1_3 <- PC1_3[PC1_3[1] != 0,]
  PC1_3$extra <- abs(as.numeric(unlist(PC1_3[1])))
  
  empty_vector <- c(empty_vector,PC1_3$KO)
  
  table <- data.frame(PC1_3$KO)
  inter <- as.data.frame(do.call(cbind, PC1_3[3]))
  
  table$contrib <- unlist(inter[1])
  colnames(table) <- c("KO", "contrib")
  table$component <- paste("PC", comp, sep="")
  
  plot <- plotLoadings(splsda.data, comp = comp, method = 'mean', contrib = 'max', 
                       name.var = var.name.short)
  
  new_data <- data.frame(plot$X$GroupContrib)
  row.names(new_data) <- row.names(plot$X)
  empty_vector_contrib <- rbind(empty_vector_contrib,new_data)
  
  table$Data <- new_data$plot.X.GroupContrib[match(table$KO,row.names(new_data))]
  empty_vector_all_2 <- rbind(empty_vector_all_2, table)
  
}

#otu table
KO_SSC=read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Caulobacteraceae_KO_no_dom.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t", row.names =1) #make the SampleID column into the row.names
colnames(samples_df)[5]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object

#Sample subsetting
samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

#At subset - Enterobacteriaceae - dom vs no-dom 
samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum == "AtSC")

OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)
samples_sub = sample_data(samples_df_sub_3)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

bray_2 <- as.matrix(beta_isolate_KO)

str(samples_df_sub_3)
str(bray_2)

#Bind metadata with distance matrix
pcoa = cmdscale(bray_2, k=10, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z", "a", "b", "c", "d", "e", "f", "g") 
eig = pcoa$eig

empty_vector_all_2$Variance <- NA
empty_vector_all_2$Variance[empty_vector_all_2$component == "PC1"] <- eig[1]/sum(eig)
empty_vector_all_2$Group <- "Caulobacteraceae_Plant"

write.table(empty_vector_all_2, paste(working_directory, "sPLS-DA/isolate_subset_data/PLSDA_Caulobacteraceae_no_dom_Plant.tsv", sep=""), sep = "\t", quote = F, row.names = F, col.names = T)

###sPLS-DA - Burkholderiaceae without dominators and only AtSC - Plants =====
#otu table
KO_SSC=read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Burkholderiaceae_KO_no_dom.tsv", sep =""), header=TRUE,sep="\t", row.names = 1)

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t", row.names =1) #make the SampleID column into the row.names
colnames(samples_df)[5]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object

#Sample subsetting
samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

#At subset - Enterobacteriaceae - dom vs no-dom 
samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum == "AtSC")

OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)
samples_sub = sample_data(samples_df_sub_3)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

bray_2 <- as.matrix(beta_isolate_KO)

str(samples_df_sub_3)
str(bray_2)

#Bind metadata with distance matrix
pcoa = cmdscale(bray_2, k=10, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z", "a", "b", "c", "d", "e", "f", "g") 
eig = pcoa$eig
points_2 <- points[order(row.names(points)), ]
samples_df_sub_6 <- samples_df_sub_3[row.names(samples_df_sub_3) %in% row.names(points),]
samples_df_sub_7 <- samples_df_sub_6[order(row.names(samples_df_sub_6)), ]
points_3 <- cbind(points_2,samples_df_sub_7)
colnames(points_3) <- c("x", "y", "z","a", "b", "c", "d", "e", "f", "g",colnames(samples_df_sub_7))


# Principal Coordinates Analysis (PCoA)
pcoa_result <- cmdscale(bray_2, eig = TRUE, k = 2)  # k = number of dimensions

# Extract coordinates for plotting
pcoa_coords <- as.data.frame(pcoa_result$points)

# Adding metadata environmental data
pcoa_coords$Condition <- samples_df_sub_7$Condition[match(row.names(pcoa_coords), row.names(samples_df_sub_7))]

KO_SSC_3 <- KO_SSC[,colnames(KO_SSC) %in% row.names(pcoa_coords)]
KO_SSC_4 <- t(t(KO_SSC_3)/rowSums(t(KO_SSC_3)))

data_plot <- data.frame(matrix(NA, ncol = 5))
colnames(data_plot) <- c("KO", "contrib", "component", "Cluster", "Study")
data_plot_2 <- data_plot[-1,]

#mixOmics
KO_SSC_5 <- t(KO_SSC_4)
KO_SSC_7 <- na.omit(KO_SSC_5[match(row.names(points_3),row.names(KO_SSC_5)),])
points_4 <- points_3[row.names(points_3) %in% row.names(KO_SSC_7),]
points_3 <- points_4

#remove columns with 0's
final.plsda <- plsda(KO_SSC_7,points_3$Condition, ncomp = 10)

set.seed(30) # For reproducibility with this handbook, remove otherwise
perf.plsda <- perf(final.plsda, validation = 'Mfold', folds = 3, 
                   progressBar = FALSE,  # Set to TRUE to track progress
                   nrepeat = 10)     

list.keepX <- c(1:10,  seq(20, 100, 10))
list.keepX

tune.splsda_data <- tune.splsda(KO_SSC_7,points_3$Condition, ncomp = 4, validation = 'Mfold', 
                                folds = 5, dist = 'max.dist', 
                                test.keepX = list.keepX, nrepeat = 10)
ncomp <- tune.splsda_data$choice.ncomp$ncomp 
select.keepX <- tune.splsda_data$choice.keepX[1:ncomp]  

splsda.data <- splsda(KO_SSC_7,points_3$Condition, ncomp = ncomp, keepX = select.keepX) 
var.name.short <- colnames(KO_SSC_7)

empty_vector <- vector()

empty_vector_all <- data.frame(matrix(NA, ncol = 4))
colnames(empty_vector_all) <- c("KO", "contrib", "component", "Data")
empty_vector_all_2 <- empty_vector_all[-1,]

empty_vector_contrib <- data_frame()

for (comp in 1:ncomp) {
  list <- as.data.frame(splsda.data$loadings$X)
  nonredun <- row.names(list)[rowSums(list)!=0]
  nonredun_2 <- list[row.names(list) %in% nonredun,]
  nonredun_3 <- as.data.frame(cbind(nonredun, nonredun_2))
  comp_2 <- comp+1
  if (comp == 1) {
    PC1 <- nonredun_3[comp_2]
    PC1$KO <- nonredun_3$nonredun
  } else {
    PC1 <- nonredun_3[comp_2]
    PC1$KO <- row.names(PC1)
  }
  colnames(PC1) <- c(paste("comp",comp,sep=""),"KO")
  PC1_2 <- PC1[order(PC1[,paste("comp",comp,sep="")]),]
  PC1_3 <- PC1_2[PC1_2[,-1] !=0,]
  PC1_3 <- PC1_3[PC1_3[1] != 0,]
  PC1_3$extra <- abs(as.numeric(unlist(PC1_3[1])))
  
  empty_vector <- c(empty_vector,PC1_3$KO)
  
  table <- data.frame(PC1_3$KO)
  inter <- as.data.frame(do.call(cbind, PC1_3[3]))
  
  table$contrib <- unlist(inter[1])
  colnames(table) <- c("KO", "contrib")
  table$component <- paste("PC", comp, sep="")
  
  plot <- plotLoadings(splsda.data, comp = comp, method = 'mean', contrib = 'max', 
                       name.var = var.name.short)
  
  new_data <- data.frame(plot$X$GroupContrib)
  row.names(new_data) <- row.names(plot$X)
  empty_vector_contrib <- rbind(empty_vector_contrib,new_data)
  
  table$Data <- new_data$plot.X.GroupContrib[match(table$KO,row.names(new_data))]
  empty_vector_all_2 <- rbind(empty_vector_all_2, table)
  
}

#otu table
KO_SSC=read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Burkholderiaceae_KO_no_dom.tsv", sep =""), header=TRUE,sep="\t", row.names = 1)

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t", row.names =1) #make the SampleID column into the row.names
colnames(samples_df)[5]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object

#Sample subsetting
samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

#At subset - Enterobacteriaceae - dom vs no-dom 
samples_df_sub_3 <- subset(samples_df_sub_2, samples_df_sub_2$Inoculum == "AtSC")

OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)
samples_sub = sample_data(samples_df_sub_3)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

bray_2 <- as.matrix(beta_isolate_KO)

str(samples_df_sub_3)
str(bray_2)

#Bind metadata with distance matrix
pcoa = cmdscale(bray_2, k=10, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z", "a", "b", "c", "d", "e", "f", "g") 
eig = pcoa$eig

empty_vector_all_2$Variance <- NA
empty_vector_all_2$Variance[empty_vector_all_2$component == "PC1"] <- eig[1]/sum(eig)
empty_vector_all_2$Variance[empty_vector_all_2$component == "PC2"] <- eig[2]/sum(eig)
empty_vector_all_2$Group <- "Burkholderiaceae_Plant"

write.table(empty_vector_all_2, paste(working_directory, "sPLS-DA/output/PLSDA_Burkholderiaceae_no_dom_Plant.tsv", sep=""), sep = "\t", quote = F, row.names = F, col.names = T)

############################ Niche replacement files ###############################
###Creating Data_with_synthetic_input.tsv =====

#Since only one input samples was sequenced for each SynCom, we will subdivide these into four samples by random subsampling at a lower sequencing depth
otu_table <- read.table(paste(working_directory,"LjSC_Family_drop_out_experiment/original_data_norm.tsv", sep = ""), sep = "\t", row.names =1, header =T)

Input_data <- otu_table[,grep("INPUT", colnames(otu_table))]
OTU = otu_table(as.matrix(Input_data),taxa_are_rows = TRUE)
phylo_sub = phyloseq(OTU)

racur_data=rarecurve(x =floor(t(phylo_sub@.Data)), step = 250, label = F, tidy = T)
racur_data$Site <- gsub("LDT1_INPUT", "Burkholderiaceae drop out", racur_data$Site)
racur_data$Site <- gsub("LDT2_INPUT", "Caulobacteraceae drop out", racur_data$Site)
racur_data$Site <- gsub("LDT3_INPUT", "Pseudomonadaceae drop out", racur_data$Site)
racur_data$Site <- gsub("LDT4_INPUT", "Rhizobiaceae drop out", racur_data$Site)
racur_data$Site <- gsub("LDT5_INPUT", "All other families drop out", racur_data$Site)
racur_data$Site <- gsub("LDT6_INPUT", "Full LjSC", racur_data$Site)

values <- c()

for (group in unique(racur_data$Site)){
  racur_sub <- racur_data[racur_data$Site == paste(group),]
  number <- tail(racur_sub$Species, n=1) * 0.95
  depth <- racur_sub$Sample[racur_sub$Species > number][1]
  
  values <- c(values, depth)
}

#Selected rarefaction depth
rarefaction_depth <- sum(values)/length(values)

#creating subsets
otu_table_no_input <- otu_table[,!grepl("INPUT", colnames(otu_table))]
Input_data <- otu_table[,grep("INPUT", colnames(otu_table))]

OTU = otu_table(as.matrix(Input_data),taxa_are_rows = TRUE)
phylo_sub = phyloseq(OTU)
phylo_sub_2 <- rarefy_even_depth(phylo_sub, sample.size = rarefaction_depth, rngseed = 1)
phylo_sub_3 <- rarefy_even_depth(phylo_sub, sample.size = rarefaction_depth, rngseed = 2)
phylo_sub_4 <- rarefy_even_depth(phylo_sub, sample.size = rarefaction_depth, rngseed = 3)
phylo_sub_5 <- rarefy_even_depth(phylo_sub, sample.size = rarefaction_depth, rngseed = 4)

otu_2 <- phylo_sub_2@.Data
otu_3 <- phylo_sub_3@.Data
otu_4 <- phylo_sub_4@.Data
otu_5 <- phylo_sub_5@.Data

colnames(otu_2) <- gsub("INPUT","INPUT_B1", colnames(otu_2))
colnames(otu_3) <- gsub("INPUT","INPUT_B2", colnames(otu_3))
colnames(otu_4) <- gsub("INPUT","INPUT_B3", colnames(otu_4))
colnames(otu_5) <- gsub("INPUT","INPUT_B4", colnames(otu_5))

reset_names_2 <- row.names(otu_table_no_input)[!row.names(otu_table_no_input) %in% row.names(otu_2) ]
for (group in reset_names_2){
  new <- data.frame(t(data.frame(c(0,0,0,0,0,0))))
  row.names(new) <- paste(group)
  colnames(new) <- colnames(otu_2)
  otu_2 <- rbind(otu_2, new)
}

reset_names_3 <- row.names(otu_table_no_input)[!row.names(otu_table_no_input) %in% row.names(otu_3) ]
for (group in reset_names_3){
  new <- data.frame(t(data.frame(c(0,0,0,0,0,0))))
  row.names(new) <- paste(group)
  colnames(new) <- colnames(otu_3)
  otu_3 <- rbind(otu_3, new)
}

reset_names_4 <- row.names(otu_table_no_input)[!row.names(otu_table_no_input) %in% row.names(otu_4) ]
for (group in reset_names_4){
  new <- data.frame(t(data.frame(c(0,0,0,0,0,0))))
  row.names(new) <- paste(group)
  colnames(new) <- colnames(otu_4)
  otu_4 <- rbind(otu_4, new)
}

reset_names_5 <- row.names(otu_table_no_input)[!row.names(otu_table_no_input) %in% row.names(otu_5) ]
for (group in reset_names_5){
  new <- data.frame(t(data.frame(c(0,0,0,0,0,0))))
  row.names(new) <- paste(group)
  colnames(new) <- colnames(otu_5)
  otu_5 <- rbind(otu_5, new)
}

otu_2_sub <- otu_2[row.names(otu_table_no_input),]
otu_3_sub <- otu_3[row.names(otu_table_no_input),]
otu_4_sub <- otu_4[row.names(otu_table_no_input),]
otu_5_sub <- otu_5[row.names(otu_table_no_input),]

new_table <- cbind(otu_table_no_input,otu_2, otu_3, otu_4, otu_5)

write.table(new_table,paste(working_directory, "LjSC_Family_drop_out_experiment/Data_with_synthetic_input.tsv", sep =""), sep = "\t", quote =F, col.names = T, row.names =T)

###Creating Figure_5d_validation_in_vivo_Fam_drop_out_with_nod.txt =====
pathway_selection <- c("ABC transporters", "Quorum sensing", "Two-component system", "Purine metabolism", "Oxidative phosphorylation", "Pentose and glucuronate interconversions", "Porphyrin metabolism", "Galactose metabolism", "Cell cycle - Caulobacter", "Exopolysaccharide biosynthesis")

#Insilico drop-out - recalculating pathways
Families <- c("Burkholderiaceae", "Caulobacteraceae", "Pseudomonadaceae", "other", "all")

#Invivo drop-out - recalculating pathways
top <- read.table(paste(working_directory, "Annotations/pathway_top.txt", sep = ""), header=F, sep="\t")
KO_to_pathway <- read.table(paste(working_directory, "Annotations/KO_to_pathway.txt", sep = ""), header=T, sep="\t")
KO_to_pathway$V3 <- top$V2[match(KO_to_pathway$V2, top$V1)]

KO_to_pathway_2 <- read.table(paste(working_directory, "Annotations/KO_to_pathway_unannotated_2.txt", sep = ""), header=F, sep="\t")
colnames(KO_to_pathway_2) <- c("KO","new_category")

for (KO in KO_to_pathway_2$KO){
  KO_to_pathway$V3[KO_to_pathway$V1 == paste(KO)] <- KO_to_pathway_2$new_category[KO_to_pathway_2$KO == paste(KO)]
}

KO_table = read.table(paste(working_directory, "KO_genome/KO_LjSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
KO_table_2 <- t(KO_table)

#Significant KOs from SSC DEseq2 analysis
input_table <- read.table(paste(working_directory, "DESeq2/Sig_KO_all.txt", sep = ""), header=T, sep="\t")
input_table_2 <- input_table$KO[input_table$Plant == "Lotus" & input_table$SynCom == "LjSC"]

SynComs <- c("LDT1","LDT2","LDT3","LDT4","LDT5","LDT6")
hop_4 <- data.frame()

for (path in pathway_selection){
  KOs <- as.vector(na.omit(KO_to_pathway$V1[KO_to_pathway$V3 == paste(path)]))
  KOs_2 <- KOs[KOs %in% input_table_2]
  if(length(KOs_2) > 0){
    hop_2 <- data.frame()
    if(length(KOs_2) > 0 ){
      for (KO in KOs_2){
        for (syncom in SynComs){
          norm_SSC =read.table(paste(working_directory, "LjSC_Family_drop_out_experiment/Data_with_synthetic_input.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
          norm_SSC_2 <- norm_SSC[,grep("ROOT", colnames(norm_SSC))]
          norm_SSC_3 <- norm_SSC_2[,grep(paste(syncom), colnames(norm_SSC_2))]
          norm_SSC_input <- norm_SSC[,grep("INPUT", colnames(norm_SSC))]
          
          KO_table_3 <- data.frame(KO_table_2[,colnames(KO_table_2) == paste(KO)])
          colnames(KO_table_3) <- "KO"
          KO_table_4 <- KO_table_3[row.names(KO_table_3) %in% row.names(norm_SSC_3),]
          names(KO_table_4) <- row.names(KO_table_3)[row.names(KO_table_3) %in% row.names(norm_SSC_3)]
          
          norm_SSC_5 <- t(t(norm_SSC_3)/rowSums(t(norm_SSC_3)))
          
          No <- names(KO_table_4)[KO_table_4 == 0]
          Yes <- names(KO_table_4)[KO_table_4 != 0]
          
          norm_SSC_yes <- norm_SSC_5[row.names(norm_SSC_5) %in% Yes,]
          norm_SSC_no <- norm_SSC_5[row.names(norm_SSC_5) %in% No,]
          
          if (length(Yes) > 1){
            Yes_sum <- sum(colSums(norm_SSC_yes))/length(colSums(norm_SSC_yes))
          } else if (length(Yes) == 1){
            Yes_sum <- sum(norm_SSC_yes)/length(norm_SSC_yes)
          } else {
            Yes_sum <- 0
          }
          
          if (length(No) > 1){
            No_sum <- sum(colSums(norm_SSC_no))/length(colSums(norm_SSC_no))
          } else if (length(No) == 1){
            No_sum <- sum(norm_SSC_no)/length(norm_SSC_no)
          } else {
            No_sum <- 0
          }
          
          norm_SSC_input_2 <- t(t(norm_SSC_input)/rowSums(t(norm_SSC_input)))
          norm_SSC_input_3 <- norm_SSC_input_2[,grep(paste(syncom), colnames(norm_SSC_input_2))]
          norm_SSC_input_yes <- norm_SSC_input_3[row.names(norm_SSC_input_3) %in% Yes,]
          norm_SSC_input_no <- norm_SSC_input_3[row.names(norm_SSC_input_3) %in% No,]
          
          if (length(Yes) > 1){
            Yes_sum_input <- sum(colSums(norm_SSC_input_yes))/length(colSums(norm_SSC_input_yes))
          } else if (length(Yes) == 1){
            Yes_sum_input <- sum(norm_SSC_input_yes)/length(norm_SSC_input_yes)
          } else {
            Yes_sum_input <- 0
          }
          
          if (length(No) > 1){
            No_sum_input <- sum(colSums(norm_SSC_input_no))/length(colSums(norm_SSC_input_no))
          } else if (length(No) == 1){
            No_sum_input <- sum(norm_SSC_input_no)/length(norm_SSC_input_no)
          } else {
            No_sum_input <- 0
          }
          
          hop <- t(data.frame(c(paste(KO), Yes_sum, No_sum, Yes_sum_input, No_sum_input, length(Yes), length(No),paste(syncom))))
          hop_2 <- rbind(hop_2, hop)
        }
      }
      
      hop_2$V2 <- as.numeric(hop_2$V2)
      hop_2$V3 <- as.numeric(hop_2$V3)
      hop_2$V4 <- as.numeric(hop_2$V4)
      hop_2$V5 <- as.numeric(hop_2$V5)
      hop_2$V6 <- as.numeric(hop_2$V6)
      hop_2$V7 <- as.numeric(hop_2$V7)
      
      int_value <- min(hop_2$V4[hop_2$V4 != 0])
      int_value_2 <- min(hop_2$V5[hop_2$V5 != 0])
      
      hop_2$V4[hop_2$V4 == 0] <- int_value
      hop_2$V5[hop_2$V5 == 0] <- int_value_2
      
      hop_2$V9 <- (hop_2$V2/hop_2$V4)
      hop_2$V10 <- (hop_2$V3/hop_2$V5)
      
      vec <- hop_2$V6/(hop_2$V6 + hop_2$V7)
      vec[is.nan(vec)] <- 0
      
      value <- sum(vec)/length(hop_2$V6)
      
      for (syncom in SynComs){
        hop_2_6 <- hop_2[hop_2$V8 == paste(syncom),]
        value_path_yes <- sum(as.numeric(hop_2_6$V9))/length(hop_2_6$V9)
        value_path_no <- sum(as.numeric(hop_2_6$V10))/length(hop_2_6$V10)
        value_RA <-  sum(as.numeric(hop_2_6$V2))/length(hop_2_6$V2)
        hop_3 <- t(data.frame(c(value_RA, value_path_yes, value_path_no, paste(syncom), paste(path), value)))
        hop_4 <- rbind(hop_4, hop_3)
      }
    }
  }
}

row.names(hop_4) <- NULL
colnames(hop_4) <- c("RA", "Present_fold_change", "Absent_fold_change", "Family", "Pathway", "Percentage")

write.table(hop_4, paste(working_directory, "LjSC_Family_drop_out_experiment/Figure_5d_validation_in_vivo_Fam_drop_out_with_nod.txt", sep = ""), quote = F, col.names = T, row.names = F, sep = "\t")

###Creating Figure_5f_validation_in_vivo_Fam_drop_out_no_nod.txt =====
pathway_selection <- c("Vitamin B6 metabolism","Transcriptional regulator","Secretion","Quorum sensing",
                       "Phenylalanine, tyrosine and tryptophan biosynthesis","Pantothenate and CoA biosynthesis",
                       "Oxidative phosphorylation","Membrane protein","Glycerophospholipid metabolism",
                       "Folate biosynthesis","Flagellar assembly","Exopolysaccharide biosynthesis","Cysteine and methionine metabolism",
                       "Arginine biosynthesis", "Methane metabolism", "Aminobenzoate degradation", "Ascorbate and aldarate metabolism")

#Insilico drop-out - recalculating pathways
Families <- c("Burkholderiaceae", "Caulobacteraceae", "Pseudomonadaceae", "Rhizobiaceae", "other", "all")

input_table <- read.table(paste(working_directory, "DESeq2/Sig_KO_all_no_nod_rhizo.txt", sep = ""), header=T, sep="\t")
input_table_2 <- input_table$KO[input_table$plant == "Lotus" & input_table$SynCom == "LjSC"]

#Invivo drop-out - recalculating pathways
top <- read.table(paste(working_directory, "Annotations/pathway_top.txt", sep = ""), header=F, sep="\t")
KO_to_pathway <- read.table(paste(working_directory, "Annotations/KO_to_pathway.txt", sep = ""), header=T, sep="\t")
KO_to_pathway$V3 <- top$V2[match(KO_to_pathway$V2, top$V1)]

KO_to_pathway_2 <- read.table(paste(working_directory, "Annotations/KO_to_pathway_unannotated_2.txt", sep = ""), header=F, sep="\t")
colnames(KO_to_pathway_2) <- c("KO","new_category")

for (KO in KO_to_pathway_2$KO){
  KO_to_pathway$V3[KO_to_pathway$V1 == paste(KO)] <- KO_to_pathway_2$new_category[KO_to_pathway_2$KO == paste(KO)]
}

KO_table = read.table(paste(working_directory, "KO_genome/KO_LjSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
KO_table_2 <- t(KO_table)

SynComs <- c("LDT1","LDT2","LDT3","LDT4","LDT5","LDT6")
hop_4 <- data.frame()

for (path in pathway_selection){
  KOs <- as.vector(na.omit(KO_to_pathway$V1[KO_to_pathway$V3 == paste(path)]))
  KOs_2 <- KOs[KOs %in% input_table_2]
  if(length(KOs_2) > 0){
    hop_2 <- data.frame()
    if(length(KOs_2) > 0 ){
      for (KO in KOs_2){
        for (syncom in SynComs){
          norm_SSC =read.table(paste(working_directory, "LjSC_Family_drop_out_experiment/Data_with_synthetic_input.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
          norm_SSC_1 <- norm_SSC[!row.names(norm_SSC) %in% c("LjNodule214", "LjRoot228", "LjRoot234"),]
          norm_SSC_2 <- norm_SSC_1[,grep("ROOT", colnames(norm_SSC_1))]
          norm_SSC_3 <- norm_SSC_2[,grep(paste(syncom), colnames(norm_SSC_2))]
          norm_SSC_input <- norm_SSC_1[,grep("INPUT", colnames(norm_SSC_1))]
          
          KO_table_3 <- data.frame(KO_table_2[,colnames(KO_table_2) == paste(KO)])
          colnames(KO_table_3) <- "KO"
          KO_table_4 <- KO_table_3[row.names(KO_table_3) %in% row.names(norm_SSC_3),]
          names(KO_table_4) <- row.names(KO_table_3)[row.names(KO_table_3) %in% row.names(norm_SSC_3)]
          
          norm_SSC_5 <- t(t(norm_SSC_3)/rowSums(t(norm_SSC_3)))
          
          No <- names(KO_table_4)[KO_table_4 == 0]
          Yes <- names(KO_table_4)[KO_table_4 != 0]
          
          norm_SSC_yes <- norm_SSC_5[row.names(norm_SSC_5) %in% Yes,]
          norm_SSC_no <- norm_SSC_5[row.names(norm_SSC_5) %in% No,]
          
          if (length(Yes) > 1){
            Yes_sum <- sum(colSums(norm_SSC_yes))/length(colSums(norm_SSC_yes))
          } else if (length(Yes) == 1){
            Yes_sum <- sum(norm_SSC_yes)/length(norm_SSC_yes)
          } else {
            Yes_sum <- 0
          }
          
          if (length(No) > 1){
            No_sum <- sum(colSums(norm_SSC_no))/length(colSums(norm_SSC_no))
          } else if (length(No) == 1){
            No_sum <- sum(norm_SSC_no)/length(norm_SSC_no)
          } else {
            No_sum <- 0
          }
          
          norm_SSC_input_2 <- t(t(norm_SSC_input)/rowSums(t(norm_SSC_input)))
          norm_SSC_input_3 <- norm_SSC_input_2[,grep(paste(syncom), colnames(norm_SSC_input_2))]
          norm_SSC_input_yes <- norm_SSC_input_3[row.names(norm_SSC_input_3) %in% Yes,]
          norm_SSC_input_no <- norm_SSC_input_3[row.names(norm_SSC_input_3) %in% No,]
          
          if (length(Yes) > 1){
            Yes_sum_input <- sum(colSums(norm_SSC_input_yes))/length(colSums(norm_SSC_input_yes))
          } else if (length(Yes) == 1){
            Yes_sum_input <- sum(norm_SSC_input_yes)/length(norm_SSC_input_yes)
          } else {
            Yes_sum_input <- 0
          }
          
          if (length(No) > 1){
            No_sum_input <- sum(colSums(norm_SSC_input_no))/length(colSums(norm_SSC_input_no))
          } else if (length(No) == 1){
            No_sum_input <- sum(norm_SSC_input_no)/length(norm_SSC_input_no)
          } else {
            No_sum_input <- 0
          }
          
          hop <- t(data.frame(c(paste(KO), Yes_sum, No_sum, Yes_sum_input, No_sum_input, length(Yes), length(No),paste(syncom))))
          hop_2 <- rbind(hop_2, hop)
        }
      }
      
      hop_2$V2 <- as.numeric(hop_2$V2)
      hop_2$V3 <- as.numeric(hop_2$V3)
      hop_2$V4 <- as.numeric(hop_2$V4)
      hop_2$V5 <- as.numeric(hop_2$V5)
      hop_2$V6 <- as.numeric(hop_2$V6)
      hop_2$V7 <- as.numeric(hop_2$V7)
      
      int_value <- min(hop_2$V4[hop_2$V4 != 0])
      int_value_2 <- min(hop_2$V5[hop_2$V5 != 0])
      
      hop_2$V4[hop_2$V4 == 0] <- int_value
      hop_2$V5[hop_2$V5 == 0] <- int_value_2
      
      hop_2$V9 <- (hop_2$V2/hop_2$V4)
      hop_2$V10 <- (hop_2$V3/hop_2$V5)
      
      vec <- hop_2$V6/(hop_2$V6 + hop_2$V7)
      vec[is.nan(vec)] <- 0
      
      value <- sum(vec)/length(hop_2$V6)
      
      for (syncom in SynComs){
        hop_2_6 <- hop_2[hop_2$V8 == paste(syncom),]
        value_path_yes <- sum(as.numeric(hop_2_6$V9))/length(hop_2_6$V9)
        value_path_no <- sum(as.numeric(hop_2_6$V10))/length(hop_2_6$V10)
        value_RA <-  sum(as.numeric(hop_2_6$V2))/length(hop_2_6$V2)
        hop_3 <- t(data.frame(c(value_RA, value_path_yes, value_path_no, paste(syncom), paste(path), value)))
        hop_4 <- rbind(hop_4, hop_3)
      }
    }
  }
}

row.names(hop_4) <- NULL
colnames(hop_4) <- c("RA", "Present_fold_change", "Absent_fold_change", "Family", "Pathway", "Percentage")

write.table(hop_4, paste(working_directory, "LjSC_Family_drop_out_experiment/Figure_5f_validation_in_vivo_Fam_drop_out_no_nod.txt", sep = ""), quote = F, col.names = T, row.names = F, sep = "\t")


#################### Files - KOs in natural soils/rhizospheres ####################
###852 and 266 KOs in natural rhizospheres - Merging data with Lopez et al. 2023 =====

ST <- read.table(paste(working_directory,"KO_genome/KO_SSC.tsv", sep = ""),header = TRUE)

res_ara <- readRDS(paste(working_directory,"Natural_soil/ANCOMBC_output/ANCOMBC2_ara.rds", sep = ""))
res_ara <- res_ara%>%dplyr::select(taxon, lfc_TreatmentSoil, p_TreatmentSoil)
colnames(res_ara) <- c("KO", "LFC", "padj")
res_ara$LFC <- -1*res_ara$LFC
res_ara_JJ <- readRDS(paste(working_directory,"Natural_soil/ANCOMBC_output/ANCOMBC2_ara_JJ.rds", sep = ""))
res_ara_JJ <- res_ara_JJ%>%dplyr::select(taxon, lfc_TreatmentSoil, p_TreatmentSoil)
colnames(res_ara_JJ) <- c("KO", "LFC", "padj")
res_ara_JJ$LFC <- -1*res_ara_JJ$LFC

res_ara$padj <- p.adjust(res_ara$padj, method = "fdr")
res_ara_JJ$padj <- p.adjust(res_ara_JJ$padj, method = "fdr")

ancombc2_araC <- res_ara
ancombc2_araJ <- res_ara_JJ

#Since the ANCOMBC analyses on ancombc2_araC led to no significant results, it was later discarded, only ancombc2_araJ was taken
ancombc2_araC_R <- filter(ancombc2_araC, padj<0.05&LFC>0)
ancombc2_araJ_R <- filter(ancombc2_araJ, padj<0.05&LFC>0)

### KO sets
stringent_KO <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/266KO.tsv", sep =""))
general_KO <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/852KO.tsv", sep =""))

Jose_KO <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/Lopez_et_al_2023_output.tsv", sep= ""))
Jose_KO_R <- filter(Jose_KO, niche_association == "Rhizosphere")

Jose_KO_R_2 <- unique(dplyr::select(Jose_KO_R, Orthogroup_Id, molecular_function))
Jose_KO_R_2$Lopez <- "YES"
Jose_852 <- right_join(Jose_KO_R_2, general_KO, by=c("Orthogroup_Id"="KO"))
Jose_266 <- inner_join(Jose_KO_R_2, stringent_KO, by=c("Orthogroup_Id"="KO"))
AraJJ_852 <- filter(general_KO, KO%in%ancombc2_araJ_R$KO)
AraJJ_852_2 <- dplyr::select(AraJJ_852, KO)
AraJJ_852_2$Sanchez <- "YES"
Jose_852_araJJ <- left_join(Jose_852, AraJJ_852_2, by=c("Orthogroup_Id"="KO"))
Jose_852_araJJ <- Jose_852_araJJ[,c(1,2,6,3,4,5)]
write.table(Jose_852_araJJ, paste(working_directory, "Natural_soil/852_sig_in_natural_rhizospheres.tsv", sep = ""),sep = "\t", row.names = FALSE)

### stringent
AraJJ_266 <- filter(stringent_KO, KO%in%ancombc2_araJ_R$KO)
AraJJ_266_2 <- dplyr::select(AraJJ_266, KO)
AraJJ_266_2$Sanchez <- "YES"
Jose_266_araJJ <- left_join(Jose_266, AraJJ_266_2, by=c("Orthogroup_Id"="KO"))
Jose_266_araJJ <- Jose_266_araJJ[,c(1,2,6,3,4,5)]
write.table(Jose_266_araJJ, paste(working_directory, "Natural_soil/266_sig_in_natural_rhizospheres.tsv",sep = ""), sep = "\t", row.names = FALSE)


###Generating KO lists in natural soil and rhizosphere metagenome samples =====
#Cucumber - Ofek Lalzar et al. (2014)

#Function to extract KOs from eggnog files
get_KO <- function(tab){
  tab_2 <- dplyr::select(tab, c("X.query", "KEGG_ko"))
  tab_2_unnested <- tab_2 %>%
    separate_rows(KEGG_ko, sep = ",")
  tab_2_unnested <- filter(tab_2_unnested, KEGG_ko != "-")
  tab_2_unnested$KEGG_ko <- gsub("ko:", "", tab_2_unnested$KEGG_ko)
  return(tab_2_unnested)
}

cuc_soil_1 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR908279.emapper.annotations", sep = ""), skip = 4)
cuc_soil_2 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR908281.emapper.annotations", sep = ""), skip = 4)                   
cuc_rhizo_1 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR908208.emapper.annotations", sep = ""),skip = 4)
cuc_rhizo_2 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR908211.emapper.annotations", sep = ""), skip = 4)
cuc_rhizo_3 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR908272.emapper.annotations", sep = ""), skip = 4)

cucS1 <- get_KO(cuc_soil_1)
cucS2 <- get_KO(cuc_soil_2)
cucR1 <- get_KO(cuc_rhizo_1)
cucR2 <- get_KO(cuc_rhizo_2)
cucR3 <- get_KO(cuc_rhizo_3)

write.table(rbind(cucS1, cucS2), paste(working_directory,"Natural_soil/KO_lists/KO_list_cucumber_soil.tsv",sep =""), row.names = F, quote = F)
write.table(rbind(cucR1, cucR2, cucR3), paste(working_directory,"Natural_soil/KO_lists/KO_list_cucumber_rhizosphere.tsv",sep =""), row.names = F, quote = F)

#Wheat - Ofek Lalzar et al. (2014)
wheat_soil_1 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR908290.emapper.annotations", sep = ""), skip = 4)
wheat_soil_2 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR908291.emapper.annotations", sep = ""), skip = 4)                   
wheat_rhizo_1 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR908273.emapper.annotations", sep = ""),skip = 4)
wheat_rhizo_2 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR908275.emapper.annotations", sep = ""), skip = 4)
wheat_rhizo_3 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR908276.emapper.annotations", sep = ""), skip = 4)

wheatS1 <- get_KO(wheat_soil_1)
wheatS2 <- get_KO(wheat_soil_2)
wheatR1 <- get_KO(wheat_rhizo_1)
wheatR2 <- get_KO(wheat_rhizo_2)
wheatR3 <- get_KO(wheat_rhizo_3)

write.table(rbind(wheatS1, wheatS2), paste(working_directory,"Natural_soil/KO_lists/KO_list_wheat_ofek_soil.tsv",sep =""), row.names = F, quote = F)
write.table(rbind(wheatR1, wheatR2, wheatR3), paste(working_directory,"Natural_soil/KO_lists/KO_list_wheat_ofek_rhizosphere.tsv",sep =""), row.names = F, quote = F)

# Arabidopsis - Stringlis et al. (2018)
ara_soil_1 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR6797242.emapper.annotations",sep =""), skip = 4)
ara_soil_2 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR6797243.emapper.annotations",sep =""), skip = 4)                   
ara_soil_3 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR6797244.emapper.annotations",sep =""), skip = 4) 
ara_rhizo_1 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR6797246.emapper.annotations",sep =""), skip = 4)
ara_rhizo_2 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR6797249.emapper.annotations",sep =""), skip = 4)
ara_rhizo_3 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR6797250.emapper.annotations",sep =""), skip = 4)

araS1 <- get_KO(ara_soil_1)
araS2 <- get_KO(ara_soil_2)
araS3 <- get_KO(ara_soil_3)
araR1 <- get_KO(ara_rhizo_1)
araR2 <- get_KO(ara_rhizo_2)
araR3 <- get_KO(ara_rhizo_3)

# to get the list of KOs for overlap
write.table(rbind(araS1, araS2, araS3), paste(working_directory,"Natural_soil/KO_lists/KO_list_stringlis_arabidopsis_soil.tsv",sep =""), row.names = F, quote = F)
write.table(rbind(araR1, araR2, araR3), paste(working_directory,"Natural_soil/KO_lists/KO_list_stringlis_arabidopsis_rhizosphere.tsv",sep =""), row.names = F, quote = F)

metadata <- read.delim(paste(working_directory,"Natural_soil/Metadata_metagenomes_GenusLevel.tsv", sep = ""), header =T, row.names =1, sep = "\t")
JJ_samples <- filter(metadata, origin=="Athal_Sanchez")

eggnog_files <- list.files(paste(working_directory,"Natural_soil/preliminary_files/eggnog_Arabidopsis_extra/", sep = ""), pattern = "annotation", full.names = TRUE)

combined_data <- do.call(rbind, lapply(eggnog_files, function(file) {
  read.delim(file, skip = 4)
}))
ko_table <- get_KO(combined_data)
soil_samples <- filter(JJ_samples, source=="Soil")
rhizo_samples <- filter(JJ_samples, source=="Root")
ara_RS <- ko_table
ara_RS$Contig <- gsub("(_[^_]+)$", "", ara_RS$X.query)
ara_RS$niche <- sub("_.*$", "", ara_RS$X.query)
ara_RS$niche <- gsub("-","_", ara_RS$niche)

soil_ara_KO <- ara_RS %>% filter(niche%in%soil_samples$id)
soil_ara_KO <- soil_ara_KO[,1:2]
rhizo_ara_KO <- ara_RS %>% filter(niche%in%rhizo_samples$id)
rhizo_ara_KO <- rhizo_ara_KO[,1:2]
write.table(soil_ara_KO, paste(working_directory,"Natural_soil/KO_lists/KO_list_juanjo_arabidopsis_soil.tsv",sep =""), row.names = F, quote = F)
write.table(rhizo_ara_KO, paste(working_directory,"Natural_soil/KO_lists/KO_list_juanjo_arabidopsis_rhizosphere.tsv",sep =""), row.names = F, quote = F)

###ANCOM-BC results - Natural Rhizosphere vs Soil - Arabidopsis datasets =====

#Data is either taken from Lopez et al. 2023, or from an ANCOMBC analysis on the Arabidopsis metagenome data in Reijerscamp soil that is shown here below

#Function to extract KOs from eggnog files
get_KO <- function(tab){
  tab_2 <- dplyr::select(tab, c("X.query", "KEGG_ko"))
  tab_2_unnested <- tab_2 %>%
    separate_rows(KEGG_ko, sep = ",")
  tab_2_unnested <- filter(tab_2_unnested, KEGG_ko != "-")
  tab_2_unnested$KEGG_ko <- gsub("ko:", "", tab_2_unnested$KEGG_ko)
  return(tab_2_unnested)
}

#Loading Arabidopsis data - Stringlis et al. (2018)
ara_soil_1 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR6797242.emapper.annotations",sep =""), skip = 4)
ara_soil_2 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR6797243.emapper.annotations",sep =""), skip = 4)                   
ara_soil_3 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR6797244.emapper.annotations",sep =""), skip = 4) 
ara_rhizo_1 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR6797246.emapper.annotations",sep =""), skip = 4)
ara_rhizo_2 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR6797249.emapper.annotations",sep =""), skip = 4)
ara_rhizo_3 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR6797250.emapper.annotations",sep =""), skip = 4)

araS1 <- get_KO(ara_soil_1)
araS2 <- get_KO(ara_soil_2)
araS3 <- get_KO(ara_soil_3)
araR1 <- get_KO(ara_rhizo_1)
araR2 <- get_KO(ara_rhizo_2)
araR3 <- get_KO(ara_rhizo_3)

metadata <- read.delim(paste(working_directory,"Natural_soil/Metadata_metagenomes_GenusLevel.tsv", sep = ""), header =T, row.names =1, sep = "\t")
Stringlis_samples <- filter(metadata, origin=="Athal_Stringlis")

ara_RS <- rbind(araS1, araS2, araS3, araR1, araR2, araR3)
ara_RS$Sample <- gsub("_.*", "", ara_RS$X.query)
colnames(Stringlis_samples)[colnames(Stringlis_samples) == "label"] <- "Sample"
ara_RS_cov <- left_join(ara_RS, Stringlis_samples)
ara_RS_cov <- ara_RS_cov %>% dplyr::select(-X.query, -origin, -method, -source)
ara_RS_cov_sum <- ara_RS_cov %>% dplyr::count(KEGG_ko, Sample)
ara_RS_cov_sum_2 <- ara_RS_cov_sum[!grepl("#", ara_RS_cov_sum$Sample),]
feature.table <- ara_RS_cov_sum_2 %>% pivot_wider(names_from = Sample, values_from = n, values_fill = 0)

colnames(feature.table) <- c("KEGG_ko", "S1", "S2", "S3", "R1", "R2", "R3") 
feature.table <- feature.table[, c("KEGG_ko", "R1", "R2", "R3", "S1", "S2", "S3")]
row_names <- feature.table$KEGG_ko
feature.table <- feature.table[-1]
rownames(feature.table) <- row_names
otu_table <- otu_table(feature.table, taxa_are_rows = TRUE)
taxonomy_table <- NULL
sample_metadata <- data.frame(
  Sample_ID = c("R1", "R2", "R3","S1", "S2", "S3"),
  Treatment = c("Rhizosphere", "Rhizosphere", "Rhizosphere",
                "Soil", "Soil", "Soil")
)
rownames(sample_metadata) <- sample_metadata$Sample_ID
sample_data_object <- sample_data(sample_metadata)
physeq <- phyloseq(otu_table,
                   taxonomy_table,
                   sample_data_object)
SSs_ancombc <- ancombc2(physeq, fix_formula = "Treatment")
results <- SSs_ancombc$res
saveRDS(results,paste(working_directory,"Natural_soil/ANCOMBC_output/ANCOMBC2_ara.rds", sep = ""))

### Arabidopsis_Juanjo_dataset
metadata <- read.delim(paste(working_directory,"Natural_soil/Metadata_metagenomes_GenusLevel.tsv", sep = ""), header =T, row.names =1, sep = "\t")
JJ_samples <- filter(metadata, origin=="Athal_Sanchez")

eggnog_files <- list.files(paste(working_directory,"Natural_soil/preliminary_files/eggnog_Arabidopsis_extra/", sep = ""), pattern = "annotation", full.names = TRUE)

combined_data <- do.call(rbind, lapply(eggnog_files, function(file) {
  read.delim(file, skip = 4)
}))
ko_table <- get_KO(combined_data)
soil_samples <- filter(JJ_samples, source=="Soil")
rhizo_samples <- filter(JJ_samples, source=="Root")
ara_RS <- ko_table
ara_RS$Contig <- gsub("(_[^_]+)$", "", ara_RS$X.query)
ara_RS$Sample <- sub("_.*$", "", ara_RS$X.query)
ara_RS$Sample <- gsub("-","_", ara_RS$Sample)

colnames(JJ_samples)[colnames(JJ_samples) == "label"] <- "Sample"
ara_RS_cov <- left_join(ara_RS, JJ_samples)
ara_RS_cov <- ara_RS_cov %>% dplyr::select(-X.query, -Contig, -origin, -method, -source)
ara_RS_cov_sum <- ara_RS_cov %>% dplyr::count(KEGG_ko, Sample)
ara_RS_cov_sum_2 <- ara_RS_cov_sum[!grepl("#", ara_RS_cov_sum$Sample),]
feature.table <- ara_RS_cov_sum_2 %>% pivot_wider(names_from = Sample, values_from = n, values_fill = 0)
feature.table <- feature.table[,c("KEGG_ko",row.names(JJ_samples))]
colnames(feature.table) <- c("KEGG_ko", "S1", "S2", "S3", "S4", "S5", "R1", "R2", "R3", 
                             "R4", "R5", "R6", "R7", "R8", "R9", "R10", "R11")
feature.table <- feature.table[, c("KEGG_ko", "R1", "R2", "R3", 
                                   "R4", "R5", "R6", "R7", "R8", "R9", "R10", "R11",
                                   "S1", "S2", "S3", "S4", "S5")]
row_names <- feature.table$KEGG_ko
feature.table <- feature.table[-1]
rownames(feature.table) <- row_names
otu_table <- otu_table(feature.table, taxa_are_rows = TRUE)
taxonomy_table <- NULL
sample_metadata <- data.frame(
  Sample_ID = c("R1", "R2", "R3", 
                "R4", "R5", "R6", "R7", "R8", "R9", "R10", "R11", "S1", "S2", "S3", "S4", "S5"),
  Treatment = c("Rhizosphere", "Rhizosphere", "Rhizosphere", "Rhizosphere", "Rhizosphere", "Rhizosphere",
                "Rhizosphere", "Rhizosphere", "Rhizosphere", "Rhizosphere", "Rhizosphere", 
                "Soil", "Soil", "Soil", "Soil", "Soil")
)
rownames(sample_metadata) <- sample_metadata$Sample_ID
sample_data_object <- sample_data(sample_metadata)
physeq <- phyloseq(otu_table,
                   taxonomy_table,
                   sample_data_object)
SSs_ancombc <- ancombc2(physeq, fix_formula = "Treatment")
results <- SSs_ancombc$res
saveRDS(results, paste(working_directory,"Natural_soil/ANCOMBC_output/ANCOMBC2_ara_JJ.rds", sep = ""))

############################# Shiny app files ######################################

#The shiny app allows the user to generate the same plots in S24-S29 for any KO from the dataset
#The generated files to make the Shiny app work were created in the following scripts

#Set output of these scripts
results_sa.dir <- paste(working_directory, "Shiny_app/", sep = "")

###Generation of ternary_KOs_av_med_all.txt - FC of isolates for each KO (split over 4 SynComs) - Original dataset =====
KO_table = read.table(paste(working_directory, "KO_genome/KO_SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
colnames(KO_table) <- gsub("X", "", colnames(KO_table))

tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)

input_table_3 <- row.names(KO_table)

SynComs <- c("AtSC","LjSC", "HvSC", "SSC")
hop_2 <- data.frame()

for (KO in input_table_3){
  for (syncom in SynComs){
    syncom_table <- read.table(paste(working_directory, "Isolate_tables/Original/", syncom,"_norm.tsv", sep= ""), sep = "\t", header =T, row.names =1)
    syncom_table_2 <- syncom_table[,grep("ES", colnames(syncom_table))]
    syncom_table_3 <- syncom_table_2[,grep(paste(syncom), colnames(syncom_table_2))]
    syncom_table_4 <- syncom_table_3[,!grepl("HL", colnames(syncom_table_3))]
    
    KO_table_sub <- KO_table[row.names(KO_table) == paste(KO),]
    
    if (syncom == "SSC"){
      KO_table_sub_2 <- KO_table_sub
    } else {
      KO_table_sub_2 <- KO_table_sub[,colnames(KO_table_sub) %in% row.names(tax_df_2)[tax_df_2$SynCom == paste(syncom)]]
    }
    
    KO_table_sub_yes <- names(KO_table_sub_2)[KO_table_sub_2 > 0]
    KO_table_sub_no <- names(KO_table_sub_2)[KO_table_sub_2 == 0]
    
    syncom_table_5 <- t(t(syncom_table_4)/rowSums(t(syncom_table_4)))
    
    syncom_table_At <- syncom_table_5[,grep("At_", colnames(syncom_table_5))]
    syncom_table_Hv <- syncom_table_5[,grep("Hv_", colnames(syncom_table_5))]
    syncom_table_Lj <- syncom_table_5[,grep("Lj_", colnames(syncom_table_5))]
    
    #Averages
    syncom_table_At_2 <- rowSums(syncom_table_At)/length(colnames(syncom_table_At))
    syncom_table_Hv_2 <- rowSums(syncom_table_Hv)/length(colnames(syncom_table_Hv))
    syncom_table_Lj_2 <- rowSums(syncom_table_Lj)/length(colnames(syncom_table_Lj))
    
    At_RA <- sum(syncom_table_At_2[names(syncom_table_At_2) %in% KO_table_sub_yes])
    Hv_RA <- sum(syncom_table_Hv_2[names(syncom_table_Hv_2) %in% KO_table_sub_yes])
    Lj_RA <- sum(syncom_table_Lj_2[names(syncom_table_Lj_2) %in% KO_table_sub_yes])
    
    syncom_table_inp <- syncom_table[,grep("Input", colnames(syncom_table))]
    syncom_table_inp_2 <- t(t(syncom_table_inp)/rowSums(t(syncom_table_inp)))
    
    syncom_table_inp_3  <- rowSums(syncom_table_inp_2)/length(colnames(syncom_table_inp_2))
    
    Input_RA <- sum(syncom_table_inp_3[names(syncom_table_inp_3) %in% KO_table_sub_yes])
    
    if(Input_RA == 0){
      At_val <- At_RA
      Hv_val <- Hv_RA
      Lj_val <- Lj_RA
    } else {
      At_val <- At_RA/Input_RA
      Hv_val <- Hv_RA/Input_RA
      Lj_val <- Lj_RA/Input_RA
    }
    
    len_val <- length(KO_table_sub_yes)/length(KO_table_sub_2)
    
    hop <- t(data.frame(c(paste(KO), At_val, Hv_val, Lj_val, len_val, paste(syncom), length(colnames(syncom_table_At)), length(colnames(syncom_table_Hv)),length(colnames(syncom_table_Lj)))))
    
    hop_2 <- rbind(hop_2, hop)
  }
}

row.names(hop_2) <- NULL
colnames(hop_2) <- c("KO", "Arabidopsis", "Barley", "Lotus", "Proportion_of_strains", "SynCom", "No_of_samples_At", "No_of_samples_Hv", "No_of_samples_Lj")

write.table(hop_2, paste(results_sa.dir, "ternary_KOs_av_med_all.txt", sep = ""), col.names=T, row.names =T, quote =F, sep = "\t")

###Generation of ternary_KOs_av_med_all_no_dom.txt - FC of isolates for each KO (split over 4 SynComs) - No dominances dataset =====
KO_table = read.table(paste(working_directory, "KO_genome/KO_SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
colnames(KO_table) <- gsub("X", "", colnames(KO_table))

tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)

input_table_3 <- row.names(KO_table)

SynComs <- c("AtSC","LjSC", "HvSC", "SSC")
hop_2 <- data.frame()

for (KO in input_table_3){
  for (syncom in SynComs){
    syncom_table <- read.table(paste(working_directory,"Isolate_tables/No_dominances/", syncom,"_norm.tsv", sep= ""), sep = "\t", header =T, row.names =1)
    syncom_table_2 <- syncom_table[,grep("ES", colnames(syncom_table))]
    syncom_table_3 <- syncom_table_2[,grep(paste(syncom), colnames(syncom_table_2))]
    syncom_table_4 <- syncom_table_3[,!grepl("HL", colnames(syncom_table_3))]
    
    KO_table_sub <- KO_table[row.names(KO_table) == paste(KO),]
    
    if (syncom == "SSC"){
      KO_table_sub_2 <- KO_table_sub
    } else {
      KO_table_sub_2 <- KO_table_sub[,colnames(KO_table_sub) %in% row.names(tax_df_2)[tax_df_2$SynCom == paste(syncom)]]
    }
    
    KO_table_sub_yes <- names(KO_table_sub_2)[KO_table_sub_2 > 0]
    KO_table_sub_no <- names(KO_table_sub_2)[KO_table_sub_2 == 0]
    
    syncom_table_5 <- t(t(syncom_table_4)/rowSums(t(syncom_table_4)))
    
    syncom_table_At <- syncom_table_5[,grep("At_", colnames(syncom_table_5))]
    syncom_table_Hv <- syncom_table_5[,grep("Hv_", colnames(syncom_table_5))]
    syncom_table_Lj <- syncom_table_5[,grep("Lj_", colnames(syncom_table_5))]
    
    #Averages
    syncom_table_At_2 <- rowSums(syncom_table_At)/length(colnames(syncom_table_At))
    syncom_table_Hv_2 <- rowSums(syncom_table_Hv)/length(colnames(syncom_table_Hv))
    syncom_table_Lj_2 <- rowSums(syncom_table_Lj)/length(colnames(syncom_table_Lj))
    
    At_RA <- sum(syncom_table_At_2[names(syncom_table_At_2) %in% KO_table_sub_yes])
    Hv_RA <- sum(syncom_table_Hv_2[names(syncom_table_Hv_2) %in% KO_table_sub_yes])
    Lj_RA <- sum(syncom_table_Lj_2[names(syncom_table_Lj_2) %in% KO_table_sub_yes])
    
    syncom_table_inp <- syncom_table[,grep("Input", colnames(syncom_table))]
    syncom_table_inp_2 <- t(t(syncom_table_inp)/rowSums(t(syncom_table_inp)))
    
    syncom_table_inp_3  <- rowSums(syncom_table_inp_2)/length(colnames(syncom_table_inp_2))
    
    Input_RA <- sum(syncom_table_inp_3[names(syncom_table_inp_3) %in% KO_table_sub_yes])
    
    if(Input_RA == 0){
      At_val <- At_RA
      Hv_val <- Hv_RA
      Lj_val <- Lj_RA
    } else {
      At_val <- At_RA/Input_RA
      Hv_val <- Hv_RA/Input_RA
      Lj_val <- Lj_RA/Input_RA
    }
    
    len_val <- length(KO_table_sub_yes)/length(KO_table_sub_2)
    
    hop <- t(data.frame(c(paste(KO), At_val, Hv_val, Lj_val, len_val, paste(syncom), length(colnames(syncom_table_At)), length(colnames(syncom_table_Hv)),length(colnames(syncom_table_Lj)))))
    
    hop_2 <- rbind(hop_2, hop)
  }
}

row.names(hop_2) <- NULL
colnames(hop_2) <- c("KO", "Arabidopsis", "Barley", "Lotus", "Proportion_of_strains", "SynCom", "No_of_samples_At", "No_of_samples_Hv", "No_of_samples_Lj")

write.table(hop_2, paste(results_sa.dir,"ternary_KOs_av_med_all_no_dom.txt",sep = ""), col.names=T, row.names =T, quote =F, sep = "\t")

###Generation of Abundances_full.txt and Family_piedonuts_full.tsv - Gene/KO, Isolate, and family abundances of KOs (on basis of S20-S24) - Original dataset =====
hop_4 <- read.table(paste(working_directory, "Shiny_app/ternary_KOs_av_med_all.txt", sep = ""), header =T, sep ="\t", row.names = 1)

tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)

samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
samples_df_2$Condition[samples_df_2$Condition == "At"] <- "Arabidopsis"
samples_df_2$Condition[samples_df_2$Condition == "Lj"] <- "Lotus"
samples_df_2$Condition[samples_df_2$Condition == "Hv"] <- "Barley"

samples_df_3 <- samples_df_2[samples_df_2$Compartment == "ES",]
samples_df_4 <- samples_df_3[samples_df_3$Inoculum != "NS",]

SynComs <- c("AtSC","HvSC","LjSC", "SSC")

KO_table = read.table(paste(working_directory, "KO_genome/KO_SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
colnames(KO_table) <- gsub("X", "", colnames(KO_table))
colnames(KO_table) <- gsub("M.1", "M-1", colnames(KO_table))
colnames(KO_table) <- gsub("M.6", "M-6", colnames(KO_table))

together_2 <- data.frame()
fams <- data.frame()

for (gene in row.names(KO_table)){
  hop_4_sub <- hop_4[hop_4$KO == paste(gene),]
  print(gene)
  
  for (syncom in SynComs) {
    norm_KO = read.table(paste(working_directory, "KO_tables/Original/", syncom, ".tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
    norm_iso = read.table(paste(working_directory, "Isolate_tables/Original/", syncom, "_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
    
    if (syncom == "SSC"){
      KO_table_2 <- KO_table
    } else {
      KO_table_2 <- KO_table[,colnames(KO_table) %in% row.names(tax_df_2)[tax_df_2$SynCom == paste(syncom)]]
    }
    
    KO_table_3 <- KO_table_2[row.names(KO_table_2) == paste(gene),]
    KO_table_4 <- colSums(KO_table_3)
    KO_table_5 <- na.omit(names(KO_table_4)[KO_table_4 != 0])
    
    samples_df_5 <- samples_df_4[samples_df_4$Inoculum == paste(syncom),]
    samples_df_8 <- samples_df_5[!grepl("HL_orig", row.names(samples_df_5)),]
    samples_df_7 <- samples_df_5[!grepl("HL_orig", row.names(samples_df_5)),]
    
    norm_KO_2 <- norm_KO[,colnames(norm_KO) %in% row.names(samples_df_7)]
    
    norm_KO_3 <- t(t(norm_KO_2)/rowSums(t(norm_KO_2)))
    norm_KO_4 <- norm_KO_3[row.names(norm_KO_3) == paste(gene),]
    
    if (length(row.names(KO_table_3)) != 0){
      norm_iso_5 <- norm_iso[,colnames(norm_iso) %in% row.names(samples_df_7)]
      norm_iso_6 <- t(t(norm_iso_5)/rowSums(t(norm_iso_5)))
      norm_iso_7 <- norm_iso_6[row.names(norm_iso_6) %in% KO_table_5,]
      families <- unique(tax_df_2$family[row.names(tax_df_2) %in% KO_table_5])
      
      norm_iso_sub <- data.frame()
      if(length(families) > 0){
        for (fam in families){
          if (length(KO_table_5) == 1){
            norm_iso_8 <- norm_iso_7
          } else {
            if (is.matrix(norm_iso_7 == TRUE)){
              norm_iso_8 <- norm_iso_7[row.names(norm_iso_7) %in% KO_table_5[KO_table_5 %in% row.names(tax_df_2)[tax_df_2$family == paste(fam)]],]
            } else {
              norm_iso_8 <- norm_iso_7
            }
          }
          
          
          if (is.matrix(norm_iso_8) == FALSE){
            norm_iso_9 <- norm_iso_8
          } else {
            norm_iso_9 <- colSums(norm_iso_8)
          }
          
          fam_value <- sum(norm_iso_9)/length(norm_iso_9)
          fam_data <- t(data.frame(c(paste(fam), fam_value, paste(syncom))))
          norm_iso_sub <- rbind(norm_iso_sub, fam_data)
        } 
        
        norm_iso_sub$V4 <- paste(gene)
        fams <- rbind(fams, norm_iso_sub)
      }
      
      
      if (length(row.names(KO_table_3)) != 0){
        norm_iso_10 <- norm_iso[,colnames(norm_iso) %in% row.names(samples_df_8)]
        norm_iso_11 <- t(t(norm_iso_10)/rowSums(t(norm_iso_10)))
        norm_iso_12 <- norm_iso_11[row.names(norm_iso_11) %in% KO_table_5,]
      } 
      
      Groups <- c("Arabidopsis", "Barley", "Lotus")
      
      norm_KO_2 <- norm_KO[,colnames(norm_KO) %in% row.names(samples_df_8)]
      
      norm_KO_3 <- t(t(norm_KO_2)/rowSums(t(norm_KO_2)))
      norm_KO_4 <- norm_KO_3[row.names(norm_KO_3) == paste(gene),]
      
      for (group in Groups){
        
        samples_df_9 <- samples_df_8[samples_df_8$Condition == paste(group),]
        
        norm_KO_5 <- norm_KO_4[names(norm_KO_4) %in% row.names(samples_df_9)]
        norm_KO_6 <- norm_KO_5
        
        if (syncom == "SSC"){
          if (length(KO_table_5) > 1){
            At_iso <- KO_table_5[KO_table_5 %in% row.names(tax_df_2)[tax_df_2$SynCom == "AtSC"]] 
            Lj_iso <- KO_table_5[KO_table_5 %in% row.names(tax_df_2)[tax_df_2$SynCom == "LjSC"]] 
            Hv_iso <- KO_table_5[KO_table_5 %in% row.names(tax_df_2)[tax_df_2$SynCom == "HvSC"]] 
            
            if (is.matrix(norm_iso_12) == TRUE){
              norm_iso_12_At <- norm_iso_12[row.names(norm_iso_12) %in% At_iso,]
              norm_iso_12_Lj <- norm_iso_12[row.names(norm_iso_12) %in% Lj_iso,]
              norm_iso_12_Hv <- norm_iso_12[row.names(norm_iso_12) %in% Hv_iso,]
            } else {
              
              if (length(At_iso) > 0){
                norm_iso_12_At <- norm_iso_12
              } else {
                norm_iso_12_At <- NULL
              }
              if (length(Lj_iso) > 0){
                norm_iso_12_Lj <- norm_iso_12
              } else {
                norm_iso_12_Lj <- NULL
              }
              if (length(Hv_iso) > 0){
                norm_iso_12_Hv <- norm_iso_12
              } else {
                norm_iso_12_Hv <- NULL
              }
            }
            
            if (length(At_iso) > 1 & is.matrix(norm_iso_12_At) == TRUE){
              norm_iso_13_At <- colSums(norm_iso_12_At)
            } else {
              norm_iso_13_At <- sum(norm_iso_12_At)
            }
            
            if (length(Hv_iso) > 1 & is.matrix(norm_iso_12_Hv) == TRUE){
              norm_iso_13_Hv <- colSums(norm_iso_12_Hv)
            } else {
              norm_iso_13_Hv <- sum(norm_iso_12_Hv)
            }
            
            if (length(Lj_iso) > 1 & is.matrix(norm_iso_12_Lj) == TRUE){
              norm_iso_13_Lj <- colSums(norm_iso_12_Lj)
            } else {
              norm_iso_13_Lj <- sum(norm_iso_12_Lj)
            }
            
            norm_iso_14_At <- norm_iso_13_At[names(norm_iso_13_At) %in% row.names(samples_df_9)]
            norm_iso_14_Lj <- norm_iso_13_Lj[names(norm_iso_13_Lj) %in% row.names(samples_df_9)]
            norm_iso_14_Hv <- norm_iso_13_Hv[names(norm_iso_13_Hv) %in% row.names(samples_df_9)]
            
            value_iso_At <- sum(norm_iso_14_At)/length(names(norm_iso_14_At))
            value_iso_Lj <- sum(norm_iso_14_Lj)/length(names(norm_iso_14_Lj))
            value_iso_Hv <- sum(norm_iso_14_Hv)/length(names(norm_iso_14_Hv))
            
          } else if (length(KO_table_5) == 1) {
            plant_sel <- tax_df_2$SynCom[row.names(tax_df_2) == paste(KO_table_5)]
            norm_iso_14 <- norm_iso_12[names(norm_iso_12) %in% row.names(samples_df_9)]
            
            if (plant_sel == "AtSC"){
              value_iso_At <- sum(norm_iso_14)/length(names(norm_iso_14))
              value_iso_Lj <- 0
              value_iso_Hv <- 0
            } else if (plant_sel == "HvSC"){
              value_iso_At <- 0
              value_iso_Lj <- 0
              value_iso_Hv <- sum(norm_iso_14)/length(names(norm_iso_14))
            } else if (plant_sel == "LjSC"){
              value_iso_At <- 0
              value_iso_Lj <- sum(norm_iso_14)/length(names(norm_iso_14))
              value_iso_Hv <- 0
            }
          } else {
            value_iso_Hv <- 0
            value_iso_At <- 0
            value_iso_Lj <- 0
          }
        } else {
          if (length(KO_table_5) > 1){
            if(is.matrix(norm_iso_12) == TRUE){
              norm_iso_13 <- colSums(norm_iso_12)
              norm_iso_14 <- norm_iso_13[names(norm_iso_13) %in% row.names(samples_df_9)]
              value_iso <- sum(norm_iso_14)/length(names(norm_iso_14))
            } else {
              norm_iso_13 <- norm_iso_12[names(norm_iso_12) %in% row.names(samples_df_9)]
              norm_iso_14 <- sum(norm_iso_13)
              value_iso <- norm_iso_14
            }
          } else if (length(KO_table_5) == 1) {
            norm_iso_14 <- norm_iso_12[names(norm_iso_12) %in% row.names(samples_df_9)]
            value_iso <- sum(norm_iso_14)/length(names(norm_iso_14))
          } else {
            value_iso <- 0
          }
        }
        
        norm_KO_7 <- norm_KO_6[names(norm_KO_6) %in% row.names(samples_df_9)]
        
        value <- sum(norm_KO_7)/length(names(norm_KO_7))
        
        if (syncom == "SSC"){
          together_3 <- t(data.frame(c(paste(gene), paste(syncom), paste(group), as.numeric(value), "AtSC", as.numeric(value_iso_At))))
          together_4 <- t(data.frame(c(paste(gene), paste(syncom), paste(group), as.numeric(value), "HvSC", as.numeric(value_iso_Hv))))
          together_5 <- t(data.frame(c(paste(gene), paste(syncom), paste(group), as.numeric(value), "LjSC", as.numeric(value_iso_Lj))))
          together <- rbind(together_3, together_4, together_5)
        } else {
          together <- t(data.frame(c(paste(gene), paste(syncom), paste(group), as.numeric(value), paste(syncom), as.numeric(value_iso))))
        }
        
        row.names(together) <- NULL
        
        together_2 <- rbind(together_2, together)
      }
    }
  }
}

row.names(together_2) <- NULL
colnames(together_2) <- c("Gene", "SynCom", "Plant", "RA_KO", "Origin", "RA_Iso")

together_2$RA_KO[together_2$RA_KO == NaN] <- 0
together_2$RA_Iso[together_2$RA_Iso == NaN] <- 0

row.names(fams) <- NULL
colnames(fams) <- c("Family", "RA", "SynCom", "Gene")

write.table(together_2, paste(results_sa.dir, "Abundances_full.tsv", col.names =T, row.names =F, sep = "\t", quote =F))
write.table(fams, paste(results_sa.dir, "Family_piedonuts_full.tsv", col.names =T, row.names =F, sep = "\t", quote =F))

###Generation of Abundances_no_dom.txt and Family_piedonuts_full.tsv - Gene/KO, Isolate, and family abundances of KOs (on basis of S20-S24) - No dominances dataset =====
hop_4 <- read.table(paste(working_directory, "Shiny_app/ternary_KOs_av_med_all_no_dom.txt", sep = ""), header =T, sep ="\t", row.names = 1)

tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)

samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
samples_df_2$Condition[samples_df_2$Condition == "At"] <- "Arabidopsis"
samples_df_2$Condition[samples_df_2$Condition == "Lj"] <- "Lotus"
samples_df_2$Condition[samples_df_2$Condition == "Hv"] <- "Barley"

samples_df_3 <- samples_df_2[samples_df_2$Compartment == "ES",]
samples_df_4 <- samples_df_3[samples_df_3$Inoculum != "NS",]

SynComs <- c("AtSC","HvSC","LjSC", "SSC")

KO_table <- read.table(paste(working_directory, "KO_genome/KO_SSC.tsv", sep = ""), header=T, sep = "\t", row.names =1)
colnames(KO_table) <- gsub("X", "", colnames(KO_table))
colnames(KO_table) <- gsub("M.1", "M-1", colnames(KO_table))
colnames(KO_table) <- gsub("M.6", "M-6", colnames(KO_table))

together_2 <- data.frame()
fams <- data.frame()

for (gene in row.names(KO_table)){
  hop_4_sub <- hop_4[hop_4$KO == paste(gene),]
  print(gene)
  
  for (syncom in SynComs) {
    norm_KO = read.table(paste(working_directory, "KO_tables/No_dominances/", syncom, ".tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
    norm_iso = read.table(paste(working_directory, "Isolate_tables/No_dominances/", syncom,"_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
    
    if (syncom == "SSC"){
      KO_table_2 <- KO_table
    } else {
      KO_table_2 <- KO_table[,colnames(KO_table) %in% row.names(tax_df_2)[tax_df_2$SynCom == paste(syncom)]]
    }
    
    KO_table_3 <- KO_table_2[row.names(KO_table_2) == paste(gene),]
    KO_table_4 <- colSums(KO_table_3)
    KO_table_5 <- na.omit(names(KO_table_4)[KO_table_4 != 0])
    
    samples_df_5 <- samples_df_4[samples_df_4$Inoculum == paste(syncom),]
    samples_df_8 <- samples_df_5[!grepl("HL_orig", row.names(samples_df_5)),]
    samples_df_7 <- samples_df_5[!grepl("HL_orig", row.names(samples_df_5)),]
    
    norm_KO_2 <- norm_KO[,colnames(norm_KO) %in% row.names(samples_df_7)]
    
    norm_KO_3 <- t(t(norm_KO_2)/rowSums(t(norm_KO_2)))
    norm_KO_4 <- norm_KO_3[row.names(norm_KO_3) == paste(gene),]
    
    if (length(row.names(KO_table_3)) != 0){
      norm_iso_5 <- norm_iso[,colnames(norm_iso) %in% row.names(samples_df_7)]
      norm_iso_6 <- t(t(norm_iso_5)/rowSums(t(norm_iso_5)))
      norm_iso_7 <- norm_iso_6[row.names(norm_iso_6) %in% KO_table_5,]
      families <- unique(tax_df_2$family[row.names(tax_df_2) %in% KO_table_5])
      
      norm_iso_sub <- data.frame()
      if(length(families) > 0){
        for (fam in families){
          if (length(KO_table_5) == 1){
            norm_iso_8 <- norm_iso_7
          } else {
            if (is.matrix(norm_iso_7 == TRUE)){
              norm_iso_8 <- norm_iso_7[row.names(norm_iso_7) %in% KO_table_5[KO_table_5 %in% row.names(tax_df_2)[tax_df_2$family == paste(fam)]],]
            } else {
              norm_iso_8 <- norm_iso_7
            }
          }
          
          
          if (is.matrix(norm_iso_8) == FALSE){
            norm_iso_9 <- norm_iso_8
          } else {
            norm_iso_9 <- colSums(norm_iso_8)
          }
          
          fam_value <- sum(norm_iso_9)/length(norm_iso_9)
          fam_data <- t(data.frame(c(paste(fam), fam_value, paste(syncom))))
          norm_iso_sub <- rbind(norm_iso_sub, fam_data)
        } 
        
        norm_iso_sub$V4 <- paste(gene)
        fams <- rbind(fams, norm_iso_sub)
      }
      
      
      if (length(row.names(KO_table_3)) != 0){
        norm_iso_10 <- norm_iso[,colnames(norm_iso) %in% row.names(samples_df_8)]
        norm_iso_11 <- t(t(norm_iso_10)/rowSums(t(norm_iso_10)))
        norm_iso_12 <- norm_iso_11[row.names(norm_iso_11) %in% KO_table_5,]
      } 
      
      Groups <- c("Arabidopsis", "Barley", "Lotus")
      
      norm_KO_2 <- norm_KO[,colnames(norm_KO) %in% row.names(samples_df_8)]
      
      norm_KO_3 <- t(t(norm_KO_2)/rowSums(t(norm_KO_2)))
      norm_KO_4 <- norm_KO_3[row.names(norm_KO_3) == paste(gene),]
      
      for (group in Groups){
        
        samples_df_9 <- samples_df_8[samples_df_8$Condition == paste(group),]
        
        norm_KO_5 <- norm_KO_4[names(norm_KO_4) %in% row.names(samples_df_9)]
        norm_KO_6 <- norm_KO_5
        
        if (syncom == "SSC"){
          if (length(KO_table_5) > 1){
            At_iso <- KO_table_5[KO_table_5 %in% row.names(tax_df_2)[tax_df_2$SynCom == "AtSC"]] 
            Lj_iso <- KO_table_5[KO_table_5 %in% row.names(tax_df_2)[tax_df_2$SynCom == "LjSC"]] 
            Hv_iso <- KO_table_5[KO_table_5 %in% row.names(tax_df_2)[tax_df_2$SynCom == "HvSC"]] 
            
            if (is.matrix(norm_iso_12) == TRUE){
              norm_iso_12_At <- norm_iso_12[row.names(norm_iso_12) %in% At_iso,]
              norm_iso_12_Lj <- norm_iso_12[row.names(norm_iso_12) %in% Lj_iso,]
              norm_iso_12_Hv <- norm_iso_12[row.names(norm_iso_12) %in% Hv_iso,]
            } else {
              
              if (length(At_iso) > 0){
                norm_iso_12_At <- norm_iso_12
              } else {
                norm_iso_12_At <- NULL
              }
              if (length(Lj_iso) > 0){
                norm_iso_12_Lj <- norm_iso_12
              } else {
                norm_iso_12_Lj <- NULL
              }
              if (length(Hv_iso) > 0){
                norm_iso_12_Hv <- norm_iso_12
              } else {
                norm_iso_12_Hv <- NULL
              }
            }
            
            if (length(At_iso) > 1 & is.matrix(norm_iso_12_At) == TRUE){
              norm_iso_13_At <- colSums(norm_iso_12_At)
            } else {
              norm_iso_13_At <- sum(norm_iso_12_At)
            }
            
            if (length(Hv_iso) > 1 & is.matrix(norm_iso_12_Hv) == TRUE){
              norm_iso_13_Hv <- colSums(norm_iso_12_Hv)
            } else {
              norm_iso_13_Hv <- sum(norm_iso_12_Hv)
            }
            
            if (length(Lj_iso) > 1 & is.matrix(norm_iso_12_Lj) == TRUE){
              norm_iso_13_Lj <- colSums(norm_iso_12_Lj)
            } else {
              norm_iso_13_Lj <- sum(norm_iso_12_Lj)
            }
            
            norm_iso_14_At <- norm_iso_13_At[names(norm_iso_13_At) %in% row.names(samples_df_9)]
            norm_iso_14_Lj <- norm_iso_13_Lj[names(norm_iso_13_Lj) %in% row.names(samples_df_9)]
            norm_iso_14_Hv <- norm_iso_13_Hv[names(norm_iso_13_Hv) %in% row.names(samples_df_9)]
            
            value_iso_At <- sum(norm_iso_14_At)/length(names(norm_iso_14_At))
            value_iso_Lj <- sum(norm_iso_14_Lj)/length(names(norm_iso_14_Lj))
            value_iso_Hv <- sum(norm_iso_14_Hv)/length(names(norm_iso_14_Hv))
            
          } else if (length(KO_table_5) == 1) {
            plant_sel <- tax_df_2$SynCom[row.names(tax_df_2) == paste(KO_table_5)]
            norm_iso_14 <- norm_iso_12[names(norm_iso_12) %in% row.names(samples_df_9)]
            
            if (plant_sel == "AtSC"){
              value_iso_At <- sum(norm_iso_14)/length(names(norm_iso_14))
              value_iso_Lj <- 0
              value_iso_Hv <- 0
            } else if (plant_sel == "HvSC"){
              value_iso_At <- 0
              value_iso_Lj <- 0
              value_iso_Hv <- sum(norm_iso_14)/length(names(norm_iso_14))
            } else if (plant_sel == "LjSC"){
              value_iso_At <- 0
              value_iso_Lj <- sum(norm_iso_14)/length(names(norm_iso_14))
              value_iso_Hv <- 0
            }
          } else {
            value_iso_Hv <- 0
            value_iso_At <- 0
            value_iso_Lj <- 0
          }
        } else {
          if (length(KO_table_5) > 1){
            if(is.matrix(norm_iso_12) == TRUE){
              norm_iso_13 <- colSums(norm_iso_12)
              norm_iso_14 <- norm_iso_13[names(norm_iso_13) %in% row.names(samples_df_9)]
              value_iso <- sum(norm_iso_14)/length(names(norm_iso_14))
            } else {
              norm_iso_13 <- norm_iso_12[names(norm_iso_12) %in% row.names(samples_df_9)]
              norm_iso_14 <- sum(norm_iso_13)
              value_iso <- norm_iso_14
            }
          } else if (length(KO_table_5) == 1) {
            norm_iso_14 <- norm_iso_12[names(norm_iso_12) %in% row.names(samples_df_9)]
            value_iso <- sum(norm_iso_14)/length(names(norm_iso_14))
          } else {
            value_iso <- 0
          }
        }
        
        norm_KO_7 <- norm_KO_6[names(norm_KO_6) %in% row.names(samples_df_9)]
        
        value <- sum(norm_KO_7)/length(names(norm_KO_7))
        
        if (syncom == "SSC"){
          together_3 <- t(data.frame(c(paste(gene), paste(syncom), paste(group), as.numeric(value), "AtSC", as.numeric(value_iso_At))))
          together_4 <- t(data.frame(c(paste(gene), paste(syncom), paste(group), as.numeric(value), "HvSC", as.numeric(value_iso_Hv))))
          together_5 <- t(data.frame(c(paste(gene), paste(syncom), paste(group), as.numeric(value), "LjSC", as.numeric(value_iso_Lj))))
          together <- rbind(together_3, together_4, together_5)
        } else {
          together <- t(data.frame(c(paste(gene), paste(syncom), paste(group), as.numeric(value), paste(syncom), as.numeric(value_iso))))
        }
        
        row.names(together) <- NULL
        
        together_2 <- rbind(together_2, together)
      }
    }
  }
}

row.names(together_2) <- NULL
colnames(together_2) <- c("Gene", "SynCom", "Plant", "RA_KO", "Origin", "RA_Iso")

together_2$RA_KO[together_2$RA_KO == NaN] <- 0
together_2$RA_Iso[together_2$RA_Iso == NaN] <- 0

row.names(fams) <- NULL
colnames(fams) <- c("Family", "RA", "SynCom", "Gene")

write.table(together_2, paste(results_sa.dir, "Abundances_no_dom.tsv", col.names =T, row.names =F, sep = "\t", quote =F))
write.table(fams, paste(results_sa.dir, "Family_piedonuts_no_dom.tsv", col.names =T, row.names =F, sep = "\t", quote =F))

###Generation of Ternary_isolates.tsv - FC of isolates for each KO (Median value across the 4 SynComs) (on basis of figure 5b - lenient selection + host specificity) - Original dataset =====
KO_table = read.table(paste(working_directory, "KO_genome/KO_SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
colnames(KO_table) <- gsub("X", "", colnames(KO_table))

tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)

SynComs <- c("AtSC","LjSC", "HvSC", "SSC")
hop_2 <- data.frame()

for (KO in row.names(KO_table)){
  for (syncom in SynComs){
    syncom_table <- read.table(paste(working_directory, "Isolate_tables/Original/", syncom,"_norm.tsv", sep= ""), sep = "\t", header =T, row.names =1)
    syncom_table_2 <- syncom_table[,grep("ES", colnames(syncom_table))]
    syncom_table_3 <- syncom_table_2[,grep(paste(syncom), colnames(syncom_table_2))]
    syncom_table_4 <- syncom_table_3[,!grepl("HL", colnames(syncom_table_3))]
    
    KO_table_sub <- KO_table[row.names(KO_table) == paste(KO),]
    
    if (syncom == "SSC"){
      KO_table_sub_2 <- KO_table_sub
    } else {
      KO_table_sub_2 <- KO_table_sub[,colnames(KO_table_sub) %in% row.names(tax_df_2)[tax_df_2$SynCom == paste(syncom)]]
    }
    
    KO_table_sub_yes <- names(KO_table_sub_2)[KO_table_sub_2 > 0]
    KO_table_sub_no <- names(KO_table_sub_2)[KO_table_sub_2 == 0]
    
    syncom_table_5 <- t(t(syncom_table_4)/rowSums(t(syncom_table_4)))
    
    syncom_table_At <- syncom_table_5[,grep("At_", colnames(syncom_table_5))]
    syncom_table_Hv <- syncom_table_5[,grep("Hv_", colnames(syncom_table_5))]
    syncom_table_Lj <- syncom_table_5[,grep("Lj_", colnames(syncom_table_5))]
    
    #Averages
    syncom_table_At_2 <- rowSums(syncom_table_At)/length(colnames(syncom_table_At))
    syncom_table_Hv_2 <- rowSums(syncom_table_Hv)/length(colnames(syncom_table_Hv))
    syncom_table_Lj_2 <- rowSums(syncom_table_Lj)/length(colnames(syncom_table_Lj))
    
    At_RA <- sum(syncom_table_At_2[names(syncom_table_At_2) %in% KO_table_sub_yes])
    Hv_RA <- sum(syncom_table_Hv_2[names(syncom_table_Hv_2) %in% KO_table_sub_yes])
    Lj_RA <- sum(syncom_table_Lj_2[names(syncom_table_Lj_2) %in% KO_table_sub_yes])
    
    syncom_table_inp <- syncom_table[,grep("Input", colnames(syncom_table))]
    syncom_table_inp_2 <- t(t(syncom_table_inp)/rowSums(t(syncom_table_inp)))
    
    syncom_table_inp_3  <- rowSums(syncom_table_inp_2)/length(colnames(syncom_table_inp_2))
    
    Input_RA <- sum(syncom_table_inp_3[names(syncom_table_inp_3) %in% KO_table_sub_yes])
    
    if(Input_RA == 0){
      At_val <- At_RA
      Hv_val <- Hv_RA
      Lj_val <- Lj_RA
    } else {
      At_val <- At_RA/Input_RA
      Hv_val <- Hv_RA/Input_RA
      Lj_val <- Lj_RA/Input_RA
    }
    
    len_val <- length(KO_table_sub_yes)/length(KO_table_sub_2)
    
    hop <- t(data.frame(c(paste(KO), At_val, Hv_val, Lj_val, len_val, paste(syncom), length(colnames(syncom_table_At)), length(colnames(syncom_table_Hv)),length(colnames(syncom_table_Lj)))))
    
    hop_2 <- rbind(hop_2, hop)
  }
}

row.names(hop_2) <- NULL
colnames(hop_2) <- c("KO", "At_val", "Hv_val", "Lj_val", "No_of_strains", "SynCom", "No_of_samples_At", "No_of_samples_Hv", "No_of_samples_Lj")

hop_4 <- data.frame()

for (KO in row.names(KO_table)){
  hop_sub <- hop_2[hop_2$KO == paste(KO),]
  
  new_3 <- data.frame()
  for (syncom in unique(hop_sub$SynCom)){
    hop_sub_2 <- hop_sub[hop_sub$SynCom == paste(syncom),]
    At_val <- as.numeric(hop_sub_2$At_val)
    Hv_val <- as.numeric(hop_sub_2$Hv_val)
    Lj_val <- as.numeric(hop_sub_2$Lj_val)
    new_2 <- data.frame(At_val,Hv_val, Lj_val,hop_sub_2$No_of_samples_At,hop_sub_2$No_of_samples_Hv, hop_sub_2$No_of_samples_Lj, paste(syncom))
    new_3 <- rbind(new_3,new_2)
  }
  
  syncom_table_inp_3 = apply(syncom_table_inp_2, 1, median, na.rm=TRUE)
  
  #Medians
  At_val <- median(as.numeric(new_3$At_val))
  Hv_val <- median(as.numeric(new_3$Hv_val))
  Lj_val <- median(as.numeric(new_3$Lj_val))
  
  No_of_strains <- sum(as.numeric(hop_sub$No_of_strains))/length(hop_sub$No_of_strains)
  
  hop_3 <- t(data.frame(c(paste(KO), At_val,Hv_val, Lj_val, No_of_strains)))
  
  hop_4 <- rbind(hop_4, hop_3)
}

row.names(hop_4) <- NULL
colnames(hop_4) <- c("KO", "Arabidopsis", "Barley", "Lotus", "Proportion_of_strains")

write.table(hop_4, paste(results_sa.dir, "Ternary_isolates.tsv", sep = ""), quote =F, col.names =T, row.names =T, sep = "\t")

###Generation of Ternary_isolates_no_dom.tsv - FC of isolates for each KO (Median value across the 4 SynComs) (on basis of figure 5e - strict selection) - No dominances dataset =====
KO_table = read.table(paste(working_directory, "KO_genome/KO_SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
colnames(KO_table) <- gsub("X", "", colnames(KO_table))

tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)

SynComs <- c("AtSC","LjSC", "HvSC", "SSC")
hop_2 <- data.frame()

for (KO in row.names(KO_table)){
  for (syncom in SynComs){
    syncom_table <- read.table(paste(working_directory, "Isolate_tables/No_dominances/", syncom,"_norm.tsv", sep= ""), sep = "\t", header =T, row.names =1)
    syncom_table_2 <- syncom_table[,grep("ES", colnames(syncom_table))]
    syncom_table_3 <- syncom_table_2[,grep(paste(syncom), colnames(syncom_table_2))]
    syncom_table_4 <- syncom_table_3[,!grepl("HL", colnames(syncom_table_3))]
    
    KO_table_sub <- KO_table[row.names(KO_table) == paste(KO),]
    
    if (syncom == "SSC"){
      KO_table_sub_2 <- KO_table_sub
    } else {
      KO_table_sub_2 <- KO_table_sub[,colnames(KO_table_sub) %in% row.names(tax_df_2)[tax_df_2$SynCom == paste(syncom)]]
    }
    
    KO_table_sub_yes <- names(KO_table_sub_2)[KO_table_sub_2 > 0]
    KO_table_sub_no <- names(KO_table_sub_2)[KO_table_sub_2 == 0]
    
    syncom_table_5 <- t(t(syncom_table_4)/rowSums(t(syncom_table_4)))
    
    syncom_table_At <- syncom_table_5[,grep("At_", colnames(syncom_table_5))]
    syncom_table_Hv <- syncom_table_5[,grep("Hv_", colnames(syncom_table_5))]
    syncom_table_Lj <- syncom_table_5[,grep("Lj_", colnames(syncom_table_5))]
    
    #Averages
    syncom_table_At_2 <- rowSums(syncom_table_At)/length(colnames(syncom_table_At))
    syncom_table_Hv_2 <- rowSums(syncom_table_Hv)/length(colnames(syncom_table_Hv))
    syncom_table_Lj_2 <- rowSums(syncom_table_Lj)/length(colnames(syncom_table_Lj))
    
    At_RA <- sum(syncom_table_At_2[names(syncom_table_At_2) %in% KO_table_sub_yes])
    Hv_RA <- sum(syncom_table_Hv_2[names(syncom_table_Hv_2) %in% KO_table_sub_yes])
    Lj_RA <- sum(syncom_table_Lj_2[names(syncom_table_Lj_2) %in% KO_table_sub_yes])
    
    syncom_table_inp <- syncom_table[,grep("Input", colnames(syncom_table))]
    syncom_table_inp_2 <- t(t(syncom_table_inp)/rowSums(t(syncom_table_inp)))
    
    syncom_table_inp_3  <- rowSums(syncom_table_inp_2)/length(colnames(syncom_table_inp_2))
    
    Input_RA <- sum(syncom_table_inp_3[names(syncom_table_inp_3) %in% KO_table_sub_yes])
    
    if(Input_RA == 0){
      At_val <- At_RA
      Hv_val <- Hv_RA
      Lj_val <- Lj_RA
    } else {
      At_val <- At_RA/Input_RA
      Hv_val <- Hv_RA/Input_RA
      Lj_val <- Lj_RA/Input_RA
    }
    
    len_val <- length(KO_table_sub_yes)/length(KO_table_sub_2)
    
    hop <- t(data.frame(c(paste(KO), At_val, Hv_val, Lj_val, len_val, paste(syncom), length(colnames(syncom_table_At)), length(colnames(syncom_table_Hv)),length(colnames(syncom_table_Lj)))))
    
    hop_2 <- rbind(hop_2, hop)
  }
}

row.names(hop_2) <- NULL
colnames(hop_2) <- c("KO", "At_val", "Hv_val", "Lj_val", "No_of_strains", "SynCom", "No_of_samples_At", "No_of_samples_Hv", "No_of_samples_Lj")

hop_4 <- data.frame()

for (KO in row.names(KO_table)){
  hop_sub <- hop_2[hop_2$KO == paste(KO),]
  
  new_3 <- data.frame()
  for (syncom in unique(hop_sub$SynCom)){
    hop_sub_2 <- hop_sub[hop_sub$SynCom == paste(syncom),]
    At_val <- as.numeric(hop_sub_2$At_val)# * as.numeric(hop_sub_2$No_of_samples_At)
    Hv_val <- as.numeric(hop_sub_2$Hv_val)# * as.numeric(hop_sub_2$No_of_samples_Hv)
    Lj_val <- as.numeric(hop_sub_2$Lj_val)# * as.numeric(hop_sub_2$No_of_samples_Lj)
    new_2 <- data.frame(At_val,Hv_val, Lj_val,hop_sub_2$No_of_samples_At,hop_sub_2$No_of_samples_Hv, hop_sub_2$No_of_samples_Lj, paste(syncom))
    new_3 <- rbind(new_3,new_2)
  }
  
  syncom_table_inp_3 = apply(syncom_table_inp_2, 1, median, na.rm=TRUE)
  
  #Medians
  At_val <- median(as.numeric(new_3$At_val))
  Hv_val <- median(as.numeric(new_3$Hv_val))
  Lj_val <- median(as.numeric(new_3$Lj_val))
  
  No_of_strains <- sum(as.numeric(hop_sub$No_of_strains))/length(hop_sub$No_of_strains)
  
  hop_3 <- t(data.frame(c(paste(KO), At_val,Hv_val, Lj_val, No_of_strains)))
  
  hop_4 <- rbind(hop_4, hop_3)
}

row.names(hop_4) <- NULL
colnames(hop_4) <- c("KO", "Arabidopsis", "Barley", "Lotus", "Proportion_of_strains")

write.table(hop_4, paste(results_sa.dir, "Ternary_isolates_no_dom.tsv", sep =""), quote =F, col.names =T, row.names =T, sep = "\t")

###Generation of boxplots_full.txt - FC of isolates with pathway - Original dataset ======
top <- read.table(paste(working_directory,"Annotations/pathway_top.txt", sep = ""), header=F, sep="\t")
KO_to_pathway <- read.table(paste(working_directory,"Annotations/KO_to_pathway.txt", sep = ""), header=T, sep="\t")
KO_to_pathway$V3 <- top$V2[match(KO_to_pathway$V2, top$V1)]

KO_to_pathway_2 <- read.table(paste(working_directory,"Annotations/KO_to_pathway_unannotated_2.txt", sep = ""), header=F, sep="\t")
colnames(KO_to_pathway_2) <- c("KO","new_category")

for (KO in KO_to_pathway_2$KO){
  KO_to_pathway$V3[KO_to_pathway$V1 == paste(KO)] <- KO_to_pathway_2$new_category[KO_to_pathway_2$KO == paste(KO)]
}

new_table <- read.table(paste(working_directory, "Shiny_app/KOs.txt", sep = ""), header=T, sep="\t",quote ="")
KO_to_pathway$V4 <- new_table$Gene[match(KO_to_pathway$V1, new_table$KO)]
KO_to_pathway_3 <- KO_to_pathway[,c(1,4,2,3)]

KO_to_pathway_3$V5 <- NA

for (KO in top$V2){
  KO_to_pathway_3$V5[KO_to_pathway_3$V3 == paste(KO)] <- top$V3[top$V2 == paste(KO)]
}

KO_to_pathway_3$V1[is.na(KO_to_pathway_3$V1)] <- "Unknown"
KO_to_pathway_3$V2[is.na(KO_to_pathway_3$V2)] <- "Unknown"
KO_to_pathway_3$V3[is.na(KO_to_pathway_3$V3)] <- "Unknown"
KO_to_pathway_3$V4[is.na(KO_to_pathway_3$V4)] <- "Unknown"
KO_to_pathway_3$V5[is.na(KO_to_pathway_3$V5)] <- "Unknown"

colnames(KO_to_pathway_3) <- c("KO", "Description", "pathway", "Pathway_Description", "Category")

write.table(KO_to_pathway_3,paste(results_sa.dir, "Gene_descriptions.txt", sep = ""), col.names =T, row.names = F, sep = "\t", quote =F)

KO_table = read.table(paste(working_directory,"KO_genome/KO_SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
colnames(KO_table) <- gsub("X", "", colnames(KO_table))

pathways <- unique(KO_to_pathway_3$Pathway_Description)
pathways_2 <- pathways[pathways != "Unknown"]

Categories <- na.omit(unique(top$V3[top$V2 %in% pathways]))

groups <- Categories
SynComs <- c("AtSC","LjSC", "HvSC", "SSC")
Plant <- c("At", "Hv", "Lj")
hop_4 <- data.frame()

for (cat in groups){
  paths <- top$V2[top$V3 == paste(cat)]
  for (path in paths){
    KOs <- na.omit(KO_to_pathway$V1[KO_to_pathway$V3 == paste(path)])
    
    hop_2 <- data.frame()
    if(length(KOs) > 0 ){
      for (KO in KOs){
        for (syncom in SynComs){
          norm_SSC =read.table(paste(working_directory, "Isolate_tables/Original/", syncom,"_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
          norm_SSC_2 <- norm_SSC[,grep("ES", colnames(norm_SSC))]
          norm_SSC_3 <- norm_SSC_2[,grep(paste(syncom), colnames(norm_SSC_2))]
          norm_SSC_input <- norm_SSC[,grep("Input", colnames(norm_SSC))]
          
          KO_table_2 <- KO_table[row.names(KO_table) == paste(KO),]
          KO_table_3 <- KO_table_2[,colnames(KO_table_2) %in% row.names(norm_SSC_3)]
          
          for (plant in Plant){
            norm_SSC_4 <- norm_SSC_3[,grep(paste(plant, "_", sep = ""), colnames(norm_SSC_3))]
            norm_SSC_5 <- t(t(norm_SSC_4)/rowSums(t(norm_SSC_4)))
            
            No <- names(KO_table_3)[KO_table_3 == 0]
            Yes <- names(KO_table_3)[KO_table_3 != 0]
            
            norm_SSC_yes <- norm_SSC_5[row.names(norm_SSC_5) %in% Yes,]
            norm_SSC_no <- norm_SSC_5[row.names(norm_SSC_5) %in% No,]
            
            if (length(Yes) > 1){
              Yes_sum <- sum(colSums(norm_SSC_yes))/length(colSums(norm_SSC_yes))
            } else if (length(Yes) == 1){
              Yes_sum <- sum(norm_SSC_yes)/length(norm_SSC_yes)
            } else {
              Yes_sum <- 0
            }
            
            
            if (length(No) > 1){
              No_sum <- sum(colSums(norm_SSC_no))/length(colSums(norm_SSC_no))
            } else if (length(No) == 1){
              No_sum <- sum(norm_SSC_no)/length(norm_SSC_no)
            } else {
              No_sum <- 0
            }
            
            norm_SSC_input_2 <- t(t(norm_SSC_input)/rowSums(t(norm_SSC_input)))
            
            norm_SSC_input_yes <- norm_SSC_input_2[row.names(norm_SSC_input_2) %in% Yes,]
            norm_SSC_input_no <- norm_SSC_input_2[row.names(norm_SSC_input_2) %in% No,]
            
            if (length(Yes) > 1){
              Yes_sum_input <- sum(colSums(norm_SSC_input_yes))/length(colSums(norm_SSC_input_yes))
            } else if (length(Yes) == 1){
              Yes_sum_input <- sum(norm_SSC_input_yes)/length(norm_SSC_input_yes)
            } else {
              Yes_sum_input <- 0
            }
            
            if (length(No) > 1){
              No_sum_input <- sum(colSums(norm_SSC_input_no))/length(colSums(norm_SSC_input_no))
            } else if (length(No) == 1){
              No_sum_input <- sum(norm_SSC_input_no)/length(norm_SSC_input_no)
            } else {
              No_sum_input <- 0
            }
            
            
            
            hop <- t(data.frame(c(paste(KO), Yes_sum, No_sum, Yes_sum_input, No_sum_input, length(Yes), length(No), paste(plant), paste(syncom))))
            hop_2 <- rbind(hop_2, hop)
          }
        }
      }
      hop_2$V2 <- as.numeric(hop_2$V2)
      hop_2$V3 <- as.numeric(hop_2$V3)
      hop_2$V4 <- as.numeric(hop_2$V4)
      hop_2$V5 <- as.numeric(hop_2$V5)
      hop_2$V6 <- as.numeric(hop_2$V6)
      hop_2$V7 <- as.numeric(hop_2$V7)
      
      int_value <- min(hop_2$V4[hop_2$V4 != 0])
      int_value_2 <- min(hop_2$V5[hop_2$V5 != 0])
      
      hop_2$V4[hop_2$V4 == 0] <- int_value
      hop_2$V5[hop_2$V5 == 0] <- int_value_2
      
      hop_2$V10 <- (hop_2$V2/hop_2$V4)
      hop_2$V11 <- (hop_2$V3/hop_2$V5)
      
      value <- sum(hop_2$V6/(hop_2$V6 + hop_2$V7))/length(hop_2$V6)
      
      for (plant in Plant){
        for (syncom in SynComs){
          hop_2_5 <- hop_2[hop_2$V8 == paste(plant),]
          hop_2_6 <- hop_2_5[hop_2_5$V9 == paste(syncom),]
          value_path_yes <- sum(as.numeric(hop_2_6$V10))/length(hop_2_6$V10)
          value_path_no <- sum(as.numeric(hop_2_6$V11))/length(hop_2_6$V11)
          value_RA <-  sum(as.numeric(hop_2_6$V2))/length(hop_2_6$V2)
          hop_3 <- t(data.frame(c(paste(plant),value_RA, value_path_yes, value_path_no, paste(syncom), paste(path), paste(cat), value)))
          hop_4 <- rbind(hop_4, hop_3)
        }
      }
    }
  }
}

row.names(hop_4) <- NULL

write.table(hop_4, paste(results_sa.dir,"boxplots_full.txt", sep = ""), quote = F, col.names = T, row.names = T, sep = "\t")

###Generation of boxplots_no_dom.txt - FC of isolates with pathway (on basis of figure 5f and S25) - No dominances dataset ======
top <- read.table(paste(working_directory,"Annotations/pathway_top.txt", sep = ""), header=F, sep="\t")
KO_to_pathway <- read.table(paste(working_directory,"Annotations/KO_to_pathway.txt", sep = ""), header=T, sep="\t")
KO_to_pathway$V3 <- top$V2[match(KO_to_pathway$V2, top$V1)]

KO_to_pathway_2 <- read.table(paste(working_directory,"Annotations/KO_to_pathway_unannotated_2.txt", sep = ""), header=F, sep="\t")
colnames(KO_to_pathway_2) <- c("KO","new_category")

for (KO in KO_to_pathway_2$KO){
  KO_to_pathway$V3[KO_to_pathway$V1 == paste(KO)] <- KO_to_pathway_2$new_category[KO_to_pathway_2$KO == paste(KO)]
}

new_table <- read.table(paste(working_directory, "Shiny_app/KOs.txt", sep = ""), header=T, sep="\t",quote ="")
KO_to_pathway$V4 <- new_table$Gene[match(KO_to_pathway$V1, new_table$KO)]
KO_to_pathway_3 <- KO_to_pathway[,c(1,4,2,3)]

KO_to_pathway_3$V5 <- NA

for (KO in top$V2){
  KO_to_pathway_3$V5[KO_to_pathway_3$V3 == paste(KO)] <- top$V3[top$V2 == paste(KO)]
}

KO_to_pathway_3$V1[is.na(KO_to_pathway_3$V1)] <- "Unknown"
KO_to_pathway_3$V2[is.na(KO_to_pathway_3$V2)] <- "Unknown"
KO_to_pathway_3$V3[is.na(KO_to_pathway_3$V3)] <- "Unknown"
KO_to_pathway_3$V4[is.na(KO_to_pathway_3$V4)] <- "Unknown"
KO_to_pathway_3$V5[is.na(KO_to_pathway_3$V5)] <- "Unknown"

colnames(KO_to_pathway_3) <- c("KO", "Description", "pathway", "Pathway_Description", "Category")

KO_table = read.table(paste(working_directory,"KO_genome/KO_SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
colnames(KO_table) <- gsub("X", "", colnames(KO_table))

pathways <- unique(KO_to_pathway_3$Pathway_Description)
pathways_2 <- pathways[pathways != "Unknown"]

Categories <- na.omit(unique(top$V3[top$V2 %in% pathways]))

groups <- Categories
SynComs <- c("AtSC","LjSC", "HvSC", "SSC")
Plant <- c("At", "Hv", "Lj")
hop_4 <- data.frame()

for (cat in groups){
  paths <- top$V2[top$V3 == paste(cat)]
  for (path in paths){
    KOs <- na.omit(KO_to_pathway$V1[KO_to_pathway$V3 == paste(path)])
    
    hop_2 <- data.frame()
    if(length(KOs) > 0 ){
      for (KO in KOs){
        for (syncom in SynComs){
          norm_SSC =read.table(paste(working_directory, "Isolate_tables/No_dominances/", syncom,"_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
          norm_SSC_2 <- norm_SSC[,grep("ES", colnames(norm_SSC))]
          norm_SSC_3 <- norm_SSC_2[,grep(paste(syncom), colnames(norm_SSC_2))]
          norm_SSC_input <- norm_SSC[,grep("Input", colnames(norm_SSC))]
          
          KO_table_2 <- KO_table[row.names(KO_table) == paste(KO),]
          KO_table_3 <- KO_table_2[,colnames(KO_table_2) %in% row.names(norm_SSC_3)]
          
          for (plant in Plant){
            norm_SSC_4 <- norm_SSC_3[,grep(paste(plant, "_", sep = ""), colnames(norm_SSC_3))]
            norm_SSC_5 <- t(t(norm_SSC_4)/rowSums(t(norm_SSC_4)))
            
            No <- names(KO_table_3)[KO_table_3 == 0]
            Yes <- names(KO_table_3)[KO_table_3 != 0]
            
            norm_SSC_yes <- norm_SSC_5[row.names(norm_SSC_5) %in% Yes,]
            norm_SSC_no <- norm_SSC_5[row.names(norm_SSC_5) %in% No,]
            
            if (length(Yes) > 1){
              Yes_sum <- sum(colSums(norm_SSC_yes))/length(colSums(norm_SSC_yes))
            } else if (length(Yes) == 1){
              Yes_sum <- sum(norm_SSC_yes)/length(norm_SSC_yes)
            } else {
              Yes_sum <- 0
            }
            
            
            if (length(No) > 1){
              No_sum <- sum(colSums(norm_SSC_no))/length(colSums(norm_SSC_no))
            } else if (length(No) == 1){
              No_sum <- sum(norm_SSC_no)/length(norm_SSC_no)
            } else {
              No_sum <- 0
            }
            
            norm_SSC_input_2 <- t(t(norm_SSC_input)/rowSums(t(norm_SSC_input)))
            
            norm_SSC_input_yes <- norm_SSC_input_2[row.names(norm_SSC_input_2) %in% Yes,]
            norm_SSC_input_no <- norm_SSC_input_2[row.names(norm_SSC_input_2) %in% No,]
            
            if (length(Yes) > 1){
              Yes_sum_input <- sum(colSums(norm_SSC_input_yes))/length(colSums(norm_SSC_input_yes))
            } else if (length(Yes) == 1){
              Yes_sum_input <- sum(norm_SSC_input_yes)/length(norm_SSC_input_yes)
            } else {
              Yes_sum_input <- 0
            }
            
            if (length(No) > 1){
              No_sum_input <- sum(colSums(norm_SSC_input_no))/length(colSums(norm_SSC_input_no))
            } else if (length(No) == 1){
              No_sum_input <- sum(norm_SSC_input_no)/length(norm_SSC_input_no)
            } else {
              No_sum_input <- 0
            }
            
            
            
            hop <- t(data.frame(c(paste(KO), Yes_sum, No_sum, Yes_sum_input, No_sum_input, length(Yes), length(No), paste(plant), paste(syncom))))
            hop_2 <- rbind(hop_2, hop)
          }
        }
      }
      hop_2$V2 <- as.numeric(hop_2$V2)
      hop_2$V3 <- as.numeric(hop_2$V3)
      hop_2$V4 <- as.numeric(hop_2$V4)
      hop_2$V5 <- as.numeric(hop_2$V5)
      hop_2$V6 <- as.numeric(hop_2$V6)
      hop_2$V7 <- as.numeric(hop_2$V7)
      
      int_value <- min(hop_2$V4[hop_2$V4 != 0])
      int_value_2 <- min(hop_2$V5[hop_2$V5 != 0])
      
      hop_2$V4[hop_2$V4 == 0] <- int_value
      hop_2$V5[hop_2$V5 == 0] <- int_value_2
      
      hop_2$V10 <- (hop_2$V2/hop_2$V4)
      hop_2$V11 <- (hop_2$V3/hop_2$V5)
      
      value <- sum(hop_2$V6/(hop_2$V6 + hop_2$V7))/length(hop_2$V6)
      
      for (plant in Plant){
        for (syncom in SynComs){
          hop_2_5 <- hop_2[hop_2$V8 == paste(plant),]
          hop_2_6 <- hop_2_5[hop_2_5$V9 == paste(syncom),]
          value_path_yes <- sum(as.numeric(hop_2_6$V10))/length(hop_2_6$V10)
          value_path_no <- sum(as.numeric(hop_2_6$V11))/length(hop_2_6$V11)
          value_RA <-  sum(as.numeric(hop_2_6$V2))/length(hop_2_6$V2)
          hop_3 <- t(data.frame(c(paste(plant),value_RA, value_path_yes, value_path_no, paste(syncom), paste(path), paste(cat), value)))
          hop_4 <- rbind(hop_4, hop_3)
        }
      }
    }
  }
}

row.names(hop_4) <- NULL

write.table(hop_4, paste(results_sa.dir,"boxplots_no_dom.txt", sep = ""), quote = F, col.names = T, row.names = T, sep = "\t")

### Notes #### ####
#Please note that these scripts can simply be changed from the full dataset to the dataset without dominators by loading in the different table (with suffix _no_dom instead of _full)
#Please set the files in a folder and name that folder here:
working_directory_SA = paste(working_directory, "Shiny_app/", sep = "")

#And if you want to make use of the scripts to produce figures, you can set the results directory here
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

###Similar to figure 5f - General enriched pathways - Box plots strains with vs strains without =====
table <- read.table(paste(working_directory_SA,"boxplots_full.txt", sep = ""),header = T, sep = "\t", row.names =1)

colnames(table) <- c("Plant", "RA", "Have", "No_Have", "SynCom", "pathway", "Category_2", "No_of_strains")
table$Category <- paste(table$pathway, " (n = ",format(round(as.numeric(table$No_of_strains),2), nsmall = 2), ")", sep = "")

table_2 <- table[,colnames(table) != "No_Have"]
table_3 <- table[,colnames(table) != "Have"]

colnames(table_2) <- c("Plant", "RA", "value","SynCom", "pathway", "Category_2","No_of_strains", "Category")
colnames(table_3) <- c("Plant", "RA", "value","SynCom", "pathway", "Category_2","No_of_strains", "Category")

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
  RA <- sum(table_sub$RA)/length(table_sub$RA)
  
  new <- t(data.frame(c(Category_2, description, round(as.numeric(No_of_strains),2), RA,Yes_RA_change, No_RA_change, stat$p.format)))
  new_2 <- rbind(new_2, new)
}

row.names(new_2) <- NULL
colnames(new_2) <- c("Category", "Pathway", "No_of_strains", "RA","RA change", "RA change w/o pathway", "Adjusted_p_value")
new_2$Adjusted_p_value <- as.numeric(new_2$Adjusted_p_value)

new_3 <- new_2[order(new_2$Adjusted_p_value, decreasing =F),]
new_3$`RA change` <- round(as.numeric(new_3$`RA change`), 2)
new_3$`RA change w/o pathway` <- round(as.numeric(new_3$`RA change w/o pathway`), 2)

new_4 <- new_3[new_3$`RA change` > new_3$`RA change w/o pathway`,]
new_5 <- new_4[new_4$Adjusted_p_value < 0.05,]

#Here pathways of interest can be put, e.g: 
pathways_of_interest <- c("Amino sugar and nucleotide sugar metabolism","Inositol phosphate metabolism", "Biofilm formation - Escherichia coli","Biosynthesis of various antibiotics","Quorum sensing", "Flagellar assembly", "Sulfur metabolism", "Cysteine and methionine metabolism", "Starch and sucrose metabolism", "ABC transporters", "Transcriptional regulator", "Biosynthesis of enediyne antibiotics")

table_6 <- table_5[table_5$pathway %in% pathways_of_interest,]
table_7 <- table_6[table_6$Present == "Present",] 
table_8 <- table_6[table_6$Present == "Absent",] 
table_7$Present <- "Present"
table_8$Present <- "Absent"

table_9 <- rbind(table_7, table_8)

hex <- hue_pal()(length(unique(table_9$pathway))) 

swr = function(string, nwrap=23) {
  paste(strwrap(string, width=nwrap), collapse="\n")
}
swr = Vectorize(swr)

# Create line breaks in Year
table_9$Category = swr(table_9$Category)

plot_top_sig <- ggplot(table_9, aes(x=Present , y=value, color = pathway)) + 
  geom_boxplot(outlier.shape = NA)+
  ggtitle("Pathways") + 
  geom_jitter(aes(shape=Plant,size = 3), cex =3)+
  theme(legend.position="top")+ 
  scale_color_manual(values = hex) +
  guides(color = FALSE) +
  guides(size = FALSE) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(size =12), axis.title = element_text(size = 14), axis.text.y = element_text(size=12), legend.title = element_text(size=16), legend.text = element_text(size=14), plot.title = element_text(size=20)) +  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  ylab("Fold change vs Input") + 
  xlab("") +
  facet_wrap(~Category,scales = "free_y") +
  theme(strip.text = element_text(size=12))
plot_top_sig_2 <- plot_top_sig + stat_compare_means(comparisons = list(c("Absent", "Present")), method = "t.test", vjust = 1.5) +
  guides(shape = guide_legend(override.aes = list(size = 5)))
plot_top_sig_2

pdf(paste(results.dir_SA,"Box_plots.pdf", sep=""), width=16, height=8)
print(plot_top_sig_2)
dev.off()

###Similar to figure 5b - Ternary plot - Median KO value of SynComs =====
table <- read.table(paste(working_directory_SA, "Ternary_isolates.tsv", sep = ""), header =T, row.names =1, sep ="\t")

table$Arabidopsis <- as.numeric(table$Arabidopsis)
table$Barley <- as.numeric(table$Barley)
table$Lotus <- as.numeric(table$Lotus)
table$Proportion_of_strains <- as.numeric(table$Proportion_of_strains)

top <- read.table(paste(working_directory_SA,"Gene_descriptions.txt", sep = ""), header=T, sep="\t")

table$pathway <- top$Pathway_Description[match(table$KO, top$KO)]
table$category <- top$Category[match(table$KO,top$KO)]

#Ternary plot
nv = 0.005
pn = position_nudge_tern(y=nv,x=-nv/2,z=-nv/2)

hex <- hue_pal()(length(unique(table$category))) 

ternary_orig <- ggplot(data=table,aes(x=Arabidopsis,y=Barley, z=Lotus, size = Proportion_of_strains, color = category)) +
  geom_point() +
  coord_tern() +
  theme_bw()+
  scale_color_manual(values = hex) +
  geom_text(position=pn,aes(label=as.character(table$KO)),check_overlap=T, size=3, angle = 25, color = 'black') +
  ggtitle("Plant-specific KOs") + theme(plot.title = element_text(hjust = 0.5, size = 20)) + 
  labs(color = "Category", size = "Proportion of strains with KO") +
  theme(text = element_text(size=18)) + theme(legend.text=element_text(size=16)) +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
ternary_orig

pdf(paste(results.dir_SA,"ternary_all.pdf", sep=""), width=28, height=19)
print(ternary_orig)
dev.off()

###Similar to figure S19-S23 - Gene abundance bar chart =====
table = read.table(paste(working_directory_SA,"Abundances_full.tsv", sep = ""), header =T, sep = "\t")
genes = table$Gene

for (gene in genes){
  table_sub <- table [table$Gene == paste(gene),]
  table_sub$RA_KO <- as.numeric(table_sub$RA_KO)
  
  table_sub$New_column <- paste(table_sub$SynCom, table_sub$Origin, sep = "-")
  table_sub$New_column <- gsub("HvSC-HvSC", "HvSC", table_sub$New_column)
  table_sub$New_column <- gsub("AtSC-AtSC", "AtSC", table_sub$New_column)
  table_sub$New_column <- gsub("LjSC-LjSC", "LjSC", table_sub$New_column)
  table_sub$New_column <- factor(table_sub$New_column, levels = c("AtSC", "HvSC", "LjSC", "SSC-AtSC", "SSC-HvSC", "SSC-LjSC"))
  
  g1 <- ggplot(table_sub %>% filter(SynCom != "SSC"),
               aes(x= Plant, weight=RA_KO, fill=New_column)) +
    theme_classic() +
    geom_bar(position = "dodge", width=0.5, just = 0.5) + 
    labs(x="Plant") +
    scale_fill_manual(values = c("#A3A500","#00B0F6","#00BF7D")) +
    ggtitle(paste(gene)) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    ylab("Relative abundance - Gene") + 
    xlab("Plant") +
    labs(fill = "SynCom") +
    theme(axis.text.x = element_text(size = 14, hjust =0.2), axis.title.y = element_text(size = 18),axis.title.x = element_blank(), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=24))
  g1
  
  g2 <- g1 + geom_bar(data=table_sub %>% filter(SynCom == "SSC"),
                      aes(x=Plant, fill=New_column),
                      position=position_stacknudge(x = 0.335), width=0.17) +
    scale_fill_manual(values = c("#A3A500","#00B0F6","#00BF7D","#fcddd9", "#fabab3","#F8766D"))
  
  pdf(paste(results.dir_SA,"RA_", gene,".pdf", sep=""), width=8, height=4)
  print(g2)
  dev.off()
}

###Similar to figure S19-S23 - Isolate abundance bar chart =====
table = read.table(paste(working_directory_SA,"Abundances_full.tsv", sep = ""), header =T, sep = "\t")
genes = table$Gene

for (gene in genes){
  table_sub <- table [table$Gene == paste(gene),]
  table_sub$RA_Iso <- as.numeric(table_sub$RA_Iso)
  
  table_sub$New_column <- paste(table_sub$SynCom, table_sub$Origin, sep = "-")
  table_sub$New_column <- gsub("HvSC-HvSC", "HvSC", table_sub$New_column)
  table_sub$New_column <- gsub("AtSC-AtSC", "AtSC", table_sub$New_column)
  table_sub$New_column <- gsub("LjSC-LjSC", "LjSC", table_sub$New_column)
  table_sub$New_column <- factor(table_sub$New_column, levels = c("AtSC", "HvSC", "LjSC", "SSC-AtSC", "SSC-HvSC", "SSC-LjSC"))
  
  g1 <- ggplot(table_sub %>% filter(SynCom != "SSC"),
               aes(x= Plant, weight=RA_Iso, fill=New_column)) +
    theme_classic() +
    geom_bar(position = "dodge", width=0.5, just = 0.5) + 
    labs(x="Plant") +
    ylim(0,1)+
    scale_fill_manual(values = c("#A3A500","#00B0F6","#00BF7D")) +
    ggtitle(paste(gene)) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    ylab("Relative abundance - Isolate") + 
    xlab("Plant") +
    labs(fill = "SynCom") +
    theme(axis.text.x = element_text(size = 14, hjust =0.2), axis.title.y = element_text(size = 18),axis.title.x = element_blank(), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=24))
  g1
  
  g2 <- g1 + geom_bar(data=table_sub %>% filter(SynCom == "SSC"),
                      aes(x=Plant, fill=New_column),position=position_stacknudge(x = 0.335), width=0.17) +
    scale_fill_manual(values = c("#A3A500","#00B0F6","#00BF7D","#fcddd9", "#fabab3","#F8766D"))
  
  pdf(paste(results.dir_SA,"RA_", gene,".pdf", sep=""), width=8, height=4)
  print(g2)
  dev.off()
}

###Similar to figure S19-S23 - Ternary plots - per gene and SynCom =====
table <- read.table(paste(working_directory_SA,"ternary_KOs_av_med_all.txt", sep = ""), header =T, sep ="\t", row.names = 1)
SynComs <- c("AtSC", "LjSC", "HvSC", "SSC")

for (gene in unique(table$KO)){
  
  table_2 <- table[table$KO == paste(gene),]
  
  table_gene_2 <-data.frame()
  for (syncom in SynComs){
    table_sub_2 <- table_2[table_2$SynCom == paste(syncom),]
    At_val <- sum(table_sub_2$Arabidopsis)/length(table_sub_2$Arabidopsis)
    Hv_val <- sum(table_sub_2$Barley)/length(table_sub_2$Barley)
    Lj_val <- sum(table_sub_2$Lotus)/length(table_sub_2$Lotus)
    Prop_val <- sum(table_sub_2$Proportion_of_strains)/length(table_sub_2$Proportion_of_strains)
    
    table_gene <- t(data.frame(c(paste(gene), At_val, Hv_val, Lj_val, Prop_val, paste(syncom))))
    table_gene_2 <- rbind(table_gene_2, table_gene)
  }
  
  row.names(table_gene_2) <- NULL
  colnames(table_gene_2) <- c("gene", "Arabidopsis", "Barley", "Lotus", "Proportion_of_strains", "SynCom")
  
  table_gene_2$Arabidopsis <- as.numeric(table_gene_2$Arabidopsis)
  table_gene_2$Barley <- as.numeric(table_gene_2$Barley)
  table_gene_2$Lotus <- as.numeric(table_gene_2$Lotus)
  table_gene_2$Proportion_of_strains <- as.numeric(table_gene_2$Proportion_of_strains)
  nv = 0.005
  pn = position_nudge_tern(y=nv,x=-nv/2,z=-nv/2)
  
  ternary <- ggtern(data=table_gene_2,aes(x=Arabidopsis,y=Barley, z=Lotus, color = SynCom)) +
    geom_point(size = 6) +
    theme_bw()+
    scale_color_manual(values =c("#A3A500","#00B0F6","#00BF7D","#F8766D") ) +
    ggtitle(paste(gene)) + 
    theme(plot.title = element_text(hjust = 0.5, size = 20)) + 
    labs(color = "Gene") +
    theme(text = element_text(size=18)) + theme(legend.text=element_text(size=16)) +
    theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
  ternary 
  
  pdf(paste(results.dir_SA,"ternary_", gene,".pdf", sep=""), width=8, height=6)
  print(ternary)
  dev.off()
}

###Similar to figure S19-S23 - Family PieDonut plots =====
#Load in tax file
tax_df = read.table(paste(working_directory_SA,"SSC_taxonomy_GTDB.tsv", sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)

#Set fixed colors for every family
Fam_colors <- data.frame(unique(tax_df_2$family))
colnames(Fam_colors) <- "Family"
hex <- hue_pal()(length(Fam_colors$Family)) 
Fam_colors$Colors <- hex

#Load in family file
fams = read.table(paste(working_directory_SA,"Family_piedonuts_full.tsv", sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)

source(paste(working_directory_SA, "PieDonutCustom_fams_GS.R", sep = ""))

for (gene in unique(fams$Gene)){
  fams_2 <- fams[fams$Gene == paste(gene),]
  fams_3 <- fams_2 %>% dplyr::select (-Gene)
  fams_4 <- fams_3[c(3,1,2)]
  fams_4$RA <- as.numeric(fams_4$RA)
  
  fams_sub_2 <- fams_4
  fams_sub_2$Color <- Fam_colors$Colors[match(fams_sub_2$Family, Fam_colors$Family)]
  
  fams_2$Rel_RA_2 <- round(as.numeric(fams_2$RA)/sum(as.numeric(fams_2$RA))*10000,0)
  fams_2$Combination <- paste(fams_2$SynCom, fams_2$Family,sep = "_")
  
  pie_data_2 <- data.frame()
  
  for (combi in fams_2$Combination){
    new <- fams_2$Rel_RA_2[fams_2$Combination == paste(combi)]
    syncom <- fams_2$SynCom[fams_2$Combination == paste(combi)]
    family <- fams_2$Family[fams_2$Combination == paste(combi)]
    
    for (i in 1:new){
      pie_data <- data.frame(paste(syncom), paste(family))
      pie_data_2 <- rbind(pie_data_2, pie_data)
    }
  }
  
  colnames(pie_data_2) <- c("SynCom", "Family")
  
  fams_sub_2$Combination <- paste(fams_sub_2$SynCom, fams_sub_2$Family, sep = "_")
  fams_sub_3 <- fams_sub_2[order(fams_sub_2$Combination),]
  colors <- fams_sub_3$Color
  
  SynCom_colors <- data.frame(c("AtSC", "HvSC", "LjSC", "SSC"),c("#A3A500","#00B0F6","#00BF7D","#F8766D"))
  colnames(SynCom_colors) <- c("SynCom", "Colour")
  SynCom_colors_2 <- SynCom_colors$Colour[SynCom_colors$SynCom %in% unique(fams_2$SynCom)]
  
  pdf(paste(results.dir_SA,"pie_", gene,".pdf", sep=""), width=6, height=6)
  print(PieDonutCustom_fams(pie_data_2,aes(pies=SynCom,donuts=Family),showRatioThreshold = 0.02))
  dev.off()
}
