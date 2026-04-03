working_directory <- ""
results.dir <- paste(working_directory,"results/", sep = "")

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

