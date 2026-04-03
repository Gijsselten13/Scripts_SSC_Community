working_directory <- ""
results.dir <- paste(working_directory,"results/", sep = "")

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

#Figure 1d and 1e together (requires running script of 1d first)
combined_plot <- plot_1d_venn + plot_1e + plot_layout(ncol = 2)

print(combined_plot)

pdf(paste(results.dir,"Figure_1D_1E_together.pdf", sep=""), width=16, height=6)
print(combined_plot)
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
