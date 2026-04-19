library("ggplot2") #Version 3.4.2
library("ggrepel") #Version 0.9.3

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

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
