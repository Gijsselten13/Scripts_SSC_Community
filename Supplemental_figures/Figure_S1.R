library("ggplot2") #Version 3.4.2
library("cowplot") #Version 1.1.3

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

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
