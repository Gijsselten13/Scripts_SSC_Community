library("ggtern") #Version 3.4.1

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

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
