library("ggtern") #Version 3.4.1
library("scales") #Version 1.2.1

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

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
