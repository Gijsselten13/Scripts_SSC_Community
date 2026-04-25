library("ggplot2") #Version 3.4.2
library("ggtern") #Version 3.4.1
library("scales") #Version 1.2.1

working_directory <- ""
working_directory_SA = paste(working_directory, "Shiny_app/", sep = "")
dir.create(paste(working_directory_SA, "results", sep = ""))
results.dir_SA <- paste(working_directory_SA,"results/", sep = "")

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
