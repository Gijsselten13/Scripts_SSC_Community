library("dplyr") #Version 1.1.2
library("ggplot2") #Version 3.4.2
library("ggrepel") #Version 0.9.3

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

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
