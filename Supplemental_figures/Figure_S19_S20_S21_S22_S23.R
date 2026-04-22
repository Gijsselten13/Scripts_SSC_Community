library("dplyr") #Version 1.1.2
library("ggplot2") #Version 3.4.2
library("ggpp") #Version 0.5.5
library("ggtern") #Version 3.4.1
library("scales") #Version 1.2.1
library("webr") #Version 0.1.6

working_directory <- ""
dir.create(paste(results.dir, "/Small_plots", sep = ""))
results.dir_2 <- paste(results.dir,"Small_plots/", sep = "")

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
