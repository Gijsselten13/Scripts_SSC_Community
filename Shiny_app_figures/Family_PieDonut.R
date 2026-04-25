library("dplyr") #Version 1.1.2
library("scales") #Version 1.2.1
library("ggplot2") #Version 3.4.2
library("webr") #Version 0.1.6

working_directory <- ""
working_directory_SA = paste(working_directory, "Shiny_app/", sep = "")
dir.create(paste(working_directory_SA, "results", sep = ""))
results.dir_SA <- paste(working_directory_SA,"results/", sep = "")

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
