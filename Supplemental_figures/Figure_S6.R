library("dplyr") #Version 1.1.2
library("phyloseq") #Version 1.44.0
library("reshape2") #Version 1.4.4
library("stringr") #Version 1.5.0
library("ggplot2") #Version 3.4.2
library("ggpubr") #Version 0.6.0

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Figure S6 - Taxonomic profile - SynCom color =====
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
samples_df = read.table(paste(working_directory,"SSC_R2_metadata.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
colnames(samples_df)[6]="Nutrient"

#Fuse and rename
samples_df$Sample=paste(samples_df$Inoculum,samples_df$Condition,samples_df$Compartment, samples_df$Nutrient, samples_df$Experiment, sep ="_")
samples_df$Plant_Inoculum_compartment=paste(samples_df$Condition, samples_df$Inoculum, samples_df$Compartment, sep ="_")
samples_df$Sample[samples_df$Sample == "SSC_Input_Input_Input_R1"] <- "SSC_Input_R1"
samples_df$Sample[samples_df$Sample == "AtSC_Input_Input_Input_R1"] <- "AtSC_Input_R1"
samples_df$Sample[samples_df$Sample == "HvSC_Input_Input_Input_R1"] <- "HvSC_Input_R1"
samples_df$Sample[samples_df$Sample == "LjSC_Input_Input_Input_R1"] <- "LjSC_Input_R1"
samples_df$Sample[samples_df$Sample == "SSC_Input_Input_Input_R2"] <- "SSC_Input_R2"
samples_df$Sample[samples_df$Sample == "AtSC_Input_Input_Input_R2"] <- "AtSC_Input_R2"
samples_df$Sample[samples_df$Sample == "HvSC_Input_Input_Input_R2"] <- "HvSC_Input_R2"
samples_df$Sample[samples_df$Sample == "LjSC_Input_Input_Input_R2"] <- "LjSC_Input_R2"

# TAX = tax_table(tax_mat)
TAX = tax_table(as.matrix(tax_df_2))
samples = sample_data(samples_df)

SynComs <- c("AtSC", "LjSC", "HvSC", "SSC")

barplot_df_4 <- data.frame()

for (inoculum in SynComs){
  #Set the OTU, TAX and sample data for making phyloseq object
  norm_otu=read.table(paste(working_directory,"Isolate_tables/Original/SSC_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
  OTU = otu_table(as.matrix(norm_otu),taxa_are_rows = TRUE)
  
  samples_df_2 <-subset(samples_df, samples_df$Inoculum == paste(inoculum))
  samples_df_inp <- subset(samples_df_2, samples_df_2$Compartment == "Input")
  samples_df_nod <- subset(samples_df_2, samples_df_2$Compartment == "NOD")
  samples_df_es <- subset(samples_df_2, samples_df_2$Compartment == "ES")
  samples_df_sub <- rbind(samples_df_inp,samples_df_nod, samples_df_es)
  
  df <- as.data.frame(lapply(sample_data(samples_df_sub),function (y) if(class(y)!="factor" ) as.factor(y) else y),stringsAsFactors=T)
  row.names(df) <- row.names(samples_df_sub)
  
  samples_sub = sample_data(df)
  phylo_inoculum = phyloseq(OTU,TAX, samples_sub)
  
  #order taxa by abundance
  top100 <- names(sort(taxa_sums(phylo_inoculum), decreasing=TRUE))[1:1000]
  phylo_fraction <- prune_taxa(taxa = top100, x = phylo_inoculum)
  
  phylo_fraction_2 <- merge_samples(phylo_fraction, group = "Sample")
  Lalala=t(as.data.frame(phylo_fraction_2@otu_table))
  
  lololo=merge.data.frame(x = Lalala, tax_df, by = 0 )
  rownames(lololo)=lololo$Row.names
  lololo=lololo[,-1]
  
  barplot_df=melt(data = lololo)
  Hank_the_normalizer <- function(df,group,amount){
    df_2 <- df %>% dplyr::group_by_at(group) %>% dplyr::summarise(total=sum(.data[[amount]]))
    df_3 <- df_2$total
    names(df_3) <- df_2[[group]]
    df$total <- df_3[as.character(df[[group]])]
    df$Rel <- df[[amount]] / df$total
    return(df)
  }
  barplot_df_2 <- Hank_the_normalizer(barplot_df,"variable","value")
  barplot_df_3 <- cbind(barplot_df_2,str_split_fixed(barplot_df_2$variable,'_',5))
  colnames(barplot_df_3)[colnames(barplot_df_3) == "1"] <- "SynCom_2"
  colnames(barplot_df_3)[colnames(barplot_df_3) == "2"] <- "Plant"
  colnames(barplot_df_3)[colnames(barplot_df_3) == "3"] <- "Compartment"
  colnames(barplot_df_3)[colnames(barplot_df_3) == "4"] <- "Nutrient"
  colnames(barplot_df_3)[colnames(barplot_df_3) == "5"] <- "Experiment"
  
  barplot_df_4 <- rbind(barplot_df_4,barplot_df_3)
}

barplot_df_5 <- barplot_df_4[barplot_df_4$Nutrient == "low",]
barplot_df_5_sub <- barplot_df_4[barplot_df_4$Plant == "Input",]
barplot_df_5 <- rbind(barplot_df_5, barplot_df_5_sub)
barplot_df_5$variable <- factor(barplot_df_5$variable, levels = c("AtSC_Input_R1", "AtSC_Input_R2", "AtSC_At_ES_low_R1", "AtSC_At_RZ_low_R2", "AtSC_Hv_ES_low_R1", "AtSC_Hv_ES_low_R2", "AtSC_Lj_ES_low_R1", "AtSC_Lj_ES_low_R2", "HvSC_Input_R1", "HvSC_Input_R2", "HvSC_At_ES_low_R1", "HvSC_At_ES_low_R2", "HvSC_Hv_ES_low_R1", "HvSC_Hv_ES_low_R2", "HvSC_Lj_ES_low_R1", "HvSC_Lj_ES_low_R2", "HvSC_Lj_NOD_low_R2", "LjSC_Input_R1", "LjSC_Input_R2", "LjSC_At_ES_low_R1", "LjSC_At_ES_low_R2", "LjSC_Hv_ES_low_R1", "LjSC_Hv_ES_low_R2", "LjSC_Lj_ES_low_R1", "LjSC_Lj_ES_low_R2", "LjSC_Lj_NOD_low_R2", "SSC_Input_R1", "SSC_Input_R2", "SSC_At_ES_low_R1", "SSC_At_ES_low_R2", "SSC_Hv_ES_low_R1", "SSC_Hv_ES_low_R2", "SSC_Lj_ES_low_R1", "SSC_Lj_ES_low_R2", "SSC_Lj_NOD_low_R2"))

barplot_df_6 <- barplot_df_4[barplot_df_4$Nutrient == "high",]
barplot_df_6_sub <- barplot_df_4[barplot_df_4$Plant == "Input",]
barplot_df_6 <- rbind(barplot_df_6, barplot_df_6_sub)
barplot_df_6$variable <- factor(barplot_df_6$variable, levels = c("AtSC_Input_R1", "AtSC_Input_R2", "AtSC_At_ES_high_R2", "AtSC_Hv_ES_high_R2", "AtSC_Lj_ES_high_R2", "HvSC_Input_R1", "HvSC_Input_R2", "HvSC_At_ES_high_R2", "HvSC_Hv_ES_high_R2", "HvSC_Lj_ES_high_R2", "HvSC_Lj_NOD_high_R2", "LjSC_Input_R1", "LjSC_Input_R2", "LjSC_At_ES_high_R2", "LjSC_Hv_ES_high_R2", "LjSC_Lj_ES_high_R2", "LjSC_Lj_NOD_high_R2", "SSC_Input_R1", "SSC_Input_R2", "SSC_At_ES_high_R2", "SSC_Hv_ES_high_R2", "SSC_Lj_ES_high_R2", "SSC_Lj_NOD_high_R2"))
barplot_df_6$Plant <- factor(barplot_df_6$Plant, levels = c("Input", "At", "Hv", "Lj"))

bar <- ggplot(barplot_df_5, aes(fill=SynCom, y=Rel, x=variable)) + 
  geom_bar(position="stack", stat="identity") +  ggtitle("Taxonomic profile") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_classic() +
  labs(x ="Sample", y = "Relative abundance - low nutrient samples", fill = "Inoculum") +
  theme(panel.background=element_blank(),panel.grid=element_blank(),axis.line.x=element_line(size=.5, colour="black"),axis.line.y=element_line(size=.5, colour="black"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black", size=7),legend.position="right",legend.background=element_blank(),legend.key=element_blank(),legend.text= element_text(size=10),text=element_text(family="sans", size=10))+
  theme(axis.text.x = element_text(size = 14, angle = 25,hjust=1),axis.title.x = element_blank(), axis.title.y = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=18))
bar_2 <- bar+facet_wrap(~SynCom_2, scales = "free", nrow = 1 ) +theme(strip.text.x = element_text(size = 18))
bar_2

bar_3 <- ggplot(barplot_df_6, aes(fill=SynCom, y=Rel, x=variable)) + 
  geom_bar(position="stack", stat="identity") +  ggtitle("Taxonomic profile") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_classic() +
  labs(x ="Sample", y = "Relative abundance - high nutrient samples", fill = "Inoculum") +
  theme(panel.background=element_blank(),panel.grid=element_blank(),axis.line.x=element_line(size=.5, colour="black"),axis.line.y=element_line(size=.5, colour="black"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black", size=7),legend.position="right",legend.background=element_blank(),legend.key=element_blank(),legend.text= element_text(size=10),text=element_text(family="sans", size=10))+
  theme(axis.text.x = element_text(size = 14, angle = 25,hjust=1), axis.title.x = element_blank(), axis.title.y = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=18))
bar_4 <- bar_3+facet_wrap(~SynCom_2, scales = "free", nrow = 1 ) + theme(strip.text.x = element_text(size = 18))
bar_4

bar_5 <- ggarrange(bar_2, bar_4, ncol =1, nrow =2, common.legend = T, legend = "right")

bar_5

pdf(paste(results.dir,"Figure_S6_Tax_profile_cont.pdf", sep=""), width=35, height=15)
print(bar_5)
dev.off()
