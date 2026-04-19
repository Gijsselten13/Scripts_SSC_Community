library("dplyr") #Version 1.1.2
library("phyloseq") #Version 1.44.0
library("vegan") #Version 2.6-4
library("ggplot2") #Version 3.4.2
library("ggpubr") #Version 0.6.0

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Figure S13b - Intrafamily diversity =====
SynComs <- c("AtSC", "LjSC", "HvSC")

fam_5 <- data.frame(matrix(NA, ncol =7))
colnames(fam_5) <- c("Isolate", "Rel", "Rel_prop", "KO", "KO_prop", "Family", "SynCom")
fam_6 <- fam_5[-1,]

for (syncom in SynComs){
  norm_table =read.table(paste(working_directory,"Isolate_tables/Original/", syncom, "_norm.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
  #KO table
  KO_table =read.table(paste(working_directory,"KO_genome/KO_", syncom, ".tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
  colnames(KO_table) <- gsub("X", "", colnames(KO_table))
  
  if (syncom == "AtSC"){
    colnames(KO_table)[grep("M.16",colnames(KO_table))] <- "M-16"
    colnames(KO_table)[grep("M.6",colnames(KO_table))] <- "M-6"
    colnames(KO_table)[grep("M.10",colnames(KO_table))] <- "M-10"
    colnames(KO_table)[grep("M.11_2",colnames(KO_table))] <- "M-11_2"
    colnames(KO_table)[394] <- "M-11"
  }
  
  #taxonomy
  tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
  rownames(tax_df) <- tax_df$isolate
  tax_df_2 <- tax_df %>% dplyr::select (-isolate)
  colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")
  
  tax_df_3 <- table(tax_df_2$Family)
  tax_df_4 <- names(tax_df_3)[tax_df_3 > 10]
  
  #Samples TABLE
  samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
  rownames(samples_df) <- samples_df$sample_id
  samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
  
  #Subset for the right SynCom
  samples_df_3 <- subset(samples_df_2, samples_df_2$Compartment == "ES")
  samples_df_4 <- subset(samples_df_3, samples_df_3$Inoculum == paste(syncom))
  
  #Subset microbiome table for the right SynCom
  norm_table_2 <- norm_table[,colnames(norm_table) %in% row.names(samples_df_4)]
  
  #Subset taxonomy accordingly 
  if (syncom == "SSC"){
    tax_df_3 <- tax_df_2[tax_df_2$Family %in% tax_df_4,]
  } else {
    tax_df_3 <- tax_df_2[tax_df_2$Family %in% tax_df_4,]
    tax_df_3 <- tax_df_3[tax_df_3$SynCom == paste(syncom),]
  }
  
  #Set the OTU, TAX and sample data for making phyloseq object
  OTU = otu_table(as.matrix(norm_table_2),taxa_are_rows = TRUE)
  #TAX = tax_table(tax_mat)
  TAX = tax_table(as.matrix(tax_df_3))
  samples_sub = sample_data(samples_df_4)
  
  phylo = phyloseq(OTU,TAX, samples_sub)
  
  phylo_RA=microbiome::transform(x = phylo, transform = "compositional" )
  ps_family <- phyloseq::tax_glom(phylo, "Family")
  phylo_RA_fam=microbiome::transform(x = ps_family, transform = "compositional" )
  
  isolate_tab <- phylo_RA@otu_table
  OTU1 = as(otu_table(phylo_RA_fam), "matrix")
  TAX1 = as.data.frame(as(tax_table(phylo_RA_fam), "matrix"))
  
  row.names(OTU1) <- TAX1$Family
  Families <- unique(tax_df_3$Family)
  
  fam_3 <- data.frame(matrix(NA, ncol =7))
  colnames(fam_3) <- c("Isolate", "Rel", "Rel_prop", "KO", "KO_prop", "Family","Family_KO")
  fam_4 <- fam_3[-1,]
  
  for (family in Families) {
    isolate_set <- row.names(tax_df_3)[tax_df_3$Family == paste(family)]
    
    isolate_set_2 <- isolate_set[isolate_set %in% row.names(isolate_tab)]
    
    fam <- data.frame(matrix(NA, ncol = 3))
    colnames(fam) <- c("Isolate", "Rel", "KO")
    fam_2 <- fam[-1,]
    
    KO_table_2 <- KO_table[,colnames(KO_table) %in% isolate_set_2]
    veg_dist <- as.matrix(vegdist(t(KO_table_2)), method = "bray", diag = T)
    veg_dist_2 <- 1-veg_dist
    
    for (isolate in isolate_set_2){
      isolate_tab_2 <-isolate_tab[row.names(isolate_tab) == paste(isolate),]
      isolate_value <- rowSums(isolate_tab_2)/length(isolate_tab_2)
      names(isolate_value) <- NULL
      
      if (length(isolate_set_2) > 1){
        KO_table_3 <- KO_table_2[, colnames(KO_table_2) == paste(isolate)]
        KO_table_4 <- KO_table_3[KO_table_3 != 0]
        KO_value <- length(KO_table_4)
      } else {
        KO_table_4 <- KO_table_2[KO_table_2 != 0]
        KO_value <- length(KO_table_4)
      }
      
      new <- t(data.frame(c(paste(isolate), isolate_value, KO_value)))
      fam_2 <- rbind(fam_2, new)
    }
    
    fam_tab_2 <- OTU1[row.names(OTU1) == paste(family),]
    fam_value <- sum(fam_tab_2)/length(fam_tab_2)
    fam_2$V2 <- as.numeric(fam_2$V2)
    fam_2$V4 <- fam_2$V2/fam_value
    
    KO_table_fam <- KO_table[,colnames(KO_table) %in% isolate_set]
    if (length(isolate_set) >1){
      KO_table_fam_2 <- rowSums(KO_table_fam)
      KO_table_fam_3 <- KO_table_fam_2[KO_table_fam_2 != 0]
      fam_KO <- length(KO_table_fam_3)
    } else {
      fam_KO <- sum(KO_table_fam)
    }
    
    fam_2$V3 <- as.numeric(fam_2$V3)
    fam_2$V5 <- fam_2$V3/fam_KO
    fam_2$V6 <- paste(family)
    fam_2$V7 <- fam_KO
    
    colnames(fam_2) <- c("Isolate", "Rel", "KO", "Rel_Prop", "KO_prop", "Family", "Family_KO")
    fam_4 <- rbind(fam_4, fam_2)
  }
  row.names(fam_4) <- NULL
  fam_4$SynCom <- paste(syncom)
  fam_6 <- rbind(fam_6, fam_4)
}

fam_6$Rel <- as.numeric(fam_6$Rel)

families <- unique(fam_6$Family)

plot_list <- list()
i <- 1

fam_6

for (family in families){
  fam_7 <- fam_6[fam_6$Family == paste(family),]
  average <- sum(fam_7$KO_prop)/length(fam_7$KO_prop)
  
  SynCom_colors <- data.frame(c("AtSC", "HvSC", "LjSC"),c("#A3A500","#00B0F6","#00BF7D"))
  colnames(SynCom_colors) <- c("SynCom", "color")            
  
  fam_8 <- fam_7[order(fam_7$SynCom),]
  fam_colors <- SynCom_colors$color[SynCom_colors$SynCom %in% fam_8$SynCom]
  
  plot <- ggplot(fam_8, aes(x = Rel_Prop, y = KO_prop, 
                            color = as.factor(SynCom))) +
    geom_point(size = 3) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) + 
    labs(x = "Relative abundance proportion",y = paste("KO Proportion (n = ", unique(fam_8$Family_KO), ")", sep = ""), color = "Inoculum") +
    ggtitle(paste(family, " - (n = ", length(fam_8$Isolate), ")", sep = "")) +
    ylim(0.4,1) +
    scale_color_manual(values = fam_colors)+
    geom_smooth(method="nls", se=FALSE, formula=y~a*log(x)+k,
                method.args=list(start=c(a=1, k=1))) +
    scale_size_continuous(limits = c(0,0.37)) +
    theme(axis.text.x = element_text(size = 10), axis.title = element_text(size = 14), axis.text.y = element_text(size=10), legend.title = element_text(size=16), legend.text = element_text(size=12), plot.title = element_text(size=18)) +
    guides(shape = guide_legend(override.aes = list(size = 5)))
  plot_list[[i]] <- plot
  i <- i + 1
}

all_plots <- ggarrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],plot_list[[6]],plot_list[[7]],plot_list[[8]],plot_list[[9]],plot_list[[10]],plot_list[[11]],plot_list[[12]],plot_list[[13]],plot_list[[14]],plot_list[[15]],plot_list[[16]],plot_list[[17]], nrow =5, ncol = 4, common.legend = T)

pdf(paste(results.dir,"Figure_S13b_Family_plot_loglinear.pdf", sep=""), width=18, height=20)
print(all_plots)
dev.off()
