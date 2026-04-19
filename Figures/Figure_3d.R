library("dplyr") #Version 1.1.2
library("phyloseq") #Version 1.44.0
library("vegan") #Version 2.6-4
library("ggplot2") #Version 3.4.2

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Figure 3d - Intrafamily diversity heatmap =====
SynComs <- c("AtSC", "LjSC", "HvSC", "SSC")

fam_5 <- data.frame(matrix(NA, ncol =9))
colnames(fam_5) <- c("Isolate", "Rel", "Rel_prop", "KO", "KO_prop", "Family", "Rel_prop_Z", "KO_prop_Z", "SynCom")
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
    fam_2$V8 <- scale(fam_2$V4)
    fam_2$V9 <- scale(fam_2$V5)
    
    colnames(fam_2) <- c("Isolate", "Rel", "KO", "Rel_Prop", "KO_prop", "Family", "Family_KO", "Rel_Prop_Z","KO_prop_Z")
    fam_4 <- rbind(fam_4, fam_2)
  }
  row.names(fam_4) <- NULL
  fam_4$SynCom <- paste(syncom)
  fam_6 <- rbind(fam_6, fam_4)
}

fam_6$Rel <- as.numeric(fam_6$Rel)

families <- unique(fam_6$Family)

table_cor_2 <- data.frame()

for (family in families){
  fam_7 <- fam_6[fam_6$Family == paste(family),]
  average <- sum(fam_7$KO_prop)/length(fam_7$KO_prop)
  
  SynCom_colors <- data.frame(c("AtSC", "HvSC", "LjSC"),c("#A3A500","#00B0F6","#00BF7D"))
  colnames(SynCom_colors) <- c("SynCom", "color")            
  
  fam_8 <- fam_7[order(fam_7$SynCom),]
  fam_colors <- SynCom_colors$color[SynCom_colors$SynCom %in% fam_8$SynCom]
  
  fam_8$Rel_log <- log2(fam_8$Rel_Prop)
  fam_8$Rel_log[fam_8$Rel_log == "-Inf"] <- 0
  
  for (syncom in SynComs){
    fam_sub <- fam_8[fam_8$SynCom == paste(syncom),]
    
    if (length(fam_sub$Isolate) < 6){
      pval <- NA
      cor <- NA
    } else {
      fam_sub_2 <- fam_sub[fam_sub$Rel_Prop != 0, ]
      val <- lm(formula = KO_prop ~ log(Rel_Prop), data = fam_sub_2)
      sum <- summary(val)
      cor <- sum$adj.r.squared
      pval_sum <- sum$coefficients
      pval <- pval_sum[8]
    }
    
    table_cor <- data.frame(t(data.frame(c(paste(syncom), cor, pval, paste(family)))))
    table_cor_2 <- rbind(table_cor_2, table_cor)
  }
}

row.names(table_cor_2) <- NULL
colnames(table_cor_2) <- c("SynCom", "R2", "Pvalue", "Family")

table_cor_2$R2 <- as.numeric(table_cor_2$R2)
table_cor_2$Pvalue <- as.numeric(table_cor_2$Pvalue)

table_cor_2$Sig <- ""
table_cor_2$Sig[table_cor_2$Pvalue < 0.05] <- "*"

new_fam <- data.frame()
#to get the order
for (family in families){
  table_cor_sub <- table_cor_2[table_cor_2$Family == paste(family),]
  table_cor_sub_2 <- table_cor_sub[!is.na(table_cor_sub$R2),]
  sum_val <- sum(table_cor_sub_2$R2)
  
  new_fam_2 <- data.frame(t(data.frame(c(paste(family), sum_val))))
  new_fam <- rbind(new_fam, new_fam_2)
}

new_fam_2 <- new_fam$X1[order(new_fam$X2, decreasing = F)]

table_cor_2$Family <- factor(table_cor_2$Family, levels = new_fam_2)

#Get the order based on family abundance in the data. 
#taxonomy
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv", sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")

tax_df_3 <- table(tax_df_2$Family)
tax_df_4 <- names(tax_df_3)[tax_df_3 > 10]

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep = ""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)

#Subset for the right SynCom
samples_df_3 <- subset(samples_df_2, samples_df_2$Compartment == "ES")
samples_df_4 <- subset(samples_df_3, samples_df_3$Condition != "NS")

#Subset microbiome table for the right SynCom
norm_table_2 <- norm_table[,colnames(norm_table) %in% row.names(samples_df_4)]

#Subset taxonomy accordingly 
tax_df_3 <- tax_df_2[tax_df_2$Family %in% tax_df_4,]

#Set the OTU, TAX and sample data for making phyloseq object
OTU = otu_table(as.matrix(norm_table_2),taxa_are_rows = TRUE)
#TAX = tax_table(tax_mat)
TAX = tax_table(as.matrix(tax_df_3))
samples_sub = sample_data(samples_df_4)

phylo = phyloseq(OTU,TAX, samples_sub)

phylo_RA=microbiome::transform(x = phylo, transform = "compositional" )
isolate_tab <- phylo_RA@otu_table

fam_order_2 <- data.frame()

for (family in families){
  tax_df_isolates <- row.names(tax_df_3)[tax_df_3$Family == paste(family)]
  isolate_tab_2 <- isolate_tab[row.names(isolate_tab) %in% tax_df_isolates,]
  average <- sum(rowSums(isolate_tab_2)/length(colnames(isolate_tab_2)))
  
  fam_order <- data.frame(t(data.frame(c(paste(family), average))))
  
  fam_order_2 <- rbind(fam_order_2, fam_order)
}

fam_order_3 <- fam_order_2$X1[order(fam_order_2$X2, decreasing =F)]
table_cor_2$Family <- factor(table_cor_2$Family, levels = fam_order_3)

fam_order_2$X3 <- "RA"
fam_order_2$X4 <- ""
fam_order_2$X5 <- "-"

fam_order_4 <- fam_order_2[,c(3,2,4,1,5)]
colnames(fam_order_4) <- colnames(table_cor_2)
table_cor_2 <- rbind(table_cor_2, fam_order_4)
row.names(table_cor_2) <- NULL
table_cor_2$R2 <- as.numeric(table_cor_2$R2)

table_cor_2$Sig[is.na(table_cor_2$Pvalue)] <- "-"
table_cor_2$Sig[table_cor_2$Pvalue < 0.05] <- "*"

table_cor_3 <- table_cor_2[table_cor_2$SynCom != "RA",]

table_cor_3$SynCom <- factor(table_cor_3$SynCom, levels = c("AtSC", "HvSC", "LjSC", "SSC"))

# Count number of isolates which belong to every family
Family_count <- data.frame(table(tax_df_2$Family))
Family_count$V2 <- paste(Family_count$Var1, " (n = ", Family_count$Freq, ")", sep ="")
# Create a named vector for mapping V1 to V2
mapping <- setNames(Family_count$V2, Family_count$Var1)
# Use the X_subset to reorder and subset the mapping to generate Y
fam_order_final <- mapping[fam_order_3]

table_cor_3$Family_2 <- Family_count$V2[match(table_cor_3$Family, Family_count$Var1)]
table_cor_3$Family_2 <- factor(table_cor_3$Family_2, levels = fam_order_final)

Plot_fam <- ggplot(table_cor_3, aes(SynCom, Family_2)) +
  geom_tile(height=0.98, mapping = aes(fill = R2)) +
  geom_text(aes(label = Sig), size =4) +
  scale_fill_gradient2(low = "#D55e00", mid = "white", high = "#56b4e9", midpoint =0, na.value = "lightgrey")+
  theme_classic() +
  labs(x ="Inoculum", y = "Family", fill = "R2") +
  theme(panel.background=element_blank(),panel.grid=element_blank(),axis.line.x=element_line(size=.5, colour="black"),axis.line.y=element_line(size=.5, colour="black"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black", size=7),legend.position="right",legend.text= element_text(size=10),text=element_text(family="sans", size=10))+
  theme(axis.text.x = element_text(size = 14, angle = 25,hjust=1),axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18), axis.text.y = element_text(size=14, face = rep("italic")), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=18)) +
  ggtitle("Root colonization versus functional diversity") +
  theme(plot.title = element_text(hjust = 0.5))
Plot_fam

pdf(paste(results.dir, "Figure_3d_family_heatmap.pdf", sep=""), width=6, height=6)
print(Plot_fam)
dev.off()

table_cor_4 <- table_cor_2[table_cor_2$SynCom == "RA",]
table_cor_4$Family <- factor(table_cor_4$Family, levels = fam_order_3)

plot_bar_fam <- ggplot(table_cor_4, aes(x=Family, y=R2)) + 
  geom_bar(stat = "identity", width = 0.98) +
  coord_flip() +
  theme(panel.background=element_blank(),panel.grid=element_blank(),axis.line.x=element_line(size=.5, colour="black"),axis.line.y=element_line(size=.5, colour="black"),axis.ticks=element_line(color="black")) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 14,angle = 25,hjust=1)) +
  ylab("Relative Abundance")
plot_bar_fam

pdf(paste(results.dir, "Figure_3d_family_bar.pdf", sep=""), width=1.5, height=6)
print(plot_bar_fam)
dev.off()

#Figure 3d - Genera heatmap
SynComs <- c("AtSC", "LjSC", "HvSC", "SSC")

genus_5 <- data.frame(matrix(NA, ncol =7))
colnames(genus_5) <- c("Isolate", "Rel", "Rel_prop", "KO", "KO_prop", "Genus", "SynCom")
genus_6 <- genus_5[-1,]

for (syncom in SynComs){
  norm_table =read.table(paste(working_directory, "Isolate_tables/Original/",syncom,"_norm.tsv",sep = ""), header=TRUE,sep="\t", row.names = 1)
  #KO table
  KO_table = read.table(paste(working_directory, "KO_genome/KO_",syncom, ".tsv", sep = ""), header = TRUE, sep = "\t", row.names =1)
  colnames(KO_table) <- gsub("X", "", colnames(KO_table))
  
  if (syncom == "AtSC"){
    colnames(KO_table)[grep("M.16",colnames(KO_table))] <- "M-16"
    colnames(KO_table)[grep("M.6",colnames(KO_table))] <- "M-6"
    colnames(KO_table)[grep("M.10",colnames(KO_table))] <- "M-10"
    colnames(KO_table)[grep("M.11_2",colnames(KO_table))] <- "M-11_2"
    colnames(KO_table)[394] <- "M-11"
  }
  
  if (syncom == "SSC"){
    colnames(KO_table)[grep("M.16",colnames(KO_table))] <- "M-16"
    colnames(KO_table)[grep("M.6",colnames(KO_table))] <- "M-6"
    colnames(KO_table)[grep("M.10",colnames(KO_table))] <- "M-10"
    colnames(KO_table)[grep("M.11_2",colnames(KO_table))] <- "M-11_2"
    colnames(KO_table)[571] <- "M-11"
  }
  
  #taxonomy
  tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv", sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
  rownames(tax_df) <- tax_df$isolate
  tax_df_2 <- tax_df %>% dplyr::select (-isolate)
  colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")
  
  tax_df_3 <- table(tax_df_2$Genus)
  tax_df_4 <- na.omit(names(tax_df_3)[tax_df_3 > 5])
  
  #Samples TABLE
  samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep = ""), header=TRUE,sep="\t") #make the SampleID column into the row.names
  rownames(samples_df) <- samples_df$sample_id
  samples_df_2 <- samples_df %>% dplyr::select (-sample_id)
  
  #Subset for the right SynCom
  samples_df_3 <- subset(samples_df_2, samples_df_2$Compartment == "ES")
  samples_df_4 <- subset(samples_df_3, samples_df_3$Inoculum == paste(syncom))
  
  #Subset microbiome table for the right SynCom
  norm_table_2 <- norm_table[,colnames(norm_table) %in% row.names(samples_df_4)]
  
  #Subset taxonomy accordingly 
  if (syncom == "SSC"){
    tax_df_3 <- tax_df_2[tax_df_2$Genus %in% tax_df_4,]
  } else {
    tax_df_3 <- tax_df_2[tax_df_2$Genus %in% tax_df_4,]
    tax_df_3 <- tax_df_3[tax_df_3$SynCom == paste(syncom),]
  }
  
  #Set the OTU, TAX and sample data for making phyloseq object
  OTU = otu_table(as.matrix(norm_table_2),taxa_are_rows = TRUE)
  #TAX = tax_table(tax_mat)
  TAX = tax_table(as.matrix(tax_df_3))
  samples_sub = sample_data(samples_df_4)
  
  phylo = phyloseq(OTU,TAX, samples_sub)
  
  phylo_RA=microbiome::transform(x = phylo, transform = "compositional" )
  ps_genus <- phyloseq::tax_glom(phylo, "Genus")
  phylo_RA_gen=microbiome::transform(x = ps_genus, transform = "compositional" )
  
  isolate_tab <- phylo_RA@otu_table
  OTU1 = as(otu_table(phylo_RA_gen), "matrix")
  TAX1 = as.data.frame(as(tax_table(phylo_RA_gen), "matrix"))
  
  row.names(OTU1) <- TAX1$Genus
  Genera <- na.omit(unique(tax_df_3$Genus))
  
  genus_3 <- data.frame(matrix(NA, ncol =7))
  colnames(genus_3) <- c("Isolate", "Rel", "Rel_prop", "KO", "KO_prop", "Genus","Family_KO")
  genus_4 <- genus_3[-1,]
  
  for (genera in Genera) {
    isolate_set <- row.names(tax_df_3)[tax_df_3$Genus == paste(genera)]
    
    isolate_set_2 <- isolate_set[isolate_set %in% row.names(isolate_tab)]
    
    genus <- data.frame(matrix(NA, ncol = 3))
    colnames(genus) <- c("Isolate", "Rel", "KO")
    genus_2 <- genus[-1,]
    
    KO_table_2 <- KO_table[,colnames(KO_table) %in% isolate_set_2]
    veg_dist <- as.matrix(vegdist(t(KO_table_2)), method = "bray", diag = T)
    veg_dist_2 <- 1-veg_dist
    
    for (isolate in isolate_set_2){
      isolate_tab_2 <-isolate_tab[row.names(isolate_tab) == paste(isolate),]
      isolate_value <- rowSums(isolate_tab_2)/length(isolate_tab_2)
      names(isolate_value) <- NULL
      
      if (length(na.omit(isolate_set_2)) > 1){
        KO_table_3 <- KO_table_2[, colnames(KO_table_2) == paste(isolate)]
        KO_table_4 <- KO_table_3[KO_table_3 != 0]
        KO_value <- length(KO_table_4)
      } else {
        KO_table_4 <- KO_table_2[KO_table_2 != 0]
        KO_value <- length(KO_table_4)
      }
      
      new <- t(data.frame(c(paste(isolate), isolate_value, KO_value)))
      genus_2 <- rbind(genus_2, new)
    }
    
    genus_tab_2 <- OTU1[row.names(OTU1) == paste(genera),]
    genus_value <- sum(genus_tab_2)/length(genus_tab_2)
    genus_2$V2 <- as.numeric(genus_2$V2)
    genus_2$V4 <- genus_2$V2/genus_value
    
    KO_table_genus <- KO_table[,colnames(KO_table) %in% isolate_set]
    if (length(na.omit(isolate_set)) >1){
      KO_table_genus_2 <- rowSums(KO_table_genus)
      KO_table_genus_3 <- KO_table_genus_2[KO_table_genus_2 != 0]
      genus_KO <- length(KO_table_genus_3)
    } else {
      genus_KO <- sum(KO_table_genus)
    }
    
    if (length(genus_2$V1) != 0){
      genus_2$V3 <- as.numeric(genus_2$V3)
      genus_2$V5 <- genus_2$V3/genus_KO
      genus_2$V6 <- paste(genera)
      genus_2$V7 <- genus_KO
      colnames(genus_2) <- c("Isolate", "Rel", "KO", "Rel_Prop", "KO_prop", "Genus", "Genus_KO")
      genus_4 <- rbind(genus_4, genus_2)
    }
  }
  
  row.names(genus_4) <- NULL
  genus_4$SynCom <- paste(syncom)
  genus_6 <- rbind(genus_6, genus_4)
}

genus_6$Rel <- as.numeric(genus_6$Rel)

Genera <- unique(genus_6$Genus)

table_cor_2 <- data.frame()

for (genera in Genera){
  genus_7 <- genus_6[genus_6$Genus == paste(genera),]
  average <- sum(genus_7$KO_prop)/length(genus_7$KO_prop)
  
  SynCom_colors <- data.frame(c("AtSC", "HvSC", "LjSC"),c("#A3A500","#00B0F6","#00BF7D"))
  colnames(SynCom_colors) <- c("SynCom", "color")            
  
  genus_8 <- genus_7[order(genus_7$SynCom),]
  genus_colors <- SynCom_colors$color[SynCom_colors$SynCom %in% genus_8$SynCom]
  
  genus_8$Rel_log <- log2(genus_8$Rel_Prop)
  genus_8$Rel_log[genus_8$Rel_log == "-Inf"] <- 0
  
  for (syncom in SynComs){
    genus_sub <- genus_8[genus_8$SynCom == paste(syncom),]
    
    if (length(genus_sub$Isolate) < 6){
      pval <- NA
      cor <- NA
    } else {
      genus_sub_2 <- genus_sub[genus_sub$Rel_Prop != 0, ]
      val <- lm(formula = KO_prop ~ log(Rel_Prop), data = genus_sub_2)
      sum <- summary(val)
      cor <- sum$adj.r.squared
      pval_sum <- sum$coefficients
      pval <- pval_sum[8]
    }
    
    table_cor <- data.frame(t(data.frame(c(paste(syncom), cor, pval, paste(genera)))))
    table_cor_2 <- rbind(table_cor_2, table_cor)
  }
}

row.names(table_cor_2) <- NULL
colnames(table_cor_2) <- c("SynCom", "R2", "Pvalue", "Genus")

table_cor_2$R2 <- as.numeric(table_cor_2$R2)
table_cor_2$Pvalue <- as.numeric(table_cor_2$Pvalue)

table_cor_2$Sig <- ""
table_cor_2$Sig[table_cor_2$Pvalue < 0.05] <- "*"

new_genus <- data.frame()
#to get the order
for (genera in Genera){
  table_cor_sub <- table_cor_2[table_cor_2$Genus == paste(genera),]
  table_cor_sub_2 <- table_cor_sub[!is.na(table_cor_sub$R2),]
  sum_val <- sum(table_cor_sub_2$R2)
  
  new_genus_2 <- data.frame(t(data.frame(c(paste(genera), sum_val))))
  new_genus <- rbind(new_genus, new_genus_2)
}

new_genus_2 <- new_genus$X1[order(new_genus$X2, decreasing = F)]

table_cor_2$Genus <- factor(table_cor_2$Genus, levels = new_genus_2)

#Get the order based on family abundance in the data. 
#taxonomy
tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv", sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")

tax_df_3 <- table(tax_df_2$Genus)
tax_df_4 <- names(tax_df_3)[tax_df_3 > 5]

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep = ""), header=TRUE,sep="\t") #make the SampleID column into the row.names
rownames(samples_df) <- samples_df$sample_id
samples_df_2 <- samples_df %>% dplyr::select (-sample_id)

#Subset for the right SynCom
samples_df_3 <- subset(samples_df_2, samples_df_2$Compartment == "ES")
samples_df_4 <- subset(samples_df_3, samples_df_3$Condition != "NS")

#Subset microbiome table for the right SynCom
norm_table_2 <- norm_table[,colnames(norm_table) %in% row.names(samples_df_4)]

#Subset taxonomy accordingly 
tax_df_3 <- tax_df_2[tax_df_2$Genus %in% tax_df_4,]

#Set the OTU, TAX and sample data for making phyloseq object
OTU = otu_table(as.matrix(norm_table_2),taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(tax_df_3))
samples_sub = sample_data(samples_df_4)

phylo = phyloseq(OTU,TAX, samples_sub)

phylo_RA=microbiome::transform(x = phylo, transform = "compositional" )
isolate_tab <- phylo_RA@otu_table

genus_order_2 <- data.frame()

for (genera in Genera){
  tax_df_isolates <- row.names(tax_df_3)[tax_df_3$Genus == paste(genera)]
  isolate_tab_2 <- isolate_tab[row.names(isolate_tab) %in% tax_df_isolates,]
  average <- sum(rowSums(isolate_tab_2)/length(colnames(isolate_tab_2)))
  
  genus_order <- data.frame(t(data.frame(c(paste(genera), average))))
  
  genus_order_2 <- rbind(genus_order_2, genus_order)
}

genus_order_3 <- genus_order_2$X1[order(genus_order_2$X2, decreasing =F)]

genus_order_2$X3 <- "RA"
genus_order_2$X4 <- ""
genus_order_2$X5 <- ""

genus_order_4 <- genus_order_2[,c(3,2,4,1,5)]
colnames(genus_order_4) <- colnames(table_cor_2)
table_cor_2 <- rbind(table_cor_2, genus_order_4)
row.names(table_cor_2) <- NULL
table_cor_2$R2 <- as.numeric(table_cor_2$R2)

table_cor_2$Sig[is.na(table_cor_2$Pvalue)] <- "-"
table_cor_2$Sig[table_cor_2$Pvalue < 0.05] <- "*"

table_cor_3 <- table_cor_2[table_cor_2$SynCom != "RA",]

table_cor_3$Genus <- factor(table_cor_3$Genus, levels = genus_order_3)
table_cor_3$SynCom <- factor(table_cor_3$SynCom, levels = c("AtSC", "HvSC", "LjSC", "SSC"))

table_cor_4 <- table_cor_3[table_cor_3$Genus %in% Genera,]
table_cor_4$Genus <- tax_df_2$Genus[match(table_cor_4$Genus, tax_df_2$Genus)]

Genera_order <- c("Acidovorax", "Rhizobacter", "Variovorax", "Polaromonas", "Cupriavidus", "Pelomonas")

# Count number of isolates which belong to every family
Genus_count <- data.frame(table(tax_df_2$Genus))
Genus_count$V2 <- paste(Genus_count$Var1, " (n = ", Genus_count$Freq, ")", sep ="")
# Create a named vector for mapping V1 to V2
mapping <- setNames(Genus_count$V2, Genus_count$Var1)
# Use the X_subset to reorder and subset the mapping to generate Y
gen_order_final <- mapping[Genera_order]

table_cor_4$Genus_2 <- Genus_count$V2[match(table_cor_3$Genus, Genus_count$Var1)]
table_cor_4 <- table_cor_4[table_cor_4$Genus_2 %in% gen_order_final,]
table_cor_4$Genus_2 <- factor(table_cor_4$Genus_2, levels = rev(gen_order_final))

Plot_genus <- ggplot(table_cor_4, aes(SynCom, Genus_2)) +
  geom_tile(aes(fill = R2)) +
  geom_text(aes(label = Sig), size =4) +
  scale_fill_gradient2(low = "#D55e00", mid = "white", high = "#56b4e9", midpoint =0, na.value = "lightgrey")+
  theme_classic() +
  labs(x ="SynCom", y = "Genus", fill = "R2") +
  theme(panel.background=element_blank(),panel.grid=element_blank(),axis.line.x=element_line(size=.5, colour="black"),axis.line.y=element_line(size=.5, colour="black"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black", size=7),legend.position="right",legend.text= element_text(size=10),text=element_text(family="sans", size=10))+
  theme(axis.text.x = element_text(size = 14, angle = 25,hjust=1),axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=18)) +
  theme(plot.title = element_text(hjust = 0.5))
Plot_genus

pdf(paste(results.dir, "Figure_3d_genus_heatmap.pdf", sep=""), width=6, height=3)
print(Plot_genus)
dev.off()

table_cor_5 <- table_cor_2[table_cor_2$SynCom == "RA",]
table_cor_5$Genus <- factor(table_cor_5$Genus, levels = genus_order_3)
table_cor_6 <- table_cor_5[table_cor_5$Genus %in% Genera_order,]
table_cor_6 <- table_cor_6[order(table_cor_6$R2, decreasing = TRUE), ]

plot_bar_genus <- ggplot(table_cor_6, aes(x=Genus, y=R2)) + 
  geom_bar(stat = "identity", width = 0.98) +
  coord_flip() +
  ylim(0,0.4) +
  theme(panel.background=element_blank(),panel.grid=element_blank(),axis.line.x=element_line(size=.5, colour="black"),axis.line.y=element_line(size=.5, colour="black"),axis.ticks=element_line(color="black")) +
  theme(axis.text.y = element_blank(),axis.title.x = element_text(size = 18), axis.title.y = element_blank()) +
  ylab("Relative Abundance")
plot_bar_genus

pdf(paste(results.dir, "Figure_3d_genus_bar.pdf", sep=""), width=1.5, height=2.117)
print(plot_bar_genus)
dev.off()


