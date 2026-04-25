library("dplyr") #Version 1.1.2
library("phyloseq") #Version 1.44.0
library("DESeq2") #Version 1.40.0

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###DESeq2 on the LjSC Family drop out experiment =====
KO_table <- read.table(paste(working_directory, "LjSC_Family_drop_out_experiment/KO_LjSC_Family_drop_out.tsv", sep = ""), header =T, row.names =1)
KO_table_2 <- KO_table[,!grepl("NOD", colnames(KO_table))]
KO_table_3 <- KO_table_2[, !grepl("LDT7", colnames(KO_table_2))]

#Metadata creation
metadata <- data.frame(colnames(KO_table_3))
colnames(metadata) <- "Sample"
metadata$Compartment <- "Root"
metadata$Compartment[grep("INPUT", metadata$Sample)] <- "Input"

metadata$Subset <- "Full LjSC"
metadata$Subset[grep("LDT1", metadata$Sample)] <- "Burkholderiaceae drop out"
metadata$Subset[grep("LDT2", metadata$Sample)] <- "Caulobacteraceae drop out"
metadata$Subset[grep("LDT3", metadata$Sample)] <- "Pseudomonadaceae drop out"
metadata$Subset[grep("LDT4", metadata$Sample)] <- "Rhizobiaceae drop out"
metadata$Subset[grep("LDT5", metadata$Sample)] <- "All other families drop out"

#DESeq2
sigtab_col_all_2 <- data.frame(matrix(NA, ncol = 7))
colnames(sigtab_col_all_2) <- c("baseMean", "log2FoldChange","lfcSE", "pvalue", "padj", "KO","Subset")
sigtab_col_all_2 <- sigtab_col_all_2[-1,]

for (inoculum in unique(metadata$Subset)) {
  KO_table_2_sub <- KO_table_3[, colnames(KO_table_3) %in% metadata$Sample[metadata$Subset == paste(inoculum)]]
  
  otu_mat <- as.matrix(KO_table_2_sub)
  
  OTU = otu_table(otu_mat,taxa_are_rows = TRUE)
  metadata_2 <- metadata[metadata$Subset == paste(inoculum),]
  row.names(metadata_2) <- metadata_2$Sample
  samples = sample_data(metadata_2)
  samples$Compartment <- relevel(factor(samples$Compartment), ref = "Input")
  
  #import as phyloseq object
  phylo <- phyloseq(OTU,samples)
  
  #make Deseq2 object, where the column 'genotype' will be used as comparison
  DEseq2_meta_col = phyloseq_to_deseq2(phylo, ~Compartment)
  
  #work around for error: every gene contains at least one zero, cannot compute log geometric means
  # calculate geometric means prior to estimate size factors
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geoMeans = apply(counts(DEseq2_meta_col), 1, gm_mean)
  
  DEseq2_meta_col = estimateSizeFactors(DEseq2_meta_col, geoMeans = geoMeans)
  
  DEseq2_meta_col = DESeq(DEseq2_meta_col, fitType="local")
  
  #check model fitting of dispersion
  dispersion_plot_col <- plotDispEsts(DEseq2_meta_col)
  
  #investigate DESeq2 results
  resultsNames(DEseq2_meta_col)
  
  #Here is the input for each comparison between mutant vs WT
  res_col <- lfcShrink(DEseq2_meta_col, coef="Compartment_Root_vs_Input")
  
  res_col$KO <- row.names(res_col)
  
  res_col$Subset <- paste(inoculum)
  sigtab_col_all_2 <- rbind(sigtab_col_all_2,as.data.frame(res_col))
}

write.table(sigtab_col_all_2 , file = paste(working_directory, "LjSC_Family_drop_out_experiment/DESeq2_Root_vs_input_Fam_drop.txt", sep = ""), quote=F, sep="\t", col.names=T, row.names=F)
