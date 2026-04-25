library("dplyr") #Version 1.1.2
library("phyloseq") #Version 1.44.0
library("DESeq2") #Version 1.40.0

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###DESeq2 on the dataset without nodulators =====
SynCom <- c("AtSC", "LjSC", "HvSC", "SSC")

sigtab_col_all_2 <- data.frame(matrix(NA, ncol = 8))
colnames(sigtab_col_all_2) <- c("baseMean", "log2FoldChange","lfcSE", "pvalue", "padj", "plant", "KO","SynCom")
sigtab_col_all_2 <- sigtab_col_all_2[-1,]

for (inoculum in SynCom) {
  #Insert OTU table
  otu_mat = read.table(paste(working_directory,"KO_tables/No_nodulators/", inoculum,".tsv", sep = ""), header=T, sep="\t")
  #Insert metadata file
  samples_df = read.table(paste(working_directory,"SSC_R2_metadata.tsv", sep =""), header=TRUE,sep="\t") #make the SampleID column into the row.names
  #make the SampleID column into the row.names
  row.names(samples_df) <- samples_df$sample_id
  samples_df <- samples_df %>% dplyr::select (-sample_id)
  row.names(otu_mat) <- otu_mat$function.
  otu_mat <- otu_mat %>% dplyr::select (-function.)
  
  ##DESeq2 Analysis
  otu_mat <- as.matrix(otu_mat)
  
  OTU = otu_table(otu_mat,taxa_are_rows = TRUE)
  samples = sample_data(samples_df)
  samples$Condition <- relevel(factor(samples$Condition), ref = "Input")
  
  #import as phyloseq object
  phylo <- phyloseq(OTU,samples)
  
  phylo_1 = subset_samples(phylo, Compartment != "RZ")
  phylo_2 = subset_samples(phylo_1, Compartment != "BS")
  phylo_3 = subset_samples(phylo_2, Inoculum == paste(inoculum))
  
  #make Deseq2 object, where the column 'genotype' will be used as comparison
  DEseq2_meta_col = phyloseq_to_deseq2(phylo_3, ~Condition)
  
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
  res_col_Hv <- lfcShrink(DEseq2_meta_col, coef="Condition_Hv_vs_Input")
  res_col_At <- lfcShrink(DEseq2_meta_col, coef="Condition_At_vs_Input")
  res_col_Lj <- lfcShrink(DEseq2_meta_col, coef="Condition_Lj_vs_Input")
  
  alpha = 0.05
  sigtab_col_Hv = res_col_Hv[which(res_col_Hv$padj < alpha), ]
  sigtab_col_At = res_col_At[which(res_col_At$padj < alpha), ]
  sigtab_col_Lj = res_col_Lj[which(res_col_Lj$padj < alpha), ]
  
  sigtab_col_Hv <- sigtab_col_Hv[sigtab_col_Hv$log2FoldChange >= 0, ]
  sigtab_col_At <- sigtab_col_At[sigtab_col_At$log2FoldChange >= 0, ]
  sigtab_col_Lj <- sigtab_col_Lj[sigtab_col_Lj$log2FoldChange >= 0, ]
  
  sigtab_col_Hv$plant <- "Barley"
  sigtab_col_At$plant <- "Arabidopsis"
  sigtab_col_Lj$plant <- "Lotus"
  
  sigtab_col_Lj$KO <- row.names(sigtab_col_Lj)
  sigtab_col_At$KO <- row.names(sigtab_col_At)
  sigtab_col_Hv$KO <- row.names(sigtab_col_Hv)
  
  sigtab_col_all <- rbind(sigtab_col_Lj, sigtab_col_Hv, sigtab_col_At)
  sigtab_col_all$SynCom <- paste(inoculum)
  sigtab_col_all_2 <- rbind(sigtab_col_all_2,as.data.frame(sigtab_col_all))
}

write.table(sigtab_col_all_2 , file = paste(working_directory, "DESeq2/Sig_KO_all_no_nod.txt", sep = ""), quote=F, sep="\t", col.names=T, row.names=F)
