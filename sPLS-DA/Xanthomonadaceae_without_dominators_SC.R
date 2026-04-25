library("phyloseq") #Version 1.44.0
library("vegan") #Version 2.6-4
library("mixOmics") #Version 6.30.0
library("dplyr") #Version 1.1.2

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###sPLS-DA - Xanthomonadaceae without dominators - SynComs =====
#otu table
KO_SSC=read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Xanthomonadaceae_KO_no_dom.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t", row.names =1) #make the SampleID column into the row.names
colnames(samples_df)[5]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object

#Sample subsetting
samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

#At subset - Enterobacteriaceae - dom vs no-dom 

OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)
samples_sub = sample_data(samples_df_sub_2)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

bray_2 <- as.matrix(beta_isolate_KO)

str(samples_df_sub_4)
str(bray_2)

#Bind metadata with distance matrix
pcoa = cmdscale(bray_2, k=10, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z", "a", "b", "c", "d", "e", "f", "g") 
eig = pcoa$eig
points_2 <- points[order(row.names(points)), ]
samples_df_sub_3 <- samples_df_sub_2[row.names(samples_df_sub_2) %in% row.names(points),]
samples_df_sub_4 <- samples_df_sub_3[order(row.names(samples_df_sub_3)), ]
points_3 <- cbind(points_2,samples_df_sub_4)
colnames(points_3) <- c("x", "y", "z","a", "b", "c", "d", "e", "f", "g",colnames(samples_df_sub_4))

# Principal Coordinates Analysis (PCoA)
pcoa_result <- cmdscale(bray_2, eig = TRUE, k = 2)  # k = number of dimensions

# Extract coordinates for plotting
pcoa_coords <- as.data.frame(pcoa_result$points)

# Adding metadata environmental data
pcoa_coords$Inoculum <- samples_df_sub_7$Inoculum[match(row.names(pcoa_coords), row.names(samples_df_sub_7))]

KO_SSC_3 <- KO_SSC[,colnames(KO_SSC) %in% row.names(pcoa_coords)]
KO_SSC_4 <- t(t(KO_SSC_3)/rowSums(t(KO_SSC_3)))

data_plot <- data.frame(matrix(NA, ncol = 5))
colnames(data_plot) <- c("KO", "contrib", "component", "Cluster", "Study")
data_plot_2 <- data_plot[-1,]

#mixOmics
KO_SSC_5 <- t(KO_SSC_4)
KO_SSC_6 <- KO_SSC_5[match(row.names(points_3),row.names(KO_SSC_5)),]
KO_SSC_7 <- KO_SSC_6[, colSums(KO_SSC_6 != 0) > 0]

#remove columns with 0's
final.plsda <- plsda(KO_SSC_7,points_3$Inoculum, ncomp = 10)

set.seed(30) # For reproducibility with this handbook, remove otherwise
perf.plsda <- perf(final.plsda, validation = 'Mfold', folds = 3, 
                   progressBar = FALSE,  # Set to TRUE to track progress
                   nrepeat = 10)     

list.keepX <- c(1:10,  seq(20, 100, 10))
list.keepX

tune.splsda_data <- tune.splsda(KO_SSC_7,points_3$Inoculum, ncomp = 4, validation = 'Mfold', 
                                folds = 5, dist = 'max.dist', 
                                test.keepX = list.keepX, nrepeat = 10)
ncomp <- tune.splsda_data$choice.ncomp$ncomp 
select.keepX <- tune.splsda_data$choice.keepX[1:ncomp]  

splsda.data <- splsda(KO_SSC_7,points_3$Inoculum, ncomp = ncomp, keepX = select.keepX) 
var.name.short <- colnames(KO_SSC_7)

empty_vector <- vector()

empty_vector_all <- data.frame(matrix(NA, ncol = 4))
colnames(empty_vector_all) <- c("KO", "contrib", "component", "Data")
empty_vector_all_2 <- empty_vector_all[-1,]

empty_vector_contrib <- data_frame()

for (comp in 1:ncomp) {
  list <- as.data.frame(splsda.data$loadings$X)
  nonredun <- row.names(list)[rowSums(list)!=0]
  nonredun_2 <- list[row.names(list) %in% nonredun,]
  nonredun_3 <- as.data.frame(cbind(nonredun, nonredun_2))
  comp_2 <- comp+1
  if (comp == 1) {
    PC1 <- nonredun_3[comp_2]
    PC1$KO <- nonredun_3$nonredun
  } else {
    PC1 <- nonredun_3[comp_2]
    PC1$KO <- row.names(PC1)
  }
  colnames(PC1) <- c(paste("comp",comp,sep=""),"KO")
  PC1_2 <- PC1[order(PC1[,paste("comp",comp,sep="")]),]
  PC1_3 <- PC1_2[PC1_2[,-1] !=0,]
  PC1_3 <- PC1_3[PC1_3[1] != 0,]
  PC1_3$extra <- abs(as.numeric(unlist(PC1_3[1])))
  
  empty_vector <- c(empty_vector,PC1_3$KO)
  
  table <- data.frame(PC1_3$KO)
  inter <- as.data.frame(do.call(cbind, PC1_3[3]))
  
  table$contrib <- unlist(inter[1])
  colnames(table) <- c("KO", "contrib")
  table$component <- paste("PC", comp, sep="")
  
  plot <- plotLoadings(splsda.data, comp = comp, method = 'mean', contrib = 'max', 
                       name.var = var.name.short)
  
  new_data <- data.frame(plot$X$GroupContrib)
  row.names(new_data) <- row.names(plot$X)
  empty_vector_contrib <- rbind(empty_vector_contrib,new_data)
  
  table$Data <- new_data$plot.X.GroupContrib[match(table$KO,row.names(new_data))]
  empty_vector_all_2 <- rbind(empty_vector_all_2, table)
}

#otu table
KO_SSC=read.table(paste(working_directory, "sPLS-DA/isolate_subset_data/Xanthomonadaceae_KO_no_dom.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)

#Samples TABLE
samples_df = read.table(paste(working_directory,"SSC_R2_metadata_no_HL.tsv", sep =""), header=TRUE,sep="\t", row.names =1) #make the SampleID column into the row.names
colnames(samples_df)[5]="Nutrient"
samples_df$Exp_Plant_compartment_inoculum_nutrient=paste(samples_df$Experiment, samples_df$Compartment, samples_df$Inoculum, samples_df$Nutrient, sep ="_")
samples_df$Plant_compartment_nutrient=paste(samples_df$Condition, samples_df$Compartment, samples_df$Nutrient, sep ="_")

#Phyloseq preparaton
#Set the OTU, TAX and sample data for making phyloseq object

#Sample subsetting
samples_df_sub <- subset(samples_df, samples_df$Compartment == "ES")
samples_df_sub_2 <- subset(samples_df_sub, samples_df_sub$Inoculum != "NS")

#At subset - Enterobacteriaceae - dom vs no-dom 
OTU_KO = otu_table(as.matrix(KO_SSC),taxa_are_rows = TRUE)
samples_sub = sample_data(samples_df_sub_2)

phylo_sub_KO = phyloseq(OTU_KO, samples_sub)
phylo_sub_KO_RA=microbiome::transform(x = phylo_sub_KO, transform = "compositional" )
beta_isolate_KO <- as.matrix(vegdist(t(phylo_sub_KO_RA@otu_table@.Data), method = "bray", diag = T))

bray_2 <- as.matrix(beta_isolate_KO)

str(samples_df_sub_2)
str(bray_2)

#Bind metadata with distance matrix
pcoa = cmdscale(bray_2, k=10, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z", "a", "b", "c", "d", "e", "f", "g") 
eig = pcoa$eig

empty_vector_all_2$Variance <- NA
empty_vector_all_2$Variance[empty_vector_all_2$component == "PC1"] <- eig[1]/sum(eig)
empty_vector_all_2$Variance[empty_vector_all_2$component == "PC2"] <- eig[2]/sum(eig)
empty_vector_all_2$Variance[empty_vector_all_2$component == "PC3"] <- eig[3]/sum(eig)
empty_vector_all_2$Variance[empty_vector_all_2$component == "PC4"] <- eig[4]/sum(eig)
empty_vector_all_2$Group <- "Xanthomonadaceae_SynCom"

write.table(empty_vector_all_2, paste(working_directory, "sPLS-DA/output/PLSDA_Xanthomonadaceae_no_dom_SynCom.tsv", sep=""), sep = "\t", quote = F, row.names = F, col.names = T)
