library("dplyr") #Version 1.1.2

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###852 and 266 KOs in natural rhizospheres - Merging data with Lopez et al. 2023 =====
ST <- read.table(paste(working_directory,"KO_genome/KO_SSC.tsv", sep = ""),header = TRUE)

res_ara <- readRDS(paste(working_directory,"Natural_soil/ANCOMBC_output/ANCOMBC2_ara.rds", sep = ""))
res_ara <- res_ara%>%dplyr::select(taxon, lfc_TreatmentSoil, p_TreatmentSoil)
colnames(res_ara) <- c("KO", "LFC", "padj")
res_ara$LFC <- -1*res_ara$LFC
res_ara_JJ <- readRDS(paste(working_directory,"Natural_soil/ANCOMBC_output/ANCOMBC2_ara_JJ.rds", sep = ""))
res_ara_JJ <- res_ara_JJ%>%dplyr::select(taxon, lfc_TreatmentSoil, p_TreatmentSoil)
colnames(res_ara_JJ) <- c("KO", "LFC", "padj")
res_ara_JJ$LFC <- -1*res_ara_JJ$LFC

res_ara$padj <- p.adjust(res_ara$padj, method = "fdr")
res_ara_JJ$padj <- p.adjust(res_ara_JJ$padj, method = "fdr")

ancombc2_araC <- res_ara
ancombc2_araJ <- res_ara_JJ

#Since the ANCOMBC analyses on ancombc2_araC led to no significant results, it was later discarded, only ancombc2_araJ was taken
ancombc2_araC_R <- filter(ancombc2_araC, padj<0.05&LFC>0)
ancombc2_araJ_R <- filter(ancombc2_araJ, padj<0.05&LFC>0)

### KO sets
stringent_KO <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/266KO.tsv", sep =""))
general_KO <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/852KO.tsv", sep =""))

Jose_KO <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/Lopez_et_al_2023_output.tsv", sep= ""))
Jose_KO_R <- filter(Jose_KO, niche_association == "Rhizosphere")

Jose_KO_R_2 <- unique(dplyr::select(Jose_KO_R, Orthogroup_Id, molecular_function))
Jose_KO_R_2$Lopez <- "YES"
Jose_852 <- right_join(Jose_KO_R_2, general_KO, by=c("Orthogroup_Id"="KO"))
Jose_266 <- inner_join(Jose_KO_R_2, stringent_KO, by=c("Orthogroup_Id"="KO"))
AraJJ_852 <- filter(general_KO, KO%in%ancombc2_araJ_R$KO)
AraJJ_852_2 <- dplyr::select(AraJJ_852, KO)
AraJJ_852_2$Sanchez <- "YES"
Jose_852_araJJ <- left_join(Jose_852, AraJJ_852_2, by=c("Orthogroup_Id"="KO"))
Jose_852_araJJ <- Jose_852_araJJ[,c(1,2,6,3,4,5)]
write.table(Jose_852_araJJ, paste(working_directory, "Natural_soil/852_sig_in_natural_rhizospheres.tsv", sep = ""),sep = "\t", row.names = FALSE)

### stringent
AraJJ_266 <- filter(stringent_KO, KO%in%ancombc2_araJ_R$KO)
AraJJ_266_2 <- dplyr::select(AraJJ_266, KO)
AraJJ_266_2$Sanchez <- "YES"
Jose_266_araJJ <- left_join(Jose_266, AraJJ_266_2, by=c("Orthogroup_Id"="KO"))
Jose_266_araJJ <- Jose_266_araJJ[,c(1,2,6,3,4,5)]
write.table(Jose_266_araJJ, paste(working_directory, "Natural_soil/266_sig_in_natural_rhizospheres.tsv",sep = ""), sep = "\t", row.names = FALSE)

