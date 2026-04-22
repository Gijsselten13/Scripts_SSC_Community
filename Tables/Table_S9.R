library("dplyr") #Version 1.1.2

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Table S9 - KO overlap in natural rhizospheres vs natural soils =====

groups <- c("cucumber_rhizosphere","cucumber_soil","juanjo_arabidopsis_soil", "juanjo_arabidopsis_rhizosphere","stringlis_arabidopsis_rhizosphere","stringlis_arabidopsis_soil","wheat_ofek_rhizosphere","wheat_ofek_soil")

compile_dat <- data.frame()

for (group in groups){
  date_set <- read.table(paste(working_directory, "Natural_soil/KO_lists/KO_list_", group, ".tsv", sep =""), sep = " ", header =T)
  date_set$V3 <- paste(group)
  compile_dat <- rbind(compile_dat, date_set)
}

Soil_ofek <- unique(c(compile_dat$KEGG_ko[compile_dat$V3 == "cucumber_soil"], compile_dat$KEGG_ko[compile_dat$V3 == "wheat_ofek_soil"]))
Soil_jj_string <- unique(c(compile_dat$KEGG_ko[compile_dat$V3 == "juanjo_arabidopsis_soil"], compile_dat$KEGG_ko[compile_dat$V3 == "stringlis_arabidopsis_soil"]))
Rhizosphere_ofek <- unique(c(compile_dat$KEGG_ko[compile_dat$V3 == "cucumber_rhizosphere"], compile_dat$KEGG_ko[compile_dat$V3 == "wheat_ofek_rhizosphere"]))
Rhizosphere_jj_string <- unique(c(compile_dat$KEGG_ko[compile_dat$V3 == "juanjo_arabidopsis_rhizosphere"], compile_dat$KEGG_ko[compile_dat$V3 == "stringlis_arabidopsis_rhizosphere"]))

Soil_overlap <- table(c(Soil_ofek, Soil_jj_string))
Expectancy <- length(names(Soil_overlap)[Soil_overlap == 2])/length(names(Soil_overlap))

Root_overlap <- table(c(Rhizosphere_ofek, Rhizosphere_jj_string))
Root_overlap_2 <- length(names(Root_overlap)[Root_overlap == 2])

together <- c(Root_overlap_2, length(names(Root_overlap)) - Root_overlap_2)

binom_out_RA <- binom.test(together,length(names(Root_overlap_2)), Expectancy)

Enrichment <- binom_out_RA$p.value

dataset <- data.frame(c(length(Soil_ofek), length(Soil_jj_string)), 
                      c(length(names(Soil_overlap)[Soil_overlap == 2]),length(names(Soil_overlap)[Soil_overlap == 2])),
                      c(Expectancy, Expectancy),
                      c(length(Rhizosphere_ofek), length(Rhizosphere_jj_string)),
                      c(Root_overlap_2,Root_overlap_2),
                      c(Root_overlap_2/length(names(Root_overlap))),
                      c(Enrichment,Enrichment))
row.names(dataset) <- c("Maon soil", "Reijerscamp soil")
colnames(dataset) <- c("Soil KOs", "Soil KO overlap", "Soil KO overlap (%)", "Rhizosphere KOs", "Rhizosphere KO overlap", "Rhizosphere KO overlap (%)", "Enrichment root vs soil")

write.table(dataset, paste(results.dir,"Table_S9_Soil_Rhizosphere_KO_overlap.tsv", sep = ""), sep = "\t", quote =F, row.names =T, col.names =T)

###Table S9 - Stats of 266 and 852 KOs in natural rhizospheres =====
Jose_266 <- read.table(paste(working_directory, "Natural_soil/266_sig_in_natural_rhizospheres.tsv", sep = ""), header =T)
Jose_852 <- read.table(paste(working_directory, "Natural_soil/852_sig_in_natural_rhizospheres.tsv", sep = ""), header =T)

#266 KOs - overlap
input_table <- read.table(paste(working_directory,"DESeq2/Sig_KO_all_no_nod_rhizo.txt", sep = ""), header=T, sep="\t")
input_table_2 <- table(input_table$KO)
input_table_3 <- names(input_table_2)[input_table_2 == 12]

Jose_both <- rbind(Jose_266, Jose_852)
Jose_both_2 <- unique(Jose_both)
Jose_266_2 <- Jose_both_2[Jose_both_2$Orthogroup_Id %in% input_table_3,]

JJ_overlap_sig_266 <- length(na.omit(Jose_266_2$Orthogroup_Id[Jose_266_2$Sanchez == "YES"]))
Jose_overlap_sig_266 <- length(na.omit(Jose_266_2$Orthogroup_Id[Jose_266_2$Lopez == "YES"]))

#852 KO overlap
JJ_overlap_sig_852 <- length(na.omit(Jose_852$Orthogroup_Id[Jose_852$Sanchez == "YES"]))
Jose_overlap_sig_852 <- length(na.omit(Jose_852$Orthogroup_Id[Jose_852$Lopez == "YES"]))

#Enrichment tests
stat_test <- function(vector1, vector2){
  # Calculate the sizes of the two vectors
  size_vector1 <- length(vector1)
  size_vector2 <- length(vector2)
  # Calculate the overlap
  overlap <- length(intersect(vector1, vector2))  # Number of shared identifiers
  # Hypergeometric test
  p_value <- phyper(
    q = overlap - 1,            # Overlap - 1 for "probability of more extreme overlaps"
    m = size_vector1,           # Size of the first group
    n = universe_size - size_vector1,  # Size of the remaining universe
    k = size_vector2,           # Size of the second group
    lower.tail = FALSE          # Use the upper tail for significance of overlap
  )
  #p_value
  fisher_test <- fisher.test(matrix(c(overlap,
                                      size_vector1 - overlap,
                                      size_vector2 - overlap,
                                      universe_size - size_vector1 - size_vector2 + overlap),
                                    nrow = 2))
  #fisher_test$p.value
  
  # Output
  #cat("Overlap:", overlap, "\n")
  #cat("P-value:", p_value, "\n")
  return(data.frame(HG_pvalue = p_value, Fisher_pvalue = fisher_test$p.value, 
                    length_overlap = overlap))
}

#Load KO lists
ST <- read.table(paste(working_directory,"KO_genome/KO_SSC.tsv", sep = ""),header = TRUE)
stringent_KO <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/266KO.tsv", sep =""))
general_KO <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/852KO.tsv", sep =""))

### Compare to Jose KOs - 266
Jose_KO <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/Lopez_et_al_2023_output.tsv", sep= ""))
Jose_KO_R <- filter(Jose_KO, niche_association == "Rhizosphere")
Jose_KO_R_2 <- unique(dplyr::select(Jose_KO_R, Orthogroup_Id, molecular_function))
Jose_KO_R_2$Lopez <- "YES"

vector1_Jose_266 <- unique(Jose_KO_R_2$Orthogroup_Id)
vector2 <- unique(stringent_KO$KO)
KO_universe <- union(Jose_KO$Orthogroup_Id, ST$sequence)
universe_size <- length(KO_universe)
universe_size_Jose_266 <- length(KO_universe)
Jose_266_stat <- stat_test(vector1_Jose_266, vector2)

### Compare to Jose KOs - 852
vector1_Jose_852 <- unique(Jose_KO_R_2$Orthogroup_Id)
vector2 <- unique(general_KO$KO)
KO_universe <- union(Jose_KO$Orthogroup_Id, ST$sequence)
universe_size <- length(KO_universe)
universe_size_Jose_852 <- length(KO_universe)
Jose_852_stat <- stat_test(vector1_Jose_852, vector2)

### Compare to JJ KOs - 266
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

vector1_JJ_266 <- unique(ancombc2_araJ_R$KO)
vector2 <- unique(stringent_KO$KO)
KO_universe <- union(ancombc2_araJ_R$KO, ST$sequence)
universe_size <- length(KO_universe)
universe_size_JJ_266 <- length(KO_universe)
JJ_266_stat <- stat_test(vector1_JJ_266, vector2)

### Compare to JJ KOs - 852
vector1_JJ_852 <- unique(ancombc2_araJ_R$KO)
vector2 <- unique(general_KO$KO)
KO_universe <- union(ancombc2_araJ_R$KO, ST$sequence)
universe <- length(KO_universe)
universe_size_JJ_852 <- length(KO_universe)
JJ_852_stat <- stat_test(vector1_JJ_852, vector2)

dataset <- data.frame(c("Lopez et al. (2023)", "Lopez et al. (2023)", "Sanchez et al. unpublished", "Sanchez et al. unpublished"),
                      c("266", "852", "266", "852"),
                      c(length(vector1_Jose_266), length(vector1_Jose_852), length(vector1_JJ_266), length(vector1_JJ_852)), 
                      c(universe_size_Jose_266,universe_size_Jose_852, universe_size_JJ_266,universe_size_JJ_852),
                      c(Jose_overlap_sig_266,Jose_overlap_sig_852,JJ_overlap_sig_266, JJ_overlap_sig_852),
                      c(Jose_266_stat$HG_pvalue,Jose_852_stat$HG_pvalue,JJ_266_stat$HG_pvalue,JJ_852_stat$HG_pvalue),
                      c(Jose_266_stat$Fisher_pvalue,Jose_852_stat$Fisher_pvalue,JJ_266_stat$Fisher_pvalue,JJ_852_stat$Fisher_pvalue))

colnames(dataset) <- c("Data", "KO set", "Number of Rhizosphere KOs", "Union of KO set and Rhizosphere KOs","Overlap KOs", "Hypergeometric test p-value", "Fisher test p-value")

write.table(dataset, paste(results.dir,"Table_S9_266_852_KOs_in_Natural_soil.tsv", sep = ""), sep = "\t", quote =F, row.names =T, col.names =T)
