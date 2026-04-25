library("dplyr") #Version 1.1.2
library("tidyr") #Version 1.3.0
library("phyloseq") #Version 1.44.0
library("nloptr") #Version 2.1.1
library("ANCOMBC") #Version 2.6.0

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###ANCOM-BC results - Natural Rhizosphere vs Soil - Arabidopsis datasets =====

#Data is either taken from Lopez et al. 2023, or from an ANCOMBC analysis on the Arabidopsis metagenome data in Reijerscamp soil that is shown here below

#Function to extract KOs from eggnog files
get_KO <- function(tab){
  tab_2 <- dplyr::select(tab, c("X.query", "KEGG_ko"))
  tab_2_unnested <- tab_2 %>%
    separate_rows(KEGG_ko, sep = ",")
  tab_2_unnested <- filter(tab_2_unnested, KEGG_ko != "-")
  tab_2_unnested$KEGG_ko <- gsub("ko:", "", tab_2_unnested$KEGG_ko)
  return(tab_2_unnested)
}

#Loading Arabidopsis data - Stringlis et al. (2018)
ara_soil_1 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR6797242.emapper.annotations",sep =""), skip = 4)
ara_soil_2 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR6797243.emapper.annotations",sep =""), skip = 4)                   
ara_soil_3 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR6797244.emapper.annotations",sep =""), skip = 4) 
ara_rhizo_1 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR6797246.emapper.annotations",sep =""), skip = 4)
ara_rhizo_2 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR6797249.emapper.annotations",sep =""), skip = 4)
ara_rhizo_3 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR6797250.emapper.annotations",sep =""), skip = 4)

araS1 <- get_KO(ara_soil_1)
araS2 <- get_KO(ara_soil_2)
araS3 <- get_KO(ara_soil_3)
araR1 <- get_KO(ara_rhizo_1)
araR2 <- get_KO(ara_rhizo_2)
araR3 <- get_KO(ara_rhizo_3)

metadata <- read.delim(paste(working_directory,"Natural_soil/Metadata_metagenomes_GenusLevel.tsv", sep = ""), header =T, row.names =1, sep = "\t")
Stringlis_samples <- filter(metadata, origin=="Athal_Stringlis")

ara_RS <- rbind(araS1, araS2, araS3, araR1, araR2, araR3)
ara_RS$Sample <- gsub("_.*", "", ara_RS$X.query)
colnames(Stringlis_samples)[colnames(Stringlis_samples) == "label"] <- "Sample"
ara_RS_cov <- left_join(ara_RS, Stringlis_samples)
ara_RS_cov <- ara_RS_cov %>% dplyr::select(-X.query, -origin, -method, -source)
ara_RS_cov_sum <- ara_RS_cov %>% dplyr::count(KEGG_ko, Sample)
ara_RS_cov_sum_2 <- ara_RS_cov_sum[!grepl("#", ara_RS_cov_sum$Sample),]
feature.table <- ara_RS_cov_sum_2 %>% pivot_wider(names_from = Sample, values_from = n, values_fill = 0)

colnames(feature.table) <- c("KEGG_ko", "S1", "S2", "S3", "R1", "R2", "R3") 
feature.table <- feature.table[, c("KEGG_ko", "R1", "R2", "R3", "S1", "S2", "S3")]
row_names <- feature.table$KEGG_ko
feature.table <- feature.table[-1]
rownames(feature.table) <- row_names
otu_table <- otu_table(feature.table, taxa_are_rows = TRUE)
taxonomy_table <- NULL
sample_metadata <- data.frame(
  Sample_ID = c("R1", "R2", "R3","S1", "S2", "S3"),
  Treatment = c("Rhizosphere", "Rhizosphere", "Rhizosphere",
                "Soil", "Soil", "Soil")
)
rownames(sample_metadata) <- sample_metadata$Sample_ID
sample_data_object <- sample_data(sample_metadata)
physeq <- phyloseq(otu_table,
                   taxonomy_table,
                   sample_data_object)
SSs_ancombc <- ancombc2(physeq, fix_formula = "Treatment")
results <- SSs_ancombc$res
saveRDS(results,paste(working_directory,"Natural_soil/ANCOMBC_output/ANCOMBC2_ara.rds", sep = ""))

### Arabidopsis_Juanjo_dataset
metadata <- read.delim(paste(working_directory,"Natural_soil/Metadata_metagenomes_GenusLevel.tsv", sep = ""), header =T, row.names =1, sep = "\t")
JJ_samples <- filter(metadata, origin=="Athal_Sanchez")

eggnog_files <- list.files(paste(working_directory,"Natural_soil/preliminary_files/eggnog_Arabidopsis_extra/", sep = ""), pattern = "annotation", full.names = TRUE)

combined_data <- do.call(rbind, lapply(eggnog_files, function(file) {
  read.delim(file, skip = 4)
}))
ko_table <- get_KO(combined_data)
soil_samples <- filter(JJ_samples, source=="Soil")
rhizo_samples <- filter(JJ_samples, source=="Root")
ara_RS <- ko_table
ara_RS$Contig <- gsub("(_[^_]+)$", "", ara_RS$X.query)
ara_RS$Sample <- sub("_.*$", "", ara_RS$X.query)
ara_RS$Sample <- gsub("-","_", ara_RS$Sample)

colnames(JJ_samples)[colnames(JJ_samples) == "label"] <- "Sample"
ara_RS_cov <- left_join(ara_RS, JJ_samples)
ara_RS_cov <- ara_RS_cov %>% dplyr::select(-X.query, -Contig, -origin, -method, -source)
ara_RS_cov_sum <- ara_RS_cov %>% dplyr::count(KEGG_ko, Sample)
ara_RS_cov_sum_2 <- ara_RS_cov_sum[!grepl("#", ara_RS_cov_sum$Sample),]
feature.table <- ara_RS_cov_sum_2 %>% pivot_wider(names_from = Sample, values_from = n, values_fill = 0)
feature.table <- feature.table[,c("KEGG_ko",row.names(JJ_samples))]
colnames(feature.table) <- c("KEGG_ko", "S1", "S2", "S3", "S4", "S5", "R1", "R2", "R3", 
                             "R4", "R5", "R6", "R7", "R8", "R9", "R10", "R11")
feature.table <- feature.table[, c("KEGG_ko", "R1", "R2", "R3", 
                                   "R4", "R5", "R6", "R7", "R8", "R9", "R10", "R11",
                                   "S1", "S2", "S3", "S4", "S5")]
row_names <- feature.table$KEGG_ko
feature.table <- feature.table[-1]
rownames(feature.table) <- row_names
otu_table <- otu_table(feature.table, taxa_are_rows = TRUE)
taxonomy_table <- NULL
sample_metadata <- data.frame(
  Sample_ID = c("R1", "R2", "R3", 
                "R4", "R5", "R6", "R7", "R8", "R9", "R10", "R11", "S1", "S2", "S3", "S4", "S5"),
  Treatment = c("Rhizosphere", "Rhizosphere", "Rhizosphere", "Rhizosphere", "Rhizosphere", "Rhizosphere",
                "Rhizosphere", "Rhizosphere", "Rhizosphere", "Rhizosphere", "Rhizosphere", 
                "Soil", "Soil", "Soil", "Soil", "Soil")
)
rownames(sample_metadata) <- sample_metadata$Sample_ID
sample_data_object <- sample_data(sample_metadata)
physeq <- phyloseq(otu_table,
                   taxonomy_table,
                   sample_data_object)
SSs_ancombc <- ancombc2(physeq, fix_formula = "Treatment")
results <- SSs_ancombc$res
saveRDS(results, paste(working_directory,"Natural_soil/ANCOMBC_output/ANCOMBC2_ara_JJ.rds", sep = ""))
