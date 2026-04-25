library("dplyr") #Version 1.1.2
library("tidyr") #Version 1.3.0

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Generating KO lists in natural soil and rhizosphere metagenome samples =====
#Cucumber - Ofek Lalzar et al. (2014)

#Function to extract KOs from eggnog files
get_KO <- function(tab){
  tab_2 <- dplyr::select(tab, c("X.query", "KEGG_ko"))
  tab_2_unnested <- tab_2 %>%
    separate_rows(KEGG_ko, sep = ",")
  tab_2_unnested <- filter(tab_2_unnested, KEGG_ko != "-")
  tab_2_unnested$KEGG_ko <- gsub("ko:", "", tab_2_unnested$KEGG_ko)
  return(tab_2_unnested)
}

cuc_soil_1 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR908279.emapper.annotations", sep = ""), skip = 4)
cuc_soil_2 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR908281.emapper.annotations", sep = ""), skip = 4)                   
cuc_rhizo_1 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR908208.emapper.annotations", sep = ""),skip = 4)
cuc_rhizo_2 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR908211.emapper.annotations", sep = ""), skip = 4)
cuc_rhizo_3 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR908272.emapper.annotations", sep = ""), skip = 4)

cucS1 <- get_KO(cuc_soil_1)
cucS2 <- get_KO(cuc_soil_2)
cucR1 <- get_KO(cuc_rhizo_1)
cucR2 <- get_KO(cuc_rhizo_2)
cucR3 <- get_KO(cuc_rhizo_3)

write.table(rbind(cucS1, cucS2), paste(working_directory,"Natural_soil/KO_lists/KO_list_cucumber_soil.tsv",sep =""), row.names = F, quote = F)
write.table(rbind(cucR1, cucR2, cucR3), paste(working_directory,"Natural_soil/KO_lists/KO_list_cucumber_rhizosphere.tsv",sep =""), row.names = F, quote = F)

#Wheat - Ofek Lalzar et al. (2014)
wheat_soil_1 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR908290.emapper.annotations", sep = ""), skip = 4)
wheat_soil_2 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR908291.emapper.annotations", sep = ""), skip = 4)                   
wheat_rhizo_1 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR908273.emapper.annotations", sep = ""),skip = 4)
wheat_rhizo_2 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR908275.emapper.annotations", sep = ""), skip = 4)
wheat_rhizo_3 <- read.delim(paste(working_directory,"Natural_soil/preliminary_files/eggnog_input/eggnog_SRR908276.emapper.annotations", sep = ""), skip = 4)

wheatS1 <- get_KO(wheat_soil_1)
wheatS2 <- get_KO(wheat_soil_2)
wheatR1 <- get_KO(wheat_rhizo_1)
wheatR2 <- get_KO(wheat_rhizo_2)
wheatR3 <- get_KO(wheat_rhizo_3)

write.table(rbind(wheatS1, wheatS2), paste(working_directory,"Natural_soil/KO_lists/KO_list_wheat_ofek_soil.tsv",sep =""), row.names = F, quote = F)
write.table(rbind(wheatR1, wheatR2, wheatR3), paste(working_directory,"Natural_soil/KO_lists/KO_list_wheat_ofek_rhizosphere.tsv",sep =""), row.names = F, quote = F)

# Arabidopsis - Stringlis et al. (2018)
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

# to get the list of KOs for overlap
write.table(rbind(araS1, araS2, araS3), paste(working_directory,"Natural_soil/KO_lists/KO_list_stringlis_arabidopsis_soil.tsv",sep =""), row.names = F, quote = F)
write.table(rbind(araR1, araR2, araR3), paste(working_directory,"Natural_soil/KO_lists/KO_list_stringlis_arabidopsis_rhizosphere.tsv",sep =""), row.names = F, quote = F)

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
ara_RS$niche <- sub("_.*$", "", ara_RS$X.query)
ara_RS$niche <- gsub("-","_", ara_RS$niche)

soil_ara_KO <- ara_RS %>% filter(niche%in%soil_samples$id)
soil_ara_KO <- soil_ara_KO[,1:2]
rhizo_ara_KO <- ara_RS %>% filter(niche%in%rhizo_samples$id)
rhizo_ara_KO <- rhizo_ara_KO[,1:2]
write.table(soil_ara_KO, paste(working_directory,"Natural_soil/KO_lists/KO_list_juanjo_arabidopsis_soil.tsv",sep =""), row.names = F, quote = F)
write.table(rhizo_ara_KO, paste(working_directory,"Natural_soil/KO_lists/KO_list_juanjo_arabidopsis_rhizosphere.tsv",sep =""), row.names = F, quote = F)
