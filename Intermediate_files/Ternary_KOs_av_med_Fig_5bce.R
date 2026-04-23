library("dplyr") #Version 1.1.2

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Script to generate ternary_KOs_av_med.txt with dom and without dom - Necessary for 5b, c, and e - lenient and strict thresholds =====

KO_table = read.table(paste(working_directory, "KO_genome/KO_SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
colnames(KO_table) <- gsub("X", "", colnames(KO_table))

tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)

rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)

input_table <- read.table(paste(working_directory, "DESeq2/Sig_KO_all.txt", sep = ""), header=T, sep="\t")
input_table_2 <- table(input_table$KO)
input_table_3 <- names(input_table_2)[input_table_2 > 6]

input_table_At <- input_table[input_table$Plant == "Arabidopsis",]
input_table_At_2 <- table(input_table_At$KO)
input_table_At_3 <- names(input_table_At_2)[input_table_At_2 > 2]

input_table_Hv <- input_table[input_table$Plant == "Barley",]
input_table_Hv_2 <- table(input_table_Hv$KO)
input_table_Hv_3 <- names(input_table_Hv_2)[input_table_Hv_2 > 2]

input_table_Lj <- input_table[input_table$Plant == "Lotus",]
input_table_Lj_2 <- table(input_table_Lj$KO)
input_table_Lj_3 <- names(input_table_Lj_2)[input_table_Lj_2 > 2]

input_table_3 <- unique(c(input_table_3,input_table_At_3, input_table_Hv_3, input_table_Lj_3))

SynComs <- c("AtSC","LjSC", "HvSC", "SSC")
hop_2 <- data.frame()

for (KO in input_table_3){
  for (syncom in SynComs){
    syncom_table <- read.table(paste(working_directory,"Isolate_tables/Original/", syncom,"_norm.tsv", sep= ""), sep = "\t", header =T, row.names =1)
    syncom_table_2 <- syncom_table[,grep("ES", colnames(syncom_table))]
    syncom_table_3 <- syncom_table_2[,grep(paste(syncom), colnames(syncom_table_2))]
    syncom_table_4 <- syncom_table_3[,!grepl("HL", colnames(syncom_table_3))]
    
    KO_table_sub <- KO_table[row.names(KO_table) == paste(KO),]
    
    if (syncom == "SSC"){
      KO_table_sub_2 <- KO_table_sub
    } else {
      KO_table_sub_2 <- KO_table_sub[,colnames(KO_table_sub) %in% row.names(tax_df_2)[tax_df_2$SynCom == paste(syncom)]]
    }
    
    KO_table_sub_yes <- names(KO_table_sub_2)[KO_table_sub_2 > 0]
    KO_table_sub_no <- names(KO_table_sub_2)[KO_table_sub_2 == 0]
    
    syncom_table_5 <- t(t(syncom_table_4)/rowSums(t(syncom_table_4)))
    
    syncom_table_At <- syncom_table_5[,grep("At_", colnames(syncom_table_5))]
    syncom_table_Hv <- syncom_table_5[,grep("Hv_", colnames(syncom_table_5))]
    syncom_table_Lj <- syncom_table_5[,grep("Lj_", colnames(syncom_table_5))]
    
    #Averages
    syncom_table_At_2 <- rowSums(syncom_table_At)/length(colnames(syncom_table_At))
    syncom_table_Hv_2 <- rowSums(syncom_table_Hv)/length(colnames(syncom_table_Hv))
    syncom_table_Lj_2 <- rowSums(syncom_table_Lj)/length(colnames(syncom_table_Lj))
    
    At_RA <- sum(syncom_table_At_2[names(syncom_table_At_2) %in% KO_table_sub_yes])
    Hv_RA <- sum(syncom_table_Hv_2[names(syncom_table_Hv_2) %in% KO_table_sub_yes])
    Lj_RA <- sum(syncom_table_Lj_2[names(syncom_table_Lj_2) %in% KO_table_sub_yes])
    
    syncom_table_inp <- syncom_table[,grep("Input", colnames(syncom_table))]
    syncom_table_inp_2 <- t(t(syncom_table_inp)/rowSums(t(syncom_table_inp)))
    
    syncom_table_inp_3  <- rowSums(syncom_table_inp_2)/length(colnames(syncom_table_inp_2))
    
    Input_RA <- sum(syncom_table_inp_3[names(syncom_table_inp_3) %in% KO_table_sub_yes])
    
    if(Input_RA == 0){
      At_val <- At_RA
      Hv_val <- Hv_RA
      Lj_val <- Lj_RA
    } else {
      At_val <- At_RA/Input_RA
      Hv_val <- Hv_RA/Input_RA
      Lj_val <- Lj_RA/Input_RA
    }
    
    len_val <- length(KO_table_sub_yes)/length(KO_table_sub_2)
    
    hop <- t(data.frame(c(paste(KO), At_val, Hv_val, Lj_val, len_val, paste(syncom), length(colnames(syncom_table_At)), length(colnames(syncom_table_Hv)),length(colnames(syncom_table_Lj)))))
    
    hop_2 <- rbind(hop_2, hop)
  }
}

row.names(hop_2) <- NULL
colnames(hop_2) <- c("KO", "At_val", "Hv_val", "Lj_val", "No_of_strains", "SynCom", "No_of_samples_At", "No_of_samples_Hv", "No_of_samples_Lj")

hop_4 <- data.frame()

for (KO in input_table_3){
  hop_sub <- hop_2[hop_2$KO == paste(KO),]
  
  new_3 <- data.frame()
  for (syncom in unique(hop_sub$SynCom)){
    hop_sub_2 <- hop_sub[hop_sub$SynCom == paste(syncom),]
    At_val <- as.numeric(hop_sub_2$At_val)
    Hv_val <- as.numeric(hop_sub_2$Hv_val)
    Lj_val <- as.numeric(hop_sub_2$Lj_val)
    new_2 <- data.frame(At_val,Hv_val, Lj_val,hop_sub_2$No_of_samples_At,hop_sub_2$No_of_samples_Hv, hop_sub_2$No_of_samples_Lj, paste(syncom))
    new_3 <- rbind(new_3,new_2)
  }
  
  syncom_table_inp_3 = apply(syncom_table_inp_2, 1, median, na.rm=TRUE)
  
  #Medians
  At_val <- median(as.numeric(new_3$At_val))
  Hv_val <- median(as.numeric(new_3$Hv_val))
  Lj_val <- median(as.numeric(new_3$Lj_val))
  
  No_of_strains <- sum(as.numeric(hop_sub$No_of_strains))/length(hop_sub$No_of_strains)
  
  hop_3 <- t(data.frame(c(paste(KO), At_val,Hv_val, Lj_val, No_of_strains)))
  
  hop_4 <- rbind(hop_4, hop_3)
}

row.names(hop_4) <- NULL
colnames(hop_4) <- c("KO", "Arabidopsis", "Barley", "Lotus", "Proportion_of_strains")

write.table(hop_4, paste(working_directory,"Functionality/852/ternary_KOs_av_med.txt", sep = ""), quote =F, col.names =T, row.names =T, sep ="\t")

#Generating file for ternary plots
KO_table = read.table(paste(working_directory, "KO_genome/KO_SSC.tsv", sep = ""), header=TRUE,sep="\t", row.names = 1)
colnames(KO_table) <- gsub("X", "", colnames(KO_table))

tax_df = read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)

input_table <- read.table(paste(working_directory, "DESeq2/Sig_KO_all_no_nod_rhizo.txt", sep = ""), header=T, sep="\t")
input_table_2 <- table(input_table$KO)
input_table_3 <- names(input_table_2)[input_table_2 == 12]

input_table_At <- input_table[input_table$Plant == "Arabidopsis",]
input_table_At_2 <- table(input_table_At$KO)
input_table_At_3 <- names(input_table_At_2)[input_table_At_2 > 2]

input_table_Hv <- input_table[input_table$Plant == "Barley",]
input_table_Hv_2 <- table(input_table_Hv$KO)
input_table_Hv_3 <- names(input_table_Hv_2)[input_table_Hv_2 > 2]

input_table_Lj <- input_table[input_table$Plant == "Lotus",]
input_table_Lj_2 <- table(input_table_Lj$KO)
input_table_Lj_3 <- names(input_table_Lj_2)[input_table_Lj_2 > 2]

input_table_3 <- unique(c(input_table_3,input_table_At_3, input_table_Hv_3, input_table_Lj_3))

SynComs <- c("AtSC","LjSC", "HvSC", "SSC")
hop_2 <- data.frame()

for (KO in input_table_3){
  for (syncom in SynComs){
    syncom_table <- read.table(paste(working_directory, "Isolate_tables/No_dominances/", syncom,"_norm.tsv", sep= ""), sep = "\t", header =T, row.names =1)
    syncom_table_2 <- syncom_table[,grep("ES", colnames(syncom_table))]
    syncom_table_3 <- syncom_table_2[,grep(paste(syncom), colnames(syncom_table_2))]
    syncom_table_4 <- syncom_table_3[,!grepl("HL", colnames(syncom_table_3))]
    
    KO_table_sub <- KO_table[row.names(KO_table) == paste(KO),]
    
    if (syncom == "SSC"){
      KO_table_sub_2 <- KO_table_sub
    } else {
      KO_table_sub_2 <- KO_table_sub[,colnames(KO_table_sub) %in% row.names(tax_df_2)[tax_df_2$SynCom == paste(syncom)]]
    }
    
    KO_table_sub_yes <- names(KO_table_sub_2)[KO_table_sub_2 > 0]
    KO_table_sub_no <- names(KO_table_sub_2)[KO_table_sub_2 == 0]
    
    syncom_table_5 <- t(t(syncom_table_4)/rowSums(t(syncom_table_4)))
    
    syncom_table_At <- syncom_table_5[,grep("At_", colnames(syncom_table_5))]
    syncom_table_Hv <- syncom_table_5[,grep("Hv_", colnames(syncom_table_5))]
    syncom_table_Lj <- syncom_table_5[,grep("Lj_", colnames(syncom_table_5))]
    
    #Averages
    syncom_table_At_2 <- rowSums(syncom_table_At)/length(colnames(syncom_table_At))
    syncom_table_Hv_2 <- rowSums(syncom_table_Hv)/length(colnames(syncom_table_Hv))
    syncom_table_Lj_2 <- rowSums(syncom_table_Lj)/length(colnames(syncom_table_Lj))
    
    At_RA <- sum(syncom_table_At_2[names(syncom_table_At_2) %in% KO_table_sub_yes])
    Hv_RA <- sum(syncom_table_Hv_2[names(syncom_table_Hv_2) %in% KO_table_sub_yes])
    Lj_RA <- sum(syncom_table_Lj_2[names(syncom_table_Lj_2) %in% KO_table_sub_yes])
    
    syncom_table_inp <- syncom_table[,grep("Input", colnames(syncom_table))]
    syncom_table_inp_2 <- t(t(syncom_table_inp)/rowSums(t(syncom_table_inp)))
    
    syncom_table_inp_3  <- rowSums(syncom_table_inp_2)/length(colnames(syncom_table_inp_2))
    
    Input_RA <- sum(syncom_table_inp_3[names(syncom_table_inp_3) %in% KO_table_sub_yes])
    
    if(Input_RA == 0){
      At_val <- At_RA
      Hv_val <- Hv_RA
      Lj_val <- Lj_RA
    } else {
      At_val <- At_RA/Input_RA
      Hv_val <- Hv_RA/Input_RA
      Lj_val <- Lj_RA/Input_RA
    }
    
    len_val <- length(KO_table_sub_yes)/length(KO_table_sub_2)
    
    hop <- t(data.frame(c(paste(KO), At_val, Hv_val, Lj_val, len_val, paste(syncom), length(colnames(syncom_table_At)), length(colnames(syncom_table_Hv)),length(colnames(syncom_table_Lj)))))
    
    hop_2 <- rbind(hop_2, hop)
  }
}

row.names(hop_2) <- NULL
colnames(hop_2) <- c("KO", "At_val", "Hv_val", "Lj_val", "No_of_strains", "SynCom", "No_of_samples_At", "No_of_samples_Hv", "No_of_samples_Lj")

hop_4 <- data.frame()

for (KO in input_table_3){
  hop_sub <- hop_2[hop_2$KO == paste(KO),]
  
  new_3 <- data.frame()
  for (syncom in unique(hop_sub$SynCom)){
    hop_sub_2 <- hop_sub[hop_sub$SynCom == paste(syncom),]
    At_val <- as.numeric(hop_sub_2$At_val)
    Hv_val <- as.numeric(hop_sub_2$Hv_val)
    Lj_val <- as.numeric(hop_sub_2$Lj_val)
    new_2 <- data.frame(At_val,Hv_val, Lj_val,hop_sub_2$No_of_samples_At,hop_sub_2$No_of_samples_Hv, hop_sub_2$No_of_samples_Lj, paste(syncom))
    new_3 <- rbind(new_3,new_2)
  }
  
  syncom_table_inp_3 = apply(syncom_table_inp_2, 1, median, na.rm=TRUE)
  
  
  #Medians
  At_val <- median(as.numeric(new_3$At_val))
  Hv_val <- median(as.numeric(new_3$Hv_val))
  Lj_val <- median(as.numeric(new_3$Lj_val))
  
  No_of_strains <- sum(as.numeric(hop_sub$No_of_strains))/length(hop_sub$No_of_strains)
  
  hop_3 <- t(data.frame(c(paste(KO), At_val,Hv_val, Lj_val, No_of_strains)))
  
  hop_4 <- rbind(hop_4, hop_3)
}

row.names(hop_4) <- NULL
colnames(hop_4) <- c("KO", "Arabidopsis", "Barley", "Lotus", "Proportion_of_strains")

write.table(hop_4, paste(working_directory, "Functionality/266/ternary_KOs_av_med_no_dom_266.txt", sep = ""), quote =F, col.names =T, row.names =T, sep ="\t")
