library("reshape2") #Version 1.4.4
library("dplyr") #Version 1.1.2
library("stringr") #Version 1.5.0
library("ggpubr") #Version 0.6.0
library("rstatix") #Version 0.7.2
library("multcompView") #Version 0.1-9
library("plyr") #Version 1.8.8

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Figure S3 - Shoot weights =====
importdat <- read.table(paste(working_directory,"SSC_R3_shoot_weights.txt", sep =""), sep = "\t", header=TRUE, dec = ".")

# Manipulate dataframe so it can be processed by ggplot2 package
good_table=t(importdat)
good_table=melt(good_table)
colnames(good_table)=c("local_condition", "replicate", "mass")

good_table=mutate(good_table,plant=word(local_condition,start = 1, sep=fixed("_")))
good_table=mutate(good_table,inoculum=word(local_condition,start = 2, sep=fixed("_")))
good_table=mutate(good_table,nutrient=word(local_condition,start = 4, sep=fixed("_")))
good_table=mutate(good_table,condition=word(local_condition,start = 1, sep=fixed("_"), end = 4))

# Remove NA values because of different number of samples being tested
good_table=na.omit(good_table)
good_table$mass=as.numeric(good_table$mass)

#Reformat results from first experiment
good_table$condition <-gsub("R1_", "",good_table$condition)
good_table_R1 <- good_table[grepl(pattern = "R1_",good_table$local_condition),]
good_table_R1$condition <- paste(good_table_R1$condition, "_Low", sep ="")
good_table_R1$plant <- good_table_R1$inoculum
good_table_R1$nutrient <- "Low"
new <- str_split(good_table_R1$local_condition, pattern = "_")
good_table_R1$inoculum <- data.table::transpose(new)[[3]]

good_table_R2 <- good_table[!grepl(pattern = "R1_",good_table$local_condition),]

good_table_2 <- good_table_R2

good_table_2$inoc <- factor(good_table_2$inoc, levels = c("SSC","AtSC","LjSC","HvSC","NS"))
colnames(good_table_2)[colnames(good_table_2) == "inoculum"] <- "SynCom"

PLANTS <- c("At", "Hv", "Lj")
FULL_NAMES <- c("Arabidopsis", "Barley", "Lotus")
NUTRIENT=c("Low", "High")
Plot_list=list()

max_scale <- data.frame(c("Arabidopsis", "Barley", "Lotus"), c(max(good_table$mass[good_table$plant == "At"])+0.1*max(good_table$mass[good_table$plant == "At"]),
                                                               max(good_table$mass[good_table$plant == "Hv"])+0.1*max(good_table$mass[good_table$plant == "Hv"]),
                                                               max(good_table$mass[good_table$plant == "Lj"])+0.1*max(good_table$mass[good_table$plant == "Lj"])))
colnames(max_scale) <- c("Host","Max")

for (i in 0:2) {
  for (j in 1:2){
    
    good_table_3=subset (x = good_table_2, subset = plant==PLANTS[i+1])
    good_table_3$SynCom <- factor(good_table_3$SynCom, levels = c("SSC","AtSC", "LjSC", "HvSC","NS"))
    good_table_3$Experiment <- "R2"
    good_table_3$Experiment[grepl("R1_", good_table_3$local_condition)] <- "R1"
    
    max_value <- max_scale$Max[max_scale$Host == FULL_NAMES[i+1]]
    good_table_3=subset (x = good_table_3, subset = nutrient==NUTRIENT[j])
    
    nutrient.labs <- paste(NUTRIENT[j],"nutrient")
    names(nutrient.labs) <-NUTRIENT[j]
    
    colnames(good_table_3)[colnames(good_table_3) == "SynCom"] <- "Inoculum"
    good_table_3$Inoculum <- factor(good_table_3$Inoculum, levels = c("AtSC", "HvSC", "LjSC", "SSC", "NS"))
    
    # Visualization the data using box plots. Plot weight by groups.
    
    bxpot <- ggboxplot(
      data = good_table_3,
      x = "Inoculum", 
      y = "mass",
      combine = FALSE,
      merge = FALSE,
      legend="left",
      fill =  "Inoculum",
      color = "black",
      palette = c("#A3A500","#00B0F6","#00BF7D","#F8766D","white"),
      title = ifelse(test = j==1,yes =  paste("R2 - ", FULL_NAMES[i+1], sep = ""),no =  ""),
      font.title=c(14,"italic"),
      xlab = "Inoculum",
      ylab = "Shoot mass (g)",
      bxp.errorbar = FALSE,
      bxp.errorbar.width = 0.4,
      facet.by = NULL,
      scales="free",
      panel.labs = NULL,
      short.panel.labs = TRUE,
      linetype = "solid",
      size = NULL,
      width = 0.8,
      notch = FALSE,
      outlier.shape = NA,
      select = NULL,
      remove = NULL,
      order = NULL,
      add = "none",
      add.params = list(),
      error.plot = "pointrange",
      label = NULL,
      font.label = list(size = 11, color = "black"),
      label.select = NULL,
      repel = FALSE,
      label.rectangle = FALSE,
      ggtheme = theme_pubr()
    ) +theme(axis.text.x = element_blank(),axis.title.x = element_blank(), title = element_text(hjust = 0.5, size = 12),plot.title = element_text(hjust = 0.5, size=12))+
      facet_wrap(facets = "nutrient", labeller = labeller(nutrient = nutrient.labs)
      ) + theme_classic() +
      coord_cartesian(ylim = c(0, max_value))      
    
    
    #  Kruskal wallis test followed by Dunn post hoc test
    # Test computation 
    res.kruskal <- good_table_3 %>% kruskal_test(mass ~ Inoculum)
    res.kruskal
    
    # Dunn posthoc Pairwise comparisons
    pwc <- good_table_3 %>% 
      # group_by(growth) %>%
      dunn_test(mass ~ Inoculum, p.adjust.method = "BH")
    
    # Generating letters for kruskal+dunn pairwise comparisons
    
    # Make fake ANOVA TUKEYHSD test only to replace the non parametric p-values
    tukey_values= data.frame()
    fit=aov(data=good_table_3,mass ~ Inoculum)
    anova(fit)
    res=TukeyHSD(fit)
    
    # Replace the adjusted pvalues from Kruskal+Dunn
    res[[1]][,4]=pwc$p.adj
    Tukey.levels <- res[[1]][,4]
    Labels_pairwise <- multcompLetters(Tukey.levels)['Letters']
    Inoculum <- names(Labels_pairwise[['Letters']])
    
    boxplot.df <- ddply(good_table_3, .(Inoculum), function (x) fivenum(max_scale$Max[max_scale$Host == FULL_NAMES[i+1]]*0.95))
    
    # Create a data frame out of the factor levels and Tukey's homogenous group letters
    plot.levels <- data.frame(Inoculum, labels = Labels_pairwise[['Letters']],
                              stringsAsFactors = FALSE)
    
    # Merge it with the labels
    labels.df <- merge(plot.levels, boxplot.df, by = "Inoculum" , sort = FALSE)
    
    p1=bxpot+
      geom_text(data = labels.df, aes(x = Inoculum, y = V1, label = labels))+
      geom_jitter(position=position_jitter(0.2)) +
      scale_shape_manual(values = c(15))
    
    Plot_list[[2*i+j]]=p1+rremove("x.ticks")
    
    
  }
}

#R1 data
Plot_list_R1=list()

for (i in 0:2) {
  
  good_table_3=subset (x = good_table_R1, subset = plant==PLANTS[i+1])
  colnames(good_table_3)[colnames(good_table_3) == "inoculum"] <- "SynCom"
  good_table_3$SynCom <- factor(good_table_3$SynCom, levels = c("SSC","AtSC", "LjSC", "HvSC","NS"))
  good_table_3$Experiment <- "R1"
  
  nutrient.labs <- paste(NUTRIENT[1],"nutrient")
  names(nutrient.labs) <-NUTRIENT[1]
  max <- max_scale$Max[max_scale$Host == FULL_NAMES[i+1]]
  
  colnames(good_table_3)[colnames(good_table_3) == "SynCom"] <- "Inoculum"
  good_table_3$Inoculum <- factor(good_table_3$Inoculum, levels = c("AtSC", "HvSC", "LjSC", "SSC", "NS"))
  
  # Visualization the data using box plots. Plot weight by groups.
  
  bxpot <- ggboxplot(
    data = good_table_3,
    x = "Inoculum", 
    y = "mass",
    combine = FALSE,
    merge = FALSE,
    legend="left",
    fill =  "Inoculum",
    color = "black",
    palette = c("#A3A500","#00B0F6","#00BF7D","#F8766D","white"),
    title = paste("R1 - ",FULL_NAMES[i+1],sep = ""),
    font.title=c(14,"italic"),
    xlab = "Inoculum",
    ylab = "Shoot mass (g)",
    bxp.errorbar = FALSE,
    bxp.errorbar.width = 0.4,
    facet.by = NULL,
    scales="free",
    panel.labs = NULL,
    short.panel.labs = TRUE,
    linetype = "solid",
    size = NULL,
    width = 0.8,
    notch = FALSE,
    outlier.shape = NA,
    select = NULL,
    remove = NULL,
    order = NULL,
    add = "none",
    add.params = list(),
    error.plot = "pointrange",
    label = NULL,
    font.label = list(size = 11, color = "black"),
    label.select = NULL,
    repel = FALSE,
    label.rectangle = FALSE,
    ggtheme = theme_pubr()
  ) +theme(axis.text.x = element_blank(),axis.title.x = element_blank(), title = element_text(hjust = 0.5, size = 12),plot.title = element_text(hjust = 0.5, size=12))+
    facet_wrap(facets = "nutrient", labeller = labeller(nutrient = nutrient.labs)
    ) + theme_classic() +
    scale_y_continuous(limits = c(0, max, na.rm = TRUE))
  
  #  Kruskal wallis test followed by Dunn post hoc test
  # Test computation 
  res.kruskal <- good_table_3 %>% kruskal_test(mass ~ Inoculum)
  res.kruskal
  
  # Dunn posthoc Pairwise comparisons
  pwc <- good_table_3 %>% 
    # group_by(growth) %>%
    dunn_test(mass ~ Inoculum, p.adjust.method = "BH")
  
  # Generating letters for kruskal+dunn pairwise comparisons
  
  # Make fake ANOVA TUKEYHSD test only to replace the non parametric p-values
  tukey_values= data.frame()
  fit=aov(data=good_table_3,mass ~ Inoculum)
  anova(fit)
  res=TukeyHSD(fit)
  
  # Replace the adjusted pvalues from Kruskal+Dunn
  res[[1]][,4]=pwc$p.adj
  Tukey.levels <- res[[1]][,4]
  Labels_pairwise <- multcompLetters(Tukey.levels)['Letters']
  Inoculum <- names(Labels_pairwise[['Letters']])
  
  boxplot.df <- ddply(good_table_3, .(Inoculum), function (x) fivenum(max_scale$Max[max_scale$Host == FULL_NAMES[i+1]]*0.95))
  
  # Create a data frame out of the factor levels and Tukey's homogenous group letters
  plot.levels <- data.frame(Inoculum, labels = Labels_pairwise[['Letters']],
                            stringsAsFactors = FALSE)
  
  # Merge it with the labels
  labels.df <- merge(plot.levels, boxplot.df, by = "Inoculum" , sort = FALSE)
  
  
  p1=bxpot+
    geom_text(data = labels.df, aes(x = Inoculum, y = V1, label = labels))+
    geom_jitter(position=position_jitter(0.2)) +
    scale_shape_manual(values = c(15))
  
  
  Plot_list_R1[[i+1]]=p1+rremove("x.ticks")
  
}

plot_shoots <- ggarrange(Plot_list[[1]],Plot_list[[3]]+rremove("ylab"),Plot_list[[5]]+rremove("ylab"),Plot_list[[2]],Plot_list[[4]]+rremove("ylab"),Plot_list[[6]]+rremove("ylab"), 
                         Plot_list_R1[[1]],Plot_list_R1[[2]]+rremove("ylab"),Plot_list_R1[[3]]+rremove("ylab"), ncol = 3,nrow = 3, common.legend = T, legend = "right" )

pdf(paste(results.dir,"Figure_S3_Shoot_weights.pdf", sep=""), width=15, height=15)
print(plot_shoots)
dev.off()
