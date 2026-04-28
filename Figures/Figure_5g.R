library("ggplot2") #Version 3.4.2
library("reshape2") #Version 1.4.4
library("tidyverse") #Version 2.0.0

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Figure 5g - Barnard statistics of enriched pathways across families =====
# Load data
taxonomy <- read.table(paste(working_directory,"SSC_taxonomy_GTDB.tsv",sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)
isos_90_266 <- read.table(paste(working_directory, "Barnard_stat/isos_90perc_266kos.tsv", sep = ""))$V1
kos266 <- read.table(paste(working_directory, "Barnard_stat/barnard_90_266.tsv", sep = ""), sep="\t", header = TRUE, check.names = FALSE)
fams <- taxonomy %>% group_by(family) %>% summarise(count = n()) %>% arrange(desc(count))
map2path <- read.table(paste(working_directory,"Annotations/pathway_top.txt", sep = ""), header=F, sep="\t")
ko2path <- read.table(paste(working_directory,"Annotations/KO_to_pathway.txt", sep = ""), header=T, sep="\t")

ko2path$V3 <- map2path$V2[match(ko2path$V2, map2path$V1)]

KO_to_pathway_2 <- read.table(paste(working_directory,"Annotations/KO_to_pathway_unannotated_2.txt", sep = ""), header=F, sep="\t")
colnames(KO_to_pathway_2) <- c("KO","new_category")

for (KO in KO_to_pathway_2$KO){
  ko2path$V3[ko2path$V1 == paste(KO)] <- KO_to_pathway_2$new_category[KO_to_pathway_2$KO == paste(KO)]
}

ko2path$V4 <- map2path$V3[match(ko2path$V3, map2path$V2)]

# Create pathway dataframe
colnames(ko2path) <- c("ko","map","pathway","category_2")

# Add pathways to Barnard tables
kos266_pway <- merge(kos266, ko2path, by="ko", all.x = T)

# Rownames: 266
fams_266_top <- taxonomy[taxonomy$isolate %in% isos_90_266,] %>% group_by(family) %>% summarise(top = n())
fams_266 <- merge(fams, fams_266_top, by="family") %>% mutate(nontop = count - top) %>% mutate(label = paste(family, " (", top, "/", nontop, ")", sep=""))

# Add label for plot
kos266 <- merge(kos266, fams_266[,c('family','count','label')])

# DotPlot: 266
kos266_pval <- kos266 %>% mutate(group = ifelse(pvalue < 1, "Keep", "Remove")) %>% group_by(group) %>% mutate(adj_pvalue = p.adjust(pvalue, method="BH"))

sig_pways_266 <- kos266_pway %>% filter(ko %in% subset(kos266_pval, adj_pvalue < 0.05)$ko)
sig_pways_266 <- sig_pways_266 %>% select(ko,pathway) %>% unique() %>% arrange(pathway)

kos266_pval <- merge(kos266_pval, sig_pways_266, by="ko", all = T)
kos266_pval$pathway[kos266_pval$pathway == "Unknown"] <- "Undefined"

roworder <- (kos266 %>% select(count, label) %>% unique() %>% arrange(count))$label
kos266_pval$label <- factor(kos266_pval$label, levels = roworder)
kos266_pval <- na.omit(kos266_pval)

p <- subset(kos266_pval, adj_pvalue < 0.05) %>% ggplot() + geom_point(aes(x=ko, y=label, fill=statistic, size=abs(statistic)/2), shape=21, color="gray15") +
  scale_fill_gradient2(low = "navy", high = "darkred", mid = "white" , midpoint = 0) +
  facet_grid(~pathway, scales = "free_x", space = "free_x") +
  labs(y="Family", x="KO", fill="Z-statistic") +
  guides(size=FALSE) +
  theme_bw() %+replace% theme(axis.title = element_text(size=24),
                              axis.text.y = element_text(size=14, hjust=1),
                              axis.text.x = element_text(size=12, angle = 90),
                              legend.text = element_text(size=18),
                              legend.title = element_text(size=20),
                              strip.background = element_blank(),
                              strip.text = element_text(angle=90, hjust=0),
                              panel.background = element_rect(fill="gray95"),
                              panel.grid.major = element_line(color="white"))
p

pdf(file = paste(results.dir, "Figure_5g_Significant.Dotplot.90_266.horizontal.pdf", sep = ""),width = 16, height = 6)
p
dev.off()
