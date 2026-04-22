library("dplyr") #Version 1.1.2
library("ggplot2") #Version 3.4.2
library("ggpubr") #Version 0.6.0
library("scales") #Version 1.2.1
library("ggrepel") #Version 0.9.3

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Figure S27b - Fold change curve =====
hop_4 <- read.table(paste(working_directory,"Functionality/266/boxplots.txt", sep = ""), header = T, sep = "\t", row.names =1)

colnames(hop_4) <- c("Plant","RA", "Have", "No_Have", "SynCom", "pathway", "Category_2", "No_of_strains")
hop_4$Category <- paste(hop_4$pathway, " (n% = ",format(round(as.numeric(hop_4$No_of_strains*100),1), nsmall = 1), ")", sep = "")

hop_4$Log <- log10(hop_4$Have)

hop_6 <- hop_4[,colnames(hop_4) != "No_Have"]
hop_7 <- hop_4[,colnames(hop_4) != "Have"]

colnames(hop_6) <- c("Plant","RA", "FC", "SynCom", "pathway", "Category_2","No_of_strains","Category","unused")
colnames(hop_7) <- c("Plant", "RA", "FC", "SynCom", "pathway", "Category_2","No_of_strains","Category","unused")

hop_6$Present <-"Present"
hop_7$Present <-"Absent"

hop_8 <- rbind(hop_6, hop_7)

hop_8$Present <- factor(hop_8$Present, levels = c("Present", "Absent"))
hop_9 <- na.omit(hop_8)

list_of_cats <- unique(hop_9$pathway)
new_scat_2 <- data.frame()

# Transform fold changes into LFC
hop_9 <- hop_9 %>%
  mutate(value = ifelse(FC == 0, 0, log2(FC)))

for (cat in list_of_cats){
  hop_sub <- hop_9[hop_9$pathway == paste(cat),]
  hop_sub$value <- as.numeric(hop_sub$value)
  
  stat <- compare_means(value~Present, hop_sub, method = "wilcox.test")
  
  # Calculate average Log2foldchange plant vs input of absent and present subpopulations
  Yes_LFC <- sum(hop_sub$value[hop_sub$Present == "Present"])/length(hop_sub$value[hop_sub$Present == "Present"])
  No_LFC <- sum(hop_sub$value[hop_sub$Present == "Absent"])/length(hop_sub$value[hop_sub$Present == "Absent"])
  
  Category_2 <- unique(hop_sub$Category_2)
  description <- unique(hop_sub$pathway)
  No_of_strains <- unique(hop_sub$No_of_strains)
  
  Average_RA <- sum(hop_sub$RA[hop_sub$Present == "Present"])/length(hop_sub$RA[hop_sub$Present == "Present"])
  
  new_scat <- t(data.frame(c(Category_2, description, round(as.numeric(No_of_strains),1),Average_RA, Yes_LFC, No_LFC, stat$p.format)))
  new_scat_2 <- rbind(new_scat_2, new_scat)
  
}

row.names(new_scat_2) <- NULL
colnames(new_scat_2) <- c("Category", "Pathway", "No_of_strains", "Average_RA", "LFC_vs_input_pathway_PRESENT_subpop", "LFC_vs_input_pathway_ABSENT_subpop", "Raw p-values")
new_scat_2$`Raw p-values` <-  as.numeric(new_scat_2$`Raw p-values`)

# Applying multiple adjustment methods for comparison
new_scat_2$padj_bonf = p.adjust(new_scat_2$`Raw p-values`, method = "bonferroni")  

#  Table manipulation, sorting by lowest pvals, rounding LFCs
new_scat_3 <- new_scat_2[order(new_scat_2$padj_bonf, decreasing =F),]
new_scat_3$LFC_vs_input_pathway_PRESENT_subpop <- round(as.numeric(new_scat_3$LFC_vs_input_pathway_PRESENT_subpop), 2)
new_scat_3$`LFC_vs_input_pathway_ABSENT_subpop` <- round(as.numeric(new_scat_3$`LFC_vs_input_pathway_ABSENT_subpop`), 2)

new_scat_4 <- new_scat_3[new_scat_3$LFC_vs_input_pathway_PRESENT_subpop > new_scat_3$`LFC_vs_input_pathway_ABSENT_subpop`,]
new_scat_5 <- new_scat_4[new_scat_4$padj_bonf < 0.05,]
new_scat_5$No_of_strains <- as.numeric(new_scat_5$No_of_strains)
new_scat_5$Average_RA <- as.numeric(new_scat_5$Average_RA)

Encoding(new_scat_5$Pathway) <- 'latin1'

# Create a "Square" metrics that select the best tradeoff between RA and fold change
new_scat_5$Square <- sqrt(((new_scat_5$Average_RA/max(new_scat_5$Average_RA)) * (new_scat_5$Average_RA/max(new_scat_5$Average_RA))) + ((new_scat_5$LFC_vs_input_pathway_PRESENT_subpop/max(new_scat_5$LFC_vs_input_pathway_PRESENT_subpop)) * (new_scat_5$LFC_vs_input_pathway_PRESENT_subpop/max(new_scat_5$LFC_vs_input_pathway_PRESENT_subpop))))

#  Create my other metric to best select pathways of interest
new_scat_5$diff=new_scat_5$LFC_vs_input_pathway_PRESENT_subpop-new_scat_5$LFC_vs_input_pathway_ABSENT_subpop
new_scat_5$diff_RA=new_scat_5$diff*new_scat_5$No_of_strains

# Define a legend theme, with custom text size and key size
text_size <- 12
key_size <- 0

legend_theme <- theme(legend.text = element_text(size = text_size),  # Apply text size
                      legend.key.size = unit(key_size, "in"), # Apply key size
                      strip.background=element_rect(colour="gray50", size=0.3), # Change 'size' for thickness
                      axis.text=element_text(color="gray50"),
                      axis.line = element_line(color="gray50", size=0.3),
                      axis.line.x = element_line(color="gray50", size=0.3), 
                      axis.line.y = element_line(color="gray50", size=0.3),
                      axis.ticks.x = element_line(color="gray50", size=0.3),
                      axis.ticks.y =element_line(color="gray50", size=0.3),
                      axis.text.x = element_text(size = 14), 
                      axis.title = element_text(size = 18), 
                      axis.text.y = element_text(size = 14), 
                      legend.title = element_text(size = 18), 
)

# Get 'Category' finite number of unique values
categories <- unique(new_scat_5$Category)
# Generate a color palette
color_palette <- scales::hue_pal()(length(categories))

# Create a named vector where names are categories and values are colors
named_colors <- setNames(color_palette, categories)

# Add a new column to new_scat_5 that maps each category to its color
new_scat_5$Color <- named_colors[new_scat_5$Category]

# Snippet for data point selection using exponential decay function
# Exponential decay function
A <- 3  
B <- 0.93
exp_decay <- function(x) {
  A * exp(-(log(A) / B) * x)
}

# Calculate decay values
new_scat_5$decay_value <- exp_decay(new_scat_5$Average_RA)

# Identify points above the exponential decay line
new_scat_5$above_decay <- new_scat_5$diff > new_scat_5$decay_value

# Compute distance from decay line for filtering
new_scat_5$distance = new_scat_5$diff - new_scat_5$decay_value
new_scat_5 <- new_scat_5[order(-new_scat_5$diff),]

# Select top points
top_points <- new_scat_5[new_scat_5$above_decay, ]
top_points <- top_points[order(-top_points$diff), ]

# Define non-informative pathways to remove
non_informative_pathways <- c("Biosynthesis", "Thermogenesis", "RNA polymerase", "RNA protein")  # Add more as needed

# Filter these out from the top_points
top_points <- top_points[!top_points$Pathway %in% non_informative_pathways, ]

# First, make sure 'alphabet_capital' is defined and has enough letters
alphabet_capital <- LETTERS  # Using predefined R variable LETTERS for capital alphabets
# Number of letters you want to assign
num_letters_to_assign <- 30  # Adjust this number based on your specific need
top_points$dist_category <- NA

if ("dist_category" %in% names(top_points) && "diff" %in% names(top_points) && "decay_value" %in% names(top_points)) {
  # Get indices where the condition is TRUE
  indices_to_change <- which(top_points$diff > top_points$decay_value)
  
  # Only take as many as are available
  num_to_assign <- min(num_letters_to_assign, length(indices_to_change))
  
  # Subset the indices safely
  selected_indices <- indices_to_change[1:num_to_assign]
  
  # Assign letters only to those positions
  top_points$dist_category[selected_indices] <- as.character(alphabet_capital[1:num_to_assign])
  
} else {
  cat("One or more specified columns do not exist in the dataframe top_points")
}

bottom_points <- new_scat_5[!new_scat_5$above_decay, ]
bottom_points$dist_category <- NA

top_points <- rbind(top_points, bottom_points)

plot_cor <- ggscatter(top_points, x = "Average_RA", y = "diff", color = "Category", size = "No_of_strains") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))+
  ylab("Log2Foldchange difference Present vs Absent") + 
  xlab("Relative abundance") +
  theme(legend.position = "right") +
  guides(color = guide_legend(ncol = 1)) +   
  labs(color = "Category", size = "Proportion of strains") +
  legend_theme +
  geom_text_repel(aes(label = as.character(dist_category), color= as.factor(Category)), show.legend = FALSE, size = 5, max.overlaps = Inf) 
plot_cor

pdf(paste(results.dir,"Figure_S27b_General_plot_exp_decay.pdf", sep=""), width=15, height=7)
print(plot_cor)
dev.off()
