library("dplyr") #Version 1.1.2
library("ggplot2") #Version 3.4.2
library("ggpubr") #Version 0.6.0

working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory,"results/", sep = "")

###Figure 4a - Family KO R2 effects =====
##Panel A subsetted no Burkho genera, and <0.05
#Contribution of the different bacterial families to the observed community-level split by inoculum
# Define hosts
Hosts <- c("At", "Hv", "Lj")
PL_colors <- c("Lj" = "#7570b3", "Hv" = "#d95f02",  "At"= "#1b9e77")

# Define SynComs
SynComs <- c("AtSC", "HvSC", "LjSC", "SSC")
SC_colors <- c("#F8766D", "#00BF7D", "#00B0F6", "#A3A500")


# Read data - With Dominators
fam_data <- read.table(paste0(working_directory, "Family_R2/SSC_Fam_R2_effects_with_dom.txt"),
                       sep="\t", header=TRUE, row.names=1)

combined_syncom_with_dom <- rbind(fam_data) %>%
  mutate(dataset = "Dominators")

# Read data - No Dominators
fam_data5 <- read.table(paste0(working_directory, "Family_R2/SSC_Fam_R2_effects_no_dom.txt"),
                        sep="\t", header=TRUE, row.names=1)


combined_syncom_no_dom <- rbind(fam_data5) %>%
  mutate(dataset = "No_Dominators")

# Merge datasets
pyramid_data_fam <- bind_rows(combined_syncom_no_dom, combined_syncom_with_dom)

family_order <- c("Chitinophagaceae", "Microbacteriaceae", "Micrococcaceae",
                  "Xanthobacteraceae", "Sphingobacteriaceae", "Rhodanobacteraceae",
                  "Sphingomonadaceae", "Flavobacteriaceae", "Devosiaceae",
                  "Beijerinckiaceae", "Enterobacteriaceae", "Caulobacteraceae",
                  "Pseudomonadaceae", "Xanthomonadaceae", "Burkholderiaceae",
                  "Rhizobiaceae")

# Filter for KO only and subset by hosts
pyramid_data_fam_host <- pyramid_data_fam %>% filter(Subset %in% Hosts)

# Factor levels
pyramid_data_fam_host$Subset <- factor(pyramid_data_fam_host$Subset, levels = c("Lj", "Hv", "At"))

pyramid_data_fam_host$Family <- factor(pyramid_data_fam_host$Family, levels = family_order)
pyramid_data_fam_host$R2_change <- abs(pyramid_data_fam_host$R2_change)
pyramid_data_fam_host$dataset <- factor(pyramid_data_fam_host$dataset, levels = c("Dominators", "No_Dominators"))

# Compute sum of R2_change per Family and dataset
sum_per_group_host <- aggregate(R2_change ~ Family + dataset, data = pyramid_data_fam_host, sum, na.rm = TRUE)

#Contribution of the different bacterial families to the observed community-level split by plant host

# Filter for KO only and subset by Inoculum
pyramid_data_fam_inoc <- pyramid_data_fam %>% filter(Subset %in% SynComs)

# Factor levels
pyramid_data_fam_inoc$Subset <- factor(pyramid_data_fam_inoc$Subset, levels = c("SSC", "LjSC", "HvSC", "AtSC"))
pyramid_data_fam_inoc$Family <- factor(pyramid_data_fam_inoc$Family, levels = family_order)
pyramid_data_fam_inoc$R2_change <- abs(pyramid_data_fam_inoc$R2_change)
pyramid_data_fam_inoc$dataset <- factor(pyramid_data_fam_inoc$dataset, levels = c("Dominators", "No_Dominators"))

# Compute sum of R2_change per Family and dataset
sum_per_group_inoc <- aggregate(R2_change ~ Family + dataset, data = pyramid_data_fam_inoc, sum, na.rm = TRUE)

# Get families with R2_change > 0.05 in either dataset
families_keep <- union(
  unique(sum_per_group_host$Family[sum_per_group_host$R2_change > 0.05]),
  unique(sum_per_group_inoc$Family[sum_per_group_inoc$R2_change > 0.05])
)

# Apply the same family filter to both datasets
pyramid_data_fam_host <- pyramid_data_fam_host %>%
  filter(Family %in% families_keep) %>%
  droplevels()

pyramid_data_fam_inoc <- pyramid_data_fam_inoc %>%
  filter(Family %in% families_keep) %>%
  droplevels()

# Define bar width
barwidth <- 0.3

# Calculate positions for Dominators and Non-Dominators
Dominators <- filter(pyramid_data_fam_host, dataset == "Dominators") %>%
  group_by(Family) %>%
  arrange(-Subset) 

No_Dominators <- filter(pyramid_data_fam_host, dataset == "No_Dominators") %>%
  group_by(Family) %>%
  arrange(-Subset) 

# Plot with the calculated positions
plot8 <- ggplot() +
  geom_bar(data = Dominators,
           mapping = aes(x = as.numeric(Family) + barwidth/2, y = R2_change, fill = Subset, alpha = "Dominators"),
           stat = "identity",
           position = 'stack',
           color = "black",
           size = 0.2,
           width = barwidth) +
  geom_bar(data = No_Dominators,
           mapping = aes(x = as.numeric(Family) - barwidth/2, y = R2_change, fill = Subset, color = Subset, alpha = "No_Dominators"),
           stat = "identity",
           position = 'stack',
           size = 0.2,
           width = barwidth) +
  ggtitle("Family R² Effects With vs. Without Dominators") +
  theme(plot.title = element_text(hjust = 0.5, size = 10)) +
  theme_classic() +
  labs(x = "Family", y = "Effect on Inoculum R²", fill = "Plant Host", alpha="Dataset") +
  scale_fill_manual(values = PL_colors) +
  scale_color_manual(values = PL_colors) +
  scale_alpha_manual(values = c("Dominators" = 1,"No_Dominators" = 0.3)) +
  scale_y_continuous(limits = c(0.32, 0), 
                     expand = expansion(mult = c(-0.05, 0)), 
                     trans = "reverse",
                     breaks = seq(0.3, 0, by = -0.1),  # Ensure 0.3 is included in breaks
                     labels = scales::label_number(accuracy = 0.01))+  # Format labels with two decimal places
  scale_x_continuous(breaks = 1:length(levels(pyramid_data_fam_host$Family)), labels = levels(pyramid_data_fam_host$Family), position = "top") + # Adjust the expand parameter here
  coord_flip() + # Flip the coordinates
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black", size = 7),
        legend.position = "right",
        legend.background = element_blank(),
        legend.key = element_blank(),
        text = element_text(family = "sans", size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(angle = 0, vjust = 0.5), # Adjust y-axis title position
        axis.text.y = element_text(face = "italic", size = 10, angle = 0, hjust = 0)) # Adjust y-axis text position

# Display the plot
print(plot8)

# Generate Dominators and Non-Dominators datasets
Dominators <- filter(pyramid_data_fam_inoc, dataset == "Dominators") %>%
  group_by(Family) %>%
  arrange(-Subset)

No_Dominators <- filter(pyramid_data_fam_inoc, dataset == "No_Dominators") %>%
  group_by(Family) %>%
  arrange(-Subset)

# Plot native
plot9 <- ggplot() +
  geom_bar(data = Dominators,
           mapping = aes(x = as.numeric(Family) + barwidth/2, y = R2_change, fill = Subset, alpha = "Dominators"),
           stat = "identity",
           position = 'stack',
           color = "black",
           size = 0.2,
           width = barwidth) +
  geom_bar(data = No_Dominators,
           mapping = aes(x = as.numeric(Family) - barwidth/2, y = R2_change, fill = Subset, alpha = "No_Dominators", color=Subset),
           stat = "identity",
           position = 'stack',
           size = 0.2,
           width = barwidth) +
  ggtitle("Family R² Effects With vs. Without Dominators") +
  theme(plot.title = element_text(hjust = 0.5, size = 10)) +
  theme_classic() +
  labs(x = "Family", y = "Effect on SynCom R²", fill = "Inoculum", alpha="Dataset") +
  scale_fill_manual(values = SC_colors) +
  scale_color_manual(values = SC_colors) +
  scale_alpha_manual(values = c("Dominators" = 1, "No_Dominators" = 0.3)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black", size = 7),
        legend.position = "right",
        legend.background = element_blank(),
        legend.key = element_blank(),
        text = element_text(family = "sans", size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = "italic", size = 10)) +
  scale_x_continuous(breaks = 1:length(levels(pyramid_data_fam_inoc$Family)), labels = levels(pyramid_data_fam_inoc$Family))+
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01), expand = expansion(mult = c(0, 0.05)), breaks = seq(0,0.3, 0.1)) +
  coord_flip()

# Display the plot
plot9

# Plot with Rhizobiaceae break
plot9_lim <- ggplot() +
  geom_bar(data = Dominators,
           mapping = aes(x = as.numeric(Family) + barwidth/2, y = R2_change, fill = Subset, alpha = "Dominators"),
           stat = "identity",
           position = 'stack',
           color = "black",
           size = 0.2,
           width = barwidth) +
  geom_bar(data = No_Dominators,
           mapping = aes(x = as.numeric(Family) - barwidth/2, y = R2_change, fill = Subset, alpha = "No_Dominators", color=Subset),
           stat = "identity",
           position = 'stack',
           size = 0.2,
           width = barwidth) +
  ggtitle("Family R² Effects With vs. Without Dominators") +
  theme(plot.title = element_text(hjust = 0.5, size = 10)) +
  theme_classic() +
  labs(x = "Family", y = "Effect on Host R²", fill = "Inoculum", alpha="Dataset") +
  scale_fill_manual(values = SC_colors) +
  scale_color_manual(values = SC_colors) +
  scale_alpha_manual(values = c("Dominators" = 1, "No_Dominators" = 0.3)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black", size = 7),
        legend.position = "right",
        legend.background = element_blank(),
        legend.key = element_blank(),
        text = element_text(family = "sans", size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = "italic", size = 10)) +
  scale_x_continuous(breaks = 1:length(levels(pyramid_data_fam_inoc$Family)), labels = levels(pyramid_data_fam_inoc$Family))+
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01), expand = expansion(mult = c(0, 0.05)), breaks = seq(0,0.3, 0.1), limits = c(0,0.32) ) +
  coord_flip()

plot9_lim




# Combine syncom and plant R² family effect (NATIVE)
pyramid_dom=ggarrange(print(plot8), print(plot9), 
                      ncol = 2, nrow = 1,    
                      common.legend = F,  
                      legend = "none")      
pyramid_dom

# Combine syncom and plant R² family effect with limits
pyramid_dom_lim=ggarrange(print(plot8), print(plot9_lim), 
                          ncol = 2, nrow = 1,    
                          common.legend = F,  
                          legend = "none")      
pyramid_dom_lim

# Extract legends from both plots
legend8 <- get_legend(plot8 + theme(legend.position = "right"))
legend9 <- get_legend(plot9 + theme(legend.position = "right"))

# Arrange legends in a single plot
legend_plot <- ggarrange(legend8, legend9, ncol = 2, nrow = 1)

# Display the legend plot
legend_plot

# # Save plot
pdf(paste0(results.dir, "Figure_4a_Pyramid_Plot_main_lim.pdf"), width = 10, height = 8)
print(pyramid_dom_lim)
dev.off()

pdf(paste0(results.dir, "Figure_4a_Pyramid_Plot_legends.pdf"), width = 10, height = 8)
print(legend_plot)
dev.off()
