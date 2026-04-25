# SSC Community scripts 

The R scripts in this repository generate the figures, supplementary figures, supplementary tables in the manuscript: *Functional capacities drive recruitment of bacteria into plant root microbiota*. These scripts can be run on the SSC data that is downloadable under this Zenodo link: https://zenodo.org/records/15656403. By downloading the data, the full path where the data is stored can be inserted in the *working_directory* line, and everything should be compatible.

In addition, some files that are used to generate these figures and tables were generated separately and can be found in any of the other folders. Lastly, Rscripts that are part of the Shiny app belonging to this manuscript are included as well (Shiny app: https://pm-bacterial-genetics-au.shinyapps.io/SSC_community_app/ ). An overview of the scripts:

## Figures

Rscripts to generate the different panels of the five main figures in the manuscript. Some supplementary figures (S7, S11, and S28) are included here as well as the code to generate the main figures is entangled with the code for these specific supplementary figures. Similarly, the code for half of Table S4 is also included here.

## Supplemental figures

Rscripts to generate the different supplementary figures in the manuscript. The script to also generate Table S6 is also included in this folder as the code to generate Figure S15 is entangled with Table S6.

## Tables

Rscripts to generate the different supplementary tables in the manuscript. The script to generate Figure S12 is also included in this folder.

## Intermediate files

Rscripts to generate intermediate files that were used to generate figures, supplementary figures and supplementary tables. These include:
- KO intravariability across the four SynComs (Fig S1A)
- KO intravariability across 1000 simulated one-strain-per-SynComs (Fig S1B)
- Order of isolates from highest functional diversity to lowest in the four SynComs - pangenome order (Fig 1d)
- Order of isolates from highest functional diversity to lowest in 1000 simulated one-strain-per-SynComs - pangenome order (Fig 1e)
- Correlation R2 value calculation from 1000 in-silico SynComs between functional diversity and root colonization (Fig S13)
- The plant-specificity of KOs calculated by the fold change of isolates with the KO in Root vs Input (Fig 5bce)
- The fold change of Root vs Input of KOs across plants (Core) (Fig 5F, S28).
- Top 70 isolates in the root microbiome (Isolates from highest abundance to lowest that together make 70% relative abundance) (Fig 4a)
- Top 90 isolates in the root microbiome (Isolates from highest abundance to lowest that together make 90% relative abundance) (Fig 5g, S30)

## Family R2

Rscripts to calculate the effect of families and Burkholderiaceae genera on the PERMANOVA R2 effect of SynComs and Hosts (Fig 4a)

## sPLS-DA

Rscripts to do Sparse Partial Least Squares Discriminant Analysis to investigate which KOs from a bacterial family play a role in the difference between hosts or SynComs subsetted for one or multiple SynComs (Fig 4bc)

## DESeq2 

Rscripts to do DESeq2 differential abundance analyses on KOs between root vs input inoculum for the SuperSynCom dataset or the LjSC family drop out SynCom dataset.

## Niche_replacement

Rscripts to calculate the niche replacement score of families in the LjSC Family drop out experiment. These scores were calculated for a combination of multiple families and plant-specific or core KO/pathways, which indicate whether the loss of a family causes other families to under/overcompensate for the loss of the pathways of interest.

## Natural_Soil

Rscripts to calculate whether the stringent and lenient core KO selections (266 and 852 KOs respectively) are also enriched in natural rhizospheres vs natural soils. 

## Shiny_app_figures

Rscripts to generate figures for every KO that are in the Shiny app, including KO abundances, isolate abundances, family contributions to the KO, plant-specificity in fold changes (root vs input)). Also versions of Figure 5b, 5f, and S27b that take all KOs are included. Please be aware that these scripts take the *original* dataset. Scripts with the dataset in which dominating strains (nodulators and *Rhizobacter* P2_G4) were removed in-silico are not included but can easily be replicated by changing the input files that are also included under the Zenodo link.

## Shiny_app_files

Rscripts to generate files on which the Shiny app runs, and are used to generate the figures in *Shiny_app_figures*. 
 
