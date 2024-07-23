rm(list=ls())

library(gplots)
library(ggplot2)
library(VennDiagram)
library(viridis)
library(readxl)
library(conflicted)
library(dplyr)
library(tidyverse) # contains ggplot for volcano plots
library(RColorBrewer) # to color plots
library(ggrepel) # annotating volcano plot
library(ggfortify) # PCA plots
library(factoextra) # PCA plot ellipses
library(MSnSet.utils) # Different PCA plot option (simpler?)
library('corrr') # PCA plot 3rd option
library(ggcorrplot) # PCA plot 3rd option
library("FactoMineR") # PCA plot 3rd option
library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot)
library("xlsx")

conflicts_prefer(base::setdiff, base::intersect, dplyr::filter())


setwd("C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics")

# Loading full ECF and WM Proteomics Datasets ----

ECF_proteomics_data <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/ECF_dataset_full.xlsx"))

WM_proteomics_data <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/WM_dataset_full.xlsx"))

-----

# Loading 2-group comparison DEP datasets ----

#Significant WM and ECF protein datasets filtered in Excel 

#WT vs mdx WM

mdxvsWT_baseline_sig_proteins_WM <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/Data files for R/mdxB_WTB_WM_sig_proteins_12_7_23.xlsx"))

mdxvsWT_scruff_sig_proteins_WM <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/Data files for R/mdxS_WTS_WM_sig_proteins_12_7_23.xlsx"))

#WT vs mdx ECF

mdxvsWT_baseline_sig_proteins_ECF <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/Data files for R/mdxB_WTB_ECF_sig_proteins_12_7_23.xlsx"))

mdxvsWT_scruff_sig_proteins_ECF <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/Data files for R/mdxS_WTS_ECF_sig_proteins_12_7_23.xlsx"))

#MTBD vs mdx WM

mdxvsMTBD_baseline_sig_proteins_WM <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/Data files for R/mdxB_MTBD_B_WM_sig_proteins_12_7_23.xlsx"))

mdxvsMTBD_scruff_sig_proteins_WM <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/Data files for R/mdxS_MTBD_S_WM_sig_proteins_12_7_23.xlsx"))

#MTBD vs mdx ECF

mdxvsMTBD_baseline_sig_proteins_ECF <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/Data files for R/mdxB_MTBD_B_ECF_sig_proteins_12_7_23.xlsx"))

mdxvsMTBD_scruff_sig_proteins_ECF <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/Data files for R/mdxS_MTBD_S_ECF_sig_proteins_12_7_23.xlsx"))

#MTBD vs WT ECF

MTBDvsWT_baseline_sig_proteins_ECF <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/Data files for R/MTBD_B_WTB_ECF_sig_proteins_12_14_23.xlsx"))

MTBDvsWT_scruff_sig_proteins_ECF <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/Data files for R/MTBD_S_WTS_ECF_sig_proteins_12_19_23.xlsx"))

----


##########################
# Venn Diagrams

#library(RColorBrewer)
#myCol <- brewer.pal(3, "Pastel2")

display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}
# Suppresses data file generation each time a Venn diagram is generated
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

#mdx vs WT WM Venn Diagram ----

mdxvsWT_WM_venn <- display_venn(mdxvsWT_WM_venn_list <- list(
  'mdx vs WT baseline WM' = mdxvsWT_baseline_sig_proteins_WM$'Accession',
  'mdx vs WT post-scruff WM' = mdxvsWT_scruff_sig_proteins_WM$'Accession'), 
  # Circles
  lwd = 2,
  lty = 'blank',
  category.names = c("DE proteins mdx vs WT baseline", 
                     "DE proteins mdx vs WT post-scruff"),
  fill = c("#4f826e", "#A95AA1"),
  # Numbers
  cex = 1.2,
  fontface = "italic",
  # Set names
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-75, 55),
  cat.dist = c(-.07, -.07),
  main = "Quadriceps whole muscle",
  main.fontface = "bold"
) ----

# mdx vs WT ECF Venn Diagram ----

mdxvsWT_ECF_venn <- display_venn(mdxvsWT_ECF_venn_list <- list(
  'mdx vs WT baseline ECF' = mdxvsWT_baseline_sig_proteins_ECF$'Accession',
  'mdx vs WT post-scruff ECF' = mdxvsWT_scruff_sig_proteins_ECF$'Accession'), 
  # Circles
  lwd = 2,
  lty = 'blank',
  category.names = c("DE proteins mdx vs WT baseline", 
                     "DE proteins mdx vs WT post-scruff"),
  fill = c("#f2cd99", "#724f82"),
  # Numbers
  cex = 1.2,
  fontface = "italic",
  # Set names
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-50, 60),
  cat.dist = c(-.32, -.32),
  main = "Quadriceps ECF",
  main.fontface = "bold",
  rotation.degree = 180
) ----

# Baseline mdx vs WT ECF vs WM Venn Diagram ----

mdxvsWT_ECFvsWM_baseline_venn <- display_venn(mdxvsWT_ECFvsWM_baseline_venn_list <- list(
  'mdx vs WT baseline ECF' = mdxvsWT_baseline_sig_proteins_ECF$'Accession',
  'mdx vs WT baseline WM' = mdxvsWT_baseline_sig_proteins_WM$'Accession'), 
  # Circles
  lwd = 2,
  lty = 'blank',
  category.names = c("ECF mdx vs WT DEP", 
                     "WM mdx vs WT DEP"),
  fill = c("#f4adad", "#e34a33"),
  # Numbers
  cex = 1,
  fontface = "italic",
  # Set names
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-30, 20),
  cat.dist = c(-.28, -.3),
  main = "Quadriceps Extracellular fluid vs Whole muscle Baseline",
  main.fontface = "bold",
  rotation.degree = 180
) ----

# Post-scruff mdx vs WT ECF vs WM Venn Diagram ----

mdxvsWT_ECFvsWM_scruff_venn <- display_venn(mdxvsWT_ECFvsWM_scruff_venn_list <- list(
  'mdx vs WT post-scruff ECF' = mdxvsWT_scruff_sig_proteins_ECF$'Accession',
  'mdx vs WT post-scruff WM' = mdxvsWT_scruff_sig_proteins_WM$'Accession'), 
  # Circles
  lwd = 2,
  lty = 'blank',
  category.names = c("ECF mdx vs WT DEP", 
                     "WM mdx vs WT DEP"),
  fill = c("#BEF7FF", "#458CFF"),
  # Numbers
  cex = 0.9,
  fontface = "italic",
  # Set names
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-20, 8),
  cat.dist = c(-.1, -.06),
  main = "Quadriceps Extracellular fluid vs Whole muscle Post-Scruff",
  main.fontface = "bold"
) ----
  
# mdx vs MTBD ECF Venn Diagram ----

mdxvsMTBD_ECF_venn <- display_venn(mdxvsMTBD_ECF_venn_list <- list(
  'mdx vs MTBD baseline ECF' = mdxvsMTBD_baseline_sig_proteins_ECF$'Accession',
  'mdx vs MTBD post-scruff ECF' = mdxvsMTBD_scruff_sig_proteins_ECF$'Accession'), 
  # Circles
  lwd = 2,
  lty = 'blank',
  category.names = c("DE proteins mdx vs MTBD baseline", 
                     "DE proteins mdx vs MTBD post-scruff"),
  fill = c("#92A072", "#8072A0"),
  # Numbers
  cex = 1.2,
  fontface = "italic",
  # Set names
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-55, 55),
  cat.dist = c(-.08, -.075),
  main = "Quadriceps ECF",
  main.fontface = "bold"
) ----
  
# mdx vs MTBD WM Venn Diagram ----

mdxvsMTBD_WM_venn <- display_venn(mdxvsMTBD_WM_venn_list <- list(
  'mdx vs MTBD baseline WM' = mdxvsMTBD_baseline_sig_proteins_WM$'Accession',
  'mdx vs MTBD post-scruff WM' = mdxvsMTBD_scruff_sig_proteins_WM$'Accession'), 
  # Circles
  lwd = 2,
  lty = 'blank',
  category.names = c("DE proteins mdx vs MTBD baseline", 
                     "DE proteins mdx vs MTBD post-scruff"),
  fill = c("#20D3D1", "#D32022"),
  # Numbers
  cex = 1.2,
  fontface = "italic",
  # Set names
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-55, 15),
  cat.dist = c(-.08, -.065),
  main = "Quadriceps whole muscle",
  main.fontface = "bold"
) ----
  
# Baseline mdx vs MTBD ECF vs WM Venn Diagram ----

mdxvsMTBD_ECFvsWM_baseline_venn <- display_venn(mdxvsMTBD_ECFvsWM_baseline_venn_list <- list(
  'mdx vs MTBD baseline ECF' = mdxvsMTBD_baseline_sig_proteins_ECF$'Accession',
  'mdx vs MTBD baseline WM' = mdxvsMTBD_baseline_sig_proteins_WM$'Accession'), 
  # Circles
  lwd = 2,
  lty = 'blank',
  category.names = c("ECF mdx vs MTBD DEP", 
                     "WM mdx vs MTBD DEP"),
  fill = c("#DCD19A", "#9AA5DC"),
  # Numbers
  cex = 1,
  fontface = "italic",
  # Set names
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(130, 50),
  cat.dist = c(-.28, -.1),
  main = "Quadriceps Extracellular fluid vs Whole muscle Baseline",
  main.fontface = "bold"
) ----
# Post-scruff mdx vs MTBD ECF vs WM Venn Diagram ----

mdxvsMTBD_ECFvsWM_scruff_venn <- display_venn(mdxvsMTBD_ECFvsWM_scruff_venn_list <- list(
  'mdx vs MTBD post-scruff ECF' = mdxvsMTBD_scruff_sig_proteins_ECF$'Accession',
  'mdx vs MTBD post-scruff WM' = mdxvsMTBD_scruff_sig_proteins_WM$'Accession'), 
  # Circles
  lwd = 2,
  lty = 'blank',
  category.names = c("ECF mdx vs MTBD DEP", 
                     "WM mdx vs MTBD DEP"),
  fill = c("#9ADCC5", "#DC9AB1"),
  # Numbers
  cex = 0.9,
  fontface = "italic",
  # Set names
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-20, 8),
  cat.dist = c(-.1, -.06),
  main = "Quadriceps Extracellular fluid vs Whole muscle Post-Scruff",
  main.fontface = "bold"
) ----
# ECF Baseline mdx upreg proteins vs WT or MTBD Venn Diagram ----

ECF_mdxvsWT_or_MTBD_baseline_venn <- display_venn(ECF_mdxvsWT_or_MTBD_baseline_venn_list <- list(
  'mdx vs WT baseline' = mdxvsWT_baseline_sig_proteins_ECF$'Accession',
  'mdx vs MTBD baseline' = mdxvsMTBD_baseline_sig_proteins_ECF$'Accession'), 
  # Circles
  lwd = 2,
  lty = 'blank',
  category.names = c("mdx vs WT DEP", 
                     "mdx vs MTBD DEP"),
  fill = c("#9ecad9","#cad99e"),
  # Numbers
  cex = 1.5,
  fontface = "italic",
  # Set names
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-65, 70),
  cat.dist = c(-.295, -.345),
  main = "Quadriceps ECF Baseline",
  main.fontface = "bold",
  rotation.degree = 180
) ----

# ECF Post-scruff mdx upreg proteins vs WT or MTBD Venn Diagram ----

ECF_mdxvsWT_or_MTBD_scruff_venn <- display_venn(ECF_mdxvsWT_or_MTBD_scruff_venn_list <- list(
  'mdx vs WT post-scruff' = mdxvsWT_scruff_sig_proteins_ECF$'Accession',
  'mdx vs MTBD post-scruff' = mdxvsMTBD_scruff_sig_proteins_ECF$'Accession'), 
  # Circles
  lwd = 2,
  lty = 'blank',
  category.names = c("mdx vs WT DEP", 
                     "mdx vs MTBD DEP"),
  fill = c("#da9be2", "#a3e29b"),
  # Numbers
  cex = 1.2,
  fontface = "italic",
  # Set names
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-80, 80),
  cat.dist = c(0.04, .05),
  main = "Quadriceps ECF Post-Scruff",
  main.fontface = "bold",
  rotation.degree = 180
) ----
  
# WM Baseline mdx upreg proteins vs WT or MTBD Venn Diagram ----

WM_mdxvsWT_or_MTBD_baseline_venn <- display_venn(WM_mdxvsWT_or_MTBD_baseline_venn_list <- list(
  'mdx vs WT baseline' = mdxvsWT_baseline_sig_proteins_WM$'Accession',
  'mdx vs MTBD baseline' = mdxvsMTBD_baseline_sig_proteins_WM$'Accession'), 
  # Circles
  lwd = 2,
  lty = 'blank',
  category.names = c("mdx vs MTBD/mdx", 
                     "mdx vs WT"),
  ext.text = FALSE,
  fill = c("#EEC584","#997070"),
  # Numbers
  cex = 1.5,
  fontface = "italic",
  # Set names
  cat.cex = 1.4,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-65, 67),
  cat.dist = c(-.34, -.36),
  #main = "Quadriceps Whole muscle Baseline",
  main.fontface = "bold",
  rotation.degree=180
) 
----
# WM Post-scruff mdx upreg proteins vs WT or MTBD Venn Diagram ----

WM_mdxvsWT_or_MTBD_scruff_venn <- display_venn(WM_mdxvsWT_or_MTBD_scruff_venn_list <- list(
  'mdx vs WT post-scruff' = mdxvsWT_scruff_sig_proteins_WM$'Accession',
  'mdx vs MTBD post-scruff' = mdxvsMTBD_scruff_sig_proteins_WM$'Accession'), 
  # Circles
  lwd = 2,
  lty = 'blank',
  category.names = c("mdx vs WT DEP", 
                     "mdx vs MTBD DEP"),
  fill = c("#A05F7C", "#5FA083"),
  # Numbers
  cex = 1.2,
  fontface = "italic",
  # Set names
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-75, 75),
  cat.dist = c(-.045, -.05),
  main = "Quadriceps Whole muscle Post-Scruff",
  main.fontface = "bold"
) ----
# ECF 4-D comparisons WT, mdx, MTBD ----

ECF_4D_venn <- display_venn(ECF_4D_venn_list <- list(
  'mdx vs WT baseline' = mdxvsWT_baseline_sig_proteins_ECF$'Accession',
  'mdx vs MTBD baseline' = mdxvsMTBD_baseline_sig_proteins_ECF$'Accession',
  'mdx vs WT post-scruff' = mdxvsWT_scruff_sig_proteins_ECF$'Accession',
  'mdx vs MTBD post-scruff' = mdxvsMTBD_scruff_sig_proteins_ECF$'Accession'), 
  # Circles
  lwd = 2,
  lty = 'blank',
  category.names = c(
    "mdx vs WT baseline DEP",
    "mdx vs MTBD baseline DEP",
    "mdx vs WT post-scruff DEP", 
    "mdx vs MTBD post-scruff DEP"),
  fill = c("#9475A7","#A77B75","#75A1A7", "#88A775"),
  # Numbers
  cex = 1.2,
  fontface = "italic",
  # Set names
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-14.5, 2, -20, -20),
  cat.dist = c(0.05, .13, 0.05, 0.05),
  main = "Quadriceps ECF Comparisons WT, mdx, and MTBD",
  main.fontface = "bold"
) 

----
  
# WM 4-D comparisons WT, mdx, MTBD ----
    
WM_4D_venn <- display_venn(WM_4D_venn_list <- list(
  'mdx vs WT baseline' = mdxvsWT_baseline_sig_proteins_WM$'Accession',
  'mdx vs MTBD baseline' = mdxvsMTBD_baseline_sig_proteins_WM$'Accession',
  'mdx vs WT post-scruff' = mdxvsWT_scruff_sig_proteins_WM$'Accession',
  'mdx vs MTBD post-scruff' = mdxvsMTBD_scruff_sig_proteins_WM$'Accession'), 
  # Circles
  lwd = 2,
  lty = 'blank',
  category.names = c(
    "mdx vs WT baseline DEP",
    "mdx vs MTBD baseline DEP",
    "mdx vs WT post-scruff DEP", 
    "mdx vs MTBD post-scruff DEP"),
  fill = c("#a39c89","#89A38F","#8990A3", "#A3899D"),
  # Numbers
  cex = 1.2,
  fontface = "italic",
  # Set names
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-15, 2, -20, -20),
  cat.dist = c(0.05, .13, 0.05, 0.05),
  main = "Quadriceps Whole muscle Comparisons WT, mdx, and MTBD",
  main.fontface = "bold"
) 
-----
  
##########################################################
# Filtering Datasets by Venn Diagram overlaps/non-overlaps
##########################################################

# WM vs ECF

  # Baseline

WM_mdxvsWT_baseline_ECFvsWM_unique <- 
  data.frame('WM_mdxvsWT_baseline_ECFvsWM_unique' = setdiff
             (mdxvsWT_baseline_sig_proteins_WM$'Accession',
               mdxvsWT_baseline_sig_proteins_ECF$'Accession')) %>% 
  'colnames<-' ("Accession")

ECF_mdxvsWT_baseline_ECFvsWM_unique <- 
  data.frame('ECF_mdxvsWT_baseline_ECFvsWM_unique' = setdiff
             (mdxvsWT_baseline_sig_proteins_ECF$'Accession',
               mdxvsWT_baseline_sig_proteins_WM$'Accession')) %>% 
  'colnames<-' ("Accession")

  # Post-scruff

WM_mdxvsWT_scruff_ECFvsWM_unique <- 
  data.frame('WM_mdxvsWT_scruff_ECFvsWM_unique' = setdiff
             (mdxvsWT_scruff_sig_proteins_WM$'Accession',
               mdxvsWT_scruff_sig_proteins_ECF$'Accession')) %>% 
  'colnames<-' ("Accession")

ECF_mdxvsWT_scruff_ECFvsWM_unique <- 
  data.frame('ECF_mdxvsWT_scruff_ECFvsWM_unique' = setdiff
             (mdxvsWT_scruff_sig_proteins_ECF$'Accession',
               mdxvsWT_scruff_sig_proteins_WM$'Accession')) %>% 
  'colnames<-' ("Accession")

# WM

# Proteins sig in mdx vs MTBD baseline, excluding mdx vs WT DEPs

mdxvsMTBD_WM_baseline_exclude_WT_baseline <- 
  data.frame('mdxvsMTBD_WM_baseline_exclude_WT_baseline' = setdiff
             (mdxvsMTBD_baseline_sig_proteins_WM$'Accession',
               mdxvsWT_baseline_sig_proteins_WM$'Accession')) %>% 
  'colnames<-' ("Accession")

mdxvsMTBD_WM_baseline_exclude_all_WT <- 
  data.frame('mdxvsMTBD_WM_baseline_exclude_all_WT' = setdiff
             (mdxvsMTBD_WM_baseline_exclude_WT_baseline$'Accession',
               mdxvsWT_scruff_sig_proteins_WM$'Accession')) %>% 
  'colnames<-' ("Accession")

mdxvsMTBD_WM_baseline_compensatory <- WM_proteomics_data %>%
  filter(Accession %in% mdxvsMTBD_WM_baseline_exclude_all_WT$Accession)

write.xlsx(mdxvsMTBD_WM_baseline_compensatory, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/mdxvsMTBD_WM_baseline_compensatory.xlsx")

# Proteins sig in mdx WM vs WT at baseline, not post-scruff

mdxvsWT_WM_baseline_unique <- 
  data.frame('unique_mdxvsWT_WM_baseline' = setdiff
             (mdxvsWT_baseline_sig_proteins_WM$'Accession',
               mdxvsWT_scruff_sig_proteins_WM$'Accession')) %>% 
  'colnames<-' ("Accession")

# Proteins sig in mdx WM vs WT after scruff, not baseline

mdxvsWT_WM_scruff_unique <- 
  data.frame('unique_mdxvsWT_WM_scruff' = base::setdiff
             (mdxvsWT_scruff_sig_proteins_WM$'Accession',
               mdxvsWT_baseline_sig_proteins_WM$'Accession')) %>% 
  'colnames<-' ("Accession")

# Overlapping sig mdx vs WT WM proteins baseline and post-scruff

mdxvsWT_WM_scruff_overlap <- 
  data.frame('overlap_mdxvsWT_WM' = intersect
             (mdxvsWT_baseline_sig_proteins_WM$'Accession',
               mdxvsWT_scruff_sig_proteins_WM$'Accession')) %>% 
  'colnames<-' ("Accession")

# Filtering main datasets by unique proteins only

mdxvsWT_baseline_filtered_sig_WM <- mdxvsWT_baseline_WM %>%
  filter(Accession %in% mdxvsWT_WM_baseline_unique$Accession)

mdxvsWT_scruff_filtered_sig_WM <- mdxvsWT_scruff_WM %>%
  filter(Accession %in% mdxvsWT_WM_scruff_unique$Accession)

## Filtering ECF ##

# Proteins sig in mdx ECF vs WT at baseline, not post-scruff

mdxvsWT_ECF_baseline_unique <- 
  data.frame('unique_mdxvsWT_ECF_baseline' = setdiff
             (mdxvsWT_baseline_sig_proteins_ECF$'Accession',
               mdxvsWT_scruff_sig_proteins_ECF$'Accession')) %>% 
  'colnames<-' ("Accession")

# Proteins sig in mdx ECF vs WT after scruff, not baseline

mdxvsWT_ECF_scruff_unique <- 
  data.frame('unique_mdxvsWT_ECF_scruff' = base::setdiff
             (mdxvsWT_scruff_sig_proteins_ECF$'Accession',
               mdxvsWT_baseline_sig_proteins_ECF$'Accession')) %>% 
  'colnames<-' ("Accession")

# Overlapping sig mdx vs WT ECF proteins at baseline and post-scruff

mdxvsWT_ECF_scruff_overlap <- 
  data.frame('overlap_mdxvsWT_ECF' = intersect
             (mdxvsWT_baseline_sig_proteins_ECF$'Accession',
               mdxvsWT_scruff_sig_proteins_ECF$'Accession')) %>% 
  'colnames<-' ("Accession")

# Filtering main datasets by unique proteins only

mdxvsWT_baseline_filtered_sig_ECF <- mdxvsWT_baseline_ECF %>%
  filter(Accession %in% mdxvsWT_ECF_baseline_unique$Accession)

mdxvsWT_scruff_filtered_sig_ECF <- mdxvsWT_scruff_ECF %>%
  filter(Accession %in% mdxvsWT_ECF_scruff_unique$Accession)

###################################################################################
# Filtering ECF and WM datasets to obtain only Venn non-overlaps

# Filtering datasets by sig in ECF mdx vs WT baseline and not other conditions ----

mdxvsWT_not_MTBD_ECF_baseline_unique <- 
  data.frame('unique_mdxvsWT_not_MTBD_baseline_ECF' = base::setdiff
             (mdxvsWT_baseline_sig_proteins_ECF$'Accession',
               mdxvsMTBD_baseline_sig_proteins_ECF$'Accession')) %>% 
  'colnames<-' ("Accession")

#Already made df for unique in mdx vs WT baseline not post-scruff above

#Could remove if want to keep DEP in mdxvsMTBD post-scruff
mdxvsWT_baseline_not_mdxvsMTBD_scruff_ECF_unique <- 
  data.frame('unique_mdxvsWT_baseline_not_mdxvsMTBD_scruff_ECF' = base::setdiff
             (mdxvsWT_baseline_sig_proteins_ECF$'Accession',
               mdxvsMTBD_scruff_sig_proteins_ECF$'Accession')) %>% 
  'colnames<-' ("Accession")

only_mdxvsWT_baseline_sig_ECF <- mdxvsWT_baseline_sig_proteins_ECF %>%
  filter(Accession %in% mdxvsWT_not_MTBD_ECF_baseline_unique$Accession &
           Accession %in% mdxvsWT_ECF_baseline_unique$Accession &
           Accession %in% mdxvsWT_baseline_not_mdxvsMTBD_scruff_ECF_unique$Accession)

full_only_mdxvsWT_baseline_sig_ECF <- ECF_proteomics_data %>%
  filter(Accession %in% mdxvsWT_not_MTBD_ECF_baseline_unique$Accession &
           Accession %in% mdxvsWT_ECF_baseline_unique$Accession &
           Accession %in% mdxvsWT_baseline_not_mdxvsMTBD_scruff_ECF_unique$Accession)

#write.xlsx(full_only_mdxvsWT_baseline_sig_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/ECF_mdxvsWT_baseline_unique_DEP_updated.xlsx",
#           sheetName="full dataset")

----

# Filtering datasets by sig in ECF mdx vs WT scruff and not other conditions ----
  
mdxvsWT_not_MTBD_ECF_scruff_unique <- 
  data.frame('unique_mdxvsWT_not_MTBD_scruff_ECF' = base::setdiff
             (mdxvsWT_scruff_sig_proteins_ECF$'Accession',
               mdxvsMTBD_scruff_sig_proteins_ECF$'Accession')) %>% 
  'colnames<-' ("Accession")

#Already made df for unique in mdx vs WT post-scruff not baseline above

#Could remove if want to keep DEP in mdxvsMTBD baseline
mdxvsWT_scruff_not_mdxvsMTBD_baseline_ECF_unique <- 
  data.frame('unique_mdxvsWT_scruff_not_mdxvsMTBD_baseline_ECF' = base::setdiff
             (mdxvsWT_scruff_sig_proteins_ECF$'Accession',
               mdxvsMTBD_baseline_sig_proteins_ECF$'Accession')) %>% 
  'colnames<-' ("Accession")

only_mdxvsWT_scruff_sig_ECF <- mdxvsWT_scruff_sig_proteins_ECF %>%
  filter(Accession %in% mdxvsWT_not_MTBD_ECF_scruff_unique$Accession &
           Accession %in% mdxvsWT_ECF_scruff_unique$Accession &
           Accession %in% mdxvsWT_scruff_not_mdxvsMTBD_baseline_ECF_unique$Accession)

full_only_mdxvsWT_scruff_sig_ECF <- ECF_proteomics_data %>%
  filter(Accession %in% mdxvsWT_not_MTBD_ECF_scruff_unique$Accession &
           Accession %in% mdxvsWT_ECF_scruff_unique$Accession &
           Accession %in% mdxvsWT_scruff_not_mdxvsMTBD_baseline_ECF_unique$Accession)

#write.xlsx(full_only_mdxvsWT_scruff_sig_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/ECF_mdxvsWT_scruff_unique_DEP.xlsx",
#           sheetName="full dataset")

----
# Filtering datasets by sig in ECF mdx vs MTBD baseline and not other conditions ----

mdxvsMTBD_not_WT_ECF_baseline_unique <- 
  data.frame('unique_mdxvsMTBD_not_WT_baseline_ECF' = base::setdiff
             (mdxvsMTBD_baseline_sig_proteins_ECF$'Accession',
               mdxvsWT_baseline_sig_proteins_ECF$'Accession')) %>% 
  'colnames<-' ("Accession")

#Could remove if want to keep all mdxvsMTBD proteins, regardless of scruff
mdxvsMTBD_baseline_not_scruff_ECF_unique <- 
  data.frame('unique_mdxvsMTBD_baseline_not_scruff_ECF' = base::setdiff
             (mdxvsMTBD_baseline_sig_proteins_ECF$'Accession',
               mdxvsMTBD_scruff_sig_proteins_ECF$'Accession')) %>%
  'colnames<-' ("Accession")

#Could remove if want to keep proteins also different in mdx vs WT post-scruff
mdxvsMTBD_baseline_not_mdxvsWTscruff_ECF_unique <- 
  data.frame('unique_mdxvsMTBD_baseline_not_mdxvsWTscruff_ECF' = base::setdiff
             (mdxvsMTBD_baseline_sig_proteins_ECF$'Accession',
               mdxvsWT_scruff_sig_proteins_ECF$'Accession')) %>%
  'colnames<-' ("Accession")

only_mdxvsMTBD_baseline_sig_ECF <- mdxvsMTBD_baseline_sig_proteins_ECF %>%
  filter(Accession %in% mdxvsMTBD_not_WT_ECF_baseline_unique$Accession &
           Accession %in% mdxvsMTBD_baseline_not_scruff_ECF_unique$Accession &
           Accession %in% mdxvsMTBD_baseline_not_mdxvsWTscruff_ECF_unique$Accession)

full_only_mdxvsMTBD_baseline_sig_ECF <- ECF_proteomics_data %>%
  filter(Accession %in% mdxvsMTBD_not_WT_ECF_baseline_unique$Accession &
           Accession %in% mdxvsMTBD_baseline_not_scruff_ECF_unique$Accession &
           Accession %in% mdxvsMTBD_baseline_not_mdxvsWTscruff_ECF_unique$Accession)

#write.xlsx(full_only_mdxvsMTBD_baseline_sig_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/ECF_mdxvsMTBD_baseline_unique_DEP.xlsx",
#           sheetName="full dataset")

----
  
# Filtering datasets by sig in ECF mdx vs MTBD scruff and not other conditions ----

mdxvsMTBD_not_WT_ECF_scruff_unique <- 
  data.frame('unique_mdxvsMTBD_not_WT_scruff_ECF' = base::setdiff
             (mdxvsMTBD_scruff_sig_proteins_ECF$'Accession',
               mdxvsWT_scruff_sig_proteins_ECF$'Accession')) %>% 
  'colnames<-' ("Accession")

#Could remove if want to keep all MTBD proteins, regardless of scruff
mdxvsMTBD_scruff_not_baseline_ECF_unique <- 
  data.frame('unique_mdxvsMTBD_scruff_not_baseline_ECF' = base::setdiff
             (mdxvsMTBD_scruff_sig_proteins_ECF$'Accession',
               mdxvsMTBD_baseline_sig_proteins_ECF$'Accession')) %>%
               'colnames<-' ("Accession")

#Could remove if want to keep DEP in mdx vs WT baseline
mdxvsMTBD_scruff_not_mdxvsWTbaseline_ECF_unique <- 
  data.frame('unique_mdxvsMTBD_scruff_not_mdxvsWTbaseline_ECF' = base::setdiff
             (mdxvsMTBD_scruff_sig_proteins_ECF$'Accession',
               mdxvsWT_baseline_sig_proteins_ECF$'Accession')) %>%
  'colnames<-' ("Accession")

only_mdxvsMTBD_scruff_sig_ECF <- mdxvsMTBD_scruff_sig_proteins_ECF %>%
  filter(Accession %in% mdxvsMTBD_not_WT_ECF_scruff_unique$Accession &
           Accession %in% mdxvsMTBD_scruff_not_baseline_ECF_unique$Accession &
           Accession %in% mdxvsMTBD_scruff_not_mdxvsWTbaseline_ECF_unique$Accession)

full_only_mdxvsMTBD_scruff_sig_ECF <- ECF_proteomics_data %>%
  filter(Accession %in% mdxvsMTBD_not_WT_ECF_scruff_unique$Accession &
           Accession %in% mdxvsMTBD_scruff_not_baseline_ECF_unique$Accession &
           Accession %in% mdxvsMTBD_scruff_not_mdxvsWTbaseline_ECF_unique$Accession)

#write.xlsx(full_only_mdxvsMTBD_scruff_sig_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/ECF_mdxvsMTBD_scruff_unique_DEP.xlsx",
#           sheetName="full dataset")

----
  
# Filtering datasets by sig in WM mdx vs WT baseline and not other conditions ----

mdxvsWT_not_MTBD_WM_baseline_unique <- 
  data.frame('unique_mdxvsWT_not_MTBD_baseline_WM' = base::setdiff
             (mdxvsWT_baseline_sig_proteins_WM$'Accession',
               mdxvsMTBD_baseline_sig_proteins_WM$'Accession')) %>% 
  'colnames<-' ("Accession")

#Already made df for unique in mdx vs WT baseline not post-scruff above

#Could remove if want to keep DEP in mdxvsMTBD post-scruff
mdxvsWT_baseline_not_mdxvsMTBD_scruff_WM_unique <- 
  data.frame('unique_mdxvsWT_baseline_not_mdxvsMTBD_scruff_WM' = base::setdiff
             (mdxvsWT_baseline_sig_proteins_WM$'Accession',
               mdxvsMTBD_scruff_sig_proteins_WM$'Accession')) %>% 
  'colnames<-' ("Accession")

only_mdxvsWT_baseline_sig_WM <- mdxvsWT_baseline_sig_proteins_WM %>%
  filter(Accession %in% mdxvsWT_not_MTBD_WM_baseline_unique$Accession &
           Accession %in% mdxvsWT_WM_baseline_unique$Accession &
           Accession %in% mdxvsWT_baseline_not_mdxvsMTBD_scruff_WM_unique$Accession)

full_only_mdxvsWT_baseline_sig_WM <- WM_proteomics_data %>%
  filter(Accession %in% mdxvsWT_not_MTBD_WM_baseline_unique$Accession &
           Accession %in% mdxvsWT_WM_baseline_unique$Accession &
           Accession %in% mdxvsWT_baseline_not_mdxvsMTBD_scruff_WM_unique$Accession)

#write.xlsx(full_only_mdxvsWT_baseline_sig_WM, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/WM_mdxvsWT_baseline_unique_DEP_updated.xlsx",
#           sheetName="full dataset")

----
  
# Filtering datasets by sig in WM mdx vs WT scruff and not other conditions ----

mdxvsWT_not_MTBD_WM_scruff_unique <- 
  data.frame('unique_mdxvsWT_not_MTBD_scruff_WM' = base::setdiff
             (mdxvsWT_scruff_sig_proteins_WM$'Accession',
               mdxvsMTBD_scruff_sig_proteins_WM$'Accession')) %>% 
  'colnames<-' ("Accession")

#Already made df for unique in mdx vs WT post-scruff not baseline above

#Could remove if want to keep DEP in mdxvsMTBD baseline
mdxvsWT_scruff_not_mdxvsMTBD_baseline_WM_unique <- 
  data.frame('unique_mdxvsWT_scruff_not_mdxvsMTBD_baseline_WM' = base::setdiff
             (mdxvsWT_scruff_sig_proteins_WM$'Accession',
               mdxvsMTBD_baseline_sig_proteins_WM$'Accession')) %>% 
  'colnames<-' ("Accession")

only_mdxvsWT_scruff_sig_WM <- mdxvsWT_scruff_sig_proteins_WM %>%
  filter(Accession %in% mdxvsWT_not_MTBD_WM_scruff_unique$Accession &
           Accession %in% mdxvsWT_WM_scruff_unique$Accession &
           Accession %in% mdxvsWT_scruff_not_mdxvsMTBD_baseline_WM_unique$Accession)

full_only_mdxvsWT_scruff_sig_WM <- WM_proteomics_data %>%
  filter(Accession %in% mdxvsWT_not_MTBD_WM_scruff_unique$Accession &
           Accession %in% mdxvsWT_WM_scruff_unique$Accession &
           Accession %in% mdxvsWT_scruff_not_mdxvsMTBD_baseline_WM_unique$Accession)

#write.xlsx(full_only_mdxvsWT_scruff_sig_WM, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/WM_mdxvsWT_scruff_unique_DEP.xlsx",
#           sheetName="full dataset")

----

# Filtering datasets by sig in WM mdx vs MTBD baseline and not other conditions ----

mdxvsMTBD_not_WT_WM_baseline_unique <- 
  data.frame('unique_mdxvsMTBD_not_WT_baseline_WM' = base::setdiff
             (mdxvsMTBD_baseline_sig_proteins_WM$'Accession',
               mdxvsWT_baseline_sig_proteins_WM$'Accession')) %>% 
  'colnames<-' ("Accession")

#Could remove if want to keep all mdxvsMTBD proteins, regardless of scruff
mdxvsMTBD_baseline_not_scruff_WM_unique <- 
  data.frame('unique_mdxvsMTBD_baseline_not_scruff_WM' = base::setdiff
             (mdxvsMTBD_baseline_sig_proteins_WM$'Accession',
               mdxvsMTBD_scruff_sig_proteins_WM$'Accession')) %>%
  'colnames<-' ("Accession")

#Could remove if want to keep proteins also different in mdx vs WT post-scruff
mdxvsMTBD_baseline_not_mdxvsWTscruff_WM_unique <- 
  data.frame('unique_mdxvsMTBD_baseline_not_mdxvsWTscruff_WM' = base::setdiff
             (mdxvsMTBD_baseline_sig_proteins_WM$'Accession',
               mdxvsWT_scruff_sig_proteins_WM$'Accession')) %>%
  'colnames<-' ("Accession")

only_mdxvsMTBD_baseline_sig_WM <- mdxvsMTBD_baseline_sig_proteins_WM %>%
  filter(Accession %in% mdxvsMTBD_not_WT_WM_baseline_unique$Accession &
           Accession %in% mdxvsMTBD_baseline_not_scruff_WM_unique$Accession &
           Accession %in% mdxvsMTBD_baseline_not_mdxvsWTscruff_WM_unique$Accession)

full_only_mdxvsMTBD_baseline_sig_WM <- WM_proteomics_data %>%
  filter(Accession %in% mdxvsMTBD_not_WT_WM_baseline_unique$Accession &
           Accession %in% mdxvsMTBD_baseline_not_scruff_WM_unique$Accession &
           Accession %in% mdxvsMTBD_baseline_not_mdxvsWTscruff_WM_unique$Accession)

#write.xlsx(full_only_mdxvsMTBD_baseline_sig_WM, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/WM_mdxvsMTBD_baseline_unique_DEP.xlsx",
#           sheetName="full dataset")

----

# Filtering datasets by sig in WM mdx vs MTBD scruff and not other conditions ----

mdxvsMTBD_not_WT_WM_scruff_unique <- 
  data.frame('unique_mdxvsMTBD_not_WT_scruff_WM' = base::setdiff
             (mdxvsMTBD_scruff_sig_proteins_WM$'Accession',
               mdxvsWT_scruff_sig_proteins_WM$'Accession')) %>% 
  'colnames<-' ("Accession")

#Could remove if want to keep all MTBD proteins, regardless of scruff
mdxvsMTBD_scruff_not_baseline_WM_unique <- 
  data.frame('unique_mdxvsMTBD_scruff_not_baseline_WM' = base::setdiff
             (mdxvsMTBD_scruff_sig_proteins_WM$'Accession',
               mdxvsMTBD_baseline_sig_proteins_WM$'Accession')) %>%
  'colnames<-' ("Accession")

#Could remove if want to keep DEP in mdx vs WT baseline
mdxvsMTBD_scruff_not_mdxvsWTbaseline_WM_unique <- 
  data.frame('unique_mdxvsMTBD_scruff_not_mdxvsWTbaseline_WM' = base::setdiff
             (mdxvsMTBD_scruff_sig_proteins_WM$'Accession',
               mdxvsWT_baseline_sig_proteins_WM$'Accession')) %>%
  'colnames<-' ("Accession")

only_mdxvsMTBD_scruff_sig_WM <- mdxvsMTBD_scruff_sig_proteins_WM %>%
  filter(Accession %in% mdxvsMTBD_not_WT_WM_scruff_unique$Accession &
           Accession %in% mdxvsMTBD_scruff_not_baseline_WM_unique$Accession &
           Accession %in% mdxvsMTBD_scruff_not_mdxvsWTbaseline_WM_unique$Accession)

full_only_mdxvsMTBD_scruff_sig_WM <- WM_proteomics_data %>%
  filter(Accession %in% mdxvsMTBD_not_WT_WM_scruff_unique$Accession &
           Accession %in% mdxvsMTBD_scruff_not_baseline_WM_unique$Accession &
           Accession %in% mdxvsMTBD_scruff_not_mdxvsWTbaseline_WM_unique$Accession)

#write.xlsx(full_only_mdxvsMTBD_scruff_sig_WM, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/WM_mdxvsMTBD_scruff_unique_DEP.xlsx",
#           sheetName="full dataset")

----
  
# Filtering datasets by overlapping sig in ECF mdx vs WT and mdx vs MTBD ----

# Baseline (mutually DEP in mdx vs WT and vs MTBD)

baseline_overlap_ECF <- 
  data.frame('baseline_overlap_ECF' = intersect
             (mdxvsWT_baseline_sig_proteins_ECF$Accession,
               mdxvsMTBD_baseline_sig_proteins_ECF$Accession)) %>% 
  'colnames<-' ("Accession")

uniquely_overlapping_baseline_ECF <-
  baseline_overlap_ECF %>% 
  base::subset(!(.$Accession %in% mdxvsWT_scruff_sig_proteins_ECF$Accession)) %>%
  as.data.frame(.) %>%
  base::subset(!(.$Accession %in% mdxvsMTBD_scruff_sig_proteins_ECF$Accession)) %>%
  as.data.frame(.)

full_uniquely_overlapping_baseline_ECF <- ECF_proteomics_data %>%
  filter(Accession %in% uniquely_overlapping_baseline_ECF$Accession)

#write.csv(full_uniquely_overlapping_baseline_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/uniquely_overlapping_baseline_ECF.csv")

# Post-scruff (mutually DEP in mdx vs WT and vs MTBD)

scruff_overlap_ECF <- 
  data.frame('scruff_overlap_ECF' = intersect
             (mdxvsWT_scruff_sig_proteins_ECF$Accession,
               mdxvsMTBD_scruff_sig_proteins_ECF$Accession)) %>% 
  'colnames<-' ("Accession")

uniquely_overlapping_scruff_ECF <-
  scruff_overlap_ECF %>% 
  base::subset(!(.$Accession %in% mdxvsWT_baseline_sig_proteins_ECF$Accession)) %>%
  as.data.frame(.) %>%
  base::subset(!(.$Accession %in% mdxvsMTBD_baseline_sig_proteins_ECF$Accession)) %>%
  as.data.frame(.)

full_uniquely_overlapping_scruff_ECF <- ECF_proteomics_data %>%
  filter(Accession %in% uniquely_overlapping_scruff_ECF$Accession)

#write.csv(full_uniquely_overlapping_scruff_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/uniquely_overlapping_scruff_ECF.csv")

#Scruff ECF secreted proteins, both overlapping and unique in mdx vs WT and mdx vs MTBD

accession_only_mdxvsWT_scruff_ECF <-
  data.frame(only_mdxvsWT_scruff_sig_ECF$Accession) %>% 
  'colnames<-' ("Accession")

accession_only_mdxvsMTBD_scruff_ECF <-
  data.frame(only_mdxvsMTBD_scruff_sig_ECF$Accession) %>% 
  'colnames<-' ("Accession")

overlap_and_unique_scruff_ECF <-
  rbind(uniquely_overlapping_scruff_ECF,
        accession_only_mdxvsWT_scruff_ECF,
        accession_only_mdxvsMTBD_scruff_ECF)

#write.csv(overlap_and_unique_scruff_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/overlapping_and_unique_scruff_ECF.csv")

-----
  

# Filtering datasets by overlapping sig in WM mdx vs WT and mdx vs MTBD ----

# Baseline (mutually DEP in mdx vs WT and vs MTBD)

baseline_overlap_WM <- 
  data.frame('baseline_overlap_WM' = intersect
             (mdxvsWT_baseline_sig_proteins_WM$Accession,
               mdxvsMTBD_baseline_sig_proteins_WM$Accession)) %>% 
  'colnames<-' ("Accession")

uniquely_overlapping_baseline_WM <-
  baseline_overlap_WM %>% 
  base::subset(!(.$Accession %in% mdxvsWT_scruff_sig_proteins_WM$Accession)) %>%
  as.data.frame(.) %>%
  base::subset(!(.$Accession %in% mdxvsMTBD_scruff_sig_proteins_WM$Accession)) %>%
  as.data.frame(.)

full_uniquely_overlapping_baseline_WM <- WM_proteomics_data %>%
  filter(Accession %in% uniquely_overlapping_baseline_WM$Accession)

#write.csv(full_uniquely_overlapping_baseline_WM, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/uniquely_overlapping_baseline_WM.csv")

# Post-scruff (mutually DEP in mdx vs WT and vs MTBD)

scruff_overlap_WM <- 
  data.frame('scruff_overlap_WM' = intersect
             (mdxvsWT_scruff_sig_proteins_WM$Accession,
               mdxvsMTBD_scruff_sig_proteins_WM$Accession)) %>% 
  'colnames<-' ("Accession")

uniquely_overlapping_scruff_WM <-
  scruff_overlap_WM %>% 
  base::subset(!(.$Accession %in% mdxvsWT_baseline_sig_proteins_WM$Accession)) %>%
  as.data.frame(.) %>%
  base::subset(!(.$Accession %in% mdxvsMTBD_baseline_sig_proteins_WM$Accession)) %>%
  as.data.frame(.)

full_uniquely_overlapping_scruff_WM <- WM_proteomics_data %>%
  filter(Accession %in% uniquely_overlapping_scruff_WM$Accession)

#write.csv(full_uniquely_overlapping_scruff_WM, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/uniquely_overlapping_scruff_WM.csv")

-----
##################################################################################
# SubcellulaRVis -----

#devtools::install_github("jowatson2011/subcellularvis")

library(subcellularvis)

#subcellularapp()

# List of Uniprot IDs for all ECF quantified proteins
ECF_all_proteins_UniprotID <- data.frame(ECF_proteomics_data$Accession) %>%
  'colnames<-' ("Uniprot ID")

#write.csv(ECF_all_proteins_UniprotID, 
#          file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/ECF_all_proteins_UniprotID.csv")

# List of Uniprot IDs for only ECF (not WM) quantified proteins

ECF_proteome_only_UniprotID <- 
  data.frame('ECF_proteome_only' = base::setdiff
             (ECF_proteomics_data$'Accession',
               WM_proteomics_data$'Accession')) %>% 
  'colnames<-' ("Uniprot ID")

#write.csv(ECF_proteome_only_UniprotID, 
#          file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/ECF_proteome_only_UniprotID.csv")

# List of Uniprot IDs for all WM quantified proteins

WM_all_proteins_UniprotID <- data.frame(WM_proteomics_data$Accession) %>%
  'colnames<-' ("Uniprot ID")  

#write.csv(WM_all_proteins_UniprotID, 
#          file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/WM_all_proteins_UniprotID.csv")

# List of Uniprot IDs for only WM (not ECF) quantified proteins

WM_proteome_only_UniprotID <- 
  data.frame('WM_proteome_only' = base::setdiff
             (WM_proteomics_data$'Accession',
               ECF_proteomics_data$'Accession')) %>% 
  'colnames<-' ("Uniprot ID")

#write.csv(WM_proteome_only_UniprotID, 
#          file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/WM_proteome_only_UniprotID.csv")

# List of Uniprot IDs for overlapping ECF and WM quantified proteins

ECF_WM_proteome_overlap_UniprotID <- 
  data.frame('ECF_WM_proteome_overlap' = intersect
             (ECF_proteomics_data$'Accession',
               WM_proteomics_data$'Accession')) %>% 
  'colnames<-' ("Uniprot ID")

#write.csv(ECF_WM_proteome_overlap_UniprotID, 
#          file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF_WM_proteome_overlap_UniprotID.csv")

# Venn Diagram ECF vs WM proteomes
ECFvsWM_whole_proteomes <- display_venn(ECFvsWM_whole_proteomes_venn_list <- list(
  'ECF proteome' = ECF_all_proteins_UniprotID$`Uniprot ID`,
  'WM proteome' = WM_all_proteins_UniprotID$`Uniprot ID`), 
  # Circles
  lwd = 2,
  lty = 'blank',
  category.names = c("ECF proteome", 
                     "QM proteome"),
  fill = c("#64839b","#9b7c64"),
  # Numbers
  cex = 2,
  fontface = "italic",
  # Set names
  cat.cex = 1.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-115, 115),
  cat.dist = c(-.045, -.05),
  main = "Quadriceps ECF vs Quad Muscle Proteome Overlap",
  main.fontface = "bold"
  #rotation.degree = 180
)

-----

########################################################################
# Subsetting 2-group comparison ECF datasets by ECF-only proteome filter
# to identify potential secreted protein candidates
########################################################################

# Secreted protein candidate datasets

#mdx vs WT baseline
secreted_cand_mdxvsWT_baseline_ECF <- 
  data.frame('secreted_cand_mdxvsWT_baseline_ECF' = intersect
             (mdxvsWT_baseline_sig_proteins_ECF$'Accession', 
               ECF_proteome_only_UniprotID$`Uniprot ID`)) %>% 
  'colnames<-' ("Accession")

#mdx vs WT post-scruff
secreted_cand_mdxvsWT_scruff_ECF <- 
  data.frame('secreted_cand_mdxvsWT_scruff_ECF' = intersect
             (mdxvsWT_scruff_sig_proteins_ECF$'Accession', 
               ECF_proteome_only_UniprotID$`Uniprot ID`)) %>% 
  'colnames<-' ("Accession")

#mdx vs MTBD baseline
secreted_cand_mdxvsMTBD_baseline_ECF <- 
  data.frame('secreted_cand_mdxvsMTBD_baseline_ECF' = intersect
             (mdxvsMTBD_baseline_sig_proteins_ECF$'Accession', 
               ECF_proteome_only_UniprotID$`Uniprot ID`)) %>% 
  'colnames<-' ("Accession")

#mdx vs MTBD post-scruff
secreted_cand_mdxvsMTBD_scruff_ECF <- 
  data.frame('secreted_cand_mdxvsMTBD_scruff_ECF' = intersect
             (mdxvsMTBD_scruff_sig_proteins_ECF$'Accession', 
               ECF_proteome_only_UniprotID$`Uniprot ID`)) %>% 
  'colnames<-' ("Accession")

#MTBD vs WT baseline
secreted_cand_MTBDvsWT_baseline_ECF <- 
  data.frame('secreted_cand_MTBDvsWT_baseline_ECF' = intersect
             (MTBDvsWT_baseline_sig_proteins_ECF$'Accession', 
               ECF_proteome_only_UniprotID$`Uniprot ID`)) %>% 
  'colnames<-' ("Accession")

#MTBD vs WT post-scruff
secreted_cand_MTBDvsWT_scruff_ECF <- 
  data.frame('secreted_cand_MTBDvsWT_scruff_ECF' = intersect
             (MTBDvsWT_scruff_sig_proteins_ECF$'Accession', 
               ECF_proteome_only_UniprotID$`Uniprot ID`)) %>% 
  'colnames<-' ("Accession")

# Secreted protein candidate non-overlap exported datasets

# mdx vs WT baseline
secreted_only_mdxvsWT_baseline_ECF <- only_mdxvsWT_baseline_sig_ECF %>%
  filter(Accession %in% ECF_proteome_only_UniprotID$`Uniprot ID`)

full_secreted_only_mdxvsWT_baseline_ECF <- ECF_proteomics_data %>%
  filter(Accession %in% secreted_only_mdxvsWT_baseline_ECF$Accession)

#write.xlsx(full_secreted_only_mdxvsWT_baseline_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/secreted_ECF_mdxvsWT_baseline.xlsx",
#           sheetName="full dataset")

# mdx vs WT post-scruff
secreted_only_mdxvsWT_scruff_ECF <- only_mdxvsWT_scruff_sig_ECF %>%
  filter(Accession %in% ECF_proteome_only_UniprotID$`Uniprot ID`)

full_secreted_only_mdxvsWT_scruff_ECF <- ECF_proteomics_data %>%
  filter(Accession %in% secreted_only_mdxvsWT_scruff_ECF$Accession)

#write.xlsx(full_secreted_only_mdxvsWT_scruff_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/secreted_ECF_mdxvsWT_scruff.xlsx",
#           sheetName="full dataset")

# mdx vs MTBD baseline
secreted_only_mdxvsMTBD_baseline_ECF <- only_mdxvsMTBD_baseline_sig_ECF %>%
  filter(Accession %in% ECF_proteome_only_UniprotID$`Uniprot ID`)

full_secreted_only_mdxvsMTBD_baseline_ECF <- ECF_proteomics_data %>%
  filter(Accession %in% secreted_only_mdxvsMTBD_baseline_ECF$Accession)

#write.xlsx(full_secreted_only_mdxvsMTBD_baseline_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/secreted_ECF_mdxvsMTBD_baseline.xlsx",
#           sheetName="full dataset")

# mdx vs MTBD post-scruff
secreted_only_mdxvsMTBD_scruff_ECF <- only_mdxvsMTBD_scruff_sig_ECF %>%
  filter(Accession %in% ECF_proteome_only_UniprotID$`Uniprot ID`)

full_secreted_only_mdxvsMTBD_scruff_ECF <- ECF_proteomics_data %>%
  filter(Accession %in% secreted_only_mdxvsMTBD_scruff_ECF$Accession)

#write.xlsx(full_secreted_only_mdxvsMTBD_scruff_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/secreted_ECF_mdxvsMTBD_scruff.xlsx",
#           sheetName="full dataset")

# Baseline (mutually DEP in mdx vs WT and vs MTBD)

secreted_baseline_overlap_ECF <- 
  data.frame('secreted_baseline_overlap_ECF' = intersect
             (secreted_cand_mdxvsWT_baseline_ECF$Accession,
               secreted_cand_mdxvsMTBD_baseline_ECF$Accession)) %>% 
  'colnames<-' ("Accession")

full_secreted_baseline_overlap_ECF <- ECF_proteomics_data %>%
  filter(Accession %in% secreted_baseline_overlap_ECF$Accession)

#write.xlsx(full_secreted_baseline_overlap_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/secreted_baseline_with_scruff_pool_ECF.xlsx",
#           sheetName="secreted_not_picky_baseline")

secreted_uniquely_overlapping_baseline_ECF <-
  secreted_baseline_overlap_ECF %>% 
  base::subset(!(.$Accession %in% secreted_cand_mdxvsWT_scruff_ECF$Accession)) %>%
  as.data.frame(.) %>%
  base::subset(!(.$Accession %in% secreted_cand_mdxvsMTBD_scruff_ECF$Accession)) %>%
  as.data.frame(.)
  
full_secreted_uniquely_overlapping_baseline_ECF <- ECF_proteomics_data %>%
  filter(Accession %in% secreted_uniquely_overlapping_baseline_ECF$Accession)

#write.xlsx(full_secreted_uniquely_overlapping_baseline_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/secreted_uniquely_overlapping_baseline_ECF.xlsx",
#           sheetName="full dataset")

# Post-scruff (mutually DEP in mdx vs WT and vs MTBD)

secreted_scruff_overlap_ECF <- 
  data.frame('secreted_scruff_overlap_ECF' = intersect
             (secreted_cand_mdxvsWT_scruff_ECF$Accession,
               secreted_cand_mdxvsMTBD_scruff_ECF$Accession)) %>% 
  'colnames<-' ("Accession")

full_secreted_scruff_overlap_ECF <- ECF_proteomics_data %>%
  filter(Accession %in% secreted_scruff_overlap_ECF$Accession)

#write.xlsx(full_secreted_scruff_overlap_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/secreted_scruff_with_baseline_pool_ECF.xlsx",
#           sheetName="secreted_not_picky_scruff")

secreted_uniquely_overlapping_scruff_ECF <-
  secreted_scruff_overlap_ECF %>% 
  base::subset(!(.$Accession %in% secreted_cand_mdxvsWT_baseline_ECF$Accession)) %>%
  as.data.frame(.) %>%
  base::subset(!(.$Accession %in% secreted_cand_mdxvsMTBD_baseline_ECF$Accession)) %>%
  as.data.frame(.)

full_secreted_uniquely_overlapping_scruff_ECF <- ECF_proteomics_data %>%
  filter(Accession %in% secreted_uniquely_overlapping_scruff_ECF$Accession)

#write.xlsx(full_secreted_uniquely_overlapping_scruff_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/secreted_uniquely_overlapping_scruff_ECF.xlsx",
#           sheetName="full dataset")

#Scruff ECF secreted proteins, both overlapping and unique in mdx vs WT and mdx vs MTBD

accession_secreted_only_mdxvsWT_scruff_ECF <-
  data.frame(secreted_only_mdxvsWT_scruff_ECF$Accession) %>% 
  'colnames<-' ("Accession")

accession_secreted_only_mdxvsMTBD_scruff_ECF <-
  data.frame(secreted_only_mdxvsMTBD_scruff_ECF$Accession) %>% 
  'colnames<-' ("Accession")

overlap_and_unique_secreted_scruff_ECF <-
  rbind(secreted_uniquely_overlapping_scruff_ECF,
        accession_secreted_only_mdxvsWT_scruff_ECF,
        accession_secreted_only_mdxvsMTBD_scruff_ECF)

#write.csv(overlap_and_unique_secreted_scruff_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/secreted_overlapping_and_unique_scruff_ECF.csv")

# Common overlap between all mdx DEPs found in ECF proteome only (12 proteins)

secreted_all_group_overlap_ECF <- 
  data.frame('secreted_all_group_overlap_ECF' = intersect
             (secreted_baseline_overlap_ECF$Accession,
               secreted_scruff_overlap_ECF$Accession)) %>% 
  'colnames<-' ("Accession")

full_secreted_all_group_overlap_ECF <- ECF_proteomics_data %>%
  filter(Accession %in% secreted_all_group_overlap_ECF$Accession)

#write.xlsx(full_secreted_all_group_overlap_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/secreted_all_group_overlap_ECF.xlsx",
#           sheetName="full dataset")

# Secreted candidate ECF protein 4-D Venn comparisons WT, mdx, MTBD ----

ECF_4D_venn <- display_venn(ECF_secreted_cand_4D_venn_list <- list(
  'mdx vs WT baseline' = secreted_cand_mdxvsWT_baseline_ECF$'Accession',
  'mdx vs MTBD baseline' = secreted_cand_mdxvsMTBD_baseline_ECF$'Accession',
  'mdx vs WT post-scruff' = secreted_cand_mdxvsWT_scruff_ECF$'Accession',
  'mdx vs MTBD post-scruff' = secreted_cand_mdxvsMTBD_scruff_ECF$'Accession'), 
  # Circles
  lwd = 2,
  lty = 'blank',
  category.names = c(
    "mdx vs WT baseline DEP",
    "mdx vs MTBD baseline DEP",
    "mdx vs WT post-scruff DEP", 
    "mdx vs MTBD post-scruff DEP"),
  fill = c("#A8D0DB","#E49273","#7180AC", "#A37A74"),
  # Numbers
  cex = 1.2,
  fontface = "italic",
  # Set names
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-6, 1, -20, 5),
  cat.dist = c(0.14, .14, 0.07, 0.06),
  main = "Quadriceps Secreted Candidate ECF Protein Comparisons WT, mdx, and MTBD",
  main.fontface = "bold"
) 


#####################################
# Volcano Plots
#####################################
# WM and ECF protein DE datasets for volcano plots ----

  #WT vs mdx WM

mdxvsWT_baseline_WM <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/Data files for R/mdxB_WTB_WM_volcano.xlsx"))

mdxvsWT_scruff_WM <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/Data files for R/mdxS_WTS_WM_volcano.xlsx"))

  #WT vs MTBD WM

MTBDvsWT_baseline_WM <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/Data files for R/MTBD_B_WTB_WM_volcano.xlsx"))

MTBDvsWT_scruff_WM <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/Data files for R/MTBD_S_WTS_WM_volcano.xlsx"))

  #WT vs mdx ECF

mdxvsWT_baseline_ECF <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/Data files for R/mdxB_WTB_ECF_volcano.xlsx"))

mdxvsWT_scruff_ECF <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/Data files for R/mdxS_WTS_ECF_volcano.xlsx"))

  #MTBD vs WT ECF

MTBDvsWT_baseline_ECF <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/Data files for R/MTBD_B_WTB_ECF_volcano.xlsx"))

MTBDvsWT_scruff_ECF <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/Data files for R/MTBD_S_WTS_ECF_volcano.xlsx"))

  #MTBD vs mdx WM

mdxvsMTBD_baseline_WM <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/Data files for R/mdxB_MTBD_B_WM_volcano.xlsx"))

mdxvsMTBD_scruff_WM <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/Data files for R/mdxS_MTBD_S_WM_volcano.xlsx"))

  #MTBD vs mdx ECF

mdxvsMTBD_baseline_ECF <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/Data files for R/mdxB_MTBD_B_ECF_volcano.xlsx"))

mdxvsMTBD_scruff_ECF <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/Data files for R/mdxS_MTBD_S_ECF_volcano.xlsx"))

  # Scruff vs baseline same genotype comparisons WM

WT_scruff_vs_baseline_WM <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/Data files for R/WTS_WTB_WM_volcano.xlsx"))

mdx_scruff_vs_baseline_WM <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/Data files for R/mdxS_mdxB_WM_volcano.xlsx"))

MTBD_scruff_vs_baseline_WM <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/Data files for R/MTBD_S_MTBD_B_WM_volcano.xlsx"))

  # Scruff vs baseline same genotype comparisons ECF

WT_scruff_vs_baseline_ECF <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/Data files for R/WTS_WTB_ECF_volcano.xlsx"))

mdx_scruff_vs_baseline_ECF <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/Data files for R/mdxS_mdxB_ECF_volcano.xlsx"))

MTBD_scruff_vs_baseline_ECF <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/Data files for R/MTBD_S_MTBD_B_ECF_volcano.xlsx"))

----

#################################
# Defining Volcano plot function 
volcano_fun <- function(volcano_df){
    volcano_df %>%
      # Add a threshold for significant observations
      mutate(threshold =
               if_else(`log_fc` >= 1 & `log_pval` >= 1.3, "A",
                       if_else(`log_fc` <= -1 & `log_pval` >= 1.3,"B","C"))) %>%
      # Plot with points coloured according to the threshold
      ggplot(aes(x=`log_fc`, y=`log_pval`, colour = threshold)) +
      geom_point(alpha = 0.5) + # Alpha sets the transparency of the points
      # Add dotted lines to indicate the threshold, semi-transparent
      geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) + 
      geom_vline(xintercept = c(-1,1), linetype = 2, alpha = 0.5) +
      coord_cartesian(ylim = c(0, 5), xlim = c(-4, 5)) +
      labs(color = 'Legend') +
      # Set the colour of the points
      scale_colour_manual(values = c("A"= "#E2635C", "B"= "#007EA7", 
                          "C" = "#8F9491"),
                          labels = c("Upregulated", "Downregulated", "Not Significant")) +
      xlab("log2(fold change)") + ylab("-log10(adjusted P-value)") + # Relabel the axes
      #ggtitle("ECF mdx vs MTBD baseline DEPs") +
      theme_minimal() + #Set the theme
      theme(legend.position = c(.98, .98),
            legend.justification = c("right", "top"),
            legend.box.just = "right",
            legend.margin = margin(6, 6, 6, 6,),
            legend.box.background = element_rect(color="#427681",size=1)) +
    theme(text = element_text(size=14))
}
volcano_fun_2 <- function(volcano_df){
  volcano_df %>%
    # Add a threshold for significant observations
    mutate(threshold =
             if_else(`log_fc` >= 1 & `log_pval` >= 1.3, "A",
                     if_else(`log_fc` <= -1 & `log_pval` >= 1.3,"B","C"))) %>%
    # Plot with points coloured according to the threshold
    ggplot(aes(x=`log_fc`, y=`log_pval`, colour = threshold)) +
    geom_point(alpha = 0.5) + # Alpha sets the transparency of the points
    # Add dotted lines to indicate the threshold, semi-transparent
    geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) + 
    geom_vline(xintercept = c(-1,1), linetype = 2, alpha = 0.5) +
    coord_cartesian(ylim = c(0, 5), xlim = c(-4, 2.5)) +
    labs(color = 'Legend') +
    # Set the colour of the points
    scale_colour_manual(values = c("A"= "#E2635C", "B"= "#007EA7", 
                                   "C" = "#8F9491"),
                        labels = c("Downregulated", "Not Significant")) +
    xlab("log2(fold change)") + ylab("-log10(adjusted P-value)") + # Relabel the axes
    #ggtitle("ECF mdx vs MTBD baseline DEPs") +
    theme_minimal() + #Set the theme
    theme(legend.position = c(.97, .97),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(6, 6, 6, 6,),
          legend.box.background = element_rect(color="#427681",size=1)) +
    theme(text = element_text(size=14))
}
volcano_fun_3 <- function(volcano_df){
  volcano_df %>%
    # Add a threshold for significant observations
    mutate(threshold =
             if_else(`log_fc` >= 1 & `log_pval` >= 1.3, "A",
                     if_else(`log_fc` <= -1 & `log_pval` >= 1.3,"B","C"))) %>%
    # Plot with points coloured according to the threshold
    ggplot(aes(x=`log_fc`, y=`log_pval`, colour = threshold)) +
    geom_point(alpha = 0.5) + # Alpha sets the transparency of the points
    # Add dotted lines to indicate the threshold, semi-transparent
    geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) + 
    geom_vline(xintercept = c(-1,1), linetype = 2, alpha = 0.5) +
    coord_cartesian(ylim = c(0, 3), xlim = c(-4, 4)) +
    labs(color = 'Legend') +
    # Set the colour of the points
    scale_colour_manual(values = c("A"= "#E2635C", "B"= "#007EA7", 
                                   "C" = "#8F9491"),
                        labels = c("Downregulated", "Not Significant")) +
    xlab("log2(fold change)") + ylab("-log10(adjusted P-value)") + # Relabel the axes
    #ggtitle("ECF mdx vs MTBD baseline DEPs") +
    theme_minimal() + #Set the theme
    theme(legend.position = c(.97, .97),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(6, 6, 6, 6,),
          legend.box.background = element_rect(color="#427681",size=1)) +
    theme(text = element_text(size=14))
}
volcano_fun_4 <- function(volcano_df){
  volcano_df %>%
    # Add a threshold for significant observations
    mutate(threshold =
             if_else(`log_fc` >= 1 & `log_pval` >= 1.3, "A",
                     if_else(`log_fc` <= -1 & `log_pval` >= 1.3,"B","C"))) %>%
    # Plot with points coloured according to the threshold
    ggplot(aes(x=`log_fc`, y=`log_pval`, colour = threshold)) +
    geom_point(alpha = 0.5) + # Alpha sets the transparency of the points
    # Add dotted lines to indicate the threshold, semi-transparent
    geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) + 
    geom_vline(xintercept = c(-1,1), linetype = 2, alpha = 0.5) +
    coord_cartesian(ylim = c(0, 2), xlim = c(-4, 5)) +
    labs(color = 'Legend') +
    # Set the colour of the points
    scale_colour_manual(values = c("A"= "#E2635C", "B"= "#007EA7", 
                                   "C" = "#8F9491"),
                        labels = "Not Significant") +
    xlab("log2(fold change)") + ylab("-log10(adjusted P-value)") + # Relabel the axes
    #ggtitle("ECF mdx vs MTBD baseline DEPs") +
    theme_minimal() + #Set the theme
    theme(legend.position = c(.97, .97),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(6, 6, 6, 6,),
          legend.box.background = element_rect(color="#427681",size=1)) +
    theme(text = element_text(size=14))
}
volcano_fun_5 <- function(volcano_df){
  volcano_df %>%
    # Add a threshold for significant observations
    mutate(threshold =
             if_else(`log_fc` >= 1 & `log_pval` >= 1.3, "A",
                     if_else(`log_fc` <= -1 & `log_pval` >= 1.3,"B","C"))) %>%
    # Plot with points coloured according to the threshold
    ggplot(aes(x=`log_fc`, y=`log_pval`, colour = threshold)) +
    geom_point(alpha = 0.5) + # Alpha sets the transparency of the points
    # Add dotted lines to indicate the threshold, semi-transparent
    geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) + 
    geom_vline(xintercept = c(-1,1), linetype = 2, alpha = 0.5) +
    coord_cartesian(ylim = c(0, 6.5), xlim = c(-3, 5)) +
    labs(color = 'Legend') +
    # Set the colour of the points
    scale_colour_manual(values = c("A"= "#E2635C", "B"= "#007EA7", 
                                   "C" = "#8F9491"),
                        labels = c("Upregulated", "Downregulated", "Not Significant")) +
    xlab("log2(fold change)") + ylab("-log10(adjusted P-value)") + # Relabel the axes
    #ggtitle("ECF mdx vs MTBD baseline DEPs") +
    theme_minimal() + #Set the theme
          theme(legend.position = c(.99, .99),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(6, 6, 6, 6,),
          legend.box.background = element_rect(color="#427681",size=1)) +
    theme(text = element_text(size=14))
}
volcano_fun_6 <- function(volcano_df){
  volcano_df %>%
    # Add a threshold for significant observations
    mutate(threshold =
             if_else(`log_fc` >= 1 & `log_pval` >= 1.3, "A",
                     if_else(`log_fc` <= -1 & `log_pval` >= 1.3,"B","C"))) %>%
    # Plot with points coloured according to the threshold
    ggplot(aes(x=`log_fc`, y=`log_pval`, colour = threshold)) +
    geom_point(alpha = 0.5) + # Alpha sets the transparency of the points
    # Add dotted lines to indicate the threshold, semi-transparent
    geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) + 
    geom_vline(xintercept = c(-1,1), linetype = 2, alpha = 0.5) +
    coord_cartesian(ylim = c(0, 6), xlim = c(-5, 5)) +
    labs(color = 'Legend') +
    # Set the colour of the points
    scale_colour_manual(values = c("A"= "#E2635C", "B"= "#007EA7", 
                                   "C" = "#8F9491"),
                        labels = c("Upregulated", "Downregulated", "Not Significant")) +
    xlab("log2(fold change)") + ylab("-log10(adjusted P-value)") + # Relabel the axes
    #ggtitle("ECF mdx vs MTBD baseline DEPs") +
    theme_minimal() + #Set the theme
    theme(legend.position = c(.97, .97),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(6, 6, 6, 6,),
          legend.box.background = element_rect(color="#427681",size=1)) +
    theme(text = element_text(size=14))
}
volcano_fun_7 <- function(volcano_df){
  volcano_df %>%
    # Add a threshold for significant observations
    mutate(threshold =
             if_else(`log_fc` >= 1 & `log_pval` >= 1.3, "A",
                     if_else(`log_fc` <= -1 & `log_pval` >= 1.3,"B","C"))) %>%
    # Plot with points coloured according to the threshold
    ggplot(aes(x=`log_fc`, y=`log_pval`, colour = threshold)) +
    geom_point(alpha = 0.5) + # Alpha sets the transparency of the points
    # Add dotted lines to indicate the threshold, semi-transparent
    geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) + 
    geom_vline(xintercept = c(-1,1), linetype = 2, alpha = 0.5) +
    coord_cartesian(ylim = c(0, 2), xlim = c(-3, 3)) +
    labs(color = 'Legend') +
    # Set the colour of the points
    scale_colour_manual(values = c("A"= "#E2635C", "B"= "#007EA7", 
                                   "C" = "#8F9491"),
                        labels = "Not Significant") +
    xlab("log2(fold change)") + ylab("-log10(adjusted P-value)") + # Relabel the axes
    #ggtitle("ECF mdx vs MTBD baseline DEPs") +
    theme_minimal() + #Set the theme
    theme(legend.position = c(.97, .97),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(6, 6, 6, 6,),
          legend.box.background = element_rect(color="#427681",size=1)) +
    theme(text = element_text(size=14))
}

# Creating volcano plots for 2-group ECF comparisons ----

## Volcano plot for mdx vs WT ECF baseline
volcano_mdxvsWT_baseline_ECF <-
  data.frame('gene ID' = mdxvsWT_baseline_ECF$Gene.Symbol,
             'log_fc' = mdxvsWT_baseline_ECF$L2FC,
             'log_pval' = mdxvsWT_baseline_ECF$X.log10pval,
             check.names = FALSE)
# Add a column to the data frame to specify if they are UP- or DOWN- regulated
volcano_mdxvsWT_baseline_ECF$diffexpressed <- "NO"
  # if log2Foldchange > 1 and -log10pval > 1.3, set as "UP"
volcano_mdxvsWT_baseline_ECF$diffexpressed[volcano_mdxvsWT_baseline_ECF$log_fc > 1 & 
  volcano_mdxvsWT_baseline_ECF$log_pval > 1.3] <- "UP"
  # if log2Foldchange < -1 and -log10pval > 1.3, set as "DOWN"
volcano_mdxvsWT_baseline_ECF$diffexpressed[volcano_mdxvsWT_baseline_ECF$log_fc < -1 & 
  volcano_mdxvsWT_baseline_ECF$log_pval > 1.3] <- "DOWN"
 
volcano_fun(volcano_mdxvsWT_baseline_ECF)

## Volcano plot for mdx vs MTBD ECF baseline
volcano_mdxvsMTBD_baseline_ECF <-
  data.frame('gene ID' = mdxvsMTBD_baseline_ECF$Gene.Symbol,
             'log_fc' = mdxvsMTBD_baseline_ECF$L2FC,
             'log_pval' = mdxvsMTBD_baseline_ECF$X.log10pval,
             check.names = FALSE)
# Add a column to the data frame to specify if they are UP- or DOWN- regulated
volcano_mdxvsMTBD_baseline_ECF$diffexpressed <- "NO"
# if log2Foldchange > 1 and -log10pval > 1.3, set as "UP"
volcano_mdxvsMTBD_baseline_ECF$diffexpressed[volcano_mdxvsMTBD_baseline_ECF$log_fc > 1 & 
                                             volcano_mdxvsMTBD_baseline_ECF$log_pval > 1.3] <- "UP"
# if log2Foldchange < -1 and -log10pval > 1.3, set as "DOWN"
volcano_mdxvsMTBD_baseline_ECF$diffexpressed[volcano_mdxvsMTBD_baseline_ECF$log_fc < -1 & 
                                             volcano_mdxvsMTBD_baseline_ECF$log_pval > 1.3] <- "DOWN"

volcano_fun(volcano_mdxvsMTBD_baseline_ECF)

#write.xlsx(volcano_mdxvsWT_baseline_ECF,file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/volcano_ECF_mdxvsWT_baseline.xlsx",
#           sheetName="mdxvsWT baseline ECF")
#write.xlsx(volcano_mdxvsMTBD_baseline_ECF,file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/volcano_ECF_mdxvsMTBD_baseline.xlsx",
#      sheetName="mdxvsMTBD baseline ECF")

## Volcano plot for MTBD vs WT ECF baseline
volcano_MTBDvsWT_baseline_ECF <-
  data.frame('gene ID' = MTBDvsWT_baseline_ECF$Gene.Symbol,
             'log_fc' = MTBDvsWT_baseline_ECF$L2FC,
             'log_pval' = MTBDvsWT_baseline_ECF$X.log10pval,
             check.names = FALSE)
# Add a column to the data frame to specify if they are UP- or DOWN- regulated
volcano_MTBDvsWT_baseline_ECF$diffexpressed <- "NO"
# if log2Foldchange > 1 and -log10pval > 1.3, set as "UP"
volcano_MTBDvsWT_baseline_ECF$diffexpressed[volcano_MTBDvsWT_baseline_ECF$log_fc > 1 & 
                                             volcano_MTBDvsWT_baseline_ECF$log_pval > 1.3] <- "UP"
# if log2Foldchange < -1 and -log10pval > 1.3, set as "DOWN"
volcano_MTBDvsWT_baseline_ECF$diffexpressed[volcano_MTBDvsWT_baseline_ECF$log_fc < -1 & 
                                             volcano_MTBDvsWT_baseline_ECF$log_pval > 1.3] <- "DOWN"

volcano_fun_2(volcano_MTBDvsWT_baseline_ECF)

## Volcano plot for mdx vs WT ECF post-scruff
volcano_mdxvsWT_scruff_ECF <-
  data.frame('gene ID' = mdxvsWT_scruff_ECF$Gene.Symbol,
             'log_fc' = mdxvsWT_scruff_ECF$L2FC,
             'log_pval' = mdxvsWT_scruff_ECF$X.log10pval,
             check.names = FALSE)
# Add a column to the data frame to specify if they are UP- or DOWN- regulated
volcano_mdxvsWT_scruff_ECF$diffexpressed <- "NO"
# if log2Foldchange > 1 and -log10pval > 1.3, set as "UP"
volcano_mdxvsWT_scruff_ECF$diffexpressed[volcano_mdxvsWT_scruff_ECF$log_fc > 1 & 
                                             volcano_mdxvsWT_scruff_ECF$log_pval > 1.3] <- "UP"
# if log2Foldchange < -1 and -log10pval > 1.3, set as "DOWN"
volcano_mdxvsWT_scruff_ECF$diffexpressed[volcano_mdxvsWT_scruff_ECF$log_fc < -1 & 
                                             volcano_mdxvsWT_scruff_ECF$log_pval > 1.3] <- "DOWN"

volcano_fun(volcano_mdxvsWT_scruff_ECF)

## Volcano plot for mdx vs MTBD ECF post-scruff
volcano_mdxvsMTBD_scruff_ECF <-
  data.frame('gene ID' = mdxvsMTBD_scruff_ECF$Gene.Symbol,
             'log_fc' = mdxvsMTBD_scruff_ECF$L2FC,
             'log_pval' = mdxvsMTBD_scruff_ECF$X.log10pval,
             check.names = FALSE)
# Add a column to the data frame to specify if they are UP- or DOWN- regulated
volcano_mdxvsMTBD_scruff_ECF$diffexpressed <- "NO"
# if log2Foldchange > 1 and -log10pval > 1.3, set as "UP"
volcano_mdxvsMTBD_scruff_ECF$diffexpressed[volcano_mdxvsMTBD_scruff_ECF$log_fc > 1 & 
                                           volcano_mdxvsMTBD_scruff_ECF$log_pval > 1.3] <- "UP"
# if log2Foldchange < -1 and -log10pval > 1.3, set as "DOWN"
volcano_mdxvsMTBD_scruff_ECF$diffexpressed[volcano_mdxvsMTBD_scruff_ECF$log_fc < -1 & 
                                           volcano_mdxvsMTBD_scruff_ECF$log_pval > 1.3] <- "DOWN"

volcano_fun(volcano_mdxvsMTBD_scruff_ECF)

# write.xlsx(volcano_mdxvsWT_scruff_ECF,file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/volcano_ECF_mdxvsWT_scruff.xlsx",
#            sheetName="mdxvsWT scruff ECF")
# write.xlsx(volcano_mdxvsMTBD_scruff_ECF,file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/volcano_ECF_mdxvsMTBD_scruff.xlsx",
#            sheetName="mdxvsMTBD scruff ECF")

## Volcano plot for MTBD vs WT ECF post-scruff
volcano_MTBDvsWT_scruff_ECF <-
  data.frame('gene ID' = MTBDvsWT_scruff_ECF$Gene.Symbol,
             'log_fc' = MTBDvsWT_scruff_ECF$L2FC,
             'log_pval' = MTBDvsWT_scruff_ECF$X.log10pval,
             check.names = FALSE)
# Add a column to the data frame to specify if they are UP- or DOWN- regulated
volcano_MTBDvsWT_scruff_ECF$diffexpressed <- "NO"
# if log2Foldchange > 1 and -log10pval > 1.3, set as "UP"
volcano_MTBDvsWT_scruff_ECF$diffexpressed[volcano_MTBDvsWT_scruff_ECF$log_fc > 1 & 
                                              volcano_MTBDvsWT_scruff_ECF$log_pval > 1.3] <- "UP"
# if log2Foldchange < -1 and -log10pval > 1.3, set as "DOWN"
volcano_MTBDvsWT_scruff_ECF$diffexpressed[volcano_MTBDvsWT_scruff_ECF$log_fc < -1 & 
                                              volcano_MTBDvsWT_scruff_ECF$log_pval > 1.3] <- "DOWN"

volcano_fun_3(volcano_MTBDvsWT_scruff_ECF)

## Volcano plot for mdx scruff vs baseline
volcano_mdxscruffvsbaseline_ECF <-
  data.frame('gene ID' = mdx_scruff_vs_baseline_ECF$Gene.Symbol,
             'log_fc' = mdx_scruff_vs_baseline_ECF$L2FC,
             'log_pval' = mdx_scruff_vs_baseline_ECF$X.log10pval,
             check.names = FALSE)
# Add a column to the data frame to specify if they are UP- or DOWN- regulated
volcano_mdxscruffvsbaseline_ECF$diffexpressed <- "NO"
# if log2Foldchange > 1 and -log10pval > 1.3, set as "UP"
volcano_mdxscruffvsbaseline_ECF$diffexpressed[volcano_mdxscruffvsbaseline_ECF$log_fc > 1 & 
                                            volcano_mdxscruffvsbaseline_ECF$log_pval > 1.3] <- "UP"
# if log2Foldchange < -1 and -log10pval > 1.3, set as "DOWN"
volcano_mdxscruffvsbaseline_ECF$diffexpressed[volcano_mdxscruffvsbaseline_ECF$log_fc < -1 & 
                                            volcano_mdxscruffvsbaseline_ECF$log_pval > 1.3] <- "DOWN"

volcano_fun_4(volcano_mdxscruffvsbaseline_ECF)

## Volcano plot for WT scruff vs baseline
volcano_WTscruffvsbaseline_ECF <-
  data.frame('gene ID' = WT_scruff_vs_baseline_ECF$Gene.Symbol,
             'log_fc' = WT_scruff_vs_baseline_ECF$L2FC,
             'log_pval' = WT_scruff_vs_baseline_ECF$X.log10pval,
             check.names = FALSE)
# Add a column to the data frame to specify if they are UP- or DOWN- regulated
volcano_WTscruffvsbaseline_ECF$diffexpressed <- "NO"
# if log2Foldchange > 1 and -log10pval > 1.3, set as "UP"
volcano_WTscruffvsbaseline_ECF$diffexpressed[volcano_WTscruffvsbaseline_ECF$log_fc > 1 & 
                                                volcano_WTscruffvsbaseline_ECF$log_pval > 1.3] <- "UP"
# if log2Foldchange < -1 and -log10pval > 1.3, set as "DOWN"
volcano_WTscruffvsbaseline_ECF$diffexpressed[volcano_WTscruffvsbaseline_ECF$log_fc < -1 & 
                                                volcano_WTscruffvsbaseline_ECF$log_pval > 1.3] <- "DOWN"

volcano_fun_4(volcano_WTscruffvsbaseline_ECF)

## Volcano plot for MTBD scruff vs baseline
volcano_MTBDscruffvsbaseline_ECF <-
  data.frame('gene ID' = MTBD_scruff_vs_baseline_ECF$Gene.Symbol,
             'log_fc' = MTBD_scruff_vs_baseline_ECF$L2FC,
             'log_pval' = MTBD_scruff_vs_baseline_ECF$X.log10pval,
             check.names = FALSE)
# Add a column to the data frame to specify if they are UP- or DOWN- regulated
volcano_MTBDscruffvsbaseline_ECF$diffexpressed <- "NO"
# if log2Foldchange > 1 and -log10pval > 1.3, set as "UP"
volcano_MTBDscruffvsbaseline_ECF$diffexpressed[volcano_MTBDscruffvsbaseline_ECF$log_fc > 1 & 
                                                volcano_MTBDscruffvsbaseline_ECF$log_pval > 1.3] <- "UP"
# if log2Foldchange < -1 and -log10pval > 1.3, set as "DOWN"
volcano_MTBDscruffvsbaseline_ECF$diffexpressed[volcano_MTBDscruffvsbaseline_ECF$log_fc < -1 & 
                                                volcano_MTBDscruffvsbaseline_ECF$log_pval > 1.3] <- "DOWN"

volcano_fun_4(volcano_MTBDscruffvsbaseline_ECF)
----

  
# Creating volcano plots for 2-group WM comparisons ----

## Volcano plot for mdx vs WT WM baseline
volcano_mdxvsWT_baseline_WM <-
  data.frame('gene ID' = mdxvsWT_baseline_WM$Gene.Symbol,
             'log_fc' = mdxvsWT_baseline_WM$L2FC,
             'log_pval' = mdxvsWT_baseline_WM$X.log10pval,
             check.names = FALSE)
# Add a column to the data frame to specify if they are UP- or DOWN- regulated
volcano_mdxvsWT_baseline_WM$diffexpressed <- "NO"
# if log2Foldchange > 1 and -log10pval > 1.3, set as "UP"
volcano_mdxvsWT_baseline_WM$diffexpressed[volcano_mdxvsWT_baseline_WM$log_fc > 1 & 
                                             volcano_mdxvsWT_baseline_WM$log_pval > 1.3] <- "UP"
# if log2Foldchange < -1 and -log10pval > 1.3, set as "DOWN"
volcano_mdxvsWT_baseline_WM$diffexpressed[volcano_mdxvsWT_baseline_WM$log_fc < -1 & 
                                             volcano_mdxvsWT_baseline_WM$log_pval > 1.3] <- "DOWN"

volcano_fun_5(volcano_mdxvsWT_baseline_WM)

## Volcano plot for mdx vs MTBD WM baseline
volcano_mdxvsMTBD_baseline_WM <-
  data.frame('gene ID' = mdxvsMTBD_baseline_WM$Gene.Symbol,
             'log_fc' = mdxvsMTBD_baseline_WM$L2FC,
             'log_pval' = mdxvsMTBD_baseline_WM$X.log10pval,
             check.names = FALSE)
# Add a column to the data frame to specify if they are UP- or DOWN- regulated
volcano_mdxvsMTBD_baseline_WM$diffexpressed <- "NO"
# if log2Foldchange > 1 and -log10pval > 1.3, set as "UP"
volcano_mdxvsMTBD_baseline_WM$diffexpressed[volcano_mdxvsMTBD_baseline_WM$log_fc > 1 & 
                                               volcano_mdxvsMTBD_baseline_WM$log_pval > 1.3] <- "UP"
# if log2Foldchange < -1 and -log10pval > 1.3, set as "DOWN"
volcano_mdxvsMTBD_baseline_WM$diffexpressed[volcano_mdxvsMTBD_baseline_WM$log_fc < -1 & 
                                               volcano_mdxvsMTBD_baseline_WM$log_pval > 1.3] <- "DOWN"

#write.xlsx(volcano_mdxvsWT_baseline_WM,file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/volcano_WM_mdxvsWT_baseline.xlsx",
#           sheetName="mdxvsWT baseline WM")
#write.xlsx(volcano_mdxvsMTBD_baseline_WM,file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/volcano_WM_mdxvsMTBD_baseline.xlsx",
#      sheetName="mdxvsMTBD baseline WM")

volcano_fun_5(volcano_mdxvsMTBD_baseline_WM)

## Volcano plot for MTBD vs WT WM baseline
volcano_MTBDvsWT_baseline_WM <-
  data.frame('gene ID' = MTBDvsWT_baseline_WM$Gene.Symbol,
             'log_fc' = MTBDvsWT_baseline_WM$L2FC,
             'log_pval' = MTBDvsWT_baseline_WM$X.log10pval,
             check.names = FALSE)
# Add a column to the data frame to specify if they are UP- or DOWN- regulated
volcano_MTBDvsWT_baseline_WM$diffexpressed <- "NO"
# if log2Foldchange > 1 and -log10pval > 1.3, set as "UP"
volcano_MTBDvsWT_baseline_WM$diffexpressed[volcano_MTBDvsWT_baseline_WM$log_fc > 1 & 
                                              volcano_MTBDvsWT_baseline_WM$log_pval > 1.3] <- "UP"
# if log2Foldchange < -1 and -log10pval > 1.3, set as "DOWN"
volcano_MTBDvsWT_baseline_WM$diffexpressed[volcano_MTBDvsWT_baseline_WM$log_fc < -1 & 
                                              volcano_MTBDvsWT_baseline_WM$log_pval > 1.3] <- "DOWN"

volcano_fun_4(volcano_MTBDvsWT_baseline_WM)

## Volcano plot for mdx vs WT WM post-scruff
volcano_mdxvsWT_scruff_WM <-
  data.frame('gene ID' = mdxvsWT_scruff_WM$Gene.Symbol,
             'log_fc' = mdxvsWT_scruff_WM$L2FC,
             'log_pval' = mdxvsWT_scruff_WM$X.log10pval,
             check.names = FALSE)
# Add a column to the data frame to specify if they are UP- or DOWN- regulated
volcano_mdxvsWT_scruff_WM$diffexpressed <- "NO"
# if log2Foldchange > 1 and -log10pval > 1.3, set as "UP"
volcano_mdxvsWT_scruff_WM$diffexpressed[volcano_mdxvsWT_scruff_WM$log_fc > 1 & 
                                           volcano_mdxvsWT_scruff_WM$log_pval > 1.3] <- "UP"
# if log2Foldchange < -1 and -log10pval > 1.3, set as "DOWN"
volcano_mdxvsWT_scruff_WM$diffexpressed[volcano_mdxvsWT_scruff_WM$log_fc < -1 & 
                                           volcano_mdxvsWT_scruff_WM$log_pval > 1.3] <- "DOWN"

volcano_fun_6(volcano_mdxvsWT_scruff_WM)

## Volcano plot for mdx vs MTBD WM post-scruff
volcano_mdxvsMTBD_scruff_WM <-
  data.frame('gene ID' = mdxvsMTBD_scruff_WM$Gene.Symbol,
             'log_fc' = mdxvsMTBD_scruff_WM$L2FC,
             'log_pval' = mdxvsMTBD_scruff_WM$X.log10pval,
             check.names = FALSE)
# Add a column to the data frame to specify if they are UP- or DOWN- regulated
volcano_mdxvsMTBD_scruff_WM$diffexpressed <- "NO"
# if log2Foldchange > 1 and -log10pval > 1.3, set as "UP"
volcano_mdxvsMTBD_scruff_WM$diffexpressed[volcano_mdxvsMTBD_scruff_WM$log_fc > 1 & 
                                             volcano_mdxvsMTBD_scruff_WM$log_pval > 1.3] <- "UP"
# if log2Foldchange < -1 and -log10pval > 1.3, set as "DOWN"
volcano_mdxvsMTBD_scruff_WM$diffexpressed[volcano_mdxvsMTBD_scruff_WM$log_fc < -1 & 
                                             volcano_mdxvsMTBD_scruff_WM$log_pval > 1.3] <- "DOWN"

#write.xlsx(volcano_mdxvsWT_scruff_WM,file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/volcano_WM_mdxvsWT_scruff.xlsx",
#           sheetName="mdxvsWT scruff WM")
#write.xlsx(volcano_mdxvsMTBD_scruff_WM,file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/volcano_WM_mdxvsMTBD_scruff.xlsx",
#      sheetName="mdxvsMTBD scruff WM")

volcano_fun(volcano_mdxvsMTBD_scruff_WM)

## Volcano plot for MTBD vs WT WM post-scruff
volcano_MTBDvsWT_scruff_WM <-
  data.frame('gene ID' = MTBDvsWT_scruff_WM$Gene.Symbol,
             'log_fc' = MTBDvsWT_scruff_WM$L2FC,
             'log_pval' = MTBDvsWT_scruff_WM$X.log10pval,
             check.names = FALSE)
# Add a column to the data frame to specify if they are UP- or DOWN- regulated
volcano_MTBDvsWT_scruff_WM$diffexpressed <- "NO"
# if log2Foldchange > 1 and -log10pval > 1.3, set as "UP"
volcano_MTBDvsWT_scruff_WM$diffexpressed[volcano_MTBDvsWT_scruff_WM$log_fc > 1 & 
                                            volcano_MTBDvsWT_scruff_WM$log_pval > 1.3] <- "UP"
# if log2Foldchange < -1 and -log10pval > 1.3, set as "DOWN"
volcano_MTBDvsWT_scruff_WM$diffexpressed[volcano_MTBDvsWT_scruff_WM$log_fc < -1 & 
                                            volcano_MTBDvsWT_scruff_WM$log_pval > 1.3] <- "DOWN"

volcano_fun(volcano_MTBDvsWT_scruff_WM)

## Volcano plot for mdx scruff vs baseline
volcano_mdxscruffvsbaseline_WM <-
  data.frame('gene ID' = mdx_scruff_vs_baseline_WM$Gene.Symbol,
             'log_fc' = mdx_scruff_vs_baseline_WM$L2FC,
             'log_pval' = mdx_scruff_vs_baseline_WM$X.log10pval,
             check.names = FALSE)
# Add a column to the data frame to specify if they are UP- or DOWN- regulated
volcano_mdxscruffvsbaseline_WM$diffexpressed <- "NO"
# if log2Foldchange > 1 and -log10pval > 1.3, set as "UP"
volcano_mdxscruffvsbaseline_WM$diffexpressed[volcano_mdxscruffvsbaseline_WM$log_fc > 1 & 
                                                volcano_mdxscruffvsbaseline_WM$log_pval > 1.3] <- "UP"
# if log2Foldchange < -1 and -log10pval > 1.3, set as "DOWN"
volcano_mdxscruffvsbaseline_WM$diffexpressed[volcano_mdxscruffvsbaseline_WM$log_fc < -1 & 
                                                volcano_mdxscruffvsbaseline_WM$log_pval > 1.3] <- "DOWN"

volcano_fun_7(volcano_mdxscruffvsbaseline_WM)

## Volcano plot for WT scruff vs baseline
volcano_WTscruffvsbaseline_WM <-
  data.frame('gene ID' = WT_scruff_vs_baseline_WM$Gene.Symbol,
             'log_fc' = WT_scruff_vs_baseline_WM$L2FC,
             'log_pval' = WT_scruff_vs_baseline_WM$X.log10pval,
             check.names = FALSE)
# Add a column to the data frame to specify if they are UP- or DOWN- regulated
volcano_WTscruffvsbaseline_WM$diffexpressed <- "NO"
# if log2Foldchange > 1 and -log10pval > 1.3, set as "UP"
volcano_WTscruffvsbaseline_WM$diffexpressed[volcano_WTscruffvsbaseline_WM$log_fc > 1 & 
                                               volcano_WTscruffvsbaseline_WM$log_pval > 1.3] <- "UP"
# if log2Foldchange < -1 and -log10pval > 1.3, set as "DOWN"
volcano_WTscruffvsbaseline_WM$diffexpressed[volcano_WTscruffvsbaseline_WM$log_fc < -1 & 
                                               volcano_WTscruffvsbaseline_WM$log_pval > 1.3] <- "DOWN"

volcano_fun_7(volcano_WTscruffvsbaseline_WM)

## Volcano plot for MTBD scruff vs baseline
volcano_MTBDscruffvsbaseline_WM <-
  data.frame('gene ID' = MTBD_scruff_vs_baseline_WM$Gene.Symbol,
             'log_fc' = MTBD_scruff_vs_baseline_WM$L2FC,
             'log_pval' = MTBD_scruff_vs_baseline_WM$X.log10pval,
             check.names = FALSE)
# Add a column to the data frame to specify if they are UP- or DOWN- regulated
volcano_MTBDscruffvsbaseline_WM$diffexpressed <- "NO"
# if log2Foldchange > 1 and -log10pval > 1.3, set as "UP"
volcano_MTBDscruffvsbaseline_WM$diffexpressed[volcano_MTBDscruffvsbaseline_WM$log_fc > 1 & 
                                                 volcano_MTBDscruffvsbaseline_WM$log_pval > 1.3] <- "UP"
# if log2Foldchange < -1 and -log10pval > 1.3, set as "DOWN"
volcano_MTBDscruffvsbaseline_WM$diffexpressed[volcano_MTBDscruffvsbaseline_WM$log_fc < -1 & 
                                                 volcano_MTBDscruffvsbaseline_WM$log_pval > 1.3] <- "DOWN"

volcano_fun_7(volcano_MTBDscruffvsbaseline_WM)
----
  
########################################################################
# Subsetting 2-group comparison WM datasets by WM-only proteome filter
# to compare WM-enriched proteins with ECF-enriched proteins
########################################################################

# WM_unique protein candidate datasets

#mdx vs WT baseline
WM_unique_mdxvsWT_baseline <- 
  data.frame('WM_unique_mdxvsWT_baseline' = intersect
             (mdxvsWT_baseline_sig_proteins_WM$'Accession', 
               WM_proteome_only_UniprotID$`Uniprot ID`)) %>% 
  'colnames<-' ("Accession")

#mdx vs WT post-scruff
WM_unique_mdxvsWT_scruff <- 
  data.frame('WM_unique_mdxvsWT_scruff' = intersect
             (mdxvsWT_scruff_sig_proteins_WM$'Accession', 
               WM_proteome_only_UniprotID$`Uniprot ID`)) %>% 
  'colnames<-' ("Accession")

#mdx vs MTBD baseline
WM_unique_mdxvsMTBD_baseline <- 
  data.frame('WM_unique_mdxvsMTBD_baseline' = intersect
             (mdxvsMTBD_baseline_sig_proteins_WM$'Accession', 
               WM_proteome_only_UniprotID$`Uniprot ID`)) %>% 
  'colnames<-' ("Accession")

#mdx vs MTBD post-scruff
WM_unique_mdxvsMTBD_scruff <- 
  data.frame('WM_unique_mdxvsMTBD_scruff' = intersect
             (mdxvsMTBD_scruff_sig_proteins_WM$'Accession', 
               WM_proteome_only_UniprotID$`Uniprot ID`)) %>% 
  'colnames<-' ("Accession")

# No significant proteins in MTBD vs WT baseline or post-scruff datasets, so omitting these

# WM_unique protein candidate non-overlap exported datasets

# mdx vs WT baseline
WM_unique_only_mdxvsWT_baseline <- only_mdxvsWT_baseline_sig_WM %>%
  filter(Accession %in% WM_proteome_only_UniprotID$`Uniprot ID`)

full_WM_unique_only_mdxvsWT_baseline <- WM_proteomics_data %>%
  filter(Accession %in% WM_unique_only_mdxvsWT_baseline$Accession)

#write.xlsx(full_WM_unique_only_mdxvsWT_baseline, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/WM_unique_mdxvsWT_baseline.xlsx",
#           sheetName="full dataset")

# mdx vs WT post-scruff
WM_unique_only_mdxvsWT_scruff <- only_mdxvsWT_scruff_sig_WM %>%
  filter(Accession %in% WM_proteome_only_UniprotID$`Uniprot ID`)

full_WM_unique_only_mdxvsWT_scruff <- WM_proteomics_data %>%
  filter(Accession %in% WM_unique_only_mdxvsWT_scruff$Accession)

#write.xlsx(full_WM_unique_only_mdxvsWT_scruff, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/WM_unique_mdxvsWT_scruff.xlsx",
#           sheetName="full dataset")

# mdx vs MTBD baseline
WM_unique_only_mdxvsMTBD_baseline <- only_mdxvsMTBD_baseline_sig_WM %>%
  filter(Accession %in% WM_proteome_only_UniprotID$`Uniprot ID`)

full_WM_unique_only_mdxvsMTBD_baseline <- WM_proteomics_data %>%
  filter(Accession %in% WM_unique_only_mdxvsMTBD_baseline$Accession)

#write.xlsx(full_WM_unique_only_mdxvsMTBD_baseline, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/WM_unique_mdxvsMTBD_baseline.xlsx",
#           sheetName="full dataset")

# mdx vs MTBD post-scruff
WM_unique_only_mdxvsMTBD_scruff <- only_mdxvsMTBD_scruff_sig_WM %>%
  filter(Accession %in% WM_proteome_only_UniprotID$`Uniprot ID`)

full_WM_unique_only_mdxvsMTBD_scruff <- WM_proteomics_data %>%
  filter(Accession %in% WM_unique_only_mdxvsMTBD_scruff$Accession)

#write.xlsx(full_WM_unique_only_mdxvsMTBD_scruff, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/WM_unique_mdxvsMTBD_scruff.xlsx",
#           sheetName="full dataset")

# Baseline (mutually DEP in mdx vs WT and vs MTBD)

WM_unique_baseline_overlap <- 
  data.frame('WM_unique_baseline_overlap' = intersect
             (WM_unique_mdxvsWT_baseline$Accession,
               WM_unique_mdxvsMTBD_baseline$Accession)) %>% 
  'colnames<-' ("Accession")

full_WM_unique_baseline_overlap <- WM_proteomics_data %>%
  filter(Accession %in% WM_unique_baseline_overlap$Accession)

#write.xlsx(full_WM_unique_baseline_overlap, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/WM_only_baseline_with_scruff_pool.xlsx",
#           sheetName="secreted_not_picky_baseline")

WM_proteome_only_uniquely_overlapping_baseline <-
  WM_unique_baseline_overlap %>% 
  base::subset(!(.$Accession %in% WM_unique_mdxvsWT_scruff$Accession)) %>%
  as.data.frame(.) %>%
  base::subset(!(.$Accession %in% WM_unique_mdxvsMTBD_scruff$Accession)) %>%
  as.data.frame(.)

full_WM_proteome_only_uniquely_overlapping_baseline <- WM_proteomics_data %>%
  filter(Accession %in% WM_proteome_only_uniquely_overlapping_baseline$Accession)

#write.xlsx(full_WM_proteome_only_uniquely_overlapping_baseline, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/WM_proteome_only_uniquely_overlapping_baseline.xlsx",
#           sheetName="full dataset")

# Post-scruff (mutually DEP in mdx vs WT and vs MTBD)

WM_unique_scruff_overlap <- 
  data.frame('WM_unique_scruff_overlap' = intersect
             (WM_unique_mdxvsWT_scruff$Accession,
               WM_unique_mdxvsMTBD_scruff$Accession)) %>% 
  'colnames<-' ("Accession")

full_WM_unique_scruff_overlap <- WM_proteomics_data %>%
  filter(Accession %in% WM_unique_scruff_overlap$Accession)

#write.xlsx(full_WM_unique_scruff_overlap, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/WM_only_scruff_with_scruff_pool.xlsx",
#           sheetName="secreted_not_picky_scruff")

WM_proteome_only_uniquely_overlapping_scruff <-
  WM_unique_scruff_overlap %>% 
  base::subset(!(.$Accession %in% WM_unique_mdxvsWT_baseline$Accession)) %>%
  as.data.frame(.) %>%
  base::subset(!(.$Accession %in% WM_unique_mdxvsMTBD_baseline$Accession)) %>%
  as.data.frame(.)

full_WM_proteome_only_uniquely_overlapping_scruff <- WM_proteomics_data %>%
  filter(Accession %in% WM_proteome_only_uniquely_overlapping_scruff$Accession)

#write.xlsx(full_WM_proteome_only_uniquely_overlapping_scruff, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/WM_proteome_only_uniquely_overlapping_scruff.xlsx",
#           sheetName="full dataset")

# All group mdx DEP overlap in WM only proteome (26 proteins)

WM_only_all_group_overlap <- 
  data.frame('WM_only_all_group_overlap' = intersect
             (WM_unique_baseline_overlap$Accession,
               WM_unique_scruff_overlap$Accession)) %>% 
  'colnames<-' ("Accession")

full_WM_only_all_group_overlap <- WM_proteomics_data %>%
  filter(Accession %in% WM_only_all_group_overlap$Accession)

#write.xlsx(full_WM_only_all_group_overlap, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/WM_only_all_group_overlap.xlsx",
#           sheetName="full dataset")

----

# WM-only protein 4-D Venn comparisons WT, mdx, MTBD ----

WM_unique_4D_venn <- display_venn(WM_unique_4D_venn_list <- list(
  'mdx vs WT baseline' = WM_unique_mdxvsWT_baseline$'Accession',
  'mdx vs MTBD baseline' = WM_unique_mdxvsMTBD_baseline$'Accession',
  'mdx vs WT post-scruff' = WM_unique_mdxvsWT_scruff$'Accession',
  'mdx vs MTBD post-scruff' = WM_unique_mdxvsMTBD_scruff$'Accession'), 
  # Circles
  lwd = 2,
  lty = 'blank',
  category.names = c(
    "mdx vs WT baseline DEP",
    "mdx vs MTBD baseline DEP",
    "mdx vs WT post-scruff DEP", 
    "mdx vs MTBD post-scruff DEP"),
  fill = c("#EFD6AC","#183A37","#C44900", "#432534"),
  # Numbers
  cex = 1.2,
  fontface = "italic",
  # Set names
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-12, 3, -20, 20),
  cat.dist = c(0.12, .12, 0.06, 0.06),
  main = "Quadriceps Whole Muscle Proteome Only Protein Comparisons WT, mdx, and MTBD",
  main.fontface = "bold"
) 

########################################################################
# Subsetting 2-group comparison ECF datasets by OVERLAPPING proteome filter
# to identify potential secreted protein candidates 
# based on differing ECF and WM abundance patterns in mdx
########################################################################

# ECF DEP candidate datasets 
# in overlapping ECF and WM pool

#mdx vs WT baseline
shared_mdxvsWT_baseline_ECF <- 
  data.frame('shared_mdxvsWT_baseline_ECF' = intersect
             (mdxvsWT_baseline_sig_proteins_ECF$'Accession', 
               ECF_WM_proteome_overlap_UniprotID$`Uniprot ID`)) %>% 
  'colnames<-' ("Accession")

#mdx vs WT post-scruff
shared_mdxvsWT_scruff_ECF <- 
  data.frame('shared_mdxvsWT_scruff_ECF' = intersect
             (mdxvsWT_scruff_sig_proteins_ECF$'Accession', 
               ECF_WM_proteome_overlap_UniprotID$`Uniprot ID`)) %>% 
  'colnames<-' ("Accession")

#mdx vs MTBD baseline
shared_mdxvsMTBD_baseline_ECF <- 
  data.frame('shared_mdxvsMTBD_baseline_ECF' = intersect
             (mdxvsMTBD_baseline_sig_proteins_ECF$'Accession', 
               ECF_WM_proteome_overlap_UniprotID$`Uniprot ID`)) %>% 
  'colnames<-' ("Accession")

#mdx vs MTBD post-scruff
shared_mdxvsMTBD_scruff_ECF <- 
  data.frame('shared_mdxvsMTBD_scruff_ECF' = intersect
             (mdxvsMTBD_scruff_sig_proteins_ECF$'Accession', 
               ECF_WM_proteome_overlap_UniprotID$`Uniprot ID`)) %>% 
  'colnames<-' ("Accession")

#MTBD vs WT baseline
shared_MTBDvsWT_baseline_ECF <- 
  data.frame('shared_MTBDvsWT_baseline_ECF' = intersect
             (MTBDvsWT_baseline_sig_proteins_ECF$'Accession', 
               ECF_WM_proteome_overlap_UniprotID$`Uniprot ID`)) %>% 
  'colnames<-' ("Accession")

#MTBD vs WT post-scruff
shared_MTBDvsWT_scruff_ECF <- 
  data.frame('shared_MTBDvsWT_scruff_ECF' = intersect
             (MTBDvsWT_scruff_sig_proteins_ECF$'Accession', 
               ECF_WM_proteome_overlap_UniprotID$`Uniprot ID`)) %>% 
  'colnames<-' ("Accession")

# ECF DEP candidate non-overlap exported datasets

# mdx vs WT baseline
shared_only_mdxvsWT_baseline_ECF <- only_mdxvsWT_baseline_sig_ECF %>%
  filter(Accession %in% ECF_WM_proteome_overlap_UniprotID$`Uniprot ID`)

full_shared_only_mdxvsWT_baseline_ECF <- ECF_proteomics_data %>%
  filter(Accession %in% shared_only_mdxvsWT_baseline_ECF$Accession)

#write.xlsx(full_shared_only_mdxvsWT_baseline_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/ECF_compartment_overlap_mdxvsWT_baseline.xlsx",
#           sheetName="full dataset")

# mdx vs WT post-scruff
shared_only_mdxvsWT_scruff_ECF <- only_mdxvsWT_scruff_sig_ECF %>%
  filter(Accession %in% ECF_WM_proteome_overlap_UniprotID$`Uniprot ID`)

full_shared_only_mdxvsWT_scruff_ECF <- ECF_proteomics_data %>%
  filter(Accession %in% shared_only_mdxvsWT_scruff_ECF$Accession)

#write.xlsx(full_shared_only_mdxvsWT_scruff_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/ECF TMTplex/ECF_compartment_overlap_mdxvsWT_scruff.xlsx",
#           sheetName="full dataset")

# mdx vs MTBD baseline
shared_only_mdxvsMTBD_baseline_ECF <- only_mdxvsMTBD_baseline_sig_ECF %>%
  filter(Accession %in% ECF_WM_proteome_overlap_UniprotID$`Uniprot ID`)

full_shared_only_mdxvsMTBD_baseline_ECF <- ECF_proteomics_data %>%
  filter(Accession %in% shared_only_mdxvsMTBD_baseline_ECF$Accession)

#write.xlsx(full_shared_only_mdxvsMTBD_baseline_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/ECF TMTplex/ECF_compartment_overlap_mdxvsMTBD_baseline.xlsx",
#           sheetName="full dataset")

# mdx vs MTBD post-scruff
shared_only_mdxvsMTBD_scruff_ECF <- only_mdxvsMTBD_scruff_sig_ECF %>%
  filter(Accession %in% ECF_WM_proteome_overlap_UniprotID$`Uniprot ID`)

full_shared_only_mdxvsMTBD_scruff_ECF <- ECF_proteomics_data %>%
  filter(Accession %in% shared_only_mdxvsMTBD_scruff_ECF$Accession)

#write.xlsx(full_secreted_shared_only_mdxvsMTBD_scruff_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/ECF TMTplex/ECF_secreted_compartment_overlap_mdxvsMTBD_scruff.xlsx",
#           sheetName="full dataset")

# Baseline (mutually DEP in mdx vs WT and vs MTBD)
# Scruff DEPs not filtered out yet

shared_baseline_overlap_ECF <- 
  data.frame('shared_baseline_overlap_ECF' = intersect
             (shared_mdxvsWT_baseline_ECF$Accession,
               shared_mdxvsMTBD_baseline_ECF$Accession)) %>% 
  'colnames<-' ("Accession")

full_shared_baseline_overlap_ECF <- ECF_proteomics_data %>%
  filter(Accession %in% shared_baseline_overlap_ECF$Accession)

#write.xlsx(full_shared_baseline_overlap_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/shared_baseline_with_scruff_pool_ECF.xlsx",
#           sheetName="shared_not_picky_baseline")

shared_uniquely_overlapping_baseline_ECF <-
  shared_baseline_overlap_ECF %>% 
  base::subset(!(.$Accession %in% shared_mdxvsWT_scruff_ECF$Accession)) %>%
  as.data.frame(.) %>%
  base::subset(!(.$Accession %in% shared_mdxvsMTBD_scruff_ECF$Accession)) %>%
  as.data.frame(.)

full_shared_uniquely_overlapping_baseline_ECF <- ECF_proteomics_data %>%
  filter(Accession %in% shared_uniquely_overlapping_baseline_ECF$Accession)

#write.xlsx(full_shared_uniquely_overlapping_baseline_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/ECF_compartment_overlap_baseline_mdx_DEPs.xlsx",
#           sheetName="full dataset")

# Post-scruff (mutually DEP in mdx vs WT and vs MTBD)

shared_scruff_overlap_ECF <- 
  data.frame('shared_scruff_overlap_ECF' = intersect
             (shared_mdxvsWT_scruff_ECF$Accession,
               shared_mdxvsMTBD_scruff_ECF$Accession)) %>% 
  'colnames<-' ("Accession")

full_shared_scruff_overlap_ECF <- ECF_proteomics_data %>%
  filter(Accession %in% shared_scruff_overlap_ECF$Accession)

#write.xlsx(full_shared_scruff_overlap_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/shared_scruff_with_baseline_pool_ECF.xlsx",
#           sheetName="shared_not_picky_scruff")

shared_uniquely_overlapping_scruff_ECF <-
  shared_scruff_overlap_ECF %>% 
  base::subset(!(.$Accession %in% shared_mdxvsWT_baseline_ECF$Accession)) %>%
  as.data.frame(.) %>%
  base::subset(!(.$Accession %in% shared_mdxvsMTBD_baseline_ECF$Accession)) %>%
  as.data.frame(.)

full_shared_uniquely_overlapping_scruff_ECF <- ECF_proteomics_data %>%
  filter(Accession %in% shared_uniquely_overlapping_scruff_ECF$Accession)

#write.xlsx(full_shared_uniquely_overlapping_scruff_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/ECF_compartment_overlap_scruff_mdx_DEPs.xlsx",
#           sheetName="full dataset")

# Scruff ECF secreted proteins, both overlapping and unique in mdx vs WT and mdx vs MTBD
# Filtered by overlapping proteome dataset
# Probably won't use this code

accession_shared_only_mdxvsWT_scruff_ECF <-
  data.frame(shared_only_mdxvsWT_scruff_ECF$Accession) %>% 
  'colnames<-' ("Accession")

accession_shared_only_mdxvsMTBD_scruff_ECF <-
  data.frame(shared_only_mdxvsMTBD_scruff_ECF$Accession) %>% 
  'colnames<-' ("Accession")

overlap_and_unique_shared_scruff_ECF <-
  rbind(shared_uniquely_overlapping_scruff_ECF,
        accession_shared_only_mdxvsWT_scruff_ECF,
        accession_shared_only_mdxvsMTBD_scruff_ECF)

#write.csv(overlap_and_unique_shared_scruff_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/shared_overlapping_and_unique_scruff_ECF.csv")

# Common overlap between all mdx DEPs

shared_all_group_overlap_ECF <- 
  data.frame('shared_all_group_overlap_ECF' = intersect
             (shared_baseline_overlap_ECF$Accession,
               shared_scruff_overlap_ECF$Accession)) %>% 
  'colnames<-' ("Accession")

full_shared_all_group_overlap_ECF <- ECF_proteomics_data %>%
  filter(Accession %in% shared_all_group_overlap_ECF$Accession)

#write.xlsx(full_shared_all_group_overlap_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/shared_all_group_overlap_ECF.xlsx",
#           sheetName="full dataset")

# ECF overlapping compartment DEP 4-D Venn comparisons WT, mdx, MTBD ----

ECF_compartment_overlap_4D_venn <- display_venn(ECF_compartment_overlap_4D_venn_list <- 
  list(
  'mdx vs WT baseline' = shared_mdxvsWT_baseline_ECF$'Accession',
  'mdx vs MTBD baseline' = shared_mdxvsMTBD_baseline_ECF$'Accession',
  'mdx vs WT post-scruff' = shared_mdxvsWT_scruff_ECF$'Accession',
  'mdx vs MTBD post-scruff' = shared_mdxvsMTBD_scruff_ECF$'Accession'), 
  # Circles
  lwd = 2,
  lty = 'blank',
  category.names = c(
    "mdx vs WT baseline DEP",
    "mdx vs MTBD baseline DEP",
    "mdx vs WT post-scruff DEP", 
    "mdx vs MTBD post-scruff DEP"),
  fill = c("#1E555C","#3A2E39","#F4D8CD", "#EDB183"),
  # Numbers
  cex = 1.5,
  fontface = "italic",
  # Set names
  cat.cex = 1.3,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-6, 1, -20, 5),
  cat.dist = c(0.14, .14, 0.07, 0.06),
  main = "Quadriceps ECF DEPs found in both compartments",
  main.fontface = "bold"
) 

-----

########################################################################
# Subsetting 2-group comparison ECF datasets by OVERLAPPING proteome filter
# to identify potential secreted protein candidates 
# based on differing ECF and WM abundance patterns in mdx
########################################################################

# ECF DEP candidate datasets 
# in overlapping ECF and WM pool

#mdx vs WT baseline
shared_mdxvsWT_baseline_ECF <- 
  data.frame('shared_mdxvsWT_baseline_ECF' = intersect
             (mdxvsWT_baseline_sig_proteins_ECF$'Accession', 
               ECF_WM_proteome_overlap_UniprotID$`Uniprot ID`)) %>% 
  'colnames<-' ("Accession")

#mdx vs WT post-scruff
shared_mdxvsWT_scruff_ECF <- 
  data.frame('shared_mdxvsWT_scruff_ECF' = intersect
             (mdxvsWT_scruff_sig_proteins_ECF$'Accession', 
               ECF_WM_proteome_overlap_UniprotID$`Uniprot ID`)) %>% 
  'colnames<-' ("Accession")

#mdx vs MTBD baseline
shared_mdxvsMTBD_baseline_ECF <- 
  data.frame('shared_mdxvsMTBD_baseline_ECF' = intersect
             (mdxvsMTBD_baseline_sig_proteins_ECF$'Accession', 
               ECF_WM_proteome_overlap_UniprotID$`Uniprot ID`)) %>% 
  'colnames<-' ("Accession")

#mdx vs MTBD post-scruff
shared_mdxvsMTBD_scruff_ECF <- 
  data.frame('shared_mdxvsMTBD_scruff_ECF' = intersect
             (mdxvsMTBD_scruff_sig_proteins_ECF$'Accession', 
               ECF_WM_proteome_overlap_UniprotID$`Uniprot ID`)) %>% 
  'colnames<-' ("Accession")

#MTBD vs WT baseline
shared_MTBDvsWT_baseline_ECF <- 
  data.frame('shared_MTBDvsWT_baseline_ECF' = intersect
             (MTBDvsWT_baseline_sig_proteins_ECF$'Accession', 
               ECF_WM_proteome_overlap_UniprotID$`Uniprot ID`)) %>% 
  'colnames<-' ("Accession")

#MTBD vs WT post-scruff
shared_MTBDvsWT_scruff_ECF <- 
  data.frame('shared_MTBDvsWT_scruff_ECF' = intersect
             (MTBDvsWT_scruff_sig_proteins_ECF$'Accession', 
               ECF_WM_proteome_overlap_UniprotID$`Uniprot ID`)) %>% 
  'colnames<-' ("Accession")

# ECF DEP candidate non-overlap exported datasets

# mdx vs WT baseline
shared_only_mdxvsWT_baseline_ECF <- only_mdxvsWT_baseline_sig_ECF %>%
  filter(Accession %in% ECF_WM_proteome_overlap_UniprotID$`Uniprot ID`)

full_shared_only_mdxvsWT_baseline_ECF <- ECF_proteomics_data %>%
  filter(Accession %in% shared_only_mdxvsWT_baseline_ECF$Accession)

#write.xlsx(full_shared_only_mdxvsWT_baseline_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/ECF_compartment_overlap_mdxvsWT_baseline.xlsx",
#           sheetName="full dataset")

# mdx vs WT post-scruff
shared_only_mdxvsWT_scruff_ECF <- only_mdxvsWT_scruff_sig_ECF %>%
  filter(Accession %in% ECF_WM_proteome_overlap_UniprotID$`Uniprot ID`)

full_shared_only_mdxvsWT_scruff_ECF <- ECF_proteomics_data %>%
  filter(Accession %in% shared_only_mdxvsWT_scruff_ECF$Accession)

#write.xlsx(full_shared_only_mdxvsWT_scruff_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/ECF TMTplex/ECF_compartment_overlap_mdxvsWT_scruff.xlsx",
#           sheetName="full dataset")

# mdx vs MTBD baseline
shared_only_mdxvsMTBD_baseline_ECF <- only_mdxvsMTBD_baseline_sig_ECF %>%
  filter(Accession %in% ECF_WM_proteome_overlap_UniprotID$`Uniprot ID`)

full_shared_only_mdxvsMTBD_baseline_ECF <- ECF_proteomics_data %>%
  filter(Accession %in% shared_only_mdxvsMTBD_baseline_ECF$Accession)

#write.xlsx(full_shared_only_mdxvsMTBD_baseline_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/ECF TMTplex/ECF_compartment_overlap_mdxvsMTBD_baseline.xlsx",
#           sheetName="full dataset")

# mdx vs MTBD post-scruff
shared_only_mdxvsMTBD_scruff_ECF <- only_mdxvsMTBD_scruff_sig_ECF %>%
  filter(Accession %in% ECF_WM_proteome_overlap_UniprotID$`Uniprot ID`)

full_shared_only_mdxvsMTBD_scruff_ECF <- ECF_proteomics_data %>%
  filter(Accession %in% shared_only_mdxvsMTBD_scruff_ECF$Accession)

#write.xlsx(full_secreted_shared_only_mdxvsMTBD_scruff_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/ECF TMTplex/ECF_secreted_compartment_overlap_mdxvsMTBD_scruff.xlsx",
#           sheetName="full dataset")

# Baseline (mutually DEP in mdx vs WT and vs MTBD)
# Scruff DEPs not filtered out yet

shared_baseline_overlap_ECF <- 
  data.frame('shared_baseline_overlap_ECF' = intersect
             (shared_mdxvsWT_baseline_ECF$Accession,
               shared_mdxvsMTBD_baseline_ECF$Accession)) %>% 
  'colnames<-' ("Accession")

full_shared_baseline_overlap_ECF <- ECF_proteomics_data %>%
  filter(Accession %in% shared_baseline_overlap_ECF$Accession)

#write.xlsx(full_shared_baseline_overlap_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/shared_baseline_with_scruff_pool_ECF.xlsx",
#           sheetName="shared_not_picky_baseline")

shared_uniquely_overlapping_baseline_ECF <-
  shared_baseline_overlap_ECF %>% 
  base::subset(!(.$Accession %in% shared_mdxvsWT_scruff_ECF$Accession)) %>%
  as.data.frame(.) %>%
  base::subset(!(.$Accession %in% shared_mdxvsMTBD_scruff_ECF$Accession)) %>%
  as.data.frame(.)

full_shared_uniquely_overlapping_baseline_ECF <- ECF_proteomics_data %>%
  filter(Accession %in% shared_uniquely_overlapping_baseline_ECF$Accession)

#write.xlsx(full_shared_uniquely_overlapping_baseline_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/ECF_compartment_overlap_baseline_mdx_DEPs.xlsx",
#           sheetName="full dataset")

# Post-scruff (mutually DEP in mdx vs WT and vs MTBD)

shared_scruff_overlap_ECF <- 
  data.frame('shared_scruff_overlap_ECF' = intersect
             (shared_mdxvsWT_scruff_ECF$Accession,
               shared_mdxvsMTBD_scruff_ECF$Accession)) %>% 
  'colnames<-' ("Accession")

full_shared_scruff_overlap_ECF <- ECF_proteomics_data %>%
  filter(Accession %in% shared_scruff_overlap_ECF$Accession)

#write.xlsx(full_shared_scruff_overlap_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/shared_scruff_with_baseline_pool_ECF.xlsx",
#           sheetName="shared_not_picky_scruff")

shared_uniquely_overlapping_scruff_ECF <-
  shared_scruff_overlap_ECF %>% 
  base::subset(!(.$Accession %in% shared_mdxvsWT_baseline_ECF$Accession)) %>%
  as.data.frame(.) %>%
  base::subset(!(.$Accession %in% shared_mdxvsMTBD_baseline_ECF$Accession)) %>%
  as.data.frame(.)

full_shared_uniquely_overlapping_scruff_ECF <- ECF_proteomics_data %>%
  filter(Accession %in% shared_uniquely_overlapping_scruff_ECF$Accession)

#write.xlsx(full_shared_uniquely_overlapping_scruff_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/ECF_compartment_overlap_scruff_mdx_DEPs.xlsx",
#           sheetName="full dataset")

# Scruff ECF secreted proteins, both overlapping and unique in mdx vs WT and mdx vs MTBD
# Filtered by overlapping proteome dataset
# Probably won't use this code

accession_shared_only_mdxvsWT_scruff_ECF <-
  data.frame(shared_only_mdxvsWT_scruff_ECF$Accession) %>% 
  'colnames<-' ("Accession")

accession_shared_only_mdxvsMTBD_scruff_ECF <-
  data.frame(shared_only_mdxvsMTBD_scruff_ECF$Accession) %>% 
  'colnames<-' ("Accession")

overlap_and_unique_shared_scruff_ECF <-
  rbind(shared_uniquely_overlapping_scruff_ECF,
        accession_shared_only_mdxvsWT_scruff_ECF,
        accession_shared_only_mdxvsMTBD_scruff_ECF)

#write.csv(overlap_and_unique_shared_scruff_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/shared_overlapping_and_unique_scruff_ECF.csv")

# Common overlap between all mdx DEPs

shared_all_group_overlap_ECF <- 
  data.frame('shared_all_group_overlap_ECF' = intersect
             (shared_baseline_overlap_ECF$Accession,
               shared_scruff_overlap_ECF$Accession)) %>% 
  'colnames<-' ("Accession")

full_shared_all_group_overlap_ECF <- ECF_proteomics_data %>%
  filter(Accession %in% shared_all_group_overlap_ECF$Accession)

#write.xlsx(full_shared_all_group_overlap_ECF, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/shared_all_group_overlap_ECF.xlsx",
#           sheetName="full dataset")

# ECF DEP overlapping proteome compartments 4-D Venn comparisons WT, mdx, MTBD ----

ECF_compartment_overlap_4D_venn <- display_venn(ECF_compartment_overlap_4D_venn_list <- 
  list(
  'mdx vs WT baseline' = shared_mdxvsWT_baseline_ECF$'Accession',
  'mdx vs MTBD baseline' = shared_mdxvsMTBD_baseline_ECF$'Accession',
  'mdx vs WT post-scruff' = shared_mdxvsWT_scruff_ECF$'Accession',
  'mdx vs MTBD post-scruff' = shared_mdxvsMTBD_scruff_ECF$'Accession'), 
  # Circles
  lwd = 2,
  lty = 'blank',
  category.names = c(
    "mdx vs WT baseline DEP",
    "mdx vs MTBD baseline DEP",
    "mdx vs WT post-scruff DEP", 
    "mdx vs MTBD post-scruff DEP"),
  fill = c("#1E555C","#3A2E39","#F4D8CD", "#EDB183"),
  # Numbers
  cex = 1.5,
  fontface = "italic",
  # Set names
  cat.cex = 1.3,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-6, 1, -20, 5),
  cat.dist = c(0.14, .14, 0.07, 0.06),
  main = "Quadriceps Secreted ECF DEPs found in both compartments",
  main.fontface = "bold"
) 

-----
  
########################################################################
# Subsetting 2-group comparison WM datasets by OVERLAPPING proteome filter
# to identify potential secreted protein candidates 
# based on differing ECF and WM abundance patterns in mdx
########################################################################

# WM DEP candidate datasets 
# in overlapping ECF and WM pool

#mdx vs WT baseline
shared_mdxvsWT_baseline_WM <- 
  data.frame('shared_mdxvsWT_baseline_WM' = intersect
             (mdxvsWT_baseline_sig_proteins_WM$'Accession', 
               ECF_WM_proteome_overlap_UniprotID$`Uniprot ID`)) %>% 
  'colnames<-' ("Accession")

#mdx vs WT post-scruff
shared_mdxvsWT_scruff_WM <- 
  data.frame('shared_mdxvsWT_scruff_WM' = intersect
             (mdxvsWT_scruff_sig_proteins_WM$'Accession', 
               ECF_WM_proteome_overlap_UniprotID$`Uniprot ID`)) %>% 
  'colnames<-' ("Accession")

#mdx vs MTBD baseline
shared_mdxvsMTBD_baseline_WM <- 
  data.frame('shared_mdxvsMTBD_baseline_WM' = intersect
             (mdxvsMTBD_baseline_sig_proteins_WM$'Accession', 
               ECF_WM_proteome_overlap_UniprotID$`Uniprot ID`)) %>% 
  'colnames<-' ("Accession")

#mdx vs MTBD post-scruff
shared_mdxvsMTBD_scruff_WM <- 
  data.frame('shared_mdxvsMTBD_scruff_WM' = intersect
             (mdxvsMTBD_scruff_sig_proteins_WM$'Accession', 
               ECF_WM_proteome_overlap_UniprotID$`Uniprot ID`)) %>% 
  'colnames<-' ("Accession")

# WM DEP candidate non-overlap exported datasets

# mdx vs WT baseline
shared_only_mdxvsWT_baseline_WM <- only_mdxvsWT_baseline_sig_WM %>%
  filter(Accession %in% ECF_WM_proteome_overlap_UniprotID$`Uniprot ID`)

full_shared_only_mdxvsWT_baseline_WM <- WM_proteomics_data %>%
  filter(Accession %in% shared_only_mdxvsWT_baseline_WM$Accession)

#write.xlsx(full_shared_only_mdxvsWT_baseline_WM, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/WM_compartment_overlap_mdxvsWT_baseline.xlsx",
#           sheetName="full dataset")

# mdx vs WT post-scruff
shared_only_mdxvsWT_scruff_WM <- only_mdxvsWT_scruff_sig_WM %>%
  filter(Accession %in% ECF_WM_proteome_overlap_UniprotID$`Uniprot ID`)

full_shared_only_mdxvsWT_scruff_WM <- WM_proteomics_data %>%
  filter(Accession %in% shared_only_mdxvsWT_scruff_WM$Accession)

#write.xlsx(full_shared_only_mdxvsWT_scruff_WM, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/WM TMTplex/WM_compartment_overlap_mdxvsWT_scruff.xlsx",
#           sheetName="full dataset")

# mdx vs MTBD baseline
shared_only_mdxvsMTBD_baseline_WM <- only_mdxvsMTBD_baseline_sig_WM %>%
  filter(Accession %in% ECF_WM_proteome_overlap_UniprotID$`Uniprot ID`)

full_shared_only_mdxvsMTBD_baseline_WM <- WM_proteomics_data %>%
  filter(Accession %in% shared_only_mdxvsMTBD_baseline_WM$Accession)

#write.xlsx(full_shared_only_mdxvsMTBD_baseline_WM, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/WM TMTplex/WM_compartment_overlap_mdxvsMTBD_baseline.xlsx",
#           sheetName="full dataset")

# mdx vs MTBD post-scruff
shared_only_mdxvsMTBD_scruff_WM <- only_mdxvsMTBD_scruff_sig_WM %>%
  filter(Accession %in% ECF_WM_proteome_overlap_UniprotID$`Uniprot ID`)

full_shared_only_mdxvsMTBD_scruff_WM <- WM_proteomics_data %>%
  filter(Accession %in% shared_only_mdxvsMTBD_scruff_WM$Accession)

#write.xlsx(full_secreted_shared_only_mdxvsMTBD_scruff_WM, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/WM TMTplex/WM_secreted_compartment_overlap_mdxvsMTBD_scruff.xlsx",
#           sheetName="full dataset")

# Baseline (mutually DEP in mdx vs WT and vs MTBD)
# Scruff DEPs not filtered out yet

shared_baseline_overlap_WM <- 
  data.frame('shared_baseline_overlap_WM' = intersect
             (shared_mdxvsWT_baseline_WM$Accession,
               shared_mdxvsMTBD_baseline_WM$Accession)) %>% 
  'colnames<-' ("Accession")

full_shared_baseline_overlap_WM <- WM_proteomics_data %>%
  filter(Accession %in% shared_baseline_overlap_WM$Accession)

#write.xlsx(full_shared_baseline_overlap_WM, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/shared_baseline_with_scruff_pool_WM.xlsx",
#           sheetName="shared_not_picky_baseline")

shared_uniquely_overlapping_baseline_WM <-
  shared_baseline_overlap_WM %>% 
  base::subset(!(.$Accession %in% shared_mdxvsWT_scruff_WM$Accession)) %>%
  as.data.frame(.) %>%
  base::subset(!(.$Accession %in% shared_mdxvsMTBD_scruff_WM$Accession)) %>%
  as.data.frame(.)

full_shared_uniquely_overlapping_baseline_WM <- WM_proteomics_data %>%
  filter(Accession %in% shared_uniquely_overlapping_baseline_WM$Accession)

#write.xlsx(full_shared_uniquely_overlapping_baseline_WM, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/WM_compartment_overlap_baseline_mdx_DEPs.xlsx",
#           sheetName="full dataset")

# Post-scruff (mutually DEP in mdx vs WT and vs MTBD)

shared_scruff_overlap_WM <- 
  data.frame('shared_scruff_overlap_WM' = intersect
             (shared_mdxvsWT_scruff_WM$Accession,
               shared_mdxvsMTBD_scruff_WM$Accession)) %>% 
  'colnames<-' ("Accession")

full_shared_scruff_overlap_WM <- WM_proteomics_data %>%
  filter(Accession %in% shared_scruff_overlap_WM$Accession)

#write.xlsx(full_shared_scruff_overlap_WM, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/shared_scruff_with_baseline_pool_WM.xlsx",
#           sheetName="shared_not_picky_scruff")

shared_uniquely_overlapping_scruff_WM <-
  shared_scruff_overlap_WM %>% 
  base::subset(!(.$Accession %in% shared_mdxvsWT_baseline_WM$Accession)) %>%
  as.data.frame(.) %>%
  base::subset(!(.$Accession %in% shared_mdxvsMTBD_baseline_WM$Accession)) %>%
  as.data.frame(.)

full_shared_uniquely_overlapping_scruff_WM <- WM_proteomics_data %>%
  filter(Accession %in% shared_uniquely_overlapping_scruff_WM$Accession)

#write.xlsx(full_shared_uniquely_overlapping_scruff_WM, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/WM_compartment_overlap_scruff_mdx_DEPs.xlsx",
#           sheetName="full dataset")

# Scruff WM secreted proteins, both overlapping and unique in mdx vs WT and mdx vs MTBD
# Filtered by overlapping proteome dataset
# Probably won't use this code

accession_shared_only_mdxvsWT_scruff_WM <-
  data.frame(shared_only_mdxvsWT_scruff_WM$Accession) %>% 
  'colnames<-' ("Accession")

accession_shared_only_mdxvsMTBD_scruff_WM <-
  data.frame(shared_only_mdxvsMTBD_scruff_WM$Accession) %>% 
  'colnames<-' ("Accession")

overlap_and_unique_shared_scruff_WM <-
  rbind(shared_uniquely_overlapping_scruff_WM,
        accession_shared_only_mdxvsWT_scruff_WM,
        accession_shared_only_mdxvsMTBD_scruff_WM)

#write.csv(overlap_and_unique_shared_scruff_WM, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/shared_overlapping_and_unique_scruff_WM.csv")

# Common overlap between all mdx DEPs

shared_all_group_overlap_WM <- 
  data.frame('shared_all_group_overlap_WM' = intersect
             (shared_baseline_overlap_WM$Accession,
               shared_scruff_overlap_WM$Accession)) %>% 
  'colnames<-' ("Accession")

full_shared_all_group_overlap_WM <- WM_proteomics_data %>%
  filter(Accession %in% shared_all_group_overlap_WM$Accession)

#write.xlsx(full_shared_all_group_overlap_WM, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/shared_all_group_overlap_WM.xlsx",
#           sheetName="full dataset")

# WM DEP overlapping proteome compartments 4-D Venn comparisons WT, mdx, MTBD ----

WM_compartment_overlap_4D_venn <- display_venn(WM_compartment_overlap_4D_venn_list <- 
  list(
    'mdx vs WT baseline' = shared_mdxvsWT_baseline_WM$'Accession',
    'mdx vs MTBD baseline' = shared_mdxvsMTBD_baseline_WM$'Accession',
    'mdx vs WT post-scruff' = shared_mdxvsWT_scruff_WM$'Accession',
    'mdx vs MTBD post-scruff' = shared_mdxvsMTBD_scruff_WM$'Accession'
      ),
  # Circles
  lwd = 2,
  lty = 'blank',
  category.names = c(
  "mdx vs WT baseline DEP",
  "mdx vs MTBD baseline DEP",
  "mdx vs WT post-scruff DEP",
  "mdx vs MTBD post-scruff DEP"
      ),
  fill = c("#5C6D70", "#484A47", "#E88873", "#A37774"),
  # Numbers
  cex = 1.5,
  fontface = "italic",
  # Set names
  cat.cex = 1.3,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-6, 1,-20, 5),
  cat.dist = c(0.14, .14, 0.07, 0.06),
  main = "Quadriceps Muscle DEPs found in both compartments",
  main.fontface = "bold"
) 
    
#####################
# PCA Plots
#####################
# ECF PCA plot
#####################

ECF_prot_transposed <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/Data files for R/ECF_proteomics_data_transposed_PCA.xlsx"))

ECF_prot_PCA <- ECF_prot_transposed[2:1571]

ECF_prot_PCA_no_missing <- ECF_prot_PCA[ , colSums(is.na(ECF_prot_PCA))==0]

# Dplyr method: ECF_prot_PCA_no_missing <- ECF_prot_PCA %>% select_if(~ !any(is.na(.)))
pca_res_ECF <- prcomp(ECF_prot_PCA_no_missing, scale. = TRUE)

# For ellipses
#viz_pca_ind(pca_res, 
#            habillage=ECF_prot_transposed$Group,
#            addEllipses=TRUE)

PCi_ECF<-data.frame(pca_res_ECF$x,Legend=ECF_prot_transposed$Legend)

# ECF PCA plot without eigenvectors (not a biplot)
ggplot(PCi_ECF,aes(x=PC1,y=PC2,col=Legend,frame=T))+
  geom_point(size=3,alpha=0.6) + #Size and alpha just for fun
                       theme_minimal() +
                       theme(text = element_text(size=15)) +
  stat_ellipse(alpha=0.02,geom = "polygon", aes(fill=Legend)) +
  scale_color_manual(values = c("#414487","#7ad151","#2a788e","#fde725","#440154","#22a884"))

#PCA<-autoplot(pca_res, data = ECF_prot_transposed, colour = 'Group',frame=T)

## Other PCA plot option using corrr and FactoMineR ##

corr_matrix_ECF <- cor(ECF_prot_PCA_no_missing)
#ggcorrplot(corr_matrix_ECF)

data_normalized_ECF <- scale(ECF_prot_PCA_no_missing)
head(data_normalized_ECF)

# For some reason, can't use princomp like example showed
# Something about needing more units than variables and R-mode vs Q-mode
data.pca_ECF <- prcomp(data_normalized_ECF)
summary(data.pca_ECF)

# Scree plot
fviz_eig(data.pca_ECF, addlabels = TRUE)

# Biplot
fviz_pca_var(data.pca_ECF, col.var = "black", select.var = list(contrib=5))

fviz_pca_ind(data.pca_ECF, col.ind = c("#414487","#414487","#414487","#414487","#414487",
                                       "#7ad151","#7ad151","#7ad151","#7ad151","#7ad151",
                                       "#2a788e","#2a788e","#2a788e","#2a788e","#2a788e",
                                       "#fde725","#fde725","#fde725","#fde725","#fde725",
                                       "#440154","#440154","#440154","#440154","#440154",
                                       "#22a884","#22a884","#22a884","#22a884","#22a884"),
             addEllipses = TRUE)

cos2_ECF50 <- fviz_cos2(data.pca_ECF, choice = "var", axes = 1:2, 
                      sort.val = "desc",
                      top=50) # Can change top to any numerical value
cos2_ECF50_PC1only <- fviz_cos2(data.pca_ECF, choice = "var", axes = 1, 
                        sort.val = "desc",
                        top=50) # Can change top to any numerical value

top_features_ECF <- as.data.frame(cos2_ECF50$data)
top_features_ECF_PC1 <- as.data.frame(cos2_ECF50_PC1only$data)

#write.xlsx(top_features_ECF, 
#           "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/ECF_PC1_PC2_top_features.xlsx")
#write.xlsx(top_features_ECF_PC1, 
#           "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/ECF TMTplex/ECF_PC1_only_top_features.xlsx")

fviz_pca_var(data.pca_ECF, col.var = "cos2",
             gradient.cols = c("black", "lightblue", "purple"),
             ggtheme = theme_minimal(),
             repel = TRUE, 
             select.var = list(name = NULL, cos2 = NULL, contrib = 20))

# Most useful and visually appealing ECF biplot
# Can change contrib number to add or subtract vectors
fviz_pca_biplot(data.pca_ECF, axes = c(1,2), label = "var", repel = TRUE,
                col.var = "black",
                geom = "point",
                pointsize = 2,
                habillage=ECF_prot_transposed$Legend,
                palette = c("#414487","#7ad151","#2a788e","#fde725","#440154","#22a884"), 
                addEllipses=TRUE, ellipse.level=0.95,
                ellipse.type = "t",
                select.var = list(contrib = 10),
                ggtheme = theme_minimal())

----

######################
# WM PCA Plot
######################

WM_prot_transposed <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/Data files for R/WM_prot_transposed.xlsx"))

WM_prot_PCA <- WM_prot_transposed[2:2408]
WM_prot_PCA_no_missing <- WM_prot_PCA[ , colSums(is.na(WM_prot_PCA))==0]
pca_res_WM <- prcomp(WM_prot_PCA_no_missing, scale. = TRUE)


# WM PCA plot without eigenvectors (not a biplot)
PCi_WM<-data.frame(pca_res_WM$x,Legend=WM_prot_transposed$Legend)

# WM PCA plot without eigenvectors (not a biplot)
ggplot(PCi_WM,aes(x=PC1,y=PC2,col=Legend,frame=T))+
  geom_point(size=3,alpha=0.6) + #Size and alpha just for fun
  theme_minimal() +
  theme(text = element_text(size=15)) +
  stat_ellipse(alpha=0.02,geom = "polygon", aes(fill=Legend)) +
  scale_color_manual(values = c("#414487","#7ad151","#2a788e","#fde725","#440154","#22a884"))

data_normalized_WM <- scale(WM_prot_PCA_no_missing)
  
data.pca_WM <- prcomp(data_normalized_WM) 

# Scree plot
fviz_eig(data.pca_WM, addlabels = TRUE)

# Biplot
fviz_pca_var(data.pca_WM, col.var = "black", select.var = list(contrib=5))

fviz_pca_ind(data.pca_WM, col.ind = c("#414487","#414487","#414487","#414487","#414487",
                                       "#7ad151","#7ad151","#7ad151","#7ad151","#7ad151",
                                       "#2a788e","#2a788e","#2a788e","#2a788e","#2a788e",
                                       "#fde725","#fde725","#fde725","#fde725","#fde725",
                                       "#440154","#440154","#440154","#440154","#440154",
                                       "#22a884","#22a884","#22a884","#22a884","#22a884"),
             addEllipses = TRUE)

cos2_WM50 <- 
fviz_cos2(data.pca_WM, choice = "var", axes = 1:2, 
                        sort.val = "desc",
                        top=50) # Can change top to any numerical value
cos2_WM50_PC1only <- 
fviz_cos2(data.pca_WM, choice = "var", axes = 1, 
                                sort.val = "desc",
                                top=50) # Can change top to any numerical value

top_features_WM <- as.data.frame(cos2_WM50$data)
top_features_WM_PC1 <- as.data.frame(cos2_WM50_PC1only$data)

#write.xlsx(top_features_WM, 
#           "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/WM_PC1_PC2_top_features.xlsx")
#write.xlsx(top_features_WM_PC1, 
#           "C:/Users/JOH18358/OneDrive/Ervasti Lab/SkM EF proteomics/WM TMTplex/WM_PC1_only_top_features.xlsx")

fviz_pca_var(data.pca_WM, col.var = "cos2",
             gradient.cols = c("black", "lightblue", "purple"),
             ggtheme = theme_minimal(),
             repel = TRUE, 
             select.var = list(name = NULL, cos2 = NULL, contrib = 20))

# Most useful and visually appealing WM biplot
# Can change contrib number to add or subtract vectors
fviz_pca_biplot(data.pca_WM, axes = c(1,2), label = "var", repel = TRUE,
                col.var = "black",
                geom = "point",
                pointsize = 2,
                habillage=WM_prot_transposed$Legend,
                palette = c("#414487","#7ad151","#2a788e","#fde725","#440154","#22a884"), 
                addEllipses=TRUE, ellipse.level=0.95,
                ellipse.type = "t",
                select.var = list(contrib = 5),
                ggtheme = theme_minimal())
