{
  library(ggplot2)
  library(RColorBrewer)
  library(lattice)
  library(reshape2)
  library(pheatmap)
  library(plyr)
  library(knitr)
  library(circlize)
  library(FactoMineR)
  library("factoextra")
  library("corrplot")
  library(readxl)
  library(ggplot2)
  library(lubridate)
  library(dplyr)
  library(tidyr)
  library(here)
  library(stringr)
  library(tidyverse)
  library(writexl)
  library(dendextend) #comparaison dendogrammes
}

#dendogramme expression
{
  setwd("/Users/pgermon/Documents/Divers_labo/Publis/2024_Gitton/Soumission_PCI_revised/scripts")
  
  table_analyse <- read_xlsx("table_analyse_complete_2024_01_30.xlsx")
  clusters_pop_sep <- read.csv("Cluster_pop_separees.csv", header=T)
  
  matrix_pca_expression <- table_analyse[,!is.element(colnames(table_analyse),
                                               c("Quartier","Name", "Date",
                                                 "Volume","ARN","Date_prel","Description",
                                                 "Echantillon","Type","duplicats",
                                                 "ACTB","GAPDH",
                                                 "Extraction","Macrophage","Lymphocyte","Neutrophile", "CCS",
                                                 "Total", "Serie",
                                                 "Macro_num", 
                                                 "Lympho_num", "Neutro_num",                                      
                                                 "Race", "SCS_cate",                                         
                                                 "SCS", "SCS_cate_2",                                       
                                                 "Cluster_3_samples", "Cluster_4_samples", "Bacterio", "Bacterio_pathoproof",                                         
                                                 "Staphylococcus sp.", "Str. uberis",                                      
                                                 "Str. dysgalactiae", "Beta-lactamase gene",                              
                                                 "E. coli", "C. bovis",                                         
                                                 "Enterococcus sp. (including faecalis and faecium)", "Klebsiella sp. (including oxytoca and pneumonia)", 
                                                 "M. bovis", "Mycoplasma sp.",                                   
                                                 "Prototheca sp.", "S. marcescens",                                    
                                                 "Staph. aureus", "Str. agalactiae",                                  
                                                 "T. pyogenes and P. indolicus", "Yeasts",                                         
                                                 "Quartier", "Date",                                             
                                                 "Volume", "Echantillon",                                      
                                                 "Type", "Extraction",                                       
                                                 "Macrophage", "Lymphocyte",                                       
                                                 "Neutrophile", "CCS",                                              
                                                 "Total", "Serie", "jours_lait"
                                               ))]
  matrix_pca_expression <- matrix_pca_expression[,c(is.element(colnames(matrix_pca_expression),clusters_pop_sep$Gene[is.element(clusters_pop_sep$Nom,c("Macro","Lc", "PMN2"))]))]
  rownames(matrix_pca_expression) <- table_analyse$Quartier
  
    
    myBreaks <- c(seq(-3, 3,by=0.1))
    paletteLength<-length(myBreaks)
    myColor <- colorRampPalette(c("blue","white","red"))(paletteLength)
    
    ann_colors = list(
      type=c(
        "tri"="#66A61E",
        "centri"="#CCEBC5",
        "lait"="purple",
        "Lympho"="blue",
        "Macro"="green",
        "Neutro"="red"),
      Cluster_genes = c("Lc"="blue",
                        "Macro"="green",
                        "PMN1"="red",
                        "Macro_PMN1"="yellow",
                        "Macro_PMN2"="yellow",
                        "PMN2"="purple"
      ),
      Cluster_samples = c("1" = "blue",
                          "2" = "cyan",
                          "3" = "burlywood1",
                          "4" = "red"
                          # "5" = "yellow",
                          # "6" = "black",
                          # "7" = "purple"
      )
    )
    
    
    annotation_row<-data.frame(
      Cluster_samples=factor(table_analyse$Cluster_4_samples),
      CCS=log2(table_analyse$CCS /100)+3
    )
    
    rownames(annotation_row) <- table_analyse$Quartier
    
    # annotation_row$Names <- rownames(annotation_row)<-paste0(table_analyse$Quartier,"_",table_analyse$Type, "_", table_analyse$Name)
    # clusters_samples <- read.csv("Clusters_4_samples_separees.csv", header = T)[,2:3]
    # annotation_row <- merge(annotation_row, clusters_samples, by="Names", all = T)
    # rownames(annotation_row) <- annotation_row$Names
    # annotation_row$Names <- NULL
    
    # annotation_row$Cluster_samples[annotation_row$type == "Macro"] <- 5
    # annotation_row$Cluster_samples[annotation_row$type == "Neutro"] <- 6
    # annotation_row$Cluster_samples[annotation_row$type == "Lympho"] <- 7
    
    
    annotation_col <- data.frame(
      Cluster_genes = factor(clusters_pop_sep$Nom)
    )
    rownames(annotation_col) <- clusters_pop_sep$Gene
    
    labs.row<-rownames(matrix_pca_expression)
    
    
    titre<-paste0("fc - gènes ")
    
    
    matrix_pca_expression <- scale(matrix_pca_expression)
    
    p_expression <- pheatmap(matrix_pca_expression,
                  # cluster_rows = F,
                  cutree_cols = 3,
                  cutree_rows = 4,
                  clustering_method="ward.D2",
                  clustering_distance_cols ="correlation",
                  clustering_distance_rows ="correlation",
                  annotation_colors = ann_colors,
                  annotation_row = annotation_row,
                  annotation_col = annotation_col,
                  breaks=myBreaks,
                  legend_breaks=c(-3,0,3),
                  
                  fontsize_row=8,
                  color=myColor,
                  # labels_row = labs.row,
                  main=titre
                  
    )
    
  
  
  dendro_expression <- p_expression[[1]] %>%
    as.dendrogram()
}

#dendogramme cytometrie
{
  matrix_pca_cytometrie <- as.matrix(table_analyse[,c("Total","CCS", "Macro_num","Lympho_num", "Neutro_num")])
  rownames(matrix_pca_cytometrie) <- table_analyse$Quartier
  
  res.pca_cytometrie <- PCA(log(matrix_pca_cytometrie[,3:5]),
                          graph = F)
  res.pca_cytometrie$ind$coord
  p_cytometrie <- pheatmap(res.pca_cytometrie$ind$coord,
                           cutree_rows = 4,
                           clustering_method="ward.D",
                           clustering_distance_cols ="euclidean",
                           clustering_distance_rows ="euclidean")
  dendro_cytometrie <- p_cytometrie[[1]] %>%
    as.dendrogram()
  
}
#création répertoire pour figures ####
{
  date_jour<-as.character(Sys.Date())
  
  mainDir<-c("/Users/pgermon/Documents/Divers_labo/Projets_2_suivi/2019_07_MASTICELLS_DSCC/2023_lait_Pin_Haras_HF_NO/2023_03_28_Fluidigm")
  subDir<-paste0("Figures_comparaison_dendrogrammes_",date_jour)
  wd_figures<-file.path(mainDir, subDir)
  dir.create(wd_figures)
  setwd(wd_figures)
  
}

# Compute distance matrix####

#voir http://talgalili.github.io/dendextend/articles/dendextend.html
# Create a list to hold dendrograms
dend_list <- dendlist(dendro_expression, dendro_cytometrie)

# Align and plot two dendrograms side by side
set.seed(1)
pdf("Figure_XX_tanglegram.pdf",12,12)

dendlist(dendro_expression, dendro_cytometrie) %>%
  untangle(method = "step2side") %>% # Find the best alignment layout
  tanglegram(main_left = "dendrogramme expression",
             main_right = "dendrogramme cytometrie",
             columns_width = c(5,3,5),
             margin_inner = 8,
             cex_main = 1.5)



dev.off()

cor_bakers_gamma(dendro_expression, dendro_cytometrie)
#Baker’s Gamma Index
# 
# Baker’s Gamma Index (see baker’s paper from 1974) is a measure of association (similarity) between two trees of Hierarchical clustering (dendrograms).
# It is defined as the rank correlation between the stages at which pairs of objects combine in each of the two trees.
# 
# Or more detailed: It is calculated by taking two items, and see what is the highest possible level of k (number of cluster groups created when cutting the tree) for which the two item still belongs to the same tree.
# That k is returned, and the same is done for these two items for the second tree.
# There are n over 2 combinations of such pairs of items from the items in the tree, and all of these numbers are calculated for each of the two trees.
# Then, these two sets of numbers (a set for the items in each tree) are paired according to the pairs of items compared, and a Spearman correlation is calculated.
# 
# The value can range between -1 to 1. With near 0 values meaning that the two trees are not statistically similar.
# For exact p-value one should use a permutation test. One such option will be to permute over the labels of one tree many times, calculating the distribution under the null hypothesis (keeping the trees topologies constant).
# 
# Notice that this measure is not affected by the height of a branch but only of its relative position compared with other branches.


set.seed(23235)
the_cor <- cor_bakers_gamma(dendro_expression, dendro_expression)
the_cor2 <- cor_bakers_gamma(dendro_expression, dendro_cytometrie)
the_cor
the_cor2


R <- 1000
cor_bakers_gamma_results <- numeric(R)
dend_mixed <- dendro_expression
for(i in 1:R) {
  dend_mixed <- sample.dendrogram(dend_mixed, replace = FALSE)
  cor_bakers_gamma_results[i] <- cor_bakers_gamma(dendro_expression, dend_mixed)
}

pdf("Figure_XX_Baker_index.pdf",12,12)
plot(density(cor_bakers_gamma_results),
     main = "Baker's gamma distribution under H0",
     xlim = c(-1,1))
abline(v = 0, lty = 2)
# abline(v = the_cor, lty = 2, col = 2)
abline(v = the_cor2, lty = 2, col = 4)
legend("topleft", legend = c("Baker's index"), fill = c(4))
round(sum(the_cor2 < cor_bakers_gamma_results)/ R, 4)

title(sub = paste("One sided p-value:", round(sum(the_cor2 < cor_bakers_gamma_results)/ R, 4))
)

dev.off()

