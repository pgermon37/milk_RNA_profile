#### chargement des libraries####
{
  library(ggplot2)
  library(RColorBrewer)
  library(lattice)
  library(reshape2)
  library(pheatmap)
  library(plyr)
  library(knitr)
  library(gplots)
  library(circlize)
  library(FactoMineR)
  library(factoextra)
  library(corrplot)
  library(readxl)
  library(ggplot2)
  library(lubridate)
  library(dplyr)
  library(tidyr)
  library(here)
  library(stringr)
  library(tidyverse)
  library(writexl)
  library(psych)
}

####Load tables####
{
  setwd("/Users/pgermon/Documents/Divers_labo/Publis/2024_Gitton/Soumission_PCI_revised/scripts")
  
  table_analyse <- read_xlsx("table_analysis_complete.xlsx")
}

####_####


#analysis by cell types####

#PCA on cell types####
{
  matrix_heatmap <- table_analyse[,!is.element(colnames(table_analyse),
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
                                                "Cluster", "Bacterio", "Bacterio_pathoproof",                                         
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

matrix_pca <- as.matrix(table_analyse[,c("Total","CCS", "Macro_num","Lympho_num", "Neutro_num", "Cluster_4_samples")])

rownames(matrix_pca) <- paste0(table_analyse$Quartier,"_",table_analyse$Type, "_", table_analyse$Name)


library(FactoMineR)

res.pca <- PCA(log(matrix_pca[,3:5]),
               graph = F)

eig.val <- get_eigenvalue(res.pca)

fviz_eig(res.pca, addlabels = TRUE, ylim = c(0,80))
ind <- get_pca_ind(res.pca)

# pdf("Figure_XX_PCA_cytometry_w_SCS_pour_publi.pdf",12,12)
tiff("Figure_XX_PCA_cytometry_w_SCS_pour_publi.tiff",width = 600, height = 600, units = "px")
color <- log2((matrix_pca[,2] *1000) / 100000) +3
fviz_pca_ind(res.pca,
             axes = c(1,2),
             # pointshape = as.factor(table_analyse$Cluster_3_samples),
             pointsize = 5,
             # col.ind = "cos2",
             col.ind=color,
             # fill.ind=color,
             gradient.cols = c("green","blue","burlywood1","yellow"),
             # palette = c("blue","green","red"),
             # palette = c("blue","cyan","burlywood1","red"),
             repel = TRUE, # Avoid text overlapping (slow if many points)
             # addEllipses = T,
             ellipse.type = "confidence",
             ellipse.level=0.99,
             mean.point = FALSE,
             label = "none",
             title=paste0(" "),
             legend.title = "SCS") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, size=10),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.title.x.top = element_text(size=25),
        strip.text = element_text(size=15),
        legend.title = element_text(size=20),
        legend.text = element_text(size=15)
  )
dev.off()

# pdf("Figure_XX_PCA_cytometry_w_cluster_pour publi.pdf",12,12)
tiff("Figure_XX_PCA_cytometry_w_cluster_pour publi.tiff",width = 600, height = 600, units = "px")
color <- as.factor(matrix_pca[,6])
fviz_pca_ind(res.pca,
             axes = c(1,2),
             # pointshape = as.factor(table_analyse$Cluster_3_samples),
             pointsize = 5,
             # col.ind = "cos2",
             col.ind=color,
             # fill.ind=color,
             # gradient.cols = c("green","blue","burlywood1","yellow"),
             # palette = c("blue","green","red"),
             palette = c("blue","cyan","burlywood1","red"),
             repel = TRUE, # Avoid text overlapping (slow if many points)
             addEllipses = T,
             ellipse.type = "confidence",
             ellipse.level=0.99,
             mean.point = FALSE,
             label = "none",
             title = paste0(" "),
             legend.title = "Cluster") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, size=10),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.title.x.top = element_text(size=25),
        strip.text = element_text(size=15),
        legend.title = element_text(size=20),
        legend.text = element_text(size=15)
  )
dev.off()

}
####fin PCA####



