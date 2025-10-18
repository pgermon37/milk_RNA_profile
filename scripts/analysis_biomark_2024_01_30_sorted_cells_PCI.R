#libraries####
{
library(ggplot2)
library(RColorBrewer)

# library(reshape2)
# library(ComplexHeatmap)
library(pheatmap)
#library(mixOmics)
library(plyr)
library(knitr)
#library(gage)
#library(gageData)
library(gplots)
#library("AnnotationHub")
library(circlize)
library(FactoMineR)
library("factoextra")
library("corrplot")
library(readxl)
library(ggplot2)
# library(sqldf) # utilise X11
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
setwd("/Users/pgermon/Documents/Divers_labo/Publis/2024_Gitton/Soumission_PCI/scripts")

#Tm table####
table_Tm<-read_xlsx("Table_tm_Fluidigm.xlsx")

#description of samples####
description <- read_xlsx("2023-03-28-Fluidigm-Plan de plaque P96x96.xlsx", sheet="P96-Echantillons")

#cell counts table####
comptages <- read_xlsx("2023_comptages.xlsx")


#raw BiomarkHD results####
export <- read.csv("2023-03-28-Fluidigm-results.csv", skip=12, sep=";", dec=",",header=F)
colnames(export) <- c("chamber",
                      "Name",
                      "type_ech",
                      "rel_conc",
                      "reference",
                      "Gene",
                      "type_gene",
                      "ref_gene",
                      "ct",
                      "ct_quality",
                      "ct_call",
                      "Ct_threshold",
                      "delta_ct_sample",
                      "delta_ct_sample_quality",
                      "delta_ct_sample_call",
                      "delta_ct_reagent",
                      "delta_ct_reagent_quality",
                      "delta_ct_reagent_call",
                      "delta_delta_ct",
                      "delta_delta_ct_fold_change",
                      "delta_delta_ct_quality",
                      "delta_delta_ct_call",
                      "tm_in",
                      "tm_out",
                      "tm_ratio")  

}
####_####
####Clean tables####
{
  #remove non relevant samples####
  export1 <- export %>%
  filter(!str_detect(Name,"ANSES"))


  #add Tm column for genes studied####
  export2  <-  export1 %>%
    left_join(table_Tm, by=c('Gene'))

  #add description column####
  export3  <-  export2 %>%
    left_join(description, by=c('Name'))
  
  #remove samples with no Tm####
  export4 <- export3 %>%
    filter(export3$tm_std!="ND")
  export4$tm_std <- as.numeric(export4$tm_std)
  
  # first filter : determine samples Tm within the range of +/- 0.7° of correct Tm####
  export6 <- export4 %>%
    mutate(tm_max=tm_std + 0.7) %>%
    filter(!tm_in > tm_max) %>%
    mutate(tm_min=tm_std - 0.7) %>%
    filter(!tm_in<tm_min)
}  
####_####

#clean table to remove genes and samples with too many missing values####
{
#identify samples corresponding to sorted cells####
{
  quartier_tri <- description$Quartier[description$Type=="tri" & !is.na(description$Type)]
  quartier_centri <- description$Quartier[description$Type=="centri" & !is.na(description$Type)]
  quartier_lait <- description$Quartier[description$Type=="lait" & !is.na(description$Type)]
  quartier_populations_sep <- description$Quartier[(description$Type=="Macro" |
                                                      description$Type=="Lympho" |
                                                      description$Type=="Neutro" ) &
                                                     !is.na(description$Type)]
  
}

#compare sorted cells
analysis <- c("Macro-Lympho-Neutro")
quarter_to_keep <- quartier_populations_sep
samples_to_keep <- description$Name[is.element(description$Quartier,quarter_to_keep) &
                                             description$Extraction=="qiagen" &
                                             description$Type!="centri"]

#restrict table to samples to keep####
export7 <- export6[is.element(export6$Name,samples_to_keep),]

nb_echantillons <- length(levels(factor(export7$Name)))

#select genes with less than 10% missing values####
list_gene_ct999 <- export7 %>%
  count(Gene, ct==999) %>%
  rename("ct999"=c(2)) %>%
  filter(!n< nb_echantillons*0.9) %>%
  mutate(ct999=as.factor(ct999)) %>%
  filter(!ct999==TRUE)

export8 <- export7 %>%
  semi_join(list_gene_ct999, by=c('Gene')) 


#simplified table containing only sample names, gene names and Ct values####
table_simplifiee <- export8[,c("Name","Gene", "ct")]

  
table_ct <- table_simplifiee %>%
  pivot_wider(names_from = Gene, values_from = ct)

#adjust missing Ct values for a gene to max(Ct) +1 ####
for (i in (2:ncol(table_ct))) {
  # i <- 23
  max <- max(table_ct[i][!is.na(table_ct[i])])
  table_ct[i][is.na(table_ct[i])] <- max +1
}
}

####calculation of ddCt####
{
  # reference for DDct https://toptipbio.com/qpcr-multiple-reference-genes/####
  moyenne <- sapply(table_ct[,2:ncol(table_ct)],mean, na.rm=T)
  control<-c("CTRL",moyenne)
  
  
  #creer matrice avec toutes les valeurs
  matrice <- data.frame(table_ct[,2:ncol(table_ct)])
  #convertir tous en nombre
  matrice[,1:ncol(matrice)] <- apply(matrice[,1:ncol(matrice)], 2, function(x) as.numeric(as.character(x)))
  
  #matrice des ∆Ct : retirer la valeur du groupe contrôle, ici geometrique_moyenne
  matrice_delta_ct <- -sweep(matrice,2,moyenne)
  
  #matrice des RQ
  matrice_RQ <- sapply(matrice_delta_ct, function(x) 2^x)
  row.names(matrice_RQ) <- table_ct$Name
  
  
  #Calculate the geometric mean of the reference genes RQ values
  genes_REF<-c("ACTB","GAPDH")
  matrice_RQ_ref <- matrice_RQ[,genes_REF]
  
  RQ_geom_mean <- geometric.mean(t(matrice_RQ_ref), na.rm=T)
  
  matrice_RQ_geom <- matrice_RQ / RQ_geom_mean
  
  table_fc <- data.frame("Name"=rownames(matrice_RQ_geom), matrice_RQ_geom)
  
  table_log_fc <- data.frame("Name"=rownames(matrice_RQ_geom), log2(matrice_RQ_geom))

}  

  


#add description to table_log_fc####
table_log_fc <- merge(table_log_fc, description, by="Name")

#add cell counts####
table_log_fc <- merge(table_log_fc, comptages, by="Quartier")

#create duplicates column
table_log_fc$duplicats <- rep(1,nrow(table_log_fc))
table_log_fc$duplicats[grep("_2",table_log_fc$Name)] <- 2
table_log_fc$Name <- gsub("_2","",table_log_fc$Name)

#definition table####

table_analyse <- table_log_fc
  
#aggregate duplicates in table_log_FC####

table_analyse_agg <- aggregate(table_log_fc[,!is.element(colnames(table_log_fc),
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
                                                           "Total", "Serie", "jours_lait",
                                                           "GAPDH","ACTB"
                                                         ))], by=list(Name=table_log_fc$Name), FUN=mean)


table_log_fc_agg <- merge(table_analyse_agg,
                          table_log_fc[table_log_fc$duplicats==1,c("Quartier","Name", "Date",
                                                                   "Volume",
                                                                   "Echantillon","Type",
                                                                   "Extraction","Macrophage","Lymphocyte","Neutrophile", "CCS",
                                                                   "Total", "Serie",
                                                                   "Macro_num", 
                                                                   "Lympho_num", "Neutro_num",                                      
                                                                   "Race", "SCS_cate",                                         
                                                                   "SCS", "SCS_cate_2", "jours_lait",                                      
                                                                   "Cluster_3_samples", "Cluster_4_samples", "Bacterio", "Bacterio_pathoproof",                                         
                                                                   "Staphylococcus sp.", "Str. uberis",                                      
                                                                   "Str. dysgalactiae", "Beta-lactamase gene",                              
                                                                   "E. coli", "C. bovis",                                         
                                                                   "Enterococcus sp. (including faecalis and faecium)", "Klebsiella sp. (including oxytoca and pneumonia)", 
                                                                   "M. bovis", "Mycoplasma sp.",                                   
                                                                   "Prototheca sp.", "S. marcescens",                                    
                                                                   "Staph. aureus", "Str. agalactiae",                                  
                                                                   "T. pyogenes and P. indolicus", "Yeasts"                                           
                                                                   )],
                          by="Name", no.dups=T)



table_analyse <- table_log_fc_agg

# write_xlsx(table_analyse, "table_analyse_cell_triee_2024_01_30.xlsx")


#cell type specific clusters####
clusters_pop_sep <- read.csv("Cluster_pop_separees.csv", header=T)

#heatmap####
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
  
 
  #restriction to specified set of genes

  
  matrix_heatmap <- matrix_heatmap[,c(is.element(colnames(matrix_heatmap),clusters_pop_sep$Gene[is.element(clusters_pop_sep$Nom,c("Macro","Lc", "PMN2", "PMN1", "Macro_PMN1", "Macro_PMN2"))]))]
  
  
  
  rownames(matrix_heatmap) <- paste0(table_analyse$Quartier,"_",table_analyse$Type)
  
  myBreaks <- c(seq(-3, 3,by=0.1))
  paletteLength<-length(myBreaks)
  myColor <- colorRampPalette(c("blue","white","red"))(paletteLength)
 
  ann_colors = list(
    type=c(
      "Lympho"="blue",
      "Macro"="green",
      "Neutro"="red"),
    Gene_clusters = c("Lc"="blue",
              "Macro"="green",
              "PMN1"="red",
              "PMN2"="purple",
              "Macro_PMN1"="yellow",
             "Macro_PMN2"="yellow"
             )
  )


  annotation_row <- data.frame(
    CCS=log2(table_analyse$CCS /100)+3
  )
  
  rownames(annotation_row) <- paste0(table_analyse$Quartier,"_",table_analyse$Type)
  
  annotation_col <- data.frame(
    Gene_clusters = factor(clusters_pop_sep$Nom)
  )
  rownames(annotation_col) <- clusters_pop_sep$Gene

  labs.row <- rownames(matrix_heatmap)

  matrix_heatmap <- scale(matrix_heatmap)

  p <- pheatmap(matrix_heatmap,
                # cluster_rows = F,
                cutree_cols = 6,
                cutree_rows = 3,
                clustering_method = "ward.D",
                clustering_distance_cols = "correlation",
                clustering_distance_rows = "correlation",
                annotation_colors = ann_colors,
                annotation_row = annotation_row,
                annotation_col = annotation_col,
                breaks = myBreaks,
                legend_breaks = c(-3,0,3),

                fontsize_row = 14,
                fontsize_col = 10,
                fontsize = 16,
                
                color = myColor,
                labels_row = labs.row,
                # main=titre
                )
  
}

# pdf("Figure_heatmap_cell_triees.pdf", 16,12)
tiff("Figure_XX_heatmap_cell_triees.tiff",width = 1200, height = 600, units = "px")
p
dev.off()



####PCA####
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
  matrix_pca <- matrix_heatmap[,c(is.element(colnames(matrix_heatmap),clusters_pop_sep$Gene[is.element(clusters_pop_sep$Nom,c("Macro","Lc", "PMN2", "PMN1", "Macro_PMN1", "Macro_PMN2"))]))]
  
  matrix_pca <- cbind(table_analyse[,c("Total","CCS", "Neutrophile","Lymphocyte","Macrophage")], matrix_pca)
  
  rownames(matrix_pca) <- paste0(table_analyse$Quartier,"_",table_analyse$Type)
  
  
  library(FactoMineR)
  
  colonne_sup_quanti <- which(colnames(matrix_pca)%in%c("Total","CCS","Neutrophile","Lymphocyte","Macrophage"))
  
  
  res.pca <- PCA(matrix_pca,
                 quanti.sup = colonne_sup_quanti,
                 graph = F)
  
  eig.val <- get_eigenvalue(res.pca)

fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 60))
ind <- get_pca_ind(res.pca)


color <- as.factor(table_analyse$Type)


pca_ind_plot <- fviz_pca_ind(res.pca,
             axes = c(1,2),
             # pointshape = as.factor(table_analyse$Cluster_4_samples),
             pointsize = 5,
             # col.ind = "cos2",
             col.ind=color,
             # fill.ind=color,
             # gradient.cols = c("green","blue","orange","yellow"),
             palette = c("blue","green","red"),
             # palette = c("blue","cyan","orange","red"),
             repel = TRUE, # Avoid text overlapping (slow if many points)
             addEllipses = T,
             ellipse.type = "confidence",
             ellipse.level=0.99,
             mean.point = FALSE,
             label = "none",
             title = " ",
             legend.title = "Cell types"
             ) +
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

tiff("Figure_XX_PCA_cell_triee.tiff",width = 600, height = 600, units = "px")
pca_ind_plot
dev.off()
}
