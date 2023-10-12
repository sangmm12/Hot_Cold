
setwd("D:/R/Hot_tumor_and_cold_tumor/GSVA/KEGG/Veen")

DEG_BLCA <- read.csv("D:/R/Hot_tumor_and_cold_tumor/GSVA/KEGG/BLCA/t_results.csv",header = T)
DEG_BLCA_up <- DEG_BLCA$X[which(DEG_BLCA$t_value>= 2.58)]
DEG_BLCA_down <- DEG_BLCA$X[which(DEG_BLCA$t_value<= -2.58)]


DEG_CESC <- read.csv("D:/R/Hot_tumor_and_cold_tumor/GSVA/KEGG/CESC/t_results.csv",header = T)
DEG_CESC_up <-  DEG_CESC$X[which(DEG_CESC$t_value>= 2.58)]
DEG_CESC_down <-  DEG_CESC$X[which(DEG_CESC$t_value<= -2.58)]


DEG_PAAD <- read.csv("D:/R/Hot_tumor_and_cold_tumor/GSVA/KEGG/PAAD/t_results.csv",header = T)
DEG_PAAD_up <-  DEG_PAAD$X[which(DEG_PAAD$t_value>= 2.58)]
DEG_PAAD_down <-  DEG_PAAD$X[which(DEG_PAAD$t_value<= -2.58)]


DEG_SARC <- read.csv("D:/R/Hot_tumor_and_cold_tumor/GSVA/KEGG/SARC/t_results.csv",header = T)
DEG_SARC_up <-  DEG_SARC$X[which(DEG_SARC$t_value>= 2.58)]
DEG_SARC_down <-  DEG_SARC$X[which(DEG_SARC$t_value<= -2.58)]

DEG_SKCM <- read.csv("D:/R/Hot_tumor_and_cold_tumor/GSVA/KEGG/SKCM/t_results.csv",header = T)
DEG_SKCM_up <-  DEG_SKCM$X[which(DEG_SKCM$t_value>= 2.58)]
DEG_SKCM_down <-  DEG_SKCM$X[which(DEG_SKCM$t_value<= -2.58)]



library(VennDiagram)

venn.diagram(list(BLCA=DEG_BLCA_up,CESC=DEG_CESC_up,PAAD=DEG_PAAD_up,SARC=DEG_SARC_up,SKCM=DEG_SKCM_up), 
             fill=c("#729ECE","#FF9E4A","#67BF5C","#ED665D","#AD8BC9"), #,"#67BF5C","#ED665D","#AD8BC9"
             alpha=c(0.5,0.5,0.5,0.5,0.5), 
             col = c("#729ECE","#FF9E4A","#67BF5C","#ED665D","#AD8BC9"),
             cex=1,
             cat.pos = c(-11,10,10,10,11),
             cat.dist = c(0.05,0.05,-0.05,-0.06,0.05),
             cat.fontface=4, 
             # fontfamily=1,
             filename ="GSVA_up.tiff",
             height = 1600, 
             width = 1600, resolution = 500)



library(VennDiagram)

venn.diagram(list(BLCA=DEG_BLCA_down,CESC=DEG_CESC_down,PAAD=DEG_PAAD_down,SARC=DEG_SARC_down,SKCM=DEG_SKCM_down), 
             fill=c("#729ECE","#FF9E4A","#67BF5C","#ED665D","#AD8BC9"), #,"#67BF5C","#ED665D","#AD8BC9"
             alpha=c(0.5,0.5,0.5,0.5,0.5), 
             col = c("#729ECE","#FF9E4A","#67BF5C","#ED665D","#AD8BC9"),
             cex=1,
             cat.pos = c(-11,10,10,10,11),
            cat.dist = c(0.05,0.05,-0.05,-0.06,0.05),
             cat.fontface=4, 
             # fontfamily=1,
             filename ="GSVA_down.tiff",
             height = 1600, 
             width = 1600, resolution = 500)



names(which(table(c(DEG_BLCA_down,DEG_CESC_down,DEG_PAAD_down,DEG_SARC_down))==4))
# KEGG_ADHERENS_JUNCTION
# EGG_ALZHEIMERS_DISEASE
# KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM

names(which(table(c(DEG_BLCA_down,DEG_CESC_down,DEG_SKCM_down,DEG_SARC_down))==4))
# KEGG_AUTOIMMUNE_THYROID_DISEASE
# KEGG_BIOSYNTHESIS_OF_UNSATURATED_FATTY_ACIDS
# KEGG_BLADDER_CANCER

names(which(table(c(DEG_PAAD_up,DEG_CESC_up,DEG_SKCM_up,DEG_SARC_up))==4))
