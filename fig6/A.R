
setwd("D:/R/Hot_tumor_and_cold_tumor/GSVA2/Veen")

DEG_BLCA <- read.csv("D:/R/Hot_tumor_and_cold_tumor/GSVA2/BLCA/t_results_hallmark.csv",header = T)
DEG_BLCA_up <- DEG_BLCA$X[which(DEG_BLCA$t_value>= 2.58)]
DEG_BLCA_down <- DEG_BLCA$X[which(DEG_BLCA$t_value<= -2.58)]


DEG_CESC <- read.csv("D:/R/Hot_tumor_and_cold_tumor/GSVA2/CESC/t_results_hallmark.csv",header = T)
DEG_CESC_up <-  DEG_CESC$X[which(DEG_CESC$t_value>= 2.58)]
DEG_CESC_down <-  DEG_CESC$X[which(DEG_CESC$t_value<= -2.58)]


DEG_PAAD <- read.csv("D:/R/Hot_tumor_and_cold_tumor/GSVA2/PAAD/t_results_hallmark.csv",header = T)
DEG_PAAD_up <-  DEG_PAAD$X[which(DEG_PAAD$t_value>= 2.58)]
DEG_PAAD_down <-  DEG_PAAD$X[which(DEG_PAAD$t_value<= -2.58)]


DEG_SARC <- read.csv("D:/R/Hot_tumor_and_cold_tumor/GSVA2/SARC/t_results_hallmark.csv",header = T)
DEG_SARC_up <-  DEG_SARC$X[which(DEG_SARC$t_value>= 2.58)]
DEG_SARC_down <-  DEG_SARC$X[which(DEG_SARC$t_value<= -2.58)]

DEG_SKCM <- read.csv("D:/R/Hot_tumor_and_cold_tumor/GSVA2/SKCM/t_results_hallmark.csv",header = T)
DEG_SKCM_up <-  DEG_SKCM$X[which(DEG_SKCM$t_value>= 2.58)]
DEG_SKCM_down <-  DEG_SKCM$X[which(DEG_SKCM$t_value<= -2.58)]




library(VennDiagram)

venn.diagram(list(BLCA=DEG_BLCA_up,CESC=DEG_CESC_up,PAAD=DEG_PAAD_up,SARC=DEG_SARC_up,SKCM=DEG_SKCM_up), 
             fill=c("#729ECE","#FF9E4A","#67BF5C","#ED665D","#AD8BC9"),
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
             fill=c("#729ECE","#FF9E4A","#67BF5C","#ED665D","#AD8BC9"), 
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
#  "HALLMARK_ANGIOGENESIS"    "HALLMARK_APICAL_JUNCTION" "HALLMARK_APICAL_SURFACE" 

names(which(table(c(DEG_BLCA_down,DEG_CESC_down,DEG_SKCM_down,DEG_SARC_down))==4))
#   "HALLMARK_APOPTOSIS"            "HALLMARK_BILE_ACID_METABOLISM"
