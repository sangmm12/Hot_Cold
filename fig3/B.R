rm(list=ls())

setwd("D:/R/Hot_tumor_and_cold_tumor/cor/cell_score/zong")

cancers <- c("BLCA","CESC","PAAD","SARC","SKCM")
#
#cancer <- c("BLCA")
for(cancer in cancers){
  
  
  library(data.table)
  
  dat <-  fread(paste("D:/R/Hot_tumor_and_cold_tumor/cluster_cell/cluster_cellzong/",cancer,"/data_",cancer,".csv",sep=''),header = T)
  
  dat <- as.data.frame(dat)
  
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]
  
  
  celllist <- c("T_cells_CD8","NK_cells_activated","T_cells_follicular_helper","Tcells_proliferation","Cytolytic_activity","Immunogenic_cell_death","Macrophages_M2","MDSCs","T_cells_regulatory_.Tregs.","T_cell_co-inhibition","MHC_class_I","Inflammation-promoting")
  
  data <- dat[,match(celllist,colnames(dat))]
  
  colnames(data) <- c("T_cells_CD8","NK_cells_activated","T_cells_follicular_helper","Tcells_proliferation","Cytolytic_activity","Immunogenic_cell_death","Macrophages_M2","MDSCs","Tregs","T_cell_co_inhibition","MHC_class_I","Inflammation-promoting")
  
  cells <- c("T_cells_CD8","NK_cells_activated","T_cells_follicular_helper","Macrophages_M2","Tregs")
  scores <- c("Tcells_proliferation","Cytolytic_activity","Immunogenic_cell_death","MDSCs","T_cell_co_inhibition","MHC_class_I","Inflammation-promoting")
  
   dat_gene <- data[,match( cells,colnames(data))]
   dat_im <- data[,match(scores,colnames(data))]
   
   library(psych)
   data.corr <- corr.test(dat_im, dat_gene, method="pearson", adjust="fdr")
   data.r <- data.corr$r  # 相关系数
   data.p <- data.corr$p  # p值
   
   paste("data.r_",cancer,".csv",sep='')
   write.csv(data.r,file= paste(cancer,"_data.r.csv",sep=''),quote=F)
   write.csv(data.p,file= paste(cancer,"_data.p.csv",sep=''),quote=F)
   
   
   library(pheatmap)
   getSig <- function(dc) {
     print(dc)
     sc <- ' '
     if (dc < 0.0001) {sc <- '****'}
     else if (dc < 0.001){sc <- '***'}
     else if (dc < 0.01){sc <- '**'}
     else if (dc < 0.05) {sc <- '*'}
     else{sc <- ''
     }
     return(sc)
   }
   sig.mat <- matrix(sapply(data.p, getSig), nrow=nrow(data.p))
   str(sig.mat)
   
   paletteLength <- 1000
   myColor <- colorRampPalette(c("#36648b", "white", "#e94644"))(paletteLength)
   
   test <- data.r
   myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1), 
                 seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
   #myBreaks <- c(seq(-0.2, 0, length.out=ceiling(paletteLength/2) + 1), 
   #seq(0.4/paletteLength, 0.4, length.out=floor(paletteLength/2)))
   
   
   #chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...
   
   pdf(paste(cancer,"_cor_gene_score.pdf",sep=''),width =4,height =4)
   pheatmap(data.r,
            color=myColor,
            breaks=myBreaks,
            clustering_method="average", cluster_rows=F, cluster_cols=F,display_numbers=sig.mat
   )
   
   dev.off()
   
  
 
  
}
