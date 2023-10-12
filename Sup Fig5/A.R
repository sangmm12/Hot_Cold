setwd("D:/R/Hot_tumor_and_cold_tumor/immune data1/cor_cell$gene/zong")


cancers <- c("BLCA","CESC","PAAD","SARC","SKCM")

out_put_list <- c()

library(dplyr)
library(patchwork)
library(ggplotify)

#cancer <- c("BLCA")
for(cancer in cancers){
  
  
  library(data.table)
  dat <- fread(paste("D:/R/Hot_tumor_and_cold_tumor/Immunecells/8cancers/",cancer,"_BIBERSORT.csv",sep=''))
  dat_immu <- as.data.frame(dat)
  rownames(dat_immu) <- dat_immu[,1]
  dat_immu <- dat_immu[,-1]
  
  
  
  dat_DDR <-  fread(paste("D:/R/Hot_tumor_and_cold_tumor/immune data1/immu_",cancer,"_zong.csv",sep=''),header = T)
  dat_DDR <- as.data.frame(dat_DDR)
  
  rownames(dat_DDR) <- dat_DDR[,1]
  dat_DDR <- dat_DDR[,-1]
  
  # colnames(dat_DDR) <- dat_DDR[1,]
  # dat_DDR <- dat_DDR[-1,]
  
  dat_DDR <- t(dat_DDR)
  #dat_DDR <- as.data.frame(dat_DDR)
  
  
  all_name <- names(which(table(c(rownames(dat_immu),rownames(dat_DDR)))==2))
  
  
  dat_gene <- dat_DDR[match(all_name,rownames(dat_DDR)),]
  
  
  #dat_gene <- dat_gene[,match(gene_DDR,colnames(dat_gene))]
  
  
  dat_imm <- dat_immu[match(all_name,rownames(dat_immu)),]
  
  
  #cell <- c("T_cells_CD8","T_cells_CD4_memory_resting","T_cells_CD4_memory_activated","T_cells_follicular_helper","T_cells_regulatory_.Tregs.","T_cells_gamma_delta")
  dat_im <- dat_imm
  
  
  colSums(dat_im)
  
  
  
  library(psych)
  data.corr <- corr.test(dat_gene, dat_im, method="pearson", adjust="fdr")
  data.r <- data.corr$r  # 相关系数
  data.p <- data.corr$p  # p值
  
  paste("data.r_",cancer,".csv",sep='')
  write.csv(data.r,file= paste("data.r_",cancer,".csv",sep=''),quote=F)
  write.csv(data.p,file= paste("data.p_",cancer,".csv",sep=''),quote=F)
  
  
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
  # myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1), 
  #               seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
  myBreaks <- c(seq(-0.6,  -0.00001, length.out=ceiling(paletteLength/2) + 1), 
                seq(0, 0.6, length.out=floor(paletteLength/2)))
  
  #chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...
  
  eval(parse(text = paste(cancer,'<- as.ggplot(pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F,cluster_cols=F,cellwidth =200,cellheight =200, display_numbers=sig.mat,fontsize=length(colnames(dat_gene))))
')))
  out_put_list <- append(out_put_list,cancer)
  
  class1 <- rep("Inhibitory", times=24)
  class2 <- rep("Stimulaotry", times=36)
  class3 <- rep("Chemokine", times=37)
  class4 <- rep("Immunoinhibitor", times=10)
  class5 <- rep("Immunostimulator", times=27)
  class6 <- rep("MHC", times=21)
  class7 <- rep("Receptor", times=17)
  
  
  
  class<- c(class1,class2,class3,class4,class5,class6,class7)
  
  annotation_row <- data.frame(class)
  
  rownames(annotation_row) <- colnames(dat_gene)
  
  
  pdf(paste("cor_",cancer,"_immu.pdf",sep=''),width =7,height = 27)
  pheatmap(data.r, 
           color=myColor,
           breaks=myBreaks,
           clustering_method="average",cluster_rows = F,annotation_row = annotation_row,cluster_cols=F,display_numbers=sig.mat)

  
  
  dev.off()
  
  
}

pdf(paste("cor_zong_immu.pdf",sep=''),width =460,height = 600)

eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=5,labels=out_put_list,label_size=0.1))',sep='')))

dev.off()