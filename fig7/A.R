setwd("D:/R/Hot_tumor_and_cold_tumor/immune data1/cor_cell$gene/zong/cor12")


cancers <- c("BLCA","CESC","PAAD","SARC","SKCM")
dat_p <- c()

#cancer <- c("BLCA")
for(cancer in cancers){
  
  library(data.table)
  
  dat <- fread(paste("D:/R/Hot_tumor_and_cold_tumor/immune data1/cor_cell$gene/zong/cor4zong/data.p_",cancer,".csv",sep=''))
  dat <- as.data.frame(dat)
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]
  
  if(length(dat_p)==0){
    dat_p <- dat
    
  }else {
    dat_p <- cbind(dat_p,dat)
  }
}

number <- rep(0, times=length(rownames(dat_p)))
dat_p <- t(dat_p)

for (i in 1:length(colnames(dat_p))) {
  aaa <- subset(dat_p, dat_p[,i] < 0.05)
  number[i] <- length(rownames(aaa))
  
}

dat_p <- t(dat_p)
dat_p <- cbind(dat_p,number)



gene <- subset(dat_p, dat_p[,21] >= 12)

genelist <- rownames(gene)


cancers <- c("BLCA","CESC","PAAD","SARC","SKCM")

out_put_list <- c()

library(dplyr)
library(patchwork)
library(ggplotify)


celllist <- c("NK_cells_activated","Macrophages_M2","T_cells_follicular_helper","T_cells_CD8")

#cancer <- c("PAAD")
for(cancer in cancers){
  
  
  library(data.table)
  dat <- fread(paste("D:/R/Hot_tumor_and_cold_tumor/Immunecells/8cancers/",cancer,"_BIBERSORT.csv",sep=''))
  dat_immu <- as.data.frame(dat)
  rownames(dat_immu) <- dat_immu[,1]
  dat_immu <- dat_immu[,-1]
  
  
  
  dat_DDR <- fread(paste("D:/R/Hot_tumor_and_cold_tumor/immune data1/immu_",cancer,"_zong.csv",sep=''),header = T)
  dat_DDR <- as.data.frame(dat_DDR)
  
  rownames(dat_DDR) <- dat_DDR[,1]
  dat_DDR <- dat_DDR[,-1]
  
  length(rownames(dat_DDR))
  
  dat_DDR <- dat_DDR[match(genelist,rownames(dat_DDR)),]
  
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
  
  data.r <- data.r[,match(celllist,colnames(data.r))]
  data.p <- data.p[,match(celllist,colnames(data.p))]
  
  
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
  
  #chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ..
  
  class1 <- rep("Inhibitory", times=7)
  class2 <- rep("Stimulaotry", times=11)
  class3 <- rep("Chemokine", times=3)
  class4 <- rep("Immunoinhibitor", times=2)
  class5 <- rep("Immunostimulator", times=3)
  class6 <- rep("MHC", times=9)
  class7 <- rep("Receptor", times=5)
  
  class<- c(class1,class2,class3,class4,class5,class6,class7)
  
  annotation_row <- data.frame(class)
  
  rownames(annotation_row) <- colnames(dat_gene)
  
  
  pdf(paste("cor_",cancer,"_immu10.pdf",sep=''),width =4,height = 10)
  pheatmap(data.r, 
           color=myColor,
           breaks=myBreaks,
           clustering_method="average",cluster_rows = F,cluster_cols=F,annotation_row = annotation_row,display_numbers=sig.mat)
  
  
  
  dev.off()
  
  
}



