setwd("D:/R/Hot_tumor_and_cold_tumor/immune data1/DGEs/cor/KEGG2")


cancers <- c("BLCA","CESC","PAAD","SARC","SKCM")

data <-  fread(paste("D:/R/Hot_tumor_and_cold_tumor/immune data1/DGEs/cor/KEGG/data.p_BLCA.csv",sep=''))
pathyways <- colnames(data)
pathyways <- pathyways[-1]

for(i in 1:length(pathyways)){
  
  pathyway <- pathyways[i]
  
  #cancer <- c("BLCA")
  for(cancer in cancers){
    
    library(data.table)
    
    datap <-  fread(paste("D:/R/Hot_tumor_and_cold_tumor/immune data1/DGEs/cor/KEGG/data.p_",cancer,".csv",sep=''))
    datap <- as.data.frame(datap)
    
    rownames(datap) <- datap[,1]
    datap <- datap[,-1]
    
    datar <-  fread(paste("D:/R/Hot_tumor_and_cold_tumor/immune data1/DGEs/cor/KEGG/data.r_",cancer,".csv",sep=''))
    datar <- as.data.frame(datar)
    
    rownames(datar) <- datar[,1]
    datar <- datar[,-1]
    
    if(cancer=="BLCA"){
      pathyway_p <- data.frame(rownames(datap),datap[,i])
      pathyway_r <- data.frame(rownames(datar),datar[,i])
      
    }else {
      pathyway_p <- data.frame(pathyway_p,datap[,i])
      pathyway_r <- data.frame(pathyway_r,datar[,i])
      
    }
  }
  
  rownames(pathyway_p) <- pathyway_p[,1]
  pathyway_p<- pathyway_p[,-1]
  colnames(pathyway_p) <- cancers
  data.p <- pathyway_p
  #data.p <- as.numeric(unlist(data.p))
  data.p <- data.frame(data.p)
  data.p <- as.matrix(data.p)
  
  rownames(pathyway_r) <- pathyway_r[,1]
  pathyway_r<- pathyway_r[,-1]
  colnames(pathyway_r) <- cancers
  data.r <- pathyway_r
  
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
  myBreaks <- c(seq(-0.3, 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(0.4/paletteLength, 0.7, length.out=floor(paletteLength/2)))
  
  class1 <- rep("up", times=22)
  class2 <- rep("down", times=1)
  
  class<- c(class1,class2)
  
  annotation_row <- data.frame(class)
  
  rownames(annotation_row) <- rownames(data.r)
  
  
  #chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...
  
  pdf(paste("cor_",pathyway,"_DGEs_GSVA.pdf",sep=''),width =3.5,height = 7)
  pheatmap(data.r, 
           color=myColor,
           breaks=myBreaks,
           clustering_method="average", cluster_rows=F, cluster_cols=F,display_numbers=sig.mat,
           annotation_row = annotation_row
  )
  
  dev.off()
  
  
  
  
}

