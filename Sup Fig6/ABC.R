setwd("D:/R/Hot_tumor_and_cold_tumor/immune data1/gene/Heterogeneity")

cancer <- c("PAAD")

hets <- c("LOH","MATH","MSI","NEO","ploidy","purity","TMB","HRD")

het <- c("NEO")

data.r <- c()
data.p <- c()

for(het in hets){
  
  library(data.table)
  dat <- fread("D:/R/Hot_tumor_and_cold_tumor/immune data1/gene/exp_NT5E&VEGFA.csv",header = T)
  #dat <- fread("D:/R/Hot_tumor_and_cold_tumor/immune data1/gene/KM/sur_PAAD_NT5E&VEGFA_hypoxia.csv",header = T)
  dat <- as.data.frame(dat)
  
  rownames(dat) <- dat[,1]
  dat_gene <- dat[,-1]

  #dat_gene <-  dat_gene[,3:4]
  
  if(het=="HRD"){
    datLOH <- read.csv(paste("D:/R/RNA_mood_PAAD/Heterogeneity/",het,".csv",sep=''))
    n1 <- grep(cancer,datLOH$X)
    dataLOH <- datLOH[n1,]
    dataLOH <- dataLOH[,-1]
    dataLOH <- dataLOH[,-2]
    dataLOH <- data.frame(dataLOH)
    
  }else {
    datLOH <- read.csv(paste("D:/R/RNA_mood_PAAD/Heterogeneity/",het,".csv",sep=''))
    n1 <- grep(cancer,datLOH$X)
    dataLOH <- datLOH[n1,]
    dataLOH <- dataLOH[,-1]
    dataLOH <- dataLOH[,-2:-36]
    dataLOH <- data.frame(dataLOH)
    
  }
  
  
  
  all_name <- names(which(table(c(rownames(dat_gene),dataLOH[,2] ))==2))
  
  data_gene <- dat_gene[match(all_name,rownames(dat_gene)),]
  
  dataLOH1<- dataLOH[match(all_name,dataLOH[,2]),]
  
  data<- data.frame(data_gene,Heterogeneity=dataLOH1[,1])
  
  
  
  dat_1 <- data[,grep("NT5E",colnames(data))]
  dat_2 <- data[,grep("Heterogeneity",colnames(data))]
  
  dat.cor <- cor.test(log(dat_1+1), log(dat_2+1), method="pearson", adjust="fdr")
  dat.r <- dat.cor$estimate  # 相关系数
  dat.p <- dat.cor$p.value  # p值
  
  r <- dat.r
  p <- dat.p
 
  
  plot(dat_1,dat_2)
  
  library(ggpubr)
  library(reshape2)
  library(ggsci)
  d1 <- data.frame(NT5E=dat_1 ,VEGFA=dat_2)
  NT5E=dat_1
  VEGFA=dat_2
  
  library(Hmisc)
  
  #ggscatterstats(data=d1,x=log(NT5E+1),y=log(VEGFA+1),bf.message = F)
  
  
  pdf(paste("cor_NT5E_",het,".pdf",sep=''),height=5,width=5)
  p1 <- ggplot(d1, aes( x = log(NT5E+1), y =log(VEGFA+1))) + geom_point( size=3) + geom_smooth(method = 'lm', formula = y ~ x, se = T)+ 
    theme_classic() +labs(x="log(NT5E+1)",y=paste("log(",het,"+1)",sep='')) + stat_cor(method = "pearson",label.x =2, label.y =1)
  print(p1)
  dev.off()
  
 p1$coordinates
  
  
  dat_1 <- data[,grep("CD276",colnames(data))]
  dat_2 <- data[,grep("Heterogeneity",colnames(data))]
  
  dat.cor <- cor.test(log(dat_1+1), log(dat_2+1), method="pearson", adjust="fdr")
  dat.r <- dat.cor$estimate  # 相关系数
  dat.p <- dat.cor$p.value  # p值
  
  r <- c(r,dat.r)
  p <- c(p,dat.p)
  
  plot(dat_1,dat_2)
  
  library(ggpubr)
  library(reshape2)
  library(ggsci)
  d1 <- data.frame(NT5E=dat_1 ,VEGFA=dat_2)
  NT5E=dat_1
  VEGFA=dat_2
  
  library(Hmisc)
  
  #ggscatterstats(data=d1,x=log(NT5E+1),y=log(VEGFA+1),bf.message = F)
  
  
  pdf(paste("cor_CD276_",het,".pdf",sep=''),height=5,width=5)
  p2 <- ggplot(d1, aes( x = log(NT5E+1), y =log(VEGFA+1))) + geom_point( size=3) + geom_smooth(method = 'lm', formula = y ~ x, se = T)+ 
    theme_classic() +labs(x="log(CD276+1)",y=paste("log(",het,"+1)",sep='')) + stat_cor(method = "pearson",label.x =min(log(NT5E+1)), label.y =max(log(VEGFA+1)))
  print(p2)
  dev.off()
  
  
  if(length(data.r)==0){
    data.r <- r
    data.p <- p
    
  }else {
    data.r <- cbind(data.r,r)
    data.p <- cbind(data.p,p)
    
  }
  

}

colnames(data.r) <- hets
rownames(data.r) <- c("NT5E","CD276")

colnames(data.p) <- hets
rownames(data.p) <- c("NT5E","CD276")



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
myBreaks <- c(seq(-max(abs(max(data.r)),abs(min(data.r))),  -0.00001, length.out=ceiling(paletteLength/2) + 1), 
              seq(0, max(abs(max(data.r)),abs(min(data.r))), length.out=floor(paletteLength/2)))

#chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...

pdf(paste("cor_",cancer,"_NT5E&CD276_Heterogeneity.pdf",sep=''),width =7,height = 2.3)
pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average",cluster_rows = F,cluster_cols=F,display_numbers=sig.mat)



dev.off()

