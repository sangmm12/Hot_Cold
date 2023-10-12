rm(list=ls())

setwd("D:/R/Hot_tumor_and_cold_tumor/Immunecells/CIBERSORT/OVER100/cor")

cancers <- c("BLCA","CESC","KIPAN","PAAD","SARC")


#cancer <- c("SKCM")
for(cancer in cancers){

  
  
  library(Hmisc)
  
  if(dir.exists(cancer)==TRUE){
    print(cancer)
  }else {
    dir.create(cancer)
  }
  
  dat <- read.csv(paste("D:/R/Hot_tumor_and_cold_tumor/Immunecells/CIBERSORT/OVER100/",cancer,"/exp_",cancer,"_BIBERSORT.csv",sep=''),row.names = 1)
 
  res <- rcorr(as.matrix(log(t(dat)+1)))
  
  P <- res$P
  
  write.csv(P,file=paste(cancer,"/P.csv",sep=''),quote=F) 
  
  
  library(corrplot)
  dat2 <- log(t(dat)+1)
  
  M <- cor(dat2)
  
  pdf(paste(cancer,"/cor1_",cancer,".pdf",sep=''),height = 12,width = 12)
  corrplot(M,method="number")
  dev.off()
  
  
  pdf(paste(cancer,"/cor2_",cancer,".pdf",sep=''),height = 10,width = 10)
  corrplot(M,type="lower")
  dev.off()
  
  write.csv(M,file=paste(cancer,"/R.csv",sep=''),quote=F)
  
  # M[order(c(M),decreasing=TRUE)[1:200]]
  # 
  # dat_1 <- dat2[,grep("RAN",colnames(dat2))]
  # dat_2 <- dat2[,grep("GPN3",colnames(dat2))]
  # 
  # 
  # plot(dat_1,dat_2)
  # 
  # library(ggpubr)
  # library(reshape2)
  # library(ggsci)
  # d1 <- data.frame(BLM=dat_1 , FANCI=dat_2)
  # 
  # 
  # #ggscatterstats(data=d1,x=log(BLM+1),y=log(FANC+1),bf.message = F)
  # 
  # 
  # pdf("cor_max.pdf",height=5,width=5)
  # ggplot(d1, aes( x = BLM, y =FANCI)) + geom_point( size=3) + geom_smooth(method = 'lm', formula = y ~ x, se = T)+ 
  #   theme_classic() +labs(x="log(BLM+1)",y="log(FANCI+1)") + stat_cor(method = "pearson",label.x =3, label.y =1.5)
  # 
  # dev.off()
  # 
  # 
  # M[order(c(M),decreasing=F)[1:200]]
  # 
  # 
  # 
  # 
  # dat_1 <- dat2[,grep("AHCY",colnames(dat2))]
  # dat_2 <- dat2[,grep("AHNAK",colnames(dat2))]
  # 
  # d1 <- data.frame(BLM=dat_1 , FANCI=dat_2)
  # 
  # pdf("cor_min.pdf",height=5,width=5)
  # ggplot(d1, aes( x = BLM, y =FANCI)) + geom_point( size=3) + geom_smooth(method = 'lm', formula = y ~ x, se = T)+ 
  #   theme_classic()+labs(x="log(CUL3+1)",y="log(POLR2L+1)") + stat_cor(method = "pearson",label.x =2.2, label.y =3)
  # 
  # dev.off()
  # 
}