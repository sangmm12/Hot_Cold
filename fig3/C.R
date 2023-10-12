rm(list=ls())

setwd("D:/R/Hot_tumor_and_cold_tumor/cor/cell/cor_single")

cancers <- c("BLCA","CESC","PAAD","SARC","SKCM")
#
#cancer <- c("BLCA")
for(cancer in cancers){
  
  library(data.table)
  dat <- fread(paste("D:/R/Hot_tumor_and_cold_tumor/Immunecells/8cancers/",cancer,"_BIBERSORT.csv",sep=''))
  data <- as.data.frame(dat)
  aa <- colnames(data)
  aa <- aa[-1]
  
  for(i in 1:length(aa)){
    for(j in 1:length(aa)){
      
      #i=1
      #j=1
      cell <- aa[i]
      score <- aa[j]
      
      dat_1 <- data[,grep(cell,colnames(data))]
      dat_2 <- data[,grep(score,colnames(data))]
      
      
      plot(dat_1,dat_2)
      
      library(ggpubr)
      library(reshape2)
      library(ggsci)
      d1 <- data.frame(NT5E=dat_1 ,VEGFA=dat_2)
      
      library(Hmisc)
      
      #ggscatterstats(data=d1,x=log(NT5E+1),y=log(VEGFA+1),bf.message = F)
      
      
      pdf( paste("cor_",cancer,"_",cell,"_",score,".pdf",sep=''),height=5,width=5)
      # p <- ggplot(d1, aes( x = NT5E, y =VEGFA)) + geom_point( size=3) + geom_smooth(method = 'lm', formula = y ~ x, se = T)+
      #   theme(axis.text = element_text(size = 25),axis.title = element_text(size = 25))+
      #   theme_classic() +labs(x="log(NT5E+1)",y="log(CD276+1)") + stat_cor(method = "pearson",size=8,,label.x =min(dat_1), label.y =max(dat_2))
      
      p <- ggplot(d1,aes(x=NT5E,y=VEGFA)) + geom_point(shape=19) +
        xlab( paste(cell,sep='')) + ylab(paste(score,sep=''))+geom_smooth(method = lm)+theme_classic()+
        theme(axis.text = element_text(size = 15),axis.title = element_text(size = 25))+
        stat_cor(method = "pearson",size=8,label.x =min(dat_1), label.y =max(dat_2))
      
      
      print(p)
      dev.off()
      
      
    }
  }
}







############################################
cancer <- c("BLCA")

library(data.table)
dat <- fread(paste("D:/R/Hot_tumor_and_cold_tumor/Immunecells/8cancers/",cancer,"_BIBERSORT.csv",sep=''))
data <- as.data.frame(dat)
 

cell <- c("T_cells_CD8")
score <- c("T_cells_CD4_memory_activated")

dat_1 <- data[,grep(cell,colnames(data))]
dat_2 <- data[,grep(score,colnames(data))]


plot(dat_1,dat_2)

library(ggpubr)
library(reshape2)
library(ggsci)
d1 <- data.frame(NT5E=dat_1 ,VEGFA=dat_2)

library(Hmisc)

#ggscatterstats(data=d1,x=log(NT5E+1),y=log(VEGFA+1),bf.message = F)


pdf( paste("cor_",cancer,"_",cell,"_",score,".pdf",sep=''),height=5,width=5)
# p <- ggplot(d1, aes( x = NT5E, y =VEGFA)) + geom_point( size=3) + geom_smooth(method = 'lm', formula = y ~ x, se = T)+
#   theme(axis.text = element_text(size = 25),axis.title = element_text(size = 25))+
#   theme_classic() +labs(x="log(NT5E+1)",y="log(CD276+1)") + stat_cor(method = "pearson",size=8,,label.x =min(dat_1), label.y =max(dat_2))

p <- ggplot(d1,aes(x=NT5E,y=VEGFA)) + geom_point(shape=19) +
  xlab( paste(cell,sep='')) + ylab(paste(score,sep=''))+geom_smooth(method = lm)+theme_classic()+
  theme(axis.text = element_text(size = 15),axis.title = element_text(size = 25))+
  stat_cor(method = "pearson",size=8,label.x =min(dat_1), label.y =max(dat_2))


print(p)
dev.off()

     
        






############################################
cancer <- c("CESC")

library(data.table)
dat <- fread(paste("D:/R/Hot_tumor_and_cold_tumor/Immunecells/8cancers/",cancer,"_BIBERSORT.csv",sep=''))
data <- as.data.frame(dat)


cell <- c("Plasma_cells")
score <- c("B_cells_naive")

dat_1 <- data[,grep(cell,colnames(data))]
dat_2 <- data[,grep(score,colnames(data))]


plot(dat_1,dat_2)

library(ggpubr)
library(reshape2)
library(ggsci)
d1 <- data.frame(NT5E=dat_1 ,VEGFA=dat_2)

library(Hmisc)

#ggscatterstats(data=d1,x=log(NT5E+1),y=log(VEGFA+1),bf.message = F)


pdf( paste("cor_",cancer,"_",cell,"_",score,".pdf",sep=''),height=5,width=5)
# p <- ggplot(d1, aes( x = NT5E, y =VEGFA)) + geom_point( size=3) + geom_smooth(method = 'lm', formula = y ~ x, se = T)+
#   theme(axis.text = element_text(size = 25),axis.title = element_text(size = 25))+
#   theme_classic() +labs(x="log(NT5E+1)",y="log(CD276+1)") + stat_cor(method = "pearson",size=8,,label.x =min(dat_1), label.y =max(dat_2))

p <- ggplot(d1,aes(x=NT5E,y=VEGFA)) + geom_point(shape=19) +
  xlab( paste(cell,sep='')) + ylab(paste(score,sep=''))+geom_smooth(method = lm)+theme_classic()+
  theme(axis.text = element_text(size = 15),axis.title = element_text(size = 25))+
  stat_cor(method = "pearson",size=8,label.x =min(dat_1), label.y =max(dat_2))


print(p)
dev.off()








############################################
cancer <- c("PAAD")

library(data.table)
dat <- fread(paste("D:/R/Hot_tumor_and_cold_tumor/Immunecells/8cancers/",cancer,"_BIBERSORT.csv",sep=''))
data <- as.data.frame(dat)


cell <- c("T_cells_CD8")
score <- c("Macrophages_M0")

dat_1 <- data[,grep(cell,colnames(data))]
dat_2 <- data[,grep(score,colnames(data))]


plot(dat_1,dat_2)

library(ggpubr)
library(reshape2)
library(ggsci)
d1 <- data.frame(NT5E=dat_1 ,VEGFA=dat_2)

library(Hmisc)

#ggscatterstats(data=d1,x=log(NT5E+1),y=log(VEGFA+1),bf.message = F)


pdf( paste("cor_",cancer,"_",cell,"_",score,".pdf",sep=''),height=5,width=5)
# p <- ggplot(d1, aes( x = NT5E, y =VEGFA)) + geom_point( size=3) + geom_smooth(method = 'lm', formula = y ~ x, se = T)+
#   theme(axis.text = element_text(size = 25),axis.title = element_text(size = 25))+
#   theme_classic() +labs(x="log(NT5E+1)",y="log(CD276+1)") + stat_cor(method = "pearson",size=8,,label.x =min(dat_1), label.y =max(dat_2))

p <- ggplot(d1,aes(x=NT5E,y=VEGFA)) + geom_point(shape=19) +
  xlab( paste(cell,sep='')) + ylab(paste(score,sep=''))+geom_smooth(method = lm)+theme_classic()+
  theme(axis.text = element_text(size = 15),axis.title = element_text(size = 25))+
  stat_cor(method = "pearson",size=8,label.x =min(dat_1), label.y =max(dat_2))


print(p)
dev.off()







############################################
cancer <- c("SARC")

library(data.table)
dat <- fread(paste("D:/R/Hot_tumor_and_cold_tumor/Immunecells/8cancers/",cancer,"_BIBERSORT.csv",sep=''))
data <- as.data.frame(dat)


cell <- c("T_cells_CD8")
score <- c("T_cells_follicular_helper")

dat_1 <- data[,grep(cell,colnames(data))]
dat_2 <- data[,grep(score,colnames(data))]


plot(dat_1,dat_2)

library(ggpubr)
library(reshape2)
library(ggsci)
d1 <- data.frame(NT5E=dat_1 ,VEGFA=dat_2)

library(Hmisc)

#ggscatterstats(data=d1,x=log(NT5E+1),y=log(VEGFA+1),bf.message = F)


pdf( paste("cor_",cancer,"_",cell,"_",score,".pdf",sep=''),height=5,width=5)
# p <- ggplot(d1, aes( x = NT5E, y =VEGFA)) + geom_point( size=3) + geom_smooth(method = 'lm', formula = y ~ x, se = T)+
#   theme(axis.text = element_text(size = 25),axis.title = element_text(size = 25))+
#   theme_classic() +labs(x="log(NT5E+1)",y="log(CD276+1)") + stat_cor(method = "pearson",size=8,,label.x =min(dat_1), label.y =max(dat_2))

p <- ggplot(d1,aes(x=NT5E,y=VEGFA)) + geom_point(shape=19) +
  xlab( paste(cell,sep='')) + ylab(paste(score,sep=''))+geom_smooth(method = lm)+theme_classic()+
  theme(axis.text = element_text(size = 15),axis.title = element_text(size = 25))+
  stat_cor(method = "pearson",size=8,label.x =min(dat_1), label.y =max(dat_2))


print(p)
dev.off()







############################################
cancer <- c("SKCM")

library(data.table)
dat <- fread(paste("D:/R/Hot_tumor_and_cold_tumor/Immunecells/8cancers/",cancer,"_BIBERSORT.csv",sep=''))
data <- as.data.frame(dat)


cell <- c("T_cells_CD8")
score <- c("NK_cells_activated")

dat_1 <- data[,grep(cell,colnames(data))]
dat_2 <- data[,grep(score,colnames(data))]


plot(dat_1,dat_2)

library(ggpubr)
library(reshape2)
library(ggsci)
d1 <- data.frame(NT5E=dat_1 ,VEGFA=dat_2)

library(Hmisc)

#ggscatterstats(data=d1,x=log(NT5E+1),y=log(VEGFA+1),bf.message = F)


pdf( paste("cor_",cancer,"_",cell,"_",score,".pdf",sep=''),height=5,width=5)
# p <- ggplot(d1, aes( x = NT5E, y =VEGFA)) + geom_point( size=3) + geom_smooth(method = 'lm', formula = y ~ x, se = T)+
#   theme(axis.text = element_text(size = 25),axis.title = element_text(size = 25))+
#   theme_classic() +labs(x="log(NT5E+1)",y="log(CD276+1)") + stat_cor(method = "pearson",size=8,,label.x =min(dat_1), label.y =max(dat_2))

p <- ggplot(d1,aes(x=NT5E,y=VEGFA)) + geom_point(shape=19) +
  xlab( paste(cell,sep='')) + ylab(paste(score,sep=''))+geom_smooth(method = lm)+theme_classic()+
  theme(axis.text = element_text(size = 15),axis.title = element_text(size = 25))+
  stat_cor(method = "pearson",size=8,label.x =min(dat_1), label.y =max(dat_2))


print(p)
dev.off()


