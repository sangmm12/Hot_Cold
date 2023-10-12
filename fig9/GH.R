setwd("D:/R/Hot_tumor_and_cold_tumor/immune data1/gene/cor")

library(data.table)
dat <- fread("D:/R/Hot_tumor_and_cold_tumor/immune data1/gene/KM/sur_PAAD_NT5E&VEGFA_hypoxia.csv",header = T)
dat <- as.data.frame(dat)

rownames(dat) <- dat[,1]
data <- dat[,-1]



#############################################

dat_1 <- data[,grep("NT5E",colnames(data))]
dat_2 <- data[,grep("CD276",colnames(data))]


plot(dat_1,dat_2)

library(ggpubr)
library(reshape2)
library(ggsci)
d1 <- data.frame(NT5E=dat_1 ,VEGFA=dat_2)

library(Hmisc)

#ggscatterstats(data=d1,x=log(NT5E+1),y=log(VEGFA+1),bf.message = F)


pdf("cor_NT5E_CD276.pdf",height=5,width=5)
ggplot(d1, aes( x = log(NT5E+1), y =log(VEGFA+1))) + geom_point( size=3) + geom_smooth(method = 'lm', formula = y ~ x, se = T)+ 
  theme_classic() +labs(x="log(NT5E+1)",y="log(CD276+1)") + stat_cor(method = "pearson",label.x =2, label.y =1)

dev.off()




#############################################

dat_1 <- data[,grep("CD276",colnames(data))]
dat_2 <- data[,grep("Hypoxia",colnames(data))]


plot(dat_1,dat_2)

library(ggpubr)
library(reshape2)
library(ggsci)
d1 <- data.frame(NT5E=dat_1 ,VEGFA=dat_2)

library(Hmisc)

#ggscatterstats(data=d1,x=log(NT5E+1),y=log(VEGFA+1),bf.message = F)


pdf("cor_CD276_hypoxia.pdf",height=5,width=5)
ggplot(d1, aes( x = log(NT5E+1), y =log(VEGFA+1))) + geom_point( size=3) + geom_smooth(method = 'lm', formula = y ~ x, se = T)+ 
  theme_classic() +labs(x="log(CD276+1)",y="log(Hypoxia+1)") + stat_cor(method = "pearson",label.x =2, label.y =1.55)

dev.off()



#############################################

dat_1 <- data[,grep("NT5E",colnames(data))]
dat_2 <- data[,grep("Hypoxia",colnames(data))]


plot(dat_1,dat_2)

library(ggpubr)
library(reshape2)
library(ggsci)
d1 <- data.frame(NT5E=dat_1 ,VEGFA=dat_2)

library(Hmisc)

#ggscatterstats(data=d1,x=log(NT5E+1),y=log(VEGFA+1),bf.message = F)


pdf("cor_NT5E_Hypoxia.pdf",height=5,width=5)
ggplot(d1, aes( x = log(NT5E+1), y =log(VEGFA+1))) + geom_point( size=3) + geom_smooth(method = 'lm', formula = y ~ x, se = T)+ 
  theme_classic() +labs(x="log(NT5E+1)",y="log(Hypoxia+1)") + stat_cor(method = "pearson",label.x =1, label.y =1.55)

dev.off()



