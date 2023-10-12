setwd("D:/R/Hot_tumor_and_cold_tumor/drug1/DEGs_updownRMA/zong/cor_gene")

cancer <- c("PAAD")

library(data.table)
dat <- fread(paste("D:/R/Hot_tumor_and_cold_tumor/DRUG1/DEGs_updownRMA/out_put_",cancer,"_DEGs_updown_zong.csv",sep=''))
dat <- as.data.frame(dat)


rownames(dat) <- dat[,1]
dat_gene <- dat[,-1]


is.numeric(dat_gene)



data <- fread(paste("D:/R/Hot_tumor_and_cold_tumor/immune data1/immu_",cancer,"_zong.csv",sep=''),header = T)
data <- as.data.frame(data)

genename <- data[,1]
data <- data[,-1]
data=as.data.frame(lapply(data,as.numeric),check.names=F)
rownames(data) <- genename




p_value <- fread(paste("D:/R/Hot_tumor_and_cold_tumor/immune data1/cluster_gene/p_value.csv",sep=''),header = T)
up_down <- fread(paste("D:/R/Hot_tumor_and_cold_tumor/immune data1/cluster_gene/updown.csv",sep=''),header = T)
p_value <- as.data.frame(p_value)
up_down <- as.data.frame(up_down)


rownames(p_value) <- p_value$V1
p_value <- p_value[,-1]

rownames(up_down) <- up_down[,1]
up_down <- up_down[,-1]

up_down <- up_down[match(rownames(p_value),rownames(up_down)),]

for(i in 1:length(rownames(up_down))){
  for(t in 1:length(colnames(up_down))){
    p_value[i,t] <- as.numeric(p_value[i,t])*as.numeric(up_down[i,t])
  }
}


DEG_BLCA_up <-  rownames(p_value)[which((p_value$BLCA>0) & (p_value$BLCA<0.05) )]
DEG_BLCA_down <-  rownames(p_value)[which((p_value$BLCA<0) & (p_value$BLCA>-0.05) )]

DEG_CESC_up <-  rownames(p_value)[which((p_value$CESC>0) & (p_value$CESC<0.05) )]
DEG_CESC_down <-  rownames(p_value)[which((p_value$CESC<0) & (p_value$CESC>-0.05) )]

DEG_PAAD_up <-  rownames(p_value)[which((p_value$PAAD>0) & (p_value$PAAD<0.05) )]
DEG_PAAD_down <-  rownames(p_value)[which((p_value$PAAD<0) & (p_value$PAAD>-0.05) )]

DEG_SARC_up <-  rownames(p_value)[which((p_value$SARC>0) & (p_value$SARC<0.05) )]
DEG_SARC_down <-  rownames(p_value)[which((p_value$SARC<0) & (p_value$SARC>-0.05) )]

DEG_SKCM_up <-  rownames(p_value)[which((p_value$SKCM>0) & (p_value$SKCM<0.05) )]
DEG_SKCM_down <-  rownames(p_value)[which((p_value$SKCM<0) & (p_value$SKCM>-0.05) )]


genelist <- c(DEG_PAAD_up,DEG_PAAD_down)
genelist <- unique(genelist)



dat_im <- data[match(genelist,rownames(data)),]
dat_im <- t(dat_im)

all_name <- names(which(table(c(rownames(dat_im),rownames(dat_gene) ))==2))


dat_gene <- dat_gene[match(all_name,rownames(dat_gene)),]
dat_im <- dat_im[match(all_name,rownames(dat_im)),]

#############################################

dat_1 <- dat_gene[,grep("Dasatinib",colnames(dat_gene))]
dat_2 <- dat_im[,grep("CCL14",colnames(dat_im))]


plot(dat_1,dat_2)

library(ggpubr)
library(reshape2)
library(ggsci)
d1 <- data.frame(NT5E=dat_1 ,VEGFA=dat_2)

library(Hmisc)

#ggscatterstats(data=d1,x=log(NT5E+1),y=log(VEGFA+1),bf.message = F)


pdf("cor_Dasatinib_CCL14.pdf",height=5,width=5)
ggplot(d1, aes( x = NT5E, y =VEGFA)) + geom_point( size=3) + geom_smooth(method = 'lm', formula = y ~ x, se = T)+ 
  theme_classic() +labs(x="Dasatinib",y="CCL14") + stat_cor(method = "pearson",label.x =min(dat_1), label.y =max(dat_2))

dev.off()




dat_1 <- dat_gene[,grep("Dasatinib",colnames(dat_gene))]
dat_2 <- dat_im[,grep("TGFB1",colnames(dat_im))]


plot(dat_1,dat_2)

library(ggpubr)
library(reshape2)
library(ggsci)
d1 <- data.frame(NT5E=dat_1 ,VEGFA=dat_2)

library(Hmisc)

#ggscatterstats(data=d1,x=log(NT5E+1),y=log(VEGFA+1),bf.message = F)


pdf("cor_Dasatinib_TGFB1.pdf",height=5,width=5)
ggplot(d1, aes( x = NT5E, y =VEGFA)) + geom_point( size=3) + geom_smooth(method = 'lm', formula = y ~ x, se = T)+ 
  theme_classic() +labs(x="Dasatinib",y="TGFB1") + stat_cor(method = "pearson",label.x =min(dat_1), label.y =max(dat_2))

dev.off()





dat_1 <- dat_gene[,grep("Tozasertib",colnames(dat_gene))]
dat_2 <- dat_im[,grep("CCL14",colnames(dat_im))]


plot(dat_1,dat_2)

library(ggpubr)
library(reshape2)
library(ggsci)
d1 <- data.frame(NT5E=dat_1 ,VEGFA=dat_2)

library(Hmisc)

#ggscatterstats(data=d1,x=log(NT5E+1),y=log(VEGFA+1),bf.message = F)


pdf("cor_Tozasertib_CCL14.pdf",height=5,width=5)
ggplot(d1, aes( x = NT5E, y =VEGFA)) + geom_point( size=3) + geom_smooth(method = 'lm', formula = y ~ x, se = T)+ 
  theme_classic() +labs(x="Tozasertib",y="CCL14") + stat_cor(method = "pearson",label.x =min(dat_1), label.y =max(dat_2))

dev.off()


dat_1 <- dat_gene[,grep("Tozasertib",colnames(dat_gene))]
dat_2 <- dat_im[,grep("TNFSF4",colnames(dat_im))]


plot(dat_1,dat_2)

library(ggpubr)
library(reshape2)
library(ggsci)
d1 <- data.frame(NT5E=dat_1 ,VEGFA=dat_2)

library(Hmisc)

#ggscatterstats(data=d1,x=log(NT5E+1),y=log(VEGFA+1),bf.message = F)


pdf("cor_Tozasertib_TNFSF4.pdf",height=5,width=5)
ggplot(d1, aes( x = NT5E, y =VEGFA)) + geom_point( size=3) + geom_smooth(method = 'lm', formula = y ~ x, se = T)+ 
  theme_classic() +labs(x="Tozasertib",y="TNFSF4") + stat_cor(method = "pearson",label.x =min(dat_1), label.y =max(dat_2))

dev.off()