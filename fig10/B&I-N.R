setwd("D:/R/Hot_tumor_and_cold_tumor/immune data1/cor_estimate$gene")

cancer <- c("PAAD")

library(data.table)
dat <- fread(paste("D:/R/Hot_tumor_and_cold_tumor/cluster_immuscore/Estimate_",cancer,".csv",sep=''))
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

# for(i in 1:length(colnames(dat_gene)))
# {
#   dat_gene[,i] = as.numeric(unlist(dat_gene[,i]))
# }
# i=1
# nrow(dat_gene)
# ncol(dat_gene)
# 

colSums(dat_im)


library(psych)
data.corr <- corr.test(dat_im, dat_gene, method="pearson", adjust="fdr")
data.r <- data.corr$r  # 相关系数
data.p <- data.corr$p  # p值

paste("data.r_",cancer,".csv",sep='')
write.csv(data.r,file= paste("data.r_xcell.csv",sep=''),quote=F)
write.csv(data.p,file= paste("data.p_xcell.csv",sep=''),quote=F)


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


pdf(paste("cor_PAAD_DEGs_estimate.pdf",sep=''),width =2.2,height = 16)
pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F, cluster_cols=F,display_numbers=sig.mat
)

dev.off()





dat_2 <- dat_im[,grep("HLA-DPB1",colnames(dat_im))]
dat_1 <- dat_gene[,grep("ESTIMATEScore",colnames(dat_gene))]


plot(dat_1,dat_2)

library(ggpubr)
library(reshape2)
library(ggsci)
d1 <- data.frame(NT5E=dat_1 ,VEGFA=dat_2)

library(ggplot2)
library(ggpubr)
library(ggExtra)

#ggscatterstats(data=d1,x=log(NT5E+1),y=log(VEGFA+1),bf.message = F)



p1 <- ggplot(d1, aes( x = NT5E, y =VEGFA)) + geom_point( size=3) + geom_smooth(method = 'lm', formula = y ~ x, se = T)+ 
      theme_classic() +labs(x="ESTIMATEScore",y="HLA-DPB1") + stat_cor(method = "pearson",label.x =min(dat_1), label.y =max(dat_2))
p <- ggMarginal(p1, type = "density", xparams = list(fill = "orange"),yparams = list(fill = "blue"))
pdf("cor_HLA-DPB1_ESTIMATEScore1.pdf",height=5,width=5)
print(p)
dev.off()




dat_2 <- dat_im[,grep("VEGFA",colnames(dat_im))]
dat_1 <- dat_gene[,grep("ESTIMATEScore",colnames(dat_gene))]


plot(dat_1,dat_2)

library(ggpubr)
library(reshape2)
library(ggsci)
d1 <- data.frame(NT5E=dat_1 ,VEGFA=dat_2)

library(ggplot2)
library(ggpubr)
library(ggExtra)

#ggscatterstats(data=d1,x=log(NT5E+1),y=log(VEGFA+1),bf.message = F)


p1 <- ggplot(d1, aes( x = NT5E, y =VEGFA)) + geom_point( size=3) + geom_smooth(method = 'lm', formula = y ~ x, se = T)+ 
  theme_classic() +labs(x="ESTIMATEScore",y="VEGFA") + stat_cor(method = "pearson",label.x =min(dat_1), label.y =max(dat_2))
p <- ggMarginal(p1, type = "density", xparams = list(fill = "orange"),yparams = list(fill = "blue"))
pdf("cor_VEGFA_ESTIMATEScore.pdf",height=5,width=5)
print(p)
dev.off()




dat_2 <- dat_im[,grep("HLA-DPB1",colnames(dat_im))]
dat_1 <- dat_gene[,grep("ImmuneScore",colnames(dat_gene))]


plot(dat_1,dat_2)

library(ggpubr)
library(reshape2)
library(ggsci)
d1 <- data.frame(NT5E=dat_1 ,VEGFA=dat_2)

library(ggplot2)
library(ggpubr)
library(ggExtra)

#ggscatterstats(data=d1,x=log(NT5E+1),y=log(VEGFA+1),bf.message = F)



p1 <- ggplot(d1, aes( x = NT5E, y =VEGFA)) + geom_point( size=3) + geom_smooth(method = 'lm', formula = y ~ x, se = T)+ 
  theme_classic() +labs(x="ImmuneScore",y="HLA-DPB1") + stat_cor(method = "pearson",label.x =min(dat_1), label.y =max(dat_2))
p <- ggMarginal(p1, type = "density", xparams = list(fill = "orange"),yparams = list(fill = "blue"))
pdf("cor_HLA-DPB1_ImmuneScore.pdf",height=5,width=5)
print(p)
dev.off()




dat_1 <- dat_im[,grep("VEGFA",colnames(dat_im))]
dat_2 <- dat_gene[,grep("ImmuneScore",colnames(dat_gene))]


plot(dat_1,dat_2)

library(ggpubr)
library(reshape2)
library(ggsci)
d1 <- data.frame(NT5E=dat_1 ,VEGFA=dat_2)

library(ggplot2)
library(ggpubr)
library(ggExtra)

#ggscatterstats(data=d1,x=log(NT5E+1),y=log(VEGFA+1),bf.message = F)


p1 <- ggplot(d1, aes( x = NT5E, y =VEGFA)) + geom_point( size=3) + geom_smooth(method = 'lm', formula = y ~ x, se = T)+ 
  theme_classic() +labs(x="ImmuneScore",y="VEGFA") + stat_cor(method = "pearson",label.x =min(dat_1), label.y =max(dat_2))
p <- ggMarginal(p1, type = "density", xparams = list(fill = "orange"),yparams = list(fill = "blue"))
pdf("cor_VEGFA_ImmuneScore.pdf",height=5,width=5)
print(p)
dev.off()





dat_2 <- dat_im[,grep("HLA-DPB1",colnames(dat_im))]
dat_1 <- dat_gene[,grep("StromalScore",colnames(dat_gene))]


plot(dat_1,dat_2)

library(ggpubr)
library(reshape2)
library(ggsci)
d1 <- data.frame(NT5E=dat_1 ,VEGFA=dat_2)

library(ggplot2)
library(ggpubr)
library(ggExtra)

#ggscatterstats(data=d1,x=log(NT5E+1),y=log(VEGFA+1),bf.message = F)



p1 <- ggplot(d1, aes( x = NT5E, y =VEGFA)) + geom_point( size=3) + geom_smooth(method = 'lm', formula = y ~ x, se = T)+ 
  theme_classic() +labs(x="StromalScore",y="HLA-DPB1") + stat_cor(method = "pearson",label.x =min(dat_1), label.y =max(dat_2))
p <- ggMarginal(p1, type = "density", xparams = list(fill = "orange"),yparams = list(fill = "blue"))
pdf("cor_HLA-DPB1_StromalScore.pdf",height=5,width=5)
print(p)
dev.off()




dat_2 <- dat_im[,grep("VEGFA",colnames(dat_im))]
dat_1 <- dat_gene[,grep("StromalScore",colnames(dat_gene))]


plot(dat_1,dat_2)

library(ggpubr)
library(reshape2)
library(ggsci)
d1 <- data.frame(NT5E=dat_1 ,VEGFA=dat_2)

library(ggplot2)
library(ggpubr)
library(ggExtra)

#ggscatterstats(data=d1,x=log(NT5E+1),y=log(VEGFA+1),bf.message = F)


p1 <- ggplot(d1, aes( x = NT5E, y =VEGFA)) + geom_point( size=3) + geom_smooth(method = 'lm', formula = y ~ x, se = T)+ 
  theme_classic() +labs(x="StromalScore",y="VEGFA") + stat_cor(method = "pearson",label.x =min(dat_1), label.y =max(dat_2))
p <- ggMarginal(p1, type = "density", xparams = list(fill = "orange"),yparams = list(fill = "blue"))
pdf("cor_VEGFA_StromalScore.pdf",height=5,width=5)
print(p)
dev.off()

