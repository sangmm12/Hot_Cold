setwd("D:/R/Hot_tumor_and_cold_tumor/drug1/DEGs_updownRMA/zong/cluster_drug")


cancer <- c("PAAD")



#for(cancer in cancers){

cluster <- read.csv( paste("D:/R/Hot_tumor_and_cold_tumor/Immunecells/CIBERSORT/OVER100/",cancer,"/",cancer,"_clusterresult/",cancer,"_cluster_hc.csv",sep=''),header = F)
rownames(cluster) <- cluster[,1]

library(data.table)
dat <- fread(paste("D:/R/Hot_tumor_and_cold_tumor/DRUG1/DEGs_updownRMA/out_put_",cancer,"_DEGs_updown_zong.csv",sep=''))
dat <- as.data.frame(dat)

rownames(dat) <- dat[,1]
dat<- dat[,-1]

dat_immu <- dat


all_name <- names(which(table(c(rownames(dat_immu),cluster[,1] ))==2))


dat_cluster <- cluster[match(all_name,cluster[,1]),]

dat_imm <- dat_immu[match(all_name,rownames(dat_immu)),]

drug <- c("Dasatinib","Tozasertib")
dat_im <- dat_imm[,match(drug,colnames(dat_imm))]

dat_cluster <- dat_cluster[c('V2')]

dat <- data.frame()
for(coln in colnames(dat_im))
{
  for(i in all_name)
  {
    dat <- rbind(dat,c(coln,dat_cluster[match(i,rownames(dat_cluster)),],dat_im[c(coln)][match(i,rownames(dat_im)),]))
  }
}

dat[,3] = as.numeric(dat[,3])

colnames(dat) <- c("Gene","Group","value")


library(ggpubr)
compare_means(value ~ Group, data = dat, group.by = "Gene",method = "anova")


pdf(paste(cancer,"_DEGs_updown_cluster_Drug2.pdf",sep=''),width=8,height = 8)
p <- ggboxplot(dat, x = "Gene", y = "value",
               color = "Group", palette = c("#2ecc71","#e74c3c"), 
               add = "jitter",x.text.angle=60) # palette可以按照期刊选择相应的配色，如"npg"等
p <- p+xlab("Drug")+ylab("Expression Value")
p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
p
#print(p)
dev.off()

#}


for(cell in unique(dat$Gene)){
  
  dat1 <- subset(dat,dat$Gene==cell)
  
  dat1$Group <- factor(dat1$Group,levels=c("hot","cold"))
  
  
  library(ggpubr)
  compare_means(value ~ Group, data = dat1, group.by = "Gene",method = "anova")
  
  # if(cell=="docetaxel:tanespimycin (2:1 mol/mol)"){
  #   cell1 <- "docetaxel:tanespimycin"
  # }else if(cell=="doxorubicin:navitoclax (2:1 mol/mol)"){
  #   cell1 <- "doxorubicin:navitoclax"
  # }else {
  #   cell1 <- cell
  # }
  # 
  
  pdf(paste("DEGs_updown_cluster_",cell,".pdf",sep=''),width=2.7,height = 6)
  p <- ggboxplot(dat1, x = "Group", y = "value",
                 color = "Group", palette = c("FA7F6F"), 
                 add = "jitter") # palette可以按照期刊选择相应的配色，如"npg"等
  p <- p+xlab(" ")+ylab(paste(cell," (IC50)",sep=''))
  p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
  p <- p + theme(axis.text = element_text(size = 15),axis.title.y = element_text(size = 15))
  # p <- p + scale_x_discrete(limits=c("hot","cold")
  p
  
  print(p)
  dev.off()
  
}
