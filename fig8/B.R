
setwd("D:/R/Hot_tumor_and_cold_tumor/immune data1/gene/fenqi")

cancer <- c("PAAD")

fenqi <- c("Age")

cluster <- read.csv(paste("D:/R/Hot_tumor_and_cold_tumor/immune data1/gene/fenqi/data/",fenqi,".csv",sep=''),header = F)
n1 <- grep(cancer,cluster$V1)
cluster <- cluster[n1,]
cluster <- cluster[,-1:-2]
rownames(cluster) <- cluster$V4


library(data.table)
dat <- fread("D:/R/Hot_tumor_and_cold_tumor/immune data1/gene/exp_NT5E&VEGFA.csv",header = T)
dat <- as.data.frame(dat)

rownames(dat) <- dat[,1]
dat <- dat[,-1:-2]
dat_immu <- as.data.frame(dat)


all_name <- names(which(table(c(rownames(dat_immu),cluster[,2] ))==2))


dat_cluster <- cluster[match(all_name,cluster[,2]),]

dat_im <- dat_immu[match(all_name,rownames(dat_immu)),]

dat_cluster <- dat_cluster[c('V3')]

class(dat_im)
class(dat_cluster)

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
xxxx <- matrix(unique(dat[,1])) #colname

temp_length <- length(xxxx)*2


library(ggpubr)
compare_means(value ~ Group, data = dat, group.by = "Gene",method = "anova")


pdf(paste(fenqi,"_NT5E&CD276.pdf",sep=''),width=temp_length/3.5+2,height = 8)
p <- ggboxplot(dat, x = "Gene", y = "value",
               color = "Group", palette = "FA7F6F", 
               add = "jitter") # palette可以按照期刊选择相应的配色，如"npg"等
p <- p+xlab("Gene")+ylab("Expression Value")
p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
p
#print(p)
dev.off()

for(gene in unique(dat$Gene)){
  
  n1 <- grep(gene,dat$Gene)
  
  dat_1 <- dat[n1,]
  
  xxxx <- matrix(unique(dat[,1])) #colname
  
  temp_length <- length(xxxx)*2
  
  
  library(ggpubr)
  compare_means(value ~ Group, data = dat, group.by = "Gene",method = "anova")
  
  
  pdf(paste(fenqi,"_",gene,".pdf",sep=''),width=temp_length/3.5+2.2,height = 8)
  p <- ggboxplot(dat_1, x = "Group", y = "value",
                 color = "Group", palette = "FA7F6F", 
                 add = "jitter") # palette可以按照期刊选择相应的配色，如"npg"等
  p <- p+xlab(fenqi)+ylab(paste("Expression Value(",fenqi,")",sep=''))
  p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
  p
  print(p)
  dev.off()
  
}







fenqilist <- c("Gender","M","N")
#,"T"

for(fenqi in fenqilist){
  
  
  cluster <- read.csv(paste("D:/R/Hot_tumor_and_cold_tumor/immune data1/gene/fenqi/data/",fenqi,".csv",sep=''),header = F)
  n1 <- grep(cancer,cluster$V1)
  cluster <- cluster[n1,]
  cluster <- cluster[,-1]
  rownames(cluster) <- cluster$V3
  
  
  library(data.table)
  dat <- fread("D:/R/Hot_tumor_and_cold_tumor/immune data1/gene/exp_NT5E&VEGFA.csv",header = T)
  dat <- as.data.frame(dat)
  
  rownames(dat) <- dat[,1]
  dat <- dat[,-1:-2]
  dat_immu <- as.data.frame(dat)
  
  
  all_name <- names(which(table(c(rownames(dat_immu),cluster[,2] ))==2))
  
  
  dat_cluster <- cluster[match(all_name,cluster[,2]),]
  
  dat_im <- dat_immu[match(all_name,rownames(dat_immu)),]
  
  dat_cluster <- dat_cluster[c('V2')]
  
  class(dat_im)
  class(dat_cluster)
  
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
  
  if(fenqi=="N"){
    dat$Group <- factor(dat$Group,levels=c("N0","N1"))
  }
  xxxx <- matrix(unique(dat[,1])) #colname
  
  temp_length <- length(xxxx)*2
  
  
  library(ggpubr)
  compare_means(value ~ Group, data = dat, group.by = "Gene",method = "anova")
  
  
  pdf(paste(fenqi,"_NT5E&VEGFA.pdf",sep=''),width=temp_length/3.5+2,height = 8)
  p <- ggboxplot(dat, x = "Gene", y = "value",
                 color = "Group", palette = "FA7F6F", 
                 add = "jitter") # palette可以按照期刊选择相应的配色，如"npg"等
  p <- p+xlab("Cell Type")+ylab("Expression Value")
  p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
  p
  print(p)
  dev.off()
  
  for(gene in unique(dat$Gene)){
    
    n1 <- grep(gene,dat$Gene)
    
    dat_1 <- dat[n1,]
    
    xxxx <- matrix(unique(dat[,1])) #colname
    
    temp_length <- length(xxxx)*2
    
    
    library(ggpubr)
    compare_means(value ~ Group, data = dat, group.by = "Gene",method = "anova")
    
    
    pdf(paste(fenqi,"_",gene,".pdf",sep=''),width=temp_length/3.5+2.2,height = 8)
    p <- ggboxplot(dat_1, x = "Group", y = "value",
                   color = "Group", palette = "FA7F6F", 
                   add = "jitter") # palette可以按照期刊选择相应的配色，如"npg"等
    p <- p+xlab(fenqi)+ylab(paste("Expression Value(",fenqi,")",sep=''))
    p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
    p
    print(p)
    dev.off()
    
  }
}





fenqilist <- c("Grade")
#,"T"

for(fenqi in fenqilist){
  
  
  cluster <- read.csv(paste("D:/R/Hot_tumor_and_cold_tumor/immune data1/gene/fenqi/data/",fenqi,".csv",sep=''),header = F)
  n1 <- grep(cancer,cluster$V1)
  cluster <- cluster[n1,]
  cluster <- cluster[,-1]
  rownames(cluster) <- cluster$V3
  
  
  library(data.table)
  dat <- fread("D:/R/Hot_tumor_and_cold_tumor/immune data1/gene/exp_NT5E&VEGFA.csv",header = T)
  dat <- as.data.frame(dat)
  
  rownames(dat) <- dat[,1]
  dat <- dat[,-1:-2]
  dat_immu <- as.data.frame(dat)
  
  
  all_name <- names(which(table(c(rownames(dat_immu),cluster[,2] ))==2))
  
  
  dat_cluster <- cluster[match(all_name,cluster[,2]),]
  
  dat_im <- dat_immu[match(all_name,rownames(dat_immu)),]
  
  dat_cluster <- dat_cluster[c('V2')]
  
  class(dat_im)
  class(dat_cluster)
  
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
  
  dat$Group <- factor(dat$Group,levels=c("G1","G2","G3"))
  
  xxxx <- matrix(unique(dat[,1])) #colname
  
  temp_length <- length(xxxx)*2
  
  
  library(ggpubr)
  compare_means(value ~ Group, data = dat, group.by = "Gene",method = "anova")
  
  
  pdf(paste(fenqi,"_NT5E&VEGFA.pdf",sep=''),width=temp_length/3.5+2,height = 8)
  p <- ggboxplot(dat, x = "Gene", y = "value",
                 color = "Group", palette = "FA7F6F", 
                 add = "jitter") # palette可以按照期刊选择相应的配色，如"npg"等
  p <- p+xlab("Cell Type")+ylab("Expression Value")
  p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
  p
  print(p)
  dev.off()
  
  for(gene in unique(dat$Gene)){
    
    n1 <- grep(gene,dat$Gene)
    
    dat_1 <- dat[n1,]
    
    xxxx <- matrix(unique(dat[,1])) #colname
    
    temp_length <- length(xxxx)*2
    
    
    library(ggpubr)
    compare_means(value ~ Group, data = dat, group.by = "Gene",method = "anova")
    
    
    pdf(paste(fenqi,"_",gene,".pdf",sep=''),width=temp_length/3.5+3,height = 8)
    p <- ggboxplot(dat_1, x = "Group", y = "value",
                   color = "Group", palette = "FA7F6F", 
                   add = "jitter") # palette可以按照期刊选择相应的配色，如"npg"等
    p <- p+xlab(fenqi)+ylab(paste("Expression Value(",fenqi,")",sep=''))
    p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
    p
    print(p)
    dev.off()
    
  }
}







fenqilist <- c("Stage")
#,"T"

for(fenqi in fenqilist){
  
  
  cluster <- read.csv(paste("D:/R/Hot_tumor_and_cold_tumor/immune data1/gene/fenqi/data/",fenqi,".csv",sep=''),header = F)
  n1 <- grep(cancer,cluster$V1)
  cluster <- cluster[n1,]
  cluster <- cluster[,-1]
  rownames(cluster) <- cluster$V3
  
  
  library(data.table)
  dat <- fread("D:/R/Hot_tumor_and_cold_tumor/immune data1/gene/exp_NT5E&VEGFA.csv",header = T)
  dat <- as.data.frame(dat)
  
  rownames(dat) <- dat[,1]
  dat <- dat[,-1:-2]
  dat_immu <- as.data.frame(dat)
  
  
  all_name <- names(which(table(c(rownames(dat_immu),cluster[,2] ))==2))
  
  
  dat_cluster <- cluster[match(all_name,cluster[,2]),]
  
  dat_im <- dat_immu[match(all_name,rownames(dat_immu)),]
  
  dat_cluster <- dat_cluster[c('V2')]
  
  class(dat_im)
  class(dat_cluster)
  
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
  
  dat$Group <- factor(dat$Group,levels=c("Stage I","Stage II","Stage III","Stage IV"))
  
  xxxx <- matrix(unique(dat[,1])) #colname
  
  temp_length <- length(xxxx)*2
  
  
  library(ggpubr)
  compare_means(value ~ Group, data = dat, group.by = "Gene",method = "anova")
  
  
  pdf(paste(fenqi,"_NT5E&VEGFA.pdf",sep=''),width=temp_length/3.5+2,height = 8)
  p <- ggboxplot(dat, x = "Gene", y = "value",
                 color = "Group", palette = "FA7F6F", 
                 add = "jitter") # palette可以按照期刊选择相应的配色，如"npg"等
  p <- p+xlab("Cell Type")+ylab("Expression Value")
  p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
  p
  print(p)
  dev.off()
  
  for(gene in unique(dat$Gene)){
    
    n1 <- grep(gene,dat$Gene)
    
    dat_1 <- dat[n1,]
    
    xxxx <- matrix(unique(dat[,1])) #colname
    
    temp_length <- length(xxxx)*2
    
    
    library(ggpubr)
    compare_means(value ~ Group, data = dat, group.by = "Gene",method = "anova")
    
    
    pdf(paste(fenqi,"_",gene,".pdf",sep=''),width=temp_length/3.5+4,height = 8)
    p <- ggboxplot(dat_1, x = "Group", y = "value",
                   color = "Group", palette = "FA7F6F", 
                   add = "jitter") # palette可以按照期刊选择相应的配色，如"npg"等
    p <- p+xlab(fenqi)+ylab(paste("Expression Value(",fenqi,")",sep=''))
    p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
    p
    print(p)
    dev.off()
    
  }
}


