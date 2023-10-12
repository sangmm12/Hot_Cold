rm(list=ls())

setwd("D:/R/Hot_tumor_and_cold_tumor/cluster_cell/cluster_cellzong/pheatmap")

p_value <- c()

cancers <- c("BLCA","CESC","KIRP","LGG","PAAD","SARC","SKCM","THYM")

#cancer <- c("BLCA")
for(cancer in cancers){
  
  cluster <- read.csv(paste("D:/R/Hot_tumor_and_cold_tumor/Immunecells/CIBERSORT/OVER100/",cancer,"/",cancer,"_clusterresult/",cancer,"_cluster_hc.csv",sep=''),header = F)
  rownames(cluster) <- cluster[,1]
  
  library(data.table)
  dat_imm <- fread(paste("D:/R/Hot_tumor_and_cold_tumor/cluster_cell/cluster_cellzong/",cancer,"/data_",cancer,".csv",sep=''),header = T)
  dat_imm <- as.data.frame(dat_imm)
  rownames(dat_imm) <- dat_imm[,1]
  dat_imm <- dat_imm[,-1]

  
  all_name <- names(which(table(c(rownames(dat_imm),cluster[,1] ))==2))
  
  dat_cluster <- cluster[match(all_name,cluster[,1]),]
  
  dat_im <- dat_imm[match(all_name,rownames(dat_imm)),]
  
  
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
  
  xxxx <- matrix(unique(dat[,1])) #colname
  
  temp_length <- length(xxxx)*2
  
  
  library(ggpubr)
  p_gene <- compare_means(value ~ Group, data = dat, group.by = "Gene",method = "anova")
  
  p_value <- cbind(p_value,cancer=p_gene$p)
  rownames(p_value) <- p_gene$Gene
  
} 

colnames(p_value) <- cancers

write.csv(p_value,file= "D:/R/Hot_tumor_and_cold_tumor/cluster_cell/cluster_cellzong/p_value.csv",quote=F)












####################################


library(data.table)
p_value <- fread("D:/R/Hot_tumor_and_cold_tumor/cluster_cell/cluster_cellzong/p_value.csv",header = T)
p_value <- as.data.frame(p_value)

rownames(p_value) <- p_value[,1]
p_value <- p_value[,-1]


#p_value[10,8] <- 0

up_down <- fread("D:/R/Hot_tumor_and_cold_tumor/cluster_cell/cluster_cellzong/updown.csv",header = T)
up_down <- as.data.frame(up_down)

up_down[up_down == "Check_point"] = "Check-point"
up_down[up_down == "Inflammation_promoting"] = "Inflammation-promoting"
up_down[up_down == "T_cell_co_inhibition"] = "T_cell_co-inhibition"
up_down[up_down == "T_cell_co_stimulation"] = "T_cell_co-stimulation"

rownames(up_down) <- up_down[,1]
up_down <- up_down[,-1]

#list <- c("T_cells_CD8","NK_cells_activated","Macrophages_M2","Cytolytic_activity","Tcells_proliferation")


df <- -log10(p_value)

#df <- t(scale(t(df)))
#df <- scale(df)

for(i in 1:length(rownames(up_down))){
  for(t in 1:8){
    df[i,t] <- as.numeric(df[i,t])*as.numeric(up_down[i,t])
  }
}

# df <- df[,-3:-4]
# df <- df[,-6]
# 
# p_value <- p_value[,-3:-4]
# p_value <- p_value[,-6]

p_value[10,8] <- 100
df[10,8] <- 0

cancerslist <- c("BLCA","CESC","PAAD","SARC","SKCM","KIRP","LGG","THYM")

df <- df[,match(cancerslist,colnames(df))]
p_value <- p_value[,match(cancerslist,colnames(p_value))]


library(pheatmap)

pdf("cluster_CIBERSORT.pdf",5,length(rownames(df))/2)


paletteLength = 1000
#immune
#myColor <- colorRampPalette(c("white", "#FF7C00"))(paletteLength)
#exp
#myColor <- colorRampPalette(c("white", "red"))(paletteLength)
#cell
#myColor <- colorRampPalette(c("white","blue"))(paletteLength)
#drug
#myColor <- colorRampPalette(c("white", "#660BAB"))(paletteLength)
#yzx_gx
#myColor <- colorRampPalette(c("white", "#C7007D"))(paletteLength)
#df[10,8] <- 0
max(df)
bk <- c(seq( -max(abs(max(df)),abs(min(df))), -0.1,by=0.01),seq(0,max(abs(max(df)),abs(min(df))),by=0.01))
myColor <- c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2))

#myBreaks <- c(seq( min(df),-min(df), length.out=floor(paletteLength/2)))

#bk <- c(seq(-9,-0.1,by=0.01),seq(0,9,by=0.01))

#######################################
getSig <- function(dc) {
  sc <- ' '
  if (dc < 0.0001) {sc <- '****'}
  else if (dc < 0.001){sc <- '***'}
  else if (dc < 0.01){sc <- '**'}
  else if (dc < 0.05) {sc <- '*'}
  else{sc <- ''}
  return(sc)
}

#p_value[10,8] <- 0
sig.mat <- matrix(sapply(as.matrix(p_value), getSig), nrow=nrow(as.matrix(p_value)))
str(sig.mat)
########################################


xx <- pheatmap(df,
               color=myColor,
               breaks=bk,
               clustering_method="average", cluster_rows=F,cluster_cols=F, cellwidth = 20,cellheight = 20,main="-log10(p)",display_numbers=sig.mat)
print(xx)
dev.off()












####################################



p_value <- fread("D:/R/Hot_tumor_and_cold_tumor/cluster_cell/cluster_cellzong/p_value.csv",header = T)
p_value <- as.data.frame(p_value)

rownames(p_value) <- p_value[,1]
p_value <- p_value[,-1]


#p_value[10,8] <- 0

up_down <- fread("D:/R/Hot_tumor_and_cold_tumor/cluster_cell/cluster_cellzong/updown.csv",header = T)
up_down <- as.data.frame(up_down)

up_down[up_down == "Check_point"] = "Check-point"
up_down[up_down == "Inflammation_promoting"] = "Inflammation-promoting"
up_down[up_down == "T_cell_co_inhibition"] = "T_cell_co-inhibition"
up_down[up_down == "T_cell_co_stimulation"] = "T_cell_co-stimulation"

rownames(up_down) <- up_down[,1]
up_down <- up_down[,-1]

list <- c("T_cells_CD8","NK_cells_activated","Macrophages_M2","Cytolytic_activity","T_cell_co-inhibition")


up_down <- up_down[match(list,rownames(up_down)),]
p_value <- p_value[match(list,rownames(p_value)),]

p_os <- c(0.00021,0.02,0.03,0.01,0.04,0.05,0.000058,0.04)
up_os <- c(1,1,-1,1,1,1,1,1)

up_down <- rbind(up_down,up_os)
p_value <- rbind(p_value,p_os)

list <- c("T_cells_CD8","NK_cells_activated","Macrophages_M2","Cytolytic_activity","T_cell_co-inhibition","Survival")
rownames(up_down) <- list
rownames(p_value) <- list

df <- -log10(p_value)

#df <- t(scale(t(df)))
#df <- scale(df)

for(i in 1:length(rownames(up_down))){
  for(t in 1:8){
    df[i,t] <- as.numeric(df[i,t])*as.numeric(up_down[i,t])
  }
}



library(pheatmap)

pdf("cluster_Inhibitory.pdf",5,length(rownames(df))/2)


paletteLength = 1000
#immune
#myColor <- colorRampPalette(c("white", "#FF7C00"))(paletteLength)
#exp
#myColor <- colorRampPalette(c("white", "red"))(paletteLength)
#cell
#myColor <- colorRampPalette(c("white","blue"))(paletteLength)
#drug
#myColor <- colorRampPalette(c("white", "#660BAB"))(paletteLength)
#yzx_gx
#myColor <- colorRampPalette(c("white", "#C7007D"))(paletteLength)
#df[10,8] <- 0
max(df)
bk <- c(seq( -max(df), -0.1,by=0.01),seq(0,max(df),by=0.01))
myColor <- c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2))

#myBreaks <- c(seq( min(df),-min(df), length.out=floor(paletteLength/2)))

#bk <- c(seq(-9,-0.1,by=0.01),seq(0,9,by=0.01))

#######################################
getSig <- function(dc) {
  sc <- ' '
  if (dc < 0.0001) {sc <- '****'}
  else if (dc < 0.001){sc <- '***'}
  else if (dc < 0.01){sc <- '**'}
  else if (dc < 0.05) {sc <- '*'}
  else{sc <- ''}
  return(sc)
}

#p_value[10,8] <- 0
sig.mat <- matrix(sapply(as.matrix(p_value), getSig), nrow=nrow(as.matrix(p_value)))
str(sig.mat)












########################################


xx <- pheatmap(df,
               color=myColor,
               breaks=bk,
               clustering_method="average", cluster_rows=F,cluster_cols=F, cellwidth = 20,cellheight = 20,main="-log10(p)",display_numbers=sig.mat)
print(xx)
dev.off()



df <- df[,-3:-4]
df <- df[,-6]

p_value <- p_value[,-3:-4]
p_value <- p_value[,-6]


pdf("cluster_Inhibitory_5.pdf",5,length(rownames(df))/2)


paletteLength = 1000
#immune
#myColor <- colorRampPalette(c("white", "#FF7C00"))(paletteLength)
#exp
#myColor <- colorRampPalette(c("white", "red"))(paletteLength)
#cell
#myColor <- colorRampPalette(c("white","blue"))(paletteLength)
#drug
#myColor <- colorRampPalette(c("white", "#660BAB"))(paletteLength)
#yzx_gx
#myColor <- colorRampPalette(c("white", "#C7007D"))(paletteLength)
#df[10,8] <- 0
max(df)
bk <- c(seq( -max(df), -0.1,by=0.01),seq(0,max(df),by=0.01))
myColor <- c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2))

#myBreaks <- c(seq( min(df),-min(df), length.out=floor(paletteLength/2)))

#bk <- c(seq(-9,-0.1,by=0.01),seq(0,9,by=0.01))

#######################################
getSig <- function(dc) {
  sc <- ' '
  if (dc < 0.0001) {sc <- '****'}
  else if (dc < 0.001){sc <- '***'}
  else if (dc < 0.01){sc <- '**'}
  else if (dc < 0.05) {sc <- '*'}
  else{sc <- ''}
  return(sc)
}

#p_value[10,8] <- 0
sig.mat <- matrix(sapply(as.matrix(p_value), getSig), nrow=nrow(as.matrix(p_value)))
str(sig.mat)
########################################


xx <- pheatmap(df,
               color=myColor,
               breaks=bk,
               clustering_method="average", cluster_rows=F,cluster_cols=F, cellwidth = 20,cellheight = 20,main="-log10(p)",display_numbers=sig.mat)
print(xx)
dev.off()

