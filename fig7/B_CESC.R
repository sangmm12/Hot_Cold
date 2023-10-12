

setwd("D:/R/Hot_tumor_and_cold_tumor/immune data1/table")



cells <- c(	"B_cells_naive","B_cells_memory","Plasma_cells","T_cells_CD8","T_cells_CD4_naive","T_cells_CD4_memory_resting","T_cells_CD4_memory_activated","T_cells_follicular_helper","T_cells_regulatory_.Tregs.","T_cells_gamma_delta","NK_cells_resting","NK_cells_activated","Monocytes","Macrophages_M0","Macrophages_M1","Macrophages_M2","Dendritic_cells_resting","Dendritic_cells_activated","Mast_cells_resting","Mast_cells_activated","Eosinophils","Neutrophils")


#cancers <- c("BLCA","CESC","PAAD","SARC")
immunames <- c("Inhibitory","Stimulaotry","chemokine","MHC","receptor","Immunoinhibitor","Immunostimulator")
#immuname <- c("Inhibitory")



cancer <- c("CESC")

#for(cancer in cancers){


data_BLCA_p <- c()
data_BLCA_r <- c()

if(dir.exists(cancer)==TRUE){
  print(cancer)
}else {
  dir.create(cancer)
}

library(data.table)

for(immuname in immunames){
  
  datp <- fread(paste("D:/R/Hot_tumor_and_cold_tumor/immune data1/cor_cell$gene/data.p_",cancer,"_",immuname,".csv",sep=''),header = T)
  datr <- fread(paste("D:/R/Hot_tumor_and_cold_tumor/immune data1/cor_cell$gene/data.r_",cancer,"_",immuname,".csv",sep=''),header = T)
  
  data_BLCA_p <- rbind(data_BLCA_p,datp)
  data_BLCA_r <- rbind(data_BLCA_r,datr)
}

data_BLCA_p <- as.matrix(data_BLCA_p)
data_BLCA_r <- as.matrix(data_BLCA_r)


data_BLCA_p <- data_BLCA_p[!duplicated(data_BLCA_p[,1]),]
data_BLCA_r <- data_BLCA_r[!duplicated(data_BLCA_r[,1]),]



write.csv(data_BLCA_p,file=paste(cancer,"/data_",cancer,"_p.csv",sep=''),row.names = F)
write.csv(data_BLCA_r,file=paste(cancer,"/data_",cancer,"_r.csv",sep=''),row.names = F)



data_BLCA_DEG <- fread(paste("D:/R/Hot_tumor_and_cold_tumor/immune data1/cluster_gene/p_value.csv",sep=''),header = T)
data_BLCA_DEG <- as.matrix(data_BLCA_DEG)


# aa <-colnames(data_BLCA_r)
# aa <- aa[-1]
cells <- c(	"B_cells_naive","B_cells_memory","Plasma_cells","T_cells_CD8","T_cells_CD4_naive","T_cells_CD4_memory_resting","T_cells_CD4_memory_activated","T_cells_follicular_helper","T_cells_regulatory_.Tregs.","T_cells_gamma_delta","NK_cells_resting","NK_cells_activated","Monocytes","Macrophages_M0","Macrophages_M1","Macrophages_M2","Dendritic_cells_resting","Dendritic_cells_activated","Mast_cells_resting","Mast_cells_activated","Eosinophils","Neutrophils")
#cells <- c(cells[4],cells[8],cells[16],cells[12])

for(i in 1:length(cells)){
  # i=1
  #i=8
  # i=12
  # i=16
  cell <- cells[i]
  
  
  lin_name <- names(which(table(c(data_BLCA_p[,1],data_BLCA_r[,1]))==2))
  
  dat_p <- data_BLCA_p[,c(1,i+1)]
  dat_r <- data_BLCA_r[,c(1,i+1)]
  
  
  dat_p <- dat_p[match(lin_name,dat_p[,1]),]
  dat_r <- dat_r[match(lin_name,dat_r[,1]),]
  
  
  dat_line <- data.frame(Gene=lin_name,
                         p=as.numeric(dat_p[,2]),
                         r=as.numeric(dat_r[,2]))
  
  
  data_positive <- dat_line[which(dat_line$p < 0.05 & dat_line$r > 0.2),]
  
  data_negative <- dat_line[which(dat_line$p < 0.05 & dat_line$r < -0.2),]
  
  data1_BLCA_DEG <- data.frame(Gene=data_BLCA_DEG[,1],
                               p=as.numeric(data_BLCA_DEG[,3]))
  
  
  
  data2_BLCA_DEG <- data1_BLCA_DEG[which(as.numeric(data1_BLCA_DEG[,2]) < 0.05),1:2]
  
  
  write.csv(data_positive,file=paste(cancer,"/",cancer,"_",cell,"_positive.csv",sep=''),row.names = F)
  write.csv(data_negative,file=paste(cancer,"/",cancer,"_",cell,"_negative.csv",sep=''),row.names = F)
  write.csv(data2_BLCA_DEG,file=paste(cancer,"/",cancer,"_",cell,"_DEGs.csv",sep=''),row.names = F)
  
  
  
  Gene_po_DEG <- names(which(table(c(data_positive$Gene,data2_BLCA_DEG$Gene))==2))
  Gene_ne_DEG <- names(which(table(c(data_negative$Gene,data2_BLCA_DEG$Gene))==2))
  
  # write.table(Gene_po_DEG,file=paste(cancer,"Gene_po_DEG.txt",sep=''),quote=F)
  # write.table(Gene_ne_DEG,file=paste(cancer,"Gene_ne_DEG.txt",sep=''),quote=F)
  
  
  
  
  data_positive1 <- data_positive[match(Gene_po_DEG,data_positive$Gene),]
  
  data_negative1 <- data_negative[match(Gene_ne_DEG,data_negative$Gene),]
  
  
  data2_BLCA_DEG_po <- data2_BLCA_DEG[match(Gene_po_DEG,data2_BLCA_DEG$Gene),]
  
  data2_BLCA_DEG_ne<- data2_BLCA_DEG[match(Gene_ne_DEG,data2_BLCA_DEG$Gene),]
  
  Gene_po_DEG <- data.frame(data_positive1,P.DEG=data2_BLCA_DEG_po[,2])
  Gene_ne_DEG <- data.frame(data_negative1,P.DEG=data2_BLCA_DEG_ne[,2])
  
  write.csv(Gene_po_DEG,file=paste(cancer,"/",cancer,"_",cell,"_positive_DEG.csv",sep=''),row.names = F)
  write.csv(Gene_ne_DEG,file=paste(cancer,"/",cancer,"_",cell,"_negative_DEG.csv",sep=''),row.names = F)
  
  
  
  veen_positive <- c(data_positive$Gene)
  veen_negative <- c(data_negative$Gene)
  veen_DEG <- c(data2_BLCA_DEG$Gene)
  
  library(VennDiagram)
  # 
  # venn.diagram(list(Positive=veen_positive,Negative=veen_negative,DEGs=veen_DEG),
  #              fill=c("#729ECE","#FF9E4A","#67BF5C"), #,"#67BF5C","#ED665D","#AD8BC9"
  #              alpha=c(0.5,0.5,0.5),
  #              col = c("#729ECE","#FF9E4A","#67BF5C"),
  #              cex=1,
  #              # cat.pos = 10,
  #              cat.pos = c(180,180,180),
  #              cat.dist = c(0.03,0.01,0.05),
  #              cat.fontface=4,
  #              # fontfamily=1,
  #              filename =paste(cancer,"/veen_",cancer,"_",cell,".tiff",sep=''),
  #              height = 1600,
  #              width = 1600, resolution = 500)
  
  
  venn.diagram(list(Positive=veen_positive,Negative=veen_negative,DEGs=veen_DEG),
               fill=c("#729ECE","#FF9E4A","#67BF5C"), #,"#67BF5C","#ED665D","#AD8BC9"
               alpha=c(0.5,0.5,0.5),
               col = c("#729ECE","#FF9E4A","#67BF5C"),
               cex=1,
               # cat.pos = 10,
               cat.dist = 0.01,
               cat.fontface=4,
               # fontfamily=1,
               filename =paste(cancer,"/veen_",cancer,"_",cell,".tiff",sep=''),
               height = 1600,
               width = 1600, resolution = 500)
  
  
  
  
}





# 
# 
# 
# data2_BLCA_DEG <- data1_BLCA_DEG[which(as.numeric(data1_BLCA_DEG[,2]) < 0.05),1:2]
# 
# 
# 
# 
# Gene_po_DEG <- names(which(table(c(data_positive$Gene,data2_BLCA_DEG$Gene))==2))
# Gene_ne_DEG <- names(which(table(c(data_negative$Gene,data2_BLCA_DEG$Gene))==2))
# 
# 
# veen_positive <- c(data_positive$Gene)
# veen_negative <- c(data_negative$Gene)
# veen_DEG <- c(data2_BLCA_DEG$Gene)
# 
# library(VennDiagram)
# 
# venn.diagram(list(Positive=veen_positive,Negative=veen_negative,DEGs=veen_DEG), 
#              fill=c("#729ECE","#FF9E4A","#67BF5C"), #,"#67BF5C","#ED665D","#AD8BC9"
#              alpha=c(0.5,0.5,0.5), 
#              col = c("#729ECE","#FF9E4A","#67BF5C"),
#              cex=1,
#              # cat.pos = 10,
#              cat.dist = 0.01,
#              cat.fontface=4, 
#              # fontfamily=1,
#              filename =paste("veen",cancer,"_",cell,".tiff",sep=''),
#              height = 1600, 
#              width = 1600, resolution = 500)
# 
# 
































