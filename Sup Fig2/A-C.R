setwd('D:/R/Hot_tumor_and_cold_tumor/fenqi/leida/leida12_2')


cancers <- c("BLCA","CESC","PAAD","SARC","SKCM","LGG","KIRP","THYM")

meanhot <- c()
meancold <- c()

#cancer <- c("BLCA")
for(cancer in cancers){
  
  library(data.table)
  
  dat <-  fread(paste("D:/R/Hot_tumor_and_cold_tumor/cluster_cell/cluster_cellzong/",cancer,"/data_",cancer,".csv",sep=''),header = T)
  
  dat <- as.data.frame(dat)
  
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]
  
  celllist <- c("T_cells_CD8","NK_cells_activated","T_cells_follicular_helper","Tcells_proliferation","Cytolytic_activity","Immunogenic_cell_death","Macrophages_M2","MDSCs","T_cells_regulatory_.Tregs.","T_cell_co-inhibition","MHC_class_I","Inflammation-promoting")
  
  dat <- dat[,match(celllist,colnames(dat))]
  
  colnames(dat) <- c("T_cells_CD8","NK_cells_activated","T_cells_follicular_helper","Tcells_proliferation","Cytolytic_activity","Immunogenic_cell_death","Macrophages_M2","MDSCs","Tregs","T_cell_co_inhibition","MHC_class_I","Inflammation-promoting")
  
  
  #aa <- scale(dat$T_cells_CD8)
  dat <- scale(dat)
  
  
  
  max_min1 <- data.frame(
    T_cells_CD8 = c(5,1), NK_cells_activated = c(5,1), T_cells_follicular_helper = c(5,1),
    Tcells_proliferation = c(5,1),Cytolytic_activity = c(5,1),Immunogenic_cell_death = c(5,1),
    Macrophages_M2 = c(5,1),MDSCs = c(5,1), Tregs = c(5,1),
    T_cell_co_inhibition= c(5,1),
    MHC_class_I = c(5,1), Inflammation_promoting = c(5,1)
  )
  
  rownames(max_min1) <- c("Max", "Min")
  
  
  cluster <- read.csv(paste("D:/R/Hot_tumor_and_cold_tumor/Immunecells/CIBERSORT/OVER100/",cancer,"/",cancer,"_clusterresult/",cancer,"_cluster_hc.csv",sep=''),header = F)
  
  
  c1 <- cluster[which(cluster$V2=="cold"),]
  c2 <- cluster[which(cluster$V2=="hot"),]
  
  
  
  dat_c1 <- as.matrix(dat)[na.omit(match(c1[,1],rownames(dat))),]
  dat_c2 <- as.matrix(dat)[na.omit(match(c2[,1],rownames(dat))),]
  
  result1 = data.frame(apply(dat_c1,2,mean))
  mean1 <- t(result1)
  mean1 <- (as.numeric(mean1))*2+3
  
  if(length(meancold)==0){
    meancold <-  mean1
  }else {
    meancold <-  rbind(meancold,mean1)
  }
  
  result2 = data.frame(apply(dat_c2,2,mean))
  mean2 <- t(result2)
  mean2 <- (as.numeric(mean2))*2+3
  
  if(length(meanhot)==0){
    meanhot <-  mean2
  }else {
    meanhot <-  rbind(meanhot,mean2)
  }
  
  
  
  df_c1c2 <- data.frame(rbind(max_min1, mean1,mean2))
  
  rownames(df_c1c2) <- c("Max", "Min","cold","hot")
  
  library(fmsb)
  pdf(paste(cancer,"_leida_hc15.pdf",sep=''),width = 10,height = 8)
  
  radarchart(
    df_c1c2, axistype = 1,
    # Customize the polygon
    pcol = c( "#70CF07", "#FC4E07"), pfcol = scales::alpha(c("#70CF07", "#FC4E07"),0.5), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "black", 
    # Variable labels
    vlcex = 0.7, vlabels = colnames(df_c1c2),
    caxislabels = c(  1,2, 3,4,5))
  # Add an horizontal legend
  legend(
    x = -1, legend = rownames(df_c1c2[-c(1,2),]), horiz = TRUE,
    bty = "n", pch = 20 , col = c( "#70CF07", "#FC4E07"),
    text.col = "black", cex = 1, pt.cex = 1.5
  )
  
  dev.off()
  
  
}





rownames(meanhot) <- cancers
rownames(meancold) <- cancers

# meancold <- meancold+1
# meanhot <- meanhot+1

write.csv(meancold,file=paste("meancold.csv",sep=''),row.names = T)
write.csv(meanhot,file=paste("meanhot.csv",sep=''),row.names = T)








for(cancer in cancers){
  
  library(data.table)
  
  dat <-  fread(paste("D:/R/Hot_tumor_and_cold_tumor/cluster_cell/cluster_cellzong/",cancer,"/data_",cancer,".csv",sep=''),header = T)
  
  dat <- as.data.frame(dat)
  
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]
  
  celllist <- c("T_cells_CD8","NK_cells_activated","T_cells_follicular_helper","Tcells_proliferation","Cytolytic_activity","Immunogenic_cell_death","Macrophages_M2","MDSCs","T_cells_regulatory_.Tregs.","T_cell_co-inhibition","MHC_class_I","Inflammation-promoting")
  
  dat <- dat[,match(celllist,colnames(dat))]
  
  colnames(dat) <- c("T_cells_CD8","NK_cells_activated","T_cells_follicular_helper","Tcells_proliferation","Cytolytic_activity","Immunogenic_cell_death","Macrophages_M2","MDSCs","Tregs","T_cell_co_inhibition","MHC_class_I","Inflammation-promoting")
  
  
  #aa <- scale(dat$T_cells_CD8)
  dat <- scale(dat)
  
  
  
  max_min1 <- data.frame(
    T_cells_CD8 = c(5,1), NK_cells_activated = c(5,1), T_cells_follicular_helper = c(5,1),
    Tcells_proliferation = c(5,1),Cytolytic_activity = c(5,1),Immunogenic_cell_death = c(5,1),
    Macrophages_M2 = c(5,1),MDSCs = c(5,1), Tregs = c(5,1),
    T_cell_co_inhibition= c(5,1),
    MHC_class_I = c(5,1), Inflammation_promoting = c(5,1)
  )
  
  rownames(max_min1) <- c("Max", "Min")
  
  
  cluster <- read.csv(paste("D:/R/Hot_tumor_and_cold_tumor/Immunecells/CIBERSORT/OVER100/",cancer,"/",cancer,"_clusterresult/",cancer,"_cluster_hc.csv",sep=''),header = F)
  
  
  c1 <- cluster[which(cluster$V2=="cold"),]
  c2 <- cluster[which(cluster$V2=="hot"),]
  
  
  
  dat_c1 <- as.matrix(dat)[na.omit(match(c1[,1],rownames(dat))),]
  dat_c2 <- as.matrix(dat)[na.omit(match(c2[,1],rownames(dat))),]
  
  result1 = data.frame(apply(dat_c1,2,mean))
  mean1 <- t(result1)
  #mean1 <- (as.numeric(mean1))*3.6+2.5
  mean1 <- (as.numeric(mean1))*2+3
  
  result2 = data.frame(apply(dat_c2,2,mean))
  mean2 <- t(result2)
  #mean2 <- (as.numeric(mean2))*3.6+2.5
  mean2 <- (as.numeric(mean2))*2+3
  
  
  df_c1c2 <- data.frame(rbind(max_min1, mean1,mean2))
  
  rownames(df_c1c2) <- c("Max", "Min","cold","hot")
  colnames(df_c1c2) <- c(1:12)
  
  library(fmsb)
  pdf(paste(cancer,"_leida.pdf",sep=''),width = 10,height = 8)
  
  radarchart(
    df_c1c2, axistype = 1,
    # Customize the polygon
    pcol = c( "#70CF07", "#FC4E07"), pfcol = scales::alpha(c("#70CF07", "#FC4E07"),0.5), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "black", 
    # Variable labels
    vlcex = 0.7, vlabels = colnames(df_c1c2),
    caxislabels = c( 1,2, 3,4,5))
  # Add an horizontal legend
  legend(
    x = -1, legend = rownames(df_c1c2[-c(1,2),]), horiz = TRUE,
    bty = "n", pch = 20 , col = c( "#70CF07", "#FC4E07"),
    text.col = "black", cex = 1, pt.cex = 1.5
  )
  
  dev.off()
  
  
}

