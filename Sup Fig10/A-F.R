setwd('D:/R/Hot_tumor_and_cold_tumor/fenqi/leida/leida_PAAD')

meanhot <- c()
meancold <- c()

cancer <- c("PAAD")



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
dat <- dat*2+3



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


for (i in 1:length(rownames(dat_c1))) {
  
  aa<- rownames(dat_c1)
  sn <- aa[i]
  df_c1c2 <- data.frame(rbind(max_min1,dat_c1[i,]))
  rownames(df_c1c2) <- c("Max", "Min","sn")
  colnames(df_c1c2) <- c(1:12)
  
  library(fmsb)
  pdf(paste("cold",sn,".pdf",sep=''),width = 10,height = 8)
  
  radarchart(
    df_c1c2, axistype = 1,
    # Customize the polygon
    pcol = c( "#70CF07"), pfcol = scales::alpha(c("#70CF07"),0.5), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "black", 
    # Variable labels
    vlcex = 0.7, vlabels = colnames(df_c1c2),
    caxislabels = c(  1,2, 3,4,5))
  # Add an horizontal legend
  dev.off()
  
}

for (i in 1:length(rownames(dat_c2))) {
  
  aa<- rownames(dat_c2)
  sn <- aa[i]
  df_c1c2 <- data.frame(rbind(max_min1,dat_c2[i,]))
  rownames(df_c1c2) <- c("Max", "Min","sn")
  colnames(df_c1c2) <- c(1:12)
  
  library(fmsb)
  pdf(paste("hot",sn,".pdf",sep=''),width = 10,height = 8)
  
  radarchart(
    df_c1c2, axistype = 1,
    # Customize the polygon
    pcol = c( "#FC4E07"), pfcol = scales::alpha(c("#FC4E07"),0.5), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "black", 
    # Variable labels
    vlcex = 0.7, vlabels = colnames(df_c1c2),
    caxislabels = c(  1,2, 3,4,5))
  # Add an horizontal legend
  dev.off()

}


 
  
