setwd('D:/R/Hot_tumor_and_cold_tumor/fenqi/leida/leida12_2/zhifang')


cancers <- c("BLCA","CESC","PAAD","SARC","SKCM","LGG","KIRP","THYM")


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
  dat_immu <- scale(dat)
  
  
  cluster <- read.csv(paste("D:/R/Hot_tumor_and_cold_tumor/Immunecells/CIBERSORT/OVER100/",cancer,"/",cancer,"_clusterresult/",cancer,"_cluster_hc.csv",sep=''),header = F)
  rownames(cluster) <- cluster[,1]
  
  
  all_name <- names(which(table(c(rownames(dat_immu),cluster[,1] ))==2))
  
  
  dat_cluster <- cluster[match(all_name,cluster[,1]),]
  
  dat_im <- dat_immu[match(all_name,rownames(dat_immu)),]
  
  dat_cluster <- dat_cluster[c('V2')]
  
  dat_im <- as.data.frame(dat_im)
  class(dat_im)
  
  dat <- data.frame()
  for(coln in colnames(dat_im))
  {
    for(i in all_name)
    {
      dat <- rbind(dat,c(coln,dat_cluster[match(i,rownames(dat_cluster)),],dat_im[c(coln)][match(i,rownames(dat_im)),]))
    }
  }
  
  dat[,3] = as.numeric(dat[,3])
  
  colnames(dat) <- c("cell","Group","value")
  
  dat$value <- dat$value*2+3
  
  dat$cell <- factor(dat$cell,levels = c("T_cells_CD8","NK_cells_activated","T_cells_follicular_helper","Tcells_proliferation","Cytolytic_activity","Immunogenic_cell_death","Macrophages_M2","MDSCs","Tregs","T_cell_co_inhibition","MHC_class_I","Inflammation-promoting"))
  
  library(ggplot2)
  
  pdf( paste(cancer,"_zhifang1.pdf",sep=''),height=5,width=8)

  p = ggplot(dat, aes(y=value, x=cell))+
    geom_boxplot(aes(fill=Group),outlier.shape = NA)+
    theme_bw()+
    scale_fill_manual(values = c("#88ada6","#ED665D"))
  p <- p + theme(axis.text.x = element_text(angle =60, hjust =1, vjust =1))
  p <- p + stat_compare_means(aes(group = Group, label = after_stat(p.signif)),method = "anova")

  print(p)
  dev.off()
  
  # 
  # pdf( paste(cancer,"_zhifang2.pdf",sep=''),height=5,width=8)
  # 
  #   p <- ggplot(dat,aes(x=cell, y=value, fill=Group))+geom_bar(stat='identity', width=0.5, position=position_dodge(0.6))#identity意味着把y当做值去输入，如果改成bin，就会计算y出现的频数。dodge意味是各组是左右分布而不是上下重叠
  # p <- p + scale_fill_manual(values=c(cold = "#88ada6", hot = "#ED665D"))
  # p <- p + theme(axis.text.x = element_text(angle =60, hjust =1, vjust =1))
  # p <- p + stat_compare_means(aes(group = Group, label = after_stat(p.signif)),method = "anova")
  # 
  # 
  # print(p)
  # dev.off()
  
}