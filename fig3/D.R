rm(list=ls())

setwd("D:/R/Hot_tumor_and_cold_tumor/cor/cell_score")

cancers <- c("BLCA","CESC","PAAD","SARC","SKCM")
#
#cancer <- c("BLCA")
for(cancer in cancers){
  
  
  library(data.table)
  
  dat <-  fread(paste("D:/R/Hot_tumor_and_cold_tumor/cluster_cell/cluster_cellzong/",cancer,"/data_",cancer,".csv",sep=''),header = T)
  
  dat <- as.data.frame(dat)
  
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]
  
  celllist <- c("T_cells_CD8","NK_cells_activated","T_cells_follicular_helper","Macrophages_M2","Inflammation-promoting","Cytolytic_activity","T_cell_co-inhibition","Tcells_proliferation","MDSCs")
  
  data <- dat[,match(celllist,colnames(dat))]
  
  colnames(data) <- c("T_cells_CD8","NK_cells_activated","T_cells_follicular_helper","Macrophages_M2","Inflammation_promoting","Cytolytic_activity","T_cell_co_inhibition","Tcells_proliferation","MDSCs")
  
  cells <- c("T_cells_CD8","NK_cells_activated","T_cells_follicular_helper","Macrophages_M2")
  cell1s <- c("T","NK","Tfh","M2")
  scores <- c("Inflammation_promoting","Cytolytic_activity","T_cell_co_inhibition","Tcells_proliferation","MDSCs")
  
  for(i in 1:length(cells)){
    for(j in 1:length(scores)){
      cell <- cells[i]
      cell1 <- cell1s[i]
      score <- scores[j]
      
      dat_1 <- data[,grep(cell,colnames(data))]
      dat_2 <- data[,grep(score,colnames(data))]
      
      
      plot(dat_1,dat_2)
      
      library(ggpubr)
      library(reshape2)
      library(ggsci)
      d1 <- data.frame(NT5E=dat_1 ,VEGFA=dat_2)
      
      library(Hmisc)
      
      #ggscatterstats(data=d1,x=log(NT5E+1),y=log(VEGFA+1),bf.message = F)
      
      
      pdf( paste("cor_",cancer,"_",cell1,"_",score,".pdf",sep=''),height=5,width=5)
      # p <- ggplot(d1, aes( x = NT5E, y =VEGFA)) + geom_point( size=3) + geom_smooth(method = 'lm', formula = y ~ x, se = T)+
      #   theme(axis.text = element_text(size = 25),axis.title = element_text(size = 25))+
      #   theme_classic() +labs(x="log(NT5E+1)",y="log(CD276+1)") + stat_cor(method = "pearson",size=8,,label.x =min(dat_1), label.y =max(dat_2))
      
      p <- ggplot(d1,aes(x=NT5E,y=VEGFA)) + geom_point(shape=19) +
        xlab( paste(cell,sep='')) + ylab(paste(score,sep=''))+geom_smooth(method = lm)+theme_classic()+
        theme(axis.text = element_text(size = 15),axis.title = element_text(size = 25))+
        stat_cor(method = "pearson",size=8,label.x =min(dat_1), label.y =max(dat_2))
      
      
      print(p)
      dev.off()
      
      
    }
  }
}
