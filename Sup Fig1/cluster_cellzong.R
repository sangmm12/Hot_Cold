rm(list=ls())

setwd("D:/R/Hot_tumor_and_cold_tumor/cluster_cell/cluster_cellzong")

updownzong <- c()

cancers <- c("BLCA","CESC","KIRP","LGG","PAAD","SARC","SKCM","THYM")
#cancers <- c("SARC","SKCM","THYM")
#cancer <- c("SARC")
for(cancer in cancers){
  
  if(dir.exists(cancer)==TRUE){
    print(cancer)
  }else {
    dir.create(cancer)
  }
  
  
  cluster <- read.csv(paste("D:/R/Hot_tumor_and_cold_tumor/Immunecells/CIBERSORT/OVER100/",cancer,"/",cancer,"_clusterresult/",cancer,"_cluster_hc.csv",sep=''),header = F)
  rownames(cluster) <- cluster[,1]
  
  
  dat <- read.csv(paste("D:/R/Hot_tumor_and_cold_tumor/Immunecells/CIBERSORT/OVER100/",cancer,"/",cancer,"_BIBERSORT.csv",sep=''))

  rownames(dat) <- dat[,1]
  dat_immu <- dat[,-1]
  
  list1 <- colnames(dat_immu)
  
  
  score <- read.csv(paste("D:/R/Hot_tumor_and_cold_tumor/cluster_immuscore/Estimate_",cancer,".csv",sep=''))
  
  rownames(score) <- score[,1]
  dat_immu2 <- score[,-1]
  
  list2 <- colnames(dat_immu2)
  
  
  library(data.table)
  other <- fread(paste("D:/R/Hot_tumor_and_cold_tumor/Immunecells/ssGSEA/ssgsea_others_Tcell_",cancer,".csv",sep=''))
  other <- as.data.frame(other)
  
  rownames(other) <- other[,1]
  dat_immu3 <- other[,-1]
  dat_immu3 <- t(dat_immu3)
  
  list3 <- colnames(dat_immu3)
  
  
  all_name <- names(which(table(c(rownames(dat_immu),rownames(dat_immu2),rownames(dat_immu3)))==3))
  
  dat_imm <- dat_immu[match(all_name,rownames(dat_immu)),]
  dat_imm2 <- dat_immu2[match(all_name,rownames(dat_immu2)),]
  dat_imm3 <- dat_immu3[match(all_name,rownames(dat_immu3)),]
  
  dat_imm <- cbind(dat_imm,dat_imm2,dat_imm3) 
  
  write.csv(dat_imm,file = paste("D:/R/Hot_tumor_and_cold_tumor/cluster_cell/cluster_cellzong/",cancer,"/data_",cancer,".csv",sep=''))
  
  
  all_name <- names(which(table(c(rownames(dat_imm),cluster[,1] ))==2))
  
  
  dat_cluster <- cluster[match(all_name,cluster[,1]),]
  

  #cell <- c("B_cells_naive","B_cells_memory","Plasma_cells","T_cells_CD8","T_cells_CD4_naive","T_cells_CD4_memory_resting","T_cells_CD4_memory_activated","T_cells_follicular_helper","T_cells_regulatory_.Tregs.","T_cells_gamma_delta","NK_cells_resting","NK_cells_activated","Monocytes","Macrophages_M0","Macrophages_M1","Macrophages_M2","Dendritic_cells_resting","Dendritic_cells_activated","Mast_cells_resting","Mast_cells_activated","Eosinophils","Neutrophils")
  #dat_im <- dat_imm[cell]
  
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
  compare_means(value ~ Group, data = dat, group.by = "Gene",method = "anova")
  
  
  pdf(paste(cancer,"/",cancer,"_cluster_ns.pdf",sep=''),width=temp_length/3.5,height = 8)
  p <- ggboxplot(dat, x = "Gene", y = "value",
                 color = "Group", palette = "FA7F6F", 
                 add = "jitter",x.text.angle=60) # palette可以按照期刊选择相应的配色，如"npg"等
  p <- p+xlab("Cell Type")+ylab("Expression Value")
  p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
  p
  print(p)
  dev.off()
  
  
  
  
  
  #####################
  dat1 <- dat[dat$Gene %in% list1,]
  
  xxxx <- matrix(unique(dat1[,1])) #colname
  
  temp_length <- length(xxxx)*2
  
  dat1$Group <- factor(dat1$Group,levels=c("hot","cold"))
  library(ggpubr)
  compare_means(value ~ Group, data = dat1, group.by = "Gene",method = "anova")
  
  
  pdf(paste(cancer,"/",cancer,"_cluster_cell_ns.pdf",sep=''),width=temp_length/3.5,height = 8)
  p <- ggboxplot(dat1, x = "Gene", y = "value",
                 color = "Group", palette = "FA7F6F", 
                 add = "jitter",x.text.angle=60) # palette可以按照期刊选择相应的配色，如"npg"等
  p <- p+xlab("Cell Type")+ylab("Expression Value")
  p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
  p
  print(p)
  dev.off()
  
  out_put_list <- c()
  mean <- c()
  
  
  cell <- c("Macrophages_M2")
  for(cell in unique(dat1$Gene)){
    
    dat1_1 <- subset(dat1,dat1$Gene==cell)
    
    dat1_1$Group <- factor(dat1_1$Group,levels=c("hot","cold"))
    
  
    
    library(ggpubr)
     compare_means(value ~ Group, data = dat1_1, group.by = "Gene",method = "anova")
    
    
    #pdf(paste(cancer,"/",cancer,"_cluster_Immuneinfiltration",cell,".pdf",sep=''),width=3,height = 8)
    p <- ggboxplot(dat1_1, x = "Group", y = "value",
                   color = "Group", palette = "FA7F6F", 
                   add = "jitter",x.text.angle=60) # palette可以按照期刊选择相应的配色，如"npg"等
    p <- p+ylab(cell)
    #p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
    # p
    # print(p)
    # dev.off()
    # print(p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova"))
    # print(p)
    # 
    
    eval(parse(text = paste(cell,'<- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")')))
    out_put_list <- append(out_put_list,cell)
    # dev.off()
    
    
    #means <- aggregate(value ~ Group, dat1_1, mean)
    
    dat1_1 <- dat1_1[,-1]
    library(dplyr)
    means=group_by(dat1_1, Group) %>% summarise_each(funs(mean))
    
    if(length(mean)==0){
      mean <- means
    }else {
      mean <- cbind(mean,means)
    }
    
    
  }
  
  pdf(paste(cancer,"/",cancer,"_cluster_cell_fen.pdf",sep=''),width = 8,height = 8*length(out_put_list)/4)
  
  eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=4))',sep='')))
  #eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=4,labels=out_put_list,label_size=40))',sep='')))
  
  dev.off()
  
  n1 <- grep("Group",colnames(mean))
  mean1 <- mean[,-n1]
  colnames(mean1) <- unique(dat1$Gene)
  rownames(mean1) <- c("hot","cold")
  
  
    

  
  #####################
  dat2 <- dat[dat$Gene %in% list2,]
  
  xxxx <- matrix(unique(dat2[,1])) #colname
  
  temp_length <- length(xxxx)*2
  
  dat2$Group <- factor(dat2$Group,levels=c("hot","cold"))
  library(ggpubr)
  compare_means(value ~ Group, data = dat2, group.by = "Gene",method = "anova")
  
  
  pdf(paste(cancer,"/",cancer,"_cluster_Immuneinscore_ns.pdf",sep=''),width=temp_length/3.5,height = 8)
  p <- ggboxplot(dat2, x = "Gene", y = "value",
                 color = "Group", palette = "FA7F6F", 
                 add = "jitter",x.text.angle=60) # palette可以按照期刊选择相应的配色，如"npg"等
  p <- p+xlab("Cell Type")+ylab("Expression Value")
  p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
  p
  print(p)
  dev.off()
  
  out_put_list <- c()
  mean <- c()
  
  
  cell <- c("Macrophages_M2")
  for(cell in unique(dat2$Gene)){
    
    dat2_1 <- subset(dat2,dat2$Gene==cell)
    
    dat2_1$Group <- factor(dat2_1$Group,levels=c("hot","cold"))
    
    
    
    library(ggpubr)
    compare_means(value ~ Group, data = dat2_1, group.by = "Gene",method = "anova")
    
    
    #pdf(paste(cancer,"/",cancer,"_cluster_Immuneinfiltration",cell,".pdf",sep=''),width=3,height = 8)
    p <- ggboxplot(dat2_1, x = "Group", y = "value",
                   color = "Group", palette = "FA7F6F", 
                   add = "jitter",x.text.angle=60) # palette可以按照期刊选择相应的配色，如"npg"等
    p <- p+ylab(cell)
    #p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
    # p
    # print(p)
    # dev.off()
    # print(p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova"))
    # print(p)
    # 
    
    eval(parse(text = paste(cell,'<- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")')))
    out_put_list <- append(out_put_list,cell)
    # dev.off()
    
    
    #means <- aggregate(value ~ Group, dat1_1, mean)
    
    dat2_1 <- dat2_1[,-1]
    library(dplyr)
    means=group_by(dat2_1, Group) %>% summarise_each(funs(mean))
    
    if(length(mean)==0){
      mean <- means
    }else {
      mean <- cbind(mean,means)
    }
    
    
  }
  
  pdf(paste(cancer,"/",cancer,"_cluster_Immuneinscore_fen.pdf",sep=''),width = 8,height = 8*length(out_put_list)/4)
  
  eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=4))',sep='')))
  #eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=4,labels=out_put_list,label_size=40))',sep='')))
  
  dev.off()
  
  n1 <- grep("Group",colnames(mean))
  mean2 <- mean[,-n1]
  colnames(mean2) <- unique(dat2$Gene)
  rownames(mean2) <- c("hot","cold")
  
 
  
  
  
   
  
  #####################
  dat3 <- dat[dat$Gene %in% list3,]
  
  
  xxxx <- matrix(unique(dat3[,1])) #colname
  
  temp_length <- length(xxxx)*2
  
  dat3$Group <- factor(dat3$Group,levels=c("hot","cold"))
  library(ggpubr)
  compare_means(value ~ Group, data = dat3, group.by = "Gene",method = "anova")
  
  
  pdf(paste(cancer,"/",cancer,"_cluster_activities_ns.pdf",sep=''),width=temp_length/3.5,height = 8)
  p <- ggboxplot(dat3, x = "Gene", y = "value",
                 color = "Group", palette = "FA7F6F", 
                 add = "jitter",x.text.angle=60) # palette可以按照期刊选择相应的配色，如"npg"等
  p <- p+xlab("Cell Type")+ylab("Expression Value")
  p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
  p
  print(p)
  dev.off()
  
  out_put_list <- c()
  mean <- c()
  
  dat3[dat3 == "Check-point"] = "Check_point"
  dat3[dat3 == "Inflammation-promoting"] = "Inflammation_promoting"
  dat3[dat3 == "T_cell_co-inhibition"] = "T_cell_co_inhibition"
  dat3[dat3 == "T_cell_co-stimulation"] = "T_cell_co_stimulation"
  
  
  
  
  cell <- c("Macrophages_M2")
  for(cell in unique(dat3$Gene)){
    
    dat3_1 <- subset(dat3,dat3$Gene==cell)
    
    dat3_1$Group <- factor(dat3_1$Group,levels=c("hot","cold"))
    
    
    
    library(ggpubr)
    compare_means(value ~ Group, data = dat3_1, group.by = "Gene",method = "anova")
    
    
    #pdf(paste(cancer,"/",cancer,"_cluster_activities_",cell,".pdf",sep=''),width=3,height = 8)
    p <- ggboxplot(dat3_1, x = "Group", y = "value",
                   color = "Group", palette = "FA7F6F", 
                   add = "jitter",x.text.angle=60) # palette可以按照期刊选择相应的配色，如"npg"等
    p <- p+ylab(cell)
    #p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
    # p
    # print(p)
    # dev.off()
    # print(p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova"))
    # print(p)
    # 
    
    eval(parse(text = paste(paste(cell),'<- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")')))
    out_put_list <- append(out_put_list,cell)
    # dev.off()
    
    
    #means <- aggregate(value ~ Group, dat1_1, mean)
    
    dat3_1 <- dat3_1[,-1]
    library(dplyr)
    means=group_by(dat3_1, Group) %>% summarise_each(funs(mean))
    
    if(length(mean)==0){
      mean <- means
    }else {
      mean <- cbind(mean,means)
    }
    
    
  }
  
  pdf(paste(cancer,"/",cancer,"_cluster_activities_fen.pdf",sep=''),width = 8,height = 8*length(out_put_list)/4)
  
  eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=4))',sep='')))
  #eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=4,labels=out_put_list,label_size=40))',sep='')))
  
  dev.off()
  
  n1 <- grep("Group",colnames(mean))
  mean3 <- mean[,-n1]
  colnames(mean3) <- unique(dat3$Gene)
  rownames(mean3) <- c("hot","cold")
  
  meanzong <- cbind(mean1,mean2,mean3)
  meanzong <-  t(meanzong) 
  
  updown <- rep(0, times=length(rownames(meanzong)))
  
  for( i in 1:length(rownames(meanzong))){
    if( meanzong[i,1] > meanzong[i,2] ){
      updown[i]= 1
    }else{
      updown[i]= -1
    }
  }
  
  if(length(updownzong)==0){
    updownzong <- updown
  }else {
    updownzong <- data.frame(updownzong,updown)
  }
 
}

colnames(updownzong) <- cancers
rownames(updownzong) <- rownames(meanzong)

write.csv(updownzong,file = "updown.csv")
