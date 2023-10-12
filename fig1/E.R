rm(list=ls())

setwd("D:/R/Hot_tumor_and_cold_tumor/Immunecells/immunecells")


cancers <- c("BLCA","CESC","KIRP","LGG","PAAD","SARC","SKCM","THYM")
#
#cancer <- c("BLCA")
for(cancer in cancers){
  
  
  if(dir.exists(cancer)==TRUE){
    print(cancer)
  }else {
    dir.create(cancer)
  }
  
  
  cluster <- read.csv(paste("D:/R/Hot_tumor_and_cold_tumor/Immunecells/CIBERSORT/OVER100/",cancer,"/",cancer,"_clusterresult/",cancer,"_cluster_hc.csv",sep=''),header = F)
  rownames(cluster) <- cluster[,1]
  
  library(data.table)
  BIBER <- fread(paste("D:/R/Hot_tumor_and_cold_tumor/Immunecells/8cancers/",cancer,"_BIBERSORT.csv",sep=''))
  BIBER <- as.data.frame(BIBER)
  
  BIBER1 <- data.frame(BIBER$V1,BIBER$T_cells_CD8,BIBER$NK_cells_resting,BIBER$NK_cells_activated,BIBER$Macrophages_M0,BIBER$Macrophages_M1,BIBER$Macrophages_M2)
  rownames(BIBER1) <- BIBER1[,1]
  BIBER1 <- BIBER1[,-1]
  colnames(BIBER1) <- c("T_cells_CD8","NK_cells_resting","NK_cells_activated","Macrophages_M0","Macrophages_M1","Macrophages_M2")
  
  all_name <- names(which(table(c(rownames(BIBER1),rownames(cluster)))==2))
  
  BIBER2 <- BIBER1[match(all_name,rownames(BIBER1)),]
  dat_cluster <- cluster[match(all_name,rownames(cluster)),]
  
  dat_cluster <- dat_cluster[c('V2')]
  
  dat <- data.frame()
  for(coln in colnames(BIBER2))
  {
    for(i in all_name)
    {
      dat <- rbind(dat,c(coln,dat_cluster[match(i,rownames(dat_cluster)),],BIBER2[c(coln)][match(i,rownames(BIBER2)),]))
    }
  }
  
  dat[,3] = as.numeric(dat[,3])
  
  colnames(dat) <- c("Gene","Group","value")
  
  dat$Group <- factor(dat$Group,levels=c("hot","cold"))
  
  xxxx <- matrix(unique(dat[,1])) #colname
  
  temp_length <- length(xxxx)*2
  
  
  library(ggpubr)
  xx <- compare_means(value ~ Group, data = dat, group.by = "Gene",method = "anova")
  
  BIBER_p <- data.frame(xx$Gene,xx$p)
  colnames(BIBER_p) <- c("CELL","BIBER_p")
  
  pdf(paste(cancer,"/",cancer,"_CIBER_ns.pdf",sep=''),width=temp_length/3.5,height = 8)
  p <- ggboxplot(dat, x = "Gene", y = "value",
                 color = "Group", palette = "FA7F6F", 
                 add = "jitter",x.text.angle=60) # palette可以按照期刊选择相应的配色，如"npg"等
  p <- p+xlab("Cell Type")+ylab("Expression Value")
  p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
  p
  print(p)
  dev.off()
  
  out_put_list <- c()
  mean <- c()
  
  
  for(cell in unique(dat$Gene)){
    
    n1 <- grep(cell,dat$Gene)
    
    dat_1 <- dat[n1,]
    
    dat_1$Group <- factor(dat_1$Group,levels=c("hot","cold"))
    
    
    
    library(ggpubr)
    compare_means(value ~ Group, data = dat_1, group.by = "Gene",method = "anova")
    
    
    #pdf(paste(cancer,"/",cancer,"_cluster_activities_",cell,".pdf",sep=''),width=3,height = 8)
    p <- ggboxplot(dat_1, x = "Group", y = "value",
                   color = "Group", palette = "FA7F6F", 
                   add = "jitter") # palette可以按照期刊选择相应的配色，如"npg"等
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
    
    dat_1 <- dat_1[,-1]
    library(dplyr)
    means=group_by(dat_1, Group) %>% summarise_each(funs(mean))
    
    if(length(mean)==0){
      mean <- means
    }else {
      mean <- cbind(mean,means)
    }
  }
  pdf(paste(cancer,"/",cancer,"_CIBER_fen.pdf",sep=''),width = 8,height = 8*length(out_put_list)/4)
  
  eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=4))',sep='')))
  #eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=4,labels=out_put_list,label_size=40))',sep='')))
  
  dev.off()
  
  n1 <- grep("Group",colnames(mean))
  mean <- mean[,-n1]
  colnames(mean) <- unique(dat$Gene)
  rownames(mean) <- c("hot","cold")
  
  mean <-  t(mean) 
  
  updown <- rep(0, times=length(rownames(mean)))
  
  for( i in 1:length(rownames(mean))){
    if( mean[i,1] > mean[i,2] ){
      updown[i]= 1
    }else{
      updown[i]= -1
    }
  }
  
  BIBER_updown <- data.frame(c("T_cells_CD8","NK_cells_resting","NK_cells_activated","Macrophages_M0","Macrophages_M1","Macrophages_M2"),updown)
  colnames(BIBER_updown) <- c("CELL","CIBERSORT_ud")
  
  
  
  
  
  
  
  ##############################
  
  
  
  TIMER <- fread(paste("D:/R/Hot_tumor_and_cold_tumor/Immunecells/8cancers/",cancer,"_TIMER.csv",sep=''))
  TIMER <- as.data.frame(TIMER)
  
  TIMER1 <- data.frame(TIMER$V1,TIMER$T_cell_CD8,TIMER$Macrophage)
  rownames(TIMER1) <- TIMER1[,1]
  TIMER1 <- TIMER1[,-1]
  colnames(TIMER1) <- c("T_cells_CD8","Macrophages")
  
  all_name <- names(which(table(c(rownames(TIMER1),rownames(cluster)))==2))
  
  TIMER2 <- TIMER1[match(all_name,rownames(TIMER1)),]
  dat_cluster <- cluster[match(all_name,rownames(cluster)),]
  
  dat_cluster <- dat_cluster[c('V2')]
  
  dat <- data.frame()
  for(coln in colnames(TIMER2))
  {
    for(i in all_name)
    {
      dat <- rbind(dat,c(coln,dat_cluster[match(i,rownames(dat_cluster)),],TIMER2[c(coln)][match(i,rownames(TIMER2)),]))
    }
  }
  
  dat[,3] = as.numeric(dat[,3])
  
  colnames(dat) <- c("Gene","Group","value")
  
  dat$Group <- factor(dat$Group,levels=c("hot","cold"))
  
  xxxx <- matrix(unique(dat[,1])) #colname
  
  temp_length <- length(xxxx)*2
  
  
  library(ggpubr)
  xx <- compare_means(value ~ Group, data = dat, group.by = "Gene",method = "anova")
  
  TIMER_p <- data.frame(xx$Gene,xx$p)
  colnames(TIMER_p) <- c("CELL","TIMER_p")
  
  pdf(paste(cancer,"/",cancer,"_TIMER_ns.pdf",sep=''),width=temp_length/3.5,height = 8)
  p <- ggboxplot(dat, x = "Gene", y = "value",
                 color = "Group", palette = "FA7F6F", 
                 add = "jitter",x.text.angle=60) # palette可以按照期刊选择相应的配色，如"npg"等
  p <- p+xlab("Cell Type")+ylab("Expression Value")
  p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
  p
  print(p)
  dev.off()
  
  out_put_list <- c()
  mean <- c()
  
  
  for(cell in unique(dat$Gene)){
    
    n1 <- grep(cell,dat$Gene)
    
    dat_1 <- dat[n1,]
    
    dat_1$Group <- factor(dat_1$Group,levels=c("hot","cold"))
    
    
    
    library(ggpubr)
    compare_means(value ~ Group, data = dat_1, group.by = "Gene",method = "anova")
    
    
    #pdf(paste(cancer,"/",cancer,"_cluster_activities_",cell,".pdf",sep=''),width=3,height = 8)
    p <- ggboxplot(dat_1, x = "Group", y = "value",
                   color = "Group", palette = "FA7F6F", 
                   add = "jitter") # palette可以按照期刊选择相应的配色，如"npg"等
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
    
    dat_1 <- dat_1[,-1]
    library(dplyr)
    means=group_by(dat_1, Group) %>% summarise_each(funs(mean))
    
    if(length(mean)==0){
      mean <- means
    }else {
      mean <- cbind(mean,means)
    }
  }
  pdf(paste(cancer,"/",cancer,"_TIMER_fen.pdf",sep=''),width = 8,height = 8*length(out_put_list)/4)
  
  eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=4))',sep='')))
  #eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=4,labels=out_put_list,label_size=40))',sep='')))
  
  dev.off()
  
  n1 <- grep("Group",colnames(mean))
  mean <- mean[,-n1]
  colnames(mean) <- unique(dat$Gene)
  rownames(mean) <- c("hot","cold")
  
  mean <-  t(mean) 
  
  updown <- rep(0, times=length(rownames(mean)))
  
  for( i in 1:length(rownames(mean))){
    if( mean[i,1] > mean[i,2] ){
      updown[i]= 1
    }else{
      updown[i]= -1
    }
  }
  
  TIMER_updown <- data.frame(c("T_cells_CD8","Macrophages"),updown)
  colnames(TIMER_updown) <- c("CELL","TIMER_ud")
  
  
  
  
  
  
  
  
  #############################################
  
  
  EPIC <- fread(paste("D:/R/Hot_tumor_and_cold_tumor/Immunecells/8cancers/",cancer,"_EPIC.csv",sep=''))
  EPIC <- as.data.frame(EPIC)
  
  EPIC1 <- data.frame(EPIC$V1,EPIC$CD8_Tcells,EPIC$NKcells,EPIC$Macrophages)
  rownames(EPIC1) <- EPIC1[,1]
  EPIC1 <- EPIC1[,-1]
  colnames(EPIC1) <- c("T_cells_CD8","NK_cells","Macrophages")
  
  
  all_name <- names(which(table(c(rownames(EPIC1),rownames(cluster)))==2))
  
  EPIC2 <- EPIC1[match(all_name,rownames(EPIC1)),]
  dat_cluster <- cluster[match(all_name,rownames(cluster)),]
  
  dat_cluster <- dat_cluster[c('V2')]
  
  dat <- data.frame()
  for(coln in colnames(EPIC2))
  {
    for(i in all_name)
    {
      dat <- rbind(dat,c(coln,dat_cluster[match(i,rownames(dat_cluster)),],EPIC2[c(coln)][match(i,rownames(EPIC2)),]))
    }
  }
  
  dat[,3] = as.numeric(dat[,3])
  
  colnames(dat) <- c("Gene","Group","value")
  
  dat$Group <- factor(dat$Group,levels=c("hot","cold"))
  
  xxxx <- matrix(unique(dat[,1])) #colname
  
  temp_length <- length(xxxx)*2
  
  
  library(ggpubr)
  xx <- compare_means(value ~ Group, data = dat, group.by = "Gene",method = "anova")
  
  EPIC_p <- data.frame(xx$Gene,xx$p)
  colnames(EPIC_p) <- c("CELL","EPIC_p")
  
  pdf(paste(cancer,"/",cancer,"_EPIC_ns.pdf",sep=''),width=temp_length/3.5,height = 8)
  p <- ggboxplot(dat, x = "Gene", y = "value",
                 color = "Group", palette = "FA7F6F", 
                 add = "jitter",x.text.angle=60) # palette可以按照期刊选择相应的配色，如"npg"等
  p <- p+xlab("Cell Type")+ylab("Expression Value")
  p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
  p
  print(p)
  dev.off()
  
  out_put_list <- c()
  mean <- c()
  
  
  for(cell in unique(dat$Gene)){
    
    n1 <- grep(cell,dat$Gene)
    
    dat_1 <- dat[n1,]
    
    dat_1$Group <- factor(dat_1$Group,levels=c("hot","cold"))
    
    
    
    library(ggpubr)
    compare_means(value ~ Group, data = dat_1, group.by = "Gene",method = "anova")
    
    
    #pdf(paste(cancer,"/",cancer,"_cluster_activities_",cell,".pdf",sep=''),width=3,height = 8)
    p <- ggboxplot(dat_1, x = "Group", y = "value",
                   color = "Group", palette = "FA7F6F", 
                   add = "jitter") # palette可以按照期刊选择相应的配色，如"npg"等
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
    
    dat_1 <- dat_1[,-1]
    library(dplyr)
    means=group_by(dat_1, Group) %>% summarise_each(funs(mean))
    
    if(length(mean)==0){
      mean <- means
    }else {
      mean <- cbind(mean,means)
    }
  }
  pdf(paste(cancer,"/",cancer,"_EPIC_fen.pdf",sep=''),width = 8,height = 8*length(out_put_list)/4)
  
  eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=4))',sep='')))
  #eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=4,labels=out_put_list,label_size=40))',sep='')))
  
  dev.off()
  
  n1 <- grep("Group",colnames(mean))
  mean <- mean[,-n1]
  colnames(mean) <- unique(dat$Gene)
  rownames(mean) <- c("hot","cold")
  
  mean <-  t(mean) 
  
  updown <- rep(0, times=length(rownames(mean)))
  
  for( i in 1:length(rownames(mean))){
    if( mean[i,1] > mean[i,2] ){
      updown[i]= 1
    }else{
      updown[i]= -1
    }
  }
  
  EPIC_updown <- data.frame(c("T_cells_CD8","NK_cells","Macrophages"),updown)
  colnames(EPIC_updown) <- c("CELL","EPIC_ud")
  
  
  
  
  
  
  
  
  
  
  
  ###############################################
  
  
  
  MCP <- fread(paste("D:/R/Hot_tumor_and_cold_tumor/Immunecells/8cancers/",cancer,"_MCPcounter.csv",sep=''))
  MCP <- as.data.frame(MCP)
  
  MCP1 <- data.frame(MCP$V1,MCP$CD8_T_cells,MCP$NK_cells)
  rownames(MCP1) <- MCP1[,1]
  MCP1 <- MCP1[,-1]
  colnames(MCP1) <- c("T_cells_CD8","NK_cells")
  
  all_name <- names(which(table(c(rownames(MCP1),rownames(cluster)))==2))
  
  MCP2 <- MCP1[match(all_name,rownames(MCP1)),]
  dat_cluster <- cluster[match(all_name,rownames(cluster)),]
  
  dat_cluster <- dat_cluster[c('V2')]
  
  dat <- data.frame()
  for(coln in colnames(MCP2))
  {
    for(i in all_name)
    {
      dat <- rbind(dat,c(coln,dat_cluster[match(i,rownames(dat_cluster)),],MCP2[c(coln)][match(i,rownames(MCP2)),]))
    }
  }
  
  dat[,3] = as.numeric(dat[,3])
  
  colnames(dat) <- c("Gene","Group","value")
  
  dat$Group <- factor(dat$Group,levels=c("hot","cold"))
  
  xxxx <- matrix(unique(dat[,1])) #colname
  
  temp_length <- length(xxxx)*2
  
  
  library(ggpubr)
  xx <- compare_means(value ~ Group, data = dat, group.by = "Gene",method = "anova")
  
  MCP_p <- data.frame(xx$Gene,xx$p)
  colnames(MCP_p) <- c("CELL","MCPcounter_p")
  
  pdf(paste(cancer,"/",cancer,"_MCP_ns.pdf",sep=''),width=temp_length/3.5,height = 8)
  p <- ggboxplot(dat, x = "Gene", y = "value",
                 color = "Group", palette = "FA7F6F", 
                 add = "jitter",x.text.angle=60) # palette可以按照期刊选择相应的配色，如"npg"等
  p <- p+xlab("Cell Type")+ylab("Expression Value")
  p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
  p
  print(p)
  dev.off()
  
  out_put_list <- c()
  mean <- c()
  
  
  for(cell in unique(dat$Gene)){
    
    n1 <- grep(cell,dat$Gene)
    
    dat_1 <- dat[n1,]
    
    dat_1$Group <- factor(dat_1$Group,levels=c("hot","cold"))
    
    
    
    library(ggpubr)
    compare_means(value ~ Group, data = dat_1, group.by = "Gene",method = "anova")
    
    
    #pdf(paste(cancer,"/",cancer,"_cluster_activities_",cell,".pdf",sep=''),width=3,height = 8)
    p <- ggboxplot(dat_1, x = "Group", y = "value",
                   color = "Group", palette = "FA7F6F", 
                   add = "jitter") # palette可以按照期刊选择相应的配色，如"npg"等
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
    
    dat_1 <- dat_1[,-1]
    library(dplyr)
    means=group_by(dat_1, Group) %>% summarise_each(funs(mean))
    
    if(length(mean)==0){
      mean <- means
    }else {
      mean <- cbind(mean,means)
    }
  }
  pdf(paste(cancer,"/",cancer,"_MCP_fen.pdf",sep=''),width = 8,height = 8*length(out_put_list)/4)
  
  eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=4))',sep='')))
  #eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=4,labels=out_put_list,label_size=40))',sep='')))
  
  dev.off()
  
  n1 <- grep("Group",colnames(mean))
  mean <- mean[,-n1]
  colnames(mean) <- unique(dat$Gene)
  rownames(mean) <- c("hot","cold")
  
  mean <-  t(mean) 
  
  updown <- rep(0, times=length(rownames(mean)))
  
  for( i in 1:length(rownames(mean))){
    if( mean[i,1] > mean[i,2] ){
      updown[i]= 1
    }else{
      updown[i]= -1
    }
  }
  
  MCP_updown <- data.frame(c("T_cells_CD8","NK_cells"),updown)
  colnames(MCP_updown) <- c("CELL","MCPcounter_ud")
  
  
  
  
  
  
  
  
  
  
  
  
  ###################################################
  
  
  
  xCELL <- fread(paste("D:/R/Hot_tumor_and_cold_tumor/Immunecells/8cancers/",cancer,"_xCELL.csv",sep=''))
  xCELL <- as.data.frame(xCELL)
  
  xCELL1 <- data.frame(xCELL$V1,xCELL[,13],xCELL$NK_cells,xCELL$Macrophages,xCELL$Macrophages_M1,xCELL$Macrophages_M2)
  rownames(xCELL1) <- xCELL1[,1]
  xCELL1 <- xCELL1[,-1]
  colnames(xCELL1) <- c("T_cells_CD8","NK_cells","Macrophages","Macrophages_M1","Macrophages_M2")
  
  all_name <- names(which(table(c(rownames(xCELL1),rownames(cluster)))==2))
  
  xCELL2 <- xCELL1[match(all_name,rownames(xCELL1)),]
  dat_cluster <- cluster[match(all_name,rownames(cluster)),]
  
  dat_cluster <- dat_cluster[c('V2')]
  
  dat <- data.frame()
  for(coln in colnames(xCELL2))
  {
    for(i in all_name)
    {
      dat <- rbind(dat,c(coln,dat_cluster[match(i,rownames(dat_cluster)),],xCELL2[c(coln)][match(i,rownames(xCELL2)),]))
    }
  }
  
  dat[,3] = as.numeric(dat[,3])
  
  colnames(dat) <- c("Gene","Group","value")
  
  dat$Group <- factor(dat$Group,levels=c("hot","cold"))
  
  xxxx <- matrix(unique(dat[,1])) #colname
  
  temp_length <- length(xxxx)*2
  
  
  library(ggpubr)
  xx <- compare_means(value ~ Group, data = dat, group.by = "Gene",method = "anova")
  
  xCELL_p <- data.frame(xx$Gene,xx$p)
  colnames(xCELL_p) <- c("CELL","xCELL_p")
  
  pdf(paste(cancer,"/",cancer,"_xCELL_ns.pdf",sep=''),width=temp_length/3.5,height = 8)
  p <- ggboxplot(dat, x = "Gene", y = "value",
                 color = "Group", palette = "FA7F6F", 
                 add = "jitter",x.text.angle=60) # palette可以按照期刊选择相应的配色，如"npg"等
  p <- p+xlab("Cell Type")+ylab("Expression Value")
  p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
  p
  print(p)
  dev.off()
  
  out_put_list <- c()
  mean <- c()
  
  
  for(cell in unique(dat$Gene)){
    
    n1 <- grep(cell,dat$Gene)
    
    dat_1 <- dat[n1,]
    
    dat_1$Group <- factor(dat_1$Group,levels=c("hot","cold"))
    
    
    
    library(ggpubr)
    compare_means(value ~ Group, data = dat_1, group.by = "Gene",method = "anova")
    
    
    #pdf(paste(cancer,"/",cancer,"_cluster_activities_",cell,".pdf",sep=''),width=3,height = 8)
    p <- ggboxplot(dat_1, x = "Group", y = "value",
                   color = "Group", palette = "FA7F6F", 
                   add = "jitter") # palette可以按照期刊选择相应的配色，如"npg"等
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
    
    dat_1 <- dat_1[,-1]
    library(dplyr)
    means=group_by(dat_1, Group) %>% summarise_each(funs(mean))
    
    if(length(mean)==0){
      mean <- means
    }else {
      mean <- cbind(mean,means)
    }
  }
  pdf(paste(cancer,"/",cancer,"_xCELL_fen.pdf",sep=''),width = 8,height = 8*length(out_put_list)/4)
  
  eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=4))',sep='')))
  #eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=4,labels=out_put_list,label_size=40))',sep='')))
  
  dev.off()
  
  n1 <- grep("Group",colnames(mean))
  mean <- mean[,-n1]
  colnames(mean) <- unique(dat$Gene)
  rownames(mean) <- c("hot","cold")
  
  mean <-  t(mean) 
  
  updown <- rep(0, times=length(rownames(mean)))
  
  for( i in 1:length(rownames(mean))){
    if( mean[i,1] > mean[i,2] ){
      updown[i]= 1
    }else{
      updown[i]= -1
    }
  }
  
  xCELL_updown <- data.frame(c("T_cells_CD8","NK_cells","Macrophages","Macrophages_M1","Macrophages_M2"),updown)
  colnames(xCELL_updown) <- c("CELL","xCELL_ud")
  
  
  
  
  
  
  
  
  
  
  
  #######################################################
  
  
  QUANTISEQ <- fread(paste("D:/R/Hot_tumor_and_cold_tumor/Immunecells/8cancers/",cancer,"_QUANTISEQ.csv",sep=''))
  QUANTISEQ <- as.data.frame(QUANTISEQ)
  
  QUANTISEQ1 <- data.frame(QUANTISEQ$V1,QUANTISEQ$T_cells_CD8,QUANTISEQ$NK_cells,QUANTISEQ$Macrophages_M1,QUANTISEQ$Macrophages_M2)
  rownames(QUANTISEQ1) <- QUANTISEQ1[,1]
  QUANTISEQ1 <- QUANTISEQ1[,-1]
  colnames(QUANTISEQ1) <- c("T_cells_CD8","NK_cells","Macrophages_M1","Macrophages_M2")
  
  all_name <- names(which(table(c(rownames(QUANTISEQ1),rownames(cluster)))==2))
  
  QUANTISEQ2 <- QUANTISEQ1[match(all_name,rownames(QUANTISEQ1)),]
  dat_cluster <- cluster[match(all_name,rownames(cluster)),]
  
  dat_cluster <- dat_cluster[c('V2')]
  
  dat <- data.frame()
  for(coln in colnames(QUANTISEQ2))
  {
    for(i in all_name)
    {
      dat <- rbind(dat,c(coln,dat_cluster[match(i,rownames(dat_cluster)),],QUANTISEQ2[c(coln)][match(i,rownames(QUANTISEQ2)),]))
    }
  }
  
  dat[,3] = as.numeric(dat[,3])
  
  colnames(dat) <- c("Gene","Group","value")
  
  dat$Group <- factor(dat$Group,levels=c("hot","cold"))
  
  xxxx <- matrix(unique(dat[,1])) #colname
  
  temp_length <- length(xxxx)*2
  
  
  library(ggpubr)
  xx <- compare_means(value ~ Group, data = dat, group.by = "Gene",method = "anova")
  
  QUANTISEQ_p <- data.frame(xx$Gene,xx$p)
  colnames(QUANTISEQ_p) <- c("CELL","QUANTISEQ_p")
  
  pdf(paste(cancer,"/",cancer,"_QUANTISEQ_ns.pdf",sep=''),width=temp_length/3.5,height = 8)
  p <- ggboxplot(dat, x = "Gene", y = "value",
                 color = "Group", palette = "FA7F6F", 
                 add = "jitter",x.text.angle=60) # palette可以按照期刊选择相应的配色，如"npg"等
  p <- p+xlab("Cell Type")+ylab("Expression Value")
  p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
  p
  print(p)
  dev.off()
  
  out_put_list <- c()
  mean <- c()
  
  
  for(cell in unique(dat$Gene)){
    
    n1 <- grep(cell,dat$Gene)
    
    dat_1 <- dat[n1,]
    
    dat_1$Group <- factor(dat_1$Group,levels=c("hot","cold"))
    
    
    
    library(ggpubr)
    compare_means(value ~ Group, data = dat_1, group.by = "Gene",method = "anova")
    
    
    #pdf(paste(cancer,"/",cancer,"_cluster_activities_",cell,".pdf",sep=''),width=3,height = 8)
    p <- ggboxplot(dat_1, x = "Group", y = "value",
                   color = "Group", palette = "FA7F6F", 
                   add = "jitter") # palette可以按照期刊选择相应的配色，如"npg"等
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
    
    dat_1 <- dat_1[,-1]
    library(dplyr)
    means=group_by(dat_1, Group) %>% summarise_each(funs(mean))
    
    if(length(mean)==0){
      mean <- means
    }else {
      mean <- cbind(mean,means)
    }
  }
  pdf(paste(cancer,"/",cancer,"_QUANTISEQ_fen.pdf",sep=''),width = 8,height = 8*length(out_put_list)/4)
  
  eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=4))',sep='')))
  #eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=4,labels=out_put_list,label_size=40))',sep='')))
  
  dev.off()
  
  n1 <- grep("Group",colnames(mean))
  mean <- mean[,-n1]
  colnames(mean) <- unique(dat$Gene)
  rownames(mean) <- c("hot","cold")
  
  mean <-  t(mean) 
  
  updown <- rep(0, times=length(rownames(mean)))
  
  for( i in 1:length(rownames(mean))){
    if( mean[i,1] > mean[i,2] ){
      updown[i]= 1
    }else{
      updown[i]= -1
    }
  }
  
  QUANTISEQ_updown <- data.frame(c("T_cells_CD8","NK_cells","Macrophages_M1","Macrophages_M2"),updown)
  colnames(QUANTISEQ_updown) <- c("CELL","QUANTISEQ")
  
  file=ls(pattern = "data")
  ALL1=list(BIBER_p,BIBER_p,TIMER_p,EPIC_p,MCP_p,xCELL_p,QUANTISEQ_p)
  #,dataGrade,dataM,dataN,dataStage,dataT
  
  multimerge<-function(dat=list(ALL1),...){
    if(length(dat)<2)return(as.data.frame(dat))
    mergedat<-rownames(dat)
    #dat[[1]]<-NULL
    for(i in dat){
      mergedat<-merge(all=TRUE,mergedat,i,...)
    }
    return(mergedat)
  }
  p_value <- multimerge(ALL1)
  rownames(p_value) <- p_value[,1]
  p_value <- p_value[,-1]
  colnames(p_value) <- c("CIBERSORT","TIMER","EPIC","MCPcounter","xCELL","QUANTISEQ")
  
  
  file=ls(pattern = "data")
  ALL1=list(BIBER_updown,BIBER_updown,TIMER_updown,EPIC_updown,MCP_updown,xCELL_updown,QUANTISEQ_updown)
  #,dataGrade,dataM,dataN,dataStage,dataT
  
  multimerge<-function(dat=list(ALL1),...){
    if(length(dat)<2)return(as.data.frame(dat))
    mergedat<-rownames(dat)
    #dat[[1]]<-NULL
    for(i in dat){
      mergedat<-merge(all=TRUE,mergedat,i,...)
    }
    return(mergedat)
  }
  up_down <- multimerge(ALL1)
  rownames(up_down) <- up_down[,1]
  up_down <- up_down[,-1]
  colnames(up_down) <- c("CIBERSORT","TIMER","EPIC","MCPcounter","xCELL","QUANTISEQ")
  
  df <- -log10(p_value)
  
  #df <- t(scale(t(df)))
  #df <- scale(df)
  
  for(i in 1:length(rownames(up_down))){
    for(t in 1:length(colnames(up_down))){
      if( is.na(df[i,t])==TRUE){
        next
      }else{
        df[i,t] <- as.numeric(df[i,t])*as.numeric(up_down[i,t])
      }
    }
  }
  
  df[is.na(df)] = 0
  p_value[is.na(p_value)] = 100
  
  library(pheatmap)
  
  pdf(paste(cancer,"/",cancer,"_immunecells.pdf",sep=''),5,length(rownames(df))/2)
  
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
  bk <- c(seq( -100, -0.1,by=0.01),seq(0,100,by=0.01))
  #bk <- c(seq( -max(abs(max(df)),abs(min(df))), -0.1,by=0.01),seq(0,max(abs(max(df)),abs(min(df))),by=0.01))
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
  
  
  
  
  
}
