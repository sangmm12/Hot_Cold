setwd("D:/R/Hot_tumor_and_cold_tumor/组织芯片2/exp/exp2")

library(ggplot2)
library(ggpubr)
library(ggstatsplot)

library(data.table)
dat <- fread(paste("D:/R/Hot_tumor_and_cold_tumor/组织芯片2/原始分析结果1.csv",sep=''))
dat <- as.data.frame(dat)
rownames(dat) <- dat[,1]
dat<- dat[,-1]

dat <- t(dat)

aa <- colnames(dat)

aa <- colnames(dat)
#aa <- gsub("[^A-Za-z0-9_]", "_", list2)
aa[16:22]
#list1 <- c("PANCK+CD163+ Cells","PANCK-CD163+ Cells","CD163+CD276+ Cells","CD163+CD276- Cells","CD163+CD73+ Cells","CD163+CD73- Cells","CD276+CD73+ Cells")
list1 <- aa[16:22]
list2 <- aa[33:39]
list3 <- aa[67:73]
list4 <- aa[84:90]
list5 <- aa[118:124]
list6 <- aa[135:141]

lists <- c(list1,list2,list3,list4,list5,list6)


data0 <- dat[,match(lists,colnames(dat))]

lists <- gsub("%", "percent", lists)
lists <- gsub(":", "_", lists)
lists <- gsub(" ", "_", lists)
lists <- gsub("-", "0", lists)
lists <- gsub("[^A-Za-z0-9_]", "_", lists)

colnames(data0) <- lists


list1 <- gsub("%", "percent", list1)
list1 <- gsub(":", "_", list1)
list1 <- gsub(" ", "_", list1)
list1 <- gsub("-", "0", list1)
list1 <- gsub("[^A-Za-z0-9_]", "_", list1)

list2 <- gsub("%", "percent", list2)
list2 <- gsub(":", "_", list2)
list2 <- gsub(" ", "_", list2)
list2 <- gsub("-", "0", list2)
list2 <- gsub("[^A-Za-z0-9_]", "_", list2)

list3 <- gsub("%", "percent", list3)
list3 <- gsub(":", "_", list3)
list3 <- gsub(" ", "_", list3)
list3 <- gsub("-", "0", list3)
list3 <- gsub("[^A-Za-z0-9_]", "_", list3)

list4 <- gsub("%", "percent", list4)
list4 <- gsub(":", "_", list4)
list4 <- gsub(" ", "_", list4)
list4 <- gsub("-", "0", list4)
list4 <- gsub("[^A-Za-z0-9_]", "_", list4)

list5 <- gsub("%", "percent", list5)
list5 <- gsub(":", "_", list5)
list5 <- gsub(" ", "_", list5)
list5 <- gsub("-", "0", list5)
list5 <- gsub("[^A-Za-z0-9_]", "_", list5)

list6 <- gsub("%", "percent", list6)
list6 <- gsub(":", "_", list6)
list6 <- gsub(" ", "_", list6)
list6 <- gsub("-", "0", list6)
list6 <- gsub("[^A-Za-z0-9_]", "_", list6)

out_put_list1 <- c()
out_put_list2 <- c()
  
for(i in 1:7){  

  #i=1
  dat1 <- data0[,match(c(list1[i],list3[i],list5[i]),colnames(data0))]
  class1 <- rep("ALL", times=length(rownames(dat1)))
  class2 <- rep("Tumor", times=length(rownames(dat1)))
  class3 <- rep("Stroma", times=length(rownames(dat1)))
  class <- c(class1,class2,class3)
  cell1  <- unlist(dat1[,1])
  cell2  <- unlist(dat1[,2])
  cell3  <- unlist(dat1[,3])
  cell <- c(cell1,cell2,cell3)
  
  dat <- data.frame(class,cell)
  
  comparisons <- list(c("ALL","Tumor"), c("Tumor","Stroma"), c("ALL","Stroma"))
  
  pdf(file=paste("violin_",list1[i],".pdf",sep=''),onefile = FALSE,width = 6.5,height =7.5)
  p <- ggviolin(dat, x = "class", y = "cell", fill = "class", palette = c("lancet"),
                add = "boxplot", add.params = list(fill = "white"), order = c("ALL","Tumor","Stroma"),
                error.plot = "errorbar") + stat_compare_means(comparisons = comparisons)+
    labs(x=list1[i],y="Expression")+
    theme(axis.text.x = element_text(size = 20, angle = 0, hjust = 1))+
    theme(axis.title = element_text(size = 25))+
    theme(axis.text.y = element_text(size = 20))# 调整字体大小和倾斜度 
  print(p)
  
  eval(parse(text = paste(list1[i],' <- p',sep='')))
  out_put_list1 <- append(out_put_list1,list1[i])
  
  dev.off()
  
}


pdf(paste("cell.pdf",sep=''),width = 70/2,height = 5*length(out_put_list1)/7)

eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list1, collapse = ","),', ncol=7))',sep='')))
#eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=4,labels=out_put_list,label_size=40))',sep='')))

dev.off()




for(i in 1:7){  
  
  
  dat1 <- data0[,match(c(list2[i],list4[i],list6[i]),colnames(data0))]
  class1 <- rep("ALL", times=length(rownames(dat1)))
  class2 <- rep("Tumor", times=length(rownames(dat1)))
  class3 <- rep("Stroma", times=length(rownames(dat1)))
  class <- c(class1,class2,class3)
  cell1  <- unlist(dat1[,1])
  cell2  <- unlist(dat1[,2])
  cell3  <- unlist(dat1[,3])
  cell <- c(cell1,cell2,cell3)
  
  dat <- data.frame(class,cell)
  
  comparisons <- list(c("ALL","Tumor"), c("Tumor","Stroma"), c("ALL","Stroma"))
  
  pdf(file=paste("violin_百分比",list1[i],".pdf",sep=''),onefile = FALSE,width = 7,height =7.5)
  p <- ggviolin(dat, x = "class", y = "cell", fill = "class", palette = c("lancet"),
                add = "boxplot", add.params = list(fill = "white"), order = c("ALL","Tumor","Stroma"),
                error.plot = "errorbar") + stat_compare_means(comparisons = comparisons)+
    labs(x=list2[i],y="Expression")+
    theme(axis.text.x = element_text(size = 20, angle = 0, hjust = 1))+
    theme(axis.title = element_text(size = 25))+
    theme(axis.text.y = element_text(size = 20))# 调整字体大小和倾斜度 
  print(p)
  
  eval(parse(text = paste(list1[i],' <- p',sep='')))
  out_put_list2 <- append(out_put_list2,list1[i])
  
  dev.off() 
  
  
}


pdf(paste("cell2.pdf",sep=''),width = 70/2,height = 5*length(out_put_list2)/7)

eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list2, collapse = ","),', ncol=7))',sep='')))
#eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=4,labels=out_put_list,label_size=40))',sep='')))

dev.off()


