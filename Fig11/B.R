setwd("D:/R/Hot_tumor_and_cold_tumor/组织芯片2/KM")


library(data.table)
dat <- fread(paste("D:/R/Hot_tumor_and_cold_tumor/组织芯片2/生存.csv",sep=''))
dat <- as.data.frame(dat)
dat<- dat[1:79,]
rownames(dat) <- dat[,1]
dat<- dat[,-1]


colnames(dat) <- c("year1","month1","day1","year2","month2","day2")
dates1 <- with(dat, as.Date(paste(year1, month1, day1, sep = "-"), format = "%Y-%m-%d"))
dates2 <- with(dat, as.Date(paste(year2, month2, day2, sep = "-"), format = "%Y-%m-%d"))

df <- data.frame(dates1,dates2)
rownames(df) <- rownames(dat)

df$days_difference <- as.integer(df$dates2 - df$dates1)



SN <- fread(paste("D:/R/Hot_tumor_and_cold_tumor/组织芯片2/SampleName1.csv",sep=''))
SN <- as.data.frame(SN)


fenqi <- fread(paste("D:/R/Hot_tumor_and_cold_tumor/组织芯片2/分期.csv",sep=''))
fenqi <- as.data.frame(fenqi)
fenqi<- fenqi[1:79,]
rownames(fenqi) <- fenqi[,1]
fenqi<- fenqi[,-1]


all_name <- names(which(table(c(rownames(fenqi),rownames(df),SN[,1]))==3))

SN1 <- SN[match(all_name,SN[,1]),]
fenqi1 <- fenqi[match(all_name,rownames(fenqi)),]
df1 <- df[match(all_name,rownames(df)),]

data <- data.frame(SN1,fenqi1,df1$days_difference)


sur <- data[,c(2,9,10)]
rownames(sur) <- sur[,1]
sur <- sur[,-1]
# sur[,2] <- as.numeric(unlist(sur[,2]))
# sur[,1] <- as.numeric(unlist(sur[,1]))



dat <- fread(paste("D:/R/Hot_tumor_and_cold_tumor/组织芯片2/原始分析结果1.csv",sep=''))
dat <- as.data.frame(dat)
rownames(dat) <- dat[,1]
dat<- dat[,-1]

dat <- t(dat)
aa <- colnames(dat)

list1 <- c("CD8 Positive Cells","PANCK Positive Cells","CD163 Positive Cells","CD276 Positive Cells","CD73 Positive Cells")
list2 <- c("% CD8 Positive Cells","% PANCK Positive Cells","% CD163 Positive Cells","% CD276 Positive Cells","% CD73 Positive Cells" )
list3 <- c("tumor: CD8 Positive Cells","tumor: PANCK Positive Cells","tumor: CD163 Positive Cells","tumor: CD276 Positive Cells","tumor: CD73 Positive Cells" )
list4 <- c("tumor: % CD8 Positive Cells","tumor: % PANCK Positive Cells","tumor: % CD163 Positive Cells","tumor: % CD276 Positive Cells","tumor: % CD73 Positive Cells")   
list5 <- c("stroma: CD8 Positive Cells","stroma: PANCK Positive Cells","stroma: CD163 Positive Cells","stroma: CD276 Positive Cells","stroma: CD73 Positive Cells")
list6 <- c("stroma: % CD8 Positive Cells","stroma: % PANCK Positive Cells","stroma: % CD163 Positive Cells","stroma: % CD276 Positive Cells","stroma: % CD73 Positive Cells")
list7 <- aa[40:44]
list8 <- aa[91:95]
list9 <- aa[142:146]
lists <- c(list1,list2,list3,list4,list5,list6,list7,list8,list9)


dat <- dat[,match(lists,colnames(dat))]

lists <- gsub("%", "percent", lists)
lists <- gsub(":", "_", lists)
lists <- gsub(" ", "_", lists)
lists <- gsub("-", "0", lists)
lists <- gsub("[^A-Za-z0-9_]", "_", lists)

colnames(dat) <- lists

all_name <- names(which(table(c(rownames(dat),rownames(sur) ))==2))


sur1 <- sur[match(all_name,rownames(sur)),]
class(sur1)
dat1 <- dat[match(all_name,rownames(dat)),]
dat1 <- as.data.frame(dat1)
class(dat1)
dataset0 <- data.frame(Time=sur1[,2],Status=sur1[,1],dat1)

write.csv(dataset0,file = paste("sur_exp.csv",sep=''))




library ("maxstat")

library("survival")

#dataset <- fread(paste("sur_exp.csv",sep=''),header = T)
dataset0 <- as.data.frame(dataset0)

connam <- colnames(dataset0)

cutofflist <- c()
aaaa <- c()

out_put_list <- c()
#i=13
for(i in 3:length(colnames(dataset0))){
  
  if(sum(is.na(dataset0[,i])==TRUE)==length(rownames(dataset0))){
    aa <- rep(NA, times=6)
    cutoff <- NA
    #if(i==14|i==25|i==142|i==158|i==185|i==191){
    # if(i==57|i==58|i==59|i==60|i==61|i==64){
    #   aa <- rep(NA, times=6)
    
  }else if(sum(dataset0[,i]==0)> length(rownames(dataset0))/1.4){
    aa <- rep(NA, times=6)
    cutoff <- NA
    #if(i==14|i==25|i==142|i==158|i==185|i==191){
    # if(i==57|i==58|i==59|i==60|i==61|i==64){
    #   aa <- rep(NA, times=6)
    
  }else{
    cutoff<- maxstat.test(Surv(Time, Status)~dataset0[,i], data=dataset0, smethod="LogRank", pmethod="Lau94",minprop=0.25, maxprop=0.75)
    
    # cutoff
    # 
    # plot(cutoff)
    
    cutoff<- as.numeric(cutoff$estimate)
    
    dataset <- data.frame(dataset0[,1:2],dataset0[,i])
    
    library("glmnet")
    library("survival")
    library("survminer")
    
    
    rt <- dataset
    
    rt$risk <- c("Low")
    rt$risk[which(rt[,3]>cutoff)] <- c("High")
    
    rt$risk <- factor(rt$risk,levels=c("Low","High"))
    
    rt$Time  =rt$Time
    #rt$Time  =rt$Time /365
    
    fit <- survfit(Surv(Time, Status) ~ risk, data = rt)
    
    ggsurvplot(fit, data = rt)
    
    
    library(survival)
    library(survminer)
    
    fit <- survfit(Surv(Time, Status) ~ risk, data = rt)
    
    
    surPlot=ggsurvplot(fit,
                       data=rt,
                       #font.title = paste(connam[i],sep=''),
                       #ggtitle = paste(connam[i],sep=''),
                       #conf.int=TRUE,
                       legend.labs=c( "L","H"),
                       legend = "top",
                       legend.title="Risk",
                       pval=TRUE,
                       #pval.method = TRUE,
                       pval.size = 5 ,
                       xlab="Time",
                       break.time.by = ceiling((max(rt$Time))/4),
                       risk.table.title="",
                       palette=c("#009E73","red"),
                       risk.table=T,
                       risk.table.height=.25,)
    
    surPlot <- surPlot +
      labs(title =paste(connam[i],sep=''))
    
    cell <- connam[i]
    
    pdf(file=paste("KM_best_",connam[i],".pdf",sep=''),onefile = FALSE,width = 6,height =7.5)
    print(surPlot)
    
    # eval(parse(text = paste(cell,'<- surPlot')))
    # out_put_list <- append(out_put_list,cell)

    dev.off()


    
    
    data.survdiff <- survdiff(Surv(Time, Status) ~ risk, data = rt)
    p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
    p.val = round(p.val, 8)
    HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
    HR = round(HR, 4)
    up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
    up95 = round(up95, 4)
    low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])) 
    low95 = round(low95, 4)
    
    aa <- c(p.val,HR,low95,up95,paste0(HR," (",low95,"-",up95,")"),cutoff)
    
  }
  
  aaaa <- rbind(aaaa,aa)
  cutofflist <- rbind(cutofflist,cutoff)
  
}


colnames(aaaa) <- c("pvalue","HR","Lower","Upper","Hazard Ratio(95%CI)","cutoff")
rownames(aaaa) <- colnames(dat)

write.csv(aaaa,file =paste("sur_exp_single.csv",sep=''),quote=T)

# 
# pdf(paste("KM1.pdf",sep=''),width = 50,height = 10*length(out_put_list)/5)
# 
# eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=5))',sep='')))
# #eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=4,labels=out_put_list,label_size=40))',sep='')))
# 
# dev.off()

