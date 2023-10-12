setwd("D:/R/Hot_tumor_and_cold_tumor/immune data1/gene/KM/KM4")



library(data.table)
dat <- fread("D:/R/Hot_tumor_and_cold_tumor/immune data1/gene/KM/sur_PAAD_NT5E&VEGFA_hypoxia.csv",header = T)
dat <- as.data.frame(dat)
rownames(dat) <- dat[,1]
dat <- dat[,-1]

cancer <- c("PAAD")
cluster <- read.csv(paste("D:/R/Hot_tumor_and_cold_tumor/Immunecells/CIBERSORT/OVER100/",cancer,"/",cancer,"_clusterresult/",cancer,"_cluster_hc.csv",sep=''),header = F)



all_name <- names(which(table(c(rownames(dat),cluster[,1] ))==2))


dat_cluster <- cluster[match(all_name,cluster[,1]),]

data<- dat[match(all_name,rownames(dat)),1:5]

dat<- data.frame(data,cluster = dat_cluster[,2])



#########################################################################
high_high <- subset(dat,dat$NT5E>=46.11087 & dat$cluster=="hot")
high_low <- subset(dat,dat$NT5E>=46.11087 & dat$cluster=="cold")
low_high <- subset(dat,dat$NT5E<46.11087 & dat$cluster=="hot")
low_low <- subset(dat,dat$NT5E<46.11087 & dat$cluster=="cold")

class1 <- rep("NT5Ehigh-hot", times=length(rownames(high_high)))
class2 <- rep("NT5Ehigh-cold", times=length(rownames(high_low)))
class3 <- rep("NT5Elow-hot", times=length(rownames(low_high)))
class4 <- rep("NT5Elow-cold", times=length(rownames(low_low)))

length(rownames(high_high))
length(rownames(high_low))
length(rownames(low_high))
length(rownames(low_low))

group<- c(class1,class2,class3,class4)


KM4 <- rbind(high_high,high_low,low_high,low_low)

KM4 <- cbind(group,KM4)

write.csv(KM4,file=paste("sur_NT5Ehc_best",".csv",sep=''),quote=F)






group2<- c(class2,class3)

KM42 <- rbind(high_low,low_high)

KM42 <- cbind(group=group2,KM42)

write.csv(KM42,file=paste("sur_NT5Ehc2_best",".csv",sep=''),quote=F)


#########################################################################
high_high <- subset(dat,dat$VEGFA>=56.15423 & dat$cluster=="hot")
high_low <- subset(dat,dat$VEGFA>=56.15423 & dat$cluster=="cold")
low_high <- subset(dat,dat$VEGFA<56.15423 & dat$cluster=="hot")
low_low <- subset(dat,dat$VEGFA<56.15423 & dat$cluster=="cold")

class1 <- rep("VEGFAhigh-hot", times=length(rownames(high_high)))
class2 <- rep("VEGFAhigh-cold", times=length(rownames(high_low)))
class3 <- rep("VEGFAlow-hot", times=length(rownames(low_high)))
class4 <- rep("VEGFAlow-cold", times=length(rownames(low_low)))

length(rownames(high_high))
length(rownames(high_low))
length(rownames(low_high))
length(rownames(low_low))



group<- c(class1,class2,class3,class4)


KM4 <- rbind(high_high,high_low,low_high,low_low)

KM4 <- cbind(group,KM4)

write.csv(KM4,file=paste("sur_VEGFAhc_best",".csv",sep=''),quote=F)




group2<- c(class2,class3)

KM42 <- rbind(high_low,low_high)

KM42 <- cbind(group=group2,KM42)

write.csv(KM42,file=paste("sur_VEGFAhc2_best",".csv",sep=''),quote=F)







#########################################################################
high_high <- subset(dat,dat$CD276>=55.79723 & dat$cluster=="hot")
high_low <- subset(dat,dat$CD276>=55.79723 & dat$cluster=="cold")
low_high <- subset(dat,dat$CD276<55.79723 & dat$cluster=="hot")
low_low <- subset(dat,dat$CD276<55.79723 & dat$cluster=="cold")

class1 <- rep("CD276high-hot", times=length(rownames(high_high)))
class2 <- rep("CD276high-cold", times=length(rownames(high_low)))
class3 <- rep("CD276low-hot", times=length(rownames(low_high)))
class4 <- rep("CD276low-cold", times=length(rownames(low_low)))

length(rownames(high_high))
length(rownames(high_low))
length(rownames(low_high))
length(rownames(low_low))



group<- c(class1,class2,class3,class4)


KM4 <- rbind(high_high,high_low,low_high,low_low)

KM4 <- cbind(group,KM4)

write.csv(KM4,file=paste("sur_CD276hc_best",".csv",sep=''),quote=F)




group2<- c(class2,class3)

KM42 <- rbind(high_low,low_high)

KM42 <- cbind(group=group2,KM42)

write.csv(KM42,file=paste("sur_CD276hc2_best",".csv",sep=''),quote=F)


