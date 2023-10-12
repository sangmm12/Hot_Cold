setwd("D:/R/Hot_tumor_and_cold_tumor/immune data1/gene/KM/KM4")


library(data.table)
dat <- fread("D:/R/Hot_tumor_and_cold_tumor/immune data1/gene/KM/sur_PAAD_NT5E&VEGFA_hypoxia.csv",header = T)
dat <- as.data.frame(dat)





#########################################################################
high_high <- subset(dat,dat$NT5E>=46.11087 & dat$VEGFA>=56.15423)
high_low <- subset(dat,dat$NT5E>=46.11087 & dat$VEGFA<56.15423)
low_high <- subset(dat,dat$NT5E<46.11087 & dat$VEGFA>=56.15423)
low_low <- subset(dat,dat$NT5E<46.11087 & dat$VEGFA<56.15423)

class1 <- rep("NT5Ehigh-VEGFAhigh", times=length(rownames(high_high)))
class2 <- rep("NT5Ehigh-VEGFAlow", times=length(rownames(high_low)))
class3 <- rep("NT5Elow-VEGFAhigh", times=length(rownames(low_high)))
class4 <- rep("NT5Elow-VEGFAlow", times=length(rownames(low_low)))


group<- c(class1,class2,class3,class4)


KM4 <- rbind(high_high,high_low,low_high,low_low)

KM4 <- cbind(group,KM4)

write.csv(KM4,file=paste("sur_NT5E&VEGFA4_best",".csv",sep=''),quote=F)




#########################################################################
high_high <- subset(dat,dat$VEGFA>=56.15423 & dat$Hypoxia>=4.454189)
high_low <- subset(dat,dat$VEGFA>=56.15423 & dat$Hypoxia<4.454189)
low_high <- subset(dat,dat$VEGFA<56.15423 & dat$Hypoxia>=4.454189)
low_low <- subset(dat,dat$VEGFA<56.15423 & dat$Hypoxia<4.454189)

class1 <- rep("VEGFAhigh-Hypoxiahigh", times=length(rownames(high_high)))
class2 <- rep("VEGFAhigh-Hypoxialow", times=length(rownames(high_low)))
class3 <- rep("VEGFAlow-Hypoxiahigh", times=length(rownames(low_high)))
class4 <- rep("VEGFAlow-Hypoxialow", times=length(rownames(low_low)))


group<- c(class1,class2,class3,class4)


KM4 <- rbind(high_high,high_low,low_high,low_low)

KM4 <- cbind(group,KM4)

write.csv(KM4,file=paste("sur_VEGFA&Hypoxia4_best",".csv",sep=''),quote=F)








#########################################################################
high_high <- subset(dat,dat$NT5E>=46.11087 & dat$Hypoxia>=4.454189)
high_low <- subset(dat,dat$NT5E>=46.11087 & dat$Hypoxia<4.454189)
low_high <- subset(dat,dat$NT5E<46.11087 & dat$Hypoxia>=4.454189)
low_low <- subset(dat,dat$NT5E<46.11087 & dat$Hypoxia<4.454189)

class1 <- rep("NT5Ehigh-Hypoxiahigh", times=length(rownames(high_high)))
class2 <- rep("NT5Ehigh-Hypoxialow", times=length(rownames(high_low)))
class3 <- rep("NT5Elow-Hypoxiahigh", times=length(rownames(low_high)))
class4 <- rep("NT5Elow-Hypoxialow", times=length(rownames(low_low)))


group<- c(class1,class2,class3,class4)


KM4 <- rbind(high_high,high_low,low_high,low_low)

KM4 <- cbind(group,KM4)

write.csv(KM4,file=paste("sur_NT5E&Hypoxia4_best",".csv",sep=''),quote=F)







#########################################################################
high_high <- subset(dat,dat$CD276>=55.79723 & dat$Hypoxia>=4.454189)
high_low <- subset(dat,dat$CD276>=55.79723 & dat$Hypoxia<4.454189)
low_high <- subset(dat,dat$CD276<55.79723 & dat$Hypoxia>=4.454189)
low_low <- subset(dat,dat$CD276<55.79723 & dat$Hypoxia<4.454189)

class1 <- rep("CD276high-Hypoxiahigh", times=length(rownames(high_high)))
class2 <- rep("CD276high-Hypoxialow", times=length(rownames(high_low)))
class3 <- rep("CD276low-Hypoxiahigh", times=length(rownames(low_high)))
class4 <- rep("CD276low-Hypoxialow", times=length(rownames(low_low)))


group<- c(class1,class2,class3,class4)


KM4 <- rbind(high_high,high_low,low_high,low_low)

KM4 <- cbind(group,KM4)

write.csv(KM4,file=paste("sur_CD276&Hypoxia4_best",".csv",sep=''),quote=F)



