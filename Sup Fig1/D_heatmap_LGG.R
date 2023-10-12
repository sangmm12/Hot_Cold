setwd('D:/R/Hot_tumor_and_cold_tumor/fenqi')


rm(list=ls())


cancer <- c("LGG")

cluster <- read.csv( paste("D:/R/Hot_tumor_and_cold_tumor/Immunecells/CIBERSORT/OVER100/",cancer,"/",cancer,"_clusterresult/",cancer,"_cluster_hc.csv",sep=''),header = F)
colnames(cluster) <- c("SampleName","Subtype")

datAge <- read.csv("D:/R/Hot_tumor_and_cold_tumor/fenqi/Age.csv")
n1 <- grep(cancer,datAge$X)
dataAge <- datAge[n1,]
n2 <- grep("GBMLGG",dataAge$X)
dataAge <- dataAge[-n2,]
dataAge <- dataAge[,-1:-2]
colnames(dataAge) <- c("Age","SampleName")

datGender <- read.csv("D:/R/RNA_mood_PAAD/fig2/Gender.csv")
n1 <- grep(cancer,datGender$X)
dataGender <- datGender[n1,]
n2 <- grep("GBMLGG",dataGender$X)
dataGender <- dataGender[-n2,]
dataGender <- dataGender[,-1]

datGrade<- read.csv("D:/R/RNA_mood_PAAD/fig2/Grade.csv")
n1 <- grep(cancer,datGrade$X)
dataGrade <- datGrade[n1,]
n2 <- grep("GBMLGG",dataGrade$X)
dataGrade <- dataGrade[-n2,]
dataGrade <- dataGrade[,-1]

# datM<- read.csv("D:/R/RNA_mood_PAAD/fig2/M.csv")
# n1 <- grep(cancer,datM$X)
# dataM <- datM[n1,]
# n2 <- grep("GBMLGG",dataM$X)
# dataM <- dataM[-n2,]
# dataM <- dataM[,-1]

# datN<- read.csv("D:/R/RNA_mood_PAAD/fig2/N.csv")
# n1 <- grep(cancer,datN$X)
# dataN <- datN[n1,]
# n2 <- grep("GBMLGG",dataN$X)
# dataN <- dataN[-n2,]
# dataN <- dataN[,-1]

# datStage<- read.csv("D:/R/RNA_mood_PAAD/fig2/Stage.csv")
# n1 <- grep(cancer,datStage$X)
# dataStage <- datStage[n1,]
# n2 <- grep("GBMLGG",dataStage$X)
# dataStage <- dataStage[-n2,]
# dataStage <- dataStage[,-1]

# datT<- read.csv("D:/R/RNA_mood_PAAD/fig2/T.csv")
# n1 <- grep(cancer,datT$X)
# dataT <- datT[n1,]
# n2 <- grep("GBMLGG",dataT$X)
# dataT <- dataT[-n2,]
# dataT <- dataT[,-1]

file=ls(pattern = "data")
ALL1=list(dataAge,dataGender,dataGrade,cluster)
#dataGrade,

multimerge<-function(dat=list(ALL1),...){
  if(length(dat)<2)return(as.data.frame(dat))
  mergedat<-dat[[1]]
  dat[[1]]<-NULL
  for(i in dat){
    mergedat<-merge(all=TRUE,mergedat,i,...)
  }
  return(mergedat)
}
annotation_col <- multimerge(ALL1)
rownames(annotation_col) <- annotation_col[,1]
annotation_col <- annotation_col[,-1]

library(data.table)

dat <-  fread(paste("D:/R/Hot_tumor_and_cold_tumor/cluster_cell/cluster_cellzong/",cancer,"/data_",cancer,".csv",sep=''),header = T)

dat <- as.data.frame(dat)

rownames(dat) <- dat[,1]
dat <- dat[,-1]



annotation_col <- annotation_col[order(annotation_col$Subtype),]

all_name <- names(which(table(c(rownames(annotation_col),rownames(dat)))==2))


annotation_col <- annotation_col[match(all_name,rownames(annotation_col)),]

annotation_col <- annotation_col[order(annotation_col$Subtype),]


RNA_data1 <- dat[match(rownames(annotation_col),rownames(dat)),]



colnames(annotation_col)
unique(annotation_col$Subtype)
unique(annotation_col$Range)
unique(annotation_col$Grade)
unique(annotation_col$M)
unique(annotation_col$N)
unique(annotation_col$Stage)
unique(annotation_col$T)



ann_colors = list(
  Range = c(">60"="#FF9D5A","<=60"= "#A1C800"),
  Gender = c(MALE = "#E7298A", FEMALE = "#7570B3"),
  #Grade = c(G1 = "#E75a8A", G2 = "#75b0B3",G3 = "#D6ECF4"),
  M = c(M0 = "#FF88B1", M1 = "#C0A6FF"),
  N = c(N0 = "#00DCAF", N1 = "#C4BF00",N2 = "#FF82FA",N3 = "#93C0DB"),
  Stage = c("Stage I" = "#00D1FF", "Stage II" = "#F7A900", "Stage III" = "#EA95FF", "Stage IV" = "#70CF07"),
  T = c(T1 = "#82B7FF", T2 = "#FF81D4",T3 = "#00C6FF",T4 = "#00D8EE"),
  Subtype = c(cold = "#1B9E77", hot = "#FF9289"),
  class = c("Immuneinfiltration" = "#FFC426", "Immunescore" = "#009999","activities"= "#E67E57")
)

class1 <- rep("Immuneinfiltration", times=22)
class2 <- rep("Immunescore", times=3)
class3 <- rep("activities", times=16)
class<- c(class1,class2,class3)

annotation_row <- data.frame(class)

rownames(annotation_row) <- colnames(RNA_data1)
#annotation_row <- annotation_row[,-1]


bk <- c(seq(-3, -0.1,by=0.01),seq(0,3,by=0.01))
myColor <- c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2))

RNA_data2 <- scale(RNA_data1)

library(pheatmap)
library(limma)
library(colorRamp2)
pdf(paste(cancer,"/heatmap_",cancer,"_hc3.pdf",sep=''),width = 10,height = 8)
pheatmap(t(RNA_data2),
         color=myColor,
         breaks=bk,
         annotation_col=annotation_col,annotation_row = annotation_row,annotation_colors = ann_colors,
         
         # col=colorRamp2(c(-2,0,2),c('#21b6af','white','#eeba4d')),
         cluster_row=F,show_colnames = F,
         cluster_col=F)

dev.off()     

