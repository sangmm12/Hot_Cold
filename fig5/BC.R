
setwd("D:/R/Hot_tumor_and_cold_tumor/immune data1/DGEs")

  p_value <- fread(paste("D:/R/Hot_tumor_and_cold_tumor/immune data1/cluster_gene/p_value.csv",sep=''),header = T)
  up_down <- fread(paste("D:/R/Hot_tumor_and_cold_tumor/immune data1/cluster_gene/updown.csv",sep=''),header = T)
  p_value <- as.data.frame(p_value)
  up_down <- as.data.frame(up_down)
  
  
  rownames(p_value) <- p_value$V1
  p_value <- p_value[,-1]
  
  rownames(up_down) <- up_down[,1]
  up_down <- up_down[,-1]
  
  up_down <- up_down[match(rownames(p_value),rownames(up_down)),]
  
  for(i in 1:length(rownames(up_down))){
    for(t in 1:length(colnames(up_down))){
      p_value[i,t] <- as.numeric(p_value[i,t])*as.numeric(up_down[i,t])
    }
  }
  
  
  DEG_BLCA_up <-  rownames(p_value)[which((p_value$BLCA>0) & (p_value$BLCA<0.05) )]
  DEG_BLCA_down <-  rownames(p_value)[which((p_value$BLCA<0) & (p_value$BLCA>-0.05) )]
  
  DEG_CESC_up <-  rownames(p_value)[which((p_value$CESC>0) & (p_value$CESC<0.05) )]
  DEG_CESC_down <-  rownames(p_value)[which((p_value$CESC<0) & (p_value$CESC>-0.05) )]
  
  DEG_PAAD_up <-  rownames(p_value)[which((p_value$PAAD>0) & (p_value$PAAD<0.05) )]
  DEG_PAAD_down <-  rownames(p_value)[which((p_value$PAAD<0) & (p_value$PAAD>-0.05) )]
  
  DEG_SARC_up <-  rownames(p_value)[which((p_value$SARC>0) & (p_value$SARC<0.05) )]
  DEG_SARC_down <-  rownames(p_value)[which((p_value$SARC<0) & (p_value$SARC>-0.05) )]
  
  DEG_SKCM_up <-  rownames(p_value)[which((p_value$SKCM>0) & (p_value$SKCM<0.05) )]
  DEG_SKCM_down <-  rownames(p_value)[which((p_value$SKCM<0) & (p_value$SKCM>-0.05) )]
  


dat <- list("BLCA"=DEG_BLCA_up,
            "CESC"=DEG_CESC_up,
            "PAAD"=DEG_PAAD_up,
            "SARC"=DEG_SARC_up,
            "SKCM"=DEG_SKCM_up
)


library(UpSetR)         #Upset图（upset 包，适用样本数 2-7）
library(VennDiagram)

pdf("DEGs_up.pdf",width=10,height=6)
upset(fromList(dat),  # fromList一个函数，用于将列表转换为与UpSetR兼容的数据形式。
      nsets = 100,     # 绘制的最大集合个数
      sets=rev(c("BLCA","CESC","PAAD","SARC","SKCM")),
      sets.bar.color=rev(c("#6AB3F4","#D459A2","#3A8F93","#FACC13","#5A7FF8")),
      #c("#6AB3F4","#D459A2","#3A8F93","#FACC13","#5A7FF8","#3A3AC8","#8747E0","#12C2C4","#FF7B78","#B27EEC","#CE1222")
      nintersects = 40, #绘制的最大交集个数，NA则全部绘制
      order.by = "freq", # 矩阵中的交点是如何排列的。 "freq"根据交集个数排序，"degree"根据
      keep.order = T, # 保持设置与使用sets参数输入的顺序一致。默认值是FALSE，它根据集合的大小排序。
      mb.ratio = c(0.6,0.4),   # 左侧和上方条形图的比例关系
      text.scale = 2, # 文字标签的大小
      mainbar.y.label = "Count of Intersection",  # y 标题
      sets.x.label = "Datasets Size",  # x 标题
)
dev.off()

dat <- list("BLCA"=DEG_BLCA_down,
            "CESC"=DEG_CESC_down,
            "PAAD"=DEG_PAAD_down,
            "SARC"=DEG_SARC_down,
            "SKCM"=DEG_SKCM_down
)


library(UpSetR)         #Upset图（upset 包，适用样本数 2-7）
library(VennDiagram)

pdf("DEGs_down.pdf",width=10,height=6)
upset(fromList(dat),  # fromList一个函数，用于将列表转换为与UpSetR兼容的数据形式。
      nsets = 100,     # 绘制的最大集合个数
      sets=rev(c("BLCA","CESC","PAAD","SARC","SKCM")),
      sets.bar.color=rev(c("#6AB3F4","#D459A2","#3A8F93","#FACC13","#5A7FF8")),
      #c("#6AB3F4","#D459A2","#3A8F93","#FACC13","#5A7FF8","#3A3AC8","#8747E0","#12C2C4","#FF7B78","#B27EEC","#CE1222")
      nintersects = 40, #绘制的最大交集个数，NA则全部绘制
      order.by = "freq", # 矩阵中的交点是如何排列的。 "freq"根据交集个数排序，"degree"根据
      keep.order = T, # 保持设置与使用sets参数输入的顺序一致。默认值是FALSE，它根据集合的大小排序。
      mb.ratio = c(0.6,0.4),   # 左侧和上方条形图的比例关系
      text.scale = 2, # 文字标签的大小
      mainbar.y.label = "Count of Intersection",  # y 标题
      sets.x.label = "Datasets Size",  # x 标题
)
dev.off()



dat <- list("BLCA_up"=DEG_BLCA_up,
            "CESC_up"=DEG_CESC_up,
            "PAAD_up"=DEG_PAAD_up,
            "SARC_up"=DEG_SARC_up,
            "SKCM_up"=DEG_SKCM_up,
            "BLCA_down"=DEG_BLCA_down,
            "CESC_down"=DEG_CESC_down,
            "PAAD_down"=DEG_PAAD_down,
            "SARC_down"=DEG_SARC_down,
            "SKCM_down"=DEG_SKCM_down
)


library(UpSetR)         #Upset图（upset 包，适用样本数 2-7）
library(VennDiagram)

pdf("DEGs.pdf",width=16,height=8)
upset(fromList(dat),  # fromList一个函数，用于将列表转换为与UpSetR兼容的数据形式。
      nsets = 100,     # 绘制的最大集合个数
      sets=rev(c("BLCA_up","CESC_up","PAAD_up","SARC_up","SKCM_up","BLCA_down","CESC_down","PAAD_down","SARC_down","SKCM_down")),
      sets.bar.color=rev(c("#6AB3F4","#D459A2","#3A8F93","#FACC13","#5A7FF8","#6AB3F4","#D459A2","#3A8F93","#FACC13","#5A7FF8")),
      #c("#6AB3F4","#D459A2","#3A8F93","#FACC13","#5A7FF8","#3A3AC8","#8747E0","#12C2C4","#FF7B78","#B27EEC","#CE1222")
      nintersects = 60, #绘制的最大交集个数，NA则全部绘制
      order.by = "freq", # 矩阵中的交点是如何排列的。 "freq"根据交集个数排序，"degree"根据
      keep.order = T, # 保持设置与使用sets参数输入的顺序一致。默认值是FALSE，它根据集合的大小排序。
      mb.ratio = c(0.6,0.4),   # 左侧和上方条形图的比例关系
      text.scale = 2, # 文字标签的大小
      mainbar.y.label = "Count of Intersection",  # y 标题
      sets.x.label = "Datasets Size",  # x 标题
)
dev.off()


names(which(table(c(DEG_BLCA_up,DEG_CESC_up,DEG_PAAD_up,DEG_SARC_up,DEG_SKCM_up))==5))
#"BTLA CCL19 CCL4 CCR5 CD160 CD27
# CD40LG CD48 CD96 CXCR3 CXCR6 HLA-DOA
# HLA-DOB ICOS KLRC1 LAG3 LTA PDCD1
# SLAMF7 TIGIT TNFRSF17 XCL2 

#"BTLA","CCL19","CCL4","CCR5","CD160","CD27","CD40LG","CD48","CD96","CXCR3",
#CXCR6","HLA-DOA","HLA-DOB","ICOS","KLRC1","LAG3","LTA","PDCD1","SLAMF7","TIGIT","TNFRSF17","XCL2"

setdiff(names(which(table(c(DEG_CESC_up,DEG_PAAD_up,DEG_SARC_up,DEG_SKCM_up))==4)),DEG_BLCA_up)
# CCL5 CCR2 CCR4 CCR7 CXCL13 GZMA HLA-DPB1 IL2 PRF1 TNFSF14 XCR1

names(which(table(c(DEG_BLCA_down,DEG_CESC_up,DEG_PAAD_up,DEG_SARC_up,DEG_SKCM_up))==5))
#XCR1


setdiff(names(which(table(c(DEG_BLCA_up,DEG_PAAD_up,DEG_SARC_up,DEG_SKCM_up))==4)),DEG_CESC_up)
#NA

setdiff(names(which(table(c(DEG_BLCA_up,DEG_CESC_up,DEG_SARC_up,DEG_SKCM_up))==4)),DEG_PAAD_up)
#B2M BTN3A1 BTN3A2 CXCL9 HLA-A HLA-B HLA-C HLA-DQA1 HLA-DQB1 HLA-DRA HLA-E HLA-F IFNG TAP1 TAP2 TNFRSF14 TNFRSF9"    

setdiff(names(which(table(c(DEG_BLCA_up,DEG_CESC_up,DEG_PAAD_up,DEG_SKCM_up))==4)),DEG_SARC_up)
# KLRK1 TNFRSF13B

setdiff(names(which(table(c(DEG_BLCA_up,DEG_CESC_up,DEG_PAAD_up,DEG_SARC_up))==4)),DEG_SKCM_up)
#CTLA4 CXCR5 TNFRSF13C


names(which(table(c(DEG_BLCA_down,DEG_CESC_down,DEG_PAAD_down,DEG_SARC_down,DEG_SKCM_down))==5))
#"CD276"

setdiff(names(which(table(c(DEG_CESC_down,DEG_PAAD_down,DEG_SARC_down,DEG_SKCM_down))==4)),DEG_BLCA_down)
"NT5E"


