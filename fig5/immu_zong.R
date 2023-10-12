setwd("D:/R/Hot_tumor_and_cold_tumor/immune data1")


Inhibitory <- c("ADORA2A","ARG1","BTLA","CD274","CD276","CTLA4","EDNRB","HAVCR2","IDO1","IL10",
                "IL13","IL4","KIR2DL1","KIR2DL3","LAG3","PDCD1","SLAMF7","TGFB1","TIGIT","VEGFA",
                "VEGFB","C10orf54","VTCN1","IL12A") 

Stimulaotry <- c("GZMA","BTN3A1","BTN3A2","CCL5","CD27","CD28","CD40","CD40LG","CD70","CD80",
                 "CX3CL1","CXCL10","CXCL9","ENTPD1","HMGB1","ICAM1","ICOS","ICOSLG","IFNA1",
                 "IFNA2","IFNG","IL1A","IL1B","IL2","IL2RA","ITGB2","PRF1","SELP","TLR4","TNF",
                 "TNFRSF14","TNFRSF18","TNFRSF4","TNFRSF9","TNFSF4","TNFSF9" ) 

chemokine <- c("CCL1","CCL2","CCL3","CCL4","CCL5","CCL7","CCL8","CCL11","CCL13","CCL14","CCL15",
               "CCL16","CCL17","CCL18","CCL19","CCL20","CCL21","CCL22","CCL23","CCL24","CCL25",
               "CCL26","CCL27","CCL28","CX3CL1","CXCL1","CXCL2","CXCL3","CXCL5","CXCL6","CXCL8",
               "CXCL9","CXCL10","CXCL11","CXCL12","CXCL13","CXCL14","CXCL16","CXCL17","XCL1","XCL2")

Immunoinhibitor <- c("ADORA2A","BTLA","CD160","CD244","CD274","CD96","CSF1R","CTLA4","HAVCR2",
                     "IDO1","IL10","IL10RB","KDR","KIR2DL1","KIR2DL3","LAG3","LGALS9","PDCD1",
                     "PDCD1LG2","PVRL2","TGFB1","TGFBR1","TIGIT","VTCN1")


Immunostimulator <- c( "BTNL2","C10orf54","CD27","CD276","CD28","CD40","CD40LG","CD48","CD70",     
                       "CD80","CD86","CXCL12","CXCR4","ENTPD1","HHLA2","ICOS","ICOSLG","IL2RA",    
                       "IL6","IL6R","KLRC1","KLRK1","LTA","MICB","NT5E", "PVR","RAET1E",   
                       "TMEM173","TMIGD2","TNFRSF13B","TNFRSF13C","TNFRSF14","TNFRSF17","TNFRSF18","TNFRSF25","TNFRSF4",  
                       "TNFRSF8","TNFRSF9","TNFSF13","TNFSF13B","TNFSF14","TNFSF15","TNFSF18","TNFSF4","TNFSF9",   
                       "ULBP1")

MHC <- c("B2M","HLA-A","HLA-B","HLA-C","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB","HLA-DPA1","HLA-DPB1",
         "HLA-DQA1","HLA-DQA2","HLA-DQB1","HLA-DRA","HLA-DRB1","HLA-E","HLA-F","HLA-G","TAP1","TAP2",    
         "TAPBP")

receptor <- c("CCR1","CCR2","CCR3","CCR4","CCR5","CCR6","CCR7","CCR8","CCR9","CCR10","CXCR1","CXCR2","CXCR3",
              "CXCR4","CXCR5","CXCR6","XCR1","CX3CR1")

allgene0 <- c(Inhibitory,Stimulaotry,chemokine,Immunoinhibitor,Immunostimulator,MHC,receptor)
allgene0 <- unique(allgene0)
allgene <- allgene0


allgene[22] <- "VSIR"
#"C10orf54"
allgene[106] <- "NECTIN2"
#"PVRL2"
allgene[122] <- "STING1"
#"TMEM173"


#cancer <- c("SARC")

cancers <- c("BLCA","CESC","PAAD","SARC","SKCM")
for(cancer in cancers){
  
  library(data.table)
  dat <-  fread(paste("D:/R/Hot_tumor_and_cold_tumor/count_TPM/",cancer,"_convert_exp.txt",sep=''))
  dat <- as.data.frame(dat)
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]
  
  
  data <- dat[match(allgene,rownames(dat)),]
  
  rownames(data) <- allgene0
  
  write.csv(data,file=paste("immu_",cancer,"_zong.csv",sep=''),quote=F)
  
  
}



which(rownames(dat)=="VSIR")

which(rownames(dat)=="NECTIN2")

which(rownames(dat)=="STING1")

