library(limma)
library(edgeR)
library(biomaRt)
library(PLIER)
library(ggplot2)
library(RColorBrewer)
library(qvalue)
library(genefilter)
library(gridExtra)
library(dplyr)
library(data.table)
library(ternvis)
library(ggpubr)
library(cowplot)
library(plyr)
library(stringr)
library(UpSetR)
library(foreign)
library(MASS)
library(sfsmisc)
library(IHW)
library("org.Rn.eg.db")
library(MCL)
library(corrplot)
library(ggrepel)
library(CellCODE)
library(tidyverse)
library(ggsci)
library(circlize)
library(grid)
library(GenomicRanges)

setwd("C:/Users/gsmit/Documents/Mount Sinai/Sealfon Laboratory/MoTrPAC/PASS1B Raw Data")

load("Figure3D_3G_S15_112624.RData")
allmethmotifs <- read.table(file = "allmethoutput50.txt",header = T,sep = "\t")
load("enstosym.RData")

####
# Next we want to look for distal relationships between DEGs and DMRs

# Let's start with gastroc. We have a 166x730 matrix of DMRs vs DEGs

####
# gastro

gastromethcortrim <- gastromethcor*(gastromethdistance < 500000 & gastromethdistance > 0)
gastromethcortestmod <- gastromethcortest*(gastromethdistance > 0 & gastromethdistance < 500000)
gastromethcortestmod[gastromethcortestmod == 0] <- 1
gastroconnectsigmeth <- rownames(gastromethcortestmod[apply(abs(gastromethcortestmod),1,min) < 0.05,])
gastroconnectsiggenes <- colnames(gastromethcortestmod[,apply(abs(gastromethcortestmod),2,min) < 0.05])
gastromethcorfin <- gastromethcortrim[gastroconnectsigmeth,gastroconnectsiggenes]

gastromethcortrimtfs <- allmethmotifs[allmethmotifs$PositionID %in% rownames(gastromethcorfin),]
gastromethcortfs <- matrix(0L,nrow = dim(gastromethcorfin)[1],ncol = length(unique(allmethmotifs[allmethmotifs$PositionID %in% rownames(gastromethcorfin),"Motif.Name"])))
rownames(gastromethcortfs) <- rownames(gastromethcorfin)
colnames(gastromethcortfs) <- unique(allmethmotifs[allmethmotifs$PositionID %in% rownames(gastromethcorfin),"Motif.Name"])
for(i in 1:dim(gastromethcortrimtfs)[1]){
  gastromethcortfs[gastromethcortrimtfs$PositionID[i],gastromethcortrimtfs$Motif.Name[i]] <- 1
}
gastromethcortfsfintrim <- gastromethcortfs[apply(gastromethcortfs,1,max) > 0,]
gastromethcortfsfin <- gastromethcortfsfintrim[,apply(gastromethcortfsfintrim,2,max) > 0]

gastromethcorfintrim <- gastromethcorfin[rownames(gastromethcortfsfin),]
gastromethcorfinal <- gastromethcorfintrim[,apply(abs(gastromethcorfintrim),2,max) > 0]

#png(filename = "gastro_DMR_DEG_correlated_neighbors.png",width = 5,height = 5,units = "in",res = 600)
#pheatmap(gastromethcorfinal,labels_col = enstosym[colnames(gastromethcorfinal),"Symbol"],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"))
#dev.off()

####
# heart

heartmethcortrim <- heartmethcor*(heartmethdistance < 500000 & heartmethdistance > 0)
heartmethcortestmod <- heartmethcortest*(heartmethdistance > 0 & heartmethdistance < 500000)
heartmethcortestmod[heartmethcortestmod == 0] <- 1
heartconnectsigmeth <- rownames(heartmethcortestmod[apply(abs(heartmethcortestmod),1,min) < 0.05,])
heartconnectsiggenes <- colnames(heartmethcortestmod[,apply(abs(heartmethcortestmod),2,min) < 0.05])
heartmethcorfin <- heartmethcortrim[heartconnectsigmeth,heartconnectsiggenes]

heartmethcortrimtfs <- allmethmotifs[allmethmotifs$PositionID %in% rownames(heartmethcorfin),]
heartmethcortfs <- matrix(0L,nrow = dim(heartmethcorfin)[1],ncol = length(unique(allmethmotifs[allmethmotifs$PositionID %in% rownames(heartmethcorfin),"Motif.Name"])))
rownames(heartmethcortfs) <- rownames(heartmethcorfin)
colnames(heartmethcortfs) <- unique(allmethmotifs[allmethmotifs$PositionID %in% rownames(heartmethcorfin),"Motif.Name"])
for(i in 1:dim(heartmethcortrimtfs)[1]){
  heartmethcortfs[heartmethcortrimtfs$PositionID[i],heartmethcortrimtfs$Motif.Name[i]] <- 1
}
heartmethcortfsfintrim <- heartmethcortfs[apply(heartmethcortfs,1,max) > 0,]
heartmethcortfsfin <- heartmethcortfsfintrim[,apply(heartmethcortfsfintrim,2,max) > 0]

heartmethcorfintrim <- heartmethcorfin[rownames(heartmethcortfsfin),]
heartmethcorfinal <- heartmethcorfintrim[,apply(abs(heartmethcorfintrim),2,max) > 0]

#png(filename = "heart_DMR_DEG_correlated_neighbors.png",width = 5,height = 5,units = "in",res = 600)
#pheatmap(heartmethcorfinal,labels_col = enstosym[colnames(heartmethcorfinal),"Symbol"],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"))
#dev.off()


####
# hippo

hippomethcortrim <- hippomethcor*(hippomethdistance < 500000 & hippomethdistance > 0)
hippomethcortestmod <- hippomethcortest*(hippomethdistance > 0 & hippomethdistance < 500000)
hippomethcortestmod[hippomethcortestmod == 0] <- 1
hippoconnectsigmeth <- rownames(hippomethcortestmod[apply(abs(hippomethcortestmod),1,min) < 0.05,])
hippoconnectsiggenes <- colnames(hippomethcortestmod[,apply(abs(hippomethcortestmod),2,min) < 0.05])
hippomethcorfin <- hippomethcortrim[hippoconnectsigmeth,hippoconnectsiggenes]

hippomethcortrimtfs <- allmethmotifs[allmethmotifs$PositionID %in% rownames(hippomethcorfin),]
hippomethcortfs <- matrix(0L,nrow = dim(hippomethcorfin)[1],ncol = length(unique(allmethmotifs[allmethmotifs$PositionID %in% rownames(hippomethcorfin),"Motif.Name"])))
rownames(hippomethcortfs) <- rownames(hippomethcorfin)
colnames(hippomethcortfs) <- unique(allmethmotifs[allmethmotifs$PositionID %in% rownames(hippomethcorfin),"Motif.Name"])
for(i in 1:dim(hippomethcortrimtfs)[1]){
  hippomethcortfs[hippomethcortrimtfs$PositionID[i],hippomethcortrimtfs$Motif.Name[i]] <- 1
}
hippomethcortfsfintrim <- hippomethcortfs[apply(hippomethcortfs,1,max) > 0,]
hippomethcortfsfin <- hippomethcortfsfintrim[,apply(hippomethcortfsfintrim,2,max) > 0]

hippomethcorfintrim <- hippomethcorfin[rownames(hippomethcortfsfin),]
hippomethcorfinal <- hippomethcorfintrim[,apply(abs(hippomethcorfintrim),2,max) > 0]

#png(filename = "hippo_DMR_DEG_correlated_neighbors.png",width = 5,height = 5,units = "in",res = 600)
#pheatmap(hippomethcorfinal,labels_col = enstosym[colnames(hippomethcorfinal),"Symbol"],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"))
#dev.off()


####
# kidney

kidneymethcortrim <- kidneymethcor*(kidneymethdistance < 500000 & kidneymethdistance > 0)
kidneymethcortestmod <- kidneymethcortest*(kidneymethdistance > 0 & kidneymethdistance < 500000)
kidneymethcortestmod[kidneymethcortestmod == 0] <- 1
kidneyconnectsigmeth <- rownames(kidneymethcortestmod[apply(abs(kidneymethcortestmod),1,min) < 0.05,])
kidneyconnectsiggenes <- colnames(kidneymethcortestmod[,apply(abs(kidneymethcortestmod),2,min) < 0.05])
kidneymethcorfin <- kidneymethcortrim[kidneyconnectsigmeth,kidneyconnectsiggenes]

kidneymethcortrimtfs <- allmethmotifs[allmethmotifs$PositionID %in% rownames(kidneymethcorfin),]
kidneymethcortfs <- matrix(0L,nrow = dim(kidneymethcorfin)[1],ncol = length(unique(allmethmotifs[allmethmotifs$PositionID %in% rownames(kidneymethcorfin),"Motif.Name"])))
rownames(kidneymethcortfs) <- rownames(kidneymethcorfin)
colnames(kidneymethcortfs) <- unique(allmethmotifs[allmethmotifs$PositionID %in% rownames(kidneymethcorfin),"Motif.Name"])
for(i in 1:dim(kidneymethcortrimtfs)[1]){
  kidneymethcortfs[kidneymethcortrimtfs$PositionID[i],kidneymethcortrimtfs$Motif.Name[i]] <- 1
}
kidneymethcortfsfintrim <- kidneymethcortfs[apply(kidneymethcortfs,1,max) > 0,]
kidneymethcortfsfin <- kidneymethcortfsfintrim[,apply(kidneymethcortfsfintrim,2,max) > 0]

kidneymethcorfintrim <- kidneymethcorfin[rownames(kidneymethcortfsfin),]
kidneymethcorfinal <- kidneymethcorfintrim[,apply(abs(kidneymethcorfintrim),2,max) > 0]

#png(filename = "kidney_DMR_DEG_correlated_neighbors.png",width = 5,height = 5,units = "in",res = 600)
#pheatmap(kidneymethcorfinal,labels_col = enstosym[colnames(kidneymethcorfinal),"Symbol"],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"))
#dev.off()


####
# liver

livermethcortrim <- livermethcor*(livermethdistance < 500000 & livermethdistance > 0)
livermethcortestmod <- livermethcortest*(livermethdistance > 0 & livermethdistance < 500000)
livermethcortestmod[livermethcortestmod == 0] <- 1
liverconnectsigmeth <- rownames(livermethcortestmod[apply(abs(livermethcortestmod),1,min) < 0.05,])
liverconnectsiggenes <- colnames(livermethcortestmod[,apply(abs(livermethcortestmod),2,min) < 0.05])
livermethcorfin <- livermethcortrim[liverconnectsigmeth,liverconnectsiggenes]

livermethcortrimtfs <- allmethmotifs[allmethmotifs$PositionID %in% rownames(livermethcorfin),]
livermethcortfs <- matrix(0L,nrow = dim(livermethcorfin)[1],ncol = length(unique(allmethmotifs[allmethmotifs$PositionID %in% rownames(livermethcorfin),"Motif.Name"])))
rownames(livermethcortfs) <- rownames(livermethcorfin)
colnames(livermethcortfs) <- unique(allmethmotifs[allmethmotifs$PositionID %in% rownames(livermethcorfin),"Motif.Name"])
for(i in 1:dim(livermethcortrimtfs)[1]){
  livermethcortfs[livermethcortrimtfs$PositionID[i],livermethcortrimtfs$Motif.Name[i]] <- 1
}
livermethcortfsfintrim <- livermethcortfs[apply(livermethcortfs,1,max) > 0,]
livermethcortfsfin <- livermethcortfsfintrim[,apply(livermethcortfsfintrim,2,max) > 0]

livermethcorfintrim <- livermethcorfin[rownames(livermethcortfsfin),]
livermethcorfinal <- livermethcorfintrim[,apply(abs(livermethcorfintrim),2,max) > 0]

#png(filename = "liver_DMR_DEG_correlated_neighbors.png",width = 5,height = 5,units = "in",res = 600)
#pheatmap(livermethcorfinal,labels_col = enstosym[colnames(livermethcorfinal),"Symbol"],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"))
#dev.off()


####
# lung

lungmethcortrim <- lungmethcor*(lungmethdistance < 500000 & lungmethdistance > 0)
lungmethcortestmod <- lungmethcortest*(lungmethdistance > 0 & lungmethdistance < 500000)
lungmethcortestmod[lungmethcortestmod == 0] <- 1
lungconnectsigmeth <- rownames(lungmethcortestmod[apply(abs(lungmethcortestmod),1,min) < 0.05,])
lungconnectsiggenes <- colnames(lungmethcortestmod[,apply(abs(lungmethcortestmod),2,min) < 0.05])
lungmethcorfin <- lungmethcortrim[lungconnectsigmeth,lungconnectsiggenes]

lungmethcortrimtfs <- allmethmotifs[allmethmotifs$PositionID %in% rownames(lungmethcorfin),]
lungmethcortfs <- matrix(0L,nrow = dim(lungmethcorfin)[1],ncol = length(unique(allmethmotifs[allmethmotifs$PositionID %in% rownames(lungmethcorfin),"Motif.Name"])))
rownames(lungmethcortfs) <- rownames(lungmethcorfin)
colnames(lungmethcortfs) <- unique(allmethmotifs[allmethmotifs$PositionID %in% rownames(lungmethcorfin),"Motif.Name"])
for(i in 1:dim(lungmethcortrimtfs)[1]){
  lungmethcortfs[lungmethcortrimtfs$PositionID[i],lungmethcortrimtfs$Motif.Name[i]] <- 1
}
lungmethcortfsfintrim <- lungmethcortfs[apply(lungmethcortfs,1,max) > 0,]
lungmethcortfsfin <- lungmethcortfsfintrim[,apply(lungmethcortfsfintrim,2,max) > 0]

lungmethcorfintrim <- lungmethcorfin[rownames(lungmethcortfsfin),]
lungmethcorfinal <- lungmethcorfintrim[,apply(abs(lungmethcorfintrim),2,max) > 0]

#png(filename = "lung_DMR_DEG_correlated_neighbors.png",width = 5,height = 5,units = "in",res = 600)
#pheatmap(lungmethcorfinal,labels_col = enstosym[colnames(lungmethcorfinal),"Symbol"],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"))
#dev.off()


####
# brown

brownmethcortrim <- brownmethcor*(brownmethdistance < 500000 & brownmethdistance > 0)
brownmethcortestmod <- brownmethcortest*(brownmethdistance > 0 & brownmethdistance < 500000)
brownmethcortestmod[brownmethcortestmod == 0] <- 1
brownconnectsigmeth <- rownames(brownmethcortestmod[apply(abs(brownmethcortestmod),1,min) < 0.05,])
brownconnectsiggenes <- colnames(brownmethcortestmod[,apply(abs(brownmethcortestmod),2,min) < 0.05])
brownmethcorfin <- brownmethcortrim[brownconnectsigmeth,brownconnectsiggenes]

brownmethcortrimtfs <- allmethmotifs[allmethmotifs$PositionID %in% rownames(brownmethcorfin),]
brownmethcortfs <- matrix(0L,nrow = dim(brownmethcorfin)[1],ncol = length(unique(allmethmotifs[allmethmotifs$PositionID %in% rownames(brownmethcorfin),"Motif.Name"])))
rownames(brownmethcortfs) <- rownames(brownmethcorfin)
colnames(brownmethcortfs) <- unique(allmethmotifs[allmethmotifs$PositionID %in% rownames(brownmethcorfin),"Motif.Name"])
for(i in 1:dim(brownmethcortrimtfs)[1]){
  brownmethcortfs[brownmethcortrimtfs$PositionID[i],brownmethcortrimtfs$Motif.Name[i]] <- 1
}
brownmethcortfsfintrim <- brownmethcortfs[apply(brownmethcortfs,1,max) > 0,]
brownmethcortfsfin <- brownmethcortfsfintrim[,apply(brownmethcortfsfintrim,2,max) > 0]

brownmethcorfintrim <- brownmethcorfin[rownames(brownmethcortfsfin),]
brownmethcorfinal <- brownmethcorfintrim[,apply(abs(brownmethcorfintrim),2,max) > 0]

#png(filename = "brown_DMR_DEG_correlated_neighbors.png",width = 5,height = 5,units = "in",res = 600)
#pheatmap(brownmethcorfinal,labels_col = enstosym[colnames(brownmethcorfinal),"Symbol"],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"))
#dev.off()


####
# white

whitemethcortrim <- whitemethcor*(whitemethdistance < 500000 & whitemethdistance > 0)
whitemethcortestmod <- whitemethcortest*(whitemethdistance > 0 & whitemethdistance < 500000)
whitemethcortestmod[whitemethcortestmod == 0] <- 1
whiteconnectsigmeth <- rownames(whitemethcortestmod[apply(abs(whitemethcortestmod),1,min) < 0.05,])
whiteconnectsiggenes <- colnames(whitemethcortestmod[,apply(abs(whitemethcortestmod),2,min) < 0.05])
whitemethcorfin <- whitemethcortrim[whiteconnectsigmeth,whiteconnectsiggenes]

whitemethcortrimtfs <- allmethmotifs[allmethmotifs$PositionID %in% rownames(whitemethcorfin),]
whitemethcortfs <- matrix(0L,nrow = dim(whitemethcorfin)[1],ncol = length(unique(allmethmotifs[allmethmotifs$PositionID %in% rownames(whitemethcorfin),"Motif.Name"])))
rownames(whitemethcortfs) <- rownames(whitemethcorfin)
colnames(whitemethcortfs) <- unique(allmethmotifs[allmethmotifs$PositionID %in% rownames(whitemethcorfin),"Motif.Name"])
for(i in 1:dim(whitemethcortrimtfs)[1]){
  whitemethcortfs[whitemethcortrimtfs$PositionID[i],whitemethcortrimtfs$Motif.Name[i]] <- 1
}
whitemethcortfsfintrim <- whitemethcortfs[apply(whitemethcortfs,1,max) > 0,]
whitemethcortfsfin <- whitemethcortfsfintrim[,apply(whitemethcortfsfintrim,2,max) > 0]

whitemethcorfintrim <- whitemethcorfin[rownames(whitemethcortfsfin),]
whitemethcorfinal <- whitemethcorfintrim[,apply(abs(whitemethcorfintrim),2,max) > 0]

#png(filename = "white_DMR_DEG_correlated_neighbors.png",width = 5,height = 5,units = "in",res = 600)
#pheatmap(whitemethcorfinal,labels_col = enstosym[colnames(whitemethcorfinal),"Symbol"],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"))
#dev.off()


###############
# Let's do TF significance- start with gastro

####
# gastro

gastrocortfsig <- matrix(0L,nrow = dim(gastromethcortfsfin)[2],ncol = 5)
rownames(gastrocortfsig) <- colnames(gastromethcortfsfin)
colnames(gastrocortfsig) <- c("Correlated_Target_Count","Total_Target_Count","Fraction_Corr","Fraction_Total","Enrichment_Significance")
gastrototcor <- length(rownames(gastromethcortfsfin))
gastrototmeth <- length(unique(allmethmotifs$PositionID))
for(i in 1:dim(gastrocortfsig)[1]){
  if(i%%20 == 0){
    print(i)
  }
  ourtf <- rownames(gastrocortfsig)[i]
  gastrocortfsig[i,"Correlated_Target_Count"] <- sum(gastromethcortfsfin[,i])
  gastrocortfsig[i,"Total_Target_Count"] <- length(unique(allmethmotifs[allmethmotifs$Motif.Name %in% ourtf,"PositionID"]))
  gastrocortfsig[i,"Fraction_Corr"] <- gastrocortfsig[i,"Correlated_Target_Count"]/gastrototcor
  gastrocortfsig[i,"Fraction_Total"] <- gastrocortfsig[i,"Total_Target_Count"]/gastrototmeth
  gastrocortfsig[i,"Enrichment_Significance"] <- binom.test(x = gastrocortfsig[i,"Correlated_Target_Count"], n = gastrototcor, p = gastrocortfsig[i,"Total_Target_Count"]/gastrototmeth,alternative = "greater")$p.value
}

####
# heart

heartcortfsig <- matrix(0L,nrow = dim(heartmethcortfsfin)[2],ncol = 5)
rownames(heartcortfsig) <- colnames(heartmethcortfsfin)
colnames(heartcortfsig) <- c("Correlated_Target_Count","Total_Target_Count","Fraction_Corr","Fraction_Total","Enrichment_Significance")
hearttotcor <- length(rownames(heartmethcortfsfin))
hearttotmeth <- length(unique(allmethmotifs$PositionID))
for(i in 1:dim(heartcortfsig)[1]){
  if(i%%20 == 0){
    print(i)
  }
  ourtf <- rownames(heartcortfsig)[i]
  heartcortfsig[i,"Correlated_Target_Count"] <- sum(heartmethcortfsfin[,i])
  heartcortfsig[i,"Total_Target_Count"] <- length(unique(allmethmotifs[allmethmotifs$Motif.Name %in% ourtf,"PositionID"]))
  heartcortfsig[i,"Fraction_Corr"] <- heartcortfsig[i,"Correlated_Target_Count"]/hearttotcor
  heartcortfsig[i,"Fraction_Total"] <- heartcortfsig[i,"Total_Target_Count"]/hearttotmeth
  heartcortfsig[i,"Enrichment_Significance"] <- binom.test(x = heartcortfsig[i,"Correlated_Target_Count"], n = hearttotcor, p = heartcortfsig[i,"Total_Target_Count"]/hearttotmeth,alternative = "greater")$p.value
}


####
# liver

livercortfsig <- matrix(0L,nrow = dim(livermethcortfsfin)[2],ncol = 5)
rownames(livercortfsig) <- colnames(livermethcortfsfin)
colnames(livercortfsig) <- c("Correlated_Target_Count","Total_Target_Count","Fraction_Corr","Fraction_Total","Enrichment_Significance")
livertotcor <- length(rownames(livermethcortfsfin))
livertotmeth <- length(unique(allmethmotifs$PositionID))
for(i in 1:dim(livercortfsig)[1]){
  if(i%%20 == 0){
    print(i)
  }
  ourtf <- rownames(livercortfsig)[i]
  livercortfsig[i,"Correlated_Target_Count"] <- sum(livermethcortfsfin[,i])
  livercortfsig[i,"Total_Target_Count"] <- length(unique(allmethmotifs[allmethmotifs$Motif.Name %in% ourtf,"PositionID"]))
  livercortfsig[i,"Fraction_Corr"] <- livercortfsig[i,"Correlated_Target_Count"]/livertotcor
  livercortfsig[i,"Fraction_Total"] <- livercortfsig[i,"Total_Target_Count"]/livertotmeth
  livercortfsig[i,"Enrichment_Significance"] <- binom.test(x = livercortfsig[i,"Correlated_Target_Count"], n = livertotcor, p = livercortfsig[i,"Total_Target_Count"]/livertotmeth,alternative = "greater")$p.value
}

####
# lung

lungcortfsig <- matrix(0L,nrow = dim(lungmethcortfsfin)[2],ncol = 5)
rownames(lungcortfsig) <- colnames(lungmethcortfsfin)
colnames(lungcortfsig) <- c("Correlated_Target_Count","Total_Target_Count","Fraction_Corr","Fraction_Total","Enrichment_Significance")
lungtotcor <- length(rownames(lungmethcortfsfin))
lungtotmeth <- length(unique(allmethmotifs$PositionID))
for(i in 1:dim(lungcortfsig)[1]){
  if(i%%20 == 0){
    print(i)
  }
  ourtf <- rownames(lungcortfsig)[i]
  lungcortfsig[i,"Correlated_Target_Count"] <- sum(lungmethcortfsfin[,i])
  lungcortfsig[i,"Total_Target_Count"] <- length(unique(allmethmotifs[allmethmotifs$Motif.Name %in% ourtf,"PositionID"]))
  lungcortfsig[i,"Fraction_Corr"] <- lungcortfsig[i,"Correlated_Target_Count"]/lungtotcor
  lungcortfsig[i,"Fraction_Total"] <- lungcortfsig[i,"Total_Target_Count"]/lungtotmeth
  lungcortfsig[i,"Enrichment_Significance"] <- binom.test(x = lungcortfsig[i,"Correlated_Target_Count"], n = lungtotcor, p = lungcortfsig[i,"Total_Target_Count"]/lungtotmeth,alternative = "greater")$p.value
}

####
# brown

browncortfsig <- matrix(0L,nrow = dim(brownmethcortfsfin)[2],ncol = 5)
rownames(browncortfsig) <- colnames(brownmethcortfsfin)
colnames(browncortfsig) <- c("Correlated_Target_Count","Total_Target_Count","Fraction_Corr","Fraction_Total","Enrichment_Significance")
browntotcor <- length(rownames(brownmethcortfsfin))
browntotmeth <- length(unique(allmethmotifs$PositionID))
for(i in 1:dim(browncortfsig)[1]){
  if(i%%20 == 0){
    print(i)
  }
  ourtf <- rownames(browncortfsig)[i]
  browncortfsig[i,"Correlated_Target_Count"] <- sum(brownmethcortfsfin[,i])
  browncortfsig[i,"Total_Target_Count"] <- length(unique(allmethmotifs[allmethmotifs$Motif.Name %in% ourtf,"PositionID"]))
  browncortfsig[i,"Fraction_Corr"] <- browncortfsig[i,"Correlated_Target_Count"]/browntotcor
  browncortfsig[i,"Fraction_Total"] <- browncortfsig[i,"Total_Target_Count"]/browntotmeth
  browncortfsig[i,"Enrichment_Significance"] <- binom.test(x = browncortfsig[i,"Correlated_Target_Count"], n = browntotcor, p = browncortfsig[i,"Total_Target_Count"]/browntotmeth,alternative = "greater")$p.value
}

####
# white

whitecortfsig <- matrix(0L,nrow = dim(whitemethcortfsfin)[2],ncol = 5)
rownames(whitecortfsig) <- colnames(whitemethcortfsfin)
colnames(whitecortfsig) <- c("Correlated_Target_Count","Total_Target_Count","Fraction_Corr","Fraction_Total","Enrichment_Significance")
whitetotcor <- length(rownames(whitemethcortfsfin))
whitetotmeth <- length(unique(allmethmotifs$PositionID))
for(i in 1:dim(whitecortfsig)[1]){
  if(i%%20 == 0){
    print(i)
  }
  ourtf <- rownames(whitecortfsig)[i]
  whitecortfsig[i,"Correlated_Target_Count"] <- sum(whitemethcortfsfin[,i])
  whitecortfsig[i,"Total_Target_Count"] <- length(unique(allmethmotifs[allmethmotifs$Motif.Name %in% ourtf,"PositionID"]))
  whitecortfsig[i,"Fraction_Corr"] <- whitecortfsig[i,"Correlated_Target_Count"]/whitetotcor
  whitecortfsig[i,"Fraction_Total"] <- whitecortfsig[i,"Total_Target_Count"]/whitetotmeth
  whitecortfsig[i,"Enrichment_Significance"] <- binom.test(x = whitecortfsig[i,"Correlated_Target_Count"], n = whitetotcor, p = whitecortfsig[i,"Total_Target_Count"]/whitetotmeth,alternative = "greater")$p.value
}

#######
# Okay I have identified the significant enrichments in TFs - now let's examine specific interesting cases

# LUNG AP-2gamma(AP2)/MCF7-TFAP2C-ChIP-Seq(GSE21234)/Homer

ourtf <- "AP-2gamma(AP2)/MCF7-TFAP2C-ChIP-Seq(GSE21234)/Homer"
ourtftargetcormat <- lungmethcorfinal[rownames(lungmethcortfsfin[lungmethcortfsfin[,ourtf] > 0,]),]
ourtftargetcormat <- ourtftargetcormat[,apply(abs(ourtftargetcormat),2,max) > 0]
png(filename = "Supplemental Figure S19C_021325.png",width = 5,height = 3,units = "in",res = 600)
pheatmap(ourtftargetcormat,labels_col = enstosym[colnames(ourtftargetcormat),"Symbol"],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),main = paste("LUNG ",gsub("\\(.*","",ourtf)," DEG-DMR Targets",sep = ""))
dev.off()
pdf(file = "Supplemental Figure S19C_021325.pdf",width = 5,height = 3)
pheatmap(ourtftargetcormat,labels_col = enstosym[colnames(ourtftargetcormat),"Symbol"],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),main = paste("LUNG ",gsub("\\(.*","",ourtf)," DEG-DMR Targets",sep = ""))
dev.off()

# LUNG NF1-halfsite(CTF)/LNCaP-NF1-ChIP-Seq(Unpublished)/Homer

ourtf <- "NF1-halfsite(CTF)/LNCaP-NF1-ChIP-Seq(Unpublished)/Homer"
ourtftargetcormat <- lungmethcorfinal[rownames(lungmethcortfsfin[lungmethcortfsfin[,ourtf] > 0,]),]
ourtftargetcormat <- ourtftargetcormat[,apply(abs(ourtftargetcormat),2,max) > 0]
png(filename = "Supplemental Figure S19B_021325.png",width = 5,height = 3,units = "in",res = 600)
pheatmap(ourtftargetcormat,labels_col = enstosym[colnames(ourtftargetcormat),"Symbol"],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),main = paste("LUNG ",gsub("\\(.*","",ourtf)," DEG-DMR Targets",sep = ""))
dev.off()
pdf(file = "Supplemental Figure S19B_021325.pdf",width = 5,height = 3)
pheatmap(ourtftargetcormat,labels_col = enstosym[colnames(ourtftargetcormat),"Symbol"],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),main = paste("LUNG ",gsub("\\(.*","",ourtf)," DEG-DMR Targets",sep = ""))
dev.off()


# WAT-SC NF1-halfsite(CTF)/LNCaP-NF1-ChIP-Seq(Unpublished)/Homer

ourtf <- "NF1-halfsite(CTF)/LNCaP-NF1-ChIP-Seq(Unpublished)/Homer"
ourtftargetcormat <- whitemethcorfinal[rownames(whitemethcortfsfin[whitemethcortfsfin[,ourtf] > 0,]),]
ourtftargetcormat <- ourtftargetcormat[,apply(abs(ourtftargetcormat),2,max) > 0]
png(filename = "Supplemental Figure S19A_021325.png",width = 12,height = 5,units = "in",res = 600)
pheatmap(ourtftargetcormat,labels_col = enstosym[colnames(ourtftargetcormat),"Symbol"],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),main = paste("WAT-SC ",gsub("\\(.*","",ourtf)," DEG-DMR Targets",sep = ""))
dev.off()
pdf(file = "Supplemental Figure S19A_021325.pdf",width = 12,height = 5)
pheatmap(ourtftargetcormat,labels_col = enstosym[colnames(ourtftargetcormat),"Symbol"],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),main = paste("WAT-SC ",gsub("\\(.*","",ourtf)," DEG-DMR Targets",sep = ""))
dev.off()


# SKM-GN Sp2(Zf)/HEK293-Sp2.eGFP-ChIP-Seq(Encode)/Homer

ourtf <- "Sp2(Zf)/HEK293-Sp2.eGFP-ChIP-Seq(Encode)/Homer"
ourtftargetcormat <- gastromethcorfinal[rownames(gastromethcortfsfin[gastromethcortfsfin[,ourtf] > 0,]),]
ourtftargetcormat <- ourtftargetcormat[,apply(abs(ourtftargetcormat),2,max) > 0]
png(file = "Figure 4N_021325.png",width = 4,height = 3,units = "in",res = 600)
pheatmap(ourtftargetcormat,labels_col = enstosym[colnames(ourtftargetcormat),"Symbol"],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),main = paste("SKM-GN ",gsub("\\(.*","",ourtf)," DEG-DMR Targets",sep = ""))
dev.off()
pdf(file = "Figure 4N_021325.pdf",width = 4,height = 3)
pheatmap(ourtftargetcormat,labels_col = enstosym[colnames(ourtftargetcormat),"Symbol"],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),main = paste("SKM-GN ",gsub("\\(.*","",ourtf)," DEG-DMR Targets",sep = ""))
dev.off()



####
# Now to generate scatter plots highlighting the correlations between individual DEG-DMR pairs

pidmetadf <- data.frame(row.names = colnames(gastrornanorm),
                        "Group" = rep("",length(colnames(gastrornanorm))),
                        "Sex" = rep("",length(colnames(gastrornanorm))))
for(i in 1:length(colnames(gastrornanorm))){
  ourpid <- rownames(pidmetadf)[i]
  pidmetadf[i,"Group"] <- pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourpid,"pass1bf0011"][1]
  pidmetadf[i,"Sex"] <- c("Female","Male")[as.numeric(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourpid,"pass1bf0027"][1])]
}
pidmetadf[pidmetadf$Group %in% c("Eight-week program Control Group"),"Group"] <- "control"
pidmetadf[pidmetadf$Group %in% c("One-week program"),"Group"] <- "1w"
pidmetadf[pidmetadf$Group %in% c("Two-week program"),"Group"] <- "2w"
pidmetadf[pidmetadf$Group %in% c("Four-week program"),"Group"] <- "4w"
pidmetadf[pidmetadf$Group %in% c("Eight-week program Training Group"),"Group"] <- "8w"


ourtf <- "NF1-halfsite(CTF)/LNCaP-NF1-ChIP-Seq(Unpublished)/Homer"
ourtftargetcormat <- whitemethcorfinal[rownames(whitemethcortfsfin[whitemethcortfsfin[,ourtf] > 0,]),]
ourtftargetcormat <- ourtftargetcormat[,apply(abs(ourtftargetcormat),2,max) > 0]

oursite <- "chr2-375476_cluster1"
oursite_white_corgenes <- colnames(ourtftargetcormat)[abs(ourtftargetcormat[oursite,]) > 0.25]

ourgene <- oursite_white_corgenes[2]
ourdf <- data.frame("RNA.L2FC" = whiternatrainingl2fcbysubject[ourgene,],
                    "Meth.L2FC" = whitemethtrainingl2fcbysubject[oursite,],
                    "Group" = pidmetadf[colnames(whiternatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(whiternatrainingl2fcbysubject),"Sex"])
png(file = "Supplemental Figure S19G_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA.L2FC,y=Meth.L2FC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("WAT-SC ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$Meth.L2FC,ourdf$RNA.L2FC),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("DMR L2FC") + xlab("DEG L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()
pdf(file = "Supplemental Figure S19G_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA.L2FC,y=Meth.L2FC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("WAT-SC ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$Meth.L2FC,ourdf$RNA.L2FC),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("DMR L2FC") + xlab("DEG L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()

ourgene <- oursite_white_corgenes[1]
ourdf <- data.frame("RNA.L2FC" = whiternatrainingl2fcbysubject[ourgene,],
                    "Meth.L2FC" = whitemethtrainingl2fcbysubject[oursite,],
                    "Group" = pidmetadf[colnames(whiternatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(whiternatrainingl2fcbysubject),"Sex"])
png(file = "Supplemental Figure S19J_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA.L2FC,y=Meth.L2FC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("WAT-SC ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$Meth.L2FC,ourdf$RNA.L2FC),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("DMR L2FC") + xlab("DEG L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()
pdf(file = "Supplemental Figure S19J_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA.L2FC,y=Meth.L2FC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("WAT-SC ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$Meth.L2FC,ourdf$RNA.L2FC),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("DMR L2FC") + xlab("DEG L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()

ourgene <- oursite_white_corgenes[3]
ourdf <- data.frame("RNA.L2FC" = whiternatrainingl2fcbysubject[ourgene,],
                    "Meth.L2FC" = whitemethtrainingl2fcbysubject[oursite,],
                    "Group" = pidmetadf[colnames(whiternatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(whiternatrainingl2fcbysubject),"Sex"])
png(file = "Supplemental Figure S19H_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA.L2FC,y=Meth.L2FC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("WAT-SC ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$Meth.L2FC,ourdf$RNA.L2FC),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("DMR L2FC") + xlab("DEG L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()
pdf(file = "Supplemental Figure S19H_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA.L2FC,y=Meth.L2FC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("WAT-SC ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$Meth.L2FC,ourdf$RNA.L2FC),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("DMR L2FC") + xlab("DEG L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()



ourgene <- oursite_white_corgenes[5]
ourdf <- data.frame("RNA.L2FC" = whiternatrainingl2fcbysubject[ourgene,],
                    "Meth.L2FC" = whitemethtrainingl2fcbysubject[oursite,],
                    "Group" = pidmetadf[colnames(whiternatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(whiternatrainingl2fcbysubject),"Sex"])
png(file = "Supplemental Figure S19I_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA.L2FC,y=Meth.L2FC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("WAT-SC ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$Meth.L2FC,ourdf$RNA.L2FC),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("DMR L2FC") + xlab("DEG L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()
pdf(file = "Supplemental Figure S19I_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA.L2FC,y=Meth.L2FC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("WAT-SC ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$Meth.L2FC,ourdf$RNA.L2FC),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("DMR L2FC") + xlab("DEG L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()


oursite <- "chr11-70924_cluster3"
oursite_white_corgenes <- colnames(ourtftargetcormat)[abs(ourtftargetcormat[oursite,]) > 0.25]

ourgene <- oursite_white_corgenes[1]
ourdf <- data.frame("RNA.L2FC" = whiternatrainingl2fcbysubject[ourgene,],
                    "Meth.L2FC" = whitemethtrainingl2fcbysubject[oursite,],
                    "Group" = pidmetadf[colnames(whiternatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(whiternatrainingl2fcbysubject),"Sex"])
png(file = "Supplemental Figure S19E_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA.L2FC,y=Meth.L2FC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("WAT-SC ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$Meth.L2FC,ourdf$RNA.L2FC),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("DMR L2FC") + xlab("DEG L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()
pdf(file = "Supplemental Figure S19E_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA.L2FC,y=Meth.L2FC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("WAT-SC ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$Meth.L2FC,ourdf$RNA.L2FC),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("DMR L2FC") + xlab("DEG L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()

ourgene <- oursite_white_corgenes[2]
ourdf <- data.frame("RNA.L2FC" = whiternatrainingl2fcbysubject[ourgene,],
                    "Meth.L2FC" = whitemethtrainingl2fcbysubject[oursite,],
                    "Group" = pidmetadf[colnames(whiternatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(whiternatrainingl2fcbysubject),"Sex"])
png(file = "Supplemental Figure S19F_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA.L2FC,y=Meth.L2FC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("WAT-SC ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$Meth.L2FC,ourdf$RNA.L2FC),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("DMR L2FC") + xlab("DEG L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()
pdf(file = "Supplemental Figure S19F_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA.L2FC,y=Meth.L2FC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("WAT-SC ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$Meth.L2FC,ourdf$RNA.L2FC),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("DMR L2FC") + xlab("DEG L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()


ourtf <- "Sp2(Zf)/HEK293-Sp2.eGFP-ChIP-Seq(Encode)/Homer"
ourtftargetcormat <- gastromethcorfinal[rownames(gastromethcortfsfin[gastromethcortfsfin[,ourtf] > 0,]),]
ourtftargetcormat <- ourtftargetcormat[,apply(abs(ourtftargetcormat),2,max) > 0]

oursite <- "chr19-75322_cluster1"
oursite_gastro_corgenes <- colnames(ourtftargetcormat)[abs(ourtftargetcormat[oursite,]) > 0.25]

ourgene <- oursite_gastro_corgenes[1]
ourdf <- data.frame("RNA.L2FC" = gastrornatrainingl2fcbysubject[ourgene,],
                    "Meth.L2FC" = gastromethtrainingl2fcbysubject[oursite,],
                    "Group" = pidmetadf[colnames(gastrornatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(gastrornatrainingl2fcbysubject),"Sex"])
png(file = "Supplemental Figure S19D_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA.L2FC,y=Meth.L2FC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$Meth.L2FC,ourdf$RNA.L2FC),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("DMR L2FC") + xlab("DEG L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()
pdf(file = "Supplemental Figure S19D_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA.L2FC,y=Meth.L2FC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$Meth.L2FC,ourdf$RNA.L2FC),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("DMR L2FC") + xlab("DEG L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()




# Lung NF1 Halfsite examples
ourtf <- "NF1-halfsite(CTF)/LNCaP-NF1-ChIP-Seq(Unpublished)/Homer"
ourtftargetcormat <- lungmethcorfinal[rownames(lungmethcortfsfin[lungmethcortfsfin[,ourtf] > 0,]),]
ourtftargetcormat <- ourtftargetcormat[,apply(abs(ourtftargetcormat),2,max) > 0]

oursite <- "chr7-140730_cluster4"
oursite_lung_corgenes <- colnames(ourtftargetcormat)[abs(ourtftargetcormat[oursite,]) > 0.25]

ourgene <- oursite_lung_corgenes[1]
ourdf <- data.frame("RNA.L2FC" = lungrnatrainingl2fcbysubject[ourgene,],
                    "Meth.L2FC" = lungmethtrainingl2fcbysubject[oursite,],
                    "Group" = pidmetadf[colnames(lungrnatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(lungrnatrainingl2fcbysubject),"Sex"])
png(file = "Supplemental Figure S19L_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA.L2FC,y=Meth.L2FC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LUNG ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$Meth.L2FC,ourdf$RNA.L2FC),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("DMR L2FC") + xlab("DEG L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()
pdf(file = "Supplemental Figure S19L_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA.L2FC,y=Meth.L2FC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LUNG ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$Meth.L2FC,ourdf$RNA.L2FC),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("DMR L2FC") + xlab("DEG L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()

ourgene <- oursite_lung_corgenes[2]
ourdf <- data.frame("RNA.L2FC" = lungrnatrainingl2fcbysubject[ourgene,],
                    "Meth.L2FC" = lungmethtrainingl2fcbysubject[oursite,],
                    "Group" = pidmetadf[colnames(lungrnatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(lungrnatrainingl2fcbysubject),"Sex"])
png(file = "Supplemental Figure S19M_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA.L2FC,y=Meth.L2FC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LUNG ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$Meth.L2FC,ourdf$RNA.L2FC),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("DMR L2FC") + xlab("DEG L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()
pdf(file = "Supplemental Figure S19M_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA.L2FC,y=Meth.L2FC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LUNG ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$Meth.L2FC,ourdf$RNA.L2FC),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("DMR L2FC") + xlab("DEG L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()

ourgene <- oursite_lung_corgenes[3]
ourdf <- data.frame("RNA.L2FC" = lungrnatrainingl2fcbysubject[ourgene,],
                    "Meth.L2FC" = lungmethtrainingl2fcbysubject[oursite,],
                    "Group" = pidmetadf[colnames(lungrnatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(lungrnatrainingl2fcbysubject),"Sex"])
png(file = "Supplemental Figure S19K_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA.L2FC,y=Meth.L2FC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LUNG ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$Meth.L2FC,ourdf$RNA.L2FC),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("DMR L2FC") + xlab("DEG L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()
pdf(file = "Supplemental Figure S19K_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA.L2FC,y=Meth.L2FC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LUNG ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$Meth.L2FC,ourdf$RNA.L2FC),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("DMR L2FC") + xlab("DEG L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()

ourgene <- oursite_lung_corgenes[4]
ourdf <- data.frame("RNA.L2FC" = lungrnatrainingl2fcbysubject[ourgene,],
                    "Meth.L2FC" = lungmethtrainingl2fcbysubject[oursite,],
                    "Group" = pidmetadf[colnames(lungrnatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(lungrnatrainingl2fcbysubject),"Sex"])
png(file = "Supplemental Figure S19N_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA.L2FC,y=Meth.L2FC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LUNG ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$Meth.L2FC,ourdf$RNA.L2FC),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("DMR L2FC") + xlab("DEG L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()
pdf(file = "Supplemental Figure S19N_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA.L2FC,y=Meth.L2FC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LUNG ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$Meth.L2FC,ourdf$RNA.L2FC),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("DMR L2FC") + xlab("DEG L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()


ourtf <- "AP-2gamma(AP2)/MCF7-TFAP2C-ChIP-Seq(GSE21234)/Homer"
ourtftargetcormat <- lungmethcorfinal[rownames(lungmethcortfsfin[lungmethcortfsfin[,ourtf] > 0,]),]
ourtftargetcormat <- ourtftargetcormat[,apply(abs(ourtftargetcormat),2,max) > 0]

oursite <- "chr4-156530_cluster11"
oursite_lung_corgenes <- colnames(ourtftargetcormat)[abs(ourtftargetcormat[oursite,]) > 0.25]

ourgene <- oursite_lung_corgenes[2]
ourdf <- data.frame("RNA.L2FC" = lungrnatrainingl2fcbysubject[ourgene,],
                    "Meth.L2FC" = lungmethtrainingl2fcbysubject[oursite,],
                    "Group" = pidmetadf[colnames(lungrnatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(lungrnatrainingl2fcbysubject),"Sex"])
png(file = "Supplemental Figure S19O_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA.L2FC,y=Meth.L2FC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LUNG ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$Meth.L2FC,ourdf$RNA.L2FC),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("DMR L2FC") + xlab("DEG L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()
pdf(file = "Supplemental Figure S19O_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA.L2FC,y=Meth.L2FC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LUNG ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$Meth.L2FC,ourdf$RNA.L2FC),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("DMR L2FC") + xlab("DEG L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()

ourgene <- oursite_lung_corgenes[1]
ourdf <- data.frame("RNA.L2FC" = lungrnatrainingl2fcbysubject[ourgene,],
                    "Meth.L2FC" = lungmethtrainingl2fcbysubject[oursite,],
                    "Group" = pidmetadf[colnames(lungrnatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(lungrnatrainingl2fcbysubject),"Sex"])
png(file = "Supplemental Figure S19P_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA.L2FC,y=Meth.L2FC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LUNG ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$Meth.L2FC,ourdf$RNA.L2FC),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("DMR L2FC") + xlab("DEG L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()
pdf(file = "Supplemental Figure S19P_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA.L2FC,y=Meth.L2FC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LUNG ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$Meth.L2FC,ourdf$RNA.L2FC),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("DMR L2FC") + xlab("DEG L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()


#####
# Figure 4M
####

cortffractiondf <- data.frame("Enrichment.Fraction" = c(gastrocortfsig[gastrocortfsig[,"Enrichment_Significance"] < 0.05 & gastrocortfsig[,"Correlated_Target_Count"] > 3,"Fraction_Corr"],
                                                        lungcortfsig[lungcortfsig[,"Enrichment_Significance"] < 0.05 & lungcortfsig[,"Correlated_Target_Count"] > 3,"Fraction_Corr"],
                                                        whitecortfsig[whitecortfsig[,"Enrichment_Significance"] < 0.05 & whitecortfsig[,"Correlated_Target_Count"] > 3 & whitecortfsig[,"Fraction_Corr"] > 0.1,"Fraction_Corr"],
                                                        gastrocortfsig[gastrocortfsig[,"Enrichment_Significance"] < 0.05 & gastrocortfsig[,"Correlated_Target_Count"] > 3,"Fraction_Total"],
                                                        lungcortfsig[lungcortfsig[,"Enrichment_Significance"] < 0.05 & lungcortfsig[,"Correlated_Target_Count"] > 3,"Fraction_Total"],
                                                        whitecortfsig[whitecortfsig[,"Enrichment_Significance"] < 0.05 & whitecortfsig[,"Correlated_Target_Count"] > 3 & whitecortfsig[,"Fraction_Corr"] > 0.1,"Fraction_Total"]),
                              "P.Value" = c(gastrocortfsig[gastrocortfsig[,"Enrichment_Significance"] < 0.05 & gastrocortfsig[,"Correlated_Target_Count"] > 3,"Enrichment_Significance"],
                                            lungcortfsig[lungcortfsig[,"Enrichment_Significance"] < 0.05 & lungcortfsig[,"Correlated_Target_Count"] > 3,"Enrichment_Significance"],
                                            whitecortfsig[whitecortfsig[,"Enrichment_Significance"] < 0.05 & whitecortfsig[,"Correlated_Target_Count"] > 3 & whitecortfsig[,"Fraction_Corr"] > 0.1,"Enrichment_Significance"],
                                            gastrocortfsig[gastrocortfsig[,"Enrichment_Significance"] < 0.05 & gastrocortfsig[,"Correlated_Target_Count"] > 3,"Enrichment_Significance"],
                                            lungcortfsig[lungcortfsig[,"Enrichment_Significance"] < 0.05 & lungcortfsig[,"Correlated_Target_Count"] > 3,"Enrichment_Significance"],
                                            whitecortfsig[whitecortfsig[,"Enrichment_Significance"] < 0.05 & whitecortfsig[,"Correlated_Target_Count"] > 3 & whitecortfsig[,"Fraction_Corr"] > 0.1,"Enrichment_Significance"]),
                              "Target.List" = c(rep("Correlated DEG-DMR Pairs",(length(gastrocortfsig[gastrocortfsig[,"Enrichment_Significance"] < 0.05 & gastrocortfsig[,"Correlated_Target_Count"] > 3,"Fraction_Corr"])+
                                                                                  length(lungcortfsig[lungcortfsig[,"Enrichment_Significance"] < 0.05 & lungcortfsig[,"Correlated_Target_Count"] > 3,"Fraction_Corr"])+
                                                                                  length(whitecortfsig[whitecortfsig[,"Enrichment_Significance"] < 0.05 & whitecortfsig[,"Correlated_Target_Count"] > 3 & whitecortfsig[,"Fraction_Corr"] > 0.1,"Fraction_Corr"]))),
                                                rep("All Methylation Sites",(length(gastrocortfsig[gastrocortfsig[,"Enrichment_Significance"] < 0.05 & gastrocortfsig[,"Correlated_Target_Count"] > 3,"Fraction_Corr"])+
                                                                               length(lungcortfsig[lungcortfsig[,"Enrichment_Significance"] < 0.05 & lungcortfsig[,"Correlated_Target_Count"] > 3,"Fraction_Corr"])+
                                                                               length(whitecortfsig[whitecortfsig[,"Enrichment_Significance"] < 0.05 & whitecortfsig[,"Correlated_Target_Count"] > 3 & whitecortfsig[,"Fraction_Corr"] > 0.1,"Fraction_Corr"])))),
                              "Tissue" = c(rep("SKM-GN",length(gastrocortfsig[gastrocortfsig[,"Enrichment_Significance"] < 0.05 & gastrocortfsig[,"Correlated_Target_Count"] > 3,"Fraction_Corr"])),
                                           rep("LUNG",length(lungcortfsig[lungcortfsig[,"Enrichment_Significance"] < 0.05 & lungcortfsig[,"Correlated_Target_Count"] > 3,"Fraction_Corr"])),
                                           rep("WAT-SC",length(whitecortfsig[whitecortfsig[,"Enrichment_Significance"] < 0.05 & whitecortfsig[,"Correlated_Target_Count"] > 3 & whitecortfsig[,"Fraction_Corr"] > 0.1,"Fraction_Corr"])),
                                           rep("SKM-GN",length(gastrocortfsig[gastrocortfsig[,"Enrichment_Significance"] < 0.05 & gastrocortfsig[,"Correlated_Target_Count"] > 3,"Fraction_Corr"])),
                                           rep("LUNG",length(lungcortfsig[lungcortfsig[,"Enrichment_Significance"] < 0.05 & lungcortfsig[,"Correlated_Target_Count"] > 3,"Fraction_Corr"])),
                                           rep("WAT-SC",length(whitecortfsig[whitecortfsig[,"Enrichment_Significance"] < 0.05 & whitecortfsig[,"Correlated_Target_Count"] > 3 & whitecortfsig[,"Fraction_Corr"] > 0.1,"Fraction_Corr"]))),
                              "TF" = c(rownames(gastrocortfsig)[gastrocortfsig[,"Enrichment_Significance"] < 0.05 & gastrocortfsig[,"Correlated_Target_Count"] > 3],
                                       rownames(lungcortfsig)[lungcortfsig[,"Enrichment_Significance"] < 0.05 & lungcortfsig[,"Correlated_Target_Count"] > 3],
                                       rownames(whitecortfsig)[whitecortfsig[,"Enrichment_Significance"] < 0.05 & whitecortfsig[,"Correlated_Target_Count"] > 3 & whitecortfsig[,"Fraction_Corr"] > 0.1],
                                       rownames(gastrocortfsig)[gastrocortfsig[,"Enrichment_Significance"] < 0.05 & gastrocortfsig[,"Correlated_Target_Count"] > 3],
                                       rownames(lungcortfsig)[lungcortfsig[,"Enrichment_Significance"] < 0.05 & lungcortfsig[,"Correlated_Target_Count"] > 3],
                                       rownames(whitecortfsig)[whitecortfsig[,"Enrichment_Significance"] < 0.05 & whitecortfsig[,"Correlated_Target_Count"] > 3 & whitecortfsig[,"Fraction_Corr"] > 0.1]))
cortffractiondf$TF <- gsub("\\(.*","",cortffractiondf$TF)
cortffractiondf$Target.List <- factor(cortffractiondf$Target.List,levels = c("Correlated DEG-DMR Pairs","All Methylation Sites"))
cortffractiondf <- cortffractiondf[c(order(cortffractiondf$Enrichment.Fraction[1:(dim(cortffractiondf)[1]/2)]),(dim(cortffractiondf)[1]/2)+order(cortffractiondf$Enrichment.Fraction[1:(dim(cortffractiondf)[1]/2)])),]
cortffractiondf$TF.Number <- as.character(rep(c(1:(dim(cortffractiondf)[1]/2)),2))
cortffractiondf$TF.Number <- factor(cortffractiondf$TF.Number,levels = c(1:(dim(cortffractiondf)[1]/2)))
cortffractiondf$TF[1] <- "EWS:FLI1"
cortffractiondf$TF[26] <- "EWS:FLI1"


cortffractiondf$Significance <- "N"
cortffractiondf[cortffractiondf$Target.List %in% "Correlated DEG-DMR Pairs" & cortffractiondf$P.Value <= 0.05 & cortffractiondf$P.Value > 0.01,"Significance"] <- "*"
cortffractiondf[cortffractiondf$Target.List %in% "Correlated DEG-DMR Pairs" & cortffractiondf$P.Value <= 0.01 & cortffractiondf$P.Value > 0.001,"Significance"] <- "**"
cortffractiondf[cortffractiondf$Target.List %in% "Correlated DEG-DMR Pairs" & cortffractiondf$P.Value <= 0.001,"Significance"] <- "***"

png(file = "Figure 4M_021325.png",width = 10,height = 5,units = "in",res = 600)
ggplot(data = cortffractiondf,aes(x = TF.Number,y = Enrichment.Fraction,group = Target.List,shape = Target.List)) + 
  geom_line(aes(linetype=Target.List),linewidth = 1) + 
  geom_point(aes(color=Tissue),size = 4) + 
  geom_point(data = cortffractiondf[cortffractiondf$Significance %in% c("*","**","***"),], aes(x=TF.Number,y= Enrichment.Fraction + 0.05,size=Significance), color="black",shape = "*") +
  theme_classic() + 
  theme(text = element_text(size = 15),axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1)) +
  scale_x_discrete(labels=cortffractiondf$TF[1:(dim(cortffractiondf)[1]/2)]) +
  scale_size_discrete(labels=c('0.01 < p <= 0.05', '0.001 < p <= 0.01','p <= 0.001')) +
  scale_color_manual(values = ann_cols$Tissue)
dev.off()

pdf(file = "Figure 4M_021325.pdf",width = 10,height = 5)
ggplot(data = cortffractiondf,aes(x = TF.Number,y = Enrichment.Fraction,group = Target.List,shape = Target.List)) + 
  geom_line(aes(linetype=Target.List),linewidth = 1) + 
  geom_point(aes(color=Tissue),size = 4) + 
  geom_point(data = cortffractiondf[cortffractiondf$Significance %in% c("*","**","***"),], aes(x=TF.Number,y= Enrichment.Fraction + 0.05,size=Significance), color="black",shape = "*") +
  theme_classic() + 
  theme(text = element_text(size = 15),axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1)) +
  scale_x_discrete(labels=cortffractiondf$TF[1:(dim(cortffractiondf)[1]/2)]) +
  scale_size_discrete(labels=c('0.01 < p <= 0.05', '0.001 < p <= 0.01','p <= 0.001')) +
  scale_color_manual(values = ann_cols$Tissue)
dev.off()

#####
# Supplemental Table S3
####

load("rnal2fcmat.RData")
load("methsigl2fcmat.RData")

dmr_deg_connection50df <- data.frame("Tissue" = "",
                                     "DEG" = "",
                                     "SYMBOL" = "",
                                     "DMR" = "",
                                     "Distance" = 0,
                                     "Correlation" = 0,
                                     "TF Motifs" = "",
                                     "RNA FW1 L2FC" = 0,
                                     "RNA FW2 L2FC" = 0,
                                     "RNA FW4 L2FC" = 0,
                                     "RNA FW8 L2FC" = 0,
                                     "RNA MW1 L2FC" = 0,
                                     "RNA MW2 L2FC" = 0,
                                     "RNA MW4 L2FC" = 0,
                                     "RNA MW8 L2FC" = 0,
                                     "METH FW1 L2FC" = 0,
                                     "METH FW2 L2FC" = 0,
                                     "METH FW4 L2FC" = 0,
                                     "METH FW8 L2FC" = 0,
                                     "METH MW1 L2FC" = 0,
                                     "METH MW2 L2FC" = 0,
                                     "METH MW4 L2FC" = 0,
                                     "METH MW8 L2FC" = 0)

for(i in 1:dim(gastromethcorfinal)[2]){
  print(i)
  ourgene <- colnames(gastromethcorfinal)[i]
  oursites <- rownames(gastromethcorfinal)[gastromethcorfinal[,ourgene] != 0]
  if(length(oursites) == 1){
    oursite <- oursites
    ourdfentry <- data.frame("Tissue" = "SKM-GN",
                             "DEG" = ourgene,
                             "SYMBOL" = enstosym[ourgene,"Symbol"],
                             "DMR" = oursite,
                             "Distance" = abs(gastrornasiganno[ourgene,"start"] - ((gastromethsiganno[oursite,"start"]+gastromethsiganno[oursite,"end"])/2)),
                             "Correlation" = gastromethcorfinal[oursite,ourgene],
                             "TF Motifs" = toString(colnames(gastromethcortfsfintrim)[gastromethcortfsfintrim[oursite,] != 0]),
                             "RNA FW1 L2FC" = gastrol2fcmat[ourgene,"F W1"],
                             "RNA FW2 L2FC" = gastrol2fcmat[ourgene,"F W2"],
                             "RNA FW4 L2FC" = gastrol2fcmat[ourgene,"F W4"],
                             "RNA FW8 L2FC" = gastrol2fcmat[ourgene,"F W8"],
                             "RNA MW1 L2FC" = gastrol2fcmat[ourgene,"M W1"],
                             "RNA MW2 L2FC" = gastrol2fcmat[ourgene,"M W2"],
                             "RNA MW4 L2FC" = gastrol2fcmat[ourgene,"M W4"],
                             "RNA MW8 L2FC" = gastrol2fcmat[ourgene,"M W8"],
                             "METH FW1 L2FC" = gastromethl2fc[oursite,"F W1"],
                             "METH FW2 L2FC" = gastromethl2fc[oursite,"F W2"],
                             "METH FW4 L2FC" = gastromethl2fc[oursite,"F W4"],
                             "METH FW8 L2FC" = gastromethl2fc[oursite,"F W8"],
                             "METH MW1 L2FC" = gastromethl2fc[oursite,"M W1"],
                             "METH MW2 L2FC" = gastromethl2fc[oursite,"M W2"],
                             "METH MW4 L2FC" = gastromethl2fc[oursite,"M W4"],
                             "METH MW8 L2FC" = gastromethl2fc[oursite,"M W8"])
    dmr_deg_connection50df <- rbind(dmr_deg_connection50df,ourdfentry)
  } else if(length(oursites) > 1){
    for(j in 1:length(oursites)){
      oursite <- oursites[j]
      ourdfentry <- data.frame("Tissue" = "SKM-GN",
                               "DEG" = ourgene,
                               "SYMBOL" = enstosym[ourgene,"Symbol"],
                               "DMR" = oursite,
                               "Distance" = abs(gastrornasiganno[ourgene,"start"] - ((gastromethsiganno[oursite,"start"]+gastromethsiganno[oursite,"end"])/2)),
                               "Correlation" = gastromethcorfinal[oursite,ourgene],
                               "TF Motifs" = toString(colnames(gastromethcortfsfintrim)[gastromethcortfsfintrim[oursite,] != 0]),
                               "RNA FW1 L2FC" = gastrol2fcmat[ourgene,"F W1"],
                               "RNA FW2 L2FC" = gastrol2fcmat[ourgene,"F W2"],
                               "RNA FW4 L2FC" = gastrol2fcmat[ourgene,"F W4"],
                               "RNA FW8 L2FC" = gastrol2fcmat[ourgene,"F W8"],
                               "RNA MW1 L2FC" = gastrol2fcmat[ourgene,"M W1"],
                               "RNA MW2 L2FC" = gastrol2fcmat[ourgene,"M W2"],
                               "RNA MW4 L2FC" = gastrol2fcmat[ourgene,"M W4"],
                               "RNA MW8 L2FC" = gastrol2fcmat[ourgene,"M W8"],
                               "METH FW1 L2FC" = gastromethl2fc[oursite,"F W1"],
                               "METH FW2 L2FC" = gastromethl2fc[oursite,"F W2"],
                               "METH FW4 L2FC" = gastromethl2fc[oursite,"F W4"],
                               "METH FW8 L2FC" = gastromethl2fc[oursite,"F W8"],
                               "METH MW1 L2FC" = gastromethl2fc[oursite,"M W1"],
                               "METH MW2 L2FC" = gastromethl2fc[oursite,"M W2"],
                               "METH MW4 L2FC" = gastromethl2fc[oursite,"M W4"],
                               "METH MW8 L2FC" = gastromethl2fc[oursite,"M W8"])
      dmr_deg_connection50df <- rbind(dmr_deg_connection50df,ourdfentry)
      
    }
  } else {
    print("found nothing")
  }
}


for(i in 1:dim(heartmethcorfinal)[2]){
  print(i)
  ourgene <- colnames(heartmethcorfinal)[i]
  oursites <- rownames(heartmethcorfinal)[heartmethcorfinal[,ourgene] != 0]
  if(length(oursites) == 1){
    oursite <- oursites
    ourdfentry <- data.frame("Tissue" = "HEART",
                             "DEG" = ourgene,
                             "SYMBOL" = enstosym[ourgene,"Symbol"],
                             "DMR" = oursite,
                             "Distance" = abs(heartrnasiganno[ourgene,"start"] - ((heartmethsiganno[oursite,"start"]+heartmethsiganno[oursite,"end"])/2)),
                             "Correlation" = heartmethcorfinal[oursite,ourgene],
                             "TF Motifs" = toString(colnames(heartmethcortfsfintrim)[heartmethcortfsfintrim[oursite,] != 0]),
                             "RNA FW1 L2FC" = heartl2fcmat[ourgene,"F W1"],
                             "RNA FW2 L2FC" = heartl2fcmat[ourgene,"F W2"],
                             "RNA FW4 L2FC" = heartl2fcmat[ourgene,"F W4"],
                             "RNA FW8 L2FC" = heartl2fcmat[ourgene,"F W8"],
                             "RNA MW1 L2FC" = heartl2fcmat[ourgene,"M W1"],
                             "RNA MW2 L2FC" = heartl2fcmat[ourgene,"M W2"],
                             "RNA MW4 L2FC" = heartl2fcmat[ourgene,"M W4"],
                             "RNA MW8 L2FC" = heartl2fcmat[ourgene,"M W8"],
                             "METH FW1 L2FC" = heartmethl2fc[oursite,"F W1"],
                             "METH FW2 L2FC" = heartmethl2fc[oursite,"F W2"],
                             "METH FW4 L2FC" = heartmethl2fc[oursite,"F W4"],
                             "METH FW8 L2FC" = heartmethl2fc[oursite,"F W8"],
                             "METH MW1 L2FC" = heartmethl2fc[oursite,"M W1"],
                             "METH MW2 L2FC" = heartmethl2fc[oursite,"M W2"],
                             "METH MW4 L2FC" = heartmethl2fc[oursite,"M W4"],
                             "METH MW8 L2FC" = heartmethl2fc[oursite,"M W8"])
    dmr_deg_connection50df <- rbind(dmr_deg_connection50df,ourdfentry)
  } else if(length(oursites) > 1){
    for(j in 1:length(oursites)){
      oursite <- oursites[j]
      ourdfentry <- data.frame("Tissue" = "HEART",
                               "DEG" = ourgene,
                               "SYMBOL" = enstosym[ourgene,"Symbol"],
                               "DMR" = oursite,
                               "Distance" = abs(heartrnasiganno[ourgene,"start"] - ((heartmethsiganno[oursite,"start"]+heartmethsiganno[oursite,"end"])/2)),
                               "Correlation" = heartmethcorfinal[oursite,ourgene],
                               "TF Motifs" = toString(colnames(heartmethcortfsfintrim)[heartmethcortfsfintrim[oursite,] != 0]),
                               "RNA FW1 L2FC" = heartl2fcmat[ourgene,"F W1"],
                               "RNA FW2 L2FC" = heartl2fcmat[ourgene,"F W2"],
                               "RNA FW4 L2FC" = heartl2fcmat[ourgene,"F W4"],
                               "RNA FW8 L2FC" = heartl2fcmat[ourgene,"F W8"],
                               "RNA MW1 L2FC" = heartl2fcmat[ourgene,"M W1"],
                               "RNA MW2 L2FC" = heartl2fcmat[ourgene,"M W2"],
                               "RNA MW4 L2FC" = heartl2fcmat[ourgene,"M W4"],
                               "RNA MW8 L2FC" = heartl2fcmat[ourgene,"M W8"],
                               "METH FW1 L2FC" = heartmethl2fc[oursite,"F W1"],
                               "METH FW2 L2FC" = heartmethl2fc[oursite,"F W2"],
                               "METH FW4 L2FC" = heartmethl2fc[oursite,"F W4"],
                               "METH FW8 L2FC" = heartmethl2fc[oursite,"F W8"],
                               "METH MW1 L2FC" = heartmethl2fc[oursite,"M W1"],
                               "METH MW2 L2FC" = heartmethl2fc[oursite,"M W2"],
                               "METH MW4 L2FC" = heartmethl2fc[oursite,"M W4"],
                               "METH MW8 L2FC" = heartmethl2fc[oursite,"M W8"])
      dmr_deg_connection50df <- rbind(dmr_deg_connection50df,ourdfentry)
      
    }
  } else {
    print("found nothing")
  }
}


for(i in 1:dim(livermethcorfinal)[2]){
  print(i)
  ourgene <- colnames(livermethcorfinal)[i]
  oursites <- rownames(livermethcorfinal)[livermethcorfinal[,ourgene] != 0]
  if(length(oursites) == 1){
    oursite <- oursites
    ourdfentry <- data.frame("Tissue" = "LIVER",
                             "DEG" = ourgene,
                             "SYMBOL" = enstosym[ourgene,"Symbol"],
                             "DMR" = oursite,
                             "Distance" = abs(liverrnasiganno[ourgene,"start"] - ((livermethsiganno[oursite,"start"]+livermethsiganno[oursite,"end"])/2)),
                             "Correlation" = livermethcorfinal[oursite,ourgene],
                             "TF Motifs" = toString(colnames(livermethcortfsfintrim)[livermethcortfsfintrim[oursite,] != 0]),
                             "RNA FW1 L2FC" = liverl2fcmat[ourgene,"F W1"],
                             "RNA FW2 L2FC" = liverl2fcmat[ourgene,"F W2"],
                             "RNA FW4 L2FC" = liverl2fcmat[ourgene,"F W4"],
                             "RNA FW8 L2FC" = liverl2fcmat[ourgene,"F W8"],
                             "RNA MW1 L2FC" = liverl2fcmat[ourgene,"M W1"],
                             "RNA MW2 L2FC" = liverl2fcmat[ourgene,"M W2"],
                             "RNA MW4 L2FC" = liverl2fcmat[ourgene,"M W4"],
                             "RNA MW8 L2FC" = liverl2fcmat[ourgene,"M W8"],
                             "METH FW1 L2FC" = livermethl2fc[oursite,"F W1"],
                             "METH FW2 L2FC" = livermethl2fc[oursite,"F W2"],
                             "METH FW4 L2FC" = livermethl2fc[oursite,"F W4"],
                             "METH FW8 L2FC" = livermethl2fc[oursite,"F W8"],
                             "METH MW1 L2FC" = livermethl2fc[oursite,"M W1"],
                             "METH MW2 L2FC" = livermethl2fc[oursite,"M W2"],
                             "METH MW4 L2FC" = livermethl2fc[oursite,"M W4"],
                             "METH MW8 L2FC" = livermethl2fc[oursite,"M W8"])
    dmr_deg_connection50df <- rbind(dmr_deg_connection50df,ourdfentry)
  } else if(length(oursites) > 1){
    for(j in 1:length(oursites)){
      oursite <- oursites[j]
      ourdfentry <- data.frame("Tissue" = "LIVER",
                               "DEG" = ourgene,
                               "SYMBOL" = enstosym[ourgene,"Symbol"],
                               "DMR" = oursite,
                               "Distance" = abs(liverrnasiganno[ourgene,"start"] - ((livermethsiganno[oursite,"start"]+livermethsiganno[oursite,"end"])/2)),
                               "Correlation" = livermethcorfinal[oursite,ourgene],
                               "TF Motifs" = toString(colnames(livermethcortfsfintrim)[livermethcortfsfintrim[oursite,] != 0]),
                               "RNA FW1 L2FC" = liverl2fcmat[ourgene,"F W1"],
                               "RNA FW2 L2FC" = liverl2fcmat[ourgene,"F W2"],
                               "RNA FW4 L2FC" = liverl2fcmat[ourgene,"F W4"],
                               "RNA FW8 L2FC" = liverl2fcmat[ourgene,"F W8"],
                               "RNA MW1 L2FC" = liverl2fcmat[ourgene,"M W1"],
                               "RNA MW2 L2FC" = liverl2fcmat[ourgene,"M W2"],
                               "RNA MW4 L2FC" = liverl2fcmat[ourgene,"M W4"],
                               "RNA MW8 L2FC" = liverl2fcmat[ourgene,"M W8"],
                               "METH FW1 L2FC" = livermethl2fc[oursite,"F W1"],
                               "METH FW2 L2FC" = livermethl2fc[oursite,"F W2"],
                               "METH FW4 L2FC" = livermethl2fc[oursite,"F W4"],
                               "METH FW8 L2FC" = livermethl2fc[oursite,"F W8"],
                               "METH MW1 L2FC" = livermethl2fc[oursite,"M W1"],
                               "METH MW2 L2FC" = livermethl2fc[oursite,"M W2"],
                               "METH MW4 L2FC" = livermethl2fc[oursite,"M W4"],
                               "METH MW8 L2FC" = livermethl2fc[oursite,"M W8"])
      dmr_deg_connection50df <- rbind(dmr_deg_connection50df,ourdfentry)
      
    }
  } else {
    print("found nothing")
  }
}


for(i in 1:dim(lungmethcorfinal)[2]){
  print(i)
  ourgene <- colnames(lungmethcorfinal)[i]
  oursites <- rownames(lungmethcorfinal)[lungmethcorfinal[,ourgene] != 0]
  if(length(oursites) == 1){
    oursite <- oursites
    ourdfentry <- data.frame("Tissue" = "LUNG",
                             "DEG" = ourgene,
                             "SYMBOL" = enstosym[ourgene,"Symbol"],
                             "DMR" = oursite,
                             "Distance" = abs(lungrnasiganno[ourgene,"start"] - ((lungmethsiganno[oursite,"start"]+lungmethsiganno[oursite,"end"])/2)),
                             "Correlation" = lungmethcorfinal[oursite,ourgene],
                             "TF Motifs" = toString(colnames(lungmethcortfsfintrim)[lungmethcortfsfintrim[oursite,] != 0]),
                             "RNA FW1 L2FC" = lungl2fcmat[ourgene,"F W1"],
                             "RNA FW2 L2FC" = lungl2fcmat[ourgene,"F W2"],
                             "RNA FW4 L2FC" = lungl2fcmat[ourgene,"F W4"],
                             "RNA FW8 L2FC" = lungl2fcmat[ourgene,"F W8"],
                             "RNA MW1 L2FC" = lungl2fcmat[ourgene,"M W1"],
                             "RNA MW2 L2FC" = lungl2fcmat[ourgene,"M W2"],
                             "RNA MW4 L2FC" = lungl2fcmat[ourgene,"M W4"],
                             "RNA MW8 L2FC" = lungl2fcmat[ourgene,"M W8"],
                             "METH FW1 L2FC" = lungmethl2fc[oursite,"F W1"],
                             "METH FW2 L2FC" = lungmethl2fc[oursite,"F W2"],
                             "METH FW4 L2FC" = lungmethl2fc[oursite,"F W4"],
                             "METH FW8 L2FC" = lungmethl2fc[oursite,"F W8"],
                             "METH MW1 L2FC" = lungmethl2fc[oursite,"M W1"],
                             "METH MW2 L2FC" = lungmethl2fc[oursite,"M W2"],
                             "METH MW4 L2FC" = lungmethl2fc[oursite,"M W4"],
                             "METH MW8 L2FC" = lungmethl2fc[oursite,"M W8"])
    dmr_deg_connection50df <- rbind(dmr_deg_connection50df,ourdfentry)
  } else if(length(oursites) > 1){
    for(j in 1:length(oursites)){
      oursite <- oursites[j]
      ourdfentry <- data.frame("Tissue" = "LUNG",
                               "DEG" = ourgene,
                               "SYMBOL" = enstosym[ourgene,"Symbol"],
                               "DMR" = oursite,
                               "Distance" = abs(lungrnasiganno[ourgene,"start"] - ((lungmethsiganno[oursite,"start"]+lungmethsiganno[oursite,"end"])/2)),
                               "Correlation" = lungmethcorfinal[oursite,ourgene],
                               "TF Motifs" = toString(colnames(lungmethcortfsfintrim)[lungmethcortfsfintrim[oursite,] != 0]),
                               "RNA FW1 L2FC" = lungl2fcmat[ourgene,"F W1"],
                               "RNA FW2 L2FC" = lungl2fcmat[ourgene,"F W2"],
                               "RNA FW4 L2FC" = lungl2fcmat[ourgene,"F W4"],
                               "RNA FW8 L2FC" = lungl2fcmat[ourgene,"F W8"],
                               "RNA MW1 L2FC" = lungl2fcmat[ourgene,"M W1"],
                               "RNA MW2 L2FC" = lungl2fcmat[ourgene,"M W2"],
                               "RNA MW4 L2FC" = lungl2fcmat[ourgene,"M W4"],
                               "RNA MW8 L2FC" = lungl2fcmat[ourgene,"M W8"],
                               "METH FW1 L2FC" = lungmethl2fc[oursite,"F W1"],
                               "METH FW2 L2FC" = lungmethl2fc[oursite,"F W2"],
                               "METH FW4 L2FC" = lungmethl2fc[oursite,"F W4"],
                               "METH FW8 L2FC" = lungmethl2fc[oursite,"F W8"],
                               "METH MW1 L2FC" = lungmethl2fc[oursite,"M W1"],
                               "METH MW2 L2FC" = lungmethl2fc[oursite,"M W2"],
                               "METH MW4 L2FC" = lungmethl2fc[oursite,"M W4"],
                               "METH MW8 L2FC" = lungmethl2fc[oursite,"M W8"])
      dmr_deg_connection50df <- rbind(dmr_deg_connection50df,ourdfentry)
      
    }
  } else {
    print("found nothing")
  }
}

for(i in 1:dim(whitemethcorfinal)[2]){
  print(i)
  ourgene <- colnames(whitemethcorfinal)[i]
  oursites <- rownames(whitemethcorfinal)[whitemethcorfinal[,ourgene] != 0]
  if(length(oursites) == 1){
    oursite <- oursites
    ourdfentry <- data.frame("Tissue" = "WAT-SC",
                             "DEG" = ourgene,
                             "SYMBOL" = enstosym[ourgene,"Symbol"],
                             "DMR" = oursite,
                             "Distance" = abs(whiternasiganno[ourgene,"start"] - ((whitemethsiganno[oursite,"start"]+whitemethsiganno[oursite,"end"])/2)),
                             "Correlation" = whitemethcorfinal[oursite,ourgene],
                             "TF Motifs" = toString(colnames(whitemethcortfsfintrim)[whitemethcortfsfintrim[oursite,] != 0]),
                             "RNA FW1 L2FC" = whitel2fcmat[ourgene,"F W1"],
                             "RNA FW2 L2FC" = whitel2fcmat[ourgene,"F W2"],
                             "RNA FW4 L2FC" = whitel2fcmat[ourgene,"F W4"],
                             "RNA FW8 L2FC" = whitel2fcmat[ourgene,"F W8"],
                             "RNA MW1 L2FC" = whitel2fcmat[ourgene,"M W1"],
                             "RNA MW2 L2FC" = whitel2fcmat[ourgene,"M W2"],
                             "RNA MW4 L2FC" = whitel2fcmat[ourgene,"M W4"],
                             "RNA MW8 L2FC" = whitel2fcmat[ourgene,"M W8"],
                             "METH FW1 L2FC" = whitemethl2fc[oursite,"F W1"],
                             "METH FW2 L2FC" = whitemethl2fc[oursite,"F W2"],
                             "METH FW4 L2FC" = whitemethl2fc[oursite,"F W4"],
                             "METH FW8 L2FC" = whitemethl2fc[oursite,"F W8"],
                             "METH MW1 L2FC" = whitemethl2fc[oursite,"M W1"],
                             "METH MW2 L2FC" = whitemethl2fc[oursite,"M W2"],
                             "METH MW4 L2FC" = whitemethl2fc[oursite,"M W4"],
                             "METH MW8 L2FC" = whitemethl2fc[oursite,"M W8"])
    dmr_deg_connection50df <- rbind(dmr_deg_connection50df,ourdfentry)
  } else if(length(oursites) > 1){
    for(j in 1:length(oursites)){
      oursite <- oursites[j]
      ourdfentry <- data.frame("Tissue" = "WAT-SC",
                               "DEG" = ourgene,
                               "SYMBOL" = enstosym[ourgene,"Symbol"],
                               "DMR" = oursite,
                               "Distance" = abs(whiternasiganno[ourgene,"start"] - ((whitemethsiganno[oursite,"start"]+whitemethsiganno[oursite,"end"])/2)),
                               "Correlation" = whitemethcorfinal[oursite,ourgene],
                               "TF Motifs" = toString(colnames(whitemethcortfsfintrim)[whitemethcortfsfintrim[oursite,] != 0]),
                               "RNA FW1 L2FC" = whitel2fcmat[ourgene,"F W1"],
                               "RNA FW2 L2FC" = whitel2fcmat[ourgene,"F W2"],
                               "RNA FW4 L2FC" = whitel2fcmat[ourgene,"F W4"],
                               "RNA FW8 L2FC" = whitel2fcmat[ourgene,"F W8"],
                               "RNA MW1 L2FC" = whitel2fcmat[ourgene,"M W1"],
                               "RNA MW2 L2FC" = whitel2fcmat[ourgene,"M W2"],
                               "RNA MW4 L2FC" = whitel2fcmat[ourgene,"M W4"],
                               "RNA MW8 L2FC" = whitel2fcmat[ourgene,"M W8"],
                               "METH FW1 L2FC" = whitemethl2fc[oursite,"F W1"],
                               "METH FW2 L2FC" = whitemethl2fc[oursite,"F W2"],
                               "METH FW4 L2FC" = whitemethl2fc[oursite,"F W4"],
                               "METH FW8 L2FC" = whitemethl2fc[oursite,"F W8"],
                               "METH MW1 L2FC" = whitemethl2fc[oursite,"M W1"],
                               "METH MW2 L2FC" = whitemethl2fc[oursite,"M W2"],
                               "METH MW4 L2FC" = whitemethl2fc[oursite,"M W4"],
                               "METH MW8 L2FC" = whitemethl2fc[oursite,"M W8"])
      dmr_deg_connection50df <- rbind(dmr_deg_connection50df,ourdfentry)
      
    }
  } else {
    print("found nothing")
  }
}

dmr_deg_connection50df <- dmr_deg_connection50df[2:dim(dmr_deg_connection50df)[1],]

# Make supplemental table - remove empty row 1 and numbered list in column 1
write.csv(dmr_deg_connection50df,file = "Supplemental Table S3_021325.csv",row.names = F)



save.image("Figure4M_to_4N_S19_SuppTableS3_021325.RData")

