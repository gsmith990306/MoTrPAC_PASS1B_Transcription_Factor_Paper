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

setwd("~/Mount Sinai/Sealfon Laboratory/MoTrPAC/PASS1B Raw Data")

load("Figure3D_3G_S9.RData")
allmethmotifs <- read.table(file = "allmethoutput50.txt",header = T,sep = "\t")
load("enstosym.RData")

####
# Next we want to look for distal relationships between DEGs and DMRs

# Let's start with gastroc. We have a 166x730 matrix of DMRs vs DEGs

####
# gastro

gastromethcortrim <- gastromethcor*(gastromethdistance < 500000 & gastromethdistance > 0)
gastromethcortrimmer <- gastromethcortrim[apply(abs(gastromethcortrim),1,max) >= 0.5,]
gastromethcorfin <- gastromethcortrimmer[,apply(abs(gastromethcortrimmer),2,max) >= 0.5]

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
heartmethcortrimmer <- heartmethcortrim[apply(abs(heartmethcortrim),1,max) >= 0.5,]
heartmethcorfin <- heartmethcortrimmer[,apply(abs(heartmethcortrimmer),2,max) >= 0.5]

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
hippomethcortrimmer <- hippomethcortrim[apply(abs(hippomethcortrim),1,max) >= 0.5,]
hippomethcorfin <- hippomethcortrimmer[,apply(abs(hippomethcortrimmer),2,max) >= 0.5]

hippomethcortrimtfs <- allmethmotifs[allmethmotifs$PositionID %in% rownames(hippomethcorfin),]
hippomethcortfs <- matrix(0L,nrow = dim(hippomethcorfin)[1],ncol = length(unique(allmethmotifs[allmethmotifs$PositionID %in% rownames(hippomethcorfin),"Motif.Name"])))
rownames(hippomethcortfs) <- rownames(hippomethcorfin)
colnames(hippomethcortfs) <- unique(allmethmotifs[allmethmotifs$PositionID %in% rownames(hippomethcorfin),"Motif.Name"])
for(i in 1:dim(hippomethcortrimtfs)[1]){
  hippomethcortfs[hippomethcortrimtfs$PositionID[i],hippomethcortrimtfs$Motif.Name[i]] <- 1
}
hippomethcortfsfintrim <- hippomethcortfs[apply(hippomethcortfs,1,max) > 0,]
# No more correlations so we end analysis in hippocampus
#hippomethcortfsfin <- hippomethcortfsfintrim[,apply(hippomethcortfsfintrim,2,max) > 0]

#hippomethcorfintrim <- hippomethcorfin[rownames(hippomethcortfsfin),]
#hippomethcorfinal <- hippomethcorfintrim[,apply(abs(hippomethcorfintrim),2,max) > 0]

#png(filename = "hippo_DMR_DEG_correlated_neighbors.png",width = 5,height = 5,units = "in",res = 600)
#pheatmap(hippomethcorfinal,labels_col = enstosym[colnames(hippomethcorfinal),"Symbol"],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"))
#dev.off()

####
# kidney

kidneymethcortrim <- kidneymethcor*(kidneymethdistance < 500000 & kidneymethdistance > 0)
kidneymethcortrimmer <- kidneymethcortrim[apply(abs(kidneymethcortrim),1,max) >= 0.5,]
kidneymethcorfin <- kidneymethcortrimmer[,apply(abs(kidneymethcortrimmer),2,max) >= 0.5]

kidneymethcortrimtfs <- allmethmotifs[allmethmotifs$PositionID %in% rownames(kidneymethcorfin),]
kidneymethcortfs <- matrix(0L,nrow = dim(kidneymethcorfin)[1],ncol = length(unique(allmethmotifs[allmethmotifs$PositionID %in% rownames(kidneymethcorfin),"Motif.Name"])))
rownames(kidneymethcortfs) <- rownames(kidneymethcorfin)
colnames(kidneymethcortfs) <- unique(allmethmotifs[allmethmotifs$PositionID %in% rownames(kidneymethcorfin),"Motif.Name"])
for(i in 1:dim(kidneymethcortrimtfs)[1]){
  kidneymethcortfs[kidneymethcortrimtfs$PositionID[i],kidneymethcortrimtfs$Motif.Name[i]] <- 1
}
kidneymethcortfsfintrim <- kidneymethcortfs[apply(kidneymethcortfs,1,max) > 0,]
# No remaining correlations so we end analysis in kidney
#kidneymethcortfsfin <- kidneymethcortfsfintrim[,apply(kidneymethcortfsfintrim,2,max) > 0]

#kidneymethcorfintrim <- kidneymethcorfin[rownames(kidneymethcortfsfin),]
#kidneymethcorfinal <- kidneymethcorfintrim[,apply(abs(kidneymethcorfintrim),2,max) > 0]

#png(filename = "kidney_DMR_DEG_correlated_neighbors.png",width = 5,height = 5,units = "in",res = 600)
#pheatmap(kidneymethcorfinal,labels_col = enstosym[colnames(kidneymethcorfinal),"Symbol"],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"))
#dev.off()

####
# liver

livermethcortrim <- livermethcor*(livermethdistance < 500000 & livermethdistance > 0)
livermethcortrimmer <- livermethcortrim[apply(abs(livermethcortrim),1,max) >= 0.5,]
livermethcorfin <- livermethcortrimmer[,apply(abs(livermethcortrimmer),2,max) >= 0.5]

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
lungmethcortrimmer <- lungmethcortrim[apply(abs(lungmethcortrim),1,max) >= 0.5,]
lungmethcorfin <- lungmethcortrimmer[,apply(abs(lungmethcortrimmer),2,max) >= 0.5]

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
brownmethcortrimmer <- brownmethcortrim[apply(abs(brownmethcortrim),1,max) >= 0.5,]
brownmethcorfin <- brownmethcortrimmer[,apply(abs(brownmethcortrimmer),2,max) >= 0.5]

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
whitemethcortrimmer <- whitemethcortrim[apply(abs(whitemethcortrim),1,max) >= 0.5,]
whitemethcorfin <- whitemethcortrimmer[,apply(abs(whitemethcortrimmer),2,max) >= 0.5]

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
png(filename = "Supplemental Figure S13H.png",width = 5,height = 5,units = "in",res = 600)
pheatmap(ourtftargetcormat,labels_col = enstosym[colnames(ourtftargetcormat),"Symbol"],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),main = paste("LUNG ",gsub("\\(.*","",ourtf)," DEG-DMR Targets",sep = ""))
dev.off()
pdf(file = "Supplemental Figure S13H.pdf",width = 5,height = 5)
pheatmap(ourtftargetcormat,labels_col = enstosym[colnames(ourtftargetcormat),"Symbol"],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),main = paste("LUNG ",gsub("\\(.*","",ourtf)," DEG-DMR Targets",sep = ""))
dev.off()

# WAT-SC NF1-halfsite(CTF)/LNCaP-NF1-ChIP-Seq(Unpublished)/Homer

ourtf <- "NF1-halfsite(CTF)/LNCaP-NF1-ChIP-Seq(Unpublished)/Homer"
ourtftargetcormat <- whitemethcorfinal[rownames(whitemethcortfsfin[whitemethcortfsfin[,ourtf] > 0,]),]
ourtftargetcormat <- ourtftargetcormat[,apply(abs(ourtftargetcormat),2,max) > 0]
png(filename = "Supplemental Figure S13A.png",width = 12,height = 5,units = "in",res = 600)
pheatmap(ourtftargetcormat,labels_col = enstosym[colnames(ourtftargetcormat),"Symbol"],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),main = paste("WAT-SC ",gsub("\\(.*","",ourtf)," DEG-DMR Targets",sep = ""))
dev.off()
pdf(file = "Supplemental Figure S13A.pdf",width = 12,height = 5)
pheatmap(ourtftargetcormat,labels_col = enstosym[colnames(ourtftargetcormat),"Symbol"],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),main = paste("WAT-SC ",gsub("\\(.*","",ourtf)," DEG-DMR Targets",sep = ""))
dev.off()


# WAT-SC AP-2gamma(AP2)/MCF7-TFAP2C-ChIP-Seq(GSE21234)/Homer

ourtf <- "AP-2gamma(AP2)/MCF7-TFAP2C-ChIP-Seq(GSE21234)/Homer"
ourtftargetcormat <- whitemethcorfinal[rownames(whitemethcortfsfin[whitemethcortfsfin[,ourtf] > 0,]),]
ourtftargetcormat <- ourtftargetcormat[,apply(abs(ourtftargetcormat),2,max) > 0]
png(file = "Supplemental Figure S13G.png",width = 12,height = 5,units = "in",res = 600)
pheatmap(ourtftargetcormat,labels_col = enstosym[colnames(ourtftargetcormat),"Symbol"],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),main = paste("WAT-SC ",gsub("\\(.*","",ourtf)," DEG-DMR Targets",sep = ""))
dev.off()
pdf(file = "Supplemental Figure S13G.pdf",width = 12,height = 5)
pheatmap(ourtftargetcormat,labels_col = enstosym[colnames(ourtftargetcormat),"Symbol"],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),main = paste("WAT-SC ",gsub("\\(.*","",ourtf)," DEG-DMR Targets",sep = ""))
dev.off()

####
# Now to generate scatter plots highlighting the correlations between individual DEG-DMR pairs
ourtf <- "NF1-halfsite(CTF)/LNCaP-NF1-ChIP-Seq(Unpublished)/Homer"
ourtftargetcormat <- whitemethcorfinal[rownames(whitemethcortfsfin[whitemethcortfsfin[,ourtf] > 0,]),]
ourtftargetcormat <- ourtftargetcormat[,apply(abs(ourtftargetcormat),2,max) > 0]

oursite <- "chr2-375476_cluster1"
oursite_white_corgenes <- colnames(ourtftargetcormat)[abs(ourtftargetcormat[oursite,]) > 0.5]

ourgene <- oursite_white_corgenes[2]
ourdf <- data.frame("RNA.L2FC" = whitel2fcmat[ourgene,],
                    "Meth.L2FC" = whitemethl2fc[oursite,],
                    "Label" = colnames(whitel2fcmat))
png(file = "New Supplemental Figure S13E.png",width = 6,height = 6,units = "in",res = 600)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=Meth.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("WAT-SC: ",enstosym[ourgene,"Symbol"]," vs ",oursite,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$Meth.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()
pdf(file = "New Supplemental Figure S13E.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=Meth.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("WAT-SC: ",enstosym[ourgene,"Symbol"]," vs ",oursite,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$Meth.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- oursite_white_corgenes[3]
ourdf <- data.frame("RNA.L2FC" = whitel2fcmat[ourgene,],
                    "Meth.L2FC" = whitemethl2fc[oursite,],
                    "Label" = colnames(whitel2fcmat))
png(file = "New Supplemental Figure S13F.png",width = 6,height = 6,units = "in",res = 600)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=Meth.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("WAT-SC: ",enstosym[ourgene,"Symbol"]," vs ",oursite,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$Meth.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()
pdf(file = "New Supplemental Figure S13F.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=Meth.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("WAT-SC: ",enstosym[ourgene,"Symbol"]," vs ",oursite,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$Meth.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- oursite_white_corgenes[4]
ourdf <- data.frame("RNA.L2FC" = whitel2fcmat[ourgene,],
                    "Meth.L2FC" = whitemethl2fc[oursite,],
                    "Label" = colnames(whitel2fcmat))
png(file = "New Supplemental Figure S13C.png",width = 6,height = 6,units = "in",res = 600)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=Meth.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("WAT-SC: ",enstosym[ourgene,"Symbol"]," vs ",oursite,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$Meth.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()
pdf(file = "New Supplemental Figure S13C.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=Meth.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("WAT-SC: ",enstosym[ourgene,"Symbol"]," vs ",oursite,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$Meth.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- oursite_white_corgenes[5]
ourdf <- data.frame("RNA.L2FC" = whitel2fcmat[ourgene,],
                    "Meth.L2FC" = whitemethl2fc[oursite,],
                    "Label" = colnames(whitel2fcmat))
png(file = "New Supplemental Figure S13B.png",width = 6,height = 6,units = "in",res = 600)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=Meth.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("WAT-SC: ",enstosym[ourgene,"Symbol"]," vs ",oursite,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$Meth.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()
pdf(file = "New Supplemental Figure S13B.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=Meth.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("WAT-SC: ",enstosym[ourgene,"Symbol"]," vs ",oursite,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$Meth.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- oursite_white_corgenes[6]
ourdf <- data.frame("RNA.L2FC" = whitel2fcmat[ourgene,],
                    "Meth.L2FC" = whitemethl2fc[oursite,],
                    "Label" = colnames(whitel2fcmat))
png(file = "New Supplemental Figure S13D.png",width = 6,height = 6,units = "in",res = 600)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=Meth.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("WAT-SC: ",enstosym[ourgene,"Symbol"]," vs ",oursite,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$Meth.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()
pdf(file = "New Supplemental Figure S13D.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=Meth.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("WAT-SC: ",enstosym[ourgene,"Symbol"]," vs ",oursite,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$Meth.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()


ourtf <- "AP-2gamma(AP2)/MCF7-TFAP2C-ChIP-Seq(GSE21234)/Homer"
ourtftargetcormat <- whitemethcorfinal[rownames(whitemethcortfsfin[whitemethcortfsfin[,ourtf] > 0,]),]
ourtftargetcormat <- ourtftargetcormat[,apply(abs(ourtftargetcormat),2,max) > 0]

ourgene <- enstosym[enstosym$Symbol %in% "B4galnt1","Ensembl"]
b4galnt1_white_corsites <- rownames(ourtftargetcormat)[abs(ourtftargetcormat[,ourgene]) > 0.5]
oursite <- b4galnt1_white_corsites[1]
ourdf <- data.frame("RNA.L2FC" = whitel2fcmat[ourgene,],
                    "Meth.L2FC" = whitemethl2fc[oursite,],
                    "Label" = colnames(whitel2fcmat))
png(file = "Supplemental Figure S13I.png",width = 6,height = 6,units = "in",res = 600)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=Meth.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("WAT-SC: ",enstosym[ourgene,"Symbol"]," vs ",oursite,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$Meth.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()
pdf(file = "Supplemental Figure S13I.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=Meth.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("WAT-SC: ",enstosym[ourgene,"Symbol"]," vs ",oursite,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$Meth.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()


ourtf <- "AP-2gamma(AP2)/MCF7-TFAP2C-ChIP-Seq(GSE21234)/Homer"
ourtftargetcormat <- lungmethcorfinal[rownames(lungmethcortfsfin[lungmethcortfsfin[,ourtf] > 0,]),]
ourtftargetcormat <- ourtftargetcormat[,apply(abs(ourtftargetcormat),2,max) > 0]


ourgene <- enstosym[enstosym$Symbol %in% "B4galnt1","Ensembl"]
b4galnt1_lung_corsites <- rownames(ourtftargetcormat)[abs(ourtftargetcormat[,ourgene]) > 0.5]
oursite <- b4galnt1_lung_corsites[1]
ourdf <- data.frame("RNA.L2FC" = lungl2fcmat[ourgene,],
                    "Meth.L2FC" = lungmethl2fc[oursite,],
                    "Label" = colnames(lungl2fcmat))
png(file = "Supplemental Figure S13J.png",width = 6,height = 6,units = "in",res = 600)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=Meth.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LUNG: ",enstosym[ourgene,"Symbol"]," vs ",oursite,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$Meth.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()
pdf(file = "Supplemental Figure S13J.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=Meth.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LUNG: ",enstosym[ourgene,"Symbol"]," vs ",oursite,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$Meth.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Arhgap9","Ensembl"]
arhgap9_lung_corsites <- rownames(ourtftargetcormat)[abs(ourtftargetcormat[,ourgene]) > 0.5]
oursite <- arhgap9_lung_corsites[1]
ourdf <- data.frame("RNA.L2FC" = lungl2fcmat[ourgene,],
                    "Meth.L2FC" = lungmethl2fc[oursite,],
                    "Label" = colnames(lungl2fcmat))
png(file = "Supplemental Figure S13K.png",width = 6,height = 6,units = "in",res = 600)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=Meth.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LUNG: ",enstosym[ourgene,"Symbol"]," vs ",oursite,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$Meth.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()
pdf(file = "Supplemental Figure S13K.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=Meth.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LUNG: ",enstosym[ourgene,"Symbol"]," vs ",oursite,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$Meth.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()


#####
# Figure 4H
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

png(file = "Figure 4H.png",width = 10,height = 5,units = "in",res = 600)
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

pdf(file = "Figure 4H.pdf",width = 10,height = 5)
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

save.image("Figure4H_S13.RData")


#####
# Supplemental Table S2
####

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
    ourdfentry <- data.frame("Tissue" = "white",
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
      ourdfentry <- data.frame("Tissue" = "white",
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
write.csv(dmr_deg_connection50df,file = "Supplemental Table 2.csv",row.names = F)

