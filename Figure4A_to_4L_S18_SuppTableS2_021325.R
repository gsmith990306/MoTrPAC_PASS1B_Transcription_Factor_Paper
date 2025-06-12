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

load("Figure3C_3F_S14_112624.RData")
#load("New_43024_Figure3C_3F_S8.RData")
allpeakmotifs <- readRDS("allpeakmotifs.RDS")
load("activepeakfiles.RData")
tfanno <- readRDS("tfanno.RDS")
load("rnal2fcmat.RData")
load("atacsigl2fcmat.RData")


gastro50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% gastroactivepeaks,]
heart50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% heartactivepeaks,]
hippo50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% hippoactivepeaks,]
kidney50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% kidneyactivepeaks,]
liver50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% liveractivepeaks,]
lung50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% lungactivepeaks,]
brown50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% brownactivepeaks,]
white50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% whiteactivepeaks,]


gastropeakcormod <- gastropeakcor*(gastropeakdistance > 0 & gastropeakdistance < 500000)
gastropeakcortestmod <- gastropeakcortest*(gastropeakdistance > 0 & gastropeakdistance < 500000)
gastropeakcortestmod[gastropeakcortestmod == 0] <- 1
gastroconnectsigpeaks <- rownames(gastropeakcortestmod[apply(abs(gastropeakcortestmod),1,min) < 0.05,])
gastroconnectsiggenes <- colnames(gastropeakcortestmod[,apply(abs(gastropeakcortestmod),2,min) < 0.05])
gastropeakcortrim <- gastropeakcormod[gastroconnectsigpeaks,gastroconnectsiggenes]

gastro50connectsigpeakstfanno <- gastroconnectsigpeaks[gastroconnectsigpeaks %in% gastro50peakmotifs$PositionID]
gastro50connectsigpeakstfs <- unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% gastro50connectsigpeakstfanno,"Motif.Name"])
gastro50expressedsigpeakstfens <- intersect(tfanno[gastro50connectsigpeakstfs,"Ensembl"],rownames(gastrornanorm))
gastro50expressedsigpeakstfs <- rownames(tfanno[tfanno$Ensembl %in% gastro50expressedsigpeakstfens,])
gastro50connectsigpeaksexprtfanno <- intersect(gastro50connectsigpeakstfanno,gastro50peakmotifs[gastro50peakmotifs$Motif.Name %in% gastro50expressedsigpeakstfs,"PositionID"])


heartpeakcormod <- heartpeakcor*(heartpeakdistance > 0 & heartpeakdistance < 500000)
heartpeakcortestmod <- heartpeakcortest*(heartpeakdistance > 0 & heartpeakdistance < 500000)
heartpeakcortestmod[heartpeakcortestmod == 0] <- 1
heartconnectsigpeaks <- rownames(heartpeakcortestmod[apply(abs(heartpeakcortestmod),1,min) < 0.05,])
heartconnectsiggenes <- colnames(heartpeakcortestmod[,apply(abs(heartpeakcortestmod),2,min) < 0.05])
heartpeakcortrim <- heartpeakcormod[heartconnectsigpeaks,heartconnectsiggenes]

heart50connectsigpeakstfanno <- heartconnectsigpeaks[heartconnectsigpeaks %in% heart50peakmotifs$PositionID]
heart50connectsigpeakstfs <- unique(heart50peakmotifs[heart50peakmotifs$PositionID %in% heart50connectsigpeakstfanno,"Motif.Name"])
heart50expressedsigpeakstfens <- intersect(tfanno[heart50connectsigpeakstfs,"Ensembl"],rownames(heartrnanorm))
heart50expressedsigpeakstfs <- rownames(tfanno[tfanno$Ensembl %in% heart50expressedsigpeakstfens,])
heart50connectsigpeaksexprtfanno <- intersect(heart50connectsigpeakstfanno,heart50peakmotifs[heart50peakmotifs$Motif.Name %in% heart50expressedsigpeakstfs,"PositionID"])


hippopeakcormod <- hippopeakcor*(hippopeakdistance > 0 & hippopeakdistance < 500000)
hippopeakcortestmod <- hippopeakcortest*(hippopeakdistance > 0 & hippopeakdistance < 500000)
hippopeakcortestmod[hippopeakcortestmod == 0] <- 1
hippoconnectsigpeaks <- rownames(hippopeakcortestmod[apply(abs(hippopeakcortestmod),1,min) < 0.05,])
hippoconnectsiggenes <- colnames(hippopeakcortestmod[,apply(abs(hippopeakcortestmod),2,min) < 0.05])
hippopeakcortrim <- hippopeakcormod[hippoconnectsigpeaks,hippoconnectsiggenes]

hippo50connectsigpeakstfanno <- hippoconnectsigpeaks[hippoconnectsigpeaks %in% hippo50peakmotifs$PositionID]
hippo50connectsigpeakstfs <- unique(hippo50peakmotifs[hippo50peakmotifs$PositionID %in% hippo50connectsigpeakstfanno,"Motif.Name"])
hippo50expressedsigpeakstfens <- intersect(tfanno[hippo50connectsigpeakstfs,"Ensembl"],rownames(hippornanorm))
hippo50expressedsigpeakstfs <- rownames(tfanno[tfanno$Ensembl %in% hippo50expressedsigpeakstfens,])
hippo50connectsigpeaksexprtfanno <- intersect(hippo50connectsigpeakstfanno,hippo50peakmotifs[hippo50peakmotifs$Motif.Name %in% hippo50expressedsigpeakstfs,"PositionID"])


kidneypeakcormod <- kidneypeakcor*(kidneypeakdistance > 0 & kidneypeakdistance < 500000)
kidneypeakcortestmod <- kidneypeakcortest*(kidneypeakdistance > 0 & kidneypeakdistance < 500000)
kidneypeakcortestmod[kidneypeakcortestmod == 0] <- 1
kidneyconnectsigpeaks <- rownames(kidneypeakcortestmod[apply(abs(kidneypeakcortestmod),1,min) < 0.05,])
kidneyconnectsiggenes <- colnames(kidneypeakcortestmod[,apply(abs(kidneypeakcortestmod),2,min) < 0.05])
kidneypeakcortrim <- kidneypeakcormod[kidneyconnectsigpeaks,kidneyconnectsiggenes]

kidney50connectsigpeakstfanno <- kidneyconnectsigpeaks[kidneyconnectsigpeaks %in% kidney50peakmotifs$PositionID]
kidney50connectsigpeakstfs <- unique(kidney50peakmotifs[kidney50peakmotifs$PositionID %in% kidney50connectsigpeakstfanno,"Motif.Name"])
kidney50expressedsigpeakstfens <- intersect(tfanno[kidney50connectsigpeakstfs,"Ensembl"],rownames(kidneyrnanorm))
kidney50expressedsigpeakstfs <- rownames(tfanno[tfanno$Ensembl %in% kidney50expressedsigpeakstfens,])
kidney50connectsigpeaksexprtfanno <- intersect(kidney50connectsigpeakstfanno,kidney50peakmotifs[kidney50peakmotifs$Motif.Name %in% kidney50expressedsigpeakstfs,"PositionID"])


liverpeakcormod <- liverpeakcor*(liverpeakdistance > 0 & liverpeakdistance < 500000)
liverpeakcortestmod <- liverpeakcortest*(liverpeakdistance > 0 & liverpeakdistance < 500000)
liverpeakcortestmod[liverpeakcortestmod == 0] <- 1
liverconnectsigpeaks <- rownames(liverpeakcortestmod[apply(abs(liverpeakcortestmod),1,min) < 0.05,])
liverconnectsiggenes <- colnames(liverpeakcortestmod[,apply(abs(liverpeakcortestmod),2,min) < 0.05])
liverpeakcortrim <- liverpeakcormod[liverconnectsigpeaks,liverconnectsiggenes]

liver50connectsigpeakstfanno <- liverconnectsigpeaks[liverconnectsigpeaks %in% liver50peakmotifs$PositionID]
liver50connectsigpeakstfs <- unique(liver50peakmotifs[liver50peakmotifs$PositionID %in% liver50connectsigpeakstfanno,"Motif.Name"])
liver50expressedsigpeakstfens <- intersect(tfanno[liver50connectsigpeakstfs,"Ensembl"],rownames(liverrnanorm))
liver50expressedsigpeakstfs <- rownames(tfanno[tfanno$Ensembl %in% liver50expressedsigpeakstfens,])
liver50connectsigpeaksexprtfanno <- intersect(liver50connectsigpeakstfanno,liver50peakmotifs[liver50peakmotifs$Motif.Name %in% liver50expressedsigpeakstfs,"PositionID"])


lungpeakcormod <- lungpeakcor*(lungpeakdistance > 0 & lungpeakdistance < 500000)
lungpeakcortestmod <- lungpeakcortest*(lungpeakdistance > 0 & lungpeakdistance < 500000)
lungpeakcortestmod[lungpeakcortestmod == 0] <- 1
lungconnectsigpeaks <- rownames(lungpeakcortestmod[apply(abs(lungpeakcortestmod),1,min) < 0.05,])
lungconnectsiggenes <- colnames(lungpeakcortestmod[,apply(abs(lungpeakcortestmod),2,min) < 0.05])
lungpeakcortrim <- lungpeakcormod[lungconnectsigpeaks,lungconnectsiggenes]

lung50connectsigpeakstfanno <- lungconnectsigpeaks[lungconnectsigpeaks %in% lung50peakmotifs$PositionID]
lung50connectsigpeakstfs <- unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% lung50connectsigpeakstfanno,"Motif.Name"])
lung50expressedsigpeakstfens <- intersect(tfanno[lung50connectsigpeakstfs,"Ensembl"],rownames(lungrnanorm))
lung50expressedsigpeakstfs <- rownames(tfanno[tfanno$Ensembl %in% lung50expressedsigpeakstfens,])
lung50connectsigpeaksexprtfanno <- intersect(lung50connectsigpeakstfanno,lung50peakmotifs[lung50peakmotifs$Motif.Name %in% lung50expressedsigpeakstfs,"PositionID"])


brownpeakcormod <- brownpeakcor*(brownpeakdistance > 0 & brownpeakdistance < 500000)
brownpeakcortestmod <- brownpeakcortest*(brownpeakdistance > 0 & brownpeakdistance < 500000)
brownpeakcortestmod[brownpeakcortestmod == 0] <- 1
brownconnectsigpeaks <- rownames(brownpeakcortestmod[apply(abs(brownpeakcortestmod),1,min) < 0.05,])
brownconnectsiggenes <- colnames(brownpeakcortestmod[,apply(abs(brownpeakcortestmod),2,min) < 0.05])
brownpeakcortrim <- brownpeakcormod[brownconnectsigpeaks,brownconnectsiggenes]

brown50connectsigpeakstfanno <- brownconnectsigpeaks[brownconnectsigpeaks %in% brown50peakmotifs$PositionID]
brown50connectsigpeakstfs <- unique(brown50peakmotifs[brown50peakmotifs$PositionID %in% brown50connectsigpeakstfanno,"Motif.Name"])
brown50expressedsigpeakstfens <- intersect(tfanno[brown50connectsigpeakstfs,"Ensembl"],rownames(brownrnanorm))
brown50expressedsigpeakstfs <- rownames(tfanno[tfanno$Ensembl %in% brown50expressedsigpeakstfens,])
brown50connectsigpeaksexprtfanno <- intersect(brown50connectsigpeakstfanno,brown50peakmotifs[brown50peakmotifs$Motif.Name %in% brown50expressedsigpeakstfs,"PositionID"])


whitepeakcormod <- whitepeakcor*(whitepeakdistance > 0 & whitepeakdistance < 500000)
whitepeakcortestmod <- whitepeakcortest*(whitepeakdistance > 0 & whitepeakdistance < 500000)
whitepeakcortestmod[whitepeakcortestmod == 0] <- 1
whiteconnectsigpeaks <- rownames(whitepeakcortestmod[apply(abs(whitepeakcortestmod),1,min) < 0.05,])
whiteconnectsiggenes <- colnames(whitepeakcortestmod[,apply(abs(whitepeakcortestmod),2,min) < 0.05])
whitepeakcortrim <- whitepeakcormod[whiteconnectsigpeaks,whiteconnectsiggenes]

white50connectsigpeakstfanno <- whiteconnectsigpeaks[whiteconnectsigpeaks %in% white50peakmotifs$PositionID]
white50connectsigpeakstfs <- unique(white50peakmotifs[white50peakmotifs$PositionID %in% white50connectsigpeakstfanno,"Motif.Name"])
white50expressedsigpeakstfens <- intersect(tfanno[white50connectsigpeakstfs,"Ensembl"],rownames(whiternanorm))
white50expressedsigpeakstfs <- rownames(tfanno[tfanno$Ensembl %in% white50expressedsigpeakstfens,])
white50connectsigpeaksexprtfanno <- intersect(white50connectsigpeakstfanno,white50peakmotifs[white50peakmotifs$Motif.Name %in% white50expressedsigpeakstfs,"PositionID"])



gastro50peakcortrimmer <- gastropeakcortrim[gastro50connectsigpeaksexprtfanno,]
gastro50peakcorfin <- gastro50peakcortrimmer[,abs(apply(gastro50peakcortrimmer,2,mean)) > 0]
heart50peakcortrimmer <- heartpeakcortrim[heart50connectsigpeaksexprtfanno,]
#heart50peakcorfin <- heart50peakcortrimmer[,abs(apply(heart50peakcortrimmer,2,mean)) > 0]
hippo50peakcortrimmer <- hippopeakcortrim[hippo50connectsigpeaksexprtfanno,]
hippo50peakcorfin <- hippo50peakcortrimmer[,abs(apply(hippo50peakcortrimmer,2,mean)) > 0]
kidney50peakcortrimmer <- kidneypeakcortrim[kidney50connectsigpeaksexprtfanno,]
kidney50peakcorfin <- kidney50peakcortrimmer[,abs(apply(kidney50peakcortrimmer,2,mean)) > 0]
liver50peakcortrimmer <- liverpeakcortrim[liver50connectsigpeaksexprtfanno,]
liver50peakcorfin <- liver50peakcortrimmer[,abs(apply(liver50peakcortrimmer,2,mean)) > 0]
lung50peakcortrimmer <- lungpeakcortrim[lung50connectsigpeaksexprtfanno,]
lung50peakcorfin <- lung50peakcortrimmer[,abs(apply(lung50peakcortrimmer,2,mean)) > 0]
brown50peakcortrimmer <- brownpeakcortrim[brown50connectsigpeaksexprtfanno,]
brown50peakcorfin <- brown50peakcortrimmer[,abs(apply(brown50peakcortrimmer,2,mean)) > 0]
white50peakcortrimmer <- whitepeakcortrim[white50connectsigpeaksexprtfanno,]
white50peakcorfin <- white50peakcortrimmer[,abs(apply(white50peakcortrimmer,2,mean)) > 0]

load("enstosym.RData")


pdf(file = "Supplemental Figure S18M_021325.pdf",width = 8,height = 5)
pheatmap(gastro50peakcorfin[unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin) & gastro50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),colnames(gastro50peakcorfin[unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin) & gastro50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),colSums(gastro50peakcorfin[unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin) & gastro50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),]) != 0])],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_col = enstosym[colnames(gastro50peakcorfin[unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin) & gastro50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),colSums(gastro50peakcorfin[unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin) & gastro50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),]) != 0]),"Symbol"],fontsize = 20,cluster_rows = F, cluster_cols = F)
dev.off()

png(file = "Supplemental Figure S18M_021325.png",width = 8,height = 5,units = "in",res = 600)
pheatmap(gastro50peakcorfin[unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin) & gastro50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),colnames(gastro50peakcorfin[unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin) & gastro50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),colSums(gastro50peakcorfin[unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin) & gastro50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),]) != 0])],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_col = enstosym[colnames(gastro50peakcorfin[unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin) & gastro50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),colSums(gastro50peakcorfin[unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin) & gastro50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),]) != 0]),"Symbol"],fontsize = 20,cluster_rows = F, cluster_cols = F)
dev.off()

pdf(file = "Supplemental Figure S18N_021325.pdf",width = 8,height = 5)
pheatmap(lung50peakcorfin[unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin) & lung50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),colnames(lung50peakcorfin[unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin) & lung50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),colSums(lung50peakcorfin[unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin) & lung50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),]) != 0])],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_col = enstosym[colnames(lung50peakcorfin[unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin) & lung50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),colSums(lung50peakcorfin[unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin) & lung50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),]) != 0]),"Symbol"],fontsize = 20,cluster_rows = F, cluster_cols = F)
dev.off()

png(file = "Supplemental Figure S18N_021325.png",width = 8,height = 5,units = "in",res = 600)
pheatmap(lung50peakcorfin[unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin) & lung50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),colnames(lung50peakcorfin[unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin) & lung50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),colSums(lung50peakcorfin[unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin) & lung50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),]) != 0])],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_col = enstosym[colnames(lung50peakcorfin[unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin) & lung50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),colSums(lung50peakcorfin[unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin) & lung50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),]) != 0]),"Symbol"],fontsize = 20,cluster_rows = F, cluster_cols = F)
dev.off()

pdf(file = "Figure 4D_021325.pdf",width = 8,height = 6)
pheatmap(liver50peakcorfin[unique(liver50peakmotifs[liver50peakmotifs$PositionID %in% rownames(liver50peakcorfin) & liver50peakmotifs$Motif.Name %in% "BMAL1(bHLH)/Liver-Bmal1-ChIP-Seq(GSE39860)/Homer","PositionID"]),colnames(liver50peakcorfin[unique(liver50peakmotifs[liver50peakmotifs$PositionID %in% rownames(liver50peakcorfin) & liver50peakmotifs$Motif.Name %in% "BMAL1(bHLH)/Liver-Bmal1-ChIP-Seq(GSE39860)/Homer","PositionID"]),colSums(liver50peakcorfin[unique(liver50peakmotifs[liver50peakmotifs$PositionID %in% rownames(liver50peakcorfin) & liver50peakmotifs$Motif.Name %in% "BMAL1(bHLH)/Liver-Bmal1-ChIP-Seq(GSE39860)/Homer","PositionID"]),]) != 0])],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_col = enstosym[colnames(liver50peakcorfin[unique(liver50peakmotifs[liver50peakmotifs$PositionID %in% rownames(liver50peakcorfin) & liver50peakmotifs$Motif.Name %in% "BMAL1(bHLH)/Liver-Bmal1-ChIP-Seq(GSE39860)/Homer","PositionID"]),colSums(liver50peakcorfin[unique(liver50peakmotifs[liver50peakmotifs$PositionID %in% rownames(liver50peakcorfin) & liver50peakmotifs$Motif.Name %in% "BMAL1(bHLH)/Liver-Bmal1-ChIP-Seq(GSE39860)/Homer","PositionID"]),]) != 0]),"Symbol"],fontsize = 20,cluster_rows = F, cluster_cols = F)
dev.off()

png(file = "Figure 4D_021325.png",width = 8,height = 6,units = "in",res = 600)
pheatmap(liver50peakcorfin[unique(liver50peakmotifs[liver50peakmotifs$PositionID %in% rownames(liver50peakcorfin) & liver50peakmotifs$Motif.Name %in% "BMAL1(bHLH)/Liver-Bmal1-ChIP-Seq(GSE39860)/Homer","PositionID"]),colnames(liver50peakcorfin[unique(liver50peakmotifs[liver50peakmotifs$PositionID %in% rownames(liver50peakcorfin) & liver50peakmotifs$Motif.Name %in% "BMAL1(bHLH)/Liver-Bmal1-ChIP-Seq(GSE39860)/Homer","PositionID"]),colSums(liver50peakcorfin[unique(liver50peakmotifs[liver50peakmotifs$PositionID %in% rownames(liver50peakcorfin) & liver50peakmotifs$Motif.Name %in% "BMAL1(bHLH)/Liver-Bmal1-ChIP-Seq(GSE39860)/Homer","PositionID"]),]) != 0])],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_col = enstosym[colnames(liver50peakcorfin[unique(liver50peakmotifs[liver50peakmotifs$PositionID %in% rownames(liver50peakcorfin) & liver50peakmotifs$Motif.Name %in% "BMAL1(bHLH)/Liver-Bmal1-ChIP-Seq(GSE39860)/Homer","PositionID"]),colSums(liver50peakcorfin[unique(liver50peakmotifs[liver50peakmotifs$PositionID %in% rownames(liver50peakcorfin) & liver50peakmotifs$Motif.Name %in% "BMAL1(bHLH)/Liver-Bmal1-ChIP-Seq(GSE39860)/Homer","PositionID"]),]) != 0]),"Symbol"],fontsize = 20,cluster_rows = F, cluster_cols = F)
dev.off()

#####
# Supplemental Table S2
####

dar_deg_connection50df <- data.frame("Tissue" = "",
                                     "DEG" = "",
                                     "SYMBOL" = "",
                                     "DAR" = "",
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
                                     "ATAC FW1 L2FC" = 0,
                                     "ATAC FW2 L2FC" = 0,
                                     "ATAC FW4 L2FC" = 0,
                                     "ATAC FW8 L2FC" = 0,
                                     "ATAC MW1 L2FC" = 0,
                                     "ATAC MW2 L2FC" = 0,
                                     "ATAC MW4 L2FC" = 0,
                                     "ATAC MW8 L2FC" = 0)

for(i in 1:dim(gastro50peakcorfin)[2]){
  print(i)
  ourgene <- colnames(gastro50peakcorfin)[i]
  ourpeaks <- rownames(gastro50peakcorfin)[gastro50peakcorfin[,ourgene] != 0]
  if(length(ourpeaks) == 1){
    ourpeak <- ourpeaks
    ourdfentry <- data.frame("Tissue" = "SKM-GN",
                             "DEG" = ourgene,
                             "SYMBOL" = enstosym[ourgene,"Symbol"],
                             "DAR" = ourpeak,
                             "Distance" = abs(gastrornasiganno[ourgene,"start"] - ((peakanno[ourpeak,"start"]+peakanno[ourpeak,"end"])/2)),
                             "Correlation" = gastro50peakcorfin[ourpeak,ourgene],
                             "TF Motifs" = toString(unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% ourpeak,"Motif.Name"])),
                             "RNA FW1 L2FC" = gastrol2fcmat[ourgene,"F W1"],
                             "RNA FW2 L2FC" = gastrol2fcmat[ourgene,"F W2"],
                             "RNA FW4 L2FC" = gastrol2fcmat[ourgene,"F W4"],
                             "RNA FW8 L2FC" = gastrol2fcmat[ourgene,"F W8"],
                             "RNA MW1 L2FC" = gastrol2fcmat[ourgene,"M W1"],
                             "RNA MW2 L2FC" = gastrol2fcmat[ourgene,"M W2"],
                             "RNA MW4 L2FC" = gastrol2fcmat[ourgene,"M W4"],
                             "RNA MW8 L2FC" = gastrol2fcmat[ourgene,"M W8"],
                             "ATAC FW1 L2FC" = gastrosigatacl2fc[ourpeak,"F W1"],
                             "ATAC FW2 L2FC" = gastrosigatacl2fc[ourpeak,"F W2"],
                             "ATAC FW4 L2FC" = gastrosigatacl2fc[ourpeak,"F W4"],
                             "ATAC FW8 L2FC" = gastrosigatacl2fc[ourpeak,"F W8"],
                             "ATAC MW1 L2FC" = gastrosigatacl2fc[ourpeak,"M W1"],
                             "ATAC MW2 L2FC" = gastrosigatacl2fc[ourpeak,"M W2"],
                             "ATAC MW4 L2FC" = gastrosigatacl2fc[ourpeak,"M W4"],
                             "ATAC MW8 L2FC" = gastrosigatacl2fc[ourpeak,"M W8"])
    dar_deg_connection50df <- rbind(dar_deg_connection50df,ourdfentry)
  } else if(length(ourpeaks) > 1){
    for(j in 1:length(ourpeaks)){
      ourpeak <- ourpeaks[j]
      ourdfentry <- data.frame("Tissue" = "SKM-GN",
                               "DEG" = ourgene,
                               "SYMBOL" = enstosym[ourgene,"Symbol"],
                               "DAR" = ourpeak,
                               "Distance" = abs(gastrornasiganno[ourgene,"start"] - ((peakanno[ourpeak,"start"]+peakanno[ourpeak,"end"])/2)),
                               "Correlation" = gastro50peakcorfin[ourpeak,ourgene],
                               "TF Motifs" = toString(unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% ourpeak,"Motif.Name"])),
                               "RNA FW1 L2FC" = gastrol2fcmat[ourgene,"F W1"],
                               "RNA FW2 L2FC" = gastrol2fcmat[ourgene,"F W2"],
                               "RNA FW4 L2FC" = gastrol2fcmat[ourgene,"F W4"],
                               "RNA FW8 L2FC" = gastrol2fcmat[ourgene,"F W8"],
                               "RNA MW1 L2FC" = gastrol2fcmat[ourgene,"M W1"],
                               "RNA MW2 L2FC" = gastrol2fcmat[ourgene,"M W2"],
                               "RNA MW4 L2FC" = gastrol2fcmat[ourgene,"M W4"],
                               "RNA MW8 L2FC" = gastrol2fcmat[ourgene,"M W8"],
                               "ATAC FW1 L2FC" = gastrosigatacl2fc[ourpeak,"F W1"],
                               "ATAC FW2 L2FC" = gastrosigatacl2fc[ourpeak,"F W2"],
                               "ATAC FW4 L2FC" = gastrosigatacl2fc[ourpeak,"F W4"],
                               "ATAC FW8 L2FC" = gastrosigatacl2fc[ourpeak,"F W8"],
                               "ATAC MW1 L2FC" = gastrosigatacl2fc[ourpeak,"M W1"],
                               "ATAC MW2 L2FC" = gastrosigatacl2fc[ourpeak,"M W2"],
                               "ATAC MW4 L2FC" = gastrosigatacl2fc[ourpeak,"M W4"],
                               "ATAC MW8 L2FC" = gastrosigatacl2fc[ourpeak,"M W8"])
      dar_deg_connection50df <- rbind(dar_deg_connection50df,ourdfentry)
      
    }
  } else {
    print("found nothing")
  }
}



for(i in 1:dim(kidney50peakcorfin)[2]){
  print(i)
  ourgene <- colnames(kidney50peakcorfin)[i]
  ourpeaks <- rownames(kidney50peakcorfin)[kidney50peakcorfin[,ourgene] != 0]
  if(length(ourpeaks) == 1){
    ourpeak <- ourpeaks
    ourdfentry <- data.frame("Tissue" = "KIDNEY",
                             "DEG" = ourgene,
                             "SYMBOL" = enstosym[ourgene,"Symbol"],
                             "DAR" = ourpeak,
                             "Distance" = abs(kidneyrnasiganno[ourgene,"start"] - ((peakanno[ourpeak,"start"]+peakanno[ourpeak,"end"])/2)),
                             "Correlation" = kidney50peakcorfin[ourpeak,ourgene],
                             "TF Motifs" = toString(unique(kidney50peakmotifs[kidney50peakmotifs$PositionID %in% ourpeak,"Motif.Name"])),
                             "RNA FW1 L2FC" = kidneyl2fcmat[ourgene,"F W1"],
                             "RNA FW2 L2FC" = kidneyl2fcmat[ourgene,"F W2"],
                             "RNA FW4 L2FC" = kidneyl2fcmat[ourgene,"F W4"],
                             "RNA FW8 L2FC" = kidneyl2fcmat[ourgene,"F W8"],
                             "RNA MW1 L2FC" = kidneyl2fcmat[ourgene,"M W1"],
                             "RNA MW2 L2FC" = kidneyl2fcmat[ourgene,"M W2"],
                             "RNA MW4 L2FC" = kidneyl2fcmat[ourgene,"M W4"],
                             "RNA MW8 L2FC" = kidneyl2fcmat[ourgene,"M W8"],
                             "ATAC FW1 L2FC" = kidneysigatacl2fc[ourpeak,"F W1"],
                             "ATAC FW2 L2FC" = kidneysigatacl2fc[ourpeak,"F W2"],
                             "ATAC FW4 L2FC" = kidneysigatacl2fc[ourpeak,"F W4"],
                             "ATAC FW8 L2FC" = kidneysigatacl2fc[ourpeak,"F W8"],
                             "ATAC MW1 L2FC" = kidneysigatacl2fc[ourpeak,"M W1"],
                             "ATAC MW2 L2FC" = kidneysigatacl2fc[ourpeak,"M W2"],
                             "ATAC MW4 L2FC" = kidneysigatacl2fc[ourpeak,"M W4"],
                             "ATAC MW8 L2FC" = kidneysigatacl2fc[ourpeak,"M W8"])
    dar_deg_connection50df <- rbind(dar_deg_connection50df,ourdfentry)
  } else if(length(ourpeaks) > 1){
    for(j in 1:length(ourpeaks)){
      ourpeak <- ourpeaks[j]
      ourdfentry <- data.frame("Tissue" = "KIDNEY",
                               "DEG" = ourgene,
                               "SYMBOL" = enstosym[ourgene,"Symbol"],
                               "DAR" = ourpeak,
                               "Distance" = abs(kidneyrnasiganno[ourgene,"start"] - ((peakanno[ourpeak,"start"]+peakanno[ourpeak,"end"])/2)),
                               "Correlation" = kidney50peakcorfin[ourpeak,ourgene],
                               "TF Motifs" = toString(unique(kidney50peakmotifs[kidney50peakmotifs$PositionID %in% ourpeak,"Motif.Name"])),
                               "RNA FW1 L2FC" = kidneyl2fcmat[ourgene,"F W1"],
                               "RNA FW2 L2FC" = kidneyl2fcmat[ourgene,"F W2"],
                               "RNA FW4 L2FC" = kidneyl2fcmat[ourgene,"F W4"],
                               "RNA FW8 L2FC" = kidneyl2fcmat[ourgene,"F W8"],
                               "RNA MW1 L2FC" = kidneyl2fcmat[ourgene,"M W1"],
                               "RNA MW2 L2FC" = kidneyl2fcmat[ourgene,"M W2"],
                               "RNA MW4 L2FC" = kidneyl2fcmat[ourgene,"M W4"],
                               "RNA MW8 L2FC" = kidneyl2fcmat[ourgene,"M W8"],
                               "ATAC FW1 L2FC" = kidneysigatacl2fc[ourpeak,"F W1"],
                               "ATAC FW2 L2FC" = kidneysigatacl2fc[ourpeak,"F W2"],
                               "ATAC FW4 L2FC" = kidneysigatacl2fc[ourpeak,"F W4"],
                               "ATAC FW8 L2FC" = kidneysigatacl2fc[ourpeak,"F W8"],
                               "ATAC MW1 L2FC" = kidneysigatacl2fc[ourpeak,"M W1"],
                               "ATAC MW2 L2FC" = kidneysigatacl2fc[ourpeak,"M W2"],
                               "ATAC MW4 L2FC" = kidneysigatacl2fc[ourpeak,"M W4"],
                               "ATAC MW8 L2FC" = kidneysigatacl2fc[ourpeak,"M W8"])
      dar_deg_connection50df <- rbind(dar_deg_connection50df,ourdfentry)
      
    }
  } else {
    print("found nothing")
  }
}


for(i in 1:dim(liver50peakcorfin)[2]){
  print(i)
  ourgene <- colnames(liver50peakcorfin)[i]
  ourpeaks <- rownames(liver50peakcorfin)[liver50peakcorfin[,ourgene] != 0]
  if(length(ourpeaks) == 1){
    ourpeak <- ourpeaks
    ourdfentry <- data.frame("Tissue" = "LIVER",
                             "DEG" = ourgene,
                             "SYMBOL" = enstosym[ourgene,"Symbol"],
                             "DAR" = ourpeak,
                             "Distance" = abs(liverrnasiganno[ourgene,"start"] - ((peakanno[ourpeak,"start"]+peakanno[ourpeak,"end"])/2)),
                             "Correlation" = liver50peakcorfin[ourpeak,ourgene],
                             "TF Motifs" = toString(unique(liver50peakmotifs[liver50peakmotifs$PositionID %in% ourpeak,"Motif.Name"])),
                             "RNA FW1 L2FC" = liverl2fcmat[ourgene,"F W1"],
                             "RNA FW2 L2FC" = liverl2fcmat[ourgene,"F W2"],
                             "RNA FW4 L2FC" = liverl2fcmat[ourgene,"F W4"],
                             "RNA FW8 L2FC" = liverl2fcmat[ourgene,"F W8"],
                             "RNA MW1 L2FC" = liverl2fcmat[ourgene,"M W1"],
                             "RNA MW2 L2FC" = liverl2fcmat[ourgene,"M W2"],
                             "RNA MW4 L2FC" = liverl2fcmat[ourgene,"M W4"],
                             "RNA MW8 L2FC" = liverl2fcmat[ourgene,"M W8"],
                             "ATAC FW1 L2FC" = liversigatacl2fc[ourpeak,"F W1"],
                             "ATAC FW2 L2FC" = liversigatacl2fc[ourpeak,"F W2"],
                             "ATAC FW4 L2FC" = liversigatacl2fc[ourpeak,"F W4"],
                             "ATAC FW8 L2FC" = liversigatacl2fc[ourpeak,"F W8"],
                             "ATAC MW1 L2FC" = liversigatacl2fc[ourpeak,"M W1"],
                             "ATAC MW2 L2FC" = liversigatacl2fc[ourpeak,"M W2"],
                             "ATAC MW4 L2FC" = liversigatacl2fc[ourpeak,"M W4"],
                             "ATAC MW8 L2FC" = liversigatacl2fc[ourpeak,"M W8"])
    dar_deg_connection50df <- rbind(dar_deg_connection50df,ourdfentry)
  } else if(length(ourpeaks) > 1){
    for(j in 1:length(ourpeaks)){
      ourpeak <- ourpeaks[j]
      ourdfentry <- data.frame("Tissue" = "LIVER",
                               "DEG" = ourgene,
                               "SYMBOL" = enstosym[ourgene,"Symbol"],
                               "DAR" = ourpeak,
                               "Distance" = abs(liverrnasiganno[ourgene,"start"] - ((peakanno[ourpeak,"start"]+peakanno[ourpeak,"end"])/2)),
                               "Correlation" = liver50peakcorfin[ourpeak,ourgene],
                               "TF Motifs" = toString(unique(liver50peakmotifs[liver50peakmotifs$PositionID %in% ourpeak,"Motif.Name"])),
                               "RNA FW1 L2FC" = liverl2fcmat[ourgene,"F W1"],
                               "RNA FW2 L2FC" = liverl2fcmat[ourgene,"F W2"],
                               "RNA FW4 L2FC" = liverl2fcmat[ourgene,"F W4"],
                               "RNA FW8 L2FC" = liverl2fcmat[ourgene,"F W8"],
                               "RNA MW1 L2FC" = liverl2fcmat[ourgene,"M W1"],
                               "RNA MW2 L2FC" = liverl2fcmat[ourgene,"M W2"],
                               "RNA MW4 L2FC" = liverl2fcmat[ourgene,"M W4"],
                               "RNA MW8 L2FC" = liverl2fcmat[ourgene,"M W8"],
                               "ATAC FW1 L2FC" = liversigatacl2fc[ourpeak,"F W1"],
                               "ATAC FW2 L2FC" = liversigatacl2fc[ourpeak,"F W2"],
                               "ATAC FW4 L2FC" = liversigatacl2fc[ourpeak,"F W4"],
                               "ATAC FW8 L2FC" = liversigatacl2fc[ourpeak,"F W8"],
                               "ATAC MW1 L2FC" = liversigatacl2fc[ourpeak,"M W1"],
                               "ATAC MW2 L2FC" = liversigatacl2fc[ourpeak,"M W2"],
                               "ATAC MW4 L2FC" = liversigatacl2fc[ourpeak,"M W4"],
                               "ATAC MW8 L2FC" = liversigatacl2fc[ourpeak,"M W8"])
      dar_deg_connection50df <- rbind(dar_deg_connection50df,ourdfentry)
      
    }
  } else {
    print("found nothing")
  }
}


for(i in 1:dim(lung50peakcorfin)[2]){
  print(i)
  ourgene <- colnames(lung50peakcorfin)[i]
  ourpeaks <- rownames(lung50peakcorfin)[lung50peakcorfin[,ourgene] != 0]
  if(length(ourpeaks) == 1){
    ourpeak <- ourpeaks
    ourdfentry <- data.frame("Tissue" = "LUNG",
                             "DEG" = ourgene,
                             "SYMBOL" = enstosym[ourgene,"Symbol"],
                             "DAR" = ourpeak,
                             "Distance" = abs(lungrnasiganno[ourgene,"start"] - ((peakanno[ourpeak,"start"]+peakanno[ourpeak,"end"])/2)),
                             "Correlation" = lung50peakcorfin[ourpeak,ourgene],
                             "TF Motifs" = toString(unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% ourpeak,"Motif.Name"])),
                             "RNA FW1 L2FC" = lungl2fcmat[ourgene,"F W1"],
                             "RNA FW2 L2FC" = lungl2fcmat[ourgene,"F W2"],
                             "RNA FW4 L2FC" = lungl2fcmat[ourgene,"F W4"],
                             "RNA FW8 L2FC" = lungl2fcmat[ourgene,"F W8"],
                             "RNA MW1 L2FC" = lungl2fcmat[ourgene,"M W1"],
                             "RNA MW2 L2FC" = lungl2fcmat[ourgene,"M W2"],
                             "RNA MW4 L2FC" = lungl2fcmat[ourgene,"M W4"],
                             "RNA MW8 L2FC" = lungl2fcmat[ourgene,"M W8"],
                             "ATAC FW1 L2FC" = lungsigatacl2fc[ourpeak,"F W1"],
                             "ATAC FW2 L2FC" = lungsigatacl2fc[ourpeak,"F W2"],
                             "ATAC FW4 L2FC" = lungsigatacl2fc[ourpeak,"F W4"],
                             "ATAC FW8 L2FC" = lungsigatacl2fc[ourpeak,"F W8"],
                             "ATAC MW1 L2FC" = lungsigatacl2fc[ourpeak,"M W1"],
                             "ATAC MW2 L2FC" = lungsigatacl2fc[ourpeak,"M W2"],
                             "ATAC MW4 L2FC" = lungsigatacl2fc[ourpeak,"M W4"],
                             "ATAC MW8 L2FC" = lungsigatacl2fc[ourpeak,"M W8"])
    dar_deg_connection50df <- rbind(dar_deg_connection50df,ourdfentry)
  } else if(length(ourpeaks) > 1){
    for(j in 1:length(ourpeaks)){
      ourpeak <- ourpeaks[j]
      ourdfentry <- data.frame("Tissue" = "LUNG",
                               "DEG" = ourgene,
                               "SYMBOL" = enstosym[ourgene,"Symbol"],
                               "DAR" = ourpeak,
                               "Distance" = abs(lungrnasiganno[ourgene,"start"] - ((peakanno[ourpeak,"start"]+peakanno[ourpeak,"end"])/2)),
                               "Correlation" = lung50peakcorfin[ourpeak,ourgene],
                               "TF Motifs" = toString(unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% ourpeak,"Motif.Name"])),
                               "RNA FW1 L2FC" = lungl2fcmat[ourgene,"F W1"],
                               "RNA FW2 L2FC" = lungl2fcmat[ourgene,"F W2"],
                               "RNA FW4 L2FC" = lungl2fcmat[ourgene,"F W4"],
                               "RNA FW8 L2FC" = lungl2fcmat[ourgene,"F W8"],
                               "RNA MW1 L2FC" = lungl2fcmat[ourgene,"M W1"],
                               "RNA MW2 L2FC" = lungl2fcmat[ourgene,"M W2"],
                               "RNA MW4 L2FC" = lungl2fcmat[ourgene,"M W4"],
                               "RNA MW8 L2FC" = lungl2fcmat[ourgene,"M W8"],
                               "ATAC FW1 L2FC" = lungsigatacl2fc[ourpeak,"F W1"],
                               "ATAC FW2 L2FC" = lungsigatacl2fc[ourpeak,"F W2"],
                               "ATAC FW4 L2FC" = lungsigatacl2fc[ourpeak,"F W4"],
                               "ATAC FW8 L2FC" = lungsigatacl2fc[ourpeak,"F W8"],
                               "ATAC MW1 L2FC" = lungsigatacl2fc[ourpeak,"M W1"],
                               "ATAC MW2 L2FC" = lungsigatacl2fc[ourpeak,"M W2"],
                               "ATAC MW4 L2FC" = lungsigatacl2fc[ourpeak,"M W4"],
                               "ATAC MW8 L2FC" = lungsigatacl2fc[ourpeak,"M W8"])
      dar_deg_connection50df <- rbind(dar_deg_connection50df,ourdfentry)
      
    }
  } else {
    print("found nothing")
  }
}

dar_deg_connection50df <- dar_deg_connection50df[2:dim(dar_deg_connection50df)[1],]

# Make supplemental table - remove empty row 1 and numbered list in column 1
write.csv(dar_deg_connection50df,file = "Supplemental Table S2_021325.csv",row.names = F)

# Generate a metadata df based on pid

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


ann_cols <- list("Tissue" = c("SKM-GN" = "#088c03",
                              "HEART" = "#f28b2f",
                              "HIPPOC" = "#bf7534",
                              "KIDNEY"= "#7553a7",
                              "LIVER" = "#da6c75",
                              "LUNG" = "#04bf8a",
                              "BAT" = "#8c5220",
                              "WAT-SC" = "#214da6"),
                 "Sex" = c("Female" = "#ff6eff",
                           "Male" = "#5555ff"),
                 "Group" = c("control" = "grey",
                             "1w" = "#F7FCB9",
                             "2w" = "#ADDD8E",
                             "4w" = "#238443",
                             "8w" = "#002612"),
                 "Region" = c("3'.UTR" = "#E377C2FF",
                              "5'.UTR" = "#D62728FF",
                              "Distal.Intergenic" = "#BCBD22FF",
                              "Downstream" = "#7F7F7FFF",
                              "Exon" = "#9467BDFF",
                              "Intron" = "#8C564BFF",
                              "Promoter.(1-2kb)" = "#2CA02CFF",
                              "Promoter.(<=1kb)" = "#FF7F0EFF",
                              "Upstream" = "#1F77B4FF"))

ourgene <- enstosym[enstosym$Symbol %in% "Srsf2","Ensembl"]
srsf2_gastro_corpeaks <- rownames(gastropeakcortrim)[abs(gastropeakcortrim[,ourgene]) > 0]
ourpeak <- srsf2_gastro_corpeaks[1]
ourdf <- data.frame("RNA" = as.vector(t(gastrornatrainingl2fcbysubject[ourgene,])),
                    "ATAC" = as.vector(t(gastroatactrainingl2fcbysubject[ourpeak,])),
                    "Group" = pidmetadf[colnames(gastrornatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(gastrornatrainingl2fcbysubject),"Sex"])
pdf(file = "Supplemental Figure S18Q_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()

png(file = "Supplemental Figure S18Q_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()



ourgene <- enstosym[enstosym$Symbol %in% "Mpeg1","Ensembl"]
mpeg1_lung_corpeaks <- rownames(lungpeakcortrim)[abs(lungpeakcortrim[,ourgene]) > 0]
ourpeak <- mpeg1_lung_corpeaks[1]
ourdf <- data.frame("RNA" = as.vector(t(lungrnatrainingl2fcbysubject[ourgene,])),
                    "ATAC" = as.vector(t(lungatactrainingl2fcbysubject[ourpeak,])),
                    "Group" = pidmetadf[colnames(lungrnatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(lungrnatrainingl2fcbysubject),"Sex"])
pdf(file = "Figure 4E_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LUNG ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()

png(file = "Figure 4E_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LUNG ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()


ourgene <- enstosym[enstosym$Symbol %in% "Abhd2","Ensembl"]
abhd2_liver_corpeaks <- rownames(liverpeakcortrim)[abs(liverpeakcortrim[,ourgene]) > 0]
ourpeak <- abhd2_liver_corpeaks[7]
ourdf <- data.frame("RNA" = as.vector(t(liverrnatrainingl2fcbysubject[ourgene,])),
                    "ATAC" = as.vector(t(liveratactrainingl2fcbysubject[ourpeak,])),
                    "Group" = pidmetadf[colnames(liverrnatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(liverrnatrainingl2fcbysubject),"Sex"])
pdf(file = "Figure 4K_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LIVER ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()

png(file = "Figure 4K_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LIVER ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Abhd2","Ensembl"]
abhd2_liver_corpeaks <- rownames(liverpeakcortrim)[abs(liverpeakcortrim[,ourgene]) > 0]
ourpeak <- abhd2_liver_corpeaks[6]
ourdf <- data.frame("RNA" = as.vector(t(liverrnatrainingl2fcbysubject[ourgene,])),
                    "ATAC" = as.vector(t(liveratactrainingl2fcbysubject[ourpeak,])),
                    "Group" = pidmetadf[colnames(liverrnatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(liverrnatrainingl2fcbysubject),"Sex"])
pdf(file = "Figure 4L_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LIVER ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()

png(file = "Figure 4L_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LIVER ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()



ourgene <- enstosym[enstosym$Symbol %in% "Onecut1","Ensembl"]
onecut1_liver_corpeaks <- rownames(liverpeakcortrim)[abs(liverpeakcortrim[,ourgene]) > 0]
ourpeak <- onecut1_liver_corpeaks[1]
ourdf <- data.frame("RNA" = as.vector(t(liverrnatrainingl2fcbysubject[ourgene,])),
                    "ATAC" = as.vector(t(liveratactrainingl2fcbysubject[ourpeak,])),
                    "Group" = pidmetadf[colnames(liverrnatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(liverrnatrainingl2fcbysubject),"Sex"])
pdf(file = "Supplemental Figure S18K_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LIVER ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()

png(file = "Supplemental Figure S18K_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LIVER ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()



ourgene <- enstosym[enstosym$Symbol %in% "Kcna7","Ensembl"]
kcna7_gastro_corpeaks <- rownames(gastropeakcortrim)[abs(gastropeakcortrim[,ourgene]) > 0]
ourpeak <- kcna7_gastro_corpeaks[1]
ourdf <- data.frame("RNA" = as.vector(t(gastrornatrainingl2fcbysubject[ourgene,])),
                    "ATAC" = as.vector(t(gastroatactrainingl2fcbysubject[ourpeak,])),
                    "Group" = pidmetadf[colnames(gastrornatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(gastrornatrainingl2fcbysubject),"Sex"])
pdf(file = "Supplemental Figure S18O_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()
png(file = "Supplemental Figure S18O_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Eif4g3","Ensembl"]
eif4g3_gastro_corpeaks <- rownames(gastropeakcortrim)[abs(gastropeakcortrim[,ourgene]) > 0]
ourpeak <- eif4g3_gastro_corpeaks[1]
ourdf <- data.frame("RNA" = as.vector(t(gastrornatrainingl2fcbysubject[ourgene,])),
                    "ATAC" = as.vector(t(gastroatactrainingl2fcbysubject[ourpeak,])),
                    "Group" = pidmetadf[colnames(gastrornatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(gastrornatrainingl2fcbysubject),"Sex"])
pdf(file = "Supplemental Figure S18P_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()
png(file = "Supplemental Figure S18P_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Oas2","Ensembl"]
oas2_lung_corpeaks <- rownames(lungpeakcortrim)[abs(lungpeakcortrim[,ourgene]) > 0]
ourpeak <- oas2_lung_corpeaks[1]
ourdf <- data.frame("RNA" = as.vector(t(lungrnatrainingl2fcbysubject[ourgene,])),
                    "ATAC" = as.vector(t(lungatactrainingl2fcbysubject[ourpeak,])),
                    "Group" = pidmetadf[colnames(lungrnatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(lungrnatrainingl2fcbysubject),"Sex"])
pdf(file = "Supplemental Figure S18B_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LUNG ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()
png(file = "Supplemental Figure S18B_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LUNG ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()


ourgene <- enstosym[enstosym$Symbol %in% "Gm2a","Ensembl"]
gm2a_lung_corpeaks <- rownames(lungpeakcortrim)[abs(lungpeakcortrim[,ourgene]) > 0]
ourpeak <- gm2a_lung_corpeaks[1]
ourdf <- data.frame("RNA" = as.vector(t(lungrnatrainingl2fcbysubject[ourgene,])),
                    "ATAC" = as.vector(t(lungatactrainingl2fcbysubject[ourpeak,])),
                    "Group" = pidmetadf[colnames(lungrnatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(lungrnatrainingl2fcbysubject),"Sex"])
pdf(file = "Figure 4F_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LUNG ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()
png(file = "Figure 4F_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LUNG ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Dennd1c","Ensembl"]
dennd1c_lung_corpeaks <- rownames(lungpeakcortrim)[abs(lungpeakcortrim[,ourgene]) > 0]
ourpeak <- dennd1c_lung_corpeaks[1]
ourdf <- data.frame("RNA" = as.vector(t(lungrnatrainingl2fcbysubject[ourgene,])),
                    "ATAC" = as.vector(t(lungatactrainingl2fcbysubject[ourpeak,])),
                    "Group" = pidmetadf[colnames(lungrnatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(lungrnatrainingl2fcbysubject),"Sex"])
pdf(file = "Supplemental Figure S18A_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LUNG ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()
png(file = "Supplemental Figure S18A_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LUNG ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()


ourgene <- enstosym[enstosym$Symbol %in% "Bcl11b","Ensembl"]
bcl11b_lung_corpeaks <- rownames(lungpeakcortrim)[abs(lungpeakcortrim[,ourgene]) > 0]
ourpeak <- bcl11b_lung_corpeaks[1]
ourdf <- data.frame("RNA" = as.vector(t(lungrnatrainingl2fcbysubject[ourgene,])),
                    "ATAC" = as.vector(t(lungatactrainingl2fcbysubject[ourpeak,])),
                    "Group" = pidmetadf[colnames(lungrnatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(lungrnatrainingl2fcbysubject),"Sex"])
pdf(file = "Figure 4G_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LUNG ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()
png(file = "Figure 4G_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LUNG ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()


ourgene <- enstosym[enstosym$Symbol %in% "Ppard","Ensembl"]
ppard_gastro_corpeaks <- rownames(gastropeakcortrim)[abs(gastropeakcortrim[,ourgene]) > 0]
ourpeak <- ppard_gastro_corpeaks[1]
ourdf <- data.frame("RNA" = as.vector(t(gastrornatrainingl2fcbysubject[ourgene,])),
                    "ATAC" = as.vector(t(gastroatactrainingl2fcbysubject[ourpeak,])),
                    "Group" = pidmetadf[colnames(gastrornatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(gastrornatrainingl2fcbysubject),"Sex"])
pdf(file = "Figure 4J_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()
png(file = "Figure 4J_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()


ourgene <- enstosym[enstosym$Symbol %in% "Hjv","Ensembl"]
hjv_gastro_corpeaks <- rownames(gastropeakcortrim)[abs(gastropeakcortrim[,ourgene]) > 0]
ourpeak <- hjv_gastro_corpeaks[3]
ourdf <- data.frame("RNA" = as.vector(t(gastrornatrainingl2fcbysubject[ourgene,])),
                    "ATAC" = as.vector(t(gastroatactrainingl2fcbysubject[ourpeak,])),
                    "Group" = pidmetadf[colnames(gastrornatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(gastrornatrainingl2fcbysubject),"Sex"])
pdf(file = "Figure 4H_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()
png(file = "Figure 4H_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()


ourgene <- enstosym[enstosym$Symbol %in% "Fkbp5","Ensembl"]
fkbp5_gastro_corpeaks <- rownames(gastropeakcortrim)[abs(gastropeakcortrim[,ourgene]) > 0]
ourpeak <- fkbp5_gastro_corpeaks[1]
ourdf <- data.frame("RNA" = as.vector(t(gastrornatrainingl2fcbysubject[ourgene,])),
                    "ATAC" = as.vector(t(gastroatactrainingl2fcbysubject[ourpeak,])),
                    "Group" = pidmetadf[colnames(gastrornatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(gastrornatrainingl2fcbysubject),"Sex"])
pdf(file = "Figure 4I_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()
png(file = "Figure 4I_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Cited4","Ensembl"]
cited4_gastro_corpeaks <- rownames(gastropeakcortrim)[abs(gastropeakcortrim[,ourgene]) > 0]
ourpeak <- cited4_gastro_corpeaks[1]
ourdf <- data.frame("RNA" = as.vector(t(gastrornatrainingl2fcbysubject[ourgene,])),
                    "ATAC" = as.vector(t(gastroatactrainingl2fcbysubject[ourpeak,])),
                    "Group" = pidmetadf[colnames(gastrornatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(gastrornatrainingl2fcbysubject),"Sex"])
pdf(file = "Supplemental Figure S18C_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()
png(file = "Supplemental Figure S18C_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()


ourgene <- enstosym[enstosym$Symbol %in% "Cab39","Ensembl"]
cab39_gastro_corpeaks <- rownames(gastropeakcortrim)[abs(gastropeakcortrim[,ourgene]) > 0]
ourpeak <- cab39_gastro_corpeaks[1]
ourdf <- data.frame("RNA" = as.vector(t(gastrornatrainingl2fcbysubject[ourgene,])),
                    "ATAC" = as.vector(t(gastroatactrainingl2fcbysubject[ourpeak,])),
                    "Group" = pidmetadf[colnames(gastrornatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(gastrornatrainingl2fcbysubject),"Sex"])
pdf(file = "Supplemental Figure S18D_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()
png(file = "Supplemental Figure S18D_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()



ourgene <- enstosym[enstosym$Symbol %in% "Sik1","Ensembl"]
sik1_gastro_corpeaks <- rownames(gastropeakcortrim)[abs(gastropeakcortrim[,ourgene]) > 0]
ourpeak <- sik1_gastro_corpeaks[1]
ourdf <- data.frame("RNA" = as.vector(t(gastrornatrainingl2fcbysubject[ourgene,])),
                    "ATAC" = as.vector(t(gastroatactrainingl2fcbysubject[ourpeak,])),
                    "Group" = pidmetadf[colnames(gastrornatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(gastrornatrainingl2fcbysubject),"Sex"])
pdf(file = "Supplemental Figure S18E_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()
png(file = "Supplemental Figure S18E_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()


ourgene <- enstosym[enstosym$Symbol %in% "Phlda3","Ensembl"]
phlda3_gastro_corpeaks <- rownames(gastropeakcortrim)[abs(gastropeakcortrim[,ourgene]) > 0]
ourpeak <- phlda3_gastro_corpeaks[1]
ourdf <- data.frame("RNA" = as.vector(t(gastrornatrainingl2fcbysubject[ourgene,])),
                    "ATAC" = as.vector(t(gastroatactrainingl2fcbysubject[ourpeak,])),
                    "Group" = pidmetadf[colnames(gastrornatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(gastrornatrainingl2fcbysubject),"Sex"])
pdf(file = "Supplemental Figure S18F_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()
png(file = "Supplemental Figure S18F_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()


ourgene <- enstosym[enstosym$Symbol %in% "Hspb6","Ensembl"]
hspb6_lung_corpeaks <- rownames(lungpeakcortrim)[abs(lungpeakcortrim[,ourgene]) > 0]
ourpeak <- hspb6_lung_corpeaks[1]
ourdf <- data.frame("RNA" = as.vector(t(lungrnatrainingl2fcbysubject[ourgene,])),
                    "ATAC" = as.vector(t(lungatactrainingl2fcbysubject[ourpeak,])),
                    "Group" = pidmetadf[colnames(lungrnatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(lungrnatrainingl2fcbysubject),"Sex"])
pdf(file = "Supplemental Figure S18R_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LUNG ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()
png(file = "Supplemental Figure S18R_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LUNG ",enstosym[ourgene,"Symbol"],": Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()



ourgene <- enstosym[enstosym$Symbol %in% "G6pc1","Ensembl"]
g6pc1_liver_corpeaks <- rownames(liverpeakcortrim)[abs(liverpeakcortrim[,ourgene]) > 0]
ourpeak <- g6pc1_liver_corpeaks[6]
ourdf <- data.frame("RNA" = as.vector(t(liverrnatrainingl2fcbysubject[ourgene,])),
                    "ATAC" = as.vector(t(liveratactrainingl2fcbysubject[ourpeak,])),
                    "Group" = pidmetadf[colnames(liverrnatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(liverrnatrainingl2fcbysubject),"Sex"])
pdf(file = "Supplemental Figure S18G_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LIVER ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()

png(file = "Supplemental Figure S18G_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LIVER ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()


ourpeak <- g6pc1_liver_corpeaks[10]
ourdf <- data.frame("RNA" = as.vector(t(liverrnatrainingl2fcbysubject[ourgene,])),
                    "ATAC" = as.vector(t(liveratactrainingl2fcbysubject[ourpeak,])),
                    "Group" = pidmetadf[colnames(liverrnatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(liverrnatrainingl2fcbysubject),"Sex"])
pdf(file = "Supplemental Figure S18H_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LIVER ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()

png(file = "Supplemental Figure S18H_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LIVER ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()

ourpeak <- g6pc1_liver_corpeaks[5]
ourdf <- data.frame("RNA" = as.vector(t(liverrnatrainingl2fcbysubject[ourgene,])),
                    "ATAC" = as.vector(t(liveratactrainingl2fcbysubject[ourpeak,])),
                    "Group" = pidmetadf[colnames(liverrnatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(liverrnatrainingl2fcbysubject),"Sex"])
pdf(file = "Supplemental Figure S18I_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LIVER ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()

png(file = "Supplemental Figure S18I_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LIVER ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()



ourgene <- enstosym[enstosym$Symbol %in% "Chp1","Ensembl"]
chp1_liver_corpeaks <- rownames(liverpeakcortrim)[abs(liverpeakcortrim[,ourgene]) > 0]
ourpeak <- chp1_liver_corpeaks
ourdf <- data.frame("RNA" = as.vector(t(liverrnatrainingl2fcbysubject[ourgene,])),
                    "ATAC" = as.vector(t(liveratactrainingl2fcbysubject[ourpeak,])),
                    "Group" = pidmetadf[colnames(liverrnatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(liverrnatrainingl2fcbysubject),"Sex"])
pdf(file = "Supplemental Figure S18J_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LIVER ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()

png(file = "Supplemental Figure S18J_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LIVER ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()


ourgene <- enstosym[enstosym$Symbol %in% "Lpar3","Ensembl"]
lpar3_liver_corpeaks <- rownames(liverpeakcortrim)[abs(liverpeakcortrim[,ourgene]) > 0]
ourpeak <- lpar3_liver_corpeaks[2]
ourdf <- data.frame("RNA" = as.vector(t(liverrnatrainingl2fcbysubject[ourgene,])),
                    "ATAC" = as.vector(t(liveratactrainingl2fcbysubject[ourpeak,])),
                    "Group" = pidmetadf[colnames(liverrnatrainingl2fcbysubject),"Group"],
                    "Sex" = pidmetadf[colnames(liverrnatrainingl2fcbysubject),"Sex"])
pdf(file = "Supplemental Figure S18L_021325.pdf",width = 6,height = 6)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LIVER ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()

png(file = "Supplemental Figure S18L_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LIVER ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC L2FC") + xlab("RNA L2FC") + scale_color_manual(values = ann_cols$Group)
dev.off()


######
# I'm going to recreate Figure 4H for the DAR-DEG connections


####
# gastro

gastropeakcortrimtfs <- gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin),]
gastropeakcortfs <- matrix(0L,nrow = dim(gastro50peakcorfin)[1],ncol = length(unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin),"Motif.Name"])))
rownames(gastropeakcortfs) <- rownames(gastro50peakcorfin)
colnames(gastropeakcortfs) <- unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin),"Motif.Name"])
for(i in 1:dim(gastropeakcortrimtfs)[1]){
  gastropeakcortfs[gastropeakcortrimtfs$PositionID[i],gastropeakcortrimtfs$Motif.Name[i]] <- 1
}

gastrocortfsig <- matrix(0L,nrow = dim(gastropeakcortfs)[2],ncol = 5)
rownames(gastrocortfsig) <- colnames(gastropeakcortfs)
colnames(gastrocortfsig) <- c("Correlated_Target_Count","Total_Target_Count","Fraction_Corr","Fraction_Total","Enrichment_Significance")
gastrototcor <- length(rownames(gastropeakcortfs))
gastrototpeaks <- length(unique(gastro50peakmotifs$PositionID))
for(i in 1:dim(gastrocortfsig)[1]){
  if(i%%20 == 0){
    print(i)
  }
  ourtf <- rownames(gastrocortfsig)[i]
  gastrocortfsig[i,"Correlated_Target_Count"] <- sum(gastropeakcortfs[,i])
  gastrocortfsig[i,"Total_Target_Count"] <- length(unique(gastro50peakmotifs[gastro50peakmotifs$Motif.Name %in% ourtf,"PositionID"]))
  gastrocortfsig[i,"Fraction_Corr"] <- gastrocortfsig[i,"Correlated_Target_Count"]/gastrototcor
  gastrocortfsig[i,"Fraction_Total"] <- gastrocortfsig[i,"Total_Target_Count"]/gastrototpeaks
  gastrocortfsig[i,"Enrichment_Significance"] <- binom.test(x = gastrocortfsig[i,"Correlated_Target_Count"], n = gastrototcor, p = gastrocortfsig[i,"Total_Target_Count"]/gastrototpeaks,alternative = "greater")$p.value
}


####
# kidney

kidneypeakcortrimtfs <- kidney50peakmotifs[kidney50peakmotifs$PositionID %in% rownames(kidney50peakcorfin),]
kidneypeakcortfs <- matrix(0L,nrow = dim(kidney50peakcorfin)[1],ncol = length(unique(kidney50peakmotifs[kidney50peakmotifs$PositionID %in% rownames(kidney50peakcorfin),"Motif.Name"])))
rownames(kidneypeakcortfs) <- rownames(kidney50peakcorfin)
colnames(kidneypeakcortfs) <- unique(kidney50peakmotifs[kidney50peakmotifs$PositionID %in% rownames(kidney50peakcorfin),"Motif.Name"])
for(i in 1:dim(kidneypeakcortrimtfs)[1]){
  kidneypeakcortfs[kidneypeakcortrimtfs$PositionID[i],kidneypeakcortrimtfs$Motif.Name[i]] <- 1
}

kidneycortfsig <- matrix(0L,nrow = dim(kidneypeakcortfs)[2],ncol = 5)
rownames(kidneycortfsig) <- colnames(kidneypeakcortfs)
colnames(kidneycortfsig) <- c("Correlated_Target_Count","Total_Target_Count","Fraction_Corr","Fraction_Total","Enrichment_Significance")
kidneytotcor <- length(rownames(kidneypeakcortfs))
kidneytotpeaks <- length(unique(kidney50peakmotifs$PositionID))
for(i in 1:dim(kidneycortfsig)[1]){
  if(i%%20 == 0){
    print(i)
  }
  ourtf <- rownames(kidneycortfsig)[i]
  kidneycortfsig[i,"Correlated_Target_Count"] <- sum(kidneypeakcortfs[,i])
  kidneycortfsig[i,"Total_Target_Count"] <- length(unique(kidney50peakmotifs[kidney50peakmotifs$Motif.Name %in% ourtf,"PositionID"]))
  kidneycortfsig[i,"Fraction_Corr"] <- kidneycortfsig[i,"Correlated_Target_Count"]/kidneytotcor
  kidneycortfsig[i,"Fraction_Total"] <- kidneycortfsig[i,"Total_Target_Count"]/kidneytotpeaks
  kidneycortfsig[i,"Enrichment_Significance"] <- binom.test(x = kidneycortfsig[i,"Correlated_Target_Count"], n = kidneytotcor, p = kidneycortfsig[i,"Total_Target_Count"]/kidneytotpeaks,alternative = "greater")$p.value
}



####
# liver

liverpeakcortrimtfs <- liver50peakmotifs[liver50peakmotifs$PositionID %in% rownames(liver50peakcorfin),]
liverpeakcortfs <- matrix(0L,nrow = dim(liver50peakcorfin)[1],ncol = length(unique(liver50peakmotifs[liver50peakmotifs$PositionID %in% rownames(liver50peakcorfin),"Motif.Name"])))
rownames(liverpeakcortfs) <- rownames(liver50peakcorfin)
colnames(liverpeakcortfs) <- unique(liver50peakmotifs[liver50peakmotifs$PositionID %in% rownames(liver50peakcorfin),"Motif.Name"])
for(i in 1:dim(liverpeakcortrimtfs)[1]){
  liverpeakcortfs[liverpeakcortrimtfs$PositionID[i],liverpeakcortrimtfs$Motif.Name[i]] <- 1
}

livercortfsig <- matrix(0L,nrow = dim(liverpeakcortfs)[2],ncol = 5)
rownames(livercortfsig) <- colnames(liverpeakcortfs)
colnames(livercortfsig) <- c("Correlated_Target_Count","Total_Target_Count","Fraction_Corr","Fraction_Total","Enrichment_Significance")
livertotcor <- length(rownames(liverpeakcortfs))
livertotpeaks <- length(unique(liver50peakmotifs$PositionID))
for(i in 1:dim(livercortfsig)[1]){
  if(i%%20 == 0){
    print(i)
  }
  ourtf <- rownames(livercortfsig)[i]
  livercortfsig[i,"Correlated_Target_Count"] <- sum(liverpeakcortfs[,i])
  livercortfsig[i,"Total_Target_Count"] <- length(unique(liver50peakmotifs[liver50peakmotifs$Motif.Name %in% ourtf,"PositionID"]))
  livercortfsig[i,"Fraction_Corr"] <- livercortfsig[i,"Correlated_Target_Count"]/livertotcor
  livercortfsig[i,"Fraction_Total"] <- livercortfsig[i,"Total_Target_Count"]/livertotpeaks
  livercortfsig[i,"Enrichment_Significance"] <- binom.test(x = livercortfsig[i,"Correlated_Target_Count"], n = livertotcor, p = livercortfsig[i,"Total_Target_Count"]/livertotpeaks,alternative = "greater")$p.value
}



####
# lung

lungpeakcortrimtfs <- lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin),]
lungpeakcortfs <- matrix(0L,nrow = dim(lung50peakcorfin)[1],ncol = length(unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin),"Motif.Name"])))
rownames(lungpeakcortfs) <- rownames(lung50peakcorfin)
colnames(lungpeakcortfs) <- unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin),"Motif.Name"])
for(i in 1:dim(lungpeakcortrimtfs)[1]){
  lungpeakcortfs[lungpeakcortrimtfs$PositionID[i],lungpeakcortrimtfs$Motif.Name[i]] <- 1
}

lungcortfsig <- matrix(0L,nrow = dim(lungpeakcortfs)[2],ncol = 5)
rownames(lungcortfsig) <- colnames(lungpeakcortfs)
colnames(lungcortfsig) <- c("Correlated_Target_Count","Total_Target_Count","Fraction_Corr","Fraction_Total","Enrichment_Significance")
lungtotcor <- length(rownames(lungpeakcortfs))
lungtotpeaks <- length(unique(lung50peakmotifs$PositionID))
for(i in 1:dim(lungcortfsig)[1]){
  if(i%%20 == 0){
    print(i)
  }
  ourtf <- rownames(lungcortfsig)[i]
  lungcortfsig[i,"Correlated_Target_Count"] <- sum(lungpeakcortfs[,i])
  lungcortfsig[i,"Total_Target_Count"] <- length(unique(lung50peakmotifs[lung50peakmotifs$Motif.Name %in% ourtf,"PositionID"]))
  lungcortfsig[i,"Fraction_Corr"] <- lungcortfsig[i,"Correlated_Target_Count"]/lungtotcor
  lungcortfsig[i,"Fraction_Total"] <- lungcortfsig[i,"Total_Target_Count"]/lungtotpeaks
  lungcortfsig[i,"Enrichment_Significance"] <- binom.test(x = lungcortfsig[i,"Correlated_Target_Count"], n = lungtotcor, p = lungcortfsig[i,"Total_Target_Count"]/lungtotpeaks,alternative = "greater")$p.value
}



####
# brown

brownpeakcortrimtfs <- brown50peakmotifs[brown50peakmotifs$PositionID %in% rownames(brown50peakcorfin),]
brownpeakcortfs <- matrix(0L,nrow = dim(brown50peakcorfin)[1],ncol = length(unique(brown50peakmotifs[brown50peakmotifs$PositionID %in% rownames(brown50peakcorfin),"Motif.Name"])))
rownames(brownpeakcortfs) <- rownames(brown50peakcorfin)
colnames(brownpeakcortfs) <- unique(brown50peakmotifs[brown50peakmotifs$PositionID %in% rownames(brown50peakcorfin),"Motif.Name"])
for(i in 1:dim(brownpeakcortrimtfs)[1]){
  brownpeakcortfs[brownpeakcortrimtfs$PositionID[i],brownpeakcortrimtfs$Motif.Name[i]] <- 1
}

browncortfsig <- matrix(0L,nrow = dim(brownpeakcortfs)[2],ncol = 5)
rownames(browncortfsig) <- colnames(brownpeakcortfs)
colnames(browncortfsig) <- c("Correlated_Target_Count","Total_Target_Count","Fraction_Corr","Fraction_Total","Enrichment_Significance")
browntotcor <- length(rownames(brownpeakcortfs))
browntotpeaks <- length(unique(brown50peakmotifs$PositionID))
for(i in 1:dim(browncortfsig)[1]){
  if(i%%20 == 0){
    print(i)
  }
  ourtf <- rownames(browncortfsig)[i]
  browncortfsig[i,"Correlated_Target_Count"] <- sum(brownpeakcortfs[,i])
  browncortfsig[i,"Total_Target_Count"] <- length(unique(brown50peakmotifs[brown50peakmotifs$Motif.Name %in% ourtf,"PositionID"]))
  browncortfsig[i,"Fraction_Corr"] <- browncortfsig[i,"Correlated_Target_Count"]/browntotcor
  browncortfsig[i,"Fraction_Total"] <- browncortfsig[i,"Total_Target_Count"]/browntotpeaks
  browncortfsig[i,"Enrichment_Significance"] <- binom.test(x = browncortfsig[i,"Correlated_Target_Count"], n = browntotcor, p = browncortfsig[i,"Total_Target_Count"]/browntotpeaks,alternative = "greater")$p.value
}




cortffractiondf <- data.frame("Enrichment.Fraction" = c(gastrocortfsig[gastrocortfsig[,"Enrichment_Significance"] < 0.05 & gastrocortfsig[,"Correlated_Target_Count"] > 3,"Fraction_Corr"],
                                                        livercortfsig[livercortfsig[,"Enrichment_Significance"] < 0.05 & livercortfsig[,"Correlated_Target_Count"] > 3,"Fraction_Corr"],
                                                        lungcortfsig[lungcortfsig[,"Enrichment_Significance"] < 0.05 & lungcortfsig[,"Correlated_Target_Count"] > 3,"Fraction_Corr"],
                                                        gastrocortfsig[gastrocortfsig[,"Enrichment_Significance"] < 0.05 & gastrocortfsig[,"Correlated_Target_Count"] > 3,"Fraction_Total"],
                                                        livercortfsig[livercortfsig[,"Enrichment_Significance"] < 0.05 & livercortfsig[,"Correlated_Target_Count"] > 3,"Fraction_Total"],
                                                        lungcortfsig[lungcortfsig[,"Enrichment_Significance"] < 0.05 & lungcortfsig[,"Correlated_Target_Count"] > 3,"Fraction_Total"]),
                              "P.Value" = c(gastrocortfsig[gastrocortfsig[,"Enrichment_Significance"] < 0.05 & gastrocortfsig[,"Correlated_Target_Count"] > 3,"Enrichment_Significance"],
                                            livercortfsig[livercortfsig[,"Enrichment_Significance"] < 0.05 & livercortfsig[,"Correlated_Target_Count"] > 3,"Enrichment_Significance"],
                                            lungcortfsig[lungcortfsig[,"Enrichment_Significance"] < 0.05 & lungcortfsig[,"Correlated_Target_Count"] > 3,"Enrichment_Significance"],
                                            gastrocortfsig[gastrocortfsig[,"Enrichment_Significance"] < 0.05 & gastrocortfsig[,"Correlated_Target_Count"] > 3,"Enrichment_Significance"],
                                            livercortfsig[livercortfsig[,"Enrichment_Significance"] < 0.05 & livercortfsig[,"Correlated_Target_Count"] > 3,"Enrichment_Significance"],
                                            lungcortfsig[lungcortfsig[,"Enrichment_Significance"] < 0.05 & lungcortfsig[,"Correlated_Target_Count"] > 3,"Enrichment_Significance"]),
                              "Target.List" = c(rep("Correlated DEG-DAR Pairs",(length(gastrocortfsig[gastrocortfsig[,"Enrichment_Significance"] < 0.05 & gastrocortfsig[,"Correlated_Target_Count"] > 3,"Fraction_Corr"])+
                                                                                  length(livercortfsig[livercortfsig[,"Enrichment_Significance"] < 0.05 & livercortfsig[,"Correlated_Target_Count"] > 3,"Fraction_Corr"])+
                                                                                  length(lungcortfsig[lungcortfsig[,"Enrichment_Significance"] < 0.05 & lungcortfsig[,"Correlated_Target_Count"] > 3,"Fraction_Corr"]))),
                                                rep("All Peaks",(length(gastrocortfsig[gastrocortfsig[,"Enrichment_Significance"] < 0.05 & gastrocortfsig[,"Correlated_Target_Count"] > 3,"Fraction_Corr"])+
                                                                               length(livercortfsig[livercortfsig[,"Enrichment_Significance"] < 0.05 & livercortfsig[,"Correlated_Target_Count"] > 3,"Fraction_Corr"])+
                                                                               length(lungcortfsig[lungcortfsig[,"Enrichment_Significance"] < 0.05 & lungcortfsig[,"Correlated_Target_Count"] > 3,"Fraction_Corr"])))),
                              "Tissue" = c(rep("SKM-GN",length(gastrocortfsig[gastrocortfsig[,"Enrichment_Significance"] < 0.05 & gastrocortfsig[,"Correlated_Target_Count"] > 3,"Fraction_Corr"])),
                                           rep("LIVER",length(livercortfsig[livercortfsig[,"Enrichment_Significance"] < 0.05 & livercortfsig[,"Correlated_Target_Count"] > 3,"Fraction_Corr"])),
                                           rep("LUNG",length(lungcortfsig[lungcortfsig[,"Enrichment_Significance"] < 0.05 & lungcortfsig[,"Correlated_Target_Count"] > 3,"Fraction_Corr"])),
                                           rep("SKM-GN",length(gastrocortfsig[gastrocortfsig[,"Enrichment_Significance"] < 0.05 & gastrocortfsig[,"Correlated_Target_Count"] > 3,"Fraction_Corr"])),
                                           rep("LIVER",length(livercortfsig[livercortfsig[,"Enrichment_Significance"] < 0.05 & livercortfsig[,"Correlated_Target_Count"] > 3,"Fraction_Corr"])),
                                           rep("LUNG",length(lungcortfsig[lungcortfsig[,"Enrichment_Significance"] < 0.05 & lungcortfsig[,"Correlated_Target_Count"] > 3,"Fraction_Corr"]))),
                              "TF" = c(rownames(gastrocortfsig)[gastrocortfsig[,"Enrichment_Significance"] < 0.05 & gastrocortfsig[,"Correlated_Target_Count"] > 3],
                                       rownames(livercortfsig)[livercortfsig[,"Enrichment_Significance"] < 0.05 & livercortfsig[,"Correlated_Target_Count"] > 3],
                                       rownames(lungcortfsig)[lungcortfsig[,"Enrichment_Significance"] < 0.05 & lungcortfsig[,"Correlated_Target_Count"] > 3],
                                       rownames(gastrocortfsig)[gastrocortfsig[,"Enrichment_Significance"] < 0.05 & gastrocortfsig[,"Correlated_Target_Count"] > 3],
                                       rownames(livercortfsig)[livercortfsig[,"Enrichment_Significance"] < 0.05 & livercortfsig[,"Correlated_Target_Count"] > 3],
                                       rownames(lungcortfsig)[lungcortfsig[,"Enrichment_Significance"] < 0.05 & lungcortfsig[,"Correlated_Target_Count"] > 3]))
cortffractiondf$TF <- gsub("\\(.*","",cortffractiondf$TF)
cortffractiondf$Target.List <- factor(cortffractiondf$Target.List,levels = c("Correlated DEG-DAR Pairs","All Peaks"))
cortffractiondf <- cortffractiondf[c(order(cortffractiondf$Enrichment.Fraction[1:(dim(cortffractiondf)[1]/2)]),(dim(cortffractiondf)[1]/2)+order(cortffractiondf$Enrichment.Fraction[1:(dim(cortffractiondf)[1]/2)])),]
cortffractiondf$TF.Number <- as.character(rep(c(1:(dim(cortffractiondf)[1]/2)),2))
cortffractiondf$TF.Number <- factor(cortffractiondf$TF.Number,levels = c(1:(dim(cortffractiondf)[1]/2)))


cortffractiondf$Significance <- "N"
cortffractiondf[cortffractiondf$Target.List %in% "Correlated DEG-DAR Pairs" & cortffractiondf$P.Value <= 0.05 & cortffractiondf$P.Value > 0.01,"Significance"] <- "*"
cortffractiondf[cortffractiondf$Target.List %in% "Correlated DEG-DAR Pairs" & cortffractiondf$P.Value <= 0.01 & cortffractiondf$P.Value > 0.001,"Significance"] <- "**"
cortffractiondf[cortffractiondf$Target.List %in% "Correlated DEG-DAR Pairs" & cortffractiondf$P.Value <= 0.001,"Significance"] <- "***"

png(file = "Figure 4A_021325.png",width = 10,height = 5,units = "in",res = 600)
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

pdf(file = "Figure 4A_021325.pdf",width = 10,height = 5)
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

png(file = "Figure 4B_021325.png",width = 8,height = 5,units = "in",res = 600)
pheatmap(lung50peakcorfin[unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin) & lung50peakmotifs$Motif.Name %in% "Sp2(Zf)/HEK293-Sp2.eGFP-ChIP-Seq(Encode)/Homer","PositionID"]),colnames(lung50peakcorfin[unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin) & lung50peakmotifs$Motif.Name %in% "Sp2(Zf)/HEK293-Sp2.eGFP-ChIP-Seq(Encode)/Homer","PositionID"]),colSums(lung50peakcorfin[unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin) & lung50peakmotifs$Motif.Name %in% "Sp2(Zf)/HEK293-Sp2.eGFP-ChIP-Seq(Encode)/Homer","PositionID"]),]) != 0])],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_col = enstosym[colnames(lung50peakcorfin[unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin) & lung50peakmotifs$Motif.Name %in% "Sp2(Zf)/HEK293-Sp2.eGFP-ChIP-Seq(Encode)/Homer","PositionID"]),colSums(lung50peakcorfin[unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin) & lung50peakmotifs$Motif.Name %in% "Sp2(Zf)/HEK293-Sp2.eGFP-ChIP-Seq(Encode)/Homer","PositionID"]),]) != 0]),"Symbol"],fontsize = 20,cluster_rows = F, cluster_cols = F)
dev.off()

pdf(file = "Figure 4B_021325.pdf",width = 8,height = 5)
pheatmap(lung50peakcorfin[unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin) & lung50peakmotifs$Motif.Name %in% "Sp2(Zf)/HEK293-Sp2.eGFP-ChIP-Seq(Encode)/Homer","PositionID"]),colnames(lung50peakcorfin[unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin) & lung50peakmotifs$Motif.Name %in% "Sp2(Zf)/HEK293-Sp2.eGFP-ChIP-Seq(Encode)/Homer","PositionID"]),colSums(lung50peakcorfin[unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin) & lung50peakmotifs$Motif.Name %in% "Sp2(Zf)/HEK293-Sp2.eGFP-ChIP-Seq(Encode)/Homer","PositionID"]),]) != 0])],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_col = enstosym[colnames(lung50peakcorfin[unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin) & lung50peakmotifs$Motif.Name %in% "Sp2(Zf)/HEK293-Sp2.eGFP-ChIP-Seq(Encode)/Homer","PositionID"]),colSums(lung50peakcorfin[unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin) & lung50peakmotifs$Motif.Name %in% "Sp2(Zf)/HEK293-Sp2.eGFP-ChIP-Seq(Encode)/Homer","PositionID"]),]) != 0]),"Symbol"],fontsize = 20,cluster_rows = F, cluster_cols = F)
dev.off()


png(file = "Figure 4C_021325.png",width = 8,height = 5,units = "in",res = 600)
pheatmap(gastro50peakcorfin[unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin) & gastro50peakmotifs$Motif.Name %in% "BMYB(HTH)/Hela-BMYB-ChIP-Seq(GSE27030)/Homer","PositionID"]),colnames(gastro50peakcorfin[unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin) & gastro50peakmotifs$Motif.Name %in% "BMYB(HTH)/Hela-BMYB-ChIP-Seq(GSE27030)/Homer","PositionID"]),colSums(gastro50peakcorfin[unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin) & gastro50peakmotifs$Motif.Name %in% "BMYB(HTH)/Hela-BMYB-ChIP-Seq(GSE27030)/Homer","PositionID"]),]) != 0])],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_col = enstosym[colnames(gastro50peakcorfin[unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin) & gastro50peakmotifs$Motif.Name %in% "BMYB(HTH)/Hela-BMYB-ChIP-Seq(GSE27030)/Homer","PositionID"]),colSums(gastro50peakcorfin[unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin) & gastro50peakmotifs$Motif.Name %in% "BMYB(HTH)/Hela-BMYB-ChIP-Seq(GSE27030)/Homer","PositionID"]),]) != 0]),"Symbol"],fontsize = 20,cluster_rows = F, cluster_cols = F)
dev.off()

pdf(file = "Figure 4C_021325.pdf",width = 8,height = 5)
pheatmap(gastro50peakcorfin[unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin) & gastro50peakmotifs$Motif.Name %in% "BMYB(HTH)/Hela-BMYB-ChIP-Seq(GSE27030)/Homer","PositionID"]),colnames(gastro50peakcorfin[unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin) & gastro50peakmotifs$Motif.Name %in% "BMYB(HTH)/Hela-BMYB-ChIP-Seq(GSE27030)/Homer","PositionID"]),colSums(gastro50peakcorfin[unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin) & gastro50peakmotifs$Motif.Name %in% "BMYB(HTH)/Hela-BMYB-ChIP-Seq(GSE27030)/Homer","PositionID"]),]) != 0])],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_col = enstosym[colnames(gastro50peakcorfin[unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin) & gastro50peakmotifs$Motif.Name %in% "BMYB(HTH)/Hela-BMYB-ChIP-Seq(GSE27030)/Homer","PositionID"]),colSums(gastro50peakcorfin[unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin) & gastro50peakmotifs$Motif.Name %in% "BMYB(HTH)/Hela-BMYB-ChIP-Seq(GSE27030)/Homer","PositionID"]),]) != 0]),"Symbol"],fontsize = 20,cluster_rows = F, cluster_cols = F)
dev.off()


save.image("Figure4A_to_4L_S18_SuppTableS2_021325.RData")

save(file = "peakgenecormod_021325.RData",list = c("gastropeakcormod","heartpeakcormod","hippopeakcormod","kidneypeakcormod","liverpeakcormod","lungpeakcormod","brownpeakcormod","whitepeakcormod"))