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

load("Figure3C_3F_S8.RData")
allpeakmotifs <- readRDS("allpeakmotifs.RDS")
load("activepeakfiles.RData")
tfanno <- readRDS("tfanno.RDS")

gastro50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% gastroactivepeaks,]
heart50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% heartactivepeaks,]
hippo50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% hippoactivepeaks,]
kidney50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% kidneyactivepeaks,]
liver50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% liveractivepeaks,]
lung50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% lungactivepeaks,]
brown50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% brownactivepeaks,]
white50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% whiteactivepeaks,]

gastropeakcormod <- gastropeakcor*(gastropeakdistance > 0 & gastropeakdistance < 500000)
gastroconnectsigpeaks <- rownames(gastropeakcormod[apply(abs(gastropeakcormod),1,max) > 0.5,])
gastroconnectsiggenes <- colnames(gastropeakcormod[,apply(abs(gastropeakcormod),2,max) > 0.5])
gastropeakcortrim <- gastropeakcormod[gastroconnectsigpeaks,gastroconnectsiggenes]

gastro50connectsigpeakstfanno <- gastroconnectsigpeaks[gastroconnectsigpeaks %in% gastro50peakmotifs$PositionID]
gastro50connectsigpeakstfs <- unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% gastro50connectsigpeakstfanno,"Motif.Name"])
gastro50expressedsigpeakstfens <- intersect(tfanno[gastro50connectsigpeakstfs,"Ensembl"],rownames(gastrol2fcmat))
gastro50expressedsigpeakstfs <- rownames(tfanno[tfanno$Ensembl %in% gastro50expressedsigpeakstfens,])
gastro50connectsigpeaksexprtfanno <- intersect(gastro50connectsigpeakstfanno,gastro50peakmotifs[gastro50peakmotifs$Motif.Name %in% gastro50expressedsigpeakstfs,"PositionID"])


heartpeakcormod <- heartpeakcor*(heartpeakdistance > 0 & heartpeakdistance < 500000)
heartconnectsigpeaks <- rownames(heartpeakcormod[apply(abs(heartpeakcormod),1,max) > 0.5,])
heartconnectsiggenes <- colnames(heartpeakcormod[,apply(abs(heartpeakcormod),2,max) > 0.5])
heartpeakcortrim <- heartpeakcormod[heartconnectsigpeaks,heartconnectsiggenes]

heart50connectsigpeakstfanno <- heartconnectsigpeaks[heartconnectsigpeaks %in% heart50peakmotifs$PositionID]
heart50connectsigpeakstfs <- unique(heart50peakmotifs[heart50peakmotifs$PositionID %in% heart50connectsigpeakstfanno,"Motif.Name"])
heart50expressedsigpeakstfens <- intersect(tfanno[heart50connectsigpeakstfs,"Ensembl"],rownames(heartl2fcmat))
heart50expressedsigpeakstfs <- rownames(tfanno[tfanno$Ensembl %in% heart50expressedsigpeakstfens,])
heart50connectsigpeaksexprtfanno <- intersect(heart50connectsigpeakstfanno,heart50peakmotifs[heart50peakmotifs$Motif.Name %in% heart50expressedsigpeakstfs,"PositionID"])


hippopeakcormod <- hippopeakcor*(hippopeakdistance > 0 & hippopeakdistance < 500000)
hippoconnectsigpeaks <- rownames(hippopeakcormod[apply(abs(hippopeakcormod),1,max) > 0.5,])
hippoconnectsiggenes <- colnames(hippopeakcormod[,apply(abs(hippopeakcormod),2,max) > 0.5])
hippopeakcortrim <- hippopeakcormod[hippoconnectsigpeaks,hippoconnectsiggenes]

hippo50connectsigpeakstfanno <- hippoconnectsigpeaks[hippoconnectsigpeaks %in% hippo50peakmotifs$PositionID]
hippo50connectsigpeakstfs <- unique(hippo50peakmotifs[hippo50peakmotifs$PositionID %in% hippo50connectsigpeakstfanno,"Motif.Name"])
hippo50expressedsigpeakstfens <- intersect(tfanno[hippo50connectsigpeakstfs,"Ensembl"],rownames(hippol2fcmat))
hippo50expressedsigpeakstfs <- rownames(tfanno[tfanno$Ensembl %in% hippo50expressedsigpeakstfens,])
hippo50connectsigpeaksexprtfanno <- intersect(hippo50connectsigpeakstfanno,hippo50peakmotifs[hippo50peakmotifs$Motif.Name %in% hippo50expressedsigpeakstfs,"PositionID"])


kidneypeakcormod <- kidneypeakcor*(kidneypeakdistance > 0 & kidneypeakdistance < 500000)
kidneyconnectsigpeaks <- rownames(kidneypeakcormod[apply(abs(kidneypeakcormod),1,max) > 0.5,])
kidneyconnectsiggenes <- colnames(kidneypeakcormod[,apply(abs(kidneypeakcormod),2,max) > 0.5])
kidneypeakcortrim <- kidneypeakcormod[kidneyconnectsigpeaks,kidneyconnectsiggenes]

kidney50connectsigpeakstfanno <- kidneyconnectsigpeaks[kidneyconnectsigpeaks %in% kidney50peakmotifs$PositionID]
kidney50connectsigpeakstfs <- unique(kidney50peakmotifs[kidney50peakmotifs$PositionID %in% kidney50connectsigpeakstfanno,"Motif.Name"])
kidney50expressedsigpeakstfens <- intersect(tfanno[kidney50connectsigpeakstfs,"Ensembl"],rownames(kidneyl2fcmat))
kidney50expressedsigpeakstfs <- rownames(tfanno[tfanno$Ensembl %in% kidney50expressedsigpeakstfens,])
kidney50connectsigpeaksexprtfanno <- intersect(kidney50connectsigpeakstfanno,kidney50peakmotifs[kidney50peakmotifs$Motif.Name %in% kidney50expressedsigpeakstfs,"PositionID"])


liverpeakcormod <- liverpeakcor*(liverpeakdistance > 0 & liverpeakdistance < 500000)
liverconnectsigpeaks <- rownames(liverpeakcormod[apply(abs(liverpeakcormod),1,max) > 0.5,])
liverconnectsiggenes <- colnames(liverpeakcormod[,apply(abs(liverpeakcormod),2,max) > 0.5])
liverpeakcortrim <- liverpeakcormod[liverconnectsigpeaks,liverconnectsiggenes]

liver50connectsigpeakstfanno <- liverconnectsigpeaks[liverconnectsigpeaks %in% liver50peakmotifs$PositionID]
liver50connectsigpeakstfs <- unique(liver50peakmotifs[liver50peakmotifs$PositionID %in% liver50connectsigpeakstfanno,"Motif.Name"])
liver50expressedsigpeakstfens <- intersect(tfanno[liver50connectsigpeakstfs,"Ensembl"],rownames(liverl2fcmat))
liver50expressedsigpeakstfs <- rownames(tfanno[tfanno$Ensembl %in% liver50expressedsigpeakstfens,])
liver50connectsigpeaksexprtfanno <- intersect(liver50connectsigpeakstfanno,liver50peakmotifs[liver50peakmotifs$Motif.Name %in% liver50expressedsigpeakstfs,"PositionID"])

lungpeakcormod <- lungpeakcor*(lungpeakdistance > 0 & lungpeakdistance < 500000)
lungconnectsigpeaks <- rownames(lungpeakcormod[apply(abs(lungpeakcormod),1,max) > 0.5,])
lungconnectsiggenes <- colnames(lungpeakcormod[,apply(abs(lungpeakcormod),2,max) > 0.5])
lungpeakcortrim <- lungpeakcormod[lungconnectsigpeaks,lungconnectsiggenes]

lung50connectsigpeakstfanno <- lungconnectsigpeaks[lungconnectsigpeaks %in% lung50peakmotifs$PositionID]
lung50connectsigpeakstfs <- unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% lung50connectsigpeakstfanno,"Motif.Name"])
lung50expressedsigpeakstfens <- intersect(tfanno[lung50connectsigpeakstfs,"Ensembl"],rownames(lungl2fcmat))
lung50expressedsigpeakstfs <- rownames(tfanno[tfanno$Ensembl %in% lung50expressedsigpeakstfens,])
lung50connectsigpeaksexprtfanno <- intersect(lung50connectsigpeakstfanno,lung50peakmotifs[lung50peakmotifs$Motif.Name %in% lung50expressedsigpeakstfs,"PositionID"])


brownpeakcormod <- brownpeakcor*(brownpeakdistance > 0 & brownpeakdistance < 500000)
brownconnectsigpeaks <- rownames(brownpeakcormod[apply(abs(brownpeakcormod),1,max) > 0.5,])
brownconnectsiggenes <- colnames(brownpeakcormod[,apply(abs(brownpeakcormod),2,max) > 0.5])
brownpeakcortrim <- brownpeakcormod[brownconnectsigpeaks,brownconnectsiggenes]

brown50connectsigpeakstfanno <- brownconnectsigpeaks[brownconnectsigpeaks %in% brown50peakmotifs$PositionID]
brown50connectsigpeakstfs <- unique(brown50peakmotifs[brown50peakmotifs$PositionID %in% brown50connectsigpeakstfanno,"Motif.Name"])
brown50expressedsigpeakstfens <- intersect(tfanno[brown50connectsigpeakstfs,"Ensembl"],rownames(brownl2fcmat))
brown50expressedsigpeakstfs <- rownames(tfanno[tfanno$Ensembl %in% brown50expressedsigpeakstfens,])
brown50connectsigpeaksexprtfanno <- intersect(brown50connectsigpeakstfanno,brown50peakmotifs[brown50peakmotifs$Motif.Name %in% brown50expressedsigpeakstfs,"PositionID"])


whitepeakcormod <- whitepeakcor*(whitepeakdistance > 0 & whitepeakdistance < 500000)
whiteconnectsigpeaks <- rownames(whitepeakcormod[apply(abs(whitepeakcormod),1,max) > 0.5,])
whiteconnectsiggenes <- colnames(whitepeakcormod[,apply(abs(whitepeakcormod),2,max) > 0.5])
whitepeakcortrim <- whitepeakcormod[whiteconnectsigpeaks,whiteconnectsiggenes]

white50connectsigpeakstfanno <- whiteconnectsigpeaks[whiteconnectsigpeaks %in% white50peakmotifs$PositionID]
white50connectsigpeakstfs <- unique(white50peakmotifs[white50peakmotifs$PositionID %in% white50connectsigpeakstfanno,"Motif.Name"])
white50expressedsigpeakstfens <- intersect(tfanno[white50connectsigpeakstfs,"Ensembl"],rownames(whitel2fcmat))
white50expressedsigpeakstfs <- rownames(tfanno[tfanno$Ensembl %in% white50expressedsigpeakstfens,])
white50connectsigpeaksexprtfanno <- intersect(white50connectsigpeakstfanno,white50peakmotifs[white50peakmotifs$Motif.Name %in% white50expressedsigpeakstfs,"PositionID"])

gastro50peakcortrimmer <- gastropeakcortrim[gastro50connectsigpeaksexprtfanno,]
gastro50peakcorfin <- gastro50peakcortrimmer[,abs(apply(gastro50peakcortrimmer,2,mean)) > 0]
heart50peakcortrimmer <- heartpeakcortrim[heart50connectsigpeaksexprtfanno,]
heart50peakcorfin <- heart50peakcortrimmer[,abs(apply(heart50peakcortrimmer,2,mean)) > 0]
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

pdf(file = "Figure 4A.pdf",width = 8,height = 5)
pheatmap(gastro50peakcorfin[unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin) & gastro50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),colnames(gastro50peakcorfin[unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin) & gastro50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),colSums(gastro50peakcorfin[unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin) & gastro50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),]) != 0])],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_col = enstosym[colnames(gastro50peakcorfin[unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin) & gastro50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),colSums(gastro50peakcorfin[unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin) & gastro50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),]) != 0]),"Symbol"],fontsize = 20,cluster_rows = F, cluster_cols = F)
dev.off()

png(file = "Figure 4A.png",width = 8,height = 5,units = "in",res = 600)
pheatmap(gastro50peakcorfin[unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin) & gastro50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),colnames(gastro50peakcorfin[unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin) & gastro50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),colSums(gastro50peakcorfin[unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin) & gastro50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),]) != 0])],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_col = enstosym[colnames(gastro50peakcorfin[unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin) & gastro50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),colSums(gastro50peakcorfin[unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin) & gastro50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),]) != 0]),"Symbol"],fontsize = 20,cluster_rows = F, cluster_cols = F)
dev.off()

pdf(file = "Figure 4B.pdf",width = 8,height = 5)
pheatmap(lung50peakcorfin[unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin) & lung50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),colnames(lung50peakcorfin[unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin) & lung50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),colSums(lung50peakcorfin[unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin) & lung50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),]) != 0])],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_col = enstosym[colnames(lung50peakcorfin[unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin) & lung50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),colSums(lung50peakcorfin[unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin) & lung50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),]) != 0]),"Symbol"],fontsize = 20,cluster_rows = F, cluster_cols = F)
dev.off()

png(file = "Figure 4B.png",width = 8,height = 5,units = "in",res = 600)
pheatmap(lung50peakcorfin[unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin) & lung50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),colnames(lung50peakcorfin[unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin) & lung50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),colSums(lung50peakcorfin[unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin) & lung50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),]) != 0])],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_col = enstosym[colnames(lung50peakcorfin[unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin) & lung50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),colSums(lung50peakcorfin[unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin) & lung50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),]) != 0]),"Symbol"],fontsize = 20,cluster_rows = F, cluster_cols = F)
dev.off()

pdf(file = "Figure 4C.pdf",width = 9,height = 11)
pheatmap(liver50peakcorfin[unique(liver50peakmotifs[liver50peakmotifs$PositionID %in% rownames(liver50peakcorfin) & liver50peakmotifs$Motif.Name %in% "Smad3(MAD)/NPC-Smad3-ChIP-Seq(GSE36673)/Homer","PositionID"]),colnames(liver50peakcorfin[unique(liver50peakmotifs[liver50peakmotifs$PositionID %in% rownames(liver50peakcorfin) & liver50peakmotifs$Motif.Name %in% "Smad3(MAD)/NPC-Smad3-ChIP-Seq(GSE36673)/Homer","PositionID"]),colSums(liver50peakcorfin[unique(liver50peakmotifs[liver50peakmotifs$PositionID %in% rownames(liver50peakcorfin) & liver50peakmotifs$Motif.Name %in% "Smad3(MAD)/NPC-Smad3-ChIP-Seq(GSE36673)/Homer","PositionID"]),]) != 0])],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_col = enstosym[colnames(liver50peakcorfin[unique(liver50peakmotifs[liver50peakmotifs$PositionID %in% rownames(liver50peakcorfin) & liver50peakmotifs$Motif.Name %in% "Smad3(MAD)/NPC-Smad3-ChIP-Seq(GSE36673)/Homer","PositionID"]),colSums(liver50peakcorfin[unique(liver50peakmotifs[liver50peakmotifs$PositionID %in% rownames(liver50peakcorfin) & liver50peakmotifs$Motif.Name %in% "Smad3(MAD)/NPC-Smad3-ChIP-Seq(GSE36673)/Homer","PositionID"]),]) != 0]),"Symbol"],fontsize = 20,cluster_rows = F, cluster_cols = F)
dev.off()

png(file = "Figure 4C.png",width = 9,height = 11,units = "in",res = 600)
pheatmap(liver50peakcorfin[unique(liver50peakmotifs[liver50peakmotifs$PositionID %in% rownames(liver50peakcorfin) & liver50peakmotifs$Motif.Name %in% "Smad3(MAD)/NPC-Smad3-ChIP-Seq(GSE36673)/Homer","PositionID"]),colnames(liver50peakcorfin[unique(liver50peakmotifs[liver50peakmotifs$PositionID %in% rownames(liver50peakcorfin) & liver50peakmotifs$Motif.Name %in% "Smad3(MAD)/NPC-Smad3-ChIP-Seq(GSE36673)/Homer","PositionID"]),colSums(liver50peakcorfin[unique(liver50peakmotifs[liver50peakmotifs$PositionID %in% rownames(liver50peakcorfin) & liver50peakmotifs$Motif.Name %in% "Smad3(MAD)/NPC-Smad3-ChIP-Seq(GSE36673)/Homer","PositionID"]),]) != 0])],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_col = enstosym[colnames(liver50peakcorfin[unique(liver50peakmotifs[liver50peakmotifs$PositionID %in% rownames(liver50peakcorfin) & liver50peakmotifs$Motif.Name %in% "Smad3(MAD)/NPC-Smad3-ChIP-Seq(GSE36673)/Homer","PositionID"]),colSums(liver50peakcorfin[unique(liver50peakmotifs[liver50peakmotifs$PositionID %in% rownames(liver50peakcorfin) & liver50peakmotifs$Motif.Name %in% "Smad3(MAD)/NPC-Smad3-ChIP-Seq(GSE36673)/Homer","PositionID"]),]) != 0]),"Symbol"],fontsize = 20,cluster_rows = F, cluster_cols = F)
dev.off()




#####
# Supplemental Table S1
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


for(i in 1:dim(heart50peakcorfin)[2]){
  print(i)
  ourgene <- colnames(heart50peakcorfin)[i]
  ourpeaks <- rownames(heart50peakcorfin)[heart50peakcorfin[,ourgene] != 0]
  if(length(ourpeaks) == 1){
    ourpeak <- ourpeaks
    ourdfentry <- data.frame("Tissue" = "HEART",
                             "DEG" = ourgene,
                             "SYMBOL" = enstosym[ourgene,"Symbol"],
                             "DAR" = ourpeak,
                             "Distance" = abs(heartrnasiganno[ourgene,"start"] - ((peakanno[ourpeak,"start"]+peakanno[ourpeak,"end"])/2)),
                             "Correlation" = heart50peakcorfin[ourpeak,ourgene],
                             "TF Motifs" = toString(unique(heart50peakmotifs[heart50peakmotifs$PositionID %in% ourpeak,"Motif.Name"])),
                             "RNA FW1 L2FC" = heartl2fcmat[ourgene,"F W1"],
                             "RNA FW2 L2FC" = heartl2fcmat[ourgene,"F W2"],
                             "RNA FW4 L2FC" = heartl2fcmat[ourgene,"F W4"],
                             "RNA FW8 L2FC" = heartl2fcmat[ourgene,"F W8"],
                             "RNA MW1 L2FC" = heartl2fcmat[ourgene,"M W1"],
                             "RNA MW2 L2FC" = heartl2fcmat[ourgene,"M W2"],
                             "RNA MW4 L2FC" = heartl2fcmat[ourgene,"M W4"],
                             "RNA MW8 L2FC" = heartl2fcmat[ourgene,"M W8"],
                             "ATAC FW1 L2FC" = heartsigatacl2fc[ourpeak,"F W1"],
                             "ATAC FW2 L2FC" = heartsigatacl2fc[ourpeak,"F W2"],
                             "ATAC FW4 L2FC" = heartsigatacl2fc[ourpeak,"F W4"],
                             "ATAC FW8 L2FC" = heartsigatacl2fc[ourpeak,"F W8"],
                             "ATAC MW1 L2FC" = heartsigatacl2fc[ourpeak,"M W1"],
                             "ATAC MW2 L2FC" = heartsigatacl2fc[ourpeak,"M W2"],
                             "ATAC MW4 L2FC" = heartsigatacl2fc[ourpeak,"M W4"],
                             "ATAC MW8 L2FC" = heartsigatacl2fc[ourpeak,"M W8"])
    dar_deg_connection50df <- rbind(dar_deg_connection50df,ourdfentry)
  } else if(length(ourpeaks) > 1){
    for(j in 1:length(ourpeaks)){
      ourpeak <- ourpeaks[j]
      ourdfentry <- data.frame("Tissue" = "HEART",
                               "DEG" = ourgene,
                               "SYMBOL" = enstosym[ourgene,"Symbol"],
                               "DAR" = ourpeak,
                               "Distance" = abs(heartrnasiganno[ourgene,"start"] - ((peakanno[ourpeak,"start"]+peakanno[ourpeak,"end"])/2)),
                               "Correlation" = heart50peakcorfin[ourpeak,ourgene],
                               "TF Motifs" = toString(unique(heart50peakmotifs[heart50peakmotifs$PositionID %in% ourpeak,"Motif.Name"])),
                               "RNA FW1 L2FC" = heartl2fcmat[ourgene,"F W1"],
                               "RNA FW2 L2FC" = heartl2fcmat[ourgene,"F W2"],
                               "RNA FW4 L2FC" = heartl2fcmat[ourgene,"F W4"],
                               "RNA FW8 L2FC" = heartl2fcmat[ourgene,"F W8"],
                               "RNA MW1 L2FC" = heartl2fcmat[ourgene,"M W1"],
                               "RNA MW2 L2FC" = heartl2fcmat[ourgene,"M W2"],
                               "RNA MW4 L2FC" = heartl2fcmat[ourgene,"M W4"],
                               "RNA MW8 L2FC" = heartl2fcmat[ourgene,"M W8"],
                               "ATAC FW1 L2FC" = heartsigatacl2fc[ourpeak,"F W1"],
                               "ATAC FW2 L2FC" = heartsigatacl2fc[ourpeak,"F W2"],
                               "ATAC FW4 L2FC" = heartsigatacl2fc[ourpeak,"F W4"],
                               "ATAC FW8 L2FC" = heartsigatacl2fc[ourpeak,"F W8"],
                               "ATAC MW1 L2FC" = heartsigatacl2fc[ourpeak,"M W1"],
                               "ATAC MW2 L2FC" = heartsigatacl2fc[ourpeak,"M W2"],
                               "ATAC MW4 L2FC" = heartsigatacl2fc[ourpeak,"M W4"],
                               "ATAC MW8 L2FC" = heartsigatacl2fc[ourpeak,"M W8"])
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
write.csv(dar_deg_connection50df,file = "Supplemental Table 1.csv",row.names = F)


#####
# Figure 4D-G
####

ourgene <- enstosym[enstosym$Symbol %in% "Igf2","Ensembl"]
igf2_gastro_corpeaks <- rownames(gastropeakcortrim)[abs(gastropeakcortrim[,ourgene]) > 0.5]
ourpeak <- igf2_gastro_corpeaks
ourdf <- data.frame("RNA.L2FC" = gastrol2fcmat[ourgene,],
                    "ATAC.L2FC" = gastrosigatacl2fc[ourpeak,],
                    "Label" = colnames(gastrol2fcmat))
pdf(file = "Figure 4D.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("SKM-GN: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

png(file = "Figure 4D.png",width = 6,height = 6,units = "in",res = 600)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("SKM-GN: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Sall2","Ensembl"]
sall2_gastro_corpeaks <- rownames(gastropeakcortrim)[abs(gastropeakcortrim[,ourgene]) > 0.5]
ourpeak <- sall2_gastro_corpeaks
ourdf <- data.frame("RNA.L2FC" = gastrol2fcmat[ourgene,],
                    "ATAC.L2FC" = gastrosigatacl2fc[ourpeak,],
                    "Label" = colnames(gastrol2fcmat))

pdf(file = "Figure 4E.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("SKM-GN: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

png(file = "Figure 4E.png",width = 6,height = 6,units = "in",res = 600)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("SKM-GN: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Nfkb2","Ensembl"]
nfkb2_lung_corpeaks <- rownames(lungpeakcortrim)[abs(lungpeakcortrim[,ourgene]) > 0.5]
ourpeak <- nfkb2_lung_corpeaks
ourdf <- data.frame("RNA.L2FC" = lungl2fcmat[ourgene,],
                    "ATAC.L2FC" = lungsigatacl2fc[ourpeak,],
                    "Label" = colnames(lungl2fcmat))
pdf(file = "Figure 4F.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LUNG: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

png(file = "Figure 4F.png",width = 6,height = 6,units = "in",res = 600)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LUNG: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Fkbp4","Ensembl"]
fkbp4_liver_corpeaks <- rownames(liverpeakcortrim)[abs(liverpeakcortrim[,ourgene]) > 0.5]
ourpeak <- fkbp4_liver_corpeaks[2]
ourdf <- data.frame("RNA.L2FC" = liverl2fcmat[ourgene,],
                    "ATAC.L2FC" = liversigatacl2fc[ourpeak,],
                    "Label" = colnames(liverl2fcmat))
pdf(file = "Figure 4G.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

png(file = "Figure 4G.png",width = 6,height = 6,units = "in",res = 600)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

#####
# Supplemental Figure S12
####

ourgene <- enstosym[enstosym$Symbol %in% "Glul","Ensembl"]
glul_liver_corpeaks <- rownames(liverpeakcortrim)[abs(liverpeakcortrim[,ourgene]) > 0.5]
ourpeak <- glul_liver_corpeaks[2]
ourdf <- data.frame("RNA.L2FC" = liverl2fcmat[ourgene,],
                    "ATAC.L2FC" = liversigatacl2fc[ourpeak,],
                    "Label" = colnames(liverl2fcmat))
pdf(file = "Supplemental Figure S12A.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()
png(file = "Supplemental Figure S12A.png",width = 6,height = 6,units = "in",res = 600)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Glul","Ensembl"]
glul_liver_corpeaks <- rownames(liverpeakcortrim)[abs(liverpeakcortrim[,ourgene]) > 0.5]
ourpeak <- glul_liver_corpeaks[1]
ourdf <- data.frame("RNA.L2FC" = liverl2fcmat[ourgene,],
                    "ATAC.L2FC" = liversigatacl2fc[ourpeak,],
                    "Label" = colnames(liverl2fcmat))
pdf(file = "Supplemental Figure S12B.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()
png(file = "Supplemental Figure S12B.png",width = 6,height = 6,units = "in",res = 600)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Lpar3","Ensembl"]
lpar3_liver_corpeaks <- rownames(liverpeakcortrim)[abs(liverpeakcortrim[,ourgene]) > 0.5]
ourpeak <- lpar3_liver_corpeaks[5]
ourdf <- data.frame("RNA.L2FC" = liverl2fcmat[ourgene,],
                    "ATAC.L2FC" = liversigatacl2fc[ourpeak,],
                    "Label" = colnames(liverl2fcmat))
pdf(file = "Supplemental Figure S12C.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()
png(file = "Supplemental Figure S12C.png",width = 6,height = 6,units = "in",res = 600)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Lpar3","Ensembl"]
lpar3_liver_corpeaks <- rownames(liverpeakcortrim)[abs(liverpeakcortrim[,ourgene]) > 0.5]
ourpeak <- lpar3_liver_corpeaks[6]
ourdf <- data.frame("RNA.L2FC" = liverl2fcmat[ourgene,],
                    "ATAC.L2FC" = liversigatacl2fc[ourpeak,],
                    "Label" = colnames(liverl2fcmat))
pdf(file = "Supplemental Figure S12D.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()
png(file = "Supplemental Figure S12D.png",width = 6,height = 6,units = "in",res = 600)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Serpina4","Ensembl"]
serpina4_liver_corpeaks <- rownames(liverpeakcortrim)[abs(liverpeakcortrim[,ourgene]) > 0.5]
ourpeak <- serpina4_liver_corpeaks[3]
ourdf <- data.frame("RNA.L2FC" = liverl2fcmat[ourgene,],
                    "ATAC.L2FC" = liversigatacl2fc[ourpeak,],
                    "Label" = colnames(liverl2fcmat))
pdf(file = "Supplemental Figure S12E.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()
png(file = "Supplemental Figure S12E.png",width = 6,height = 6,units = "in",res = 600)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Serpina4","Ensembl"]
serpina4_liver_corpeaks <- rownames(liverpeakcortrim)[abs(liverpeakcortrim[,ourgene]) > 0.5]
ourpeak <- serpina4_liver_corpeaks[4]
ourdf <- data.frame("RNA.L2FC" = liverl2fcmat[ourgene,],
                    "ATAC.L2FC" = liversigatacl2fc[ourpeak,],
                    "Label" = colnames(liverl2fcmat))
pdf(file = "Supplemental Figure S12F.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()
png(file = "Supplemental Figure S12F.png",width = 6,height = 6,units = "in",res = 600)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Abhd2","Ensembl"]
abhd2_liver_corpeaks <- rownames(liverpeakcortrim)[abs(liverpeakcortrim[,ourgene]) > 0.5]
ourpeak <- abhd2_liver_corpeaks[5]
ourdf <- data.frame("RNA.L2FC" = liverl2fcmat[ourgene,],
                    "ATAC.L2FC" = liversigatacl2fc[ourpeak,],
                    "Label" = colnames(liverl2fcmat))
pdf(file = "Supplemental Figure S12G.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()
png(file = "Supplemental Figure S12G.png",width = 6,height = 6,units = "in",res = 600)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Onecut1","Ensembl"]
onecut1_liver_corpeaks <- rownames(liverpeakcortrim)[abs(liverpeakcortrim[,ourgene]) > 0.5]
ourpeak <- onecut1_liver_corpeaks
ourdf <- data.frame("RNA.L2FC" = liverl2fcmat[ourgene,],
                    "ATAC.L2FC" = liversigatacl2fc[ourpeak,],
                    "Label" = colnames(liverl2fcmat))
pdf(file = "Supplemental Figure S12H.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()
png(file = "Supplemental Figure S12H.png",width = 6,height = 6,units = "in",res = 600)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Ccnd1","Ensembl"]
ccnd1_liver_corpeaks <- rownames(liverpeakcortrim)[abs(liverpeakcortrim[,ourgene]) > 0.5]
ourpeak <- ccnd1_liver_corpeaks[2]
ourdf <- data.frame("RNA.L2FC" = liverl2fcmat[ourgene,],
                    "ATAC.L2FC" = liversigatacl2fc[ourpeak,],
                    "Label" = colnames(liverl2fcmat))
pdf(file = "Supplemental Figure S12I.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()
png(file = "Supplemental Figure S12I.png",width = 6,height = 6,units = "in",res = 600)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Xbp1","Ensembl"]
xbp1_liver_corpeaks <- rownames(liverpeakcortrim)[abs(liverpeakcortrim[,ourgene]) > 0.5]
ourpeak <- xbp1_liver_corpeaks[1]
ourdf <- data.frame("RNA.L2FC" = liverl2fcmat[ourgene,],
                    "ATAC.L2FC" = liversigatacl2fc[ourpeak,],
                    "Label" = colnames(liverl2fcmat))
pdf(file = "Supplemental Figure S12J.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()
png(file = "Supplemental Figure S12J.png",width = 6,height = 6,units = "in",res = 600)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Egfr","Ensembl"]
egfr_liver_corpeaks <- rownames(liverpeakcortrim)[abs(liverpeakcortrim[,ourgene]) > 0.5]
ourpeak <- egfr_liver_corpeaks[3]
ourdf <- data.frame("RNA.L2FC" = liverl2fcmat[ourgene,],
                    "ATAC.L2FC" = liversigatacl2fc[ourpeak,],
                    "Label" = colnames(liverl2fcmat))
pdf(file = "Supplemental Figure S12K.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()
png(file = "Supplemental Figure S12K.png",width = 6,height = 6,units = "in",res = 600)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Fkbp5","Ensembl"]
fkbp5_gastro_corpeaks <- rownames(gastropeakcortrim)[abs(gastropeakcortrim[,ourgene]) > 0.5]
ourpeak <- fkbp5_gastro_corpeaks
ourdf <- data.frame("RNA.L2FC" = gastrol2fcmat[ourgene,],
                    "ATAC.L2FC" = gastrosigatacl2fc[ourpeak,],
                    "Label" = colnames(gastrol2fcmat))
pdf(file = "Supplemental Figure S12L.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("SKM-GN: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()
png(file = "Supplemental Figure S12L.png",width = 6,height = 6,units = "in",res = 600)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("SKM-GN: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Ppard","Ensembl"]
ppard_gastro_corpeaks <- rownames(gastropeakcortrim)[abs(gastropeakcortrim[,ourgene]) > 0.5]
ourpeak <- ppard_gastro_corpeaks
ourdf <- data.frame("RNA.L2FC" = gastrol2fcmat[ourgene,],
                    "ATAC.L2FC" = gastrosigatacl2fc[ourpeak,],
                    "Label" = colnames(gastrol2fcmat))
pdf(file = "Supplemental Figure S12M.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("SKM-GN: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()
png(file = "Supplemental Figure S12M.png",width = 6,height = 6,units = "in",res = 600)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("SKM-GN: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Sik1","Ensembl"]
sik1_gastro_corpeaks <- rownames(gastropeakcortrim)[abs(gastropeakcortrim[,ourgene]) > 0.5]
ourpeak <- sik1_gastro_corpeaks
ourdf <- data.frame("RNA.L2FC" = gastrol2fcmat[ourgene,],
                    "ATAC.L2FC" = gastrosigatacl2fc[ourpeak,],
                    "Label" = colnames(gastrol2fcmat))
pdf(file = "Supplemental Figure S12N.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("SKM-GN: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()
png(file = "Supplemental Figure S12N.png",width = 6,height = 6,units = "in",res = 600)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("SKM-GN: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Igf1","Ensembl"]
igf1_liver_corpeaks <- rownames(liverpeakcortrim)[abs(liverpeakcortrim[,ourgene]) > 0.5]
ourpeak <- igf1_liver_corpeaks[2]
ourdf <- data.frame("RNA.L2FC" = liverl2fcmat[ourgene,],
                    "ATAC.L2FC" = liversigatacl2fc[ourpeak,],
                    "Label" = colnames(liverl2fcmat))
pdf(file = "Supplemental Figure S12O.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()
png(file = "Supplemental Figure S12O.png",width = 6,height = 6,units = "in",res = 600)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

# Supplemental Figure S12P - 500x600
pdf(file = "Supplemental Figure S12P.pdf",width = 5,height = 6)
pheatmap(liverl2fcmat[intersect(tfanno[liver50peakmotifs[liver50peakmotifs$PositionID %in% igf1_liver_corpeaks,"Motif.Name"],"Ensembl"],rownames(liverl2fcmat)),],cluster_rows = F,cluster_cols = F,labels_row = enstosym[intersect(tfanno[liver50peakmotifs[liver50peakmotifs$PositionID %in% igf1_liver_corpeaks,"Motif.Name"],"Ensembl"],rownames(liverl2fcmat)),"Symbol"],breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = 0,display_numbers = T)
dev.off()
png(file = "Supplemental Figure S12P.png",width = 5,height = 6,units = "in",res = 600)
pheatmap(liverl2fcmat[intersect(tfanno[liver50peakmotifs[liver50peakmotifs$PositionID %in% igf1_liver_corpeaks,"Motif.Name"],"Ensembl"],rownames(liverl2fcmat)),],cluster_rows = F,cluster_cols = F,labels_row = enstosym[intersect(tfanno[liver50peakmotifs[liver50peakmotifs$PositionID %in% igf1_liver_corpeaks,"Motif.Name"],"Ensembl"],rownames(liverl2fcmat)),"Symbol"],breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = 0,display_numbers = T)
dev.off()


save.image("Figure4A_to_4G_S12_SuppTableS1.RData")

