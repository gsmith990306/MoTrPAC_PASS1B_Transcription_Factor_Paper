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

#setwd("~/Mount Sinai/Sealfon Laboratory/MoTrPAC/PASS1B Raw Data")
setwd("C:/Users/gsmit/Documents/Mount Sinai/Sealfon Laboratory/MoTrPAC/PASS1B Raw Data")

load("rnasigl2fcmat.RData")

#####
# Supplemental Figure S2
####

load("rnal2fcmat.RData")

pdf(file = "Supplemental Figure S2Aii_112624.pdf",width = 6,height = 5)
pheatmap(gastrosigrnal2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S2Ai_112624.pdf",width = 6,height = 5)
hist(gastrosigrnal2fc[abs(gastrosigrnal2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "SKM-GN DEG L2FC")
dev.off()

pdf(file = "Supplemental Figure S2Aiii_112624.pdf",width = 6,height = 5)
hist(gastrol2fcmat[abs(gastrol2fcmat) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "SKM-GN All Gene L2FC")
dev.off()

pdf(file = "Supplemental Figure S2Bii_112624.pdf",width = 6,height = 5)
pheatmap(heartsigrnal2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S2Bi_112624.pdf",width = 6,height = 5)
hist(heartsigrnal2fc[abs(heartsigrnal2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "HEART DEG L2FC")
dev.off()

pdf(file = "Supplemental Figure S2Biii_112624.pdf",width = 6,height = 5)
hist(heartl2fcmat[abs(heartl2fcmat) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "HEART All Gene L2FC")
dev.off()

pdf(file = "Supplemental Figure S2Cii_112624.pdf",width = 6,height = 5)
pheatmap(hipposigrnal2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S2Ci_112624.pdf",width = 6,height = 5)
hist(hipposigrnal2fc[abs(hipposigrnal2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "HIPPOC DEG L2FC")
dev.off()

pdf(file = "Supplemental Figure S2Ciii_112624.pdf",width = 6,height = 5)
hist(hippol2fcmat[abs(hippol2fcmat) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "HIPPOC All Gene L2FC")
dev.off()

pdf(file = "Supplemental Figure S2Dii_112624.pdf",width = 6,height = 5)
pheatmap(kidneysigrnal2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S2Di_112624.pdf",width = 6,height = 5)
hist(kidneysigrnal2fc[abs(kidneysigrnal2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "KIDNEY DEG L2FC")
dev.off()

pdf(file = "Supplemental Figure S2Diii_112624.pdf",width = 6,height = 5)
hist(kidneyl2fcmat[abs(kidneyl2fcmat) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "KIDNEY All Gene L2FC")
dev.off()

pdf(file = "Supplemental Figure S2Eii_112624.pdf",width = 6,height = 5)
pheatmap(liversigrnal2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S2Ei_112624.pdf",width = 6,height = 5)
hist(liversigrnal2fc[abs(liversigrnal2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "LIVER DEG L2FC")
dev.off()

pdf(file = "Supplemental Figure S2Eiii_112624.pdf",width = 6,height = 5)
hist(liverl2fcmat[abs(liverl2fcmat) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "LIVER All Gene L2FC")
dev.off()

pdf(file = "Supplemental Figure S2Fii_112624.pdf",width = 6,height = 5)
pheatmap(lungsigrnal2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S2Fi_112624.pdf",width = 6,height = 5)
hist(lungsigrnal2fc[abs(lungsigrnal2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "LUNG DEG L2FC")
dev.off()

pdf(file = "Supplemental Figure S2Fiii_112624.pdf",width = 6,height = 5)
hist(lungl2fcmat[abs(lungl2fcmat) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "LUNG All Gene L2FC")
dev.off()

pdf(file = "Supplemental Figure S2Gii_112624.pdf",width = 6,height = 5)
pheatmap(brownsigrnal2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S2Gi_112624.pdf",width = 6,height = 5)
hist(brownsigrnal2fc[abs(brownsigrnal2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "BAT DEG L2FC")
dev.off()

pdf(file = "Supplemental Figure S2Giii_112624.pdf",width = 6,height = 5)
hist(brownl2fcmat[abs(brownl2fcmat) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "BAT All Gene L2FC")
dev.off()

pdf(file = "Supplemental Figure S2Hii_112624.pdf",width = 6,height = 5)
pheatmap(whitesigrnal2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S2Hi_112624.pdf",width = 6,height = 5)
hist(whitesigrnal2fc[abs(whitesigrnal2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "WAT-SC DEG L2FC")
dev.off()

pdf(file = "Supplemental Figure S2Hiii_112624.pdf",width = 6,height = 5)
hist(whitel2fcmat[abs(whitel2fcmat) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "WAT-SC All Gene L2FC")
dev.off()



png(file = "Supplemental Figure S2Aii_112624.png",width = 6,height = 5,units = "in",res = 600)
pheatmap(gastrosigrnal2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

png(file = "Supplemental Figure S2Ai_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(gastrosigrnal2fc[abs(gastrosigrnal2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "SKM-GN DEG L2FC")
dev.off()

png(file = "Supplemental Figure S2Aiii_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(gastrol2fcmat[abs(gastrol2fcmat) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "SKM-GN All Gene L2FC")
dev.off()

png(file = "Supplemental Figure S2Bii_112624.png",width = 6,height = 5,units = "in",res = 600)
pheatmap(heartsigrnal2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

png(file = "Supplemental Figure S2Bi_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(heartsigrnal2fc[abs(heartsigrnal2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "HEART DEG L2FC")
dev.off()

png(file = "Supplemental Figure S2Biii_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(heartl2fcmat[abs(heartl2fcmat) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "HEART All Gene L2FC")
dev.off()

png(file = "Supplemental Figure S2Cii_112624.png",width = 6,height = 5,units = "in",res = 600)
pheatmap(hipposigrnal2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

png(file = "Supplemental Figure S2Ci_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(hipposigrnal2fc[abs(hipposigrnal2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "HIPPOC DEG L2FC")
dev.off()

png(file = "Supplemental Figure S2Ciii_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(hippol2fcmat[abs(hippol2fcmat) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "HIPPOC All Gene L2FC")
dev.off()

png(file = "Supplemental Figure S2Dii_112624.png",width = 6,height = 5,units = "in",res = 600)
pheatmap(kidneysigrnal2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

png(file = "Supplemental Figure S2Di_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(kidneysigrnal2fc[abs(kidneysigrnal2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "KIDNEY DEG L2FC")
dev.off()

png(file = "Supplemental Figure S2Diii_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(kidneyl2fcmat[abs(kidneyl2fcmat) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "KIDNEY All Gene L2FC")
dev.off()

png(file = "Supplemental Figure S2Eii_112624.png",width = 6,height = 5,units = "in",res = 600)
pheatmap(liversigrnal2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

png(file = "Supplemental Figure S2Ei_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(liversigrnal2fc[abs(liversigrnal2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "LIVER DEG L2FC")
dev.off()

png(file = "Supplemental Figure S2Eiii_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(liverl2fcmat[abs(liverl2fcmat) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "LIVER All Gene L2FC")
dev.off()

png(file = "Supplemental Figure S2Fii_112624.png",width = 6,height = 5,units = "in",res = 600)
pheatmap(lungsigrnal2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

png(file = "Supplemental Figure S2Fi_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(lungsigrnal2fc[abs(lungsigrnal2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "LUNG DEG L2FC")
dev.off()

png(file = "Supplemental Figure S2Fiii_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(lungl2fcmat[abs(lungl2fcmat) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "LUNG All Gene L2FC")
dev.off()

png(file = "Supplemental Figure S2Gii_112624.png",width = 6,height = 5,units = "in",res = 600)
pheatmap(brownsigrnal2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

png(file = "Supplemental Figure S2Gi_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(brownsigrnal2fc[abs(brownsigrnal2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "BAT DEG L2FC")
dev.off()

png(file = "Supplemental Figure S2Giii_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(brownl2fcmat[abs(brownl2fcmat) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "BAT All Gene L2FC")
dev.off()

png(file = "Supplemental Figure S2Hii_112624.png",width = 6,height = 5,units = "in",res = 600)
pheatmap(whitesigrnal2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

png(file = "Supplemental Figure S2Hi_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(whitesigrnal2fc[abs(whitesigrnal2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "WAT-SC DEG L2FC")
dev.off()

png(file = "Supplemental Figure S2Hiii_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(whitel2fcmat[abs(whitel2fcmat) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "WAT-SC All Gene L2FC")
dev.off()

#####
# Supplemental Figure S3
####

load("atacsigl2fcmat.RData")
load("atacl2fcmat.RData")


pdf(file = "Supplemental Figure S3Aii_112624.pdf",width = 6,height = 5)
pheatmap(gastrosigatacl2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S3Ai_112624.pdf",width = 6,height = 5)
hist(gastrosigatacl2fc[abs(gastrosigatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "SKM-GN DAR L2FC")
dev.off()

pdf(file = "Supplemental Figure S3Aiii_112624.pdf",width = 6,height = 5)
hist(gastroatacl2fc[abs(gastroatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "SKM-GN All Peaks L2FC")
dev.off()


pdf(file = "Supplemental Figure S3Bii_112624.pdf",width = 6,height = 5)
pheatmap(heartsigatacl2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S3Bi_112624.pdf",width = 6,height = 5)
hist(heartsigatacl2fc[abs(heartsigatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "HEART DAR L2FC")
dev.off()

pdf(file = "Supplemental Figure S3Biii_112624.pdf",width = 6,height = 5)
hist(heartatacl2fc[abs(heartatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "HEART All Peaks L2FC")
dev.off()

pdf(file = "Supplemental Figure S3Cii_112624.pdf",width = 6,height = 5)
pheatmap(hipposigatacl2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S3Ci_112624.pdf",width = 6,height = 5)
hist(hipposigatacl2fc[abs(hipposigatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "HIPPOC DAR L2FC")
dev.off()

pdf(file = "Supplemental Figure S3Ciii_112624.pdf",width = 6,height = 5)
hist(hippoatacl2fc[abs(hippoatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "HIPPOC All Peaks L2FC")
dev.off()

pdf(file = "Supplemental Figure S3Dii_112624.pdf",width = 6,height = 5)
pheatmap(kidneysigatacl2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S3Di_112624.pdf",width = 6,height = 5)
hist(kidneysigatacl2fc[abs(kidneysigatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "KIDNEY DAR L2FC")
dev.off()

pdf(file = "Supplemental Figure S3Diii_112624.pdf",width = 6,height = 5)
hist(kidneyatacl2fc[abs(kidneyatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "KIDNEY All Peaks L2FC")
dev.off()

pdf(file = "Supplemental Figure S3Eii_112624.pdf",width = 6,height = 5)
pheatmap(liversigatacl2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S3Ei_112624.pdf",width = 6,height = 5)
hist(liversigatacl2fc[abs(liversigatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "LIVER DAR L2FC")
dev.off()

pdf(file = "Supplemental Figure S3Eiii_112624.pdf",width = 6,height = 5)
hist(liveratacl2fc[abs(liveratacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "LIVER All Peaks L2FC")
dev.off()

pdf(file = "Supplemental Figure S3Fii_112624.pdf",width = 6,height = 5)
pheatmap(lungsigatacl2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S3Fi_112624.pdf",width = 6,height = 5)
hist(lungsigatacl2fc[abs(lungsigatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "LUNG DAR L2FC")
dev.off()

pdf(file = "Supplemental Figure S3Fiii_112624.pdf",width = 6,height = 5)
hist(lungatacl2fc[abs(lungatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "LUNG All Peaks L2FC")
dev.off()

pdf(file = "Supplemental Figure S3Gii_112624.pdf",width = 6,height = 5)
pheatmap(brownsigatacl2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S3Gi_112624.pdf",width = 6,height = 5)
hist(brownsigatacl2fc[abs(brownsigatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "BAT DAR L2FC")
dev.off()

pdf(file = "Supplemental Figure S3Giii_112624.pdf",width = 6,height = 5)
hist(brownatacl2fc[abs(brownatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "BAT All Peaks L2FC")
dev.off()

pdf(file = "Supplemental Figure S3Hii_112624.pdf",width = 6,height = 5)
pheatmap(whitesigatacl2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S3Hi_112624.pdf",width = 6,height = 5)
hist(whitesigatacl2fc[abs(whitesigatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "WAT-SC DAR L2FC")
dev.off()

pdf(file = "Supplemental Figure S3Hiii_112624.pdf",width = 6,height = 5)
hist(whiteatacl2fc[abs(whiteatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "WAT-SC All Peaks L2FC")
dev.off()

png(file = "Supplemental Figure S3Aii_112624.png",width = 6,height = 5,units = "in",res = 600)
pheatmap(gastrosigatacl2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

png(file = "Supplemental Figure S3Ai_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(gastrosigatacl2fc[abs(gastrosigatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "SKM-GN DAR L2FC")
dev.off()

png(file = "Supplemental Figure S3Aiii_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(gastroatacl2fc[abs(gastroatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "SKM-GN All Peaks L2FC")
dev.off()


png(file = "Supplemental Figure S3Bii_112624.png",width = 6,height = 5,units = "in",res = 600)
pheatmap(heartsigatacl2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

png(file = "Supplemental Figure S3Bi_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(heartsigatacl2fc[abs(heartsigatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "HEART DAR L2FC")
dev.off()

png(file = "Supplemental Figure S3Biii_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(heartatacl2fc[abs(heartatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "HEART All Peaks L2FC")
dev.off()

png(file = "Supplemental Figure S3Cii_112624.png",width = 6,height = 5,units = "in",res = 600)
pheatmap(hipposigatacl2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

png(file = "Supplemental Figure S3Ci_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(hipposigatacl2fc[abs(hipposigatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "HIPPOC DAR L2FC")
dev.off()

png(file = "Supplemental Figure S3Ciii_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(hippoatacl2fc[abs(hippoatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "HIPPOC All Peaks L2FC")
dev.off()

png(file = "Supplemental Figure S3Dii_112624.png",width = 6,height = 5,units = "in",res = 600)
pheatmap(kidneysigatacl2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

png(file = "Supplemental Figure S3Di_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(kidneysigatacl2fc[abs(kidneysigatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "KIDNEY DAR L2FC")
dev.off()

png(file = "Supplemental Figure S3Diii_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(kidneyatacl2fc[abs(kidneyatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "KIDNEY All Peaks L2FC")
dev.off()

png(file = "Supplemental Figure S3Eii_112624.png",width = 6,height = 5,units = "in",res = 600)
pheatmap(liversigatacl2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

png(file = "Supplemental Figure S3Ei_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(liversigatacl2fc[abs(liversigatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "LIVER DAR L2FC")
dev.off()

png(file = "Supplemental Figure S3Eiii_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(liveratacl2fc[abs(liveratacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "LIVER All Peaks L2FC")
dev.off()

png(file = "Supplemental Figure S3Fii_112624.png",width = 6,height = 5,units = "in",res = 600)
pheatmap(lungsigatacl2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

png(file = "Supplemental Figure S3Fi_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(lungsigatacl2fc[abs(lungsigatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "LUNG DAR L2FC")
dev.off()

png(file = "Supplemental Figure S3Fiii_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(lungatacl2fc[abs(lungatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "LUNG All Peaks L2FC")
dev.off()

png(file = "Supplemental Figure S3Gii_112624.png",width = 6,height = 5,units = "in",res = 600)
pheatmap(brownsigatacl2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

png(file = "Supplemental Figure S3Gi_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(brownsigatacl2fc[abs(brownsigatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "BAT DAR L2FC")
dev.off()

png(file = "Supplemental Figure S3Giii_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(brownatacl2fc[abs(brownatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "BAT All Peaks L2FC")
dev.off()

png(file = "Supplemental Figure S3Hii_112624.png",width = 6,height = 5,units = "in",res = 600)
pheatmap(whitesigatacl2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

png(file = "Supplemental Figure S3Hi_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(whitesigatacl2fc[abs(whitesigatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "WAT-SC DAR L2FC")
dev.off()

png(file = "Supplemental Figure S3Hiii_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(whiteatacl2fc[abs(whiteatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "WAT-SC All Peaks L2FC")
dev.off()


#####
# Supplemental Figure S4
####

load("new_methsigl2fcmat_51624.RData")
load("methl2fcmat.RData")

pdf(file = "Supplemental Figure S4Ai_112624.pdf",width = 6,height = 5)
hist(gastromethsigl2fc[abs(gastromethsigl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "SKM-GN DMR L2FC")
dev.off()

pdf(file = "Supplemental Figure S4Aii_112624.pdf",width = 6,height = 5)
pheatmap(gastromethsigl2fc,breaks = seq(-4,4,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S4Aiii_112624.pdf",width = 6,height = 5)
hist(gastromethl2fc[abs(gastromethl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "SKM-GN All Sites L2FC")
dev.off()

pdf(file = "Supplemental Figure S4Bi_112624.pdf",width = 6,height = 5)
hist(heartmethsigl2fc[abs(heartmethsigl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "HEART DMR L2FC")
dev.off()

pdf(file = "Supplemental Figure S4Bii_112624.pdf",width = 6,height = 5)
pheatmap(heartmethsigl2fc,breaks = seq(-4,4,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S4Biii_112624.pdf",width = 6,height = 5)
hist(heartmethl2fc[abs(heartmethl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "HEART All Sites L2FC")
dev.off()


pdf(file = "Supplemental Figure S4Ci_112624.pdf",width = 6,height = 5)
hist(hippomethsigl2fc[abs(hippomethsigl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "HIPPOC DMR L2FC")
dev.off()

pdf(file = "Supplemental Figure S4Cii_112624.pdf",width = 6,height = 5)
pheatmap(hippomethsigl2fc,breaks = seq(-4,4,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S4Ciii_112624.pdf",width = 6,height = 5)
hist(hippomethl2fc[abs(hippomethl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "HIPPOC All Sites L2FC")
dev.off()

pdf(file = "Supplemental Figure S4Di_112624.pdf",width = 6,height = 5)
hist(kidneymethsigl2fc[abs(kidneymethsigl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "KIDNEY DMR L2FC")
dev.off()

pdf(file = "Supplemental Figure S4Dii_112624.pdf",width = 6,height = 5)
pheatmap(kidneymethsigl2fc,breaks = seq(-4,4,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S4Diii_112624.pdf",width = 6,height = 5)
hist(kidneymethl2fc[abs(kidneymethl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "KIDNEY All Sites L2FC")
dev.off()

pdf(file = "Supplemental Figure S4Ei_112624.pdf",width = 6,height = 5)
hist(livermethsigl2fc[abs(livermethsigl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "LIVER DMR L2FC")
dev.off()

pdf(file = "Supplemental Figure S4Eii_112624.pdf",width = 6,height = 5)
pheatmap(livermethsigl2fc,breaks = seq(-4,4,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S4Eiii_112624.pdf",width = 6,height = 5)
hist(livermethl2fc[abs(livermethl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "LIVER All Sites L2FC")
dev.off()

pdf(file = "Supplemental Figure S4Fi_112624.pdf",width = 6,height = 5)
hist(lungmethsigl2fc[abs(lungmethsigl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "LUNG DMR L2FC")
dev.off()

pdf(file = "Supplemental Figure S4Fii_112624.pdf",width = 6,height = 5)
pheatmap(lungmethsigl2fc,breaks = seq(-4,4,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S4Fiii_112624.pdf",width = 6,height = 5)
hist(lungmethl2fc[abs(lungmethl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "LUNG All Sites L2FC")
dev.off()

pdf(file = "Supplemental Figure S4Gi_112624.pdf",width = 6,height = 5)
hist(brownmethsigl2fc[abs(brownmethsigl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "BAT DMR L2FC")
dev.off()

pdf(file = "Supplemental Figure S4Gii_112624.pdf",width = 6,height = 5)
pheatmap(brownmethsigl2fc,breaks = seq(-4,4,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S4Giii_112624.pdf",width = 6,height = 5)
hist(brownmethl2fc[abs(brownmethl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "BAT All Sites L2FC")
dev.off()

pdf(file = "Supplemental Figure S4Hi_112624.pdf",width = 6,height = 5)
hist(whitemethsigl2fc[abs(whitemethsigl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "WAT-SC DMR L2FC")
dev.off()

pdf(file = "Supplemental Figure S4Hii_112624.pdf",width = 6,height = 5)
pheatmap(whitemethsigl2fc,breaks = seq(-4,4,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S4Hiii_112624.pdf",width = 6,height = 5)
hist(whitemethl2fc[abs(whitemethl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "WAT-SC All Sites L2FC")
dev.off()

png(file = "Supplemental Figure S4Ai_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(gastromethsigl2fc[abs(gastromethsigl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "SKM-GN DMR L2FC")
dev.off()

png(file = "Supplemental Figure S4Aii_112624.png",width = 6,height = 5,units = "in",res = 600)
pheatmap(gastromethsigl2fc,breaks = seq(-4,4,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

png(file = "Supplemental Figure S4Aiii_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(gastromethl2fc[abs(gastromethl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "SKM-GN All Sites L2FC")
dev.off()

png(file = "Supplemental Figure S4Bi_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(heartmethsigl2fc[abs(heartmethsigl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "HEART DMR L2FC")
dev.off()

png(file = "Supplemental Figure S4Bii_112624.png",width = 6,height = 5,units = "in",res = 600)
pheatmap(heartmethsigl2fc,breaks = seq(-4,4,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

png(file = "Supplemental Figure S4Biii_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(heartmethl2fc[abs(heartmethl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "HEART All Sites L2FC")
dev.off()


png(file = "Supplemental Figure S4Ci_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(hippomethsigl2fc[abs(hippomethsigl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "HIPPOC DMR L2FC")
dev.off()

png(file = "Supplemental Figure S4Cii_112624.png",width = 6,height = 5,units = "in",res = 600)
pheatmap(hippomethsigl2fc,breaks = seq(-4,4,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

png(file = "Supplemental Figure S4Ciii_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(hippomethl2fc[abs(hippomethl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "HIPPOC All Sites L2FC")
dev.off()

png(file = "Supplemental Figure S4Di_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(kidneymethsigl2fc[abs(kidneymethsigl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "KIDNEY DMR L2FC")
dev.off()

png(file = "Supplemental Figure S4Dii_112624.png",width = 6,height = 5,units = "in",res = 600)
pheatmap(kidneymethsigl2fc,breaks = seq(-4,4,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

png(file = "Supplemental Figure S4Diii_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(kidneymethl2fc[abs(kidneymethl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "KIDNEY All Sites L2FC")
dev.off()

png(file = "Supplemental Figure S4Ei_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(livermethsigl2fc[abs(livermethsigl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "LIVER DMR L2FC")
dev.off()

png(file = "Supplemental Figure S4Eii_112624.png",width = 6,height = 5,units = "in",res = 600)
pheatmap(livermethsigl2fc,breaks = seq(-4,4,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

png(file = "Supplemental Figure S4Eiii_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(livermethl2fc[abs(livermethl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "LIVER All Sites L2FC")
dev.off()

png(file = "Supplemental Figure S4Fi_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(lungmethsigl2fc[abs(lungmethsigl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "LUNG DMR L2FC")
dev.off()

png(file = "Supplemental Figure S4Fii_112624.png",width = 6,height = 5,units = "in",res = 600)
pheatmap(lungmethsigl2fc,breaks = seq(-4,4,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

png(file = "Supplemental Figure S4Fiii_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(lungmethl2fc[abs(lungmethl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "LUNG All Sites L2FC")
dev.off()

png(file = "Supplemental Figure S4Gi_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(brownmethsigl2fc[abs(brownmethsigl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "BAT DMR L2FC")
dev.off()

png(file = "Supplemental Figure S4Gii_112624.png",width = 6,height = 5,units = "in",res = 600)
pheatmap(brownmethsigl2fc,breaks = seq(-4,4,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

png(file = "Supplemental Figure S4Giii_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(brownmethl2fc[abs(brownmethl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "BAT All Sites L2FC")
dev.off()

png(file = "Supplemental Figure S4Hi_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(whitemethsigl2fc[abs(whitemethsigl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "WAT-SC DMR L2FC")
dev.off()

png(file = "Supplemental Figure S4Hii_112624.png",width = 6,height = 5,units = "in",res = 600)
pheatmap(whitemethsigl2fc,breaks = seq(-4,4,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12,border_color = NA)
dev.off()

png(file = "Supplemental Figure S4Hiii_112624.png",width = 6,height = 5,units = "in",res = 600)
hist(whitemethl2fc[abs(whitemethl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "WAT-SC All Sites L2FC")
dev.off()


save.image("FigureS2_S3_S4_112624.RData")