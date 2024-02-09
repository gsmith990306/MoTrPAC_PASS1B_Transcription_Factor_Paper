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

load("rnasigl2fcmat.RData")

#####
# Supplemental Figure S2
####

pdf(file = "Figure S2Aii.pdf",width = 6,height = 5)
pheatmap(gastrosigrnal2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

pdf(file = "Figure S2Ai.pdf",width = 6,height = 5)
hist(gastrosigrnal2fc[abs(gastrosigrnal2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "SKM-GN DEG L2FC")
dev.off()

pdf(file = "Figure S2Bii.pdf",width = 6,height = 5)
pheatmap(heartsigrnal2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

pdf(file = "Figure S2Bi.pdf",width = 6,height = 5)
hist(heartsigrnal2fc[abs(heartsigrnal2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "HEART DEG L2FC")
dev.off()

pdf(file = "Figure SCii.pdf",width = 6,height = 5)
pheatmap(hipposigrnal2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

pdf(file = "Figure S2Ci.pdf",width = 6,height = 5)
hist(hipposigrnal2fc[abs(hipposigrnal2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "HIPPOC DEG L2FC")
dev.off()

pdf(file = "Figure S2Dii.pdf",width = 6,height = 5)
pheatmap(kidneysigrnal2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

pdf(file = "Figure S2Di.pdf",width = 6,height = 5)
hist(kidneysigrnal2fc[abs(kidneysigrnal2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "KIDNEY DEG L2FC")
dev.off()

pdf(file = "Figure S2Eii.pdf",width = 6,height = 5)
pheatmap(liversigrnal2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

pdf(file = "Figure S2Ei.pdf",width = 6,height = 5)
hist(liversigrnal2fc[abs(liversigrnal2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "LIVER DEG L2FC")
dev.off()

pdf(file = "Figure S2Fii.pdf",width = 6,height = 5)
pheatmap(lungsigrnal2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

pdf(file = "Figure S2Fi.pdf",width = 6,height = 5)
hist(lungsigrnal2fc[abs(lungsigrnal2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "LUNG DEG L2FC")
dev.off()

pdf(file = "Figure S2Gii.pdf",width = 6,height = 5)
pheatmap(brownsigrnal2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

pdf(file = "Figure S2Gi.pdf",width = 6,height = 5)
hist(brownsigrnal2fc[abs(brownsigrnal2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "BAT DEG L2FC")
dev.off()

pdf(file = "Figure S2Hii.pdf",width = 6,height = 5)
pheatmap(whitesigrnal2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

pdf(file = "Figure S2Hi.pdf",width = 6,height = 5)
hist(whitesigrnal2fc[abs(whitesigrnal2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "WAT-SC DEG L2FC")
dev.off()


#####
# Supplemental Figure S3
####

load("atacsigl2fcmat.RData")

pdf(file = "Figure S3Aii.pdf",width = 6,height = 5)
pheatmap(gastrosigatacl2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

pdf(file = "Figure S3Ai.pdf",width = 6,height = 5)
hist(gastrosigatacl2fc[abs(gastrosigatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "SKM-GN DAR L2FC")
dev.off()

pdf(file = "Figure S3Bii.pdf",width = 6,height = 5)
pheatmap(heartsigatacl2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

pdf(file = "Figure S3Bi.pdf",width = 6,height = 5)
hist(heartsigatacl2fc[abs(heartsigatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "HEART DAR L2FC")
dev.off()

pdf(file = "Figure S3Cii.pdf",width = 6,height = 5)
pheatmap(hipposigatacl2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

pdf(file = "Figure S3Ci.pdf",width = 6,height = 5)
hist(hipposigatacl2fc[abs(hipposigatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "HIPPOC DAR L2FC")
dev.off()

pdf(file = "Figure S3Dii.pdf",width = 6,height = 5)
pheatmap(kidneysigatacl2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

pdf(file = "Figure S3Di.pdf",width = 6,height = 5)
hist(kidneysigatacl2fc[abs(kidneysigatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "KIDNEY DAR L2FC")
dev.off()

pdf(file = "Figure S3Eii.pdf",width = 6,height = 5)
pheatmap(liversigatacl2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

pdf(file = "Figure S3Ei.pdf",width = 6,height = 5)
hist(liversigatacl2fc[abs(liversigatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "LIVER DAR L2FC")
dev.off()

pdf(file = "Figure S3Fii.pdf",width = 6,height = 5)
pheatmap(lungsigatacl2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

pdf(file = "Figure S3Fi.pdf",width = 6,height = 5)
hist(lungsigatacl2fc[abs(lungsigatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "LUNG DAR L2FC")
dev.off()

pdf(file = "Figure S3Gii.pdf",width = 6,height = 5)
pheatmap(brownsigatacl2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

pdf(file = "Figure S3Gi.pdf",width = 6,height = 5)
hist(brownsigatacl2fc[abs(brownsigatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "BAT DAR L2FC")
dev.off()

pdf(file = "Figure S3Hii.pdf",width = 6,height = 5)
pheatmap(whitesigatacl2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

pdf(file = "Figure S3Hi.pdf",width = 6,height = 5)
hist(whitesigatacl2fc[abs(whitesigatacl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "WAT-SC DAR L2FC")
dev.off()

#####
# Supplemental Figure S4
####

load("methsigl2fcmat.RData")

pdf(file = "Figure S4Ai.pdf",width = 6,height = 5)
hist(gastromethl2fc[abs(gastromethl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "SKM-GN DMR L2FC")
dev.off()

pdf(file = "Figure S4Aii.pdf",width = 6,height = 5)
pheatmap(gastromethl2fc,breaks = seq(-4,4,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

png(file = "Figure S4Ai.png",width = 6,height = 5,units = "in",res = 600)
hist(gastromethl2fc[abs(gastromethl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "SKM-GN DMR L2FC")
dev.off()

png(file = "Figure S4Aii.png",width = 6,height = 5,units = "in",res = 600)
pheatmap(gastromethl2fc,breaks = seq(-4,4,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

pdf(file = "Figure S4Bi.pdf",width = 6,height = 5)
hist(heartmethl2fc[abs(heartmethl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "HEART DMR L2FC")
dev.off()

pdf(file = "Figure S4Bii.pdf",width = 6,height = 5)
pheatmap(heartmethl2fc,breaks = seq(-4,4,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

png(file = "Figure S4Bi.png",width = 6,height = 5,units = "in",res = 600)
hist(heartmethl2fc[abs(heartmethl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "HEART DMR L2FC")
dev.off()

png(file = "Figure S4Bii.png",width = 6,height = 5,units = "in",res = 600)
pheatmap(heartmethl2fc,breaks = seq(-4,4,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

pdf(file = "Figure S4Ci.pdf",width = 6,height = 5)
hist(hippomethl2fc[abs(hippomethl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "HIPPOC DMR L2FC")
dev.off()

pdf(file = "Figure S4Cii.pdf",width = 6,height = 5)
pheatmap(hippomethl2fc,breaks = seq(-4,4,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

png(file = "Figure S4Ci.png",width = 6,height = 5,units = "in",res = 600)
hist(hippomethl2fc[abs(hippomethl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "HIPPOC DMR L2FC")
dev.off()

png(file = "Figure S4Cii.png",width = 6,height = 5,units = "in",res = 600)
pheatmap(hippomethl2fc,breaks = seq(-4,4,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

pdf(file = "Figure S4Di.pdf",width = 6,height = 5)
hist(kidneymethl2fc[abs(kidneymethl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "KIDNEY DMR L2FC")
dev.off()

pdf(file = "Figure S4Dii.pdf",width = 6,height = 5)
pheatmap(kidneymethl2fc,breaks = seq(-4,4,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

png(file = "Figure S4Di.png",width = 6,height = 5,units = "in",res = 600)
hist(kidneymethl2fc[abs(kidneymethl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "KIDNEY DMR L2FC")
dev.off()

png(file = "Figure S4Dii.png",width = 6,height = 5,units = "in",res = 600)
pheatmap(kidneymethl2fc,breaks = seq(-4,4,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

pdf(file = "Figure S4Ei.pdf",width = 6,height = 5)
hist(livermethl2fc[abs(livermethl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "LIVER DMR L2FC")
dev.off()

pdf(file = "Figure S4Eii.pdf",width = 6,height = 5)
pheatmap(livermethl2fc,breaks = seq(-4,4,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

png(file = "Figure S4Ei.png",width = 6,height = 5,units = "in",res = 600)
hist(livermethl2fc[abs(livermethl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "LIVER DMR L2FC")
dev.off()

png(file = "Figure S4Eii.png",width = 6,height = 5,units = "in",res = 600)
pheatmap(livermethl2fc,breaks = seq(-4,4,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

pdf(file = "Figure S4Fi.pdf",width = 6,height = 5)
hist(lungmethl2fc[abs(lungmethl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "LUNG DMR L2FC")
dev.off()

pdf(file = "Figure S4Fii.pdf",width = 6,height = 5)
pheatmap(lungmethl2fc,breaks = seq(-4,4,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

png(file = "Figure S4Fi.png",width = 6,height = 5,units = "in",res = 600)
hist(lungmethl2fc[abs(lungmethl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "LUNG DMR L2FC")
dev.off()

png(file = "Figure S4Fii.png",width = 6,height = 5,units = "in",res = 600)
pheatmap(lungmethl2fc,breaks = seq(-4,4,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

pdf(file = "Figure S4Gi.pdf",width = 6,height = 5)
hist(brownmethl2fc[abs(brownmethl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "BAT DMR L2FC")
dev.off()

pdf(file = "Figure S4Gii.pdf",width = 6,height = 5)
pheatmap(brownmethl2fc,breaks = seq(-4,4,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

png(file = "Figure S4Gi.png",width = 6,height = 5,units = "in",res = 600)
hist(brownmethl2fc[abs(brownmethl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "BAT DMR L2FC")
dev.off()

png(file = "Figure S4Gii.png",width = 6,height = 5,units = "in",res = 600)
pheatmap(brownmethl2fc,breaks = seq(-4,4,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

pdf(file = "Figure S4Hi.pdf",width = 6,height = 5)
hist(whitemethl2fc[abs(whitemethl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "WAT-SC DMR L2FC")
dev.off()

pdf(file = "Figure S4Hii.pdf",width = 6,height = 5)
pheatmap(whitemethl2fc,breaks = seq(-4,4,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

png(file = "Figure S4Hi.png",width = 6,height = 5,units = "in",res = 600)
hist(whitemethl2fc[abs(whitemethl2fc) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "WAT-SC DMR L2FC")
dev.off()

png(file = "Figure S4Hii.png",width = 6,height = 5,units = "in",res = 600)
pheatmap(whitemethl2fc,breaks = seq(-4,4,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

save.image("FigureS2_S3_S4.RData")