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
load("omesigdata.RData")

# Supplemental Figure S1
tissuedegintersectmat <- matrix(0L,nrow = 8,ncol = 8)
rownames(tissuedegintersectmat) <- c("SKM-GN","HEART","HIPPOC","KIDNEY",
                                     "LIVER","LUNG","BAT","WAT-SC")
colnames(tissuedegintersectmat) <- c("SKM-GN","HEART","HIPPOC","KIDNEY",
                                     "LIVER","LUNG","BAT","WAT-SC")

tissuedarintersectmat <- matrix(0L,nrow = 8,ncol = 8)
rownames(tissuedarintersectmat) <- c("SKM-GN","HEART","HIPPOC","KIDNEY",
                                     "LIVER","LUNG","BAT","WAT-SC")
colnames(tissuedarintersectmat) <- c("SKM-GN","HEART","HIPPOC","KIDNEY",
                                     "LIVER","LUNG","BAT","WAT-SC")

tissuedegintersectmatfrac <- matrix(0L,nrow = 8,ncol = 8)
rownames(tissuedegintersectmatfrac) <- c("SKM-GN","HEART","HIPPOC","KIDNEY",
                                         "LIVER","LUNG","BAT","WAT-SC")
colnames(tissuedegintersectmatfrac) <- c("SKM-GN","HEART","HIPPOC","KIDNEY",
                                         "LIVER","LUNG","BAT","WAT-SC")

tissuedarintersectmatfrac <- matrix(0L,nrow = 8,ncol = 8)
rownames(tissuedarintersectmatfrac) <- c("SKM-GN","HEART","HIPPOC","KIDNEY",
                                         "LIVER","LUNG","BAT","WAT-SC")
colnames(tissuedarintersectmatfrac) <- c("SKM-GN","HEART","HIPPOC","KIDNEY",
                                         "LIVER","LUNG","BAT","WAT-SC")


tissuedeglist <- list(gastrornasig,
                      heartrnasig,
                      hippornasig,
                      kidneyrnasig,
                      liverrnasig,
                      lungrnasig,
                      brownrnasig,
                      whiternasig)
tissuedarlist <- list(gastroatacsig,
                      heartatacsig,
                      hippoatacsig,
                      kidneyatacsig,
                      liveratacsig,
                      lungatacsig,
                      brownatacsig,
                      whiteatacsig)

for(i in 1:8){
  for(j in 1:8){
    tissuedegintersectmat[i,j] <- length(intersect(tissuedeglist[[i]],tissuedeglist[[j]]))
    tissuedarintersectmat[i,j] <- length(intersect(tissuedarlist[[i]],tissuedarlist[[j]]))
    
    tissuedegintersectmatfrac[i,j] <- length(intersect(tissuedeglist[[i]],tissuedeglist[[j]]))/length(union(tissuedeglist[[i]],tissuedeglist[[j]]))
    tissuedarintersectmatfrac[i,j] <- length(intersect(tissuedarlist[[i]],tissuedarlist[[j]]))/length(union(tissuedarlist[[i]],tissuedarlist[[j]]))
  }
}

tissuedmrintersectmatfrac <- matrix(0L,nrow = 8,ncol = 8)
rownames(tissuedmrintersectmatfrac) <- c("SKM-GN","HEART","HIPPOC","KIDNEY",
                                         "LIVER","LUNG","BAT","WAT-SC")
colnames(tissuedmrintersectmatfrac) <- c("SKM-GN","HEART","HIPPOC","KIDNEY",
                                         "LIVER","LUNG","BAT","WAT-SC")

tissuedmrlist <- list(gastromethsig,
                      heartmethsig,
                      hippomethsig,
                      kidneymethsig,
                      livermethsig,
                      lungmethsig,
                      brownmethsig,
                      whitemethsig)

for(i in 1:8){
  for(j in 1:8){
    tissuedmrintersectmatfrac[i,j] <- length(intersect(tissuedmrlist[[i]],tissuedmrlist[[j]]))/length(union(tissuedmrlist[[i]],tissuedmrlist[[j]]))
  }
}

tissuemeta <- data.frame(row.names = rownames(tissuedmrintersectmatfrac),
                         "Tissue" = rownames(tissuedmrintersectmatfrac))

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
                 "Group" = c("control" = "white",
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

pdf(file = "Figure S1A.pdf", width=7, height=6)
pheatmap(tissuedegintersectmatfrac,breaks = seq(0,1,length.out = 101),color = colorpanel(101,"white","firebrick","firebrick"),display_numbers = T,number_color = "black",annotation_row = tissuemeta,annotation_col = tissuemeta,annotation_colors = ann_cols,angle_col = 315,annotation_names_col = F,annotation_names_row = F,fontsize = 15,annotation_legend = F)
dev.off()

pdf(file = "Figure S1B.pdf", width=7, height=6)
pheatmap(tissuedarintersectmatfrac,breaks = seq(0,1,length.out = 101),color = colorpanel(101,"white","firebrick","firebrick"),display_numbers = T,number_color = "black",annotation_row = tissuemeta,annotation_col = tissuemeta,annotation_colors = ann_cols,angle_col = 315,annotation_names_col = F,annotation_names_row = F,fontsize = 15,annotation_legend = F)
dev.off()

png(file = "Figure S1A.png", width=7, height=6,units = "in",res = 600)
pheatmap(tissuedegintersectmatfrac,breaks = seq(0,1,length.out = 101),color = colorpanel(101,"white","firebrick","firebrick"),display_numbers = T,number_color = "black",annotation_row = tissuemeta,annotation_col = tissuemeta,annotation_colors = ann_cols,angle_col = 315,annotation_names_col = F,annotation_names_row = F,fontsize = 15,annotation_legend = F)
dev.off()

png(file = "Figure S1B.png", width=7, height=6,units = "in",res = 600)
pheatmap(tissuedarintersectmatfrac,breaks = seq(0,1,length.out = 101),color = colorpanel(101,"white","firebrick","firebrick"),display_numbers = T,number_color = "black",annotation_row = tissuemeta,annotation_col = tissuemeta,annotation_colors = ann_cols,angle_col = 315,annotation_names_col = F,annotation_names_row = F,fontsize = 15,annotation_legend = F)
dev.off()

pdf(file = "Figure S1C.pdf", width=7, height=6)
pheatmap(tissuedmrintersectmatfrac,breaks = seq(0,1,length.out = 101),color = colorpanel(101,"white","firebrick","firebrick"),display_numbers = T,number_color = "black",annotation_row = tissuemeta,annotation_col = tissuemeta,annotation_colors = ann_cols,angle_col = 315,annotation_names_col = F,annotation_names_row = F,fontsize = 15,annotation_legend = F)
dev.off()

png(file = "Figure S1C.png", width=7, height=6,units = "in",res = 600)
pheatmap(tissuedmrintersectmatfrac,breaks = seq(0,1,length.out = 101),color = colorpanel(101,"white","firebrick","firebrick"),display_numbers = T,number_color = "black",annotation_row = tissuemeta,annotation_col = tissuemeta,annotation_colors = ann_cols,angle_col = 315,annotation_names_col = F,annotation_names_row = F,fontsize = 15,annotation_legend = F)
dev.off()


save.image("FigureS1.RData")