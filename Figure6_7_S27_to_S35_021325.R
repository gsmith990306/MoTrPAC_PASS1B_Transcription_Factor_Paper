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

tfanno <- readRDS("tfanno.RDS")
load("rnacontrolnormandmeta.RData")

#####
# Figure 6A
####

# Homer TF enrichment among DARs in each tissue
gastroatacsigtf50 <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/gastroatacsigfilejun2022out/knownResults.csv",header = T,row.names = 1)
heartatacsigtf50 <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/heartatacsigfilejun2022out/knownResults.csv",header = T,row.names = 1)
kidneyatacsigtf50 <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/kidneyatacsigfilejun2022out/knownResults.csv",header = T,row.names = 1)
liveratacsigtf50 <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/liveratacsigfilejun2022out/knownResults.csv",header = T,row.names = 1)
lungatacsigtf50 <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/lungatacsigfilejun2022out/knownResults.csv",header = T,row.names = 1)
brownatacsigtf50 <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/brownatacsigfilejun2022out/knownResults.csv",header = T,row.names = 1)

atacsigtflist <- Reduce(intersect,list(rownames(gastroatacsigtf50),rownames(heartatacsigtf50),rownames(kidneyatacsigtf50)))
atacsigtflabel <- gsub("\\(.*","",atacsigtflist)


atacsigtoptfs50 <- Reduce(union,list(rownames(gastroatacsigtf50)[1:10],
                                     rownames(heartatacsigtf50)[1:10],
                                     rownames(kidneyatacsigtf50)[1:10],
                                     rownames(liveratacsigtf50)[1:10],
                                     rownames(lungatacsigtf50)[1:10],
                                     rownames(brownatacsigtf50)[1:10]))

trimatacsigtoptfs50ens <- unique(tfanno[atacsigtoptfs50,"Ensembl"])[!is.na(unique(tfanno[atacsigtoptfs50,"Ensembl"]))]
trimatacsigtoptfs <- intersect(atacsigtoptfs50,rownames(tfanno[tfanno$Ensembl %in% trimatacsigtoptfs50ens,]))

tfatacsigpmat <- -1*cbind(gastroatacsigtf50[atacsigtflist,"Log.P.value"],
                          heartatacsigtf50[atacsigtflist,"Log.P.value"],
                          kidneyatacsigtf50[atacsigtflist,"Log.P.value"],
                          liveratacsigtf50[atacsigtflist,"Log.P.value"],
                          lungatacsigtf50[atacsigtflist,"Log.P.value"],
                          brownatacsigtf50[atacsigtflist,"Log.P.value"])
rownames(tfatacsigpmat) <- atacsigtflist
colnames(tfatacsigpmat) <- c("SKM-GN","HEART","KIDNEY","LIVER","LUNG","BAT")


tfatacsigpctmat <- cbind(as.numeric(gsub("%","",gastroatacsigtf50[atacsigtflist,"X..of.Target.Sequences.with.Motif"])),
                         as.numeric(gsub("%","",heartatacsigtf50[atacsigtflist,"X..of.Target.Sequences.with.Motif"])),
                         as.numeric(gsub("%","",kidneyatacsigtf50[atacsigtflist,"X..of.Target.Sequences.with.Motif"])),
                         as.numeric(gsub("%","",liveratacsigtf50[atacsigtflist,"X..of.Target.Sequences.with.Motif"])),
                         as.numeric(gsub("%","",lungatacsigtf50[atacsigtflist,"X..of.Target.Sequences.with.Motif"])),
                         as.numeric(gsub("%","",brownatacsigtf50[atacsigtflist,"X..of.Target.Sequences.with.Motif"])))
rownames(tfatacsigpctmat) <- atacsigtflist
colnames(tfatacsigpctmat) <- c("SKM-GN","HEART","KIDNEY","LIVER","LUNG","BAT")

finalatacsigtfgenes <- intersect(tfanno[trimatacsigtoptfs,"Ensembl"],rownames(rnacontrolnorm))
finalatacsigtfs <- intersect(trimatacsigtoptfs,rownames(tfanno[tfanno$Ensembl %in% finalatacsigtfgenes,]))

tissuemeta <- data.frame(row.names = c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC"),
                         "Tissue" = c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC"))

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


# Figure 6A
pdf(file = "Figure 6A_021325.pdf",width = 6,height = 8)
pheatmap(tfatacsigpmat[finalatacsigtfs,],cluster_rows = F,cluster_cols = F,labels_row = gsub("\\(.*","",finalatacsigtfs),annotation_col = tissuemeta,annotation_colors = ann_cols,show_colnames = F,color = colorpanel(101,"white","firebrick"),border_color = NA)
dev.off()

png(file = "Figure 6A_021325.png",width = 6,height = 8,units = "in",res = 600)
pheatmap(tfatacsigpmat[finalatacsigtfs,],cluster_rows = F,cluster_cols = F,labels_row = gsub("\\(.*","",finalatacsigtfs),annotation_col = tissuemeta,annotation_colors = ann_cols,show_colnames = F,color = colorpanel(101,"white","firebrick"),border_color = NA)
dev.off()

#####
# Figure 6B
####

gastromethsigtf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/gastromethsigannopeakout/knownResults.csv",header = T,row.names = 1)
gastromethsigpromproxtf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/gastropromproxmethsigannopeakout/knownResults.csv",header = T,row.names = 1)
gastromethsiginttf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/gastrointmethsigannopeakout/knownResults.csv",header = T,row.names = 1)
gastromethsigdisttf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/gastrodistmethsigannopeakout/knownResults.csv",header = T,row.names = 1)
gastromethsigdowntf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/gastrodownmethsigannopeakout/knownResults.csv",header = T,row.names = 1)
gastromethsigextf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/gastroexmethsigannopeakout/knownResults.csv",header = T,row.names = 1)

heartmethsigtf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/heartmethsigannopeakout/knownResults.csv",header = T,row.names = 1)
heartmethsigpromproxtf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/heartpromproxmethsigannopeakout/knownResults.csv",header = T,row.names = 1)
heartmethsiginttf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/heartintmethsigannopeakout/knownResults.csv",header = T,row.names = 1)
heartmethsigdisttf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/heartdistmethsigannopeakout/knownResults.csv",header = T,row.names = 1)
heartmethsigdowntf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/heartdownmethsigannopeakout/knownResults.csv",header = T,row.names = 1)
heartmethsigextf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/heartexmethsigannopeakout/knownResults.csv",header = T,row.names = 1)

hippomethsigtf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/hippomethsigannopeakout/knownResults.csv",header = T,row.names = 1)
hippomethsigpromproxtf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/hippopromproxmethsigannopeakout/knownResults.csv",header = T,row.names = 1)
hippomethsiginttf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/hippointmethsigannopeakout/knownResults.csv",header = T,row.names = 1)
hippomethsigdisttf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/hippodistmethsigannopeakout/knownResults.csv",header = T,row.names = 1)
hippomethsigdowntf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/hippodownmethsigannopeakout/knownResults.csv",header = T,row.names = 1)
hippomethsigextf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/hippoexmethsigannopeakout/knownResults.csv",header = T,row.names = 1)

kidneymethsigtf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/kidneymethsigannopeakout/knownResults.csv",header = T,row.names = 1)
kidneymethsigpromproxtf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/kidneypromproxmethsigannopeakout/knownResults.csv",header = T,row.names = 1)
kidneymethsiginttf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/kidneyintmethsigannopeakout/knownResults.csv",header = T,row.names = 1)
kidneymethsigdisttf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/kidneydistmethsigannopeakout/knownResults.csv",header = T,row.names = 1)
kidneymethsigdowntf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/kidneydownmethsigannopeakout/knownResults.csv",header = T,row.names = 1)
kidneymethsigextf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/kidneyexmethsigannopeakout/knownResults.csv",header = T,row.names = 1)

livermethsigtf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/livermethsigannopeakout/knownResults.csv",header = T,row.names = 1)
livermethsigpromproxtf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/liverpromproxmethsigannopeakout/knownResults.csv",header = T,row.names = 1)
livermethsiginttf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/liverintmethsigannopeakout/knownResults.csv",header = T,row.names = 1)
livermethsigdisttf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/liverdistmethsigannopeakout/knownResults.csv",header = T,row.names = 1)
livermethsigdowntf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/liverdownmethsigannopeakout/knownResults.csv",header = T,row.names = 1)
livermethsigextf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/liverexmethsigannopeakout/knownResults.csv",header = T,row.names = 1)

lungmethsigtf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/lungmethsigannopeakout/knownResults.csv",header = T,row.names = 1)
lungmethsigpromproxtf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/lungpromproxmethsigannopeakout/knownResults.csv",header = T,row.names = 1)
lungmethsiginttf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/lungintmethsigannopeakout/knownResults.csv",header = T,row.names = 1)
lungmethsigdisttf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/lungdistmethsigannopeakout/knownResults.csv",header = T,row.names = 1)
lungmethsigdowntf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/lungdownmethsigannopeakout/knownResults.csv",header = T,row.names = 1)
lungmethsigextf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/lungexmethsigannopeakout/knownResults.csv",header = T,row.names = 1)

brownmethsigtf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/brownmethsigannopeakout/knownResults.csv",header = T,row.names = 1)
brownmethsigpromproxtf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/brownpromproxmethsigannopeakout/knownResults.csv",header = T,row.names = 1)
brownmethsiginttf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/brownintmethsigannopeakout/knownResults.csv",header = T,row.names = 1)
brownmethsigdisttf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/browndistmethsigannopeakout/knownResults.csv",header = T,row.names = 1)
brownmethsigdowntf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/browndownmethsigannopeakout/knownResults.csv",header = T,row.names = 1)
brownmethsigextf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/brownexmethsigannopeakout/knownResults.csv",header = T,row.names = 1)

whitemethsigtf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/whitemethsigannopeakout/knownResults.csv",header = T,row.names = 1)
whitemethsigpromproxtf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/whitepromproxmethsigannopeakout/knownResults.csv",header = T,row.names = 1)
whitemethsiginttf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/whiteintmethsigannopeakout/knownResults.csv",header = T,row.names = 1)
whitemethsigdisttf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/whitedistmethsigannopeakout/knownResults.csv",header = T,row.names = 1)
whitemethsigdowntf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/whitedownmethsigannopeakout/knownResults.csv",header = T,row.names = 1)
whitemethsigextf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/whiteexmethsigannopeakout/knownResults.csv",header = T,row.names = 1)

methsigtflist <- Reduce(intersect,list(rownames(gastromethsigtf),rownames(heartmethsigtf),rownames(kidneymethsigtf)))
methsigtflabel <- gsub("\\(.*","",methsigtflist)

methsigtoptfs <- Reduce(union,list(rownames(gastromethsigtf)[1:10],
                                   rownames(heartmethsigtf)[1:10],
                                   rownames(hippomethsigtf)[1:10],
                                   rownames(kidneymethsigtf)[1:10],
                                   rownames(livermethsigtf)[1:10],
                                   rownames(lungmethsigtf)[1:10],
                                   rownames(brownmethsigtf)[1:10],
                                   rownames(whitemethsigtf)[1:10]))

tfanno <- read.csv("PASS1B Transcription Factor Paper Data/TF Data/tflist.csv",header = T,row.names = 1)
tfanno$Ensembl <- gsub(" ","",tfanno$Ensembl)
tfanno <- tfanno[!tfanno$Ensembl == "",]


trimmethsigtoptfsens <- unique(tfanno[methsigtoptfs,"Ensembl"])[!is.na(unique(tfanno[methsigtoptfs,"Ensembl"]))]
trimmethsigtoptfs <- intersect(methsigtoptfs,rownames(tfanno[tfanno$Ensembl %in% trimmethsigtoptfsens,]))

tfmethsigpmat <- -1*cbind(gastromethsigtf[methsigtflist,"Log.P.value"],
                          heartmethsigtf[methsigtflist,"Log.P.value"],
                          hippomethsigtf[methsigtflist,"Log.P.value"],
                          kidneymethsigtf[methsigtflist,"Log.P.value"],
                          livermethsigtf[methsigtflist,"Log.P.value"],
                          lungmethsigtf[methsigtflist,"Log.P.value"],
                          brownmethsigtf[methsigtflist,"Log.P.value"],
                          whitemethsigtf[methsigtflist,"Log.P.value"])
rownames(tfmethsigpmat) <- methsigtflist
colnames(tfmethsigpmat) <- c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC")


tfmethsigpctmat <- cbind(as.numeric(gsub("%","",gastromethsigtf[methsigtflist,"X..of.Target.Sequences.with.Motif"])),
                         as.numeric(gsub("%","",heartmethsigtf[methsigtflist,"X..of.Target.Sequences.with.Motif"])),
                         as.numeric(gsub("%","",hippomethsigtf[methsigtflist,"X..of.Target.Sequences.with.Motif"])),
                         as.numeric(gsub("%","",kidneymethsigtf[methsigtflist,"X..of.Target.Sequences.with.Motif"])),
                         as.numeric(gsub("%","",livermethsigtf[methsigtflist,"X..of.Target.Sequences.with.Motif"])),
                         as.numeric(gsub("%","",lungmethsigtf[methsigtflist,"X..of.Target.Sequences.with.Motif"])),
                         as.numeric(gsub("%","",brownmethsigtf[methsigtflist,"X..of.Target.Sequences.with.Motif"])),
                         as.numeric(gsub("%","",whitemethsigtf[methsigtflist,"X..of.Target.Sequences.with.Motif"])))
rownames(tfmethsigpctmat) <- methsigtflist
colnames(tfmethsigpctmat) <- c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC")

finalmethsigtfgenes <- intersect(tfanno[trimmethsigtoptfs,"Ensembl"],rownames(rnacontrolnorm))
finalmethsigtfs <- intersect(trimmethsigtoptfs,rownames(tfanno[tfanno$Ensembl %in% finalmethsigtfgenes,]))

# Figure 6B
pdf(file = "Figure 6B_021325.pdf",width = 6,height = 8)
pheatmap(tfmethsigpmat[finalmethsigtfs,],cluster_rows = F,cluster_cols = F,labels_row = gsub("\\(.*","",finalmethsigtfs),annotation_col = tissuemeta,annotation_colors = ann_cols,show_colnames = F,color = colorpanel(101,"white","firebrick"),border_color = NA)
dev.off()

png(file = "Figure 6B_021325.png",width = 6,height = 8,units = "in",res = 600)
pheatmap(tfmethsigpmat[finalmethsigtfs,],cluster_rows = F,cluster_cols = F,labels_row = gsub("\\(.*","",finalmethsigtfs),annotation_col = tissuemeta,annotation_colors = ann_cols,show_colnames = F,color = colorpanel(101,"white","firebrick"),border_color = NA)
dev.off()



#####
# Figure 6C
####

gastrotf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/gastrodegpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
gastropromtf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/gastrodegprompeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
gastrodisttf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/gastrodegdistpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
gastrointtf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/gastrodegintpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
gastroextf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/gastrodegexpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
gastrouptf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/gastrodeguppeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
gastrodowntf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/gastrodegdownpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)

hearttf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/heartdegpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
heartpromtf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/heartdegprompeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
heartdisttf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/heartdegdistpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
heartinttf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/heartdegintpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
heartextf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/heartdegexpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
heartuptf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/heartdeguppeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
heartdowntf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/heartdegdownpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)

hippotf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/hippodegpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
hippopromtf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/hippodegprompeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
hippodisttf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/hippodegdistpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
hippointtf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/hippodegintpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
hippoextf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/hippodegexpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
hippouptf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/hippodeguppeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
hippodowntf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/hippodegdownpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)

kidneytf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/kidneydegpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
kidneypromtf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/kidneydegprompeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
kidneydisttf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/kidneydegdistpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
kidneyinttf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/kidneydegintpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
kidneyextf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/kidneydegexpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
kidneyuptf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/kidneydeguppeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
kidneydowntf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/kidneydegdownpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)

livertf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/liverdegpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
liverpromtf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/liverdegprompeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
liverdisttf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/liverdegdistpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
liverinttf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/liverdegintpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
liverextf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/liverdegexpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
liveruptf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/liverdeguppeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
liverdowntf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/liverdegdownpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)

lungtf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/lungdegpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
lungpromtf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/lungdegprompeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
lungdisttf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/lungdegdistpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
lunginttf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/lungdegintpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
lungextf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/lungdegexpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
lunguptf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/lungdeguppeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
lungdowntf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/lungdegdownpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)

browntf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/browndegpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
brownpromtf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/browndegprompeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
browndisttf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/browndegdistpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
browninttf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/browndegintpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
brownextf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/browndegexpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
brownuptf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/browndeguppeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
browndowntf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/browndegdownpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)

whitetf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/whitedegpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
whitepromtf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/whitedegprompeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
whitedisttf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/whitedegdistpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
whiteinttf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/whitedegintpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
whiteextf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/whitedegexpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
whiteuptf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/whitedeguppeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
whitedowntf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/whitedegdownpeakfilejun2022out/knownResults.csv",header = T,row.names = 1)

gastro_activetf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/gastro_activepeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
heart_activetf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/heart_activepeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
hippo_activetf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/hippo_activepeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
kidney_activetf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/kidney_activepeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
liver_activetf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/liver_activepeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
lung_activetf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/lung_activepeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
brown_activetf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/brown_activepeakfilejun2022out/knownResults.csv",header = T,row.names = 1)
white_activetf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/white_activepeakfilejun2022out/knownResults.csv",header = T,row.names = 1)


tflist <- Reduce(intersect,list(rownames(gastro_activetf),rownames(heart_activetf),rownames(hippo_activetf)))
tflabel <- gsub("\\(.*","",tflist)

gastrosigpromproxtf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/tfgastrosigpromproxjun2022out/knownResults.csv",header = T,row.names = 1)
gastrosigpromfartf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/tfgastrosigpromfarjun2022out/knownResults.csv",header = T,row.names = 1)
heartsigpromproxtf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/tfheartsigpromproxjun2022out/knownResults.csv",header = T,row.names = 1)
heartsigpromfartf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/tfheartsigpromfarjun2022out/knownResults.csv",header = T,row.names = 1)
hipposigpromproxtf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/tfhipposigpromproxjun2022out/knownResults.csv",header = T,row.names = 1)
hipposigpromfartf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/tfhipposigpromfarjun2022out/knownResults.csv",header = T,row.names = 1)
kidneysigpromproxtf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/tfkidneysigpromproxjun2022out/knownResults.csv",header = T,row.names = 1)
kidneysigpromfartf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/tfkidneysigpromfarjun2022out/knownResults.csv",header = T,row.names = 1)
liversigpromproxtf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/tfliversigpromproxjun2022out/knownResults.csv",header = T,row.names = 1)
liversigpromfartf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/tfliversigpromfarjun2022out/knownResults.csv",header = T,row.names = 1)
lungsigpromproxtf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/tflungsigpromproxjun2022out/knownResults.csv",header = T,row.names = 1)
lungsigpromfartf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/tflungsigpromfarjun2022out/knownResults.csv",header = T,row.names = 1)
whitesigpromproxtf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/tfwhitesigpromproxjun2022out/knownResults.csv",header = T,row.names = 1)
whitesigpromfartf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/tfwhitesigpromfarjun2022out/knownResults.csv",header = T,row.names = 1)
brownsigpromproxtf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/tfbrownsigpromproxjun2022out/knownResults.csv",header = T,row.names = 1)
brownsigpromfartf <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/tfbrownsigpromfarjun2022out/knownResults.csv",header = T,row.names = 1)



tf50logpmat <- -1*cbind(gastrotf[tflist,"Log.P.value"],
                        hearttf[tflist,"Log.P.value"],
                        hippotf[tflist,"Log.P.value"],
                        kidneytf[tflist,"Log.P.value"],
                        livertf[tflist,"Log.P.value"],
                        lungtf[tflist,"Log.P.value"],
                        browntf[tflist,"Log.P.value"],
                        whitetf[tflist,"Log.P.value"])
rownames(tf50logpmat) <- tflist
colnames(tf50logpmat) <- c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC")

# We select the TFs most highly enriched in each of the tissues
sigtflist <- intersect(Reduce(union,list(rownames(gastrotf)[1:10],
                                         rownames(hearttf)[1:10],
                                         rownames(hippotf)[1:10],
                                         rownames(kidneytf)[1:10],
                                         rownames(livertf)[1:10],
                                         rownames(lungtf)[1:10],
                                         rownames(browntf)[1:10],
                                         rownames(whitetf)[1:10])),tflist)
sigtflabel <- gsub("\\(.*","",sigtflist)
# We remove TFs without associated ensembl ids in the data and a duplicate FOXA1 entry
finalsigtflist <- c("Six2(Homeobox)/NephronProgenitor-Six2-ChIP-Seq(GSE39837)/Homer",
                    "Mef2b(MADS)/HEK293-Mef2b.V5-ChIP-Seq(GSE67450)/Homer",
                    "Mef2c(MADS)/GM12878-Mef2c-ChIP-Seq(GSE32465)/Homer",          
                    "Mef2a(MADS)/HL1-Mef2a.biotin-ChIP-Seq(GSE21529)/Homer",         
                    "NF1-halfsite(CTF)/LNCaP-NF1-ChIP-Seq(Unpublished)/Homer",       
                    "Sp2(Zf)/HEK293-Sp2.eGFP-ChIP-Seq(Encode)/Homer",                
                    "Jun-AP1(bZIP)/K562-cJun-ChIP-Seq(GSE31477)/Homer",              
                    "Fosl2(bZIP)/3T3L1-Fosl2-ChIP-Seq(GSE56872)/Homer",              
                    "Fos(bZIP)/TSC-Fos-ChIP-Seq(GSE110950)/Homer",                   
                    "JunB(bZIP)/DendriticCells-Junb-ChIP-Seq(GSE36099)/Homer",       
                    "Fra2(bZIP)/Striatum-Fra2-ChIP-Seq(GSE43429)/Homer",             
                    "KLF14(Zf)/HEK293-KLF14.GFP-ChIP-Seq(GSE58341)/Homer",           
                    "Sox9(HMG)/Limb-SOX9-ChIP-Seq(GSE73225)/Homer",                  
                    "GLIS3(Zf)/Thyroid-Glis3.GFP-ChIP-Seq(GSE103297)/Homer",         
                    "BATF(bZIP)/Th17-BATF-ChIP-Seq(GSE39756)/Homer",                 
                    "Atf3(bZIP)/GBM-ATF3-ChIP-Seq(GSE33912)/Homer",                
                    "Fra1(bZIP)/BT549-Fra1-ChIP-Seq(GSE46166)/Homer",                
                    "Nur77(NR)/K562-NR4A1-ChIP-Seq(GSE31363)/Homer",                 
                    "HNF1b(Homeobox)/PDAC-HNF1B-ChIP-Seq(GSE64557)/Homer",           
                    "Hnf1(Homeobox)/Liver-Foxa2-Chip-Seq(GSE25694)/Homer",           
                    "ERRg(NR)/Kidney-ESRRG-ChIP-Seq(GSE104905)/Homer",               
                    "Zic2(Zf)/ESC-Zic2-ChIP-Seq(SRP197560)/Homer",                   
                    "Zic(Zf)/Cerebellum-ZIC1.2-ChIP-Seq(GSE60731)/Homer",            
                    "Foxa2(Forkhead)/Liver-Foxa2-ChIP-Seq(GSE25694)/Homer",          
                    "Sox10(HMG)/SciaticNerve-Sox3-ChIP-Seq(GSE35132)/Homer",         
                    "Sox3(HMG)/NPC-Sox3-ChIP-Seq(GSE33059)/Homer",                   
                    "Foxo3(Forkhead)/U2OS-Foxo3-ChIP-Seq(E-MTAB-2701)/Homer",        
                    "FoxL2(Forkhead)/Ovary-FoxL2-ChIP-Seq(GSE60858)/Homer",          
                    "Foxf1(Forkhead)/Lung-Foxf1-ChIP-Seq(GSE77951)/Homer",           
                    "FOXA1(Forkhead)/MCF7-FOXA1-ChIP-Seq(GSE26831)/Homer",           
                    "RXR(NR),DR1/3T3L1-RXR-ChIP-Seq(GSE13511)/Homer",                
                    "Fox:Ebox(Forkhead,bHLH)/Panc1-Foxa2-ChIP-Seq(GSE47459)/Homer",  
                    "FOXK1(Forkhead)/HEK293-FOXK1-ChIP-Seq(GSE51673)/Homer",         
                    "FOXM1(Forkhead)/MCF7-FOXM1-ChIP-Seq(GSE72977)/Homer",           
                    "Erra(NR)/HepG2-Erra-ChIP-Seq(GSE31477)/Homer",                  
                    "Foxo1(Forkhead)/RAW-Foxo1-ChIP-Seq(Fan_et_al.)/Homer",          
                    "NFIL3(bZIP)/HepG2-NFIL3-ChIP-Seq(Encode)/Homer",                
                    "SpiB(ETS)/OCILY3-SPIB-ChIP-Seq(GSE56857)/Homer",                
                    "KLF5(Zf)/LoVo-KLF5-ChIP-Seq(GSE49402)/Homer",                   
                    "TEAD3(TEA)/HepG2-TEAD3-ChIP-Seq(Encode)/Homer",                 
                    "Elf4(ETS)/BMDM-Elf4-ChIP-Seq(GSE88699)/Homer",                  
                    "ELF5(ETS)/T47D-ELF5-ChIP-Seq(GSE30407)/Homer",                  
                    "Etv2(ETS)/ES-ER71-ChIP-Seq(GSE59402)/Homer",                    
                    "ETV1(ETS)/GIST48-ETV1-ChIP-Seq(GSE22441)/Homer",                
                    "ERG(ETS)/VCaP-ERG-ChIP-Seq(GSE14097)/Homer",                    
                    "BORIS(Zf)/K562-CTCFL-ChIP-Seq(GSE32465)/Homer",                 
                    "EBF2(EBF)/BrownAdipose-EBF2-ChIP-Seq(GSE97114)/Homer",          
                    "PU.1(ETS)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer",
                    "ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer")

ann_colsheatmaptrim <- list("Tissue" = c("SKM-GN" = "#088c03",
                                         "HEART" = "#f28b2f",
                                         "HIPPOC" = "#bf7534",
                                         "KIDNEY"= "#7553a7",
                                         "LIVER" = "#da6c75",
                                         "LUNG" = "#04bf8a",
                                         "BAT" = "#8c5220",
                                         "WAT-SC" = "#214da6"),
                            "Sex" = c("Female" = "#ff6eff",
                                      "Male" = "#5555ff"),
                            "Week" = c("control" = "white",
                                       "1w" = "#F7FCB9",
                                       "2w" = "#ADDD8E",
                                       "4w" = "#238443",
                                       "8w" = "#002612",
                                       "All" = "purple",
                                       "Background" = "grey"),
                            "Region" = c("Distal Intergenic" = "#BCBD22FF",
                                         "Intron" = "#8C564BFF",
                                         "Promoter" = "#FF7F0EFF",
                                         "All" = "navy"),
                            "Training.Response" = c("Background" = "grey",
                                                    "Down.Reg" = "skyblue2",
                                                    "All.Sig" = "lightgoldenrod2",
                                                    "Up.Reg" = "indianred2"))


# Figure 6C - 592x580
pdf(file = "Figure 6C_021325.pdf",width = 6,height = 8)
pheatmap(tf50logpmat[finalsigtflist,],labels_row = gsub("\\(.*","",finalsigtflist),angle_col = 0,color = colorpanel(101,"white","firebrick"),breaks = seq(0,40,length.out = 101),cluster_cols = F,cluster_rows = F,annotation_col = tissuemeta,annotation_colors = ann_colsheatmaptrim,show_colnames = F,border_color = NA)
dev.off()

png(file = "Figure 6C_021325.png",width = 6,height = 8, units = "in",res = 600)
pheatmap(tf50logpmat[finalsigtflist,],labels_row = gsub("\\(.*","",finalsigtflist),angle_col = 0,color = colorpanel(101,"white","firebrick"),breaks = seq(0,40,length.out = 101),cluster_cols = F,cluster_rows = F,annotation_col = tissuemeta,annotation_colors = ann_colsheatmaptrim,show_colnames = F,border_color = NA)
dev.off()


#####
# Figure 6D
####

allpeakmotifs <- readRDS("allpeakmotifs.RDS")
load("PASS1B Transcription Factor Paper Data/DEA Analysis/pass1b-06_epigen-atac-seq_dea.RData")
peakanno <- readRDS("peakanno.RDS")

atactraining <- epigen_atac_seq$training_dea
atactraining$custom_annotation <- ""
atactraining$custom_annotation <- peakanno[atactraining$feature_ID,"custom_annotation"]

rm(epigen_atac_seq)
gc()

homerPositionAnnotation = allpeakmotifs %>% dplyr::mutate(motif_name=`Motif.Name`) %>% dplyr::select(PositionID, motif_name)
homerPositionAnnotation = unique(homerPositionAnnotation)
atacPeaksAnnotated = atactraining %>% dplyr::left_join(homerPositionAnnotation, by=c("feature_ID"="PositionID"))
pcutoff = 0.05
#atacPeaksAnnotatedWithMotif = atactraining[atactraining$feature_ID %in% allpeakmotifs$PositionID,]
atacPeaksAnnotatedWithMotif = atacPeaksAnnotated %>% dplyr::filter(!is.na(motif_name)) %>% dplyr::mutate(is_sig=adj_p_value<pcutoff)

rvByMotif = lapply(unique(atacPeaksAnnotatedWithMotif$motif_name), function(thisMotif){
  ## go by motif
  message(thisMotif)
  subdat = atacPeaksAnnotatedWithMotif %>% dplyr::filter(motif_name==thisMotif)
  
  ## by tissue
  ans = lapply(unique(subdat$tissue), function(thisTissue){
    message(thisTissue)
    
    subdatByTissue = subdat %>% dplyr::filter(tissue==thisTissue)
    mytab = table(subdatByTissue$custom_annotation, factor(subdatByTissue$is_sig,c(TRUE,FALSE)))
    n_peaks = nrow(subdatByTissue)
    ## 
    rvByMotifByTissueByRegion = lapply( rownames(mytab),  function(thisRegion){
      tibble(motif_name=thisMotif,
             tissue = thisTissue, 
             region = thisRegion,
             n_peaks = n_peaks,
             n_sig_in_region = mytab[thisRegion,"TRUE"],
             n_nosig_in_region = mytab[thisRegion, "FALSE"],
             n_sig_outside_region = sum(mytab[setdiff(rownames(mytab),thisRegion),"TRUE"]),
             n_nosig_outisde_region = sum(mytab[setdiff(rownames(mytab),thisRegion),"FALSE"]))
    }
    )
    rvByMotifByTissueByRegion = do.call(rbind, rvByMotifByTissueByRegion)
    
    return(rvByMotifByTissueByRegion)
  })
  rvByMotifByTissue = do.call(rbind, ans)
})

motifCountByRegion = do.call(rbind, rvByMotif)

stopifnot(all(rowSums(motifCountByRegion%>% dplyr::select(n_sig_in_region:n_nosig_outisde_region))  == motifCountByRegion$n_peaks))

motifCountByRegionHasSig =  motifCountByRegion %>% dplyr::filter(n_sig_in_region>0)

do_fisher = function(x, alternative = "two.sided") {
  fisher.test(matrix(x,ncol=2), alternative = alternative)$p.value
}

motifCountByRegionFin <- motifCountByRegion

motifCountByRegionFin$pvalue = apply(motifCountByRegionFin[, c("n_sig_in_region", "n_sig_outside_region","n_nosig_in_region","n_nosig_outisde_region")], 1 , do_fisher, alternative="greater")

motifCountByRegionFin = motifCountByRegionFin %>% dplyr::mutate(padj=p.adjust(pvalue, method="BH")) %>% arrange(padj)

motifCountByRegionFin%>% dplyr::filter(padj<0.05) %>% dplyr::select(tissue,region) %>% ftable()

motifCountByRegionFin %>% dplyr::mutate(is_sig=padj<0.05) %>% dplyr::select(region,is_sig) %>% table()

motifCountByRegionHasSig$pvalue = apply(motifCountByRegionHasSig[, c("n_sig_in_region", "n_sig_outside_region","n_nosig_in_region","n_nosig_outisde_region")], 1 , do_fisher, alternative="greater")

motifCountByRegionHasSig = motifCountByRegionHasSig %>% dplyr::mutate(padj=p.adjust(pvalue, method="BH")) %>% arrange(padj)

motifCountByRegionHasSig%>% dplyr::filter(padj<0.05) %>% dplyr::select(tissue,region) %>% ftable()

motifCountByRegionHasSig %>% dplyr::mutate(is_sig=padj<0.05) %>% dplyr::select(region,is_sig) %>% table()

write_tsv(motifCountByRegionHasSig,file = "homer_motif_region_enrichment_by_tissue.txt")

currsigtflist <- unique(motifCountByRegionHasSig[motifCountByRegionHasSig$padj < 0.05,]$motif_name)
#currsigtflist <- unique(motifCountByRegionFin[motifCountByRegionFin$padj < 0.1,]$motif_name)

motifCountByRegioncurrsig <- matrix(1L,nrow = length(currsigtflist),ncol = 80)
rownames(motifCountByRegioncurrsig) <- currsigtflist
colnames(motifCountByRegioncurrsig) <- paste(rep(unique(peakanno$custom_annotation),8),
                                             c(rep("t52-hippocampus",10),
                                               rep("t55-gastrocnemius",10),
                                               rep("t58-heart",10),
                                               rep("t59-kidney",10),
                                               rep("t66-lung",10),
                                               rep("t68-liver",10),
                                               rep("t69-brown-adipose",10),
                                               rep("t70-white-adipose",10)),
                                             sep = "_")
motifCountByRegioncurrsigmeta <- data.frame(row.names = colnames(motifCountByRegioncurrsig),
                                            "Tissue" = c(rep("HIPPOC",10),
                                                         rep("SKM-GN",10),
                                                         rep("HEART",10),
                                                         rep("KIDNEY",10),
                                                         rep("LUNG",10),
                                                         rep("LIVER",10),
                                                         rep("BAT",10),
                                                         rep("WAT-SC",10)),
                                            "Region" = rep(unique(peakanno$custom_annotation),8))
for(i in 1:dim(motifCountByRegioncurrsig)[1]){
  currsig <- rownames(motifCountByRegioncurrsig)[i]
  currsigdata <- motifCountByRegionFin[motifCountByRegionFin$motif_name %in% currsig,]
  for(j in 1:dim(motifCountByRegioncurrsig)[2]){
    #motifCountByTissuecurrsig[i,j] <- currsigdata[currsigdata$tissue %in% gsub(".*_","",colnames(motifCountByTissuecurrsig))[j] & currsigdata$region %in% gsub("_.*","",colnames(motifCountByTissuecurrsig))[j],"padj"]
    if(length(currsigdata[currsigdata$tissue %in% gsub(".*_","",colnames(motifCountByRegioncurrsig))[j] & currsigdata$region %in% gsub("_.*","",colnames(motifCountByRegioncurrsig))[j],]$pvalue) > 0){
      motifCountByRegioncurrsig[i,j] <- currsigdata[currsigdata$tissue %in% gsub(".*_","",colnames(motifCountByRegioncurrsig))[j] & currsigdata$region %in% gsub("_.*","",colnames(motifCountByRegioncurrsig))[j],]$pvalue
    }
    
  }
}
motifCountByRegioncurrsigmeta <- motifCountByRegioncurrsigmeta[order(motifCountByRegioncurrsigmeta$Region),]
motifCountByRegioncurrsig <- motifCountByRegioncurrsig[,rownames(motifCountByRegioncurrsigmeta)]
motifCountByRegioncurrsigmeta$Region <- gsub(" ",".",motifCountByRegioncurrsigmeta$Region)
motifCountByRegioncurrsigmeta$Region <- gsub(".\\(<5kb\\)","",motifCountByRegioncurrsigmeta$Region)
motifCountByRegioncurrsigmetatrim <- motifCountByRegioncurrsigmeta[(!motifCountByRegioncurrsigmeta$Region %in% "Overlaps Gene") & (!motifCountByRegioncurrsigmeta$Tissue %in% c("HIPPOC","WAT-SC")),]

motifCountByRegioncurrsigtrim <- motifCountByRegioncurrsig[,rownames(motifCountByRegioncurrsigmetatrim)]
motifCountByRegioncurrsigtrim <- motifCountByRegioncurrsigtrim[apply(motifCountByRegioncurrsigtrim,1,min) < 0.05,]

motifCountByRegioncurrsigmetatrim$Region <- gsub(" ",".",motifCountByRegioncurrsigmetatrim$Region)
motifCountByRegioncurrsigmetatrim[motifCountByRegioncurrsigmetatrim$Region %in% "Upstream.(<5kb)","Region"] <- "Upstream"
motifCountByRegioncurrsigmetatrim[motifCountByRegioncurrsigmetatrim$Region %in% "Downstream.(<5kb)","Region"] <- "Downstream"

logpcutoff <- 2.5
darsigcolumns <- colnames(-log10(motifCountByRegioncurrsigtrim[,rownames(motifCountByRegioncurrsigmetatrim)]))[apply(-log10(motifCountByRegioncurrsigtrim[,rownames(motifCountByRegioncurrsigmetatrim)]),2,max) > logpcutoff]
ourcolumn <- darsigcolumns[1]
oursort <- rownames(motifCountByRegioncurrsigtrim)[order(-log10(motifCountByRegioncurrsigtrim[,ourcolumn]),decreasing = TRUE)]
ourvals <- rownames(motifCountByRegioncurrsigtrim)[-log10(motifCountByRegioncurrsigtrim[,ourcolumn]) > logpcutoff]
if(length(ourvals) < logpcutoff){
  oursigset <- oursort[1:length(ourvals)]
} else {
  oursigset <- oursort[1:5]
}

for(i in 2:length(darsigcolumns)){
  ourcolumn <- darsigcolumns[i]
  oursort <- rownames(motifCountByRegioncurrsigtrim)[order(-log10(motifCountByRegioncurrsigtrim[,ourcolumn]),decreasing = TRUE)]
  ourvals <- rownames(motifCountByRegioncurrsigtrim)[-log10(motifCountByRegioncurrsigtrim[,ourcolumn]) > logpcutoff]
  if(length(ourvals) < logpcutoff){
    oursigset <- c(oursigset,oursort[1:length(ourvals)])
  } else {
    oursigset <- c(oursigset,oursort[1:5])
  } 
}
oursigset <- unique(oursigset)


pdf("Figure 6D_021325.pdf",width = 7,height = 7)
pheatmap(-log10(motifCountByRegioncurrsigtrim[oursigset,rownames(motifCountByRegioncurrsigmetatrim[!motifCountByRegioncurrsigmetatrim$Region %in% "Overlaps.Gene",])]),labels_row = gsub("\\/.*","",gsub("\\(.*","",oursigset)),show_colnames = F,annotation_col = motifCountByRegioncurrsigmetatrim[,c("Region","Tissue")],breaks = seq(0,6,length.out = 101),color = colorpanel(101,"white","firebrick"),cluster_cols = F,annotation_colors = ann_cols,border_color = NA)
dev.off()

png("Figure 6D_021325.png",width = 7,height = 7,units = "in",res = 600)
pheatmap(-log10(motifCountByRegioncurrsigtrim[oursigset,rownames(motifCountByRegioncurrsigmetatrim[!motifCountByRegioncurrsigmetatrim$Region %in% "Overlaps.Gene",])]),labels_row = gsub("\\/.*","",gsub("\\(.*","",oursigset)),show_colnames = F,annotation_col = motifCountByRegioncurrsigmetatrim[,c("Region","Tissue")],breaks = seq(0,6,length.out = 101),color = colorpanel(101,"white","firebrick"),cluster_cols = F,annotation_colors = ann_cols,border_color = NA)
dev.off()


#####
# Figure 6F
####

load("atacnormmatrices.RData")
load("activepeakfiles.RData")
load("omesigdata.RData")

atactraining$isdegap <- 0

gastroatactraining <- atactraining[atactraining$tissue_abbreviation %in% "SKM-GN",]
heartatactraining <- atactraining[atactraining$tissue_abbreviation %in% "HEART",]
hippoatactraining <- atactraining[atactraining$tissue_abbreviation %in% "HIPPOC",]
kidneyatactraining <- atactraining[atactraining$tissue_abbreviation %in% "KIDNEY",]
liveratactraining <- atactraining[atactraining$tissue_abbreviation %in% "LIVER",]
lungatactraining <- atactraining[atactraining$tissue_abbreviation %in% "LUNG",]
brownatactraining <- atactraining[atactraining$tissue_abbreviation %in% "BAT",]
whiteatactraining <- atactraining[atactraining$tissue_abbreviation %in% "WAT-SC",]

gastroatactraining$isdegap <- (gastroatactraining$feature_ID %in% gastroactivepeaks) * (peakanno[gastroatactraining$feature_ID,"ensembl_gene"] %in% gastrornasig)
heartatactraining$isdegap <- (heartatactraining$feature_ID %in% heartactivepeaks) * (peakanno[heartatactraining$feature_ID,"ensembl_gene"] %in% heartrnasig)
hippoatactraining$isdegap <- (hippoatactraining$feature_ID %in% hippoactivepeaks) * (peakanno[hippoatactraining$feature_ID,"ensembl_gene"] %in% hippornasig)
kidneyatactraining$isdegap <- (kidneyatactraining$feature_ID %in% kidneyactivepeaks) * (peakanno[kidneyatactraining$feature_ID,"ensembl_gene"] %in% kidneyrnasig)
liveratactraining$isdegap <- (liveratactraining$feature_ID %in% liveractivepeaks) * (peakanno[liveratactraining$feature_ID,"ensembl_gene"] %in% liverrnasig)
lungatactraining$isdegap <- (lungatactraining$feature_ID %in% lungactivepeaks) * (peakanno[lungatactraining$feature_ID,"ensembl_gene"] %in% lungrnasig)
brownatactraining$isdegap <- (brownatactraining$feature_ID %in% brownactivepeaks) * (peakanno[brownatactraining$feature_ID,"ensembl_gene"] %in% brownrnasig)
whiteatactraining$isdegap <- (whiteatactraining$feature_ID %in% whiteactivepeaks) * (peakanno[whiteatactraining$feature_ID,"ensembl_gene"] %in% whiternasig)

atactraining <- rbind(gastroatactraining,heartatactraining,hippoatactraining,kidneyatactraining,
                      liveratactraining,lungatactraining,brownatactraining,whiteatactraining)
gc()

homerPositionAnnotation = allpeakmotifs %>% dplyr::mutate(motif_name=`Motif.Name`) %>% dplyr::select(PositionID, motif_name)
homerPositionAnnotation = unique(homerPositionAnnotation)
atacPeaksAnnotated = atactraining %>% dplyr::left_join(homerPositionAnnotation, by=c("feature_ID"="PositionID"))
pcutoff = 0.05
#atacPeaksAnnotatedWithMotif = atactraining[atactraining$feature_ID %in% allpeakmotifs$PositionID,]
atacPeaksAnnotatedWithMotif = atacPeaksAnnotated %>% dplyr::filter(!is.na(motif_name)) %>% dplyr::mutate(is_sig=isdegap>0)

gc()

rvByMotif = lapply(unique(atacPeaksAnnotatedWithMotif$motif_name), function(thisMotif){
  ## go by motif
  message(thisMotif)
  subdat = atacPeaksAnnotatedWithMotif %>% dplyr::filter(motif_name==thisMotif)
  
  ## by tissue
  ans = lapply(unique(subdat$tissue), function(thisTissue){
    message(thisTissue)
    
    subdatByTissue = subdat %>% dplyr::filter(tissue==thisTissue)
    mytab = table(subdatByTissue$custom_annotation, factor(subdatByTissue$is_sig,c(TRUE,FALSE)))
    n_peaks = nrow(subdatByTissue)
    ## 
    rvByMotifByTissueByRegion = lapply( rownames(mytab),  function(thisRegion){
      tibble(motif_name=thisMotif,
             tissue = thisTissue, 
             region = thisRegion,
             n_peaks = n_peaks,
             n_sig_in_region = mytab[thisRegion,"TRUE"],
             n_nosig_in_region = mytab[thisRegion, "FALSE"],
             n_sig_outside_region = sum(mytab[setdiff(rownames(mytab),thisRegion),"TRUE"]),
             n_nosig_outisde_region = sum(mytab[setdiff(rownames(mytab),thisRegion),"FALSE"]))
    }
    )
    rvByMotifByTissueByRegion = do.call(rbind, rvByMotifByTissueByRegion)
    
    return(rvByMotifByTissueByRegion)
  })
  rvByMotifByTissue = do.call(rbind, ans)
})

motifCountByRegion = do.call(rbind, rvByMotif)

stopifnot(all(rowSums(motifCountByRegion%>% dplyr::select(n_sig_in_region:n_nosig_outisde_region))  == motifCountByRegion$n_peaks))

motifCountByRegionHasSig =  motifCountByRegion %>% dplyr::filter(n_sig_in_region>0)

do_fisher = function(x, alternative = alternative) {
  fisher.test(matrix(unlist(x),ncol=2), alternative = alternative)$p.value
}

motifCountByRegionFin <- motifCountByRegion

motifCountByRegionFin$pvalue = apply(motifCountByRegionFin[, c("n_sig_in_region", "n_sig_outside_region","n_nosig_in_region","n_nosig_outisde_region")], 1 , do_fisher, alternative="greater")

motifCountByRegionFin = motifCountByRegionFin %>% dplyr::mutate(padj=p.adjust(pvalue, method="BH")) %>% arrange(padj)

motifCountByRegionFin%>% dplyr::filter(padj<0.05) %>% dplyr::select(tissue,region) %>% ftable()

motifCountByRegionFin %>% dplyr::mutate(is_sig=padj<0.05) %>% dplyr::select(region,is_sig) %>% table()

motifCountByRegionHasSig$pvalue = apply(motifCountByRegionHasSig[, c("n_sig_in_region", "n_sig_outside_region","n_nosig_in_region","n_nosig_outisde_region")], 1 , do_fisher, alternative="greater")

motifCountByRegionHasSig = motifCountByRegionHasSig %>% dplyr::mutate(padj=p.adjust(pvalue, method="BH")) %>% arrange(padj)

motifCountByRegionHasSig%>% dplyr::filter(padj<0.05) %>% dplyr::select(tissue,region) %>% ftable()

motifCountByRegionHasSig %>% dplyr::mutate(is_sig=padj<0.05) %>% dplyr::select(region,is_sig) %>% table()

write_tsv(motifCountByRegionHasSig,file = "homer_degap_motif_region_enrichment_by_tissue.txt")

currsigtflist <- unique(motifCountByRegionHasSig[motifCountByRegionHasSig$padj < 0.05,]$motif_name)
#currsigtflist <- unique(motifCountByRegionFin[motifCountByRegionFin$padj < 0.1,]$motif_name)

motifCountByRegioncurrsig <- matrix(1L,nrow = length(currsigtflist),ncol = 80)
rownames(motifCountByRegioncurrsig) <- currsigtflist
colnames(motifCountByRegioncurrsig) <- paste(rep(unique(peakanno$custom_annotation),8),
                                             c(rep("t52-hippocampus",10),
                                               rep("t55-gastrocnemius",10),
                                               rep("t58-heart",10),
                                               rep("t59-kidney",10),
                                               rep("t66-lung",10),
                                               rep("t68-liver",10),
                                               rep("t69-brown-adipose",10),
                                               rep("t70-white-adipose",10)),
                                             sep = "_")
motifCountByRegioncurrsigmeta <- data.frame(row.names = colnames(motifCountByRegioncurrsig),
                                            "Tissue" = c(rep("HIPPOC",10),
                                                         rep("SKM-GN",10),
                                                         rep("HEART",10),
                                                         rep("KIDNEY",10),
                                                         rep("LUNG",10),
                                                         rep("LIVER",10),
                                                         rep("BAT",10),
                                                         rep("WAT-SC",10)),
                                            "Region" = rep(unique(peakanno$custom_annotation),8))
for(i in 1:dim(motifCountByRegioncurrsig)[1]){
  currsig <- rownames(motifCountByRegioncurrsig)[i]
  currsigdata <- motifCountByRegionFin[motifCountByRegionFin$motif_name %in% currsig,]
  for(j in 1:dim(motifCountByRegioncurrsig)[2]){
    #motifCountByTissuecurrsig[i,j] <- currsigdata[currsigdata$tissue %in% gsub(".*_","",colnames(motifCountByTissuecurrsig))[j] & currsigdata$region %in% gsub("_.*","",colnames(motifCountByTissuecurrsig))[j],"padj"]
    if(length(currsigdata[currsigdata$tissue %in% gsub(".*_","",colnames(motifCountByRegioncurrsig))[j] & currsigdata$region %in% gsub("_.*","",colnames(motifCountByRegioncurrsig))[j],]$pvalue) > 0){
      motifCountByRegioncurrsig[i,j] <- currsigdata[currsigdata$tissue %in% gsub(".*_","",colnames(motifCountByRegioncurrsig))[j] & currsigdata$region %in% gsub("_.*","",colnames(motifCountByRegioncurrsig))[j],]$pvalue
    }
    
  }
}
motifCountByRegioncurrsigmeta <- motifCountByRegioncurrsigmeta[order(motifCountByRegioncurrsigmeta$Region),]
motifCountByRegioncurrsig <- motifCountByRegioncurrsig[,rownames(motifCountByRegioncurrsigmeta)]
motifCountByRegioncurrsigmeta$Region <- gsub(" ",".",motifCountByRegioncurrsigmeta$Region)
motifCountByRegioncurrsigmeta$Region <- gsub(".\\(<5kb\\)","",motifCountByRegioncurrsigmeta$Region)
motifCountByRegioncurrsigmetatrim <- motifCountByRegioncurrsigmeta[(!motifCountByRegioncurrsigmeta$Region %in% "Overlaps.Gene"),]

motifCountByRegioncurrsigtrim <- motifCountByRegioncurrsig[,rownames(motifCountByRegioncurrsigmetatrim)]
motifCountByRegioncurrsigtrim <- motifCountByRegioncurrsigtrim[apply(motifCountByRegioncurrsigtrim,1,min) < 0.01,]

degapsigcolumns <- colnames(-log10(motifCountByRegioncurrsigtrim[,rownames(motifCountByRegioncurrsigmetatrim)]))[apply(-log10(motifCountByRegioncurrsigtrim[,rownames(motifCountByRegioncurrsigmetatrim)]),2,max) > 5]
ourcolumn <- degapsigcolumns[1]
oursort <- rownames(motifCountByRegioncurrsigtrim)[order(-log10(motifCountByRegioncurrsigtrim[,ourcolumn]),decreasing = TRUE)]
ourvals <- rownames(motifCountByRegioncurrsigtrim)[-log10(motifCountByRegioncurrsigtrim[,ourcolumn]) > 5]
if(length(ourvals) < 5){
  oursigset <- oursort[1:length(ourvals)]
} else {
  oursigset <- oursort[1:5]
}

for(i in 2:length(degapsigcolumns)){
  ourcolumn <- degapsigcolumns[i]
  oursort <- rownames(motifCountByRegioncurrsigtrim)[order(-log10(motifCountByRegioncurrsigtrim[,ourcolumn]),decreasing = TRUE)]
  ourvals <- rownames(motifCountByRegioncurrsigtrim)[-log10(motifCountByRegioncurrsigtrim[,ourcolumn]) > 5]
  if(length(ourvals) < 5){
    oursigset <- c(oursigset,oursort[1:length(ourvals)])
  } else {
    oursigset <- c(oursigset,oursort[1:5])
  } 
}
oursigset <- unique(oursigset)

pdf("Figure 6F_021325.pdf",width = 7,height = 7)
pheatmap(-log10(motifCountByRegioncurrsigtrim[oursigset,rownames(motifCountByRegioncurrsigmetatrim)]),labels_row = gsub("\\/.*","",gsub("\\(.*","",oursigset)),show_colnames = F,annotation_col = motifCountByRegioncurrsigmetatrim[,c("Region","Tissue")],breaks = seq(0,6,length.out = 101),color = colorpanel(101,"white","firebrick"),cluster_cols = F,annotation_colors = ann_cols,border_color = NA)
dev.off()

png("Figure 6F_021325.png",width = 7,height = 7,units = "in",res = 600)
pheatmap(-log10(motifCountByRegioncurrsigtrim[oursigset,rownames(motifCountByRegioncurrsigmetatrim)]),labels_row = gsub("\\/.*","",gsub("\\(.*","",oursigset)),show_colnames = F,annotation_col = motifCountByRegioncurrsigmetatrim[,c("Region","Tissue")],breaks = seq(0,10,length.out = 101),color = colorpanel(101,"white","firebrick"),cluster_cols = F,annotation_colors = ann_cols,border_color = NA)
dev.off()

rm(atactraining)
rm(gastroatactraining)
rm(heartatactraining)
rm(hippoatactraining)
rm(kidneyatactraining)
rm(liveratactraining)
rm(lungatactraining)
rm(brownatactraining)
rm(whiteatactraining)
rm(atacPeaksAnnotated)
rm(atacPeaksAnnotatedWithMotif)
gc()

#####
# Figure 6E
####

load("PASS1B Transcription Factor Paper Data/DEA Analysis/pass1b-06_epigen-rrbs_dea.RData")

methtraining <- epigen_rrbs$training_dea
rm(epigen_rrbs)
gc()

#reloading this
trimmedmethanno <- readRDS("trimmedmethanno.RDS")

methtraining$custom_annotation <- ""
methtraining$custom_annotation <- trimmedmethanno[methtraining$feature_ID,"custom_annotation"]

allmethmotifs <- read.table(file = "allmethoutput50.txt",header = T,sep = "\t")

homerPositionAnnotation = allmethmotifs %>% dplyr::mutate(motif_name=`Motif.Name`) %>% dplyr::select(PositionID, motif_name)
homerPositionAnnotation = unique(homerPositionAnnotation)
methSitesAnnotated = methtraining %>% dplyr::left_join(homerPositionAnnotation, by=c("feature_ID"="PositionID"))
pcutoff = 0.1
#atacPeaksAnnotatedWithMotif = atactraining[atactraining$feature_ID %in% allpeakmotifs$PositionID,]
methSitesAnnotatedWithMotif = methSitesAnnotated %>% dplyr::filter(!is.na(motif_name)) %>% dplyr::mutate(is_sig=adj_p_value<pcutoff)

gc()

rvByMotif = lapply(unique(methSitesAnnotatedWithMotif$motif_name), function(thisMotif){
  ## go by motif
  message(thisMotif)
  subdat = methSitesAnnotatedWithMotif %>% dplyr::filter(motif_name==thisMotif)
  
  ## by tissue
  ans = lapply(unique(subdat$tissue), function(thisTissue){
    message(thisTissue)
    
    subdatByTissue = subdat %>% dplyr::filter(tissue==thisTissue)
    mytab = table(subdatByTissue$custom_annotation, factor(subdatByTissue$is_sig,c(TRUE,FALSE)))
    n_peaks = nrow(subdatByTissue)
    ## 
    rvByMotifByTissueByRegion = lapply( rownames(mytab),  function(thisRegion){
      tibble(motif_name=thisMotif,
             tissue = thisTissue, 
             region = thisRegion,
             n_peaks = n_peaks,
             n_sig_in_region = mytab[thisRegion,"TRUE"],
             n_nosig_in_region = mytab[thisRegion, "FALSE"],
             n_sig_outside_region = sum(mytab[setdiff(rownames(mytab),thisRegion),"TRUE"]),
             n_nosig_outisde_region = sum(mytab[setdiff(rownames(mytab),thisRegion),"FALSE"]))
    }
    )
    rvByMotifByTissueByRegion = do.call(rbind, rvByMotifByTissueByRegion)
    
    return(rvByMotifByTissueByRegion)
  })
  rvByMotifByTissue = do.call(rbind, ans)
})

motifCountByRegion = do.call(rbind, rvByMotif)

stopifnot(all(rowSums(motifCountByRegion%>% dplyr::select(n_sig_in_region:n_nosig_outisde_region))  == motifCountByRegion$n_peaks))

motifCountByRegionHasSig =  motifCountByRegion %>% dplyr::filter(n_sig_in_region>0)

do_fisher = function(x, alternative = "two.sided") {
  fisher.test(matrix(x,ncol=2), alternative = alternative)$p.value
}

motifCountByRegionFin <- motifCountByRegion

motifCountByRegionFin$pvalue = apply(motifCountByRegionFin[, c("n_sig_in_region", "n_sig_outside_region","n_nosig_in_region","n_nosig_outisde_region")], 1 , do_fisher, alternative="greater")

motifCountByRegionFin = motifCountByRegionFin %>% dplyr::mutate(padj=p.adjust(pvalue, method="BH")) %>% arrange(padj)

motifCountByRegionFin%>% dplyr::filter(padj<0.1) %>% dplyr::select(tissue,region) %>% ftable()

motifCountByRegionFin %>% dplyr::mutate(is_sig=padj<0.1) %>% dplyr::select(region,is_sig) %>% table()

motifCountByRegionHasSig$pvalue = apply(motifCountByRegionHasSig[, c("n_sig_in_region", "n_sig_outside_region","n_nosig_in_region","n_nosig_outisde_region")], 1 , do_fisher, alternative="greater")

motifCountByRegionHasSig = motifCountByRegionHasSig %>% dplyr::mutate(padj=p.adjust(pvalue, method="BH")) %>% arrange(padj)

motifCountByRegionHasSig%>% dplyr::filter(padj<0.1) %>% dplyr::select(tissue,region) %>% ftable()

motifCountByRegionHasSig %>% dplyr::mutate(is_sig=padj<0.1) %>% dplyr::select(region,is_sig) %>% table()

write_tsv(motifCountByRegionHasSig,file = "homer_methSites_motif_region_enrichment_by_tissue.txt")

currmethsigtflist <- unique(motifCountByRegionHasSig[motifCountByRegionHasSig$padj < 0.1,]$motif_name)
#currmethsigtflist <- unique(motifCountByRegionFin[motifCountByRegionFin$padj < 0.1,]$motif_name)

motifCountByRegioncurrsig <- matrix(1L,nrow = length(currmethsigtflist),ncol = 80)
rownames(motifCountByRegioncurrsig) <- currmethsigtflist
colnames(motifCountByRegioncurrsig) <- paste(rep(unique(peakanno$custom_annotation),8),
                                             c(rep("t52-hippocampus",10),
                                               rep("t55-gastrocnemius",10),
                                               rep("t58-heart",10),
                                               rep("t59-kidney",10),
                                               rep("t66-lung",10),
                                               rep("t68-liver",10),
                                               rep("t69-brown-adipose",10),
                                               rep("t70-white-adipose",10)),
                                             sep = "_")
motifCountByRegioncurrsigmeta <- data.frame(row.names = colnames(motifCountByRegioncurrsig),
                                            "Tissue" = c(rep("HIPPOC",10),
                                                         rep("SKM-GN",10),
                                                         rep("HEART",10),
                                                         rep("KIDNEY",10),
                                                         rep("LUNG",10),
                                                         rep("LIVER",10),
                                                         rep("BAT",10),
                                                         rep("WAT-SC",10)),
                                            "Region" = rep(unique(peakanno$custom_annotation),8))
for(i in 1:dim(motifCountByRegioncurrsig)[1]){
  currsig <- rownames(motifCountByRegioncurrsig)[i]
  currsigdata <- motifCountByRegionFin[motifCountByRegionFin$motif_name %in% currsig,]
  for(j in 1:dim(motifCountByRegioncurrsig)[2]){
    #motifCountByTissuecurrsig[i,j] <- currsigdata[currsigdata$tissue %in% gsub(".*_","",colnames(motifCountByTissuecurrsig))[j] & currsigdata$region %in% gsub("_.*","",colnames(motifCountByTissuecurrsig))[j],"padj"]
    if(length(currsigdata[currsigdata$tissue %in% gsub(".*_","",colnames(motifCountByRegioncurrsig))[j] & currsigdata$region %in% gsub("_.*","",colnames(motifCountByRegioncurrsig))[j],]$pvalue) > 0){
      motifCountByRegioncurrsig[i,j] <- currsigdata[currsigdata$tissue %in% gsub(".*_","",colnames(motifCountByRegioncurrsig))[j] & currsigdata$region %in% gsub("_.*","",colnames(motifCountByRegioncurrsig))[j],]$pvalue
    }
    
  }
}
motifCountByRegioncurrsigmeta <- motifCountByRegioncurrsigmeta[order(motifCountByRegioncurrsigmeta$Region),]
motifCountByRegioncurrsig <- motifCountByRegioncurrsig[,rownames(motifCountByRegioncurrsigmeta)]
motifCountByRegioncurrsigmeta$Region <- gsub(" ",".",motifCountByRegioncurrsigmeta$Region)
motifCountByRegioncurrsigmeta$Region <- gsub(".\\(<5kb\\)","",motifCountByRegioncurrsigmeta$Region)
motifCountByRegioncurrsigmetatrim <- motifCountByRegioncurrsigmeta[(!motifCountByRegioncurrsigmeta$Region %in% "Overlaps.Gene"),]

motifCountByRegioncurrsigtrim <- motifCountByRegioncurrsig[,rownames(motifCountByRegioncurrsigmetatrim)]
motifCountByRegioncurrsigtrim <- motifCountByRegioncurrsigtrim[apply(motifCountByRegioncurrsigtrim,1,min) < 0.001,]

logpcutoff <- 3.5
dmrsigcolumns <- colnames(-log10(motifCountByRegioncurrsigtrim[,rownames(motifCountByRegioncurrsigmetatrim)]))[apply(-log10(motifCountByRegioncurrsigtrim[,rownames(motifCountByRegioncurrsigmetatrim)]),2,max) > logpcutoff]
ourcolumn <- dmrsigcolumns[1]
oursort <- rownames(motifCountByRegioncurrsigtrim)[order(-log10(motifCountByRegioncurrsigtrim[,ourcolumn]),decreasing = TRUE)]
ourvals <- rownames(motifCountByRegioncurrsigtrim)[-log10(motifCountByRegioncurrsigtrim[,ourcolumn]) > logpcutoff]
if(length(ourvals) < logpcutoff){
  oursigset <- oursort[1:length(ourvals)]
} else {
  oursigset <- oursort[1:5]
}

for(i in 2:length(dmrsigcolumns)){
  ourcolumn <- dmrsigcolumns[i]
  oursort <- rownames(motifCountByRegioncurrsigtrim)[order(-log10(motifCountByRegioncurrsigtrim[,ourcolumn]),decreasing = TRUE)]
  ourvals <- rownames(motifCountByRegioncurrsigtrim)[-log10(motifCountByRegioncurrsigtrim[,ourcolumn]) > logpcutoff]
  if(length(ourvals) < logpcutoff){
    oursigset <- c(oursigset,oursort[1:length(ourvals)])
  } else {
    oursigset <- c(oursigset,oursort[1:5])
  } 
}
oursigset <- unique(oursigset)

pdf("Figure 6E_021325.pdf",width = 7,height = 7)
pheatmap(-log10(motifCountByRegioncurrsigtrim[oursigset,rownames(motifCountByRegioncurrsigmetatrim)]),labels_row = gsub("\\/.*","",gsub("\\(.*","",oursigset)),show_colnames = F,annotation_col = motifCountByRegioncurrsigmetatrim[,c("Region","Tissue")],breaks = seq(0,6,length.out = 101),color = colorpanel(101,"white","firebrick"),cluster_cols = F,annotation_colors = ann_cols,border_color = NA)
dev.off()

png("Figure 6E_021325.png",width = 7,height = 7,units = "in",res = 600)
pheatmap(-log10(motifCountByRegioncurrsigtrim[oursigset,rownames(motifCountByRegioncurrsigmetatrim)]),labels_row = gsub("\\/.*","",gsub("\\(.*","",oursigset)),show_colnames = F,annotation_col = motifCountByRegioncurrsigmetatrim[,c("Region","Tissue")],breaks = seq(0,6,length.out = 101),color = colorpanel(101,"white","firebrick"),cluster_cols = F,annotation_colors = ann_cols,border_color = NA)
dev.off()

rm(methSitesAnnotated)
rm(methSitesAnnotatedWithMotif)
gc()

#####
# Supplemental Figures - DARs
####

load("enstosym.RData")
load("rnacontrolnormandmeta.RData")

tfatacsigpctmat <- cbind(as.numeric(gsub("%","",gastroatacsigtf50[atacsigtflist,"X..of.Target.Sequences.with.Motif"])),
                         as.numeric(gsub("%","",heartatacsigtf50[atacsigtflist,"X..of.Target.Sequences.with.Motif"])),
                         as.numeric(gsub("%","",kidneyatacsigtf50[atacsigtflist,"X..of.Target.Sequences.with.Motif"])),
                         as.numeric(gsub("%","",liveratacsigtf50[atacsigtflist,"X..of.Target.Sequences.with.Motif"])),
                         as.numeric(gsub("%","",lungatacsigtf50[atacsigtflist,"X..of.Target.Sequences.with.Motif"])),
                         as.numeric(gsub("%","",brownatacsigtf50[atacsigtflist,"X..of.Target.Sequences.with.Motif"])))
rownames(tfatacsigpctmat) <- atacsigtflist
colnames(tfatacsigpctmat) <- c("SKM-GN","HEART","KIDNEY","LIVER","LUNG","BAT")

finalatacsigtfgenes <- intersect(tfanno[trimatacsigtoptfs,"Ensembl"],rownames(rnacontrolnorm))
finalatacsigtfs <- intersect(trimatacsigtoptfs,rownames(tfanno[tfanno$Ensembl %in% finalatacsigtfgenes,]))

pdf(file = "Supplemental Figure S27A_021325.pdf",width = 6,height = 6)
pheatmap(t(scale(t(tfatacsigpctmat[finalatacsigtfs,]))),cluster_rows = F,cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_row = gsub("\\(.*","",finalatacsigtfs),annotation_col = tissuemeta,annotation_colors = ann_cols,show_colnames = F,border_color = NA)
dev.off()

png(file = "Supplemental Figure S27A_021325.png",width = 6,height = 6,units = "in",res = 600)
pheatmap(t(scale(t(tfatacsigpctmat[finalatacsigtfs,]))),cluster_rows = F,cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_row = gsub("\\(.*","",finalatacsigtfs),annotation_col = tissuemeta,annotation_colors = ann_cols,show_colnames = F,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S31A_021325.pdf",width = 6,height = 6)
pheatmap(t(scale(t(rnacontrolnorm[tfanno[finalatacsigtfs,"Ensembl"],controlcolumns[c(1:20,31:70)]]))),cluster_rows = F,cluster_cols = F,breaks = seq(-3,3,length.out = 101),color = colorpanel(101,"blue","white","red"),annotation_col = rnacontrolmeta[,c("Sex","Tissue")],annotation_colors = ann_cols,show_colnames = F,labels_row = enstosym[tfanno[finalatacsigtfs,"Ensembl"],"Symbol"],cellwidth = 4.0,border_color = NA)
dev.off()

png(file = "Supplemental Figure S31A_021325.png",width = 6,height = 6,units = "in",res = 600)
pheatmap(t(scale(t(rnacontrolnorm[tfanno[finalatacsigtfs,"Ensembl"],controlcolumns[c(1:20,31:70)]]))),cluster_rows = F,cluster_cols = F,breaks = seq(-3,3,length.out = 101),color = colorpanel(101,"blue","white","red"),annotation_col = rnacontrolmeta[,c("Sex","Tissue")],annotation_colors = ann_cols,show_colnames = F,labels_row = enstosym[tfanno[finalatacsigtfs,"Ensembl"],"Symbol"],cellwidth = 4.0,border_color = NA)
dev.off()

load("rnal2fcmat.RData")

atacsigtfl2fcmat <- matrix(0L,nrow = length(finalatacsigtfs),ncol = 48)
rownames(atacsigtfl2fcmat) <- tfanno[finalatacsigtfs,"Ensembl"]
colnames(atacsigtfl2fcmat) <- c("SKM-GN_female_w1","SKM-GN_female_w2","SKM-GN_female_w4","SKM-GN_female_w8",
                                "SKM-GN_male_w1","SKM-GN_male_w2","SKM-GN_male_w4","SKM-GN_male_w8",
                                "HEART_female_w1","HEART_female_w2","HEART_female_w4","HEART_female_w8",
                                "HEART_male_w1","HEART_male_w2","HEART_male_w4","HEART_male_w8",
                                "KIDNEY_female_w1","KIDNEY_female_w2","KIDNEY_female_w4","KIDNEY_female_w8",
                                "KIDNEY_male_w1","KIDNEY_male_w2","KIDNEY_male_w4","KIDNEY_male_w8",
                                "LIVER_female_w1","LIVER_female_w2","LIVER_female_w4","LIVER_female_w8",
                                "LIVER_male_w1","LIVER_male_w2","LIVER_male_w4","LIVER_male_w8",
                                "LUNG_female_w1","LUNG_female_w2","LUNG_female_w4","LUNG_female_w8",
                                "LUNG_male_w1","LUNG_male_w2","LUNG_male_w4","LUNG_male_w8",
                                "BAT_female_w1","BAT_female_w2","BAT_female_w4","BAT_female_w8",
                                "BAT_male_w1","BAT_male_w2","BAT_male_w4","BAT_male_w8")

for(i in 1:length(finalatacsigtfs)){
  ourens <- tfanno[finalatacsigtfs[i],"Ensembl"]
  if(ourens %in% rownames(gastrol2fcmat)){
    atacsigtfl2fcmat[i,c(1:8)] <- gastrol2fcmat[ourens,]
  }
  if(ourens %in% rownames(heartl2fcmat)){
    atacsigtfl2fcmat[i,c(9:16)] <- heartl2fcmat[ourens,]
  }
  if(ourens %in% rownames(kidneyl2fcmat)){
    atacsigtfl2fcmat[i,c(17:24)] <- kidneyl2fcmat[ourens,]
  }
  if(ourens %in% rownames(liverl2fcmat)){
    atacsigtfl2fcmat[i,c(25:32)] <- liverl2fcmat[ourens,]
  }
  if(ourens %in% rownames(lungl2fcmat)){
    atacsigtfl2fcmat[i,c(33:40)] <- lungl2fcmat[ourens,]
  }
  if(ourens %in% rownames(brownl2fcmat)){
    atacsigtfl2fcmat[i,c(41:48)] <- brownl2fcmat[ourens,]
  }
  
}

atacsigtfl2fcmeta <- data.frame(row.names = colnames(atacsigtfl2fcmat),
                                "Tissue" = c(rep("SKM-GN",8),
                                             rep("HEART",8),
                                             rep("KIDNEY",8),
                                             rep("LIVER",8),
                                             rep("LUNG",8),
                                             rep("BAT",8)),
                                "Sex" = rep(c(rep("Female",4),rep("Male",4)),6),
                                "Group" = rep(c("1w","2w","4w","8w"),12))

pdf(file = "Supplemental Figure S32A_021325.pdf",width = 6.75,height = 6)
pheatmap(atacsigtfl2fcmat,cluster_rows = F,cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_row = enstosym[rownames(atacsigtfl2fcmat),"Symbol"],annotation_col = atacsigtfl2fcmeta[,c("Group","Sex","Tissue")],annotation_colors = ann_cols,cellwidth = 6,show_colnames = F,border_color = NA)
dev.off()

png(file = "Supplemental Figure S32A_021325.png",width = 6.75,height = 6,units = "in",res = 600)
pheatmap(atacsigtfl2fcmat,cluster_rows = F,cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_row = enstosym[rownames(atacsigtfl2fcmat),"Symbol"],annotation_col = atacsigtfl2fcmeta[,c("Group","Sex","Tissue")],annotation_colors = ann_cols,cellwidth = 6,show_colnames = F,border_color = NA)
dev.off()


tfatacsigpctmatz <- t(scale(t(tfatacsigpctmat)))

atactfcontrolmat <- rnacontrolnorm[tfanno[finalatacsigtfs,"Ensembl"],controlcolumns]
atactfcontrolmatz <- t(scale(t(atactfcontrolmat)))
atactfcontrolmeanmatz <- cbind(rowMeans(atactfcontrolmatz[,c(1:10)]),
                               rowMeans(atactfcontrolmatz[,c(11:20)]),
                               rowMeans(atactfcontrolmatz[,c(31:40)]),
                               rowMeans(atactfcontrolmatz[,c(41:50)]),
                               rowMeans(atactfcontrolmatz[,c(51:60)]),
                               rowMeans(atactfcontrolmatz[,c(61:70)]))
colnames(atactfcontrolmeanmatz) <- c("SKM-GN","HEART","KIDNEY","LIVER","LUNG","BAT")
atactfcontrolmeanmat <- cbind(rowMeans(atactfcontrolmatz[,c(1:10)]),
                              rowMeans(atactfcontrolmat[,c(11:20)]),
                              rowMeans(atactfcontrolmat[,c(31:40)]),
                              rowMeans(atactfcontrolmat[,c(41:50)]),
                              rowMeans(atactfcontrolmat[,c(51:60)]),
                              rowMeans(atactfcontrolmat[,c(61:70)]))
colnames(atactfcontrolmeanmat) <- c("SKM-GN","HEART","KIDNEY","LIVER","LUNG","BAT")

atactfexprenrichcor <- rbind(c(cor(tfatacsigpctmatz[finalatacsigtfs,"SKM-GN"],atactfcontrolmeanmatz[,"SKM-GN"]),
                               cor(tfatacsigpctmatz[finalatacsigtfs,"HEART"],atactfcontrolmeanmatz[,"HEART"]),
                               cor(tfatacsigpctmatz[finalatacsigtfs,"KIDNEY"],atactfcontrolmeanmatz[,"KIDNEY"]),
                               cor(tfatacsigpctmatz[finalatacsigtfs,"LIVER"],atactfcontrolmeanmatz[,"LIVER"]),
                               cor(tfatacsigpctmatz[finalatacsigtfs,"LUNG"],atactfcontrolmeanmatz[,"LUNG"]),
                               cor(tfatacsigpctmatz[finalatacsigtfs,"BAT"],atactfcontrolmeanmatz[,"BAT"])),
                             c(cor(tfatacsigpctmatz[finalatacsigtfs,"SKM-GN"],atactfcontrolmeanmatz[,"SKM-GN"]),
                               cor(tfatacsigpctmatz[finalatacsigtfs,"HEART"],atactfcontrolmeanmatz[,"HEART"]),
                               cor(tfatacsigpctmatz[finalatacsigtfs,"KIDNEY"],atactfcontrolmeanmatz[,"KIDNEY"]),
                               cor(tfatacsigpctmatz[finalatacsigtfs,"LIVER"],atactfcontrolmeanmatz[,"LIVER"]),
                               cor(tfatacsigpctmatz[finalatacsigtfs,"LUNG"],atactfcontrolmeanmatz[,"LUNG"]),
                               cor(tfatacsigpctmatz[finalatacsigtfs,"BAT"],atactfcontrolmeanmatz[,"BAT"])))
colnames(atactfexprenrichcor) <- c("SKM-GN",
                                   "HEART",
                                   "KIDNEY",
                                   "LIVER",
                                   "LUNG",
                                   "BAT")
rownames(atactfexprenrichcor) <- c("Correlation","Duplicate")

pdf(file = "Supplemental Figure S33A_021325.pdf",width = 3.5,height = 5.5)
pheatmap(atactfexprenrichcor["Correlation",],cluster_rows = F,cluster_cols = F,display_numbers = T,number_color = "black",breaks = seq(-0.6,0.6,length.out = 101),color = colorpanel(101,"royalblue4","white","firebrick"),labels_row = colnames(atactfexprenrichcor),angle_col = 0,fontsize = 20,show_rownames = T,cellheight = 60,border_color = NA)
dev.off()

png(file = "Supplemental Figure S33A_021325.png",width = 3.5,height = 5.5,units = "in",res = 600)
pheatmap(atactfexprenrichcor["Correlation",],cluster_rows = F,cluster_cols = F,display_numbers = T,number_color = "black",breaks = seq(-0.6,0.6,length.out = 101),color = colorpanel(101,"royalblue4","white","firebrick"),labels_row = colnames(atactfexprenrichcor),angle_col = 0,fontsize = 20,show_rownames = T,cellheight = 60,border_color = NA)
dev.off()

heartatacsigexprenrichcordf <- data.frame("TF" = gsub("\\(.*","",finalatacsigtfs),
                                          "Enrichment Z Score" = tfatacsigpctmatz[finalatacsigtfs,"HEART"],
                                          "Expression Z Score" = atactfcontrolmeanmatz[,"HEART"])

lungatacsigexprenrichcordf <- data.frame("TF" = gsub("\\(.*","",finalatacsigtfs),
                                         "Enrichment Z Score" = tfatacsigpctmatz[finalatacsigtfs,"LUNG"],
                                         "Expression Z Score" = atactfcontrolmeanmatz[,"LUNG"])

pdf(file = "Supplemental Figure S33B_021325.pdf",width = 6,height = 6)
ggplot(heartatacsigexprenrichcordf,aes(x=Enrichment.Z.Score,y=Expression.Z.Score,label=TF)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_text_repel(max.overlaps = 15)+geom_smooth(method=lm,  linetype="dashed",color="darkred", fill="blue") + theme_classic() + ggtitle(paste("HEART Correlation: ",toString(cor(tfatacsigpctmatz[finalatacsigtfs,"HEART"],atactfcontrolmeanmatz[,"HEART"])),sep = "")) + theme(text = element_text(size = 15))
dev.off()

png(file = "Supplemental Figure S33B_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(heartatacsigexprenrichcordf,aes(x=Enrichment.Z.Score,y=Expression.Z.Score,label=TF)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_text_repel(max.overlaps = 15)+geom_smooth(method=lm,  linetype="dashed",color="darkred", fill="blue") + theme_classic() + ggtitle(paste("HEART Correlation: ",toString(cor(tfatacsigpctmatz[finalatacsigtfs,"HEART"],atactfcontrolmeanmatz[,"HEART"])),sep = "")) + theme(text = element_text(size = 15))
dev.off()

pdf(file = "Supplemental Figure S33C_021325.pdf",width = 6,height = 6)
ggplot(lungatacsigexprenrichcordf,aes(x=Enrichment.Z.Score,y=Expression.Z.Score,label=TF)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_text_repel(max.overlaps = 15)+geom_smooth(method=lm,  linetype="dashed",color="darkred", fill="blue") + theme_classic() + ggtitle(paste("LUNG Correlation: ",toString(cor(tfatacsigpctmatz[finalatacsigtfs,"LUNG"],atactfcontrolmeanmatz[,"LUNG"])),sep = "")) + theme(text = element_text(size = 15))
dev.off()

png(file = "Supplemental Figure S33C_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(lungatacsigexprenrichcordf,aes(x=Enrichment.Z.Score,y=Expression.Z.Score,label=TF)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_text_repel(max.overlaps = 15)+geom_smooth(method=lm,  linetype="dashed",color="darkred", fill="blue") + theme_classic() + ggtitle(paste("LUNG Correlation: ",toString(cor(tfatacsigpctmatz[finalatacsigtfs,"LUNG"],atactfcontrolmeanmatz[,"LUNG"])),sep = "")) + theme(text = element_text(size = 15))
dev.off()

#####
# Supplemental Figures - DMRs
####

pdf(file = "Supplemental Figure S27B_021325.pdf",width = 6,height = 8)
pheatmap(t(scale(t(tfmethsigpctmat[finalmethsigtfs,]))),cluster_rows = F,cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_row = gsub("\\(.*","",finalmethsigtfs),annotation_col = tissuemeta,annotation_colors = ann_cols,show_colnames = F,border_color = NA)
dev.off()

png(file = "Supplemental Figure S27B_021325.png",width = 6,height = 8,units = "in",res = 600)
pheatmap(t(scale(t(tfmethsigpctmat[finalmethsigtfs,]))),cluster_rows = F,cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_row = gsub("\\(.*","",finalmethsigtfs),annotation_col = tissuemeta,annotation_colors = ann_cols,show_colnames = F,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S31B_021325.pdf",width = 6,height = 8)
pheatmap(t(scale(t(rnacontrolnorm[tfanno[finalmethsigtfs,"Ensembl"],controlcolumns[c(1:20,31:70)]]))),cluster_rows = F,cluster_cols = F,breaks = seq(-3,3,length.out = 101),color = colorpanel(101,"blue","white","red"),annotation_col = rnacontrolmeta[,c("Sex","Tissue")],annotation_colors = ann_cols,show_colnames = F,labels_row = enstosym[tfanno[finalmethsigtfs,"Ensembl"],"Symbol"],cellwidth = 4.0,border_color = NA)
dev.off()

png(file = "Supplemental Figure S31B_021325.png",width = 6,height = 8,units = "in",res = 600)
pheatmap(t(scale(t(rnacontrolnorm[tfanno[finalmethsigtfs,"Ensembl"],controlcolumns[c(1:20,31:70)]]))),cluster_rows = F,cluster_cols = F,breaks = seq(-3,3,length.out = 101),color = colorpanel(101,"blue","white","red"),annotation_col = rnacontrolmeta[,c("Sex","Tissue")],annotation_colors = ann_cols,show_colnames = F,labels_row = enstosym[tfanno[finalmethsigtfs,"Ensembl"],"Symbol"],cellwidth = 4.0,border_color = NA)
dev.off()

methsigtfl2fcmat <- matrix(0L,nrow = length(finalmethsigtfs),ncol = 64)
rownames(methsigtfl2fcmat) <- tfanno[finalmethsigtfs,"Ensembl"]
colnames(methsigtfl2fcmat) <- c("SKM-GN_female_w1","SKM-GN_female_w2","SKM-GN_female_w4","SKM-GN_female_w8",
                                "SKM-GN_male_w1","SKM-GN_male_w2","SKM-GN_male_w4","SKM-GN_male_w8",
                                "HEART_female_w1","HEART_female_w2","HEART_female_w4","HEART_female_w8",
                                "HEART_male_w1","HEART_male_w2","HEART_male_w4","HEART_male_w8",
                                "HIPPOC_female_w1","HIPPOC_female_w2","HIPPOC_female_w4","HIPPOC_female_w8",
                                "HIPPOC_male_w1","HIPPOC_male_w2","HIPPOC_male_w4","HIPPOC_male_w8",
                                "KIDNEY_female_w1","KIDNEY_female_w2","KIDNEY_female_w4","KIDNEY_female_w8",
                                "KIDNEY_male_w1","KIDNEY_male_w2","KIDNEY_male_w4","KIDNEY_male_w8",
                                "LIVER_female_w1","LIVER_female_w2","LIVER_female_w4","LIVER_female_w8",
                                "LIVER_male_w1","LIVER_male_w2","LIVER_male_w4","LIVER_male_w8",
                                "LUNG_female_w1","LUNG_female_w2","LUNG_female_w4","LUNG_female_w8",
                                "LUNG_male_w1","LUNG_male_w2","LUNG_male_w4","LUNG_male_w8",
                                "BAT_female_w1","BAT_female_w2","BAT_female_w4","BAT_female_w8",
                                "BAT_male_w1","BAT_male_w2","BAT_male_w4","BAT_male_w8",
                                "WAT-SC_female_w1","WAT-SC_female_w2","WAT-SC_female_w4","WAT-SC_female_w8",
                                "WAT-SC_male_w1","WAT-SC_male_w2","WAT-SC_male_w4","WAT-SC_male_w8")

for(i in 1:length(finalmethsigtfs)){
  ourens <- tfanno[finalmethsigtfs[i],"Ensembl"]
  if(ourens %in% rownames(gastrol2fcmat)){
    methsigtfl2fcmat[i,c(1:8)] <- gastrol2fcmat[ourens,]
  }
  if(ourens %in% rownames(heartl2fcmat)){
    methsigtfl2fcmat[i,c(9:16)] <- heartl2fcmat[ourens,]
  }
  if(ourens %in% rownames(hippol2fcmat)){
    methsigtfl2fcmat[i,c(17:24)] <- hippol2fcmat[ourens,]
  }
  if(ourens %in% rownames(kidneyl2fcmat)){
    methsigtfl2fcmat[i,c(25:32)] <- kidneyl2fcmat[ourens,]
  }
  if(ourens %in% rownames(liverl2fcmat)){
    methsigtfl2fcmat[i,c(33:40)] <- liverl2fcmat[ourens,]
  }
  if(ourens %in% rownames(lungl2fcmat)){
    methsigtfl2fcmat[i,c(41:48)] <- lungl2fcmat[ourens,]
  }
  if(ourens %in% rownames(brownl2fcmat)){
    methsigtfl2fcmat[i,c(49:56)] <- brownl2fcmat[ourens,]
  }
  if(ourens %in% rownames(whitel2fcmat)){
    methsigtfl2fcmat[i,c(57:64)] <- whitel2fcmat[ourens,]
  }
}

methsigtfl2fcmeta <- data.frame(row.names = colnames(methsigtfl2fcmat),
                                "Tissue" = c(rep("SKM-GN",8),
                                             rep("HEART",8),
                                             rep("HIPPOC",8),
                                             rep("KIDNEY",8),
                                             rep("LIVER",8),
                                             rep("LUNG",8),
                                             rep("BAT",8),
                                             rep("WAT-SC",8)),
                                "Sex" = rep(c(rep("Female",4),rep("Male",4)),8),
                                "Group" = rep(c("1w","2w","4w","8w"),16))

pdf(file = "Supplemental Figure S32B_021325.pdf",width = 6.75,height = 8)
pheatmap(methsigtfl2fcmat,cluster_rows = F,cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_row = enstosym[rownames(methsigtfl2fcmat),"Symbol"],annotation_col = methsigtfl2fcmeta[,c("Group","Sex","Tissue")],annotation_colors = ann_cols,cellwidth = 6,show_colnames = F,border_color = NA)
dev.off()
png(file = "Supplemental Figure S32B_021325.png",width = 6.75,height = 8,units = "in",res = 600)
pheatmap(methsigtfl2fcmat,cluster_rows = F,cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_row = enstosym[rownames(methsigtfl2fcmat),"Symbol"],annotation_col = methsigtfl2fcmeta[,c("Group","Sex","Tissue")],annotation_colors = ann_cols,cellwidth = 6,show_colnames = F,border_color = NA)
dev.off()

tfmethsigpctmatz <- t(scale(t(tfmethsigpctmat)))

methtfcontrolmat <- rnacontrolnorm[tfanno[finalmethsigtfs,"Ensembl"],controlcolumns]
methtfcontrolmatz <- t(scale(t(methtfcontrolmat)))
methtfcontrolmeanmatz <- cbind(rowMeans(methtfcontrolmatz[,c(1:10)]),
                               rowMeans(methtfcontrolmatz[,c(11:20)]),
                               rowMeans(methtfcontrolmatz[,c(21:30)]),
                               rowMeans(methtfcontrolmatz[,c(31:40)]),
                               rowMeans(methtfcontrolmatz[,c(41:50)]),
                               rowMeans(methtfcontrolmatz[,c(51:60)]),
                               rowMeans(methtfcontrolmatz[,c(61:70)]),
                               rowMeans(methtfcontrolmatz[,c(71:80)]))
colnames(methtfcontrolmeanmatz) <- c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC")
methtfcontrolmeanmat <- cbind(rowMeans(methtfcontrolmatz[,c(1:10)]),
                              rowMeans(methtfcontrolmat[,c(11:20)]),
                              rowMeans(methtfcontrolmat[,c(21:30)]),
                              rowMeans(methtfcontrolmat[,c(31:40)]),
                              rowMeans(methtfcontrolmat[,c(41:50)]),
                              rowMeans(methtfcontrolmat[,c(51:60)]),
                              rowMeans(methtfcontrolmat[,c(61:70)]),
                              rowMeans(methtfcontrolmat[,c(71:80)]))
colnames(methtfcontrolmeanmat) <- c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC")

methtfexprenrichcor <- rbind(c(cor(tfmethsigpctmatz[finalmethsigtfs,"SKM-GN"],methtfcontrolmeanmatz[,"SKM-GN"]),
                               cor(tfmethsigpctmatz[finalmethsigtfs,"HEART"],methtfcontrolmeanmatz[,"HEART"]),
                               cor(tfmethsigpctmatz[finalmethsigtfs,"HIPPOC"],methtfcontrolmeanmatz[,"HIPPOC"]),
                               cor(tfmethsigpctmatz[finalmethsigtfs,"KIDNEY"],methtfcontrolmeanmatz[,"KIDNEY"]),
                               cor(tfmethsigpctmatz[finalmethsigtfs,"LIVER"],methtfcontrolmeanmatz[,"LIVER"]),
                               cor(tfmethsigpctmatz[finalmethsigtfs,"LUNG"],methtfcontrolmeanmatz[,"LUNG"]),
                               cor(tfmethsigpctmatz[finalmethsigtfs,"BAT"],methtfcontrolmeanmatz[,"BAT"]),
                               cor(tfmethsigpctmatz[finalmethsigtfs,"WAT-SC"],methtfcontrolmeanmatz[,"WAT-SC"])),
                             c(cor(tfmethsigpctmatz[finalmethsigtfs,"SKM-GN"],methtfcontrolmeanmatz[,"SKM-GN"]),
                               cor(tfmethsigpctmatz[finalmethsigtfs,"HEART"],methtfcontrolmeanmatz[,"HEART"]),
                               cor(tfmethsigpctmatz[finalmethsigtfs,"HIPPOC"],methtfcontrolmeanmatz[,"HIPPOC"]),
                               cor(tfmethsigpctmatz[finalmethsigtfs,"KIDNEY"],methtfcontrolmeanmatz[,"KIDNEY"]),
                               cor(tfmethsigpctmatz[finalmethsigtfs,"LIVER"],methtfcontrolmeanmatz[,"LIVER"]),
                               cor(tfmethsigpctmatz[finalmethsigtfs,"LUNG"],methtfcontrolmeanmatz[,"LUNG"]),
                               cor(tfmethsigpctmatz[finalmethsigtfs,"BAT"],methtfcontrolmeanmatz[,"BAT"]),
                               cor(tfmethsigpctmatz[finalmethsigtfs,"WAT-SC"],methtfcontrolmeanmatz[,"WAT-SC"])))
colnames(methtfexprenrichcor) <- c("SKM-GN",
                                   "HEART",
                                   "HIPPOC",
                                   "KIDNEY",
                                   "LIVER",
                                   "LUNG",
                                   "BAT",
                                   "WAT-SC")
rownames(methtfexprenrichcor) <- c("Correlation","Duplicate")

pdf(file = "Supplemental Figure S33D_021325.pdf",width = 3.5,height = 5.5)
pheatmap(methtfexprenrichcor["Correlation",],cluster_rows = F,cluster_cols = F,display_numbers = T,number_color = "black",breaks = seq(-0.6,0.6,length.out = 101),color = colorpanel(101,"royalblue4","white","firebrick"),labels_row = colnames(methtfexprenrichcor),angle_col = 0,fontsize = 20,show_rownames = T,cellheight = 60,border_color = NA)
dev.off()

png(file = "Supplemental Figure S33D_021325.png",width = 3.5,height = 5.5,units = "in",res = 600)
pheatmap(methtfexprenrichcor["Correlation",],cluster_rows = F,cluster_cols = F,display_numbers = T,number_color = "black",breaks = seq(-0.6,0.6,length.out = 101),color = colorpanel(101,"royalblue4","white","firebrick"),labels_row = colnames(methtfexprenrichcor),angle_col = 0,fontsize = 20,show_rownames = T,cellheight = 43,border_color = NA)
dev.off()

hippomethsigexprenrichcordf <- data.frame("TF" = gsub("\\(.*","",finalmethsigtfs),
                                          "Enrichment Z Score" = tfmethsigpctmatz[finalmethsigtfs,"HIPPOC"],
                                          "Expression Z Score" = methtfcontrolmeanmatz[,"HIPPOC"])

brownmethsigexprenrichcordf <- data.frame("TF" = gsub("\\(.*","",finalmethsigtfs),
                                          "Enrichment Z Score" = tfmethsigpctmatz[finalmethsigtfs,"BAT"],
                                          "Expression Z Score" = methtfcontrolmeanmatz[,"BAT"])

pdf(file = "Supplemental Figure S33E_021325.pdf",width = 6,height = 6)
ggplot(hippomethsigexprenrichcordf,aes(x=Enrichment.Z.Score,y=Expression.Z.Score,label=TF)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_text_repel(max.overlaps = 15)+geom_smooth(method=lm,  linetype="dashed",color="darkred", fill="blue") + theme_classic() + ggtitle(paste("HIPPOC Correlation: ",toString(cor(tfmethsigpctmatz[finalmethsigtfs,"HIPPOC"],methtfcontrolmeanmatz[,"HIPPOC"])),sep = "")) + theme(text = element_text(size = 15))
dev.off()

png(file = "Supplemental Figure S33E_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(hippomethsigexprenrichcordf,aes(x=Enrichment.Z.Score,y=Expression.Z.Score,label=TF)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_text_repel(max.overlaps = 15)+geom_smooth(method=lm,  linetype="dashed",color="darkred", fill="blue") + theme_classic() + ggtitle(paste("HIPPOC Correlation: ",toString(cor(tfmethsigpctmatz[finalmethsigtfs,"HIPPOC"],methtfcontrolmeanmatz[,"HIPPOC"])),sep = "")) + theme(text = element_text(size = 15))
dev.off()

pdf(file = "Supplemental Figure S33F_021325.pdf",width = 6,height = 6)
ggplot(brownmethsigexprenrichcordf,aes(x=Enrichment.Z.Score,y=Expression.Z.Score,label=TF)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_text_repel(max.overlaps = 15)+geom_smooth(method=lm,  linetype="dashed",color="darkred", fill="blue") + theme_classic() + ggtitle(paste("BAT Correlation: ",toString(cor(tfmethsigpctmatz[finalmethsigtfs,"BAT"],methtfcontrolmeanmatz[,"BAT"])),sep = "")) + theme(text = element_text(size = 15))
dev.off()

png(file = "Supplemental Figure S33F_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(brownmethsigexprenrichcordf,aes(x=Enrichment.Z.Score,y=Expression.Z.Score,label=TF)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_text_repel(max.overlaps = 15)+geom_smooth(method=lm,  linetype="dashed",color="darkred", fill="blue") + theme_classic() + ggtitle(paste("BAT Correlation: ",toString(cor(tfmethsigpctmatz[finalmethsigtfs,"BAT"],methtfcontrolmeanmatz[,"BAT"])),sep = "")) + theme(text = element_text(size = 15))
dev.off()


#####
# Supplemental Figures - DEGaPs
####

tf50pctmat <- cbind(as.numeric(gsub("%","",gastrotf[tflist,"X..of.Target.Sequences.with.Motif"])),
                    as.numeric(gsub("%","",hearttf[tflist,"X..of.Target.Sequences.with.Motif"])),
                    as.numeric(gsub("%","",hippotf[tflist,"X..of.Target.Sequences.with.Motif"])),
                    as.numeric(gsub("%","",kidneytf[tflist,"X..of.Target.Sequences.with.Motif"])),
                    as.numeric(gsub("%","",livertf[tflist,"X..of.Target.Sequences.with.Motif"])),
                    as.numeric(gsub("%","",lungtf[tflist,"X..of.Target.Sequences.with.Motif"])),
                    as.numeric(gsub("%","",browntf[tflist,"X..of.Target.Sequences.with.Motif"])),
                    as.numeric(gsub("%","",whitetf[tflist,"X..of.Target.Sequences.with.Motif"])))
rownames(tf50pctmat) <- tflist
colnames(tf50pctmat) <- c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC")


# Supplemental Figure S27C
pdf(file = "Supplemental Figure S27C_021325.pdf",width = 6,height = 6)
pheatmap(t(scale(t(tf50pctmat[finalsigtflist,]))),labels_row = gsub("\\(.*","",finalsigtflist),angle_col = 0,color = colorpanel(101,"blue","white","red"),breaks = seq(-2,2,length.out = 101),cluster_cols = F,cluster_rows = F,annotation_col = tissuemeta,annotation_colors = ann_colsheatmaptrim,show_colnames = F,border_color = NA)
dev.off()

png(file = "Supplemental Figure S27C_021325.png",width = 6,height = 6,units = "in",res = 600)
pheatmap(t(scale(t(tf50pctmat[finalsigtflist,]))),labels_row = gsub("\\(.*","",finalsigtflist),angle_col = 0,color = colorpanel(101,"blue","white","red"),breaks = seq(-2,2,length.out = 101),cluster_cols = F,cluster_rows = F,annotation_col = tissuemeta,annotation_colors = ann_colsheatmaptrim,show_colnames = F,border_color = NA)
dev.off()

sigtfensembl <- tfanno[finalsigtflist,"Ensembl"]

ourcontrolmat <- rnacontrolnorm[sigtfensembl,controlcolumns]
ourcontrolmatz <- t(scale(t(ourcontrolmat)))
ourcontrolmeanmatz <- cbind(rowMeans(ourcontrolmatz[,c(1:10)]),
                            rowMeans(ourcontrolmatz[,c(11:20)]),
                            rowMeans(ourcontrolmatz[,c(21:30)]),
                            rowMeans(ourcontrolmatz[,c(31:40)]),
                            rowMeans(ourcontrolmatz[,c(41:50)]),
                            rowMeans(ourcontrolmatz[,c(51:60)]),
                            rowMeans(ourcontrolmatz[,c(61:70)]),
                            rowMeans(ourcontrolmatz[,c(71:80)]))
colnames(ourcontrolmeanmatz) <- c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC")
ourcontrolmeanmat <- cbind(rowMeans(ourcontrolmat[,c(1:10)]),
                           rowMeans(ourcontrolmat[,c(11:20)]),
                           rowMeans(ourcontrolmat[,c(21:30)]),
                           rowMeans(ourcontrolmat[,c(31:40)]),
                           rowMeans(ourcontrolmat[,c(41:50)]),
                           rowMeans(ourcontrolmat[,c(51:60)]),
                           rowMeans(ourcontrolmat[,c(61:70)]),
                           rowMeans(ourcontrolmat[,c(71:80)]))
colnames(ourcontrolmeanmat) <- c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC")

rnacontrolmetacolumns <- rnacontrolmeta[controlcolumns,]

tissue_cols <- list("Tissue" = c("SKM-GN" = "#088c03",
                                 "HEART" = "#f28b2f",
                                 "HIPPOC" = "#bf7534",
                                 "KIDNEY" = "#7553a7",
                                 "LIVER" = "#da6c75",
                                 "LUNG" = "#04bf8a",
                                 "BAT" = "#8c5220",
                                 "WAT-SC" = "#214da6"))

pdf(file = "Supplemental Figure S31C_021325.pdf",width = 6,height = 6)
pheatmap(t(scale(t(rnacontrolnorm[sigtfensembl,controlcolumns]))),cluster_rows = F,cluster_cols = F,labels_row = tfanno[finalsigtflist,"Gene.Name"],show_colnames = F,annotation_col = rnacontrolmetacolumns[,c("Sex","Tissue")],breaks = seq(-3,3,length.out = 101),color = colorpanel(101,"blue","white","red"),annotation_colors = ann_cols,cellwidth = 3.5,border_color = NA)
dev.off()

png(file = "Supplemental Figure S31C_021325.png",width = 6,height = 6,units = "in",res = 600)
pheatmap(t(scale(t(rnacontrolnorm[sigtfensembl,controlcolumns]))),cluster_rows = F,cluster_cols = F,labels_row = tfanno[finalsigtflist,"Gene.Name"],show_colnames = F,annotation_col = rnacontrolmetacolumns[,c("Sex","Tissue")],breaks = seq(-3,3,length.out = 101),color = colorpanel(101,"blue","white","red"),annotation_colors = ann_cols,cellwidth = 3.5,border_color = NA)
dev.off()

rnasigtfl2fcmat <- matrix(0L,nrow = length(sigtfensembl),ncol = 64)
rownames(rnasigtfl2fcmat) <- sigtfensembl
colnames(rnasigtfl2fcmat) <- c("SKM-GN F 1W","SKM-GN F 2W","SKM-GN F 4W","SKM-GN F 8W",
                               "SKM-GN M 1W","SKM-GN M 2W","SKM-GN M 4W","SKM-GN M 8W",
                               "HEART F 1W","HEART F 2W","HEART F 4W","HEART F 8W",
                               "HEART M 1W","HEART M 2W","HEART M 4W","HEART M 8W",
                               "HIPPOC F 1W","HIPPOC F 2W","HIPPOC F 4W","HIPPOC F 8W",
                               "HIPPOC M 1W","HIPPOC M 2W","HIPPOC M 4W","HIPPOC M 8W",
                               "KIDNEY F 1W","KIDNEY F 2W","KIDNEY F 4W","KIDNEY F 8W",
                               "KIDNEY M 1W","KIDNEY M 2W","KIDNEY M 4W","KIDNEY M 8W",
                               "LIVER F 1W","LIVER F 2W","LIVER F 4W","LIVER F 8W",
                               "LIVER M 1W","LIVER M 2W","LIVER M 4W","LIVER M 8W",
                               "LUNG F 1W","LUNG F 2W","LUNG F 4W","LUNG F 8W",
                               "LUNG M 1W","LUNG M 2W","LUNG M 4W","LUNG M 8W",
                               "BAT F 1W","BAT F 2W","BAT F 4W","BAT F 8W",
                               "BAT M 1W","BAT M 2W","BAT M 4W","BAT M 8W",
                               "WAT-SC F 1W","WAT-SC F 2W","WAT-SC F 4W","WAT-SC F 8W",
                               "WAT-SC M 1W","WAT-SC M 2W","WAT-SC M 4W","WAT-SC M 8W")

for(i in 1:length(sigtfensembl)){
  ourens <- sigtfensembl[i]
  if(ourens %in% rownames(gastrol2fcmat)){
    rnasigtfl2fcmat[i,c(1:8)] <- gastrol2fcmat[ourens,]
  }
  if(ourens %in% rownames(heartl2fcmat)){
    rnasigtfl2fcmat[i,c(9:16)] <- heartl2fcmat[ourens,]
  }
  if(ourens %in% rownames(hippol2fcmat)){
    rnasigtfl2fcmat[i,c(17:24)] <- hippol2fcmat[ourens,]
  }
  if(ourens %in% rownames(kidneyl2fcmat)){
    rnasigtfl2fcmat[i,c(25:32)] <- kidneyl2fcmat[ourens,]
  }
  if(ourens %in% rownames(liverl2fcmat)){
    rnasigtfl2fcmat[i,c(33:40)] <- liverl2fcmat[ourens,]
  }
  if(ourens %in% rownames(lungl2fcmat)){
    rnasigtfl2fcmat[i,c(41:48)] <- lungl2fcmat[ourens,]
  }
  if(ourens %in% rownames(brownl2fcmat)){
    rnasigtfl2fcmat[i,c(49:56)] <- brownl2fcmat[ourens,]
  }
  if(ourens %in% rownames(whitel2fcmat)){
    rnasigtfl2fcmat[i,c(57:64)] <- whitel2fcmat[ourens,]
  }
}

columnmetadf <- data.frame(row.names = colnames(rnasigtfl2fcmat),
                           "Tissue" = c(rep("SKM-GN",8),
                                        rep("HEART",8),
                                        rep("HIPPOC",8),
                                        rep("KIDNEY",8),
                                        rep("LIVER",8),
                                        rep("LUNG",8),
                                        rep("BAT",8),
                                        rep("WAT-SC",8)),
                           "Sex" = rep(c(rep("Female",4),rep("Male",4)),8),
                           "Group" = rep(c("1w","2w","4w","8w"),16))

pdf(file = "Supplemental Figure S32C_021325.pdf",width = 7,height = 6)
pheatmap(rnasigtfl2fcmat,cluster_rows = F,angle_col = 315,labels_row = tfanno[finalsigtflist,"Gene.Name"],breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,annotation_col = columnmetadf[,c("Group","Sex","Tissue")],show_colnames = F,annotation_colors = ann_cols,cellwidth = 5,border_color = NA)
dev.off()

png(file = "Supplemental Figure S32C_021325.png",width = 7,height = 6,units = "in",res = 600)
pheatmap(rnasigtfl2fcmat,cluster_rows = F,angle_col = 315,labels_row = tfanno[finalsigtflist,"Gene.Name"],breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,annotation_col = columnmetadf[,c("Group","Sex","Tissue")],show_colnames = F,annotation_colors = ann_cols,cellwidth = 5,border_color = NA)
dev.off()


tf50pctmatz <- t(scale(t(tf50pctmat)))

tfexprenrichcor <- rbind(c(cor(tf50pctmatz[finalsigtflist,"SKM-GN"],ourcontrolmeanmatz[,"SKM-GN"]),
                           cor(tf50pctmatz[finalsigtflist,"HEART"],ourcontrolmeanmatz[,"HEART"]),
                           cor(tf50pctmatz[finalsigtflist,"HIPPOC"],ourcontrolmeanmatz[,"HIPPOC"]),
                           cor(tf50pctmatz[finalsigtflist,"KIDNEY"],ourcontrolmeanmatz[,"KIDNEY"]),
                           cor(tf50pctmatz[finalsigtflist,"LIVER"],ourcontrolmeanmatz[,"LIVER"]),
                           cor(tf50pctmatz[finalsigtflist,"LUNG"],ourcontrolmeanmatz[,"LUNG"]),
                           cor(tf50pctmatz[finalsigtflist,"BAT"],ourcontrolmeanmatz[,"BAT"]),
                           cor(tf50pctmatz[finalsigtflist,"WAT-SC"],ourcontrolmeanmatz[,"WAT-SC"])),
                         c(cor(tf50pctmatz[finalsigtflist,"SKM-GN"],ourcontrolmeanmatz[,"SKM-GN"]),
                           cor(tf50pctmatz[finalsigtflist,"HEART"],ourcontrolmeanmatz[,"HEART"]),
                           cor(tf50pctmatz[finalsigtflist,"HIPPOC"],ourcontrolmeanmatz[,"HIPPOC"]),
                           cor(tf50pctmatz[finalsigtflist,"KIDNEY"],ourcontrolmeanmatz[,"KIDNEY"]),
                           cor(tf50pctmatz[finalsigtflist,"LIVER"],ourcontrolmeanmatz[,"LIVER"]),
                           cor(tf50pctmatz[finalsigtflist,"LUNG"],ourcontrolmeanmatz[,"LUNG"]),
                           cor(tf50pctmatz[finalsigtflist,"BAT"],ourcontrolmeanmatz[,"BAT"]),
                           cor(tf50pctmatz[finalsigtflist,"WAT-SC"],ourcontrolmeanmatz[,"WAT-SC"])))
colnames(tfexprenrichcor) <- c("SKM-GN",
                               "HEART",
                               "HIPPOC",
                               "KIDNEY",
                               "LIVER",
                               "LUNG",
                               "BAT",
                               "WAT-SC")
rownames(tfexprenrichcor) <- c("Correlation","Duplicate")

# Supplemental Figure S33G
pdf(file = "Supplemental Figure S33G_021325.pdf",width = 3.5,height = 5.5)
pheatmap(tfexprenrichcor["Correlation",],cluster_rows = F,cluster_cols = F,display_numbers = T,number_color = "black",breaks = seq(0,0.6,length.out = 101),color = colorpanel(101,"white","firebrick"),labels_row = colnames(tfexprenrichcor),angle_col = 0,fontsize = 20,show_rownames = T,cellheight = 43,border_color = NA)
dev.off()

png(file = "Supplemental Figure S33G_021325.png",width = 3.5,height = 5.5,units = "in",res = 600)
pheatmap(tfexprenrichcor["Correlation",],cluster_rows = F,cluster_cols = F,display_numbers = T,number_color = "black",breaks = seq(0,0.6,length.out = 101),color = colorpanel(101,"white","firebrick"),labels_row = colnames(tfexprenrichcor),angle_col = 0,fontsize = 20,show_rownames = T,cellheight = 43,border_color = NA)
dev.off()


liverexprenrichcordf <- data.frame(row.names = gsub("\\(.*","",finalsigtflist),
                                   "TF" = gsub("\\(.*","",finalsigtflist),
                                   "Enrichment Z Score" = tf50pctmatz[finalsigtflist,"LIVER"],
                                   "Expression Z Score" = ourcontrolmeanmatz[,"LIVER"])
lungexprenrichcordf <- data.frame(row.names = gsub("\\(.*","",finalsigtflist),
                                  "TF" = gsub("\\(.*","",finalsigtflist),
                                  "Enrichment Z Score" = tf50pctmatz[finalsigtflist,"LUNG"],
                                  "Expression Z Score" = ourcontrolmeanmatz[,"LUNG"])

# Supplemental Figure S33H
pdf(file = "Supplemental Figure S33H_021325.pdf",width = 6,height = 6)
ggplot(liverexprenrichcordf,aes(x=Enrichment.Z.Score,y=Expression.Z.Score,label=TF)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_text_repel(max.overlaps = 15)+geom_smooth(method=lm,  linetype="dashed",color="darkred", fill="blue") + theme_classic() + ggtitle(paste("LIVER Correlation: ",toString(cor(tf50pctmatz[finalsigtflist,"LIVER"],ourcontrolmeanmatz[,"LIVER"])),sep = "")) + theme(text = element_text(size = 15))
dev.off()

png(file = "Supplemental Figure S33H_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(liverexprenrichcordf,aes(x=Enrichment.Z.Score,y=Expression.Z.Score,label=TF)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_text_repel(max.overlaps = 15)+geom_smooth(method=lm,  linetype="dashed",color="darkred", fill="blue") + theme_classic() + ggtitle(paste("LIVER Correlation: ",toString(cor(tf50pctmatz[finalsigtflist,"LIVER"],ourcontrolmeanmatz[,"LIVER"])),sep = "")) + theme(text = element_text(size = 15))
dev.off()

# Supplemental Figure S33I
pdf(file = "Supplemental Figure S33I_021325.pdf",width = 6,height = 6)
ggplot(lungexprenrichcordf,aes(x=Enrichment.Z.Score,y=Expression.Z.Score,label=TF)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_text_repel(max.overlaps = 15)+geom_smooth(method=lm,  linetype="dashed",color="darkred", fill="blue") + theme_classic() + ggtitle(paste("LUNG Correlation: ",toString(cor(tf50pctmatz[finalsigtflist,"LUNG"],ourcontrolmeanmatz[,"LUNG"])),sep = "")) + theme(text = element_text(size = 15))
dev.off()

png(file = "Supplemental Figure S33I_021325.png",width = 6,height = 6,units = "in",res = 600)
ggplot(lungexprenrichcordf,aes(x=Enrichment.Z.Score,y=Expression.Z.Score,label=TF)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_text_repel(max.overlaps = 15)+geom_smooth(method=lm,  linetype="dashed",color="darkred", fill="blue") + theme_classic() + ggtitle(paste("LUNG Correlation: ",toString(cor(tf50pctmatz[finalsigtflist,"LUNG"],ourcontrolmeanmatz[,"LUNG"])),sep = "")) + theme(text = element_text(size = 15))
dev.off()


#####
# Supplemental Figure S28_S29_S30
####

ourtissue <- "SKM-GN"
ourdf <- data.frame(row.names = atacsigtflist,
                    "DEGaP Enrichment" = tf50logpmat[atacsigtflist,ourtissue],
                    "DAR Enrichment" = tfatacsigpmat[atacsigtflist,ourtissue],
                    "DMR Enrichment" = tfmethsigpmat[atacsigtflist,ourtissue],
                    "TF" = gsub("\\(.*","",atacsigtflist))
pdf(file = "Supplemental Figure S28A_021325.pdf",width = 6.5,height = 5.5)
ggplot(ourdf,aes(x=DEGaP.Enrichment,y=DAR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DEGaP.Enrichment,ourdf$DAR.Enrichment)),sep = ""))
dev.off()

pdf(file = "Supplemental Figure S29A_021325.pdf",width = 6.5,height = 5.5)
ggplot(ourdf,aes(x=DEGaP.Enrichment,y=DMR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DEGaP.Enrichment,ourdf$DMR.Enrichment)),sep = ""))
dev.off()

pdf(file = "Supplemental Figure S30A_021325.pdf",width = 6.5,height = 5.5)
ggplot(ourdf,aes(x=DAR.Enrichment,y=DMR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DAR.Enrichment,ourdf$DMR.Enrichment)),sep = ""))
dev.off()

png(file = "Supplemental Figure S28A_021325.png",width = 6.5,height = 5.5,units = "in",res = 600)
ggplot(ourdf,aes(x=DEGaP.Enrichment,y=DAR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DEGaP.Enrichment,ourdf$DAR.Enrichment)),sep = ""))
dev.off()

png(file = "Supplemental Figure S29A_021325.png",width = 6.5,height = 5.5,units = "in",res = 600)
ggplot(ourdf,aes(x=DEGaP.Enrichment,y=DMR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DEGaP.Enrichment,ourdf$DMR.Enrichment)),sep = ""))
dev.off()

png(file = "Supplemental Figure S30A_021325.png",width = 6.5,height = 5.5,units = "in",res = 600)
ggplot(ourdf,aes(x=DAR.Enrichment,y=DMR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DAR.Enrichment,ourdf$DMR.Enrichment)),sep = ""))
dev.off()


ourtissue <- "HEART"
ourdf <- data.frame(row.names = atacsigtflist,
                    "DEGaP Enrichment" = tf50logpmat[atacsigtflist,ourtissue],
                    "DAR Enrichment" = tfatacsigpmat[atacsigtflist,ourtissue],
                    "DMR Enrichment" = tfmethsigpmat[atacsigtflist,ourtissue],
                    "TF" = gsub("\\(.*","",atacsigtflist))
pdf(file = "Supplemental Figure S28B_021325.pdf",width = 6.5,height = 5.5)
ggplot(ourdf,aes(x=DEGaP.Enrichment,y=DAR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DEGaP.Enrichment,ourdf$DAR.Enrichment)),sep = ""))
dev.off()

pdf(file = "Supplemental Figure S29B_021325.pdf",width = 6.5,height = 5.5)
ggplot(ourdf,aes(x=DEGaP.Enrichment,y=DMR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DEGaP.Enrichment,ourdf$DMR.Enrichment)),sep = ""))
dev.off()

pdf(file = "Supplemental Figure S30B_021325.pdf",width = 6.5,height = 5.5)
ggplot(ourdf,aes(x=DAR.Enrichment,y=DMR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DAR.Enrichment,ourdf$DMR.Enrichment)),sep = ""))
dev.off()

png(file = "Supplemental Figure S28B_021325.png",width = 6.5,height = 5.5,units = "in",res = 600)
ggplot(ourdf,aes(x=DEGaP.Enrichment,y=DAR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DEGaP.Enrichment,ourdf$DAR.Enrichment)),sep = ""))
dev.off()

png(file = "Supplemental Figure S29B_021325.png",width = 6.5,height = 5.5,units = "in",res = 600)
ggplot(ourdf,aes(x=DEGaP.Enrichment,y=DMR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DEGaP.Enrichment,ourdf$DMR.Enrichment)),sep = ""))
dev.off()

png(file = "Supplemental Figure S30B_021325.png",width = 6.5,height = 5.5,units = "in",res = 600)
ggplot(ourdf,aes(x=DAR.Enrichment,y=DMR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DAR.Enrichment,ourdf$DMR.Enrichment)),sep = ""))
dev.off()

ourtissue <- "HIPPOC"
ourdf <- data.frame(row.names = atacsigtflist,
                    "DEGaP Enrichment" = tf50logpmat[atacsigtflist,ourtissue],
                    "DMR Enrichment" = tfmethsigpmat[atacsigtflist,ourtissue],
                    "TF" = gsub("\\(.*","",atacsigtflist))
pdf(file = "Supplemental Figure S29C_021325.pdf",width = 6.5,height = 5.5)
ggplot(ourdf,aes(x=DEGaP.Enrichment,y=DMR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DEGaP.Enrichment,ourdf$DMR.Enrichment)),sep = ""))
dev.off()

png(file = "Supplemental Figure S29C_021325.png",width = 6.5,height = 5.5,units = "in",res = 600)
ggplot(ourdf,aes(x=DEGaP.Enrichment,y=DMR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DEGaP.Enrichment,ourdf$DMR.Enrichment)),sep = ""))
dev.off()


ourtissue <- "KIDNEY"
ourdf <- data.frame(row.names = atacsigtflist,
                    "DEGaP Enrichment" = tf50logpmat[atacsigtflist,ourtissue],
                    "DAR Enrichment" = tfatacsigpmat[atacsigtflist,ourtissue],
                    "DMR Enrichment" = tfmethsigpmat[atacsigtflist,ourtissue],
                    "TF" = gsub("\\(.*","",atacsigtflist))
pdf(file = "Supplemental Figure S28C_021325.pdf",width = 6.5,height = 5.5)
ggplot(ourdf,aes(x=DEGaP.Enrichment,y=DAR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DEGaP.Enrichment,ourdf$DAR.Enrichment)),sep = ""))
dev.off()

pdf(file = "Supplemental Figure S29D_021325.pdf",width = 6.5,height = 5.5)
ggplot(ourdf,aes(x=DEGaP.Enrichment,y=DMR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DEGaP.Enrichment,ourdf$DMR.Enrichment)),sep = ""))
dev.off()

pdf(file = "Supplemental Figure S30C_021325.pdf",width = 6.5,height = 5.5)
ggplot(ourdf,aes(x=DAR.Enrichment,y=DMR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DAR.Enrichment,ourdf$DMR.Enrichment)),sep = ""))
dev.off()

png(file = "Supplemental Figure S28C_021325.png",width = 6.5,height = 5.5,units = "in",res = 600)
ggplot(ourdf,aes(x=DEGaP.Enrichment,y=DAR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DEGaP.Enrichment,ourdf$DAR.Enrichment)),sep = ""))
dev.off()

png(file = "Supplemental Figure S29D_021325.png",width = 6.5,height = 5.5,units = "in",res = 600)
ggplot(ourdf,aes(x=DEGaP.Enrichment,y=DMR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DEGaP.Enrichment,ourdf$DMR.Enrichment)),sep = ""))
dev.off()

png(file = "Supplemental Figure S30C_021325.png",width = 6.5,height = 5.5,units = "in",res = 600)
ggplot(ourdf,aes(x=DAR.Enrichment,y=DMR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DAR.Enrichment,ourdf$DMR.Enrichment)),sep = ""))
dev.off()


ourtissue <- "LIVER"
ourdf <- data.frame(row.names = atacsigtflist,
                    "DEGaP Enrichment" = tf50logpmat[atacsigtflist,ourtissue],
                    "DAR Enrichment" = tfatacsigpmat[atacsigtflist,ourtissue],
                    "DMR Enrichment" = tfmethsigpmat[atacsigtflist,ourtissue],
                    "TF" = gsub("\\(.*","",atacsigtflist))
pdf(file = "Supplemental Figure S28D_021325.pdf",width = 6.5,height = 5.5)
ggplot(ourdf,aes(x=DEGaP.Enrichment,y=DAR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DEGaP.Enrichment,ourdf$DAR.Enrichment)),sep = ""))
dev.off()

pdf(file = "Supplemental Figure S29E_021325.pdf",width = 6.5,height = 5.5)
ggplot(ourdf,aes(x=DEGaP.Enrichment,y=DMR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DEGaP.Enrichment,ourdf$DMR.Enrichment)),sep = ""))
dev.off()

pdf(file = "Supplemental Figure S30D_021325.pdf",width = 6.5,height = 5.5)
ggplot(ourdf,aes(x=DAR.Enrichment,y=DMR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DAR.Enrichment,ourdf$DMR.Enrichment)),sep = ""))
dev.off()

png(file = "Supplemental Figure S28D_021325.png",width = 6.5,height = 5.5,units = "in",res = 600)
ggplot(ourdf,aes(x=DEGaP.Enrichment,y=DAR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DEGaP.Enrichment,ourdf$DAR.Enrichment)),sep = ""))
dev.off()

png(file = "Supplemental Figure S29E_021325.png",width = 6.5,height = 5.5,units = "in",res = 600)
ggplot(ourdf,aes(x=DEGaP.Enrichment,y=DMR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DEGaP.Enrichment,ourdf$DMR.Enrichment)),sep = ""))
dev.off()

png(file = "Supplemental Figure S30D_021325.png",width = 6.5,height = 5.5,units = "in",res = 600)
ggplot(ourdf,aes(x=DAR.Enrichment,y=DMR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DAR.Enrichment,ourdf$DMR.Enrichment)),sep = ""))
dev.off()


ourtissue <- "LUNG"
ourdf <- data.frame(row.names = atacsigtflist,
                    "DEGaP Enrichment" = tf50logpmat[atacsigtflist,ourtissue],
                    "DAR Enrichment" = tfatacsigpmat[atacsigtflist,ourtissue],
                    "DMR Enrichment" = tfmethsigpmat[atacsigtflist,ourtissue],
                    "TF" = gsub("\\(.*","",atacsigtflist))
pdf(file = "Supplemental Figure S28E_021325.pdf",width = 6.5,height = 5.5)
ggplot(ourdf,aes(x=DEGaP.Enrichment,y=DAR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DEGaP.Enrichment,ourdf$DAR.Enrichment)),sep = ""))
dev.off()

pdf(file = "Supplemental Figure S29F_021325.pdf",width = 6.5,height = 5.5)
ggplot(ourdf,aes(x=DEGaP.Enrichment,y=DMR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DEGaP.Enrichment,ourdf$DMR.Enrichment)),sep = ""))
dev.off()

pdf(file = "Supplemental Figure S30E_021325.pdf",width = 6.5,height = 5.5)
ggplot(ourdf,aes(x=DAR.Enrichment,y=DMR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DAR.Enrichment,ourdf$DMR.Enrichment)),sep = ""))
dev.off()

png(file = "Supplemental Figure S28E_021325.png",width = 6.5,height = 5.5,units = "in",res = 600)
ggplot(ourdf,aes(x=DEGaP.Enrichment,y=DAR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DEGaP.Enrichment,ourdf$DAR.Enrichment)),sep = ""))
dev.off()

png(file = "Supplemental Figure S29F_021325.png",width = 6.5,height = 5.5,units = "in",res = 600)
ggplot(ourdf,aes(x=DEGaP.Enrichment,y=DMR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DEGaP.Enrichment,ourdf$DMR.Enrichment)),sep = ""))
dev.off()

png(file = "Supplemental Figure S30E_021325.png",width = 6.5,height = 5.5,units = "in",res = 600)
ggplot(ourdf,aes(x=DAR.Enrichment,y=DMR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DAR.Enrichment,ourdf$DMR.Enrichment)),sep = ""))
dev.off()


ourtissue <- "BAT"
ourdf <- data.frame(row.names = atacsigtflist,
                    "DEGaP Enrichment" = tf50logpmat[atacsigtflist,ourtissue],
                    "DAR Enrichment" = tfatacsigpmat[atacsigtflist,ourtissue],
                    "DMR Enrichment" = tfmethsigpmat[atacsigtflist,ourtissue],
                    "TF" = gsub("\\(.*","",atacsigtflist))
pdf(file = "Supplemental Figure S28F_021325.pdf",width = 6.5,height = 5.5)
ggplot(ourdf,aes(x=DEGaP.Enrichment,y=DAR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DEGaP.Enrichment,ourdf$DAR.Enrichment)),sep = ""))
dev.off()

pdf(file = "Supplemental Figure S29G_021325.pdf",width = 6.5,height = 5.5)
ggplot(ourdf,aes(x=DEGaP.Enrichment,y=DMR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DEGaP.Enrichment,ourdf$DMR.Enrichment)),sep = ""))
dev.off()

pdf(file = "Supplemental Figure S30F_021325.pdf",width = 6.5,height = 5.5)
ggplot(ourdf,aes(x=DAR.Enrichment,y=DMR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DAR.Enrichment,ourdf$DMR.Enrichment)),sep = ""))
dev.off()

png(file = "Supplemental Figure S28F_021325.png",width = 6.5,height = 5.5,units = "in",res = 600)
ggplot(ourdf,aes(x=DEGaP.Enrichment,y=DAR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DEGaP.Enrichment,ourdf$DAR.Enrichment)),sep = ""))
dev.off()

png(file = "Supplemental Figure S29G_021325.png",width = 6.5,height = 5.5,units = "in",res = 600)
ggplot(ourdf,aes(x=DEGaP.Enrichment,y=DMR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DEGaP.Enrichment,ourdf$DMR.Enrichment)),sep = ""))
dev.off()

png(file = "Supplemental Figure S30F_021325.png",width = 6.5,height = 5.5,units = "in",res = 600)
ggplot(ourdf,aes(x=DAR.Enrichment,y=DMR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DAR.Enrichment,ourdf$DMR.Enrichment)),sep = ""))
dev.off()

ourtissue <- "WAT-SC"
ourdf <- data.frame(row.names = atacsigtflist,
                    "DEGaP Enrichment" = tf50logpmat[atacsigtflist,ourtissue],
                    "DMR Enrichment" = tfmethsigpmat[atacsigtflist,ourtissue],
                    "TF" = gsub("\\(.*","",atacsigtflist))
pdf(file = "Supplemental Figure S29H_021325.pdf",width = 6.5,height = 5.5)
ggplot(ourdf,aes(x=DEGaP.Enrichment,y=DMR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DEGaP.Enrichment,ourdf$DMR.Enrichment)),sep = ""))
dev.off()

png(file = "Supplemental Figure S29H_021325.png",width = 6.5,height = 5.5,units = "in",res = 600)
ggplot(ourdf,aes(x=DEGaP.Enrichment,y=DMR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DEGaP.Enrichment,ourdf$DMR.Enrichment)),sep = ""))
dev.off()

gc()

save.image("Figure6_S27_to_S33_021325.RData")

#####
# Supplemental Figure S34
####

gastro50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% gastroactivepeaks,]
heart50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% heartactivepeaks,]
hippo50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% hippoactivepeaks,]
kidney50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% kidneyactivepeaks,]
liver50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% liveractivepeaks,]
lung50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% lungactivepeaks,]
brown50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% brownactivepeaks,]
white50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% whiteactivepeaks,]

####
# Generation of peak sets for significant genes (DEGaPs) for each tissue
#####

load("rnasigl2fcmat.RData")

gastrosigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrornasig,]),gastroactivepeaks)
gastrosigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrornasig & peakanno$trimanno %in% "Promoter",]),gastroactivepeaks)
gastrosigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrornasig & peakanno$trimanno %in% "Intron",]),gastroactivepeaks)
gastrosigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrornasig & peakanno$trimanno %in% "Distal Intergenic",]),gastroactivepeaks)
gastrosiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrornasig & peakanno$trimanno %in% "Upstream",]),gastroactivepeaks)
gastrosigdownpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrornasig & peakanno$trimanno %in% "Downstream",]),gastroactivepeaks)
gastrosigexpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrornasig & peakanno$trimanno %in% "Exon",]),gastroactivepeaks)
gastrosigpromproxpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrornasig & peakanno$custom_annotation %in% "Promoter (<=1kb)",]),gastroactivepeaks)
gastrosigpromfarpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrornasig & peakanno$custom_annotation %in% "Promoter (1-2kb)",]),gastroactivepeaks)

gastrow8upsig <- rownames(gastrosigrnal2fc[gastrosigrnal2fc[,"F W8"] > 0 & gastrosigrnal2fc[,"M W8"] > 0,])
gastrow8upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow8upsig,]),gastroactivepeaks)
gastrow8upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow8upsig & peakanno$trimanno %in% "Promoter",]),gastroactivepeaks)
gastrow8upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow8upsig & peakanno$trimanno %in% "Intron",]),gastroactivepeaks)
gastrow8upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow8upsig & peakanno$trimanno %in% "Distal Intergenic",]),gastroactivepeaks)
gastrow8upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow8upsig & peakanno$trimanno %in% "Upstream",]),gastroactivepeaks)

gastrow8downsig <- rownames(gastrosigrnal2fc[gastrosigrnal2fc[,"F W8"] < 0 & gastrosigrnal2fc[,"M W8"] < 0,])
gastrow8downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow8downsig,]),gastroactivepeaks)
gastrow8downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow8downsig & peakanno$trimanno %in% "Promoter",]),gastroactivepeaks)
gastrow8downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow8downsig & peakanno$trimanno %in% "Intron",]),gastroactivepeaks)
gastrow8downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow8downsig & peakanno$trimanno %in% "Distal Intergenic",]),gastroactivepeaks)
gastrow8downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow8downsig & peakanno$trimanno %in% "Upstream",]),gastroactivepeaks)

gastrow1upsig <- rownames(gastrosigrnal2fc[gastrosigrnal2fc[,"F W1"] > 0 & gastrosigrnal2fc[,"M W1"] > 0,])
gastrow1upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow1upsig,]),gastroactivepeaks)
gastrow1upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow1upsig & peakanno$trimanno %in% "Promoter",]),gastroactivepeaks)
gastrow1upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow1upsig & peakanno$trimanno %in% "Intron",]),gastroactivepeaks)
gastrow1upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow1upsig & peakanno$trimanno %in% "Distal Intergenic",]),gastroactivepeaks)
gastrow1upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow1upsig & peakanno$trimanno %in% "Upstream",]),gastroactivepeaks)

gastrow1downsig <- rownames(gastrosigrnal2fc[gastrosigrnal2fc[,"F W1"] < 0 & gastrosigrnal2fc[,"M W1"] < 0,])
gastrow1downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow1downsig,]),gastroactivepeaks)
gastrow1downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow1downsig & peakanno$trimanno %in% "Promoter",]),gastroactivepeaks)
gastrow1downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow1downsig & peakanno$trimanno %in% "Intron",]),gastroactivepeaks)
gastrow1downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow1downsig & peakanno$trimanno %in% "Distal Intergenic",]),gastroactivepeaks)
gastrow1downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow1downsig & peakanno$trimanno %in% "Upstream",]),gastroactivepeaks)


gastrow2upsig <- rownames(gastrosigrnal2fc[gastrosigrnal2fc[,"F W2"] > 0 & gastrosigrnal2fc[,"M W1"] > 0,])
gastrow2upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow2upsig,]),gastroactivepeaks)
gastrow2upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow2upsig & peakanno$trimanno %in% "Promoter",]),gastroactivepeaks)
gastrow2upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow2upsig & peakanno$trimanno %in% "Intron",]),gastroactivepeaks)
gastrow2upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow2upsig & peakanno$trimanno %in% "Distal Intergenic",]),gastroactivepeaks)
gastrow2upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow2upsig & peakanno$trimanno %in% "Upstream",]),gastroactivepeaks)

gastrow2downsig <- rownames(gastrosigrnal2fc[gastrosigrnal2fc[,"F W2"] < 0 & gastrosigrnal2fc[,"M W1"] < 0,])
gastrow2downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow2downsig,]),gastroactivepeaks)
gastrow2downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow2downsig & peakanno$trimanno %in% "Promoter",]),gastroactivepeaks)
gastrow2downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow2downsig & peakanno$trimanno %in% "Intron",]),gastroactivepeaks)
gastrow2downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow2downsig & peakanno$trimanno %in% "Distal Intergenic",]),gastroactivepeaks)
gastrow2downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow2downsig & peakanno$trimanno %in% "Upstream",]),gastroactivepeaks)


gastrow4upsig <- rownames(gastrosigrnal2fc[gastrosigrnal2fc[,"F W4"] > 0 & gastrosigrnal2fc[,"M W4"] > 0,])
gastrow4upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow4upsig,]),gastroactivepeaks)
gastrow4upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow4upsig & peakanno$trimanno %in% "Promoter",]),gastroactivepeaks)
gastrow4upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow4upsig & peakanno$trimanno %in% "Intron",]),gastroactivepeaks)
gastrow4upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow4upsig & peakanno$trimanno %in% "Distal Intergenic",]),gastroactivepeaks)
gastrow4upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow4upsig & peakanno$trimanno %in% "Upstream",]),gastroactivepeaks)

gastrow4downsig <- rownames(gastrosigrnal2fc[gastrosigrnal2fc[,"F W4"] < 0 & gastrosigrnal2fc[,"M W4"] < 0,])
gastrow4downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow4downsig,]),gastroactivepeaks)
gastrow4downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow4downsig & peakanno$trimanno %in% "Promoter",]),gastroactivepeaks)
gastrow4downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow4downsig & peakanno$trimanno %in% "Intron",]),gastroactivepeaks)
gastrow4downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow4downsig & peakanno$trimanno %in% "Distal Intergenic",]),gastroactivepeaks)
gastrow4downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow4downsig & peakanno$trimanno %in% "Upstream",]),gastroactivepeaks)


heartsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartrnasig,]),heartactivepeaks)
heartsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartrnasig & peakanno$trimanno %in% "Promoter",]),heartactivepeaks)
heartsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartrnasig & peakanno$trimanno %in% "Intron",]),heartactivepeaks)
heartsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartrnasig & peakanno$trimanno %in% "Distal Intergenic",]),heartactivepeaks)
heartsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartrnasig & peakanno$trimanno %in% "Upstream",]),heartactivepeaks)
heartsigdownpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartrnasig & peakanno$trimanno %in% "Downstream",]),heartactivepeaks)
heartsigexpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartrnasig & peakanno$trimanno %in% "Exon",]),heartactivepeaks)
heartsigpromproxpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartrnasig & peakanno$custom_annotation %in% "Promoter (<=1kb)",]),heartactivepeaks)
heartsigpromfarpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartrnasig & peakanno$custom_annotation %in% "Promoter (1-2kb)",]),heartactivepeaks)

heartw8upsig <- rownames(heartsigrnal2fc[heartsigrnal2fc[,"F W8"] > 0 & heartsigrnal2fc[,"M W8"] > 0,])
heartw8upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw8upsig,]),heartactivepeaks)
heartw8upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw8upsig & peakanno$trimanno %in% "Promoter",]),heartactivepeaks)
heartw8upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw8upsig & peakanno$trimanno %in% "Intron",]),heartactivepeaks)
heartw8upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw8upsig & peakanno$trimanno %in% "Distal Intergenic",]),heartactivepeaks)
heartw8upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw8upsig & peakanno$trimanno %in% "Upstream",]),heartactivepeaks)

heartw8downsig <- rownames(heartsigrnal2fc[heartsigrnal2fc[,"F W8"] < 0 & heartsigrnal2fc[,"M W8"] < 0,])
heartw8downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw8downsig,]),heartactivepeaks)
heartw8downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw8downsig & peakanno$trimanno %in% "Promoter",]),heartactivepeaks)
heartw8downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw8downsig & peakanno$trimanno %in% "Intron",]),heartactivepeaks)
heartw8downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw8downsig & peakanno$trimanno %in% "Distal Intergenic",]),heartactivepeaks)
heartw8downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw8downsig & peakanno$trimanno %in% "Upstream",]),heartactivepeaks)

heartw1upsig <- rownames(heartsigrnal2fc[heartsigrnal2fc[,"F W1"] > 0 & heartsigrnal2fc[,"M W1"] > 0,])
heartw1upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw1upsig,]),heartactivepeaks)
heartw1upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw1upsig & peakanno$trimanno %in% "Promoter",]),heartactivepeaks)
heartw1upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw1upsig & peakanno$trimanno %in% "Intron",]),heartactivepeaks)
heartw1upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw1upsig & peakanno$trimanno %in% "Distal Intergenic",]),heartactivepeaks)
heartw1upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw1upsig & peakanno$trimanno %in% "Upstream",]),heartactivepeaks)

heartw1downsig <- rownames(heartsigrnal2fc[heartsigrnal2fc[,"F W1"] < 0 & heartsigrnal2fc[,"M W1"] < 0,])
heartw1downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw1downsig,]),heartactivepeaks)
heartw1downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw1downsig & peakanno$trimanno %in% "Promoter",]),heartactivepeaks)
heartw1downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw1downsig & peakanno$trimanno %in% "Intron",]),heartactivepeaks)
heartw1downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw1downsig & peakanno$trimanno %in% "Distal Intergenic",]),heartactivepeaks)
heartw1downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw1downsig & peakanno$trimanno %in% "Upstream",]),heartactivepeaks)


heartw2upsig <- rownames(heartsigrnal2fc[heartsigrnal2fc[,"F W2"] > 0 & heartsigrnal2fc[,"M W1"] > 0,])
heartw2upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw2upsig,]),heartactivepeaks)
heartw2upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw2upsig & peakanno$trimanno %in% "Promoter",]),heartactivepeaks)
heartw2upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw2upsig & peakanno$trimanno %in% "Intron",]),heartactivepeaks)
heartw2upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw2upsig & peakanno$trimanno %in% "Distal Intergenic",]),heartactivepeaks)
heartw2upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw2upsig & peakanno$trimanno %in% "Upstream",]),heartactivepeaks)

heartw2downsig <- rownames(heartsigrnal2fc[heartsigrnal2fc[,"F W2"] < 0 & heartsigrnal2fc[,"M W1"] < 0,])
heartw2downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw2downsig,]),heartactivepeaks)
heartw2downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw2downsig & peakanno$trimanno %in% "Promoter",]),heartactivepeaks)
heartw2downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw2downsig & peakanno$trimanno %in% "Intron",]),heartactivepeaks)
heartw2downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw2downsig & peakanno$trimanno %in% "Distal Intergenic",]),heartactivepeaks)
heartw2downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw2downsig & peakanno$trimanno %in% "Upstream",]),heartactivepeaks)


heartw4upsig <- rownames(heartsigrnal2fc[heartsigrnal2fc[,"F W4"] > 0 & heartsigrnal2fc[,"M W4"] > 0,])
heartw4upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw4upsig,]),heartactivepeaks)
heartw4upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw4upsig & peakanno$trimanno %in% "Promoter",]),heartactivepeaks)
heartw4upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw4upsig & peakanno$trimanno %in% "Intron",]),heartactivepeaks)
heartw4upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw4upsig & peakanno$trimanno %in% "Distal Intergenic",]),heartactivepeaks)
heartw4upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw4upsig & peakanno$trimanno %in% "Upstream",]),heartactivepeaks)

heartw4downsig <- rownames(heartsigrnal2fc[heartsigrnal2fc[,"F W4"] < 0 & heartsigrnal2fc[,"M W4"] < 0,])
heartw4downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw4downsig,]),heartactivepeaks)
heartw4downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw4downsig & peakanno$trimanno %in% "Promoter",]),heartactivepeaks)
heartw4downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw4downsig & peakanno$trimanno %in% "Intron",]),heartactivepeaks)
heartw4downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw4downsig & peakanno$trimanno %in% "Distal Intergenic",]),heartactivepeaks)
heartw4downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw4downsig & peakanno$trimanno %in% "Upstream",]),heartactivepeaks)


hipposigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippornasig,]),hippoactivepeaks)
hipposigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippornasig & peakanno$trimanno %in% "Promoter",]),hippoactivepeaks)
hipposigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippornasig & peakanno$trimanno %in% "Intron",]),hippoactivepeaks)
hipposigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippornasig & peakanno$trimanno %in% "Distal Intergenic",]),hippoactivepeaks)
hipposiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippornasig & peakanno$trimanno %in% "Upstream",]),hippoactivepeaks)
hipposigdownpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippornasig & peakanno$trimanno %in% "Downstream",]),hippoactivepeaks)
hipposigexpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippornasig & peakanno$trimanno %in% "Exon",]),hippoactivepeaks)
hipposigpromproxpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippornasig & peakanno$custom_annotation %in% "Promoter (<=1kb)",]),hippoactivepeaks)
hipposigpromfarpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippornasig & peakanno$custom_annotation %in% "Promoter (1-2kb)",]),hippoactivepeaks)

hippow8upsig <- rownames(hipposigrnal2fc[hipposigrnal2fc[,"F W8"] > 0 & hipposigrnal2fc[,"M W8"] > 0,])
hippow8upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow8upsig,]),hippoactivepeaks)
hippow8upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow8upsig & peakanno$trimanno %in% "Promoter",]),hippoactivepeaks)
hippow8upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow8upsig & peakanno$trimanno %in% "Intron",]),hippoactivepeaks)
hippow8upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow8upsig & peakanno$trimanno %in% "Distal Intergenic",]),hippoactivepeaks)
hippow8upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow8upsig & peakanno$trimanno %in% "Upstream",]),hippoactivepeaks)

hippow8downsig <- rownames(hipposigrnal2fc[hipposigrnal2fc[,"F W8"] < 0 & hipposigrnal2fc[,"M W8"] < 0,])
hippow8downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow8downsig,]),hippoactivepeaks)
hippow8downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow8downsig & peakanno$trimanno %in% "Promoter",]),hippoactivepeaks)
hippow8downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow8downsig & peakanno$trimanno %in% "Intron",]),hippoactivepeaks)
hippow8downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow8downsig & peakanno$trimanno %in% "Distal Intergenic",]),hippoactivepeaks)
hippow8downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow8downsig & peakanno$trimanno %in% "Upstream",]),hippoactivepeaks)

hippow1upsig <- rownames(hipposigrnal2fc[hipposigrnal2fc[,"F W1"] > 0 & hipposigrnal2fc[,"M W1"] > 0,])
hippow1upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow1upsig,]),hippoactivepeaks)
hippow1upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow1upsig & peakanno$trimanno %in% "Promoter",]),hippoactivepeaks)
hippow1upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow1upsig & peakanno$trimanno %in% "Intron",]),hippoactivepeaks)
hippow1upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow1upsig & peakanno$trimanno %in% "Distal Intergenic",]),hippoactivepeaks)
hippow1upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow1upsig & peakanno$trimanno %in% "Upstream",]),hippoactivepeaks)

hippow1downsig <- rownames(hipposigrnal2fc[hipposigrnal2fc[,"F W1"] < 0 & hipposigrnal2fc[,"M W1"] < 0,])
hippow1downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow1downsig,]),hippoactivepeaks)
hippow1downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow1downsig & peakanno$trimanno %in% "Promoter",]),hippoactivepeaks)
hippow1downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow1downsig & peakanno$trimanno %in% "Intron",]),hippoactivepeaks)
hippow1downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow1downsig & peakanno$trimanno %in% "Distal Intergenic",]),hippoactivepeaks)
hippow1downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow1downsig & peakanno$trimanno %in% "Upstream",]),hippoactivepeaks)


hippow2upsig <- rownames(hipposigrnal2fc[hipposigrnal2fc[,"F W2"] > 0 & hipposigrnal2fc[,"M W1"] > 0,])
hippow2upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow2upsig,]),hippoactivepeaks)
hippow2upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow2upsig & peakanno$trimanno %in% "Promoter",]),hippoactivepeaks)
hippow2upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow2upsig & peakanno$trimanno %in% "Intron",]),hippoactivepeaks)
hippow2upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow2upsig & peakanno$trimanno %in% "Distal Intergenic",]),hippoactivepeaks)
hippow2upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow2upsig & peakanno$trimanno %in% "Upstream",]),hippoactivepeaks)

hippow2downsig <- rownames(hipposigrnal2fc[hipposigrnal2fc[,"F W2"] < 0 & hipposigrnal2fc[,"M W1"] < 0,])
hippow2downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow2downsig,]),hippoactivepeaks)
hippow2downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow2downsig & peakanno$trimanno %in% "Promoter",]),hippoactivepeaks)
hippow2downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow2downsig & peakanno$trimanno %in% "Intron",]),hippoactivepeaks)
hippow2downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow2downsig & peakanno$trimanno %in% "Distal Intergenic",]),hippoactivepeaks)
hippow2downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow2downsig & peakanno$trimanno %in% "Upstream",]),hippoactivepeaks)


hippow4upsig <- rownames(hipposigrnal2fc[hipposigrnal2fc[,"F W4"] > 0 & hipposigrnal2fc[,"M W4"] > 0,])
hippow4upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow4upsig,]),hippoactivepeaks)
hippow4upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow4upsig & peakanno$trimanno %in% "Promoter",]),hippoactivepeaks)
hippow4upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow4upsig & peakanno$trimanno %in% "Intron",]),hippoactivepeaks)
hippow4upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow4upsig & peakanno$trimanno %in% "Distal Intergenic",]),hippoactivepeaks)
hippow4upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow4upsig & peakanno$trimanno %in% "Upstream",]),hippoactivepeaks)

hippow4downsig <- rownames(hipposigrnal2fc[hipposigrnal2fc[,"F W4"] < 0 & hipposigrnal2fc[,"M W4"] < 0,])
hippow4downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow4downsig,]),hippoactivepeaks)
hippow4downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow4downsig & peakanno$trimanno %in% "Promoter",]),hippoactivepeaks)
hippow4downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow4downsig & peakanno$trimanno %in% "Intron",]),hippoactivepeaks)
hippow4downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow4downsig & peakanno$trimanno %in% "Distal Intergenic",]),hippoactivepeaks)
hippow4downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow4downsig & peakanno$trimanno %in% "Upstream",]),hippoactivepeaks)


kidneysigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyrnasig,]),kidneyactivepeaks)
kidneysigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyrnasig & peakanno$trimanno %in% "Promoter",]),kidneyactivepeaks)
kidneysigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyrnasig & peakanno$trimanno %in% "Intron",]),kidneyactivepeaks)
kidneysigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyrnasig & peakanno$trimanno %in% "Distal Intergenic",]),kidneyactivepeaks)
kidneysiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyrnasig & peakanno$trimanno %in% "Upstream",]),kidneyactivepeaks)
kidneysigdownpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyrnasig & peakanno$trimanno %in% "Downstream",]),kidneyactivepeaks)
kidneysigexpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyrnasig & peakanno$trimanno %in% "Exon",]),kidneyactivepeaks)
kidneysigpromproxpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyrnasig & peakanno$custom_annotation %in% "Promoter (<=1kb)",]),kidneyactivepeaks)
kidneysigpromfarpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyrnasig & peakanno$custom_annotation %in% "Promoter (1-2kb)",]),kidneyactivepeaks)

kidneyw8upsig <- rownames(kidneysigrnal2fc[kidneysigrnal2fc[,"F W8"] > 0 & kidneysigrnal2fc[,"M W8"] > 0,])
kidneyw8upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw8upsig,]),kidneyactivepeaks)
kidneyw8upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw8upsig & peakanno$trimanno %in% "Promoter",]),kidneyactivepeaks)
kidneyw8upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw8upsig & peakanno$trimanno %in% "Intron",]),kidneyactivepeaks)
kidneyw8upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw8upsig & peakanno$trimanno %in% "Distal Intergenic",]),kidneyactivepeaks)
kidneyw8upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw8upsig & peakanno$trimanno %in% "Upstream",]),kidneyactivepeaks)

kidneyw8downsig <- rownames(kidneysigrnal2fc[kidneysigrnal2fc[,"F W8"] < 0 & kidneysigrnal2fc[,"M W8"] < 0,])
kidneyw8downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw8downsig,]),kidneyactivepeaks)
kidneyw8downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw8downsig & peakanno$trimanno %in% "Promoter",]),kidneyactivepeaks)
kidneyw8downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw8downsig & peakanno$trimanno %in% "Intron",]),kidneyactivepeaks)
kidneyw8downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw8downsig & peakanno$trimanno %in% "Distal Intergenic",]),kidneyactivepeaks)
kidneyw8downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw8downsig & peakanno$trimanno %in% "Upstream",]),kidneyactivepeaks)

kidneyw1upsig <- rownames(kidneysigrnal2fc[kidneysigrnal2fc[,"F W1"] > 0 & kidneysigrnal2fc[,"M W1"] > 0,])
kidneyw1upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw1upsig,]),kidneyactivepeaks)
kidneyw1upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw1upsig & peakanno$trimanno %in% "Promoter",]),kidneyactivepeaks)
kidneyw1upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw1upsig & peakanno$trimanno %in% "Intron",]),kidneyactivepeaks)
kidneyw1upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw1upsig & peakanno$trimanno %in% "Distal Intergenic",]),kidneyactivepeaks)
kidneyw1upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw1upsig & peakanno$trimanno %in% "Upstream",]),kidneyactivepeaks)

kidneyw1downsig <- rownames(kidneysigrnal2fc[kidneysigrnal2fc[,"F W1"] < 0 & kidneysigrnal2fc[,"M W1"] < 0,])
kidneyw1downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw1downsig,]),kidneyactivepeaks)
kidneyw1downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw1downsig & peakanno$trimanno %in% "Promoter",]),kidneyactivepeaks)
kidneyw1downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw1downsig & peakanno$trimanno %in% "Intron",]),kidneyactivepeaks)
kidneyw1downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw1downsig & peakanno$trimanno %in% "Distal Intergenic",]),kidneyactivepeaks)
kidneyw1downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw1downsig & peakanno$trimanno %in% "Upstream",]),kidneyactivepeaks)


kidneyw2upsig <- rownames(kidneysigrnal2fc[kidneysigrnal2fc[,"F W2"] > 0 & kidneysigrnal2fc[,"M W1"] > 0,])
kidneyw2upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw2upsig,]),kidneyactivepeaks)
kidneyw2upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw2upsig & peakanno$trimanno %in% "Promoter",]),kidneyactivepeaks)
kidneyw2upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw2upsig & peakanno$trimanno %in% "Intron",]),kidneyactivepeaks)
kidneyw2upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw2upsig & peakanno$trimanno %in% "Distal Intergenic",]),kidneyactivepeaks)
kidneyw2upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw2upsig & peakanno$trimanno %in% "Upstream",]),kidneyactivepeaks)

kidneyw2downsig <- rownames(kidneysigrnal2fc[kidneysigrnal2fc[,"F W2"] < 0 & kidneysigrnal2fc[,"M W1"] < 0,])
kidneyw2downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw2downsig,]),kidneyactivepeaks)
kidneyw2downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw2downsig & peakanno$trimanno %in% "Promoter",]),kidneyactivepeaks)
kidneyw2downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw2downsig & peakanno$trimanno %in% "Intron",]),kidneyactivepeaks)
kidneyw2downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw2downsig & peakanno$trimanno %in% "Distal Intergenic",]),kidneyactivepeaks)
kidneyw2downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw2downsig & peakanno$trimanno %in% "Upstream",]),kidneyactivepeaks)


kidneyw4upsig <- rownames(kidneysigrnal2fc[kidneysigrnal2fc[,"F W4"] > 0 & kidneysigrnal2fc[,"M W4"] > 0,])
kidneyw4upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw4upsig,]),kidneyactivepeaks)
kidneyw4upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw4upsig & peakanno$trimanno %in% "Promoter",]),kidneyactivepeaks)
kidneyw4upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw4upsig & peakanno$trimanno %in% "Intron",]),kidneyactivepeaks)
kidneyw4upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw4upsig & peakanno$trimanno %in% "Distal Intergenic",]),kidneyactivepeaks)
kidneyw4upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw4upsig & peakanno$trimanno %in% "Upstream",]),kidneyactivepeaks)

kidneyw4downsig <- rownames(kidneysigrnal2fc[kidneysigrnal2fc[,"F W4"] < 0 & kidneysigrnal2fc[,"M W4"] < 0,])
kidneyw4downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw4downsig,]),kidneyactivepeaks)
kidneyw4downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw4downsig & peakanno$trimanno %in% "Promoter",]),kidneyactivepeaks)
kidneyw4downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw4downsig & peakanno$trimanno %in% "Intron",]),kidneyactivepeaks)
kidneyw4downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw4downsig & peakanno$trimanno %in% "Distal Intergenic",]),kidneyactivepeaks)
kidneyw4downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw4downsig & peakanno$trimanno %in% "Upstream",]),kidneyactivepeaks)

liversigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverrnasig,]),liveractivepeaks)
liversigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverrnasig & peakanno$trimanno %in% "Promoter",]),liveractivepeaks)
liversigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverrnasig & peakanno$trimanno %in% "Intron",]),liveractivepeaks)
liversigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverrnasig & peakanno$trimanno %in% "Distal Intergenic",]),liveractivepeaks)
liversiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverrnasig & peakanno$trimanno %in% "Upstream",]),liveractivepeaks)
liversigdownpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverrnasig & peakanno$trimanno %in% "Downstream",]),liveractivepeaks)
liversigexpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverrnasig & peakanno$trimanno %in% "Exon",]),liveractivepeaks)
liversigpromproxpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverrnasig & peakanno$custom_annotation %in% "Promoter (<=1kb)",]),liveractivepeaks)
liversigpromfarpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverrnasig & peakanno$custom_annotation %in% "Promoter (1-2kb)",]),liveractivepeaks)

liverw8upsig <- rownames(liversigrnal2fc[liversigrnal2fc[,"F W8"] > 0 & liversigrnal2fc[,"M W8"] > 0,])
liverw8upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw8upsig,]),liveractivepeaks)
liverw8upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw8upsig & peakanno$trimanno %in% "Promoter",]),liveractivepeaks)
liverw8upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw8upsig & peakanno$trimanno %in% "Intron",]),liveractivepeaks)
liverw8upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw8upsig & peakanno$trimanno %in% "Distal Intergenic",]),liveractivepeaks)
liverw8upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw8upsig & peakanno$trimanno %in% "Upstream",]),liveractivepeaks)

liverw8downsig <- rownames(liversigrnal2fc[liversigrnal2fc[,"F W8"] < 0 & liversigrnal2fc[,"M W8"] < 0,])
liverw8downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw8downsig,]),liveractivepeaks)
liverw8downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw8downsig & peakanno$trimanno %in% "Promoter",]),liveractivepeaks)
liverw8downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw8downsig & peakanno$trimanno %in% "Intron",]),liveractivepeaks)
liverw8downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw8downsig & peakanno$trimanno %in% "Distal Intergenic",]),liveractivepeaks)
liverw8downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw8downsig & peakanno$trimanno %in% "Upstream",]),liveractivepeaks)

liverw1upsig <- rownames(liversigrnal2fc[liversigrnal2fc[,"F W1"] > 0 & liversigrnal2fc[,"M W1"] > 0,])
liverw1upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw1upsig,]),liveractivepeaks)
liverw1upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw1upsig & peakanno$trimanno %in% "Promoter",]),liveractivepeaks)
liverw1upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw1upsig & peakanno$trimanno %in% "Intron",]),liveractivepeaks)
liverw1upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw1upsig & peakanno$trimanno %in% "Distal Intergenic",]),liveractivepeaks)
liverw1upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw1upsig & peakanno$trimanno %in% "Upstream",]),liveractivepeaks)

liverw1downsig <- rownames(liversigrnal2fc[liversigrnal2fc[,"F W1"] < 0 & liversigrnal2fc[,"M W1"] < 0,])
liverw1downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw1downsig,]),liveractivepeaks)
liverw1downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw1downsig & peakanno$trimanno %in% "Promoter",]),liveractivepeaks)
liverw1downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw1downsig & peakanno$trimanno %in% "Intron",]),liveractivepeaks)
liverw1downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw1downsig & peakanno$trimanno %in% "Distal Intergenic",]),liveractivepeaks)
liverw1downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw1downsig & peakanno$trimanno %in% "Upstream",]),liveractivepeaks)


liverw2upsig <- rownames(liversigrnal2fc[liversigrnal2fc[,"F W2"] > 0 & liversigrnal2fc[,"M W1"] > 0,])
liverw2upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw2upsig,]),liveractivepeaks)
liverw2upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw2upsig & peakanno$trimanno %in% "Promoter",]),liveractivepeaks)
liverw2upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw2upsig & peakanno$trimanno %in% "Intron",]),liveractivepeaks)
liverw2upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw2upsig & peakanno$trimanno %in% "Distal Intergenic",]),liveractivepeaks)
liverw2upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw2upsig & peakanno$trimanno %in% "Upstream",]),liveractivepeaks)

liverw2downsig <- rownames(liversigrnal2fc[liversigrnal2fc[,"F W2"] < 0 & liversigrnal2fc[,"M W1"] < 0,])
liverw2downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw2downsig,]),liveractivepeaks)
liverw2downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw2downsig & peakanno$trimanno %in% "Promoter",]),liveractivepeaks)
liverw2downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw2downsig & peakanno$trimanno %in% "Intron",]),liveractivepeaks)
liverw2downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw2downsig & peakanno$trimanno %in% "Distal Intergenic",]),liveractivepeaks)
liverw2downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw2downsig & peakanno$trimanno %in% "Upstream",]),liveractivepeaks)


liverw4upsig <- rownames(liversigrnal2fc[liversigrnal2fc[,"F W4"] > 0 & liversigrnal2fc[,"M W4"] > 0,])
liverw4upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw4upsig,]),liveractivepeaks)
liverw4upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw4upsig & peakanno$trimanno %in% "Promoter",]),liveractivepeaks)
liverw4upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw4upsig & peakanno$trimanno %in% "Intron",]),liveractivepeaks)
liverw4upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw4upsig & peakanno$trimanno %in% "Distal Intergenic",]),liveractivepeaks)
liverw4upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw4upsig & peakanno$trimanno %in% "Upstream",]),liveractivepeaks)

liverw4downsig <- rownames(liversigrnal2fc[liversigrnal2fc[,"F W4"] < 0 & liversigrnal2fc[,"M W4"] < 0,])
liverw4downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw4downsig,]),liveractivepeaks)
liverw4downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw4downsig & peakanno$trimanno %in% "Promoter",]),liveractivepeaks)
liverw4downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw4downsig & peakanno$trimanno %in% "Intron",]),liveractivepeaks)
liverw4downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw4downsig & peakanno$trimanno %in% "Distal Intergenic",]),liveractivepeaks)
liverw4downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw4downsig & peakanno$trimanno %in% "Upstream",]),liveractivepeaks)

lungsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungrnasig,]),lungactivepeaks)
lungsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungrnasig & peakanno$trimanno %in% "Promoter",]),lungactivepeaks)
lungsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungrnasig & peakanno$trimanno %in% "Intron",]),lungactivepeaks)
lungsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungrnasig & peakanno$trimanno %in% "Distal Intergenic",]),lungactivepeaks)
lungsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungrnasig & peakanno$trimanno %in% "Upstream",]),lungactivepeaks)
lungsigdownpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungrnasig & peakanno$trimanno %in% "Downstream",]),lungactivepeaks)
lungsigexpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungrnasig & peakanno$trimanno %in% "Exon",]),lungactivepeaks)
lungsigpromproxpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungrnasig & peakanno$custom_annotation %in% "Promoter (<=1kb)",]),lungactivepeaks)
lungsigpromfarpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungrnasig & peakanno$custom_annotation %in% "Promoter (1-2kb)",]),lungactivepeaks)

lungw8upsig <- rownames(lungsigrnal2fc[lungsigrnal2fc[,"F W8"] > 0 & lungsigrnal2fc[,"M W8"] > 0,])
lungw8upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw8upsig,]),lungactivepeaks)
lungw8upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw8upsig & peakanno$trimanno %in% "Promoter",]),lungactivepeaks)
lungw8upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw8upsig & peakanno$trimanno %in% "Intron",]),lungactivepeaks)
lungw8upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw8upsig & peakanno$trimanno %in% "Distal Intergenic",]),lungactivepeaks)
lungw8upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw8upsig & peakanno$trimanno %in% "Upstream",]),lungactivepeaks)

lungw8downsig <- rownames(lungsigrnal2fc[lungsigrnal2fc[,"F W8"] < 0 & lungsigrnal2fc[,"M W8"] < 0,])
lungw8downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw8downsig,]),lungactivepeaks)
lungw8downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw8downsig & peakanno$trimanno %in% "Promoter",]),lungactivepeaks)
lungw8downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw8downsig & peakanno$trimanno %in% "Intron",]),lungactivepeaks)
lungw8downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw8downsig & peakanno$trimanno %in% "Distal Intergenic",]),lungactivepeaks)
lungw8downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw8downsig & peakanno$trimanno %in% "Upstream",]),lungactivepeaks)

lungw1upsig <- rownames(lungsigrnal2fc[lungsigrnal2fc[,"F W1"] > 0 & lungsigrnal2fc[,"M W1"] > 0,])
lungw1upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw1upsig,]),lungactivepeaks)
lungw1upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw1upsig & peakanno$trimanno %in% "Promoter",]),lungactivepeaks)
lungw1upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw1upsig & peakanno$trimanno %in% "Intron",]),lungactivepeaks)
lungw1upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw1upsig & peakanno$trimanno %in% "Distal Intergenic",]),lungactivepeaks)
lungw1upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw1upsig & peakanno$trimanno %in% "Upstream",]),lungactivepeaks)

lungw1downsig <- rownames(lungsigrnal2fc[lungsigrnal2fc[,"F W1"] < 0 & lungsigrnal2fc[,"M W1"] < 0,])
lungw1downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw1downsig,]),lungactivepeaks)
lungw1downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw1downsig & peakanno$trimanno %in% "Promoter",]),lungactivepeaks)
lungw1downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw1downsig & peakanno$trimanno %in% "Intron",]),lungactivepeaks)
lungw1downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw1downsig & peakanno$trimanno %in% "Distal Intergenic",]),lungactivepeaks)
lungw1downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw1downsig & peakanno$trimanno %in% "Upstream",]),lungactivepeaks)


lungw2upsig <- rownames(lungsigrnal2fc[lungsigrnal2fc[,"F W2"] > 0 & lungsigrnal2fc[,"M W1"] > 0,])
lungw2upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw2upsig,]),lungactivepeaks)
lungw2upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw2upsig & peakanno$trimanno %in% "Promoter",]),lungactivepeaks)
lungw2upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw2upsig & peakanno$trimanno %in% "Intron",]),lungactivepeaks)
lungw2upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw2upsig & peakanno$trimanno %in% "Distal Intergenic",]),lungactivepeaks)
lungw2upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw2upsig & peakanno$trimanno %in% "Upstream",]),lungactivepeaks)

lungw2downsig <- rownames(lungsigrnal2fc[lungsigrnal2fc[,"F W2"] < 0 & lungsigrnal2fc[,"M W1"] < 0,])
lungw2downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw2downsig,]),lungactivepeaks)
lungw2downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw2downsig & peakanno$trimanno %in% "Promoter",]),lungactivepeaks)
lungw2downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw2downsig & peakanno$trimanno %in% "Intron",]),lungactivepeaks)
lungw2downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw2downsig & peakanno$trimanno %in% "Distal Intergenic",]),lungactivepeaks)
lungw2downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw2downsig & peakanno$trimanno %in% "Upstream",]),lungactivepeaks)


lungw4upsig <- rownames(lungsigrnal2fc[lungsigrnal2fc[,"F W4"] > 0 & lungsigrnal2fc[,"M W4"] > 0,])
lungw4upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw4upsig,]),lungactivepeaks)
lungw4upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw4upsig & peakanno$trimanno %in% "Promoter",]),lungactivepeaks)
lungw4upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw4upsig & peakanno$trimanno %in% "Intron",]),lungactivepeaks)
lungw4upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw4upsig & peakanno$trimanno %in% "Distal Intergenic",]),lungactivepeaks)
lungw4upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw4upsig & peakanno$trimanno %in% "Upstream",]),lungactivepeaks)

lungw4downsig <- rownames(lungsigrnal2fc[lungsigrnal2fc[,"F W4"] < 0 & lungsigrnal2fc[,"M W4"] < 0,])
lungw4downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw4downsig,]),lungactivepeaks)
lungw4downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw4downsig & peakanno$trimanno %in% "Promoter",]),lungactivepeaks)
lungw4downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw4downsig & peakanno$trimanno %in% "Intron",]),lungactivepeaks)
lungw4downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw4downsig & peakanno$trimanno %in% "Distal Intergenic",]),lungactivepeaks)
lungw4downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw4downsig & peakanno$trimanno %in% "Upstream",]),lungactivepeaks)


brownsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownrnasig,]),brownactivepeaks)
brownsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownrnasig & peakanno$trimanno %in% "Promoter",]),brownactivepeaks)
brownsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownrnasig & peakanno$trimanno %in% "Intron",]),brownactivepeaks)
brownsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownrnasig & peakanno$trimanno %in% "Distal Intergenic",]),brownactivepeaks)
brownsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownrnasig & peakanno$trimanno %in% "Upstream",]),brownactivepeaks)
brownsigdownpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownrnasig & peakanno$trimanno %in% "Downstream",]),brownactivepeaks)
brownsigexpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownrnasig & peakanno$trimanno %in% "Exon",]),brownactivepeaks)
brownsigpromproxpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownrnasig & peakanno$custom_annotation %in% "Promoter (<=1kb)",]),brownactivepeaks)
brownsigpromfarpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownrnasig & peakanno$custom_annotation %in% "Promoter (1-2kb)",]),brownactivepeaks)

brownw8upsig <- rownames(brownsigrnal2fc[brownsigrnal2fc[,"F W8"] > 0 & brownsigrnal2fc[,"M W8"] > 0,])
brownw8upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw8upsig,]),brownactivepeaks)
brownw8upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw8upsig & peakanno$trimanno %in% "Promoter",]),brownactivepeaks)
brownw8upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw8upsig & peakanno$trimanno %in% "Intron",]),brownactivepeaks)
brownw8upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw8upsig & peakanno$trimanno %in% "Distal Intergenic",]),brownactivepeaks)
brownw8upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw8upsig & peakanno$trimanno %in% "Upstream",]),brownactivepeaks)

brownw8downsig <- rownames(brownsigrnal2fc[brownsigrnal2fc[,"F W8"] < 0 & brownsigrnal2fc[,"M W8"] < 0,])
brownw8downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw8downsig,]),brownactivepeaks)
brownw8downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw8downsig & peakanno$trimanno %in% "Promoter",]),brownactivepeaks)
brownw8downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw8downsig & peakanno$trimanno %in% "Intron",]),brownactivepeaks)
brownw8downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw8downsig & peakanno$trimanno %in% "Distal Intergenic",]),brownactivepeaks)
brownw8downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw8downsig & peakanno$trimanno %in% "Upstream",]),brownactivepeaks)

brownw1upsig <- rownames(brownsigrnal2fc[brownsigrnal2fc[,"F W1"] > 0 & brownsigrnal2fc[,"M W1"] > 0,])
brownw1upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw1upsig,]),brownactivepeaks)
brownw1upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw1upsig & peakanno$trimanno %in% "Promoter",]),brownactivepeaks)
brownw1upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw1upsig & peakanno$trimanno %in% "Intron",]),brownactivepeaks)
brownw1upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw1upsig & peakanno$trimanno %in% "Distal Intergenic",]),brownactivepeaks)
brownw1upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw1upsig & peakanno$trimanno %in% "Upstream",]),brownactivepeaks)

brownw1downsig <- rownames(brownsigrnal2fc[brownsigrnal2fc[,"F W1"] < 0 & brownsigrnal2fc[,"M W1"] < 0,])
brownw1downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw1downsig,]),brownactivepeaks)
brownw1downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw1downsig & peakanno$trimanno %in% "Promoter",]),brownactivepeaks)
brownw1downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw1downsig & peakanno$trimanno %in% "Intron",]),brownactivepeaks)
brownw1downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw1downsig & peakanno$trimanno %in% "Distal Intergenic",]),brownactivepeaks)
brownw1downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw1downsig & peakanno$trimanno %in% "Upstream",]),brownactivepeaks)


brownw2upsig <- rownames(brownsigrnal2fc[brownsigrnal2fc[,"F W2"] > 0 & brownsigrnal2fc[,"M W1"] > 0,])
brownw2upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw2upsig,]),brownactivepeaks)
brownw2upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw2upsig & peakanno$trimanno %in% "Promoter",]),brownactivepeaks)
brownw2upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw2upsig & peakanno$trimanno %in% "Intron",]),brownactivepeaks)
brownw2upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw2upsig & peakanno$trimanno %in% "Distal Intergenic",]),brownactivepeaks)
brownw2upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw2upsig & peakanno$trimanno %in% "Upstream",]),brownactivepeaks)

brownw2downsig <- rownames(brownsigrnal2fc[brownsigrnal2fc[,"F W2"] < 0 & brownsigrnal2fc[,"M W1"] < 0,])
brownw2downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw2downsig,]),brownactivepeaks)
brownw2downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw2downsig & peakanno$trimanno %in% "Promoter",]),brownactivepeaks)
brownw2downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw2downsig & peakanno$trimanno %in% "Intron",]),brownactivepeaks)
brownw2downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw2downsig & peakanno$trimanno %in% "Distal Intergenic",]),brownactivepeaks)
brownw2downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw2downsig & peakanno$trimanno %in% "Upstream",]),brownactivepeaks)

brownw4upsig <- rownames(brownsigrnal2fc[brownsigrnal2fc[,"F W4"] > 0 & brownsigrnal2fc[,"M W4"] > 0,])
brownw4upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw4upsig,]),brownactivepeaks)
brownw4upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw4upsig & peakanno$trimanno %in% "Promoter",]),brownactivepeaks)
brownw4upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw4upsig & peakanno$trimanno %in% "Intron",]),brownactivepeaks)
brownw4upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw4upsig & peakanno$trimanno %in% "Distal Intergenic",]),brownactivepeaks)
brownw4upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw4upsig & peakanno$trimanno %in% "Upstream",]),brownactivepeaks)

brownw4downsig <- rownames(brownsigrnal2fc[brownsigrnal2fc[,"F W4"] < 0 & brownsigrnal2fc[,"M W4"] < 0,])
brownw4downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw4downsig,]),brownactivepeaks)
brownw4downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw4downsig & peakanno$trimanno %in% "Promoter",]),brownactivepeaks)
brownw4downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw4downsig & peakanno$trimanno %in% "Intron",]),brownactivepeaks)
brownw4downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw4downsig & peakanno$trimanno %in% "Distal Intergenic",]),brownactivepeaks)
brownw4downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw4downsig & peakanno$trimanno %in% "Upstream",]),brownactivepeaks)

whitesigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whiternasig,]),whiteactivepeaks)
whitesigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whiternasig & peakanno$trimanno %in% "Promoter",]),whiteactivepeaks)
whitesigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whiternasig & peakanno$trimanno %in% "Intron",]),whiteactivepeaks)
whitesigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whiternasig & peakanno$trimanno %in% "Distal Intergenic",]),whiteactivepeaks)
whitesiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whiternasig & peakanno$trimanno %in% "Upstream",]),whiteactivepeaks)
whitesigdownpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whiternasig & peakanno$trimanno %in% "Downstream",]),whiteactivepeaks)
whitesigexpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whiternasig & peakanno$trimanno %in% "Exon",]),whiteactivepeaks)
whitesigpromproxpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whiternasig & peakanno$custom_annotation %in% "Promoter (<=1kb)",]),whiteactivepeaks)
whitesigpromfarpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whiternasig & peakanno$custom_annotation %in% "Promoter (1-2kb)",]),whiteactivepeaks)

whitew8upsig <- rownames(whitesigrnal2fc[whitesigrnal2fc[,"F W8"] > 0 & whitesigrnal2fc[,"M W8"] > 0,])
whitew8upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew8upsig,]),whiteactivepeaks)
whitew8upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew8upsig & peakanno$trimanno %in% "Promoter",]),whiteactivepeaks)
whitew8upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew8upsig & peakanno$trimanno %in% "Intron",]),whiteactivepeaks)
whitew8upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew8upsig & peakanno$trimanno %in% "Distal Intergenic",]),whiteactivepeaks)
whitew8upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew8upsig & peakanno$trimanno %in% "Upstream",]),whiteactivepeaks)

whitew8downsig <- rownames(whitesigrnal2fc[whitesigrnal2fc[,"F W8"] < 0 & whitesigrnal2fc[,"M W8"] < 0,])
whitew8downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew8downsig,]),whiteactivepeaks)
whitew8downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew8downsig & peakanno$trimanno %in% "Promoter",]),whiteactivepeaks)
whitew8downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew8downsig & peakanno$trimanno %in% "Intron",]),whiteactivepeaks)
whitew8downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew8downsig & peakanno$trimanno %in% "Distal Intergenic",]),whiteactivepeaks)
whitew8downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew8downsig & peakanno$trimanno %in% "Upstream",]),whiteactivepeaks)

whitew1upsig <- rownames(whitesigrnal2fc[whitesigrnal2fc[,"F W1"] > 0 & whitesigrnal2fc[,"M W1"] > 0,])
whitew1upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew1upsig,]),whiteactivepeaks)
whitew1upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew1upsig & peakanno$trimanno %in% "Promoter",]),whiteactivepeaks)
whitew1upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew1upsig & peakanno$trimanno %in% "Intron",]),whiteactivepeaks)
whitew1upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew1upsig & peakanno$trimanno %in% "Distal Intergenic",]),whiteactivepeaks)
whitew1upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew1upsig & peakanno$trimanno %in% "Upstream",]),whiteactivepeaks)

whitew1downsig <- rownames(whitesigrnal2fc[whitesigrnal2fc[,"F W1"] < 0 & whitesigrnal2fc[,"M W1"] < 0,])
whitew1downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew1downsig,]),whiteactivepeaks)
whitew1downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew1downsig & peakanno$trimanno %in% "Promoter",]),whiteactivepeaks)
whitew1downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew1downsig & peakanno$trimanno %in% "Intron",]),whiteactivepeaks)
whitew1downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew1downsig & peakanno$trimanno %in% "Distal Intergenic",]),whiteactivepeaks)
whitew1downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew1downsig & peakanno$trimanno %in% "Upstream",]),whiteactivepeaks)


whitew2upsig <- rownames(whitesigrnal2fc[whitesigrnal2fc[,"F W2"] > 0 & whitesigrnal2fc[,"M W1"] > 0,])
whitew2upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew2upsig,]),whiteactivepeaks)
whitew2upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew2upsig & peakanno$trimanno %in% "Promoter",]),whiteactivepeaks)
whitew2upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew2upsig & peakanno$trimanno %in% "Intron",]),whiteactivepeaks)
whitew2upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew2upsig & peakanno$trimanno %in% "Distal Intergenic",]),whiteactivepeaks)
whitew2upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew2upsig & peakanno$trimanno %in% "Upstream",]),whiteactivepeaks)

whitew2downsig <- rownames(whitesigrnal2fc[whitesigrnal2fc[,"F W2"] < 0 & whitesigrnal2fc[,"M W1"] < 0,])
whitew2downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew2downsig,]),whiteactivepeaks)
whitew2downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew2downsig & peakanno$trimanno %in% "Promoter",]),whiteactivepeaks)
whitew2downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew2downsig & peakanno$trimanno %in% "Intron",]),whiteactivepeaks)
whitew2downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew2downsig & peakanno$trimanno %in% "Distal Intergenic",]),whiteactivepeaks)
whitew2downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew2downsig & peakanno$trimanno %in% "Upstream",]),whiteactivepeaks)


whitew4upsig <- rownames(whitesigrnal2fc[whitesigrnal2fc[,"F W4"] > 0 & whitesigrnal2fc[,"M W4"] > 0,])
whitew4upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew4upsig,]),whiteactivepeaks)
whitew4upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew4upsig & peakanno$trimanno %in% "Promoter",]),whiteactivepeaks)
whitew4upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew4upsig & peakanno$trimanno %in% "Intron",]),whiteactivepeaks)
whitew4upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew4upsig & peakanno$trimanno %in% "Distal Intergenic",]),whiteactivepeaks)
whitew4upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew4upsig & peakanno$trimanno %in% "Upstream",]),whiteactivepeaks)

whitew4downsig <- rownames(whitesigrnal2fc[whitesigrnal2fc[,"F W4"] < 0 & whitesigrnal2fc[,"M W4"] < 0,])
whitew4downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew4downsig,]),whiteactivepeaks)
whitew4downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew4downsig & peakanno$trimanno %in% "Promoter",]),whiteactivepeaks)
whitew4downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew4downsig & peakanno$trimanno %in% "Intron",]),whiteactivepeaks)
whitew4downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew4downsig & peakanno$trimanno %in% "Distal Intergenic",]),whiteactivepeaks)
whitew4downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew4downsig & peakanno$trimanno %in% "Upstream",]),whiteactivepeaks)


gastro50motifdf <- data.frame(row.names = unique(gastro50peakmotifs$Motif.Name),
                              "Enrichment in Active Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in Sig Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in Sig Prom Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in Sig Int Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in Sig Dist Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in W1 UpReg Sig Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in W1 UpReg Sig Prom Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in W1 UpReg Sig Int Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in W1 UpReg Sig Dist Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in W1 DownReg Sig Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in W1 DownReg Sig Prom Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in W1 DownReg Sig Int Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in W1 DownReg Sig Dist Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in W2 UpReg Sig Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in W2 UpReg Sig Prom Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in W2 UpReg Sig Int Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in W2 UpReg Sig Dist Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in W2 DownReg Sig Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in W2 DownReg Sig Prom Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in W2 DownReg Sig Int Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in W2 DownReg Sig Dist Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in W4 UpReg Sig Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in W4 UpReg Sig Prom Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in W4 UpReg Sig Int Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in W4 UpReg Sig Dist Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in W4 DownReg Sig Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in W4 DownReg Sig Prom Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in W4 DownReg Sig Int Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in W4 DownReg Sig Dist Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in W8 UpReg Sig Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in W8 UpReg Sig Prom Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in W8 UpReg Sig Int Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in W8 UpReg Sig Dist Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in W8 DownReg Sig Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in W8 DownReg Sig Prom Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in W8 DownReg Sig Int Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))),
                              "Enrichment in W8 DownReg Sig Dist Peaks" = rep(0,length(unique(gastro50peakmotifs$Motif.Name))))

for(i in 1:length(rownames(gastro50motifdf))){
  
  if(i%%20 == 0){
    print(i)
  }
  
  ourmotif <- rownames(gastro50motifdf)[i]
  ourmotifdf <- gastro50peakmotifs[gastro50peakmotifs$Motif.Name %in% ourmotif,]
  
  gastro50motifdf[i,"Enrichment.in.Active.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastroactivepeaks,])[1]/length(gastroactivepeaks)
  gastro50motifdf[i,"Enrichment.in.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrosigpeak,])[1]/length(gastrosigpeak)
  gastro50motifdf[i,"Enrichment.in.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrosigprompeak,])[1]/length(gastrosigprompeak)
  gastro50motifdf[i,"Enrichment.in.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrosigintpeak,])[1]/length(gastrosigintpeak)
  gastro50motifdf[i,"Enrichment.in.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrosigdistpeak,])[1]/length(gastrosigdistpeak)
  gastro50motifdf[i,"Enrichment.in.W1.UpReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrow1upsigpeak,])[1]/length(gastrow1upsigpeak)
  gastro50motifdf[i,"Enrichment.in.W1.UpReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrow1upsigprompeak,])[1]/length(gastrow1upsigprompeak)
  gastro50motifdf[i,"Enrichment.in.W1.UpReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrow1upsigintpeak,])[1]/length(gastrow1upsigintpeak)
  gastro50motifdf[i,"Enrichment.in.W1.UpReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrow1upsigdistpeak,])[1]/length(gastrow1upsigdistpeak)
  gastro50motifdf[i,"Enrichment.in.W1.DownReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrow1downsigpeak,])[1]/length(gastrow1downsigpeak)
  gastro50motifdf[i,"Enrichment.in.W1.DownReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrow1downsigprompeak,])[1]/length(gastrow1downsigprompeak)
  gastro50motifdf[i,"Enrichment.in.W1.DownReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrow1downsigintpeak,])[1]/length(gastrow1downsigintpeak)
  gastro50motifdf[i,"Enrichment.in.W1.DownReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrow1downsigdistpeak,])[1]/length(gastrow1downsigdistpeak)
  gastro50motifdf[i,"Enrichment.in.W2.UpReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrow2upsigpeak,])[1]/length(gastrow2upsigpeak)
  gastro50motifdf[i,"Enrichment.in.W2.UpReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrow2upsigprompeak,])[1]/length(gastrow2upsigprompeak)
  gastro50motifdf[i,"Enrichment.in.W2.UpReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrow2upsigintpeak,])[1]/length(gastrow2upsigintpeak)
  gastro50motifdf[i,"Enrichment.in.W2.UpReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrow2upsigdistpeak,])[1]/length(gastrow2upsigdistpeak)
  gastro50motifdf[i,"Enrichment.in.W2.DownReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrow2downsigpeak,])[1]/length(gastrow2downsigpeak)
  gastro50motifdf[i,"Enrichment.in.W2.DownReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrow2downsigprompeak,])[1]/length(gastrow2downsigprompeak)
  gastro50motifdf[i,"Enrichment.in.W2.DownReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrow2downsigintpeak,])[1]/length(gastrow2downsigintpeak)
  gastro50motifdf[i,"Enrichment.in.W2.DownReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrow2downsigdistpeak,])[1]/length(gastrow2downsigdistpeak)
  gastro50motifdf[i,"Enrichment.in.W4.UpReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrow4upsigpeak,])[1]/length(gastrow4upsigpeak)
  gastro50motifdf[i,"Enrichment.in.W4.UpReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrow4upsigprompeak,])[1]/length(gastrow4upsigprompeak)
  gastro50motifdf[i,"Enrichment.in.W4.UpReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrow4upsigintpeak,])[1]/length(gastrow4upsigintpeak)
  gastro50motifdf[i,"Enrichment.in.W4.UpReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrow4upsigdistpeak,])[1]/length(gastrow4upsigdistpeak)
  gastro50motifdf[i,"Enrichment.in.W4.DownReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrow4downsigpeak,])[1]/length(gastrow4downsigpeak)
  gastro50motifdf[i,"Enrichment.in.W4.DownReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrow4downsigprompeak,])[1]/length(gastrow4downsigprompeak)
  gastro50motifdf[i,"Enrichment.in.W4.DownReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrow4downsigintpeak,])[1]/length(gastrow4downsigintpeak)
  gastro50motifdf[i,"Enrichment.in.W4.DownReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrow4downsigdistpeak,])[1]/length(gastrow4downsigdistpeak)
  gastro50motifdf[i,"Enrichment.in.W8.UpReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrow8upsigpeak,])[1]/length(gastrow8upsigpeak)
  gastro50motifdf[i,"Enrichment.in.W8.UpReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrow8upsigprompeak,])[1]/length(gastrow8upsigprompeak)
  gastro50motifdf[i,"Enrichment.in.W8.UpReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrow8upsigintpeak,])[1]/length(gastrow8upsigintpeak)
  gastro50motifdf[i,"Enrichment.in.W8.UpReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrow8upsigdistpeak,])[1]/length(gastrow8upsigdistpeak)
  gastro50motifdf[i,"Enrichment.in.W8.DownReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrow8downsigpeak,])[1]/length(gastrow8downsigpeak)
  gastro50motifdf[i,"Enrichment.in.W8.DownReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrow8downsigprompeak,])[1]/length(gastrow8downsigprompeak)
  gastro50motifdf[i,"Enrichment.in.W8.DownReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrow8downsigintpeak,])[1]/length(gastrow8downsigintpeak)
  gastro50motifdf[i,"Enrichment.in.W8.DownReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% gastrow8downsigdistpeak,])[1]/length(gastrow8downsigdistpeak)
  
}



heart50motifdf <- data.frame(row.names = unique(heart50peakmotifs$Motif.Name),
                             "Enrichment in Active Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in Sig Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in Sig Prom Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in Sig Int Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in Sig Dist Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in W1 UpReg Sig Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in W1 UpReg Sig Prom Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in W1 UpReg Sig Int Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in W1 UpReg Sig Dist Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in W1 DownReg Sig Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in W1 DownReg Sig Prom Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in W1 DownReg Sig Int Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in W1 DownReg Sig Dist Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in W2 UpReg Sig Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in W2 UpReg Sig Prom Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in W2 UpReg Sig Int Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in W2 UpReg Sig Dist Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in W2 DownReg Sig Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in W2 DownReg Sig Prom Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in W2 DownReg Sig Int Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in W2 DownReg Sig Dist Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in W4 UpReg Sig Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in W4 UpReg Sig Prom Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in W4 UpReg Sig Int Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in W4 UpReg Sig Dist Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in W4 DownReg Sig Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in W4 DownReg Sig Prom Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in W4 DownReg Sig Int Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in W4 DownReg Sig Dist Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in W8 UpReg Sig Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in W8 UpReg Sig Prom Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in W8 UpReg Sig Int Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in W8 UpReg Sig Dist Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in W8 DownReg Sig Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in W8 DownReg Sig Prom Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in W8 DownReg Sig Int Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))),
                             "Enrichment in W8 DownReg Sig Dist Peaks" = rep(0,length(unique(heart50peakmotifs$Motif.Name))))

for(i in 1:length(rownames(heart50motifdf))){
  
  if(i%%20 == 0){
    print(i)
  }
  
  ourmotif <- rownames(heart50motifdf)[i]
  ourmotifdf <- heart50peakmotifs[heart50peakmotifs$Motif.Name %in% ourmotif,]
  
  heart50motifdf[i,"Enrichment.in.Active.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartactivepeaks,])[1]/length(heartactivepeaks)
  heart50motifdf[i,"Enrichment.in.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartsigpeak,])[1]/length(heartsigpeak)
  heart50motifdf[i,"Enrichment.in.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartsigprompeak,])[1]/length(heartsigprompeak)
  heart50motifdf[i,"Enrichment.in.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartsigintpeak,])[1]/length(heartsigintpeak)
  heart50motifdf[i,"Enrichment.in.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartsigdistpeak,])[1]/length(heartsigdistpeak)
  heart50motifdf[i,"Enrichment.in.W1.UpReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartw1upsigpeak,])[1]/length(heartw1upsigpeak)
  heart50motifdf[i,"Enrichment.in.W1.UpReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartw1upsigprompeak,])[1]/length(heartw1upsigprompeak)
  heart50motifdf[i,"Enrichment.in.W1.UpReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartw1upsigintpeak,])[1]/length(heartw1upsigintpeak)
  heart50motifdf[i,"Enrichment.in.W1.UpReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartw1upsigdistpeak,])[1]/length(heartw1upsigdistpeak)
  heart50motifdf[i,"Enrichment.in.W1.DownReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartw1downsigpeak,])[1]/length(heartw1downsigpeak)
  heart50motifdf[i,"Enrichment.in.W1.DownReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartw1downsigprompeak,])[1]/length(heartw1downsigprompeak)
  heart50motifdf[i,"Enrichment.in.W1.DownReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartw1downsigintpeak,])[1]/length(heartw1downsigintpeak)
  heart50motifdf[i,"Enrichment.in.W1.DownReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartw1downsigdistpeak,])[1]/length(heartw1downsigdistpeak)
  heart50motifdf[i,"Enrichment.in.W2.UpReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartw2upsigpeak,])[1]/length(heartw2upsigpeak)
  heart50motifdf[i,"Enrichment.in.W2.UpReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartw2upsigprompeak,])[1]/length(heartw2upsigprompeak)
  heart50motifdf[i,"Enrichment.in.W2.UpReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartw2upsigintpeak,])[1]/length(heartw2upsigintpeak)
  heart50motifdf[i,"Enrichment.in.W2.UpReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartw2upsigdistpeak,])[1]/length(heartw2upsigdistpeak)
  heart50motifdf[i,"Enrichment.in.W2.DownReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartw2downsigpeak,])[1]/length(heartw2downsigpeak)
  heart50motifdf[i,"Enrichment.in.W2.DownReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartw2downsigprompeak,])[1]/length(heartw2downsigprompeak)
  heart50motifdf[i,"Enrichment.in.W2.DownReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartw2downsigintpeak,])[1]/length(heartw2downsigintpeak)
  heart50motifdf[i,"Enrichment.in.W2.DownReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartw2downsigdistpeak,])[1]/length(heartw2downsigdistpeak)
  heart50motifdf[i,"Enrichment.in.W4.UpReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartw4upsigpeak,])[1]/length(heartw4upsigpeak)
  heart50motifdf[i,"Enrichment.in.W4.UpReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartw4upsigprompeak,])[1]/length(heartw4upsigprompeak)
  heart50motifdf[i,"Enrichment.in.W4.UpReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartw4upsigintpeak,])[1]/length(heartw4upsigintpeak)
  heart50motifdf[i,"Enrichment.in.W4.UpReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartw4upsigdistpeak,])[1]/length(heartw4upsigdistpeak)
  heart50motifdf[i,"Enrichment.in.W4.DownReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartw4downsigpeak,])[1]/length(heartw4downsigpeak)
  heart50motifdf[i,"Enrichment.in.W4.DownReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartw4downsigprompeak,])[1]/length(heartw4downsigprompeak)
  heart50motifdf[i,"Enrichment.in.W4.DownReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartw4downsigintpeak,])[1]/length(heartw4downsigintpeak)
  heart50motifdf[i,"Enrichment.in.W4.DownReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartw4downsigdistpeak,])[1]/length(heartw4downsigdistpeak)
  heart50motifdf[i,"Enrichment.in.W8.UpReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartw8upsigpeak,])[1]/length(heartw8upsigpeak)
  heart50motifdf[i,"Enrichment.in.W8.UpReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartw8upsigprompeak,])[1]/length(heartw8upsigprompeak)
  heart50motifdf[i,"Enrichment.in.W8.UpReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartw8upsigintpeak,])[1]/length(heartw8upsigintpeak)
  heart50motifdf[i,"Enrichment.in.W8.UpReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartw8upsigdistpeak,])[1]/length(heartw8upsigdistpeak)
  heart50motifdf[i,"Enrichment.in.W8.DownReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartw8downsigpeak,])[1]/length(heartw8downsigpeak)
  heart50motifdf[i,"Enrichment.in.W8.DownReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartw8downsigprompeak,])[1]/length(heartw8downsigprompeak)
  heart50motifdf[i,"Enrichment.in.W8.DownReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartw8downsigintpeak,])[1]/length(heartw8downsigintpeak)
  heart50motifdf[i,"Enrichment.in.W8.DownReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% heartw8downsigdistpeak,])[1]/length(heartw8downsigdistpeak)
  
}


hippo50motifdf <- data.frame(row.names = unique(hippo50peakmotifs$Motif.Name),
                             "Enrichment in Active Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in Sig Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in Sig Prom Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in Sig Int Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in Sig Dist Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in W1 UpReg Sig Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in W1 UpReg Sig Prom Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in W1 UpReg Sig Int Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in W1 UpReg Sig Dist Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in W1 DownReg Sig Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in W1 DownReg Sig Prom Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in W1 DownReg Sig Int Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in W1 DownReg Sig Dist Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in W2 UpReg Sig Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in W2 UpReg Sig Prom Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in W2 UpReg Sig Int Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in W2 UpReg Sig Dist Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in W2 DownReg Sig Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in W2 DownReg Sig Prom Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in W2 DownReg Sig Int Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in W2 DownReg Sig Dist Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in W4 UpReg Sig Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in W4 UpReg Sig Prom Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in W4 UpReg Sig Int Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in W4 UpReg Sig Dist Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in W4 DownReg Sig Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in W4 DownReg Sig Prom Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in W4 DownReg Sig Int Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in W4 DownReg Sig Dist Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in W8 UpReg Sig Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in W8 UpReg Sig Prom Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in W8 UpReg Sig Int Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in W8 UpReg Sig Dist Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in W8 DownReg Sig Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in W8 DownReg Sig Prom Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in W8 DownReg Sig Int Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))),
                             "Enrichment in W8 DownReg Sig Dist Peaks" = rep(0,length(unique(hippo50peakmotifs$Motif.Name))))

for(i in 1:length(rownames(hippo50motifdf))){
  
  if(i%%20 == 0){
    print(i)
  }
  
  ourmotif <- rownames(hippo50motifdf)[i]
  ourmotifdf <- hippo50peakmotifs[hippo50peakmotifs$Motif.Name %in% ourmotif,]
  
  hippo50motifdf[i,"Enrichment.in.Active.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hippoactivepeaks,])[1]/length(hippoactivepeaks)
  hippo50motifdf[i,"Enrichment.in.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hipposigpeak,])[1]/length(hipposigpeak)
  hippo50motifdf[i,"Enrichment.in.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hipposigprompeak,])[1]/length(hipposigprompeak)
  hippo50motifdf[i,"Enrichment.in.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hipposigintpeak,])[1]/length(hipposigintpeak)
  hippo50motifdf[i,"Enrichment.in.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hipposigdistpeak,])[1]/length(hipposigdistpeak)
  hippo50motifdf[i,"Enrichment.in.W1.UpReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hippow1upsigpeak,])[1]/length(hippow1upsigpeak)
  hippo50motifdf[i,"Enrichment.in.W1.UpReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hippow1upsigprompeak,])[1]/length(hippow1upsigprompeak)
  hippo50motifdf[i,"Enrichment.in.W1.UpReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hippow1upsigintpeak,])[1]/length(hippow1upsigintpeak)
  hippo50motifdf[i,"Enrichment.in.W1.UpReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hippow1upsigdistpeak,])[1]/length(hippow1upsigdistpeak)
  hippo50motifdf[i,"Enrichment.in.W1.DownReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hippow1downsigpeak,])[1]/length(hippow1downsigpeak)
  hippo50motifdf[i,"Enrichment.in.W1.DownReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hippow1downsigprompeak,])[1]/length(hippow1downsigprompeak)
  hippo50motifdf[i,"Enrichment.in.W1.DownReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hippow1downsigintpeak,])[1]/length(hippow1downsigintpeak)
  hippo50motifdf[i,"Enrichment.in.W1.DownReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hippow1downsigdistpeak,])[1]/length(hippow1downsigdistpeak)
  hippo50motifdf[i,"Enrichment.in.W2.UpReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hippow2upsigpeak,])[1]/length(hippow2upsigpeak)
  hippo50motifdf[i,"Enrichment.in.W2.UpReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hippow2upsigprompeak,])[1]/length(hippow2upsigprompeak)
  hippo50motifdf[i,"Enrichment.in.W2.UpReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hippow2upsigintpeak,])[1]/length(hippow2upsigintpeak)
  hippo50motifdf[i,"Enrichment.in.W2.UpReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hippow2upsigdistpeak,])[1]/length(hippow2upsigdistpeak)
  hippo50motifdf[i,"Enrichment.in.W2.DownReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hippow2downsigpeak,])[1]/length(hippow2downsigpeak)
  hippo50motifdf[i,"Enrichment.in.W2.DownReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hippow2downsigprompeak,])[1]/length(hippow2downsigprompeak)
  hippo50motifdf[i,"Enrichment.in.W2.DownReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hippow2downsigintpeak,])[1]/length(hippow2downsigintpeak)
  hippo50motifdf[i,"Enrichment.in.W2.DownReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hippow2downsigdistpeak,])[1]/length(hippow2downsigdistpeak)
  hippo50motifdf[i,"Enrichment.in.W4.UpReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hippow4upsigpeak,])[1]/length(hippow4upsigpeak)
  hippo50motifdf[i,"Enrichment.in.W4.UpReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hippow4upsigprompeak,])[1]/length(hippow4upsigprompeak)
  hippo50motifdf[i,"Enrichment.in.W4.UpReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hippow4upsigintpeak,])[1]/length(hippow4upsigintpeak)
  hippo50motifdf[i,"Enrichment.in.W4.UpReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hippow4upsigdistpeak,])[1]/length(hippow4upsigdistpeak)
  hippo50motifdf[i,"Enrichment.in.W4.DownReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hippow4downsigpeak,])[1]/length(hippow4downsigpeak)
  hippo50motifdf[i,"Enrichment.in.W4.DownReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hippow4downsigprompeak,])[1]/length(hippow4downsigprompeak)
  hippo50motifdf[i,"Enrichment.in.W4.DownReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hippow4downsigintpeak,])[1]/length(hippow4downsigintpeak)
  hippo50motifdf[i,"Enrichment.in.W4.DownReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hippow4downsigdistpeak,])[1]/length(hippow4downsigdistpeak)
  hippo50motifdf[i,"Enrichment.in.W8.UpReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hippow8upsigpeak,])[1]/length(hippow8upsigpeak)
  hippo50motifdf[i,"Enrichment.in.W8.UpReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hippow8upsigprompeak,])[1]/length(hippow8upsigprompeak)
  hippo50motifdf[i,"Enrichment.in.W8.UpReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hippow8upsigintpeak,])[1]/length(hippow8upsigintpeak)
  hippo50motifdf[i,"Enrichment.in.W8.UpReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hippow8upsigdistpeak,])[1]/length(hippow8upsigdistpeak)
  hippo50motifdf[i,"Enrichment.in.W8.DownReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hippow8downsigpeak,])[1]/length(hippow8downsigpeak)
  hippo50motifdf[i,"Enrichment.in.W8.DownReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hippow8downsigprompeak,])[1]/length(hippow8downsigprompeak)
  hippo50motifdf[i,"Enrichment.in.W8.DownReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hippow8downsigintpeak,])[1]/length(hippow8downsigintpeak)
  hippo50motifdf[i,"Enrichment.in.W8.DownReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% hippow8downsigdistpeak,])[1]/length(hippow8downsigdistpeak)
  
}

kidney50motifdf <- data.frame(row.names = unique(kidney50peakmotifs$Motif.Name),
                              "Enrichment in Active Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in Sig Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in Sig Prom Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in Sig Int Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in Sig Dist Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in W1 UpReg Sig Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in W1 UpReg Sig Prom Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in W1 UpReg Sig Int Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in W1 UpReg Sig Dist Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in W1 DownReg Sig Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in W1 DownReg Sig Prom Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in W1 DownReg Sig Int Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in W1 DownReg Sig Dist Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in W2 UpReg Sig Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in W2 UpReg Sig Prom Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in W2 UpReg Sig Int Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in W2 UpReg Sig Dist Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in W2 DownReg Sig Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in W2 DownReg Sig Prom Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in W2 DownReg Sig Int Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in W2 DownReg Sig Dist Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in W4 UpReg Sig Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in W4 UpReg Sig Prom Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in W4 UpReg Sig Int Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in W4 UpReg Sig Dist Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in W4 DownReg Sig Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in W4 DownReg Sig Prom Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in W4 DownReg Sig Int Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in W4 DownReg Sig Dist Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in W8 UpReg Sig Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in W8 UpReg Sig Prom Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in W8 UpReg Sig Int Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in W8 UpReg Sig Dist Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in W8 DownReg Sig Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in W8 DownReg Sig Prom Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in W8 DownReg Sig Int Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))),
                              "Enrichment in W8 DownReg Sig Dist Peaks" = rep(0,length(unique(kidney50peakmotifs$Motif.Name))))

for(i in 1:length(rownames(kidney50motifdf))){
  
  if(i%%20 == 0){
    print(i)
  }
  
  ourmotif <- rownames(kidney50motifdf)[i]
  ourmotifdf <- kidney50peakmotifs[kidney50peakmotifs$Motif.Name %in% ourmotif,]
  
  kidney50motifdf[i,"Enrichment.in.Active.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneyactivepeaks,])[1]/length(kidneyactivepeaks)
  kidney50motifdf[i,"Enrichment.in.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneysigpeak,])[1]/length(kidneysigpeak)
  kidney50motifdf[i,"Enrichment.in.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneysigprompeak,])[1]/length(kidneysigprompeak)
  kidney50motifdf[i,"Enrichment.in.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneysigintpeak,])[1]/length(kidneysigintpeak)
  kidney50motifdf[i,"Enrichment.in.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneysigdistpeak,])[1]/length(kidneysigdistpeak)
  kidney50motifdf[i,"Enrichment.in.W1.UpReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneyw1upsigpeak,])[1]/length(kidneyw1upsigpeak)
  kidney50motifdf[i,"Enrichment.in.W1.UpReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneyw1upsigprompeak,])[1]/length(kidneyw1upsigprompeak)
  kidney50motifdf[i,"Enrichment.in.W1.UpReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneyw1upsigintpeak,])[1]/length(kidneyw1upsigintpeak)
  kidney50motifdf[i,"Enrichment.in.W1.UpReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneyw1upsigdistpeak,])[1]/length(kidneyw1upsigdistpeak)
  kidney50motifdf[i,"Enrichment.in.W1.DownReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneyw1downsigpeak,])[1]/length(kidneyw1downsigpeak)
  kidney50motifdf[i,"Enrichment.in.W1.DownReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneyw1downsigprompeak,])[1]/length(kidneyw1downsigprompeak)
  kidney50motifdf[i,"Enrichment.in.W1.DownReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneyw1downsigintpeak,])[1]/length(kidneyw1downsigintpeak)
  kidney50motifdf[i,"Enrichment.in.W1.DownReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneyw1downsigdistpeak,])[1]/length(kidneyw1downsigdistpeak)
  kidney50motifdf[i,"Enrichment.in.W2.UpReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneyw2upsigpeak,])[1]/length(kidneyw2upsigpeak)
  kidney50motifdf[i,"Enrichment.in.W2.UpReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneyw2upsigprompeak,])[1]/length(kidneyw2upsigprompeak)
  kidney50motifdf[i,"Enrichment.in.W2.UpReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneyw2upsigintpeak,])[1]/length(kidneyw2upsigintpeak)
  kidney50motifdf[i,"Enrichment.in.W2.UpReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneyw2upsigdistpeak,])[1]/length(kidneyw2upsigdistpeak)
  kidney50motifdf[i,"Enrichment.in.W2.DownReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneyw2downsigpeak,])[1]/length(kidneyw2downsigpeak)
  kidney50motifdf[i,"Enrichment.in.W2.DownReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneyw2downsigprompeak,])[1]/length(kidneyw2downsigprompeak)
  kidney50motifdf[i,"Enrichment.in.W2.DownReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneyw2downsigintpeak,])[1]/length(kidneyw2downsigintpeak)
  kidney50motifdf[i,"Enrichment.in.W2.DownReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneyw2downsigdistpeak,])[1]/length(kidneyw2downsigdistpeak)
  kidney50motifdf[i,"Enrichment.in.W4.UpReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneyw4upsigpeak,])[1]/length(kidneyw4upsigpeak)
  kidney50motifdf[i,"Enrichment.in.W4.UpReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneyw4upsigprompeak,])[1]/length(kidneyw4upsigprompeak)
  kidney50motifdf[i,"Enrichment.in.W4.UpReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneyw4upsigintpeak,])[1]/length(kidneyw4upsigintpeak)
  kidney50motifdf[i,"Enrichment.in.W4.UpReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneyw4upsigdistpeak,])[1]/length(kidneyw4upsigdistpeak)
  kidney50motifdf[i,"Enrichment.in.W4.DownReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneyw4downsigpeak,])[1]/length(kidneyw4downsigpeak)
  kidney50motifdf[i,"Enrichment.in.W4.DownReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneyw4downsigprompeak,])[1]/length(kidneyw4downsigprompeak)
  kidney50motifdf[i,"Enrichment.in.W4.DownReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneyw4downsigintpeak,])[1]/length(kidneyw4downsigintpeak)
  kidney50motifdf[i,"Enrichment.in.W4.DownReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneyw4downsigdistpeak,])[1]/length(kidneyw4downsigdistpeak)
  kidney50motifdf[i,"Enrichment.in.W8.UpReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneyw8upsigpeak,])[1]/length(kidneyw8upsigpeak)
  kidney50motifdf[i,"Enrichment.in.W8.UpReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneyw8upsigprompeak,])[1]/length(kidneyw8upsigprompeak)
  kidney50motifdf[i,"Enrichment.in.W8.UpReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneyw8upsigintpeak,])[1]/length(kidneyw8upsigintpeak)
  kidney50motifdf[i,"Enrichment.in.W8.UpReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneyw8upsigdistpeak,])[1]/length(kidneyw8upsigdistpeak)
  kidney50motifdf[i,"Enrichment.in.W8.DownReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneyw8downsigpeak,])[1]/length(kidneyw8downsigpeak)
  kidney50motifdf[i,"Enrichment.in.W8.DownReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneyw8downsigprompeak,])[1]/length(kidneyw8downsigprompeak)
  kidney50motifdf[i,"Enrichment.in.W8.DownReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneyw8downsigintpeak,])[1]/length(kidneyw8downsigintpeak)
  kidney50motifdf[i,"Enrichment.in.W8.DownReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% kidneyw8downsigdistpeak,])[1]/length(kidneyw8downsigdistpeak)
  
}


liver50motifdf <- data.frame(row.names = unique(liver50peakmotifs$Motif.Name),
                             "Enrichment in Active Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in Sig Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in Sig Prom Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in Sig Int Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in Sig Dist Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in W1 UpReg Sig Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in W1 UpReg Sig Prom Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in W1 UpReg Sig Int Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in W1 UpReg Sig Dist Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in W1 DownReg Sig Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in W1 DownReg Sig Prom Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in W1 DownReg Sig Int Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in W1 DownReg Sig Dist Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in W2 UpReg Sig Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in W2 UpReg Sig Prom Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in W2 UpReg Sig Int Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in W2 UpReg Sig Dist Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in W2 DownReg Sig Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in W2 DownReg Sig Prom Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in W2 DownReg Sig Int Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in W2 DownReg Sig Dist Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in W4 UpReg Sig Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in W4 UpReg Sig Prom Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in W4 UpReg Sig Int Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in W4 UpReg Sig Dist Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in W4 DownReg Sig Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in W4 DownReg Sig Prom Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in W4 DownReg Sig Int Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in W4 DownReg Sig Dist Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in W8 UpReg Sig Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in W8 UpReg Sig Prom Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in W8 UpReg Sig Int Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in W8 UpReg Sig Dist Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in W8 DownReg Sig Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in W8 DownReg Sig Prom Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in W8 DownReg Sig Int Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))),
                             "Enrichment in W8 DownReg Sig Dist Peaks" = rep(0,length(unique(liver50peakmotifs$Motif.Name))))

for(i in 1:length(rownames(liver50motifdf))){
  
  if(i%%20 == 0){
    print(i)
  }
  
  ourmotif <- rownames(liver50motifdf)[i]
  ourmotifdf <- liver50peakmotifs[liver50peakmotifs$Motif.Name %in% ourmotif,]
  
  liver50motifdf[i,"Enrichment.in.Active.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liveractivepeaks,])[1]/length(liveractivepeaks)
  liver50motifdf[i,"Enrichment.in.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liversigpeak,])[1]/length(liversigpeak)
  liver50motifdf[i,"Enrichment.in.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liversigprompeak,])[1]/length(liversigprompeak)
  liver50motifdf[i,"Enrichment.in.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liversigintpeak,])[1]/length(liversigintpeak)
  liver50motifdf[i,"Enrichment.in.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liversigdistpeak,])[1]/length(liversigdistpeak)
  liver50motifdf[i,"Enrichment.in.W1.UpReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liverw1upsigpeak,])[1]/length(liverw1upsigpeak)
  liver50motifdf[i,"Enrichment.in.W1.UpReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liverw1upsigprompeak,])[1]/length(liverw1upsigprompeak)
  liver50motifdf[i,"Enrichment.in.W1.UpReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liverw1upsigintpeak,])[1]/length(liverw1upsigintpeak)
  liver50motifdf[i,"Enrichment.in.W1.UpReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liverw1upsigdistpeak,])[1]/length(liverw1upsigdistpeak)
  liver50motifdf[i,"Enrichment.in.W1.DownReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liverw1downsigpeak,])[1]/length(liverw1downsigpeak)
  liver50motifdf[i,"Enrichment.in.W1.DownReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liverw1downsigprompeak,])[1]/length(liverw1downsigprompeak)
  liver50motifdf[i,"Enrichment.in.W1.DownReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liverw1downsigintpeak,])[1]/length(liverw1downsigintpeak)
  liver50motifdf[i,"Enrichment.in.W1.DownReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liverw1downsigdistpeak,])[1]/length(liverw1downsigdistpeak)
  liver50motifdf[i,"Enrichment.in.W2.UpReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liverw2upsigpeak,])[1]/length(liverw2upsigpeak)
  liver50motifdf[i,"Enrichment.in.W2.UpReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liverw2upsigprompeak,])[1]/length(liverw2upsigprompeak)
  liver50motifdf[i,"Enrichment.in.W2.UpReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liverw2upsigintpeak,])[1]/length(liverw2upsigintpeak)
  liver50motifdf[i,"Enrichment.in.W2.UpReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liverw2upsigdistpeak,])[1]/length(liverw2upsigdistpeak)
  liver50motifdf[i,"Enrichment.in.W2.DownReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liverw2downsigpeak,])[1]/length(liverw2downsigpeak)
  liver50motifdf[i,"Enrichment.in.W2.DownReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liverw2downsigprompeak,])[1]/length(liverw2downsigprompeak)
  liver50motifdf[i,"Enrichment.in.W2.DownReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liverw2downsigintpeak,])[1]/length(liverw2downsigintpeak)
  liver50motifdf[i,"Enrichment.in.W2.DownReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liverw2downsigdistpeak,])[1]/length(liverw2downsigdistpeak)
  liver50motifdf[i,"Enrichment.in.W4.UpReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liverw4upsigpeak,])[1]/length(liverw4upsigpeak)
  liver50motifdf[i,"Enrichment.in.W4.UpReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liverw4upsigprompeak,])[1]/length(liverw4upsigprompeak)
  liver50motifdf[i,"Enrichment.in.W4.UpReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liverw4upsigintpeak,])[1]/length(liverw4upsigintpeak)
  liver50motifdf[i,"Enrichment.in.W4.UpReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liverw4upsigdistpeak,])[1]/length(liverw4upsigdistpeak)
  liver50motifdf[i,"Enrichment.in.W4.DownReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liverw4downsigpeak,])[1]/length(liverw4downsigpeak)
  liver50motifdf[i,"Enrichment.in.W4.DownReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liverw4downsigprompeak,])[1]/length(liverw4downsigprompeak)
  liver50motifdf[i,"Enrichment.in.W4.DownReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liverw4downsigintpeak,])[1]/length(liverw4downsigintpeak)
  liver50motifdf[i,"Enrichment.in.W4.DownReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liverw4downsigdistpeak,])[1]/length(liverw4downsigdistpeak)
  liver50motifdf[i,"Enrichment.in.W8.UpReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liverw8upsigpeak,])[1]/length(liverw8upsigpeak)
  liver50motifdf[i,"Enrichment.in.W8.UpReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liverw8upsigprompeak,])[1]/length(liverw8upsigprompeak)
  liver50motifdf[i,"Enrichment.in.W8.UpReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liverw8upsigintpeak,])[1]/length(liverw8upsigintpeak)
  liver50motifdf[i,"Enrichment.in.W8.UpReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liverw8upsigdistpeak,])[1]/length(liverw8upsigdistpeak)
  liver50motifdf[i,"Enrichment.in.W8.DownReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liverw8downsigpeak,])[1]/length(liverw8downsigpeak)
  liver50motifdf[i,"Enrichment.in.W8.DownReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liverw8downsigprompeak,])[1]/length(liverw8downsigprompeak)
  liver50motifdf[i,"Enrichment.in.W8.DownReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liverw8downsigintpeak,])[1]/length(liverw8downsigintpeak)
  liver50motifdf[i,"Enrichment.in.W8.DownReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% liverw8downsigdistpeak,])[1]/length(liverw8downsigdistpeak)
  
}


lung50motifdf <- data.frame(row.names = unique(lung50peakmotifs$Motif.Name),
                            "Enrichment in Active Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in Sig Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in Sig Prom Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in Sig Int Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in Sig Dist Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in W1 UpReg Sig Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in W1 UpReg Sig Prom Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in W1 UpReg Sig Int Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in W1 UpReg Sig Dist Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in W1 DownReg Sig Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in W1 DownReg Sig Prom Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in W1 DownReg Sig Int Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in W1 DownReg Sig Dist Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in W2 UpReg Sig Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in W2 UpReg Sig Prom Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in W2 UpReg Sig Int Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in W2 UpReg Sig Dist Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in W2 DownReg Sig Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in W2 DownReg Sig Prom Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in W2 DownReg Sig Int Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in W2 DownReg Sig Dist Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in W4 UpReg Sig Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in W4 UpReg Sig Prom Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in W4 UpReg Sig Int Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in W4 UpReg Sig Dist Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in W4 DownReg Sig Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in W4 DownReg Sig Prom Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in W4 DownReg Sig Int Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in W4 DownReg Sig Dist Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in W8 UpReg Sig Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in W8 UpReg Sig Prom Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in W8 UpReg Sig Int Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in W8 UpReg Sig Dist Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in W8 DownReg Sig Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in W8 DownReg Sig Prom Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in W8 DownReg Sig Int Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))),
                            "Enrichment in W8 DownReg Sig Dist Peaks" = rep(0,length(unique(lung50peakmotifs$Motif.Name))))

for(i in 1:length(rownames(lung50motifdf))){
  
  if(i%%20 == 0){
    print(i)
  }
  
  ourmotif <- rownames(lung50motifdf)[i]
  ourmotifdf <- lung50peakmotifs[lung50peakmotifs$Motif.Name %in% ourmotif,]
  
  lung50motifdf[i,"Enrichment.in.Active.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungactivepeaks,])[1]/length(lungactivepeaks)
  lung50motifdf[i,"Enrichment.in.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungsigpeak,])[1]/length(lungsigpeak)
  lung50motifdf[i,"Enrichment.in.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungsigprompeak,])[1]/length(lungsigprompeak)
  lung50motifdf[i,"Enrichment.in.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungsigintpeak,])[1]/length(lungsigintpeak)
  lung50motifdf[i,"Enrichment.in.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungsigdistpeak,])[1]/length(lungsigdistpeak)
  lung50motifdf[i,"Enrichment.in.W1.UpReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungw1upsigpeak,])[1]/length(lungw1upsigpeak)
  lung50motifdf[i,"Enrichment.in.W1.UpReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungw1upsigprompeak,])[1]/length(lungw1upsigprompeak)
  lung50motifdf[i,"Enrichment.in.W1.UpReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungw1upsigintpeak,])[1]/length(lungw1upsigintpeak)
  lung50motifdf[i,"Enrichment.in.W1.UpReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungw1upsigdistpeak,])[1]/length(lungw1upsigdistpeak)
  lung50motifdf[i,"Enrichment.in.W1.DownReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungw1downsigpeak,])[1]/length(lungw1downsigpeak)
  lung50motifdf[i,"Enrichment.in.W1.DownReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungw1downsigprompeak,])[1]/length(lungw1downsigprompeak)
  lung50motifdf[i,"Enrichment.in.W1.DownReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungw1downsigintpeak,])[1]/length(lungw1downsigintpeak)
  lung50motifdf[i,"Enrichment.in.W1.DownReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungw1downsigdistpeak,])[1]/length(lungw1downsigdistpeak)
  lung50motifdf[i,"Enrichment.in.W2.UpReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungw2upsigpeak,])[1]/length(lungw2upsigpeak)
  lung50motifdf[i,"Enrichment.in.W2.UpReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungw2upsigprompeak,])[1]/length(lungw2upsigprompeak)
  lung50motifdf[i,"Enrichment.in.W2.UpReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungw2upsigintpeak,])[1]/length(lungw2upsigintpeak)
  lung50motifdf[i,"Enrichment.in.W2.UpReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungw2upsigdistpeak,])[1]/length(lungw2upsigdistpeak)
  lung50motifdf[i,"Enrichment.in.W2.DownReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungw2downsigpeak,])[1]/length(lungw2downsigpeak)
  lung50motifdf[i,"Enrichment.in.W2.DownReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungw2downsigprompeak,])[1]/length(lungw2downsigprompeak)
  lung50motifdf[i,"Enrichment.in.W2.DownReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungw2downsigintpeak,])[1]/length(lungw2downsigintpeak)
  lung50motifdf[i,"Enrichment.in.W2.DownReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungw2downsigdistpeak,])[1]/length(lungw2downsigdistpeak)
  lung50motifdf[i,"Enrichment.in.W4.UpReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungw4upsigpeak,])[1]/length(lungw4upsigpeak)
  lung50motifdf[i,"Enrichment.in.W4.UpReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungw4upsigprompeak,])[1]/length(lungw4upsigprompeak)
  lung50motifdf[i,"Enrichment.in.W4.UpReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungw4upsigintpeak,])[1]/length(lungw4upsigintpeak)
  lung50motifdf[i,"Enrichment.in.W4.UpReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungw4upsigdistpeak,])[1]/length(lungw4upsigdistpeak)
  lung50motifdf[i,"Enrichment.in.W4.DownReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungw4downsigpeak,])[1]/length(lungw4downsigpeak)
  lung50motifdf[i,"Enrichment.in.W4.DownReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungw4downsigprompeak,])[1]/length(lungw4downsigprompeak)
  lung50motifdf[i,"Enrichment.in.W4.DownReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungw4downsigintpeak,])[1]/length(lungw4downsigintpeak)
  lung50motifdf[i,"Enrichment.in.W4.DownReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungw4downsigdistpeak,])[1]/length(lungw4downsigdistpeak)
  lung50motifdf[i,"Enrichment.in.W8.UpReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungw8upsigpeak,])[1]/length(lungw8upsigpeak)
  lung50motifdf[i,"Enrichment.in.W8.UpReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungw8upsigprompeak,])[1]/length(lungw8upsigprompeak)
  lung50motifdf[i,"Enrichment.in.W8.UpReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungw8upsigintpeak,])[1]/length(lungw8upsigintpeak)
  lung50motifdf[i,"Enrichment.in.W8.UpReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungw8upsigdistpeak,])[1]/length(lungw8upsigdistpeak)
  lung50motifdf[i,"Enrichment.in.W8.DownReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungw8downsigpeak,])[1]/length(lungw8downsigpeak)
  lung50motifdf[i,"Enrichment.in.W8.DownReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungw8downsigprompeak,])[1]/length(lungw8downsigprompeak)
  lung50motifdf[i,"Enrichment.in.W8.DownReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungw8downsigintpeak,])[1]/length(lungw8downsigintpeak)
  lung50motifdf[i,"Enrichment.in.W8.DownReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% lungw8downsigdistpeak,])[1]/length(lungw8downsigdistpeak)
  
}


brown50motifdf <- data.frame(row.names = unique(brown50peakmotifs$Motif.Name),
                             "Enrichment in Active Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in Sig Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in Sig Prom Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in Sig Int Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in Sig Dist Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in W1 UpReg Sig Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in W1 UpReg Sig Prom Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in W1 UpReg Sig Int Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in W1 UpReg Sig Dist Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in W1 DownReg Sig Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in W1 DownReg Sig Prom Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in W1 DownReg Sig Int Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in W1 DownReg Sig Dist Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in W2 UpReg Sig Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in W2 UpReg Sig Prom Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in W2 UpReg Sig Int Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in W2 UpReg Sig Dist Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in W2 DownReg Sig Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in W2 DownReg Sig Prom Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in W2 DownReg Sig Int Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in W2 DownReg Sig Dist Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in W4 UpReg Sig Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in W4 UpReg Sig Prom Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in W4 UpReg Sig Int Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in W4 UpReg Sig Dist Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in W4 DownReg Sig Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in W4 DownReg Sig Prom Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in W4 DownReg Sig Int Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in W4 DownReg Sig Dist Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in W8 UpReg Sig Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in W8 UpReg Sig Prom Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in W8 UpReg Sig Int Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in W8 UpReg Sig Dist Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in W8 DownReg Sig Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in W8 DownReg Sig Prom Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in W8 DownReg Sig Int Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))),
                             "Enrichment in W8 DownReg Sig Dist Peaks" = rep(0,length(unique(brown50peakmotifs$Motif.Name))))

for(i in 1:length(rownames(brown50motifdf))){
  
  if(i%%20 == 0){
    print(i)
  }
  
  ourmotif <- rownames(brown50motifdf)[i]
  ourmotifdf <- brown50peakmotifs[brown50peakmotifs$Motif.Name %in% ourmotif,]
  
  brown50motifdf[i,"Enrichment.in.Active.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownactivepeaks,])[1]/length(brownactivepeaks)
  brown50motifdf[i,"Enrichment.in.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownsigpeak,])[1]/length(brownsigpeak)
  brown50motifdf[i,"Enrichment.in.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownsigprompeak,])[1]/length(brownsigprompeak)
  brown50motifdf[i,"Enrichment.in.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownsigintpeak,])[1]/length(brownsigintpeak)
  brown50motifdf[i,"Enrichment.in.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownsigdistpeak,])[1]/length(brownsigdistpeak)
  brown50motifdf[i,"Enrichment.in.W1.UpReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownw1upsigpeak,])[1]/length(brownw1upsigpeak)
  brown50motifdf[i,"Enrichment.in.W1.UpReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownw1upsigprompeak,])[1]/length(brownw1upsigprompeak)
  brown50motifdf[i,"Enrichment.in.W1.UpReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownw1upsigintpeak,])[1]/length(brownw1upsigintpeak)
  brown50motifdf[i,"Enrichment.in.W1.UpReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownw1upsigdistpeak,])[1]/length(brownw1upsigdistpeak)
  brown50motifdf[i,"Enrichment.in.W1.DownReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownw1downsigpeak,])[1]/length(brownw1downsigpeak)
  brown50motifdf[i,"Enrichment.in.W1.DownReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownw1downsigprompeak,])[1]/length(brownw1downsigprompeak)
  brown50motifdf[i,"Enrichment.in.W1.DownReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownw1downsigintpeak,])[1]/length(brownw1downsigintpeak)
  brown50motifdf[i,"Enrichment.in.W1.DownReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownw1downsigdistpeak,])[1]/length(brownw1downsigdistpeak)
  brown50motifdf[i,"Enrichment.in.W2.UpReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownw2upsigpeak,])[1]/length(brownw2upsigpeak)
  brown50motifdf[i,"Enrichment.in.W2.UpReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownw2upsigprompeak,])[1]/length(brownw2upsigprompeak)
  brown50motifdf[i,"Enrichment.in.W2.UpReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownw2upsigintpeak,])[1]/length(brownw2upsigintpeak)
  brown50motifdf[i,"Enrichment.in.W2.UpReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownw2upsigdistpeak,])[1]/length(brownw2upsigdistpeak)
  brown50motifdf[i,"Enrichment.in.W2.DownReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownw2downsigpeak,])[1]/length(brownw2downsigpeak)
  brown50motifdf[i,"Enrichment.in.W2.DownReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownw2downsigprompeak,])[1]/length(brownw2downsigprompeak)
  brown50motifdf[i,"Enrichment.in.W2.DownReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownw2downsigintpeak,])[1]/length(brownw2downsigintpeak)
  brown50motifdf[i,"Enrichment.in.W2.DownReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownw2downsigdistpeak,])[1]/length(brownw2downsigdistpeak)
  brown50motifdf[i,"Enrichment.in.W4.UpReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownw4upsigpeak,])[1]/length(brownw4upsigpeak)
  brown50motifdf[i,"Enrichment.in.W4.UpReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownw4upsigprompeak,])[1]/length(brownw4upsigprompeak)
  brown50motifdf[i,"Enrichment.in.W4.UpReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownw4upsigintpeak,])[1]/length(brownw4upsigintpeak)
  brown50motifdf[i,"Enrichment.in.W4.UpReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownw4upsigdistpeak,])[1]/length(brownw4upsigdistpeak)
  brown50motifdf[i,"Enrichment.in.W4.DownReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownw4downsigpeak,])[1]/length(brownw4downsigpeak)
  brown50motifdf[i,"Enrichment.in.W4.DownReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownw4downsigprompeak,])[1]/length(brownw4downsigprompeak)
  brown50motifdf[i,"Enrichment.in.W4.DownReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownw4downsigintpeak,])[1]/length(brownw4downsigintpeak)
  brown50motifdf[i,"Enrichment.in.W4.DownReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownw4downsigdistpeak,])[1]/length(brownw4downsigdistpeak)
  brown50motifdf[i,"Enrichment.in.W8.UpReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownw8upsigpeak,])[1]/length(brownw8upsigpeak)
  brown50motifdf[i,"Enrichment.in.W8.UpReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownw8upsigprompeak,])[1]/length(brownw8upsigprompeak)
  brown50motifdf[i,"Enrichment.in.W8.UpReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownw8upsigintpeak,])[1]/length(brownw8upsigintpeak)
  brown50motifdf[i,"Enrichment.in.W8.UpReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownw8upsigdistpeak,])[1]/length(brownw8upsigdistpeak)
  brown50motifdf[i,"Enrichment.in.W8.DownReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownw8downsigpeak,])[1]/length(brownw8downsigpeak)
  brown50motifdf[i,"Enrichment.in.W8.DownReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownw8downsigprompeak,])[1]/length(brownw8downsigprompeak)
  brown50motifdf[i,"Enrichment.in.W8.DownReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownw8downsigintpeak,])[1]/length(brownw8downsigintpeak)
  brown50motifdf[i,"Enrichment.in.W8.DownReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% brownw8downsigdistpeak,])[1]/length(brownw8downsigdistpeak)
  
}


white50motifdf <- data.frame(row.names = unique(white50peakmotifs$Motif.Name),
                             "Enrichment in Active Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in Sig Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in Sig Prom Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in Sig Int Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in Sig Dist Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in W1 UpReg Sig Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in W1 UpReg Sig Prom Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in W1 UpReg Sig Int Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in W1 UpReg Sig Dist Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in W1 DownReg Sig Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in W1 DownReg Sig Prom Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in W1 DownReg Sig Int Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in W1 DownReg Sig Dist Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in W2 UpReg Sig Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in W2 UpReg Sig Prom Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in W2 UpReg Sig Int Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in W2 UpReg Sig Dist Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in W2 DownReg Sig Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in W2 DownReg Sig Prom Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in W2 DownReg Sig Int Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in W2 DownReg Sig Dist Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in W4 UpReg Sig Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in W4 UpReg Sig Prom Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in W4 UpReg Sig Int Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in W4 UpReg Sig Dist Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in W4 DownReg Sig Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in W4 DownReg Sig Prom Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in W4 DownReg Sig Int Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in W4 DownReg Sig Dist Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in W8 UpReg Sig Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in W8 UpReg Sig Prom Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in W8 UpReg Sig Int Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in W8 UpReg Sig Dist Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in W8 DownReg Sig Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in W8 DownReg Sig Prom Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in W8 DownReg Sig Int Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))),
                             "Enrichment in W8 DownReg Sig Dist Peaks" = rep(0,length(unique(white50peakmotifs$Motif.Name))))

for(i in 1:length(rownames(white50motifdf))){
  
  if(i%%20 == 0){
    print(i)
  }
  
  ourmotif <- rownames(white50motifdf)[i]
  ourmotifdf <- white50peakmotifs[white50peakmotifs$Motif.Name %in% ourmotif,]
  
  white50motifdf[i,"Enrichment.in.Active.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whiteactivepeaks,])[1]/length(whiteactivepeaks)
  white50motifdf[i,"Enrichment.in.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitesigpeak,])[1]/length(whitesigpeak)
  white50motifdf[i,"Enrichment.in.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitesigprompeak,])[1]/length(whitesigprompeak)
  white50motifdf[i,"Enrichment.in.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitesigintpeak,])[1]/length(whitesigintpeak)
  white50motifdf[i,"Enrichment.in.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitesigdistpeak,])[1]/length(whitesigdistpeak)
  white50motifdf[i,"Enrichment.in.W1.UpReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitew1upsigpeak,])[1]/length(whitew1upsigpeak)
  white50motifdf[i,"Enrichment.in.W1.UpReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitew1upsigprompeak,])[1]/length(whitew1upsigprompeak)
  white50motifdf[i,"Enrichment.in.W1.UpReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitew1upsigintpeak,])[1]/length(whitew1upsigintpeak)
  white50motifdf[i,"Enrichment.in.W1.UpReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitew1upsigdistpeak,])[1]/length(whitew1upsigdistpeak)
  white50motifdf[i,"Enrichment.in.W1.DownReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitew1downsigpeak,])[1]/length(whitew1downsigpeak)
  white50motifdf[i,"Enrichment.in.W1.DownReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitew1downsigprompeak,])[1]/length(whitew1downsigprompeak)
  white50motifdf[i,"Enrichment.in.W1.DownReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitew1downsigintpeak,])[1]/length(whitew1downsigintpeak)
  white50motifdf[i,"Enrichment.in.W1.DownReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitew1downsigdistpeak,])[1]/length(whitew1downsigdistpeak)
  white50motifdf[i,"Enrichment.in.W2.UpReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitew2upsigpeak,])[1]/length(whitew2upsigpeak)
  white50motifdf[i,"Enrichment.in.W2.UpReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitew2upsigprompeak,])[1]/length(whitew2upsigprompeak)
  white50motifdf[i,"Enrichment.in.W2.UpReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitew2upsigintpeak,])[1]/length(whitew2upsigintpeak)
  white50motifdf[i,"Enrichment.in.W2.UpReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitew2upsigdistpeak,])[1]/length(whitew2upsigdistpeak)
  white50motifdf[i,"Enrichment.in.W2.DownReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitew2downsigpeak,])[1]/length(whitew2downsigpeak)
  white50motifdf[i,"Enrichment.in.W2.DownReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitew2downsigprompeak,])[1]/length(whitew2downsigprompeak)
  white50motifdf[i,"Enrichment.in.W2.DownReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitew2downsigintpeak,])[1]/length(whitew2downsigintpeak)
  white50motifdf[i,"Enrichment.in.W2.DownReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitew2downsigdistpeak,])[1]/length(whitew2downsigdistpeak)
  white50motifdf[i,"Enrichment.in.W4.UpReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitew4upsigpeak,])[1]/length(whitew4upsigpeak)
  white50motifdf[i,"Enrichment.in.W4.UpReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitew4upsigprompeak,])[1]/length(whitew4upsigprompeak)
  white50motifdf[i,"Enrichment.in.W4.UpReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitew4upsigintpeak,])[1]/length(whitew4upsigintpeak)
  white50motifdf[i,"Enrichment.in.W4.UpReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitew4upsigdistpeak,])[1]/length(whitew4upsigdistpeak)
  white50motifdf[i,"Enrichment.in.W4.DownReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitew4downsigpeak,])[1]/length(whitew4downsigpeak)
  white50motifdf[i,"Enrichment.in.W4.DownReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitew4downsigprompeak,])[1]/length(whitew4downsigprompeak)
  white50motifdf[i,"Enrichment.in.W4.DownReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitew4downsigintpeak,])[1]/length(whitew4downsigintpeak)
  white50motifdf[i,"Enrichment.in.W4.DownReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitew4downsigdistpeak,])[1]/length(whitew4downsigdistpeak)
  white50motifdf[i,"Enrichment.in.W8.UpReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitew8upsigpeak,])[1]/length(whitew8upsigpeak)
  white50motifdf[i,"Enrichment.in.W8.UpReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitew8upsigprompeak,])[1]/length(whitew8upsigprompeak)
  white50motifdf[i,"Enrichment.in.W8.UpReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitew8upsigintpeak,])[1]/length(whitew8upsigintpeak)
  white50motifdf[i,"Enrichment.in.W8.UpReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitew8upsigdistpeak,])[1]/length(whitew8upsigdistpeak)
  white50motifdf[i,"Enrichment.in.W8.DownReg.Sig.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitew8downsigpeak,])[1]/length(whitew8downsigpeak)
  white50motifdf[i,"Enrichment.in.W8.DownReg.Sig.Prom.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitew8downsigprompeak,])[1]/length(whitew8downsigprompeak)
  white50motifdf[i,"Enrichment.in.W8.DownReg.Sig.Int.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitew8downsigintpeak,])[1]/length(whitew8downsigintpeak)
  white50motifdf[i,"Enrichment.in.W8.DownReg.Sig.Dist.Peaks"] <- dim(ourmotifdf[ourmotifdf$PositionID %in% whitew8downsigdistpeak,])[1]/length(whitew8downsigdistpeak)
  
}

motifdfmetadata <- data.frame(row.names = colnames(gastro50motifdf),
                              "Region" = c("All","All","Promoter","Intron","Distal Intergenic",
                                           "All","Promoter","Intron","Distal Intergenic",
                                           "All","Promoter","Intron","Distal Intergenic",
                                           "All","Promoter","Intron","Distal Intergenic",
                                           "All","Promoter","Intron","Distal Intergenic",
                                           "All","Promoter","Intron","Distal Intergenic",
                                           "All","Promoter","Intron","Distal Intergenic",
                                           "All","Promoter","Intron","Distal Intergenic",
                                           "All","Promoter","Intron","Distal Intergenic"),
                              "Training.Response" = c("Background","All.Sig","All.Sig","All.Sig","All.Sig",
                                                      "Up.Reg","Up.Reg","Up.Reg","Up.Reg",
                                                      "Down.Reg","Down.Reg","Down.Reg","Down.Reg",
                                                      "Up.Reg","Up.Reg","Up.Reg","Up.Reg",
                                                      "Down.Reg","Down.Reg","Down.Reg","Down.Reg",
                                                      "Up.Reg","Up.Reg","Up.Reg","Up.Reg",
                                                      "Down.Reg","Down.Reg","Down.Reg","Down.Reg",
                                                      "Up.Reg","Up.Reg","Up.Reg","Up.Reg",
                                                      "Down.Reg","Down.Reg","Down.Reg","Down.Reg"),
                              "Week" = c("Background","All","All","All","All",
                                         "1w","1w","1w","1w","1w","1w","1w","1w",
                                         "2w","2w","2w","2w","2w","2w","2w","2w",
                                         "4w","4w","4w","4w","4w","4w","4w","4w",
                                         "8w","8w","8w","8w","8w","8w","8w","8w"))

altflist <- Reduce(intersect,list(rownames(gastro50motifdf),rownames(heart50motifdf),rownames(hippo50motifdf)))
fintflist <- intersect(tflist,altflist)
fintflabel <- gsub("\\(.*","",fintflist)

#activevsdegcorrmat <- matrix(0L,nrow = 3,ncol = 8)
#rownames(activevsdegcorrmat) <- c("Active vs DEG","Active Prom vs DEG Prom","Up-Reg Prom vs Down-Reg Prom")
#colnames(activevsdegcorrmat) <- c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC")
#activevsdegcorrmat["Active vs DEG","SKM-GN"] <- cor(as.numeric(gsub("%","",gastro_activetf[fintflist,"X..of.Target.Sequences.with.Motif"])),as.numeric(gsub("%","",gastrotf[fintflist,"X..of.Target.Sequences.with.Motif"])))
#activevsdegcorrmat["Active vs DEG","HEART"] <- cor(as.numeric(gsub("%","",heart_activetf[fintflist,"X..of.Target.Sequences.with.Motif"])),as.numeric(gsub("%","",hearttf[fintflist,"X..of.Target.Sequences.with.Motif"])))
#activevsdegcorrmat["Active vs DEG","HIPPOC"] <- cor(as.numeric(gsub("%","",hippo_activetf[fintflist,"X..of.Target.Sequences.with.Motif"])),as.numeric(gsub("%","",hippotf[fintflist,"X..of.Target.Sequences.with.Motif"])))
#activevsdegcorrmat["Active vs DEG","KIDNEY"] <- cor(as.numeric(gsub("%","",kidney_activetf[fintflist,"X..of.Target.Sequences.with.Motif"])),as.numeric(gsub("%","",kidneytf[fintflist,"X..of.Target.Sequences.with.Motif"])))
#activevsdegcorrmat["Active vs DEG","LIVER"] <- cor(as.numeric(gsub("%","",liver_activetf[fintflist,"X..of.Target.Sequences.with.Motif"])),as.numeric(gsub("%","",livertf[fintflist,"X..of.Target.Sequences.with.Motif"])))
#activevsdegcorrmat["Active vs DEG","LUNG"] <- cor(as.numeric(gsub("%","",lung_activetf[fintflist,"X..of.Target.Sequences.with.Motif"])),as.numeric(gsub("%","",lungtf[fintflist,"X..of.Target.Sequences.with.Motif"])))
#activevsdegcorrmat["Active vs DEG","BAT"] <- cor(as.numeric(gsub("%","",brown_activetf[fintflist,"X..of.Target.Sequences.with.Motif"])),as.numeric(gsub("%","",browntf[fintflist,"X..of.Target.Sequences.with.Motif"])))
#activevsdegcorrmat["Active vs DEG","WAT-SC"] <- cor(as.numeric(gsub("%","",white_activetf[fintflist,"X..of.Target.Sequences.with.Motif"])),as.numeric(gsub("%","",whitetf[fintflist,"X..of.Target.Sequences.with.Motif"])))

#activevsdegcorrmat["Active Prom vs DEG Prom","SKM-GN"] <- cor(tf50activeprompeaks[fintflist,c("SKM.GN.Enrichment.in.Active.Prom.Peaks")],100*gastro50motifdf[fintflist,"Enrichment.in.Sig.Prom.Peaks"])
#activevsdegcorrmat["Active Prom vs DEG Prom","HEART"] <- cor(tf50activeprompeaks[fintflist,c("HEART.Enrichment.in.Active.Prom.Peaks")],100*as.numeric(gsub("%","",heart50motifdf[fintflist,"Enrichment.in.Sig.Prom.Peaks"])))
#activevsdegcorrmat["Active Prom vs DEG Prom","HIPPOC"] <- cor(tf50activeprompeaks[fintflist,c("HIPPOC.Enrichment.in.Active.Prom.Peaks")],100*as.numeric(gsub("%","",hippo50motifdf[fintflist,"Enrichment.in.Sig.Prom.Peaks"])))
#activevsdegcorrmat["Active Prom vs DEG Prom","KIDNEY"] <- cor(tf50activeprompeaks[fintflist,c("KIDNEY.Enrichment.in.Active.Prom.Peaks")],100*as.numeric(gsub("%","",kidney50motifdf[fintflist,"Enrichment.in.Sig.Prom.Peaks"])))
#activevsdegcorrmat["Active Prom vs DEG Prom","LIVER"] <- cor(tf50activeprompeaks[fintflist,c("LIVER.Enrichment.in.Active.Prom.Peaks")],100*as.numeric(gsub("%","",liver50motifdf[fintflist,"Enrichment.in.Sig.Prom.Peaks"])))
#activevsdegcorrmat["Active Prom vs DEG Prom","LUNG"] <- cor(tf50activeprompeaks[fintflist,c("LUNG.Enrichment.in.Active.Prom.Peaks")],100*as.numeric(gsub("%","",lung50motifdf[fintflist,"Enrichment.in.Sig.Prom.Peaks"])))
#activevsdegcorrmat["Active Prom vs DEG Prom","BAT"] <- cor(tf50activeprompeaks[fintflist,c("BAT.Enrichment.in.Active.Prom.Peaks")],100*as.numeric(gsub("%","",brown50motifdf[fintflist,"Enrichment.in.Sig.Prom.Peaks"])))
#activevsdegcorrmat["Active Prom vs DEG Prom","WAT-SC"] <- cor(tf50activeprompeaks[fintflist,c("WAT.SC.Enrichment.in.Active.Prom.Peaks")],100*as.numeric(gsub("%","",white50motifdf[fintflist,"Enrichment.in.Sig.Prom.Peaks"])))

#activevsdegcorrmat["Up-Reg Prom vs Down-Reg Prom","SKM-GN"] <- cor(gastro50motifdf[,c("Enrichment.in.W8.UpReg.Sig.Prom.Peaks")],gastro50motifdf[,c("Enrichment.in.W8.DownReg.Sig.Prom.Peaks")])
#activevsdegcorrmat["Up-Reg Prom vs Down-Reg Prom","HEART"] <- cor(heart50motifdf[,c("Enrichment.in.W8.UpReg.Sig.Prom.Peaks")],heart50motifdf[,c("Enrichment.in.W8.DownReg.Sig.Prom.Peaks")])
#activevsdegcorrmat["Up-Reg Prom vs Down-Reg Prom","HIPPOC"] <- cor(hippo50motifdf[,c("Enrichment.in.W8.UpReg.Sig.Prom.Peaks")],hippo50motifdf[,c("Enrichment.in.W8.DownReg.Sig.Prom.Peaks")])
#activevsdegcorrmat["Up-Reg Prom vs Down-Reg Prom","KIDNEY"] <- cor(kidney50motifdf[,c("Enrichment.in.W8.UpReg.Sig.Prom.Peaks")],kidney50motifdf[,c("Enrichment.in.W8.DownReg.Sig.Prom.Peaks")])
#activevsdegcorrmat["Up-Reg Prom vs Down-Reg Prom","LIVER"] <- cor(liver50motifdf[,c("Enrichment.in.W8.UpReg.Sig.Prom.Peaks")],liver50motifdf[,c("Enrichment.in.W8.DownReg.Sig.Prom.Peaks")])
#activevsdegcorrmat["Up-Reg Prom vs Down-Reg Prom","LUNG"] <- cor(lung50motifdf[,c("Enrichment.in.W8.UpReg.Sig.Prom.Peaks")],lung50motifdf[,c("Enrichment.in.W8.DownReg.Sig.Prom.Peaks")])
#activevsdegcorrmat["Up-Reg Prom vs Down-Reg Prom","BAT"] <- cor(brown50motifdf[,c("Enrichment.in.W8.UpReg.Sig.Prom.Peaks")],brown50motifdf[,c("Enrichment.in.W8.DownReg.Sig.Prom.Peaks")])
#activevsdegcorrmat["Up-Reg Prom vs Down-Reg Prom","WAT-SC"] <- cor(white50motifdf[,c("Enrichment.in.W8.UpReg.Sig.Prom.Peaks")],white50motifdf[,c("Enrichment.in.W8.DownReg.Sig.Prom.Peaks")])

#activevsdegcorrdf <- data.frame("Correlation" = as.vector(t(activevsdegcorrmat)),
#                                "Comparison" = c(rep("Active vs DEG",8),
#                                                 rep("Active Prom vs DEG Prom",8),
#                                                 rep("Up-Reg Prom vs Down-Reg Prom",8)),
#                                "Tissue" = rep(c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC"),3))
#activevsdegcorrdf$Comparison <- factor(activevsdegcorrdf$Comparison,levels = c("Active vs DEG",
#                                                                               "Active Prom vs DEG Prom",
#                                                                               "Up-Reg Prom vs Down-Reg Prom"))


# Supplemental Figure S34A
pdf(file = "Supplemental Figure S34A_021325.pdf",width = 7,height = 5.5)
pheatmap(t(scale(t(gastro50motifdf[apply(gastro50motifdf,1,max) > 0.01,order(motifdfmetadata$Region,motifdfmetadata$Training.Response)]))),labels_row = gsub("\\(.*","",rownames(gastro50motifdf))[apply(gastro50motifdf,1,max) > 0.01],annotation_col = motifdfmetadata[,c("Week","Training.Response","Region")],annotation_colors = ann_colsheatmaptrim,show_colnames = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,cluster_cols = F,border_color = NA)
dev.off()
# Supplemental Figure S34B
pdf(file = "Supplemental Figure S34B_021325.pdf",width = 7,height = 5.5)
pheatmap(t(scale(t(heart50motifdf[apply(heart50motifdf,1,max) > 0.01,order(motifdfmetadata$Region,motifdfmetadata$Training.Response)]))),labels_row = gsub("\\(.*","",rownames(heart50motifdf))[apply(heart50motifdf,1,max) > 0.01],annotation_col = motifdfmetadata[,c("Week","Training.Response","Region")],annotation_colors = ann_colsheatmaptrim,show_colnames = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,cluster_cols = F,border_color = NA)
dev.off()
# Supplemental Figure S34C
pdf(file = "Supplemental Figure S34C_021325.pdf",width = 7,height = 5.5)
pheatmap(t(scale(t(hippo50motifdf[apply(hippo50motifdf,1,max) > 0.01,order(motifdfmetadata$Region,motifdfmetadata$Training.Response)]))),labels_row = gsub("\\(.*","",rownames(hippo50motifdf))[apply(hippo50motifdf,1,max) > 0.01],annotation_col = motifdfmetadata[,c("Week","Training.Response","Region")],annotation_colors = ann_colsheatmaptrim,show_colnames = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,cluster_cols = F,border_color = NA)
dev.off()
# Supplemental Figure S34D
pdf(file = "Supplemental Figure S34D_021325.pdf",width = 7,height = 5.5)
pheatmap(t(scale(t(kidney50motifdf[apply(kidney50motifdf,1,max) > 0.01,order(motifdfmetadata$Region,motifdfmetadata$Training.Response)]))),labels_row = gsub("\\(.*","",rownames(kidney50motifdf))[apply(kidney50motifdf,1,max) > 0.01],annotation_col = motifdfmetadata[,c("Week","Training.Response","Region")],annotation_colors = ann_colsheatmaptrim,show_colnames = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,cluster_cols = F,border_color = NA)
dev.off()
# Supplemental Figure S34E
pdf(file = "Supplemental Figure S34E_021325.pdf",width = 7,height = 5.5)
pheatmap(t(scale(t(liver50motifdf[apply(liver50motifdf,1,max) > 0.01,order(motifdfmetadata$Region,motifdfmetadata$Training.Response)]))),labels_row = gsub("\\(.*","",rownames(liver50motifdf))[apply(liver50motifdf,1,max) > 0.01],annotation_col = motifdfmetadata[,c("Week","Training.Response","Region")],annotation_colors = ann_colsheatmaptrim,show_colnames = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,cluster_cols = F,border_color = NA)
dev.off()
# Supplemental Figure S34F
pdf(file = "Supplemental Figure S34F_021325.pdf",width = 7,height = 5.5)
pheatmap(t(scale(t(lung50motifdf[apply(lung50motifdf,1,max) > 0.01,order(motifdfmetadata$Region,motifdfmetadata$Training.Response)]))),labels_row = gsub("\\(.*","",rownames(lung50motifdf))[apply(lung50motifdf,1,max) > 0.01],annotation_col = motifdfmetadata[,c("Week","Training.Response","Region")],annotation_colors = ann_colsheatmaptrim,show_colnames = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,cluster_cols = F,border_color = NA)
dev.off()
# Supplemental Figure S34G
pdf(file = "Supplemental Figure S34G_021325.pdf",width = 7,height = 5.5)
pheatmap(t(scale(t(brown50motifdf[apply(brown50motifdf,1,max) > 0.01,order(motifdfmetadata$Region,motifdfmetadata$Training.Response)]))),labels_row = gsub("\\(.*","",rownames(brown50motifdf))[apply(brown50motifdf,1,max) > 0.01],annotation_col = motifdfmetadata[,c("Week","Training.Response","Region")],annotation_colors = ann_colsheatmaptrim,show_colnames = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,cluster_cols = F,border_color = NA)
dev.off()
# Supplemental Figure S34H
pdf(file = "Supplemental Figure S34H_021325.pdf",width = 7,height = 5.5)
pheatmap(t(scale(t(white50motifdf[apply(white50motifdf,1,max) > 0.01,order(motifdfmetadata$Region,motifdfmetadata$Training.Response)]))),labels_row = gsub("\\(.*","",rownames(white50motifdf))[apply(white50motifdf,1,max) > 0.01],annotation_col = motifdfmetadata[,c("Week","Training.Response","Region")],annotation_colors = ann_colsheatmaptrim,show_colnames = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,cluster_cols = F,border_color = NA)
dev.off()

png(file = "Supplemental Figure S34A_021325.png",width = 7,height = 5.5,units = "in",res = 600)
pheatmap(t(scale(t(gastro50motifdf[apply(gastro50motifdf,1,max) > 0.01,order(motifdfmetadata$Region,motifdfmetadata$Training.Response)]))),labels_row = gsub("\\(.*","",rownames(gastro50motifdf))[apply(gastro50motifdf,1,max) > 0.01],annotation_col = motifdfmetadata[,c("Week","Training.Response","Region")],annotation_colors = ann_colsheatmaptrim,show_colnames = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,cluster_cols = F,border_color = NA)
dev.off()
# Supplemental Figure S34B
png(file = "Supplemental Figure S34B_021325.png",width = 7,height = 5.5,units = "in",res = 600)
pheatmap(t(scale(t(heart50motifdf[apply(heart50motifdf,1,max) > 0.01,order(motifdfmetadata$Region,motifdfmetadata$Training.Response)]))),labels_row = gsub("\\(.*","",rownames(heart50motifdf))[apply(heart50motifdf,1,max) > 0.01],annotation_col = motifdfmetadata[,c("Week","Training.Response","Region")],annotation_colors = ann_colsheatmaptrim,show_colnames = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,cluster_cols = F,border_color = NA)
dev.off()
# Supplemental Figure S34C
png(file = "Supplemental Figure S34C_021325.png",width = 7,height = 5.5,units = "in",res = 600)
pheatmap(t(scale(t(hippo50motifdf[apply(hippo50motifdf,1,max) > 0.01,order(motifdfmetadata$Region,motifdfmetadata$Training.Response)]))),labels_row = gsub("\\(.*","",rownames(hippo50motifdf))[apply(hippo50motifdf,1,max) > 0.01],annotation_col = motifdfmetadata[,c("Week","Training.Response","Region")],annotation_colors = ann_colsheatmaptrim,show_colnames = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,cluster_cols = F,border_color = NA)
dev.off()
# Supplemental Figure S34D
png(file = "Supplemental Figure S34D_021325.png",width = 7,height = 5.5,units = "in",res = 600)
pheatmap(t(scale(t(kidney50motifdf[apply(kidney50motifdf,1,max) > 0.01,order(motifdfmetadata$Region,motifdfmetadata$Training.Response)]))),labels_row = gsub("\\(.*","",rownames(kidney50motifdf))[apply(kidney50motifdf,1,max) > 0.01],annotation_col = motifdfmetadata[,c("Week","Training.Response","Region")],annotation_colors = ann_colsheatmaptrim,show_colnames = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,cluster_cols = F,border_color = NA)
dev.off()
# Supplemental Figure S34E
png(file = "Supplemental Figure S34E_021325.png",width = 7,height = 5.5,units = "in",res = 600)
pheatmap(t(scale(t(liver50motifdf[apply(liver50motifdf,1,max) > 0.01,order(motifdfmetadata$Region,motifdfmetadata$Training.Response)]))),labels_row = gsub("\\(.*","",rownames(liver50motifdf))[apply(liver50motifdf,1,max) > 0.01],annotation_col = motifdfmetadata[,c("Week","Training.Response","Region")],annotation_colors = ann_colsheatmaptrim,show_colnames = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,cluster_cols = F,border_color = NA)
dev.off()
# Supplemental Figure S34F
png(file = "Supplemental Figure S34F_021325.png",width = 7,height = 5.5,units = "in",res = 600)
pheatmap(t(scale(t(lung50motifdf[apply(lung50motifdf,1,max) > 0.01,order(motifdfmetadata$Region,motifdfmetadata$Training.Response)]))),labels_row = gsub("\\(.*","",rownames(lung50motifdf))[apply(lung50motifdf,1,max) > 0.01],annotation_col = motifdfmetadata[,c("Week","Training.Response","Region")],annotation_colors = ann_colsheatmaptrim,show_colnames = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,cluster_cols = F,border_color = NA)
dev.off()
# Supplemental Figure S34G
png(file = "Supplemental Figure S34G_021325.png",width = 7,height = 5.5,units = "in",res = 600)
pheatmap(t(scale(t(brown50motifdf[apply(brown50motifdf,1,max) > 0.01,order(motifdfmetadata$Region,motifdfmetadata$Training.Response)]))),labels_row = gsub("\\(.*","",rownames(brown50motifdf))[apply(brown50motifdf,1,max) > 0.01],annotation_col = motifdfmetadata[,c("Week","Training.Response","Region")],annotation_colors = ann_colsheatmaptrim,show_colnames = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,cluster_cols = F,border_color = NA)
dev.off()
# Supplemental Figure S34H
png(file = "Supplemental Figure S34H_021325.png",width = 7,height = 5.5,units = "in",res = 600)
pheatmap(t(scale(t(white50motifdf[apply(white50motifdf,1,max) > 0.01,order(motifdfmetadata$Region,motifdfmetadata$Training.Response)]))),labels_row = gsub("\\(.*","",rownames(white50motifdf))[apply(white50motifdf,1,max) > 0.01],annotation_col = motifdfmetadata[,c("Week","Training.Response","Region")],annotation_colors = ann_colsheatmaptrim,show_colnames = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,cluster_cols = F,border_color = NA)
dev.off()

#####
# Supplemental Figure S35
####

scalegastro50motifdf <- t(scale(t(gastro50motifdf[apply(gastro50motifdf,1,max) > 0.02,])))
gastromotif50promupreg <- Reduce(union,list(rownames(scalegastro50motifdf[scalegastro50motifdf[,"Enrichment.in.W1.UpReg.Sig.Prom.Peaks"] > 2,]),
                                            rownames(scalegastro50motifdf[scalegastro50motifdf[,"Enrichment.in.W2.UpReg.Sig.Prom.Peaks"] > 2,]),
                                            rownames(scalegastro50motifdf[scalegastro50motifdf[,"Enrichment.in.W4.UpReg.Sig.Prom.Peaks"] > 2,]),
                                            rownames(scalegastro50motifdf[scalegastro50motifdf[,"Enrichment.in.W8.UpReg.Sig.Prom.Peaks"] > 2,])))
gastromotif50promdownreg <- Reduce(union,list(rownames(scalegastro50motifdf[scalegastro50motifdf[,"Enrichment.in.W1.DownReg.Sig.Prom.Peaks"] > 2,]),
                                              rownames(scalegastro50motifdf[scalegastro50motifdf[,"Enrichment.in.W2.DownReg.Sig.Prom.Peaks"] > 2,]),
                                              rownames(scalegastro50motifdf[scalegastro50motifdf[,"Enrichment.in.W4.DownReg.Sig.Prom.Peaks"] > 2,]),
                                              rownames(scalegastro50motifdf[scalegastro50motifdf[,"Enrichment.in.W8.DownReg.Sig.Prom.Peaks"] > 2,])))
gastromotif50promupregspec <- gastromotif50promupreg[!gastromotif50promupreg %in% gastromotif50promdownreg]
gastromotif50promdownregspec <- gastromotif50promdownreg[!gastromotif50promdownreg %in% gastromotif50promupreg]


scaleheart50motifdf <- t(scale(t(heart50motifdf[apply(heart50motifdf,1,max) > 0.02,])))
heartmotif50promupreg <- Reduce(union,list(rownames(scaleheart50motifdf[scaleheart50motifdf[,"Enrichment.in.W1.UpReg.Sig.Prom.Peaks"] > 2,]),
                                           rownames(scaleheart50motifdf[scaleheart50motifdf[,"Enrichment.in.W2.UpReg.Sig.Prom.Peaks"] > 2,]),
                                           rownames(scaleheart50motifdf[scaleheart50motifdf[,"Enrichment.in.W4.UpReg.Sig.Prom.Peaks"] > 2,]),
                                           rownames(scaleheart50motifdf[scaleheart50motifdf[,"Enrichment.in.W8.UpReg.Sig.Prom.Peaks"] > 2,])))
heartmotif50promdownreg <- Reduce(union,list(rownames(scaleheart50motifdf[scaleheart50motifdf[,"Enrichment.in.W1.DownReg.Sig.Prom.Peaks"] > 2,]),
                                             rownames(scaleheart50motifdf[scaleheart50motifdf[,"Enrichment.in.W2.DownReg.Sig.Prom.Peaks"] > 2,]),
                                             rownames(scaleheart50motifdf[scaleheart50motifdf[,"Enrichment.in.W4.DownReg.Sig.Prom.Peaks"] > 2,]),
                                             rownames(scaleheart50motifdf[scaleheart50motifdf[,"Enrichment.in.W8.DownReg.Sig.Prom.Peaks"] > 2,])))
heartmotif50promupregspec <- heartmotif50promupreg[!heartmotif50promupreg %in% heartmotif50promdownreg]
heartmotif50promdownregspec <- heartmotif50promdownreg[!heartmotif50promdownreg %in% heartmotif50promupreg]


scalehippo50motifdf <- t(scale(t(hippo50motifdf[apply(hippo50motifdf,1,max) > 0.02,])))
hippomotif50promupreg <- Reduce(union,list(rownames(scalehippo50motifdf[scalehippo50motifdf[,"Enrichment.in.W1.UpReg.Sig.Prom.Peaks"] > 2,]),
                                           rownames(scalehippo50motifdf[scalehippo50motifdf[,"Enrichment.in.W2.UpReg.Sig.Prom.Peaks"] > 2,]),
                                           rownames(scalehippo50motifdf[scalehippo50motifdf[,"Enrichment.in.W4.UpReg.Sig.Prom.Peaks"] > 2,]),
                                           rownames(scalehippo50motifdf[scalehippo50motifdf[,"Enrichment.in.W8.UpReg.Sig.Prom.Peaks"] > 2,])))
hippomotif50promdownreg <- Reduce(union,list(rownames(scalehippo50motifdf[scalehippo50motifdf[,"Enrichment.in.W1.DownReg.Sig.Prom.Peaks"] > 2,]),
                                             rownames(scalehippo50motifdf[scalehippo50motifdf[,"Enrichment.in.W2.DownReg.Sig.Prom.Peaks"] > 2,]),
                                             rownames(scalehippo50motifdf[scalehippo50motifdf[,"Enrichment.in.W4.DownReg.Sig.Prom.Peaks"] > 2,]),
                                             rownames(scalehippo50motifdf[scalehippo50motifdf[,"Enrichment.in.W8.DownReg.Sig.Prom.Peaks"] > 2,])))
hippomotif50promupregspec <- hippomotif50promupreg[!hippomotif50promupreg %in% hippomotif50promdownreg]
hippomotif50promdownregspec <- hippomotif50promdownreg[!hippomotif50promdownreg %in% hippomotif50promupreg]


scalekidney50motifdf <- t(scale(t(kidney50motifdf[apply(kidney50motifdf,1,max) > 0.02,])))
kidneymotif50promupreg <- Reduce(union,list(rownames(scalekidney50motifdf[scalekidney50motifdf[,"Enrichment.in.W1.UpReg.Sig.Prom.Peaks"] > 2,]),
                                            rownames(scalekidney50motifdf[scalekidney50motifdf[,"Enrichment.in.W2.UpReg.Sig.Prom.Peaks"] > 2,]),
                                            rownames(scalekidney50motifdf[scalekidney50motifdf[,"Enrichment.in.W4.UpReg.Sig.Prom.Peaks"] > 2,]),
                                            rownames(scalekidney50motifdf[scalekidney50motifdf[,"Enrichment.in.W8.UpReg.Sig.Prom.Peaks"] > 2,])))
kidneymotif50promdownreg <- Reduce(union,list(rownames(scalekidney50motifdf[scalekidney50motifdf[,"Enrichment.in.W1.DownReg.Sig.Prom.Peaks"] > 2,]),
                                              rownames(scalekidney50motifdf[scalekidney50motifdf[,"Enrichment.in.W2.DownReg.Sig.Prom.Peaks"] > 2,]),
                                              rownames(scalekidney50motifdf[scalekidney50motifdf[,"Enrichment.in.W4.DownReg.Sig.Prom.Peaks"] > 2,]),
                                              rownames(scalekidney50motifdf[scalekidney50motifdf[,"Enrichment.in.W8.DownReg.Sig.Prom.Peaks"] > 2,])))
kidneymotif50promupregspec <- kidneymotif50promupreg[!kidneymotif50promupreg %in% kidneymotif50promdownreg]
kidneymotif50promdownregspec <- kidneymotif50promdownreg[!kidneymotif50promdownreg %in% kidneymotif50promupreg]


scaleliver50motifdf <- t(scale(t(liver50motifdf[apply(liver50motifdf,1,max) > 0.02,])))
livermotif50promupreg <- Reduce(union,list(rownames(scaleliver50motifdf[scaleliver50motifdf[,"Enrichment.in.W1.UpReg.Sig.Prom.Peaks"] > 2,]),
                                           rownames(scaleliver50motifdf[scaleliver50motifdf[,"Enrichment.in.W2.UpReg.Sig.Prom.Peaks"] > 2,]),
                                           rownames(scaleliver50motifdf[scaleliver50motifdf[,"Enrichment.in.W4.UpReg.Sig.Prom.Peaks"] > 2,]),
                                           rownames(scaleliver50motifdf[scaleliver50motifdf[,"Enrichment.in.W8.UpReg.Sig.Prom.Peaks"] > 2,])))
livermotif50promdownreg <- Reduce(union,list(rownames(scaleliver50motifdf[scaleliver50motifdf[,"Enrichment.in.W1.DownReg.Sig.Prom.Peaks"] > 2,]),
                                             rownames(scaleliver50motifdf[scaleliver50motifdf[,"Enrichment.in.W2.DownReg.Sig.Prom.Peaks"] > 2,]),
                                             rownames(scaleliver50motifdf[scaleliver50motifdf[,"Enrichment.in.W4.DownReg.Sig.Prom.Peaks"] > 2,]),
                                             rownames(scaleliver50motifdf[scaleliver50motifdf[,"Enrichment.in.W8.DownReg.Sig.Prom.Peaks"] > 2,])))
livermotif50promupregspec <- livermotif50promupreg[!livermotif50promupreg %in% livermotif50promdownreg]
livermotif50promdownregspec <- livermotif50promdownreg[!livermotif50promdownreg %in% livermotif50promupreg]


scalelung50motifdf <- t(scale(t(lung50motifdf[apply(lung50motifdf,1,max) > 0.02,])))
lungmotif50promupreg <- Reduce(union,list(rownames(scalelung50motifdf[scalelung50motifdf[,"Enrichment.in.W1.UpReg.Sig.Prom.Peaks"] > 2,]),
                                          rownames(scalelung50motifdf[scalelung50motifdf[,"Enrichment.in.W2.UpReg.Sig.Prom.Peaks"] > 2,]),
                                          rownames(scalelung50motifdf[scalelung50motifdf[,"Enrichment.in.W4.UpReg.Sig.Prom.Peaks"] > 2,]),
                                          rownames(scalelung50motifdf[scalelung50motifdf[,"Enrichment.in.W8.UpReg.Sig.Prom.Peaks"] > 2,])))
lungmotif50promdownreg <- Reduce(union,list(rownames(scalelung50motifdf[scalelung50motifdf[,"Enrichment.in.W1.DownReg.Sig.Prom.Peaks"] > 2,]),
                                            rownames(scalelung50motifdf[scalelung50motifdf[,"Enrichment.in.W2.DownReg.Sig.Prom.Peaks"] > 2,]),
                                            rownames(scalelung50motifdf[scalelung50motifdf[,"Enrichment.in.W4.DownReg.Sig.Prom.Peaks"] > 2,]),
                                            rownames(scalelung50motifdf[scalelung50motifdf[,"Enrichment.in.W8.DownReg.Sig.Prom.Peaks"] > 2,])))

lungmotif50promupregspec <- lungmotif50promupreg[!lungmotif50promupreg %in% lungmotif50promdownreg]
lungmotif50promdownregspec <- lungmotif50promdownreg[!lungmotif50promdownreg %in% lungmotif50promupreg]


scalebrown50motifdf <- t(scale(t(brown50motifdf[apply(brown50motifdf,1,max) > 0.02,])))
brownmotif50promupreg <- Reduce(union,list(rownames(scalebrown50motifdf[scalebrown50motifdf[,"Enrichment.in.W1.UpReg.Sig.Prom.Peaks"] > 2,]),
                                           rownames(scalebrown50motifdf[scalebrown50motifdf[,"Enrichment.in.W2.UpReg.Sig.Prom.Peaks"] > 2,]),
                                           rownames(scalebrown50motifdf[scalebrown50motifdf[,"Enrichment.in.W4.UpReg.Sig.Prom.Peaks"] > 2,]),
                                           rownames(scalebrown50motifdf[scalebrown50motifdf[,"Enrichment.in.W8.UpReg.Sig.Prom.Peaks"] > 2,])))
brownmotif50promdownreg <- Reduce(union,list(rownames(scalebrown50motifdf[scalebrown50motifdf[,"Enrichment.in.W1.DownReg.Sig.Prom.Peaks"] > 2,]),
                                             rownames(scalebrown50motifdf[scalebrown50motifdf[,"Enrichment.in.W2.DownReg.Sig.Prom.Peaks"] > 2,]),
                                             rownames(scalebrown50motifdf[scalebrown50motifdf[,"Enrichment.in.W4.DownReg.Sig.Prom.Peaks"] > 2,]),
                                             rownames(scalebrown50motifdf[scalebrown50motifdf[,"Enrichment.in.W8.DownReg.Sig.Prom.Peaks"] > 2,])))
brownmotif50promupregspec <- brownmotif50promupreg[!brownmotif50promupreg %in% brownmotif50promdownreg]
brownmotif50promdownregspec <- brownmotif50promdownreg[!brownmotif50promdownreg %in% brownmotif50promupreg]


scalewhite50motifdf <- t(scale(t(white50motifdf[apply(white50motifdf,1,max) > 0.02,])))
whitemotif50promupreg <- Reduce(union,list(rownames(scalewhite50motifdf[scalewhite50motifdf[,"Enrichment.in.W1.UpReg.Sig.Prom.Peaks"] > 2,]),
                                           rownames(scalewhite50motifdf[scalewhite50motifdf[,"Enrichment.in.W2.UpReg.Sig.Prom.Peaks"] > 2,]),
                                           rownames(scalewhite50motifdf[scalewhite50motifdf[,"Enrichment.in.W4.UpReg.Sig.Prom.Peaks"] > 2,]),
                                           rownames(scalewhite50motifdf[scalewhite50motifdf[,"Enrichment.in.W8.UpReg.Sig.Prom.Peaks"] > 2,])))
whitemotif50promdownreg <- Reduce(union,list(rownames(scalewhite50motifdf[scalewhite50motifdf[,"Enrichment.in.W1.DownReg.Sig.Prom.Peaks"] > 2,]),
                                             rownames(scalewhite50motifdf[scalewhite50motifdf[,"Enrichment.in.W2.DownReg.Sig.Prom.Peaks"] > 2,]),
                                             rownames(scalewhite50motifdf[scalewhite50motifdf[,"Enrichment.in.W4.DownReg.Sig.Prom.Peaks"] > 2,]),
                                             rownames(scalewhite50motifdf[scalewhite50motifdf[,"Enrichment.in.W8.DownReg.Sig.Prom.Peaks"] > 2,])))
whitemotif50promupregspec <- whitemotif50promupreg[!whitemotif50promupreg %in% whitemotif50promdownreg]
whitemotif50promdownregspec <- whitemotif50promdownreg[!whitemotif50promdownreg %in% whitemotif50promupreg]

# Fig S35A - 750x1500
pdf(file = "Supplemental Figure S35A_021325.pdf",width = 7.5,height = 15)
pheatmap(t(scale(t(gastro50motifdf[c(gastromotif50promupregspec,gastromotif50promdownregspec),
                                   c("Enrichment.in.W1.UpReg.Sig.Prom.Peaks",
                                     "Enrichment.in.W2.UpReg.Sig.Prom.Peaks",
                                     "Enrichment.in.W4.UpReg.Sig.Prom.Peaks",
                                     "Enrichment.in.W8.UpReg.Sig.Prom.Peaks",
                                     "Enrichment.in.W1.DownReg.Sig.Prom.Peaks",
                                     "Enrichment.in.W2.DownReg.Sig.Prom.Peaks",
                                     "Enrichment.in.W4.DownReg.Sig.Prom.Peaks",
                                     "Enrichment.in.W8.DownReg.Sig.Prom.Peaks")]))),
         labels_row = gsub("\\(.*","",c(gastromotif50promupregspec,gastromotif50promdownregspec)),
         breaks = seq(-2,2,length.out = 101),
         color = colorpanel(101,"blue","white","red"),
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = motifdfmetadata[,c("Week","Training.Response")],
         annotation_colors = ann_colsheatmaptrim,
         show_colnames = F,
         fontsize = 15,
         border_color = NA)
dev.off()
png(file = "Supplemental Figure S35A_021325.png",width = 7.5,height = 15,units = "in",res = 600)
pheatmap(t(scale(t(gastro50motifdf[c(gastromotif50promupregspec,gastromotif50promdownregspec),
                                   c("Enrichment.in.W1.UpReg.Sig.Prom.Peaks",
                                     "Enrichment.in.W2.UpReg.Sig.Prom.Peaks",
                                     "Enrichment.in.W4.UpReg.Sig.Prom.Peaks",
                                     "Enrichment.in.W8.UpReg.Sig.Prom.Peaks",
                                     "Enrichment.in.W1.DownReg.Sig.Prom.Peaks",
                                     "Enrichment.in.W2.DownReg.Sig.Prom.Peaks",
                                     "Enrichment.in.W4.DownReg.Sig.Prom.Peaks",
                                     "Enrichment.in.W8.DownReg.Sig.Prom.Peaks")]))),
         labels_row = gsub("\\(.*","",c(gastromotif50promupregspec,gastromotif50promdownregspec)),
         breaks = seq(-2,2,length.out = 101),
         color = colorpanel(101,"blue","white","red"),
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = motifdfmetadata[,c("Week","Training.Response")],
         annotation_colors = ann_colsheatmaptrim,
         show_colnames = F,
         fontsize = 15,
         border_color = NA)
dev.off()

# Fig S35B - 750x1500
pdf(file = "Supplemental Figure S35B_021325.pdf",width = 7.5,height = 15)
pheatmap(t(scale(t(heart50motifdf[c(heartmotif50promupregspec,heartmotif50promdownregspec),
                                  c("Enrichment.in.W1.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W2.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W4.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W8.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W1.DownReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W2.DownReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W4.DownReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W8.DownReg.Sig.Prom.Peaks")]))),
         labels_row = gsub("\\(.*","",c(heartmotif50promupregspec,heartmotif50promdownregspec)),
         breaks = seq(-2,2,length.out = 101),
         color = colorpanel(101,"blue","white","red"),
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = motifdfmetadata[,c("Week","Training.Response")],
         annotation_colors = ann_colsheatmaptrim,
         show_colnames = F,
         fontsize = 15,
         border_color = NA)
dev.off()
png(file = "Supplemental Figure S35B_021325.png",width = 7.5,height = 15,units = "in",res = 600)
pheatmap(t(scale(t(heart50motifdf[c(heartmotif50promupregspec,heartmotif50promdownregspec),
                                  c("Enrichment.in.W1.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W2.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W4.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W8.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W1.DownReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W2.DownReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W4.DownReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W8.DownReg.Sig.Prom.Peaks")]))),
         labels_row = gsub("\\(.*","",c(heartmotif50promupregspec,heartmotif50promdownregspec)),
         breaks = seq(-2,2,length.out = 101),
         color = colorpanel(101,"blue","white","red"),
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = motifdfmetadata[,c("Week","Training.Response")],
         annotation_colors = ann_colsheatmaptrim,
         show_colnames = F,
         fontsize = 15,
         border_color = NA)
dev.off()

# Fig S35C - 750x1500
pdf(file = "Supplemental Figure S35C_021325.pdf",width = 7.5,height = 15)
pheatmap(t(scale(t(hippo50motifdf[c(hippomotif50promupregspec,hippomotif50promdownregspec),
                                  c("Enrichment.in.W1.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W2.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W4.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W8.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W1.DownReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W2.DownReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W4.DownReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W8.DownReg.Sig.Prom.Peaks")]))),
         labels_row = gsub("\\(.*","",c(hippomotif50promupregspec,hippomotif50promdownregspec)),
         breaks = seq(-2,2,length.out = 101),
         color = colorpanel(101,"blue","white","red"),
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = motifdfmetadata[,c("Week","Training.Response")],
         annotation_colors = ann_colsheatmaptrim,
         show_colnames = F,
         fontsize = 15,
         border_color = NA)
dev.off()
png(file = "Supplemental Figure S35C_021325.png",width = 7.5,height = 15,units = "in",res = 600)
pheatmap(t(scale(t(hippo50motifdf[c(hippomotif50promupregspec,hippomotif50promdownregspec),
                                  c("Enrichment.in.W1.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W2.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W4.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W8.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W1.DownReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W2.DownReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W4.DownReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W8.DownReg.Sig.Prom.Peaks")]))),
         labels_row = gsub("\\(.*","",c(hippomotif50promupregspec,hippomotif50promdownregspec)),
         breaks = seq(-2,2,length.out = 101),
         color = colorpanel(101,"blue","white","red"),
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = motifdfmetadata[,c("Week","Training.Response")],
         annotation_colors = ann_colsheatmaptrim,
         show_colnames = F,
         fontsize = 15,
         border_color = NA)
dev.off()

# Fig S35D - 750x1500
pdf(file = "Supplemental Figure S35D_021325.pdf",width = 7.5,height = 15)
pheatmap(t(scale(t(kidney50motifdf[c(kidneymotif50promupregspec,kidneymotif50promdownregspec),
                                   c("Enrichment.in.W1.UpReg.Sig.Prom.Peaks",
                                     "Enrichment.in.W2.UpReg.Sig.Prom.Peaks",
                                     "Enrichment.in.W4.UpReg.Sig.Prom.Peaks",
                                     "Enrichment.in.W8.UpReg.Sig.Prom.Peaks",
                                     "Enrichment.in.W1.DownReg.Sig.Prom.Peaks",
                                     "Enrichment.in.W2.DownReg.Sig.Prom.Peaks",
                                     "Enrichment.in.W4.DownReg.Sig.Prom.Peaks",
                                     "Enrichment.in.W8.DownReg.Sig.Prom.Peaks")]))),
         labels_row = gsub("\\(.*","",c(kidneymotif50promupregspec,kidneymotif50promdownregspec)),
         breaks = seq(-2,2,length.out = 101),
         color = colorpanel(101,"blue","white","red"),
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = motifdfmetadata[,c("Week","Training.Response")],
         annotation_colors = ann_colsheatmaptrim,
         show_colnames = F,
         fontsize = 15,
         border_color = NA)
dev.off()
png(file = "Supplemental Figure S35D_021325.png",width = 7.5,height = 15,units = "in",res = 600)
pheatmap(t(scale(t(kidney50motifdf[c(kidneymotif50promupregspec,kidneymotif50promdownregspec),
                                   c("Enrichment.in.W1.UpReg.Sig.Prom.Peaks",
                                     "Enrichment.in.W2.UpReg.Sig.Prom.Peaks",
                                     "Enrichment.in.W4.UpReg.Sig.Prom.Peaks",
                                     "Enrichment.in.W8.UpReg.Sig.Prom.Peaks",
                                     "Enrichment.in.W1.DownReg.Sig.Prom.Peaks",
                                     "Enrichment.in.W2.DownReg.Sig.Prom.Peaks",
                                     "Enrichment.in.W4.DownReg.Sig.Prom.Peaks",
                                     "Enrichment.in.W8.DownReg.Sig.Prom.Peaks")]))),
         labels_row = gsub("\\(.*","",c(kidneymotif50promupregspec,kidneymotif50promdownregspec)),
         breaks = seq(-2,2,length.out = 101),
         color = colorpanel(101,"blue","white","red"),
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = motifdfmetadata[,c("Week","Training.Response")],
         annotation_colors = ann_colsheatmaptrim,
         show_colnames = F,
         fontsize = 15,
         border_color = NA)
dev.off()

# Fig S35E - 750x1500
pdf(file = "Supplemental Figure S35E_021325.pdf",width = 7.5,height = 15)
pheatmap(t(scale(t(liver50motifdf[c(livermotif50promupregspec,livermotif50promdownregspec),
                                  c("Enrichment.in.W1.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W2.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W4.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W8.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W1.DownReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W2.DownReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W4.DownReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W8.DownReg.Sig.Prom.Peaks")]))),
         labels_row = gsub("\\(.*","",c(livermotif50promupregspec,livermotif50promdownregspec)),
         breaks = seq(-2,2,length.out = 101),
         color = colorpanel(101,"blue","white","red"),
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = motifdfmetadata[,c("Week","Training.Response")],
         annotation_colors = ann_colsheatmaptrim,
         show_colnames = F,
         fontsize = 15,
         border_color = NA)
dev.off()
png(file = "Supplemental Figure S35E_021325.png",width = 7.5,height = 15,units = "in",res = 600)
pheatmap(t(scale(t(liver50motifdf[c(livermotif50promupregspec,livermotif50promdownregspec),
                                  c("Enrichment.in.W1.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W2.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W4.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W8.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W1.DownReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W2.DownReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W4.DownReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W8.DownReg.Sig.Prom.Peaks")]))),
         labels_row = gsub("\\(.*","",c(livermotif50promupregspec,livermotif50promdownregspec)),
         breaks = seq(-2,2,length.out = 101),
         color = colorpanel(101,"blue","white","red"),
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = motifdfmetadata[,c("Week","Training.Response")],
         annotation_colors = ann_colsheatmaptrim,
         show_colnames = F,
         fontsize = 15,
         border_color = NA)
dev.off()

# Fig S35F - 750x1500
pdf(file = "Supplemental Figure S35F_021325.pdf",width = 7.5,height = 15)
pheatmap(t(scale(t(lung50motifdf[c(lungmotif50promupregspec,lungmotif50promdownregspec),
                                 c("Enrichment.in.W1.UpReg.Sig.Prom.Peaks",
                                   "Enrichment.in.W2.UpReg.Sig.Prom.Peaks",
                                   "Enrichment.in.W4.UpReg.Sig.Prom.Peaks",
                                   "Enrichment.in.W8.UpReg.Sig.Prom.Peaks",
                                   "Enrichment.in.W1.DownReg.Sig.Prom.Peaks",
                                   "Enrichment.in.W2.DownReg.Sig.Prom.Peaks",
                                   "Enrichment.in.W4.DownReg.Sig.Prom.Peaks",
                                   "Enrichment.in.W8.DownReg.Sig.Prom.Peaks")]))),
         labels_row = gsub("\\(.*","",c(lungmotif50promupregspec,lungmotif50promdownregspec)),
         breaks = seq(-2,2,length.out = 101),
         color = colorpanel(101,"blue","white","red"),
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = motifdfmetadata[,c("Week","Training.Response")],
         annotation_colors = ann_colsheatmaptrim,
         show_colnames = F,
         fontsize = 15,
         border_color = NA)
dev.off()
png(file = "Supplemental Figure S35F_021325.png",width = 7.5,height = 15,units = "in",res = 600)
pheatmap(t(scale(t(lung50motifdf[c(lungmotif50promupregspec,lungmotif50promdownregspec),
                                 c("Enrichment.in.W1.UpReg.Sig.Prom.Peaks",
                                   "Enrichment.in.W2.UpReg.Sig.Prom.Peaks",
                                   "Enrichment.in.W4.UpReg.Sig.Prom.Peaks",
                                   "Enrichment.in.W8.UpReg.Sig.Prom.Peaks",
                                   "Enrichment.in.W1.DownReg.Sig.Prom.Peaks",
                                   "Enrichment.in.W2.DownReg.Sig.Prom.Peaks",
                                   "Enrichment.in.W4.DownReg.Sig.Prom.Peaks",
                                   "Enrichment.in.W8.DownReg.Sig.Prom.Peaks")]))),
         labels_row = gsub("\\(.*","",c(lungmotif50promupregspec,lungmotif50promdownregspec)),
         breaks = seq(-2,2,length.out = 101),
         color = colorpanel(101,"blue","white","red"),
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = motifdfmetadata[,c("Week","Training.Response")],
         annotation_colors = ann_colsheatmaptrim,
         show_colnames = F,
         fontsize = 15,
         border_color = NA)
dev.off()

# Fig S35G - 750x1500
pdf(file = "Supplemental Figure S35G_021325.pdf",width = 7.5,height = 15)
pheatmap(t(scale(t(brown50motifdf[c(brownmotif50promupregspec,brownmotif50promdownregspec),
                                  c("Enrichment.in.W1.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W2.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W4.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W8.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W1.DownReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W2.DownReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W4.DownReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W8.DownReg.Sig.Prom.Peaks")]))),
         labels_row = gsub("\\(.*","",c(brownmotif50promupregspec,brownmotif50promdownregspec)),
         breaks = seq(-2,2,length.out = 101),
         color = colorpanel(101,"blue","white","red"),
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = motifdfmetadata[,c("Week","Training.Response")],
         annotation_colors = ann_colsheatmaptrim,
         show_colnames = F,
         fontsize = 15,
         border_color = NA)
dev.off()
png(file = "Supplemental Figure S35G_021325.png",width = 7.5,height = 15,units = "in",res = 600)
pheatmap(t(scale(t(brown50motifdf[c(brownmotif50promupregspec,brownmotif50promdownregspec),
                                  c("Enrichment.in.W1.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W2.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W4.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W8.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W1.DownReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W2.DownReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W4.DownReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W8.DownReg.Sig.Prom.Peaks")]))),
         labels_row = gsub("\\(.*","",c(brownmotif50promupregspec,brownmotif50promdownregspec)),
         breaks = seq(-2,2,length.out = 101),
         color = colorpanel(101,"blue","white","red"),
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = motifdfmetadata[,c("Week","Training.Response")],
         annotation_colors = ann_colsheatmaptrim,
         show_colnames = F,
         fontsize = 15,
         border_color = NA)
dev.off()

# Fig S35H - 750x1500
pdf(file = "Supplemental Figure S35H_021325.pdf",width = 7.5,height = 15)
pheatmap(t(scale(t(white50motifdf[c(whitemotif50promupregspec,whitemotif50promdownregspec),
                                  c("Enrichment.in.W1.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W2.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W4.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W8.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W1.DownReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W2.DownReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W4.DownReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W8.DownReg.Sig.Prom.Peaks")]))),
         labels_row = gsub("\\(.*","",c(whitemotif50promupregspec,whitemotif50promdownregspec)),
         breaks = seq(-2,2,length.out = 101),
         color = colorpanel(101,"blue","white","red"),
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = motifdfmetadata[,c("Week","Training.Response")],
         annotation_colors = ann_colsheatmaptrim,
         show_colnames = F,
         fontsize = 15,
         border_color = NA)
dev.off()
png(file = "Supplemental Figure S35H_021325.png",width = 7.5,height = 15,units = "in",res = 600)
pheatmap(t(scale(t(white50motifdf[c(whitemotif50promupregspec,whitemotif50promdownregspec),
                                  c("Enrichment.in.W1.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W2.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W4.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W8.UpReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W1.DownReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W2.DownReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W4.DownReg.Sig.Prom.Peaks",
                                    "Enrichment.in.W8.DownReg.Sig.Prom.Peaks")]))),
         labels_row = gsub("\\(.*","",c(whitemotif50promupregspec,whitemotif50promdownregspec)),
         breaks = seq(-2,2,length.out = 101),
         color = colorpanel(101,"blue","white","red"),
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = motifdfmetadata[,c("Week","Training.Response")],
         annotation_colors = ann_colsheatmaptrim,
         show_colnames = F,
         fontsize = 15,
         border_color = NA)
dev.off()

#####
# Figure 7
####

upreg50list <- list("SKM-GN" = gastromotif50promupregspec,
                    "HEART" = heartmotif50promupregspec,
                    "HIPPOC" = hippomotif50promupregspec,
                    "KIDNEY" = kidneymotif50promupregspec,
                    "LIVER" = livermotif50promupregspec,
                    "LUNG" = lungmotif50promupregspec,
                    "BAT" = brownmotif50promupregspec,
                    "WAT-SC" = whitemotif50promupregspec)

allupreg50tf <- Reduce(union,upreg50list)
upregupset50mat <- matrix(0L,nrow = length(allupreg50tf),ncol = 8)
rownames(upregupset50mat) <- allupreg50tf
colnames(upregupset50mat) <- names(upreg50list)

for(i in 1:length(allupreg50tf)){
  for(j in 1:length(names(upreg50list))){
    ourtf <- allupreg50tf[i]
    ourtiss <- names(upreg50list)[j]
    if(ourtf %in% upreg50list[[ourtiss]]){
      upregupset50mat[i,j] <- 1
    }
  }
}


downreg50list <- list("SKM-GN" = gastromotif50promdownregspec,
                      "HEART" = heartmotif50promdownregspec,
                      "HIPPOC" = hippomotif50promdownregspec,
                      "KIDNEY" = kidneymotif50promdownregspec,
                      "LIVER" = livermotif50promdownregspec,
                      "LUNG" = lungmotif50promdownregspec,
                      "BAT" = brownmotif50promdownregspec,
                      "WAT-SC" = whitemotif50promdownregspec)

alldownreg50tf <- Reduce(union,downreg50list)
downregdownset50mat <- matrix(0L,nrow = length(alldownreg50tf),ncol = 8)
rownames(downregdownset50mat) <- alldownreg50tf
colnames(downregdownset50mat) <- names(downreg50list)

for(i in 1:length(alldownreg50tf)){
  for(j in 1:length(names(downreg50list))){
    ourtf <- alldownreg50tf[i]
    ourtiss <- names(downreg50list)[j]
    if(ourtf %in% downreg50list[[ourtiss]]){
      downregdownset50mat[i,j] <- 1
    }
  }
}

gastromotif50promdf <- gastro50motifdf[,c("Enrichment.in.W1.UpReg.Sig.Prom.Peaks",
                                          "Enrichment.in.W2.UpReg.Sig.Prom.Peaks",
                                          "Enrichment.in.W4.UpReg.Sig.Prom.Peaks",
                                          "Enrichment.in.W8.UpReg.Sig.Prom.Peaks",
                                          "Enrichment.in.W1.DownReg.Sig.Prom.Peaks",
                                          "Enrichment.in.W2.DownReg.Sig.Prom.Peaks",
                                          "Enrichment.in.W4.DownReg.Sig.Prom.Peaks",
                                          "Enrichment.in.W8.DownReg.Sig.Prom.Peaks")]
gastromotif50promz <- t(scale(t(gastromotif50promdf)))
gastromotif50promz[is.na(gastromotif50promz)] <- 0

heartmotif50promdf <- heart50motifdf[,c("Enrichment.in.W1.UpReg.Sig.Prom.Peaks",
                                        "Enrichment.in.W2.UpReg.Sig.Prom.Peaks",
                                        "Enrichment.in.W4.UpReg.Sig.Prom.Peaks",
                                        "Enrichment.in.W8.UpReg.Sig.Prom.Peaks",
                                        "Enrichment.in.W1.DownReg.Sig.Prom.Peaks",
                                        "Enrichment.in.W2.DownReg.Sig.Prom.Peaks",
                                        "Enrichment.in.W4.DownReg.Sig.Prom.Peaks",
                                        "Enrichment.in.W8.DownReg.Sig.Prom.Peaks")]
heartmotif50promz <- t(scale(t(heartmotif50promdf)))
heartmotif50promz[is.na(heartmotif50promz)] <- 0

hippomotif50promdf <- hippo50motifdf[,c("Enrichment.in.W1.UpReg.Sig.Prom.Peaks",
                                        "Enrichment.in.W2.UpReg.Sig.Prom.Peaks",
                                        "Enrichment.in.W4.UpReg.Sig.Prom.Peaks",
                                        "Enrichment.in.W8.UpReg.Sig.Prom.Peaks",
                                        "Enrichment.in.W1.DownReg.Sig.Prom.Peaks",
                                        "Enrichment.in.W2.DownReg.Sig.Prom.Peaks",
                                        "Enrichment.in.W4.DownReg.Sig.Prom.Peaks",
                                        "Enrichment.in.W8.DownReg.Sig.Prom.Peaks")]
hippomotif50promz <- t(scale(t(hippomotif50promdf)))
hippomotif50promz[is.na(hippomotif50promz)] <- 0

kidneymotif50promdf <- kidney50motifdf[,c("Enrichment.in.W1.UpReg.Sig.Prom.Peaks",
                                          "Enrichment.in.W2.UpReg.Sig.Prom.Peaks",
                                          "Enrichment.in.W4.UpReg.Sig.Prom.Peaks",
                                          "Enrichment.in.W8.UpReg.Sig.Prom.Peaks",
                                          "Enrichment.in.W1.DownReg.Sig.Prom.Peaks",
                                          "Enrichment.in.W2.DownReg.Sig.Prom.Peaks",
                                          "Enrichment.in.W4.DownReg.Sig.Prom.Peaks",
                                          "Enrichment.in.W8.DownReg.Sig.Prom.Peaks")]
kidneymotif50promz <- t(scale(t(kidneymotif50promdf)))
kidneymotif50promz[is.na(kidneymotif50promz)] <- 0

livermotif50promdf <- liver50motifdf[,c("Enrichment.in.W1.UpReg.Sig.Prom.Peaks",
                                        "Enrichment.in.W2.UpReg.Sig.Prom.Peaks",
                                        "Enrichment.in.W4.UpReg.Sig.Prom.Peaks",
                                        "Enrichment.in.W8.UpReg.Sig.Prom.Peaks",
                                        "Enrichment.in.W1.DownReg.Sig.Prom.Peaks",
                                        "Enrichment.in.W2.DownReg.Sig.Prom.Peaks",
                                        "Enrichment.in.W4.DownReg.Sig.Prom.Peaks",
                                        "Enrichment.in.W8.DownReg.Sig.Prom.Peaks")]
livermotif50promz <- t(scale(t(livermotif50promdf)))
livermotif50promz[is.na(livermotif50promz)] <- 0

lungmotif50promdf <- lung50motifdf[,c("Enrichment.in.W1.UpReg.Sig.Prom.Peaks",
                                      "Enrichment.in.W2.UpReg.Sig.Prom.Peaks",
                                      "Enrichment.in.W4.UpReg.Sig.Prom.Peaks",
                                      "Enrichment.in.W8.UpReg.Sig.Prom.Peaks",
                                      "Enrichment.in.W1.DownReg.Sig.Prom.Peaks",
                                      "Enrichment.in.W2.DownReg.Sig.Prom.Peaks",
                                      "Enrichment.in.W4.DownReg.Sig.Prom.Peaks",
                                      "Enrichment.in.W8.DownReg.Sig.Prom.Peaks")]
lungmotif50promz <- t(scale(t(lungmotif50promdf)))
lungmotif50promz[is.na(lungmotif50promz)] <- 0

whitemotif50promdf <- white50motifdf[,c("Enrichment.in.W1.UpReg.Sig.Prom.Peaks",
                                        "Enrichment.in.W2.UpReg.Sig.Prom.Peaks",
                                        "Enrichment.in.W4.UpReg.Sig.Prom.Peaks",
                                        "Enrichment.in.W8.UpReg.Sig.Prom.Peaks",
                                        "Enrichment.in.W1.DownReg.Sig.Prom.Peaks",
                                        "Enrichment.in.W2.DownReg.Sig.Prom.Peaks",
                                        "Enrichment.in.W4.DownReg.Sig.Prom.Peaks",
                                        "Enrichment.in.W8.DownReg.Sig.Prom.Peaks")]
whitemotif50promz <- t(scale(t(whitemotif50promdf)))
whitemotif50promz[is.na(whitemotif50promz)] <- 0

brownmotif50promdf <- brown50motifdf[,c("Enrichment.in.W1.UpReg.Sig.Prom.Peaks",
                                        "Enrichment.in.W2.UpReg.Sig.Prom.Peaks",
                                        "Enrichment.in.W4.UpReg.Sig.Prom.Peaks",
                                        "Enrichment.in.W8.UpReg.Sig.Prom.Peaks",
                                        "Enrichment.in.W1.DownReg.Sig.Prom.Peaks",
                                        "Enrichment.in.W2.DownReg.Sig.Prom.Peaks",
                                        "Enrichment.in.W4.DownReg.Sig.Prom.Peaks",
                                        "Enrichment.in.W8.DownReg.Sig.Prom.Peaks")]
brownmotif50promz <- t(scale(t(brownmotif50promdf)))
brownmotif50promz[is.na(brownmotif50promz)] <- 0

tfmotif50prommeanzmat <- cbind(apply(gastromotif50promz[,c("Enrichment.in.W1.UpReg.Sig.Prom.Peaks",
                                                           "Enrichment.in.W2.UpReg.Sig.Prom.Peaks",
                                                           "Enrichment.in.W4.UpReg.Sig.Prom.Peaks",
                                                           "Enrichment.in.W8.UpReg.Sig.Prom.Peaks")],1,mean),
                               apply(heartmotif50promz[,c("Enrichment.in.W1.UpReg.Sig.Prom.Peaks",
                                                          "Enrichment.in.W2.UpReg.Sig.Prom.Peaks",
                                                          "Enrichment.in.W4.UpReg.Sig.Prom.Peaks",
                                                          "Enrichment.in.W8.UpReg.Sig.Prom.Peaks")],1,mean),
                               apply(hippomotif50promz[,c("Enrichment.in.W1.UpReg.Sig.Prom.Peaks",
                                                          "Enrichment.in.W2.UpReg.Sig.Prom.Peaks",
                                                          "Enrichment.in.W4.UpReg.Sig.Prom.Peaks",
                                                          "Enrichment.in.W8.UpReg.Sig.Prom.Peaks")],1,mean),
                               apply(kidneymotif50promz[,c("Enrichment.in.W1.UpReg.Sig.Prom.Peaks",
                                                           "Enrichment.in.W2.UpReg.Sig.Prom.Peaks",
                                                           "Enrichment.in.W4.UpReg.Sig.Prom.Peaks",
                                                           "Enrichment.in.W8.UpReg.Sig.Prom.Peaks")],1,mean),
                               apply(livermotif50promz[,c("Enrichment.in.W1.UpReg.Sig.Prom.Peaks",
                                                          "Enrichment.in.W2.UpReg.Sig.Prom.Peaks",
                                                          "Enrichment.in.W4.UpReg.Sig.Prom.Peaks",
                                                          "Enrichment.in.W8.UpReg.Sig.Prom.Peaks")],1,mean),
                               apply(lungmotif50promz[,c("Enrichment.in.W1.UpReg.Sig.Prom.Peaks",
                                                         "Enrichment.in.W2.UpReg.Sig.Prom.Peaks",
                                                         "Enrichment.in.W4.UpReg.Sig.Prom.Peaks",
                                                         "Enrichment.in.W8.UpReg.Sig.Prom.Peaks")],1,mean),
                               apply(brownmotif50promz[,c("Enrichment.in.W1.UpReg.Sig.Prom.Peaks",
                                                          "Enrichment.in.W2.UpReg.Sig.Prom.Peaks",
                                                          "Enrichment.in.W4.UpReg.Sig.Prom.Peaks",
                                                          "Enrichment.in.W8.UpReg.Sig.Prom.Peaks")],1,mean),
                               apply(whitemotif50promz[,c("Enrichment.in.W1.UpReg.Sig.Prom.Peaks",
                                                          "Enrichment.in.W2.UpReg.Sig.Prom.Peaks",
                                                          "Enrichment.in.W4.UpReg.Sig.Prom.Peaks",
                                                          "Enrichment.in.W8.UpReg.Sig.Prom.Peaks")],1,mean),
                               apply(gastromotif50promz[,c("Enrichment.in.W1.DownReg.Sig.Prom.Peaks",
                                                           "Enrichment.in.W2.DownReg.Sig.Prom.Peaks",
                                                           "Enrichment.in.W4.DownReg.Sig.Prom.Peaks",
                                                           "Enrichment.in.W8.DownReg.Sig.Prom.Peaks")],1,mean),
                               apply(heartmotif50promz[,c("Enrichment.in.W1.DownReg.Sig.Prom.Peaks",
                                                          "Enrichment.in.W2.DownReg.Sig.Prom.Peaks",
                                                          "Enrichment.in.W4.DownReg.Sig.Prom.Peaks",
                                                          "Enrichment.in.W8.DownReg.Sig.Prom.Peaks")],1,mean),
                               apply(hippomotif50promz[,c("Enrichment.in.W1.DownReg.Sig.Prom.Peaks",
                                                          "Enrichment.in.W2.DownReg.Sig.Prom.Peaks",
                                                          "Enrichment.in.W4.DownReg.Sig.Prom.Peaks",
                                                          "Enrichment.in.W8.DownReg.Sig.Prom.Peaks")],1,mean),
                               apply(kidneymotif50promz[,c("Enrichment.in.W1.DownReg.Sig.Prom.Peaks",
                                                           "Enrichment.in.W2.DownReg.Sig.Prom.Peaks",
                                                           "Enrichment.in.W4.DownReg.Sig.Prom.Peaks",
                                                           "Enrichment.in.W8.DownReg.Sig.Prom.Peaks")],1,mean),
                               apply(livermotif50promz[,c("Enrichment.in.W1.DownReg.Sig.Prom.Peaks",
                                                          "Enrichment.in.W2.DownReg.Sig.Prom.Peaks",
                                                          "Enrichment.in.W4.DownReg.Sig.Prom.Peaks",
                                                          "Enrichment.in.W8.DownReg.Sig.Prom.Peaks")],1,mean),
                               apply(lungmotif50promz[,c("Enrichment.in.W1.DownReg.Sig.Prom.Peaks",
                                                         "Enrichment.in.W2.DownReg.Sig.Prom.Peaks",
                                                         "Enrichment.in.W4.DownReg.Sig.Prom.Peaks",
                                                         "Enrichment.in.W8.DownReg.Sig.Prom.Peaks")],1,mean),
                               apply(brownmotif50promz[,c("Enrichment.in.W1.DownReg.Sig.Prom.Peaks",
                                                          "Enrichment.in.W2.DownReg.Sig.Prom.Peaks",
                                                          "Enrichment.in.W4.DownReg.Sig.Prom.Peaks",
                                                          "Enrichment.in.W8.DownReg.Sig.Prom.Peaks")],1,mean),
                               apply(whitemotif50promz[,c("Enrichment.in.W1.DownReg.Sig.Prom.Peaks",
                                                          "Enrichment.in.W2.DownReg.Sig.Prom.Peaks",
                                                          "Enrichment.in.W4.DownReg.Sig.Prom.Peaks",
                                                          "Enrichment.in.W8.DownReg.Sig.Prom.Peaks")],1,mean))
colnames(tfmotif50prommeanzmat) <- c("SKM-GN UP","HEART UP","HIPPOC UP","KIDNEY UP",
                                     "LIVER UP","LUNG UP","BAT UP","WAT-SC UP",
                                     "SKM-GN DOWN","HEART DOWN","HIPPOC DOWN","KIDNEY DOWN",
                                     "LIVER DOWN","LUNG DOWN","BAT DOWN","WAT-SC DOWN")

tfmotif50prommeanzmatupreg <- tfmotif50prommeanzmat[Reduce(union,upreg50list),]
tfmotif50prommeanzmatdownreg <- tfmotif50prommeanzmat[Reduce(union,downreg50list),]

promzmatmetadf <- data.frame(row.names = colnames(tfmotif50prommeanzmatupreg),
                             "Tissue" = gsub(" .*","",colnames(tfmotif50prommeanzmatupreg)))

# Figure 7A
pdf(file = "Figure 7A_021325.pdf",width = 5,height = 6)
pheatmap(tfmotif50prommeanzmatupreg[apply(abs(tfmotif50prommeanzmatupreg[,c(3,7,5,8,2,6,1,4)]),1,max) > 0.88,c(3,7,5,8,2,6,1,4)],show_rownames = T,cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = 0,labels_row = gsub("\\(.*","",rownames(tfmotif50prommeanzmatupreg[apply(abs(tfmotif50prommeanzmatupreg[,c(3,7,5,8,2,6,1,4)]),1,max) > 0.88,])),annotation_col = promzmatmetadf,annotation_colors = tissue_cols,show_colnames = F,cellwidth = 18,border_color = NA)
dev.off()

# Figure 7B
pdf(file = "Figure 7B_021325.pdf",width = 5,height = 6)
pheatmap(tfmotif50prommeanzmatdownreg[apply(abs(tfmotif50prommeanzmatdownreg[,c(11,15,10,16,9,12,13,14)]),1,max) > 0.88,c(11,15,10,16,9,12,13,14)],show_rownames = T,cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = 0,labels_row = gsub("\\(.*","",rownames(tfmotif50prommeanzmatdownreg[apply(abs(tfmotif50prommeanzmatdownreg[,c(11,15,10,16,9,12,13,14)]),1,max) > 0.88,])),annotation_col = promzmatmetadf,annotation_colors = tissue_cols,show_colnames = F,cellwidth = 18,border_color = NA)
dev.off()

# Figure 7A
png(file = "Figure 7A_021325.png",width = 5,height = 6,units = "in",res = 600)
pheatmap(tfmotif50prommeanzmatupreg[apply(abs(tfmotif50prommeanzmatupreg[,c(3,7,5,8,2,6,1,4)]),1,max) > 0.88,c(3,7,5,8,2,6,1,4)],show_rownames = T,cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = 0,labels_row = gsub("\\(.*","",rownames(tfmotif50prommeanzmatupreg[apply(abs(tfmotif50prommeanzmatupreg[,c(3,7,5,8,2,6,1,4)]),1,max) > 0.88,])),annotation_col = promzmatmetadf,annotation_colors = tissue_cols,show_colnames = F,cellwidth = 18,border_color = NA)
dev.off()

# Figure 7B
png(file = "Figure 7B_021325.png",width = 5,height = 6,units = "in",res = 600)
pheatmap(tfmotif50prommeanzmatdownreg[apply(abs(tfmotif50prommeanzmatdownreg[,c(11,15,10,16,9,12,13,14)]),1,max) > 0.88,c(11,15,10,16,9,12,13,14)],show_rownames = T,cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = 0,labels_row = gsub("\\(.*","",rownames(tfmotif50prommeanzmatdownreg[apply(abs(tfmotif50prommeanzmatdownreg[,c(11,15,10,16,9,12,13,14)]),1,max) > 0.88,])),annotation_col = promzmatmetadf,annotation_colors = tissue_cols,show_colnames = F,cellwidth = 18,border_color = NA)
dev.off()

tfmotifprommeanupregsign50 <- sign(tfmotif50prommeanzmatupreg[apply(abs(tfmotif50prommeanzmatupreg[,c(1:8)]),1,max) > 0.88,c(1:8)]) * (abs(tfmotif50prommeanzmatupreg[apply(abs(tfmotif50prommeanzmatupreg[,c(1:8)]),1,max) > 0.88,c(1:8)]) > 0.5)
tfmotifprommeandownregsign50 <- sign(tfmotif50prommeanzmatdownreg[apply(abs(tfmotif50prommeanzmatdownreg[,c(9:16)]),1,max) > 0.88,c(9:16)]) * (abs(tfmotif50prommeanzmatdownreg[apply(abs(tfmotif50prommeanzmatdownreg[,c(9:16)]),1,max) > 0.88,c(9:16)]) > 0.5)


prommeanupregtissuecomp50 <- matrix(0L,nrow = 8,ncol = 8)
rownames(prommeanupregtissuecomp50) <- c("SKM-GN","HEART","HIPPOC","KIDNEY",
                                         "LIVER","LUNG","BAT","WAT-SC")
colnames(prommeanupregtissuecomp50) <- c("SKM-GN","HEART","HIPPOC","KIDNEY",
                                         "LIVER","LUNG","BAT","WAT-SC")
for(i in 1:8){
  for(j in 1:8){
    prommeanupregtissuecomp50[i,j] <- sum(tfmotifprommeanupregsign50[,i] == 1 & tfmotifprommeanupregsign50[,j] == 1) #+ sum(tfmotifprommeanupregsign[,i] == -1 & tfmotifprommeanupregsign[,j] == -1)
  }
}


prommeandownregtissuecomp50 <- matrix(0L,nrow = 8,ncol = 8)
rownames(prommeandownregtissuecomp50) <- c("SKM-GN","HEART","HIPPOC","KIDNEY",
                                           "LIVER","LUNG","BAT","WAT-SC")
colnames(prommeandownregtissuecomp50) <- c("SKM-GN","HEART","HIPPOC","KIDNEY",
                                           "LIVER","LUNG","BAT","WAT-SC")
for(i in 1:8){
  for(j in 1:8){
    prommeandownregtissuecomp50[i,j] <- sum(tfmotifprommeandownregsign50[,i] == 1 & tfmotifprommeandownregsign50[,j] == 1) #+ sum(tfmotifprommeandownregsign[,i] == -1 & tfmotifprommeandownregsign[,j] == -1)
  }
}

# Figure 7C
pdf(file = "Figure 7C_021325.pdf",width = 5,height = 4)
pheatmap(prommeanupregtissuecomp50,display_numbers = T,angle_col = 0,color = colorpanel(101,"white","firebrick"),annotation_row = tissuemeta,annotation_col = tissuemeta,annotation_colors = ann_cols,show_rownames = F,show_colnames = F,fontsize = 15,number_color = "black",number_format = "%d",breaks = seq(0,30,length.out = 101),annotation_names_col = F,annotation_names_row = F,annotation_legend = F,border_color = NA)
dev.off()

png(file = "Figure 7C_021325.png",width = 5,height = 4,units = "in",res = 600)
pheatmap(prommeanupregtissuecomp50,display_numbers = T,angle_col = 0,color = colorpanel(101,"white","firebrick"),annotation_row = tissuemeta,annotation_col = tissuemeta,annotation_colors = ann_cols,show_rownames = F,show_colnames = F,fontsize = 15,number_color = "black",number_format = "%d",breaks = seq(0,30,length.out = 101),annotation_names_col = F,annotation_names_row = F,annotation_legend = F,border_color = NA)
dev.off()

# Figure 7D
pdf(file = "Figure 7D_021325.pdf",width = 5,height = 4)
pheatmap(prommeandownregtissuecomp50,display_numbers = T,angle_col = 0,color = colorpanel(101,"white","firebrick"),annotation_row = tissuemeta,annotation_col = tissuemeta,annotation_colors = ann_cols,show_rownames = F,show_colnames = F,fontsize = 15,number_color = "black",number_format = "%d",breaks = seq(0,30,length.out = 101),annotation_names_col = F,annotation_names_row = F,annotation_legend = F,border_color = NA)
dev.off()

png(file = "Figure 7D_021325.png",width = 5,height = 4,units = "in",res = 600)
pheatmap(prommeandownregtissuecomp50,display_numbers = T,angle_col = 0,color = colorpanel(101,"white","firebrick"),annotation_row = tissuemeta,annotation_col = tissuemeta,annotation_colors = ann_cols,show_rownames = F,show_colnames = F,fontsize = 15,number_color = "black",number_format = "%d",breaks = seq(0,30,length.out = 101),annotation_names_col = F,annotation_names_row = F,annotation_legend = F,border_color = NA)
dev.off()

save.image("Figure7_021325.RData")