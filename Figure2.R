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

#load("PASS1B Transcription Factor Paper Data/DEA Analysis/pass1b-06_transcript-rna-seq_dea.RData")
load("omesigdata.RData")

# Figure 2A

load("rnacontrolnormandmeta.RData")

combodeglist <- c(gastrornasig,heartrnasig,hippornasig,kidneyrnasig,liverrnasig,lungrnasig,brownrnasig,whiternasig)

degbaselineexpressionmat <- matrix(0L,nrow = length(combodeglist),ncol = 8)
rownames(degbaselineexpressionmat) <- combodeglist
colnames(degbaselineexpressionmat) <- c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC")

combodegmeta <- data.frame(row.names = unique(unique(combodeglist)),"SKM-GN" = rep(0,length(unique(combodeglist))),"HEART" = rep(0,length(unique(combodeglist))),
                           "HIPPOC" = rep(0,length(unique(combodeglist))),"KIDNEY" = rep(0,length(unique(combodeglist))),"LIVER" = rep(0,length(unique(combodeglist))),"LUNG" = rep(0,length(unique(combodeglist))),
                           "BAT" = rep(0,length(unique(combodeglist))),"WAT-SC" = rep(0,length(unique(combodeglist))))
for(i in 1:length(rownames(combodegmeta))){
  ourgene <- rownames(combodegmeta)[i]
  if(ourgene %in% gastrornasig){
    combodegmeta[i,"SKM.GN"] <- 1
  }
  if(ourgene %in% heartrnasig){
    combodegmeta[i,"HEART"] <- 1
  }
  if(ourgene %in% hippornasig){
    combodegmeta[i,"HIPPOC"] <- 1
  }
  if(ourgene %in% kidneyrnasig){
    combodegmeta[i,"KIDNEY"] <- 1
  }
  if(ourgene %in% liverrnasig){
    combodegmeta[i,"LIVER"] <- 1
  }
  if(ourgene %in% lungrnasig){
    combodegmeta[i,"LUNG"] <- 1
  }
  if(ourgene %in% brownrnasig){
    combodegmeta[i,"BAT"] <- 1
  }
  if(ourgene %in% whiternasig){
    combodegmeta[i,"WAT.SC"] <- 1
  }
}


rnacontrolnormmeans <- cbind(rowMeans(rnacontrolnorm[,c(61:70)]),
                             rowMeans(rnacontrolnorm[,c(71:80)]),
                             rowMeans(rnacontrolnorm[,c(81:90)]),
                             rowMeans(rnacontrolnorm[,c(101:110)]),
                             rowMeans(rnacontrolnorm[,c(111:120)]),
                             rowMeans(rnacontrolnorm[,c(121:130)]),
                             rowMeans(rnacontrolnorm[,c(31:40)]),
                             rowMeans(rnacontrolnorm[,c(171:180)]))
rownames(rnacontrolnormmeans) <- rownames(rnacontrolnorm)
colnames(rnacontrolnormmeans) <- c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC")

tissuemeta <- data.frame(row.names = colnames(rnacontrolnormmeans),
                         "Tissue" = colnames(rnacontrolnormmeans))
ann_new_cols <- list("Tissue" = c("SKM-GN" = "#088c03",
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
                                  "Upstream" = "#1F77B4FF"),
                     "SKM.GN" = c("0" = "grey90","1" = "#088c03"),
                     "HEART" = c("0" = "grey90","1" = "#f28b2f"),
                     "HIPPOC" = c("0" = "grey90","1" = "#bf7534"),
                     "KIDNEY" = c("0" = "grey90","1" = "#7553a7"),
                     "LIVER" = c("0" = "grey90","1" = "#da6c75"),
                     "LUNG" = c("0" = "grey90","1" = "#04bf8a"),
                     "BAT" = c("0" = "grey90","1" = "#8c5220"),
                     "WAT.SC" = c("0" = "grey90","1" = "#214da6"))


pdf(file = "Figure 2A.pdf",width = 5,height = 6)
pheatmap(t(scale(t(rnacontrolnormmeans[combodeglist[combodeglist %in% rownames(rnacontrolnorm)],]))),breaks = seq(-3,3,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_rows = F,cluster_cols = F,annotation_col = tissuemeta,show_rownames = F,show_colnames = F,annotation_row = combodegmeta,annotation_colors = ann_new_cols)
dev.off()

png(file = "Figure 2A.png",width = 5,height = 6,units = "in",res = 600)
pheatmap(t(scale(t(rnacontrolnormmeans[combodeglist[combodeglist %in% rownames(rnacontrolnorm)],]))),breaks = seq(-3,3,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_rows = F,cluster_cols = F,annotation_col = tissuemeta,show_rownames = F,show_colnames = F,annotation_row = combodegmeta,annotation_colors = ann_new_cols)
dev.off()

# Figure 2B

#load("PASS1B Transcription Factor Paper Data/DEA Analysis/pass1b-06_epigen-atac-seq_dea.RData")

atacseqcontrolnorm <- readRDS(file = "atacseqcontrolnorm.RDS")

ataccontrolnormmeans <- cbind(rowMeans(atacseqcontrolnorm[,c(11:20)]),
                              rowMeans(atacseqcontrolnorm[,c(21:30)]),
                              rowMeans(atacseqcontrolnorm[,c(31:40)]),
                              rowMeans(atacseqcontrolnorm[,c(41:50)]),
                              rowMeans(atacseqcontrolnorm[,c(51:60)]),
                              rowMeans(atacseqcontrolnorm[,c(61:70)]),
                              rowMeans(atacseqcontrolnorm[,c(1:10)]),
                              rowMeans(atacseqcontrolnorm[,c(71:80)]))
rownames(ataccontrolnormmeans) <- rownames(atacseqcontrolnorm)
colnames(ataccontrolnormmeans) <- c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC")

combodarlist <- c(gastroatacsig,heartatacsig,hippoatacsig,kidneyatacsig,liveratacsig,lungatacsig,brownatacsig,whiteatacsig)

combodarmeta <- data.frame(row.names = unique(unique(combodarlist)),"SKM-GN" = rep(0,length(unique(combodarlist))),"HEART" = rep(0,length(unique(combodarlist))),
                           "HIPPOC" = rep(0,length(unique(combodarlist))),"KIDNEY" = rep(0,length(unique(combodarlist))),"LIVER" = rep(0,length(unique(combodarlist))),"LUNG" = rep(0,length(unique(combodarlist))),
                           "BAT" = rep(0,length(unique(combodarlist))),"WAT-SC" = rep(0,length(unique(combodarlist))))
for(i in 1:length(rownames(combodarmeta))){
  ourgene <- rownames(combodarmeta)[i]
  if(ourgene %in% gastroatacsig){
    combodarmeta[i,"SKM.GN"] <- 1
  }
  if(ourgene %in% heartatacsig){
    combodarmeta[i,"HEART"] <- 1
  }
  if(ourgene %in% hippoatacsig){
    combodarmeta[i,"HIPPOC"] <- 1
  }
  if(ourgene %in% kidneyatacsig){
    combodarmeta[i,"KIDNEY"] <- 1
  }
  if(ourgene %in% liveratacsig){
    combodarmeta[i,"LIVER"] <- 1
  }
  if(ourgene %in% lungatacsig){
    combodarmeta[i,"LUNG"] <- 1
  }
  if(ourgene %in% brownatacsig){
    combodarmeta[i,"BAT"] <- 1
  }
  if(ourgene %in% whiteatacsig){
    combodarmeta[i,"WAT.SC"] <- 1
  }
}

pdf(file = "Figure 2B.pdf",width = 5,height = 6)
pheatmap(t(scale(t(ataccontrolnormmeans[combodarlist[combodarlist %in% rownames(ataccontrolnormmeans)],]))),breaks = seq(-3,3,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_rows = F,cluster_cols = F,annotation_col = tissuemeta,show_rownames = F,show_colnames = F,annotation_row = combodarmeta,annotation_colors = ann_new_cols)
dev.off()

png(file = "Figure 2B.png",width = 5,height = 6,units = "in",res = 600)
pheatmap(t(scale(t(ataccontrolnormmeans[combodarlist[combodarlist %in% rownames(ataccontrolnormmeans)],]))),breaks = seq(-3,3,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_rows = F,cluster_cols = F,annotation_col = tissuemeta,show_rownames = F,show_colnames = F,annotation_row = combodarmeta,annotation_colors = ann_new_cols)
dev.off()

rm(transcript_rna_seq)
rm(epigen_atac_seq)
gc()

# Figure 2C

pass1bphenodata <- readRDS("pass1bphenodata.rds")

gastromethM <- read.delim(file = "PASS1B Transcription Factor Paper Data/RRBS Data/gastrocnemius_M_matrix.txt",header = TRUE, row.names = 1,sep = "\t")
colnames(gastromethM) <- gsub("X","",colnames(gastromethM))

gastromethmeta <- data.frame(row.names = colnames(gastromethM),
                             "Sex" = rep("Female",length(colnames(gastromethM))),
                             "Group" = rep("",length(colnames(gastromethM))))
for(i in 1:length(colnames(gastromethM))){
  ourid <- colnames(gastromethM)[i]
  if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"] == 2){
    gastromethmeta[i,"Sex"] <- "Male"
  }
  gastromethmeta[i,"Group"] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"]
}
gastromethmeta$Cohort <- paste(gastromethmeta$Sex,gastromethmeta$Group,sep = "_")

gastromethmeta[gastromethmeta$Group %in% "Eight-week program Control Group","Group"] <- "control"
gastromethmeta[gastromethmeta$Group %in% "One-week program","Group"] <- "1w"
gastromethmeta[gastromethmeta$Group %in% "Two-week program","Group"] <- "2w"
gastromethmeta[gastromethmeta$Group %in% "Four-week program","Group"] <- "4w"
gastromethmeta[gastromethmeta$Group %in% "Eight-week program Training Group","Group"] <- "8w"

gastromethmeta$is_outlier <- t(gastromethM["is_outlier",rownames(gastromethmeta)])

heartmethM <- read.delim(file = "PASS1B Transcription Factor Paper Data/RRBS Data/heart_M_matrix.txt",header = TRUE, row.names = 1,sep = "\t")
colnames(heartmethM) <- gsub("X","",colnames(heartmethM))

heartmethmeta <- data.frame(row.names = colnames(heartmethM),
                            "Sex" = rep("Female",length(colnames(heartmethM))),
                            "Group" = rep("",length(colnames(heartmethM))))
for(i in 1:length(colnames(heartmethM))){
  ourid <- colnames(heartmethM)[i]
  if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"] == 2){
    heartmethmeta[i,"Sex"] <- "Male"
  }
  heartmethmeta[i,"Group"] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"]
}
heartmethmeta$Cohort <- paste(heartmethmeta$Sex,heartmethmeta$Group,sep = "_")

heartmethmeta[heartmethmeta$Group %in% "Eight-week program Control Group","Group"] <- "control"
heartmethmeta[heartmethmeta$Group %in% "One-week program","Group"] <- "1w"
heartmethmeta[heartmethmeta$Group %in% "Two-week program","Group"] <- "2w"
heartmethmeta[heartmethmeta$Group %in% "Four-week program","Group"] <- "4w"
heartmethmeta[heartmethmeta$Group %in% "Eight-week program Training Group","Group"] <- "8w"

heartmethmeta$is_outlier <- t(heartmethM["is_outlier",rownames(heartmethmeta)])

hippomethM <- read.delim(file = "PASS1B Transcription Factor Paper Data/RRBS Data/hippocampus_M_matrix.txt",header = TRUE, row.names = 1,sep = "\t")
colnames(hippomethM) <- gsub("X","",colnames(hippomethM))

hippomethmeta <- data.frame(row.names = colnames(hippomethM),
                            "Sex" = rep("Female",length(colnames(hippomethM))),
                            "Group" = rep("",length(colnames(hippomethM))))
for(i in 1:length(colnames(hippomethM))){
  ourid <- colnames(hippomethM)[i]
  if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"] == 2){
    hippomethmeta[i,"Sex"] <- "Male"
  }
  hippomethmeta[i,"Group"] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"]
}
hippomethmeta$Cohort <- paste(hippomethmeta$Sex,hippomethmeta$Group,sep = "_")

hippomethmeta[hippomethmeta$Group %in% "Eight-week program Control Group","Group"] <- "control"
hippomethmeta[hippomethmeta$Group %in% "One-week program","Group"] <- "1w"
hippomethmeta[hippomethmeta$Group %in% "Two-week program","Group"] <- "2w"
hippomethmeta[hippomethmeta$Group %in% "Four-week program","Group"] <- "4w"
hippomethmeta[hippomethmeta$Group %in% "Eight-week program Training Group","Group"] <- "8w"

hippomethmeta$is_outlier <- t(hippomethM["is_outlier",rownames(hippomethmeta)])

kidneymethM <- read.delim(file = "PASS1B Transcription Factor Paper Data/RRBS Data/kidney_M_matrix.txt",header = TRUE, row.names = 1,sep = "\t")
colnames(kidneymethM) <- gsub("X","",colnames(kidneymethM))

kidneymethmeta <- data.frame(row.names = colnames(kidneymethM),
                             "Sex" = rep("Female",length(colnames(kidneymethM))),
                             "Group" = rep("",length(colnames(kidneymethM))))
for(i in 1:length(colnames(kidneymethM))){
  ourid <- colnames(kidneymethM)[i]
  if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"] == 2){
    kidneymethmeta[i,"Sex"] <- "Male"
  }
  kidneymethmeta[i,"Group"] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"]
}
kidneymethmeta$Cohort <- paste(kidneymethmeta$Sex,kidneymethmeta$Group,sep = "_")

kidneymethmeta[kidneymethmeta$Group %in% "Eight-week program Control Group","Group"] <- "control"
kidneymethmeta[kidneymethmeta$Group %in% "One-week program","Group"] <- "1w"
kidneymethmeta[kidneymethmeta$Group %in% "Two-week program","Group"] <- "2w"
kidneymethmeta[kidneymethmeta$Group %in% "Four-week program","Group"] <- "4w"
kidneymethmeta[kidneymethmeta$Group %in% "Eight-week program Training Group","Group"] <- "8w"

kidneymethmeta$is_outlier <- t(kidneymethM["is_outlier",rownames(kidneymethmeta)])

livermethM <- read.delim(file = "PASS1B Transcription Factor Paper Data/RRBS Data/liver_M_matrix.txt",header = TRUE, row.names = 1,sep = "\t")
colnames(livermethM) <- gsub("X","",colnames(livermethM))

livermethmeta <- data.frame(row.names = colnames(livermethM),
                            "Sex" = rep("Female",length(colnames(livermethM))),
                            "Group" = rep("",length(colnames(livermethM))))
for(i in 1:length(colnames(livermethM))){
  ourid <- colnames(livermethM)[i]
  if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"] == 2){
    livermethmeta[i,"Sex"] <- "Male"
  }
  livermethmeta[i,"Group"] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"]
}
livermethmeta$Cohort <- paste(livermethmeta$Sex,livermethmeta$Group,sep = "_")

livermethmeta[livermethmeta$Group %in% "Eight-week program Control Group","Group"] <- "control"
livermethmeta[livermethmeta$Group %in% "One-week program","Group"] <- "1w"
livermethmeta[livermethmeta$Group %in% "Two-week program","Group"] <- "2w"
livermethmeta[livermethmeta$Group %in% "Four-week program","Group"] <- "4w"
livermethmeta[livermethmeta$Group %in% "Eight-week program Training Group","Group"] <- "8w"

livermethmeta$is_outlier <- t(livermethM["is_outlier",rownames(livermethmeta)])

lungmethM <- read.delim(file = "PASS1B Transcription Factor Paper Data/RRBS Data/lung_M_matrix.txt",header = TRUE, row.names = 1,sep = "\t")
colnames(lungmethM) <- gsub("X","",colnames(lungmethM))

lungmethmeta <- data.frame(row.names = colnames(lungmethM),
                           "Sex" = rep("Female",length(colnames(lungmethM))),
                           "Group" = rep("",length(colnames(lungmethM))))
for(i in 1:length(colnames(lungmethM))){
  ourid <- colnames(lungmethM)[i]
  if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"] == 2){
    lungmethmeta[i,"Sex"] <- "Male"
  }
  lungmethmeta[i,"Group"] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"]
}
lungmethmeta$Cohort <- paste(lungmethmeta$Sex,lungmethmeta$Group,sep = "_")

lungmethmeta[lungmethmeta$Group %in% "Eight-week program Control Group","Group"] <- "control"
lungmethmeta[lungmethmeta$Group %in% "One-week program","Group"] <- "1w"
lungmethmeta[lungmethmeta$Group %in% "Two-week program","Group"] <- "2w"
lungmethmeta[lungmethmeta$Group %in% "Four-week program","Group"] <- "4w"
lungmethmeta[lungmethmeta$Group %in% "Eight-week program Training Group","Group"] <- "8w"

lungmethmeta$is_outlier <- t(lungmethM["is_outlier",rownames(lungmethmeta)])

brownmethM <- read.delim(file = "PASS1B Transcription Factor Paper Data/RRBS Data/brown-adipose_M_matrix.txt",header = TRUE, row.names = 1,sep = "\t")
colnames(brownmethM) <- gsub("X","",colnames(brownmethM))

brownmethmeta <- data.frame(row.names = colnames(brownmethM),
                            "Sex" = rep("Female",length(colnames(brownmethM))),
                            "Group" = rep("",length(colnames(brownmethM))))
for(i in 1:length(colnames(brownmethM))){
  ourid <- colnames(brownmethM)[i]
  if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"] == 2){
    brownmethmeta[i,"Sex"] <- "Male"
  }
  brownmethmeta[i,"Group"] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"]
}
brownmethmeta$Cohort <- paste(brownmethmeta$Sex,brownmethmeta$Group,sep = "_")

brownmethmeta[brownmethmeta$Group %in% "Eight-week program Control Group","Group"] <- "control"
brownmethmeta[brownmethmeta$Group %in% "One-week program","Group"] <- "1w"
brownmethmeta[brownmethmeta$Group %in% "Two-week program","Group"] <- "2w"
brownmethmeta[brownmethmeta$Group %in% "Four-week program","Group"] <- "4w"
brownmethmeta[brownmethmeta$Group %in% "Eight-week program Training Group","Group"] <- "8w"

brownmethmeta$is_outlier <- t(brownmethM["is_outlier",rownames(brownmethmeta)])

whitemethM <- read.delim(file = "PASS1B Transcription Factor Paper Data/RRBS Data/white-adipose_M_matrix.txt",header = TRUE, row.names = 1,sep = "\t")
colnames(whitemethM) <- gsub("X","",colnames(whitemethM))

whitemethmeta <- data.frame(row.names = colnames(whitemethM),
                            "Sex" = rep("Female",length(colnames(whitemethM))),
                            "Group" = rep("",length(colnames(whitemethM))))
for(i in 1:length(colnames(whitemethM))){
  ourid <- colnames(whitemethM)[i]
  if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"] == 2){
    whitemethmeta[i,"Sex"] <- "Male"
  }
  whitemethmeta[i,"Group"] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"]
}
whitemethmeta$Cohort <- paste(whitemethmeta$Sex,whitemethmeta$Group,sep = "_")

whitemethmeta[whitemethmeta$Group %in% "Eight-week program Control Group","Group"] <- "control"
whitemethmeta[whitemethmeta$Group %in% "One-week program","Group"] <- "1w"
whitemethmeta[whitemethmeta$Group %in% "Two-week program","Group"] <- "2w"
whitemethmeta[whitemethmeta$Group %in% "Four-week program","Group"] <- "4w"
whitemethmeta[whitemethmeta$Group %in% "Eight-week program Training Group","Group"] <- "8w"

whitemethmeta$is_outlier <- t(whitemethM["is_outlier",rownames(whitemethmeta)])

combodmrlist <- c(gastromethsig,heartmethsig,hippomethsig,kidneymethsig,livermethsig,lungmethsig,brownmethsig,whitemethsig)

dmrbaselineexpressionmat <- matrix(0L,nrow = length(combodmrlist),ncol = 8)
rownames(dmrbaselineexpressionmat) <- combodmrlist
colnames(dmrbaselineexpressionmat) <- c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC")

for(i in 1:length(combodmrlist)){
  if(combodmrlist[i] %in% rownames(gastromethM)){
    dmrbaselineexpressionmat[i,"SKM-GN"] <- mean(as.numeric(gastromethM[combodmrlist[i],rownames(gastromethmeta[gastromethmeta$Group %in% "control" & gastromethmeta$is_outlier == "FALSE",])]))
  }
  if(combodmrlist[i] %in% rownames(heartmethM)){
    dmrbaselineexpressionmat[i,"HEART"] <- mean(as.numeric(heartmethM[combodmrlist[i],rownames(heartmethmeta[heartmethmeta$Group %in% "control" & heartmethmeta$is_outlier == "FALSE",])]))
  }
  if(combodmrlist[i] %in% rownames(hippomethM)){
    dmrbaselineexpressionmat[i,"HIPPOC"] <- mean(as.numeric(hippomethM[combodmrlist[i],rownames(hippomethmeta[hippomethmeta$Group %in% "control" & hippomethmeta$is_outlier == "FALSE",])]))
  }
  if(combodmrlist[i] %in% rownames(kidneymethM)){
    dmrbaselineexpressionmat[i,"KIDNEY"] <- mean(as.numeric(kidneymethM[combodmrlist[i],rownames(kidneymethmeta[kidneymethmeta$Group %in% "control" & kidneymethmeta$is_outlier == "FALSE",])]))
  }
  if(combodmrlist[i] %in% rownames(livermethM)){
    dmrbaselineexpressionmat[i,"LIVER"] <- mean(as.numeric(livermethM[combodmrlist[i],rownames(livermethmeta[livermethmeta$Group %in% "control" & livermethmeta$is_outlier == "FALSE",])]))
  }
  if(combodmrlist[i] %in% rownames(lungmethM)){
    dmrbaselineexpressionmat[i,"LUNG"] <- mean(as.numeric(lungmethM[combodmrlist[i],rownames(lungmethmeta[lungmethmeta$Group %in% "control" & lungmethmeta$is_outlier == "FALSE",])]))
  }
  if(combodmrlist[i] %in% rownames(brownmethM)){
    dmrbaselineexpressionmat[i,"BAT"] <- mean(as.numeric(brownmethM[combodmrlist[i],rownames(brownmethmeta[brownmethmeta$Group %in% "control" & brownmethmeta$is_outlier == "FALSE",])]))
  }
  if(combodmrlist[i] %in% rownames(whitemethM)){
    dmrbaselineexpressionmat[i,"WAT-SC"] <- mean(as.numeric(whitemethM[combodmrlist[i],rownames(whitemethmeta[whitemethmeta$Group %in% "control" & whitemethmeta$is_outlier == "FALSE",])]))
  }
}

combodmrmeta <- data.frame(row.names = unique(unique(combodmrlist)),"SKM-GN" = rep(0,length(unique(combodmrlist))),"HEART" = rep(0,length(unique(combodmrlist))),
                           "HIPPOC" = rep(0,length(unique(combodmrlist))),"KIDNEY" = rep(0,length(unique(combodmrlist))),"LIVER" = rep(0,length(unique(combodmrlist))),"LUNG" = rep(0,length(unique(combodmrlist))),
                           "BAT" = rep(0,length(unique(combodmrlist))),"WAT-SC" = rep(0,length(unique(combodmrlist))))
for(i in 1:length(rownames(combodmrmeta))){
  ourdmr <- rownames(combodmrmeta)[i]
  if(ourdmr %in% gastromethsig){
    combodmrmeta[i,"SKM.GN"] <- 1
  }
  if(ourdmr %in% heartmethsig){
    combodmrmeta[i,"HEART"] <- 1
  }
  if(ourdmr %in% hippomethsig){
    combodmrmeta[i,"HIPPOC"] <- 1
  }
  if(ourdmr %in% kidneymethsig){
    combodmrmeta[i,"KIDNEY"] <- 1
  }
  if(ourdmr %in% livermethsig){
    combodmrmeta[i,"LIVER"] <- 1
  }
  if(ourdmr %in% lungmethsig){
    combodmrmeta[i,"LUNG"] <- 1
  }
  if(ourdmr %in% brownmethsig){
    combodmrmeta[i,"BAT"] <- 1
  }
  if(ourdmr %in% whitemethsig){
    combodmrmeta[i,"WAT.SC"] <- 1
  }
}


pdf(file = "Figure 2C.pdf",width = 5,height = 6)
pheatmap(t(scale(t(dmrbaselineexpressionmat))),show_rownames = F,annotation_col = tissuemeta,show_colnames = F,breaks = seq(-3,3,length.out = 101),color = colorpanel(101,"blue","white","red"),annotation_row = combodmrmeta,cluster_cols = F,cluster_rows = F,annotation_colors = ann_new_cols)
dev.off()

png(file = "Figure 2C.png",width = 5,height = 6,units = "in",res = 600)
pheatmap(t(scale(t(dmrbaselineexpressionmat))),show_rownames = F,annotation_col = tissuemeta,show_colnames = F,breaks = seq(-3,3,length.out = 101),color = colorpanel(101,"blue","white","red"),annotation_row = combodmrmeta,cluster_cols = F,cluster_rows = F,annotation_colors = ann_new_cols)
dev.off()

# Figure 2D

rnasigtissuez <- t(scale(t(rnacontrolnormmeans[combodeglist[combodeglist %in% rownames(rnacontrolnorm)],])))
atacsigtissuez <- t(scale(t(ataccontrolnormmeans[combodarlist[combodarlist %in% rownames(ataccontrolnormmeans)],])))
methsigtissuez <- t(scale(t(dmrbaselineexpressionmat)))

rnasiglineplotdf <- data.frame("Control.Expression.Z.Score" = c(colMeans(rnasigtissuez[gastrornasig,]),
                                                                colMeans(rnasigtissuez[heartrnasig,]),
                                                                colMeans(rnasigtissuez[hippornasig,]),
                                                                colMeans(rnasigtissuez[kidneyrnasig,]),
                                                                colMeans(rnasigtissuez[liverrnasig,]),
                                                                colMeans(rnasigtissuez[lungrnasig,]),
                                                                colMeans(rnasigtissuez[brownrnasig,]),
                                                                colMeans(rnasigtissuez[intersect(whiternasig,rownames(rnasigtissuez)),])),
                               "Control.Tissue" = rep(c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC"),8),
                               "Tissue.DEGs" = c(rep("SKM-GN",8),rep("HEART",8),rep("HIPPOC",8),rep("KIDNEY",8),
                                                 rep("LIVER",8),rep("LUNG",8),rep("BAT",8),rep("WAT-SC",8)))

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


png("Figure 2D.png",width = 7,height = 3,units = "in",res = 600)
ggline(rnasiglineplotdf, "Control.Tissue", "Control.Expression.Z.Score",
       color = "Tissue.DEGs", palette = ann_cols$Tissue,size = 1)
dev.off()

pdf("Figure 2D.pdf",width = 7,height = 3)
ggline(rnasiglineplotdf, "Control.Tissue", "Control.Expression.Z.Score",
       color = "Tissue.DEGs", palette = ann_cols$Tissue,size = 1)
dev.off()

# Figure 2E

atacsiglineplotdf <- data.frame("Control.Expression.Z.Score" = c(colMeans(atacsigtissuez[intersect(gastroatacsig,rownames(atacsigtissuez)),]),
                                                                 colMeans(atacsigtissuez[intersect(heartatacsig,rownames(atacsigtissuez)),]),
                                                                 colMeans(atacsigtissuez[hippoatacsig,]),
                                                                 colMeans(atacsigtissuez[kidneyatacsig,]),
                                                                 colMeans(atacsigtissuez[liveratacsig,]),
                                                                 colMeans(atacsigtissuez[lungatacsig,]),
                                                                 colMeans(atacsigtissuez[brownatacsig,]),
                                                                 colMeans(atacsigtissuez[intersect(whiteatacsig,rownames(atacsigtissuez)),])),
                                "Control.Tissue" = rep(c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC"),8),
                                "Tissue.DARs" = c(rep("SKM-GN",8),rep("HEART",8),rep("HIPPOC",8),rep("KIDNEY",8),
                                                  rep("LIVER",8),rep("LUNG",8),rep("BAT",8),rep("WAT-SC",8)))

png("Figure 2E.png",width = 7,height = 3,units = "in",res = 600)
ggline(atacsiglineplotdf, "Control.Tissue", "Control.Expression.Z.Score",
       color = "Tissue.DARs", palette = ann_cols$Tissue,size = 1)
dev.off()

pdf("Figure 2E.pdf",width = 7,height = 3)
ggline(atacsiglineplotdf, "Control.Tissue", "Control.Expression.Z.Score",
       color = "Tissue.DARs", palette = ann_cols$Tissue,size = 1)
dev.off()

# Figure 2F

methsiglineplotdf <- data.frame("Control.Expression.Z.Score" = c(colMeans(methsigtissuez[intersect(gastromethsig,rownames(methsigtissuez)),]),
                                                                 colMeans(methsigtissuez[intersect(heartmethsig,rownames(methsigtissuez)),]),
                                                                 colMeans(methsigtissuez[hippomethsig,]),
                                                                 colMeans(methsigtissuez[kidneymethsig,]),
                                                                 colMeans(methsigtissuez[livermethsig,]),
                                                                 colMeans(methsigtissuez[lungmethsig,]),
                                                                 colMeans(methsigtissuez[brownmethsig,]),
                                                                 colMeans(methsigtissuez[intersect(whitemethsig,rownames(methsigtissuez)),])),
                                "Control.Tissue" = rep(c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC"),8),
                                "Tissue.DMRs" = c(rep("SKM-GN",8),rep("HEART",8),rep("HIPPOC",8),rep("KIDNEY",8),
                                                  rep("LIVER",8),rep("LUNG",8),rep("BAT",8),rep("WAT-SC",8)))


png("Figure 2F.png",width = 7,height = 3,units = "in",res = 600)
ggline(methsiglineplotdf, "Control.Tissue", "Control.Expression.Z.Score",
       color = "Tissue.DMRs", palette = ann_cols$Tissue,size = 1)
dev.off()

pdf("Figure 2F.pdf",width = 7,height = 3)
ggline(methsiglineplotdf, "Control.Tissue", "Control.Expression.Z.Score",
       color = "Tissue.DMRs", palette = ann_cols$Tissue,size = 1)
dev.off()

save.image("Figure2.RData")