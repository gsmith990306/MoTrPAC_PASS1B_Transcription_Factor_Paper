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

load(file = "PASS1B Transcription Factor Paper Data/DEA Analysis/pass1b-06_transcript-rna-seq_dea.RData")
load(file = "PASS1B Transcription Factor Paper Data/DEA Analysis/pass1b-06_epigen-atac-seq_dea.RData")

combodeglist <- c(gastrornasig,heartrnasig,hippornasig,kidneyrnasig,liverrnasig,lungrnasig,brownrnasig,whiternasig)
combodarlist <- c(gastroatacsig,heartatacsig,hippoatacsig,kidneyatacsig,liveratacsig,lungatacsig,brownatacsig,whiteatacsig)

combogenelist <- c(transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue_abbreviation %in% "SKM-GN","feature_ID"],
                   transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue_abbreviation %in% "HEART","feature_ID"],
                   transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue_abbreviation %in% "HIPPOC","feature_ID"],
                   transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue_abbreviation %in% "KIDNEY","feature_ID"],
                   transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue_abbreviation %in% "LIVER","feature_ID"],
                   transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue_abbreviation %in% "LUNG","feature_ID"],
                   transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue_abbreviation %in% "BAT","feature_ID"],
                   transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue_abbreviation %in% "WAT-SC","feature_ID"])

combopeaklist <- c(epigen_atac_seq$training_dea[epigen_atac_seq$training_dea$tissue_abbreviation %in% "SKM-GN","feature_ID"],
                   epigen_atac_seq$training_dea[epigen_atac_seq$training_dea$tissue_abbreviation %in% "HEART","feature_ID"],
                   epigen_atac_seq$training_dea[epigen_atac_seq$training_dea$tissue_abbreviation %in% "HIPPOC","feature_ID"],
                   epigen_atac_seq$training_dea[epigen_atac_seq$training_dea$tissue_abbreviation %in% "KIDNEY","feature_ID"],
                   epigen_atac_seq$training_dea[epigen_atac_seq$training_dea$tissue_abbreviation %in% "LIVER","feature_ID"],
                   epigen_atac_seq$training_dea[epigen_atac_seq$training_dea$tissue_abbreviation %in% "LUNG","feature_ID"],
                   epigen_atac_seq$training_dea[epigen_atac_seq$training_dea$tissue_abbreviation %in% "BAT","feature_ID"],
                   epigen_atac_seq$training_dea[epigen_atac_seq$training_dea$tissue_abbreviation %in% "WAT-SC","feature_ID"])

rm(transcript_rna_seq)
rm(epigen_atac_seq)
gc()

load("PASS1B Transcription Factor Paper Data/DEA Analysis/pass1b-06_epigen-rrbs_dea.RData")

combodmrlist <- c(gastromethsig,heartmethsig,hippomethsig,kidneymethsig,livermethsig,lungmethsig,brownmethsig,whitemethsig)

combomethlist <- c(epigen_rrbs$training_dea[epigen_rrbs$training_dea$tissue_abbreviation %in% "SKM-GN","feature_ID"],
                   epigen_rrbs$training_dea[epigen_rrbs$training_dea$tissue_abbreviation %in% "HEART","feature_ID"],
                   epigen_rrbs$training_dea[epigen_rrbs$training_dea$tissue_abbreviation %in% "HIPPOC","feature_ID"],
                   epigen_rrbs$training_dea[epigen_rrbs$training_dea$tissue_abbreviation %in% "KIDNEY","feature_ID"],
                   epigen_rrbs$training_dea[epigen_rrbs$training_dea$tissue_abbreviation %in% "LIVER","feature_ID"],
                   epigen_rrbs$training_dea[epigen_rrbs$training_dea$tissue_abbreviation %in% "LUNG","feature_ID"],
                   epigen_rrbs$training_dea[epigen_rrbs$training_dea$tissue_abbreviation %in% "BAT","feature_ID"],
                   epigen_rrbs$training_dea[epigen_rrbs$training_dea$tissue_abbreviation %in% "WAT-SC","feature_ID"])


methcomboquantmat <- cbind(c(table(table(combopeaklist))[1],
                             table(table(combopeaklist))[2],
                             table(table(combopeaklist))[3],
                             table(table(combopeaklist))[4],
                             table(table(combopeaklist))[5],
                             table(table(combopeaklist))[6],
                             table(table(combopeaklist))[7],
                             table(table(combopeaklist))[8]),
                           c(table(table(combogenelist))[1],
                             table(table(combogenelist))[2],
                             table(table(combogenelist))[3],
                             table(table(combogenelist))[4],
                             table(table(combogenelist))[5],
                             table(table(combogenelist))[6],
                             table(table(combogenelist))[7],
                             table(table(combogenelist))[8]),
                           c(table(table(combomethlist))[1],
                             table(table(combomethlist))[2],
                             table(table(combomethlist))[3],
                             table(table(combomethlist))[4],
                             table(table(combomethlist))[5],
                             table(table(combomethlist))[6],
                             table(table(combomethlist))[7],
                             table(table(combomethlist))[8]),
                           c(table(table(combodarlist))[1],
                             table(table(combodarlist))[2],
                             table(table(combodarlist))[3],
                             table(table(combodarlist))[4],
                             table(table(combodarlist))[5],
                             table(table(combodarlist))[6],
                             table(table(combodarlist))[7],
                             table(table(combodarlist))[8]),
                           c(table(table(combodeglist))[1],
                             table(table(combodeglist))[2],
                             table(table(combodeglist))[3],
                             table(table(combodeglist))[4],
                             table(table(combodeglist))[5],
                             table(table(combodeglist))[6],
                             table(table(combodeglist))[7],
                             table(table(combodeglist))[8]),
                           c(table(table(combodmrlist))[1],
                             table(table(combodmrlist))[2],
                             table(table(combodmrlist))[3],
                             table(table(combodmrlist))[4],
                             table(table(combodmrlist))[5],
                             table(table(combodmrlist))[6],
                             table(table(combodmrlist))[7],
                             table(table(combodmrlist))[8]))

methcomboquantmat[is.na(methcomboquantmat)] <- 0
colnames(methcomboquantmat) <- c("ATAC Peaks","Genes","Methyl Sites","DARs","DEGs","DMRs")
rownames(methcomboquantmat) <- c("Tissue Specific",
                                 "Shared by 2 tissues",
                                 "Shared by 3 tissues",
                                 "Shared by 4 tissues",
                                 "Shared by 5 tissues",
                                 "Shared by 6 tissues",
                                 "Shared by 7 tissues",
                                 "Shared by 8 tissues")

methcomboquantmatquo <- methcomboquantmat
for(i in 1:6){
  methcomboquantmatquo[,i] <- methcomboquantmat[,i]/sum(methcomboquantmat[,i])
}

methcomboquantdf <- data.frame("Percentage" = as.vector(methcomboquantmatquo),
                               "Specificity" = rep(c("Tissue Specific",
                                                     "Shared by 2 tissues",
                                                     "Shared by 3 tissues",
                                                     "Shared by 4 tissues",
                                                     "Shared by 5 tissues",
                                                     "Shared by 6 tissues",
                                                     "Shared by 7 tissues",
                                                     "Shared by 8 tissues"),6),
                               "Comparison" = c(rep("ATAC Peaks",8),
                                                rep("Genes",8),
                                                rep("Methyl Sites",8),
                                                rep("DARs",8),
                                                rep("DEGs",8),
                                                rep("DMRs",8)))
methcomboquantdf$Specificity <- factor(methcomboquantdf$Specificity,
                                       levels = c("Tissue Specific",
                                                  "Shared by 2 tissues",
                                                  "Shared by 3 tissues",
                                                  "Shared by 4 tissues",
                                                  "Shared by 5 tissues",
                                                  "Shared by 6 tissues",
                                                  "Shared by 7 tissues",
                                                  "Shared by 8 tissues"))
methcomboquantdf$Comparison <- factor(methcomboquantdf$Comparison,
                                      levels = c("ATAC Peaks",
                                                 "Genes",
                                                 "Methyl Sites",
                                                 "DARs",
                                                 "DEGs",
                                                 "DMRs"))

pdf(file = "Figure 1C.pdf", width=7, height=8)
ggplot(methcomboquantdf, aes(fill=Specificity, y=Percentage, x=Comparison)) + 
  geom_bar(position="fill", stat="identity") + scale_fill_brewer(palette = "Spectral") + theme_classic() + theme(text = element_text(size = 20),axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1)) + xlab(label = "")
dev.off()

png(file = "Figure 1C.png", width=7, height=8,units = "in",res = 600)
ggplot(methcomboquantdf, aes(fill=Specificity, y=Percentage, x=Comparison)) + 
  geom_bar(position="fill", stat="identity") + scale_fill_brewer(palette = "Spectral") + theme_classic() + theme(text = element_text(size = 20),axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1)) + xlab(label = "")
dev.off()

rm(epigen_rrbs)
gc()
save.image("Figure1C.RData")