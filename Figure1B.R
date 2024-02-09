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

#####
# Figure 1B
####

# Here we measure the frequency of differential analytes in the RNAseq, ATACseq and RRBS data

# To be edited to your directory containing the folder with all necessary paper data files
setwd("~/Mount Sinai/Sealfon Laboratory/MoTrPAC/PASS1B Raw Data")

load(file = "PASS1B Transcription Factor Paper Data/DEA Analysis/pass1b-06_transcript-rna-seq_dea.RData")
load(file = "PASS1B Transcription Factor Paper Data/DEA Analysis/pass1b-06_epigen-atac-seq_dea.RData")

gastrornasig <- unique(transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue_abbreviation %in% "SKM-GN" &
                                                         transcript_rna_seq$training_dea$adj_p_value < 0.1,"feature_ID"])
heartrnasig <- unique(transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue_abbreviation %in% "HEART" &
                                                        transcript_rna_seq$training_dea$adj_p_value < 0.1,"feature_ID"])
hippornasig <- unique(transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue_abbreviation %in% "HIPPOC" &
                                                        transcript_rna_seq$training_dea$adj_p_value < 0.1,"feature_ID"])
kidneyrnasig <- unique(transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue_abbreviation %in% "KIDNEY" &
                                                         transcript_rna_seq$training_dea$adj_p_value < 0.1,"feature_ID"])
liverrnasig <- unique(transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue_abbreviation %in% "LIVER" &
                                                        transcript_rna_seq$training_dea$adj_p_value < 0.1,"feature_ID"])
lungrnasig <- unique(transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue_abbreviation %in% "LUNG" &
                                                       transcript_rna_seq$training_dea$adj_p_value < 0.1,"feature_ID"])
brownrnasig <- unique(transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue_abbreviation %in% "BAT" &
                                                        transcript_rna_seq$training_dea$adj_p_value < 0.1,"feature_ID"])
whiternasig <- unique(transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue_abbreviation %in% "WAT-SC" &
                                                        transcript_rna_seq$training_dea$adj_p_value < 0.1,"feature_ID"])

gastroatacsig <- unique(epigen_atac_seq$training_dea[epigen_atac_seq$training_dea$tissue_abbreviation %in% "SKM-GN" &
                                                       epigen_atac_seq$training_dea$adj_p_value < 0.1,"feature_ID"])
heartatacsig <- unique(epigen_atac_seq$training_dea[epigen_atac_seq$training_dea$tissue_abbreviation %in% "HEART" &
                                                      epigen_atac_seq$training_dea$adj_p_value < 0.1,"feature_ID"])
hippoatacsig <- unique(epigen_atac_seq$training_dea[epigen_atac_seq$training_dea$tissue_abbreviation %in% "HIPPOC" &
                                                      epigen_atac_seq$training_dea$adj_p_value < 0.1,"feature_ID"])
kidneyatacsig <- unique(epigen_atac_seq$training_dea[epigen_atac_seq$training_dea$tissue_abbreviation %in% "KIDNEY" &
                                                       epigen_atac_seq$training_dea$adj_p_value < 0.1,"feature_ID"])
liveratacsig <- unique(epigen_atac_seq$training_dea[epigen_atac_seq$training_dea$tissue_abbreviation %in% "LIVER" &
                                                      epigen_atac_seq$training_dea$adj_p_value < 0.1,"feature_ID"])
lungatacsig <- unique(epigen_atac_seq$training_dea[epigen_atac_seq$training_dea$tissue_abbreviation %in% "LUNG" &
                                                     epigen_atac_seq$training_dea$adj_p_value < 0.1,"feature_ID"])
brownatacsig <- unique(epigen_atac_seq$training_dea[epigen_atac_seq$training_dea$tissue_abbreviation %in% "BAT" &
                                                      epigen_atac_seq$training_dea$adj_p_value < 0.1,"feature_ID"])
whiteatacsig <- unique(epigen_atac_seq$training_dea[epigen_atac_seq$training_dea$tissue_abbreviation %in% "WAT-SC" &
                                                      epigen_atac_seq$training_dea$adj_p_value < 0.1,"feature_ID"])

tissuetotalcountmat <- rbind(c(length(unique(transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue_abbreviation %in% "SKM-GN","feature_ID"])),
                               length(unique(transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue_abbreviation %in% "HEART","feature_ID"])),
                               length(unique(transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue_abbreviation %in% "HIPPOC","feature_ID"])),
                               length(unique(transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue_abbreviation %in% "KIDNEY","feature_ID"])),
                               length(unique(transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue_abbreviation %in% "LIVER","feature_ID"])),
                               length(unique(transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue_abbreviation %in% "LUNG","feature_ID"])),
                               length(unique(transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue_abbreviation %in% "BAT","feature_ID"])),
                               length(unique(transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue_abbreviation %in% "WAT-SC","feature_ID"]))),
                             c(length(unique(epigen_atac_seq$training_dea[epigen_atac_seq$training_dea$tissue_abbreviation %in% "SKM-GN","feature_ID"])),
                               length(unique(epigen_atac_seq$training_dea[epigen_atac_seq$training_dea$tissue_abbreviation %in% "HEART","feature_ID"])),
                               length(unique(epigen_atac_seq$training_dea[epigen_atac_seq$training_dea$tissue_abbreviation %in% "HIPPOC","feature_ID"])),
                               length(unique(epigen_atac_seq$training_dea[epigen_atac_seq$training_dea$tissue_abbreviation %in% "KIDNEY","feature_ID"])),
                               length(unique(epigen_atac_seq$training_dea[epigen_atac_seq$training_dea$tissue_abbreviation %in% "LIVER","feature_ID"])),
                               length(unique(epigen_atac_seq$training_dea[epigen_atac_seq$training_dea$tissue_abbreviation %in% "LUNG","feature_ID"])),
                               length(unique(epigen_atac_seq$training_dea[epigen_atac_seq$training_dea$tissue_abbreviation %in% "BAT","feature_ID"])),
                               length(unique(epigen_atac_seq$training_dea[epigen_atac_seq$training_dea$tissue_abbreviation %in% "WAT-SC","feature_ID"]))))
rownames(tissuetotalcountmat) <- c("DEGs","DARs")

colnames(tissuetotalcountmat) <- c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC")

tissuemeta <- data.frame(row.names = colnames(tissuetotalcountmat),
                         "Tissue" = colnames(tissuetotalcountmat))

rm(transcript_rna_seq)
rm(epigen_atac_seq)
gc()

load("PASS1B Transcription Factor Paper Data/DEA Analysis/pass1b-06_epigen-rrbs_dea.RData")

gastromethsig <- epigen_rrbs$training_dea[epigen_rrbs$training_dea$adj_p_value < 0.1 & epigen_rrbs$training_dea$tissue_abbreviation %in% "SKM-GN","feature_ID"]
heartmethsig <- epigen_rrbs$training_dea[epigen_rrbs$training_dea$adj_p_value < 0.1 & epigen_rrbs$training_dea$tissue_abbreviation %in% "HEART","feature_ID"]
hippomethsig <- epigen_rrbs$training_dea[epigen_rrbs$training_dea$adj_p_value < 0.1 & epigen_rrbs$training_dea$tissue_abbreviation %in% "HIPPOC","feature_ID"]
kidneymethsig <- epigen_rrbs$training_dea[epigen_rrbs$training_dea$adj_p_value < 0.1 & epigen_rrbs$training_dea$tissue_abbreviation %in% "KIDNEY","feature_ID"]
livermethsig <- epigen_rrbs$training_dea[epigen_rrbs$training_dea$adj_p_value < 0.1 & epigen_rrbs$training_dea$tissue_abbreviation %in% "LIVER","feature_ID"]
lungmethsig <- epigen_rrbs$training_dea[epigen_rrbs$training_dea$adj_p_value < 0.1 & epigen_rrbs$training_dea$tissue_abbreviation %in% "LUNG","feature_ID"]
brownmethsig <- epigen_rrbs$training_dea[epigen_rrbs$training_dea$adj_p_value < 0.1 & epigen_rrbs$training_dea$tissue_abbreviation %in% "BAT","feature_ID"]
whitemethsig <- epigen_rrbs$training_dea[epigen_rrbs$training_dea$adj_p_value < 0.1 & epigen_rrbs$training_dea$tissue_abbreviation %in% "WAT-SC","feature_ID"]


methtissuesigcountmat <- rbind(c(length(gastrornasig),
                                 length(heartrnasig),
                                 length(hippornasig),
                                 length(kidneyrnasig),
                                 length(liverrnasig),
                                 length(lungrnasig),
                                 length(brownrnasig),
                                 length(whiternasig)),
                               c(length(gastroatacsig),
                                 length(heartatacsig),
                                 length(hippoatacsig),
                                 length(kidneyatacsig),
                                 length(liveratacsig),
                                 length(lungatacsig),
                                 length(brownatacsig),
                                 length(whiteatacsig)),
                               c(length(gastromethsig),
                                 length(heartmethsig),
                                 length(hippomethsig),
                                 length(kidneymethsig),
                                 length(livermethsig),
                                 length(lungmethsig),
                                 length(brownmethsig),
                                 length(whitemethsig)))

methtissuetotalcountmat <- rbind(tissuetotalcountmat,
                                 c(length(unique(epigen_rrbs$training_dea[epigen_rrbs$training_dea$tissue_abbreviation %in% "SKM-GN","feature_ID"])),
                                   length(unique(epigen_rrbs$training_dea[epigen_rrbs$training_dea$tissue_abbreviation %in% "HEART","feature_ID"])),
                                   length(unique(epigen_rrbs$training_dea[epigen_rrbs$training_dea$tissue_abbreviation %in% "HIPPOC","feature_ID"])),
                                   length(unique(epigen_rrbs$training_dea[epigen_rrbs$training_dea$tissue_abbreviation %in% "KIDNEY","feature_ID"])),
                                   length(unique(epigen_rrbs$training_dea[epigen_rrbs$training_dea$tissue_abbreviation %in% "LIVER","feature_ID"])),
                                   length(unique(epigen_rrbs$training_dea[epigen_rrbs$training_dea$tissue_abbreviation %in% "LUNG","feature_ID"])),
                                   length(unique(epigen_rrbs$training_dea[epigen_rrbs$training_dea$tissue_abbreviation %in% "BAT","feature_ID"])),
                                   length(unique(epigen_rrbs$training_dea[epigen_rrbs$training_dea$tissue_abbreviation %in% "WAT-SC","feature_ID"]))))

rownames(methtissuesigcountmat) <- c("DEGs","DARs","DMRs")
rownames(methtissuetotalcountmat) <- c("DEGs","DARs","DMRs")

colnames(methtissuesigcountmat) <- c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC")
colnames(methtissuetotalcountmat) <- c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC")

methtissuesigquomat <- methtissuesigcountmat/methtissuetotalcountmat

tissuemeta <- data.frame(row.names = colnames(methtissuesigquomat),
                         "Tissue" = colnames(methtissuesigquomat))

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

pdf(file = "Figure 1B.pdf", width=13, height=5.5)
pheatmap(methtissuesigquomat,cluster_cols = F,cluster_rows = F,angle_col = 0,display_numbers = methtissuesigcountmat,breaks = seq(0,0.2,length.out = 9),color = brewer.pal(9,"Reds"),annotation_col = tissuemeta,annotation_colors = ann_cols,fontsize = 20,number_color = "black",cellwidth = 80,cellheight = 110)
dev.off()

png(file = "Figure 1B.png", width=13, height=5.5,units = "in",res = 600)
pheatmap(methtissuesigquomat,cluster_cols = F,cluster_rows = F,angle_col = 0,display_numbers = methtissuesigcountmat,breaks = seq(0,0.2,length.out = 9),color = brewer.pal(9,"Reds"),annotation_col = tissuemeta,annotation_colors = ann_cols,fontsize = 20,number_color = "black",cellwidth = 80,cellheight = 110)
dev.off()

rm(epigen_rrbs)
gc()

save.image("Figure1B.RData")