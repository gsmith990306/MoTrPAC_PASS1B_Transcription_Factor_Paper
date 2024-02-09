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
tfproanno <- readRDS("tfproanno.RDS")
tfanno <- readRDS("tfanno.RDS")
load("rnal2fcmat.RData")
load("enstosym.RData")

tfgastrornasig <- intersect(gastrornasig,tfproanno$Ensembl)
tfheartrnasig <- intersect(heartrnasig,tfproanno$Ensembl)
tfhippornasig <- intersect(hippornasig,tfproanno$Ensembl)
tfkidneyrnasig <- intersect(kidneyrnasig,tfproanno$Ensembl)
tfliverrnasig <- intersect(liverrnasig,tfproanno$Ensembl)
tflungrnasig <- intersect(lungrnasig,tfproanno$Ensembl)
tfbrownrnasig <- intersect(brownrnasig,tfproanno$Ensembl)
tfwhiternasig <- intersect(whiternasig,tfproanno$Ensembl)

tfrnasigmat <- rbind(gastrol2fcmat[tfgastrornasig,],
                     heartl2fcmat[tfheartrnasig,],
                     hippol2fcmat[tfhippornasig,],
                     kidneyl2fcmat[tfkidneyrnasig,],
                     liverl2fcmat[tfliverrnasig,],
                     lungl2fcmat[tflungrnasig,],
                     brownl2fcmat[tfbrownrnasig,],
                     whitel2fcmat[tfwhiternasig,])

rownames(tfrnasigmat)[1:length(tfgastrornasig)] <- paste("SKM-GN",rownames(tfrnasigmat)[1:length(tfgastrornasig)],sep = "_")
rownames(tfrnasigmat)[(length(tfgastrornasig)+1):(length(tfgastrornasig)+length(tfheartrnasig))] <- paste("HEART",rownames(tfrnasigmat)[(length(tfgastrornasig)+1):(length(tfgastrornasig)+length(tfheartrnasig))],sep = "_")
rownames(tfrnasigmat)[(length(tfgastrornasig)+length(tfheartrnasig)+1):(length(tfgastrornasig)+length(tfheartrnasig)+length(tfhippornasig))] <- paste("HIPPOC",rownames(tfrnasigmat)[(length(tfgastrornasig)+length(tfheartrnasig)+1):(length(tfgastrornasig)+length(tfheartrnasig)+length(tfhippornasig))],sep = "_")
rownames(tfrnasigmat)[(length(tfgastrornasig)+length(tfheartrnasig)+length(tfhippornasig)+1):(length(tfgastrornasig)+length(tfheartrnasig)+length(tfhippornasig)+length(tfkidneyrnasig))] <- paste("KIDNEY",rownames(tfrnasigmat)[(length(tfgastrornasig)+length(tfheartrnasig)+length(tfhippornasig)+1):(length(tfgastrornasig)+length(tfheartrnasig)+length(tfhippornasig)+length(tfkidneyrnasig))],sep = "_")
rownames(tfrnasigmat)[(length(tfgastrornasig)+length(tfheartrnasig)+length(tfhippornasig)+length(tfkidneyrnasig)+1):(length(tfgastrornasig)+length(tfheartrnasig)+length(tfhippornasig)+length(tfkidneyrnasig)+length(tfliverrnasig))] <- paste("LIVER",rownames(tfrnasigmat)[(length(tfgastrornasig)+length(tfheartrnasig)+length(tfhippornasig)+length(tfkidneyrnasig)+1):(length(tfgastrornasig)+length(tfheartrnasig)+length(tfhippornasig)+length(tfkidneyrnasig)+length(tfliverrnasig))],sep = "_")
rownames(tfrnasigmat)[(length(tfgastrornasig)+length(tfheartrnasig)+length(tfhippornasig)+length(tfkidneyrnasig)+length(tfliverrnasig)+1):(length(tfgastrornasig)+length(tfheartrnasig)+length(tfhippornasig)+length(tfkidneyrnasig)+length(tfliverrnasig)+length(tflungrnasig))] <- paste("LUNG",rownames(tfrnasigmat)[(length(tfgastrornasig)+length(tfheartrnasig)+length(tfhippornasig)+length(tfkidneyrnasig)+length(tfliverrnasig)+1):(length(tfgastrornasig)+length(tfheartrnasig)+length(tfhippornasig)+length(tfkidneyrnasig)+length(tfliverrnasig)+length(tflungrnasig))],sep = "_")
rownames(tfrnasigmat)[(length(tfgastrornasig)+length(tfheartrnasig)+length(tfhippornasig)+length(tfkidneyrnasig)+length(tfliverrnasig)+length(tflungrnasig)+1):(length(tfgastrornasig)+length(tfheartrnasig)+length(tfhippornasig)+length(tfkidneyrnasig)+length(tfliverrnasig)+length(tflungrnasig)+length(tfbrownrnasig))] <- paste("BAT",rownames(tfrnasigmat)[(length(tfgastrornasig)+length(tfheartrnasig)+length(tfhippornasig)+length(tfkidneyrnasig)+length(tfliverrnasig)+length(tflungrnasig)+1):(length(tfgastrornasig)+length(tfheartrnasig)+length(tfhippornasig)+length(tfkidneyrnasig)+length(tfliverrnasig)+length(tflungrnasig)+length(tfbrownrnasig))],sep = "_")
rownames(tfrnasigmat)[(length(tfgastrornasig)+length(tfheartrnasig)+length(tfhippornasig)+length(tfkidneyrnasig)+length(tfliverrnasig)+length(tflungrnasig)+length(tfbrownrnasig)+1):(length(tfgastrornasig)+length(tfheartrnasig)+length(tfhippornasig)+length(tfkidneyrnasig)+length(tfliverrnasig)+length(tflungrnasig)+length(tfbrownrnasig)+length(tfwhiternasig))] <- paste("WAT-SC",rownames(tfrnasigmat)[(length(tfgastrornasig)+length(tfheartrnasig)+length(tfhippornasig)+length(tfkidneyrnasig)+length(tfliverrnasig)+length(tflungrnasig)+length(tfbrownrnasig)+1):(length(tfgastrornasig)+length(tfheartrnasig)+length(tfhippornasig)+length(tfkidneyrnasig)+length(tfliverrnasig)+length(tflungrnasig)+length(tfbrownrnasig)+length(tfwhiternasig))],sep = "_")

tfsigmeta <- data.frame(row.names = colnames(tfrnasigmat),
                        "Group" = c("1w","2w","4w","8w","1w","2w","4w","8w"),
                        "Sex" = c(rep("Female",4),rep("Male",4)))

tfrnasigmeta <- data.frame(row.names = rownames(tfrnasigmat),"Tissue" = gsub("_.*","",rownames(tfrnasigmat)))

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

png("Figure 5A.png",width = 6,height = 8,units = "in",res = 600)
pheatmap(tfrnasigmat[apply(abs(tfrnasigmat),1,max) > 0.75,],cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),annotation_col = tfsigmeta[,c("Group","Sex")],annotation_row = tfrnasigmeta,annotation_colors = ann_cols,show_colnames = F,cluster_rows = F,labels_row = enstosym[gsub(".*_","",rownames(tfrnasigmat[apply(abs(tfrnasigmat),1,max) > 0.75,])),"Symbol"],display_numbers = T,number_color = "black")
dev.off()

pdf("Figure 5A.pdf",width = 6,height = 8)
pheatmap(tfrnasigmat[apply(abs(tfrnasigmat),1,max) > 0.75,],cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),annotation_col = tfsigmeta[,c("Group","Sex")],annotation_row = tfrnasigmeta,annotation_colors = ann_cols,show_colnames = F,cluster_rows = F,labels_row = enstosym[gsub(".*_","",rownames(tfrnasigmat[apply(abs(tfrnasigmat),1,max) > 0.75,])),"Symbol"],display_numbers = T,number_color = "black")
dev.off

#####
# Figure 5B
####

load("PASS1B Transcription Factor Paper Data/DEA Analysis/pass1b-06_proteomics-prot-pr_dea.RData")
load("prol2fcmat.RData")

tfgastroprosig <- prot_pr$training_dea[prot_pr$training_dea$feature_ID %in% tfproanno$Gastro.Pro.ID & prot_pr$training_dea$tissue_abbreviation %in% "SKM-GN" & prot_pr$training_dea$p_value < 0.06,"feature_ID"]
tfheartprosig <- prot_pr$training_dea[prot_pr$training_dea$feature_ID %in% tfproanno$Heart.Pro.ID & prot_pr$training_dea$tissue_abbreviation %in% "HEART" & prot_pr$training_dea$p_value < 0.05,"feature_ID"]
tfkidneyprosig <- prot_pr$training_dea[prot_pr$training_dea$feature_ID %in% tfproanno$Kidney.Pro.ID & prot_pr$training_dea$tissue_abbreviation %in% "KIDNEY" & prot_pr$training_dea$p_value < 0.05,"feature_ID"]
tfliverprosig <- prot_pr$training_dea[prot_pr$training_dea$feature_ID %in% tfproanno$Liver.Pro.ID & prot_pr$training_dea$tissue_abbreviation %in% "LIVER" & prot_pr$training_dea$p_value < 0.05,"feature_ID"]
tflungprosig <- prot_pr$training_dea[prot_pr$training_dea$feature_ID %in% tfproanno$Lung.Pro.ID & prot_pr$training_dea$tissue_abbreviation %in% "LUNG" & prot_pr$training_dea$p_value < 0.05,"feature_ID"]
tfwhiteprosig <- prot_pr$training_dea[prot_pr$training_dea$feature_ID %in% tfproanno$WhiteAd.Pro.ID & prot_pr$training_dea$tissue_abbreviation %in% "WAT-SC" & prot_pr$training_dea$p_value < 0.05,"feature_ID"]

# We want a heatmap where each row is a sig protein in its tissue of interest

tfproannogastroselect <- tfproanno[tfproanno$Gastro.Pro.ID %in% tfgastroprosig,c("Gene.Name","Gastro.Pro.ID")]
tfproannogastroselect <- tfproannogastroselect[!duplicated(tfproannogastroselect$Gastro.Pro.ID),]

tfproannoheartselect <- tfproanno[tfproanno$Heart.Pro.ID %in% tfheartprosig,c("Gene.Name","Heart.Pro.ID")]
tfproannoheartselect <- tfproannoheartselect[!duplicated(tfproannoheartselect$Heart.Pro.ID),]

tfproannokidneyselect <- tfproanno[tfproanno$Kidney.Pro.ID %in% tfkidneyprosig,c("Gene.Name","Kidney.Pro.ID")]
tfproannokidneyselect <- tfproannokidneyselect[!duplicated(tfproannokidneyselect$Kidney.Pro.ID),]

tfproannoliverselect <- tfproanno[tfproanno$Liver.Pro.ID %in% tfliverprosig,c("Gene.Name","Liver.Pro.ID")]
tfproannoliverselect <- tfproannoliverselect[!duplicated(tfproannoliverselect$Liver.Pro.ID),]

tfproannolungselect <- tfproanno[tfproanno$Lung.Pro.ID %in% tflungprosig,c("Gene.Name","Lung.Pro.ID")]
tfproannolungselect <- tfproannolungselect[!duplicated(tfproannolungselect$Lung.Pro.ID),]

tfproannowhiteselect <- tfproanno[tfproanno$WhiteAd.Pro.ID %in% tfwhiteprosig,c("Gene.Name","WhiteAd.Pro.ID")]
tfproannowhiteselect <- tfproannowhiteselect[!duplicated(tfproannowhiteselect$WhiteAd.Pro.ID),]

tfgastroprosig <- tfproannogastroselect$Gastro.Pro.ID
tfheartprosig <- tfproannoheartselect$Heart.Pro.ID
tfkidneyprosig <- tfproannokidneyselect$Kidney.Pro.ID
tfliverprosig <- tfproannoliverselect$Liver.Pro.ID
tflungprosig <- tfproannolungselect$Lung.Pro.ID
tfwhiteprosig <- tfproannowhiteselect$WhiteAd.Pro.ID

tfprosigmat <- rbind(gastroprol2fc[tfgastroprosig,],
                     heartprol2fc[tfheartprosig,],
                     kidneyprol2fc[tfkidneyprosig,],
                     liverprol2fc[tfliverprosig,],
                     lungprol2fc[tflungprosig,],
                     whiteprol2fc[tfwhiteprosig,])
rownames(tfprosigmat)[1:length(tfgastroprosig)] <- paste("SKM-GN",rownames(tfprosigmat)[1:length(tfgastroprosig)],sep = "_")
rownames(tfprosigmat)[(length(tfgastroprosig)+1):(length(tfgastroprosig)+length(tfheartprosig))] <- paste("HEART",rownames(tfprosigmat)[(length(tfgastroprosig)+1):(length(tfgastroprosig)+length(tfheartprosig))],sep = "_")
rownames(tfprosigmat)[(length(tfgastroprosig)+length(tfheartprosig)+1):(length(tfgastroprosig)+length(tfheartprosig)+length(tfkidneyprosig))] <- paste("KIDNEY",rownames(tfprosigmat)[(length(tfgastroprosig)+length(tfheartprosig)+1):(length(tfgastroprosig)+length(tfheartprosig)+length(tfkidneyprosig))],sep = "_")
rownames(tfprosigmat)[(length(tfgastroprosig)+length(tfheartprosig)+length(tfkidneyprosig)+1):(length(tfgastroprosig)+length(tfheartprosig)+length(tfkidneyprosig)+length(tfliverprosig))] <- paste("LIVER",rownames(tfprosigmat)[(length(tfgastroprosig)+length(tfheartprosig)+length(tfkidneyprosig)+1):(length(tfgastroprosig)+length(tfheartprosig)+length(tfkidneyprosig)+length(tfliverprosig))],sep = "_")
rownames(tfprosigmat)[(length(tfgastroprosig)+length(tfheartprosig)+length(tfkidneyprosig)+length(tfliverprosig)+1):(length(tfgastroprosig)+length(tfheartprosig)+length(tfkidneyprosig)+length(tfliverprosig)+length(tflungprosig))] <- paste("LUNG",rownames(tfprosigmat)[(length(tfgastroprosig)+length(tfheartprosig)+length(tfkidneyprosig)+length(tfliverprosig)+1):(length(tfgastroprosig)+length(tfheartprosig)+length(tfkidneyprosig)+length(tfliverprosig)+length(tflungprosig))],sep = "_")
rownames(tfprosigmat)[(length(tfgastroprosig)+length(tfheartprosig)+length(tfkidneyprosig)+length(tfliverprosig)+length(tflungprosig)+1):(length(tfgastroprosig)+length(tfheartprosig)+length(tfkidneyprosig)+length(tfliverprosig)+length(tflungprosig)+length(tfwhiteprosig))] <- paste("WAT-SC",rownames(tfprosigmat)[(length(tfgastroprosig)+length(tfheartprosig)+length(tfkidneyprosig)+length(tfliverprosig)+length(tflungprosig)+1):(length(tfgastroprosig)+length(tfheartprosig)+length(tfkidneyprosig)+length(tfliverprosig)+length(tflungprosig)+length(tfwhiteprosig))],sep = "_")

tfsigmeta <- data.frame(row.names = colnames(tfprosigmat),
                        "Group" = c("1w","2w","4w","8w","1w","2w","4w","8w"),
                        "Sex" = c(rep("Female",4),rep("Male",4)))

tfprosigmeta <- data.frame(row.names = rownames(tfprosigmat),
                           "Tissue" = gsub("_.*","",rownames(tfprosigmat)))

png("Figure 5B.png",width = 6,height = 8,units = "in",res = 600)
pheatmap(tfprosigmat[apply(abs(tfprosigmat),1,max) > 0.25,],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),annotation_col = tfsigmeta[,c("Group","Sex")],annotation_row = tfprosigmeta,annotation_colors = ann_cols,show_colnames = F,cluster_rows = F,labels_row = toupper(c(tfproannogastroselect$Gene.Name,tfproannoheartselect$Gene.Name,tfproannokidneyselect$Gene.Name,tfproannoliverselect$Gene.Name,tfproannolungselect$Gene.Name,tfproannowhiteselect$Gene.Name))[apply(abs(tfprosigmat),1,max) > 0.25],display_numbers = T,number_color = "black")
dev.off()

pdf("Figure 5B.pdf",width = 6,height = 8)
pheatmap(tfprosigmat[apply(abs(tfprosigmat),1,max) > 0.25,],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),annotation_col = tfsigmeta[,c("Group","Sex")],annotation_row = tfprosigmeta,annotation_colors = ann_cols,show_colnames = F,cluster_rows = F,labels_row = toupper(c(tfproannogastroselect$Gene.Name,tfproannoheartselect$Gene.Name,tfproannokidneyselect$Gene.Name,tfproannoliverselect$Gene.Name,tfproannolungselect$Gene.Name,tfproannowhiteselect$Gene.Name))[apply(abs(tfprosigmat),1,max) > 0.25],display_numbers = T,number_color = "black")
dev.off()

#####
# Figure 5C
####

load("PASS1B Transcription Factor Paper Data/DEA Analysis/pass1b-06_proteomics-prot-ph_dea.RData")
load("prophl2fcmat.RData")

prot_ph$training_dea$pr_feature_ID <- sapply(strsplit(prot_ph$training_dea$feature_ID, "_"), function(x) {
  paste(x[1],x[2],sep = "_")
})
prot_ph$timewise_dea$pr_feature_ID <- sapply(strsplit(prot_ph$timewise_dea$feature_ID, "_"), function(x) {
  paste(x[1],x[2],sep = "_")
})

tfgastroprophsig <- prot_ph$training_dea[prot_ph$training_dea$pr_feature_ID %in% tfproanno$Gastro.Pro.ID & prot_ph$training_dea$tissue_abbreviation %in% "SKM-GN" & prot_ph$training_dea$p_value < 0.01,"feature_ID"]
tfheartprophsig <- prot_ph$training_dea[prot_ph$training_dea$pr_feature_ID %in% tfproanno$Heart.Pro.ID & prot_ph$training_dea$tissue_abbreviation %in% "HEART" & prot_ph$training_dea$p_value < 0.01,"feature_ID"]
tfkidneyprophsig <- prot_ph$training_dea[prot_ph$training_dea$pr_feature_ID %in% tfproanno$Kidney.Pro.ID & prot_ph$training_dea$tissue_abbreviation %in% "KIDNEY" & prot_ph$training_dea$p_value < 0.01,"feature_ID"]
tfliverprophsig <- prot_ph$training_dea[prot_ph$training_dea$pr_feature_ID %in% tfproanno$Liver.Pro.ID & prot_ph$training_dea$tissue_abbreviation %in% "LIVER" & prot_ph$training_dea$p_value < 0.01,"feature_ID"]
tflungprophsig <- prot_ph$training_dea[prot_ph$training_dea$pr_feature_ID %in% tfproanno$Lung.Pro.ID & prot_ph$training_dea$tissue_abbreviation %in% "LUNG" & prot_ph$training_dea$p_value < 0.01,"feature_ID"]
tfwhiteprophsig <- prot_ph$training_dea[prot_ph$training_dea$pr_feature_ID %in% tfproanno$WhiteAd.Pro.ID & prot_ph$training_dea$tissue_abbreviation %in% "WAT-SC" & prot_ph$training_dea$p_value < 0.01,"feature_ID"]


sapply(strsplit(tfgastroprophsig, "_"), function(x) {
  paste(x[1],x[2],sep = "_")
})


tfprophannogastroselect <- data.frame(row.names = tfgastroprophsig,"Phospho" = tfgastroprophsig,"Protein_ID" = sapply(strsplit(tfgastroprophsig, "_"), function(x) {
  paste(x[1],x[2],sep = "_")
}))
tfprophannogastroselect$Symbol = ""
for(i in 1:length(tfgastroprophsig)){
  ourpro <- tfprophannogastroselect[i,"Protein_ID"]
  tfprophannogastroselect[i,"Symbol"] <- toupper(tfproanno[tfproanno$Gastro.Pro.ID %in% ourpro,"Gene.Name"][1])
}

tfprophannoheartselect <- data.frame(row.names = tfheartprophsig,"Phospho" = tfheartprophsig,"Protein_ID" = sapply(strsplit(tfheartprophsig, "_"), function(x) {
  paste(x[1],x[2],sep = "_")
}))
tfprophannoheartselect$Symbol = ""
for(i in 1:length(tfheartprophsig)){
  ourpro <- tfprophannoheartselect[i,"Protein_ID"]
  tfprophannoheartselect[i,"Symbol"] <- toupper(tfproanno[tfproanno$Heart.Pro.ID %in% ourpro,"Gene.Name"][1])
}

tfprophannokidneyselect <- data.frame(row.names = tfkidneyprophsig,"Phospho" = tfkidneyprophsig,"Protein_ID" = sapply(strsplit(tfkidneyprophsig, "_"), function(x) {
  paste(x[1],x[2],sep = "_")
}))
tfprophannokidneyselect$Symbol = ""
for(i in 1:length(tfkidneyprophsig)){
  ourpro <- tfprophannokidneyselect[i,"Protein_ID"]
  tfprophannokidneyselect[i,"Symbol"] <- toupper(tfproanno[tfproanno$Kidney.Pro.ID %in% ourpro,"Gene.Name"][1])
}

tfprophannoliverselect <- data.frame(row.names = tfliverprophsig,"Phospho" = tfliverprophsig,"Protein_ID" = sapply(strsplit(tfliverprophsig, "_"), function(x) {
  paste(x[1],x[2],sep = "_")
}))
tfprophannoliverselect$Symbol = ""
for(i in 1:length(tfliverprophsig)){
  ourpro <- tfprophannoliverselect[i,"Protein_ID"]
  tfprophannoliverselect[i,"Symbol"] <- toupper(tfproanno[tfproanno$Liver.Pro.ID %in% ourpro,"Gene.Name"][1])
}

tfprophannolungselect <- data.frame(row.names = tflungprophsig,"Phospho" = tflungprophsig,"Protein_ID" = sapply(strsplit(tflungprophsig, "_"), function(x) {
  paste(x[1],x[2],sep = "_")
}))
tfprophannolungselect$Symbol = ""
for(i in 1:length(tflungprophsig)){
  ourpro <- tfprophannolungselect[i,"Protein_ID"]
  tfprophannolungselect[i,"Symbol"] <- toupper(tfproanno[tfproanno$Lung.Pro.ID %in% ourpro,"Gene.Name"][1])
}

tfprophannowhiteselect <- data.frame(row.names = tfwhiteprophsig,"Phospho" = tfwhiteprophsig,"Protein_ID" = sapply(strsplit(tfwhiteprophsig, "_"), function(x) {
  paste(x[1],x[2],sep = "_")
}))
tfprophannowhiteselect$Symbol = ""
for(i in 1:length(tfwhiteprophsig)){
  ourpro <- tfprophannowhiteselect[i,"Protein_ID"]
  tfprophannowhiteselect[i,"Symbol"] <- toupper(tfproanno[tfproanno$WhiteAd.Pro.ID %in% ourpro,"Gene.Name"][1])
}

tfprophsigmat <- rbind(gastroprophl2fc[tfgastroprophsig,],
                       heartprophl2fc[tfheartprophsig,],
                       kidneyprophl2fc[tfkidneyprophsig,],
                       liverprophl2fc[tfliverprophsig,],
                       lungprophl2fc[tflungprophsig,],
                       whiteprophl2fc[tfwhiteprophsig,])
rownames(tfprophsigmat)[1:length(tfgastroprophsig)] <- paste("SKM-GN",rownames(tfprophsigmat)[1:length(tfgastroprophsig)],sep = "_")
rownames(tfprophsigmat)[(length(tfgastroprophsig)+1):(length(tfgastroprophsig)+length(tfheartprophsig))] <- paste("HEART",rownames(tfprophsigmat)[(length(tfgastroprophsig)+1):(length(tfgastroprophsig)+length(tfheartprophsig))],sep = "_")
rownames(tfprophsigmat)[(length(tfgastroprophsig)+length(tfheartprophsig)+1):(length(tfgastroprophsig)+length(tfheartprophsig)+length(tfkidneyprophsig))] <- paste("KIDNEY",rownames(tfprophsigmat)[(length(tfgastroprophsig)+length(tfheartprophsig)+1):(length(tfgastroprophsig)+length(tfheartprophsig)+length(tfkidneyprophsig))],sep = "_")
rownames(tfprophsigmat)[(length(tfgastroprophsig)+length(tfheartprophsig)+length(tfkidneyprophsig)+1):(length(tfgastroprophsig)+length(tfheartprophsig)+length(tfkidneyprophsig)+length(tfliverprophsig))] <- paste("LIVER",rownames(tfprophsigmat)[(length(tfgastroprophsig)+length(tfheartprophsig)+length(tfkidneyprophsig)+1):(length(tfgastroprophsig)+length(tfheartprophsig)+length(tfkidneyprophsig)+length(tfliverprophsig))],sep = "_")
rownames(tfprophsigmat)[(length(tfgastroprophsig)+length(tfheartprophsig)+length(tfkidneyprophsig)+length(tfliverprophsig)+1):(length(tfgastroprophsig)+length(tfheartprophsig)+length(tfkidneyprophsig)+length(tfliverprophsig)+length(tflungprophsig))] <- paste("LUNG",rownames(tfprophsigmat)[(length(tfgastroprophsig)+length(tfheartprophsig)+length(tfkidneyprophsig)+length(tfliverprophsig)+1):(length(tfgastroprophsig)+length(tfheartprophsig)+length(tfkidneyprophsig)+length(tfliverprophsig)+length(tflungprophsig))],sep = "_")
rownames(tfprophsigmat)[(length(tfgastroprophsig)+length(tfheartprophsig)+length(tfkidneyprophsig)+length(tfliverprophsig)+length(tflungprophsig)+1):(length(tfgastroprophsig)+length(tfheartprophsig)+length(tfkidneyprophsig)+length(tfliverprophsig)+length(tflungprophsig)+length(tfwhiteprophsig))] <- paste("WAT-SC",rownames(tfprophsigmat)[(length(tfgastroprophsig)+length(tfheartprophsig)+length(tfkidneyprophsig)+length(tfliverprophsig)+length(tflungprophsig)+1):(length(tfgastroprophsig)+length(tfheartprophsig)+length(tfkidneyprophsig)+length(tfliverprophsig)+length(tflungprophsig)+length(tfwhiteprophsig))],sep = "_")

tfprophsigmeta <- data.frame(row.names = rownames(tfprophsigmat),
                             "Tissue" = gsub("_.*","",rownames(tfprophsigmat)))



png("Figure 5C.png",width = 6,height = 8,units = "in",res = 600)
pheatmap(tfprophsigmat[apply(abs(tfprophsigmat),1,max) > 0.4,],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),annotation_col = tfsigmeta[,c("Group","Sex")],annotation_row = tfprophsigmeta,annotation_colors = ann_cols,show_colnames = F,cluster_rows = F,
         labels_row = c(paste(toupper(tfprophannogastroselect$Symbol),gsub(".*_","",tfprophannogastroselect$Phospho),sep = "_"),
                        paste(toupper(tfprophannoheartselect$Symbol),gsub(".*_","",tfprophannoheartselect$Phospho),sep = "_"),
                        paste(toupper(tfprophannokidneyselect$Symbol),gsub(".*_","",tfprophannokidneyselect$Phospho),sep = "_"),
                        paste(toupper(tfprophannoliverselect$Symbol),gsub(".*_","",tfprophannoliverselect$Phospho),sep = "_"),
                        paste(toupper(tfprophannolungselect$Symbol),gsub(".*_","",tfprophannolungselect$Phospho),sep = "_"),
                        paste(toupper(tfprophannowhiteselect$Symbol),gsub(".*_","",tfprophannowhiteselect$Phospho),sep = "_"))[apply(abs(tfprophsigmat),1,max) > 0.4],
         display_numbers = T,number_color = "black")
dev.off()

pdf("Figure 5C.pdf",width = 6,height = 8)
pheatmap(tfprophsigmat[apply(abs(tfprophsigmat),1,max) > 0.4,],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),annotation_col = tfsigmeta[,c("Group","Sex")],annotation_row = tfprophsigmeta,annotation_colors = ann_cols,show_colnames = F,cluster_rows = F,
         labels_row = c(paste(toupper(tfprophannogastroselect$Symbol),gsub(".*_","",tfprophannogastroselect$Phospho),sep = "_"),
                        paste(toupper(tfprophannoheartselect$Symbol),gsub(".*_","",tfprophannoheartselect$Phospho),sep = "_"),
                        paste(toupper(tfprophannokidneyselect$Symbol),gsub(".*_","",tfprophannokidneyselect$Phospho),sep = "_"),
                        paste(toupper(tfprophannoliverselect$Symbol),gsub(".*_","",tfprophannoliverselect$Phospho),sep = "_"),
                        paste(toupper(tfprophannolungselect$Symbol),gsub(".*_","",tfprophannolungselect$Phospho),sep = "_"),
                        paste(toupper(tfprophannowhiteselect$Symbol),gsub(".*_","",tfprophannowhiteselect$Phospho),sep = "_"))[apply(abs(tfprophsigmat),1,max) > 0.4],
         display_numbers = T,number_color = "black")
dev.off()

save.image("Figure5A_5B_5C.RData")

#####
# Supplemental Figure S14
####

# We need to investigate what are the compelling relationships between changing TFs at the
# protein and phosphoprotein level and their enrichment among DEGs

oursigtfprolist <- unique(c(rownames(tfproanno[tfproanno$Gastro.Pro.ID %in% gsub(".*XP","XP",gsub(".*NP","NP",rownames(tfprosigmat[apply(abs(tfprosigmat),1,max) > 0.25,]))),]),
                            rownames(tfproanno[tfproanno$Heart.Pro.ID %in% gsub(".*XP","XP",gsub(".*NP","NP",rownames(tfprosigmat[apply(abs(tfprosigmat),1,max) > 0.25,]))),]),
                            rownames(tfproanno[tfproanno$Kidney.Pro.ID %in% gsub(".*XP","XP",gsub(".*NP","NP",rownames(tfprosigmat[apply(abs(tfprosigmat),1,max) > 0.25,]))),]),
                            rownames(tfproanno[tfproanno$Liver.Pro.ID %in% gsub(".*XP","XP",gsub(".*NP","NP",rownames(tfprosigmat[apply(abs(tfprosigmat),1,max) > 0.25,]))),]),
                            rownames(tfproanno[tfproanno$Lung.Pro.ID %in% gsub(".*XP","XP",gsub(".*NP","NP",rownames(tfprosigmat[apply(abs(tfprosigmat),1,max) > 0.25,]))),]),
                            rownames(tfproanno[tfproanno$WhiteAd.Pro.ID %in% gsub(".*XP","XP",gsub(".*NP","NP",rownames(tfprosigmat[apply(abs(tfprosigmat),1,max) > 0.25,]))),])))

oursigtfprophlist <- unique(c(rownames(tfproanno[tfproanno$Gastro.Pro.ID %in% gsub(".*XP","XP",gsub(".*NP","NP",sub("_[^_]+$","",gsub(".*XP","XP",gsub(".*NP","NP",rownames(tfprophsigmat[apply(abs(tfprophsigmat),1,max) > 0.4,])))))),]),
                              rownames(tfproanno[tfproanno$Heart.Pro.ID %in% gsub(".*XP","XP",gsub(".*NP","NP",sub("_[^_]+$","",gsub(".*XP","XP",gsub(".*NP","NP",rownames(tfprophsigmat[apply(abs(tfprophsigmat),1,max) > 0.4,])))))),]),
                              rownames(tfproanno[tfproanno$Kidney.Pro.ID %in% gsub(".*XP","XP",gsub(".*NP","NP",sub("_[^_]+$","",gsub(".*XP","XP",gsub(".*NP","NP",rownames(tfprophsigmat[apply(abs(tfprophsigmat),1,max) > 0.4,])))))),]),
                              rownames(tfproanno[tfproanno$Liver.Pro.ID %in% gsub(".*XP","XP",gsub(".*NP","NP",sub("_[^_]+$","",gsub(".*XP","XP",gsub(".*NP","NP",rownames(tfprophsigmat[apply(abs(tfprophsigmat),1,max) > 0.4,])))))),]),
                              rownames(tfproanno[tfproanno$Lung.Pro.ID %in% gsub(".*XP","XP",gsub(".*NP","NP",sub("_[^_]+$","",gsub(".*XP","XP",gsub(".*NP","NP",rownames(tfprophsigmat[apply(abs(tfprophsigmat),1,max) > 0.4,])))))),]),
                              rownames(tfproanno[tfproanno$WhiteAd.Pro.ID %in% gsub(".*XP","XP",gsub(".*NP","NP",sub("_[^_]+$","",gsub(".*XP","XP",gsub(".*NP","NP",rownames(tfprophsigmat[apply(abs(tfprophsigmat),1,max) > 0.4,])))))),])))

gastrotfproandprophlist <- unique(c(rownames(tfproanno[tfproanno$Gastro.Pro.ID %in% tfgastroprosig,]),
                                    rownames(tfproanno[tfproanno$Gastro.Pro.ID %in% sub("_[^_]+$","",tfgastroprophsig),])))
gastrotfproandprophmeta <- data.frame(row.names = gastrotfproandprophlist,"ID" = gastrotfproandprophlist)
gastrotfproandprophmeta$Pro.Sig <- 0
gastrotfproandprophmeta$Proph.Sig <- 0
gastrotfproandprophmeta[gastrotfproandprophmeta$ID %in% rownames(tfproanno[tfproanno$Gastro.Pro.ID %in% tfgastroprosig,]),"Pro.Sig"] <- 1
gastrotfproandprophmeta[gastrotfproandprophmeta$ID %in% rownames(tfproanno[tfproanno$Gastro.Pro.ID %in% sub("_[^_]+$","",tfgastroprophsig),]),"Proph.Sig"] <- 1

hearttfproandprophlist <- unique(c(rownames(tfproanno[tfproanno$Heart.Pro.ID %in% tfheartprosig,]),
                                   rownames(tfproanno[tfproanno$Heart.Pro.ID %in% sub("_[^_]+$","",tfheartprophsig),])))
hearttfproandprophmeta <- data.frame(row.names = hearttfproandprophlist,"ID" = hearttfproandprophlist)
hearttfproandprophmeta$Pro.Sig <- 0
hearttfproandprophmeta$Proph.Sig <- 0
hearttfproandprophmeta[hearttfproandprophmeta$ID %in% rownames(tfproanno[tfproanno$Heart.Pro.ID %in% tfheartprosig,]),"Pro.Sig"] <- 1
hearttfproandprophmeta[hearttfproandprophmeta$ID %in% rownames(tfproanno[tfproanno$Heart.Pro.ID %in% sub("_[^_]+$","",tfheartprophsig),]),"Proph.Sig"] <- 1

kidneytfproandprophlist <- unique(c(rownames(tfproanno[tfproanno$Kidney.Pro.ID %in% tfkidneyprosig,]),
                                    rownames(tfproanno[tfproanno$Kidney.Pro.ID %in% sub("_[^_]+$","",tfkidneyprophsig),])))
kidneytfproandprophmeta <- data.frame(row.names = kidneytfproandprophlist,"ID" = kidneytfproandprophlist)
kidneytfproandprophmeta$Pro.Sig <- 0
kidneytfproandprophmeta$Proph.Sig <- 0
kidneytfproandprophmeta[kidneytfproandprophmeta$ID %in% rownames(tfproanno[tfproanno$Kidney.Pro.ID %in% tfkidneyprosig,]),"Pro.Sig"] <- 1
kidneytfproandprophmeta[kidneytfproandprophmeta$ID %in% rownames(tfproanno[tfproanno$Kidney.Pro.ID %in% sub("_[^_]+$","",tfkidneyprophsig),]),"Proph.Sig"] <- 1

livertfproandprophlist <- unique(c(rownames(tfproanno[tfproanno$Liver.Pro.ID %in% tfliverprosig,]),
                                   rownames(tfproanno[tfproanno$Liver.Pro.ID %in% sub("_[^_]+$","",tfliverprophsig),])))
livertfproandprophmeta <- data.frame(row.names = livertfproandprophlist,"ID" = livertfproandprophlist)
livertfproandprophmeta$Pro.Sig <- 0
livertfproandprophmeta$Proph.Sig <- 0
livertfproandprophmeta[livertfproandprophmeta$ID %in% rownames(tfproanno[tfproanno$Liver.Pro.ID %in% tfliverprosig,]),"Pro.Sig"] <- 1
livertfproandprophmeta[livertfproandprophmeta$ID %in% rownames(tfproanno[tfproanno$Liver.Pro.ID %in% sub("_[^_]+$","",tfliverprophsig),]),"Proph.Sig"] <- 1

lungtfproandprophlist <- unique(c(rownames(tfproanno[tfproanno$Lung.Pro.ID %in% tflungprosig,]),
                                  rownames(tfproanno[tfproanno$Lung.Pro.ID %in% sub("_[^_]+$","",tflungprophsig),])))
lungtfproandprophmeta <- data.frame(row.names = lungtfproandprophlist,"ID" = lungtfproandprophlist)
lungtfproandprophmeta$Pro.Sig <- 0
lungtfproandprophmeta$Proph.Sig <- 0
lungtfproandprophmeta[lungtfproandprophmeta$ID %in% rownames(tfproanno[tfproanno$Lung.Pro.ID %in% tflungprosig,]),"Pro.Sig"] <- 1
lungtfproandprophmeta[lungtfproandprophmeta$ID %in% rownames(tfproanno[tfproanno$Lung.Pro.ID %in% sub("_[^_]+$","",tflungprophsig),]),"Proph.Sig"] <- 1

whitetfproandprophlist <- unique(c(rownames(tfproanno[tfproanno$WhiteAd.Pro.ID %in% tfwhiteprosig,]),
                                   rownames(tfproanno[tfproanno$WhiteAd.Pro.ID %in% sub("_[^_]+$","",tfwhiteprophsig),])))
whitetfproandprophmeta <- data.frame(row.names = whitetfproandprophlist,"ID" = whitetfproandprophlist)
whitetfproandprophmeta$Pro.Sig <- 0
whitetfproandprophmeta$Proph.Sig <- 0
whitetfproandprophmeta[whitetfproandprophmeta$ID %in% rownames(tfproanno[tfproanno$WhiteAd.Pro.ID %in% tfwhiteprosig,]),"Pro.Sig"] <- 1
whitetfproandprophmeta[whitetfproandprophmeta$ID %in% rownames(tfproanno[tfproanno$WhiteAd.Pro.ID %in% sub("_[^_]+$","",tfwhiteprophsig),]),"Proph.Sig"] <- 1


allpeakmotifs <- readRDS("allpeakmotifs.RDS")

load("activepeakfiles.RData")
peakanno <- readRDS("peakanno.RDS")

gastro50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% gastroactivepeaks,]
heart50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% heartactivepeaks,]
hippo50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% hippoactivepeaks,]
kidney50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% kidneyactivepeaks,]
liver50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% liveractivepeaks,]
lung50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% lungactivepeaks,]
brown50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% brownactivepeaks,]
white50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% whiteactivepeaks,]


# SKM-GN 

gastrotftargetpromproxgenesigpct <- matrix(0L,nrow = length(gastrotfproandprophlist),ncol = 1)
rownames(gastrotftargetpromproxgenesigpct) <- gastrotfproandprophlist
colnames(gastrotftargetpromproxgenesigpct) <- "Percent.Significant"

gastrotftargetintrongenesigpct <- matrix(0L,nrow = length(gastrotfproandprophlist),ncol = 1)
rownames(gastrotftargetintrongenesigpct) <- gastrotfproandprophlist
colnames(gastrotftargetintrongenesigpct) <- "Percent.Significant"

gastrotftargetexongenesigpct <- matrix(0L,nrow = length(gastrotfproandprophlist),ncol = 1)
rownames(gastrotftargetexongenesigpct) <- gastrotfproandprophlist
colnames(gastrotftargetexongenesigpct) <- "Percent.Significant"

gastrotftargetintergenicgenesigpct <- matrix(0L,nrow = length(gastrotfproandprophlist),ncol = 1)
rownames(gastrotftargetintergenicgenesigpct) <- gastrotfproandprophlist
colnames(gastrotftargetintergenicgenesigpct) <- "Percent.Significant"



gastrotftargetpromproxgenesigcount <- matrix(0L,nrow = length(gastrotfproandprophlist),ncol = 1)
rownames(gastrotftargetpromproxgenesigcount) <- gastrotfproandprophlist
colnames(gastrotftargetpromproxgenesigcount) <- "Count.Significant"

gastrotftargetintrongenesigcount <- matrix(0L,nrow = length(gastrotfproandprophlist),ncol = 1)
rownames(gastrotftargetintrongenesigcount) <- gastrotfproandprophlist
colnames(gastrotftargetintrongenesigcount) <- "Count.Significant"

gastrotftargetexongenesigcount <- matrix(0L,nrow = length(gastrotfproandprophlist),ncol = 1)
rownames(gastrotftargetexongenesigcount) <- gastrotfproandprophlist
colnames(gastrotftargetexongenesigcount) <- "Count.Significant"

gastrotftargetintergenicgenesigcount <- matrix(0L,nrow = length(gastrotfproandprophlist),ncol = 1)
rownames(gastrotftargetintergenicgenesigcount) <- gastrotfproandprophlist
colnames(gastrotftargetintergenicgenesigcount) <- "Count.Significant"


for(i in 1:length(gastrotfproandprophlist)){
  ourtf <- gastrotfproandprophlist[i]
  ourtftargetpeaks <- gastro50peakmotifs[gastro50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
  ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
  ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
  ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
  ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
  ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])
  gastrotftargetpromproxgenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetpromproxgenes,gastrornasig))
  gastrotftargetintrongenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetintrongenes,gastrornasig))
  gastrotftargetexongenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetexongenes,gastrornasig))
  gastrotftargetintergenicgenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetintergenicgenes,gastrornasig))
  gastrotftargetpromproxgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetpromproxgenes,gastrornasig))/length(ourtftargetpromproxgenes)
  gastrotftargetintrongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintrongenes,gastrornasig))/length(ourtftargetintrongenes)
  gastrotftargetexongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetexongenes,gastrornasig))/length(ourtftargetexongenes)
  gastrotftargetintergenicgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintergenicgenes,gastrornasig))/length(ourtftargetintergenicgenes)
}


# HEART 

hearttftargetpromproxgenesigpct <- matrix(0L,nrow = length(hearttfproandprophlist),ncol = 1)
rownames(hearttftargetpromproxgenesigpct) <- hearttfproandprophlist
colnames(hearttftargetpromproxgenesigpct) <- "Percent.Significant"

hearttftargetintrongenesigpct <- matrix(0L,nrow = length(hearttfproandprophlist),ncol = 1)
rownames(hearttftargetintrongenesigpct) <- hearttfproandprophlist
colnames(hearttftargetintrongenesigpct) <- "Percent.Significant"

hearttftargetexongenesigpct <- matrix(0L,nrow = length(hearttfproandprophlist),ncol = 1)
rownames(hearttftargetexongenesigpct) <- hearttfproandprophlist
colnames(hearttftargetexongenesigpct) <- "Percent.Significant"

hearttftargetintergenicgenesigpct <- matrix(0L,nrow = length(hearttfproandprophlist),ncol = 1)
rownames(hearttftargetintergenicgenesigpct) <- hearttfproandprophlist
colnames(hearttftargetintergenicgenesigpct) <- "Percent.Significant"



hearttftargetpromproxgenesigcount <- matrix(0L,nrow = length(hearttfproandprophlist),ncol = 1)
rownames(hearttftargetpromproxgenesigcount) <- hearttfproandprophlist
colnames(hearttftargetpromproxgenesigcount) <- "Count.Significant"

hearttftargetintrongenesigcount <- matrix(0L,nrow = length(hearttfproandprophlist),ncol = 1)
rownames(hearttftargetintrongenesigcount) <- hearttfproandprophlist
colnames(hearttftargetintrongenesigcount) <- "Count.Significant"

hearttftargetexongenesigcount <- matrix(0L,nrow = length(hearttfproandprophlist),ncol = 1)
rownames(hearttftargetexongenesigcount) <- hearttfproandprophlist
colnames(hearttftargetexongenesigcount) <- "Count.Significant"

hearttftargetintergenicgenesigcount <- matrix(0L,nrow = length(hearttfproandprophlist),ncol = 1)
rownames(hearttftargetintergenicgenesigcount) <- hearttfproandprophlist
colnames(hearttftargetintergenicgenesigcount) <- "Count.Significant"


for(i in 1:length(hearttfproandprophlist)){
  ourtf <- hearttfproandprophlist[i]
  ourtftargetpeaks <- heart50peakmotifs[heart50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
  ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
  ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
  ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
  ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
  ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])
  hearttftargetpromproxgenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetpromproxgenes,heartrnasig))
  hearttftargetintrongenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetintrongenes,heartrnasig))
  hearttftargetexongenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetexongenes,heartrnasig))
  hearttftargetintergenicgenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetintergenicgenes,heartrnasig))
  hearttftargetpromproxgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetpromproxgenes,heartrnasig))/length(ourtftargetpromproxgenes)
  hearttftargetintrongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintrongenes,heartrnasig))/length(ourtftargetintrongenes)
  hearttftargetexongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetexongenes,heartrnasig))/length(ourtftargetexongenes)
  hearttftargetintergenicgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintergenicgenes,heartrnasig))/length(ourtftargetintergenicgenes)
}


# KIDNEY

kidneytftargetpromproxgenesigpct <- matrix(0L,nrow = length(kidneytfproandprophlist),ncol = 1)
rownames(kidneytftargetpromproxgenesigpct) <- kidneytfproandprophlist
colnames(kidneytftargetpromproxgenesigpct) <- "Percent.Significant"

kidneytftargetintrongenesigpct <- matrix(0L,nrow = length(kidneytfproandprophlist),ncol = 1)
rownames(kidneytftargetintrongenesigpct) <- kidneytfproandprophlist
colnames(kidneytftargetintrongenesigpct) <- "Percent.Significant"

kidneytftargetexongenesigpct <- matrix(0L,nrow = length(kidneytfproandprophlist),ncol = 1)
rownames(kidneytftargetexongenesigpct) <- kidneytfproandprophlist
colnames(kidneytftargetexongenesigpct) <- "Percent.Significant"

kidneytftargetintergenicgenesigpct <- matrix(0L,nrow = length(kidneytfproandprophlist),ncol = 1)
rownames(kidneytftargetintergenicgenesigpct) <- kidneytfproandprophlist
colnames(kidneytftargetintergenicgenesigpct) <- "Percent.Significant"



kidneytftargetpromproxgenesigcount <- matrix(0L,nrow = length(kidneytfproandprophlist),ncol = 1)
rownames(kidneytftargetpromproxgenesigcount) <- kidneytfproandprophlist
colnames(kidneytftargetpromproxgenesigcount) <- "Count.Significant"

kidneytftargetintrongenesigcount <- matrix(0L,nrow = length(kidneytfproandprophlist),ncol = 1)
rownames(kidneytftargetintrongenesigcount) <- kidneytfproandprophlist
colnames(kidneytftargetintrongenesigcount) <- "Count.Significant"

kidneytftargetexongenesigcount <- matrix(0L,nrow = length(kidneytfproandprophlist),ncol = 1)
rownames(kidneytftargetexongenesigcount) <- kidneytfproandprophlist
colnames(kidneytftargetexongenesigcount) <- "Count.Significant"

kidneytftargetintergenicgenesigcount <- matrix(0L,nrow = length(kidneytfproandprophlist),ncol = 1)
rownames(kidneytftargetintergenicgenesigcount) <- kidneytfproandprophlist
colnames(kidneytftargetintergenicgenesigcount) <- "Count.Significant"


for(i in 1:length(kidneytfproandprophlist)){
  ourtf <- kidneytfproandprophlist[i]
  ourtftargetpeaks <- kidney50peakmotifs[kidney50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
  ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
  ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
  ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
  ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
  ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])
  kidneytftargetpromproxgenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetpromproxgenes,kidneyrnasig))
  kidneytftargetintrongenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetintrongenes,kidneyrnasig))
  kidneytftargetexongenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetexongenes,kidneyrnasig))
  kidneytftargetintergenicgenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetintergenicgenes,kidneyrnasig))
  kidneytftargetpromproxgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetpromproxgenes,kidneyrnasig))/length(ourtftargetpromproxgenes)
  kidneytftargetintrongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintrongenes,kidneyrnasig))/length(ourtftargetintrongenes)
  kidneytftargetexongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetexongenes,kidneyrnasig))/length(ourtftargetexongenes)
  kidneytftargetintergenicgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintergenicgenes,kidneyrnasig))/length(ourtftargetintergenicgenes)
}


# LIVER

livertftargetpromproxgenesigpct <- matrix(0L,nrow = length(livertfproandprophlist),ncol = 1)
rownames(livertftargetpromproxgenesigpct) <- livertfproandprophlist
colnames(livertftargetpromproxgenesigpct) <- "Percent.Significant"

livertftargetintrongenesigpct <- matrix(0L,nrow = length(livertfproandprophlist),ncol = 1)
rownames(livertftargetintrongenesigpct) <- livertfproandprophlist
colnames(livertftargetintrongenesigpct) <- "Percent.Significant"

livertftargetexongenesigpct <- matrix(0L,nrow = length(livertfproandprophlist),ncol = 1)
rownames(livertftargetexongenesigpct) <- livertfproandprophlist
colnames(livertftargetexongenesigpct) <- "Percent.Significant"

livertftargetintergenicgenesigpct <- matrix(0L,nrow = length(livertfproandprophlist),ncol = 1)
rownames(livertftargetintergenicgenesigpct) <- livertfproandprophlist
colnames(livertftargetintergenicgenesigpct) <- "Percent.Significant"



livertftargetpromproxgenesigcount <- matrix(0L,nrow = length(livertfproandprophlist),ncol = 1)
rownames(livertftargetpromproxgenesigcount) <- livertfproandprophlist
colnames(livertftargetpromproxgenesigcount) <- "Count.Significant"

livertftargetintrongenesigcount <- matrix(0L,nrow = length(livertfproandprophlist),ncol = 1)
rownames(livertftargetintrongenesigcount) <- livertfproandprophlist
colnames(livertftargetintrongenesigcount) <- "Count.Significant"

livertftargetexongenesigcount <- matrix(0L,nrow = length(livertfproandprophlist),ncol = 1)
rownames(livertftargetexongenesigcount) <- livertfproandprophlist
colnames(livertftargetexongenesigcount) <- "Count.Significant"

livertftargetintergenicgenesigcount <- matrix(0L,nrow = length(livertfproandprophlist),ncol = 1)
rownames(livertftargetintergenicgenesigcount) <- livertfproandprophlist
colnames(livertftargetintergenicgenesigcount) <- "Count.Significant"


for(i in 1:length(livertfproandprophlist)){
  ourtf <- livertfproandprophlist[i]
  ourtftargetpeaks <- liver50peakmotifs[liver50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
  ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
  ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
  ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
  ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
  ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])
  livertftargetpromproxgenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetpromproxgenes,liverrnasig))
  livertftargetintrongenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetintrongenes,liverrnasig))
  livertftargetexongenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetexongenes,liverrnasig))
  livertftargetintergenicgenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetintergenicgenes,liverrnasig))
  livertftargetpromproxgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetpromproxgenes,liverrnasig))/length(ourtftargetpromproxgenes)
  livertftargetintrongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintrongenes,liverrnasig))/length(ourtftargetintrongenes)
  livertftargetexongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetexongenes,liverrnasig))/length(ourtftargetexongenes)
  livertftargetintergenicgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintergenicgenes,liverrnasig))/length(ourtftargetintergenicgenes)
}


# LUNG

lungtftargetpromproxgenesigpct <- matrix(0L,nrow = length(lungtfproandprophlist),ncol = 1)
rownames(lungtftargetpromproxgenesigpct) <- lungtfproandprophlist
colnames(lungtftargetpromproxgenesigpct) <- "Percent.Significant"

lungtftargetintrongenesigpct <- matrix(0L,nrow = length(lungtfproandprophlist),ncol = 1)
rownames(lungtftargetintrongenesigpct) <- lungtfproandprophlist
colnames(lungtftargetintrongenesigpct) <- "Percent.Significant"

lungtftargetexongenesigpct <- matrix(0L,nrow = length(lungtfproandprophlist),ncol = 1)
rownames(lungtftargetexongenesigpct) <- lungtfproandprophlist
colnames(lungtftargetexongenesigpct) <- "Percent.Significant"

lungtftargetintergenicgenesigpct <- matrix(0L,nrow = length(lungtfproandprophlist),ncol = 1)
rownames(lungtftargetintergenicgenesigpct) <- lungtfproandprophlist
colnames(lungtftargetintergenicgenesigpct) <- "Percent.Significant"



lungtftargetpromproxgenesigcount <- matrix(0L,nrow = length(lungtfproandprophlist),ncol = 1)
rownames(lungtftargetpromproxgenesigcount) <- lungtfproandprophlist
colnames(lungtftargetpromproxgenesigcount) <- "Count.Significant"

lungtftargetintrongenesigcount <- matrix(0L,nrow = length(lungtfproandprophlist),ncol = 1)
rownames(lungtftargetintrongenesigcount) <- lungtfproandprophlist
colnames(lungtftargetintrongenesigcount) <- "Count.Significant"

lungtftargetexongenesigcount <- matrix(0L,nrow = length(lungtfproandprophlist),ncol = 1)
rownames(lungtftargetexongenesigcount) <- lungtfproandprophlist
colnames(lungtftargetexongenesigcount) <- "Count.Significant"

lungtftargetintergenicgenesigcount <- matrix(0L,nrow = length(lungtfproandprophlist),ncol = 1)
rownames(lungtftargetintergenicgenesigcount) <- lungtfproandprophlist
colnames(lungtftargetintergenicgenesigcount) <- "Count.Significant"


for(i in 1:length(lungtfproandprophlist)){
  ourtf <- lungtfproandprophlist[i]
  ourtftargetpeaks <- lung50peakmotifs[lung50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
  ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
  ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
  ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
  ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
  ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])
  lungtftargetpromproxgenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetpromproxgenes,lungrnasig))
  lungtftargetintrongenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetintrongenes,lungrnasig))
  lungtftargetexongenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetexongenes,lungrnasig))
  lungtftargetintergenicgenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetintergenicgenes,lungrnasig))
  lungtftargetpromproxgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetpromproxgenes,lungrnasig))/length(ourtftargetpromproxgenes)
  lungtftargetintrongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintrongenes,lungrnasig))/length(ourtftargetintrongenes)
  lungtftargetexongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetexongenes,lungrnasig))/length(ourtftargetexongenes)
  lungtftargetintergenicgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintergenicgenes,lungrnasig))/length(ourtftargetintergenicgenes)
}


# WAT-SC

whitetftargetpromproxgenesigpct <- matrix(0L,nrow = length(whitetfproandprophlist),ncol = 1)
rownames(whitetftargetpromproxgenesigpct) <- whitetfproandprophlist
colnames(whitetftargetpromproxgenesigpct) <- "Percent.Significant"

whitetftargetintrongenesigpct <- matrix(0L,nrow = length(whitetfproandprophlist),ncol = 1)
rownames(whitetftargetintrongenesigpct) <- whitetfproandprophlist
colnames(whitetftargetintrongenesigpct) <- "Percent.Significant"

whitetftargetexongenesigpct <- matrix(0L,nrow = length(whitetfproandprophlist),ncol = 1)
rownames(whitetftargetexongenesigpct) <- whitetfproandprophlist
colnames(whitetftargetexongenesigpct) <- "Percent.Significant"

whitetftargetintergenicgenesigpct <- matrix(0L,nrow = length(whitetfproandprophlist),ncol = 1)
rownames(whitetftargetintergenicgenesigpct) <- whitetfproandprophlist
colnames(whitetftargetintergenicgenesigpct) <- "Percent.Significant"



whitetftargetpromproxgenesigcount <- matrix(0L,nrow = length(whitetfproandprophlist),ncol = 1)
rownames(whitetftargetpromproxgenesigcount) <- whitetfproandprophlist
colnames(whitetftargetpromproxgenesigcount) <- "Count.Significant"

whitetftargetintrongenesigcount <- matrix(0L,nrow = length(whitetfproandprophlist),ncol = 1)
rownames(whitetftargetintrongenesigcount) <- whitetfproandprophlist
colnames(whitetftargetintrongenesigcount) <- "Count.Significant"

whitetftargetexongenesigcount <- matrix(0L,nrow = length(whitetfproandprophlist),ncol = 1)
rownames(whitetftargetexongenesigcount) <- whitetfproandprophlist
colnames(whitetftargetexongenesigcount) <- "Count.Significant"

whitetftargetintergenicgenesigcount <- matrix(0L,nrow = length(whitetfproandprophlist),ncol = 1)
rownames(whitetftargetintergenicgenesigcount) <- whitetfproandprophlist
colnames(whitetftargetintergenicgenesigcount) <- "Count.Significant"


for(i in 1:length(whitetfproandprophlist)){
  ourtf <- whitetfproandprophlist[i]
  ourtftargetpeaks <- white50peakmotifs[white50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
  ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
  ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
  ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
  ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
  ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])
  whitetftargetpromproxgenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetpromproxgenes,whiternasig))
  whitetftargetintrongenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetintrongenes,whiternasig))
  whitetftargetexongenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetexongenes,whiternasig))
  whitetftargetintergenicgenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetintergenicgenes,whiternasig))
  whitetftargetpromproxgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetpromproxgenes,whiternasig))/length(ourtftargetpromproxgenes)
  whitetftargetintrongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintrongenes,whiternasig))/length(ourtftargetintrongenes)
  whitetftargetexongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetexongenes,whiternasig))/length(ourtftargetexongenes)
  whitetftargetintergenicgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintergenicgenes,whiternasig))/length(ourtftargetintergenicgenes)
}

gastroprotargetsigdf <- data.frame("Percent.Significant" = c(gastrotftargetpromproxgenesigpct[gastrotftargetpromproxgenesigcount > 3],
                                                             gastrotftargetintrongenesigpct[gastrotftargetpromproxgenesigcount > 3],
                                                             gastrotftargetexongenesigpct[gastrotftargetpromproxgenesigcount > 3],
                                                             gastrotftargetintergenicgenesigpct[gastrotftargetpromproxgenesigcount > 3]),
                                   "Region" = c(rep("Promoter (<=1kb)",length(gastrotftargetpromproxgenesigpct[gastrotftargetpromproxgenesigcount > 3])),
                                                rep("Intron",length(gastrotftargetpromproxgenesigpct[gastrotftargetpromproxgenesigcount > 3])),
                                                rep("Exon",length(gastrotftargetpromproxgenesigpct[gastrotftargetpromproxgenesigcount > 3])),
                                                rep("Distal Intergenic",length(gastrotftargetpromproxgenesigpct[gastrotftargetpromproxgenesigcount > 3]))),
                                   "Transcription Factor" = rep(gsub("\\(.*","",rownames(gastrotftargetpromproxgenesigpct)[gastrotftargetpromproxgenesigcount > 3]),4))

heartprotargetsigdf <- data.frame("Percent.Significant" = c(hearttftargetpromproxgenesigpct[hearttftargetpromproxgenesigcount > 3],
                                                            hearttftargetintrongenesigpct[hearttftargetpromproxgenesigcount > 3],
                                                            hearttftargetexongenesigpct[hearttftargetpromproxgenesigcount > 3],
                                                            hearttftargetintergenicgenesigpct[hearttftargetpromproxgenesigcount > 3]),
                                  "Region" = c(rep("Promoter (<=1kb)",length(hearttftargetpromproxgenesigpct[hearttftargetpromproxgenesigcount > 3])),
                                               rep("Intron",length(hearttftargetpromproxgenesigpct[hearttftargetpromproxgenesigcount > 3])),
                                               rep("Exon",length(hearttftargetpromproxgenesigpct[hearttftargetpromproxgenesigcount > 3])),
                                               rep("Distal Intergenic",length(hearttftargetpromproxgenesigpct[hearttftargetpromproxgenesigcount > 3]))),
                                  "Transcription Factor" = rep(gsub("\\(.*","",rownames(hearttftargetpromproxgenesigpct)[hearttftargetpromproxgenesigcount > 3]),4))


kidneyprotargetsigdf <- data.frame("Percent.Significant" = c(kidneytftargetpromproxgenesigpct[kidneytftargetpromproxgenesigcount > 3],
                                                             kidneytftargetintrongenesigpct[kidneytftargetpromproxgenesigcount > 3],
                                                             kidneytftargetexongenesigpct[kidneytftargetpromproxgenesigcount > 3],
                                                             kidneytftargetintergenicgenesigpct[kidneytftargetpromproxgenesigcount > 3]),
                                   "Region" = c(rep("Promoter (<=1kb)",length(kidneytftargetpromproxgenesigpct[kidneytftargetpromproxgenesigcount > 3])),
                                                rep("Intron",length(kidneytftargetpromproxgenesigpct[kidneytftargetpromproxgenesigcount > 3])),
                                                rep("Exon",length(kidneytftargetpromproxgenesigpct[kidneytftargetpromproxgenesigcount > 3])),
                                                rep("Distal Intergenic",length(kidneytftargetpromproxgenesigpct[kidneytftargetpromproxgenesigcount > 3]))),
                                   "Transcription Factor" = rep(gsub("\\(.*","",rownames(kidneytftargetpromproxgenesigpct)[kidneytftargetpromproxgenesigcount > 3]),4))


liverprotargetsigdf <- data.frame("Percent.Significant" = c(livertftargetpromproxgenesigpct[livertftargetpromproxgenesigcount > 3],
                                                            livertftargetintrongenesigpct[livertftargetpromproxgenesigcount > 3],
                                                            livertftargetexongenesigpct[livertftargetpromproxgenesigcount > 3],
                                                            livertftargetintergenicgenesigpct[livertftargetpromproxgenesigcount > 3]),
                                  "Region" = c(rep("Promoter (<=1kb)",length(livertftargetpromproxgenesigpct[livertftargetpromproxgenesigcount > 3])),
                                               rep("Intron",length(livertftargetpromproxgenesigpct[livertftargetpromproxgenesigcount > 3])),
                                               rep("Exon",length(livertftargetpromproxgenesigpct[livertftargetpromproxgenesigcount > 3])),
                                               rep("Distal Intergenic",length(livertftargetpromproxgenesigpct[livertftargetpromproxgenesigcount > 3]))),
                                  "Transcription Factor" = rep(gsub("\\(.*","",rownames(livertftargetpromproxgenesigpct)[livertftargetpromproxgenesigcount > 3]),4))


lungprotargetsigdf <- data.frame("Percent.Significant" = c(lungtftargetpromproxgenesigpct[lungtftargetpromproxgenesigcount > 3],
                                                           lungtftargetintrongenesigpct[lungtftargetpromproxgenesigcount > 3],
                                                           lungtftargetexongenesigpct[lungtftargetpromproxgenesigcount > 3],
                                                           lungtftargetintergenicgenesigpct[lungtftargetpromproxgenesigcount > 3]),
                                 "Region" = c(rep("Promoter (<=1kb)",length(lungtftargetpromproxgenesigpct[lungtftargetpromproxgenesigcount > 3])),
                                              rep("Intron",length(lungtftargetpromproxgenesigpct[lungtftargetpromproxgenesigcount > 3])),
                                              rep("Exon",length(lungtftargetpromproxgenesigpct[lungtftargetpromproxgenesigcount > 3])),
                                              rep("Distal Intergenic",length(lungtftargetpromproxgenesigpct[lungtftargetpromproxgenesigcount > 3]))),
                                 "Transcription Factor" = rep(gsub("\\(.*","",rownames(lungtftargetpromproxgenesigpct)[lungtftargetpromproxgenesigcount > 3]),4))


whiteprotargetsigdf <- data.frame("Percent.Significant" = c(whitetftargetpromproxgenesigpct[whitetftargetpromproxgenesigcount > 3],
                                                            whitetftargetintrongenesigpct[whitetftargetpromproxgenesigcount > 3],
                                                            whitetftargetexongenesigpct[whitetftargetpromproxgenesigcount > 3],
                                                            whitetftargetintergenicgenesigpct[whitetftargetpromproxgenesigcount > 3]),
                                  "Region" = c(rep("Promoter (<=1kb)",length(whitetftargetpromproxgenesigpct[whitetftargetpromproxgenesigcount > 3])),
                                               rep("Intron",length(whitetftargetpromproxgenesigpct[whitetftargetpromproxgenesigcount > 3])),
                                               rep("Exon",length(whitetftargetpromproxgenesigpct[whitetftargetpromproxgenesigcount > 3])),
                                               rep("Distal Intergenic",length(whitetftargetpromproxgenesigpct[whitetftargetpromproxgenesigcount > 3]))),
                                  "Transcription Factor" = rep(gsub("\\(.*","",rownames(whitetftargetpromproxgenesigpct)[whitetftargetpromproxgenesigcount > 3]),4))


whiteprotargetsigdf$Transcription.Factor <- gsub("\\/.*","",whiteprotargetsigdf$Transcription.Factor)


# SKM-GN

gastrotfproandprophmeta$Proph.Sig <- gastrotfproandprophmeta$Proph.Sig*3
gastrotfproandprophmeta$Sig.Group <- abs(gastrotfproandprophmeta$Proph.Sig - gastrotfproandprophmeta$Pro.Sig)

gastrotfproandprophmetatrim <- gastrotfproandprophmeta[gastrotftargetpromproxgenesigcount > 3,]

gastroprotargetsigdf$Transcription.Factor <- factor(gastroprotargetsigdf$Transcription.Factor,levels = c(gsub("\\(.*","",rownames(gastrotfproandprophmetatrim[order(gastrotfproandprophmetatrim$Sig.Group),]))))

gastrotable <- data.frame(start = c(0.5,sum(gastrotfproandprophmetatrim$Sig.Group <= 1)+0.5,sum(gastrotfproandprophmetatrim$Sig.Group <= 2)+0.5),
                          end = c(sum(gastrotfproandprophmetatrim$Sig.Group <= 1)+0.5,sum(gastrotfproandprophmetatrim$Sig.Group <= 2)+0.5,sum(gastrotfproandprophmetatrim$Sig.Group <= 3)+0.5),
                          Sig.Group = factor(c("Pro","Pro+Phospho","Phospho"),levels = c("Pro","Pro+Phospho","Phospho")))


gastroprotargetsigdf$Transcription.Factor.Label = gastroprotargetsigdf$Transcription.Factor
gastroprotargetsigdf$Transcription.Factor <- rep(order(gastrotfproandprophmetatrim$Sig.Group),4)

png(file = "Supplemental Figure S14A.png",width = 9,height = 4.5,units = "in",res = 600)
ggplot(gastroprotargetsigdf[!gastroprotargetsigdf$Region %in% "Exon",],
       aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + 
  geom_hline(yintercept = length(gastrornasig)/dim(gastrol2fcmat)[1],linetype = "dashed",size = 1) + 
  geom_point(aes(color=Region),size = 3) + 
  geom_line(aes(color=Region),size = 2) + 
  geom_rect(data = gastrotable,aes(xmin = start,xmax = end,ymin = -Inf,ymax = 0,fill = Sig.Group),alpha = 0.5,inherit.aes = F) + 
  theme_classic() + 
  scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + 
  scale_x_continuous(breaks = c(1:max(gastroprotargetsigdf$Transcription.Factor)),
                     labels = levels(gastroprotargetsigdf$Transcription.Factor.Label)) + 
  theme(legend.position = "top",
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.box = "vertical")
dev.off()

pdf(file = "Supplemental Figure S14A.pdf",width = 9,height = 4.5)
ggplot(gastroprotargetsigdf[!gastroprotargetsigdf$Region %in% "Exon",],
       aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + 
  geom_hline(yintercept = length(gastrornasig)/dim(gastrol2fcmat)[1],linetype = "dashed",size = 1) + 
  geom_point(aes(color=Region),size = 3) + 
  geom_line(aes(color=Region),size = 2) + 
  geom_rect(data = gastrotable,aes(xmin = start,xmax = end,ymin = -Inf,ymax = 0,fill = Sig.Group),alpha = 0.5,inherit.aes = F) + 
  theme_classic() + 
  scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + 
  scale_x_continuous(breaks = c(1:max(gastroprotargetsigdf$Transcription.Factor)),
                     labels = levels(gastroprotargetsigdf$Transcription.Factor.Label)) + 
  theme(legend.position = "top",
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.box = "vertical")
dev.off()


# HEART

hearttfproandprophmeta$Proph.Sig <- hearttfproandprophmeta$Proph.Sig*3
hearttfproandprophmeta$Sig.Group <- abs(hearttfproandprophmeta$Proph.Sig - hearttfproandprophmeta$Pro.Sig)

hearttfproandprophmetatrim <- hearttfproandprophmeta[hearttftargetpromproxgenesigcount > 3,]

heartprotargetsigdf$Transcription.Factor <- factor(heartprotargetsigdf$Transcription.Factor,levels = c(gsub("\\(.*","",rownames(hearttfproandprophmetatrim[order(hearttfproandprophmetatrim$Sig.Group),]))))


hearttable <- data.frame(start = c(0.5,sum(hearttfproandprophmetatrim$Sig.Group <= 1)+0.5,sum(hearttfproandprophmetatrim$Sig.Group <= 2)+0.5),
                         end = c(sum(hearttfproandprophmetatrim$Sig.Group <= 1)+0.5,sum(hearttfproandprophmetatrim$Sig.Group <= 2)+0.5,sum(hearttfproandprophmetatrim$Sig.Group <= 3)+0.5),
                         Sig.Group = factor(c("Pro","Pro+Phospho","Phospho"),levels = c("Pro","Pro+Phospho","Phospho")))


heartprotargetsigdf$Transcription.Factor.Label = heartprotargetsigdf$Transcription.Factor
heartprotargetsigdf$Transcription.Factor <- rep(order(hearttfproandprophmetatrim$Sig.Group),4)

png(file = "Supplemental Figure S14B.png",width = 9,height = 4.5,units = "in",res = 600)
ggplot(heartprotargetsigdf[!heartprotargetsigdf$Region %in% "Exon",],
       aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + 
  geom_hline(yintercept = length(heartrnasig)/dim(heartl2fcmat)[1],linetype = "dashed",size = 1) + 
  geom_point(aes(color=Region),size = 3) + 
  geom_line(aes(color=Region),size = 2) + 
  geom_rect(data = hearttable,aes(xmin = start,xmax = end,ymin = -Inf,ymax = 0,fill = Sig.Group),alpha = 0.5,inherit.aes = F) + 
  theme_classic() + 
  scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + 
  scale_x_continuous(breaks = c(1:max(heartprotargetsigdf$Transcription.Factor)),
                     labels = levels(heartprotargetsigdf$Transcription.Factor.Label)) + 
  theme(legend.position = "top",
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.box = "vertical")
dev.off()

pdf(file = "Supplemental Figure S14B.pdf",width = 9,height = 4.5)
ggplot(heartprotargetsigdf[!heartprotargetsigdf$Region %in% "Exon",],
       aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + 
  geom_hline(yintercept = length(heartrnasig)/dim(heartl2fcmat)[1],linetype = "dashed",size = 1) + 
  geom_point(aes(color=Region),size = 3) + 
  geom_line(aes(color=Region),size = 2) + 
  geom_rect(data = hearttable,aes(xmin = start,xmax = end,ymin = -Inf,ymax = 0,fill = Sig.Group),alpha = 0.5,inherit.aes = F) + 
  theme_classic() + 
  scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + 
  scale_x_continuous(breaks = c(1:max(heartprotargetsigdf$Transcription.Factor)),
                     labels = levels(heartprotargetsigdf$Transcription.Factor.Label)) + 
  theme(legend.position = "top",
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.box = "vertical")
dev.off()


# KIDNEY

kidneytfproandprophmeta$Proph.Sig <- kidneytfproandprophmeta$Proph.Sig*3
kidneytfproandprophmeta$Sig.Group <- abs(kidneytfproandprophmeta$Proph.Sig - kidneytfproandprophmeta$Pro.Sig)

kidneytfproandprophmetatrim <- kidneytfproandprophmeta[kidneytftargetpromproxgenesigcount > 3,]

kidneyprotargetsigdf$Transcription.Factor <- factor(kidneyprotargetsigdf$Transcription.Factor,levels = c(gsub("\\(.*","",rownames(kidneytfproandprophmetatrim[order(kidneytfproandprophmetatrim$Sig.Group),]))))


kidneytable <- data.frame(start = c(0.5,sum(kidneytfproandprophmetatrim$Sig.Group <= 1)+0.5,sum(kidneytfproandprophmetatrim$Sig.Group <= 2)+0.5),
                          end = c(sum(kidneytfproandprophmetatrim$Sig.Group <= 1)+0.5,sum(kidneytfproandprophmetatrim$Sig.Group <= 2)+0.5,sum(kidneytfproandprophmetatrim$Sig.Group <= 3)+0.5),
                          Sig.Group = factor(c("Pro","Pro+Phospho","Phospho"),levels = c("Pro","Pro+Phospho","Phospho")))


kidneyprotargetsigdf$Transcription.Factor.Label = kidneyprotargetsigdf$Transcription.Factor
kidneyprotargetsigdf$Transcription.Factor <- rep(order(kidneytfproandprophmetatrim$Sig.Group),4)

png(file = "Supplemental Figure S14C.png",width = 9,height = 4.5,units = "in",res = 600)
ggplot(kidneyprotargetsigdf[!kidneyprotargetsigdf$Region %in% "Exon",],
       aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + 
  geom_hline(yintercept = length(kidneyrnasig)/dim(kidneyl2fcmat)[1],linetype = "dashed",size = 1) + 
  geom_point(aes(color=Region),size = 3) + 
  geom_line(aes(color=Region),size = 2) + 
  geom_rect(data = kidneytable,aes(xmin = start,xmax = end,ymin = -Inf,ymax = 0,fill = Sig.Group),alpha = 0.5,inherit.aes = F) + 
  theme_classic() + 
  scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + 
  scale_x_continuous(breaks = c(1:max(kidneyprotargetsigdf$Transcription.Factor)),
                     labels = levels(kidneyprotargetsigdf$Transcription.Factor.Label)) + 
  theme(legend.position = "top",
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.box = "vertical")
dev.off()

pdf(file = "Supplemental Figure S14C.pdf",width = 9,height = 4.5)
ggplot(kidneyprotargetsigdf[!kidneyprotargetsigdf$Region %in% "Exon",],
       aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + 
  geom_hline(yintercept = length(kidneyrnasig)/dim(kidneyl2fcmat)[1],linetype = "dashed",size = 1) + 
  geom_point(aes(color=Region),size = 3) + 
  geom_line(aes(color=Region),size = 2) + 
  geom_rect(data = kidneytable,aes(xmin = start,xmax = end,ymin = -Inf,ymax = 0,fill = Sig.Group),alpha = 0.5,inherit.aes = F) + 
  theme_classic() + 
  scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + 
  scale_x_continuous(breaks = c(1:max(kidneyprotargetsigdf$Transcription.Factor)),
                     labels = levels(kidneyprotargetsigdf$Transcription.Factor.Label)) + 
  theme(legend.position = "top",
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.box = "vertical")
dev.off()


# LIVER

livertfproandprophmeta$Proph.Sig <- livertfproandprophmeta$Proph.Sig*3
livertfproandprophmeta$Sig.Group <- abs(livertfproandprophmeta$Proph.Sig - livertfproandprophmeta$Pro.Sig)

livertfproandprophmetatrim <- livertfproandprophmeta[livertftargetpromproxgenesigcount > 3,]

liverprotargetsigdf$Transcription.Factor <- factor(liverprotargetsigdf$Transcription.Factor,levels = c(gsub("\\(.*","",rownames(livertfproandprophmetatrim[order(livertfproandprophmetatrim$Sig.Group),]))))


livertable <- data.frame(start = c(0.5,sum(livertfproandprophmetatrim$Sig.Group <= 1)+0.5,sum(livertfproandprophmetatrim$Sig.Group <= 2)+0.5),
                         end = c(sum(livertfproandprophmetatrim$Sig.Group <= 1)+0.5,sum(livertfproandprophmetatrim$Sig.Group <= 2)+0.5,sum(livertfproandprophmetatrim$Sig.Group <= 3)+0.5),
                         Sig.Group = factor(c("Pro","Pro+Phospho","Phospho"),levels = c("Pro","Pro+Phospho","Phospho")))


liverprotargetsigdf$Transcription.Factor.Label = liverprotargetsigdf$Transcription.Factor
liverprotargetsigdf$Transcription.Factor <- rep(order(livertfproandprophmetatrim$Sig.Group),4)

png(file = "Supplemental Figure S14D.png",width = 9,height = 4.5,units = "in",res = 600)
ggplot(liverprotargetsigdf[!liverprotargetsigdf$Region %in% "Exon",],
       aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + 
  geom_hline(yintercept = length(liverrnasig)/dim(liverl2fcmat)[1],linetype = "dashed",size = 1) + 
  geom_point(aes(color=Region),size = 3) + 
  geom_line(aes(color=Region),size = 2) + 
  geom_rect(data = livertable,aes(xmin = start,xmax = end,ymin = -Inf,ymax = 0,fill = Sig.Group),alpha = 0.5,inherit.aes = F) + 
  theme_classic() + 
  scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + 
  scale_x_continuous(breaks = c(1:max(liverprotargetsigdf$Transcription.Factor)),
                     labels = levels(liverprotargetsigdf$Transcription.Factor.Label)) + 
  theme(legend.position = "top",
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.box = "vertical")
dev.off()

pdf(file = "Supplemental Figure S14D.pdf",width = 9,height = 4.5)
ggplot(liverprotargetsigdf[!liverprotargetsigdf$Region %in% "Exon",],
       aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + 
  geom_hline(yintercept = length(liverrnasig)/dim(liverl2fcmat)[1],linetype = "dashed",size = 1) + 
  geom_point(aes(color=Region),size = 3) + 
  geom_line(aes(color=Region),size = 2) + 
  geom_rect(data = livertable,aes(xmin = start,xmax = end,ymin = -Inf,ymax = 0,fill = Sig.Group),alpha = 0.5,inherit.aes = F) + 
  theme_classic() + 
  scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + 
  scale_x_continuous(breaks = c(1:max(liverprotargetsigdf$Transcription.Factor)),
                     labels = levels(liverprotargetsigdf$Transcription.Factor.Label)) + 
  theme(legend.position = "top",
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.box = "vertical")
dev.off()


# LUNG

lungtfproandprophmeta$Proph.Sig <- lungtfproandprophmeta$Proph.Sig*3
lungtfproandprophmeta$Sig.Group <- abs(lungtfproandprophmeta$Proph.Sig - lungtfproandprophmeta$Pro.Sig)

lungtfproandprophmetatrim <- lungtfproandprophmeta[lungtftargetpromproxgenesigcount > 3,]

lungprotargetsigdf$Transcription.Factor <- factor(lungprotargetsigdf$Transcription.Factor,levels = c(gsub("\\(.*","",rownames(lungtfproandprophmetatrim[order(lungtfproandprophmetatrim$Sig.Group),]))))

lungtable <- data.frame(start = c(0.5,sum(lungtfproandprophmetatrim$Sig.Group <= 1)+0.5,sum(lungtfproandprophmetatrim$Sig.Group <= 2)+0.5),
                        end = c(sum(lungtfproandprophmetatrim$Sig.Group <= 1)+0.5,sum(lungtfproandprophmetatrim$Sig.Group <= 2)+0.5,sum(lungtfproandprophmetatrim$Sig.Group <= 3)+0.5),
                        Sig.Group = factor(c("Pro","Pro+Phospho","Phospho"),levels = c("Pro","Pro+Phospho","Phospho")))


lungprotargetsigdf$Transcription.Factor.Label = lungprotargetsigdf$Transcription.Factor
lungprotargetsigdf$Transcription.Factor <- rep(order(lungtfproandprophmetatrim$Sig.Group),4)

png(file = "Supplemental Figure S14E.png",width = 9,height = 4.5,units = "in",res = 600)
ggplot(lungprotargetsigdf[!lungprotargetsigdf$Region %in% "Exon",],
       aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + 
  geom_hline(yintercept = length(lungrnasig)/dim(lungl2fcmat)[1],linetype = "dashed",size = 1) + 
  geom_point(aes(color=Region),size = 3) + 
  geom_line(aes(color=Region),size = 2) + 
  geom_rect(data = lungtable,aes(xmin = start,xmax = end,ymin = -Inf,ymax = 0,fill = Sig.Group),alpha = 0.5,inherit.aes = F) + 
  theme_classic() + 
  scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + 
  scale_x_continuous(breaks = c(1:max(lungprotargetsigdf$Transcription.Factor)),
                     labels = levels(lungprotargetsigdf$Transcription.Factor.Label)) + 
  theme(legend.position = "top",
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.box = "vertical")
dev.off()

pdf(file = "Supplemental Figure S14E.pdf",width = 9,height = 4.5)
ggplot(lungprotargetsigdf[!lungprotargetsigdf$Region %in% "Exon",],
       aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + 
  geom_hline(yintercept = length(lungrnasig)/dim(lungl2fcmat)[1],linetype = "dashed",size = 1) + 
  geom_point(aes(color=Region),size = 3) + 
  geom_line(aes(color=Region),size = 2) + 
  geom_rect(data = lungtable,aes(xmin = start,xmax = end,ymin = -Inf,ymax = 0,fill = Sig.Group),alpha = 0.5,inherit.aes = F) + 
  theme_classic() + 
  scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + 
  scale_x_continuous(breaks = c(1:max(lungprotargetsigdf$Transcription.Factor)),
                     labels = levels(lungprotargetsigdf$Transcription.Factor.Label)) + 
  theme(legend.position = "top",
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.box = "vertical")
dev.off()



# WAT-SC

whitetfproandprophmeta$Proph.Sig <- whitetfproandprophmeta$Proph.Sig*3
whitetfproandprophmeta$Sig.Group <- abs(whitetfproandprophmeta$Proph.Sig - whitetfproandprophmeta$Pro.Sig)

whitetfproandprophmetatrim <- whitetfproandprophmeta[whitetftargetpromproxgenesigcount > 3,]

whiteprotargetsigdf$Transcription.Factor <- factor(whiteprotargetsigdf$Transcription.Factor,levels = c(gsub("\\/.*","",gsub("\\(.*","",rownames(whitetfproandprophmetatrim[order(whitetfproandprophmetatrim$Sig.Group),])))))

whitetable <- data.frame(start = c(0.5,sum(whitetfproandprophmetatrim$Sig.Group <= 1)+0.5,sum(whitetfproandprophmetatrim$Sig.Group <= 2)+0.5),
                         end = c(sum(whitetfproandprophmetatrim$Sig.Group <= 1)+0.5,sum(whitetfproandprophmetatrim$Sig.Group <= 2)+0.5,sum(whitetfproandprophmetatrim$Sig.Group <= 3)+0.5),
                         Sig.Group = factor(c("Pro","Pro+Phospho","Phospho"),levels = c("Pro","Pro+Phospho","Phospho")))


whiteprotargetsigdf$Transcription.Factor.Label = whiteprotargetsigdf$Transcription.Factor
whiteprotargetsigdf$Transcription.Factor <- rep(order(whitetfproandprophmetatrim$Sig.Group),4)

png(file = "Supplemental Figure S14F.png",width = 9,height = 4.5,units = "in",res = 600)
ggplot(whiteprotargetsigdf[!whiteprotargetsigdf$Region %in% "Exon",],
       aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + 
  geom_hline(yintercept = length(whiternasig)/dim(whitel2fcmat)[1],linetype = "dashed",size = 1) + 
  geom_point(aes(color=Region),size = 3) + 
  geom_line(aes(color=Region),size = 2) + 
  geom_rect(data = whitetable,aes(xmin = start,xmax = end,ymin = -Inf,ymax = 0,fill = Sig.Group),alpha = 0.5,inherit.aes = F) + 
  theme_classic() + 
  scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + 
  scale_x_continuous(breaks = c(1:max(whiteprotargetsigdf$Transcription.Factor)),
                     labels = levels(whiteprotargetsigdf$Transcription.Factor.Label)) + 
  theme(legend.position = "top",
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.box = "vertical")
dev.off()

pdf(file = "Supplemental Figure S14F.pdf",width = 9,height = 4.5)
ggplot(whiteprotargetsigdf[!whiteprotargetsigdf$Region %in% "Exon",],
       aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + 
  geom_hline(yintercept = length(whiternasig)/dim(whitel2fcmat)[1],linetype = "dashed",size = 1) + 
  geom_point(aes(color=Region),size = 3) + 
  geom_line(aes(color=Region),size = 2) + 
  geom_rect(data = whitetable,aes(xmin = start,xmax = end,ymin = -Inf,ymax = 0,fill = Sig.Group),alpha = 0.5,inherit.aes = F) + 
  theme_classic() + 
  scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + 
  scale_x_continuous(breaks = c(1:max(whiteprotargetsigdf$Transcription.Factor)),
                     labels = levels(whiteprotargetsigdf$Transcription.Factor.Label)) + 
  theme(legend.position = "top",
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.box = "vertical")
dev.off()


#####
# Figure 5F, 5G, and Supplemental Figure S17
####

ann_cols_procor <- list("Tissue" = c("SKM-GN" = "#088c03",
                                     "HEART" = "#f28b2f",
                                     "KIDNEY"= "#7553a7",
                                     "LIVER" = "#da6c75",
                                     "LUNG" = "#04bf8a",
                                     "WAT-SC" = "#214da6"),
                        "Sex" = c("Female" = "#ff6eff",
                                  "Male" = "#5555ff"),
                        "Week" = c("control" = "white",
                                   "1w" = "#F7FCB9",
                                   "2w" = "#ADDD8E",
                                   "4w" = "#238443",
                                   "8w" = "#002612"),
                        "Region" = c("Distal Intergenic" = "#BCBD22FF",
                                     "Intron" = "#8C564BFF",
                                     "Promoter" = "#FF7F0EFF",
                                     "All" = "navy"),
                        "Training.Response" = c("Down.Reg" = "skyblue2",
                                                "Up.Reg" = "indianred2"),
                        "Correlation" = c("1" = "red4",
                                          "0.8" = "red3",
                                          "0.6" = "red1",
                                          "0.4" = "pink2",
                                          "0.2" = "pink",
                                          "0" = "white",
                                          "-0.2" = "lightblue",
                                          "-0.4" = "lightblue3",
                                          "-0.6" = "blue1",
                                          "-0.8" = "blue3",
                                          "-1" = "blue4"),
                        "TF.L2FC" = c("1" = "#FF0000",
                                      "0.9" = "#FF1919",
                                      "0.8" = "#FF3232",
                                      "0.7" = "#FF4C4C",
                                      "0.6" = "#FF6565",
                                      "0.5" = "#FF7F7F",
                                      "0.4" = "#FF9898",
                                      "0.3" = "#FFB2B2",
                                      "0.2" = "#FFCBCB",
                                      "0.1" = "#FFE5E5",
                                      "0" = "#FFFFFF",
                                      "-0.1" = "#E5E5FF",
                                      "-0.2" = "#CCCCFF",
                                      "-0.3" = "#B2B2FF",
                                      "-0.4" = "#9999FF",
                                      "-0.5" = "#7F7FFF",
                                      "-0.6" = "#6666FF",
                                      "-0.7" = "#4C4CFF",
                                      "-0.8" = "#3333FF",
                                      "-0.9" = "#1919FF",
                                      "-1" = "#0000FF"))

ourtf <- "Mef2c(MADS)/GM12878-Mef2c-ChIP-Seq(GSE32465)/Homer"
ourtftargetpeaks <- gastro50peakmotifs[gastro50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])

sextimemeta <- data.frame(row.names = colnames(gastrol2fcmat),
                          "Sex" = c(rep("Female",4),rep("Male",4)),
                          "Week" = rep(c("1w","2w","4w","8w"),2),
                          "TF.L2FC" = gastroprol2fc[tfproanno[ourtf,"Gastro.Pro.ID"],])
sextimemeta$TF.L2FC <- round(sextimemeta$TF.L2FC/0.1)*0.1
sextimemeta$TF.L2FC[sextimemeta$TF.L2FC > 1] <- 1
sextimemeta$TF.L2FC[sextimemeta$TF.L2FC < -1] <- -1
sextimemeta$TF.L2FC <- factor(sextimemeta$TF.L2FC,levels = c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8,-0.9,-1))
ourtftargetpromproxcormeta <- data.frame(row.names = intersect(ourtftargetpromproxgenes,gastrornasig),
                                         "Correlation" = rep(0,length(intersect(ourtftargetpromproxgenes,gastrornasig))))
ourtftargetexoncormeta <- data.frame(row.names = intersect(ourtftargetexongenes,gastrornasig),
                                     "Correlation" = rep(0,length(intersect(ourtftargetexongenes,gastrornasig))))
ourtftargetintroncormeta <- data.frame(row.names = intersect(ourtftargetintrongenes,gastrornasig),
                                       "Correlation" = rep(0,length(intersect(ourtftargetintrongenes,gastrornasig))))
for(i in 1:dim(ourtftargetpromproxcormeta)[1]){
  ourtftargetpromproxcormeta[i,"Correlation"] <- cor(gastroprol2fc[tfproanno[ourtf,"Gastro.Pro.ID"],],gastrol2fcmat[intersect(ourtftargetpromproxgenes,gastrornasig)[i],])
}
ourtftargetpromproxcormeta$Correlation <- round(ourtftargetpromproxcormeta$Correlation/0.2)*0.2
ourtftargetpromproxcormeta$Correlation <- factor(ourtftargetpromproxcormeta$Correlation,levels = c(1,0.8,0.6,0.4,0.2,0,-0.2,-0.4,-0.6,-0.8,-1))

png(file = "Figure 5F_Supplemental Figure S17A.png",width = 7,height = 4,units = "in",res = 600)
pheatmap(gastrol2fcmat[intersect(ourtftargetpromproxgenes,gastrornasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,gastrornasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F)
dev.off()

pdf(file = "Figure 5F_Supplemental Figure S17A.pdf",width = 7,height = 4)
pheatmap(gastrol2fcmat[intersect(ourtftargetpromproxgenes,gastrornasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,gastrornasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F)
dev.off()


ourtf <- "Nur77(NR)/K562-NR4A1-ChIP-Seq(GSE31363)/Homer"
ourtftargetpeaks <- gastro50peakmotifs[gastro50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])

sextimemeta <- data.frame(row.names = colnames(gastrol2fcmat),
                          "Sex" = c(rep("Female",4),rep("Male",4)),
                          "Week" = rep(c("1w","2w","4w","8w"),2),
                          "TF.L2FC" = gastroprol2fc[tfproanno[ourtf,"Gastro.Pro.ID"],])
sextimemeta$TF.L2FC <- round(sextimemeta$TF.L2FC/0.1)*0.1
sextimemeta$TF.L2FC[sextimemeta$TF.L2FC > 1] <- 1
sextimemeta$TF.L2FC[sextimemeta$TF.L2FC < -1] <- -1
sextimemeta$TF.L2FC <- factor(sextimemeta$TF.L2FC,levels = c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8,-0.9,-1))
ourtftargetpromproxcormeta <- data.frame(row.names = intersect(ourtftargetpromproxgenes,gastrornasig),
                                         "Correlation" = rep(0,length(intersect(ourtftargetpromproxgenes,gastrornasig))))
ourtftargetexoncormeta <- data.frame(row.names = intersect(ourtftargetexongenes,gastrornasig),
                                     "Correlation" = rep(0,length(intersect(ourtftargetexongenes,gastrornasig))))
ourtftargetintroncormeta <- data.frame(row.names = intersect(ourtftargetintrongenes,gastrornasig),
                                       "Correlation" = rep(0,length(intersect(ourtftargetintrongenes,gastrornasig))))
for(i in 1:dim(ourtftargetpromproxcormeta)[1]){
  ourtftargetpromproxcormeta[i,"Correlation"] <- cor(gastroprol2fc[tfproanno[ourtf,"Gastro.Pro.ID"],],gastrol2fcmat[intersect(ourtftargetpromproxgenes,gastrornasig)[i],])
}
ourtftargetpromproxcormeta$Correlation <- round(ourtftargetpromproxcormeta$Correlation/0.2)*0.2
ourtftargetpromproxcormeta$Correlation <- factor(ourtftargetpromproxcormeta$Correlation,levels = c(1,0.8,0.6,0.4,0.2,0,-0.2,-0.4,-0.6,-0.8,-1))

png(file = "Supplemental Figure S17B.png",width = 7,height = 4,units = "in",res = 600)
pheatmap(gastrol2fcmat[intersect(ourtftargetpromproxgenes,gastrornasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,gastrornasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F)
dev.off()

pdf(file = "Supplemental Figure S17B.pdf",width = 7,height = 4)
pheatmap(gastrol2fcmat[intersect(ourtftargetpromproxgenes,gastrornasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,gastrornasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F)
dev.off()


ourtf <- "IRF:BATF(IRF:bZIP)/pDC-Irf8-ChIP-Seq(GSE66899)/Homer"
ourtftargetpeaks <- lung50peakmotifs[lung50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])

sextimemeta <- data.frame(row.names = colnames(lungl2fcmat),
                          "Sex" = c(rep("Female",4),rep("Male",4)),
                          "Week" = rep(c("1w","2w","4w","8w"),2),
                          "TF.L2FC" = lungprol2fc[tfproanno[ourtf,"Lung.Pro.ID"],])
sextimemeta$TF.L2FC <- round(sextimemeta$TF.L2FC/0.1)*0.1
sextimemeta$TF.L2FC[sextimemeta$TF.L2FC > 1] <- 1
sextimemeta$TF.L2FC[sextimemeta$TF.L2FC < -1] <- -1
sextimemeta$TF.L2FC <- factor(sextimemeta$TF.L2FC,levels = c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8,-0.9,-1))
ourtftargetpromproxcormeta <- data.frame(row.names = intersect(ourtftargetpromproxgenes,lungrnasig),
                                         "Correlation" = rep(0,length(intersect(ourtftargetpromproxgenes,lungrnasig))))
ourtftargetexoncormeta <- data.frame(row.names = intersect(ourtftargetexongenes,lungrnasig),
                                     "Correlation" = rep(0,length(intersect(ourtftargetexongenes,lungrnasig))))
ourtftargetintroncormeta <- data.frame(row.names = intersect(ourtftargetintrongenes,lungrnasig),
                                       "Correlation" = rep(0,length(intersect(ourtftargetintrongenes,lungrnasig))))
for(i in 1:dim(ourtftargetpromproxcormeta)[1]){
  ourtftargetpromproxcormeta[i,"Correlation"] <- cor(lungprol2fc[tfproanno[ourtf,"Lung.Pro.ID"],],lungl2fcmat[intersect(ourtftargetpromproxgenes,lungrnasig)[i],])
}
ourtftargetpromproxcormeta$Correlation <- round(ourtftargetpromproxcormeta$Correlation/0.2)*0.2
ourtftargetpromproxcormeta$Correlation <- factor(ourtftargetpromproxcormeta$Correlation,levels = c(1,0.8,0.6,0.4,0.2,0,-0.2,-0.4,-0.6,-0.8,-1))

png(file = "Figure 5G_Supplemental Figure S17C.png",width = 7,height = 4,units = "in",res = 600)
pheatmap(lungl2fcmat[intersect(ourtftargetpromproxgenes,lungrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,lungrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F)
dev.off()

pdf(file = "Figure 5G_Supplemental Figure S17C.pdf",width = 7,height = 4)
pheatmap(lungl2fcmat[intersect(ourtftargetpromproxgenes,lungrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,lungrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F)
dev.off()


ourtf <- "PBX1(Homeobox)/MCF7-PBX1-ChIP-Seq(GSE28007)/Homer"
ourtftargetpeaks <- white50peakmotifs[white50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])

sextimemeta <- data.frame(row.names = colnames(whitel2fcmat),
                          "Sex" = c(rep("Female",4),rep("Male",4)),
                          "Week" = rep(c("1w","2w","4w","8w"),2),
                          "TF.L2FC" = lungprol2fc[tfproanno[ourtf,"WhiteAd.Pro.ID"],])
sextimemeta$TF.L2FC <- round(sextimemeta$TF.L2FC/0.1)*0.1
sextimemeta$TF.L2FC[sextimemeta$TF.L2FC > 1] <- 1
sextimemeta$TF.L2FC[sextimemeta$TF.L2FC < -1] <- -1
sextimemeta$TF.L2FC <- factor(sextimemeta$TF.L2FC,levels = c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8,-0.9,-1))

ourtftargetpromproxcormeta <- data.frame(row.names = intersect(ourtftargetpromproxgenes,whiternasig),
                                         "Correlation" = rep(0,length(intersect(ourtftargetpromproxgenes,whiternasig))))
ourtftargetexoncormeta <- data.frame(row.names = intersect(ourtftargetexongenes,whiternasig),
                                     "Correlation" = rep(0,length(intersect(ourtftargetexongenes,whiternasig))))
ourtftargetintroncormeta <- data.frame(row.names = intersect(ourtftargetintrongenes,whiternasig),
                                       "Correlation" = rep(0,length(intersect(ourtftargetintrongenes,whiternasig))))
for(i in 1:dim(ourtftargetpromproxcormeta)[1]){
  ourtftargetpromproxcormeta[i,"Correlation"] <- cor(whiteprol2fc[tfproanno[ourtf,"WhiteAd.Pro.ID"],],whitel2fcmat[intersect(ourtftargetpromproxgenes,whiternasig)[i],])
}
ourtftargetpromproxcormeta$Correlation <- round(ourtftargetpromproxcormeta$Correlation/0.2)*0.2
ourtftargetpromproxcormeta$Correlation <- factor(ourtftargetpromproxcormeta$Correlation,levels = c(1,0.8,0.6,0.4,0.2,0,-0.2,-0.4,-0.6,-0.8,-1))

png(file = "Supplemental Figure S17D.png",width = 7,height = 4,units = "in",res = 600)
pheatmap(whitel2fcmat[intersect(ourtftargetpromproxgenes,whiternasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,whiternasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F)
dev.off()

pdf(file = "Supplemental Figure S17D.pdf",width = 7,height = 4)
pheatmap(whitel2fcmat[intersect(ourtftargetpromproxgenes,whiternasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,whiternasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F)
dev.off()


ourtf <- "Mef2c(MADS)/GM12878-Mef2c-ChIP-Seq(GSE32465)/Homer"
ourtftargetpeaks <- gastro50peakmotifs[gastro50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])

sextimemeta <- data.frame(row.names = colnames(gastrol2fcmat),
                          "Sex" = c(rep("Female",4),rep("Male",4)),
                          "Week" = rep(c("1w","2w","4w","8w"),2),
                          "TF.L2FC" = tfprophsigmat[rownames(tfprophsigmat)[grep(tfproanno[ourtf,"Gastro.Pro.ID"],rownames(tfprophsigmat))],])
sextimemeta$TF.L2FC <- round(sextimemeta$TF.L2FC/0.1)*0.1
sextimemeta$TF.L2FC[sextimemeta$TF.L2FC > 1] <- 1
sextimemeta$TF.L2FC[sextimemeta$TF.L2FC < -1] <- -1
sextimemeta$TF.L2FC <- factor(sextimemeta$TF.L2FC,levels = c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8,-0.9,-1))

ourtftargetpromproxcormeta <- data.frame(row.names = intersect(ourtftargetpromproxgenes,gastrornasig),
                                         "Correlation" = rep(0,length(intersect(ourtftargetpromproxgenes,gastrornasig))))
ourtftargetexoncormeta <- data.frame(row.names = intersect(ourtftargetexongenes,gastrornasig),
                                     "Correlation" = rep(0,length(intersect(ourtftargetexongenes,gastrornasig))))
ourtftargetintroncormeta <- data.frame(row.names = intersect(ourtftargetintrongenes,gastrornasig),
                                       "Correlation" = rep(0,length(intersect(ourtftargetintrongenes,gastrornasig))))
for(i in 1:dim(ourtftargetpromproxcormeta)[1]){
  ourtftargetpromproxcormeta[i,"Correlation"] <- cor(tfprophsigmat[rownames(tfprophsigmat)[grep(tfproanno[ourtf,"Gastro.Pro.ID"],rownames(tfprophsigmat))],],gastrol2fcmat[intersect(ourtftargetpromproxgenes,gastrornasig)[i],])
}
ourtftargetpromproxcormeta$Correlation <- round(ourtftargetpromproxcormeta$Correlation/0.2)*0.2
ourtftargetpromproxcormeta$Correlation <- factor(ourtftargetpromproxcormeta$Correlation,levels = c(1,0.8,0.6,0.4,0.2,0,-0.2,-0.4,-0.6,-0.8,-1))

png(file = "Supplemental Figure S17E.png",width = 7,height = 4,units = "in",res = 600)
pheatmap(gastrol2fcmat[intersect(ourtftargetpromproxgenes,gastrornasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,gastrornasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F)
dev.off()

pdf(file = "Supplemental Figure S17E.pdf",width = 7,height = 4)
pheatmap(gastrol2fcmat[intersect(ourtftargetpromproxgenes,gastrornasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,gastrornasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F)
dev.off()


ourtf <- "Mef2a(MADS)/HL1-Mef2a.biotin-ChIP-Seq(GSE21529)/Homer"
ourtftargetpeaks <- heart50peakmotifs[heart50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])

sextimemeta <- data.frame(row.names = colnames(heartl2fcmat),
                          "Sex" = c(rep("Female",4),rep("Male",4)),
                          "Week" = rep(c("1w","2w","4w","8w"),2),
                          "TF.L2FC" = tfprophsigmat[rownames(tfprophsigmat)[grep(tfproanno[ourtf,"Heart.Pro.ID"],rownames(tfprophsigmat))][1],])
sextimemeta$TF.L2FC <- round(sextimemeta$TF.L2FC/0.1)*0.1
sextimemeta$TF.L2FC[sextimemeta$TF.L2FC > 1] <- 1
sextimemeta$TF.L2FC[sextimemeta$TF.L2FC < -1] <- -1
sextimemeta$TF.L2FC <- factor(sextimemeta$TF.L2FC,levels = c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8,-0.9,-1))
ourtftargetpromproxcormeta <- data.frame(row.names = intersect(ourtftargetpromproxgenes,heartrnasig),
                                         "Correlation" = rep(0,length(intersect(ourtftargetpromproxgenes,heartrnasig))))
ourtftargetexoncormeta <- data.frame(row.names = intersect(ourtftargetexongenes,heartrnasig),
                                     "Correlation" = rep(0,length(intersect(ourtftargetexongenes,heartrnasig))))
ourtftargetintroncormeta <- data.frame(row.names = intersect(ourtftargetintrongenes,heartrnasig),
                                       "Correlation" = rep(0,length(intersect(ourtftargetintrongenes,heartrnasig))))
for(i in 1:dim(ourtftargetpromproxcormeta)[1]){
  ourtftargetpromproxcormeta[i,"Correlation"] <- cor(tfprophsigmat[rownames(tfprophsigmat)[grep(tfproanno[ourtf,"Heart.Pro.ID"],rownames(tfprophsigmat))][1],],heartl2fcmat[intersect(ourtftargetpromproxgenes,heartrnasig)[i],])
}
ourtftargetpromproxcormeta$Correlation <- round(ourtftargetpromproxcormeta$Correlation/0.2)*0.2
ourtftargetpromproxcormeta$Correlation <- factor(ourtftargetpromproxcormeta$Correlation,levels = c(1,0.8,0.6,0.4,0.2,0,-0.2,-0.4,-0.6,-0.8,-1))

png(file = "Supplemental Figure S17F.png",width = 7,height = 4,units = "in",res = 600)
pheatmap(heartl2fcmat[intersect(ourtftargetpromproxgenes,heartrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,heartrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F)
dev.off()

pdf(file = "Supplemental Figure S17F.pdf",width = 7,height = 4)
pheatmap(heartl2fcmat[intersect(ourtftargetpromproxgenes,heartrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,heartrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F)
dev.off()


ourtf <- "Atf7(bZIP)/3T3L1-Atf7-ChIP-Seq(GSE56872)/Homer"
ourtftargetpeaks <- kidney50peakmotifs[kidney50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])

sextimemeta <- data.frame(row.names = colnames(kidneyl2fcmat),
                          "Sex" = c(rep("Female",4),rep("Male",4)),
                          "Week" = rep(c("1w","2w","4w","8w"),2),
                          "TF.L2FC" = tfprophsigmat[rownames(tfprophsigmat)[grep(tfproanno[ourtf,"Kidney.Pro.ID"],rownames(tfprophsigmat))],])
sextimemeta$TF.L2FC <- round(sextimemeta$TF.L2FC/0.1)*0.1
sextimemeta$TF.L2FC[sextimemeta$TF.L2FC > 1] <- 1
sextimemeta$TF.L2FC[sextimemeta$TF.L2FC < -1] <- -1
sextimemeta$TF.L2FC <- factor(sextimemeta$TF.L2FC,levels = c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8,-0.9,-1))
ourtftargetpromproxcormeta <- data.frame(row.names = intersect(ourtftargetpromproxgenes,kidneyrnasig),
                                         "Correlation" = rep(0,length(intersect(ourtftargetpromproxgenes,kidneyrnasig))))
ourtftargetexoncormeta <- data.frame(row.names = intersect(ourtftargetexongenes,kidneyrnasig),
                                     "Correlation" = rep(0,length(intersect(ourtftargetexongenes,kidneyrnasig))))
ourtftargetintroncormeta <- data.frame(row.names = intersect(ourtftargetintrongenes,kidneyrnasig),
                                       "Correlation" = rep(0,length(intersect(ourtftargetintrongenes,kidneyrnasig))))
for(i in 1:dim(ourtftargetpromproxcormeta)[1]){
  ourtftargetpromproxcormeta[i,"Correlation"] <- cor(tfprophsigmat[rownames(tfprophsigmat)[grep(tfproanno[ourtf,"Kidney.Pro.ID"],rownames(tfprophsigmat))],],kidneyl2fcmat[intersect(ourtftargetpromproxgenes,kidneyrnasig)[i],])
}
ourtftargetpromproxcormeta$Correlation <- round(ourtftargetpromproxcormeta$Correlation/0.2)*0.2
ourtftargetpromproxcormeta$Correlation <- factor(ourtftargetpromproxcormeta$Correlation,levels = c(1,0.8,0.6,0.4,0.2,0,-0.2,-0.4,-0.6,-0.8,-1))

png(file = "Supplemental Figure S17G.png",width = 7,height = 4,units = "in",res = 600)
pheatmap(kidneyl2fcmat[intersect(ourtftargetpromproxgenes,kidneyrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,kidneyrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F)
dev.off()

pdf(file = "Supplemental Figure S17G.pdf",width = 7,height = 4)
pheatmap(kidneyl2fcmat[intersect(ourtftargetpromproxgenes,kidneyrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,kidneyrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F)
dev.off()


ourtf <- "STAT1(Stat)/HelaS3-STAT1-ChIP-Seq(GSE12782)/Homer"
ourtftargetpeaks <- liver50peakmotifs[liver50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])

sextimemeta <- data.frame(row.names = colnames(liverl2fcmat),
                          "Sex" = c(rep("Female",4),rep("Male",4)),
                          "Week" = rep(c("1w","2w","4w","8w"),2),
                          "TF.L2FC" = tfprophsigmat[rownames(tfprophsigmat)[grep(tfproanno[ourtf,"Liver.Pro.ID"],rownames(tfprophsigmat))],])
sextimemeta$TF.L2FC <- round(sextimemeta$TF.L2FC/0.1)*0.1
sextimemeta$TF.L2FC[sextimemeta$TF.L2FC > 1] <- 1
sextimemeta$TF.L2FC[sextimemeta$TF.L2FC < -1] <- -1
sextimemeta$TF.L2FC <- factor(sextimemeta$TF.L2FC,levels = c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8,-0.9,-1))
ourtftargetpromproxcormeta <- data.frame(row.names = intersect(ourtftargetpromproxgenes,liverrnasig),
                                         "Correlation" = rep(0,length(intersect(ourtftargetpromproxgenes,liverrnasig))))
ourtftargetexoncormeta <- data.frame(row.names = intersect(ourtftargetexongenes,liverrnasig),
                                     "Correlation" = rep(0,length(intersect(ourtftargetexongenes,liverrnasig))))
ourtftargetintroncormeta <- data.frame(row.names = intersect(ourtftargetintrongenes,liverrnasig),
                                       "Correlation" = rep(0,length(intersect(ourtftargetintrongenes,liverrnasig))))
for(i in 1:dim(ourtftargetpromproxcormeta)[1]){
  ourtftargetpromproxcormeta[i,"Correlation"] <- cor(tfprophsigmat[rownames(tfprophsigmat)[grep(tfproanno[ourtf,"Liver.Pro.ID"],rownames(tfprophsigmat))],],liverl2fcmat[intersect(ourtftargetpromproxgenes,liverrnasig)[i],])
}
ourtftargetpromproxcormeta$Correlation <- round(ourtftargetpromproxcormeta$Correlation/0.2)*0.2
ourtftargetpromproxcormeta$Correlation <- factor(ourtftargetpromproxcormeta$Correlation,levels = c(1,0.8,0.6,0.4,0.2,0,-0.2,-0.4,-0.6,-0.8,-1))

png(file = "Supplemental Figure S17H.png",width = 7,height = 4,units = "in",res = 600)
pheatmap(liverl2fcmat[intersect(ourtftargetpromproxgenes,liverrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,liverrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F)
dev.off()

pdf(file = "Supplemental Figure S17H.pdf",width = 7,height = 4)
pheatmap(liverl2fcmat[intersect(ourtftargetpromproxgenes,liverrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,liverrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F)
dev.off()


# We need to investigate what are the compelling relationships between changing TFs at the
# RNA level and their enrichment among DEGs

oursigtfrnalist <- unique(c(rownames(tfanno)[tfanno$Ensembl %in% tfgastrornasig],
                            rownames(tfanno)[tfanno$Ensembl %in% tfheartrnasig],
                            rownames(tfanno)[tfanno$Ensembl %in% tfhippornasig],
                            rownames(tfanno)[tfanno$Ensembl %in% tfkidneyrnasig],
                            rownames(tfanno)[tfanno$Ensembl %in% tfliverrnasig],
                            rownames(tfanno)[tfanno$Ensembl %in% tflungrnasig],
                            rownames(tfanno)[tfanno$Ensembl %in% tfbrownrnasig],
                            rownames(tfanno)[tfanno$Ensembl %in% tfwhiternasig]))

# SKM-GN 

gastrotfrnatargetpromproxgenesigpct <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfgastrornasig]),ncol = 1)
rownames(gastrotfrnatargetpromproxgenesigpct) <- rownames(tfanno)[tfanno$Ensembl %in% tfgastrornasig]
colnames(gastrotfrnatargetpromproxgenesigpct) <- "Percent.Significant"

gastrotfrnatargetintrongenesigpct <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfgastrornasig]),ncol = 1)
rownames(gastrotfrnatargetintrongenesigpct) <- rownames(tfanno)[tfanno$Ensembl %in% tfgastrornasig]
colnames(gastrotfrnatargetintrongenesigpct) <- "Percent.Significant"

gastrotfrnatargetexongenesigpct <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfgastrornasig]),ncol = 1)
rownames(gastrotfrnatargetexongenesigpct) <- rownames(tfanno)[tfanno$Ensembl %in% tfgastrornasig]
colnames(gastrotfrnatargetexongenesigpct) <- "Percent.Significant"

gastrotfrnatargetintergenicgenesigpct <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfgastrornasig]),ncol = 1)
rownames(gastrotfrnatargetintergenicgenesigpct) <- rownames(tfanno)[tfanno$Ensembl %in% tfgastrornasig]
colnames(gastrotfrnatargetintergenicgenesigpct) <- "Percent.Significant"


gastrotfrnatargetpromproxgenesigcount <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfgastrornasig]),ncol = 1)
rownames(gastrotfrnatargetpromproxgenesigcount) <- rownames(tfanno)[tfanno$Ensembl %in% tfgastrornasig]
colnames(gastrotfrnatargetpromproxgenesigcount) <- "Count.Significant"

gastrotfrnatargetintrongenesigcount <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfgastrornasig]),ncol = 1)
rownames(gastrotfrnatargetintrongenesigcount) <- rownames(tfanno)[tfanno$Ensembl %in% tfgastrornasig]
colnames(gastrotfrnatargetintrongenesigcount) <- "Count.Significant"

gastrotfrnatargetexongenesigcount <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfgastrornasig]),ncol = 1)
rownames(gastrotfrnatargetexongenesigcount) <- rownames(tfanno)[tfanno$Ensembl %in% tfgastrornasig]
colnames(gastrotfrnatargetexongenesigcount) <- "Count.Significant"

gastrotfrnatargetintergenicgenesigcount <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfgastrornasig]),ncol = 1)
rownames(gastrotfrnatargetintergenicgenesigcount) <- rownames(tfanno)[tfanno$Ensembl %in% tfgastrornasig]
colnames(gastrotfrnatargetintergenicgenesigcount) <- "Count.Significant"


for(i in 1:length(rownames(tfanno)[tfanno$Ensembl %in% tfgastrornasig])){
  ourtf <- rownames(tfanno)[tfanno$Ensembl %in% tfgastrornasig][i]
  ourtftargetpeaks <- gastro50peakmotifs[gastro50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
  ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
  ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
  ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
  ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
  ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])
  gastrotfrnatargetpromproxgenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetpromproxgenes,gastrornasig))
  gastrotfrnatargetintrongenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetintrongenes,gastrornasig))
  gastrotfrnatargetexongenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetexongenes,gastrornasig))
  gastrotfrnatargetintergenicgenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetintergenicgenes,gastrornasig))
  gastrotfrnatargetpromproxgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetpromproxgenes,gastrornasig))/length(ourtftargetpromproxgenes)
  gastrotfrnatargetintrongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintrongenes,gastrornasig))/length(ourtftargetintrongenes)
  gastrotfrnatargetexongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetexongenes,gastrornasig))/length(ourtftargetexongenes)
  gastrotfrnatargetintergenicgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintergenicgenes,gastrornasig))/length(ourtftargetintergenicgenes)
}

gastrornatargetsigdf <- data.frame("Percent.Significant" = c(gastrotfrnatargetpromproxgenesigpct[gastrotfrnatargetpromproxgenesigcount > 3],
                                                             gastrotfrnatargetintrongenesigpct[gastrotfrnatargetpromproxgenesigcount > 3],
                                                             gastrotfrnatargetexongenesigpct[gastrotfrnatargetpromproxgenesigcount > 3],
                                                             gastrotfrnatargetintergenicgenesigpct[gastrotfrnatargetpromproxgenesigcount > 3]),
                                   "Region" = c(rep("Promoter (<=1kb)",length(gastrotfrnatargetpromproxgenesigpct[gastrotfrnatargetpromproxgenesigcount > 3])),
                                                rep("Intron",length(gastrotfrnatargetpromproxgenesigpct[gastrotfrnatargetpromproxgenesigcount > 3])),
                                                rep("Exon",length(gastrotfrnatargetpromproxgenesigpct[gastrotfrnatargetpromproxgenesigcount > 3])),
                                                rep("Distal Intergenic",length(gastrotfrnatargetpromproxgenesigpct[gastrotfrnatargetpromproxgenesigcount > 3]))),
                                   "Transcription Factor" = rep(gsub("\\(.*","",rownames(gastrotfrnatargetpromproxgenesigpct)[gastrotfrnatargetpromproxgenesigcount > 3]),4))

# HEART

hearttfrnatargetpromproxgenesigpct <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfheartrnasig]),ncol = 1)
rownames(hearttfrnatargetpromproxgenesigpct) <- rownames(tfanno)[tfanno$Ensembl %in% tfheartrnasig]
colnames(hearttfrnatargetpromproxgenesigpct) <- "Percent.Significant"

hearttfrnatargetintrongenesigpct <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfheartrnasig]),ncol = 1)
rownames(hearttfrnatargetintrongenesigpct) <- rownames(tfanno)[tfanno$Ensembl %in% tfheartrnasig]
colnames(hearttfrnatargetintrongenesigpct) <- "Percent.Significant"

hearttfrnatargetexongenesigpct <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfheartrnasig]),ncol = 1)
rownames(hearttfrnatargetexongenesigpct) <- rownames(tfanno)[tfanno$Ensembl %in% tfheartrnasig]
colnames(hearttfrnatargetexongenesigpct) <- "Percent.Significant"

hearttfrnatargetintergenicgenesigpct <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfheartrnasig]),ncol = 1)
rownames(hearttfrnatargetintergenicgenesigpct) <- rownames(tfanno)[tfanno$Ensembl %in% tfheartrnasig]
colnames(hearttfrnatargetintergenicgenesigpct) <- "Percent.Significant"


hearttfrnatargetpromproxgenesigcount <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfheartrnasig]),ncol = 1)
rownames(hearttfrnatargetpromproxgenesigcount) <- rownames(tfanno)[tfanno$Ensembl %in% tfheartrnasig]
colnames(hearttfrnatargetpromproxgenesigcount) <- "Count.Significant"

hearttfrnatargetintrongenesigcount <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfheartrnasig]),ncol = 1)
rownames(hearttfrnatargetintrongenesigcount) <- rownames(tfanno)[tfanno$Ensembl %in% tfheartrnasig]
colnames(hearttfrnatargetintrongenesigcount) <- "Count.Significant"

hearttfrnatargetexongenesigcount <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfheartrnasig]),ncol = 1)
rownames(hearttfrnatargetexongenesigcount) <- rownames(tfanno)[tfanno$Ensembl %in% tfheartrnasig]
colnames(hearttfrnatargetexongenesigcount) <- "Count.Significant"

hearttfrnatargetintergenicgenesigcount <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfheartrnasig]),ncol = 1)
rownames(hearttfrnatargetintergenicgenesigcount) <- rownames(tfanno)[tfanno$Ensembl %in% tfheartrnasig]
colnames(hearttfrnatargetintergenicgenesigcount) <- "Count.Significant"


for(i in 1:length(rownames(tfanno)[tfanno$Ensembl %in% tfheartrnasig])){
  ourtf <- rownames(tfanno)[tfanno$Ensembl %in% tfheartrnasig][i]
  ourtftargetpeaks <- heart50peakmotifs[heart50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
  ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
  ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
  ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
  ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
  ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])
  hearttfrnatargetpromproxgenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetpromproxgenes,heartrnasig))
  hearttfrnatargetintrongenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetintrongenes,heartrnasig))
  hearttfrnatargetexongenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetexongenes,heartrnasig))
  hearttfrnatargetintergenicgenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetintergenicgenes,heartrnasig))
  hearttfrnatargetpromproxgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetpromproxgenes,heartrnasig))/length(ourtftargetpromproxgenes)
  hearttfrnatargetintrongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintrongenes,heartrnasig))/length(ourtftargetintrongenes)
  hearttfrnatargetexongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetexongenes,heartrnasig))/length(ourtftargetexongenes)
  hearttfrnatargetintergenicgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintergenicgenes,heartrnasig))/length(ourtftargetintergenicgenes)
}

heartrnatargetsigdf <- data.frame("Percent.Significant" = c(hearttfrnatargetpromproxgenesigpct[hearttfrnatargetpromproxgenesigcount > 3],
                                                            hearttfrnatargetintrongenesigpct[hearttfrnatargetpromproxgenesigcount > 3],
                                                            hearttfrnatargetexongenesigpct[hearttfrnatargetpromproxgenesigcount > 3],
                                                            hearttfrnatargetintergenicgenesigpct[hearttfrnatargetpromproxgenesigcount > 3]),
                                  "Region" = c(rep("Promoter (<=1kb)",length(hearttfrnatargetpromproxgenesigpct[hearttfrnatargetpromproxgenesigcount > 3])),
                                               rep("Intron",length(hearttfrnatargetpromproxgenesigpct[hearttfrnatargetpromproxgenesigcount > 3])),
                                               rep("Exon",length(hearttfrnatargetpromproxgenesigpct[hearttfrnatargetpromproxgenesigcount > 3])),
                                               rep("Distal Intergenic",length(hearttfrnatargetpromproxgenesigpct[hearttfrnatargetpromproxgenesigcount > 3]))),
                                  "Transcription Factor" = rep(gsub("\\(.*","",rownames(hearttfrnatargetpromproxgenesigpct)[hearttfrnatargetpromproxgenesigcount > 3]),4))



# HIPPOC

hippotfrnatargetpromproxgenesigpct <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfhippornasig]),ncol = 1)
rownames(hippotfrnatargetpromproxgenesigpct) <- rownames(tfanno)[tfanno$Ensembl %in% tfhippornasig]
colnames(hippotfrnatargetpromproxgenesigpct) <- "Percent.Significant"

hippotfrnatargetintrongenesigpct <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfhippornasig]),ncol = 1)
rownames(hippotfrnatargetintrongenesigpct) <- rownames(tfanno)[tfanno$Ensembl %in% tfhippornasig]
colnames(hippotfrnatargetintrongenesigpct) <- "Percent.Significant"

hippotfrnatargetexongenesigpct <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfhippornasig]),ncol = 1)
rownames(hippotfrnatargetexongenesigpct) <- rownames(tfanno)[tfanno$Ensembl %in% tfhippornasig]
colnames(hippotfrnatargetexongenesigpct) <- "Percent.Significant"

hippotfrnatargetintergenicgenesigpct <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfhippornasig]),ncol = 1)
rownames(hippotfrnatargetintergenicgenesigpct) <- rownames(tfanno)[tfanno$Ensembl %in% tfhippornasig]
colnames(hippotfrnatargetintergenicgenesigpct) <- "Percent.Significant"


hippotfrnatargetpromproxgenesigcount <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfhippornasig]),ncol = 1)
rownames(hippotfrnatargetpromproxgenesigcount) <- rownames(tfanno)[tfanno$Ensembl %in% tfhippornasig]
colnames(hippotfrnatargetpromproxgenesigcount) <- "Count.Significant"

hippotfrnatargetintrongenesigcount <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfhippornasig]),ncol = 1)
rownames(hippotfrnatargetintrongenesigcount) <- rownames(tfanno)[tfanno$Ensembl %in% tfhippornasig]
colnames(hippotfrnatargetintrongenesigcount) <- "Count.Significant"

hippotfrnatargetexongenesigcount <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfhippornasig]),ncol = 1)
rownames(hippotfrnatargetexongenesigcount) <- rownames(tfanno)[tfanno$Ensembl %in% tfhippornasig]
colnames(hippotfrnatargetexongenesigcount) <- "Count.Significant"

hippotfrnatargetintergenicgenesigcount <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfhippornasig]),ncol = 1)
rownames(hippotfrnatargetintergenicgenesigcount) <- rownames(tfanno)[tfanno$Ensembl %in% tfhippornasig]
colnames(hippotfrnatargetintergenicgenesigcount) <- "Count.Significant"


for(i in 1:length(rownames(tfanno)[tfanno$Ensembl %in% tfhippornasig])){
  ourtf <- rownames(tfanno)[tfanno$Ensembl %in% tfhippornasig][i]
  ourtftargetpeaks <- hippo50peakmotifs[hippo50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
  ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
  ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
  ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
  ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
  ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])
  hippotfrnatargetpromproxgenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetpromproxgenes,hippornasig))
  hippotfrnatargetintrongenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetintrongenes,hippornasig))
  hippotfrnatargetexongenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetexongenes,hippornasig))
  hippotfrnatargetintergenicgenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetintergenicgenes,hippornasig))
  hippotfrnatargetpromproxgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetpromproxgenes,hippornasig))/length(ourtftargetpromproxgenes)
  hippotfrnatargetintrongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintrongenes,hippornasig))/length(ourtftargetintrongenes)
  hippotfrnatargetexongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetexongenes,hippornasig))/length(ourtftargetexongenes)
  hippotfrnatargetintergenicgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintergenicgenes,hippornasig))/length(ourtftargetintergenicgenes)
}

hippornatargetsigdf <- data.frame("Percent.Significant" = c(hippotfrnatargetpromproxgenesigpct[hippotfrnatargetpromproxgenesigcount > 3],
                                                            hippotfrnatargetintrongenesigpct[hippotfrnatargetpromproxgenesigcount > 3],
                                                            hippotfrnatargetexongenesigpct[hippotfrnatargetpromproxgenesigcount > 3],
                                                            hippotfrnatargetintergenicgenesigpct[hippotfrnatargetpromproxgenesigcount > 3]),
                                  "Region" = c(rep("Promoter (<=1kb)",length(hippotfrnatargetpromproxgenesigpct[hippotfrnatargetpromproxgenesigcount > 3])),
                                               rep("Intron",length(hippotfrnatargetpromproxgenesigpct[hippotfrnatargetpromproxgenesigcount > 3])),
                                               rep("Exon",length(hippotfrnatargetpromproxgenesigpct[hippotfrnatargetpromproxgenesigcount > 3])),
                                               rep("Distal Intergenic",length(hippotfrnatargetpromproxgenesigpct[hippotfrnatargetpromproxgenesigcount > 3]))),
                                  "Transcription Factor" = rep(gsub("\\(.*","",rownames(hippotfrnatargetpromproxgenesigpct)[hippotfrnatargetpromproxgenesigcount > 3]),4))


# KIDNEY

kidneytfrnatargetpromproxgenesigpct <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfkidneyrnasig]),ncol = 1)
rownames(kidneytfrnatargetpromproxgenesigpct) <- rownames(tfanno)[tfanno$Ensembl %in% tfkidneyrnasig]
colnames(kidneytfrnatargetpromproxgenesigpct) <- "Percent.Significant"

kidneytfrnatargetintrongenesigpct <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfkidneyrnasig]),ncol = 1)
rownames(kidneytfrnatargetintrongenesigpct) <- rownames(tfanno)[tfanno$Ensembl %in% tfkidneyrnasig]
colnames(kidneytfrnatargetintrongenesigpct) <- "Percent.Significant"

kidneytfrnatargetexongenesigpct <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfkidneyrnasig]),ncol = 1)
rownames(kidneytfrnatargetexongenesigpct) <- rownames(tfanno)[tfanno$Ensembl %in% tfkidneyrnasig]
colnames(kidneytfrnatargetexongenesigpct) <- "Percent.Significant"

kidneytfrnatargetintergenicgenesigpct <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfkidneyrnasig]),ncol = 1)
rownames(kidneytfrnatargetintergenicgenesigpct) <- rownames(tfanno)[tfanno$Ensembl %in% tfkidneyrnasig]
colnames(kidneytfrnatargetintergenicgenesigpct) <- "Percent.Significant"


kidneytfrnatargetpromproxgenesigcount <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfkidneyrnasig]),ncol = 1)
rownames(kidneytfrnatargetpromproxgenesigcount) <- rownames(tfanno)[tfanno$Ensembl %in% tfkidneyrnasig]
colnames(kidneytfrnatargetpromproxgenesigcount) <- "Count.Significant"

kidneytfrnatargetintrongenesigcount <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfkidneyrnasig]),ncol = 1)
rownames(kidneytfrnatargetintrongenesigcount) <- rownames(tfanno)[tfanno$Ensembl %in% tfkidneyrnasig]
colnames(kidneytfrnatargetintrongenesigcount) <- "Count.Significant"

kidneytfrnatargetexongenesigcount <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfkidneyrnasig]),ncol = 1)
rownames(kidneytfrnatargetexongenesigcount) <- rownames(tfanno)[tfanno$Ensembl %in% tfkidneyrnasig]
colnames(kidneytfrnatargetexongenesigcount) <- "Count.Significant"

kidneytfrnatargetintergenicgenesigcount <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfkidneyrnasig]),ncol = 1)
rownames(kidneytfrnatargetintergenicgenesigcount) <- rownames(tfanno)[tfanno$Ensembl %in% tfkidneyrnasig]
colnames(kidneytfrnatargetintergenicgenesigcount) <- "Count.Significant"


for(i in 1:length(rownames(tfanno)[tfanno$Ensembl %in% tfkidneyrnasig])){
  ourtf <- rownames(tfanno)[tfanno$Ensembl %in% tfkidneyrnasig][i]
  ourtftargetpeaks <- kidney50peakmotifs[kidney50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
  ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
  ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
  ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
  ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
  ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])
  kidneytfrnatargetpromproxgenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetpromproxgenes,kidneyrnasig))
  kidneytfrnatargetintrongenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetintrongenes,kidneyrnasig))
  kidneytfrnatargetexongenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetexongenes,kidneyrnasig))
  kidneytfrnatargetintergenicgenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetintergenicgenes,kidneyrnasig))
  kidneytfrnatargetpromproxgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetpromproxgenes,kidneyrnasig))/length(ourtftargetpromproxgenes)
  kidneytfrnatargetintrongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintrongenes,kidneyrnasig))/length(ourtftargetintrongenes)
  kidneytfrnatargetexongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetexongenes,kidneyrnasig))/length(ourtftargetexongenes)
  kidneytfrnatargetintergenicgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintergenicgenes,kidneyrnasig))/length(ourtftargetintergenicgenes)
}

kidneyrnatargetsigdf <- data.frame("Percent.Significant" = c(kidneytfrnatargetpromproxgenesigpct[kidneytfrnatargetpromproxgenesigcount > 3],
                                                             kidneytfrnatargetintrongenesigpct[kidneytfrnatargetpromproxgenesigcount > 3],
                                                             kidneytfrnatargetexongenesigpct[kidneytfrnatargetpromproxgenesigcount > 3],
                                                             kidneytfrnatargetintergenicgenesigpct[kidneytfrnatargetpromproxgenesigcount > 3]),
                                   "Region" = c(rep("Promoter (<=1kb)",length(kidneytfrnatargetpromproxgenesigpct[kidneytfrnatargetpromproxgenesigcount > 3])),
                                                rep("Intron",length(kidneytfrnatargetpromproxgenesigpct[kidneytfrnatargetpromproxgenesigcount > 3])),
                                                rep("Exon",length(kidneytfrnatargetpromproxgenesigpct[kidneytfrnatargetpromproxgenesigcount > 3])),
                                                rep("Distal Intergenic",length(kidneytfrnatargetpromproxgenesigpct[kidneytfrnatargetpromproxgenesigcount > 3]))),
                                   "Transcription Factor" = rep(gsub("\\(.*","",rownames(kidneytfrnatargetpromproxgenesigpct)[kidneytfrnatargetpromproxgenesigcount > 3]),4))



# LIVER

livertfrnatargetpromproxgenesigpct <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfliverrnasig]),ncol = 1)
rownames(livertfrnatargetpromproxgenesigpct) <- rownames(tfanno)[tfanno$Ensembl %in% tfliverrnasig]
colnames(livertfrnatargetpromproxgenesigpct) <- "Percent.Significant"

livertfrnatargetintrongenesigpct <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfliverrnasig]),ncol = 1)
rownames(livertfrnatargetintrongenesigpct) <- rownames(tfanno)[tfanno$Ensembl %in% tfliverrnasig]
colnames(livertfrnatargetintrongenesigpct) <- "Percent.Significant"

livertfrnatargetexongenesigpct <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfliverrnasig]),ncol = 1)
rownames(livertfrnatargetexongenesigpct) <- rownames(tfanno)[tfanno$Ensembl %in% tfliverrnasig]
colnames(livertfrnatargetexongenesigpct) <- "Percent.Significant"

livertfrnatargetintergenicgenesigpct <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfliverrnasig]),ncol = 1)
rownames(livertfrnatargetintergenicgenesigpct) <- rownames(tfanno)[tfanno$Ensembl %in% tfliverrnasig]
colnames(livertfrnatargetintergenicgenesigpct) <- "Percent.Significant"


livertfrnatargetpromproxgenesigcount <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfliverrnasig]),ncol = 1)
rownames(livertfrnatargetpromproxgenesigcount) <- rownames(tfanno)[tfanno$Ensembl %in% tfliverrnasig]
colnames(livertfrnatargetpromproxgenesigcount) <- "Count.Significant"

livertfrnatargetintrongenesigcount <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfliverrnasig]),ncol = 1)
rownames(livertfrnatargetintrongenesigcount) <- rownames(tfanno)[tfanno$Ensembl %in% tfliverrnasig]
colnames(livertfrnatargetintrongenesigcount) <- "Count.Significant"

livertfrnatargetexongenesigcount <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfliverrnasig]),ncol = 1)
rownames(livertfrnatargetexongenesigcount) <- rownames(tfanno)[tfanno$Ensembl %in% tfliverrnasig]
colnames(livertfrnatargetexongenesigcount) <- "Count.Significant"

livertfrnatargetintergenicgenesigcount <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfliverrnasig]),ncol = 1)
rownames(livertfrnatargetintergenicgenesigcount) <- rownames(tfanno)[tfanno$Ensembl %in% tfliverrnasig]
colnames(livertfrnatargetintergenicgenesigcount) <- "Count.Significant"


for(i in 1:length(rownames(tfanno)[tfanno$Ensembl %in% tfliverrnasig])){
  ourtf <- rownames(tfanno)[tfanno$Ensembl %in% tfliverrnasig][i]
  ourtftargetpeaks <- liver50peakmotifs[liver50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
  ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
  ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
  ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
  ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
  ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])
  livertfrnatargetpromproxgenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetpromproxgenes,liverrnasig))
  livertfrnatargetintrongenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetintrongenes,liverrnasig))
  livertfrnatargetexongenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetexongenes,liverrnasig))
  livertfrnatargetintergenicgenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetintergenicgenes,liverrnasig))
  livertfrnatargetpromproxgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetpromproxgenes,liverrnasig))/length(ourtftargetpromproxgenes)
  livertfrnatargetintrongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintrongenes,liverrnasig))/length(ourtftargetintrongenes)
  livertfrnatargetexongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetexongenes,liverrnasig))/length(ourtftargetexongenes)
  livertfrnatargetintergenicgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintergenicgenes,liverrnasig))/length(ourtftargetintergenicgenes)
}

liverrnatargetsigdf <- data.frame("Percent.Significant" = c(livertfrnatargetpromproxgenesigpct[livertfrnatargetpromproxgenesigcount > 3],
                                                            livertfrnatargetintrongenesigpct[livertfrnatargetpromproxgenesigcount > 3],
                                                            livertfrnatargetexongenesigpct[livertfrnatargetpromproxgenesigcount > 3],
                                                            livertfrnatargetintergenicgenesigpct[livertfrnatargetpromproxgenesigcount > 3]),
                                  "Region" = c(rep("Promoter (<=1kb)",length(livertfrnatargetpromproxgenesigpct[livertfrnatargetpromproxgenesigcount > 3])),
                                               rep("Intron",length(livertfrnatargetpromproxgenesigpct[livertfrnatargetpromproxgenesigcount > 3])),
                                               rep("Exon",length(livertfrnatargetpromproxgenesigpct[livertfrnatargetpromproxgenesigcount > 3])),
                                               rep("Distal Intergenic",length(livertfrnatargetpromproxgenesigpct[livertfrnatargetpromproxgenesigcount > 3]))),
                                  "Transcription Factor" = rep(gsub("\\(.*","",rownames(livertfrnatargetpromproxgenesigpct)[livertfrnatargetpromproxgenesigcount > 3]),4))



# LUNG

lungtfrnatargetpromproxgenesigpct <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tflungrnasig]),ncol = 1)
rownames(lungtfrnatargetpromproxgenesigpct) <- rownames(tfanno)[tfanno$Ensembl %in% tflungrnasig]
colnames(lungtfrnatargetpromproxgenesigpct) <- "Percent.Significant"

lungtfrnatargetintrongenesigpct <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tflungrnasig]),ncol = 1)
rownames(lungtfrnatargetintrongenesigpct) <- rownames(tfanno)[tfanno$Ensembl %in% tflungrnasig]
colnames(lungtfrnatargetintrongenesigpct) <- "Percent.Significant"

lungtfrnatargetexongenesigpct <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tflungrnasig]),ncol = 1)
rownames(lungtfrnatargetexongenesigpct) <- rownames(tfanno)[tfanno$Ensembl %in% tflungrnasig]
colnames(lungtfrnatargetexongenesigpct) <- "Percent.Significant"

lungtfrnatargetintergenicgenesigpct <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tflungrnasig]),ncol = 1)
rownames(lungtfrnatargetintergenicgenesigpct) <- rownames(tfanno)[tfanno$Ensembl %in% tflungrnasig]
colnames(lungtfrnatargetintergenicgenesigpct) <- "Percent.Significant"


lungtfrnatargetpromproxgenesigcount <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tflungrnasig]),ncol = 1)
rownames(lungtfrnatargetpromproxgenesigcount) <- rownames(tfanno)[tfanno$Ensembl %in% tflungrnasig]
colnames(lungtfrnatargetpromproxgenesigcount) <- "Count.Significant"

lungtfrnatargetintrongenesigcount <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tflungrnasig]),ncol = 1)
rownames(lungtfrnatargetintrongenesigcount) <- rownames(tfanno)[tfanno$Ensembl %in% tflungrnasig]
colnames(lungtfrnatargetintrongenesigcount) <- "Count.Significant"

lungtfrnatargetexongenesigcount <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tflungrnasig]),ncol = 1)
rownames(lungtfrnatargetexongenesigcount) <- rownames(tfanno)[tfanno$Ensembl %in% tflungrnasig]
colnames(lungtfrnatargetexongenesigcount) <- "Count.Significant"

lungtfrnatargetintergenicgenesigcount <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tflungrnasig]),ncol = 1)
rownames(lungtfrnatargetintergenicgenesigcount) <- rownames(tfanno)[tfanno$Ensembl %in% tflungrnasig]
colnames(lungtfrnatargetintergenicgenesigcount) <- "Count.Significant"


for(i in 1:length(rownames(tfanno)[tfanno$Ensembl %in% tflungrnasig])){
  ourtf <- rownames(tfanno)[tfanno$Ensembl %in% tflungrnasig][i]
  ourtftargetpeaks <- lung50peakmotifs[lung50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
  ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
  ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
  ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
  ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
  ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])
  lungtfrnatargetpromproxgenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetpromproxgenes,lungrnasig))
  lungtfrnatargetintrongenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetintrongenes,lungrnasig))
  lungtfrnatargetexongenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetexongenes,lungrnasig))
  lungtfrnatargetintergenicgenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetintergenicgenes,lungrnasig))
  lungtfrnatargetpromproxgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetpromproxgenes,lungrnasig))/length(ourtftargetpromproxgenes)
  lungtfrnatargetintrongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintrongenes,lungrnasig))/length(ourtftargetintrongenes)
  lungtfrnatargetexongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetexongenes,lungrnasig))/length(ourtftargetexongenes)
  lungtfrnatargetintergenicgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintergenicgenes,lungrnasig))/length(ourtftargetintergenicgenes)
}

lungrnatargetsigdf <- data.frame("Percent.Significant" = c(lungtfrnatargetpromproxgenesigpct[lungtfrnatargetpromproxgenesigcount > 3],
                                                           lungtfrnatargetintrongenesigpct[lungtfrnatargetpromproxgenesigcount > 3],
                                                           lungtfrnatargetexongenesigpct[lungtfrnatargetpromproxgenesigcount > 3],
                                                           lungtfrnatargetintergenicgenesigpct[lungtfrnatargetpromproxgenesigcount > 3]),
                                 "Region" = c(rep("Promoter (<=1kb)",length(lungtfrnatargetpromproxgenesigpct[lungtfrnatargetpromproxgenesigcount > 3])),
                                              rep("Intron",length(lungtfrnatargetpromproxgenesigpct[lungtfrnatargetpromproxgenesigcount > 3])),
                                              rep("Exon",length(lungtfrnatargetpromproxgenesigpct[lungtfrnatargetpromproxgenesigcount > 3])),
                                              rep("Distal Intergenic",length(lungtfrnatargetpromproxgenesigpct[lungtfrnatargetpromproxgenesigcount > 3]))),
                                 "Transcription Factor" = rep(gsub("\\(.*","",rownames(lungtfrnatargetpromproxgenesigpct)[lungtfrnatargetpromproxgenesigcount > 3]),4))



# BAT

browntfrnatargetpromproxgenesigpct <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfbrownrnasig]),ncol = 1)
rownames(browntfrnatargetpromproxgenesigpct) <- rownames(tfanno)[tfanno$Ensembl %in% tfbrownrnasig]
colnames(browntfrnatargetpromproxgenesigpct) <- "Percent.Significant"

browntfrnatargetintrongenesigpct <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfbrownrnasig]),ncol = 1)
rownames(browntfrnatargetintrongenesigpct) <- rownames(tfanno)[tfanno$Ensembl %in% tfbrownrnasig]
colnames(browntfrnatargetintrongenesigpct) <- "Percent.Significant"

browntfrnatargetexongenesigpct <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfbrownrnasig]),ncol = 1)
rownames(browntfrnatargetexongenesigpct) <- rownames(tfanno)[tfanno$Ensembl %in% tfbrownrnasig]
colnames(browntfrnatargetexongenesigpct) <- "Percent.Significant"

browntfrnatargetintergenicgenesigpct <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfbrownrnasig]),ncol = 1)
rownames(browntfrnatargetintergenicgenesigpct) <- rownames(tfanno)[tfanno$Ensembl %in% tfbrownrnasig]
colnames(browntfrnatargetintergenicgenesigpct) <- "Percent.Significant"


browntfrnatargetpromproxgenesigcount <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfbrownrnasig]),ncol = 1)
rownames(browntfrnatargetpromproxgenesigcount) <- rownames(tfanno)[tfanno$Ensembl %in% tfbrownrnasig]
colnames(browntfrnatargetpromproxgenesigcount) <- "Count.Significant"

browntfrnatargetintrongenesigcount <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfbrownrnasig]),ncol = 1)
rownames(browntfrnatargetintrongenesigcount) <- rownames(tfanno)[tfanno$Ensembl %in% tfbrownrnasig]
colnames(browntfrnatargetintrongenesigcount) <- "Count.Significant"

browntfrnatargetexongenesigcount <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfbrownrnasig]),ncol = 1)
rownames(browntfrnatargetexongenesigcount) <- rownames(tfanno)[tfanno$Ensembl %in% tfbrownrnasig]
colnames(browntfrnatargetexongenesigcount) <- "Count.Significant"

browntfrnatargetintergenicgenesigcount <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfbrownrnasig]),ncol = 1)
rownames(browntfrnatargetintergenicgenesigcount) <- rownames(tfanno)[tfanno$Ensembl %in% tfbrownrnasig]
colnames(browntfrnatargetintergenicgenesigcount) <- "Count.Significant"


for(i in 1:length(rownames(tfanno)[tfanno$Ensembl %in% tfbrownrnasig])){
  ourtf <- rownames(tfanno)[tfanno$Ensembl %in% tfbrownrnasig][i]
  ourtftargetpeaks <- brown50peakmotifs[brown50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
  ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
  ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
  ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
  ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
  ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])
  browntfrnatargetpromproxgenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetpromproxgenes,brownrnasig))
  browntfrnatargetintrongenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetintrongenes,brownrnasig))
  browntfrnatargetexongenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetexongenes,brownrnasig))
  browntfrnatargetintergenicgenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetintergenicgenes,brownrnasig))
  browntfrnatargetpromproxgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetpromproxgenes,brownrnasig))/length(ourtftargetpromproxgenes)
  browntfrnatargetintrongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintrongenes,brownrnasig))/length(ourtftargetintrongenes)
  browntfrnatargetexongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetexongenes,brownrnasig))/length(ourtftargetexongenes)
  browntfrnatargetintergenicgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintergenicgenes,brownrnasig))/length(ourtftargetintergenicgenes)
}

brownrnatargetsigdf <- data.frame("Percent.Significant" = c(browntfrnatargetpromproxgenesigpct[browntfrnatargetpromproxgenesigcount > 3],
                                                            browntfrnatargetintrongenesigpct[browntfrnatargetpromproxgenesigcount > 3],
                                                            browntfrnatargetexongenesigpct[browntfrnatargetpromproxgenesigcount > 3],
                                                            browntfrnatargetintergenicgenesigpct[browntfrnatargetpromproxgenesigcount > 3]),
                                  "Region" = c(rep("Promoter (<=1kb)",length(browntfrnatargetpromproxgenesigpct[browntfrnatargetpromproxgenesigcount > 3])),
                                               rep("Intron",length(browntfrnatargetpromproxgenesigpct[browntfrnatargetpromproxgenesigcount > 3])),
                                               rep("Exon",length(browntfrnatargetpromproxgenesigpct[browntfrnatargetpromproxgenesigcount > 3])),
                                               rep("Distal Intergenic",length(browntfrnatargetpromproxgenesigpct[browntfrnatargetpromproxgenesigcount > 3]))),
                                  "Transcription Factor" = rep(gsub("\\(.*","",rownames(browntfrnatargetpromproxgenesigpct)[browntfrnatargetpromproxgenesigcount > 3]),4))



# WAT-SC

whitetfrnatargetpromproxgenesigpct <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfwhiternasig]),ncol = 1)
rownames(whitetfrnatargetpromproxgenesigpct) <- rownames(tfanno)[tfanno$Ensembl %in% tfwhiternasig]
colnames(whitetfrnatargetpromproxgenesigpct) <- "Percent.Significant"

whitetfrnatargetintrongenesigpct <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfwhiternasig]),ncol = 1)
rownames(whitetfrnatargetintrongenesigpct) <- rownames(tfanno)[tfanno$Ensembl %in% tfwhiternasig]
colnames(whitetfrnatargetintrongenesigpct) <- "Percent.Significant"

whitetfrnatargetexongenesigpct <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfwhiternasig]),ncol = 1)
rownames(whitetfrnatargetexongenesigpct) <- rownames(tfanno)[tfanno$Ensembl %in% tfwhiternasig]
colnames(whitetfrnatargetexongenesigpct) <- "Percent.Significant"

whitetfrnatargetintergenicgenesigpct <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfwhiternasig]),ncol = 1)
rownames(whitetfrnatargetintergenicgenesigpct) <- rownames(tfanno)[tfanno$Ensembl %in% tfwhiternasig]
colnames(whitetfrnatargetintergenicgenesigpct) <- "Percent.Significant"


whitetfrnatargetpromproxgenesigcount <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfwhiternasig]),ncol = 1)
rownames(whitetfrnatargetpromproxgenesigcount) <- rownames(tfanno)[tfanno$Ensembl %in% tfwhiternasig]
colnames(whitetfrnatargetpromproxgenesigcount) <- "Count.Significant"

whitetfrnatargetintrongenesigcount <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfwhiternasig]),ncol = 1)
rownames(whitetfrnatargetintrongenesigcount) <- rownames(tfanno)[tfanno$Ensembl %in% tfwhiternasig]
colnames(whitetfrnatargetintrongenesigcount) <- "Count.Significant"

whitetfrnatargetexongenesigcount <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfwhiternasig]),ncol = 1)
rownames(whitetfrnatargetexongenesigcount) <- rownames(tfanno)[tfanno$Ensembl %in% tfwhiternasig]
colnames(whitetfrnatargetexongenesigcount) <- "Count.Significant"

whitetfrnatargetintergenicgenesigcount <- matrix(0L,nrow = length(rownames(tfanno)[tfanno$Ensembl %in% tfwhiternasig]),ncol = 1)
rownames(whitetfrnatargetintergenicgenesigcount) <- rownames(tfanno)[tfanno$Ensembl %in% tfwhiternasig]
colnames(whitetfrnatargetintergenicgenesigcount) <- "Count.Significant"


for(i in 1:length(rownames(tfanno)[tfanno$Ensembl %in% tfwhiternasig])){
  ourtf <- rownames(tfanno)[tfanno$Ensembl %in% tfwhiternasig][i]
  ourtftargetpeaks <- white50peakmotifs[white50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
  ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
  ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
  ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
  ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
  ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])
  whitetfrnatargetpromproxgenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetpromproxgenes,whiternasig))
  whitetfrnatargetintrongenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetintrongenes,whiternasig))
  whitetfrnatargetexongenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetexongenes,whiternasig))
  whitetfrnatargetintergenicgenesigcount[i,"Count.Significant"] <- length(intersect(ourtftargetintergenicgenes,whiternasig))
  whitetfrnatargetpromproxgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetpromproxgenes,whiternasig))/length(ourtftargetpromproxgenes)
  whitetfrnatargetintrongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintrongenes,whiternasig))/length(ourtftargetintrongenes)
  whitetfrnatargetexongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetexongenes,whiternasig))/length(ourtftargetexongenes)
  whitetfrnatargetintergenicgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintergenicgenes,whiternasig))/length(ourtftargetintergenicgenes)
}

whiternatargetsigdf <- data.frame("Percent.Significant" = c(whitetfrnatargetpromproxgenesigpct[whitetfrnatargetpromproxgenesigcount > 3],
                                                            whitetfrnatargetintrongenesigpct[whitetfrnatargetpromproxgenesigcount > 3],
                                                            whitetfrnatargetexongenesigpct[whitetfrnatargetpromproxgenesigcount > 3],
                                                            whitetfrnatargetintergenicgenesigpct[whitetfrnatargetpromproxgenesigcount > 3]),
                                  "Region" = c(rep("Promoter (<=1kb)",length(whitetfrnatargetpromproxgenesigpct[whitetfrnatargetpromproxgenesigcount > 3])),
                                               rep("Intron",length(whitetfrnatargetpromproxgenesigpct[whitetfrnatargetpromproxgenesigcount > 3])),
                                               rep("Exon",length(whitetfrnatargetpromproxgenesigpct[whitetfrnatargetpromproxgenesigcount > 3])),
                                               rep("Distal Intergenic",length(whitetfrnatargetpromproxgenesigpct[whitetfrnatargetpromproxgenesigcount > 3]))),
                                  "Transcription Factor" = rep(gsub("\\(.*","",rownames(whitetfrnatargetpromproxgenesigpct)[whitetfrnatargetpromproxgenesigcount > 3]),4))



pdf(file = "Supplemental Figure S15A.pdf",width = 12,height = 4.5)
ggplot(gastrornatargetsigdf[!gastrornatargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(gastrornasig)/dim(gastrol2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),axis.title = element_text(size = 13),legend.text = element_text(size = 12))
dev.off()

pdf(file = "Supplemental Figure S15B.pdf",width = 12,height = 4.5)
ggplot(heartrnatargetsigdf[!heartrnatargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(heartrnasig)/dim(heartl2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),axis.title = element_text(size = 13),legend.text = element_text(size = 12))
dev.off()

pdf(file = "Supplemental Figure S15C.pdf",width = 12,height = 4.5)
ggplot(hippornatargetsigdf[!hippornatargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(hippornasig)/dim(hippol2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),axis.title = element_text(size = 13),legend.text = element_text(size = 12))
dev.off()

pdf(file = "Supplemental Figure S15D.pdf",width = 12,height = 4.5)
ggplot(kidneyrnatargetsigdf[!kidneyrnatargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(kidneyrnasig)/dim(kidneyl2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),axis.title = element_text(size = 13),legend.text = element_text(size = 12))
dev.off()

pdf(file = "Supplemental Figure S15E.pdf",width = 12,height = 4.5)
ggplot(liverrnatargetsigdf[!liverrnatargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(liverrnasig)/dim(liverl2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),axis.title = element_text(size = 13),legend.text = element_text(size = 12))
dev.off()

pdf(file = "Supplemental Figure S15F.pdf",width = 12,height = 4.5)
ggplot(lungrnatargetsigdf[!lungrnatargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(lungrnasig)/dim(lungl2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),axis.title = element_text(size = 13),legend.text = element_text(size = 12))
dev.off()

pdf(file = "Supplemental Figure S15G.pdf",width = 12,height = 4.5)
ggplot(brownrnatargetsigdf[!brownrnatargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(brownrnasig)/dim(brownl2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),axis.title = element_text(size = 13),legend.text = element_text(size = 12))
dev.off()

pdf(file = "Supplemental Figure S15H.pdf",width = 12,height = 4.5)
ggplot(whiternatargetsigdf[!whiternatargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(whiternasig)/dim(whitel2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),axis.title = element_text(size = 13),legend.text = element_text(size = 12))
dev.off()

png(file = "Supplemental Figure S15A.png",width = 12,height = 4.5,units = "in",res = 600)
ggplot(gastrornatargetsigdf[!gastrornatargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(gastrornasig)/dim(gastrol2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),axis.title = element_text(size = 13),legend.text = element_text(size = 12))
dev.off()

png(file = "Supplemental Figure S15B.png",width = 12,height = 4.5,units = "in",res = 600)
ggplot(heartrnatargetsigdf[!heartrnatargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(heartrnasig)/dim(heartl2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),axis.title = element_text(size = 13),legend.text = element_text(size = 12))
dev.off()

png(file = "Supplemental Figure S15C.png",width = 12,height = 4.5,units = "in",res = 600)
ggplot(hippornatargetsigdf[!hippornatargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(hippornasig)/dim(hippol2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),axis.title = element_text(size = 13),legend.text = element_text(size = 12))
dev.off()

png(file = "Supplemental Figure S15D.png",width = 12,height = 4.5,units = "in",res = 600)
ggplot(kidneyrnatargetsigdf[!kidneyrnatargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(kidneyrnasig)/dim(kidneyl2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),axis.title = element_text(size = 13),legend.text = element_text(size = 12))
dev.off()

png(file = "Supplemental Figure S15E.png",width = 12,height = 4.5,units = "in",res = 600)
ggplot(liverrnatargetsigdf[!liverrnatargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(liverrnasig)/dim(liverl2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),axis.title = element_text(size = 13),legend.text = element_text(size = 12))
dev.off()

png(file = "Supplemental Figure S15F.png",width = 12,height = 4.5,units = "in",res = 600)
ggplot(lungrnatargetsigdf[!lungrnatargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(lungrnasig)/dim(lungl2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),axis.title = element_text(size = 13),legend.text = element_text(size = 12))
dev.off()

png(file = "Supplemental Figure S15G.png",width = 12,height = 4.5,units = "in",res = 600)
ggplot(brownrnatargetsigdf[!brownrnatargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(brownrnasig)/dim(brownl2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),axis.title = element_text(size = 13),legend.text = element_text(size = 12))
dev.off()

png(file = "Supplemental Figure S15H.png",width = 12,height = 4.5,units = "in",res = 600)
ggplot(whiternatargetsigdf[!whiternatargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(whiternasig)/dim(whitel2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),axis.title = element_text(size = 13),legend.text = element_text(size = 12))
dev.off()



ourtf <- "Six1(Homeobox)/Myoblast-Six1-ChIP-Chip(GSE20150)/Homer"
ourtftargetpeaks <- gastro50peakmotifs[gastro50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])

sextimemeta <- data.frame(row.names = colnames(gastrol2fcmat),
                          "Sex" = c(rep("Female",4),rep("Male",4)),
                          "Week" = rep(c("1w","2w","4w","8w"),2),
                          "TF.L2FC" = gastrol2fcmat[tfanno[ourtf,"Ensembl"],])
sextimemeta$TF.L2FC <- round(sextimemeta$TF.L2FC/0.1)*0.1
sextimemeta$TF.L2FC[sextimemeta$TF.L2FC > 1] <- 1
sextimemeta$TF.L2FC[sextimemeta$TF.L2FC < -1] <- -1
sextimemeta$TF.L2FC <- factor(sextimemeta$TF.L2FC,levels = c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8,-0.9,-1))

ourtftargetpromproxcormeta <- data.frame(row.names = intersect(ourtftargetpromproxgenes,gastrornasig),
                                         "Correlation" = rep(0,length(intersect(ourtftargetpromproxgenes,gastrornasig))))
ourtftargetexoncormeta <- data.frame(row.names = intersect(ourtftargetexongenes,gastrornasig),
                                     "Correlation" = rep(0,length(intersect(ourtftargetexongenes,gastrornasig))))
ourtftargetintroncormeta <- data.frame(row.names = intersect(ourtftargetintrongenes,gastrornasig),
                                       "Correlation" = rep(0,length(intersect(ourtftargetintrongenes,gastrornasig))))
for(i in 1:dim(ourtftargetpromproxcormeta)[1]){
  ourtftargetpromproxcormeta[i,"Correlation"] <- cor(gastrol2fcmat[tfanno[ourtf,"Ensembl"],],gastrol2fcmat[intersect(ourtftargetpromproxgenes,gastrornasig)[i],])
}
ourtftargetpromproxcormeta$Correlation <- round(ourtftargetpromproxcormeta$Correlation/0.2)*0.2
ourtftargetpromproxcormeta$Correlation <- factor(ourtftargetpromproxcormeta$Correlation,levels = c(1,0.8,0.6,0.4,0.2,0,-0.2,-0.4,-0.6,-0.8,-1))

png(file = "Figure 5E_Supplemental Figure S16A.png",width = 7,height = 4,units = "in",res = 600)
pheatmap(gastrol2fcmat[intersect(ourtftargetpromproxgenes,gastrornasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,gastrornasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F)
dev.off()

pdf(file = "Figure 5E_Supplemental Figure S16A.pdf",width = 7,height = 4)
pheatmap(gastrol2fcmat[intersect(ourtftargetpromproxgenes,gastrornasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,gastrornasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F)
dev.off()

ourtf <- "JunD(bZIP)/K562-JunD-ChIP-Seq/Homer"
ourtftargetpeaks <- gastro50peakmotifs[gastro50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])

sextimemeta <- data.frame(row.names = colnames(gastrol2fcmat),
                          "Sex" = c(rep("Female",4),rep("Male",4)),
                          "Week" = rep(c("1w","2w","4w","8w"),2),
                          "TF.L2FC" = gastrol2fcmat[tfanno[ourtf,"Ensembl"],])
sextimemeta$TF.L2FC <- round(sextimemeta$TF.L2FC/0.1)*0.1
sextimemeta$TF.L2FC[sextimemeta$TF.L2FC > 1] <- 1
sextimemeta$TF.L2FC[sextimemeta$TF.L2FC < -1] <- -1
sextimemeta$TF.L2FC <- factor(sextimemeta$TF.L2FC,levels = c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8,-0.9,-1))
ourtftargetpromproxcormeta <- data.frame(row.names = intersect(ourtftargetpromproxgenes,gastrornasig),
                                         "Correlation" = rep(0,length(intersect(ourtftargetpromproxgenes,gastrornasig))))
ourtftargetexoncormeta <- data.frame(row.names = intersect(ourtftargetexongenes,gastrornasig),
                                     "Correlation" = rep(0,length(intersect(ourtftargetexongenes,gastrornasig))))
ourtftargetintroncormeta <- data.frame(row.names = intersect(ourtftargetintrongenes,gastrornasig),
                                       "Correlation" = rep(0,length(intersect(ourtftargetintrongenes,gastrornasig))))
for(i in 1:dim(ourtftargetpromproxcormeta)[1]){
  ourtftargetpromproxcormeta[i,"Correlation"] <- cor(gastrol2fcmat[tfanno[ourtf,"Ensembl"],],gastrol2fcmat[intersect(ourtftargetpromproxgenes,gastrornasig)[i],])
}
ourtftargetpromproxcormeta$Correlation <- round(ourtftargetpromproxcormeta$Correlation/0.2)*0.2
ourtftargetpromproxcormeta$Correlation <- factor(ourtftargetpromproxcormeta$Correlation,levels = c(1,0.8,0.6,0.4,0.2,0,-0.2,-0.4,-0.6,-0.8,-1))

png(file = "Supplemental Figure S16B.png",width = 7,height = 4,units = "in",res = 600)
pheatmap(gastrol2fcmat[intersect(ourtftargetpromproxgenes,gastrornasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,gastrornasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F)
dev.off()

pdf(file = "Supplemental Figure S16B.pdf",width = 7,height = 4)
pheatmap(gastrol2fcmat[intersect(ourtftargetpromproxgenes,gastrornasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,gastrornasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F)
dev.off()



ourtf <- "SF1(NR)/H295R-Nr5a1-ChIP-Seq(GSE44220)/Homer"
ourtftargetpeaks <- gastro50peakmotifs[gastro50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])

sextimemeta <- data.frame(row.names = colnames(gastrol2fcmat),
                          "Sex" = c(rep("Female",4),rep("Male",4)),
                          "Week" = rep(c("1w","2w","4w","8w"),2),
                          "TF.L2FC" = gastrol2fcmat[tfanno[ourtf,"Ensembl"],])
sextimemeta$TF.L2FC <- round(sextimemeta$TF.L2FC/0.1)*0.1
sextimemeta$TF.L2FC[sextimemeta$TF.L2FC > 1] <- 1
sextimemeta$TF.L2FC[sextimemeta$TF.L2FC < -1] <- -1
sextimemeta$TF.L2FC <- factor(sextimemeta$TF.L2FC,levels = c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8,-0.9,-1))

ourtftargetpromproxcormeta <- data.frame(row.names = intersect(ourtftargetpromproxgenes,gastrornasig),
                                         "Correlation" = rep(0,length(intersect(ourtftargetpromproxgenes,gastrornasig))))
ourtftargetexoncormeta <- data.frame(row.names = intersect(ourtftargetexongenes,gastrornasig),
                                     "Correlation" = rep(0,length(intersect(ourtftargetexongenes,gastrornasig))))
ourtftargetintroncormeta <- data.frame(row.names = intersect(ourtftargetintrongenes,gastrornasig),
                                       "Correlation" = rep(0,length(intersect(ourtftargetintrongenes,gastrornasig))))
for(i in 1:dim(ourtftargetpromproxcormeta)[1]){
  ourtftargetpromproxcormeta[i,"Correlation"] <- cor(gastrol2fcmat[tfanno[ourtf,"Ensembl"],],gastrol2fcmat[intersect(ourtftargetpromproxgenes,gastrornasig)[i],])
}
ourtftargetpromproxcormeta$Correlation <- round(ourtftargetpromproxcormeta$Correlation/0.2)*0.2
ourtftargetpromproxcormeta$Correlation <- factor(ourtftargetpromproxcormeta$Correlation,levels = c(1,0.8,0.6,0.4,0.2,0,-0.2,-0.4,-0.6,-0.8,-1))

png(file = "Supplemental Figure S16C.png",width = 7,height = 6,units = "in",res = 600)
pheatmap(gastrol2fcmat[intersect(ourtftargetpromproxgenes,gastrornasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,gastrornasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F)
dev.off()

pdf(file = "Supplemental Figure S16C.pdf",width = 7,height = 6)
pheatmap(gastrol2fcmat[intersect(ourtftargetpromproxgenes,gastrornasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,gastrornasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F)
dev.off()



ourtf <- "JunD(bZIP)/K562-JunD-ChIP-Seq/Homer"
ourtftargetpeaks <- heart50peakmotifs[heart50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])

sextimemeta <- data.frame(row.names = colnames(gastrol2fcmat),
                          "Sex" = c(rep("Female",4),rep("Male",4)),
                          "Week" = rep(c("1w","2w","4w","8w"),2),
                          "TF.L2FC" = gastrol2fcmat[tfanno[ourtf,"Ensembl"],])
sextimemeta$TF.L2FC <- round(sextimemeta$TF.L2FC/0.1)*0.1
sextimemeta$TF.L2FC[sextimemeta$TF.L2FC > 1] <- 1
sextimemeta$TF.L2FC[sextimemeta$TF.L2FC < -1] <- -1
sextimemeta$TF.L2FC <- factor(sextimemeta$TF.L2FC,levels = c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8,-0.9,-1))
ourtftargetpromproxcormeta <- data.frame(row.names = intersect(ourtftargetpromproxgenes,heartrnasig),
                                         "Correlation" = rep(0,length(intersect(ourtftargetpromproxgenes,heartrnasig))))
ourtftargetexoncormeta <- data.frame(row.names = intersect(ourtftargetexongenes,heartrnasig),
                                     "Correlation" = rep(0,length(intersect(ourtftargetexongenes,heartrnasig))))
ourtftargetintroncormeta <- data.frame(row.names = intersect(ourtftargetintrongenes,heartrnasig),
                                       "Correlation" = rep(0,length(intersect(ourtftargetintrongenes,heartrnasig))))
for(i in 1:dim(ourtftargetpromproxcormeta)[1]){
  ourtftargetpromproxcormeta[i,"Correlation"] <- cor(heartl2fcmat[tfanno[ourtf,"Ensembl"],],heartl2fcmat[intersect(ourtftargetpromproxgenes,heartrnasig)[i],])
}
ourtftargetpromproxcormeta$Correlation <- round(ourtftargetpromproxcormeta$Correlation/0.2)*0.2
ourtftargetpromproxcormeta$Correlation <- factor(ourtftargetpromproxcormeta$Correlation,levels = c(1,0.8,0.6,0.4,0.2,0,-0.2,-0.4,-0.6,-0.8,-1))

png(file = "Supplemental Figure S16D.png",width = 7,height = 4,units = "in",res = 600)
pheatmap(heartl2fcmat[intersect(ourtftargetpromproxgenes,heartrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,heartrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F)
dev.off()

pdf(file = "Supplemental Figure S16D.pdf",width = 7,height = 4)
pheatmap(heartl2fcmat[intersect(ourtftargetpromproxgenes,heartrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,heartrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F)
dev.off()


ourtf <- "Fos(bZIP)/TSC-Fos-ChIP-Seq(GSE110950)/Homer" 
ourtftargetpeaks <- kidney50peakmotifs[kidney50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])

sextimemeta <- data.frame(row.names = colnames(gastrol2fcmat),
                          "Sex" = c(rep("Female",4),rep("Male",4)),
                          "Week" = rep(c("1w","2w","4w","8w"),2),
                          "TF.L2FC" = gastrol2fcmat[tfanno[ourtf,"Ensembl"],])
sextimemeta$TF.L2FC <- round(sextimemeta$TF.L2FC/0.1)*0.1
sextimemeta$TF.L2FC[sextimemeta$TF.L2FC > 1] <- 1
sextimemeta$TF.L2FC[sextimemeta$TF.L2FC < -1] <- -1
sextimemeta$TF.L2FC <- factor(sextimemeta$TF.L2FC,levels = c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8,-0.9,-1))
ourtftargetpromproxcormeta <- data.frame(row.names = intersect(ourtftargetpromproxgenes,kidneyrnasig),
                                         "Correlation" = rep(0,length(intersect(ourtftargetpromproxgenes,kidneyrnasig))))
ourtftargetexoncormeta <- data.frame(row.names = intersect(ourtftargetexongenes,kidneyrnasig),
                                     "Correlation" = rep(0,length(intersect(ourtftargetexongenes,kidneyrnasig))))
ourtftargetintroncormeta <- data.frame(row.names = intersect(ourtftargetintrongenes,kidneyrnasig),
                                       "Correlation" = rep(0,length(intersect(ourtftargetintrongenes,kidneyrnasig))))
for(i in 1:dim(ourtftargetpromproxcormeta)[1]){
  ourtftargetpromproxcormeta[i,"Correlation"] <- cor(kidneyl2fcmat[tfanno[ourtf,"Ensembl"],],kidneyl2fcmat[intersect(ourtftargetpromproxgenes,kidneyrnasig)[i],])
}
ourtftargetpromproxcormeta$Correlation <- round(ourtftargetpromproxcormeta$Correlation/0.2)*0.2
ourtftargetpromproxcormeta$Correlation <- factor(ourtftargetpromproxcormeta$Correlation,levels = c(1,0.8,0.6,0.4,0.2,0,-0.2,-0.4,-0.6,-0.8,-1))

png(file = "Supplemental Figure S16E.png",width = 7,height = 4,units = "in",res = 600)
pheatmap(kidneyl2fcmat[intersect(ourtftargetpromproxgenes,kidneyrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,kidneyrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F)
dev.off()

pdf(file = "Supplemental Figure S16E.pdf",width = 7,height = 4)
pheatmap(kidneyl2fcmat[intersect(ourtftargetpromproxgenes,kidneyrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,kidneyrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F)
dev.off()



ourtf <- "JunD(bZIP)/K562-JunD-ChIP-Seq/Homer"
ourtftargetpeaks <- white50peakmotifs[white50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])

sextimemeta <- data.frame(row.names = colnames(gastrol2fcmat),
                          "Sex" = c(rep("Female",4),rep("Male",4)),
                          "Week" = rep(c("1w","2w","4w","8w"),2),
                          "TF.L2FC" = gastrol2fcmat[tfanno[ourtf,"Ensembl"],])
sextimemeta$TF.L2FC <- round(sextimemeta$TF.L2FC/0.1)*0.1
sextimemeta$TF.L2FC[sextimemeta$TF.L2FC > 1] <- 1
sextimemeta$TF.L2FC[sextimemeta$TF.L2FC < -1] <- -1
sextimemeta$TF.L2FC <- factor(sextimemeta$TF.L2FC,levels = c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8,-0.9,-1))
ourtftargetpromproxcormeta <- data.frame(row.names = intersect(ourtftargetpromproxgenes,whiternasig),
                                         "Correlation" = rep(0,length(intersect(ourtftargetpromproxgenes,whiternasig))))
ourtftargetexoncormeta <- data.frame(row.names = intersect(ourtftargetexongenes,whiternasig),
                                     "Correlation" = rep(0,length(intersect(ourtftargetexongenes,whiternasig))))
ourtftargetintroncormeta <- data.frame(row.names = intersect(ourtftargetintrongenes,whiternasig),
                                       "Correlation" = rep(0,length(intersect(ourtftargetintrongenes,whiternasig))))
for(i in 1:dim(ourtftargetpromproxcormeta)[1]){
  ourtftargetpromproxcormeta[i,"Correlation"] <- cor(whitel2fcmat[tfanno[ourtf,"Ensembl"],],whitel2fcmat[intersect(ourtftargetpromproxgenes,whiternasig)[i],])
}
ourtftargetpromproxcormeta$Correlation <- round(ourtftargetpromproxcormeta$Correlation/0.2)*0.2
ourtftargetpromproxcormeta$Correlation <- factor(ourtftargetpromproxcormeta$Correlation,levels = c(1,0.8,0.6,0.4,0.2,0,-0.2,-0.4,-0.6,-0.8,-1))

png(file = "Supplemental Figure S16F.png",width = 7,height = 4,units = "in",res = 600)
pheatmap(whitel2fcmat[intersect(ourtftargetpromproxgenes,whiternasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,whiternasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F)
dev.off()

pdf(file = "Supplemental Figure S16F.pdf",width = 7,height = 4)
pheatmap(whitel2fcmat[intersect(ourtftargetpromproxgenes,whiternasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,whiternasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F)
dev.off()


ourtf <- "PBX1(Homeobox)/MCF7-PBX1-ChIP-Seq(GSE28007)/Homer"
ourtftargetpeaks <- lung50peakmotifs[lung50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])

sextimemeta <- data.frame(row.names = colnames(gastrol2fcmat),
                          "Sex" = c(rep("Female",4),rep("Male",4)),
                          "Week" = rep(c("1w","2w","4w","8w"),2),
                          "TF.L2FC" = gastrol2fcmat[tfanno[ourtf,"Ensembl"],])
sextimemeta$TF.L2FC <- round(sextimemeta$TF.L2FC/0.1)*0.1
sextimemeta$TF.L2FC[sextimemeta$TF.L2FC > 1] <- 1
sextimemeta$TF.L2FC[sextimemeta$TF.L2FC < -1] <- -1
sextimemeta$TF.L2FC <- factor(sextimemeta$TF.L2FC,levels = c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8,-0.9,-1))
ourtftargetpromproxcormeta <- data.frame(row.names = intersect(ourtftargetpromproxgenes,lungrnasig),
                                         "Correlation" = rep(0,length(intersect(ourtftargetpromproxgenes,lungrnasig))))
ourtftargetexoncormeta <- data.frame(row.names = intersect(ourtftargetexongenes,lungrnasig),
                                     "Correlation" = rep(0,length(intersect(ourtftargetexongenes,lungrnasig))))
ourtftargetintroncormeta <- data.frame(row.names = intersect(ourtftargetintrongenes,lungrnasig),
                                       "Correlation" = rep(0,length(intersect(ourtftargetintrongenes,lungrnasig))))
for(i in 1:dim(ourtftargetpromproxcormeta)[1]){
  ourtftargetpromproxcormeta[i,"Correlation"] <- cor(lungl2fcmat[tfanno[ourtf,"Ensembl"],],lungl2fcmat[intersect(ourtftargetpromproxgenes,lungrnasig)[i],])
}
ourtftargetpromproxcormeta$Correlation <- round(ourtftargetpromproxcormeta$Correlation/0.2)*0.2
ourtftargetpromproxcormeta$Correlation <- factor(ourtftargetpromproxcormeta$Correlation,levels = c(1,0.8,0.6,0.4,0.2,0,-0.2,-0.4,-0.6,-0.8,-1))

png(file = "Supplemental Figure S16G.png",width = 7,height = 4,units = "in",res = 600)
pheatmap(lungl2fcmat[intersect(ourtftargetpromproxgenes,lungrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,lungrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F)
dev.off()

pdf(file = "Supplemental Figure S16G.pdf",width = 7,height = 4)
pheatmap(lungl2fcmat[intersect(ourtftargetpromproxgenes,lungrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,lungrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F)
dev.off()


ourtf <- "RAR:RXR(NR),DR0/ES-RAR-ChIP-Seq(GSE56893)/Homer"
ourtftargetpeaks <- brown50peakmotifs[brown50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])

sextimemeta <- data.frame(row.names = colnames(gastrol2fcmat),
                          "Sex" = c(rep("Female",4),rep("Male",4)),
                          "Week" = rep(c("1w","2w","4w","8w"),2),
                          "TF.L2FC" = gastrol2fcmat[tfanno[ourtf,"Ensembl"],])
sextimemeta$TF.L2FC <- round(sextimemeta$TF.L2FC/0.1)*0.1
sextimemeta$TF.L2FC[sextimemeta$TF.L2FC > 1] <- 1
sextimemeta$TF.L2FC[sextimemeta$TF.L2FC < -1] <- -1
sextimemeta$TF.L2FC <- factor(sextimemeta$TF.L2FC,levels = c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8,-0.9,-1))
ourtftargetpromproxcormeta <- data.frame(row.names = intersect(ourtftargetpromproxgenes,brownrnasig),
                                         "Correlation" = rep(0,length(intersect(ourtftargetpromproxgenes,brownrnasig))))
ourtftargetexoncormeta <- data.frame(row.names = intersect(ourtftargetexongenes,brownrnasig),
                                     "Correlation" = rep(0,length(intersect(ourtftargetexongenes,brownrnasig))))
ourtftargetintroncormeta <- data.frame(row.names = intersect(ourtftargetintrongenes,brownrnasig),
                                       "Correlation" = rep(0,length(intersect(ourtftargetintrongenes,brownrnasig))))
for(i in 1:dim(ourtftargetpromproxcormeta)[1]){
  ourtftargetpromproxcormeta[i,"Correlation"] <- cor(brownl2fcmat[tfanno[ourtf,"Ensembl"],],brownl2fcmat[intersect(ourtftargetpromproxgenes,brownrnasig)[i],])
}
ourtftargetpromproxcormeta$Correlation <- round(ourtftargetpromproxcormeta$Correlation/0.2)*0.2
ourtftargetpromproxcormeta$Correlation <- factor(ourtftargetpromproxcormeta$Correlation,levels = c(1,0.8,0.6,0.4,0.2,0,-0.2,-0.4,-0.6,-0.8,-1))

png(file = "Supplemental Figure S16H.png",width = 7,height = 4,units = "in",res = 600)
pheatmap(brownl2fcmat[intersect(ourtftargetpromproxgenes,brownrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,brownrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F)
dev.off()

pdf(file = "Supplemental Figure S16H.pdf",width = 7,height = 4)
pheatmap(brownl2fcmat[intersect(ourtftargetpromproxgenes,brownrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,brownrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F)
dev.off()


#####
# Figure 5D
####

# Calculate the exact binomial test for all of our differential TFs across tissues

# SKM-GN
gastrotfrnatargetbinomtest <- matrix(1L,nrow = length(gastrotfrnatargetpromproxgenesigcount),ncol = 1)
rownames(gastrotfrnatargetbinomtest) <- rownames(gastrotfrnatargetpromproxgenesigcount)
colnames(gastrotfrnatargetbinomtest) <- "Binom.Test.pval"
for(i in 1:length(gastrotfrnatargetpromproxgenesigcount)){
  ourtf <- rownames(gastrotfrnatargetbinomtest)[i]
  if(gastrotfrnatargetpromproxgenesigcount[i] > 0){
    ourtest <- binom.test(x = gastrotfrnatargetpromproxgenesigcount[i],n = gastrotfrnatargetpromproxgenesigcount[i]/gastrotfrnatargetpromproxgenesigpct[i],p = length(gastrornasig)/dim(gastrol2fcmat)[1],alternative = "greater")
    gastrotfrnatargetbinomtest[i,"Binom.Test.pval"] <- ourtest$p.value
  }
}

gastrotfprotargetbinomtest <- matrix(1L,nrow = length(gastrotftargetpromproxgenesigcount),ncol = 1)
rownames(gastrotfprotargetbinomtest) <- rownames(gastrotftargetpromproxgenesigcount)
colnames(gastrotfprotargetbinomtest) <- "Binom.Test.pval"
for(i in 1:length(gastrotftargetpromproxgenesigcount)){
  ourtf <- rownames(gastrotfprotargetbinomtest)[i]
  if(gastrotftargetpromproxgenesigcount[i] > 0){
    ourtest <- binom.test(x = gastrotftargetpromproxgenesigcount[i],n = gastrotftargetpromproxgenesigcount[i]/gastrotftargetpromproxgenesigpct[i],p = length(gastrornasig)/dim(gastrol2fcmat)[1],alternative = "greater")
    gastrotfprotargetbinomtest[i,"Binom.Test.pval"] <- ourtest$p.value
  }
}


# HEART
hearttfrnatargetbinomtest <- matrix(1L,nrow = length(hearttfrnatargetpromproxgenesigcount),ncol = 1)
rownames(hearttfrnatargetbinomtest) <- rownames(hearttfrnatargetpromproxgenesigcount)
colnames(hearttfrnatargetbinomtest) <- "Binom.Test.pval"
for(i in 1:length(hearttfrnatargetpromproxgenesigcount)){
  ourtf <- rownames(hearttfrnatargetbinomtest)[i]
  if(hearttfrnatargetpromproxgenesigcount[i] > 0){
    ourtest <- binom.test(x = hearttfrnatargetpromproxgenesigcount[i],n = hearttfrnatargetpromproxgenesigcount[i]/hearttfrnatargetpromproxgenesigpct[i],p = length(heartrnasig)/dim(heartl2fcmat)[1],alternative = "greater")
    hearttfrnatargetbinomtest[i,"Binom.Test.pval"] <- ourtest$p.value
  }
}

hearttfprotargetbinomtest <- matrix(1L,nrow = length(hearttftargetpromproxgenesigcount),ncol = 1)
rownames(hearttfprotargetbinomtest) <- rownames(hearttftargetpromproxgenesigcount)
colnames(hearttfprotargetbinomtest) <- "Binom.Test.pval"
for(i in 1:length(hearttftargetpromproxgenesigcount)){
  ourtf <- rownames(hearttfprotargetbinomtest)[i]
  if(hearttftargetpromproxgenesigcount[i] > 0){
    ourtest <- binom.test(x = hearttftargetpromproxgenesigcount[i],n = hearttftargetpromproxgenesigcount[i]/hearttftargetpromproxgenesigpct[i],p = length(heartrnasig)/dim(heartl2fcmat)[1],alternative = "greater")
    hearttfprotargetbinomtest[i,"Binom.Test.pval"] <- ourtest$p.value
  }
}


# HIPPOC
hippotfrnatargetbinomtest <- matrix(1L,nrow = length(hippotfrnatargetpromproxgenesigcount),ncol = 1)
rownames(hippotfrnatargetbinomtest) <- rownames(hippotfrnatargetpromproxgenesigcount)
colnames(hippotfrnatargetbinomtest) <- "Binom.Test.pval"
for(i in 1:length(hippotfrnatargetpromproxgenesigcount)){
  ourtf <- rownames(hippotfrnatargetbinomtest)[i]
  if(hippotfrnatargetpromproxgenesigcount[i] > 0){
    ourtest <- binom.test(x = hippotfrnatargetpromproxgenesigcount[i],n = hippotfrnatargetpromproxgenesigcount[i]/hippotfrnatargetpromproxgenesigpct[i],p = length(hippornasig)/dim(hippol2fcmat)[1],alternative = "greater")
    hippotfrnatargetbinomtest[i,"Binom.Test.pval"] <- ourtest$p.value
  }
}

# KIDNEY
kidneytfrnatargetbinomtest <- matrix(1L,nrow = length(kidneytfrnatargetpromproxgenesigcount),ncol = 1)
rownames(kidneytfrnatargetbinomtest) <- rownames(kidneytfrnatargetpromproxgenesigcount)
colnames(kidneytfrnatargetbinomtest) <- "Binom.Test.pval"
for(i in 1:length(kidneytfrnatargetpromproxgenesigcount)){
  ourtf <- rownames(kidneytfrnatargetbinomtest)[i]
  if(kidneytfrnatargetpromproxgenesigcount[i] > 0){
    ourtest <- binom.test(x = kidneytfrnatargetpromproxgenesigcount[i],n = kidneytfrnatargetpromproxgenesigcount[i]/kidneytfrnatargetpromproxgenesigpct[i],p = length(kidneyrnasig)/dim(kidneyl2fcmat)[1],alternative = "greater")
    kidneytfrnatargetbinomtest[i,"Binom.Test.pval"] <- ourtest$p.value
  }
}

kidneytfprotargetbinomtest <- matrix(1L,nrow = length(kidneytftargetpromproxgenesigcount),ncol = 1)
rownames(kidneytfprotargetbinomtest) <- rownames(kidneytftargetpromproxgenesigcount)
colnames(kidneytfprotargetbinomtest) <- "Binom.Test.pval"
for(i in 1:length(kidneytftargetpromproxgenesigcount)){
  ourtf <- rownames(kidneytfprotargetbinomtest)[i]
  if(kidneytftargetpromproxgenesigcount[i] > 0){
    ourtest <- binom.test(x = kidneytftargetpromproxgenesigcount[i],n = kidneytftargetpromproxgenesigcount[i]/kidneytftargetpromproxgenesigpct[i],p = length(kidneyrnasig)/dim(kidneyl2fcmat)[1],alternative = "greater")
    kidneytfprotargetbinomtest[i,"Binom.Test.pval"] <- ourtest$p.value
  }
}


# LIVER
livertfrnatargetbinomtest <- matrix(1L,nrow = length(livertfrnatargetpromproxgenesigcount),ncol = 1)
rownames(livertfrnatargetbinomtest) <- rownames(livertfrnatargetpromproxgenesigcount)
colnames(livertfrnatargetbinomtest) <- "Binom.Test.pval"
for(i in 1:length(livertfrnatargetpromproxgenesigcount)){
  ourtf <- rownames(livertfrnatargetbinomtest)[i]
  if(livertfrnatargetpromproxgenesigcount[i] > 0){
    ourtest <- binom.test(x = livertfrnatargetpromproxgenesigcount[i],n = livertfrnatargetpromproxgenesigcount[i]/livertfrnatargetpromproxgenesigpct[i],p = length(liverrnasig)/dim(liverl2fcmat)[1],alternative = "greater")
    livertfrnatargetbinomtest[i,"Binom.Test.pval"] <- ourtest$p.value
  }
}

livertfprotargetbinomtest <- matrix(1L,nrow = length(livertftargetpromproxgenesigcount),ncol = 1)
rownames(livertfprotargetbinomtest) <- rownames(livertftargetpromproxgenesigcount)
colnames(livertfprotargetbinomtest) <- "Binom.Test.pval"
for(i in 1:length(livertftargetpromproxgenesigcount)){
  ourtf <- rownames(livertfprotargetbinomtest)[i]
  if(livertftargetpromproxgenesigcount[i] > 0){
    ourtest <- binom.test(x = livertftargetpromproxgenesigcount[i],n = livertftargetpromproxgenesigcount[i]/livertftargetpromproxgenesigpct[i],p = length(liverrnasig)/dim(liverl2fcmat)[1],alternative = "greater")
    livertfprotargetbinomtest[i,"Binom.Test.pval"] <- ourtest$p.value
  }
}


# LUNG
lungtfrnatargetbinomtest <- matrix(1L,nrow = length(lungtfrnatargetpromproxgenesigcount),ncol = 1)
rownames(lungtfrnatargetbinomtest) <- rownames(lungtfrnatargetpromproxgenesigcount)
colnames(lungtfrnatargetbinomtest) <- "Binom.Test.pval"
for(i in 1:length(lungtfrnatargetpromproxgenesigcount)){
  ourtf <- rownames(lungtfrnatargetbinomtest)[i]
  if(lungtfrnatargetpromproxgenesigcount[i] > 0){
    ourtest <- binom.test(x = lungtfrnatargetpromproxgenesigcount[i],n = lungtfrnatargetpromproxgenesigcount[i]/lungtfrnatargetpromproxgenesigpct[i],p = length(lungrnasig)/dim(lungl2fcmat)[1],alternative = "greater")
    lungtfrnatargetbinomtest[i,"Binom.Test.pval"] <- ourtest$p.value
  }
}

lungtfprotargetbinomtest <- matrix(1L,nrow = length(lungtftargetpromproxgenesigcount),ncol = 1)
rownames(lungtfprotargetbinomtest) <- rownames(lungtftargetpromproxgenesigcount)
colnames(lungtfprotargetbinomtest) <- "Binom.Test.pval"
for(i in 1:length(lungtftargetpromproxgenesigcount)){
  ourtf <- rownames(lungtfprotargetbinomtest)[i]
  if(lungtftargetpromproxgenesigcount[i] > 0){
    ourtest <- binom.test(x = lungtftargetpromproxgenesigcount[i],n = lungtftargetpromproxgenesigcount[i]/lungtftargetpromproxgenesigpct[i],p = length(lungrnasig)/dim(lungl2fcmat)[1],alternative = "greater")
    lungtfprotargetbinomtest[i,"Binom.Test.pval"] <- ourtest$p.value
  }
}


# BAT
browntfrnatargetbinomtest <- matrix(1L,nrow = length(browntfrnatargetpromproxgenesigcount),ncol = 1)
rownames(browntfrnatargetbinomtest) <- rownames(browntfrnatargetpromproxgenesigcount)
colnames(browntfrnatargetbinomtest) <- "Binom.Test.pval"
for(i in 1:length(browntfrnatargetpromproxgenesigcount)){
  ourtf <- rownames(browntfrnatargetbinomtest)[i]
  if(browntfrnatargetpromproxgenesigcount[i] > 0){
    ourtest <- binom.test(x = browntfrnatargetpromproxgenesigcount[i],n = browntfrnatargetpromproxgenesigcount[i]/browntfrnatargetpromproxgenesigpct[i],p = length(brownrnasig)/dim(brownl2fcmat)[1],alternative = "greater")
    browntfrnatargetbinomtest[i,"Binom.Test.pval"] <- ourtest$p.value
  }
}



# WAT-SC
whitetfrnatargetbinomtest <- matrix(1L,nrow = length(whitetfrnatargetpromproxgenesigcount),ncol = 1)
rownames(whitetfrnatargetbinomtest) <- rownames(whitetfrnatargetpromproxgenesigcount)
colnames(whitetfrnatargetbinomtest) <- "Binom.Test.pval"
for(i in 1:length(whitetfrnatargetpromproxgenesigcount)){
  ourtf <- rownames(whitetfrnatargetbinomtest)[i]
  if(whitetfrnatargetpromproxgenesigcount[i] > 0){
    ourtest <- binom.test(x = whitetfrnatargetpromproxgenesigcount[i],n = whitetfrnatargetpromproxgenesigcount[i]/whitetfrnatargetpromproxgenesigpct[i],p = length(whiternasig)/dim(whitel2fcmat)[1],alternative = "greater")
    whitetfrnatargetbinomtest[i,"Binom.Test.pval"] <- ourtest$p.value
  }
}

whitetfprotargetbinomtest <- matrix(1L,nrow = length(whitetftargetpromproxgenesigcount),ncol = 1)
rownames(whitetfprotargetbinomtest) <- rownames(whitetftargetpromproxgenesigcount)
colnames(whitetfprotargetbinomtest) <- "Binom.Test.pval"
for(i in 1:length(whitetftargetpromproxgenesigcount)){
  ourtf <- rownames(whitetfprotargetbinomtest)[i]
  if(whitetftargetpromproxgenesigcount[i] > 0){
    ourtest <- binom.test(x = whitetftargetpromproxgenesigcount[i],n = whitetftargetpromproxgenesigcount[i]/whitetftargetpromproxgenesigpct[i],p = length(whiternasig)/dim(whitel2fcmat)[1],alternative = "greater")
    whitetfprotargetbinomtest[i,"Binom.Test.pval"] <- ourtest$p.value
  }
}



####

gastrornatargetsigdf$Transcription.Factor.Name <- gastrornatargetsigdf$Transcription.Factor
heartrnatargetsigdf$Transcription.Factor.Name <- heartrnatargetsigdf$Transcription.Factor
hippornatargetsigdf$Transcription.Factor.Name <- hippornatargetsigdf$Transcription.Factor
kidneyrnatargetsigdf$Transcription.Factor.Name <- kidneyrnatargetsigdf$Transcription.Factor
liverrnatargetsigdf$Transcription.Factor.Name <- liverrnatargetsigdf$Transcription.Factor
lungrnatargetsigdf$Transcription.Factor.Name <- lungrnatargetsigdf$Transcription.Factor
brownrnatargetsigdf$Transcription.Factor.Name <- brownrnatargetsigdf$Transcription.Factor
whiternatargetsigdf$Transcription.Factor.Name <- whiternatargetsigdf$Transcription.Factor
gastroprotargetsigdf$Transcription.Factor.Name <- gastroprotargetsigdf$Transcription.Factor.Label
heartprotargetsigdf$Transcription.Factor.Name <- heartprotargetsigdf$Transcription.Factor.Label
kidneyprotargetsigdf$Transcription.Factor.Name <- kidneyprotargetsigdf$Transcription.Factor.Label
liverprotargetsigdf$Transcription.Factor.Name <- liverprotargetsigdf$Transcription.Factor.Label
lungprotargetsigdf$Transcription.Factor.Name <- lungprotargetsigdf$Transcription.Factor.Label
whiteprotargetsigdf$Transcription.Factor.Name <- whiteprotargetsigdf$Transcription.Factor.Label



mergedtftargetsigdf <- rbind(gastrornatargetsigdf[gastrornatargetsigdf$Region %in% "Promoter (<=1kb)" & 
                                                    gastrornatargetsigdf$Transcription.Factor.Name %in% c("SF1","Six1","Six2","JunD","HIF2a"),c("Percent.Significant","Transcription.Factor.Name")],
                             gastroprotargetsigdf[gastroprotargetsigdf$Region %in% "Promoter (<=1kb)" & 
                                                    gastroprotargetsigdf$Transcription.Factor.Name %in% c("Mef2c","Mef2d","Nur77","NFAT"),c("Percent.Significant","Transcription.Factor.Name")],
                             heartrnatargetsigdf[heartrnatargetsigdf$Region %in% "Promoter (<=1kb)" & 
                                                   heartrnatargetsigdf$Transcription.Factor.Name %in% c("JunD","IRF1","HIF2a"),c("Percent.Significant","Transcription.Factor.Name")],
                             heartprotargetsigdf[heartprotargetsigdf$Region %in% "Promoter (<=1kb)" & 
                                                   heartprotargetsigdf$Transcription.Factor.Name %in% c("Mef2a","ZNF692"),c("Percent.Significant","Transcription.Factor.Name")],
                             kidneyrnatargetsigdf[kidneyrnatargetsigdf$Region %in% "Promoter (<=1kb)" & 
                                                    kidneyrnatargetsigdf$Transcription.Factor.Name %in% c("Egr1","Fos"),c("Percent.Significant","Transcription.Factor.Name")],
                             kidneyprotargetsigdf[kidneyprotargetsigdf$Region %in% "Promoter (<=1kb)" & 
                                                    kidneyprotargetsigdf$Transcription.Factor.Name %in% c("Atf7","NFAT"),c("Percent.Significant","Transcription.Factor.Name")],
                             liverprotargetsigdf[liverprotargetsigdf$Region %in% "Promoter (<=1kb)" & 
                                                   liverprotargetsigdf$Transcription.Factor.Name %in% c("Atf2","STAT1"),c("Percent.Significant","Transcription.Factor.Name")],
                             lungrnatargetsigdf[lungrnatargetsigdf$Region %in% "Promoter (<=1kb)" & 
                                                  lungrnatargetsigdf$Transcription.Factor.Name %in% c("PBX1","PU.1","IRF4"),c("Percent.Significant","Transcription.Factor.Name")],
                             lungprotargetsigdf[lungprotargetsigdf$Region %in% "Promoter (<=1kb)" & 
                                                  lungprotargetsigdf$Transcription.Factor.Name %in% c("RUNX","IRF1","IRF8","IRF:BATF","IRF4"),c("Percent.Significant","Transcription.Factor.Name")],
                             brownrnatargetsigdf[brownrnatargetsigdf$Region %in% "Promoter (<=1kb)" & 
                                                   brownrnatargetsigdf$Transcription.Factor.Name %in% c("EAR2","RAR:RXR","IRF1"),c("Percent.Significant","Transcription.Factor.Name")],
                             whiternatargetsigdf[whiternatargetsigdf$Region %in% "Promoter (<=1kb)" & 
                                                   whiternatargetsigdf$Transcription.Factor.Name %in% c("JunD","HRE"),c("Percent.Significant","Transcription.Factor.Name")],
                             whiteprotargetsigdf[whiteprotargetsigdf$Region %in% "Promoter (<=1kb)" & 
                                                   whiteprotargetsigdf$Transcription.Factor.Name %in% c("PBX1"),c("Percent.Significant","Transcription.Factor.Name")])
mergedtftargetsigdf <- mergedtftargetsigdf[c(1:31,33:35),]
mergedtftargetsigdf$Tissue <- c(rep("SKM-GN",9),
                                rep("HEART",5),
                                rep("KIDNEY",4),
                                rep("LIVER",2),
                                rep("LUNG",8),
                                rep("BAT",3),
                                rep("WAT-SC",3))
mergedtftargetsigdf$DEG.Enrichment <- c(rep(length(gastrornasig)/dim(gastrol2fcmat)[1],9),
                                        rep(length(heartrnasig)/dim(heartl2fcmat)[1],5),
                                        rep(length(kidneyrnasig)/dim(kidneyl2fcmat)[1],4),
                                        rep(length(liverrnasig)/dim(liverl2fcmat)[1],2),
                                        rep(length(lungrnasig)/dim(lungl2fcmat)[1],8),
                                        rep(length(brownrnasig)/dim(brownl2fcmat)[1],3),
                                        rep(length(whiternasig)/dim(whitel2fcmat)[1],3))
mergedtftargetsigdf$Ome <- c("RNA","RNA","RNA","RNA","RNA","Prot","Prot","Phos","Phos",
                             "RNA","RNA","RNA","Prot","Phos",
                             "RNA","RNA","Phos","Phos",
                             "Phos","Phos",
                             "RNA","RNA","RNA","Prot","Prot","Prot","Prot","Prot",
                             "RNA","RNA","RNA",
                             "RNA","RNA","Prot")
mergedtftargetsigdf$Row <- factor(as.character(1:34),levels = as.character(1:34))
mergedtftargetsigdf$Tissue <- factor(mergedtftargetsigdf$Tissue,levels = c("SKM-GN","HEART","KIDNEY","LIVER","LUNG","BAT","WAT-SC"))

mergedtftargetsigdf <- rbind(mergedtftargetsigdf,mergedtftargetsigdf)
mergedtftargetsigdf[c(35:68),"Percent.Significant"] <- mergedtftargetsigdf[c(1:34),"DEG.Enrichment"]
mergedtftargetsigdf$DEG.Frequency <- c(rep("TF.Targets",34),
                                       rep("Whole.Tissue",34))
mergedtftargetsigdf$DEG.Frequency <- factor(mergedtftargetsigdf$DEG.Frequency,levels = c("TF.Targets","Whole.Tissue"))

generaltable <- data.frame(start = c(0.5,5.5,6.5,7.5,9.5,12.5,13.5,14.5,16.5,20.5,23.5,28.5,33.5),
                           end = c(5.5,6.5,7.5,9.5,12.5,13.5,14.5,16.5,20.5,23.5,28.5,33.5,34.5),
                           Sig.Group = factor(c("RNA","Prot","Prot+Phos","Phos","RNA","Prot","Phos","RNA",
                                                "Phos","RNA","Prot","RNA","Prot"),
                                              levels = c("RNA","Prot","Prot+Phos","Phos")))


mergedtftargetsigdf$Significance <- "N"
mergedtftargetsigdf[c(2,4,5,7,9,14,22,27),"Significance"] <- "Y"


png(file = "Figure 5D.png",width = 12,height = 4.5,units = "in",res = 600)
ggplot(mergedtftargetsigdf,aes(x=Row,y=Percent.Significant,group=interaction(DEG.Frequency,Tissue))) + 
  geom_point(aes(color=Tissue),size = 3) + 
  geom_line(aes(linetype = DEG.Frequency,color=Tissue),size = 2) + 
  geom_rect(data = generaltable,aes(xmin = start,xmax = end,ymin = -Inf,ymax = 0,fill = Sig.Group),alpha = 0.5,inherit.aes = F) +
  geom_point(data = mergedtftargetsigdf[mergedtftargetsigdf$Significance == "Y", ], aes(Row, Percent.Significant + 0.05), shape = "*", size=10, color="black") +
  theme_classic() + 
  theme(legend.position = "top",
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.box="vertical",
        legend.margin=margin()) + 
  scale_x_discrete(labels=mergedtftargetsigdf$Transcription.Factor.Name) +
  scale_color_manual(values = ann_cols$Tissue[unique(mergedtftargetsigdf$Tissue)])+
  scale_linetype_manual(values = c("TF.Targets" = "solid", "Whole.Tissue" = "dotted"))+
  xlab("Transcription Factor")
dev.off()

pdf(file = "Figure 5D.pdf",width = 12,height = 4.5)
ggplot(mergedtftargetsigdf,aes(x=Row,y=Percent.Significant,group=interaction(DEG.Frequency,Tissue))) + 
  geom_point(aes(color=Tissue),size = 3) + 
  geom_line(aes(linetype = DEG.Frequency,color=Tissue),size = 2) + 
  geom_rect(data = generaltable,aes(xmin = start,xmax = end,ymin = -Inf,ymax = 0,fill = Sig.Group),alpha = 0.5,inherit.aes = F) +
  geom_point(data = mergedtftargetsigdf[mergedtftargetsigdf$Significance == "Y", ], aes(Row, Percent.Significant + 0.05), shape = "*", size=10, color="black") +
  theme_classic() + 
  theme(legend.position = "top",
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.box="vertical",
        legend.margin=margin()) + 
  scale_x_discrete(labels=mergedtftargetsigdf$Transcription.Factor.Name) +
  scale_color_manual(values = ann_cols$Tissue[unique(mergedtftargetsigdf$Tissue)])+
  scale_linetype_manual(values = c("TF.Targets" = "solid", "Whole.Tissue" = "dotted"))+
  xlab("Transcription Factor")
dev.off()

save.image("Figure5_S14_S15_S16_S17.RData")