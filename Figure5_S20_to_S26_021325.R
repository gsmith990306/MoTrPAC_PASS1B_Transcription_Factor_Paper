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
                 "Group" = c("control" = "grey50",
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

tfrnasigmattrim <- tfrnasigmat[apply(abs(tfrnasigmat),1,max) > 0.75,]
tfrnasigmatfinal <- rbind(tfrnasigmattrim[1:3,],
                          tfrnasigmat["SKM-GN_ENSRNOG00000022777",],
                          tfrnasigmat["SKM-GN_ENSRNOG00000019568",],
                          tfrnasigmat["SKM-GN_ENSRNOG00000021085",],
                          tfrnasigmat["HEART_ENSRNOG00000019568",],
                          tfrnasigmattrim[4:14,],
                          tfrnasigmat["LUNG_ENSRNOG00000004693",],
                          tfrnasigmattrim[15:18,],
                          tfrnasigmat["BAT_ENSRNOG00000009972",],
                          tfrnasigmattrim[19:31,],
                          tfrnasigmat["WAT-SC_ENSRNOG00000019568",],
                          tfrnasigmattrim[32:53,]
                          )
rownames(tfrnasigmatfinal)[4:7] <- c("SKM-GN_ENSRNOG00000022777","SKM-GN_ENSRNOG00000019568","SKM-GN_ENSRNOG00000021085","HEART_ENSRNOG00000019568")
rownames(tfrnasigmatfinal)[19] <- "LUNG_ENSRNOG00000004693"
rownames(tfrnasigmatfinal)[24] <- "BAT_ENSRNOG00000009972"
rownames(tfrnasigmatfinal)[38] <- "WAT-SC_ENSRNOG00000019568"

png("Figure 5A_021325.png",width = 6,height = 8,units = "in",res = 600)
pheatmap(tfrnasigmatfinal,cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),annotation_col = tfsigmeta[,c("Group","Sex")],annotation_row = tfrnasigmeta,annotation_colors = ann_cols,show_colnames = F,cluster_rows = F,labels_row = enstosym[gsub(".*_","",rownames(tfrnasigmatfinal)),"Symbol"],display_numbers = T,number_color = "black",border_color = NA)
dev.off()

pdf("Figure 5A_021325.pdf",width = 6,height = 8)
pheatmap(tfrnasigmatfinal,cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),annotation_col = tfsigmeta[,c("Group","Sex")],annotation_row = tfrnasigmeta,annotation_colors = ann_cols,show_colnames = F,cluster_rows = F,labels_row = enstosym[gsub(".*_","",rownames(tfrnasigmatfinal)),"Symbol"],display_numbers = T,number_color = "black",border_color = NA)
dev.off()

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

tfprosigmattrimlist <- rownames(tfprosigmat)[apply(abs(tfprosigmat),1,max) > 0.25]
tfprosigmattrimlist <- c(tfprosigmattrimlist[1:13],"LUNG_NP_059021.1",tfprosigmattrimlist[14:35])
tfprosigmattrimlistids <- gsub(".*NP","NP",gsub(".*XP","XP",tfprosigmattrimlist))
rownamedf <- data.frame(row.names = tfprosigmattrimlistids,
                        "ID_list" = tfprosigmattrimlistids)
rownamedf$proname <- ""
for(i in 1:dim(rownamedf)[1]){
  rownamedf[i,"proname"] <- toupper(tfproanno[rownames(which(tfproanno == tfprosigmattrimlistids[i],arr.ind = TRUE))[1],"Gene.Name"])
}

png("Figure 5B_021325.png",width = 6,height = 8,units = "in",res = 600)
pheatmap(tfprosigmat[tfprosigmattrimlist,],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),annotation_col = tfsigmeta[,c("Group","Sex")],annotation_row = tfprosigmeta,annotation_colors = ann_cols,show_colnames = F,cluster_rows = F,labels_row = rownamedf$proname,display_numbers = T,number_color = "black",border_color = NA)
dev.off()

pdf("Figure 5B_021325.pdf",width = 6,height = 8)
pheatmap(tfprosigmat[tfprosigmattrimlist,],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),annotation_col = tfsigmeta[,c("Group","Sex")],annotation_row = tfprosigmeta,annotation_colors = ann_cols,show_colnames = F,cluster_rows = F,labels_row = rownamedf$proname,display_numbers = T,number_color = "black",border_color = NA)
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

tfprophsigmattrim <- tfprophsigmat[apply(abs(tfprophsigmat),1,max) > 0.4,]
tfprophannoselectcombotrim <- c(paste(toupper(tfprophannogastroselect$Symbol),gsub(".*_","",tfprophannogastroselect$Phospho),sep = "_"),
                            paste(toupper(tfprophannoheartselect$Symbol),gsub(".*_","",tfprophannoheartselect$Phospho),sep = "_"),
                            paste(toupper(tfprophannokidneyselect$Symbol),gsub(".*_","",tfprophannokidneyselect$Phospho),sep = "_"),
                            paste(toupper(tfprophannoliverselect$Symbol),gsub(".*_","",tfprophannoliverselect$Phospho),sep = "_"),
                            paste(toupper(tfprophannolungselect$Symbol),gsub(".*_","",tfprophannolungselect$Phospho),sep = "_"),
                            paste(toupper(tfprophannowhiteselect$Symbol),gsub(".*_","",tfprophannowhiteselect$Phospho),sep = "_"))[apply(abs(tfprophsigmat),1,max) > 0.4]
tfprophsigmatfinal <- rbind(tfprophsigmattrim[1:5,],
                            tfprophsigmat["HEART_XP_008757763.1_S465sS480s",],
                            tfprophsigmattrim[6:9,],
                            tfprophsigmat["KIDNEY_NP_001101585.1_T53t",],
                            tfprophsigmattrim[10:42,])
rownames(tfprophsigmatfinal)[6] <- "HEART_XP_008757763.1_S465sS480s"
rownames(tfprophsigmatfinal)[11] <- "KIDNEY_NP_001101585.1_T53t"

tfprophannoselectcombofinal <- c(tfprophannoselectcombotrim[1:5],"MEF2A_S465sS480s",tfprophannoselectcombotrim[6:9],
                                 "ATF7_T53t",tfprophannoselectcombotrim[10:42])

png("Figure 5C_021325.png",width = 6,height = 8,units = "in",res = 600)
pheatmap(tfprophsigmatfinal,cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),annotation_col = tfsigmeta[,c("Group","Sex")],annotation_row = tfprophsigmeta,annotation_colors = ann_cols,show_colnames = F,cluster_rows = F,
         labels_row = tfprophannoselectcombofinal,display_numbers = T,number_color = "black",border_color = NA)
dev.off()

pdf("Figure 5C_021325.pdf",width = 6,height = 8)
pheatmap(tfprophsigmatfinal,cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),annotation_col = tfsigmeta[,c("Group","Sex")],annotation_row = tfprophsigmeta,annotation_colors = ann_cols,show_colnames = F,cluster_rows = F,
         labels_row = tfprophannoselectcombofinal,display_numbers = T,number_color = "black",border_color = NA)
dev.off()

save.image("Figure5A_5B_5C_021325.RData")

#####
# Supplemental Figure S22
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

png(file = "Supplemental Figure S22A_021325.png",width = 9,height = 4.5,units = "in",res = 600)
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

pdf(file = "Supplemental Figure S22A_021325.pdf",width = 9,height = 4.5)
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

png(file = "Supplemental Figure S22B_021325.png",width = 9,height = 4.5,units = "in",res = 600)
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

pdf(file = "Supplemental Figure S22B_021325.pdf",width = 9,height = 4.5)
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

png(file = "Supplemental Figure S22C_021325.png",width = 9,height = 4.5,units = "in",res = 600)
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

pdf(file = "Supplemental Figure S22C_021325.pdf",width = 9,height = 4.5)
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

png(file = "Supplemental Figure S22D_021325.png",width = 9,height = 4.5,units = "in",res = 600)
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

pdf(file = "Supplemental Figure S22D_021325.pdf",width = 9,height = 4.5)
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

png(file = "Supplemental Figure S22E_021325.png",width = 9,height = 4.5,units = "in",res = 600)
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

pdf(file = "Supplemental Figure S22E_021325.pdf",width = 9,height = 4.5)
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

png(file = "Supplemental Figure S22F_021325.png",width = 9,height = 4.5,units = "in",res = 600)
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

pdf(file = "Supplemental Figure S22F_021325.pdf",width = 9,height = 4.5)
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
# Supplemental Figure S25
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

png(file = "Supplemental Figure S25A_021325.png",width = 7,height = 4,units = "in",res = 600)
pheatmap(gastrol2fcmat[intersect(ourtftargetpromproxgenes,gastrornasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,gastrornasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S25A_021325.pdf",width = 7,height = 4)
pheatmap(gastrol2fcmat[intersect(ourtftargetpromproxgenes,gastrornasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,gastrornasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
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

png(file = "Supplemental Figure S25C_021325.png",width = 7,height = 4,units = "in",res = 600)
pheatmap(gastrol2fcmat[intersect(ourtftargetpromproxgenes,gastrornasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,gastrornasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S25C_021325.pdf",width = 7,height = 4)
pheatmap(gastrol2fcmat[intersect(ourtftargetpromproxgenes,gastrornasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,gastrornasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
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

png(file = "Supplemental Figure S25B_021325.png",width = 7,height = 4,units = "in",res = 600)
pheatmap(lungl2fcmat[intersect(ourtftargetpromproxgenes,lungrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,lungrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S25B_021325.pdf",width = 7,height = 4)
pheatmap(lungl2fcmat[intersect(ourtftargetpromproxgenes,lungrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,lungrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
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
                          "TF.L2FC" = whiteprol2fc[tfproanno[ourtf,"WhiteAd.Pro.ID"],])
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

png(file = "Supplemental Figure S25D_021325.png",width = 7,height = 4,units = "in",res = 600)
pheatmap(whitel2fcmat[intersect(ourtftargetpromproxgenes,whiternasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,whiternasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S25D_021325.pdf",width = 7,height = 4)
pheatmap(whitel2fcmat[intersect(ourtftargetpromproxgenes,whiternasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,whiternasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
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

png(file = "Supplemental Figure S25E_021325.png",width = 7,height = 4,units = "in",res = 600)
pheatmap(gastrol2fcmat[intersect(ourtftargetpromproxgenes,gastrornasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,gastrornasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S25E_021325.pdf",width = 7,height = 4)
pheatmap(gastrol2fcmat[intersect(ourtftargetpromproxgenes,gastrornasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,gastrornasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
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

png(file = "Supplemental Figure S25F_021325.png",width = 7,height = 4,units = "in",res = 600)
pheatmap(heartl2fcmat[intersect(ourtftargetpromproxgenes,heartrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,heartrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S25F_021325.pdf",width = 7,height = 4)
pheatmap(heartl2fcmat[intersect(ourtftargetpromproxgenes,heartrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,heartrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
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

png(file = "Supplemental Figure S25G_021325.png",width = 7,height = 4,units = "in",res = 600)
pheatmap(kidneyl2fcmat[intersect(ourtftargetpromproxgenes,kidneyrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,kidneyrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S25G_021325.pdf",width = 7,height = 4)
pheatmap(kidneyl2fcmat[intersect(ourtftargetpromproxgenes,kidneyrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,kidneyrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
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

png(file = "Supplemental Figure S25H_021325.png",width = 7,height = 4,units = "in",res = 600)
pheatmap(liverl2fcmat[intersect(ourtftargetpromproxgenes,liverrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,liverrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S25H_021325.pdf",width = 7,height = 4)
pheatmap(liverl2fcmat[intersect(ourtftargetpromproxgenes,liverrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,liverrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
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



pdf(file = "Supplemental Figure S23A_021325.pdf",width = 12,height = 4.5)
ggplot(gastrornatargetsigdf[!gastrornatargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(gastrornasig)/dim(gastrol2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),axis.title = element_text(size = 13),legend.text = element_text(size = 12))
dev.off()

pdf(file = "Supplemental Figure S23B_021325.pdf",width = 12,height = 4.5)
ggplot(heartrnatargetsigdf[!heartrnatargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(heartrnasig)/dim(heartl2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),axis.title = element_text(size = 13),legend.text = element_text(size = 12))
dev.off()

pdf(file = "Supplemental Figure S23C_021325.pdf",width = 12,height = 4.5)
ggplot(hippornatargetsigdf[!hippornatargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(hippornasig)/dim(hippol2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),axis.title = element_text(size = 13),legend.text = element_text(size = 12))
dev.off()

pdf(file = "Supplemental Figure S23D_021325.pdf",width = 12,height = 4.5)
ggplot(kidneyrnatargetsigdf[!kidneyrnatargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(kidneyrnasig)/dim(kidneyl2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),axis.title = element_text(size = 13),legend.text = element_text(size = 12))
dev.off()

pdf(file = "Supplemental Figure S23E_021325.pdf",width = 12,height = 4.5)
ggplot(liverrnatargetsigdf[!liverrnatargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(liverrnasig)/dim(liverl2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),axis.title = element_text(size = 13),legend.text = element_text(size = 12))
dev.off()

pdf(file = "Supplemental Figure S23F_021325.pdf",width = 12,height = 4.5)
ggplot(lungrnatargetsigdf[!lungrnatargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(lungrnasig)/dim(lungl2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),axis.title = element_text(size = 13),legend.text = element_text(size = 12))
dev.off()

pdf(file = "Supplemental Figure S23G_021325.pdf",width = 12,height = 4.5)
ggplot(brownrnatargetsigdf[!brownrnatargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(brownrnasig)/dim(brownl2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),axis.title = element_text(size = 13),legend.text = element_text(size = 12))
dev.off()

pdf(file = "Supplemental Figure S23H_021325.pdf",width = 12,height = 4.5)
ggplot(whiternatargetsigdf[!whiternatargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(whiternasig)/dim(whitel2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),axis.title = element_text(size = 13),legend.text = element_text(size = 12))
dev.off()

png(file = "Supplemental Figure S23A_021325.png",width = 12,height = 4.5,units = "in",res = 600)
ggplot(gastrornatargetsigdf[!gastrornatargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(gastrornasig)/dim(gastrol2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),axis.title = element_text(size = 13),legend.text = element_text(size = 12))
dev.off()

png(file = "Supplemental Figure S23B_021325.png",width = 12,height = 4.5,units = "in",res = 600)
ggplot(heartrnatargetsigdf[!heartrnatargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(heartrnasig)/dim(heartl2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),axis.title = element_text(size = 13),legend.text = element_text(size = 12))
dev.off()

png(file = "Supplemental Figure S23C_021325.png",width = 12,height = 4.5,units = "in",res = 600)
ggplot(hippornatargetsigdf[!hippornatargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(hippornasig)/dim(hippol2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),axis.title = element_text(size = 13),legend.text = element_text(size = 12))
dev.off()

png(file = "Supplemental Figure S23D_021325.png",width = 12,height = 4.5,units = "in",res = 600)
ggplot(kidneyrnatargetsigdf[!kidneyrnatargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(kidneyrnasig)/dim(kidneyl2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),axis.title = element_text(size = 13),legend.text = element_text(size = 12))
dev.off()

png(file = "Supplemental Figure S23E_021325.png",width = 12,height = 4.5,units = "in",res = 600)
ggplot(liverrnatargetsigdf[!liverrnatargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(liverrnasig)/dim(liverl2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),axis.title = element_text(size = 13),legend.text = element_text(size = 12))
dev.off()

png(file = "Supplemental Figure S23F_021325.png",width = 12,height = 4.5,units = "in",res = 600)
ggplot(lungrnatargetsigdf[!lungrnatargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(lungrnasig)/dim(lungl2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),axis.title = element_text(size = 13),legend.text = element_text(size = 12))
dev.off()

png(file = "Supplemental Figure S23G_021325.png",width = 12,height = 4.5,units = "in",res = 600)
ggplot(brownrnatargetsigdf[!brownrnatargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(brownrnasig)/dim(brownl2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),axis.title = element_text(size = 13),legend.text = element_text(size = 12))
dev.off()

png(file = "Supplemental Figure S23H_021325.png",width = 12,height = 4.5,units = "in",res = 600)
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

png(file = "Supplemental Figure S24B_021325.png",width = 7,height = 4,units = "in",res = 600)
pheatmap(gastrol2fcmat[intersect(ourtftargetpromproxgenes,gastrornasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,gastrornasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S24B_021325.pdf",width = 7,height = 4)
pheatmap(gastrol2fcmat[intersect(ourtftargetpromproxgenes,gastrornasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,gastrornasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
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

png(file = "Supplemental Figure S24D_021325.png",width = 7,height = 4,units = "in",res = 600)
pheatmap(gastrol2fcmat[intersect(ourtftargetpromproxgenes,gastrornasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,gastrornasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S24D_021325.pdf",width = 7,height = 4)
pheatmap(gastrol2fcmat[intersect(ourtftargetpromproxgenes,gastrornasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,gastrornasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
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

png(file = "Supplemental Figure S24A_021325.png",width = 7,height = 6,units = "in",res = 600)
pheatmap(gastrol2fcmat[intersect(ourtftargetpromproxgenes,gastrornasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,gastrornasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S24A_021325.pdf",width = 7,height = 6)
pheatmap(gastrol2fcmat[intersect(ourtftargetpromproxgenes,gastrornasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,gastrornasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
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
                          "TF.L2FC" = heartl2fcmat[tfanno[ourtf,"Ensembl"],])
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

png(file = "Supplemental Figure S24E_021325.png",width = 7,height = 4,units = "in",res = 600)
pheatmap(heartl2fcmat[intersect(ourtftargetpromproxgenes,heartrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,heartrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S24E_021325.pdf",width = 7,height = 4)
pheatmap(heartl2fcmat[intersect(ourtftargetpromproxgenes,heartrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,heartrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
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
                          "TF.L2FC" = kidneyl2fcmat[tfanno[ourtf,"Ensembl"],])
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

png(file = "Supplemental Figure S24C_021325.png",width = 7,height = 4,units = "in",res = 600)
pheatmap(kidneyl2fcmat[intersect(ourtftargetpromproxgenes,kidneyrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,kidneyrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S24C_021325.pdf",width = 7,height = 4)
pheatmap(kidneyl2fcmat[intersect(ourtftargetpromproxgenes,kidneyrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,kidneyrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
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
                          "TF.L2FC" = whitel2fcmat[tfanno[ourtf,"Ensembl"],])
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

png(file = "Supplemental Figure S24F_021325.png",width = 7,height = 4,units = "in",res = 600)
pheatmap(whitel2fcmat[intersect(ourtftargetpromproxgenes,whiternasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,whiternasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S24F_021325.pdf",width = 7,height = 4)
pheatmap(whitel2fcmat[intersect(ourtftargetpromproxgenes,whiternasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,whiternasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
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
                          "TF.L2FC" = lungl2fcmat[tfanno[ourtf,"Ensembl"],])
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

png(file = "Supplemental Figure S24G_021325.png",width = 7,height = 4,units = "in",res = 600)
pheatmap(lungl2fcmat[intersect(ourtftargetpromproxgenes,lungrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,lungrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S24G_021325.pdf",width = 7,height = 4)
pheatmap(lungl2fcmat[intersect(ourtftargetpromproxgenes,lungrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,lungrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
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
                          "TF.L2FC" = brownl2fcmat[tfanno[ourtf,"Ensembl"],])
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

png(file = "Supplemental Figure S24H_021325.png",width = 7,height = 4,units = "in",res = 600)
pheatmap(brownl2fcmat[intersect(ourtftargetpromproxgenes,brownrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,brownrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S24H_021325.pdf",width = 7,height = 4)
pheatmap(brownl2fcmat[intersect(ourtftargetpromproxgenes,brownrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,brownrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
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


png(file = "Figure 5D_021325.png",width = 12,height = 4.5,units = "in",res = 600)
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
  scale_color_manual(values = ann_cols$Tissue)+
  scale_linetype_manual(values = c("TF.Targets" = "solid", "Whole.Tissue" = "dotted"))+
  xlab("Transcription Factor")
dev.off()

pdf(file = "Figure 5D_021325.pdf",width = 12,height = 4.5)
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
  scale_color_manual(values = ann_cols$Tissue)+
  scale_linetype_manual(values = c("TF.Targets" = "solid", "Whole.Tissue" = "dotted"))+
  xlab("Transcription Factor")
dev.off()

save.image("Figure5A_to_5D_S22_to_S25_021325.RData")


#####
# We are going to investigate the general relationship between TFs at the RNAseq, Proteomic and Phosphoproteomic levels
####

load("rnal2fcmat.RData")
load("new_rnanormmatrices.RData")
load("proprnormmatrices.RData")
load("prophnormmatrices.RData")
pass1bphenodata <- readRDS("pass1bphenodata.rds")


####
# gastro

tfpro_crossomeannogastroselect <- tfproanno[nchar(tfproanno$Gastro.Pro.ID) > 0,c("Gene.Name","Gastro.Pro.ID","Ensembl")]
tfpro_crossomeannogastroselect <- tfpro_crossomeannogastroselect[!duplicated(tfpro_crossomeannogastroselect$Gastro.Pro.ID),]

tfpro_crossomeannogastroselect <- tfpro_crossomeannogastroselect[tfpro_crossomeannogastroselect$Ensembl %in% rownames(gastrornanorm),]

tfgastroprol2fc <- gastroprol2fc[tfpro_crossomeannogastroselect$Gastro.Pro.ID,]

tfgastropro_prol2fcmat <- gastroprol2fc[tfpro_crossomeannogastroselect$Gastro.Pro.ID,]
tfgastropro_rnal2fcmat <- gastrol2fcmat[tfpro_crossomeannogastroselect$Ensembl,]

tfgastropro_rnanorm <- gastrornanorm[rownames(tfgastropro_rnal2fcmat),]
for(i in 1:dim(tfgastropro_rnanorm)[2]){
  ourid <- colnames(tfgastropro_rnanorm)[i]
  colnames(tfgastropro_rnanorm)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0001"]
}
tfgastropro_pronorm <- gastroproprnorm[rownames(tfgastropro_prol2fcmat),]

gastroprosamplesubset <- intersect(colnames(tfgastropro_rnanorm),colnames(tfgastropro_pronorm))
tfgastropro_rnanorm <- tfgastropro_rnanorm[,gastroprosamplesubset]
tfgastropro_pronorm <- tfgastropro_pronorm[,gastroprosamplesubset]

tfgastropro_crossomecor <- matrix(0L,nrow = dim(tfgastropro_pronorm)[1],ncol = 1)
rownames(tfgastropro_crossomecor) <- rownames(tfgastropro_pronorm)
colnames(tfgastropro_crossomecor) <- c("RNA-Pro")
tfgastropro_crossomecortest <- matrix(0L,nrow = dim(tfgastropro_pronorm)[1],ncol = 1)
rownames(tfgastropro_crossomecortest) <- rownames(tfgastropro_pronorm)
colnames(tfgastropro_crossomecortest) <- c("RNA-Pro")
for(i in 1:dim(tfgastropro_crossomecor)[1]){
  ourtest <- cor.test(t(tfgastropro_rnanorm[i,]),t(tfgastropro_pronorm[i,]))
  tfgastropro_crossomecor[i,"RNA-Pro"] <- ourtest$estimate
  tfgastropro_crossomecortest[i,"RNA-Pro"] <- ourtest$p.value
}

#png("gastro_tfpro_crossomecor.png",width = 2.5,height = 12,units = "in",res = 600)
#pheatmap(tfgastropro_crossomecor,cluster_cols = F,angle_col = 0,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_row = tfpro_crossomeannogastroselect$Gene.Name,display_numbers = T,number_color = "black")
#dev.off()

tfgastroproph <- prot_ph$training_dea[prot_ph$training_dea$pr_feature_ID %in% tfproanno$Gastro.Pro.ID & prot_ph$training_dea$tissue_abbreviation %in% "SKM-GN","feature_ID"]
tfgastroproph_prophl2fcmat <- gastroprophl2fc[tfgastroproph,]
tfgastroproph_prol2fcmat <- matrix(0L,nrow = dim(tfgastroproph_prophl2fcmat)[1],ncol = 8)
colnames(tfgastroproph_prol2fcmat) <- colnames(tfgastroproph_prophl2fcmat)
rownames(tfgastroproph_prol2fcmat) <- rownames(tfgastroproph_prophl2fcmat)
for(i in 1:dim(tfgastroproph_prol2fcmat)[1]){
  ourid <- rownames(tfgastroproph_prophl2fcmat)[i]
  rownames(tfgastroproph_prol2fcmat)[i] <- prot_ph$training_dea[prot_ph$training_dea$feature_ID %in% ourid,"pr_feature_ID"][1]
}
tfgastroproph_prol2fcmat <- gastroprol2fc[rownames(tfgastroproph_prol2fcmat),]

tfgastroproph_rnal2fcmat <- matrix(0L,nrow = dim(tfgastroproph_prophl2fcmat)[1],ncol = 8)
colnames(tfgastroproph_rnal2fcmat) <- colnames(tfgastroproph_prophl2fcmat)
rownames(tfgastroproph_rnal2fcmat) <- rownames(tfgastroproph_prol2fcmat)
for(i in 1:dim(tfgastroproph_rnal2fcmat)[1]){
  ourid <- rownames(tfgastroproph_rnal2fcmat)[i]
  rownames(tfgastroproph_rnal2fcmat)[i] <- tfproanno[grep(ourid,tfproanno$Gastro.Pro.ID),"Ensembl"][1]
}

tfgastroproph_prol2fcmat <- tfgastroproph_prol2fcmat[rownames(tfgastroproph_rnal2fcmat) %in% rownames(gastrol2fcmat),]
tfgastroproph_prophl2fcmat <- tfgastroproph_prophl2fcmat[rownames(tfgastroproph_rnal2fcmat) %in% rownames(gastrol2fcmat),]
tfgastroproph_rnal2fcmat <- tfgastroproph_rnal2fcmat[rownames(tfgastroproph_rnal2fcmat) %in% rownames(gastrol2fcmat),]

tfgastroproph_rnal2fcmat <- gastrol2fcmat[rownames(tfgastroproph_rnal2fcmat),]

tfgastroproph_rnanorm <- gastrornanorm[rownames(tfgastroproph_rnal2fcmat),]
for(i in 1:dim(tfgastroproph_rnanorm)[2]){
  ourid <- colnames(tfgastroproph_rnanorm)[i]
  colnames(tfgastroproph_rnanorm)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0001"]
}

tfgastroproph_pronorm <- gastroproprnorm[rownames(tfgastroproph_prol2fcmat),]
tfgastroproph_prophnorm <- gastroprophnorm[rownames(tfgastroproph_prophl2fcmat),]
for(i in 1:dim(tfgastroproph_prophnorm)[2]){
  ourid <- colnames(tfgastroproph_prophnorm)[i]
  colnames(tfgastroproph_prophnorm)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0001"]
}

gastroprophsamplesubset <- Reduce(intersect,list(colnames(tfgastroproph_rnanorm),colnames(tfgastroproph_pronorm),colnames(tfgastroproph_prophnorm)))
tfgastroproph_rnanorm <- tfgastroproph_rnanorm[,gastroprophsamplesubset]
tfgastroproph_pronorm <- tfgastroproph_pronorm[,gastroprophsamplesubset]
tfgastroproph_prophnorm <- tfgastroproph_prophnorm[,gastroprophsamplesubset]

tfgastroproph_crossomecor <- matrix(0L,nrow = dim(tfgastroproph_prophnorm)[1],ncol = 3)
rownames(tfgastroproph_crossomecor) <- rownames(tfgastroproph_prophnorm)
colnames(tfgastroproph_crossomecor) <- c("RNA-Pro","RNA-Phos","Pro-Phos")
tfgastroproph_crossomecortest <- matrix(0L,nrow = dim(tfgastroproph_prophnorm)[1],ncol = 3)
rownames(tfgastroproph_crossomecortest) <- rownames(tfgastroproph_prophnorm)
colnames(tfgastroproph_crossomecortest) <- c("RNA-Pro","RNA-Phos","Pro-Phos")

for(i in 1:dim(tfgastroproph_crossomecor)[1]){
  ourtest <- cor.test(t(tfgastroproph_rnanorm[i,]),t(tfgastroproph_pronorm[i,]))
  tfgastroproph_crossomecor[i,"RNA-Pro"] <- ourtest$estimate
  tfgastroproph_crossomecortest[i,"RNA-Pro"] <- ourtest$p.value
  ourtest <- cor.test(t(tfgastroproph_rnanorm[i,]),t(tfgastroproph_prophnorm[i,]))
  tfgastroproph_crossomecor[i,"RNA-Phos"] <- ourtest$estimate
  tfgastroproph_crossomecortest[i,"RNA-Phos"] <- ourtest$p.value
  ourtest <- cor.test(t(tfgastroproph_pronorm[i,]),t(tfgastroproph_prophnorm[i,]))
  tfgastroproph_crossomecor[i,"Pro-Phos"] <- ourtest$estimate
  tfgastroproph_crossomecortest[i,"Pro-Phos"] <- ourtest$p.value
}

tfproph_crossomeannogastroselect <- data.frame(row.names = rownames(tfgastroproph_crossomecor),"Phospho" = rownames(tfgastroproph_crossomecor),"Protein_ID" = sapply(strsplit(rownames(tfgastroproph_crossomecor), "_"), function(x) {
  paste(x[1],x[2],sep = "_")
}))
tfproph_crossomeannogastroselect$Symbol = ""
for(i in 1:length(rownames(tfgastroproph_crossomecor))){
  ourpro <- tfproph_crossomeannogastroselect[i,"Protein_ID"]
  tfproph_crossomeannogastroselect[i,"Symbol"] <- toupper(tfproanno[tfproanno$Gastro.Pro.ID %in% ourpro,"Gene.Name"][1])
}

tfproph_crossomeannogastroselect$Sym_ID <- paste(toupper(tfproph_crossomeannogastroselect$Symbol),gsub(".*_","",tfproph_crossomeannogastroselect$Phospho),sep = "_")

#png("gastro_tfproph_crossomecor.png",width = 7.5,height = 12,units = "in",res = 600)
#pheatmap(tfgastroproph_crossomecor,cluster_cols = F,angle_col = 0,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_row = tfproph_crossomeannogastroselect$Sym_ID,display_numbers = T,number_color = "black")
#dev.off()


####
# heart

tfpro_crossomeannoheartselect <- tfproanno[nchar(tfproanno$Heart.Pro.ID) > 0,c("Gene.Name","Heart.Pro.ID","Ensembl")]
tfpro_crossomeannoheartselect <- tfpro_crossomeannoheartselect[!duplicated(tfpro_crossomeannoheartselect$Heart.Pro.ID),]

tfpro_crossomeannoheartselect <- tfpro_crossomeannoheartselect[tfpro_crossomeannoheartselect$Ensembl %in% rownames(heartrnanorm),]

tfheartprol2fc <- heartprol2fc[tfpro_crossomeannoheartselect$Heart.Pro.ID,]

tfheartpro_prol2fcmat <- heartprol2fc[tfpro_crossomeannoheartselect$Heart.Pro.ID,]
tfheartpro_rnal2fcmat <- heartl2fcmat[tfpro_crossomeannoheartselect$Ensembl,]

tfheartpro_rnanorm <- heartrnanorm[rownames(tfheartpro_rnal2fcmat),]
for(i in 1:dim(tfheartpro_rnanorm)[2]){
  ourid <- colnames(tfheartpro_rnanorm)[i]
  colnames(tfheartpro_rnanorm)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0001"]
}
tfheartpro_pronorm <- heartproprnorm[rownames(tfheartpro_prol2fcmat),]

heartprosamplesubset <- intersect(colnames(tfheartpro_rnanorm),colnames(tfheartpro_pronorm))
tfheartpro_rnanorm <- tfheartpro_rnanorm[,heartprosamplesubset]
tfheartpro_pronorm <- tfheartpro_pronorm[,heartprosamplesubset]

tfheartpro_crossomecor <- matrix(0L,nrow = dim(tfheartpro_pronorm)[1],ncol = 1)
rownames(tfheartpro_crossomecor) <- rownames(tfheartpro_pronorm)
colnames(tfheartpro_crossomecor) <- c("RNA-Pro")
tfheartpro_crossomecortest <- matrix(0L,nrow = dim(tfheartpro_pronorm)[1],ncol = 1)
rownames(tfheartpro_crossomecortest) <- rownames(tfheartpro_pronorm)
colnames(tfheartpro_crossomecortest) <- c("RNA-Pro")
for(i in 1:dim(tfheartpro_crossomecor)[1]){
  ourtest <- cor.test(t(tfheartpro_rnanorm[i,]),t(tfheartpro_pronorm[i,]))
  tfheartpro_crossomecor[i,"RNA-Pro"] <- ourtest$estimate
  tfheartpro_crossomecortest[i,"RNA-Pro"] <- ourtest$p.value
}

#png("heart_tfpro_crossomecor.png",width = 2.5,height = 12,units = "in",res = 600)
#pheatmap(tfheartpro_crossomecor,cluster_cols = F,angle_col = 0,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_row = tfpro_crossomeannoheartselect$Gene.Name,display_numbers = T,number_color = "black")
#dev.off()

tfheartproph <- prot_ph$training_dea[prot_ph$training_dea$pr_feature_ID %in% tfproanno$Heart.Pro.ID & prot_ph$training_dea$tissue_abbreviation %in% "HEART","feature_ID"]
tfheartproph_prophl2fcmat <- heartprophl2fc[tfheartproph,]
tfheartproph_prol2fcmat <- matrix(0L,nrow = dim(tfheartproph_prophl2fcmat)[1],ncol = 8)
colnames(tfheartproph_prol2fcmat) <- colnames(tfheartproph_prophl2fcmat)
rownames(tfheartproph_prol2fcmat) <- rownames(tfheartproph_prophl2fcmat)
for(i in 1:dim(tfheartproph_prol2fcmat)[1]){
  ourid <- rownames(tfheartproph_prophl2fcmat)[i]
  rownames(tfheartproph_prol2fcmat)[i] <- prot_ph$training_dea[prot_ph$training_dea$feature_ID %in% ourid,"pr_feature_ID"][1]
}
tfheartproph_prol2fcmat <- heartprol2fc[rownames(tfheartproph_prol2fcmat),]

tfheartproph_rnal2fcmat <- matrix(0L,nrow = dim(tfheartproph_prophl2fcmat)[1],ncol = 8)
colnames(tfheartproph_rnal2fcmat) <- colnames(tfheartproph_prophl2fcmat)
rownames(tfheartproph_rnal2fcmat) <- rownames(tfheartproph_prol2fcmat)
for(i in 1:dim(tfheartproph_rnal2fcmat)[1]){
  ourid <- rownames(tfheartproph_rnal2fcmat)[i]
  rownames(tfheartproph_rnal2fcmat)[i] <- tfproanno[grep(ourid,tfproanno$Heart.Pro.ID),"Ensembl"][1]
}

tfheartproph_prol2fcmat <- tfheartproph_prol2fcmat[rownames(tfheartproph_rnal2fcmat) %in% rownames(heartl2fcmat),]
tfheartproph_prophl2fcmat <- tfheartproph_prophl2fcmat[rownames(tfheartproph_rnal2fcmat) %in% rownames(heartl2fcmat),]
tfheartproph_rnal2fcmat <- tfheartproph_rnal2fcmat[rownames(tfheartproph_rnal2fcmat) %in% rownames(heartl2fcmat),]

tfheartproph_rnal2fcmat <- heartl2fcmat[rownames(tfheartproph_rnal2fcmat),]

tfheartproph_rnanorm <- heartrnanorm[rownames(tfheartproph_rnal2fcmat),]
for(i in 1:dim(tfheartproph_rnanorm)[2]){
  ourid <- colnames(tfheartproph_rnanorm)[i]
  colnames(tfheartproph_rnanorm)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0001"]
}

tfheartproph_pronorm <- heartproprnorm[rownames(tfheartproph_prol2fcmat),]
tfheartproph_prophnorm <- heartprophnorm[rownames(tfheartproph_prophl2fcmat),]
for(i in 1:dim(tfheartproph_prophnorm)[2]){
  ourid <- colnames(tfheartproph_prophnorm)[i]
  colnames(tfheartproph_prophnorm)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0001"]
}

heartprophsamplesubset <- Reduce(intersect,list(colnames(tfheartproph_rnanorm),colnames(tfheartproph_pronorm),colnames(tfheartproph_prophnorm)))
tfheartproph_rnanorm <- tfheartproph_rnanorm[,heartprophsamplesubset]
tfheartproph_pronorm <- tfheartproph_pronorm[,heartprophsamplesubset]
tfheartproph_prophnorm <- tfheartproph_prophnorm[,heartprophsamplesubset]

tfheartproph_crossomecor <- matrix(0L,nrow = dim(tfheartproph_prophnorm)[1],ncol = 3)
rownames(tfheartproph_crossomecor) <- rownames(tfheartproph_prophnorm)
colnames(tfheartproph_crossomecor) <- c("RNA-Pro","RNA-Phos","Pro-Phos")
tfheartproph_crossomecortest <- matrix(0L,nrow = dim(tfheartproph_prophnorm)[1],ncol = 3)
rownames(tfheartproph_crossomecortest) <- rownames(tfheartproph_prophnorm)
colnames(tfheartproph_crossomecortest) <- c("RNA-Pro","RNA-Phos","Pro-Phos")

for(i in 1:dim(tfheartproph_crossomecor)[1]){
  ourtest <- cor.test(t(tfheartproph_rnanorm[i,]),t(tfheartproph_pronorm[i,]))
  tfheartproph_crossomecor[i,"RNA-Pro"] <- ourtest$estimate
  tfheartproph_crossomecortest[i,"RNA-Pro"] <- ourtest$p.value
  ourtest <- cor.test(t(tfheartproph_rnanorm[i,]),t(tfheartproph_prophnorm[i,]))
  tfheartproph_crossomecor[i,"RNA-Phos"] <- ourtest$estimate
  tfheartproph_crossomecortest[i,"RNA-Phos"] <- ourtest$p.value
  ourtest <- cor.test(t(tfheartproph_pronorm[i,]),t(tfheartproph_prophnorm[i,]))
  tfheartproph_crossomecor[i,"Pro-Phos"] <- ourtest$estimate
  tfheartproph_crossomecortest[i,"Pro-Phos"] <- ourtest$p.value
}

tfproph_crossomeannoheartselect <- data.frame(row.names = rownames(tfheartproph_crossomecor),"Phospho" = rownames(tfheartproph_crossomecor),"Protein_ID" = sapply(strsplit(rownames(tfheartproph_crossomecor), "_"), function(x) {
  paste(x[1],x[2],sep = "_")
}))
tfproph_crossomeannoheartselect$Symbol = ""
for(i in 1:length(rownames(tfheartproph_crossomecor))){
  ourpro <- tfproph_crossomeannoheartselect[i,"Protein_ID"]
  tfproph_crossomeannoheartselect[i,"Symbol"] <- toupper(tfproanno[tfproanno$Heart.Pro.ID %in% ourpro,"Gene.Name"][1])
}

tfproph_crossomeannoheartselect$Sym_ID <- paste(toupper(tfproph_crossomeannoheartselect$Symbol),gsub(".*_","",tfproph_crossomeannoheartselect$Phospho),sep = "_")

#png("heart_tfproph_crossomecor.png",width = 7.5,height = 12,units = "in",res = 600)
#pheatmap(tfheartproph_crossomecor,cluster_cols = F,angle_col = 0,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_row = tfproph_crossomeannoheartselect$Sym_ID,display_numbers = T,number_color = "black")
#dev.off()


####
# kidney

tfpro_crossomeannokidneyselect <- tfproanno[nchar(tfproanno$Kidney.Pro.ID) > 0,c("Gene.Name","Kidney.Pro.ID","Ensembl")]
tfpro_crossomeannokidneyselect <- tfpro_crossomeannokidneyselect[!duplicated(tfpro_crossomeannokidneyselect$Kidney.Pro.ID),]

tfpro_crossomeannokidneyselect <- tfpro_crossomeannokidneyselect[tfpro_crossomeannokidneyselect$Ensembl %in% rownames(kidneyrnanorm),]

tfkidneyprol2fc <- kidneyprol2fc[tfpro_crossomeannokidneyselect$Kidney.Pro.ID,]

tfkidneypro_prol2fcmat <- kidneyprol2fc[tfpro_crossomeannokidneyselect$Kidney.Pro.ID,]
tfkidneypro_rnal2fcmat <- kidneyl2fcmat[tfpro_crossomeannokidneyselect$Ensembl,]

tfkidneypro_rnanorm <- kidneyrnanorm[rownames(tfkidneypro_rnal2fcmat),]
for(i in 1:dim(tfkidneypro_rnanorm)[2]){
  ourid <- colnames(tfkidneypro_rnanorm)[i]
  colnames(tfkidneypro_rnanorm)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0001"]
}
tfkidneypro_pronorm <- kidneyproprnorm[rownames(tfkidneypro_prol2fcmat),]

kidneyprosamplesubset <- intersect(colnames(tfkidneypro_rnanorm),colnames(tfkidneypro_pronorm))
tfkidneypro_rnanorm <- tfkidneypro_rnanorm[,kidneyprosamplesubset]
tfkidneypro_pronorm <- tfkidneypro_pronorm[,kidneyprosamplesubset]

tfkidneypro_crossomecor <- matrix(0L,nrow = dim(tfkidneypro_pronorm)[1],ncol = 1)
rownames(tfkidneypro_crossomecor) <- rownames(tfkidneypro_pronorm)
colnames(tfkidneypro_crossomecor) <- c("RNA-Pro")
tfkidneypro_crossomecortest <- matrix(0L,nrow = dim(tfkidneypro_pronorm)[1],ncol = 1)
rownames(tfkidneypro_crossomecortest) <- rownames(tfkidneypro_pronorm)
colnames(tfkidneypro_crossomecortest) <- c("RNA-Pro")
for(i in 1:dim(tfkidneypro_crossomecor)[1]){
  ourtest <- cor.test(t(tfkidneypro_rnanorm[i,]),t(tfkidneypro_pronorm[i,]))
  tfkidneypro_crossomecor[i,"RNA-Pro"] <- ourtest$estimate
  tfkidneypro_crossomecortest[i,"RNA-Pro"] <- ourtest$p.value
}

#png("kidney_tfpro_crossomecor.png",width = 2.5,height = 12,units = "in",res = 600)
#pheatmap(tfkidneypro_crossomecor,cluster_cols = F,angle_col = 0,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_row = tfpro_crossomeannokidneyselect$Gene.Name,display_numbers = T,number_color = "black")
#dev.off()

tfkidneyproph <- prot_ph$training_dea[prot_ph$training_dea$pr_feature_ID %in% tfproanno$Kidney.Pro.ID & prot_ph$training_dea$tissue_abbreviation %in% "KIDNEY","feature_ID"]
tfkidneyproph_prophl2fcmat <- kidneyprophl2fc[tfkidneyproph,]
tfkidneyproph_prol2fcmat <- matrix(0L,nrow = dim(tfkidneyproph_prophl2fcmat)[1],ncol = 8)
colnames(tfkidneyproph_prol2fcmat) <- colnames(tfkidneyproph_prophl2fcmat)
rownames(tfkidneyproph_prol2fcmat) <- rownames(tfkidneyproph_prophl2fcmat)
for(i in 1:dim(tfkidneyproph_prol2fcmat)[1]){
  ourid <- rownames(tfkidneyproph_prophl2fcmat)[i]
  rownames(tfkidneyproph_prol2fcmat)[i] <- prot_ph$training_dea[prot_ph$training_dea$feature_ID %in% ourid,"pr_feature_ID"][1]
}
tfkidneyproph_prol2fcmat <- kidneyprol2fc[rownames(tfkidneyproph_prol2fcmat),]

tfkidneyproph_rnal2fcmat <- matrix(0L,nrow = dim(tfkidneyproph_prophl2fcmat)[1],ncol = 8)
colnames(tfkidneyproph_rnal2fcmat) <- colnames(tfkidneyproph_prophl2fcmat)
rownames(tfkidneyproph_rnal2fcmat) <- rownames(tfkidneyproph_prol2fcmat)
for(i in 1:dim(tfkidneyproph_rnal2fcmat)[1]){
  ourid <- rownames(tfkidneyproph_rnal2fcmat)[i]
  rownames(tfkidneyproph_rnal2fcmat)[i] <- tfproanno[grep(ourid,tfproanno$Kidney.Pro.ID),"Ensembl"][1]
}

tfkidneyproph_prol2fcmat <- tfkidneyproph_prol2fcmat[rownames(tfkidneyproph_rnal2fcmat) %in% rownames(kidneyl2fcmat),]
tfkidneyproph_prophl2fcmat <- tfkidneyproph_prophl2fcmat[rownames(tfkidneyproph_rnal2fcmat) %in% rownames(kidneyl2fcmat),]
tfkidneyproph_rnal2fcmat <- tfkidneyproph_rnal2fcmat[rownames(tfkidneyproph_rnal2fcmat) %in% rownames(kidneyl2fcmat),]

tfkidneyproph_rnal2fcmat <- kidneyl2fcmat[rownames(tfkidneyproph_rnal2fcmat),]

tfkidneyproph_rnanorm <- kidneyrnanorm[rownames(tfkidneyproph_rnal2fcmat),]
for(i in 1:dim(tfkidneyproph_rnanorm)[2]){
  ourid <- colnames(tfkidneyproph_rnanorm)[i]
  colnames(tfkidneyproph_rnanorm)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0001"]
}

tfkidneyproph_pronorm <- kidneyproprnorm[rownames(tfkidneyproph_prol2fcmat),]
tfkidneyproph_prophnorm <- kidneyprophnorm[rownames(tfkidneyproph_prophl2fcmat),]
for(i in 1:dim(tfkidneyproph_prophnorm)[2]){
  ourid <- colnames(tfkidneyproph_prophnorm)[i]
  colnames(tfkidneyproph_prophnorm)[i] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$viallabel %in% ourid,"pid"]
}


kidneyprophsamplesubset <- Reduce(intersect,list(colnames(tfkidneyproph_rnanorm),colnames(tfkidneyproph_pronorm),colnames(tfkidneyproph_prophnorm)))
tfkidneyproph_rnanorm <- tfkidneyproph_rnanorm[,kidneyprophsamplesubset]
tfkidneyproph_pronorm <- tfkidneyproph_pronorm[,kidneyprophsamplesubset]
tfkidneyproph_prophnorm <- tfkidneyproph_prophnorm[,kidneyprophsamplesubset]

tfkidneyproph_crossomecor <- matrix(0L,nrow = dim(tfkidneyproph_prophnorm)[1],ncol = 3)
rownames(tfkidneyproph_crossomecor) <- rownames(tfkidneyproph_prophnorm)
colnames(tfkidneyproph_crossomecor) <- c("RNA-Pro","RNA-Phos","Pro-Phos")
tfkidneyproph_crossomecortest <- matrix(0L,nrow = dim(tfkidneyproph_prophnorm)[1],ncol = 3)
rownames(tfkidneyproph_crossomecortest) <- rownames(tfkidneyproph_prophnorm)
colnames(tfkidneyproph_crossomecortest) <- c("RNA-Pro","RNA-Phos","Pro-Phos")

for(i in 1:dim(tfkidneyproph_crossomecor)[1]){
  ourtest <- cor.test(t(tfkidneyproph_rnanorm[i,]),t(tfkidneyproph_pronorm[i,]))
  tfkidneyproph_crossomecor[i,"RNA-Pro"] <- ourtest$estimate
  tfkidneyproph_crossomecortest[i,"RNA-Pro"] <- ourtest$p.value
  ourtest <- cor.test(t(tfkidneyproph_rnanorm[i,]),t(tfkidneyproph_prophnorm[i,]))
  tfkidneyproph_crossomecor[i,"RNA-Phos"] <- ourtest$estimate
  tfkidneyproph_crossomecortest[i,"RNA-Phos"] <- ourtest$p.value
  ourtest <- cor.test(t(tfkidneyproph_pronorm[i,]),t(tfkidneyproph_prophnorm[i,]))
  tfkidneyproph_crossomecor[i,"Pro-Phos"] <- ourtest$estimate
  tfkidneyproph_crossomecortest[i,"Pro-Phos"] <- ourtest$p.value
}

tfproph_crossomeannokidneyselect <- data.frame(row.names = rownames(tfkidneyproph_crossomecor),"Phospho" = rownames(tfkidneyproph_crossomecor),"Protein_ID" = sapply(strsplit(rownames(tfkidneyproph_crossomecor), "_"), function(x) {
  paste(x[1],x[2],sep = "_")
}))
tfproph_crossomeannokidneyselect$Symbol = ""
for(i in 1:length(rownames(tfkidneyproph_crossomecor))){
  ourpro <- tfproph_crossomeannokidneyselect[i,"Protein_ID"]
  tfproph_crossomeannokidneyselect[i,"Symbol"] <- toupper(tfproanno[tfproanno$Kidney.Pro.ID %in% ourpro,"Gene.Name"][1])
}

tfproph_crossomeannokidneyselect$Sym_ID <- paste(toupper(tfproph_crossomeannokidneyselect$Symbol),gsub(".*_","",tfproph_crossomeannokidneyselect$Phospho),sep = "_")

#png("kidney_tfproph_crossomecor.png",width = 7.5,height = 12,units = "in",res = 600)
#pheatmap(tfkidneyproph_crossomecor,cluster_cols = F,angle_col = 0,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_row = tfproph_crossomeannokidneyselect$Sym_ID,display_numbers = T,number_color = "black")
#dev.off()


####
# liver

tfpro_crossomeannoliverselect <- tfproanno[nchar(tfproanno$Liver.Pro.ID) > 0,c("Gene.Name","Liver.Pro.ID","Ensembl")]
tfpro_crossomeannoliverselect <- tfpro_crossomeannoliverselect[!duplicated(tfpro_crossomeannoliverselect$Liver.Pro.ID),]

tfpro_crossomeannoliverselect <- tfpro_crossomeannoliverselect[tfpro_crossomeannoliverselect$Ensembl %in% rownames(liverrnanorm),]

tfliverprol2fc <- liverprol2fc[tfpro_crossomeannoliverselect$Liver.Pro.ID,]

tfliverpro_prol2fcmat <- liverprol2fc[tfpro_crossomeannoliverselect$Liver.Pro.ID,]
tfliverpro_rnal2fcmat <- liverl2fcmat[tfpro_crossomeannoliverselect$Ensembl,]

tfliverpro_rnanorm <- liverrnanorm[rownames(tfliverpro_rnal2fcmat),]
for(i in 1:dim(tfliverpro_rnanorm)[2]){
  ourid <- colnames(tfliverpro_rnanorm)[i]
  colnames(tfliverpro_rnanorm)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0001"]
}
tfliverpro_pronorm <- liverproprnorm[rownames(tfliverpro_prol2fcmat),]

liverprosamplesubset <- intersect(colnames(tfliverpro_rnanorm),colnames(tfliverpro_pronorm))
tfliverpro_rnanorm <- tfliverpro_rnanorm[,liverprosamplesubset]
tfliverpro_pronorm <- tfliverpro_pronorm[,liverprosamplesubset]

tfliverpro_crossomecor <- matrix(0L,nrow = dim(tfliverpro_pronorm)[1],ncol = 1)
rownames(tfliverpro_crossomecor) <- rownames(tfliverpro_pronorm)
colnames(tfliverpro_crossomecor) <- c("RNA-Pro")
tfliverpro_crossomecortest <- matrix(0L,nrow = dim(tfliverpro_pronorm)[1],ncol = 1)
rownames(tfliverpro_crossomecortest) <- rownames(tfliverpro_pronorm)
colnames(tfliverpro_crossomecortest) <- c("RNA-Pro")
for(i in 1:dim(tfliverpro_crossomecor)[1]){
  ourtest <- cor.test(t(tfliverpro_rnanorm[i,]),t(tfliverpro_pronorm[i,]))
  tfliverpro_crossomecor[i,"RNA-Pro"] <- ourtest$estimate
  tfliverpro_crossomecortest[i,"RNA-Pro"] <- ourtest$p.value
}

#png("liver_tfpro_crossomecor.png",width = 2.5,height = 12,units = "in",res = 600)
#pheatmap(tfliverpro_crossomecor,cluster_cols = F,angle_col = 0,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_row = tfpro_crossomeannoliverselect$Gene.Name,display_numbers = T,number_color = "black")
#dev.off()

tfliverproph <- prot_ph$training_dea[prot_ph$training_dea$pr_feature_ID %in% tfproanno$Liver.Pro.ID & prot_ph$training_dea$tissue_abbreviation %in% "LIVER","feature_ID"]
tfliverproph_prophl2fcmat <- liverprophl2fc[tfliverproph,]
tfliverproph_prol2fcmat <- matrix(0L,nrow = dim(tfliverproph_prophl2fcmat)[1],ncol = 8)
colnames(tfliverproph_prol2fcmat) <- colnames(tfliverproph_prophl2fcmat)
rownames(tfliverproph_prol2fcmat) <- rownames(tfliverproph_prophl2fcmat)
for(i in 1:dim(tfliverproph_prol2fcmat)[1]){
  ourid <- rownames(tfliverproph_prophl2fcmat)[i]
  rownames(tfliverproph_prol2fcmat)[i] <- prot_ph$training_dea[prot_ph$training_dea$feature_ID %in% ourid,"pr_feature_ID"][1]
}
tfliverproph_prol2fcmat <- liverprol2fc[rownames(tfliverproph_prol2fcmat),]

tfliverproph_rnal2fcmat <- matrix(0L,nrow = dim(tfliverproph_prophl2fcmat)[1],ncol = 8)
colnames(tfliverproph_rnal2fcmat) <- colnames(tfliverproph_prophl2fcmat)
rownames(tfliverproph_rnal2fcmat) <- rownames(tfliverproph_prol2fcmat)
for(i in 1:dim(tfliverproph_rnal2fcmat)[1]){
  ourid <- rownames(tfliverproph_rnal2fcmat)[i]
  rownames(tfliverproph_rnal2fcmat)[i] <- tfproanno[grep(ourid,tfproanno$Liver.Pro.ID),"Ensembl"][1]
}

tfliverproph_prol2fcmat <- tfliverproph_prol2fcmat[rownames(tfliverproph_rnal2fcmat) %in% rownames(liverl2fcmat),]
tfliverproph_prophl2fcmat <- tfliverproph_prophl2fcmat[rownames(tfliverproph_rnal2fcmat) %in% rownames(liverl2fcmat),]
tfliverproph_rnal2fcmat <- tfliverproph_rnal2fcmat[rownames(tfliverproph_rnal2fcmat) %in% rownames(liverl2fcmat),]

tfliverproph_rnal2fcmat <- liverl2fcmat[rownames(tfliverproph_rnal2fcmat),]

tfliverproph_rnanorm <- liverrnanorm[rownames(tfliverproph_rnal2fcmat),]
for(i in 1:dim(tfliverproph_rnanorm)[2]){
  ourid <- colnames(tfliverproph_rnanorm)[i]
  colnames(tfliverproph_rnanorm)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0001"]
}

tfliverproph_pronorm <- liverproprnorm[rownames(tfliverproph_prol2fcmat),]
tfliverproph_prophnorm <- liverprophnorm[rownames(tfliverproph_prophl2fcmat),]
for(i in 1:dim(tfliverproph_prophnorm)[2]){
  ourid <- colnames(tfliverproph_prophnorm)[i]
  colnames(tfliverproph_prophnorm)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0001"]
}

liverprophsamplesubset <- Reduce(intersect,list(colnames(tfliverproph_rnanorm),colnames(tfliverproph_pronorm),colnames(tfliverproph_prophnorm)))
tfliverproph_rnanorm <- tfliverproph_rnanorm[,liverprophsamplesubset]
tfliverproph_pronorm <- tfliverproph_pronorm[,liverprophsamplesubset]
tfliverproph_prophnorm <- tfliverproph_prophnorm[,liverprophsamplesubset]

tfliverproph_crossomecor <- matrix(0L,nrow = dim(tfliverproph_prophnorm)[1],ncol = 3)
rownames(tfliverproph_crossomecor) <- rownames(tfliverproph_prophnorm)
colnames(tfliverproph_crossomecor) <- c("RNA-Pro","RNA-Phos","Pro-Phos")
tfliverproph_crossomecortest <- matrix(0L,nrow = dim(tfliverproph_prophnorm)[1],ncol = 3)
rownames(tfliverproph_crossomecortest) <- rownames(tfliverproph_prophnorm)
colnames(tfliverproph_crossomecortest) <- c("RNA-Pro","RNA-Phos","Pro-Phos")

for(i in 1:dim(tfliverproph_crossomecor)[1]){
  ourtest <- cor.test(t(tfliverproph_rnanorm[i,]),t(tfliverproph_pronorm[i,]))
  tfliverproph_crossomecor[i,"RNA-Pro"] <- ourtest$estimate
  tfliverproph_crossomecortest[i,"RNA-Pro"] <- ourtest$p.value
  ourtest <- cor.test(t(tfliverproph_rnanorm[i,]),t(tfliverproph_prophnorm[i,]))
  tfliverproph_crossomecor[i,"RNA-Phos"] <- ourtest$estimate
  tfliverproph_crossomecortest[i,"RNA-Phos"] <- ourtest$p.value
  ourtest <- cor.test(t(tfliverproph_pronorm[i,]),t(tfliverproph_prophnorm[i,]))
  tfliverproph_crossomecor[i,"Pro-Phos"] <- ourtest$estimate
  tfliverproph_crossomecortest[i,"Pro-Phos"] <- ourtest$p.value
}

tfproph_crossomeannoliverselect <- data.frame(row.names = rownames(tfliverproph_crossomecor),"Phospho" = rownames(tfliverproph_crossomecor),"Protein_ID" = sapply(strsplit(rownames(tfliverproph_crossomecor), "_"), function(x) {
  paste(x[1],x[2],sep = "_")
}))
tfproph_crossomeannoliverselect$Symbol = ""
for(i in 1:length(rownames(tfliverproph_crossomecor))){
  ourpro <- tfproph_crossomeannoliverselect[i,"Protein_ID"]
  tfproph_crossomeannoliverselect[i,"Symbol"] <- toupper(tfproanno[tfproanno$Liver.Pro.ID %in% ourpro,"Gene.Name"][1])
}

tfproph_crossomeannoliverselect$Sym_ID <- paste(toupper(tfproph_crossomeannoliverselect$Symbol),gsub(".*_","",tfproph_crossomeannoliverselect$Phospho),sep = "_")

#png("liver_tfproph_crossomecor.png",width = 7.5,height = 12,units = "in",res = 600)
#pheatmap(tfliverproph_crossomecor,cluster_cols = F,angle_col = 0,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_row = tfproph_crossomeannoliverselect$Sym_ID,display_numbers = T,number_color = "black")
#dev.off()


####
# lung

tfpro_crossomeannolungselect <- tfproanno[nchar(tfproanno$Lung.Pro.ID) > 0,c("Gene.Name","Lung.Pro.ID","Ensembl")]
tfpro_crossomeannolungselect <- tfpro_crossomeannolungselect[!duplicated(tfpro_crossomeannolungselect$Lung.Pro.ID),]

tfpro_crossomeannolungselect <- tfpro_crossomeannolungselect[tfpro_crossomeannolungselect$Ensembl %in% rownames(lungrnanorm),]

tflungprol2fc <- lungprol2fc[tfpro_crossomeannolungselect$Lung.Pro.ID,]

tflungpro_prol2fcmat <- lungprol2fc[tfpro_crossomeannolungselect$Lung.Pro.ID,]
tflungpro_rnal2fcmat <- lungl2fcmat[tfpro_crossomeannolungselect$Ensembl,]

tflungpro_rnanorm <- lungrnanorm[rownames(tflungpro_rnal2fcmat),]
for(i in 1:dim(tflungpro_rnanorm)[2]){
  ourid <- colnames(tflungpro_rnanorm)[i]
  colnames(tflungpro_rnanorm)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0001"]
}
tflungpro_pronorm <- lungproprnorm[rownames(tflungpro_prol2fcmat),]

lungprosamplesubset <- intersect(colnames(tflungpro_rnanorm),colnames(tflungpro_pronorm))
tflungpro_rnanorm <- tflungpro_rnanorm[,lungprosamplesubset]
tflungpro_pronorm <- tflungpro_pronorm[,lungprosamplesubset]

tflungpro_crossomecor <- matrix(0L,nrow = dim(tflungpro_pronorm)[1],ncol = 1)
rownames(tflungpro_crossomecor) <- rownames(tflungpro_pronorm)
colnames(tflungpro_crossomecor) <- c("RNA-Pro")
tflungpro_crossomecortest <- matrix(0L,nrow = dim(tflungpro_pronorm)[1],ncol = 1)
rownames(tflungpro_crossomecortest) <- rownames(tflungpro_pronorm)
colnames(tflungpro_crossomecortest) <- c("RNA-Pro")
for(i in 1:dim(tflungpro_crossomecor)[1]){
  ourtest <- cor.test(t(tflungpro_rnanorm[i,]),t(tflungpro_pronorm[i,]))
  tflungpro_crossomecor[i,"RNA-Pro"] <- ourtest$estimate
  tflungpro_crossomecortest[i,"RNA-Pro"] <- ourtest$p.value
}

#png("lung_tfpro_crossomecor.png",width = 2.5,height = 12,units = "in",res = 600)
#pheatmap(tflungpro_crossomecor,cluster_cols = F,angle_col = 0,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_row = tfpro_crossomeannolungselect$Gene.Name,display_numbers = T,number_color = "black")
#dev.off()

tflungproph <- prot_ph$training_dea[prot_ph$training_dea$pr_feature_ID %in% tfproanno$Lung.Pro.ID & prot_ph$training_dea$tissue_abbreviation %in% "LUNG","feature_ID"]
tflungproph_prophl2fcmat <- lungprophl2fc[tflungproph,]
tflungproph_prol2fcmat <- matrix(0L,nrow = dim(tflungproph_prophl2fcmat)[1],ncol = 8)
colnames(tflungproph_prol2fcmat) <- colnames(tflungproph_prophl2fcmat)
rownames(tflungproph_prol2fcmat) <- rownames(tflungproph_prophl2fcmat)
for(i in 1:dim(tflungproph_prol2fcmat)[1]){
  ourid <- rownames(tflungproph_prophl2fcmat)[i]
  rownames(tflungproph_prol2fcmat)[i] <- prot_ph$training_dea[prot_ph$training_dea$feature_ID %in% ourid,"pr_feature_ID"][1]
}
tflungproph_prol2fcmat <- lungprol2fc[rownames(tflungproph_prol2fcmat),]

tflungproph_rnal2fcmat <- matrix(0L,nrow = dim(tflungproph_prophl2fcmat)[1],ncol = 8)
colnames(tflungproph_rnal2fcmat) <- colnames(tflungproph_prophl2fcmat)
rownames(tflungproph_rnal2fcmat) <- rownames(tflungproph_prol2fcmat)
for(i in 1:dim(tflungproph_rnal2fcmat)[1]){
  ourid <- rownames(tflungproph_rnal2fcmat)[i]
  rownames(tflungproph_rnal2fcmat)[i] <- tfproanno[grep(ourid,tfproanno$Lung.Pro.ID),"Ensembl"][1]
}

tflungproph_prol2fcmat <- tflungproph_prol2fcmat[rownames(tflungproph_rnal2fcmat) %in% rownames(lungl2fcmat),]
tflungproph_prophl2fcmat <- tflungproph_prophl2fcmat[rownames(tflungproph_rnal2fcmat) %in% rownames(lungl2fcmat),]
tflungproph_rnal2fcmat <- tflungproph_rnal2fcmat[rownames(tflungproph_rnal2fcmat) %in% rownames(lungl2fcmat),]

tflungproph_rnal2fcmat <- lungl2fcmat[rownames(tflungproph_rnal2fcmat),]

tflungproph_rnanorm <- lungrnanorm[rownames(tflungproph_rnal2fcmat),]
for(i in 1:dim(tflungproph_rnanorm)[2]){
  ourid <- colnames(tflungproph_rnanorm)[i]
  colnames(tflungproph_rnanorm)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0001"]
}

tflungproph_pronorm <- lungproprnorm[rownames(tflungproph_prol2fcmat),]
tflungproph_prophnorm <- lungprophnorm[rownames(tflungproph_prophl2fcmat),]
for(i in 1:dim(tflungproph_prophnorm)[2]){
  ourid <- colnames(tflungproph_prophnorm)[i]
  colnames(tflungproph_prophnorm)[i] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$viallabel %in% ourid,"pid"]
}

lungprophsamplesubset <- Reduce(intersect,list(colnames(tflungproph_rnanorm),colnames(tflungproph_pronorm),colnames(tflungproph_prophnorm)))
tflungproph_rnanorm <- tflungproph_rnanorm[,lungprophsamplesubset]
tflungproph_pronorm <- tflungproph_pronorm[,lungprophsamplesubset]
tflungproph_prophnorm <- tflungproph_prophnorm[,lungprophsamplesubset]

tflungproph_crossomecor <- matrix(0L,nrow = dim(tflungproph_prophnorm)[1],ncol = 3)
rownames(tflungproph_crossomecor) <- rownames(tflungproph_prophnorm)
colnames(tflungproph_crossomecor) <- c("RNA-Pro","RNA-Phos","Pro-Phos")
tflungproph_crossomecortest <- matrix(0L,nrow = dim(tflungproph_prophnorm)[1],ncol = 3)
rownames(tflungproph_crossomecortest) <- rownames(tflungproph_prophnorm)
colnames(tflungproph_crossomecortest) <- c("RNA-Pro","RNA-Phos","Pro-Phos")

for(i in 1:dim(tflungproph_crossomecor)[1]){
  ourtest <- cor.test(t(tflungproph_rnanorm[i,]),t(tflungproph_pronorm[i,]))
  tflungproph_crossomecor[i,"RNA-Pro"] <- ourtest$estimate
  tflungproph_crossomecortest[i,"RNA-Pro"] <- ourtest$p.value
  ourtest <- cor.test(t(tflungproph_rnanorm[i,]),t(tflungproph_prophnorm[i,]))
  tflungproph_crossomecor[i,"RNA-Phos"] <- ourtest$estimate
  tflungproph_crossomecortest[i,"RNA-Phos"] <- ourtest$p.value
  ourtest <- cor.test(t(tflungproph_pronorm[i,]),t(tflungproph_prophnorm[i,]))
  tflungproph_crossomecor[i,"Pro-Phos"] <- ourtest$estimate
  tflungproph_crossomecortest[i,"Pro-Phos"] <- ourtest$p.value
}

tfproph_crossomeannolungselect <- data.frame(row.names = rownames(tflungproph_crossomecor),"Phospho" = rownames(tflungproph_crossomecor),"Protein_ID" = sapply(strsplit(rownames(tflungproph_crossomecor), "_"), function(x) {
  paste(x[1],x[2],sep = "_")
}))
tfproph_crossomeannolungselect$Symbol = ""
for(i in 1:length(rownames(tflungproph_crossomecor))){
  ourpro <- tfproph_crossomeannolungselect[i,"Protein_ID"]
  tfproph_crossomeannolungselect[i,"Symbol"] <- toupper(tfproanno[tfproanno$Lung.Pro.ID %in% ourpro,"Gene.Name"][1])
}

tfproph_crossomeannolungselect$Sym_ID <- paste(toupper(tfproph_crossomeannolungselect$Symbol),gsub(".*_","",tfproph_crossomeannolungselect$Phospho),sep = "_")

#png("lung_tfproph_crossomecor.png",width = 7.5,height = 12,units = "in",res = 600)
#pheatmap(tflungproph_crossomecor,cluster_cols = F,angle_col = 0,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_row = tfproph_crossomeannolungselect$Sym_ID,display_numbers = T,number_color = "black")
#dev.off()


####
# white

tfpro_crossomeannowhiteselect <- tfproanno[nchar(tfproanno$WhiteAd.Pro.ID) > 0,c("Gene.Name","WhiteAd.Pro.ID","Ensembl")]
tfpro_crossomeannowhiteselect <- tfpro_crossomeannowhiteselect[!duplicated(tfpro_crossomeannowhiteselect$WhiteAd.Pro.ID),]

tfpro_crossomeannowhiteselect <- tfpro_crossomeannowhiteselect[tfpro_crossomeannowhiteselect$Ensembl %in% rownames(whiternanorm),]

tfwhiteprol2fc <- whiteprol2fc[tfpro_crossomeannowhiteselect$WhiteAd.Pro.ID,]

tfwhitepro_prol2fcmat <- whiteprol2fc[tfpro_crossomeannowhiteselect$WhiteAd.Pro.ID,]
tfwhitepro_rnal2fcmat <- whitel2fcmat[tfpro_crossomeannowhiteselect$Ensembl,]

tfwhitepro_rnanorm <- whiternanorm[rownames(tfwhitepro_rnal2fcmat),]
for(i in 1:dim(tfwhitepro_rnanorm)[2]){
  ourid <- colnames(tfwhitepro_rnanorm)[i]
  colnames(tfwhitepro_rnanorm)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0001"]
}
tfwhitepro_pronorm <- whiteproprnorm[rownames(tfwhitepro_prol2fcmat),]

whiteprosamplesubset <- intersect(colnames(tfwhitepro_rnanorm),colnames(tfwhitepro_pronorm))
tfwhitepro_rnanorm <- tfwhitepro_rnanorm[,whiteprosamplesubset]
tfwhitepro_pronorm <- tfwhitepro_pronorm[,whiteprosamplesubset]

tfwhitepro_crossomecor <- matrix(0L,nrow = dim(tfwhitepro_pronorm)[1],ncol = 1)
rownames(tfwhitepro_crossomecor) <- rownames(tfwhitepro_pronorm)
colnames(tfwhitepro_crossomecor) <- c("RNA-Pro")
tfwhitepro_crossomecortest <- matrix(0L,nrow = dim(tfwhitepro_pronorm)[1],ncol = 1)
rownames(tfwhitepro_crossomecortest) <- rownames(tfwhitepro_pronorm)
colnames(tfwhitepro_crossomecortest) <- c("RNA-Pro")
for(i in 1:dim(tfwhitepro_crossomecor)[1]){
  ourtest <- cor.test(t(tfwhitepro_rnanorm[i,]),t(tfwhitepro_pronorm[i,]))
  tfwhitepro_crossomecor[i,"RNA-Pro"] <- ourtest$estimate
  tfwhitepro_crossomecortest[i,"RNA-Pro"] <- ourtest$p.value
}

#png("white_tfpro_crossomecor.png",width = 2.5,height = 12,units = "in",res = 600)
#pheatmap(tfwhitepro_crossomecor,cluster_cols = F,angle_col = 0,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_row = tfpro_crossomeannowhiteselect$Gene.Name,display_numbers = T,number_color = "black")
#dev.off()

tfwhiteproph <- prot_ph$training_dea[prot_ph$training_dea$pr_feature_ID %in% tfproanno$WhiteAd.Pro.ID & prot_ph$training_dea$tissue_abbreviation %in% "WAT-SC","feature_ID"]
tfwhiteproph_prophl2fcmat <- whiteprophl2fc[tfwhiteproph,]
tfwhiteproph_prol2fcmat <- matrix(0L,nrow = dim(tfwhiteproph_prophl2fcmat)[1],ncol = 8)
colnames(tfwhiteproph_prol2fcmat) <- colnames(tfwhiteproph_prophl2fcmat)
rownames(tfwhiteproph_prol2fcmat) <- rownames(tfwhiteproph_prophl2fcmat)
for(i in 1:dim(tfwhiteproph_prol2fcmat)[1]){
  ourid <- rownames(tfwhiteproph_prophl2fcmat)[i]
  rownames(tfwhiteproph_prol2fcmat)[i] <- prot_ph$training_dea[prot_ph$training_dea$feature_ID %in% ourid,"pr_feature_ID"][1]
}
tfwhiteproph_prol2fcmat <- whiteprol2fc[rownames(tfwhiteproph_prol2fcmat),]

tfwhiteproph_rnal2fcmat <- matrix(0L,nrow = dim(tfwhiteproph_prophl2fcmat)[1],ncol = 8)
colnames(tfwhiteproph_rnal2fcmat) <- colnames(tfwhiteproph_prophl2fcmat)
rownames(tfwhiteproph_rnal2fcmat) <- rownames(tfwhiteproph_prol2fcmat)
for(i in 1:dim(tfwhiteproph_rnal2fcmat)[1]){
  ourid <- rownames(tfwhiteproph_rnal2fcmat)[i]
  rownames(tfwhiteproph_rnal2fcmat)[i] <- tfproanno[grep(ourid,tfproanno$WhiteAd.Pro.ID),"Ensembl"][1]
}

tfwhiteproph_prol2fcmat <- tfwhiteproph_prol2fcmat[rownames(tfwhiteproph_rnal2fcmat) %in% rownames(whitel2fcmat),]
tfwhiteproph_prophl2fcmat <- tfwhiteproph_prophl2fcmat[rownames(tfwhiteproph_rnal2fcmat) %in% rownames(whitel2fcmat),]
tfwhiteproph_rnal2fcmat <- tfwhiteproph_rnal2fcmat[rownames(tfwhiteproph_rnal2fcmat) %in% rownames(whitel2fcmat),]

tfwhiteproph_rnal2fcmat <- whitel2fcmat[rownames(tfwhiteproph_rnal2fcmat),]

tfwhiteproph_rnanorm <- whiternanorm[rownames(tfwhiteproph_rnal2fcmat),]
for(i in 1:dim(tfwhiteproph_rnanorm)[2]){
  ourid <- colnames(tfwhiteproph_rnanorm)[i]
  colnames(tfwhiteproph_rnanorm)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0001"]
}

tfwhiteproph_pronorm <- whiteproprnorm[rownames(tfwhiteproph_prol2fcmat),]
tfwhiteproph_prophnorm <- whiteprophnorm[rownames(tfwhiteproph_prophl2fcmat),]
for(i in 1:dim(tfwhiteproph_prophnorm)[2]){
  ourid <- colnames(tfwhiteproph_prophnorm)[i]
  colnames(tfwhiteproph_prophnorm)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0001"]
}

whiteprophsamplesubset <- Reduce(intersect,list(colnames(tfwhiteproph_rnanorm),colnames(tfwhiteproph_pronorm),colnames(tfwhiteproph_prophnorm)))
tfwhiteproph_rnanorm <- tfwhiteproph_rnanorm[,whiteprophsamplesubset]
tfwhiteproph_pronorm <- tfwhiteproph_pronorm[,whiteprophsamplesubset]
tfwhiteproph_prophnorm <- tfwhiteproph_prophnorm[,whiteprophsamplesubset]

tfwhiteproph_crossomecor <- matrix(0L,nrow = dim(tfwhiteproph_prophnorm)[1],ncol = 3)
rownames(tfwhiteproph_crossomecor) <- rownames(tfwhiteproph_prophnorm)
colnames(tfwhiteproph_crossomecor) <- c("RNA-Pro","RNA-Phos","Pro-Phos")
tfwhiteproph_crossomecortest <- matrix(0L,nrow = dim(tfwhiteproph_prophnorm)[1],ncol = 3)
rownames(tfwhiteproph_crossomecortest) <- rownames(tfwhiteproph_prophnorm)
colnames(tfwhiteproph_crossomecortest) <- c("RNA-Pro","RNA-Phos","Pro-Phos")

for(i in 1:dim(tfwhiteproph_crossomecor)[1]){
  ourtest <- cor.test(t(tfwhiteproph_rnanorm[i,]),t(tfwhiteproph_pronorm[i,]))
  tfwhiteproph_crossomecor[i,"RNA-Pro"] <- ourtest$estimate
  tfwhiteproph_crossomecortest[i,"RNA-Pro"] <- ourtest$p.value
  ourtest <- cor.test(t(tfwhiteproph_rnanorm[i,]),t(tfwhiteproph_prophnorm[i,]))
  tfwhiteproph_crossomecor[i,"RNA-Phos"] <- ourtest$estimate
  tfwhiteproph_crossomecortest[i,"RNA-Phos"] <- ourtest$p.value
  ourtest <- cor.test(t(tfwhiteproph_pronorm[i,]),t(tfwhiteproph_prophnorm[i,]))
  tfwhiteproph_crossomecor[i,"Pro-Phos"] <- ourtest$estimate
  tfwhiteproph_crossomecortest[i,"Pro-Phos"] <- ourtest$p.value
}

tfproph_crossomeannowhiteselect <- data.frame(row.names = rownames(tfwhiteproph_crossomecor),"Phospho" = rownames(tfwhiteproph_crossomecor),"Protein_ID" = sapply(strsplit(rownames(tfwhiteproph_crossomecor), "_"), function(x) {
  paste(x[1],x[2],sep = "_")
}))
tfproph_crossomeannowhiteselect$Symbol = ""
for(i in 1:length(rownames(tfwhiteproph_crossomecor))){
  ourpro <- tfproph_crossomeannowhiteselect[i,"Protein_ID"]
  tfproph_crossomeannowhiteselect[i,"Symbol"] <- toupper(tfproanno[tfproanno$WhiteAd.Pro.ID %in% ourpro,"Gene.Name"][1])
}

tfproph_crossomeannowhiteselect$Sym_ID <- paste(toupper(tfproph_crossomeannowhiteselect$Symbol),gsub(".*_","",tfproph_crossomeannowhiteselect$Phospho),sep = "_")

#png("white_tfproph_crossomecor.png",width = 7.5,height = 12,units = "in",res = 600)
#pheatmap(tfwhiteproph_crossomecor,cluster_cols = F,angle_col = 0,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_row = tfproph_crossomeannowhiteselect$Sym_ID,display_numbers = T,number_color = "black")
#dev.off()


allproph_cross_pvals <- c(as.vector(tfgastroproph_crossomecortest),
              as.vector(tfheartproph_crossomecortest),
              as.vector(tfkidneyproph_crossomecortest),
              as.vector(tfliverproph_crossomecortest),
              as.vector(tflungproph_crossomecortest),
              as.vector(tfwhiteproph_crossomecortest))

allproph_cross_adjpvals <- p.adjust(allproph_cross_pvals, method = "BH", n = length(allproph_cross_pvals))

allpro_cross_pvals <- c(tfgastropro_crossomecortest,
                        tfheartpro_crossomecortest,
                        tfkidneypro_crossomecortest,
                        tfliverpro_crossomecortest,
                        tflungpro_crossomecortest,
                        tfwhitepro_crossomecortest)

allpro_cross_adjpvals <- p.adjust(allpro_cross_pvals, method = "BH", n = length(allpro_cross_pvals))

index = 1
tfgastroproph_crossomecortest_adj <- matrix(allproph_cross_adjpvals[index:(index+length(tfgastroproph_crossomecortest)-1)],ncol = 3)
index = index + length(tfgastroproph_crossomecortest)
tfheartproph_crossomecortest_adj <- matrix(allproph_cross_adjpvals[index:(index+length(tfheartproph_crossomecortest)-1)],ncol = 3)
index = index + length(tfheartproph_crossomecortest)
tfkidneyproph_crossomecortest_adj <- matrix(allproph_cross_adjpvals[index:(index+length(tfkidneyproph_crossomecortest)-1)],ncol = 3)
index = index + length(tfkidneyproph_crossomecortest)
tfliverproph_crossomecortest_adj <- matrix(allproph_cross_adjpvals[index:(index+length(tfliverproph_crossomecortest)-1)],ncol = 3)
index = index + length(tfliverproph_crossomecortest)
tflungproph_crossomecortest_adj <- matrix(allproph_cross_adjpvals[index:(index+length(tflungproph_crossomecortest)-1)],ncol = 3)
index = index + length(tflungproph_crossomecortest)
tfwhiteproph_crossomecortest_adj <- matrix(allproph_cross_adjpvals[index:(index+length(tfwhiteproph_crossomecortest)-1)],ncol = 3)

rownames(tfgastroproph_crossomecortest_adj) <- rownames(tfgastroproph_crossomecortest)
colnames(tfgastroproph_crossomecortest_adj) <- colnames(tfgastroproph_crossomecortest)
rownames(tfheartproph_crossomecortest_adj) <- rownames(tfheartproph_crossomecortest)
colnames(tfheartproph_crossomecortest_adj) <- colnames(tfheartproph_crossomecortest)
rownames(tfkidneyproph_crossomecortest_adj) <- rownames(tfkidneyproph_crossomecortest)
colnames(tfkidneyproph_crossomecortest_adj) <- colnames(tfkidneyproph_crossomecortest)
rownames(tfliverproph_crossomecortest_adj) <- rownames(tfliverproph_crossomecortest)
colnames(tfliverproph_crossomecortest_adj) <- colnames(tfliverproph_crossomecortest)
rownames(tflungproph_crossomecortest_adj) <- rownames(tflungproph_crossomecortest)
colnames(tflungproph_crossomecortest_adj) <- colnames(tflungproph_crossomecortest)
rownames(tfwhiteproph_crossomecortest_adj) <- rownames(tfwhiteproph_crossomecortest)
colnames(tfwhiteproph_crossomecortest_adj) <- colnames(tfwhiteproph_crossomecortest)

index = 1
tfgastropro_crossomecortest_adj <- matrix(allpro_cross_adjpvals[index:(index+length(tfgastropro_crossomecortest)-1)],ncol = 1)
index = index + length(tfgastropro_crossomecortest)
tfheartpro_crossomecortest_adj <- matrix(allpro_cross_adjpvals[index:(index+length(tfheartpro_crossomecortest)-1)],ncol = 1)
index = index + length(tfheartpro_crossomecortest)
tfkidneypro_crossomecortest_adj <- matrix(allpro_cross_adjpvals[index:(index+length(tfkidneypro_crossomecortest)-1)],ncol = 1)
index = index + length(tfkidneypro_crossomecortest)
tfliverpro_crossomecortest_adj <- matrix(allpro_cross_adjpvals[index:(index+length(tfliverpro_crossomecortest)-1)],ncol = 1)
index = index + length(tfliverpro_crossomecortest)
tflungpro_crossomecortest_adj <- matrix(allpro_cross_adjpvals[index:(index+length(tflungpro_crossomecortest)-1)],ncol = 1)
index = index + length(tflungpro_crossomecortest)
tfwhitepro_crossomecortest_adj <- matrix(allpro_cross_adjpvals[index:(index+length(tfwhitepro_crossomecortest)-1)],ncol = 1)

rownames(tfgastropro_crossomecortest_adj) <- rownames(tfgastropro_crossomecortest)
colnames(tfgastropro_crossomecortest_adj) <- colnames(tfgastropro_crossomecortest)
rownames(tfheartpro_crossomecortest_adj) <- rownames(tfheartpro_crossomecortest)
colnames(tfheartpro_crossomecortest_adj) <- colnames(tfheartpro_crossomecortest)
rownames(tfkidneypro_crossomecortest_adj) <- rownames(tfkidneypro_crossomecortest)
colnames(tfkidneypro_crossomecortest_adj) <- colnames(tfkidneypro_crossomecortest)
rownames(tfliverpro_crossomecortest_adj) <- rownames(tfliverpro_crossomecortest)
colnames(tfliverpro_crossomecortest_adj) <- colnames(tfliverpro_crossomecortest)
rownames(tflungpro_crossomecortest_adj) <- rownames(tflungpro_crossomecortest)
colnames(tflungpro_crossomecortest_adj) <- colnames(tflungpro_crossomecortest)
rownames(tfwhitepro_crossomecortest_adj) <- rownames(tfwhitepro_crossomecortest)
colnames(tfwhitepro_crossomecortest_adj) <- colnames(tfwhitepro_crossomecortest)

#gastro
tfgastroproph_crossomecor_df <- data.frame("Correlation" = c(tfgastroproph_crossomecor[,"RNA-Pro"],
                                                             tfgastroproph_crossomecor[,"RNA-Phos"],
                                                             tfgastroproph_crossomecor[,"Pro-Phos"]),
                                           "P.Value" = c(tfgastroproph_crossomecortest_adj[,"RNA-Pro"],
                                                         tfgastroproph_crossomecortest_adj[,"RNA-Phos"],
                                                         tfgastroproph_crossomecortest_adj[,"Pro-Phos"]),
                                           "Comparison" = c(rep("RNA-Pro",dim(tfgastroproph_crossomecortest_adj)[1]),
                                                            rep("RNA-Phos",dim(tfgastroproph_crossomecortest_adj)[1]),
                                                            rep("Pro-Phos",dim(tfgastroproph_crossomecortest_adj)[1])))
#tfgastroproph_crossomecor_df$Log2.P.Value <- -log2(tfgastroproph_crossomecor_df$Log2.P.Value)
tfgastroproph_crossomecor_df$Label <- rep("",dim(tfgastroproph_crossomecor_df)[1])
for(i in 1:dim(tfgastroproph_crossomecor)[1]){
  ourrna <- gsub("\\..*","",rownames(tfgastroproph_rnanorm))[i]
  ourpro <- sub('^([^\\.]+\\.[^\\.]+).*', '\\1', rownames(tfgastroproph_pronorm))[i]
  ourproph <- rownames(tfgastroproph_prophnorm)[i]
  ourrnasymbol <- enstosym[ourrna,"Symbol"]
  ourprosymbol <- toupper(tfproanno[tfproanno$Gastro.Pro.ID %in% ourpro,"Gene.Name"])[1]
  ourprophsymbol <- paste(ourprosymbol,gsub(".*_","",ourproph),sep = "_")
  #tfgastroproph_crossomecor_df[i,"Label"] <- paste(ourrnasymbol,":",ourprosymbol,sep = "")
  #tfgastroproph_crossomecor_df[(i+dim(tfgastroproph_crossomecor)[1]),"Label"] <- paste(ourrnasymbol,":",ourprophsymbol,sep = "")
  #tfgastroproph_crossomecor_df[(i+2*dim(tfgastroproph_crossomecor)[1]),"Label"] <- gsub("_",":",ourprophsymbol)
  tfgastroproph_crossomecor_df[i,"Label"] <- ourprosymbol
  tfgastroproph_crossomecor_df[(i+dim(tfgastroproph_crossomecor)[1]),"Label"] <- ourprosymbol
  tfgastroproph_crossomecor_df[(i+2*dim(tfgastroproph_crossomecor)[1]),"Label"] <- ourprosymbol
  
}

tfgastroproph_omeids <- c(rownames(tfgastroproph_rnanorm),
                          rownames(tfgastroproph_pronorm),
                          rownames(tfgastroproph_prophnorm))
tfgastroproph_crossomecor_df$SigTrain <- 0
for(i in 1:dim(tfgastroproph_rnanorm)[1]){
  tfgastroproph_crossomecor_df[i,"SigTrain"] <- (rownames(tfgastroproph_rnanorm)[i] %in% gastrornasig ) + (rownames(tfgastroproph_pronorm)[i] %in% gastroprosig)
  tfgastroproph_crossomecor_df[dim(tfgastroproph_rnanorm)[1]+i,"SigTrain"] <- (rownames(tfgastroproph_rnanorm)[i] %in% gastrornasig ) + (rownames(tfgastroproph_prophnorm)[i] %in% gastroprophsig)
  tfgastroproph_crossomecor_df[dim(tfgastroproph_rnanorm)[1]+i,"SigTrain"] <- (rownames(tfgastroproph_pronorm)[i] %in% gastroprosig ) + (rownames(tfgastroproph_prophnorm)[i] %in% gastroprophsig)
}


tfgastroproph_crossomecor_dffull <- rbind(tfgastroproph_crossomecor_df[tfgastroproph_crossomecor_df$Comparison %in% "RNA-Pro",],
                                          tfgastroproph_crossomecor_df[tfgastroproph_crossomecor_df$Comparison %in% "RNA-Phos",],
                                          tfgastroproph_crossomecor_df[tfgastroproph_crossomecor_df$Comparison %in% "Pro-Phos",])
tfgastroproph_crossomecor_dffull[tfgastroproph_crossomecor_df$P.Value > 0.1,"Label"] <- ""


tfgastroproph_crossomecor_df <- rbind(distinct(tfgastroproph_crossomecor_df[tfgastroproph_crossomecor_df$Comparison %in% "RNA-Pro",],Label,.keep_all = TRUE),
                                      tfgastroproph_crossomecor_df[tfgastroproph_crossomecor_df$Comparison %in% "RNA-Phos",],
                                      tfgastroproph_crossomecor_df[tfgastroproph_crossomecor_df$Comparison %in% "Pro-Phos",])
tfgastroproph_crossomecor_df[tfgastroproph_crossomecor_df$P.Value > 0.1,"Label"] <- ""


tfgastroproph_crossomecor_df$SigTrain <- factor(tfgastroproph_crossomecor_df$SigTrain,levels = c("0","1","2"))


png(file = "Supplemental Figure S21_021325.png",width = 20, height = 12,units = "in",res = 600)
ggscatter(tfgastroproph_crossomecor_df, x = "Correlation", y = "P.Value",facet.by = "Comparison") + scale_y_sqrt(breaks = c(0.00,0.05,0.10,0.25,0.50,0.75,1.00)) + geom_text_repel(aes(label=Label,color = SigTrain),max.overlaps = 30) + theme(strip.text = element_text(size = 25),axis.text = element_text(size = 22),axis.title = element_text(size = 25),legend.position = "none") + scale_color_manual(values = c("black","orange","red"))
dev.off()

pdf(file = "Supplemental Figure S21_021325.pdf",width = 20, height = 12)
ggscatter(tfgastroproph_crossomecor_df, x = "Correlation", y = "P.Value",facet.by = "Comparison") + scale_y_sqrt(breaks = c(0.00,0.05,0.10,0.25,0.50,0.75,1.00)) + geom_text_repel(aes(label=Label,color = SigTrain),max.overlaps = 30) + theme(strip.text = element_text(size = 25),axis.text = element_text(size = 22),axis.title = element_text(size = 25),legend.position = "none") + scale_color_manual(values = c("black","orange","red"))
dev.off()

#heart
tfheartproph_crossomecor_df <- data.frame("Correlation" = c(tfheartproph_crossomecor[,"RNA-Pro"],
                                                             tfheartproph_crossomecor[,"RNA-Phos"],
                                                             tfheartproph_crossomecor[,"Pro-Phos"]),
                                           "P.Value" = c(tfheartproph_crossomecortest_adj[,"RNA-Pro"],
                                                         tfheartproph_crossomecortest_adj[,"RNA-Phos"],
                                                         tfheartproph_crossomecortest_adj[,"Pro-Phos"]),
                                           "Comparison" = c(rep("RNA-Pro",dim(tfheartproph_crossomecortest_adj)[1]),
                                                            rep("RNA-Phos",dim(tfheartproph_crossomecortest_adj)[1]),
                                                            rep("Pro-Phos",dim(tfheartproph_crossomecortest_adj)[1])))
#tfheartproph_crossomecor_df$Log2.P.Value <- -log2(tfheartproph_crossomecor_df$Log2.P.Value)
tfheartproph_crossomecor_df$Label <- rep("",dim(tfheartproph_crossomecor_df)[1])
for(i in 1:dim(tfheartproph_crossomecor)[1]){
  ourrna <- gsub("\\..*","",rownames(tfheartproph_rnanorm))[i]
  ourpro <- sub('^([^\\.]+\\.[^\\.]+).*', '\\1', rownames(tfheartproph_pronorm))[i]
  ourproph <- rownames(tfheartproph_prophnorm)[i]
  ourrnasymbol <- enstosym[ourrna,"Symbol"]
  ourprosymbol <- toupper(tfproanno[tfproanno$Heart.Pro.ID %in% ourpro,"Gene.Name"])[1]
  ourprophsymbol <- paste(ourprosymbol,gsub(".*_","",ourproph),sep = "_")
  #tfheartproph_crossomecor_df[i,"Label"] <- paste(ourrnasymbol,":",ourprosymbol,sep = "")
  #tfheartproph_crossomecor_df[(i+dim(tfheartproph_crossomecor)[1]),"Label"] <- paste(ourrnasymbol,":",ourprophsymbol,sep = "")
  #tfheartproph_crossomecor_df[(i+2*dim(tfheartproph_crossomecor)[1]),"Label"] <- gsub("_",":",ourprophsymbol)
  tfheartproph_crossomecor_df[i,"Label"] <- ourprosymbol
  tfheartproph_crossomecor_df[(i+dim(tfheartproph_crossomecor)[1]),"Label"] <- ourprosymbol
  tfheartproph_crossomecor_df[(i+2*dim(tfheartproph_crossomecor)[1]),"Label"] <- ourprosymbol
  
}

tfheartproph_omeids <- c(rownames(tfheartproph_rnanorm),
                          rownames(tfheartproph_pronorm),
                          rownames(tfheartproph_prophnorm))
tfheartproph_crossomecor_df$SigTrain <- 0
for(i in 1:dim(tfheartproph_rnanorm)[1]){
  tfheartproph_crossomecor_df[i,"SigTrain"] <- (rownames(tfheartproph_rnanorm)[i] %in% heartrnasig ) + (rownames(tfheartproph_pronorm)[i] %in% heartprosig)
  tfheartproph_crossomecor_df[dim(tfheartproph_rnanorm)[1]+i,"SigTrain"] <- (rownames(tfheartproph_rnanorm)[i] %in% heartrnasig ) + (rownames(tfheartproph_prophnorm)[i] %in% heartprophsig)
  tfheartproph_crossomecor_df[dim(tfheartproph_rnanorm)[1]+i,"SigTrain"] <- (rownames(tfheartproph_pronorm)[i] %in% heartprosig ) + (rownames(tfheartproph_prophnorm)[i] %in% heartprophsig)
}


tfheartproph_crossomecor_dffull <- rbind(tfheartproph_crossomecor_df[tfheartproph_crossomecor_df$Comparison %in% "RNA-Pro",],
                                          tfheartproph_crossomecor_df[tfheartproph_crossomecor_df$Comparison %in% "RNA-Phos",],
                                          tfheartproph_crossomecor_df[tfheartproph_crossomecor_df$Comparison %in% "Pro-Phos",])
tfheartproph_crossomecor_dffull[tfheartproph_crossomecor_df$P.Value > 0.1,"Label"] <- ""


tfheartproph_crossomecor_df <- rbind(distinct(tfheartproph_crossomecor_df[tfheartproph_crossomecor_df$Comparison %in% "RNA-Pro",],Label,.keep_all = TRUE),
                                      tfheartproph_crossomecor_df[tfheartproph_crossomecor_df$Comparison %in% "RNA-Phos",],
                                      tfheartproph_crossomecor_df[tfheartproph_crossomecor_df$Comparison %in% "Pro-Phos",])
tfheartproph_crossomecor_df[tfheartproph_crossomecor_df$P.Value > 0.1,"Label"] <- ""


tfheartproph_crossomecor_df$SigTrain <- factor(tfheartproph_crossomecor_df$SigTrain,levels = c("0","1","2"))


png(file = "Supplemental Figure S21B_021325.png",width = 20, height = 12,units = "in",res = 600)
ggscatter(tfheartproph_crossomecor_df, x = "Correlation", y = "P.Value",facet.by = "Comparison") + scale_y_sqrt(breaks = c(0.00,0.05,0.10,0.25,0.50,0.75,1.00)) + geom_text_repel(aes(label=Label,color = SigTrain),max.overlaps = 30) + theme(strip.text = element_text(size = 25),axis.text = element_text(size = 22),axis.title = element_text(size = 25),legend.position = "none") + scale_color_manual(values = c("black","orange","red"))
dev.off()

pdf(file = "Supplemental Figure S21B_021325.pdf",width = 20, height = 12)
ggscatter(tfheartproph_crossomecor_df, x = "Correlation", y = "P.Value",facet.by = "Comparison") + scale_y_sqrt(breaks = c(0.00,0.05,0.10,0.25,0.50,0.75,1.00)) + geom_text_repel(aes(label=Label,color = SigTrain),max.overlaps = 30) + theme(strip.text = element_text(size = 25),axis.text = element_text(size = 22),axis.title = element_text(size = 25),legend.position = "none") + scale_color_manual(values = c("black","orange","red"))
dev.off()

#kidney
tfkidneyproph_crossomecor_df <- data.frame("Correlation" = c(tfkidneyproph_crossomecor[,"RNA-Pro"],
                                                             tfkidneyproph_crossomecor[,"RNA-Phos"],
                                                             tfkidneyproph_crossomecor[,"Pro-Phos"]),
                                           "P.Value" = c(tfkidneyproph_crossomecortest_adj[,"RNA-Pro"],
                                                         tfkidneyproph_crossomecortest_adj[,"RNA-Phos"],
                                                         tfkidneyproph_crossomecortest_adj[,"Pro-Phos"]),
                                           "Comparison" = c(rep("RNA-Pro",dim(tfkidneyproph_crossomecortest_adj)[1]),
                                                            rep("RNA-Phos",dim(tfkidneyproph_crossomecortest_adj)[1]),
                                                            rep("Pro-Phos",dim(tfkidneyproph_crossomecortest_adj)[1])))
#tfkidneyproph_crossomecor_df$Log2.P.Value <- -log2(tfkidneyproph_crossomecor_df$Log2.P.Value)
tfkidneyproph_crossomecor_df$Label <- rep("",dim(tfkidneyproph_crossomecor_df)[1])
for(i in 1:dim(tfkidneyproph_crossomecor)[1]){
  ourrna <- gsub("\\..*","",rownames(tfkidneyproph_rnanorm))[i]
  ourpro <- sub('^([^\\.]+\\.[^\\.]+).*', '\\1', rownames(tfkidneyproph_pronorm))[i]
  ourproph <- rownames(tfkidneyproph_prophnorm)[i]
  ourrnasymbol <- enstosym[ourrna,"Symbol"]
  ourprosymbol <- toupper(tfproanno[tfproanno$Kidney.Pro.ID %in% ourpro,"Gene.Name"])[1]
  ourprophsymbol <- paste(ourprosymbol,gsub(".*_","",ourproph),sep = "_")
  #tfkidneyproph_crossomecor_df[i,"Label"] <- paste(ourrnasymbol,":",ourprosymbol,sep = "")
  #tfkidneyproph_crossomecor_df[(i+dim(tfkidneyproph_crossomecor)[1]),"Label"] <- paste(ourrnasymbol,":",ourprophsymbol,sep = "")
  #tfkidneyproph_crossomecor_df[(i+2*dim(tfkidneyproph_crossomecor)[1]),"Label"] <- gsub("_",":",ourprophsymbol)
  tfkidneyproph_crossomecor_df[i,"Label"] <- ourprosymbol
  tfkidneyproph_crossomecor_df[(i+dim(tfkidneyproph_crossomecor)[1]),"Label"] <- ourprosymbol
  tfkidneyproph_crossomecor_df[(i+2*dim(tfkidneyproph_crossomecor)[1]),"Label"] <- ourprosymbol
  
}

tfkidneyproph_omeids <- c(rownames(tfkidneyproph_rnanorm),
                          rownames(tfkidneyproph_pronorm),
                          rownames(tfkidneyproph_prophnorm))
tfkidneyproph_crossomecor_df$SigTrain <- 0
for(i in 1:dim(tfkidneyproph_rnanorm)[1]){
  tfkidneyproph_crossomecor_df[i,"SigTrain"] <- (rownames(tfkidneyproph_rnanorm)[i] %in% kidneyrnasig ) + (rownames(tfkidneyproph_pronorm)[i] %in% kidneyprosig)
  tfkidneyproph_crossomecor_df[dim(tfkidneyproph_rnanorm)[1]+i,"SigTrain"] <- (rownames(tfkidneyproph_rnanorm)[i] %in% kidneyrnasig ) + (rownames(tfkidneyproph_prophnorm)[i] %in% kidneyprophsig)
  tfkidneyproph_crossomecor_df[dim(tfkidneyproph_rnanorm)[1]+i,"SigTrain"] <- (rownames(tfkidneyproph_pronorm)[i] %in% kidneyprosig ) + (rownames(tfkidneyproph_prophnorm)[i] %in% kidneyprophsig)
}


tfkidneyproph_crossomecor_dffull <- rbind(tfkidneyproph_crossomecor_df[tfkidneyproph_crossomecor_df$Comparison %in% "RNA-Pro",],
                                          tfkidneyproph_crossomecor_df[tfkidneyproph_crossomecor_df$Comparison %in% "RNA-Phos",],
                                          tfkidneyproph_crossomecor_df[tfkidneyproph_crossomecor_df$Comparison %in% "Pro-Phos",])
tfkidneyproph_crossomecor_dffull[tfkidneyproph_crossomecor_df$P.Value > 0.1,"Label"] <- ""


tfkidneyproph_crossomecor_df <- rbind(distinct(tfkidneyproph_crossomecor_df[tfkidneyproph_crossomecor_df$Comparison %in% "RNA-Pro",],Label,.keep_all = TRUE),
                                      tfkidneyproph_crossomecor_df[tfkidneyproph_crossomecor_df$Comparison %in% "RNA-Phos",],
                                      tfkidneyproph_crossomecor_df[tfkidneyproph_crossomecor_df$Comparison %in% "Pro-Phos",])
tfkidneyproph_crossomecor_df[tfkidneyproph_crossomecor_df$P.Value > 0.1,"Label"] <- ""


tfkidneyproph_crossomecor_df$SigTrain <- factor(tfkidneyproph_crossomecor_df$SigTrain,levels = c("0","1","2"))


png(file = "Supplemental Figure S21C_021325.png",width = 20, height = 12,units = "in",res = 600)
ggscatter(tfkidneyproph_crossomecor_df, x = "Correlation", y = "P.Value",facet.by = "Comparison") + scale_y_sqrt(breaks = c(0.00,0.05,0.10,0.25,0.50,0.75,1.00)) + geom_text_repel(aes(label=Label,color = SigTrain),max.overlaps = 30) + theme(strip.text = element_text(size = 25),axis.text = element_text(size = 22),axis.title = element_text(size = 25),legend.position = "none") + scale_color_manual(values = c("black","orange","red"))
dev.off()

pdf(file = "Supplemental Figure S21C_021325.pdf",width = 20, height = 12)
ggscatter(tfkidneyproph_crossomecor_df, x = "Correlation", y = "P.Value",facet.by = "Comparison") + scale_y_sqrt(breaks = c(0.00,0.05,0.10,0.25,0.50,0.75,1.00)) + geom_text_repel(aes(label=Label,color = SigTrain),max.overlaps = 30) + theme(strip.text = element_text(size = 25),axis.text = element_text(size = 22),axis.title = element_text(size = 25),legend.position = "none") + scale_color_manual(values = c("black","orange","red"))
dev.off()

#liver
tfliverproph_crossomecor_df <- data.frame("Correlation" = c(tfliverproph_crossomecor[,"RNA-Pro"],
                                                             tfliverproph_crossomecor[,"RNA-Phos"],
                                                             tfliverproph_crossomecor[,"Pro-Phos"]),
                                           "P.Value" = c(tfliverproph_crossomecortest_adj[,"RNA-Pro"],
                                                         tfliverproph_crossomecortest_adj[,"RNA-Phos"],
                                                         tfliverproph_crossomecortest_adj[,"Pro-Phos"]),
                                           "Comparison" = c(rep("RNA-Pro",dim(tfliverproph_crossomecortest_adj)[1]),
                                                            rep("RNA-Phos",dim(tfliverproph_crossomecortest_adj)[1]),
                                                            rep("Pro-Phos",dim(tfliverproph_crossomecortest_adj)[1])))
#tfliverproph_crossomecor_df$Log2.P.Value <- -log2(tfliverproph_crossomecor_df$Log2.P.Value)
tfliverproph_crossomecor_df$Label <- rep("",dim(tfliverproph_crossomecor_df)[1])
for(i in 1:dim(tfliverproph_crossomecor)[1]){
  ourrna <- gsub("\\..*","",rownames(tfliverproph_rnanorm))[i]
  ourpro <- sub('^([^\\.]+\\.[^\\.]+).*', '\\1', rownames(tfliverproph_pronorm))[i]
  ourproph <- rownames(tfliverproph_prophnorm)[i]
  ourrnasymbol <- enstosym[ourrna,"Symbol"]
  ourprosymbol <- toupper(tfproanno[tfproanno$Liver.Pro.ID %in% ourpro,"Gene.Name"])[1]
  ourprophsymbol <- paste(ourprosymbol,gsub(".*_","",ourproph),sep = "_")
  #tfliverproph_crossomecor_df[i,"Label"] <- paste(ourrnasymbol,":",ourprosymbol,sep = "")
  #tfliverproph_crossomecor_df[(i+dim(tfliverproph_crossomecor)[1]),"Label"] <- paste(ourrnasymbol,":",ourprophsymbol,sep = "")
  #tfliverproph_crossomecor_df[(i+2*dim(tfliverproph_crossomecor)[1]),"Label"] <- gsub("_",":",ourprophsymbol)
  tfliverproph_crossomecor_df[i,"Label"] <- ourprosymbol
  tfliverproph_crossomecor_df[(i+dim(tfliverproph_crossomecor)[1]),"Label"] <- ourprosymbol
  tfliverproph_crossomecor_df[(i+2*dim(tfliverproph_crossomecor)[1]),"Label"] <- ourprosymbol
  
}

tfliverproph_omeids <- c(rownames(tfliverproph_rnanorm),
                          rownames(tfliverproph_pronorm),
                          rownames(tfliverproph_prophnorm))
tfliverproph_crossomecor_df$SigTrain <- 0
for(i in 1:dim(tfliverproph_rnanorm)[1]){
  tfliverproph_crossomecor_df[i,"SigTrain"] <- (rownames(tfliverproph_rnanorm)[i] %in% liverrnasig ) + (rownames(tfliverproph_pronorm)[i] %in% liverprosig)
  tfliverproph_crossomecor_df[dim(tfliverproph_rnanorm)[1]+i,"SigTrain"] <- (rownames(tfliverproph_rnanorm)[i] %in% liverrnasig ) + (rownames(tfliverproph_prophnorm)[i] %in% liverprophsig)
  tfliverproph_crossomecor_df[dim(tfliverproph_rnanorm)[1]+i,"SigTrain"] <- (rownames(tfliverproph_pronorm)[i] %in% liverprosig ) + (rownames(tfliverproph_prophnorm)[i] %in% liverprophsig)
}


tfliverproph_crossomecor_dffull <- rbind(tfliverproph_crossomecor_df[tfliverproph_crossomecor_df$Comparison %in% "RNA-Pro",],
                                          tfliverproph_crossomecor_df[tfliverproph_crossomecor_df$Comparison %in% "RNA-Phos",],
                                          tfliverproph_crossomecor_df[tfliverproph_crossomecor_df$Comparison %in% "Pro-Phos",])
tfliverproph_crossomecor_dffull[tfliverproph_crossomecor_df$P.Value > 0.1,"Label"] <- ""


tfliverproph_crossomecor_df <- rbind(distinct(tfliverproph_crossomecor_df[tfliverproph_crossomecor_df$Comparison %in% "RNA-Pro",],Label,.keep_all = TRUE),
                                      tfliverproph_crossomecor_df[tfliverproph_crossomecor_df$Comparison %in% "RNA-Phos",],
                                      tfliverproph_crossomecor_df[tfliverproph_crossomecor_df$Comparison %in% "Pro-Phos",])
tfliverproph_crossomecor_df[tfliverproph_crossomecor_df$P.Value > 0.1,"Label"] <- ""


tfliverproph_crossomecor_df$SigTrain <- factor(tfliverproph_crossomecor_df$SigTrain,levels = c("0","1","2"))


png(file = "Supplemental Figure S21D_021325.png",width = 20, height = 12,units = "in",res = 600)
ggscatter(tfliverproph_crossomecor_df, x = "Correlation", y = "P.Value",facet.by = "Comparison") + scale_y_sqrt(breaks = c(0.00,0.05,0.10,0.25,0.50,0.75,1.00)) + geom_text_repel(aes(label=Label,color = SigTrain),max.overlaps = 30) + theme(strip.text = element_text(size = 25),axis.text = element_text(size = 22),axis.title = element_text(size = 25),legend.position = "none") + scale_color_manual(values = c("black","orange","red"))
dev.off()

pdf(file = "Supplemental Figure S21D_021325.pdf",width = 20, height = 12)
ggscatter(tfliverproph_crossomecor_df, x = "Correlation", y = "P.Value",facet.by = "Comparison") + scale_y_sqrt(breaks = c(0.00,0.05,0.10,0.25,0.50,0.75,1.00)) + geom_text_repel(aes(label=Label,color = SigTrain),max.overlaps = 30) + theme(strip.text = element_text(size = 25),axis.text = element_text(size = 22),axis.title = element_text(size = 25),legend.position = "none") + scale_color_manual(values = c("black","orange","red"))
dev.off()


#lung
tflungproph_crossomecor_df <- data.frame("Correlation" = c(tflungproph_crossomecor[,"RNA-Pro"],
                                                             tflungproph_crossomecor[,"RNA-Phos"],
                                                             tflungproph_crossomecor[,"Pro-Phos"]),
                                           "P.Value" = c(tflungproph_crossomecortest_adj[,"RNA-Pro"],
                                                         tflungproph_crossomecortest_adj[,"RNA-Phos"],
                                                         tflungproph_crossomecortest_adj[,"Pro-Phos"]),
                                           "Comparison" = c(rep("RNA-Pro",dim(tflungproph_crossomecortest_adj)[1]),
                                                            rep("RNA-Phos",dim(tflungproph_crossomecortest_adj)[1]),
                                                            rep("Pro-Phos",dim(tflungproph_crossomecortest_adj)[1])))
#tflungproph_crossomecor_df$Log2.P.Value <- -log2(tflungproph_crossomecor_df$Log2.P.Value)
tflungproph_crossomecor_df$Label <- rep("",dim(tflungproph_crossomecor_df)[1])
for(i in 1:dim(tflungproph_crossomecor)[1]){
  ourrna <- gsub("\\..*","",rownames(tflungproph_rnanorm))[i]
  ourpro <- sub('^([^\\.]+\\.[^\\.]+).*', '\\1', rownames(tflungproph_pronorm))[i]
  ourproph <- rownames(tflungproph_prophnorm)[i]
  ourrnasymbol <- enstosym[ourrna,"Symbol"]
  ourprosymbol <- toupper(tfproanno[tfproanno$Lung.Pro.ID %in% ourpro,"Gene.Name"])[1]
  ourprophsymbol <- paste(ourprosymbol,gsub(".*_","",ourproph),sep = "_")
  #tflungproph_crossomecor_df[i,"Label"] <- paste(ourrnasymbol,":",ourprosymbol,sep = "")
  #tflungproph_crossomecor_df[(i+dim(tflungproph_crossomecor)[1]),"Label"] <- paste(ourrnasymbol,":",ourprophsymbol,sep = "")
  #tflungproph_crossomecor_df[(i+2*dim(tflungproph_crossomecor)[1]),"Label"] <- gsub("_",":",ourprophsymbol)
  tflungproph_crossomecor_df[i,"Label"] <- ourprosymbol
  tflungproph_crossomecor_df[(i+dim(tflungproph_crossomecor)[1]),"Label"] <- ourprosymbol
  tflungproph_crossomecor_df[(i+2*dim(tflungproph_crossomecor)[1]),"Label"] <- ourprosymbol
  
}

tflungproph_omeids <- c(rownames(tflungproph_rnanorm),
                          rownames(tflungproph_pronorm),
                          rownames(tflungproph_prophnorm))
tflungproph_crossomecor_df$SigTrain <- 0
for(i in 1:dim(tflungproph_rnanorm)[1]){
  tflungproph_crossomecor_df[i,"SigTrain"] <- (rownames(tflungproph_rnanorm)[i] %in% lungrnasig ) + (rownames(tflungproph_pronorm)[i] %in% lungprosig)
  tflungproph_crossomecor_df[dim(tflungproph_rnanorm)[1]+i,"SigTrain"] <- (rownames(tflungproph_rnanorm)[i] %in% lungrnasig ) + (rownames(tflungproph_prophnorm)[i] %in% lungprophsig)
  tflungproph_crossomecor_df[dim(tflungproph_rnanorm)[1]+i,"SigTrain"] <- (rownames(tflungproph_pronorm)[i] %in% lungprosig ) + (rownames(tflungproph_prophnorm)[i] %in% lungprophsig)
}


tflungproph_crossomecor_dffull <- rbind(tflungproph_crossomecor_df[tflungproph_crossomecor_df$Comparison %in% "RNA-Pro",],
                                          tflungproph_crossomecor_df[tflungproph_crossomecor_df$Comparison %in% "RNA-Phos",],
                                          tflungproph_crossomecor_df[tflungproph_crossomecor_df$Comparison %in% "Pro-Phos",])
tflungproph_crossomecor_dffull[tflungproph_crossomecor_df$P.Value > 0.1,"Label"] <- ""


tflungproph_crossomecor_df <- rbind(distinct(tflungproph_crossomecor_df[tflungproph_crossomecor_df$Comparison %in% "RNA-Pro",],Label,.keep_all = TRUE),
                                      tflungproph_crossomecor_df[tflungproph_crossomecor_df$Comparison %in% "RNA-Phos",],
                                      tflungproph_crossomecor_df[tflungproph_crossomecor_df$Comparison %in% "Pro-Phos",])
tflungproph_crossomecor_df[tflungproph_crossomecor_df$P.Value > 0.1,"Label"] <- ""


tflungproph_crossomecor_df$SigTrain <- factor(tflungproph_crossomecor_df$SigTrain,levels = c("0","1","2"))


png(file = "Supplemental Figure S21E_021325.png",width = 20, height = 12,units = "in",res = 600)
ggscatter(tflungproph_crossomecor_df, x = "Correlation", y = "P.Value",facet.by = "Comparison") + scale_y_sqrt(breaks = c(0.00,0.05,0.10,0.25,0.50,0.75,1.00)) + geom_text_repel(aes(label=Label,color = SigTrain),max.overlaps = 30) + theme(strip.text = element_text(size = 25),axis.text = element_text(size = 22),axis.title = element_text(size = 25),legend.position = "none") + scale_color_manual(values = c("black","orange","red"))
dev.off()

pdf(file = "Supplemental Figure S21E_021325.pdf",width = 20, height = 12)
ggscatter(tflungproph_crossomecor_df, x = "Correlation", y = "P.Value",facet.by = "Comparison") + scale_y_sqrt(breaks = c(0.00,0.05,0.10,0.25,0.50,0.75,1.00)) + geom_text_repel(aes(label=Label,color = SigTrain),max.overlaps = 30) + theme(strip.text = element_text(size = 25),axis.text = element_text(size = 22),axis.title = element_text(size = 25),legend.position = "none") + scale_color_manual(values = c("black","orange","red"))
dev.off()

#white
tfwhiteproph_crossomecor_df <- data.frame("Correlation" = c(tfwhiteproph_crossomecor[,"RNA-Pro"],
                                                             tfwhiteproph_crossomecor[,"RNA-Phos"],
                                                             tfwhiteproph_crossomecor[,"Pro-Phos"]),
                                           "P.Value" = c(tfwhiteproph_crossomecortest_adj[,"RNA-Pro"],
                                                         tfwhiteproph_crossomecortest_adj[,"RNA-Phos"],
                                                         tfwhiteproph_crossomecortest_adj[,"Pro-Phos"]),
                                           "Comparison" = c(rep("RNA-Pro",dim(tfwhiteproph_crossomecortest_adj)[1]),
                                                            rep("RNA-Phos",dim(tfwhiteproph_crossomecortest_adj)[1]),
                                                            rep("Pro-Phos",dim(tfwhiteproph_crossomecortest_adj)[1])))
#tfwhiteproph_crossomecor_df$Log2.P.Value <- -log2(tfwhiteproph_crossomecor_df$Log2.P.Value)
tfwhiteproph_crossomecor_df$Label <- rep("",dim(tfwhiteproph_crossomecor_df)[1])
for(i in 1:dim(tfwhiteproph_crossomecor)[1]){
  ourrna <- gsub("\\..*","",rownames(tfwhiteproph_rnanorm))[i]
  ourpro <- sub('^([^\\.]+\\.[^\\.]+).*', '\\1', rownames(tfwhiteproph_pronorm))[i]
  ourproph <- rownames(tfwhiteproph_prophnorm)[i]
  ourrnasymbol <- enstosym[ourrna,"Symbol"]
  ourprosymbol <- toupper(tfproanno[tfproanno$WhiteAd.Pro.ID %in% ourpro,"Gene.Name"])[1]
  ourprophsymbol <- paste(ourprosymbol,gsub(".*_","",ourproph),sep = "_")
  #tfwhiteproph_crossomecor_df[i,"Label"] <- paste(ourrnasymbol,":",ourprosymbol,sep = "")
  #tfwhiteproph_crossomecor_df[(i+dim(tfwhiteproph_crossomecor)[1]),"Label"] <- paste(ourrnasymbol,":",ourprophsymbol,sep = "")
  #tfwhiteproph_crossomecor_df[(i+2*dim(tfwhiteproph_crossomecor)[1]),"Label"] <- gsub("_",":",ourprophsymbol)
  tfwhiteproph_crossomecor_df[i,"Label"] <- ourprosymbol
  tfwhiteproph_crossomecor_df[(i+dim(tfwhiteproph_crossomecor)[1]),"Label"] <- ourprosymbol
  tfwhiteproph_crossomecor_df[(i+2*dim(tfwhiteproph_crossomecor)[1]),"Label"] <- ourprosymbol
  
}

tfwhiteproph_omeids <- c(rownames(tfwhiteproph_rnanorm),
                          rownames(tfwhiteproph_pronorm),
                          rownames(tfwhiteproph_prophnorm))
tfwhiteproph_crossomecor_df$SigTrain <- 0
for(i in 1:dim(tfwhiteproph_rnanorm)[1]){
  tfwhiteproph_crossomecor_df[i,"SigTrain"] <- (rownames(tfwhiteproph_rnanorm)[i] %in% whiternasig ) + (rownames(tfwhiteproph_pronorm)[i] %in% whiteprosig)
  tfwhiteproph_crossomecor_df[dim(tfwhiteproph_rnanorm)[1]+i,"SigTrain"] <- (rownames(tfwhiteproph_rnanorm)[i] %in% whiternasig ) + (rownames(tfwhiteproph_prophnorm)[i] %in% whiteprophsig)
  tfwhiteproph_crossomecor_df[dim(tfwhiteproph_rnanorm)[1]+i,"SigTrain"] <- (rownames(tfwhiteproph_pronorm)[i] %in% whiteprosig ) + (rownames(tfwhiteproph_prophnorm)[i] %in% whiteprophsig)
}


tfwhiteproph_crossomecor_dffull <- rbind(tfwhiteproph_crossomecor_df[tfwhiteproph_crossomecor_df$Comparison %in% "RNA-Pro",],
                                          tfwhiteproph_crossomecor_df[tfwhiteproph_crossomecor_df$Comparison %in% "RNA-Phos",],
                                          tfwhiteproph_crossomecor_df[tfwhiteproph_crossomecor_df$Comparison %in% "Pro-Phos",])
tfwhiteproph_crossomecor_dffull[tfwhiteproph_crossomecor_df$P.Value > 0.1,"Label"] <- ""


tfwhiteproph_crossomecor_df <- rbind(distinct(tfwhiteproph_crossomecor_df[tfwhiteproph_crossomecor_df$Comparison %in% "RNA-Pro",],Label,.keep_all = TRUE),
                                      tfwhiteproph_crossomecor_df[tfwhiteproph_crossomecor_df$Comparison %in% "RNA-Phos",],
                                      tfwhiteproph_crossomecor_df[tfwhiteproph_crossomecor_df$Comparison %in% "Pro-Phos",])
tfwhiteproph_crossomecor_df[tfwhiteproph_crossomecor_df$P.Value > 0.1,"Label"] <- ""


tfwhiteproph_crossomecor_df$SigTrain <- factor(tfwhiteproph_crossomecor_df$SigTrain,levels = c("0","1","2"))


png(file = "Supplemental Figure S21F_021325.png",width = 20, height = 12,units = "in",res = 600)
ggscatter(tfwhiteproph_crossomecor_df, x = "Correlation", y = "P.Value",facet.by = "Comparison") + scale_y_sqrt(breaks = c(0.00,0.05,0.10,0.25,0.50,0.75,1.00)) + geom_text_repel(aes(label=Label,color = SigTrain),max.overlaps = 30) + theme(strip.text = element_text(size = 25),axis.text = element_text(size = 22),axis.title = element_text(size = 25),legend.position = "none") + scale_color_manual(values = c("black","orange","red"))
dev.off()

pdf(file = "Supplemental Figure S21F_021325.pdf",width = 20, height = 12)
ggscatter(tfwhiteproph_crossomecor_df, x = "Correlation", y = "P.Value",facet.by = "Comparison") + scale_y_sqrt(breaks = c(0.00,0.05,0.10,0.25,0.50,0.75,1.00)) + geom_text_repel(aes(label=Label,color = SigTrain),max.overlaps = 30) + theme(strip.text = element_text(size = 25),axis.text = element_text(size = 22),axis.title = element_text(size = 25),legend.position = "none") + scale_color_manual(values = c("black","orange","red"))
dev.off()


temp1 <- tfgastroproph_rnanorm
colnames(temp1) <- paste("RNA",colnames(temp1),sep = "_")
temp1 <- temp1[apply(tfgastroproph_crossomecortest_adj,1,min) < 0.1,]
temp2 <- tfgastroproph_pronorm
colnames(temp2) <- paste("Pro",colnames(temp2),sep = "_")
temp2 <- temp2[apply(tfgastroproph_crossomecortest_adj,1,min) < 0.1,]
temp3 <- tfgastroproph_prophnorm
colnames(temp3) <- paste("Phos",colnames(temp3),sep = "_")
temp3 <- temp3[apply(tfgastroproph_crossomecortest_adj,1,min) < 0.1,]
tfgastroproph_sigcore_allomenorm <- cbind(t(scale(t(temp1))),t(scale(t(temp2))),t(scale(t(temp3))))
tfgastroproph_sigcore_allomenorm_metadf <- data.frame(row.names = colnames(tfgastroproph_sigcore_allomenorm),
                                                      "Ome" = gsub("_.*","",colnames(tfgastroproph_sigcore_allomenorm)),
                                                      "pid" = gsub(".*_","",colnames(tfgastroproph_sigcore_allomenorm)),
                                                      "Sex" = rep("Female",length(colnames(tfgastroproph_sigcore_allomenorm))),
                                                      "Group" = rep("control",length(colnames(tfgastroproph_sigcore_allomenorm))))
for(i in 1:dim(tfgastroproph_sigcore_allomenorm_metadf)[1]){
  ourpid <- tfgastroproph_sigcore_allomenorm_metadf$pid[i]
  tfgastroproph_sigcore_allomenorm_metadf[i,"Group"] <- pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourpid,"pass1bf0011"][1]
  tfgastroproph_sigcore_allomenorm_metadf[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourpid,"pass1bf0027"][1]]
}
tfgastroproph_sigcore_allomenorm_metadf[tfgastroproph_sigcore_allomenorm_metadf$Group %in% "Eight-week program Control Group","Group"] <- "control"
tfgastroproph_sigcore_allomenorm_metadf[tfgastroproph_sigcore_allomenorm_metadf$Group %in% "One-week program","Group"] <- "1w"
tfgastroproph_sigcore_allomenorm_metadf[tfgastroproph_sigcore_allomenorm_metadf$Group %in% "Two-week program","Group"] <- "2w"
tfgastroproph_sigcore_allomenorm_metadf[tfgastroproph_sigcore_allomenorm_metadf$Group %in% "Four-week program","Group"] <- "4w"
tfgastroproph_sigcore_allomenorm_metadf[tfgastroproph_sigcore_allomenorm_metadf$Group %in% "Eight-week program Training Group","Group"] <- "8w"
tfgastroproph_sigcore_allomenorm_metadf$Ome <- factor(tfgastroproph_sigcore_allomenorm_metadf$Ome,levels = c("RNA","Pro","Phos"))
tfgastroproph_sigcore_allomenorm_metadf$Group <- factor(tfgastroproph_sigcore_allomenorm_metadf$Group,levels = c("control","1w","2w","4w","8w"))

tfgastroproph_sigcore_allomenorm_metadf <- tfgastroproph_sigcore_allomenorm_metadf[order(tfgastroproph_sigcore_allomenorm_metadf$Ome,tfgastroproph_sigcore_allomenorm_metadf$Sex,tfgastroproph_sigcore_allomenorm_metadf$Group),]
tfgastroproph_sigcore_allomenorm <- tfgastroproph_sigcore_allomenorm[,rownames(tfgastroproph_sigcore_allomenorm_metadf)]

save.image("Figure5_021325.RData")

# We want to identify the TFs that show coordinated changes in RNA-Pro-Phos - then we want to select the
# ones that are significantly changing in response to exercise. Then we want to look at their gene targets

# gastro
tfgastrosigrows <- as.numeric(rownames(tfgastroproph_crossomecor_dffull[tfgastroproph_crossomecor_dffull$P.Value < 0.1,])) %% length(rownames(tfgastroproph_prophnorm))
tfgastrosigrows[tfgastrosigrows == 0] <- length(rownames(tfgastroproph_prophnorm))

tfgastroproph_sigtfphossites <- unique(rownames(tfgastroproph_prophnorm)[tfgastrosigrows])

tfgastroproph_sigtfphossites_df <- data.frame(row.names = tfgastroproph_sigtfphossites,
                                              "Phos_ID" = tfgastroproph_sigtfphossites,
                                              "Pro_ID" = rep("",length(tfgastroproph_sigtfphossites)),
                                              "RNA_ID" = rep("",length(tfgastroproph_sigtfphossites)))
for(i in 1:length(tfgastroproph_sigtfphossites)){
  ourphos <- tfgastroproph_sigtfphossites[i]
  ourpro <- sub('^([^_]+_[^_]+).*', '\\1',ourphos)
  ourgene <- tfproanno[tfproanno$Gastro.Pro.ID %in% ourpro,"Ensembl"][1]
  tfgastroproph_sigtfphossites_df[i,"Pro_ID"] <- ourpro
  tfgastroproph_sigtfphossites_df[i,"RNA_ID"] <- ourgene
}

gastropidmeta <- data.frame(row.names = colnames(tfgastroproph_prophnorm),
                            "pid" = colnames(tfgastroproph_prophnorm),
                            "Sex" = rep("Female",length(colnames(tfgastroproph_prophnorm))),
                            "Group" = rep("control",length(colnames(tfgastroproph_prophnorm))))
for(i in 1:dim(gastropidmeta)[1]){
  ourpid <- rownames(gastropidmeta)[i]
  gastropidmeta[i,"Group"] <- pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourpid,"pass1bf0011"][1]
  gastropidmeta[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourpid,"pass1bf0027"][1]]
}
gastropidmeta[gastropidmeta$Group %in% "One-week program","Group"] <- "1w"
gastropidmeta[gastropidmeta$Group %in% "Two-week program","Group"] <- "2w"
gastropidmeta[gastropidmeta$Group %in% "Four-week program","Group"] <- "4w"
gastropidmeta[gastropidmeta$Group %in% "Eight-week program Training Group","Group"] <- "8w"
gastropidmeta[gastropidmeta$Group %in% "Eight-week program Control Group","Group"] <- "control"

ourrna <- tfgastroproph_sigtfphossites_df$RNA_ID[1]
ourpro <- tfgastroproph_sigtfphossites_df$Pro_ID[1]
ourdf <- data.frame(row.names = colnames(tfgastroproph_rnanorm),
                    "RNA" = t(tfgastroproph_rnanorm[ourrna,]),
                    "Pro" = t(tfgastroproph_pronorm[ourpro,]),
                    "Sex" = gastropidmeta[colnames(tfgastroproph_rnanorm),"Sex"],
                    "Group" = gastropidmeta[colnames(tfgastroproph_rnanorm),"Group"])
colnames(ourdf) <- c("RNA","Pro","Sex","Group")


# heart
tfheartsigrows <- as.numeric(rownames(tfheartproph_crossomecor_dffull[tfheartproph_crossomecor_dffull$P.Value < 0.1,])) %% length(rownames(tfheartproph_prophnorm))
tfheartsigrows[tfheartsigrows == 0] <- length(rownames(tfheartproph_prophnorm))

tfheartproph_sigtfphossites <- unique(rownames(tfheartproph_prophnorm)[tfheartsigrows])
tfheartproph_sigtfphossites_df <- data.frame(row.names = tfheartproph_sigtfphossites,
                                        "Phos_ID" = tfheartproph_sigtfphossites,
                                        "Pro_ID" = rep("",length(tfheartproph_sigtfphossites)),
                                        "RNA_ID" = rep("",length(tfheartproph_sigtfphossites)))
for(i in 1:length(tfheartproph_sigtfphossites)){
  ourphos <- tfheartproph_sigtfphossites[i]
  ourpro <- sub('^([^_]+_[^_]+).*', '\\1',ourphos)
  ourgene <- tfproanno[tfproanno$Heart.Pro.ID %in% ourpro,"Ensembl"][1]
  tfheartproph_sigtfphossites_df[i,"Pro_ID"] <- ourpro
  tfheartproph_sigtfphossites_df[i,"RNA_ID"] <- ourgene
}

heartpidmeta <- data.frame(row.names = colnames(tfheartproph_prophnorm),
                            "pid" = colnames(tfheartproph_prophnorm),
                            "Sex" = rep("Female",length(colnames(tfheartproph_prophnorm))),
                            "Group" = rep("control",length(colnames(tfheartproph_prophnorm))))
for(i in 1:dim(heartpidmeta)[1]){
  ourpid <- rownames(heartpidmeta)[i]
  heartpidmeta[i,"Group"] <- pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourpid,"pass1bf0011"][1]
  heartpidmeta[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourpid,"pass1bf0027"][1]]
}
heartpidmeta[heartpidmeta$Group %in% "One-week program","Group"] <- "1w"
heartpidmeta[heartpidmeta$Group %in% "Two-week program","Group"] <- "2w"
heartpidmeta[heartpidmeta$Group %in% "Four-week program","Group"] <- "4w"
heartpidmeta[heartpidmeta$Group %in% "Eight-week program Training Group","Group"] <- "8w"
heartpidmeta[heartpidmeta$Group %in% "Eight-week program Control Group","Group"] <- "control"

ourrna <- tfheartproph_sigtfphossites_df$RNA_ID[18]
#ourpro <- tfheartproph_sigtfphossites_df$Pro_ID[18]
ourproph <- tfheartproph_sigtfphossites_df$Phos_ID[18]
ourdf <- data.frame(row.names = colnames(tfheartproph_rnanorm),
                    "RNA" = t(tfheartproph_rnanorm[ourrna,]),
                    "Phos" = t(tfheartproph_prophnorm[ourproph,]),
                    "Sex" = heartpidmeta[colnames(tfheartproph_rnanorm),"Sex"],
                    "Group" = heartpidmeta[colnames(tfheartproph_rnanorm),"Group"])
colnames(ourdf) <- c("RNA","Phos","Sex","Group")

#png(file = paste("HEART_",enstosym[ourrna,"Symbol"],"_",paste(toupper(enstosym[ourrna,"Symbol"]),"_",gsub(".*_","",ourproph),sep = ""),".png",sep = ""),width = 6,height = 6,units = "in",res = 600)
#ggscatter(ourdf, x = "RNA", y = "Phos",color = "Group",shape = "Sex",size = 4) + scale_color_manual(values = ann_cols$Group) + ggtitle(paste("HEART: ",enstosym[ourrna,"Symbol"],":",paste(toupper(enstosym[ourrna,"Symbol"]),"_",gsub(".*_","",ourproph),sep = "")," Correlation: ",substr(toString(cor(ourdf$RNA,ourdf$Phos)),1,5),sep = "")) + xlab("RNA") + ylab("Phosphoproteome") + theme(legend.box = "vertical",legend.margin = margin())
#dev.off()

# kidney
tfkidneysigrows <- as.numeric(rownames(tfkidneyproph_crossomecor_dffull[tfkidneyproph_crossomecor_dffull$P.Value < 0.1,])) %% length(rownames(tfkidneyproph_prophnorm))
tfkidneysigrows[tfkidneysigrows == 0] <- length(rownames(tfkidneyproph_prophnorm))

tfkidneyproph_sigtfphossites <- unique(rownames(tfkidneyproph_prophnorm)[tfkidneysigrows])

tfkidneyproph_sigtfphossites_df <- data.frame(row.names = tfkidneyproph_sigtfphossites,
                                              "Phos_ID" = tfkidneyproph_sigtfphossites,
                                              "Pro_ID" = rep("",length(tfkidneyproph_sigtfphossites)),
                                              "RNA_ID" = rep("",length(tfkidneyproph_sigtfphossites)))
for(i in 1:length(tfkidneyproph_sigtfphossites)){
  ourphos <- tfkidneyproph_sigtfphossites[i]
  ourpro <- sub('^([^_]+_[^_]+).*', '\\1',ourphos)
  ourgene <- tfproanno[tfproanno$Kidney.Pro.ID %in% ourpro,"Ensembl"][1]
  tfkidneyproph_sigtfphossites_df[i,"Pro_ID"] <- ourpro
  tfkidneyproph_sigtfphossites_df[i,"RNA_ID"] <- ourgene
}

kidneypidmeta <- data.frame(row.names = colnames(tfkidneyproph_prophnorm),
                            "pid" = colnames(tfkidneyproph_prophnorm),
                            "Sex" = rep("Female",length(colnames(tfkidneyproph_prophnorm))),
                            "Group" = rep("control",length(colnames(tfkidneyproph_prophnorm))))
for(i in 1:dim(kidneypidmeta)[1]){
  ourpid <- rownames(kidneypidmeta)[i]
  kidneypidmeta[i,"Group"] <- pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourpid,"pass1bf0011"][1]
  kidneypidmeta[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourpid,"pass1bf0027"][1]]
}
kidneypidmeta[kidneypidmeta$Group %in% "One-week program","Group"] <- "1w"
kidneypidmeta[kidneypidmeta$Group %in% "Two-week program","Group"] <- "2w"
kidneypidmeta[kidneypidmeta$Group %in% "Four-week program","Group"] <- "4w"
kidneypidmeta[kidneypidmeta$Group %in% "Eight-week program Training Group","Group"] <- "8w"
kidneypidmeta[kidneypidmeta$Group %in% "Eight-week program Control Group","Group"] <- "control"

ourrna <- tfkidneyproph_sigtfphossites_df$RNA_ID[1]
ourpro <- tfkidneyproph_sigtfphossites_df$Pro_ID[1]
ourdf <- data.frame(row.names = colnames(tfkidneyproph_rnanorm),
                    "RNA" = t(tfkidneyproph_rnanorm[ourrna,]),
                    "Pro" = t(tfkidneyproph_pronorm[ourpro,]),
                    "Sex" = kidneypidmeta[colnames(tfkidneyproph_rnanorm),"Sex"],
                    "Group" = kidneypidmeta[colnames(tfkidneyproph_rnanorm),"Group"])
colnames(ourdf) <- c("RNA","Pro","Sex","Group")


# liver
tfliversigrows <- as.numeric(rownames(tfliverproph_crossomecor_dffull[tfliverproph_crossomecor_dffull$P.Value < 0.1,])) %% length(rownames(tfliverproph_prophnorm))
tfliversigrows[tfliversigrows == 0] <- length(rownames(tfliverproph_prophnorm))

tfliverproph_sigtfphossites <- unique(rownames(tfliverproph_prophnorm)[tfliversigrows])

tfliverproph_sigtfphossites_df <- data.frame(row.names = tfliverproph_sigtfphossites,
                                              "Phos_ID" = tfliverproph_sigtfphossites,
                                              "Pro_ID" = rep("",length(tfliverproph_sigtfphossites)),
                                              "RNA_ID" = rep("",length(tfliverproph_sigtfphossites)))
for(i in 1:length(tfliverproph_sigtfphossites)){
  ourphos <- tfliverproph_sigtfphossites[i]
  ourpro <- sub('^([^_]+_[^_]+).*', '\\1',ourphos)
  ourgene <- tfproanno[tfproanno$Liver.Pro.ID %in% ourpro,"Ensembl"][1]
  tfliverproph_sigtfphossites_df[i,"Pro_ID"] <- ourpro
  tfliverproph_sigtfphossites_df[i,"RNA_ID"] <- ourgene
}

liverpidmeta <- data.frame(row.names = colnames(tfliverproph_prophnorm),
                            "pid" = colnames(tfliverproph_prophnorm),
                            "Sex" = rep("Female",length(colnames(tfliverproph_prophnorm))),
                            "Group" = rep("control",length(colnames(tfliverproph_prophnorm))))
for(i in 1:dim(liverpidmeta)[1]){
  ourpid <- rownames(liverpidmeta)[i]
  liverpidmeta[i,"Group"] <- pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourpid,"pass1bf0011"][1]
  liverpidmeta[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourpid,"pass1bf0027"][1]]
}
liverpidmeta[liverpidmeta$Group %in% "One-week program","Group"] <- "1w"
liverpidmeta[liverpidmeta$Group %in% "Two-week program","Group"] <- "2w"
liverpidmeta[liverpidmeta$Group %in% "Four-week program","Group"] <- "4w"
liverpidmeta[liverpidmeta$Group %in% "Eight-week program Training Group","Group"] <- "8w"
liverpidmeta[liverpidmeta$Group %in% "Eight-week program Control Group","Group"] <- "control"

ourrna <- tfliverproph_sigtfphossites_df$RNA_ID[29]
ourpro <- tfliverproph_sigtfphossites_df$Pro_ID[29]
#ourproph <- tfliverproph_sigtfphossites_df$Phos_ID[29]
ourdf <- data.frame(row.names = colnames(tfliverproph_rnanorm),
                    "RNA" = t(tfliverproph_rnanorm[ourrna,]),
                    "Pro" = t(tfliverproph_pronorm[ourpro,]),
                    "Sex" = liverpidmeta[colnames(tfliverproph_rnanorm),"Sex"],
                    "Group" = liverpidmeta[colnames(tfliverproph_rnanorm),"Group"])
colnames(ourdf) <- c("RNA","Pro","Sex","Group")

#png(file = paste("LIVER_",enstosym[ourrna,"Symbol"],"_",toupper(enstosym[ourrna,"Symbol"]),".png",sep = ""),width = 6,height = 6,units = "in",res = 600)
#ggscatter(ourdf, x = "RNA", y = "Pro",color = "Group",shape = "Sex",size = 4) + scale_color_manual(values = ann_cols$Group) + ggtitle(paste("LIVER: ",enstosym[ourrna,"Symbol"],":",toupper(enstosym[ourrna,"Symbol"])," Correlation: ",substr(toString(cor(ourdf$RNA,ourdf$Pro)),1,5),sep = "")) + xlab("RNA") + ylab("Proteome") + theme(legend.box = "vertical",legend.margin = margin())
#dev.off()


ourrna <- tfliverproph_sigtfphossites_df$RNA_ID[25]
ourpro <- tfliverproph_sigtfphossites_df$Pro_ID[25]
#ourproph <- tfliverproph_sigtfphossites_df$Phos_ID[25]
ourdf <- data.frame(row.names = colnames(tfliverproph_rnanorm),
                    "RNA" = t(tfliverproph_rnanorm[ourrna,]),
                    "Pro" = t(tfliverproph_pronorm[ourpro,]),
                    "Sex" = liverpidmeta[colnames(tfliverproph_rnanorm),"Sex"],
                    "Group" = liverpidmeta[colnames(tfliverproph_rnanorm),"Group"])
colnames(ourdf) <- c("RNA","Pro","Sex","Group")

#png(file = paste("LIVER_",enstosym[ourrna,"Symbol"],"_",toupper(enstosym[ourrna,"Symbol"]),".png",sep = ""),width = 6,height = 6,units = "in",res = 600)
#ggscatter(ourdf, x = "RNA", y = "Pro",color = "Group",shape = "Sex",size = 4) + scale_color_manual(values = ann_cols$Group) + ggtitle(paste("LIVER: ",enstosym[ourrna,"Symbol"],":",toupper(enstosym[ourrna,"Symbol"])," Correlation: ",substr(toString(cor(ourdf$RNA,ourdf$Pro)),1,5),sep = "")) + xlab("RNA") + ylab("Proteome") + theme(legend.box = "vertical",legend.margin = margin())
#dev.off()


# lung
tflungsigrows <- as.numeric(rownames(tflungproph_crossomecor_dffull[tflungproph_crossomecor_dffull$P.Value < 0.1,])) %% length(rownames(tflungproph_prophnorm))
tflungsigrows[tflungsigrows == 0] <- length(rownames(tflungproph_prophnorm))

tflungproph_sigtfphossites <- unique(rownames(tflungproph_prophnorm)[tflungsigrows])

tflungproph_sigtfphossites_df <- data.frame(row.names = tflungproph_sigtfphossites,
                                              "Phos_ID" = tflungproph_sigtfphossites,
                                              "Pro_ID" = rep("",length(tflungproph_sigtfphossites)),
                                              "RNA_ID" = rep("",length(tflungproph_sigtfphossites)))
for(i in 1:length(tflungproph_sigtfphossites)){
  ourphos <- tflungproph_sigtfphossites[i]
  ourpro <- sub('^([^_]+_[^_]+).*', '\\1',ourphos)
  ourgene <- tfproanno[tfproanno$Lung.Pro.ID %in% ourpro,"Ensembl"][1]
  tflungproph_sigtfphossites_df[i,"Pro_ID"] <- ourpro
  tflungproph_sigtfphossites_df[i,"RNA_ID"] <- ourgene
}

lungpidmeta <- data.frame(row.names = colnames(tflungproph_prophnorm),
                            "pid" = colnames(tflungproph_prophnorm),
                            "Sex" = rep("Female",length(colnames(tflungproph_prophnorm))),
                            "Group" = rep("control",length(colnames(tflungproph_prophnorm))))
for(i in 1:dim(lungpidmeta)[1]){
  ourpid <- rownames(lungpidmeta)[i]
  lungpidmeta[i,"Group"] <- pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourpid,"pass1bf0011"][1]
  lungpidmeta[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourpid,"pass1bf0027"][1]]
}
lungpidmeta[lungpidmeta$Group %in% "One-week program","Group"] <- "1w"
lungpidmeta[lungpidmeta$Group %in% "Two-week program","Group"] <- "2w"
lungpidmeta[lungpidmeta$Group %in% "Four-week program","Group"] <- "4w"
lungpidmeta[lungpidmeta$Group %in% "Eight-week program Training Group","Group"] <- "8w"
lungpidmeta[lungpidmeta$Group %in% "Eight-week program Control Group","Group"] <- "control"

ourrna <- "ENSRNOG00000014079"
ourpro <- "NP_001012226.1"
#ourproph <- tflungproph_sigtfphossites_df$Phos_ID[25]
ourdf <- data.frame(row.names = colnames(tflungproph_rnanorm),
                    "RNA" = t(tflungproph_rnanorm[ourrna,]),
                    "Pro" = t(tflungproph_pronorm[ourpro,]),
                    "Sex" = lungpidmeta[colnames(tflungproph_rnanorm),"Sex"],
                    "Group" = lungpidmeta[colnames(tflungproph_rnanorm),"Group"])
colnames(ourdf) <- c("RNA","Pro","Sex","Group")

#png(file = "LUNG_Stat4_STAT4.png",width = 6,height = 6,units = "in",res = 600)
#ggscatter(ourdf, x = "RNA", y = "Pro",color = "Group",shape = "Sex",size = 4) + scale_color_manual(values = ann_cols$Group) + ggtitle(paste("LUNG: Stat4:STAT4 Correlation: ",substr(toString(cor(ourdf$RNA,ourdf$Pro)),1,5),sep = "")) + xlab("RNA") + ylab("Proteome") + theme(legend.box = "vertical",legend.margin = margin())
#dev.off()

# white
tfwhitesigrows <- as.numeric(rownames(tfwhiteproph_crossomecor_dffull[tfwhiteproph_crossomecor_dffull$P.Value < 0.1,])) %% length(rownames(tfwhiteproph_prophnorm))
tfwhitesigrows[tfwhitesigrows == 0] <- length(rownames(tfwhiteproph_prophnorm))

tfwhiteproph_sigtfphossites <- unique(rownames(tfwhiteproph_prophnorm)[tfwhitesigrows])

tfwhiteproph_sigtfphossites_df <- data.frame(row.names = tfwhiteproph_sigtfphossites,
                                              "Phos_ID" = tfwhiteproph_sigtfphossites,
                                              "Pro_ID" = rep("",length(tfwhiteproph_sigtfphossites)),
                                              "RNA_ID" = rep("",length(tfwhiteproph_sigtfphossites)))
for(i in 1:length(tfwhiteproph_sigtfphossites)){
  ourphos <- tfwhiteproph_sigtfphossites[i]
  ourpro <- sub('^([^_]+_[^_]+).*', '\\1',ourphos)
  ourgene <- tfproanno[tfproanno$WhiteAd.Pro.ID %in% ourpro,"Ensembl"][1]
  tfwhiteproph_sigtfphossites_df[i,"Pro_ID"] <- ourpro
  tfwhiteproph_sigtfphossites_df[i,"RNA_ID"] <- ourgene
}

whitepidmeta <- data.frame(row.names = colnames(tfwhiteproph_prophnorm),
                            "pid" = colnames(tfwhiteproph_prophnorm),
                            "Sex" = rep("Female",length(colnames(tfwhiteproph_prophnorm))),
                            "Group" = rep("control",length(colnames(tfwhiteproph_prophnorm))))
for(i in 1:dim(whitepidmeta)[1]){
  ourpid <- rownames(whitepidmeta)[i]
  whitepidmeta[i,"Group"] <- pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourpid,"pass1bf0011"][1]
  whitepidmeta[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourpid,"pass1bf0027"][1]]
}
whitepidmeta[whitepidmeta$Group %in% "One-week program","Group"] <- "1w"
whitepidmeta[whitepidmeta$Group %in% "Two-week program","Group"] <- "2w"
whitepidmeta[whitepidmeta$Group %in% "Four-week program","Group"] <- "4w"
whitepidmeta[whitepidmeta$Group %in% "Eight-week program Training Group","Group"] <- "8w"
whitepidmeta[whitepidmeta$Group %in% "Eight-week program Control Group","Group"] <- "control"

ourrna <- tfwhiteproph_sigtfphossites_df$RNA_ID[79]
ourpro <- tfwhiteproph_sigtfphossites_df$Pro_ID[79]
#ourproph <- tfwhiteproph_sigtfphossites_df$Phos_ID[79]
ourdf <- data.frame(row.names = colnames(tfwhiteproph_rnanorm),
                    "RNA" = t(tfwhiteproph_rnanorm[ourrna,]),
                    "Pro" = t(tfwhiteproph_pronorm[ourpro,]),
                    "Sex" = whitepidmeta[colnames(tfwhiteproph_rnanorm),"Sex"],
                    "Group" = whitepidmeta[colnames(tfwhiteproph_rnanorm),"Group"])
colnames(ourdf) <- c("RNA","Pro","Sex","Group")

#png(file = paste("WAT-SC_",enstosym[ourrna,"Symbol"],"_",toupper(enstosym[ourrna,"Symbol"]),".png",sep = ""),width = 6,height = 6,units = "in",res = 600)
#ggscatter(ourdf, x = "RNA", y = "Pro",color = "Group",shape = "Sex",size = 4) + scale_color_manual(values = ann_cols$Group) + ggtitle(paste("WAT-SC: ",enstosym[ourrna,"Symbol"],":",toupper(enstosym[ourrna,"Symbol"])," Correlation: ",substr(toString(cor(ourdf$RNA,ourdf$Pro)),1,5),sep = "")) + xlab("RNA") + ylab("Proteome") + theme(legend.box = "vertical",legend.margin = margin())
#dev.off()

ourrna <- tfwhiteproph_sigtfphossites_df$RNA_ID[2]
ourpro <- tfwhiteproph_sigtfphossites_df$Pro_ID[2]
#ourproph <- tfwhiteproph_sigtfphossites_df$Phos_ID[2]
ourdf <- data.frame(row.names = colnames(tfwhiteproph_rnanorm),
                    "RNA" = t(tfwhiteproph_rnanorm[ourrna,]),
                    "Pro" = t(tfwhiteproph_pronorm[ourpro,]),
                    "Sex" = whitepidmeta[colnames(tfwhiteproph_rnanorm),"Sex"],
                    "Group" = whitepidmeta[colnames(tfwhiteproph_rnanorm),"Group"])
colnames(ourdf) <- c("RNA","Pro","Sex","Group")

#png(file = paste("WAT-SC_",enstosym[ourrna,"Symbol"],"_",toupper(enstosym[ourrna,"Symbol"]),".png",sep = ""),width = 6,height = 6,units = "in",res = 600)
#ggscatter(ourdf, x = "RNA", y = "Pro",color = "Group",shape = "Sex",size = 4) + scale_color_manual(values = ann_cols$Group) + ggtitle(paste("WAT-SC: ",enstosym[ourrna,"Symbol"],":",toupper(enstosym[ourrna,"Symbol"])," Correlation: ",substr(toString(cor(ourdf$RNA,ourdf$Pro)),1,5),sep = "")) + xlab("RNA") + ylab("Proteome") + theme(legend.box = "vertical",legend.margin = margin())
#dev.off()


ourrna <- tfwhiteproph_sigtfphossites_df$RNA_ID[43]
ourpro <- tfwhiteproph_sigtfphossites_df$Pro_ID[43]
#ourproph <- tfwhiteproph_sigtfphossites_df$Phos_ID[43]
ourdf <- data.frame(row.names = colnames(tfwhiteproph_rnanorm),
                    "RNA" = t(tfwhiteproph_rnanorm[ourrna,]),
                    "Pro" = t(tfwhiteproph_pronorm[ourpro,]),
                    "Sex" = whitepidmeta[colnames(tfwhiteproph_rnanorm),"Sex"],
                    "Group" = whitepidmeta[colnames(tfwhiteproph_rnanorm),"Group"])
colnames(ourdf) <- c("RNA","Pro","Sex","Group")

#png(file = paste("WAT-SC_",enstosym[ourrna,"Symbol"],"_",toupper(enstosym[ourrna,"Symbol"]),".png",sep = ""),width = 6,height = 6,units = "in",res = 600)
#ggscatter(ourdf, x = "RNA", y = "Pro",color = "Group",shape = "Sex",size = 4) + scale_color_manual(values = ann_cols$Group) + ggtitle(paste("WAT-SC: ",enstosym[ourrna,"Symbol"],":",toupper(enstosym[ourrna,"Symbol"])," Correlation: ",substr(toString(cor(ourdf$RNA,ourdf$Pro)),1,5),sep = "")) + xlab("RNA") + ylab("Proteome") + theme(legend.box = "vertical",legend.margin = margin())
#dev.off()

ourrna <- tfwhiteproph_sigtfphossites_df$RNA_ID[8]
ourpro <- tfwhiteproph_sigtfphossites_df$Pro_ID[8]
#ourproph <- tfwhiteproph_sigtfphossites_df$Phos_ID[8]
ourdf <- data.frame(row.names = colnames(tfwhiteproph_rnanorm),
                    "RNA" = t(tfwhiteproph_rnanorm[ourrna,]),
                    "Pro" = t(tfwhiteproph_pronorm[ourpro,]),
                    "Sex" = whitepidmeta[colnames(tfwhiteproph_rnanorm),"Sex"],
                    "Group" = whitepidmeta[colnames(tfwhiteproph_rnanorm),"Group"])
colnames(ourdf) <- c("RNA","Pro","Sex","Group")

#png(file = paste("WAT-SC_",enstosym[ourrna,"Symbol"],"_",toupper(enstosym[ourrna,"Symbol"]),".png",sep = ""),width = 6,height = 6,units = "in",res = 600)
#ggscatter(ourdf, x = "RNA", y = "Pro",color = "Group",shape = "Sex",size = 4) + scale_color_manual(values = ann_cols$Group) + ggtitle(paste("WAT-SC: ",enstosym[ourrna,"Symbol"],":",toupper(enstosym[ourrna,"Symbol"])," Correlation: ",substr(toString(cor(ourdf$RNA,ourdf$Pro)),1,5),sep = "")) + xlab("RNA") + ylab("Proteome") + theme(legend.box = "vertical",legend.margin = margin())
#dev.off()


ourrna <- tfwhiteproph_sigtfphossites_df$RNA_ID[73]
ourpro <- tfwhiteproph_sigtfphossites_df$Pro_ID[73]
#ourproph <- tfwhiteproph_sigtfphossites_df$Phos_ID[73]
ourdf <- data.frame(row.names = colnames(tfwhiteproph_rnanorm),
                    "RNA" = t(tfwhiteproph_rnanorm[ourrna,]),
                    "Pro" = t(tfwhiteproph_pronorm[ourpro,]),
                    "Sex" = whitepidmeta[colnames(tfwhiteproph_rnanorm),"Sex"],
                    "Group" = whitepidmeta[colnames(tfwhiteproph_rnanorm),"Group"])
colnames(ourdf) <- c("RNA","Pro","Sex","Group")

#png(file = paste("WAT-SC_",enstosym[ourrna,"Symbol"],"_",toupper(enstosym[ourrna,"Symbol"]),".png",sep = ""),width = 6,height = 6,units = "in",res = 600)
#ggscatter(ourdf, x = "RNA", y = "Pro",color = "Group",shape = "Sex",size = 4) + scale_color_manual(values = ann_cols$Group) + ggtitle(paste("WAT-SC: ",enstosym[ourrna,"Symbol"],":",toupper(enstosym[ourrna,"Symbol"])," Correlation: ",substr(toString(cor(ourdf$RNA,ourdf$Pro)),1,5),sep = "")) + xlab("RNA") + ylab("Proteome") + theme(legend.box = "vertical",legend.margin = margin())
#dev.off()


####
# Can we connect TFs with significant responses to training with correlated target DEGs with significant responses to training?

for(i in 1:dim(gastrornanorm)[2]){
  ourid <- colnames(gastrornanorm)[i]
  colnames(gastrornanorm)[i] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$viallabel %in% ourid,"pid"][1]
}
for(i in 1:dim(gastroprophnorm)[2]){
  ourid <- colnames(gastroprophnorm)[i]
  colnames(gastroprophnorm)[i] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$viallabel %in% ourid,"pid"][1]
}

for(i in 1:dim(heartrnanorm)[2]){
  ourid <- colnames(heartrnanorm)[i]
  colnames(heartrnanorm)[i] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$viallabel %in% ourid,"pid"][1]
}
for(i in 1:dim(heartprophnorm)[2]){
  ourid <- colnames(heartprophnorm)[i]
  colnames(heartprophnorm)[i] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$viallabel %in% ourid,"pid"][1]
}

for(i in 1:dim(kidneyrnanorm)[2]){
  ourid <- colnames(kidneyrnanorm)[i]
  colnames(kidneyrnanorm)[i] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$viallabel %in% ourid,"pid"][1]
}
for(i in 1:dim(kidneyprophnorm)[2]){
  ourid <- colnames(kidneyprophnorm)[i]
  colnames(kidneyprophnorm)[i] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$viallabel %in% ourid,"pid"][1]
}

for(i in 1:dim(liverrnanorm)[2]){
  ourid <- colnames(liverrnanorm)[i]
  colnames(liverrnanorm)[i] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$viallabel %in% ourid,"pid"][1]
}
for(i in 1:dim(liverprophnorm)[2]){
  ourid <- colnames(liverprophnorm)[i]
  colnames(liverprophnorm)[i] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$viallabel %in% ourid,"pid"][1]
}

for(i in 1:dim(lungrnanorm)[2]){
  ourid <- colnames(lungrnanorm)[i]
  colnames(lungrnanorm)[i] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$viallabel %in% ourid,"pid"][1]
}
for(i in 1:dim(lungprophnorm)[2]){
  ourid <- colnames(lungprophnorm)[i]
  colnames(lungprophnorm)[i] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$viallabel %in% ourid,"pid"][1]
}

for(i in 1:dim(whiternanorm)[2]){
  ourid <- colnames(whiternanorm)[i]
  colnames(whiternanorm)[i] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$viallabel %in% ourid,"pid"][1]
}
for(i in 1:dim(whiteprophnorm)[2]){
  ourid <- colnames(whiteprophnorm)[i]
  colnames(whiteprophnorm)[i] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$viallabel %in% ourid,"pid"][1]
}

mergegastrocols <- Reduce(intersect,list(colnames(gastrornanorm),colnames(gastroproprnorm),colnames(gastroprophnorm)))
gastrornanorm <- gastrornanorm[,mergegastrocols]
gastroproprnorm <- gastroproprnorm[,mergegastrocols]
gastroprophnorm <- gastroprophnorm[,mergegastrocols]

mergeheartcols <- Reduce(intersect,list(colnames(heartrnanorm),colnames(heartproprnorm),colnames(heartprophnorm)))
heartrnanorm <- heartrnanorm[,mergeheartcols]
heartproprnorm <- heartproprnorm[,mergeheartcols]
heartprophnorm <- heartprophnorm[,mergeheartcols]

mergekidneycols <- Reduce(intersect,list(colnames(kidneyrnanorm),colnames(kidneyproprnorm),colnames(kidneyprophnorm)))
kidneyrnanorm <- kidneyrnanorm[,mergekidneycols]
kidneyproprnorm <- kidneyproprnorm[,mergekidneycols]
kidneyprophnorm <- kidneyprophnorm[,mergekidneycols]

mergelivercols <- Reduce(intersect,list(colnames(liverrnanorm),colnames(liverproprnorm),colnames(liverprophnorm)))
liverrnanorm <- liverrnanorm[,mergelivercols]
liverproprnorm <- liverproprnorm[,mergelivercols]
liverprophnorm <- liverprophnorm[,mergelivercols]

mergelungcols <- Reduce(intersect,list(colnames(lungrnanorm),colnames(lungproprnorm),colnames(lungprophnorm)))
lungrnanorm <- lungrnanorm[,mergelungcols]
lungproprnorm <- lungproprnorm[,mergelungcols]
lungprophnorm <- lungprophnorm[,mergelungcols]

mergewhitecols <- Reduce(intersect,list(colnames(whiternanorm),colnames(whiteproprnorm),colnames(whiteprophnorm)))
whiternanorm <- whiternanorm[,mergewhitecols]
whiteproprnorm <- whiteproprnorm[,mergewhitecols]
whiteprophnorm <- whiteprophnorm[,mergewhitecols]

mergedtftargetsigdf$FullName <- ""

for(i in 1:dim(mergedtftargetsigdf)[1]){
  ourtf <- mergedtftargetsigdf$Transcription.Factor.Name[i]
  mergedtftargetsigdf[i,"FullName"] <- rownames(tfanno)[grep(ourtf,rownames(tfanno))][1]
}

# Define a structure to contain all the correlation mappings
tftargetcorrelationstructure <- list()

for(i in 1:(dim(mergedtftargetsigdf)[1]/2)){
  
  print(i)
  ourtf <- mergedtftargetsigdf[i,"FullName"]
  ourtissue <- mergedtftargetsigdf[i,"Tissue"]
  ourome <- mergedtftargetsigdf[i,"Ome"]
  
  
  if(ourtissue %in% "SKM-GN"){
    
    ourtftargetpeaks <- gastro50peakmotifs[gastro50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
    ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
    ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
    
    if(ourome %in% "RNA"){
      ourtftargetpromproxcormeta <- data.frame(row.names = intersect(ourtftargetpromproxgenes,rownames(gastrol2fcmat)),
                                               "Correlation" = rep(0,length(intersect(ourtftargetpromproxgenes,rownames(gastrol2fcmat)))))
      for(j in 1:dim(ourtftargetpromproxcormeta)[1]){
        ourtftargetpromproxcormeta[j,"Correlation"] <- cor(t(gastrornanorm[tfanno[ourtf,"Ensembl"],]),t(gastrornanorm[intersect(ourtftargetpromproxgenes,rownames(gastrol2fcmat))[j],]))
      }
    }
    if(ourome %in% "Prot"){
      ourtftargetpromproxcormeta <- data.frame(row.names = intersect(ourtftargetpromproxgenes,rownames(gastrol2fcmat)),
                                               "Correlation" = rep(0,length(intersect(ourtftargetpromproxgenes,rownames(gastrol2fcmat)))))
      for(j in 1:dim(ourtftargetpromproxcormeta)[1]){
        ourtftargetpromproxcormeta[j,"Correlation"] <- cor(t(gastroproprnorm[tfproanno[ourtf,"Gastro.Pro.ID"],]),t(gastrornanorm[intersect(ourtftargetpromproxgenes,rownames(gastrol2fcmat))[j],]))
      }
    }
    if(ourome %in% "Phos"){
      ourtftargetpromproxcormeta <- data.frame(row.names = intersect(ourtftargetpromproxgenes,rownames(gastrol2fcmat)),
                                               "Correlation" = rep(0,length(intersect(ourtftargetpromproxgenes,rownames(gastrol2fcmat)))))
      for(j in 1:dim(ourtftargetpromproxcormeta)[1]){
        ourtftargetpromproxcormeta[j,"Correlation"] <- cor(t(gastroprophnorm[tfgastroprophsig[grep(tfproanno[ourtf,"Gastro.Pro.ID"],tfgastroprophsig)],]),t(gastrornanorm[intersect(ourtftargetpromproxgenes,rownames(gastrol2fcmat))[j],]))
      }
    }
    
  }
 
  if(ourtissue %in% "HEART"){
    
    ourtftargetpeaks <- heart50peakmotifs[heart50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
    ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
    ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
    
    if(ourome %in% "RNA"){
      ourtftargetpromproxcormeta <- data.frame(row.names = intersect(ourtftargetpromproxgenes,rownames(heartl2fcmat)),
                                               "Correlation" = rep(0,length(intersect(ourtftargetpromproxgenes,rownames(heartl2fcmat)))))
      for(j in 1:dim(ourtftargetpromproxcormeta)[1]){
        ourtftargetpromproxcormeta[j,"Correlation"] <- cor(t(heartrnanorm[tfanno[ourtf,"Ensembl"],]),t(heartrnanorm[intersect(ourtftargetpromproxgenes,rownames(heartl2fcmat))[j],]))
      }
    }
    if(ourome %in% "Prot"){
      ourtftargetpromproxcormeta <- data.frame(row.names = intersect(ourtftargetpromproxgenes,rownames(heartl2fcmat)),
                                               "Correlation" = rep(0,length(intersect(ourtftargetpromproxgenes,rownames(heartl2fcmat)))))
      for(j in 1:dim(ourtftargetpromproxcormeta)[1]){
        ourtftargetpromproxcormeta[j,"Correlation"] <- cor(t(heartproprnorm[tfproanno[ourtf,"Heart.Pro.ID"],]),t(heartrnanorm[intersect(ourtftargetpromproxgenes,rownames(heartl2fcmat))[j],]))
      }
    }
    if(ourome %in% "Phos"){
      ourtftargetpromproxcormeta <- data.frame(row.names = intersect(ourtftargetpromproxgenes,rownames(heartl2fcmat)),
                                               "Correlation" = rep(0,length(intersect(ourtftargetpromproxgenes,rownames(heartl2fcmat)))))
      for(j in 1:dim(ourtftargetpromproxcormeta)[1]){
        ourtftargetpromproxcormeta[j,"Correlation"] <- cor(t(heartprophnorm[tfheartprophsig[grep(tfproanno[ourtf,"Heart.Pro.ID"],tfheartprophsig)],]),t(heartrnanorm[intersect(ourtftargetpromproxgenes,rownames(heartl2fcmat))[j],]))
      }
    }
    
  }
  
  if(ourtissue %in% "KIDNEY"){
    
    ourtftargetpeaks <- kidney50peakmotifs[kidney50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
    ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
    ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
    
    if(ourome %in% "RNA"){
      ourtftargetpromproxcormeta <- data.frame(row.names = intersect(ourtftargetpromproxgenes,rownames(kidneyl2fcmat)),
                                               "Correlation" = rep(0,length(intersect(ourtftargetpromproxgenes,rownames(kidneyl2fcmat)))))
      for(j in 1:dim(ourtftargetpromproxcormeta)[1]){
        ourtftargetpromproxcormeta[j,"Correlation"] <- cor(t(kidneyrnanorm[tfanno[ourtf,"Ensembl"],]),t(kidneyrnanorm[intersect(ourtftargetpromproxgenes,rownames(kidneyl2fcmat))[j],]))
      }
    }
    if(ourome %in% "Prot"){
      ourtftargetpromproxcormeta <- data.frame(row.names = intersect(ourtftargetpromproxgenes,rownames(kidneyl2fcmat)),
                                               "Correlation" = rep(0,length(intersect(ourtftargetpromproxgenes,rownames(kidneyl2fcmat)))))
      for(j in 1:dim(ourtftargetpromproxcormeta)[1]){
        ourtftargetpromproxcormeta[j,"Correlation"] <- cor(t(kidneyproprnorm[tfproanno[ourtf,"Kidney.Pro.ID"],]),t(kidneyrnanorm[intersect(ourtftargetpromproxgenes,rownames(kidneyl2fcmat))[j],]))
      }
    }
    if(ourome %in% "Phos"){
      ourtftargetpromproxcormeta <- data.frame(row.names = intersect(ourtftargetpromproxgenes,rownames(kidneyl2fcmat)),
                                               "Correlation" = rep(0,length(intersect(ourtftargetpromproxgenes,rownames(kidneyl2fcmat)))))
      for(j in 1:dim(ourtftargetpromproxcormeta)[1]){
        ourtftargetpromproxcormeta[j,"Correlation"] <- cor(t(kidneyprophnorm[tfkidneyprophsig[grep(tfproanno[ourtf,"Kidney.Pro.ID"],tfkidneyprophsig)],]),t(kidneyrnanorm[intersect(ourtftargetpromproxgenes,rownames(kidneyl2fcmat))[j],]))
      }
    }
    
  }
  
  if(ourtissue %in% "LIVER"){
    
    ourtftargetpeaks <- liver50peakmotifs[liver50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
    ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
    ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
    
    if(ourome %in% "RNA"){
      ourtftargetpromproxcormeta <- data.frame(row.names = intersect(ourtftargetpromproxgenes,rownames(liverl2fcmat)),
                                               "Correlation" = rep(0,length(intersect(ourtftargetpromproxgenes,rownames(liverl2fcmat)))))
      for(j in 1:dim(ourtftargetpromproxcormeta)[1]){
        ourtftargetpromproxcormeta[j,"Correlation"] <- cor(t(liverrnanorm[tfanno[ourtf,"Ensembl"],]),t(liverrnanorm[intersect(ourtftargetpromproxgenes,rownames(liverl2fcmat))[j],]))
      }
    }
    if(ourome %in% "Prot"){
      ourtftargetpromproxcormeta <- data.frame(row.names = intersect(ourtftargetpromproxgenes,rownames(liverl2fcmat)),
                                               "Correlation" = rep(0,length(intersect(ourtftargetpromproxgenes,rownames(liverl2fcmat)))))
      for(j in 1:dim(ourtftargetpromproxcormeta)[1]){
        ourtftargetpromproxcormeta[j,"Correlation"] <- cor(t(liverproprnorm[tfproanno[ourtf,"Liver.Pro.ID"],]),t(liverrnanorm[intersect(ourtftargetpromproxgenes,rownames(liverl2fcmat))[j],]))
      }
    }
    if(ourome %in% "Phos"){
      ourtftargetpromproxcormeta <- data.frame(row.names = intersect(ourtftargetpromproxgenes,rownames(liverl2fcmat)),
                                               "Correlation" = rep(0,length(intersect(ourtftargetpromproxgenes,rownames(liverl2fcmat)))))
      for(j in 1:dim(ourtftargetpromproxcormeta)[1]){
        ourtftargetpromproxcormeta[j,"Correlation"] <- cor(t(liverprophnorm[tfliverprophsig[grep(tfproanno[ourtf,"Liver.Pro.ID"],tfliverprophsig)][1],]),t(liverrnanorm[intersect(ourtftargetpromproxgenes,rownames(liverl2fcmat))[j],]))
      }
    }
    
  }
  
  if(ourtissue %in% "LUNG"){
    
    ourtftargetpeaks <- lung50peakmotifs[lung50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
    ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
    ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
    
    if(ourome %in% "RNA"){
      ourtftargetpromproxcormeta <- data.frame(row.names = intersect(ourtftargetpromproxgenes,rownames(lungl2fcmat)),
                                               "Correlation" = rep(0,length(intersect(ourtftargetpromproxgenes,rownames(lungl2fcmat)))))
      for(j in 1:dim(ourtftargetpromproxcormeta)[1]){
        ourtftargetpromproxcormeta[j,"Correlation"] <- cor(t(lungrnanorm[tfanno[ourtf,"Ensembl"],]),t(lungrnanorm[intersect(ourtftargetpromproxgenes,rownames(lungl2fcmat))[j],]))
      }
    }
    if(ourome %in% "Prot"){
      ourtftargetpromproxcormeta <- data.frame(row.names = intersect(ourtftargetpromproxgenes,rownames(lungl2fcmat)),
                                               "Correlation" = rep(0,length(intersect(ourtftargetpromproxgenes,rownames(lungl2fcmat)))))
      for(j in 1:dim(ourtftargetpromproxcormeta)[1]){
        ourtftargetpromproxcormeta[j,"Correlation"] <- cor(t(lungproprnorm[tfproanno[ourtf,"Lung.Pro.ID"],]),t(lungrnanorm[intersect(ourtftargetpromproxgenes,rownames(lungl2fcmat))[j],]))
      }
    }
    if(ourome %in% "Phos"){
      ourtftargetpromproxcormeta <- data.frame(row.names = intersect(ourtftargetpromproxgenes,rownames(lungl2fcmat)),
                                               "Correlation" = rep(0,length(intersect(ourtftargetpromproxgenes,rownames(lungl2fcmat)))))
      for(j in 1:dim(ourtftargetpromproxcormeta)[1]){
        ourtftargetpromproxcormeta[j,"Correlation"] <- cor(t(lungprophnorm[tflungprophsig[grep(tfproanno[ourtf,"Lung.Pro.ID"],tflungprophsig)][1],]),t(lungrnanorm[intersect(ourtftargetpromproxgenes,rownames(lungl2fcmat))[j],]))
      }
    }
    
  }
  
  if(ourtissue %in% "WAT-SC"){
    
    ourtftargetpeaks <- white50peakmotifs[white50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
    ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
    ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
    
    if(ourome %in% "RNA"){
      ourtftargetpromproxcormeta <- data.frame(row.names = intersect(ourtftargetpromproxgenes,rownames(whitel2fcmat)),
                                               "Correlation" = rep(0,length(intersect(ourtftargetpromproxgenes,rownames(whitel2fcmat)))))
      for(j in 1:dim(ourtftargetpromproxcormeta)[1]){
        ourtftargetpromproxcormeta[j,"Correlation"] <- cor(t(whiternanorm[tfanno[ourtf,"Ensembl"],]),t(whiternanorm[intersect(ourtftargetpromproxgenes,rownames(whitel2fcmat))[j],]))
      }
    }
    if(ourome %in% "Prot"){
      ourtftargetpromproxcormeta <- data.frame(row.names = intersect(ourtftargetpromproxgenes,rownames(whitel2fcmat)),
                                               "Correlation" = rep(0,length(intersect(ourtftargetpromproxgenes,rownames(whitel2fcmat)))))
      for(j in 1:dim(ourtftargetpromproxcormeta)[1]){
        ourtftargetpromproxcormeta[j,"Correlation"] <- cor(t(whiteproprnorm[tfproanno[ourtf,"WhiteAd.Pro.ID"],]),t(whiternanorm[intersect(ourtftargetpromproxgenes,rownames(whitel2fcmat))[j],]))
      }
    }
    if(ourome %in% "Phos"){
      ourtftargetpromproxcormeta <- data.frame(row.names = intersect(ourtftargetpromproxgenes,rownames(whitel2fcmat)),
                                               "Correlation" = rep(0,length(intersect(ourtftargetpromproxgenes,rownames(whitel2fcmat)))))
      for(j in 1:dim(ourtftargetpromproxcormeta)[1]){
        ourtftargetpromproxcormeta[j,"Correlation"] <- cor(t(whiteprophnorm[tfwhiteprophsig[grep(tfproanno[ourtf,"WhiteAd.Pro.ID"],tfwhiteprophsig)][1],]),t(whiternanorm[intersect(ourtftargetpromproxgenes,rownames(whitel2fcmat))[j],]))
      }
    }
    
  }
  
  tftargetcorrelationstructure[[i]] <- ourtftargetpromproxcormeta
  
}

for(i in 1:(dim(mergedtftargetsigdf)[1]/2)){
  
  tftargetcorrelationstructure[[i]]$TF <- mergedtftargetsigdf[i,"Transcription.Factor.Name"]
  tftargetcorrelationstructure[[i]]$Tissue <- mergedtftargetsigdf[i,"Tissue"]
  tftargetcorrelationstructure[[i]]$Ome <- mergedtftargetsigdf[i,"Ome"]
  tftargetcorrelationstructure[[i]]$TF_Tissue_Ome <- paste(tftargetcorrelationstructure[[i]]$TF,"_", mergedtftargetsigdf[i,"Tissue"],"_",tftargetcorrelationstructure[[i]]$Ome,sep = "")
  tftargetcorrelationstructure[[i]]$GeneEns <- rownames(tftargetcorrelationstructure[[i]])
}

tftargetcorrelationstructurecombo <- rbind(tftargetcorrelationstructure[[1]],tftargetcorrelationstructure[[2]])
for(i in 3:length(tftargetcorrelationstructure)){
  tftargetcorrelationstructurecombo <- rbind(tftargetcorrelationstructurecombo,tftargetcorrelationstructure[[i]])
}

tftargetcorrelationstructurecombo$Training.Response <- "Non.DEG"
for(i in 1:dim(tftargetcorrelationstructurecombo)[1]){
  ourgene <- tftargetcorrelationstructurecombo[i,"GeneEns"]
  if(tftargetcorrelationstructurecombo[i,"Tissue"] %in% "SKM-GN"){
    if(ourgene %in% gastrornasig){
      tftargetcorrelationstructurecombo[i,"Training.Response"] <- "DEG"
    }
  }
  if(tftargetcorrelationstructurecombo[i,"Tissue"] %in% "HEART"){
    if(ourgene %in% heartrnasig){
      tftargetcorrelationstructurecombo[i,"Training.Response"] <- "DEG"
    }
  }
  if(tftargetcorrelationstructurecombo[i,"Tissue"] %in% "KIDNEY"){
    if(ourgene %in% kidneyrnasig){
      tftargetcorrelationstructurecombo[i,"Training.Response"] <- "DEG"
    }
  }
  if(tftargetcorrelationstructurecombo[i,"Tissue"] %in% "LIVER"){
    if(ourgene %in% liverrnasig){
      tftargetcorrelationstructurecombo[i,"Training.Response"] <- "DEG"
    }
  }
  if(tftargetcorrelationstructurecombo[i,"Tissue"] %in% "LUNG"){
    if(ourgene %in% lungrnasig){
      tftargetcorrelationstructurecombo[i,"Training.Response"] <- "DEG"
    }
  }
  if(tftargetcorrelationstructurecombo[i,"Tissue"] %in% "WAT-SC"){
    if(ourgene %in% whiternasig){
      tftargetcorrelationstructurecombo[i,"Training.Response"] <- "DEG"
    }
  }
}
tftargetcorrelationstructurecombo$Training.Response <- factor(tftargetcorrelationstructurecombo$Training.Response,levels = c("Non.DEG","DEG"))
tftargetcorrelationstructurecombo$TF_Ome <- paste(tftargetcorrelationstructurecombo$TF,tftargetcorrelationstructurecombo$Ome,sep = "_")

#png(file = "Supplemental Figure S26A_021315.png",width = 10,height = 6,units = "in",res = 600)
#ggdotplot(tftargetcorrelationstructurecombo[tftargetcorrelationstructurecombo$Tissue %in% "SKM-GN",],x = "TF_Ome",y = "Correlation",size = 0.5, color = "Training.Response",fill = "Training.Response",binwidth = 0.03) + geom_hline(yintercept = 0) + ggtitle("SKM-GN EET Significant TFs")
#dev.off()

#png(file = "Supplemental Figure S26B_021315.png",width = 10,height = 6,units = "in",res = 600)
#ggdotplot(tftargetcorrelationstructurecombo[tftargetcorrelationstructurecombo$Tissue %in% "HEART",],x = "TF_Ome",y = "Correlation",size = 0.5, color = "Training.Response",fill = "Training.Response",binwidth = 0.03) + geom_hline(yintercept = 0) + ggtitle("HEART EET Significant TFs")
#dev.off()

#png(file = "Supplemental Figure S26C_021315.png",width = 10,height = 6,units = "in",res = 600)
#ggdotplot(tftargetcorrelationstructurecombo[tftargetcorrelationstructurecombo$Tissue %in% "KIDNEY",],x = "TF_Ome",y = "Correlation",size = 0.5, color = "Training.Response",fill = "Training.Response",binwidth = 0.0175) + geom_hline(yintercept = 0) + ggtitle("KIDNEY EET Significant TFs")
#dev.off()

#png(file = "Supplemental Figure S26D_021315.png",width = 10,height = 6,units = "in",res = 600)
#ggdotplot(tftargetcorrelationstructurecombo[tftargetcorrelationstructurecombo$Tissue %in% "LIVER",],x = "TF_Ome",y = "Correlation",size = 0.5, color = "Training.Response",fill = "Training.Response",binwidth = 0.03) + geom_hline(yintercept = 0) + ggtitle("LIVER EET Significant TFs")
#dev.off()

#png(file = "Supplemental Figure S26E_021315.png",width = 10,height = 6,units = "in",res = 600)
#ggdotplot(tftargetcorrelationstructurecombo[tftargetcorrelationstructurecombo$Tissue %in% "LUNG",],x = "TF_Ome",y = "Correlation",size = 0.5, color = "Training.Response",fill = "Training.Response",binwidth = 0.03) + geom_hline(yintercept = 0) + ggtitle("LUNG EET Significant TFs")
#dev.off()

#png(file = "Supplemental Figure S26F_021315.png",width = 10,height = 6,units = "in",res = 600)
#ggdotplot(tftargetcorrelationstructurecombo[tftargetcorrelationstructurecombo$Tissue %in% "WAT-SC",],x = "TF_Ome",y = "Correlation",size = 1, color = "Training.Response",fill = "Training.Response",binwidth = 0.03) + geom_hline(yintercept = 0) + ggtitle("WAT-SC EET Significant TFs")
#dev.off()

#pdf(file = "Supplemental Figure S26A_021315.pdf",width = 10,height = 6)
#ggdotplot(tftargetcorrelationstructurecombo[tftargetcorrelationstructurecombo$Tissue %in% "SKM-GN",],x = "TF_Ome",y = "Correlation",size = 0.5, color = "Training.Response",fill = "Training.Response",binwidth = 0.03) + geom_hline(yintercept = 0) + ggtitle("SKM-GN EET Significant TFs")
#dev.off()

#pdf(file = "Supplemental Figure S26B_021315.pdf",width = 10,height = 6)
#ggdotplot(tftargetcorrelationstructurecombo[tftargetcorrelationstructurecombo$Tissue %in% "HEART",],x = "TF_Ome",y = "Correlation",size = 0.5, color = "Training.Response",fill = "Training.Response",binwidth = 0.03) + geom_hline(yintercept = 0) + ggtitle("HEART EET Significant TFs")
#dev.off()

#pdf(file = "Supplemental Figure S26C_021315.pdf",width = 10,height = 6)
#ggdotplot(tftargetcorrelationstructurecombo[tftargetcorrelationstructurecombo$Tissue %in% "KIDNEY",],x = "TF_Ome",y = "Correlation",size = 0.5, color = "Training.Response",fill = "Training.Response",binwidth = 0.0175) + geom_hline(yintercept = 0) + ggtitle("KIDNEY EET Significant TFs")
#dev.off()

#pdf(file = "Supplemental Figure S26D_021315.pdf",width = 10,height = 6)
#ggdotplot(tftargetcorrelationstructurecombo[tftargetcorrelationstructurecombo$Tissue %in% "LIVER",],x = "TF_Ome",y = "Correlation",size = 0.5, color = "Training.Response",fill = "Training.Response",binwidth = 0.03) + geom_hline(yintercept = 0) + ggtitle("LIVER EET Significant TFs")
#dev.off()

#pdf(file = "Supplemental Figure S26E_021315.pdf",width = 10,height = 6)
#ggdotplot(tftargetcorrelationstructurecombo[tftargetcorrelationstructurecombo$Tissue %in% "LUNG",],x = "TF_Ome",y = "Correlation",size = 0.5, color = "Training.Response",fill = "Training.Response",binwidth = 0.03) + geom_hline(yintercept = 0) + ggtitle("LUNG EET Significant TFs")
#dev.off()

#pdf(file = "Supplemental Figure S26F_021315.pdf",width = 10,height = 6)
#ggdotplot(tftargetcorrelationstructurecombo[tftargetcorrelationstructurecombo$Tissue %in% "WAT-SC",],x = "TF_Ome",y = "Correlation",size = 1, color = "Training.Response",fill = "Training.Response",binwidth = 0.03) + geom_hline(yintercept = 0) + ggtitle("WAT-SC EET Significant TFs")
#dev.off()


save.image("Figure5_021325.RData")

ourtf <- "PU.1(ETS)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer"
ourtftargetpeaks <- lung50peakmotifs[lung50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])

sextimemeta <- data.frame(row.names = colnames(gastrol2fcmat),
                          "Sex" = c(rep("Female",4),rep("Male",4)),
                          "Week" = rep(c("1w","2w","4w","8w"),2),
                          "TF.L2FC" = lungl2fcmat[tfanno[ourtf,"Ensembl"],])
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

png(file = "Figure_5J_021325.png",width = 7,height = 5,units = "in",res = 600)
pheatmap(lungl2fcmat[intersect(ourtftargetpromproxgenes,lungrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,lungrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
dev.off()

pdf(file = "Figure_5J_021325.pdf",width = 7,height = 5)
pheatmap(lungl2fcmat[intersect(ourtftargetpromproxgenes,lungrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,lungrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
dev.off()



ourtf <- "IRF8(IRF)/BMDM-IRF8-ChIP-Seq(GSE77884)/Homer"
ourtftargetpeaks <- lung50peakmotifs[lung50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])

sextimemeta <- data.frame(row.names = colnames(gastrol2fcmat),
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
  ourtftargetpromproxcormeta[i,"Correlation"] <- cor(lungl2fcmat[tfanno[ourtf,"Ensembl"],],lungl2fcmat[intersect(ourtftargetpromproxgenes,lungrnasig)[i],])
}
ourtftargetpromproxcormeta$Correlation <- round(ourtftargetpromproxcormeta$Correlation/0.2)*0.2
ourtftargetpromproxcormeta$Correlation <- factor(ourtftargetpromproxcormeta$Correlation,levels = c(1,0.8,0.6,0.4,0.2,0,-0.2,-0.4,-0.6,-0.8,-1))

png(file = "Figure 5K_021325.png",width = 7,height = 5,units = "in",res = 600)
pheatmap(lungl2fcmat[intersect(ourtftargetpromproxgenes,lungrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,lungrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
dev.off()

pdf(file = "Figure 5K_021325.pdf",width = 7,height = 5)
pheatmap(lungl2fcmat[intersect(ourtftargetpromproxgenes,lungrnasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,lungrnasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("TF.L2FC","Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F,border_color = NA)
dev.off()


#####
# Now let us look for shared patterns of TF response across tissues. This will accompany Figure 5
####

tfgastrornasigtfs <- rownames(tfanno)[tfanno$Ensembl %in% tfgastrornasig]
tfgastroprosigtfs <- rownames(tfproanno)[tfproanno$Gastro.Pro.ID %in% tfgastroprosig]
tfgastroprophsigtfs <- rownames(tfproanno)[tfproanno$Gastro.Pro.ID %in% sub('^([^_]+_[^_]+).*', '\\1',tfgastroprophsig)]
tfgastrosigtfs <- Reduce(union,list(tfgastrornasigtfs,tfgastroprosigtfs,tfgastroprophsigtfs))

tfheartrnasigtfs <- rownames(tfanno)[tfanno$Ensembl %in% tfheartrnasig]
tfheartprosigtfs <- rownames(tfproanno)[tfproanno$Heart.Pro.ID %in% tfheartprosig]
tfheartprophsigtfs <- rownames(tfproanno)[tfproanno$Heart.Pro.ID %in% sub('^([^_]+_[^_]+).*', '\\1',tfheartprophsig)]
tfheartsigtfs <- Reduce(union,list(tfheartrnasigtfs,tfheartprosigtfs,tfheartprophsigtfs))

tfkidneyrnasigtfs <- rownames(tfanno)[tfanno$Ensembl %in% tfkidneyrnasig]
tfkidneyprosigtfs <- rownames(tfproanno)[tfproanno$Kidney.Pro.ID %in% tfkidneyprosig]
tfkidneyprophsigtfs <- rownames(tfproanno)[tfproanno$Kidney.Pro.ID %in% sub('^([^_]+_[^_]+).*', '\\1',tfkidneyprophsig)]
tfkidneysigtfs <- Reduce(union,list(tfkidneyrnasigtfs,tfkidneyprosigtfs,tfkidneyprophsigtfs))

tfliverrnasigtfs <- rownames(tfanno)[tfanno$Ensembl %in% tfliverrnasig]
tfliverprosigtfs <- rownames(tfproanno)[tfproanno$Liver.Pro.ID %in% tfliverprosig]
tfliverprophsigtfs <- rownames(tfproanno)[tfproanno$Liver.Pro.ID %in% sub('^([^_]+_[^_]+).*', '\\1',tfliverprophsig)]
tfliversigtfs <- Reduce(union,list(tfliverrnasigtfs,tfliverprosigtfs,tfliverprophsigtfs))

tflungrnasigtfs <- rownames(tfanno)[tfanno$Ensembl %in% tflungrnasig]
tflungprosigtfs <- rownames(tfproanno)[tfproanno$Lung.Pro.ID %in% tflungprosig]
tflungprophsigtfs <- rownames(tfproanno)[tfproanno$Lung.Pro.ID %in% sub('^([^_]+_[^_]+).*', '\\1',tflungprophsig)]
tflungsigtfs <- Reduce(union,list(tflungrnasigtfs,tflungprosigtfs,tflungprophsigtfs))

tfwhiternasigtfs <- rownames(tfanno)[tfanno$Ensembl %in% tfwhiternasig]
tfwhiteprosigtfs <- rownames(tfproanno)[tfproanno$WhiteAd.Pro.ID %in% tfwhiteprosig]
tfwhiteprophsigtfs <- rownames(tfproanno)[tfproanno$WhiteAd.Pro.ID %in% sub('^([^_]+_[^_]+).*', '\\1',tfwhiteprophsig)]
tfwhitesigtfs <- Reduce(union,list(tfwhiternasigtfs,tfwhiteprosigtfs,tfwhiteprophsigtfs))

tfhippornasigtfs <- rownames(tfanno)[tfanno$Ensembl %in% tfhippornasig]
tfbrownrnasigtfs <- rownames(tfanno)[tfanno$Ensembl %in% tfbrownrnasig]


tfallsigtfs <- Reduce(union,list(tfgastrosigtfs,tfheartsigtfs,tfkidneysigtfs,tfliversigtfs,tflungsigtfs,tfwhitesigtfs))
tfallrnasigtfs <- Reduce(union,list(tfgastrornasigtfs,tfheartrnasigtfs,tfhippornasigtfs,tfkidneyrnasigtfs,tfliverrnasigtfs,tflungrnasigtfs,tfbrownrnasigtfs,tfwhiternasigtfs))
tfallprosigtfs <- Reduce(union,list(tfgastroprosigtfs,tfheartprosigtfs,tfkidneyprosigtfs,tfliverprosigtfs,tflungprosigtfs,tfwhiteprosigtfs))
tfallprophsigtfs <- Reduce(union,list(tfgastroprophsigtfs,tfheartprophsigtfs,tfkidneyprophsigtfs,tfliverprophsigtfs,tflungprophsigtfs,tfwhiteprophsigtfs))


tfallsigtfsdf <- data.frame(row.names = tfallsigtfs,
                            "SKM-GN" = rep(0,length(tfallsigtfs)),
                            "HEART" = rep(0,length(tfallsigtfs)),
                            "KIDNEY" = rep(0,length(tfallsigtfs)),
                            "LIVER" = rep(0,length(tfallsigtfs)),
                            "LUNG" = rep(0,length(tfallsigtfs)),
                            "WAT-SC" = rep(0,length(tfallsigtfs)))

tfallrnasigtfsdf <- data.frame(row.names = tfallrnasigtfs,
                            "SKM-GN" = rep(0,length(tfallrnasigtfs)),
                            "HEART" = rep(0,length(tfallrnasigtfs)),
                            "HIPPOC" = rep(0,length(tfallrnasigtfs)),
                            "KIDNEY" = rep(0,length(tfallrnasigtfs)),
                            "LIVER" = rep(0,length(tfallrnasigtfs)),
                            "LUNG" = rep(0,length(tfallrnasigtfs)),
                            "BAT" = rep(0,length(tfallrnasigtfs)),
                            "WAT-SC" = rep(0,length(tfallrnasigtfs)))

tfallprosigtfsdf <- data.frame(row.names = tfallprosigtfs,
                            "SKM-GN" = rep(0,length(tfallprosigtfs)),
                            "HEART" = rep(0,length(tfallprosigtfs)),
                            "KIDNEY" = rep(0,length(tfallprosigtfs)),
                            "LIVER" = rep(0,length(tfallprosigtfs)),
                            "LUNG" = rep(0,length(tfallprosigtfs)),
                            "WAT-SC" = rep(0,length(tfallprosigtfs)))

tfallprophsigtfsdf <- data.frame(row.names = tfallprophsigtfs,
                            "SKM-GN" = rep(0,length(tfallprophsigtfs)),
                            "HEART" = rep(0,length(tfallprophsigtfs)),
                            "KIDNEY" = rep(0,length(tfallprophsigtfs)),
                            "LIVER" = rep(0,length(tfallprophsigtfs)),
                            "LUNG" = rep(0,length(tfallprophsigtfs)),
                            "WAT-SC" = rep(0,length(tfallprophsigtfs)))

for(i in 1:length(tfallsigtfs)){
  ourtf <- tfallsigtfs[i]
  if(ourtf %in% tfgastrosigtfs){
    tfallsigtfsdf[i,"SKM.GN"] <- 1
  }
  if(ourtf %in% tfheartsigtfs){
    tfallsigtfsdf[i,"HEART"] <- 1
  }
  if(ourtf %in% tfkidneysigtfs){
    tfallsigtfsdf[i,"KIDNEY"] <- 1
  }
  if(ourtf %in% tfliversigtfs){
    tfallsigtfsdf[i,"LIVER"] <- 1
  }
  if(ourtf %in% tflungsigtfs){
    tfallsigtfsdf[i,"LUNG"] <- 1
  }
  if(ourtf %in% tfwhitesigtfs){
    tfallsigtfsdf[i,"WAT.SC"] <- 1
  }
}

for(i in 1:length(tfallrnasigtfs)){
  ourtf <- tfallrnasigtfs[i]
  if(ourtf %in% tfgastrornasigtfs){
    tfallrnasigtfsdf[i,"SKM.GN"] <- 1
  }
  if(ourtf %in% tfheartrnasigtfs){
    tfallrnasigtfsdf[i,"HEART"] <- 1
  }
  if(ourtf %in% tfhippornasigtfs){
    tfallrnasigtfsdf[i,"HIPPOC"] <- 1
  }
  if(ourtf %in% tfkidneyrnasigtfs){
    tfallrnasigtfsdf[i,"KIDNEY"] <- 1
  }
  if(ourtf %in% tfliverrnasigtfs){
    tfallrnasigtfsdf[i,"LIVER"] <- 1
  }
  if(ourtf %in% tflungrnasigtfs){
    tfallrnasigtfsdf[i,"LUNG"] <- 1
  }
  if(ourtf %in% tfbrownrnasigtfs){
    tfallrnasigtfsdf[i,"BAT"] <- 1
  }
  if(ourtf %in% tfwhiternasigtfs){
    tfallrnasigtfsdf[i,"WAT.SC"] <- 1
  }
}

for(i in 1:length(tfallprosigtfs)){
  ourtf <- tfallprosigtfs[i]
  if(ourtf %in% tfgastroprosigtfs){
    tfallprosigtfsdf[i,"SKM.GN"] <- 1
  }
  if(ourtf %in% tfheartprosigtfs){
    tfallprosigtfsdf[i,"HEART"] <- 1
  }
  if(ourtf %in% tfkidneyprosigtfs){
    tfallprosigtfsdf[i,"KIDNEY"] <- 1
  }
  if(ourtf %in% tfliverprosigtfs){
    tfallprosigtfsdf[i,"LIVER"] <- 1
  }
  if(ourtf %in% tflungprosigtfs){
    tfallprosigtfsdf[i,"LUNG"] <- 1
  }
  if(ourtf %in% tfwhiteprosigtfs){
    tfallprosigtfsdf[i,"WAT.SC"] <- 1
  }
}

for(i in 1:length(tfallprophsigtfs)){
  ourtf <- tfallprophsigtfs[i]
  if(ourtf %in% tfgastroprophsigtfs){
    tfallprophsigtfsdf[i,"SKM.GN"] <- 1
  }
  if(ourtf %in% tfheartprophsigtfs){
    tfallprophsigtfsdf[i,"HEART"] <- 1
  }
  if(ourtf %in% tfkidneyprophsigtfs){
    tfallprophsigtfsdf[i,"KIDNEY"] <- 1
  }
  if(ourtf %in% tfliverprophsigtfs){
    tfallprophsigtfsdf[i,"LIVER"] <- 1
  }
  if(ourtf %in% tflungprophsigtfs){
    tfallprophsigtfsdf[i,"LUNG"] <- 1
  }
  if(ourtf %in% tfwhiteprophsigtfs){
    tfallprophsigtfsdf[i,"WAT.SC"] <- 1
  }
}

protissuemeta <- data.frame(row.names = colnames(tfallsigtfsdf),
                            "Tissue" = gsub("\\.","-",colnames(tfallsigtfsdf)))

rnatissuemeta <- data.frame(row.names = colnames(tfallrnasigtfsdf),
                            "Tissue" = gsub("\\.","-",colnames(tfallrnasigtfsdf)))

png(file = "Supplemental Figure S20A_021325.png",width = 6,height = 18,units = "in",res = 600)
pheatmap(tfallrnasigtfsdf,labels_row = gsub("\\/.*","",gsub("\\(.*","",tfallrnasigtfs)),angle_col = 0,breaks = seq(0,1,length.out = 3),color = colorpanel(3,"white","red"),labels_col = gsub("\\.","-",colnames(tfallrnasigtfsdf)),annotation_col = rnatissuemeta,annotation_colors = ann_cols,show_colnames = F,cluster_cols = F,cluster_rows = F,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S20A_021325.pdf",width = 6,height = 18)
pheatmap(tfallrnasigtfsdf,labels_row = gsub("\\/.*","",gsub("\\(.*","",tfallrnasigtfs)),angle_col = 0,breaks = seq(0,1,length.out = 3),color = colorpanel(3,"white","red"),labels_col = gsub("\\.","-",colnames(tfallrnasigtfsdf)),annotation_col = rnatissuemeta,annotation_colors = ann_cols,show_colnames = F,cluster_cols = F,cluster_rows = F,border_color = NA)
dev.off()

png(file = "Supplemental Figure S20B_021325.png",width = 6,height = 10,units = "in",res = 600)
pheatmap(tfallprosigtfsdf,labels_row = gsub("\\/.*","",gsub("\\(.*","",tfallprosigtfs)),angle_col = 0,breaks = seq(0,1,length.out = 3),color = colorpanel(3,"white","red"),labels_col = gsub("\\.","-",colnames(tfallprosigtfsdf)),annotation_col = protissuemeta,annotation_colors = ann_cols,show_colnames = F,cluster_cols = F,cluster_rows = F,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S20B_021325.pdf",width = 6,height = 10)
pheatmap(tfallprosigtfsdf,labels_row = gsub("\\/.*","",gsub("\\(.*","",tfallprosigtfs)),angle_col = 0,breaks = seq(0,1,length.out = 3),color = colorpanel(3,"white","red"),labels_col = gsub("\\.","-",colnames(tfallprosigtfsdf)),annotation_col = protissuemeta,annotation_colors = ann_cols,show_colnames = F,cluster_cols = F,cluster_rows = F,border_color = NA)
dev.off()

png(file = "Supplemental Figure S20C_021325.png",width = 6,height = 8,units = "in",res = 600)
pheatmap(tfallprophsigtfsdf,labels_row = gsub("\\/.*","",gsub("\\(.*","",tfallprophsigtfs)),angle_col = 0,breaks = seq(0,1,length.out = 3),color = colorpanel(3,"white","red"),labels_col = gsub("\\.","-",colnames(tfallprophsigtfsdf)),annotation_col = protissuemeta,annotation_colors = list("Tissue" = ann_cols$Tissue[gsub("\\.","-",colnames(tfallsigtfsdf))]),show_colnames = F,cluster_cols = F,cluster_rows = F,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S20C_021325.pdf",width = 6,height = 8)
pheatmap(tfallprophsigtfsdf,labels_row = gsub("\\/.*","",gsub("\\(.*","",tfallprophsigtfs)),angle_col = 0,breaks = seq(0,1,length.out = 3),color = colorpanel(3,"white","red"),labels_col = gsub("\\.","-",colnames(tfallprophsigtfsdf)),annotation_col = protissuemeta,annotation_colors = list("Tissue" = ann_cols$Tissue[gsub("\\.","-",colnames(tfallsigtfsdf))]),show_colnames = F,cluster_cols = F,cluster_rows = F,border_color = NA)
dev.off()

tfallrnasigtfsoverlap <- matrix(0L,nrow = 8,ncol = 8)
rownames(tfallrnasigtfsoverlap) <- rnatissuemeta$Tissue
colnames(tfallrnasigtfsoverlap) <- rnatissuemeta$Tissue

tfallprosigtfsoverlap <- matrix(0L,nrow = 6,ncol = 6)
rownames(tfallprosigtfsoverlap) <- protissuemeta$Tissue
colnames(tfallprosigtfsoverlap) <- protissuemeta$Tissue

tfallprophsigtfsoverlap <- matrix(0L,nrow = 6,ncol = 6)
rownames(tfallprophsigtfsoverlap) <- protissuemeta$Tissue
colnames(tfallprophsigtfsoverlap) <- protissuemeta$Tissue

tfrnasiglist <- list("SKM-GN" = tfgastrornasigtfs,
                     "HEART" = tfheartrnasigtfs,
                     "HIPPOC" = tfhippornasigtfs,
                     "KIDNEY" = tfkidneyrnasigtfs,
                     "LIVER" = tfliverrnasigtfs,
                     "LUNG" = tflungrnasigtfs,
                     "BAT" = tfbrownrnasigtfs,
                     "WAT-SC" = tfwhiternasigtfs)

tfprosiglist <- list("SKM-GN" = tfgastroprosigtfs,
                     "HEART" = tfheartprosigtfs,
                     "KIDNEY" = tfkidneyprosigtfs,
                     "LIVER" = tfliverprosigtfs,
                     "LUNG" = tflungprosigtfs,
                     "WAT-SC" = tfwhiteprosigtfs)

tfprophsiglist <- list("SKM-GN" = tfgastroprophsigtfs,
                     "HEART" = tfheartprophsigtfs,
                     "KIDNEY" = tfkidneyprophsigtfs,
                     "LIVER" = tfliverprophsigtfs,
                     "LUNG" = tflungprophsigtfs,
                     "WAT-SC" = tfwhiteprophsigtfs)

newrnatissuemeta <- data.frame(row.names = rownames(tfallrnasigtfsoverlap),
                               "Tissue" = rownames(tfallrnasigtfsoverlap))
newprotissuemeta <- data.frame(row.names = rownames(tfallprosigtfsoverlap),
                               "Tissue" = rownames(tfallprosigtfsoverlap))

for(i in 1:dim(tfallrnasigtfsoverlap)[1]){
  for(j in 1:dim(tfallrnasigtfsoverlap)[2]){
    tfallrnasigtfsoverlap[i,j] <- length(intersect(tfrnasiglist[[i]],tfrnasiglist[[j]]))
  }
}

for(i in 1:dim(tfallprosigtfsoverlap)[1]){
  for(j in 1:dim(tfallprosigtfsoverlap)[2]){
    tfallprosigtfsoverlap[i,j] <- length(intersect(tfprosiglist[[i]],tfprosiglist[[j]]))
  }
}

for(i in 1:dim(tfallprophsigtfsoverlap)[1]){
  for(j in 1:dim(tfallprophsigtfsoverlap)[2]){
    tfallprophsigtfsoverlap[i,j] <- length(intersect(tfprophsiglist[[i]],tfprophsiglist[[j]]))
  }
}

png(file = "Supplemental Figure S20D_021325.png",width = 6.5, height = 3,units = "in",res = 600)
pheatmap(tfallrnasigtfsoverlap,cluster_rows = F,cluster_cols = F,color = colorpanel(101,"white","firebrick"),display_numbers = T,number_color = "black",number_format = "%.0f",angle_col = 0,fontsize = 15,show_colnames = F,annotation_col = newrnatissuemeta,annotation_colors = ann_cols,border_color = NA)
dev.off()

png(file = "Supplemental Figure S20E_021325.png",width = 6.5, height = 3,units = "in",res = 600)
pheatmap(tfallprosigtfsoverlap,cluster_rows = F,cluster_cols = F,color = colorpanel(101,"white","firebrick"),display_numbers = T,number_color = "black",number_format = "%.0f",angle_col = 0,fontsize = 15,show_colnames = F,annotation_col = newprotissuemeta,annotation_colors = ann_cols,border_color = NA)
dev.off()

png(file = "Supplemental Figure S20F_021325.png",width = 6.5, height = 3,units = "in",res = 600)
pheatmap(tfallprophsigtfsoverlap,cluster_rows = F,cluster_cols = F,color = colorpanel(101,"white","firebrick"),display_numbers = T,number_color = "black",number_format = "%.0f",angle_col = 0,fontsize = 15,show_colnames = F,annotation_col = newprotissuemeta,annotation_colors = ann_cols,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S20D_021325.pdf",width = 6.5, height = 3)
pheatmap(tfallrnasigtfsoverlap,cluster_rows = F,cluster_cols = F,color = colorpanel(101,"white","firebrick"),display_numbers = T,number_color = "black",number_format = "%.0f",angle_col = 0,fontsize = 15,show_colnames = F,annotation_col = newrnatissuemeta,annotation_colors = ann_cols,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S20E_021325.pdf",width = 6.5, height = 3)
pheatmap(tfallprosigtfsoverlap,cluster_rows = F,cluster_cols = F,color = colorpanel(101,"white","firebrick"),display_numbers = T,number_color = "black",number_format = "%.0f",angle_col = 0,fontsize = 15,show_colnames = F,annotation_col = newprotissuemeta,annotation_colors = ann_cols,border_color = NA)
dev.off()

pdf(file = "Supplemental Figure S20F_021325.pdf",width = 6.5, height = 3)
pheatmap(tfallprophsigtfsoverlap,cluster_rows = F,cluster_cols = F,color = colorpanel(101,"white","firebrick"),display_numbers = T,number_color = "black",number_format = "%.0f",angle_col = 0,fontsize = 15,show_colnames = F,annotation_col = newprotissuemeta,annotation_colors = ann_cols,border_color = NA)
dev.off()



save.image("Figure5_021325.RData")


#####
# Now we are adding significance tests for tftargetcorrelationstructurecombo
####

gastrotftargetcorrelationstructurecombo <- tftargetcorrelationstructurecombo[tftargetcorrelationstructurecombo$Tissue %in% "SKM-GN",]
hearttftargetcorrelationstructurecombo <- tftargetcorrelationstructurecombo[tftargetcorrelationstructurecombo$Tissue %in% "HEART",]
kidneytftargetcorrelationstructurecombo <- tftargetcorrelationstructurecombo[tftargetcorrelationstructurecombo$Tissue %in% "KIDNEY",]
livertftargetcorrelationstructurecombo <- tftargetcorrelationstructurecombo[tftargetcorrelationstructurecombo$Tissue %in% "LIVER",]
lungtftargetcorrelationstructurecombo <- tftargetcorrelationstructurecombo[tftargetcorrelationstructurecombo$Tissue %in% "LUNG",]
whitetftargetcorrelationstructurecombo <- tftargetcorrelationstructurecombo[tftargetcorrelationstructurecombo$Tissue %in% "WAT-SC",]

gastrotftargetcorrelationstructurecombo$AbsCorrelation <- abs(gastrotftargetcorrelationstructurecombo$Correlation)
gastrotftargetcorrelationstructurecombo$p.adj <- 1
gastrotftargetcorrelationstructurecombo$p.signif <- "ns"
for(i in 1:dim(gastrotftargetcorrelationstructurecombo)[1]){
  ourtf_ome <- gastrotftargetcorrelationstructurecombo$TF_Ome[i]
  ourtest <- compare_means(AbsCorrelation ~ Training.Response,data = gastrotftargetcorrelationstructurecombo[gastrotftargetcorrelationstructurecombo$TF_Ome %in% ourtf_ome,])
  gastrotftargetcorrelationstructurecombo[i,"p.adj"] <- ourtest$p.adj
  gastrotftargetcorrelationstructurecombo[i,"p.signif"] <- ourtest$p.signif
}

hearttftargetcorrelationstructurecombo$AbsCorrelation <- abs(hearttftargetcorrelationstructurecombo$Correlation)
hearttftargetcorrelationstructurecombo$p.adj <- 1
hearttftargetcorrelationstructurecombo$p.signif <- "ns"
for(i in 1:dim(hearttftargetcorrelationstructurecombo)[1]){
  ourtf_ome <- hearttftargetcorrelationstructurecombo$TF_Ome[i]
  ourtest <- compare_means(AbsCorrelation ~ Training.Response,data = hearttftargetcorrelationstructurecombo[hearttftargetcorrelationstructurecombo$TF_Ome %in% ourtf_ome,])
  hearttftargetcorrelationstructurecombo[i,"p.adj"] <- ourtest$p.adj
  hearttftargetcorrelationstructurecombo[i,"p.signif"] <- ourtest$p.signif
}

kidneytftargetcorrelationstructurecombo$AbsCorrelation <- abs(kidneytftargetcorrelationstructurecombo$Correlation)
kidneytftargetcorrelationstructurecombo$p.adj <- 1
kidneytftargetcorrelationstructurecombo$p.signif <- "ns"
for(i in 1:dim(kidneytftargetcorrelationstructurecombo)[1]){
  ourtf_ome <- kidneytftargetcorrelationstructurecombo$TF_Ome[i]
  if(sum(kidneytftargetcorrelationstructurecombo$TF_Ome %in% ourtf_ome & kidneytftargetcorrelationstructurecombo$Training.Response %in% "DEG") > 0){
    ourtest <- compare_means(AbsCorrelation ~ Training.Response,data = kidneytftargetcorrelationstructurecombo[kidneytftargetcorrelationstructurecombo$TF_Ome %in% ourtf_ome,])
    kidneytftargetcorrelationstructurecombo[i,"p.adj"] <- ourtest$p.adj
    kidneytftargetcorrelationstructurecombo[i,"p.signif"] <- ourtest$p.signif 
  }
}

livertftargetcorrelationstructurecombo$AbsCorrelation <- abs(livertftargetcorrelationstructurecombo$Correlation)
livertftargetcorrelationstructurecombo$p.adj <- 1
livertftargetcorrelationstructurecombo$p.signif <- "ns"
for(i in 1:dim(livertftargetcorrelationstructurecombo)[1]){
  ourtf_ome <- livertftargetcorrelationstructurecombo$TF_Ome[i]
  ourtest <- compare_means(AbsCorrelation ~ Training.Response,data = livertftargetcorrelationstructurecombo[livertftargetcorrelationstructurecombo$TF_Ome %in% ourtf_ome,])
  livertftargetcorrelationstructurecombo[i,"p.adj"] <- ourtest$p.adj
  livertftargetcorrelationstructurecombo[i,"p.signif"] <- ourtest$p.signif
}

lungtftargetcorrelationstructurecombo$AbsCorrelation <- abs(lungtftargetcorrelationstructurecombo$Correlation)
lungtftargetcorrelationstructurecombo$p.adj <- 1
lungtftargetcorrelationstructurecombo$p.signif <- "ns"
for(i in 1:dim(lungtftargetcorrelationstructurecombo)[1]){
  ourtf_ome <- lungtftargetcorrelationstructurecombo$TF_Ome[i]
  ourtest <- compare_means(AbsCorrelation ~ Training.Response,data = lungtftargetcorrelationstructurecombo[lungtftargetcorrelationstructurecombo$TF_Ome %in% ourtf_ome,])
  lungtftargetcorrelationstructurecombo[i,"p.adj"] <- ourtest$p.adj
  lungtftargetcorrelationstructurecombo[i,"p.signif"] <- ourtest$p.signif
}

whitetftargetcorrelationstructurecombo$AbsCorrelation <- abs(whitetftargetcorrelationstructurecombo$Correlation)
whitetftargetcorrelationstructurecombo$p.adj <- 1
whitetftargetcorrelationstructurecombo$p.signif <- "ns"
for(i in 1:dim(whitetftargetcorrelationstructurecombo)[1]){
  ourtf_ome <- whitetftargetcorrelationstructurecombo$TF_Ome[i]
  ourtest <- compare_means(AbsCorrelation ~ Training.Response,data = whitetftargetcorrelationstructurecombo[whitetftargetcorrelationstructurecombo$TF_Ome %in% ourtf_ome,])
  whitetftargetcorrelationstructurecombo[i,"p.adj"] <- ourtest$p.adj
  whitetftargetcorrelationstructurecombo[i,"p.signif"] <- ourtest$p.signif
}

gastrostat.test <- compare_means(
  AbsCorrelation ~ Training.Response, data = gastrotftargetcorrelationstructurecombo, group.by = "TF_Ome",
  method = "t.test"
)
png(file = "Supplemental Figure S26A_021325.png",width = 10,height = 6,units = "in",res = 600)
ggdotplot(gastrotftargetcorrelationstructurecombo,x = "TF_Ome",y = "Correlation",size = 0.5, color = "Training.Response",fill = "Training.Response",binwidth = 0.03) + geom_hline(yintercept = 0) + ggtitle("SKM-GN EET Significant TFs") + stat_pvalue_manual(gastrostat.test, x = "TF_Ome",y.position = 1, label = "p.signif")
dev.off()
pdf(file = "Supplemental Figure S26A_021325.pdf",width = 10,height = 6)
ggdotplot(gastrotftargetcorrelationstructurecombo,x = "TF_Ome",y = "Correlation",size = 0.5, color = "Training.Response",fill = "Training.Response",binwidth = 0.03) + geom_hline(yintercept = 0) + ggtitle("SKM-GN EET Significant TFs") + stat_pvalue_manual(gastrostat.test, x = "TF_Ome",y.position = 1, label = "p.signif")
dev.off()

heartstat.test <- compare_means(
  AbsCorrelation ~ Training.Response, data = hearttftargetcorrelationstructurecombo, group.by = "TF_Ome",
  method = "t.test"
)
png(file = "Supplemental Figure S26B_021325.png",width = 10,height = 6,units = "in",res = 600)
ggdotplot(hearttftargetcorrelationstructurecombo,x = "TF_Ome",y = "Correlation",size = 0.5, color = "Training.Response",fill = "Training.Response",binwidth = 0.03) + geom_hline(yintercept = 0) + ggtitle("HEART EET Significant TFs") + stat_pvalue_manual(heartstat.test, x = "TF_Ome",y.position = 1, label = "p.signif")
dev.off()
pdf(file = "Supplemental Figure S26B_021325.pdf",width = 10,height = 6)
ggdotplot(hearttftargetcorrelationstructurecombo,x = "TF_Ome",y = "Correlation",size = 0.5, color = "Training.Response",fill = "Training.Response",binwidth = 0.03) + geom_hline(yintercept = 0) + ggtitle("HEART EET Significant TFs") + stat_pvalue_manual(heartstat.test, x = "TF_Ome",y.position = 1, label = "p.signif")
dev.off()

kidneystat.test <- compare_means(
  AbsCorrelation ~ Training.Response, data = kidneytftargetcorrelationstructurecombo, group.by = "TF_Ome",
  method = "t.test"
)
png(file = "Supplemental Figure S26C_021325.png",width = 10,height = 6,units = "in",res = 600)
ggdotplot(kidneytftargetcorrelationstructurecombo,x = "TF_Ome",y = "Correlation",size = 0.5, color = "Training.Response",fill = "Training.Response",binwidth = 0.03) + geom_hline(yintercept = 0) + ggtitle("KIDNEY EET Significant TFs") + stat_pvalue_manual(kidneystat.test, x = "TF_Ome",y.position = 1, label = "p.signif")
dev.off()
pdf(file = "Supplemental Figure S26C_021325.pdf",width = 10,height = 6)
ggdotplot(kidneytftargetcorrelationstructurecombo,x = "TF_Ome",y = "Correlation",size = 0.5, color = "Training.Response",fill = "Training.Response",binwidth = 0.03) + geom_hline(yintercept = 0) + ggtitle("KIDNEY EET Significant TFs") + stat_pvalue_manual(kidneystat.test, x = "TF_Ome",y.position = 1, label = "p.signif")
dev.off()

liverstat.test <- compare_means(
  AbsCorrelation ~ Training.Response, data = livertftargetcorrelationstructurecombo, group.by = "TF_Ome",
  method = "t.test"
)
png(file = "Supplemental Figure S26D_021325.png",width = 10,height = 6,units = "in",res = 600)
ggdotplot(livertftargetcorrelationstructurecombo,x = "TF_Ome",y = "Correlation",size = 0.5, color = "Training.Response",fill = "Training.Response",binwidth = 0.03) + geom_hline(yintercept = 0) + ggtitle("LIVER EET Significant TFs") + stat_pvalue_manual(liverstat.test, x = "TF_Ome",y.position = 1, label = "p.signif")
dev.off()
pdf(file = "Supplemental Figure S26D_021325.pdf",width = 10,height = 6)
ggdotplot(livertftargetcorrelationstructurecombo,x = "TF_Ome",y = "Correlation",size = 0.5, color = "Training.Response",fill = "Training.Response",binwidth = 0.03) + geom_hline(yintercept = 0) + ggtitle("LIVER EET Significant TFs") + stat_pvalue_manual(liverstat.test, x = "TF_Ome",y.position = 1, label = "p.signif")
dev.off()


lungstat.test <- compare_means(
  AbsCorrelation ~ Training.Response, data = lungtftargetcorrelationstructurecombo, group.by = "TF_Ome",
  method = "t.test"
)
png(file = "Supplemental Figure S26E_021325.png",width = 10,height = 6,units = "in",res = 600)
ggdotplot(lungtftargetcorrelationstructurecombo,x = "TF_Ome",y = "Correlation",size = 0.5, color = "Training.Response",fill = "Training.Response",binwidth = 0.03) + geom_hline(yintercept = 0) + ggtitle("LUNG EET Significant TFs") + stat_pvalue_manual(lungstat.test, x = "TF_Ome",y.position = 1, label = "p.signif")
dev.off()
pdf(file = "Supplemental Figure S26E_021325.pdf",width = 10,height = 6)
ggdotplot(lungtftargetcorrelationstructurecombo,x = "TF_Ome",y = "Correlation",size = 0.5, color = "Training.Response",fill = "Training.Response",binwidth = 0.03) + geom_hline(yintercept = 0) + ggtitle("LUNG EET Significant TFs") + stat_pvalue_manual(lungstat.test, x = "TF_Ome",y.position = 1, label = "p.signif")
dev.off()

png(file = "Figure 5E_021325.png",width = 12,height = 6,units = "in",res = 600)
ggdotplot(lungtftargetcorrelationstructurecombo,x = "TF_Ome",y = "Correlation",size = 0.5, color = "Training.Response",fill = "Training.Response",binwidth = 0.03) + geom_hline(yintercept = 0) + ggtitle("LUNG EET Significant TFs") + stat_pvalue_manual(lungstat.test, x = "TF_Ome",y.position = 1, label = "p.signif") + theme(text = element_text(size = 15))
dev.off()
pdf(file = "Figure 5E_021325.pdf",width = 12,height = 6)
ggdotplot(lungtftargetcorrelationstructurecombo,x = "TF_Ome",y = "Correlation",size = 0.5, color = "Training.Response",fill = "Training.Response",binwidth = 0.03) + geom_hline(yintercept = 0) + ggtitle("LUNG EET Significant TFs") + stat_pvalue_manual(lungstat.test, x = "TF_Ome",y.position = 1, label = "p.signif") + theme(text = element_text(size = 15))
dev.off()

whitestat.test <- compare_means(
  AbsCorrelation ~ Training.Response, data = whitetftargetcorrelationstructurecombo, group.by = "TF_Ome",
  method = "t.test"
)
png(file = "Supplemental Figure S26F_021325.png",width = 10,height = 6,units = "in",res = 600)
ggdotplot(whitetftargetcorrelationstructurecombo,x = "TF_Ome",y = "Correlation",size = 0.5, color = "Training.Response",fill = "Training.Response",binwidth = 0.03) + geom_hline(yintercept = 0) + ggtitle("WAT-SC EET Significant TFs") + stat_pvalue_manual(whitestat.test, x = "TF_Ome",y.position = 1, label = "p.signif")
dev.off()
pdf(file = "Supplemental Figure S26F_021325.pdf",width = 10,height = 6)
ggdotplot(whitetftargetcorrelationstructurecombo,x = "TF_Ome",y = "Correlation",size = 0.5, color = "Training.Response",fill = "Training.Response",binwidth = 0.03) + geom_hline(yintercept = 0) + ggtitle("WAT-SC EET Significant TFs") + stat_pvalue_manual(whitestat.test, x = "TF_Ome",y.position = 1, label = "p.signif")
dev.off()


#####
# An aside with the immune data, specifically in lung
####

lungimmunel2fc <- matrix(0L,nrow = length(unique(MotrpacRatTraining6moData::IMMUNO_LUNG_DA$feature_ID)),ncol = 8)
rownames(lungimmunel2fc) <- unique(MotrpacRatTraining6moData::IMMUNO_LUNG_DA$feature_ID)
colnames(lungimmunel2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
for(i in 1:dim(lungimmunel2fc)[1]){
  ourimmune <- rownames(lungimmunel2fc)[i]
  lungimmunel2fc[i,"F W1"] <- MotrpacRatTraining6moData::IMMUNO_LUNG_DA[MotrpacRatTraining6moData::IMMUNO_LUNG_DA$feature_ID %in% ourimmune & MotrpacRatTraining6moData::IMMUNO_LUNG_DA$sex %in% "female" & MotrpacRatTraining6moData::IMMUNO_LUNG_DA$comparison_group %in% "1w","logFC"]
  lungimmunel2fc[i,"F W2"] <- MotrpacRatTraining6moData::IMMUNO_LUNG_DA[MotrpacRatTraining6moData::IMMUNO_LUNG_DA$feature_ID %in% ourimmune & MotrpacRatTraining6moData::IMMUNO_LUNG_DA$sex %in% "female" & MotrpacRatTraining6moData::IMMUNO_LUNG_DA$comparison_group %in% "2w","logFC"]
  lungimmunel2fc[i,"F W4"] <- MotrpacRatTraining6moData::IMMUNO_LUNG_DA[MotrpacRatTraining6moData::IMMUNO_LUNG_DA$feature_ID %in% ourimmune & MotrpacRatTraining6moData::IMMUNO_LUNG_DA$sex %in% "female" & MotrpacRatTraining6moData::IMMUNO_LUNG_DA$comparison_group %in% "4w","logFC"]
  lungimmunel2fc[i,"F W8"] <- MotrpacRatTraining6moData::IMMUNO_LUNG_DA[MotrpacRatTraining6moData::IMMUNO_LUNG_DA$feature_ID %in% ourimmune & MotrpacRatTraining6moData::IMMUNO_LUNG_DA$sex %in% "female" & MotrpacRatTraining6moData::IMMUNO_LUNG_DA$comparison_group %in% "8w","logFC"]
  lungimmunel2fc[i,"M W1"] <- MotrpacRatTraining6moData::IMMUNO_LUNG_DA[MotrpacRatTraining6moData::IMMUNO_LUNG_DA$feature_ID %in% ourimmune & MotrpacRatTraining6moData::IMMUNO_LUNG_DA$sex %in% "male" & MotrpacRatTraining6moData::IMMUNO_LUNG_DA$comparison_group %in% "1w","logFC"]
  lungimmunel2fc[i,"M W2"] <- MotrpacRatTraining6moData::IMMUNO_LUNG_DA[MotrpacRatTraining6moData::IMMUNO_LUNG_DA$feature_ID %in% ourimmune & MotrpacRatTraining6moData::IMMUNO_LUNG_DA$sex %in% "male" & MotrpacRatTraining6moData::IMMUNO_LUNG_DA$comparison_group %in% "2w","logFC"]
  lungimmunel2fc[i,"M W4"] <- MotrpacRatTraining6moData::IMMUNO_LUNG_DA[MotrpacRatTraining6moData::IMMUNO_LUNG_DA$feature_ID %in% ourimmune & MotrpacRatTraining6moData::IMMUNO_LUNG_DA$sex %in% "male" & MotrpacRatTraining6moData::IMMUNO_LUNG_DA$comparison_group %in% "4w","logFC"]
  lungimmunel2fc[i,"M W8"] <- MotrpacRatTraining6moData::IMMUNO_LUNG_DA[MotrpacRatTraining6moData::IMMUNO_LUNG_DA$feature_ID %in% ourimmune & MotrpacRatTraining6moData::IMMUNO_LUNG_DA$sex %in% "male" & MotrpacRatTraining6moData::IMMUNO_LUNG_DA$comparison_group %in% "8w","logFC"]
}


#png(file = "lung_immunel2fcheatmap.png",width = 6,height = 6,units = "in",res = 600)
#pheatmap(lungimmunel2fc,cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = 0,display_numbers = T,number_color = "black")
#dev.off()

save.image("Figure5_021325.RData")


immunoplasmada <- MotrpacRatTraining6moData::IMMUNO_PLASMA_DA
immunoplasmada$feature_ID_dataset <- paste(immunoplasmada$feature_ID,"_",immunoplasmada$dataset,sep = "")

plasmaimmunel2fc <- matrix(0L,nrow = length(unique(immunoplasmada$feature_ID_dataset)),ncol = 8)
rownames(plasmaimmunel2fc) <- unique(immunoplasmada$feature_ID_dataset)
colnames(plasmaimmunel2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
for(i in 1:dim(plasmaimmunel2fc)[1]){
  ourimmune <- rownames(plasmaimmunel2fc)[i]
  plasmaimmunel2fc[i,"F W1"] <- immunoplasmada[immunoplasmada$feature_ID_dataset %in% ourimmune & immunoplasmada$sex %in% "female" & immunoplasmada$comparison_group %in% "1w","logFC"]
  plasmaimmunel2fc[i,"F W2"] <- immunoplasmada[immunoplasmada$feature_ID_dataset %in% ourimmune & immunoplasmada$sex %in% "female" & immunoplasmada$comparison_group %in% "2w","logFC"]
  plasmaimmunel2fc[i,"F W4"] <- immunoplasmada[immunoplasmada$feature_ID_dataset %in% ourimmune & immunoplasmada$sex %in% "female" & immunoplasmada$comparison_group %in% "4w","logFC"]
  plasmaimmunel2fc[i,"F W8"] <- immunoplasmada[immunoplasmada$feature_ID_dataset %in% ourimmune & immunoplasmada$sex %in% "female" & immunoplasmada$comparison_group %in% "8w","logFC"]
  plasmaimmunel2fc[i,"M W1"] <- immunoplasmada[immunoplasmada$feature_ID_dataset %in% ourimmune & immunoplasmada$sex %in% "male" & immunoplasmada$comparison_group %in% "1w","logFC"]
  plasmaimmunel2fc[i,"M W2"] <- immunoplasmada[immunoplasmada$feature_ID_dataset %in% ourimmune & immunoplasmada$sex %in% "male" & immunoplasmada$comparison_group %in% "2w","logFC"]
  plasmaimmunel2fc[i,"M W4"] <- immunoplasmada[immunoplasmada$feature_ID_dataset %in% ourimmune & immunoplasmada$sex %in% "male" & immunoplasmada$comparison_group %in% "4w","logFC"]
  plasmaimmunel2fc[i,"M W8"] <- immunoplasmada[immunoplasmada$feature_ID_dataset %in% ourimmune & immunoplasmada$sex %in% "male" & immunoplasmada$comparison_group %in% "8w","logFC"]
}


plasmaimmunepval <- matrix(0L,nrow = length(unique(immunoplasmada$feature_ID_dataset)),ncol = 8)
rownames(plasmaimmunepval) <- unique(immunoplasmada$feature_ID_dataset)
colnames(plasmaimmunepval) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
for(i in 1:dim(plasmaimmunepval)[1]){
  ourimmune <- rownames(plasmaimmunepval)[i]
  plasmaimmunepval[i,"F W1"] <- immunoplasmada[immunoplasmada$feature_ID_dataset %in% ourimmune & immunoplasmada$sex %in% "female" & immunoplasmada$comparison_group %in% "1w","p_value"]
  plasmaimmunepval[i,"F W2"] <- immunoplasmada[immunoplasmada$feature_ID_dataset %in% ourimmune & immunoplasmada$sex %in% "female" & immunoplasmada$comparison_group %in% "2w","p_value"]
  plasmaimmunepval[i,"F W4"] <- immunoplasmada[immunoplasmada$feature_ID_dataset %in% ourimmune & immunoplasmada$sex %in% "female" & immunoplasmada$comparison_group %in% "4w","p_value"]
  plasmaimmunepval[i,"F W8"] <- immunoplasmada[immunoplasmada$feature_ID_dataset %in% ourimmune & immunoplasmada$sex %in% "female" & immunoplasmada$comparison_group %in% "8w","p_value"]
  plasmaimmunepval[i,"M W1"] <- immunoplasmada[immunoplasmada$feature_ID_dataset %in% ourimmune & immunoplasmada$sex %in% "male" & immunoplasmada$comparison_group %in% "1w","p_value"]
  plasmaimmunepval[i,"M W2"] <- immunoplasmada[immunoplasmada$feature_ID_dataset %in% ourimmune & immunoplasmada$sex %in% "male" & immunoplasmada$comparison_group %in% "2w","p_value"]
  plasmaimmunepval[i,"M W4"] <- immunoplasmada[immunoplasmada$feature_ID_dataset %in% ourimmune & immunoplasmada$sex %in% "male" & immunoplasmada$comparison_group %in% "4w","p_value"]
  plasmaimmunepval[i,"M W8"] <- immunoplasmada[immunoplasmada$feature_ID_dataset %in% ourimmune & immunoplasmada$sex %in% "male" & immunoplasmada$comparison_group %in% "8w","p_value"]
}


#png(file = "plasma_immunel2fcheatmap.png",width = 8,height = 8,units = "in",res = 600)
#pheatmap(plasmaimmunel2fc,cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = 0,display_numbers = T,number_color = "black",cluster_rows = F)
#dev.off()

#png(file = "plasma_immunepvalheatmap.png",width = 8,height = 8,units = "in",res = 600)
#pheatmap(-log10(plasmaimmunepval),cluster_cols = F,color = colorpanel(101,"white","firebrick"),angle_col = 0,display_numbers = T,number_color = "black",cluster_rows = F)
#dev.off()

# gastro
immunogastroda <- MotrpacRatTraining6moData::IMMUNO_SKMGN_DA
immunogastroda$feature_ID_dataset <- paste(immunogastroda$feature_ID,"_",immunogastroda$dataset,sep = "")

gastroimmunel2fc <- matrix(0L,nrow = length(unique(immunogastroda$feature_ID_dataset)),ncol = 8)
rownames(gastroimmunel2fc) <- unique(immunogastroda$feature_ID_dataset)
colnames(gastroimmunel2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
for(i in 1:dim(gastroimmunel2fc)[1]){
  ourimmune <- rownames(gastroimmunel2fc)[i]
  gastroimmunel2fc[i,"F W1"] <- immunogastroda[immunogastroda$feature_ID_dataset %in% ourimmune & immunogastroda$sex %in% "female" & immunogastroda$comparison_group %in% "1w","logFC"]
  gastroimmunel2fc[i,"F W2"] <- immunogastroda[immunogastroda$feature_ID_dataset %in% ourimmune & immunogastroda$sex %in% "female" & immunogastroda$comparison_group %in% "2w","logFC"]
  gastroimmunel2fc[i,"F W4"] <- immunogastroda[immunogastroda$feature_ID_dataset %in% ourimmune & immunogastroda$sex %in% "female" & immunogastroda$comparison_group %in% "4w","logFC"]
  gastroimmunel2fc[i,"F W8"] <- immunogastroda[immunogastroda$feature_ID_dataset %in% ourimmune & immunogastroda$sex %in% "female" & immunogastroda$comparison_group %in% "8w","logFC"]
  gastroimmunel2fc[i,"M W1"] <- immunogastroda[immunogastroda$feature_ID_dataset %in% ourimmune & immunogastroda$sex %in% "male" & immunogastroda$comparison_group %in% "1w","logFC"]
  gastroimmunel2fc[i,"M W2"] <- immunogastroda[immunogastroda$feature_ID_dataset %in% ourimmune & immunogastroda$sex %in% "male" & immunogastroda$comparison_group %in% "2w","logFC"]
  gastroimmunel2fc[i,"M W4"] <- immunogastroda[immunogastroda$feature_ID_dataset %in% ourimmune & immunogastroda$sex %in% "male" & immunogastroda$comparison_group %in% "4w","logFC"]
  gastroimmunel2fc[i,"M W8"] <- immunogastroda[immunogastroda$feature_ID_dataset %in% ourimmune & immunogastroda$sex %in% "male" & immunogastroda$comparison_group %in% "8w","logFC"]
}


gastroimmunepval <- matrix(0L,nrow = length(unique(immunogastroda$feature_ID_dataset)),ncol = 8)
rownames(gastroimmunepval) <- unique(immunogastroda$feature_ID_dataset)
colnames(gastroimmunepval) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
for(i in 1:dim(gastroimmunepval)[1]){
  ourimmune <- rownames(gastroimmunepval)[i]
  gastroimmunepval[i,"F W1"] <- immunogastroda[immunogastroda$feature_ID_dataset %in% ourimmune & immunogastroda$sex %in% "female" & immunogastroda$comparison_group %in% "1w","p_value"]
  gastroimmunepval[i,"F W2"] <- immunogastroda[immunogastroda$feature_ID_dataset %in% ourimmune & immunogastroda$sex %in% "female" & immunogastroda$comparison_group %in% "2w","p_value"]
  gastroimmunepval[i,"F W4"] <- immunogastroda[immunogastroda$feature_ID_dataset %in% ourimmune & immunogastroda$sex %in% "female" & immunogastroda$comparison_group %in% "4w","p_value"]
  gastroimmunepval[i,"F W8"] <- immunogastroda[immunogastroda$feature_ID_dataset %in% ourimmune & immunogastroda$sex %in% "female" & immunogastroda$comparison_group %in% "8w","p_value"]
  gastroimmunepval[i,"M W1"] <- immunogastroda[immunogastroda$feature_ID_dataset %in% ourimmune & immunogastroda$sex %in% "male" & immunogastroda$comparison_group %in% "1w","p_value"]
  gastroimmunepval[i,"M W2"] <- immunogastroda[immunogastroda$feature_ID_dataset %in% ourimmune & immunogastroda$sex %in% "male" & immunogastroda$comparison_group %in% "2w","p_value"]
  gastroimmunepval[i,"M W4"] <- immunogastroda[immunogastroda$feature_ID_dataset %in% ourimmune & immunogastroda$sex %in% "male" & immunogastroda$comparison_group %in% "4w","p_value"]
  gastroimmunepval[i,"M W8"] <- immunogastroda[immunogastroda$feature_ID_dataset %in% ourimmune & immunogastroda$sex %in% "male" & immunogastroda$comparison_group %in% "8w","p_value"]
}


#png(file = "gastro_immunel2fcheatmap.png",width = 8,height = 8,units = "in",res = 600)
#pheatmap(gastroimmunel2fc,cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = 0,display_numbers = T,number_color = "black",cluster_rows = F)
#dev.off()

#png(file = "gastro_immunepvalheatmap.png",width = 8,height = 8,units = "in",res = 600)
#pheatmap(-log10(gastroimmunepval),cluster_cols = F,color = colorpanel(101,"white","firebrick"),angle_col = 0,display_numbers = T,number_color = "black",cluster_rows = F)
#dev.off()


# lung
immunolungda <- MotrpacRatTraining6moData::IMMUNO_LUNG_DA
immunolungda$feature_ID_dataset <- paste(immunolungda$feature_ID,"_",immunolungda$dataset,sep = "")

lungimmunel2fc <- matrix(0L,nrow = length(unique(immunolungda$feature_ID_dataset)),ncol = 8)
rownames(lungimmunel2fc) <- unique(immunolungda$feature_ID_dataset)
colnames(lungimmunel2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
for(i in 1:dim(lungimmunel2fc)[1]){
  ourimmune <- rownames(lungimmunel2fc)[i]
  lungimmunel2fc[i,"F W1"] <- immunolungda[immunolungda$feature_ID_dataset %in% ourimmune & immunolungda$sex %in% "female" & immunolungda$comparison_group %in% "1w","logFC"]
  lungimmunel2fc[i,"F W2"] <- immunolungda[immunolungda$feature_ID_dataset %in% ourimmune & immunolungda$sex %in% "female" & immunolungda$comparison_group %in% "2w","logFC"]
  lungimmunel2fc[i,"F W4"] <- immunolungda[immunolungda$feature_ID_dataset %in% ourimmune & immunolungda$sex %in% "female" & immunolungda$comparison_group %in% "4w","logFC"]
  lungimmunel2fc[i,"F W8"] <- immunolungda[immunolungda$feature_ID_dataset %in% ourimmune & immunolungda$sex %in% "female" & immunolungda$comparison_group %in% "8w","logFC"]
  lungimmunel2fc[i,"M W1"] <- immunolungda[immunolungda$feature_ID_dataset %in% ourimmune & immunolungda$sex %in% "male" & immunolungda$comparison_group %in% "1w","logFC"]
  lungimmunel2fc[i,"M W2"] <- immunolungda[immunolungda$feature_ID_dataset %in% ourimmune & immunolungda$sex %in% "male" & immunolungda$comparison_group %in% "2w","logFC"]
  lungimmunel2fc[i,"M W4"] <- immunolungda[immunolungda$feature_ID_dataset %in% ourimmune & immunolungda$sex %in% "male" & immunolungda$comparison_group %in% "4w","logFC"]
  lungimmunel2fc[i,"M W8"] <- immunolungda[immunolungda$feature_ID_dataset %in% ourimmune & immunolungda$sex %in% "male" & immunolungda$comparison_group %in% "8w","logFC"]
}


lungimmunepval <- matrix(0L,nrow = length(unique(immunolungda$feature_ID_dataset)),ncol = 8)
rownames(lungimmunepval) <- unique(immunolungda$feature_ID_dataset)
colnames(lungimmunepval) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
for(i in 1:dim(lungimmunepval)[1]){
  ourimmune <- rownames(lungimmunepval)[i]
  lungimmunepval[i,"F W1"] <- immunolungda[immunolungda$feature_ID_dataset %in% ourimmune & immunolungda$sex %in% "female" & immunolungda$comparison_group %in% "1w","p_value"]
  lungimmunepval[i,"F W2"] <- immunolungda[immunolungda$feature_ID_dataset %in% ourimmune & immunolungda$sex %in% "female" & immunolungda$comparison_group %in% "2w","p_value"]
  lungimmunepval[i,"F W4"] <- immunolungda[immunolungda$feature_ID_dataset %in% ourimmune & immunolungda$sex %in% "female" & immunolungda$comparison_group %in% "4w","p_value"]
  lungimmunepval[i,"F W8"] <- immunolungda[immunolungda$feature_ID_dataset %in% ourimmune & immunolungda$sex %in% "female" & immunolungda$comparison_group %in% "8w","p_value"]
  lungimmunepval[i,"M W1"] <- immunolungda[immunolungda$feature_ID_dataset %in% ourimmune & immunolungda$sex %in% "male" & immunolungda$comparison_group %in% "1w","p_value"]
  lungimmunepval[i,"M W2"] <- immunolungda[immunolungda$feature_ID_dataset %in% ourimmune & immunolungda$sex %in% "male" & immunolungda$comparison_group %in% "2w","p_value"]
  lungimmunepval[i,"M W4"] <- immunolungda[immunolungda$feature_ID_dataset %in% ourimmune & immunolungda$sex %in% "male" & immunolungda$comparison_group %in% "4w","p_value"]
  lungimmunepval[i,"M W8"] <- immunolungda[immunolungda$feature_ID_dataset %in% ourimmune & immunolungda$sex %in% "male" & immunolungda$comparison_group %in% "8w","p_value"]
}


#png(file = "lung_immunel2fcheatmap.png",width = 8,height = 8,units = "in",res = 600)
#pheatmap(lungimmunel2fc,cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = 0,display_numbers = T,number_color = "black",cluster_rows = F)
#dev.off()

#png(file = "lung_immunepvalheatmap.png",width = 8,height = 8,units = "in",res = 600)
#pheatmap(-log10(lungimmunepval),cluster_cols = F,color = colorpanel(101,"white","firebrick"),angle_col = 0,display_numbers = T,number_color = "black",cluster_rows = F)
#dev.off()



# heart
immunoheartda <- MotrpacRatTraining6moData::IMMUNO_HEART_DA
immunoheartda$feature_ID_dataset <- paste(immunoheartda$feature_ID,"_",immunoheartda$dataset,sep = "")

heartimmunel2fc <- matrix(0L,nrow = length(unique(immunoheartda$feature_ID_dataset)),ncol = 8)
rownames(heartimmunel2fc) <- unique(immunoheartda$feature_ID_dataset)
colnames(heartimmunel2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
for(i in 1:dim(heartimmunel2fc)[1]){
  ourimmune <- rownames(heartimmunel2fc)[i]
  heartimmunel2fc[i,"F W1"] <- immunoheartda[immunoheartda$feature_ID_dataset %in% ourimmune & immunoheartda$sex %in% "female" & immunoheartda$comparison_group %in% "1w","logFC"]
  heartimmunel2fc[i,"F W2"] <- immunoheartda[immunoheartda$feature_ID_dataset %in% ourimmune & immunoheartda$sex %in% "female" & immunoheartda$comparison_group %in% "2w","logFC"]
  heartimmunel2fc[i,"F W4"] <- immunoheartda[immunoheartda$feature_ID_dataset %in% ourimmune & immunoheartda$sex %in% "female" & immunoheartda$comparison_group %in% "4w","logFC"]
  heartimmunel2fc[i,"F W8"] <- immunoheartda[immunoheartda$feature_ID_dataset %in% ourimmune & immunoheartda$sex %in% "female" & immunoheartda$comparison_group %in% "8w","logFC"]
  heartimmunel2fc[i,"M W1"] <- immunoheartda[immunoheartda$feature_ID_dataset %in% ourimmune & immunoheartda$sex %in% "male" & immunoheartda$comparison_group %in% "1w","logFC"]
  heartimmunel2fc[i,"M W2"] <- immunoheartda[immunoheartda$feature_ID_dataset %in% ourimmune & immunoheartda$sex %in% "male" & immunoheartda$comparison_group %in% "2w","logFC"]
  heartimmunel2fc[i,"M W4"] <- immunoheartda[immunoheartda$feature_ID_dataset %in% ourimmune & immunoheartda$sex %in% "male" & immunoheartda$comparison_group %in% "4w","logFC"]
  heartimmunel2fc[i,"M W8"] <- immunoheartda[immunoheartda$feature_ID_dataset %in% ourimmune & immunoheartda$sex %in% "male" & immunoheartda$comparison_group %in% "8w","logFC"]
}


heartimmunepval <- matrix(0L,nrow = length(unique(immunoheartda$feature_ID_dataset)),ncol = 8)
rownames(heartimmunepval) <- unique(immunoheartda$feature_ID_dataset)
colnames(heartimmunepval) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
for(i in 1:dim(heartimmunepval)[1]){
  ourimmune <- rownames(heartimmunepval)[i]
  heartimmunepval[i,"F W1"] <- immunoheartda[immunoheartda$feature_ID_dataset %in% ourimmune & immunoheartda$sex %in% "female" & immunoheartda$comparison_group %in% "1w","p_value"]
  heartimmunepval[i,"F W2"] <- immunoheartda[immunoheartda$feature_ID_dataset %in% ourimmune & immunoheartda$sex %in% "female" & immunoheartda$comparison_group %in% "2w","p_value"]
  heartimmunepval[i,"F W4"] <- immunoheartda[immunoheartda$feature_ID_dataset %in% ourimmune & immunoheartda$sex %in% "female" & immunoheartda$comparison_group %in% "4w","p_value"]
  heartimmunepval[i,"F W8"] <- immunoheartda[immunoheartda$feature_ID_dataset %in% ourimmune & immunoheartda$sex %in% "female" & immunoheartda$comparison_group %in% "8w","p_value"]
  heartimmunepval[i,"M W1"] <- immunoheartda[immunoheartda$feature_ID_dataset %in% ourimmune & immunoheartda$sex %in% "male" & immunoheartda$comparison_group %in% "1w","p_value"]
  heartimmunepval[i,"M W2"] <- immunoheartda[immunoheartda$feature_ID_dataset %in% ourimmune & immunoheartda$sex %in% "male" & immunoheartda$comparison_group %in% "2w","p_value"]
  heartimmunepval[i,"M W4"] <- immunoheartda[immunoheartda$feature_ID_dataset %in% ourimmune & immunoheartda$sex %in% "male" & immunoheartda$comparison_group %in% "4w","p_value"]
  heartimmunepval[i,"M W8"] <- immunoheartda[immunoheartda$feature_ID_dataset %in% ourimmune & immunoheartda$sex %in% "male" & immunoheartda$comparison_group %in% "8w","p_value"]
}


#png(file = "heart_immunel2fcheatmap.png",width = 8,height = 8,units = "in",res = 600)
#pheatmap(heartimmunel2fc,cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = 0,display_numbers = T,number_color = "black",cluster_rows = F)
#dev.off()

#png(file = "heart_immunepvalheatmap.png",width = 8,height = 8,units = "in",res = 600)
#pheatmap(-log10(heartimmunepval),cluster_cols = F,color = colorpanel(101,"white","firebrick"),angle_col = 0,display_numbers = T,number_color = "black",cluster_rows = F)
#dev.off()


# hippo
immunohippoda <- MotrpacRatTraining6moData::IMMUNO_HIPPOC_DA
immunohippoda$feature_ID_dataset <- paste(immunohippoda$feature_ID,"_",immunohippoda$dataset,sep = "")

hippoimmunel2fc <- matrix(0L,nrow = length(unique(immunohippoda$feature_ID_dataset)),ncol = 8)
rownames(hippoimmunel2fc) <- unique(immunohippoda$feature_ID_dataset)
colnames(hippoimmunel2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
for(i in 1:dim(hippoimmunel2fc)[1]){
  ourimmune <- rownames(hippoimmunel2fc)[i]
  hippoimmunel2fc[i,"F W1"] <- immunohippoda[immunohippoda$feature_ID_dataset %in% ourimmune & immunohippoda$sex %in% "female" & immunohippoda$comparison_group %in% "1w","logFC"]
  hippoimmunel2fc[i,"F W2"] <- immunohippoda[immunohippoda$feature_ID_dataset %in% ourimmune & immunohippoda$sex %in% "female" & immunohippoda$comparison_group %in% "2w","logFC"]
  hippoimmunel2fc[i,"F W4"] <- immunohippoda[immunohippoda$feature_ID_dataset %in% ourimmune & immunohippoda$sex %in% "female" & immunohippoda$comparison_group %in% "4w","logFC"]
  hippoimmunel2fc[i,"F W8"] <- immunohippoda[immunohippoda$feature_ID_dataset %in% ourimmune & immunohippoda$sex %in% "female" & immunohippoda$comparison_group %in% "8w","logFC"]
  hippoimmunel2fc[i,"M W1"] <- immunohippoda[immunohippoda$feature_ID_dataset %in% ourimmune & immunohippoda$sex %in% "male" & immunohippoda$comparison_group %in% "1w","logFC"]
  hippoimmunel2fc[i,"M W2"] <- immunohippoda[immunohippoda$feature_ID_dataset %in% ourimmune & immunohippoda$sex %in% "male" & immunohippoda$comparison_group %in% "2w","logFC"]
  hippoimmunel2fc[i,"M W4"] <- immunohippoda[immunohippoda$feature_ID_dataset %in% ourimmune & immunohippoda$sex %in% "male" & immunohippoda$comparison_group %in% "4w","logFC"]
  hippoimmunel2fc[i,"M W8"] <- immunohippoda[immunohippoda$feature_ID_dataset %in% ourimmune & immunohippoda$sex %in% "male" & immunohippoda$comparison_group %in% "8w","logFC"]
}


hippoimmunepval <- matrix(0L,nrow = length(unique(immunohippoda$feature_ID_dataset)),ncol = 8)
rownames(hippoimmunepval) <- unique(immunohippoda$feature_ID_dataset)
colnames(hippoimmunepval) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
for(i in 1:dim(hippoimmunepval)[1]){
  ourimmune <- rownames(hippoimmunepval)[i]
  hippoimmunepval[i,"F W1"] <- immunohippoda[immunohippoda$feature_ID_dataset %in% ourimmune & immunohippoda$sex %in% "female" & immunohippoda$comparison_group %in% "1w","p_value"]
  hippoimmunepval[i,"F W2"] <- immunohippoda[immunohippoda$feature_ID_dataset %in% ourimmune & immunohippoda$sex %in% "female" & immunohippoda$comparison_group %in% "2w","p_value"]
  hippoimmunepval[i,"F W4"] <- immunohippoda[immunohippoda$feature_ID_dataset %in% ourimmune & immunohippoda$sex %in% "female" & immunohippoda$comparison_group %in% "4w","p_value"]
  hippoimmunepval[i,"F W8"] <- immunohippoda[immunohippoda$feature_ID_dataset %in% ourimmune & immunohippoda$sex %in% "female" & immunohippoda$comparison_group %in% "8w","p_value"]
  hippoimmunepval[i,"M W1"] <- immunohippoda[immunohippoda$feature_ID_dataset %in% ourimmune & immunohippoda$sex %in% "male" & immunohippoda$comparison_group %in% "1w","p_value"]
  hippoimmunepval[i,"M W2"] <- immunohippoda[immunohippoda$feature_ID_dataset %in% ourimmune & immunohippoda$sex %in% "male" & immunohippoda$comparison_group %in% "2w","p_value"]
  hippoimmunepval[i,"M W4"] <- immunohippoda[immunohippoda$feature_ID_dataset %in% ourimmune & immunohippoda$sex %in% "male" & immunohippoda$comparison_group %in% "4w","p_value"]
  hippoimmunepval[i,"M W8"] <- immunohippoda[immunohippoda$feature_ID_dataset %in% ourimmune & immunohippoda$sex %in% "male" & immunohippoda$comparison_group %in% "8w","p_value"]
}


#png(file = "hippo_immunel2fcheatmap.png",width = 8,height = 8,units = "in",res = 600)
#pheatmap(hippoimmunel2fc,cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = 0,display_numbers = T,number_color = "black",cluster_rows = F)
#dev.off()

#png(file = "hippo_immunepvalheatmap.png",width = 8,height = 8,units = "in",res = 600)
#pheatmap(-log10(hippoimmunepval),cluster_cols = F,color = colorpanel(101,"white","firebrick"),angle_col = 0,display_numbers = T,number_color = "black",cluster_rows = F)
#dev.off()



# kidney
immunokidneyda <- MotrpacRatTraining6moData::IMMUNO_KIDNEY_DA
immunokidneyda$feature_ID_dataset <- paste(immunokidneyda$feature_ID,"_",immunokidneyda$dataset,sep = "")

kidneyimmunel2fc <- matrix(0L,nrow = length(unique(immunokidneyda$feature_ID_dataset)),ncol = 8)
rownames(kidneyimmunel2fc) <- unique(immunokidneyda$feature_ID_dataset)
colnames(kidneyimmunel2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
for(i in 1:dim(kidneyimmunel2fc)[1]){
  ourimmune <- rownames(kidneyimmunel2fc)[i]
  kidneyimmunel2fc[i,"F W1"] <- immunokidneyda[immunokidneyda$feature_ID_dataset %in% ourimmune & immunokidneyda$sex %in% "female" & immunokidneyda$comparison_group %in% "1w","logFC"]
  kidneyimmunel2fc[i,"F W2"] <- immunokidneyda[immunokidneyda$feature_ID_dataset %in% ourimmune & immunokidneyda$sex %in% "female" & immunokidneyda$comparison_group %in% "2w","logFC"]
  kidneyimmunel2fc[i,"F W4"] <- immunokidneyda[immunokidneyda$feature_ID_dataset %in% ourimmune & immunokidneyda$sex %in% "female" & immunokidneyda$comparison_group %in% "4w","logFC"]
  kidneyimmunel2fc[i,"F W8"] <- immunokidneyda[immunokidneyda$feature_ID_dataset %in% ourimmune & immunokidneyda$sex %in% "female" & immunokidneyda$comparison_group %in% "8w","logFC"]
  kidneyimmunel2fc[i,"M W1"] <- immunokidneyda[immunokidneyda$feature_ID_dataset %in% ourimmune & immunokidneyda$sex %in% "male" & immunokidneyda$comparison_group %in% "1w","logFC"]
  kidneyimmunel2fc[i,"M W2"] <- immunokidneyda[immunokidneyda$feature_ID_dataset %in% ourimmune & immunokidneyda$sex %in% "male" & immunokidneyda$comparison_group %in% "2w","logFC"]
  kidneyimmunel2fc[i,"M W4"] <- immunokidneyda[immunokidneyda$feature_ID_dataset %in% ourimmune & immunokidneyda$sex %in% "male" & immunokidneyda$comparison_group %in% "4w","logFC"]
  kidneyimmunel2fc[i,"M W8"] <- immunokidneyda[immunokidneyda$feature_ID_dataset %in% ourimmune & immunokidneyda$sex %in% "male" & immunokidneyda$comparison_group %in% "8w","logFC"]
}


kidneyimmunepval <- matrix(0L,nrow = length(unique(immunokidneyda$feature_ID_dataset)),ncol = 8)
rownames(kidneyimmunepval) <- unique(immunokidneyda$feature_ID_dataset)
colnames(kidneyimmunepval) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
for(i in 1:dim(kidneyimmunepval)[1]){
  ourimmune <- rownames(kidneyimmunepval)[i]
  kidneyimmunepval[i,"F W1"] <- immunokidneyda[immunokidneyda$feature_ID_dataset %in% ourimmune & immunokidneyda$sex %in% "female" & immunokidneyda$comparison_group %in% "1w","p_value"]
  kidneyimmunepval[i,"F W2"] <- immunokidneyda[immunokidneyda$feature_ID_dataset %in% ourimmune & immunokidneyda$sex %in% "female" & immunokidneyda$comparison_group %in% "2w","p_value"]
  kidneyimmunepval[i,"F W4"] <- immunokidneyda[immunokidneyda$feature_ID_dataset %in% ourimmune & immunokidneyda$sex %in% "female" & immunokidneyda$comparison_group %in% "4w","p_value"]
  kidneyimmunepval[i,"F W8"] <- immunokidneyda[immunokidneyda$feature_ID_dataset %in% ourimmune & immunokidneyda$sex %in% "female" & immunokidneyda$comparison_group %in% "8w","p_value"]
  kidneyimmunepval[i,"M W1"] <- immunokidneyda[immunokidneyda$feature_ID_dataset %in% ourimmune & immunokidneyda$sex %in% "male" & immunokidneyda$comparison_group %in% "1w","p_value"]
  kidneyimmunepval[i,"M W2"] <- immunokidneyda[immunokidneyda$feature_ID_dataset %in% ourimmune & immunokidneyda$sex %in% "male" & immunokidneyda$comparison_group %in% "2w","p_value"]
  kidneyimmunepval[i,"M W4"] <- immunokidneyda[immunokidneyda$feature_ID_dataset %in% ourimmune & immunokidneyda$sex %in% "male" & immunokidneyda$comparison_group %in% "4w","p_value"]
  kidneyimmunepval[i,"M W8"] <- immunokidneyda[immunokidneyda$feature_ID_dataset %in% ourimmune & immunokidneyda$sex %in% "male" & immunokidneyda$comparison_group %in% "8w","p_value"]
}


#png(file = "kidney_immunel2fcheatmap.png",width = 8,height = 8,units = "in",res = 600)
#pheatmap(kidneyimmunel2fc,cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = 0,display_numbers = T,number_color = "black",cluster_rows = F)
#dev.off()

#png(file = "kidney_immunepvalheatmap.png",width = 8,height = 8,units = "in",res = 600)
#pheatmap(-log10(kidneyimmunepval),cluster_cols = F,color = colorpanel(101,"white","firebrick"),angle_col = 0,display_numbers = T,number_color = "black",cluster_rows = F)
#dev.off()



# liver
immunoliverda <- MotrpacRatTraining6moData::IMMUNO_LIVER_DA
immunoliverda$feature_ID_dataset <- paste(immunoliverda$feature_ID,"_",immunoliverda$dataset,sep = "")

liverimmunel2fc <- matrix(0L,nrow = length(unique(immunoliverda$feature_ID_dataset)),ncol = 8)
rownames(liverimmunel2fc) <- unique(immunoliverda$feature_ID_dataset)
colnames(liverimmunel2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
for(i in 1:dim(liverimmunel2fc)[1]){
  ourimmune <- rownames(liverimmunel2fc)[i]
  liverimmunel2fc[i,"F W1"] <- immunoliverda[immunoliverda$feature_ID_dataset %in% ourimmune & immunoliverda$sex %in% "female" & immunoliverda$comparison_group %in% "1w","logFC"]
  liverimmunel2fc[i,"F W2"] <- immunoliverda[immunoliverda$feature_ID_dataset %in% ourimmune & immunoliverda$sex %in% "female" & immunoliverda$comparison_group %in% "2w","logFC"]
  liverimmunel2fc[i,"F W4"] <- immunoliverda[immunoliverda$feature_ID_dataset %in% ourimmune & immunoliverda$sex %in% "female" & immunoliverda$comparison_group %in% "4w","logFC"]
  liverimmunel2fc[i,"F W8"] <- immunoliverda[immunoliverda$feature_ID_dataset %in% ourimmune & immunoliverda$sex %in% "female" & immunoliverda$comparison_group %in% "8w","logFC"]
  liverimmunel2fc[i,"M W1"] <- immunoliverda[immunoliverda$feature_ID_dataset %in% ourimmune & immunoliverda$sex %in% "male" & immunoliverda$comparison_group %in% "1w","logFC"]
  liverimmunel2fc[i,"M W2"] <- immunoliverda[immunoliverda$feature_ID_dataset %in% ourimmune & immunoliverda$sex %in% "male" & immunoliverda$comparison_group %in% "2w","logFC"]
  liverimmunel2fc[i,"M W4"] <- immunoliverda[immunoliverda$feature_ID_dataset %in% ourimmune & immunoliverda$sex %in% "male" & immunoliverda$comparison_group %in% "4w","logFC"]
  liverimmunel2fc[i,"M W8"] <- immunoliverda[immunoliverda$feature_ID_dataset %in% ourimmune & immunoliverda$sex %in% "male" & immunoliverda$comparison_group %in% "8w","logFC"]
}


liverimmunepval <- matrix(0L,nrow = length(unique(immunoliverda$feature_ID_dataset)),ncol = 8)
rownames(liverimmunepval) <- unique(immunoliverda$feature_ID_dataset)
colnames(liverimmunepval) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
for(i in 1:dim(liverimmunepval)[1]){
  ourimmune <- rownames(liverimmunepval)[i]
  liverimmunepval[i,"F W1"] <- immunoliverda[immunoliverda$feature_ID_dataset %in% ourimmune & immunoliverda$sex %in% "female" & immunoliverda$comparison_group %in% "1w","p_value"]
  liverimmunepval[i,"F W2"] <- immunoliverda[immunoliverda$feature_ID_dataset %in% ourimmune & immunoliverda$sex %in% "female" & immunoliverda$comparison_group %in% "2w","p_value"]
  liverimmunepval[i,"F W4"] <- immunoliverda[immunoliverda$feature_ID_dataset %in% ourimmune & immunoliverda$sex %in% "female" & immunoliverda$comparison_group %in% "4w","p_value"]
  liverimmunepval[i,"F W8"] <- immunoliverda[immunoliverda$feature_ID_dataset %in% ourimmune & immunoliverda$sex %in% "female" & immunoliverda$comparison_group %in% "8w","p_value"]
  liverimmunepval[i,"M W1"] <- immunoliverda[immunoliverda$feature_ID_dataset %in% ourimmune & immunoliverda$sex %in% "male" & immunoliverda$comparison_group %in% "1w","p_value"]
  liverimmunepval[i,"M W2"] <- immunoliverda[immunoliverda$feature_ID_dataset %in% ourimmune & immunoliverda$sex %in% "male" & immunoliverda$comparison_group %in% "2w","p_value"]
  liverimmunepval[i,"M W4"] <- immunoliverda[immunoliverda$feature_ID_dataset %in% ourimmune & immunoliverda$sex %in% "male" & immunoliverda$comparison_group %in% "4w","p_value"]
  liverimmunepval[i,"M W8"] <- immunoliverda[immunoliverda$feature_ID_dataset %in% ourimmune & immunoliverda$sex %in% "male" & immunoliverda$comparison_group %in% "8w","p_value"]
}


#png(file = "liver_immunel2fcheatmap.png",width = 8,height = 8,units = "in",res = 600)
#pheatmap(liverimmunel2fc,cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = 0,display_numbers = T,number_color = "black",cluster_rows = F)
#dev.off()

#png(file = "liver_immunepvalheatmap.png",width = 8,height = 8,units = "in",res = 600)
#pheatmap(-log10(liverimmunepval),cluster_cols = F,color = colorpanel(101,"white","firebrick"),angle_col = 0,display_numbers = T,number_color = "black",cluster_rows = F)
#dev.off()



# brown
immunobrownda <- MotrpacRatTraining6moData::IMMUNO_BAT_DA
immunobrownda$feature_ID_dataset <- paste(immunobrownda$feature_ID,"_",immunobrownda$dataset,sep = "")

brownimmunel2fc <- matrix(0L,nrow = length(unique(immunobrownda$feature_ID_dataset)),ncol = 8)
rownames(brownimmunel2fc) <- unique(immunobrownda$feature_ID_dataset)
colnames(brownimmunel2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
for(i in 1:dim(brownimmunel2fc)[1]){
  ourimmune <- rownames(brownimmunel2fc)[i]
  brownimmunel2fc[i,"F W1"] <- immunobrownda[immunobrownda$feature_ID_dataset %in% ourimmune & immunobrownda$sex %in% "female" & immunobrownda$comparison_group %in% "1w","logFC"]
  brownimmunel2fc[i,"F W2"] <- immunobrownda[immunobrownda$feature_ID_dataset %in% ourimmune & immunobrownda$sex %in% "female" & immunobrownda$comparison_group %in% "2w","logFC"]
  brownimmunel2fc[i,"F W4"] <- immunobrownda[immunobrownda$feature_ID_dataset %in% ourimmune & immunobrownda$sex %in% "female" & immunobrownda$comparison_group %in% "4w","logFC"]
  brownimmunel2fc[i,"F W8"] <- immunobrownda[immunobrownda$feature_ID_dataset %in% ourimmune & immunobrownda$sex %in% "female" & immunobrownda$comparison_group %in% "8w","logFC"]
  brownimmunel2fc[i,"M W1"] <- immunobrownda[immunobrownda$feature_ID_dataset %in% ourimmune & immunobrownda$sex %in% "male" & immunobrownda$comparison_group %in% "1w","logFC"]
  brownimmunel2fc[i,"M W2"] <- immunobrownda[immunobrownda$feature_ID_dataset %in% ourimmune & immunobrownda$sex %in% "male" & immunobrownda$comparison_group %in% "2w","logFC"]
  brownimmunel2fc[i,"M W4"] <- immunobrownda[immunobrownda$feature_ID_dataset %in% ourimmune & immunobrownda$sex %in% "male" & immunobrownda$comparison_group %in% "4w","logFC"]
  brownimmunel2fc[i,"M W8"] <- immunobrownda[immunobrownda$feature_ID_dataset %in% ourimmune & immunobrownda$sex %in% "male" & immunobrownda$comparison_group %in% "8w","logFC"]
}


brownimmunepval <- matrix(0L,nrow = length(unique(immunobrownda$feature_ID_dataset)),ncol = 8)
rownames(brownimmunepval) <- unique(immunobrownda$feature_ID_dataset)
colnames(brownimmunepval) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
for(i in 1:dim(brownimmunepval)[1]){
  ourimmune <- rownames(brownimmunepval)[i]
  brownimmunepval[i,"F W1"] <- immunobrownda[immunobrownda$feature_ID_dataset %in% ourimmune & immunobrownda$sex %in% "female" & immunobrownda$comparison_group %in% "1w","p_value"]
  brownimmunepval[i,"F W2"] <- immunobrownda[immunobrownda$feature_ID_dataset %in% ourimmune & immunobrownda$sex %in% "female" & immunobrownda$comparison_group %in% "2w","p_value"]
  brownimmunepval[i,"F W4"] <- immunobrownda[immunobrownda$feature_ID_dataset %in% ourimmune & immunobrownda$sex %in% "female" & immunobrownda$comparison_group %in% "4w","p_value"]
  brownimmunepval[i,"F W8"] <- immunobrownda[immunobrownda$feature_ID_dataset %in% ourimmune & immunobrownda$sex %in% "female" & immunobrownda$comparison_group %in% "8w","p_value"]
  brownimmunepval[i,"M W1"] <- immunobrownda[immunobrownda$feature_ID_dataset %in% ourimmune & immunobrownda$sex %in% "male" & immunobrownda$comparison_group %in% "1w","p_value"]
  brownimmunepval[i,"M W2"] <- immunobrownda[immunobrownda$feature_ID_dataset %in% ourimmune & immunobrownda$sex %in% "male" & immunobrownda$comparison_group %in% "2w","p_value"]
  brownimmunepval[i,"M W4"] <- immunobrownda[immunobrownda$feature_ID_dataset %in% ourimmune & immunobrownda$sex %in% "male" & immunobrownda$comparison_group %in% "4w","p_value"]
  brownimmunepval[i,"M W8"] <- immunobrownda[immunobrownda$feature_ID_dataset %in% ourimmune & immunobrownda$sex %in% "male" & immunobrownda$comparison_group %in% "8w","p_value"]
}


#png(file = "brown_immunel2fcheatmap.png",width = 8,height = 8,units = "in",res = 600)
#pheatmap(brownimmunel2fc,cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = 0,display_numbers = T,number_color = "black",cluster_rows = F)
#dev.off()

#png(file = "brown_immunepvalheatmap.png",width = 8,height = 8,units = "in",res = 600)
#pheatmap(-log10(brownimmunepval),cluster_cols = F,color = colorpanel(101,"white","firebrick"),angle_col = 0,display_numbers = T,number_color = "black",cluster_rows = F)
#dev.off()



# white
immunowhiteda <- MotrpacRatTraining6moData::IMMUNO_WATSC_DA
immunowhiteda$feature_ID_dataset <- paste(immunowhiteda$feature_ID,"_",immunowhiteda$dataset,sep = "")

whiteimmunel2fc <- matrix(0L,nrow = length(unique(immunowhiteda$feature_ID_dataset)),ncol = 8)
rownames(whiteimmunel2fc) <- unique(immunowhiteda$feature_ID_dataset)
colnames(whiteimmunel2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
for(i in 1:dim(whiteimmunel2fc)[1]){
  ourimmune <- rownames(whiteimmunel2fc)[i]
  whiteimmunel2fc[i,"F W1"] <- immunowhiteda[immunowhiteda$feature_ID_dataset %in% ourimmune & immunowhiteda$sex %in% "female" & immunowhiteda$comparison_group %in% "1w","logFC"]
  whiteimmunel2fc[i,"F W2"] <- immunowhiteda[immunowhiteda$feature_ID_dataset %in% ourimmune & immunowhiteda$sex %in% "female" & immunowhiteda$comparison_group %in% "2w","logFC"]
  whiteimmunel2fc[i,"F W4"] <- immunowhiteda[immunowhiteda$feature_ID_dataset %in% ourimmune & immunowhiteda$sex %in% "female" & immunowhiteda$comparison_group %in% "4w","logFC"]
  whiteimmunel2fc[i,"F W8"] <- immunowhiteda[immunowhiteda$feature_ID_dataset %in% ourimmune & immunowhiteda$sex %in% "female" & immunowhiteda$comparison_group %in% "8w","logFC"]
  whiteimmunel2fc[i,"M W1"] <- immunowhiteda[immunowhiteda$feature_ID_dataset %in% ourimmune & immunowhiteda$sex %in% "male" & immunowhiteda$comparison_group %in% "1w","logFC"]
  whiteimmunel2fc[i,"M W2"] <- immunowhiteda[immunowhiteda$feature_ID_dataset %in% ourimmune & immunowhiteda$sex %in% "male" & immunowhiteda$comparison_group %in% "2w","logFC"]
  whiteimmunel2fc[i,"M W4"] <- immunowhiteda[immunowhiteda$feature_ID_dataset %in% ourimmune & immunowhiteda$sex %in% "male" & immunowhiteda$comparison_group %in% "4w","logFC"]
  whiteimmunel2fc[i,"M W8"] <- immunowhiteda[immunowhiteda$feature_ID_dataset %in% ourimmune & immunowhiteda$sex %in% "male" & immunowhiteda$comparison_group %in% "8w","logFC"]
}


whiteimmunepval <- matrix(0L,nrow = length(unique(immunowhiteda$feature_ID_dataset)),ncol = 8)
rownames(whiteimmunepval) <- unique(immunowhiteda$feature_ID_dataset)
colnames(whiteimmunepval) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
for(i in 1:dim(whiteimmunepval)[1]){
  ourimmune <- rownames(whiteimmunepval)[i]
  whiteimmunepval[i,"F W1"] <- immunowhiteda[immunowhiteda$feature_ID_dataset %in% ourimmune & immunowhiteda$sex %in% "female" & immunowhiteda$comparison_group %in% "1w","p_value"]
  whiteimmunepval[i,"F W2"] <- immunowhiteda[immunowhiteda$feature_ID_dataset %in% ourimmune & immunowhiteda$sex %in% "female" & immunowhiteda$comparison_group %in% "2w","p_value"]
  whiteimmunepval[i,"F W4"] <- immunowhiteda[immunowhiteda$feature_ID_dataset %in% ourimmune & immunowhiteda$sex %in% "female" & immunowhiteda$comparison_group %in% "4w","p_value"]
  whiteimmunepval[i,"F W8"] <- immunowhiteda[immunowhiteda$feature_ID_dataset %in% ourimmune & immunowhiteda$sex %in% "female" & immunowhiteda$comparison_group %in% "8w","p_value"]
  whiteimmunepval[i,"M W1"] <- immunowhiteda[immunowhiteda$feature_ID_dataset %in% ourimmune & immunowhiteda$sex %in% "male" & immunowhiteda$comparison_group %in% "1w","p_value"]
  whiteimmunepval[i,"M W2"] <- immunowhiteda[immunowhiteda$feature_ID_dataset %in% ourimmune & immunowhiteda$sex %in% "male" & immunowhiteda$comparison_group %in% "2w","p_value"]
  whiteimmunepval[i,"M W4"] <- immunowhiteda[immunowhiteda$feature_ID_dataset %in% ourimmune & immunowhiteda$sex %in% "male" & immunowhiteda$comparison_group %in% "4w","p_value"]
  whiteimmunepval[i,"M W8"] <- immunowhiteda[immunowhiteda$feature_ID_dataset %in% ourimmune & immunowhiteda$sex %in% "male" & immunowhiteda$comparison_group %in% "8w","p_value"]
}


#png(file = "white_immunel2fcheatmap.png",width = 8,height = 8,units = "in",res = 600)
#pheatmap(whiteimmunel2fc,cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = 0,display_numbers = T,number_color = "black",cluster_rows = F)
#dev.off()

#png(file = "white_immunepvalheatmap.png",width = 8,height = 8,units = "in",res = 600)
#pheatmap(-log10(whiteimmunepval),cluster_cols = F,color = colorpanel(101,"white","firebrick"),angle_col = 0,display_numbers = T,number_color = "black",cluster_rows = F)
#dev.off()



immunonormdata <- MotrpacRatTraining6moData::IMMUNO_NORM_DATA_FLAT
plasmaimmunonormdata <- immunonormdata[immunonormdata$tissue %in% "PLASMA",]
rownames(plasmaimmunonormdata) <- paste(plasmaimmunonormdata$feature_ID,"_",plasmaimmunonormdata$dataset,sep = "")
plasmaimmunonormdata <- plasmaimmunonormdata[,intersect(colnames(lungrnanorm),colnames(immunonormdata))]

lungimmunonormdata <- immunonormdata[immunonormdata$tissue %in% "LUNG",]
rownames(lungimmunonormdata) <- paste(lungimmunonormdata$feature_ID,"_",lungimmunonormdata$dataset,sep = "")
lungimmunonormdata <- lungimmunonormdata[,intersect(colnames(lungrnanorm),colnames(immunonormdata))]

gastroimmunonormdata <- immunonormdata[immunonormdata$tissue %in% "SKM-GN",]
rownames(gastroimmunonormdata) <- paste(gastroimmunonormdata$feature_ID,"_",gastroimmunonormdata$dataset,sep = "")
gastroimmunonormdata <- gastroimmunonormdata[,intersect(colnames(lungrnanorm),colnames(immunonormdata))]


immunonormdatameta <- immunonormdata[,c("feature","feature_ID","tissue","assay","dataset")]
immunonormdatameta$feature_ID_dataset <- paste(immunonormdatameta$feature_ID,"_",immunonormdatameta$dataset,sep = "")

rownames(immunonormdatameta) <- paste(immunonormdatameta$feature_ID,"_",immunonormdatameta$dataset,"_",immunonormdatameta$tissue,sep = "")

immunonormdatatrim <- immunonormdata[,intersect(colnames(lungrnanorm),colnames(immunonormdata))]
immunonormdatameta$PU1LungRNAcor <- 0
immunonormdatameta$IRF8LungProcor <- 0

pu1rnaid <- tfanno["PU.1(ETS)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer","Ensembl"]
irf8proid <- tfproanno["IRF8(IRF)/BMDM-IRF8-ChIP-Seq(GSE77884)/Homer","Lung.Pro.ID"]

for(i in 1:dim(immunonormdatatrim)[1]){

  immunonormdatameta[i,"PU1LungRNAcor"] <- cor(t(immunonormdatatrim[i,colnames(immunonormdatatrim)[(!is.na(immunonormdatatrim[i,]))]]),t(lungrnanorm[pu1rnaid,colnames(immunonormdatatrim)[(!is.na(immunonormdatatrim[i,]))]]))
  immunonormdatameta[i,"IRF8LungProcor"] <- cor(t(immunonormdatatrim[i,colnames(immunonormdatatrim)[(!is.na(immunonormdatatrim[i,]))]]),t(lungproprnorm[irf8proid,colnames(immunonormdatatrim)[(!is.na(immunonormdatatrim[i,]))]]))
  
}

save.image("Figure5_021325.RData")

crosstissueIFNGl2fc <- rbind(gastroimmunel2fc["IFNY_rat-mag27plex",],
                         heartimmunel2fc["IFNY_rat-mag27plex",],
                         hippoimmunel2fc["IFNY_rat-mag27plex",],
                         kidneyimmunel2fc["IFNY_rat-mag27plex",],
                         liverimmunel2fc["IFNY_rat-mag27plex",],
                         lungimmunel2fc["IFNY_rat-mag27plex",],
                         plasmaimmunel2fc["IFNY_rat-mag27plex",],
                         brownimmunel2fc["IFNY_rat-mag27plex",],
                         whiteimmunel2fc["IFNY_rat-mag27plex",])

rownames(crosstissueIFNGl2fc) <- c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","PLASMA","BAT","WAT-SC")

#png("IFNG L2FC across tissues.png",width = 5,height = 5,units = "in",res = 600)
#pheatmap(crosstissueIFNGl2fc,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_rows = F,cluster_cols = F,display_numbers = T,number_color = "black",angle_col = 0,main = "IFNG L2FC Across Tissues")
#dev.off()

crosstissueIFNGpval <- rbind(gastroimmunepval["IFNY_rat-mag27plex",],
                             heartimmunepval["IFNY_rat-mag27plex",],
                             hippoimmunepval["IFNY_rat-mag27plex",],
                             kidneyimmunepval["IFNY_rat-mag27plex",],
                             liverimmunepval["IFNY_rat-mag27plex",],
                             lungimmunepval["IFNY_rat-mag27plex",],
                             plasmaimmunepval["IFNY_rat-mag27plex",],
                             brownimmunepval["IFNY_rat-mag27plex",],
                             whiteimmunepval["IFNY_rat-mag27plex",])

rownames(crosstissueIFNGpval) <- c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","PLASMA","BAT","WAT-SC")

#png("IFNG pval across tissues.png",width = 5,height = 5,units = "in",res = 600)
#pheatmap(-log10(crosstissueIFNGpval),color = colorpanel(101,"white","firebrick"),cluster_rows = F,cluster_cols = F,display_numbers = T,number_color = "black",angle_col = 0,main = "IFNG -log10(pval) Across Tissues")
#dev.off()


plasmasigmarkerenssymbols <- c("ENSRNOG00000001821",
                               "ENSRNOG00000012417",
                               "ENSRNOG00000047466",
                               "ENSRNOG00000047466",
                               "ENSRNOG00000012052",
                               "ENSRNOG00000053979",
                               "ENSRNOG00000007335",
                               "ENSRNOG00000010349",
                               "ENSRNOG00000005498",
                               "ENSRNOG00000005498",
                               "ENSRNOG00000007468",
                               "ENSRNOG00000004647",
                               "ENSRNOG00000007652",
                               "ENSRNOG00000012467",
                               "ENSRNOG00000008111",
                               "ENSRNOG00000022256",
                               "ENSRNOG00000045797",
                               "ENSRNOG00000045797",
                               "ENSRNOG00000047040",
                               "ENSRNOG00000002843",
                               "ENSRNOG00000011205",
                               "ENSRNOG00000017374",
                               "ENSRNOG00000030462",
                               "ENSRNOG00000020874",
                               "ENSRNOG00000020877",
                               "ENSRNOG00000010906",
                               "ENSRNOG00000012840",
                               "ENSRNOG00000055156")


plasmasigmarkerdf <- data.frame(row.names = rownames(plasmaimmunepval[apply(plasmaimmunepval,1,min) < 0.05,]),
                                "Gene.Ens" = plasmasigmarkerenssymbols)
plasmasigmarkerdf$plasmal2fc <- plasmaimmunel2fc[rownames(plasmasigmarkerdf),]

plasmaimmunonormdatameta <- immunonormdatameta[immunonormdatameta$tissue %in% "PLASMA",]
rownames(plasmaimmunonormdatameta) <- plasmaimmunonormdatameta$feature_ID_dataset

plasmasigmarkerdf$PU1LungRNAcor <- plasmaimmunonormdatameta[rownames(plasmasigmarkerdf),"PU1LungRNAcor"]
plasmasigmarkerdf$IRF8LungProcor <- plasmaimmunonormdatameta[rownames(plasmasigmarkerdf),"IRF8LungProcor"]

bloodrnanorm <- MotrpacRatTraining6moData::TRNSCRPT_BLOOD_NORM_DATA
rownames(bloodrnanorm) <- bloodrnanorm$feature_ID
bloodrnanorm <- bloodrnanorm[,c(5:dim(bloodrnanorm)[2])]

for(i in 1:dim(bloodrnanorm)[2]){
  ourid <- colnames(bloodrnanorm)[i]
  colnames(bloodrnanorm)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0001"]
}

for(i in 1:dim(hippornanorm)[2]){
  ourid <- colnames(hippornanorm)[i]
  colnames(hippornanorm)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0001"]
}

for(i in 1:dim(brownrnanorm)[2]){
  ourid <- colnames(brownrnanorm)[i]
  colnames(brownrnanorm)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0001"]
}


plasmasigmarkerdf$SKMGNRNAcor <- 0
plasmasigmarkerdf$HEARTRNAcor <- 0
plasmasigmarkerdf$HIPPOCRNAcor <- 0
plasmasigmarkerdf$KIDNEYRNAcor <- 0
plasmasigmarkerdf$LIVERRNAcor <- 0
plasmasigmarkerdf$LUNGRNAcor <- 0
plasmasigmarkerdf$BLOODRNAcor <- 0
plasmasigmarkerdf$BATRNAcor <- 0
plasmasigmarkerdf$WATSCRNAcor <- 0

for(i in 1:dim(plasmasigmarkerdf)[1]){
  ourmarker <- plasmasigmarkerdf$Gene.Ens[i]
  if(ourmarker %in% rownames(gastrornanorm)){
    ourrnanorm <- t(gastrornanorm[ourmarker,intersect(colnames(plasmaimmunonormdata),colnames(gastrornanorm))])
    ourplasmanorm <- t(plasmaimmunonormdata[rownames(plasmasigmarkerdf)[i],intersect(colnames(plasmaimmunonormdata),colnames(gastrornanorm))])
    plasmasigmarkerdf[i,"SKMGNRNAcor"] <- cor(ourplasmanorm[!is.na(ourplasmanorm),],ourrnanorm[!is.na(ourplasmanorm),])
  }
  
  if(ourmarker %in% rownames(heartrnanorm)){
    ourrnanorm <- t(heartrnanorm[ourmarker,intersect(colnames(plasmaimmunonormdata),colnames(heartrnanorm))])
    ourplasmanorm <- t(plasmaimmunonormdata[rownames(plasmasigmarkerdf)[i],intersect(colnames(plasmaimmunonormdata),colnames(heartrnanorm))])
    plasmasigmarkerdf[i,"HEARTRNAcor"] <- cor(ourplasmanorm[!is.na(ourplasmanorm),],ourrnanorm[!is.na(ourplasmanorm),])
  }
  
  if(ourmarker %in% rownames(hippornanorm)){
    ourrnanorm <- t(hippornanorm[ourmarker,intersect(colnames(plasmaimmunonormdata),colnames(hippornanorm))])
    ourplasmanorm <- t(plasmaimmunonormdata[rownames(plasmasigmarkerdf)[i],intersect(colnames(plasmaimmunonormdata),colnames(hippornanorm))])
    plasmasigmarkerdf[i,"HIPPOCRNAcor"] <- cor(ourplasmanorm[!is.na(ourplasmanorm),],ourrnanorm[!is.na(ourplasmanorm),])
  }
  
  if(ourmarker %in% rownames(kidneyrnanorm)){
    ourrnanorm <- t(kidneyrnanorm[ourmarker,intersect(colnames(plasmaimmunonormdata),colnames(kidneyrnanorm))])
    ourplasmanorm <- t(plasmaimmunonormdata[rownames(plasmasigmarkerdf)[i],intersect(colnames(plasmaimmunonormdata),colnames(kidneyrnanorm))])
    plasmasigmarkerdf[i,"KIDNEYRNAcor"] <- cor(ourplasmanorm[!is.na(ourplasmanorm),],ourrnanorm[!is.na(ourplasmanorm),])
  }
  
  if(ourmarker %in% rownames(liverrnanorm)){
    ourrnanorm <- t(liverrnanorm[ourmarker,intersect(colnames(plasmaimmunonormdata),colnames(liverrnanorm))])
    ourplasmanorm <- t(plasmaimmunonormdata[rownames(plasmasigmarkerdf)[i],intersect(colnames(plasmaimmunonormdata),colnames(liverrnanorm))])
    plasmasigmarkerdf[i,"LIVERRNAcor"] <- cor(ourplasmanorm[!is.na(ourplasmanorm),],ourrnanorm[!is.na(ourplasmanorm),])
  }
  
  if(ourmarker %in% rownames(lungrnanorm)){
    ourrnanorm <- t(lungrnanorm[ourmarker,intersect(colnames(plasmaimmunonormdata),colnames(lungrnanorm))])
    ourplasmanorm <- t(plasmaimmunonormdata[rownames(plasmasigmarkerdf)[i],intersect(colnames(plasmaimmunonormdata),colnames(lungrnanorm))])
    plasmasigmarkerdf[i,"LUNGRNAcor"] <- cor(ourplasmanorm[!is.na(ourplasmanorm),],ourrnanorm[!is.na(ourplasmanorm),])
  }
  
  if(ourmarker %in% rownames(bloodrnanorm)){
    ourrnanorm <- t(bloodrnanorm[ourmarker,intersect(colnames(plasmaimmunonormdata),colnames(bloodrnanorm))])
    ourplasmanorm <- t(plasmaimmunonormdata[rownames(plasmasigmarkerdf)[i],intersect(colnames(plasmaimmunonormdata),colnames(bloodrnanorm))])
    plasmasigmarkerdf[i,"BLOODRNAcor"] <- cor(ourplasmanorm[!is.na(ourplasmanorm),],ourrnanorm[!is.na(ourplasmanorm),])
  }
  
  if(ourmarker %in% rownames(brownrnanorm)){
    ourrnanorm <- t(brownrnanorm[ourmarker,intersect(colnames(plasmaimmunonormdata),colnames(brownrnanorm))])
    ourplasmanorm <- t(plasmaimmunonormdata[rownames(plasmasigmarkerdf)[i],intersect(colnames(plasmaimmunonormdata),colnames(brownrnanorm))])
    plasmasigmarkerdf[i,"BATRNAcor"] <- cor(ourplasmanorm[!is.na(ourplasmanorm),],ourrnanorm[!is.na(ourplasmanorm),])
  }
  
  if(ourmarker %in% rownames(whiternanorm)){
    ourrnanorm <- t(whiternanorm[ourmarker,intersect(colnames(plasmaimmunonormdata),colnames(whiternanorm))])
    ourplasmanorm <- t(plasmaimmunonormdata[rownames(plasmasigmarkerdf)[i],intersect(colnames(plasmaimmunonormdata),colnames(whiternanorm))])
    plasmasigmarkerdf[i,"WATSCRNAcor"] <- cor(ourplasmanorm[!is.na(ourplasmanorm),],ourrnanorm[!is.na(ourplasmanorm),])
  }
}


#png("PlasmaSigMarkersCorrelatedwithTFsandTheirRNALevels.png",width = 8.5,height = 6,units = "in",res = 600)
#pheatmap(plasmasigmarkerdf[,c("PU1LungRNAcor","IRF8LungProcor","SKMGNRNAcor","HEARTRNAcor","HIPPOCRNAcor","KIDNEYRNAcor","LIVERRNAcor","LUNGRNAcor","BLOODRNAcor","BATRNAcor","WATSCRNAcor")],cluster_rows = F,cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = 315,display_numbers = T,number_color = "black")
#dev.off()


#####
# Let's investigate cell type deconvolution in Lung related to IRF8
####

#load("New_CellTypeComp_Figure1E_1F_S6_6624.RData")
load("CellTypeComp_Figure1E_1F_S10_S11_112624.RData")

for(i in 1:dim(lungrnanorm)[2]){
  ourid <- colnames(lungrnanorm)[i]
  colnames(lungrnanorm)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0001"]
}

for(i in 1:dim(lungcomboSPVs)[1]){
  ourid <- rownames(lungcomboSPVs)[i]
  rownames(lungcomboSPVs)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0001"]
}

ourirf8df <- data.frame(row.names = intersect(colnames(lungproprnorm),rownames(lungcomboSPVs)),
                    "IRF8 Protein" = t(lungproprnorm[tfproanno["IRF8(IRF)/BMDM-IRF8-ChIP-Seq(GSE77884)/Homer","Lung.Pro.ID"],intersect(colnames(lungproprnorm),rownames(lungcomboSPVs))]),
                    "Granulocyte" = lungcomboSPVs[intersect(colnames(lungproprnorm),rownames(lungcomboSPVs)),"Granulocyte"],
                    "Megakaryocyte" = lungcomboSPVs[intersect(colnames(lungproprnorm),rownames(lungcomboSPVs)),"Megakaryocyte"],
                    "Erythrocyte" = lungcomboSPVs[intersect(colnames(lungproprnorm),rownames(lungcomboSPVs)),"Erythrocyte"],
                    "Monocyte" = lungcomboSPVs[intersect(colnames(lungproprnorm),rownames(lungcomboSPVs)),"Monocyte"],
                    "Dendritic.Cell" = lungcomboSPVs[intersect(colnames(lungproprnorm),rownames(lungcomboSPVs)),"Dendritic.Cell"],
                    "CD4.T.Cell" = lungcomboSPVs[intersect(colnames(lungproprnorm),rownames(lungcomboSPVs)),"CD4.T.Cell"],
                    "GMP.Cell" = lungcomboSPVs[intersect(colnames(lungproprnorm),rownames(lungcomboSPVs)),"GMP.Cell"],
                    "CD56pCD16p.NK.Cell" = lungcomboSPVs[intersect(colnames(lungproprnorm),rownames(lungcomboSPVs)),"CD56pCD16p.NK.Cell"],
                    "CD56nCD16n.NK.Cell" = lungcomboSPVs[intersect(colnames(lungproprnorm),rownames(lungcomboSPVs)),"CD56nCD16n.NK.Cell"],
                    "B.Cell" = lungcomboSPVs[intersect(colnames(lungproprnorm),rownames(lungcomboSPVs)),"B.Cell"],
                    "CD8.T.Cell" = lungcomboSPVs[intersect(colnames(lungproprnorm),rownames(lungcomboSPVs)),"CD8.T.Cell"],
                    "Endothelial.Cells" = lungcomboSPVs[intersect(colnames(lungproprnorm),rownames(lungcomboSPVs)),"Endothelial.Cells"],
                    "PCV.Endothelial.Cells" = lungcomboSPVs[intersect(colnames(lungproprnorm),rownames(lungcomboSPVs)),"PCV.Endothelial.Cells"],
                    "FAP.Cells" = lungcomboSPVs[intersect(colnames(lungproprnorm),rownames(lungcomboSPVs)),"FAP.Cells"],
                    "Pericytes" = lungcomboSPVs[intersect(colnames(lungproprnorm),rownames(lungcomboSPVs)),"Pericytes"])
ourirf8df$Sex <- "Female"
ourirf8df$Group <- "control"
for(i in 1:dim(ourirf8df)[1]){
  ourpid <- rownames(ourirf8df)[i]
  ourirf8df[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourpid,"pass1bf0027"][1]]
  ourirf8df[i,"Group"] <- pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourpid,"pass1bf0011"][1]
}
ourirf8df[ourirf8df$Group %in% "One-week program","Group"] <- "1w"
ourirf8df[ourirf8df$Group %in% "Two-week program","Group"] <- "2w"
ourirf8df[ourirf8df$Group %in% "Four-week program","Group"] <- "4w"
ourirf8df[ourirf8df$Group %in% "Eight-week program Training Group","Group"] <- "8w"
ourirf8df[ourirf8df$Group %in% "Eight-week program Control Group","Group"] <- "control"
colnames(ourirf8df)[1] <- "IRF8"

#png("IRF8 Protein vs Granulocyte Subject Plot_61124.png",width = 5,height = 5,units = "in",res = 600)
#ggscatter(ourirf8df, x = "IRF8", y = "Granulocyte",color = "Group",shape = "Sex",size = 4) + scale_color_manual(values = ann_cols$Group) + ggtitle(paste("LUNG: IRF8 vs Granulocyte Correlation: ",substr(toString(cor(ourirf8df$IRF8,ourirf8df$Granulocyte)),1,5),sep = "")) + xlab("IRF8 Protein") + ylab("Granulocyte Cell Change") + theme(legend.box = "vertical",legend.margin = margin())
#dev.off()

png("Figure 5F_021325.png",width = 5,height = 5,units = "in",res = 600)
ggscatter(ourirf8df, x = "IRF8", y = "Monocyte",color = "Group",shape = "Sex",size = 4) + scale_color_manual(values = ann_cols$Group) + ggtitle(paste("LUNG: IRF8 vs Monocyte Correlation: ",substr(toString(cor(ourirf8df$IRF8,ourirf8df$Monocyte)),1,5),sep = "")) + xlab("IRF8 Protein") + ylab("Monocyte Cell Change") + theme(legend.box = "vertical",legend.margin = margin())
dev.off()

pdf("Figure 5F_021325.pdf",width = 5,height = 5)
ggscatter(ourirf8df, x = "IRF8", y = "Monocyte",color = "Group",shape = "Sex",size = 4) + scale_color_manual(values = ann_cols$Group) + ggtitle(paste("LUNG: IRF8 vs Monocyte Correlation: ",substr(toString(cor(ourirf8df$IRF8,ourirf8df$Monocyte)),1,5),sep = "")) + xlab("IRF8 Protein") + ylab("Monocyte Cell Change") + theme(legend.box = "vertical",legend.margin = margin())
dev.off()

#png("IRF8 Protein vs Megakaryocyte Subject Plot_61124.png",width = 5,height = 5,units = "in",res = 600)
#ggscatter(ourirf8df, x = "IRF8", y = "Megakaryocyte",color = "Group",shape = "Sex",size = 4) + scale_color_manual(values = ann_cols$Group) + ggtitle(paste("LUNG: IRF8 vs Megakaryocyte Correlation: ",substr(toString(cor(ourirf8df$IRF8,ourirf8df$Megakaryocyte)),1,5),sep = "")) + xlab("IRF8 Protein") + ylab("Megakaryocyte Cell Change") + theme(legend.box = "vertical",legend.margin = margin())
#dev.off()

#png("IRF8 Protein vs Erythrocyte Subject Plot_61124.png",width = 5,height = 5,units = "in",res = 600)
#ggscatter(ourirf8df, x = "IRF8", y = "Erythrocyte",color = "Group",shape = "Sex",size = 4) + scale_color_manual(values = ann_cols$Group) + ggtitle(paste("LUNG: IRF8 vs Erythrocyte Correlation: ",substr(toString(cor(ourirf8df$IRF8,ourirf8df$Erythrocyte)),1,5),sep = "")) + xlab("IRF8 Protein") + ylab("Erythrocyte Cell Change") + theme(legend.box = "vertical",legend.margin = margin())
#dev.off()

#png("IRF8 Protein vs Dendritic Cell Subject Plot_61124.png",width = 5,height = 5,units = "in",res = 600)
#ggscatter(ourirf8df, x = "IRF8", y = "Dendritic.Cell",color = "Group",shape = "Sex",size = 4) + scale_color_manual(values = ann_cols$Group) + ggtitle(paste("LUNG: IRF8 vs Dendritic.Cell Correlation: ",substr(toString(cor(ourirf8df$IRF8,ourirf8df$Dendritic.Cell)),1,5),sep = "")) + xlab("IRF8 Protein") + ylab("Dendritic.Cell Cell Change") + theme(legend.box = "vertical",legend.margin = margin())
#dev.off()

#png("IRF8 Protein vs CD4 T Cell Subject Plot_61124.png",width = 5,height = 5,units = "in",res = 600)
#ggscatter(ourirf8df, x = "IRF8", y = "CD4.T.Cell",color = "Group",shape = "Sex",size = 4) + scale_color_manual(values = ann_cols$Group) + ggtitle(paste("LUNG: IRF8 vs CD4.T.Cell Correlation: ",substr(toString(cor(ourirf8df$IRF8,ourirf8df$CD4.T.Cell)),1,5),sep = "")) + xlab("IRF8 Protein") + ylab("CD4.T.Cell Cell Change") + theme(legend.box = "vertical",legend.margin = margin())
#dev.off()

#png("IRF8 Protein vs GMP Cell Subject Plot_61124.png",width = 5,height = 5,units = "in",res = 600)
#ggscatter(ourirf8df, x = "IRF8", y = "GMP.Cell",color = "Group",shape = "Sex",size = 4) + scale_color_manual(values = ann_cols$Group) + ggtitle(paste("LUNG: IRF8 vs GMP.Cell Correlation: ",substr(toString(cor(ourirf8df$IRF8,ourirf8df$GMP.Cell)),1,5),sep = "")) + xlab("IRF8 Protein") + ylab("GMP.Cell Cell Change") + theme(legend.box = "vertical",legend.margin = margin())
#dev.off()

#png("IRF8 Protein vs CD56pCD16p NK Cell Subject Plot_61124.png",width = 5,height = 5,units = "in",res = 600)
#ggscatter(ourirf8df, x = "IRF8", y = "CD56pCD16p.NK.Cell",color = "Group",shape = "Sex",size = 4) + scale_color_manual(values = ann_cols$Group) + ggtitle(paste("LUNG: IRF8 vs CD56pCD16p.NK.Cell Correlation: ",substr(toString(cor(ourirf8df$IRF8,ourirf8df$CD56pCD16p.NK.Cell)),1,5),sep = "")) + xlab("IRF8 Protein") + ylab("CD56pCD16p.NK.Cell Cell Change") + theme(legend.box = "vertical",legend.margin = margin())
#dev.off()

#png("IRF8 Protein vs CD56nCD16n NK Cell Subject Plot_61124.png",width = 5,height = 5,units = "in",res = 600)
#ggscatter(ourirf8df, x = "IRF8", y = "CD56nCD16n.NK.Cell",color = "Group",shape = "Sex",size = 4) + scale_color_manual(values = ann_cols$Group) + ggtitle(paste("LUNG: IRF8 vs CD56nCD16n.NK.Cell Correlation: ",substr(toString(cor(ourirf8df$IRF8,ourirf8df$CD56nCD16n.NK.Cell)),1,5),sep = "")) + xlab("IRF8 Protein") + ylab("CD56nCD16n.NK.Cell Cell Change") + theme(legend.box = "vertical",legend.margin = margin())
#dev.off()

#png("IRF8 Protein vs B Cell Subject Plot_61124.png",width = 5,height = 5,units = "in",res = 600)
#ggscatter(ourirf8df, x = "IRF8", y = "B.Cell",color = "Group",shape = "Sex",size = 4) + scale_color_manual(values = ann_cols$Group) + ggtitle(paste("LUNG: IRF8 vs B.Cell Correlation: ",substr(toString(cor(ourirf8df$IRF8,ourirf8df$B.Cell)),1,5),sep = "")) + xlab("IRF8 Protein") + ylab("B.Cell Cell Change") + theme(legend.box = "vertical",legend.margin = margin())
#dev.off()

#png("IRF8 Protein vs CD8 T Cell Subject Plot_61124.png",width = 5,height = 5,units = "in",res = 600)
#ggscatter(ourirf8df, x = "IRF8", y = "CD8.T.Cell",color = "Group",shape = "Sex",size = 4) + scale_color_manual(values = ann_cols$Group) + ggtitle(paste("LUNG: IRF8 vs CD8.T.Cell Correlation: ",substr(toString(cor(ourirf8df$IRF8,ourirf8df$CD8.T.Cell)),1,5),sep = "")) + xlab("IRF8 Protein") + ylab("CD8.T.Cell Cell Change") + theme(legend.box = "vertical",legend.margin = margin())
#dev.off()

#png("IRF8 Protein vs Endothelial Cells Subject Plot_61124.png",width = 5,height = 5,units = "in",res = 600)
#ggscatter(ourirf8df, x = "IRF8", y = "Endothelial.Cells",color = "Group",shape = "Sex",size = 4) + scale_color_manual(values = ann_cols$Group) + ggtitle(paste("LUNG: IRF8 vs Endothelial.Cells Correlation: ",substr(toString(cor(ourirf8df$IRF8,ourirf8df$Endothelial.Cells)),1,5),sep = "")) + xlab("IRF8 Protein") + ylab("Endothelial.Cells Cell Change") + theme(legend.box = "vertical",legend.margin = margin())
#dev.off()

#png("IRF8 Protein vs PCV Endothelial Cells Subject Plot_61124.png",width = 5,height = 5,units = "in",res = 600)
#ggscatter(ourirf8df, x = "IRF8", y = "PCV.Endothelial.Cells",color = "Group",shape = "Sex",size = 4) + scale_color_manual(values = ann_cols$Group) + ggtitle(paste("LUNG: IRF8 vs PCV.Endothelial.Cells Correlation: ",substr(toString(cor(ourirf8df$IRF8,ourirf8df$PCV.Endothelial.Cells)),1,5),sep = "")) + xlab("IRF8 Protein") + ylab("PCV.Endothelial.Cells Cell Change") + theme(legend.box = "vertical",legend.margin = margin())
#dev.off()

#png("IRF8 Protein vs FAP Cells Subject Plot_61124.png",width = 5,height = 5,units = "in",res = 600)
#ggscatter(ourirf8df, x = "IRF8", y = "FAP.Cells",color = "Group",shape = "Sex",size = 4) + scale_color_manual(values = ann_cols$Group) + ggtitle(paste("LUNG: IRF8 vs FAP.Cells Correlation: ",substr(toString(cor(ourirf8df$IRF8,ourirf8df$FAP.Cells)),1,5),sep = "")) + xlab("IRF8 Protein") + ylab("FAP.Cells Cell Change") + theme(legend.box = "vertical",legend.margin = margin())
#dev.off()

#png("IRF8 Protein vs Pericytes Subject Plot_61124.png",width = 5,height = 5,units = "in",res = 600)
#ggscatter(ourirf8df, x = "IRF8", y = "Pericytes",color = "Group",shape = "Sex",size = 4) + scale_color_manual(values = ann_cols$Group) + ggtitle(paste("LUNG: IRF8 vs Pericytes Correlation: ",substr(toString(cor(ourirf8df$IRF8,ourirf8df$Pericytes)),1,5),sep = "")) + xlab("IRF8 Protein") + ylab("Pericytes Cell Change") + theme(legend.box = "vertical",legend.margin = margin())
#dev.off()


ourpu1df <- data.frame(row.names = intersect(colnames(lungrnanorm),rownames(lungcomboSPVs)),
                        "PU1 RNA" = t(lungrnanorm[tfanno["PU.1(ETS)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer","Ensembl"],intersect(colnames(lungrnanorm),rownames(lungcomboSPVs))]),
                        "Granulocyte" = lungcomboSPVs[intersect(colnames(lungrnanorm),rownames(lungcomboSPVs)),"Granulocyte"])
ourpu1df$Sex <- "Female"
ourpu1df$Group <- "control"
for(i in 1:dim(ourpu1df)[1]){
  ourpid <- rownames(ourpu1df)[i]
  ourpu1df[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourpid,"pass1bf0027"][1]]
  ourpu1df[i,"Group"] <- pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourpid,"pass1bf0011"][1]
}
ourpu1df[ourpu1df$Group %in% "One-week program","Group"] <- "1w"
ourpu1df[ourpu1df$Group %in% "Two-week program","Group"] <- "2w"
ourpu1df[ourpu1df$Group %in% "Four-week program","Group"] <- "4w"
ourpu1df[ourpu1df$Group %in% "Eight-week program Training Group","Group"] <- "8w"
ourpu1df[ourpu1df$Group %in% "Eight-week program Control Group","Group"] <- "control"
colnames(ourpu1df)[1] <- "PU1"

#png("PU1 RNA vs Granulocyte Subject Plot.png",width = 5,height = 5,units = "in",res = 600)
#ggscatter(ourpu1df, x = "PU1", y = "Granulocyte",color = "Group",shape = "Sex",size = 4) + scale_color_manual(values = ann_cols$Group) + ggtitle(paste("LUNG: PU1 vs Granulocyte Correlation: ",substr(toString(cor(ourpu1df$PU1,ourpu1df$Granulocyte)),1,5),sep = "")) + xlab("PU1 RNA") + ylab("Granulocyte Cell Change") + theme(legend.box = "vertical",legend.margin = margin())
#dev.off()

save.image("Figure5_021325.RData")

#####
# Let's look at IRF8 vs CBFB and RUNX1
####

runx1ens <- enstosym[enstosym$Symbol %in% "Runx1","Ensembl"][1]
cbfbens <- enstosym[enstosym$Symbol %in% "Cbfb","Ensembl"][1]


ourdf <- data.frame(row.names = colnames(lungproprnorm),
                    "IRF8" = t(lungproprnorm["NP_001008722.1",]),
                    "CBFB" = t(lungproprnorm["NP_001013209.1",]),
                    "RUNX1" = t(lungproprnorm["NP_059021.1",]),
                    "PU1_Protein" = t(lungproprnorm["NP_001005892.1",]),
                    "PU1_RNA" = t(lungrnanorm[tfanno["PU.1(ETS)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer","Ensembl"],colnames(lungproprnorm)]),
                    "CBFB_RNA" = t(lungrnanorm[cbfbens,colnames(lungproprnorm)]),
                    "RUNX1_RNA" = t(lungrnanorm[runx1ens,colnames(lungproprnorm)]))
ourdf$Sex <- "Female"
ourdf$Group <- "control"
for(i in 1:dim(ourdf)[1]){
  ourpid <- rownames(ourdf)[i]
  ourdf[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourpid,"pass1bf0027"][1]]
  ourdf[i,"Group"] <- pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourpid,"pass1bf0011"][1]
}
ourdf[ourdf$Group %in% "One-week program","Group"] <- "1w"
ourdf[ourdf$Group %in% "Two-week program","Group"] <- "2w"
ourdf[ourdf$Group %in% "Four-week program","Group"] <- "4w"
ourdf[ourdf$Group %in% "Eight-week program Training Group","Group"] <- "8w"
ourdf[ourdf$Group %in% "Eight-week program Control Group","Group"] <- "control"
colnames(ourdf)[1:7] <- c("IRF8","CBFB","RUNX1","PU.1 Protein","PU.1 RNA","CBFB RNA","RUNX1 RNA")

png("Figure 5G_021325.png",width = 5,height = 5,units = "in",res = 600)
ggscatter(ourdf, x = "IRF8", y = "CBFB",color = "Group",shape = "Sex",size = 4) + scale_color_manual(values = ann_cols$Group) + ggtitle(paste("LUNG: IRF8 vs CBFB Correlation: ",substr(toString(cor(ourdf$IRF8,ourdf$CBFB)),1,5),sep = "")) + xlab("IRF8 Protein") + ylab("CBFB Protein") + theme(legend.box = "vertical",legend.margin = margin())
dev.off()

pdf("Figure 5G_021325.pdf",width = 5,height = 5)
ggscatter(ourdf, x = "IRF8", y = "CBFB",color = "Group",shape = "Sex",size = 4) + scale_color_manual(values = ann_cols$Group) + ggtitle(paste("LUNG: IRF8 vs CBFB Correlation: ",substr(toString(cor(ourdf$IRF8,ourdf$CBFB)),1,5),sep = "")) + xlab("IRF8 Protein") + ylab("CBFB Protein") + theme(legend.box = "vertical",legend.margin = margin())
dev.off()

png("Figure 5H_021325.png",width = 5,height = 5,units = "in",res = 600)
ggscatter(ourdf, x = "IRF8", y = "RUNX1",color = "Group",shape = "Sex",size = 4) + scale_color_manual(values = ann_cols$Group) + ggtitle(paste("LUNG: IRF8 vs RUNX1 Correlation: ",substr(toString(cor(ourdf$IRF8,ourdf$RUNX1)),1,5),sep = "")) + xlab("IRF8 Protein") + ylab("RUNX1 Protein") + theme(legend.box = "vertical",legend.margin = margin())
dev.off()

pdf("Figure 5H_021325.pdf",width = 5,height = 5)
ggscatter(ourdf, x = "IRF8", y = "RUNX1",color = "Group",shape = "Sex",size = 4) + scale_color_manual(values = ann_cols$Group) + ggtitle(paste("LUNG: IRF8 vs RUNX1 Correlation: ",substr(toString(cor(ourdf$IRF8,ourdf$RUNX1)),1,5),sep = "")) + xlab("IRF8 Protein") + ylab("RUNX1 Protein") + theme(legend.box = "vertical",legend.margin = margin())
dev.off()

png("Figure 5I_021325.png",width = 5,height = 5,units = "in",res = 600)
ggscatter(ourdf, x = "RUNX1", y = "PU.1 RNA",color = "Group",shape = "Sex",size = 4) + scale_color_manual(values = ann_cols$Group) + ggtitle(paste("LUNG: RUNX1 vs PU1 Correlation: ",substr(toString(cor(ourdf$RUNX1,ourdf$`PU.1 RNA`)),1,5),sep = "")) + xlab("RUNX1 Protein") + ylab("PU.1 RNA") + theme(legend.box = "vertical",legend.margin = margin())
dev.off()

pdf("Figure 5I_021325.pdf",width = 5,height = 5)
ggscatter(ourdf, x = "RUNX1", y = "PU.1 RNA",color = "Group",shape = "Sex",size = 4) + scale_color_manual(values = ann_cols$Group) + ggtitle(paste("LUNG: RUNX1 vs PU1 Correlation: ",substr(toString(cor(ourdf$RUNX1,ourdf$`PU.1 RNA`)),1,5),sep = "")) + xlab("RUNX1 Protein") + ylab("PU.1 RNA") + theme(legend.box = "vertical",legend.margin = margin())
dev.off()

save.image("Figure5_S20_to_S26_021325.RData")