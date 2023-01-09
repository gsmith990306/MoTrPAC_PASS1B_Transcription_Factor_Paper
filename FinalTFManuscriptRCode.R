# This file contains R code to generate the figures contained within 
# "Multiomic identification of key transcriptional regulatory programs during endurance exercise training"
# Authors: Gregory R Smith, Bingqing Zhao, Malene E. Lindholm, Archana Raja
# Table of Contents
# 1. Library loading:
#    Lines 32-63 
# 2. Loading data files
#    Lines 66-222
# 3. Generating L2FC matrices and significant gene sets
#    Lines 230-642
# 4. Generating color annotation sets
#    Lines 649-780
# 5. Generating DEGaP sets
#    Lines 786-1693
# 6. Load homer results of TF enrichment from DEGaP sets
#    Lines 1703-1806
# 7. Figure 1 and Supplemental Figures S1-S5
#    Lines 1814-3205
# 8. Figure 2 and Supplemental Figure S6
#    Lines 3212-3949
# 9. Figure 3 and Supplemental Figure S7
#    Lines 3955-4809
# 10. Figure 4 and Supplemental Figure S8
#    Lines 4818-5697
# 11. Figure 5 and Supplemental Figures S9-S13
#    Lines 5704-7817
# 12. Figure 6 and Supplemental Figures S14-S15
#    Lines 7824-8586
# 13. Figure 7 and Supplemental Figure S16
#    Lines 8593-9548

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
library(ComplexHeatmap)
library(GenomicRanges)

######################
## Loading initial important data files - file names used below have been simplified from MoTrPAC data releases

# To be edited to your directory containing the folder with all necessary paper data files
setwd("~/Mount Sinai/Sealfon Laboratory/MoTrPAC/PASS1B Raw Data")

pass1bphenodata <- read.delim(file = "PASS1B Transcription Factor Paper Data/Phenotype Data/motrpac_pass1b-06_pheno_viallabel-data.txt",header = T,sep = "\t")
phenomeasuredata <- read.delim(file = "PASS1B Transcription Factor Paper Data/Phenotype Data/pass1b-06_phenotype_motrpac_pass1b-06_pheno_viallabel_data_merged_v2.0.txt",header = T,sep = "\t",row.names = 4)
# To correct for a missing entry
phenomeasuredata["90412017702",c("calculated_variables___wgt_gain_after_train","calculated_variables___pct_body_fat_change","calculated_variables___pct_body_lean_change","calculated_variables___pct_body_fluid_change","calculated_variables___lactate_change_dueto_train","calculated_variables___vo_2_max_change")] <- phenomeasuredata["90412017014",c("calculated_variables___wgt_gain_after_train","calculated_variables___pct_body_fat_change","calculated_variables___pct_body_lean_change","calculated_variables___pct_body_fluid_change","calculated_variables___lactate_change_dueto_train","calculated_variables___vo_2_max_change")]

load(file = "PASS1B Transcription Factor Paper Data/DEA Analysis/pass1b-06_transcript-rna-seq_dea.RData")
load(file = "PASS1B Transcription Factor Paper Data/DEA Analysis/pass1b-06_epigen-atac-seq_dea.RData")
load(file = "PASS1B Transcription Factor Paper Data/DEA Analysis/pass1b-06_proteomics-prot-pr_dea.RData")

peakanno <- read.delim(file = "PASS1B Transcription Factor Paper Data/ATACSeq Data/pass1b-06_epigen-atac-seq_feature-mapping.txt",header = T,row.names = 2)
peakanno$trimanno <- gsub(" \\(.*","",peakanno$custom_annotation)

tfanno <- read.csv("PASS1B Transcription Factor Paper Data/TF Data/tflist.csv",header = T,row.names = 1)
tfanno$Ensembl <- gsub(" ","",tfanno$Ensembl)
tfanno <- tfanno[!tfanno$Ensembl == "",]

genemapping <- read.delim(file = "PASS1B Transcription Factor Paper Data/RNASeq Data/pass1b-06_transcript-rna-seq_feature-mapping.txt",header = T,sep = "\t")
enstosym <- data.frame(row.names = unique(genemapping$ensembl_gene),"Symbol" = rep("test",length(unique(genemapping$ensembl_gene))))
for(i in 1:length(unique(genemapping$ensembl_gene))){
  if(i%%1000 == 0){
    print(i)
  }
  ourens <- rownames(enstosym)[i]
  enstosym[i,"Symbol"] <- genemapping[genemapping$ensembl_gene %in% ourens,"gene_symbol"][1]
}
enstosym$Ensembl <- rownames(enstosym)
for(i in 1:dim(enstosym)[1]){
  if(is.na(enstosym[i,"Symbol"])){
    enstosym[i,"Symbol"] <- enstosym[i,"Ensembl"]
  }
}
enstosym["ENSRNOG00000055858","Symbol"] <- "Myb"

allpeakmotifs <- read.table(file = "PASS1B Transcription Factor Paper Data/TF Data/allpeaksmotifs.txt",header = T,sep = "\t")


# Generate annotated list of protein id's associated with transcription factors
gastroproanno <- read.delim(file = "PASS1B Transcription Factor Paper Data/Prot_PR Data/t55-gastrocnemius_feature_annotation.txt",header = T,row.names = 1,sep = "\t")
heartproanno <- read.delim(file = "PASS1B Transcription Factor Paper Data/Prot_PR Data/t58-heart_feature_annotations.txt",header = T,row.names = 1,sep = "\t")
kidneyproanno <- read.delim(file = "PASS1B Transcription Factor Paper Data/Prot_PR Data/t59-kidney_feature_annotations.txt",header = T,row.names = 1,sep = "\t")
liverproanno <- read.delim(file = "PASS1B Transcription Factor Paper Data/Prot_PR Data/t68-liver_feature_annotations.txt",header = T,row.names = 1,sep = "\t")
lungproanno <- read.delim(file = "PASS1B Transcription Factor Paper Data/Prot_PR Data/t66-lung_feature_annotations.txt",header = T,row.names = 1,sep = "\t")
whiteproanno <- read.delim(file = "PASS1B Transcription Factor Paper Data/Prot_PR Data/t70-white-adipose_feature_annotations.txt",header = T,row.names = 1,sep = "\t")

tfproanno <- data.frame(row.names = rownames(tfanno),
                        "Gene.Name" = tfanno$Gene.Name,
                        "Ensembl" = tfanno$Ensembl,
                        "Gastro.Pro.ID" = "",
                        "Heart.Pro.ID" = "",
                        "Kidney.Pro.ID" = "",
                        "Liver.Pro.ID" = "",
                        "Lung.Pro.ID" = "",
                        "WhiteAd.Pro.ID" = "")

for(i in 1:dim(tfanno)[1]){
  
  if(i%%20 == 0){
    print(i)
  }
  
  ourgene <- tfanno$Gene.Name[i]
  if(ourgene %in% gastroproanno$gene_symbol){
    ourpro <- gastroproanno[gastroproanno$gene_symbol %in% ourgene,]
    tfproanno[i,"Gastro.Pro.ID"] <- ourpro$protein_id[1]
  }
  
  if(ourgene %in% heartproanno$gene_symbol){
    ourpro <- heartproanno[heartproanno$gene_symbol %in% ourgene,]
    tfproanno[i,"Heart.Pro.ID"] <- ourpro$protein_id[1]
  }
  
  if(ourgene %in% kidneyproanno$gene_symbol){
    ourpro <- kidneyproanno[kidneyproanno$gene_symbol %in% ourgene,]
    tfproanno[i,"Kidney.Pro.ID"] <- ourpro$protein_id[1]
  }
  
  if(ourgene %in% liverproanno$gene_symbol){
    ourpro <- liverproanno[liverproanno$gene_symbol %in% ourgene,]
    tfproanno[i,"Liver.Pro.ID"] <- ourpro$protein_id[1]
  }
  
  if(ourgene %in% lungproanno$gene_symbol){
    ourpro <- lungproanno[lungproanno$gene_symbol %in% ourgene,]
    tfproanno[i,"Lung.Pro.ID"] <- ourpro$protein_id[1]
  }
  
  if(ourgene %in% whiteproanno$gene_symbol){
    ourpro <- whiteproanno[whiteproanno$gene_symbol %in% ourgene,]
    tfproanno[i,"WhiteAd.Pro.ID"] <- ourpro$protein_id[1]
  }
  
}


####
# Load normalized data
#####

gastrornanorm <- read.delim(file = "PASS1B Transcription Factor Paper Data/RNASeq Data/pass1b-06_transcript-rna-seq_normalized-data_t55-gastrocnemius.txt",header = T,row.names = 1)
colnames(gastrornanorm) <- gsub("X","",colnames(gastrornanorm))
heartrnanorm <- read.delim(file = "PASS1B Transcription Factor Paper Data/RNASeq Data/pass1b-06_transcript-rna-seq_normalized-data_t58-heart.txt",header = T,row.names = 1)
colnames(heartrnanorm) <- gsub("X","",colnames(heartrnanorm))
hippornanorm <- read.delim(file = "PASS1B Transcription Factor Paper Data/RNASeq Data/pass1b-06_transcript-rna-seq_normalized-data_t52-hippocampus.txt",header = T,row.names = 1)
colnames(hippornanorm) <- gsub("X","",colnames(hippornanorm))
kidneyrnanorm <- read.delim(file = "PASS1B Transcription Factor Paper Data/RNASeq Data/pass1b-06_transcript-rna-seq_normalized-data_t59-kidney.txt",header = T,row.names = 1)
colnames(kidneyrnanorm) <- gsub("X","",colnames(kidneyrnanorm))
liverrnanorm <- read.delim(file = "PASS1B Transcription Factor Paper Data/RNASeq Data/pass1b-06_transcript-rna-seq_normalized-data_t68-liver.txt",header = T,row.names = 1)
colnames(liverrnanorm) <- gsub("X","",colnames(liverrnanorm))
lungrnanorm <- read.delim(file = "PASS1B Transcription Factor Paper Data/RNASeq Data/pass1b-06_transcript-rna-seq_normalized-data_t66-lung.txt",header = T,row.names = 1)
colnames(lungrnanorm) <- gsub("X","",colnames(lungrnanorm))
brownrnanorm <- read.delim(file = "PASS1B Transcription Factor Paper Data/RNASeq Data/pass1b-06_transcript-rna-seq_normalized-data_t69-brown-adipose.txt",header = T,row.names = 1)
colnames(brownrnanorm) <- gsub("X","",colnames(brownrnanorm))
whiternanorm <- read.delim(file = "PASS1B Transcription Factor Paper Data/RNASeq Data/pass1b-06_transcript-rna-seq_normalized-data_t70-white-adipose.txt",header = T,row.names = 1)
colnames(whiternanorm) <- gsub("X","",colnames(whiternanorm))

gastroatacnorm <- read.delim(file = "PASS1B Transcription Factor Paper Data/ATACSeq Data/pass1b-06_epigen-atac-seq_quant-norm_t55-gastrocnemius.tsv",header = T,row.names = 1)
colnames(gastroatacnorm) <- gsub("X","",colnames(gastroatacnorm))
heartatacnorm <- read.delim(file = "PASS1B Transcription Factor Paper Data/ATACSeq Data/pass1b-06_epigen-atac-seq_quant-norm_t58-heart.tsv",header = T,row.names = 1)
colnames(heartatacnorm) <- gsub("X","",colnames(heartatacnorm))
hippoatacnorm <- read.delim(file = "PASS1B Transcription Factor Paper Data/ATACSeq Data/pass1b-06_epigen-atac-seq_quant-norm_t52-hippocampus.tsv",header = T,row.names = 1)
colnames(hippoatacnorm) <- gsub("X","",colnames(hippoatacnorm))
kidneyatacnorm <- read.delim(file = "PASS1B Transcription Factor Paper Data/ATACSeq Data/pass1b-06_epigen-atac-seq_quant-norm_t59-kidney.tsv",header = T,row.names = 1)
colnames(kidneyatacnorm) <- gsub("X","",colnames(kidneyatacnorm))
liveratacnorm <- read.delim(file = "PASS1B Transcription Factor Paper Data/ATACSeq Data/pass1b-06_epigen-atac-seq_quant-norm_t68-liver.tsv",header = T,row.names = 1)
colnames(liveratacnorm) <- gsub("X","",colnames(liveratacnorm))
lungatacnorm <- read.delim(file = "PASS1B Transcription Factor Paper Data/ATACSeq Data/pass1b-06_epigen-atac-seq_quant-norm_t66-lung.tsv",header = T,row.names = 1)
colnames(lungatacnorm) <- gsub("X","",colnames(lungatacnorm))
brownatacnorm <- read.delim(file = "PASS1B Transcription Factor Paper Data/ATACSeq Data/pass1b-06_epigen-atac-seq_quant-norm_t69-brown-adipose.tsv",header = T,row.names = 1)
colnames(brownatacnorm) <- gsub("X","",colnames(brownatacnorm))
whiteatacnorm <- read.delim(file = "PASS1B Transcription Factor Paper Data/ATACSeq Data/pass1b-06_epigen-atac-seq_quant-norm_t70-white-adipose.tsv",header = T,row.names = 1)
colnames(whiteatacnorm) <- gsub("X","",colnames(whiteatacnorm))

gastroatacnorm <- gastroatacnorm[3:dim(gastroatacnorm)[1],]
heartatacnorm <- heartatacnorm[3:dim(heartatacnorm)[1],]
hippoatacnorm <- hippoatacnorm[3:dim(hippoatacnorm)[1],]
kidneyatacnorm <- kidneyatacnorm[3:dim(kidneyatacnorm)[1],]
liveratacnorm <- liveratacnorm[3:dim(liveratacnorm)[1],]
lungatacnorm <- lungatacnorm[3:dim(lungatacnorm)[1],]
whiteatacnorm <- whiteatacnorm[3:dim(whiteatacnorm)[1],]
brownatacnorm <- brownatacnorm[3:dim(brownatacnorm)[1],]

# Specify active peaks in each tissue for downstream analysis
activepeakcutoff <- -1

gastroactivepeaks <- rownames(gastroatacnorm[apply(gastroatacnorm,1,median) > activepeakcutoff,])
heartactivepeaks <- rownames(heartatacnorm[apply(heartatacnorm,1,median) > activepeakcutoff,])
hippoactivepeaks <- rownames(hippoatacnorm[apply(hippoatacnorm,1,median) > activepeakcutoff,])
kidneyactivepeaks <- rownames(kidneyatacnorm[apply(kidneyatacnorm,1,median) > activepeakcutoff,])
liveractivepeaks <- rownames(liveratacnorm[apply(liveratacnorm,1,median) > activepeakcutoff,])
lungactivepeaks <- rownames(lungatacnorm[apply(lungatacnorm,1,median) > activepeakcutoff,])
brownactivepeaks <- rownames(brownatacnorm[apply(brownatacnorm,1,median) > activepeakcutoff,])
whiteactivepeaks <- rownames(whiteatacnorm[apply(whiteatacnorm,1,median) > activepeakcutoff,])



####
# Generate L2FC Matrices for RNA, ATAC and Proteomics Data
#####

gastrotimewise <- transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue %in% "t55-gastrocnemius",]
hearttimewise <- transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue %in% "t58-heart",]
hippotimewise <- transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue %in% "t52-hippocampus",]
kidneytimewise <- transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue %in% "t59-kidney",]
livertimewise <- transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue %in% "t68-liver",]
lungtimewise <- transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue %in% "t66-lung",]
browntimewise <- transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue %in% "t69-brown-adipose",]
whitetimewise <- transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue %in% "t70-white-adipose",]

gastrol2fcmat <- matrix(0L,nrow = length(unique(gastrotimewise$feature_ID)),ncol = 8)
rownames(gastrol2fcmat) <- unique(gastrotimewise$feature_ID)
colnames(gastrol2fcmat) <- c("F W1","F W2","F W4","F W8",
                             "M W1","M W2","M W4","M W8")
for(i in 1:dim(gastrol2fcmat)[1]){
  
  if(i%%100 == 0){
    print(i)
  }
  ourgene <- rownames(gastrol2fcmat)[i]
  ourgenemat <- gastrotimewise[gastrotimewise$feature_ID %in% ourgene,]
  gastrol2fcmat[i,"F W1"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "1w","logFC"]
  gastrol2fcmat[i,"F W2"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "2w","logFC"]
  gastrol2fcmat[i,"F W4"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "4w","logFC"]
  gastrol2fcmat[i,"F W8"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "8w","logFC"]
  gastrol2fcmat[i,"M W1"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "1w","logFC"]
  gastrol2fcmat[i,"M W2"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "2w","logFC"]
  gastrol2fcmat[i,"M W4"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "4w","logFC"]
  gastrol2fcmat[i,"M W8"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "8w","logFC"]
  
}

heartl2fcmat <- matrix(0L,nrow = length(unique(hearttimewise$feature_ID)),ncol = 8)
rownames(heartl2fcmat) <- unique(hearttimewise$feature_ID)
colnames(heartl2fcmat) <- c("F W1","F W2","F W4","F W8",
                            "M W1","M W2","M W4","M W8")
for(i in 1:dim(heartl2fcmat)[1]){
  
  if(i%%100 == 0){
    print(i)
  }
  
  ourgene <- rownames(heartl2fcmat)[i]
  ourgenemat <- hearttimewise[hearttimewise$feature_ID %in% ourgene,]
  heartl2fcmat[i,"F W1"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "1w","logFC"]
  heartl2fcmat[i,"F W2"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "2w","logFC"]
  heartl2fcmat[i,"F W4"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "4w","logFC"]
  heartl2fcmat[i,"F W8"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "8w","logFC"]
  heartl2fcmat[i,"M W1"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "1w","logFC"]
  heartl2fcmat[i,"M W2"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "2w","logFC"]
  heartl2fcmat[i,"M W4"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "4w","logFC"]
  heartl2fcmat[i,"M W8"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "8w","logFC"]
  
}

hippol2fcmat <- matrix(0L,nrow = length(unique(hippotimewise$feature_ID)),ncol = 8)
rownames(hippol2fcmat) <- unique(hippotimewise$feature_ID)
colnames(hippol2fcmat) <- c("F W1","F W2","F W4","F W8",
                            "M W1","M W2","M W4","M W8")
for(i in 1:dim(hippol2fcmat)[1]){
  
  if(i%%100 == 0){
    print(i)
  }
  
  ourgene <- rownames(hippol2fcmat)[i]
  ourgenemat <- hippotimewise[hippotimewise$feature_ID %in% ourgene,]
  hippol2fcmat[i,"F W1"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "1w","logFC"]
  hippol2fcmat[i,"F W2"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "2w","logFC"]
  hippol2fcmat[i,"F W4"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "4w","logFC"]
  hippol2fcmat[i,"F W8"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "8w","logFC"]
  hippol2fcmat[i,"M W1"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "1w","logFC"]
  hippol2fcmat[i,"M W2"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "2w","logFC"]
  hippol2fcmat[i,"M W4"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "4w","logFC"]
  hippol2fcmat[i,"M W8"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "8w","logFC"]
  
}

kidneyl2fcmat <- matrix(0L,nrow = length(unique(kidneytimewise$feature_ID)),ncol = 8)
rownames(kidneyl2fcmat) <- unique(kidneytimewise$feature_ID)
colnames(kidneyl2fcmat) <- c("F W1","F W2","F W4","F W8",
                             "M W1","M W2","M W4","M W8")
for(i in 1:dim(kidneyl2fcmat)[1]){
  
  if(i%%100 == 0){
    print(i)
  }
  
  ourgene <- rownames(kidneyl2fcmat)[i]
  ourgenemat <- kidneytimewise[kidneytimewise$feature_ID %in% ourgene,]
  kidneyl2fcmat[i,"F W1"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "1w","logFC"]
  kidneyl2fcmat[i,"F W2"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "2w","logFC"]
  kidneyl2fcmat[i,"F W4"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "4w","logFC"]
  kidneyl2fcmat[i,"F W8"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "8w","logFC"]
  kidneyl2fcmat[i,"M W1"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "1w","logFC"]
  kidneyl2fcmat[i,"M W2"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "2w","logFC"]
  kidneyl2fcmat[i,"M W4"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "4w","logFC"]
  kidneyl2fcmat[i,"M W8"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "8w","logFC"]
  
}


liverl2fcmat <- matrix(0L,nrow = length(unique(livertimewise$feature_ID)),ncol = 8)
rownames(liverl2fcmat) <- unique(livertimewise$feature_ID)
colnames(liverl2fcmat) <- c("F W1","F W2","F W4","F W8",
                            "M W1","M W2","M W4","M W8")
for(i in 1:dim(liverl2fcmat)[1]){
  
  if(i%%100 == 0){
    print(i)
  }
  
  ourgene <- rownames(liverl2fcmat)[i]
  ourgenemat <- livertimewise[livertimewise$feature_ID %in% ourgene,]
  liverl2fcmat[i,"F W1"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "1w","logFC"]
  liverl2fcmat[i,"F W2"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "2w","logFC"]
  liverl2fcmat[i,"F W4"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "4w","logFC"]
  liverl2fcmat[i,"F W8"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "8w","logFC"]
  liverl2fcmat[i,"M W1"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "1w","logFC"]
  liverl2fcmat[i,"M W2"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "2w","logFC"]
  liverl2fcmat[i,"M W4"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "4w","logFC"]
  liverl2fcmat[i,"M W8"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "8w","logFC"]
  
}

lungl2fcmat <- matrix(0L,nrow = length(unique(lungtimewise$feature_ID)),ncol = 8)
rownames(lungl2fcmat) <- unique(lungtimewise$feature_ID)
colnames(lungl2fcmat) <- c("F W1","F W2","F W4","F W8",
                           "M W1","M W2","M W4","M W8")
for(i in 1:dim(lungl2fcmat)[1]){
  
  if(i%%100 == 0){
    print(i)
  }
  
  ourgene <- rownames(lungl2fcmat)[i]
  ourgenemat <- lungtimewise[lungtimewise$feature_ID %in% ourgene,]
  lungl2fcmat[i,"F W1"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "1w","logFC"]
  lungl2fcmat[i,"F W2"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "2w","logFC"]
  lungl2fcmat[i,"F W4"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "4w","logFC"]
  lungl2fcmat[i,"F W8"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "8w","logFC"]
  lungl2fcmat[i,"M W1"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "1w","logFC"]
  lungl2fcmat[i,"M W2"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "2w","logFC"]
  lungl2fcmat[i,"M W4"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "4w","logFC"]
  lungl2fcmat[i,"M W8"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "8w","logFC"]
  
}

brownl2fcmat <- matrix(0L,nrow = length(unique(browntimewise$feature_ID)),ncol = 8)
rownames(brownl2fcmat) <- unique(browntimewise$feature_ID)
colnames(brownl2fcmat) <- c("F W1","F W2","F W4","F W8",
                            "M W1","M W2","M W4","M W8")
for(i in 1:dim(brownl2fcmat)[1]){
  
  if(i%%100 == 0){
    print(i)
  }
  
  ourgene <- rownames(brownl2fcmat)[i]
  ourgenemat <- browntimewise[browntimewise$feature_ID %in% ourgene,]
  brownl2fcmat[i,"F W1"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "1w","logFC"]
  brownl2fcmat[i,"F W2"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "2w","logFC"]
  brownl2fcmat[i,"F W4"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "4w","logFC"]
  brownl2fcmat[i,"F W8"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "8w","logFC"]
  brownl2fcmat[i,"M W1"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "1w","logFC"]
  brownl2fcmat[i,"M W2"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "2w","logFC"]
  brownl2fcmat[i,"M W4"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "4w","logFC"]
  brownl2fcmat[i,"M W8"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "8w","logFC"]
  
}


whitel2fcmat <- matrix(0L,nrow = length(unique(whitetimewise$feature_ID)),ncol = 8)
rownames(whitel2fcmat) <- unique(whitetimewise$feature_ID)
colnames(whitel2fcmat) <- c("F W1","F W2","F W4","F W8",
                            "M W1","M W2","M W4","M W8")
for(i in 1:dim(whitel2fcmat)[1]){
  
  if(i%%100 == 0){
    print(i)
  }
  
  ourgene <- rownames(whitel2fcmat)[i]
  ourgenemat <- whitetimewise[whitetimewise$feature_ID %in% ourgene,]
  whitel2fcmat[i,"F W1"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "1w","logFC"]
  whitel2fcmat[i,"F W2"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "2w","logFC"]
  whitel2fcmat[i,"F W4"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "4w","logFC"]
  whitel2fcmat[i,"F W8"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "8w","logFC"]
  whitel2fcmat[i,"M W1"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "1w","logFC"]
  whitel2fcmat[i,"M W2"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "2w","logFC"]
  whitel2fcmat[i,"M W4"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "4w","logFC"]
  whitel2fcmat[i,"M W8"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "8w","logFC"]
  
}

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


gastrosigl2fcmat <- gastrol2fcmat[gastrornasig,]
heartsigl2fcmat <- heartl2fcmat[heartrnasig,]
hipposigl2fcmat <- hippol2fcmat[hippornasig,]
kidneysigl2fcmat <- kidneyl2fcmat[kidneyrnasig,]
liversigl2fcmat <- liverl2fcmat[liverrnasig,]
lungsigl2fcmat <- lungl2fcmat[lungrnasig,]
brownsigl2fcmat <- brownl2fcmat[brownrnasig,]
whitesigl2fcmat <- whitel2fcmat[whiternasig,]

gastroatactimewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "SKM-GN" & epigen_atac_seq$timewise_dea$feature_ID %in% gastroatacsig,]
heartatactimewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "HEART" & epigen_atac_seq$timewise_dea$feature_ID %in% heartatacsig,]
hippoatactimewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "HIPPOC" & epigen_atac_seq$timewise_dea$feature_ID %in% hippoatacsig,]
kidneyatactimewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "KIDNEY" & epigen_atac_seq$timewise_dea$feature_ID %in% kidneyatacsig,]
liveratactimewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "LIVER" & epigen_atac_seq$timewise_dea$feature_ID %in% liveratacsig,]
lungatactimewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "LUNG" & epigen_atac_seq$timewise_dea$feature_ID %in% lungatacsig,]
brownatactimewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "BAT" & epigen_atac_seq$timewise_dea$feature_ID %in% brownatacsig,]
whiteatactimewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "WAT-SC" & epigen_atac_seq$timewise_dea$feature_ID %in% whiteatacsig,]

gastrosigatacl2fc <- matrix(0L,nrow = length(gastroatacsig),ncol = 8)
rownames(gastrosigatacl2fc) <- gastroatacsig
colnames(gastrosigatacl2fc) <- c("F W1","F W2","F W4","F W8",
                                 "M W1","M W2","M W4","M W8")
for(i in 1:length(gastroatacsig)){
  if(i%%20 == 0){
    print(paste("gastro_",toString(i)))
  }
  ourpeak <- gastroatacsig[i]
  ourpeaktimewise <- gastroatactimewise[gastroatactimewise$feature_ID %in% ourpeak,]
  gastrosigatacl2fc[i,"F W1"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "female" & ourpeaktimewise$comparison_group %in% "1w","logFC"]
  gastrosigatacl2fc[i,"F W2"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "female" & ourpeaktimewise$comparison_group %in% "2w","logFC"]
  gastrosigatacl2fc[i,"F W4"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "female" & ourpeaktimewise$comparison_group %in% "4w","logFC"]
  gastrosigatacl2fc[i,"F W8"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "female" & ourpeaktimewise$comparison_group %in% "8w","logFC"]
  gastrosigatacl2fc[i,"M W1"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "male" & ourpeaktimewise$comparison_group %in% "1w","logFC"]
  gastrosigatacl2fc[i,"M W2"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "male" & ourpeaktimewise$comparison_group %in% "2w","logFC"]
  gastrosigatacl2fc[i,"M W4"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "male" & ourpeaktimewise$comparison_group %in% "4w","logFC"]
  gastrosigatacl2fc[i,"M W8"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "male" & ourpeaktimewise$comparison_group %in% "8w","logFC"]
}


heartsigatacl2fc <- matrix(0L,nrow = length(heartatacsig),ncol = 8)
rownames(heartsigatacl2fc) <- heartatacsig
colnames(heartsigatacl2fc) <- c("F W1","F W2","F W4","F W8",
                                "M W1","M W2","M W4","M W8")
for(i in 1:length(heartatacsig)){
  if(i%%20 == 0){
    print(paste("heart_",toString(i)))
  }
  ourpeak <- heartatacsig[i]
  ourpeaktimewise <- heartatactimewise[heartatactimewise$feature_ID %in% ourpeak,]
  heartsigatacl2fc[i,"F W1"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "female" & ourpeaktimewise$comparison_group %in% "1w","logFC"]
  heartsigatacl2fc[i,"F W2"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "female" & ourpeaktimewise$comparison_group %in% "2w","logFC"]
  heartsigatacl2fc[i,"F W4"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "female" & ourpeaktimewise$comparison_group %in% "4w","logFC"]
  heartsigatacl2fc[i,"F W8"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "female" & ourpeaktimewise$comparison_group %in% "8w","logFC"]
  heartsigatacl2fc[i,"M W1"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "male" & ourpeaktimewise$comparison_group %in% "1w","logFC"]
  heartsigatacl2fc[i,"M W2"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "male" & ourpeaktimewise$comparison_group %in% "2w","logFC"]
  heartsigatacl2fc[i,"M W4"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "male" & ourpeaktimewise$comparison_group %in% "4w","logFC"]
  heartsigatacl2fc[i,"M W8"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "male" & ourpeaktimewise$comparison_group %in% "8w","logFC"]
}


hipposigatacl2fc <- matrix(0L,nrow = length(hippoatacsig),ncol = 8)
rownames(hipposigatacl2fc) <- hippoatacsig
colnames(hipposigatacl2fc) <- c("F W1","F W2","F W4","F W8",
                                "M W1","M W2","M W4","M W8")
for(i in 1:length(hippoatacsig)){
  if(i%%20 == 0){
    print(paste("hippo_",toString(i)))
  }
  ourpeak <- hippoatacsig[i]
  ourpeaktimewise <- hippoatactimewise[hippoatactimewise$feature_ID %in% ourpeak,]
  hipposigatacl2fc[i,"F W1"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "female" & ourpeaktimewise$comparison_group %in% "1w","logFC"]
  hipposigatacl2fc[i,"F W2"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "female" & ourpeaktimewise$comparison_group %in% "2w","logFC"]
  hipposigatacl2fc[i,"F W4"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "female" & ourpeaktimewise$comparison_group %in% "4w","logFC"]
  hipposigatacl2fc[i,"F W8"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "female" & ourpeaktimewise$comparison_group %in% "8w","logFC"]
  hipposigatacl2fc[i,"M W1"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "male" & ourpeaktimewise$comparison_group %in% "1w","logFC"]
  hipposigatacl2fc[i,"M W2"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "male" & ourpeaktimewise$comparison_group %in% "2w","logFC"]
  hipposigatacl2fc[i,"M W4"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "male" & ourpeaktimewise$comparison_group %in% "4w","logFC"]
  hipposigatacl2fc[i,"M W8"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "male" & ourpeaktimewise$comparison_group %in% "8w","logFC"]
}


kidneysigatacl2fc <- matrix(0L,nrow = length(kidneyatacsig),ncol = 8)
rownames(kidneysigatacl2fc) <- kidneyatacsig
colnames(kidneysigatacl2fc) <- c("F W1","F W2","F W4","F W8",
                                 "M W1","M W2","M W4","M W8")
for(i in 1:length(kidneyatacsig)){
  if(i%%20 == 0){
    print(paste("kidney_",toString(i)))
  }
  ourpeak <- kidneyatacsig[i]
  ourpeaktimewise <- kidneyatactimewise[kidneyatactimewise$feature_ID %in% ourpeak,]
  kidneysigatacl2fc[i,"F W1"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "female" & ourpeaktimewise$comparison_group %in% "1w","logFC"]
  kidneysigatacl2fc[i,"F W2"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "female" & ourpeaktimewise$comparison_group %in% "2w","logFC"]
  kidneysigatacl2fc[i,"F W4"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "female" & ourpeaktimewise$comparison_group %in% "4w","logFC"]
  kidneysigatacl2fc[i,"F W8"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "female" & ourpeaktimewise$comparison_group %in% "8w","logFC"]
  kidneysigatacl2fc[i,"M W1"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "male" & ourpeaktimewise$comparison_group %in% "1w","logFC"]
  kidneysigatacl2fc[i,"M W2"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "male" & ourpeaktimewise$comparison_group %in% "2w","logFC"]
  kidneysigatacl2fc[i,"M W4"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "male" & ourpeaktimewise$comparison_group %in% "4w","logFC"]
  kidneysigatacl2fc[i,"M W8"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "male" & ourpeaktimewise$comparison_group %in% "8w","logFC"]
}


liversigatacl2fc <- matrix(0L,nrow = length(liveratacsig),ncol = 8)
rownames(liversigatacl2fc) <- liveratacsig
colnames(liversigatacl2fc) <- c("F W1","F W2","F W4","F W8",
                                "M W1","M W2","M W4","M W8")
for(i in 1:length(liveratacsig)){
  if(i%%20 == 0){
    print(paste("liver_",toString(i)))
  }
  ourpeak <- liveratacsig[i]
  ourpeaktimewise <- liveratactimewise[liveratactimewise$feature_ID %in% ourpeak,]
  liversigatacl2fc[i,"F W1"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "female" & ourpeaktimewise$comparison_group %in% "1w","logFC"]
  liversigatacl2fc[i,"F W2"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "female" & ourpeaktimewise$comparison_group %in% "2w","logFC"]
  liversigatacl2fc[i,"F W4"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "female" & ourpeaktimewise$comparison_group %in% "4w","logFC"]
  liversigatacl2fc[i,"F W8"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "female" & ourpeaktimewise$comparison_group %in% "8w","logFC"]
  liversigatacl2fc[i,"M W1"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "male" & ourpeaktimewise$comparison_group %in% "1w","logFC"]
  liversigatacl2fc[i,"M W2"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "male" & ourpeaktimewise$comparison_group %in% "2w","logFC"]
  liversigatacl2fc[i,"M W4"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "male" & ourpeaktimewise$comparison_group %in% "4w","logFC"]
  liversigatacl2fc[i,"M W8"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "male" & ourpeaktimewise$comparison_group %in% "8w","logFC"]
}


lungsigatacl2fc <- matrix(0L,nrow = length(lungatacsig),ncol = 8)
rownames(lungsigatacl2fc) <- lungatacsig
colnames(lungsigatacl2fc) <- c("F W1","F W2","F W4","F W8",
                               "M W1","M W2","M W4","M W8")
for(i in 1:length(lungatacsig)){
  if(i%%20 == 0){
    print(paste("lung_",toString(i)))
  }
  ourpeak <- lungatacsig[i]
  ourpeaktimewise <- lungatactimewise[lungatactimewise$feature_ID %in% ourpeak,]
  lungsigatacl2fc[i,"F W1"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "female" & ourpeaktimewise$comparison_group %in% "1w","logFC"]
  lungsigatacl2fc[i,"F W2"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "female" & ourpeaktimewise$comparison_group %in% "2w","logFC"]
  lungsigatacl2fc[i,"F W4"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "female" & ourpeaktimewise$comparison_group %in% "4w","logFC"]
  lungsigatacl2fc[i,"F W8"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "female" & ourpeaktimewise$comparison_group %in% "8w","logFC"]
  lungsigatacl2fc[i,"M W1"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "male" & ourpeaktimewise$comparison_group %in% "1w","logFC"]
  lungsigatacl2fc[i,"M W2"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "male" & ourpeaktimewise$comparison_group %in% "2w","logFC"]
  lungsigatacl2fc[i,"M W4"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "male" & ourpeaktimewise$comparison_group %in% "4w","logFC"]
  lungsigatacl2fc[i,"M W8"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "male" & ourpeaktimewise$comparison_group %in% "8w","logFC"]
}


brownsigatacl2fc <- matrix(0L,nrow = length(brownatacsig),ncol = 8)
rownames(brownsigatacl2fc) <- brownatacsig
colnames(brownsigatacl2fc) <- c("F W1","F W2","F W4","F W8",
                                "M W1","M W2","M W4","M W8")
for(i in 1:length(brownatacsig)){
  if(i%%20 == 0){
    print(paste("brown_",toString(i)))
  }
  ourpeak <- brownatacsig[i]
  ourpeaktimewise <- brownatactimewise[brownatactimewise$feature_ID %in% ourpeak,]
  brownsigatacl2fc[i,"F W1"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "female" & ourpeaktimewise$comparison_group %in% "1w","logFC"]
  brownsigatacl2fc[i,"F W2"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "female" & ourpeaktimewise$comparison_group %in% "2w","logFC"]
  brownsigatacl2fc[i,"F W4"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "female" & ourpeaktimewise$comparison_group %in% "4w","logFC"]
  brownsigatacl2fc[i,"F W8"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "female" & ourpeaktimewise$comparison_group %in% "8w","logFC"]
  brownsigatacl2fc[i,"M W1"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "male" & ourpeaktimewise$comparison_group %in% "1w","logFC"]
  brownsigatacl2fc[i,"M W2"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "male" & ourpeaktimewise$comparison_group %in% "2w","logFC"]
  brownsigatacl2fc[i,"M W4"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "male" & ourpeaktimewise$comparison_group %in% "4w","logFC"]
  brownsigatacl2fc[i,"M W8"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "male" & ourpeaktimewise$comparison_group %in% "8w","logFC"]
}


whitesigatacl2fc <- matrix(0L,nrow = length(whiteatacsig),ncol = 8)
rownames(whitesigatacl2fc) <- whiteatacsig
colnames(whitesigatacl2fc) <- c("F W1","F W2","F W4","F W8",
                                "M W1","M W2","M W4","M W8")
for(i in 1:length(whiteatacsig)){
  if(i%%20 == 0){
    print(paste("white_",toString(i)))
  }
  ourpeak <- whiteatacsig[i]
  ourpeaktimewise <- whiteatactimewise[whiteatactimewise$feature_ID %in% ourpeak,]
  whitesigatacl2fc[i,"F W1"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "female" & ourpeaktimewise$comparison_group %in% "1w","logFC"]
  whitesigatacl2fc[i,"F W2"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "female" & ourpeaktimewise$comparison_group %in% "2w","logFC"]
  whitesigatacl2fc[i,"F W4"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "female" & ourpeaktimewise$comparison_group %in% "4w","logFC"]
  whitesigatacl2fc[i,"F W8"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "female" & ourpeaktimewise$comparison_group %in% "8w","logFC"]
  whitesigatacl2fc[i,"M W1"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "male" & ourpeaktimewise$comparison_group %in% "1w","logFC"]
  whitesigatacl2fc[i,"M W2"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "male" & ourpeaktimewise$comparison_group %in% "2w","logFC"]
  whitesigatacl2fc[i,"M W4"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "male" & ourpeaktimewise$comparison_group %in% "4w","logFC"]
  whitesigatacl2fc[i,"M W8"] <- ourpeaktimewise[ourpeaktimewise$sex %in% "male" & ourpeaktimewise$comparison_group %in% "8w","logFC"]
}


####
# Color Annotation sets for future plots
#####

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

ann_cols_combocor <- list("Correlation" = c("Positive" = "#indianred2",
                                            "Negative" = "#skyblue2"))

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
                                          "-1" = "blue4"))


ann_cols_sigphenocor <- list("Measure" = c("Body Fat" = "#e41a1c",
                                             "Body Water" = "#377eb8",
                                             "Body Lean" = "#4daf4a",
                                             "Lactate" = "#984ea3",
                                             "VO2 Max" = "#ff7f00",
                                             "Body Weight" = "#ffff33"),
                               "Tissue" = c("SKM-GN" = "#088c03",
                                            "HEART" = "#f28b2f",
                                            "HIPPOC" = "#bf7534",
                                            "KIDNEY"= "#7553a7",
                                            "LIVER" = "#da6c75",
                                            "LUNG" = "#04bf8a",
                                            "BAT" = "#8c5220",
                                            "WAT-SC" = "#214da6"))

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

ann_colstrim <- list("Tissue" = c("SKM-GN" = "#088c03",
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
                     "Region" = c("Distal.Intergenic" = "#BCBD22FF",
                                  "Exon" = "#9467BDFF",
                                  "Intron" = "#8C564BFF",
                                  "Promoter.(1-2kb)" = "#2CA02CFF",
                                  "Promoter.(<=1kb)" = "#FF7F0EFF",
                                  "Upstream" = "#1F77B4FF"))
genomic_features_col <- list("Region" = c("3' UTR" = "#E377C2FF",
                                          "5' UTR" = "#D62728FF",
                                          "Distal Intergenic" = "#BCBD22FF",
                                          "Downstream (<5kb)" = "#7F7F7FFF",
                                          "Exon" = "#9467BDFF",
                                          "Intron" = "#8C564BFF",
                                          "Promoter (1-2kb)" = "#2CA02CFF",
                                          "Promoter (<=1kb)" = "#FF7F0EFF",
                                          "Upstream (<5kb)" = "#1F77B4FF",
                                          "Overlaps Gene" = "black"))

####
# Generation of peak sets for significant genes (DEGaPs) for each tissue
#####

gastrosigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrornasig,]),gastroactivepeaks)
gastrosigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrornasig & peakanno$trimanno %in% "Promoter",]),gastroactivepeaks)
gastrosigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrornasig & peakanno$trimanno %in% "Intron",]),gastroactivepeaks)
gastrosigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrornasig & peakanno$trimanno %in% "Distal Intergenic",]),gastroactivepeaks)
gastrosiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrornasig & peakanno$trimanno %in% "Upstream",]),gastroactivepeaks)
gastrosigdownpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrornasig & peakanno$trimanno %in% "Downstream",]),gastroactivepeaks)
gastrosigexpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrornasig & peakanno$trimanno %in% "Exon",]),gastroactivepeaks)
gastrosigpromproxpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrornasig & peakanno$custom_annotation %in% "Promoter (<=1kb)",]),gastroactivepeaks)
gastrosigpromfarpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrornasig & peakanno$custom_annotation %in% "Promoter (1-2kb)",]),gastroactivepeaks)

gastrow8upsig <- rownames(gastrosigl2fcmat[gastrosigl2fcmat[,"F W8"] > 0 & gastrosigl2fcmat[,"M W8"] > 0,])
gastrow8upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow8upsig,]),gastroactivepeaks)
gastrow8upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow8upsig & peakanno$trimanno %in% "Promoter",]),gastroactivepeaks)
gastrow8upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow8upsig & peakanno$trimanno %in% "Intron",]),gastroactivepeaks)
gastrow8upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow8upsig & peakanno$trimanno %in% "Distal Intergenic",]),gastroactivepeaks)
gastrow8upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow8upsig & peakanno$trimanno %in% "Upstream",]),gastroactivepeaks)

gastrow8downsig <- rownames(gastrosigl2fcmat[gastrosigl2fcmat[,"F W8"] < 0 & gastrosigl2fcmat[,"M W8"] < 0,])
gastrow8downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow8downsig,]),gastroactivepeaks)
gastrow8downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow8downsig & peakanno$trimanno %in% "Promoter",]),gastroactivepeaks)
gastrow8downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow8downsig & peakanno$trimanno %in% "Intron",]),gastroactivepeaks)
gastrow8downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow8downsig & peakanno$trimanno %in% "Distal Intergenic",]),gastroactivepeaks)
gastrow8downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow8downsig & peakanno$trimanno %in% "Upstream",]),gastroactivepeaks)

gastrow1upsig <- rownames(gastrosigl2fcmat[gastrosigl2fcmat[,"F W1"] > 0 & gastrosigl2fcmat[,"M W1"] > 0,])
gastrow1upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow1upsig,]),gastroactivepeaks)
gastrow1upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow1upsig & peakanno$trimanno %in% "Promoter",]),gastroactivepeaks)
gastrow1upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow1upsig & peakanno$trimanno %in% "Intron",]),gastroactivepeaks)
gastrow1upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow1upsig & peakanno$trimanno %in% "Distal Intergenic",]),gastroactivepeaks)
gastrow1upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow1upsig & peakanno$trimanno %in% "Upstream",]),gastroactivepeaks)

gastrow1downsig <- rownames(gastrosigl2fcmat[gastrosigl2fcmat[,"F W1"] < 0 & gastrosigl2fcmat[,"M W1"] < 0,])
gastrow1downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow1downsig,]),gastroactivepeaks)
gastrow1downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow1downsig & peakanno$trimanno %in% "Promoter",]),gastroactivepeaks)
gastrow1downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow1downsig & peakanno$trimanno %in% "Intron",]),gastroactivepeaks)
gastrow1downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow1downsig & peakanno$trimanno %in% "Distal Intergenic",]),gastroactivepeaks)
gastrow1downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow1downsig & peakanno$trimanno %in% "Upstream",]),gastroactivepeaks)


gastrow2upsig <- rownames(gastrosigl2fcmat[gastrosigl2fcmat[,"F W2"] > 0 & gastrosigl2fcmat[,"M W1"] > 0,])
gastrow2upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow2upsig,]),gastroactivepeaks)
gastrow2upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow2upsig & peakanno$trimanno %in% "Promoter",]),gastroactivepeaks)
gastrow2upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow2upsig & peakanno$trimanno %in% "Intron",]),gastroactivepeaks)
gastrow2upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow2upsig & peakanno$trimanno %in% "Distal Intergenic",]),gastroactivepeaks)
gastrow2upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow2upsig & peakanno$trimanno %in% "Upstream",]),gastroactivepeaks)

gastrow2downsig <- rownames(gastrosigl2fcmat[gastrosigl2fcmat[,"F W2"] < 0 & gastrosigl2fcmat[,"M W1"] < 0,])
gastrow2downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow2downsig,]),gastroactivepeaks)
gastrow2downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow2downsig & peakanno$trimanno %in% "Promoter",]),gastroactivepeaks)
gastrow2downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow2downsig & peakanno$trimanno %in% "Intron",]),gastroactivepeaks)
gastrow2downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow2downsig & peakanno$trimanno %in% "Distal Intergenic",]),gastroactivepeaks)
gastrow2downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow2downsig & peakanno$trimanno %in% "Upstream",]),gastroactivepeaks)


gastrow4upsig <- rownames(gastrosigl2fcmat[gastrosigl2fcmat[,"F W4"] > 0 & gastrosigl2fcmat[,"M W4"] > 0,])
gastrow4upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow4upsig,]),gastroactivepeaks)
gastrow4upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow4upsig & peakanno$trimanno %in% "Promoter",]),gastroactivepeaks)
gastrow4upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow4upsig & peakanno$trimanno %in% "Intron",]),gastroactivepeaks)
gastrow4upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow4upsig & peakanno$trimanno %in% "Distal Intergenic",]),gastroactivepeaks)
gastrow4upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrow4upsig & peakanno$trimanno %in% "Upstream",]),gastroactivepeaks)

gastrow4downsig <- rownames(gastrosigl2fcmat[gastrosigl2fcmat[,"F W4"] < 0 & gastrosigl2fcmat[,"M W4"] < 0,])
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

heartw8upsig <- rownames(heartsigl2fcmat[heartsigl2fcmat[,"F W8"] > 0 & heartsigl2fcmat[,"M W8"] > 0,])
heartw8upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw8upsig,]),heartactivepeaks)
heartw8upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw8upsig & peakanno$trimanno %in% "Promoter",]),heartactivepeaks)
heartw8upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw8upsig & peakanno$trimanno %in% "Intron",]),heartactivepeaks)
heartw8upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw8upsig & peakanno$trimanno %in% "Distal Intergenic",]),heartactivepeaks)
heartw8upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw8upsig & peakanno$trimanno %in% "Upstream",]),heartactivepeaks)

heartw8downsig <- rownames(heartsigl2fcmat[heartsigl2fcmat[,"F W8"] < 0 & heartsigl2fcmat[,"M W8"] < 0,])
heartw8downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw8downsig,]),heartactivepeaks)
heartw8downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw8downsig & peakanno$trimanno %in% "Promoter",]),heartactivepeaks)
heartw8downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw8downsig & peakanno$trimanno %in% "Intron",]),heartactivepeaks)
heartw8downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw8downsig & peakanno$trimanno %in% "Distal Intergenic",]),heartactivepeaks)
heartw8downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw8downsig & peakanno$trimanno %in% "Upstream",]),heartactivepeaks)

heartw1upsig <- rownames(heartsigl2fcmat[heartsigl2fcmat[,"F W1"] > 0 & heartsigl2fcmat[,"M W1"] > 0,])
heartw1upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw1upsig,]),heartactivepeaks)
heartw1upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw1upsig & peakanno$trimanno %in% "Promoter",]),heartactivepeaks)
heartw1upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw1upsig & peakanno$trimanno %in% "Intron",]),heartactivepeaks)
heartw1upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw1upsig & peakanno$trimanno %in% "Distal Intergenic",]),heartactivepeaks)
heartw1upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw1upsig & peakanno$trimanno %in% "Upstream",]),heartactivepeaks)

heartw1downsig <- rownames(heartsigl2fcmat[heartsigl2fcmat[,"F W1"] < 0 & heartsigl2fcmat[,"M W1"] < 0,])
heartw1downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw1downsig,]),heartactivepeaks)
heartw1downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw1downsig & peakanno$trimanno %in% "Promoter",]),heartactivepeaks)
heartw1downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw1downsig & peakanno$trimanno %in% "Intron",]),heartactivepeaks)
heartw1downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw1downsig & peakanno$trimanno %in% "Distal Intergenic",]),heartactivepeaks)
heartw1downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw1downsig & peakanno$trimanno %in% "Upstream",]),heartactivepeaks)


heartw2upsig <- rownames(heartsigl2fcmat[heartsigl2fcmat[,"F W2"] > 0 & heartsigl2fcmat[,"M W1"] > 0,])
heartw2upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw2upsig,]),heartactivepeaks)
heartw2upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw2upsig & peakanno$trimanno %in% "Promoter",]),heartactivepeaks)
heartw2upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw2upsig & peakanno$trimanno %in% "Intron",]),heartactivepeaks)
heartw2upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw2upsig & peakanno$trimanno %in% "Distal Intergenic",]),heartactivepeaks)
heartw2upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw2upsig & peakanno$trimanno %in% "Upstream",]),heartactivepeaks)

heartw2downsig <- rownames(heartsigl2fcmat[heartsigl2fcmat[,"F W2"] < 0 & heartsigl2fcmat[,"M W1"] < 0,])
heartw2downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw2downsig,]),heartactivepeaks)
heartw2downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw2downsig & peakanno$trimanno %in% "Promoter",]),heartactivepeaks)
heartw2downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw2downsig & peakanno$trimanno %in% "Intron",]),heartactivepeaks)
heartw2downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw2downsig & peakanno$trimanno %in% "Distal Intergenic",]),heartactivepeaks)
heartw2downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw2downsig & peakanno$trimanno %in% "Upstream",]),heartactivepeaks)


heartw4upsig <- rownames(heartsigl2fcmat[heartsigl2fcmat[,"F W4"] > 0 & heartsigl2fcmat[,"M W4"] > 0,])
heartw4upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw4upsig,]),heartactivepeaks)
heartw4upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw4upsig & peakanno$trimanno %in% "Promoter",]),heartactivepeaks)
heartw4upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw4upsig & peakanno$trimanno %in% "Intron",]),heartactivepeaks)
heartw4upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw4upsig & peakanno$trimanno %in% "Distal Intergenic",]),heartactivepeaks)
heartw4upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartw4upsig & peakanno$trimanno %in% "Upstream",]),heartactivepeaks)

heartw4downsig <- rownames(heartsigl2fcmat[heartsigl2fcmat[,"F W4"] < 0 & heartsigl2fcmat[,"M W4"] < 0,])
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

hippow8upsig <- rownames(hipposigl2fcmat[hipposigl2fcmat[,"F W8"] > 0 & hipposigl2fcmat[,"M W8"] > 0,])
hippow8upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow8upsig,]),hippoactivepeaks)
hippow8upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow8upsig & peakanno$trimanno %in% "Promoter",]),hippoactivepeaks)
hippow8upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow8upsig & peakanno$trimanno %in% "Intron",]),hippoactivepeaks)
hippow8upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow8upsig & peakanno$trimanno %in% "Distal Intergenic",]),hippoactivepeaks)
hippow8upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow8upsig & peakanno$trimanno %in% "Upstream",]),hippoactivepeaks)

hippow8downsig <- rownames(hipposigl2fcmat[hipposigl2fcmat[,"F W8"] < 0 & hipposigl2fcmat[,"M W8"] < 0,])
hippow8downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow8downsig,]),hippoactivepeaks)
hippow8downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow8downsig & peakanno$trimanno %in% "Promoter",]),hippoactivepeaks)
hippow8downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow8downsig & peakanno$trimanno %in% "Intron",]),hippoactivepeaks)
hippow8downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow8downsig & peakanno$trimanno %in% "Distal Intergenic",]),hippoactivepeaks)
hippow8downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow8downsig & peakanno$trimanno %in% "Upstream",]),hippoactivepeaks)

hippow1upsig <- rownames(hipposigl2fcmat[hipposigl2fcmat[,"F W1"] > 0 & hipposigl2fcmat[,"M W1"] > 0,])
hippow1upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow1upsig,]),hippoactivepeaks)
hippow1upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow1upsig & peakanno$trimanno %in% "Promoter",]),hippoactivepeaks)
hippow1upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow1upsig & peakanno$trimanno %in% "Intron",]),hippoactivepeaks)
hippow1upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow1upsig & peakanno$trimanno %in% "Distal Intergenic",]),hippoactivepeaks)
hippow1upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow1upsig & peakanno$trimanno %in% "Upstream",]),hippoactivepeaks)

hippow1downsig <- rownames(hipposigl2fcmat[hipposigl2fcmat[,"F W1"] < 0 & hipposigl2fcmat[,"M W1"] < 0,])
hippow1downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow1downsig,]),hippoactivepeaks)
hippow1downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow1downsig & peakanno$trimanno %in% "Promoter",]),hippoactivepeaks)
hippow1downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow1downsig & peakanno$trimanno %in% "Intron",]),hippoactivepeaks)
hippow1downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow1downsig & peakanno$trimanno %in% "Distal Intergenic",]),hippoactivepeaks)
hippow1downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow1downsig & peakanno$trimanno %in% "Upstream",]),hippoactivepeaks)


hippow2upsig <- rownames(hipposigl2fcmat[hipposigl2fcmat[,"F W2"] > 0 & hipposigl2fcmat[,"M W1"] > 0,])
hippow2upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow2upsig,]),hippoactivepeaks)
hippow2upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow2upsig & peakanno$trimanno %in% "Promoter",]),hippoactivepeaks)
hippow2upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow2upsig & peakanno$trimanno %in% "Intron",]),hippoactivepeaks)
hippow2upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow2upsig & peakanno$trimanno %in% "Distal Intergenic",]),hippoactivepeaks)
hippow2upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow2upsig & peakanno$trimanno %in% "Upstream",]),hippoactivepeaks)

hippow2downsig <- rownames(hipposigl2fcmat[hipposigl2fcmat[,"F W2"] < 0 & hipposigl2fcmat[,"M W1"] < 0,])
hippow2downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow2downsig,]),hippoactivepeaks)
hippow2downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow2downsig & peakanno$trimanno %in% "Promoter",]),hippoactivepeaks)
hippow2downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow2downsig & peakanno$trimanno %in% "Intron",]),hippoactivepeaks)
hippow2downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow2downsig & peakanno$trimanno %in% "Distal Intergenic",]),hippoactivepeaks)
hippow2downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow2downsig & peakanno$trimanno %in% "Upstream",]),hippoactivepeaks)


hippow4upsig <- rownames(hipposigl2fcmat[hipposigl2fcmat[,"F W4"] > 0 & hipposigl2fcmat[,"M W4"] > 0,])
hippow4upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow4upsig,]),hippoactivepeaks)
hippow4upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow4upsig & peakanno$trimanno %in% "Promoter",]),hippoactivepeaks)
hippow4upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow4upsig & peakanno$trimanno %in% "Intron",]),hippoactivepeaks)
hippow4upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow4upsig & peakanno$trimanno %in% "Distal Intergenic",]),hippoactivepeaks)
hippow4upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippow4upsig & peakanno$trimanno %in% "Upstream",]),hippoactivepeaks)

hippow4downsig <- rownames(hipposigl2fcmat[hipposigl2fcmat[,"F W4"] < 0 & hipposigl2fcmat[,"M W4"] < 0,])
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

kidneyw8upsig <- rownames(kidneysigl2fcmat[kidneysigl2fcmat[,"F W8"] > 0 & kidneysigl2fcmat[,"M W8"] > 0,])
kidneyw8upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw8upsig,]),kidneyactivepeaks)
kidneyw8upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw8upsig & peakanno$trimanno %in% "Promoter",]),kidneyactivepeaks)
kidneyw8upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw8upsig & peakanno$trimanno %in% "Intron",]),kidneyactivepeaks)
kidneyw8upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw8upsig & peakanno$trimanno %in% "Distal Intergenic",]),kidneyactivepeaks)
kidneyw8upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw8upsig & peakanno$trimanno %in% "Upstream",]),kidneyactivepeaks)

kidneyw8downsig <- rownames(kidneysigl2fcmat[kidneysigl2fcmat[,"F W8"] < 0 & kidneysigl2fcmat[,"M W8"] < 0,])
kidneyw8downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw8downsig,]),kidneyactivepeaks)
kidneyw8downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw8downsig & peakanno$trimanno %in% "Promoter",]),kidneyactivepeaks)
kidneyw8downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw8downsig & peakanno$trimanno %in% "Intron",]),kidneyactivepeaks)
kidneyw8downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw8downsig & peakanno$trimanno %in% "Distal Intergenic",]),kidneyactivepeaks)
kidneyw8downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw8downsig & peakanno$trimanno %in% "Upstream",]),kidneyactivepeaks)

kidneyw1upsig <- rownames(kidneysigl2fcmat[kidneysigl2fcmat[,"F W1"] > 0 & kidneysigl2fcmat[,"M W1"] > 0,])
kidneyw1upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw1upsig,]),kidneyactivepeaks)
kidneyw1upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw1upsig & peakanno$trimanno %in% "Promoter",]),kidneyactivepeaks)
kidneyw1upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw1upsig & peakanno$trimanno %in% "Intron",]),kidneyactivepeaks)
kidneyw1upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw1upsig & peakanno$trimanno %in% "Distal Intergenic",]),kidneyactivepeaks)
kidneyw1upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw1upsig & peakanno$trimanno %in% "Upstream",]),kidneyactivepeaks)

kidneyw1downsig <- rownames(kidneysigl2fcmat[kidneysigl2fcmat[,"F W1"] < 0 & kidneysigl2fcmat[,"M W1"] < 0,])
kidneyw1downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw1downsig,]),kidneyactivepeaks)
kidneyw1downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw1downsig & peakanno$trimanno %in% "Promoter",]),kidneyactivepeaks)
kidneyw1downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw1downsig & peakanno$trimanno %in% "Intron",]),kidneyactivepeaks)
kidneyw1downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw1downsig & peakanno$trimanno %in% "Distal Intergenic",]),kidneyactivepeaks)
kidneyw1downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw1downsig & peakanno$trimanno %in% "Upstream",]),kidneyactivepeaks)


kidneyw2upsig <- rownames(kidneysigl2fcmat[kidneysigl2fcmat[,"F W2"] > 0 & kidneysigl2fcmat[,"M W1"] > 0,])
kidneyw2upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw2upsig,]),kidneyactivepeaks)
kidneyw2upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw2upsig & peakanno$trimanno %in% "Promoter",]),kidneyactivepeaks)
kidneyw2upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw2upsig & peakanno$trimanno %in% "Intron",]),kidneyactivepeaks)
kidneyw2upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw2upsig & peakanno$trimanno %in% "Distal Intergenic",]),kidneyactivepeaks)
kidneyw2upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw2upsig & peakanno$trimanno %in% "Upstream",]),kidneyactivepeaks)

kidneyw2downsig <- rownames(kidneysigl2fcmat[kidneysigl2fcmat[,"F W2"] < 0 & kidneysigl2fcmat[,"M W1"] < 0,])
kidneyw2downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw2downsig,]),kidneyactivepeaks)
kidneyw2downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw2downsig & peakanno$trimanno %in% "Promoter",]),kidneyactivepeaks)
kidneyw2downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw2downsig & peakanno$trimanno %in% "Intron",]),kidneyactivepeaks)
kidneyw2downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw2downsig & peakanno$trimanno %in% "Distal Intergenic",]),kidneyactivepeaks)
kidneyw2downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw2downsig & peakanno$trimanno %in% "Upstream",]),kidneyactivepeaks)


kidneyw4upsig <- rownames(kidneysigl2fcmat[kidneysigl2fcmat[,"F W4"] > 0 & kidneysigl2fcmat[,"M W4"] > 0,])
kidneyw4upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw4upsig,]),kidneyactivepeaks)
kidneyw4upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw4upsig & peakanno$trimanno %in% "Promoter",]),kidneyactivepeaks)
kidneyw4upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw4upsig & peakanno$trimanno %in% "Intron",]),kidneyactivepeaks)
kidneyw4upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw4upsig & peakanno$trimanno %in% "Distal Intergenic",]),kidneyactivepeaks)
kidneyw4upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyw4upsig & peakanno$trimanno %in% "Upstream",]),kidneyactivepeaks)

kidneyw4downsig <- rownames(kidneysigl2fcmat[kidneysigl2fcmat[,"F W4"] < 0 & kidneysigl2fcmat[,"M W4"] < 0,])
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

liverw8upsig <- rownames(liversigl2fcmat[liversigl2fcmat[,"F W8"] > 0 & liversigl2fcmat[,"M W8"] > 0,])
liverw8upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw8upsig,]),liveractivepeaks)
liverw8upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw8upsig & peakanno$trimanno %in% "Promoter",]),liveractivepeaks)
liverw8upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw8upsig & peakanno$trimanno %in% "Intron",]),liveractivepeaks)
liverw8upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw8upsig & peakanno$trimanno %in% "Distal Intergenic",]),liveractivepeaks)
liverw8upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw8upsig & peakanno$trimanno %in% "Upstream",]),liveractivepeaks)

liverw8downsig <- rownames(liversigl2fcmat[liversigl2fcmat[,"F W8"] < 0 & liversigl2fcmat[,"M W8"] < 0,])
liverw8downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw8downsig,]),liveractivepeaks)
liverw8downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw8downsig & peakanno$trimanno %in% "Promoter",]),liveractivepeaks)
liverw8downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw8downsig & peakanno$trimanno %in% "Intron",]),liveractivepeaks)
liverw8downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw8downsig & peakanno$trimanno %in% "Distal Intergenic",]),liveractivepeaks)
liverw8downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw8downsig & peakanno$trimanno %in% "Upstream",]),liveractivepeaks)

liverw1upsig <- rownames(liversigl2fcmat[liversigl2fcmat[,"F W1"] > 0 & liversigl2fcmat[,"M W1"] > 0,])
liverw1upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw1upsig,]),liveractivepeaks)
liverw1upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw1upsig & peakanno$trimanno %in% "Promoter",]),liveractivepeaks)
liverw1upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw1upsig & peakanno$trimanno %in% "Intron",]),liveractivepeaks)
liverw1upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw1upsig & peakanno$trimanno %in% "Distal Intergenic",]),liveractivepeaks)
liverw1upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw1upsig & peakanno$trimanno %in% "Upstream",]),liveractivepeaks)

liverw1downsig <- rownames(liversigl2fcmat[liversigl2fcmat[,"F W1"] < 0 & liversigl2fcmat[,"M W1"] < 0,])
liverw1downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw1downsig,]),liveractivepeaks)
liverw1downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw1downsig & peakanno$trimanno %in% "Promoter",]),liveractivepeaks)
liverw1downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw1downsig & peakanno$trimanno %in% "Intron",]),liveractivepeaks)
liverw1downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw1downsig & peakanno$trimanno %in% "Distal Intergenic",]),liveractivepeaks)
liverw1downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw1downsig & peakanno$trimanno %in% "Upstream",]),liveractivepeaks)


liverw2upsig <- rownames(liversigl2fcmat[liversigl2fcmat[,"F W2"] > 0 & liversigl2fcmat[,"M W1"] > 0,])
liverw2upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw2upsig,]),liveractivepeaks)
liverw2upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw2upsig & peakanno$trimanno %in% "Promoter",]),liveractivepeaks)
liverw2upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw2upsig & peakanno$trimanno %in% "Intron",]),liveractivepeaks)
liverw2upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw2upsig & peakanno$trimanno %in% "Distal Intergenic",]),liveractivepeaks)
liverw2upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw2upsig & peakanno$trimanno %in% "Upstream",]),liveractivepeaks)

liverw2downsig <- rownames(liversigl2fcmat[liversigl2fcmat[,"F W2"] < 0 & liversigl2fcmat[,"M W1"] < 0,])
liverw2downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw2downsig,]),liveractivepeaks)
liverw2downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw2downsig & peakanno$trimanno %in% "Promoter",]),liveractivepeaks)
liverw2downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw2downsig & peakanno$trimanno %in% "Intron",]),liveractivepeaks)
liverw2downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw2downsig & peakanno$trimanno %in% "Distal Intergenic",]),liveractivepeaks)
liverw2downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw2downsig & peakanno$trimanno %in% "Upstream",]),liveractivepeaks)


liverw4upsig <- rownames(liversigl2fcmat[liversigl2fcmat[,"F W4"] > 0 & liversigl2fcmat[,"M W4"] > 0,])
liverw4upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw4upsig,]),liveractivepeaks)
liverw4upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw4upsig & peakanno$trimanno %in% "Promoter",]),liveractivepeaks)
liverw4upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw4upsig & peakanno$trimanno %in% "Intron",]),liveractivepeaks)
liverw4upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw4upsig & peakanno$trimanno %in% "Distal Intergenic",]),liveractivepeaks)
liverw4upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverw4upsig & peakanno$trimanno %in% "Upstream",]),liveractivepeaks)

liverw4downsig <- rownames(liversigl2fcmat[liversigl2fcmat[,"F W4"] < 0 & liversigl2fcmat[,"M W4"] < 0,])
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

lungw8upsig <- rownames(lungsigl2fcmat[lungsigl2fcmat[,"F W8"] > 0 & lungsigl2fcmat[,"M W8"] > 0,])
lungw8upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw8upsig,]),lungactivepeaks)
lungw8upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw8upsig & peakanno$trimanno %in% "Promoter",]),lungactivepeaks)
lungw8upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw8upsig & peakanno$trimanno %in% "Intron",]),lungactivepeaks)
lungw8upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw8upsig & peakanno$trimanno %in% "Distal Intergenic",]),lungactivepeaks)
lungw8upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw8upsig & peakanno$trimanno %in% "Upstream",]),lungactivepeaks)

lungw8downsig <- rownames(lungsigl2fcmat[lungsigl2fcmat[,"F W8"] < 0 & lungsigl2fcmat[,"M W8"] < 0,])
lungw8downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw8downsig,]),lungactivepeaks)
lungw8downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw8downsig & peakanno$trimanno %in% "Promoter",]),lungactivepeaks)
lungw8downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw8downsig & peakanno$trimanno %in% "Intron",]),lungactivepeaks)
lungw8downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw8downsig & peakanno$trimanno %in% "Distal Intergenic",]),lungactivepeaks)
lungw8downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw8downsig & peakanno$trimanno %in% "Upstream",]),lungactivepeaks)

lungw1upsig <- rownames(lungsigl2fcmat[lungsigl2fcmat[,"F W1"] > 0 & lungsigl2fcmat[,"M W1"] > 0,])
lungw1upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw1upsig,]),lungactivepeaks)
lungw1upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw1upsig & peakanno$trimanno %in% "Promoter",]),lungactivepeaks)
lungw1upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw1upsig & peakanno$trimanno %in% "Intron",]),lungactivepeaks)
lungw1upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw1upsig & peakanno$trimanno %in% "Distal Intergenic",]),lungactivepeaks)
lungw1upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw1upsig & peakanno$trimanno %in% "Upstream",]),lungactivepeaks)

lungw1downsig <- rownames(lungsigl2fcmat[lungsigl2fcmat[,"F W1"] < 0 & lungsigl2fcmat[,"M W1"] < 0,])
lungw1downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw1downsig,]),lungactivepeaks)
lungw1downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw1downsig & peakanno$trimanno %in% "Promoter",]),lungactivepeaks)
lungw1downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw1downsig & peakanno$trimanno %in% "Intron",]),lungactivepeaks)
lungw1downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw1downsig & peakanno$trimanno %in% "Distal Intergenic",]),lungactivepeaks)
lungw1downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw1downsig & peakanno$trimanno %in% "Upstream",]),lungactivepeaks)


lungw2upsig <- rownames(lungsigl2fcmat[lungsigl2fcmat[,"F W2"] > 0 & lungsigl2fcmat[,"M W1"] > 0,])
lungw2upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw2upsig,]),lungactivepeaks)
lungw2upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw2upsig & peakanno$trimanno %in% "Promoter",]),lungactivepeaks)
lungw2upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw2upsig & peakanno$trimanno %in% "Intron",]),lungactivepeaks)
lungw2upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw2upsig & peakanno$trimanno %in% "Distal Intergenic",]),lungactivepeaks)
lungw2upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw2upsig & peakanno$trimanno %in% "Upstream",]),lungactivepeaks)

lungw2downsig <- rownames(lungsigl2fcmat[lungsigl2fcmat[,"F W2"] < 0 & lungsigl2fcmat[,"M W1"] < 0,])
lungw2downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw2downsig,]),lungactivepeaks)
lungw2downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw2downsig & peakanno$trimanno %in% "Promoter",]),lungactivepeaks)
lungw2downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw2downsig & peakanno$trimanno %in% "Intron",]),lungactivepeaks)
lungw2downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw2downsig & peakanno$trimanno %in% "Distal Intergenic",]),lungactivepeaks)
lungw2downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw2downsig & peakanno$trimanno %in% "Upstream",]),lungactivepeaks)


lungw4upsig <- rownames(lungsigl2fcmat[lungsigl2fcmat[,"F W4"] > 0 & lungsigl2fcmat[,"M W4"] > 0,])
lungw4upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw4upsig,]),lungactivepeaks)
lungw4upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw4upsig & peakanno$trimanno %in% "Promoter",]),lungactivepeaks)
lungw4upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw4upsig & peakanno$trimanno %in% "Intron",]),lungactivepeaks)
lungw4upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw4upsig & peakanno$trimanno %in% "Distal Intergenic",]),lungactivepeaks)
lungw4upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungw4upsig & peakanno$trimanno %in% "Upstream",]),lungactivepeaks)

lungw4downsig <- rownames(lungsigl2fcmat[lungsigl2fcmat[,"F W4"] < 0 & lungsigl2fcmat[,"M W4"] < 0,])
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

brownw8upsig <- rownames(brownsigl2fcmat[brownsigl2fcmat[,"F W8"] > 0 & brownsigl2fcmat[,"M W8"] > 0,])
brownw8upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw8upsig,]),brownactivepeaks)
brownw8upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw8upsig & peakanno$trimanno %in% "Promoter",]),brownactivepeaks)
brownw8upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw8upsig & peakanno$trimanno %in% "Intron",]),brownactivepeaks)
brownw8upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw8upsig & peakanno$trimanno %in% "Distal Intergenic",]),brownactivepeaks)
brownw8upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw8upsig & peakanno$trimanno %in% "Upstream",]),brownactivepeaks)

brownw8downsig <- rownames(brownsigl2fcmat[brownsigl2fcmat[,"F W8"] < 0 & brownsigl2fcmat[,"M W8"] < 0,])
brownw8downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw8downsig,]),brownactivepeaks)
brownw8downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw8downsig & peakanno$trimanno %in% "Promoter",]),brownactivepeaks)
brownw8downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw8downsig & peakanno$trimanno %in% "Intron",]),brownactivepeaks)
brownw8downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw8downsig & peakanno$trimanno %in% "Distal Intergenic",]),brownactivepeaks)
brownw8downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw8downsig & peakanno$trimanno %in% "Upstream",]),brownactivepeaks)

brownw1upsig <- rownames(brownsigl2fcmat[brownsigl2fcmat[,"F W1"] > 0 & brownsigl2fcmat[,"M W1"] > 0,])
brownw1upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw1upsig,]),brownactivepeaks)
brownw1upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw1upsig & peakanno$trimanno %in% "Promoter",]),brownactivepeaks)
brownw1upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw1upsig & peakanno$trimanno %in% "Intron",]),brownactivepeaks)
brownw1upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw1upsig & peakanno$trimanno %in% "Distal Intergenic",]),brownactivepeaks)
brownw1upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw1upsig & peakanno$trimanno %in% "Upstream",]),brownactivepeaks)

brownw1downsig <- rownames(brownsigl2fcmat[brownsigl2fcmat[,"F W1"] < 0 & brownsigl2fcmat[,"M W1"] < 0,])
brownw1downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw1downsig,]),brownactivepeaks)
brownw1downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw1downsig & peakanno$trimanno %in% "Promoter",]),brownactivepeaks)
brownw1downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw1downsig & peakanno$trimanno %in% "Intron",]),brownactivepeaks)
brownw1downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw1downsig & peakanno$trimanno %in% "Distal Intergenic",]),brownactivepeaks)
brownw1downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw1downsig & peakanno$trimanno %in% "Upstream",]),brownactivepeaks)


brownw2upsig <- rownames(brownsigl2fcmat[brownsigl2fcmat[,"F W2"] > 0 & brownsigl2fcmat[,"M W1"] > 0,])
brownw2upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw2upsig,]),brownactivepeaks)
brownw2upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw2upsig & peakanno$trimanno %in% "Promoter",]),brownactivepeaks)
brownw2upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw2upsig & peakanno$trimanno %in% "Intron",]),brownactivepeaks)
brownw2upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw2upsig & peakanno$trimanno %in% "Distal Intergenic",]),brownactivepeaks)
brownw2upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw2upsig & peakanno$trimanno %in% "Upstream",]),brownactivepeaks)

brownw2downsig <- rownames(brownsigl2fcmat[brownsigl2fcmat[,"F W2"] < 0 & brownsigl2fcmat[,"M W1"] < 0,])
brownw2downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw2downsig,]),brownactivepeaks)
brownw2downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw2downsig & peakanno$trimanno %in% "Promoter",]),brownactivepeaks)
brownw2downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw2downsig & peakanno$trimanno %in% "Intron",]),brownactivepeaks)
brownw2downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw2downsig & peakanno$trimanno %in% "Distal Intergenic",]),brownactivepeaks)
brownw2downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw2downsig & peakanno$trimanno %in% "Upstream",]),brownactivepeaks)

brownw4upsig <- rownames(brownsigl2fcmat[brownsigl2fcmat[,"F W4"] > 0 & brownsigl2fcmat[,"M W4"] > 0,])
brownw4upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw4upsig,]),brownactivepeaks)
brownw4upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw4upsig & peakanno$trimanno %in% "Promoter",]),brownactivepeaks)
brownw4upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw4upsig & peakanno$trimanno %in% "Intron",]),brownactivepeaks)
brownw4upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw4upsig & peakanno$trimanno %in% "Distal Intergenic",]),brownactivepeaks)
brownw4upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownw4upsig & peakanno$trimanno %in% "Upstream",]),brownactivepeaks)

brownw4downsig <- rownames(brownsigl2fcmat[brownsigl2fcmat[,"F W4"] < 0 & brownsigl2fcmat[,"M W4"] < 0,])
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

whitew8upsig <- rownames(whitesigl2fcmat[whitesigl2fcmat[,"F W8"] > 0 & whitesigl2fcmat[,"M W8"] > 0,])
whitew8upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew8upsig,]),whiteactivepeaks)
whitew8upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew8upsig & peakanno$trimanno %in% "Promoter",]),whiteactivepeaks)
whitew8upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew8upsig & peakanno$trimanno %in% "Intron",]),whiteactivepeaks)
whitew8upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew8upsig & peakanno$trimanno %in% "Distal Intergenic",]),whiteactivepeaks)
whitew8upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew8upsig & peakanno$trimanno %in% "Upstream",]),whiteactivepeaks)

whitew8downsig <- rownames(whitesigl2fcmat[whitesigl2fcmat[,"F W8"] < 0 & whitesigl2fcmat[,"M W8"] < 0,])
whitew8downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew8downsig,]),whiteactivepeaks)
whitew8downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew8downsig & peakanno$trimanno %in% "Promoter",]),whiteactivepeaks)
whitew8downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew8downsig & peakanno$trimanno %in% "Intron",]),whiteactivepeaks)
whitew8downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew8downsig & peakanno$trimanno %in% "Distal Intergenic",]),whiteactivepeaks)
whitew8downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew8downsig & peakanno$trimanno %in% "Upstream",]),whiteactivepeaks)

whitew1upsig <- rownames(whitesigl2fcmat[whitesigl2fcmat[,"F W1"] > 0 & whitesigl2fcmat[,"M W1"] > 0,])
whitew1upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew1upsig,]),whiteactivepeaks)
whitew1upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew1upsig & peakanno$trimanno %in% "Promoter",]),whiteactivepeaks)
whitew1upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew1upsig & peakanno$trimanno %in% "Intron",]),whiteactivepeaks)
whitew1upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew1upsig & peakanno$trimanno %in% "Distal Intergenic",]),whiteactivepeaks)
whitew1upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew1upsig & peakanno$trimanno %in% "Upstream",]),whiteactivepeaks)

whitew1downsig <- rownames(whitesigl2fcmat[whitesigl2fcmat[,"F W1"] < 0 & whitesigl2fcmat[,"M W1"] < 0,])
whitew1downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew1downsig,]),whiteactivepeaks)
whitew1downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew1downsig & peakanno$trimanno %in% "Promoter",]),whiteactivepeaks)
whitew1downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew1downsig & peakanno$trimanno %in% "Intron",]),whiteactivepeaks)
whitew1downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew1downsig & peakanno$trimanno %in% "Distal Intergenic",]),whiteactivepeaks)
whitew1downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew1downsig & peakanno$trimanno %in% "Upstream",]),whiteactivepeaks)


whitew2upsig <- rownames(whitesigl2fcmat[whitesigl2fcmat[,"F W2"] > 0 & whitesigl2fcmat[,"M W1"] > 0,])
whitew2upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew2upsig,]),whiteactivepeaks)
whitew2upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew2upsig & peakanno$trimanno %in% "Promoter",]),whiteactivepeaks)
whitew2upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew2upsig & peakanno$trimanno %in% "Intron",]),whiteactivepeaks)
whitew2upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew2upsig & peakanno$trimanno %in% "Distal Intergenic",]),whiteactivepeaks)
whitew2upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew2upsig & peakanno$trimanno %in% "Upstream",]),whiteactivepeaks)

whitew2downsig <- rownames(whitesigl2fcmat[whitesigl2fcmat[,"F W2"] < 0 & whitesigl2fcmat[,"M W1"] < 0,])
whitew2downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew2downsig,]),whiteactivepeaks)
whitew2downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew2downsig & peakanno$trimanno %in% "Promoter",]),whiteactivepeaks)
whitew2downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew2downsig & peakanno$trimanno %in% "Intron",]),whiteactivepeaks)
whitew2downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew2downsig & peakanno$trimanno %in% "Distal Intergenic",]),whiteactivepeaks)
whitew2downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew2downsig & peakanno$trimanno %in% "Upstream",]),whiteactivepeaks)


whitew4upsig <- rownames(whitesigl2fcmat[whitesigl2fcmat[,"F W4"] > 0 & whitesigl2fcmat[,"M W4"] > 0,])
whitew4upsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew4upsig,]),whiteactivepeaks)
whitew4upsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew4upsig & peakanno$trimanno %in% "Promoter",]),whiteactivepeaks)
whitew4upsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew4upsig & peakanno$trimanno %in% "Intron",]),whiteactivepeaks)
whitew4upsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew4upsig & peakanno$trimanno %in% "Distal Intergenic",]),whiteactivepeaks)
whitew4upsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew4upsig & peakanno$trimanno %in% "Upstream",]),whiteactivepeaks)

whitew4downsig <- rownames(whitesigl2fcmat[whitesigl2fcmat[,"F W4"] < 0 & whitesigl2fcmat[,"M W4"] < 0,])
whitew4downsigpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew4downsig,]),whiteactivepeaks)
whitew4downsigprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew4downsig & peakanno$trimanno %in% "Promoter",]),whiteactivepeaks)
whitew4downsigintpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew4downsig & peakanno$trimanno %in% "Intron",]),whiteactivepeaks)
whitew4downsigdistpeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew4downsig & peakanno$trimanno %in% "Distal Intergenic",]),whiteactivepeaks)
whitew4downsiguppeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% whitew4downsig & peakanno$trimanno %in% "Upstream",]),whiteactivepeaks)

####
# Table generation
#####

# Making Files of all active peaks in each tissue
gastroatacactiveanno <- peakanno[gastroactivepeaks,c("chrom","start","end","geneStrand")]
gastroatacactiveanno$chrom <- paste("chr",gastroatacactiveanno$chrom,sep = "")
write.table(gastroatacactiveanno,file = "gastroatacactiveanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

heartatacactiveanno <- peakanno[heartactivepeaks,c("chrom","start","end","geneStrand")]
heartatacactiveanno$chrom <- paste("chr",heartatacactiveanno$chrom,sep = "")
write.table(heartatacactiveanno,file = "heartatacactiveanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

hippoatacactiveanno <- peakanno[hippoactivepeaks,c("chrom","start","end","geneStrand")]
hippoatacactiveanno$chrom <- paste("chr",hippoatacactiveanno$chrom,sep = "")
write.table(hippoatacactiveanno,file = "hippoatacactiveanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

kidneyatacactiveanno <- peakanno[kidneyactivepeaks,c("chrom","start","end","geneStrand")]
kidneyatacactiveanno$chrom <- paste("chr",kidneyatacactiveanno$chrom,sep = "")
write.table(kidneyatacactiveanno,file = "kidneyatacactiveanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

liveratacactiveanno <- peakanno[liveractivepeaks,c("chrom","start","end","geneStrand")]
liveratacactiveanno$chrom <- paste("chr",liveratacactiveanno$chrom,sep = "")
write.table(liveratacactiveanno,file = "liveratacactiveanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

lungatacactiveanno <- peakanno[lungactivepeaks,c("chrom","start","end","geneStrand")]
lungatacactiveanno$chrom <- paste("chr",lungatacactiveanno$chrom,sep = "")
write.table(lungatacactiveanno,file = "lungatacactiveanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

brownatacactiveanno <- peakanno[brownactivepeaks,c("chrom","start","end","geneStrand")]
brownatacactiveanno$chrom <- paste("chr",brownatacactiveanno$chrom,sep = "")
write.table(brownatacactiveanno,file = "brownatacactiveanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

whiteatacactiveanno <- peakanno[whiteactivepeaks,c("chrom","start","end","geneStrand")]
whiteatacactiveanno$chrom <- paste("chr",whiteatacactiveanno$chrom,sep = "")
write.table(whiteatacactiveanno,file = "whiteatacactiveanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

# Making DAR narrowpeak files

gastroatacsiganno <- peakanno[gastroatacsig,c("chrom","start","end","geneStrand")]
gastroatacsiganno$chrom <- paste("chr",gastroatacsiganno$chrom,sep = "")
write.table(gastroatacsiganno,file = "gastroatacsiganno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

heartatacsiganno <- peakanno[heartatacsig,c("chrom","start","end","geneStrand")]
heartatacsiganno$chrom <- paste("chr",heartatacsiganno$chrom,sep = "")
write.table(heartatacsiganno,file = "heartatacsiganno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

kidneyatacsiganno <- peakanno[kidneyatacsig,c("chrom","start","end","geneStrand")]
kidneyatacsiganno$chrom <- paste("chr",kidneyatacsiganno$chrom,sep = "")
write.table(kidneyatacsiganno,file = "kidneyatacsiganno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

liveratacsiganno <- peakanno[liveratacsig,c("chrom","start","end","geneStrand")]
liveratacsiganno$chrom <- paste("chr",liveratacsiganno$chrom,sep = "")
write.table(liveratacsiganno,file = "liveratacsiganno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

lungatacsiganno <- peakanno[lungatacsig,c("chrom","start","end","geneStrand")]
lungatacsiganno$chrom <- paste("chr",lungatacsiganno$chrom,sep = "")
write.table(lungatacsiganno,file = "lungatacsiganno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

brownatacsiganno <- peakanno[brownatacsig,c("chrom","start","end","geneStrand")]
brownatacsiganno$chrom <- paste("chr",brownatacsiganno$chrom,sep = "")
write.table(brownatacsiganno,file = "brownatacsiganno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

# Making DEGaP narrowpeak files

# SKM-GN
gastrodegpeakanno <- peakanno[gastrosigpeak,c("chrom","start","end","geneStrand")]
gastrodegpeakanno$chrom <- paste("chr",gastrodegpeakanno$chrom,sep = "")
write.table(gastrodegpeakanno,file = "gastrodegpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

gastrodegprompeakanno <- peakanno[gastrosigprompeak,c("chrom","start","end","geneStrand")]
gastrodegprompeakanno$chrom <- paste("chr",gastrodegprompeakanno$chrom,sep = "")
write.table(gastrodegprompeakanno,file = "gastrodegprompeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

gastrodegpromproxpeakanno <- peakanno[gastrosigpromproxpeak,c("chrom","start","end","geneStrand")]
gastrodegpromproxpeakanno$chrom <- paste("chr",gastrodegpromproxpeakanno$chrom,sep = "")
write.table(gastrodegpromproxpeakanno,file = "gastrodegpromproxpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

gastrodegpromfarpeakanno <- peakanno[gastrosigpromfarpeak,c("chrom","start","end","geneStrand")]
gastrodegpromfarpeakanno$chrom <- paste("chr",gastrodegpromfarpeakanno$chrom,sep = "")
write.table(gastrodegpromfarpeakanno,file = "gastrodegpromfarpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

gastrodegdistpeakanno <- peakanno[gastrosigdistpeak,c("chrom","start","end","geneStrand")]
gastrodegdistpeakanno$chrom <- paste("chr",gastrodegdistpeakanno$chrom,sep = "")
write.table(gastrodegdistpeakanno,file = "gastrodegdistpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

gastrodegintpeakanno <- peakanno[gastrosigintpeak,c("chrom","start","end","geneStrand")]
gastrodegintpeakanno$chrom <- paste("chr",gastrodegintpeakanno$chrom,sep = "")
write.table(gastrodegintpeakanno,file = "gastrodegintpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

gastrodegexpeakanno <- peakanno[gastrosigexpeak,c("chrom","start","end","geneStrand")]
gastrodegexpeakanno$chrom <- paste("chr",gastrodegexpeakanno$chrom,sep = "")
write.table(gastrodegexpeakanno,file = "gastrodegexpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

gastrodeguppeakanno <- peakanno[gastrosiguppeak,c("chrom","start","end","geneStrand")]
gastrodeguppeakanno$chrom <- paste("chr",gastrodeguppeakanno$chrom,sep = "")
write.table(gastrodeguppeakanno,file = "gastrodeguppeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

gastrodegdownpeakanno <- peakanno[gastrosigdownpeak,c("chrom","start","end","geneStrand")]
gastrodegdownpeakanno$chrom <- paste("chr",gastrodegdownpeakanno$chrom,sep = "")
write.table(gastrodegdownpeakanno,file = "gastrodegdownpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

# HEART
heartdegpeakanno <- peakanno[heartsigpeak,c("chrom","start","end","geneStrand")]
heartdegpeakanno$chrom <- paste("chr",heartdegpeakanno$chrom,sep = "")
write.table(heartdegpeakanno,file = "heartdegpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

heartdegprompeakanno <- peakanno[heartsigprompeak,c("chrom","start","end","geneStrand")]
heartdegprompeakanno$chrom <- paste("chr",heartdegprompeakanno$chrom,sep = "")
write.table(heartdegprompeakanno,file = "heartdegprompeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

heartdegpromproxpeakanno <- peakanno[heartsigpromproxpeak,c("chrom","start","end","geneStrand")]
heartdegpromproxpeakanno$chrom <- paste("chr",heartdegpromproxpeakanno$chrom,sep = "")
write.table(heartdegpromproxpeakanno,file = "heartdegpromproxpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

heartdegpromfarpeakanno <- peakanno[heartsigpromfarpeak,c("chrom","start","end","geneStrand")]
heartdegpromfarpeakanno$chrom <- paste("chr",heartdegpromfarpeakanno$chrom,sep = "")
write.table(heartdegpromfarpeakanno,file = "heartdegpromfarpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

heartdegdistpeakanno <- peakanno[heartsigdistpeak,c("chrom","start","end","geneStrand")]
heartdegdistpeakanno$chrom <- paste("chr",heartdegdistpeakanno$chrom,sep = "")
write.table(heartdegdistpeakanno,file = "heartdegdistpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

heartdegintpeakanno <- peakanno[heartsigintpeak,c("chrom","start","end","geneStrand")]
heartdegintpeakanno$chrom <- paste("chr",heartdegintpeakanno$chrom,sep = "")
write.table(heartdegintpeakanno,file = "heartdegintpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

heartdegexpeakanno <- peakanno[heartsigexpeak,c("chrom","start","end","geneStrand")]
heartdegexpeakanno$chrom <- paste("chr",heartdegexpeakanno$chrom,sep = "")
write.table(heartdegexpeakanno,file = "heartdegexpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

heartdeguppeakanno <- peakanno[heartsiguppeak,c("chrom","start","end","geneStrand")]
heartdeguppeakanno$chrom <- paste("chr",heartdeguppeakanno$chrom,sep = "")
write.table(heartdeguppeakanno,file = "heartdeguppeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

heartdegdownpeakanno <- peakanno[heartsigdownpeak,c("chrom","start","end","geneStrand")]
heartdegdownpeakanno$chrom <- paste("chr",heartdegdownpeakanno$chrom,sep = "")
write.table(heartdegdownpeakanno,file = "heartdegdownpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")


# HIPPOC
hippodegpeakanno <- peakanno[hipposigpeak,c("chrom","start","end","geneStrand")]
hippodegpeakanno$chrom <- paste("chr",hippodegpeakanno$chrom,sep = "")
write.table(hippodegpeakanno,file = "hippodegpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

hippodegprompeakanno <- peakanno[hipposigprompeak,c("chrom","start","end","geneStrand")]
hippodegprompeakanno$chrom <- paste("chr",hippodegprompeakanno$chrom,sep = "")
write.table(hippodegprompeakanno,file = "hippodegprompeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

hippodegpromproxpeakanno <- peakanno[hipposigpromproxpeak,c("chrom","start","end","geneStrand")]
hippodegpromproxpeakanno$chrom <- paste("chr",hippodegpromproxpeakanno$chrom,sep = "")
write.table(hippodegpromproxpeakanno,file = "hippodegpromproxpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

hippodegpromfarpeakanno <- peakanno[hipposigpromfarpeak,c("chrom","start","end","geneStrand")]
hippodegpromfarpeakanno$chrom <- paste("chr",hippodegpromfarpeakanno$chrom,sep = "")
write.table(hippodegpromfarpeakanno,file = "hippodegpromfarpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

hippodegdistpeakanno <- peakanno[hipposigdistpeak,c("chrom","start","end","geneStrand")]
hippodegdistpeakanno$chrom <- paste("chr",hippodegdistpeakanno$chrom,sep = "")
write.table(hippodegdistpeakanno,file = "hippodegdistpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

hippodegintpeakanno <- peakanno[hipposigintpeak,c("chrom","start","end","geneStrand")]
hippodegintpeakanno$chrom <- paste("chr",hippodegintpeakanno$chrom,sep = "")
write.table(hippodegintpeakanno,file = "hippodegintpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

hippodegexpeakanno <- peakanno[hipposigexpeak,c("chrom","start","end","geneStrand")]
hippodegexpeakanno$chrom <- paste("chr",hippodegexpeakanno$chrom,sep = "")
write.table(hippodegexpeakanno,file = "hippodegexpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

hippodeguppeakanno <- peakanno[hipposiguppeak,c("chrom","start","end","geneStrand")]
hippodeguppeakanno$chrom <- paste("chr",hippodeguppeakanno$chrom,sep = "")
write.table(hippodeguppeakanno,file = "hippodeguppeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

hippodegdownpeakanno <- peakanno[hipposigdownpeak,c("chrom","start","end","geneStrand")]
hippodegdownpeakanno$chrom <- paste("chr",hippodegdownpeakanno$chrom,sep = "")
write.table(hippodegdownpeakanno,file = "hippodegdownpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

# KIDNEY
kidneydegpeakanno <- peakanno[kidneysigpeak,c("chrom","start","end","geneStrand")]
kidneydegpeakanno$chrom <- paste("chr",kidneydegpeakanno$chrom,sep = "")
write.table(kidneydegpeakanno,file = "kidneydegpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

kidneydegprompeakanno <- peakanno[kidneysigprompeak,c("chrom","start","end","geneStrand")]
kidneydegprompeakanno$chrom <- paste("chr",kidneydegprompeakanno$chrom,sep = "")
write.table(kidneydegprompeakanno,file = "kidneydegprompeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

kidneydegpromproxpeakanno <- peakanno[kidneysigpromproxpeak,c("chrom","start","end","geneStrand")]
kidneydegpromproxpeakanno$chrom <- paste("chr",kidneydegpromproxpeakanno$chrom,sep = "")
write.table(kidneydegpromproxpeakanno,file = "kidneydegpromproxpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

kidneydegpromfarpeakanno <- peakanno[kidneysigpromfarpeak,c("chrom","start","end","geneStrand")]
kidneydegpromfarpeakanno$chrom <- paste("chr",kidneydegpromfarpeakanno$chrom,sep = "")
write.table(kidneydegpromfarpeakanno,file = "kidneydegpromfarpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

kidneydegdistpeakanno <- peakanno[kidneysigdistpeak,c("chrom","start","end","geneStrand")]
kidneydegdistpeakanno$chrom <- paste("chr",kidneydegdistpeakanno$chrom,sep = "")
write.table(kidneydegdistpeakanno,file = "kidneydegdistpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

kidneydegintpeakanno <- peakanno[kidneysigintpeak,c("chrom","start","end","geneStrand")]
kidneydegintpeakanno$chrom <- paste("chr",kidneydegintpeakanno$chrom,sep = "")
write.table(kidneydegintpeakanno,file = "kidneydegintpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

kidneydegexpeakanno <- peakanno[kidneysigexpeak,c("chrom","start","end","geneStrand")]
kidneydegexpeakanno$chrom <- paste("chr",kidneydegexpeakanno$chrom,sep = "")
write.table(kidneydegexpeakanno,file = "kidneydegexpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

kidneydeguppeakanno <- peakanno[kidneysiguppeak,c("chrom","start","end","geneStrand")]
kidneydeguppeakanno$chrom <- paste("chr",kidneydeguppeakanno$chrom,sep = "")
write.table(kidneydeguppeakanno,file = "kidneydeguppeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

kidneydegdownpeakanno <- peakanno[kidneysigdownpeak,c("chrom","start","end","geneStrand")]
kidneydegdownpeakanno$chrom <- paste("chr",kidneydegdownpeakanno$chrom,sep = "")
write.table(kidneydegdownpeakanno,file = "kidneydegdownpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

# LIVER
liverdegpeakanno <- peakanno[liversigpeak,c("chrom","start","end","geneStrand")]
liverdegpeakanno$chrom <- paste("chr",liverdegpeakanno$chrom,sep = "")
write.table(liverdegpeakanno,file = "liverdegpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

liverdegprompeakanno <- peakanno[liversigprompeak,c("chrom","start","end","geneStrand")]
liverdegprompeakanno$chrom <- paste("chr",liverdegprompeakanno$chrom,sep = "")
write.table(liverdegprompeakanno,file = "liverdegprompeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

liverdegpromproxpeakanno <- peakanno[liversigpromproxpeak,c("chrom","start","end","geneStrand")]
liverdegpromproxpeakanno$chrom <- paste("chr",liverdegpromproxpeakanno$chrom,sep = "")
write.table(liverdegpromproxpeakanno,file = "liverdegpromproxpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

liverdegpromfarpeakanno <- peakanno[liversigpromfarpeak,c("chrom","start","end","geneStrand")]
liverdegpromfarpeakanno$chrom <- paste("chr",liverdegpromfarpeakanno$chrom,sep = "")
write.table(liverdegpromfarpeakanno,file = "liverdegpromfarpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

liverdegdistpeakanno <- peakanno[liversigdistpeak,c("chrom","start","end","geneStrand")]
liverdegdistpeakanno$chrom <- paste("chr",liverdegdistpeakanno$chrom,sep = "")
write.table(liverdegdistpeakanno,file = "liverdegdistpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

liverdegintpeakanno <- peakanno[liversigintpeak,c("chrom","start","end","geneStrand")]
liverdegintpeakanno$chrom <- paste("chr",liverdegintpeakanno$chrom,sep = "")
write.table(liverdegintpeakanno,file = "liverdegintpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

liverdegexpeakanno <- peakanno[liversigexpeak,c("chrom","start","end","geneStrand")]
liverdegexpeakanno$chrom <- paste("chr",liverdegexpeakanno$chrom,sep = "")
write.table(liverdegexpeakanno,file = "liverdegexpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

liverdeguppeakanno <- peakanno[liversiguppeak,c("chrom","start","end","geneStrand")]
liverdeguppeakanno$chrom <- paste("chr",liverdeguppeakanno$chrom,sep = "")
write.table(liverdeguppeakanno,file = "liverdeguppeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

liverdegdownpeakanno <- peakanno[liversigdownpeak,c("chrom","start","end","geneStrand")]
liverdegdownpeakanno$chrom <- paste("chr",liverdegdownpeakanno$chrom,sep = "")
write.table(liverdegdownpeakanno,file = "liverdegdownpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

# LUNG
lungdegpeakanno <- peakanno[lungsigpeak,c("chrom","start","end","geneStrand")]
lungdegpeakanno$chrom <- paste("chr",lungdegpeakanno$chrom,sep = "")
write.table(lungdegpeakanno,file = "lungdegpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

lungdegprompeakanno <- peakanno[lungsigprompeak,c("chrom","start","end","geneStrand")]
lungdegprompeakanno$chrom <- paste("chr",lungdegprompeakanno$chrom,sep = "")
write.table(lungdegprompeakanno,file = "lungdegprompeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

lungdegpromproxpeakanno <- peakanno[lungsigpromproxpeak,c("chrom","start","end","geneStrand")]
lungdegpromproxpeakanno$chrom <- paste("chr",lungdegpromproxpeakanno$chrom,sep = "")
write.table(lungdegpromproxpeakanno,file = "lungdegpromproxpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

lungdegpromfarpeakanno <- peakanno[lungsigpromfarpeak,c("chrom","start","end","geneStrand")]
lungdegpromfarpeakanno$chrom <- paste("chr",lungdegpromfarpeakanno$chrom,sep = "")
write.table(lungdegpromfarpeakanno,file = "lungdegpromfarpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

lungdegdistpeakanno <- peakanno[lungsigdistpeak,c("chrom","start","end","geneStrand")]
lungdegdistpeakanno$chrom <- paste("chr",lungdegdistpeakanno$chrom,sep = "")
write.table(lungdegdistpeakanno,file = "lungdegdistpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

lungdegintpeakanno <- peakanno[lungsigintpeak,c("chrom","start","end","geneStrand")]
lungdegintpeakanno$chrom <- paste("chr",lungdegintpeakanno$chrom,sep = "")
write.table(lungdegintpeakanno,file = "lungdegintpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

lungdegexpeakanno <- peakanno[lungsigexpeak,c("chrom","start","end","geneStrand")]
lungdegexpeakanno$chrom <- paste("chr",lungdegexpeakanno$chrom,sep = "")
write.table(lungdegexpeakanno,file = "lungdegexpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

lungdeguppeakanno <- peakanno[lungsiguppeak,c("chrom","start","end","geneStrand")]
lungdeguppeakanno$chrom <- paste("chr",lungdeguppeakanno$chrom,sep = "")
write.table(lungdeguppeakanno,file = "lungdeguppeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

lungdegdownpeakanno <- peakanno[lungsigdownpeak,c("chrom","start","end","geneStrand")]
lungdegdownpeakanno$chrom <- paste("chr",lungdegdownpeakanno$chrom,sep = "")
write.table(lungdegdownpeakanno,file = "lungdegdownpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

# BAT
browndegpeakanno <- peakanno[brownsigpeak,c("chrom","start","end","geneStrand")]
browndegpeakanno$chrom <- paste("chr",browndegpeakanno$chrom,sep = "")
write.table(browndegpeakanno,file = "browndegpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

browndegprompeakanno <- peakanno[brownsigprompeak,c("chrom","start","end","geneStrand")]
browndegprompeakanno$chrom <- paste("chr",browndegprompeakanno$chrom,sep = "")
write.table(browndegprompeakanno,file = "browndegprompeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

browndegpromproxpeakanno <- peakanno[brownsigpromproxpeak,c("chrom","start","end","geneStrand")]
browndegpromproxpeakanno$chrom <- paste("chr",browndegpromproxpeakanno$chrom,sep = "")
write.table(browndegpromproxpeakanno,file = "browndegpromproxpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

browndegpromfarpeakanno <- peakanno[brownsigpromfarpeak,c("chrom","start","end","geneStrand")]
browndegpromfarpeakanno$chrom <- paste("chr",browndegpromfarpeakanno$chrom,sep = "")
write.table(browndegpromfarpeakanno,file = "browndegpromfarpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

browndegdistpeakanno <- peakanno[brownsigdistpeak,c("chrom","start","end","geneStrand")]
browndegdistpeakanno$chrom <- paste("chr",browndegdistpeakanno$chrom,sep = "")
write.table(browndegdistpeakanno,file = "browndegdistpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

browndegintpeakanno <- peakanno[brownsigintpeak,c("chrom","start","end","geneStrand")]
browndegintpeakanno$chrom <- paste("chr",browndegintpeakanno$chrom,sep = "")
write.table(browndegintpeakanno,file = "browndegintpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

browndegexpeakanno <- peakanno[brownsigexpeak,c("chrom","start","end","geneStrand")]
browndegexpeakanno$chrom <- paste("chr",browndegexpeakanno$chrom,sep = "")
write.table(browndegexpeakanno,file = "browndegexpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

browndeguppeakanno <- peakanno[brownsiguppeak,c("chrom","start","end","geneStrand")]
browndeguppeakanno$chrom <- paste("chr",browndeguppeakanno$chrom,sep = "")
write.table(browndeguppeakanno,file = "browndeguppeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

browndegdownpeakanno <- peakanno[brownsigdownpeak,c("chrom","start","end","geneStrand")]
browndegdownpeakanno$chrom <- paste("chr",browndegdownpeakanno$chrom,sep = "")
write.table(browndegdownpeakanno,file = "browndegdownpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

# WAT-SC
whitedegpeakanno <- peakanno[whitesigpeak,c("chrom","start","end","geneStrand")]
whitedegpeakanno$chrom <- paste("chr",whitedegpeakanno$chrom,sep = "")
write.table(whitedegpeakanno,file = "whitedegpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

whitedegprompeakanno <- peakanno[whitesigprompeak,c("chrom","start","end","geneStrand")]
whitedegprompeakanno$chrom <- paste("chr",whitedegprompeakanno$chrom,sep = "")
write.table(whitedegprompeakanno,file = "whitedegprompeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

whitedegpromproxpeakanno <- peakanno[whitesigpromproxpeak,c("chrom","start","end","geneStrand")]
whitedegpromproxpeakanno$chrom <- paste("chr",whitedegpromproxpeakanno$chrom,sep = "")
write.table(whitedegpromproxpeakanno,file = "whitedegpromproxpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

whitedegpromfarpeakanno <- peakanno[whitesigpromfarpeak,c("chrom","start","end","geneStrand")]
whitedegpromfarpeakanno$chrom <- paste("chr",whitedegpromfarpeakanno$chrom,sep = "")
write.table(whitedegpromfarpeakanno,file = "whitedegpromfarpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

whitedegdistpeakanno <- peakanno[whitesigdistpeak,c("chrom","start","end","geneStrand")]
whitedegdistpeakanno$chrom <- paste("chr",whitedegdistpeakanno$chrom,sep = "")
write.table(whitedegdistpeakanno,file = "whitedegdistpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

whitedegintpeakanno <- peakanno[whitesigintpeak,c("chrom","start","end","geneStrand")]
whitedegintpeakanno$chrom <- paste("chr",whitedegintpeakanno$chrom,sep = "")
write.table(whitedegintpeakanno,file = "whitedegintpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

whitedegexpeakanno <- peakanno[whitesigexpeak,c("chrom","start","end","geneStrand")]
whitedegexpeakanno$chrom <- paste("chr",whitedegexpeakanno$chrom,sep = "")
write.table(whitedegexpeakanno,file = "whitedegexpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

whitedeguppeakanno <- peakanno[whitesiguppeak,c("chrom","start","end","geneStrand")]
whitedeguppeakanno$chrom <- paste("chr",whitedeguppeakanno$chrom,sep = "")
write.table(whitedeguppeakanno,file = "whitedeguppeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

whitedegdownpeakanno <- peakanno[whitesigdownpeak,c("chrom","start","end","geneStrand")]
whitedegdownpeakanno$chrom <- paste("chr",whitedegdownpeakanno$chrom,sep = "")
write.table(whitedegdownpeakanno,file = "whitedegdownpeakanno.narrowPeak",col.names = F,row.names = T,quote = F,sep = "\t")

####
# Use tables as input for homer findMotifsGenome.pl command with flags -size 50 -p 10 -preparse
####

####
## Loading the results of Homer TF enrichment analysis on gene region-specific DEGaPs
#####

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

# Homer TF enrichment among DARs in each tissue
gastroatacsigtf50 <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/gastroatacsigfilejun2022out/knownResults.csv",header = T,row.names = 1)
heartatacsigtf50 <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/heartatacsigfilejun2022out/knownResults.csv",header = T,row.names = 1)
kidneyatacsigtf50 <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/kidneyatacsigfilejun2022out/knownResults.csv",header = T,row.names = 1)
liveratacsigtf50 <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/liveratacsigfilejun2022out/knownResults.csv",header = T,row.names = 1)
lungatacsigtf50 <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/lungatacsigfilejun2022out/knownResults.csv",header = T,row.names = 1)
brownatacsigtf50 <- read.csv(file = "PASS1B Transcription Factor Paper Data/TF Data/brownatacsigfilejun2022out/knownResults.csv",header = T,row.names = 1)

atacsigtflist <- Reduce(intersect,list(rownames(gastroatacsigtf50),rownames(heartatacsigtf50),rownames(kidneyatacsigtf50)))
atacsigtflabel <- gsub("\\(.*","",atacsigtflist)

####
# Figure 1A-C
#####

# Figure 1B

tissuesigcountmat <- rbind(c(length(gastrornasig),
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
                             length(whiteatacsig)))
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
rownames(tissuesigcountmat) <- c("DEGs","DARs")
rownames(tissuetotalcountmat) <- c("DEGs","DARs")

colnames(tissuesigcountmat) <- c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC")
colnames(tissuetotalcountmat) <- c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC")

tissuesigquomat <- tissuesigcountmat/tissuetotalcountmat

tissuemeta <- data.frame(row.names = colnames(tissuesigquomat),
                         "Tissue" = colnames(tissuesigquomat))

pdf(file = "Figure 1B.pdf", width=13, height=5.5)
pheatmap(tissuesigquomat,cluster_cols = F,cluster_rows = F,angle_col = 0,display_numbers = tissuesigcountmat,breaks = seq(0,0.2,length.out = 9),color = brewer.pal(9,"Reds"),annotation_col = tissuemeta,annotation_colors = ann_cols,fontsize = 20,number_color = "black",cellwidth = 80,cellheight = 160)
dev.off()

# Figure 1C

combodeglist <- c(gastrornasig,heartrnasig,hippornasig,kidneyrnasig,liverrnasig,lungrnasig,brownrnasig,whiternasig)
combodarlist <- c(gastroatacsig,heartatacsig,hippoatacsig,kidneyatacsig,liveratacsig,lungatacsig,brownatacsig,whiteatacsig)

combogenelist <- c(rownames(gastrornanorm),
                                   rownames(heartrnanorm),
                                   rownames(hippornanorm),
                                   rownames(kidneyrnanorm),
                                   rownames(liverrnanorm),
                                   rownames(lungrnanorm),
                                   rownames(brownrnanorm),
                                   rownames(whiternanorm))

combopeaklist <- c(rownames(gastroatacnorm),
                                   rownames(heartatacnorm),
                                   rownames(hippoatacnorm),
                                   rownames(kidneyatacnorm),
                                   rownames(liveratacnorm),
                                   rownames(lungatacnorm),
                                   rownames(brownatacnorm),
                                   rownames(whiteatacnorm))

comboquantmat <- cbind(c(table(table(combopeaklist))[1],
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
                         table(table(combodeglist))[8]))
comboquantmat[is.na(comboquantmat)] <- 0
colnames(comboquantmat) <- c("All Regions","DARs","All Genes","DEGs")
rownames(comboquantmat) <- c("Tissue Specific",
                             "Shared by 2 tissues",
                             "Shared by 3 tissues",
                             "Shared by 4 tissues",
                             "Shared by 5 tissues",
                             "Shared by 6 tissues",
                             "Shared by 7 tissues",
                             "Shared by 8 tissues")

comboquantmatquo <- comboquantmat
for(i in 1:4){
  comboquantmatquo[,i] <- comboquantmat[,i]/sum(comboquantmat[,i])
}

comboquantdf <- data.frame("Percentage" = as.vector(comboquantmatquo),
                           "Specificity" = rep(c("Tissue Specific",
                                                 "Shared by 2 tissues",
                                                 "Shared by 3 tissues",
                                                 "Shared by 4 tissues",
                                                 "Shared by 5 tissues",
                                                 "Shared by 6 tissues",
                                                 "Shared by 7 tissues",
                                                 "Shared by 8 tissues"),4),
                           "Comparison" = c(rep("All Regions",8),
                                            rep("DARs",8),
                                            rep("All Genes",8),
                                            rep("DEGs",8)))
comboquantdf$Specificity <- factor(comboquantdf$Specificity,
                                   levels = c("Tissue Specific",
                                              "Shared by 2 tissues",
                                              "Shared by 3 tissues",
                                              "Shared by 4 tissues",
                                              "Shared by 5 tissues",
                                              "Shared by 6 tissues",
                                              "Shared by 7 tissues",
                                              "Shared by 8 tissues"))
comboquantdf$Comparison <- factor(comboquantdf$Comparison,
                                  levels = c("All Regions",
                                             "DARs",
                                             "All Genes",
                                             "DEGs"))
pdf(file = "Figure 1C.pdf", width=7, height=8)
ggplot(comboquantdf, aes(fill=Specificity, y=Percentage, x=Comparison)) + 
  geom_bar(position="fill", stat="identity") + scale_fill_brewer(palette = "Spectral") + theme_classic() + theme(text = element_text(size = 20),axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1)) + xlab(label = "")
dev.off()

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

pdf(file = "Figure S1A.pdf", width=7, height=6)
pheatmap(tissuedegintersectmatfrac,breaks = seq(0,1,length.out = 101),color = colorpanel(101,"white","firebrick","firebrick"),display_numbers = T,number_color = "black",annotation_row = tissuemeta,annotation_col = tissuemeta,annotation_colors = ann_cols,angle_col = 315,annotation_names_col = F,annotation_names_row = F,fontsize = 15,annotation_legend = F)
dev.off()

pdf(file = "Figure S1B.pdf", width=7, height=6)
pheatmap(tissuedarintersectmatfrac,breaks = seq(0,1,length.out = 101),color = colorpanel(101,"white","firebrick","firebrick"),display_numbers = T,number_color = "black",annotation_row = tissuemeta,annotation_col = tissuemeta,annotation_colors = ann_cols,angle_col = 315,annotation_names_col = F,annotation_names_row = F,fontsize = 15,annotation_legend = F)
dev.off()

####
# Supplemental Figure S2
#####

pdf(file = "Figure S2Aii.pdf",width = 6,height = 5)
pheatmap(gastrosigl2fcmat,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

pdf(file = "Figure S2Ai.pdf",width = 6,height = 5)
hist(gastrosigl2fcmat[abs(gastrosigl2fcmat) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "SKM-GN DEG L2FC")
dev.off()

pdf(file = "Figure S2Bii.pdf",width = 6,height = 5)
pheatmap(heartsigl2fcmat,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

pdf(file = "Figure S2Bi.pdf",width = 6,height = 5)
hist(heartsigl2fcmat[abs(heartsigl2fcmat) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "HEART DEG L2FC")
dev.off()

pdf(file = "Figure SCii.pdf",width = 6,height = 5)
pheatmap(hipposigl2fcmat,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

pdf(file = "Figure S2Ci.pdf",width = 6,height = 5)
hist(hipposigl2fcmat[abs(hipposigl2fcmat) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "HIPPOC DEG L2FC")
dev.off()

pdf(file = "Figure S2Dii.pdf",width = 6,height = 5)
pheatmap(kidneysigl2fcmat,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

pdf(file = "Figure S2Di.pdf",width = 6,height = 5)
hist(kidneysigl2fcmat[abs(kidneysigl2fcmat) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "KIDNEY DEG L2FC")
dev.off()

pdf(file = "Figure S2Eii.pdf",width = 6,height = 5)
pheatmap(liversigl2fcmat,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

pdf(file = "Figure S2Ei.pdf",width = 6,height = 5)
hist(liversigl2fcmat[abs(liversigl2fcmat) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "LIVER DEG L2FC")
dev.off()

pdf(file = "Figure S2Fii.pdf",width = 6,height = 5)
pheatmap(lungsigl2fcmat,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

pdf(file = "Figure S2Fi.pdf",width = 6,height = 5)
hist(lungsigl2fcmat[abs(lungsigl2fcmat) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "LUNG DEG L2FC")
dev.off()

pdf(file = "Figure S2Gii.pdf",width = 6,height = 5)
pheatmap(brownsigl2fcmat,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

pdf(file = "Figure S2Gi.pdf",width = 6,height = 5)
hist(brownsigl2fcmat[abs(brownsigl2fcmat) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "BAT DEG L2FC")
dev.off()

pdf(file = "Figure S2Hii.pdf",width = 6,height = 5)
pheatmap(whitesigl2fcmat,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,angle_col = 0,show_rownames = F,fontsize = 12)
dev.off()

pdf(file = "Figure S2Hi.pdf",width = 6,height = 5)
hist(whitesigl2fcmat[abs(whitesigl2fcmat) <= 4],breaks = seq(-4,4,0.2),xlab = "L2FC",main = "WAT-SC DEG L2FC")
dev.off()

####
# Supplemental Figure S3
#####

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


####
# Supplemental Figure S4
#####

# TO BE ADDED FROM MALENE AND ARCHANA

####
# Figure 1D-F
#####

# Fix two y-linked genes that have NA L2FC values for female time points
lungsigl2fcmat[is.na(lungsigl2fcmat)] <- 0
whitesigl2fcmat[is.na(whitesigl2fcmat)] <- 0

tissuel2fcdf <- data.frame(row.names = c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC"),
                           "DEG.Count" = c(dim(gastrosigl2fcmat)[1],
                                           dim(heartsigl2fcmat)[1],
                                           dim(hipposigl2fcmat)[1],
                                           dim(kidneysigl2fcmat)[1],
                                           dim(liversigl2fcmat)[1],
                                           dim(lungsigl2fcmat)[1],
                                           dim(brownsigl2fcmat)[1],
                                           dim(whitesigl2fcmat)[1]),
                           "DAR.Count" = c(dim(gastrosigatacl2fc)[1],
                                           dim(heartsigatacl2fc)[1],
                                           dim(hipposigatacl2fc)[1],
                                           dim(kidneysigatacl2fc)[1],
                                           dim(liversigatacl2fc)[1],
                                           dim(lungsigatacl2fc)[1],
                                           dim(brownsigatacl2fc)[1],
                                           dim(whitesigatacl2fc)[1]),
                           "DEG.L2FC.All.Agree" = c(sum(abs(rowSums(sign(gastrosigl2fcmat))) >= 8),
                                                    sum(abs(rowSums(sign(heartsigl2fcmat))) >= 8),
                                                    sum(abs(rowSums(sign(hipposigl2fcmat))) >= 8),
                                                    sum(abs(rowSums(sign(kidneysigl2fcmat))) >= 8),
                                                    sum(abs(rowSums(sign(liversigl2fcmat))) >= 8),
                                                    sum(abs(rowSums(sign(lungsigl2fcmat))) >= 8),
                                                    sum(abs(rowSums(sign(brownsigl2fcmat))) >= 8),
                                                    sum(abs(rowSums(sign(whitesigl2fcmat))) >= 8)),
                           "DAR.L2FC.All.Agree" = c(sum(abs(rowSums(sign(gastrosigatacl2fc))) >= 8),
                                                    sum(abs(rowSums(sign(heartsigatacl2fc))) >= 8),
                                                    sum(abs(rowSums(sign(hipposigatacl2fc))) >= 8),
                                                    sum(abs(rowSums(sign(kidneysigatacl2fc))) >= 8),
                                                    sum(abs(rowSums(sign(liversigatacl2fc))) >= 8),
                                                    sum(abs(rowSums(sign(lungsigatacl2fc))) >= 8),
                                                    sum(abs(rowSums(sign(brownsigatacl2fc))) >= 8),
                                                    sum(abs(rowSums(sign(whitesigatacl2fc))) >= 8)),
                           "DEG.L2FC.Most.Agree" = c(sum(abs(rowSums(sign(gastrosigl2fcmat))) >= 4),
                                                     sum(abs(rowSums(sign(heartsigl2fcmat))) >= 4),
                                                     sum(abs(rowSums(sign(hipposigl2fcmat))) >= 4),
                                                     sum(abs(rowSums(sign(kidneysigl2fcmat))) >= 4),
                                                     sum(abs(rowSums(sign(liversigl2fcmat))) >= 4),
                                                     sum(abs(rowSums(sign(lungsigl2fcmat))) >= 4),
                                                     sum(abs(rowSums(sign(brownsigl2fcmat))) >= 4),
                                                     sum(abs(rowSums(sign(whitesigl2fcmat))) >= 4)),
                           "DAR.L2FC.Most.Agree" = c(sum(abs(rowSums(sign(gastrosigatacl2fc))) >= 4),
                                                     sum(abs(rowSums(sign(heartsigatacl2fc))) >= 4),
                                                     sum(abs(rowSums(sign(hipposigatacl2fc))) >= 4),
                                                     sum(abs(rowSums(sign(kidneysigatacl2fc))) >= 4),
                                                     sum(abs(rowSums(sign(liversigatacl2fc))) >= 4),
                                                     sum(abs(rowSums(sign(lungsigatacl2fc))) >= 4),
                                                     sum(abs(rowSums(sign(brownsigatacl2fc))) >= 4),
                                                     sum(abs(rowSums(sign(whitesigatacl2fc))) >= 4)))

tissuequo <- data.frame(row.names = rownames(tissuel2fcdf),
                        "DEG.L2FC.All.Agree" = tissuel2fcdf$DEG.L2FC.All.Agree/tissuel2fcdf$DEG.Count,
                        "DEG.L2FC.Most.Agree" = tissuel2fcdf$DEG.L2FC.Most.Agree/tissuel2fcdf$DEG.Count,
                        "DAR.L2FC.All.Agree" = tissuel2fcdf$DAR.L2FC.All.Agree/tissuel2fcdf$DAR.Count,
                        "DAR.L2FC.Most.Agree" = tissuel2fcdf$DAR.L2FC.Most.Agree/tissuel2fcdf$DAR.Count)
tissuequodf <- data.frame(row.names = rownames(tissuequo),"Tissue" = rownames(tissuequo))


pdf(file = "Figure 1D.pdf",width = 8,height = 6)
pheatmap(t(tissuequo),cluster_rows = F,annotation_col = tissuequodf,show_colnames = F,annotation_colors = ann_cols,breaks = seq(0,1,length.out = 9),color = brewer.pal(9,"Reds"),display_numbers = T,number_color = "black",labels_row = c("DEG 100% Consistent","DEG 75% Consistent","DAR 100% Consistent","DAR 75% Consistent"),fontsize = 18,annotation_legend = F)
dev.off()


####
# Cellular Deconvolution Analysis for Figures 1E, 1F
#####

data("GSE20300")
data("GSE20300.grp")
data("cellCounts")
data("Garvan")
data("IRIS")
data("DMAP")

brownrnasymb <- unique(enstosym[rownames(brownrnanorm),"Symbol"])
brownrnatrim <- matrix(0L,nrow = length(brownrnasymb),ncol = dim(brownrnanorm)[2])
rownames(brownrnatrim) <- brownrnasymb
colnames(brownrnatrim) <- colnames(brownrnanorm)

for(i in 1:length(brownrnasymb)){
  if(i%%100 == 0){
    print(i)
  }
  oursymb <- brownrnasymb[i]
  ourens <- enstosym[enstosym$Symbol %in% oursymb,"Ensembl"]
  brownrnatrim[i,] <- colMeans(brownrnanorm[intersect(ourens,rownames(brownrnanorm)),])
}

brownrnatrim <- brownrnatrim[2:dim(brownrnatrim)[1],]
rownames(brownrnatrim) <- toupper(rownames(brownrnatrim))

adiabdz <- read.csv(file = "PASS1B Transcription Factor Paper Data/Deconvolution Data/adiabdcellz.csv",header = T,row.names = 1)
adiabdztrim <- adiabdz[intersect(rownames(brownrnatrim),rownames(adiabdz)),]
adiabdzfinal <- adiabdztrim[rowSums(is.na(adiabdztrim)) == 0,]
colnames(adiabdzfinal)[5] <- "Myeloid.Cells"

# SKM-GN

gastrornasymb <- unique(enstosym[rownames(gastrornanorm),"Symbol"])
gastrornatrim <- matrix(0L,nrow = length(gastrornasymb),ncol = dim(gastrornanorm)[2])
rownames(gastrornatrim) <- gastrornasymb
colnames(gastrornatrim) <- colnames(gastrornanorm)

for(i in 1:length(gastrornasymb)){
  if(i%%100 == 0){
    print(i)
  }
  oursymb <- gastrornasymb[i]
  ourens <- enstosym[enstosym$Symbol %in% oursymb,"Ensembl"]
  gastrornatrim[i,] <- colMeans(gastrornanorm[intersect(ourens,rownames(gastrornanorm)),])
}

gastrornatrim <- gastrornatrim[2:dim(gastrornatrim)[1],]
rownames(gastrornatrim) <- toupper(rownames(gastrornatrim))

gastrornameta <- data.frame(row.names = colnames(gastrornatrim),
                            "Sex" = rep("Female",length(colnames(gastrornatrim))),
                            "Group" = rep("",length(colnames(gastrornatrim))))
for(i in 1:length(colnames(gastrornatrim))){
  ourid <- colnames(gastrornatrim)[i]
  if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"] == 2){
    gastrornameta[i,"Sex"] <- "Male"
  }
  gastrornameta[i,"Group"] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"]
}
gastrornameta$Cohort <- paste(gastrornameta$Sex,gastrornameta$Group,sep = "_")
dmapTag=tagData(DMAP[, c("MEGA2", "ERY5", "MONO2", "GRAN3",
                         "DENDA1", "TCELLA7",
                         "NKA2", "BCELLA2")], 0.7, max=15,
                ref=gastrornatrim, ref.mean=F)
colnames(dmapTag)[1:2]=c("Megakaryocyte", "Erythrocyte")
colnames(dmapTag)[4] <- "Granulocyte"
colnames(dmapTag)[5] <- "Dendritic.Cell"
irisTag=tagData(IRIS[,c("Neutrophil-Resting","CD4Tcell-N0",
                        "Monocyte-Day0", "Bcell-naïve","NKcell-control" )],2, 
                max=15,ref=gastrornatrim, ref.mean=F)
colnames(irisTag)=c("Neutrophil","T.Cell", "Monocyte",
                    "B.Cell", "NK.Cell" )
tmpTag=combineTags(irisTag, dmapTag[,c(1,2,4,5)])


bloodimmgastroSPVs=getAllSPVs(gastrornatrim, grp=as.factor(gastrornameta$Cohort), tmpTag, method="mixed",plot=T, mix.par=0.3)

adiabdtag <- tagData(adiabdzfinal,max = 15,ref = gastrornatrim,ref.mean = F,cutoff = 0.75)
adiabdannogastroSPVs=getAllSPVs(gastrornatrim, grp=as.factor(gastrornameta$Cohort), adiabdtag, method="mixed",plot=T, mix.par=0.3)

rownames(bloodimmgastroSPVs) <- colnames(gastrornatrim)

gastrornameta[gastrornameta$Group %in% "Eight-week program Control Group","Group"] <- "control"
gastrornameta[gastrornameta$Group %in% "One-week program","Group"] <- "1w"
gastrornameta[gastrornameta$Group %in% "Two-week program","Group"] <- "2w"
gastrornameta[gastrornameta$Group %in% "Four-week program","Group"] <- "4w"
gastrornameta[gastrornameta$Group %in% "Eight-week program Training Group","Group"] <- "8w"

gastrornameta$Group <- factor(gastrornameta$Group,levels = c("control","1w","2w","4w","8w"))
gastrornameta <- gastrornameta[order(gastrornameta$Sex,gastrornameta$Group),]
bloodimmgastroSPVs <- bloodimmgastroSPVs[rownames(gastrornameta),]

rownames(adiabdannogastroSPVs) <- colnames(gastrornatrim)
adiabdannogastroSPVs <- adiabdannogastroSPVs[rownames(gastrornameta),]

# HEART

heartrnasymb <- unique(enstosym[rownames(heartrnanorm),"Symbol"])
heartrnatrim <- matrix(0L,nrow = length(heartrnasymb),ncol = dim(heartrnanorm)[2])
rownames(heartrnatrim) <- heartrnasymb
colnames(heartrnatrim) <- colnames(heartrnanorm)

for(i in 1:length(heartrnasymb)){
  if(i%%100 == 0){
    print(i)
  }
  oursymb <- heartrnasymb[i]
  ourens <- enstosym[enstosym$Symbol %in% oursymb,"Ensembl"]
  heartrnatrim[i,] <- colMeans(heartrnanorm[intersect(ourens,rownames(heartrnanorm)),])
}

heartrnatrim <- heartrnatrim[2:dim(heartrnatrim)[1],]
rownames(heartrnatrim) <- toupper(rownames(heartrnatrim))

heartrnameta <- data.frame(row.names = colnames(heartrnatrim),
                           "Sex" = rep("Female",length(colnames(heartrnatrim))),
                           "Group" = rep("",length(colnames(heartrnatrim))))
for(i in 1:length(colnames(heartrnatrim))){
  ourid <- colnames(heartrnatrim)[i]
  if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"] == 2){
    heartrnameta[i,"Sex"] <- "Male"
  }
  heartrnameta[i,"Group"] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"]
}
heartrnameta$Cohort <- paste(heartrnameta$Sex,heartrnameta$Group,sep = "_")
dmapTag=tagData(DMAP[, c("MEGA2", "ERY5", "MONO2", "GRAN3",
                         "DENDA1", "TCELLA7",
                         "NKA2", "BCELLA2")], 0.7, max=15,
                ref=heartrnatrim, ref.mean=F)
colnames(dmapTag)[1:2]=c("Megakaryocyte", "Erythrocyte")
colnames(dmapTag)[4] <- "Granulocyte"
colnames(dmapTag)[5] <- "Dendritic.Cell"
irisTag=tagData(IRIS[,c("Neutrophil-Resting","CD4Tcell-N0",
                        "Monocyte-Day0", "Bcell-naïve","NKcell-control" )],2, 
                max=15,ref=heartrnatrim, ref.mean=F)
colnames(irisTag)=c("Neutrophil","T.Cell", "Monocyte",
                    "B.Cell", "NK.Cell" )
tmpTag=combineTags(irisTag, dmapTag[,c(1,2,4,5)])

bloodimmheartSPVs=getAllSPVs(heartrnatrim, grp=as.factor(heartrnameta$Cohort), tmpTag, method="mixed",plot=T, mix.par=0.3)

adiabdtag <- tagData(adiabdzfinal,max = 15,ref = heartrnatrim,ref.mean = F,cutoff = 0.75)
adiabdannoheartSPVs=getAllSPVs(heartrnatrim, grp=as.factor(heartrnameta$Cohort), adiabdtag, method="mixed",plot=T, mix.par=0.3)

rownames(bloodimmheartSPVs) <- colnames(heartrnatrim)

heartrnameta[heartrnameta$Group %in% "Eight-week program Control Group","Group"] <- "control"
heartrnameta[heartrnameta$Group %in% "One-week program","Group"] <- "1w"
heartrnameta[heartrnameta$Group %in% "Two-week program","Group"] <- "2w"
heartrnameta[heartrnameta$Group %in% "Four-week program","Group"] <- "4w"
heartrnameta[heartrnameta$Group %in% "Eight-week program Training Group","Group"] <- "8w"

heartrnameta$Group <- factor(heartrnameta$Group,levels = c("control","1w","2w","4w","8w"))
heartrnameta <- heartrnameta[order(heartrnameta$Sex,heartrnameta$Group),]
bloodimmheartSPVs <- bloodimmheartSPVs[rownames(heartrnameta),]

rownames(adiabdannoheartSPVs) <- colnames(heartrnatrim)
adiabdannoheartSPVs <- adiabdannoheartSPVs[rownames(heartrnameta),]

# HIPPOC

hippornasymb <- unique(enstosym[rownames(hippornanorm),"Symbol"])
hippornatrim <- matrix(0L,nrow = length(hippornasymb),ncol = dim(hippornanorm)[2])
rownames(hippornatrim) <- hippornasymb
colnames(hippornatrim) <- colnames(hippornanorm)

for(i in 1:length(hippornasymb)){
  if(i%%100 == 0){
    print(i)
  }
  oursymb <- hippornasymb[i]
  ourens <- enstosym[enstosym$Symbol %in% oursymb,"Ensembl"]
  hippornatrim[i,] <- colMeans(hippornanorm[intersect(ourens,rownames(hippornanorm)),])
}

hippornatrim <- hippornatrim[2:dim(hippornatrim)[1],]
rownames(hippornatrim) <- toupper(rownames(hippornatrim))

hippornameta <- data.frame(row.names = colnames(hippornatrim),
                           "Sex" = rep("Female",length(colnames(hippornatrim))),
                           "Group" = rep("",length(colnames(hippornatrim))))
for(i in 1:length(colnames(hippornatrim))){
  ourid <- colnames(hippornatrim)[i]
  if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"] == 2){
    hippornameta[i,"Sex"] <- "Male"
  }
  hippornameta[i,"Group"] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"]
}
hippornameta$Cohort <- paste(hippornameta$Sex,hippornameta$Group,sep = "_")
dmapTag=tagData(DMAP[, c("MEGA2", "ERY5", "MONO2", "GRAN3",
                         "DENDA1", "TCELLA7",
                         "NKA2", "BCELLA2")], 0.7, max=15,
                ref=hippornatrim, ref.mean=F)
colnames(dmapTag)[1:2]=c("Megakaryocyte", "Erythrocyte")
colnames(dmapTag)[4] <- "Granulocyte"
colnames(dmapTag)[5] <- "Dendritic.Cell"
irisTag=tagData(IRIS[,c("Neutrophil-Resting","CD4Tcell-N0",
                        "Monocyte-Day0", "Bcell-naïve","NKcell-control" )],2, 
                max=15,ref=hippornatrim, ref.mean=F)
colnames(irisTag)=c("Neutrophil","T.Cell", "Monocyte",
                    "B.Cell", "NK.Cell" )
tmpTag=combineTags(irisTag, dmapTag[,c(1,2,4,5)])

bloodimmhippoSPVs=getAllSPVs(hippornatrim, grp=as.factor(hippornameta$Cohort), tmpTag, method="mixed",plot=T, mix.par=0.3)

adiabdtag <- tagData(adiabdzfinal,max = 15,ref = hippornatrim,ref.mean = F,cutoff = 0.55)
adiabdannohippoSPVs=getAllSPVs(hippornatrim, grp=as.factor(hippornameta$Cohort), adiabdtag, method="mixed",plot=T, mix.par=0.3)

rownames(bloodimmhippoSPVs) <- colnames(hippornatrim)

hippornameta[hippornameta$Group %in% "Eight-week program Control Group","Group"] <- "control"
hippornameta[hippornameta$Group %in% "One-week program","Group"] <- "1w"
hippornameta[hippornameta$Group %in% "Two-week program","Group"] <- "2w"
hippornameta[hippornameta$Group %in% "Four-week program","Group"] <- "4w"
hippornameta[hippornameta$Group %in% "Eight-week program Training Group","Group"] <- "8w"

hippornameta$Group <- factor(hippornameta$Group,levels = c("control","1w","2w","4w","8w"))
hippornameta <- hippornameta[order(hippornameta$Sex,hippornameta$Group),]
bloodimmhippoSPVs <- bloodimmhippoSPVs[rownames(hippornameta),]

rownames(adiabdannohippoSPVs) <- colnames(hippornatrim)
adiabdannohippoSPVs <- adiabdannohippoSPVs[rownames(hippornameta),]

# KIDNEY

kidneyrnasymb <- unique(enstosym[rownames(kidneyrnanorm),"Symbol"])
kidneyrnatrim <- matrix(0L,nrow = length(kidneyrnasymb),ncol = dim(kidneyrnanorm)[2])
rownames(kidneyrnatrim) <- kidneyrnasymb
colnames(kidneyrnatrim) <- colnames(kidneyrnanorm)

for(i in 1:length(kidneyrnasymb)){
  if(i%%100 == 0){
    print(i)
  }
  oursymb <- kidneyrnasymb[i]
  ourens <- enstosym[enstosym$Symbol %in% oursymb,"Ensembl"]
  kidneyrnatrim[i,] <- colMeans(kidneyrnanorm[intersect(ourens,rownames(kidneyrnanorm)),])
}

kidneyrnatrim <- kidneyrnatrim[2:dim(kidneyrnatrim)[1],]
rownames(kidneyrnatrim) <- toupper(rownames(kidneyrnatrim))

kidneyrnameta <- data.frame(row.names = colnames(kidneyrnatrim),
                            "Sex" = rep("Female",length(colnames(kidneyrnatrim))),
                            "Group" = rep("",length(colnames(kidneyrnatrim))))
for(i in 1:length(colnames(kidneyrnatrim))){
  ourid <- colnames(kidneyrnatrim)[i]
  if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"] == 2){
    kidneyrnameta[i,"Sex"] <- "Male"
  }
  kidneyrnameta[i,"Group"] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"]
}
kidneyrnameta$Cohort <- paste(kidneyrnameta$Sex,kidneyrnameta$Group,sep = "_")
dmapTag=tagData(DMAP[, c("MEGA2", "ERY5", "MONO2", "GRAN3",
                         "DENDA1", "TCELLA7",
                         "NKA2", "BCELLA2")], 0.7, max=15,
                ref=kidneyrnatrim, ref.mean=F)
colnames(dmapTag)[1:2]=c("Megakaryocyte", "Erythrocyte")
colnames(dmapTag)[4] <- "Granulocyte"
colnames(dmapTag)[5] <- "Dendritic.Cell"
irisTag=tagData(IRIS[,c("Neutrophil-Resting","CD4Tcell-N0",
                        "Monocyte-Day0", "Bcell-naïve","NKcell-control" )],2, 
                max=15,ref=kidneyrnatrim, ref.mean=F)
colnames(irisTag)=c("Neutrophil","T.Cell", "Monocyte",
                    "B.Cell", "NK.Cell" )
tmpTag=combineTags(irisTag, dmapTag[,c(1,2,4,5)])

bloodimmkidneySPVs=getAllSPVs(kidneyrnatrim, grp=as.factor(kidneyrnameta$Cohort), tmpTag, method="mixed",plot=T, mix.par=0.3)

adiabdtag <- tagData(adiabdzfinal,max = 15,ref = kidneyrnatrim,ref.mean = F,cutoff = 0.75)
adiabdannokidneySPVs=getAllSPVs(kidneyrnatrim, grp=as.factor(kidneyrnameta$Cohort), adiabdtag, method="mixed",plot=T, mix.par=0.3)

rownames(bloodimmkidneySPVs) <- colnames(kidneyrnatrim)

kidneyrnameta[kidneyrnameta$Group %in% "Eight-week program Control Group","Group"] <- "control"
kidneyrnameta[kidneyrnameta$Group %in% "One-week program","Group"] <- "1w"
kidneyrnameta[kidneyrnameta$Group %in% "Two-week program","Group"] <- "2w"
kidneyrnameta[kidneyrnameta$Group %in% "Four-week program","Group"] <- "4w"
kidneyrnameta[kidneyrnameta$Group %in% "Eight-week program Training Group","Group"] <- "8w"

kidneyrnameta$Group <- factor(kidneyrnameta$Group,levels = c("control","1w","2w","4w","8w"))
kidneyrnameta <- kidneyrnameta[order(kidneyrnameta$Sex,kidneyrnameta$Group),]
bloodimmkidneySPVs <- bloodimmkidneySPVs[rownames(kidneyrnameta),]

rownames(adiabdannokidneySPVs) <- colnames(kidneyrnatrim)
adiabdannokidneySPVs <- adiabdannokidneySPVs[rownames(kidneyrnameta),]

# LIVER

liverrnasymb <- unique(enstosym[rownames(liverrnanorm),"Symbol"])
liverrnatrim <- matrix(0L,nrow = length(liverrnasymb),ncol = dim(liverrnanorm)[2])
rownames(liverrnatrim) <- liverrnasymb
colnames(liverrnatrim) <- colnames(liverrnanorm)

for(i in 1:length(liverrnasymb)){
  if(i%%100 == 0){
    print(i)
  }
  oursymb <- liverrnasymb[i]
  ourens <- enstosym[enstosym$Symbol %in% oursymb,"Ensembl"]
  liverrnatrim[i,] <- colMeans(liverrnanorm[intersect(ourens,rownames(liverrnanorm)),])
}

liverrnatrim <- liverrnatrim[2:dim(liverrnatrim)[1],]
rownames(liverrnatrim) <- toupper(rownames(liverrnatrim))

liverrnameta <- data.frame(row.names = colnames(liverrnatrim),
                           "Sex" = rep("Female",length(colnames(liverrnatrim))),
                           "Group" = rep("",length(colnames(liverrnatrim))))
for(i in 1:length(colnames(liverrnatrim))){
  ourid <- colnames(liverrnatrim)[i]
  if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"] == 2){
    liverrnameta[i,"Sex"] <- "Male"
  }
  liverrnameta[i,"Group"] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"]
}
liverrnameta$Cohort <- paste(liverrnameta$Sex,liverrnameta$Group,sep = "_")
dmapTag=tagData(DMAP[, c("MEGA2", "ERY5", "MONO2", "GRAN3",
                         "DENDA1", "TCELLA7",
                         "NKA2", "BCELLA2")], 0.7, max=15,
                ref=liverrnatrim, ref.mean=F)
colnames(dmapTag)[1:2]=c("Megakaryocyte", "Erythrocyte")
colnames(dmapTag)[4] <- "Granulocyte"
colnames(dmapTag)[5] <- "Dendritic.Cell"
irisTag=tagData(IRIS[,c("Neutrophil-Resting","CD4Tcell-N0",
                        "Monocyte-Day0", "Bcell-naïve","NKcell-control" )],2, 
                max=15,ref=liverrnatrim, ref.mean=F)
colnames(irisTag)=c("Neutrophil","T.Cell", "Monocyte",
                    "B.Cell", "NK.Cell" )
tmpTag=combineTags(irisTag, dmapTag[,c(1,2,4,5)])

bloodimmliverSPVs=getAllSPVs(liverrnatrim, grp=as.factor(liverrnameta$Cohort), tmpTag, method="mixed",plot=T, mix.par=0.3)

adiabdtag <- tagData(adiabdzfinal,max = 15,ref = liverrnatrim,ref.mean = F,cutoff = 0.70)
adiabdannoliverSPVs=getAllSPVs(liverrnatrim, grp=as.factor(liverrnameta$Cohort), adiabdtag, method="mixed",plot=T, mix.par=0.3)

rownames(bloodimmliverSPVs) <- colnames(liverrnatrim)

liverrnameta[liverrnameta$Group %in% "Eight-week program Control Group","Group"] <- "control"
liverrnameta[liverrnameta$Group %in% "One-week program","Group"] <- "1w"
liverrnameta[liverrnameta$Group %in% "Two-week program","Group"] <- "2w"
liverrnameta[liverrnameta$Group %in% "Four-week program","Group"] <- "4w"
liverrnameta[liverrnameta$Group %in% "Eight-week program Training Group","Group"] <- "8w"

liverrnameta$Group <- factor(liverrnameta$Group,levels = c("control","1w","2w","4w","8w"))
liverrnameta <- liverrnameta[order(liverrnameta$Sex,liverrnameta$Group),]
bloodimmliverSPVs <- bloodimmliverSPVs[rownames(liverrnameta),]

rownames(adiabdannoliverSPVs) <- colnames(liverrnatrim)
adiabdannoliverSPVs <- adiabdannoliverSPVs[rownames(liverrnameta),]


# LUNG

lungrnasymb <- unique(enstosym[rownames(lungrnanorm),"Symbol"])
lungrnatrim <- matrix(0L,nrow = length(lungrnasymb),ncol = dim(lungrnanorm)[2])
rownames(lungrnatrim) <- lungrnasymb
colnames(lungrnatrim) <- colnames(lungrnanorm)

for(i in 1:length(lungrnasymb)){
  if(i%%100 == 0){
    print(i)
  }
  oursymb <- lungrnasymb[i]
  ourens <- enstosym[enstosym$Symbol %in% oursymb,"Ensembl"]
  lungrnatrim[i,] <- colMeans(lungrnanorm[intersect(ourens,rownames(lungrnanorm)),])
}

lungrnatrim <- lungrnatrim[2:dim(lungrnatrim)[1],]
rownames(lungrnatrim) <- toupper(rownames(lungrnatrim))

lungrnameta <- data.frame(row.names = colnames(lungrnatrim),
                          "Sex" = rep("Female",length(colnames(lungrnatrim))),
                          "Group" = rep("",length(colnames(lungrnatrim))))
for(i in 1:length(colnames(lungrnatrim))){
  ourid <- colnames(lungrnatrim)[i]
  if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"] == 2){
    lungrnameta[i,"Sex"] <- "Male"
  }
  lungrnameta[i,"Group"] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"]
}
lungrnameta$Cohort <- paste(lungrnameta$Sex,lungrnameta$Group,sep = "_")
dmapTag=tagData(DMAP[, c("MEGA2", "ERY5", "MONO2", "GRAN3",
                         "DENDA1", "TCELLA7",
                         "NKA2", "BCELLA2")], 0.7, max=15,
                ref=lungrnatrim, ref.mean=F)
colnames(dmapTag)[1:2]=c("Megakaryocyte", "Erythrocyte")
colnames(dmapTag)[4] <- "Granulocyte"
colnames(dmapTag)[5] <- "Dendritic.Cell"
irisTag=tagData(IRIS[,c("Neutrophil-Resting","CD4Tcell-N0",
                        "Monocyte-Day0", "Bcell-naïve","NKcell-control" )],2, 
                max=15,ref=lungrnatrim, ref.mean=F)
colnames(irisTag)=c("Neutrophil","T.Cell", "Monocyte",
                    "B.Cell", "NK.Cell" )
tmpTag=combineTags(irisTag, dmapTag[,c(1,2,4,5)])

bloodimmlungSPVs=getAllSPVs(lungrnatrim, grp=as.factor(lungrnameta$Cohort), tmpTag, method="mixed",plot=T, mix.par=0.3)

adiabdtag <- tagData(adiabdzfinal,max = 15,ref = lungrnatrim,ref.mean = F,cutoff = 0.75)
adiabdannolungSPVs=getAllSPVs(lungrnatrim, grp=as.factor(lungrnameta$Cohort), adiabdtag, method="mixed",plot=T, mix.par=0.3)

rownames(bloodimmlungSPVs) <- colnames(lungrnatrim)

lungrnameta[lungrnameta$Group %in% "Eight-week program Control Group","Group"] <- "control"
lungrnameta[lungrnameta$Group %in% "One-week program","Group"] <- "1w"
lungrnameta[lungrnameta$Group %in% "Two-week program","Group"] <- "2w"
lungrnameta[lungrnameta$Group %in% "Four-week program","Group"] <- "4w"
lungrnameta[lungrnameta$Group %in% "Eight-week program Training Group","Group"] <- "8w"

lungrnameta$Group <- factor(lungrnameta$Group,levels = c("control","1w","2w","4w","8w"))
lungrnameta <- lungrnameta[order(lungrnameta$Sex,lungrnameta$Group),]
bloodimmlungSPVs <- bloodimmlungSPVs[rownames(lungrnameta),]

rownames(adiabdannolungSPVs) <- colnames(lungrnatrim)
adiabdannolungSPVs <- adiabdannolungSPVs[rownames(lungrnameta),]

# BAT

brownrnasymb <- unique(enstosym[rownames(brownrnanorm),"Symbol"])
brownrnatrim <- matrix(0L,nrow = length(brownrnasymb),ncol = dim(brownrnanorm)[2])
rownames(brownrnatrim) <- brownrnasymb
colnames(brownrnatrim) <- colnames(brownrnanorm)

for(i in 1:length(brownrnasymb)){
  if(i%%100 == 0){
    print(i)
  }
  oursymb <- brownrnasymb[i]
  ourens <- enstosym[enstosym$Symbol %in% oursymb,"Ensembl"]
  brownrnatrim[i,] <- colMeans(brownrnanorm[intersect(ourens,rownames(brownrnanorm)),])
}

brownrnatrim <- brownrnatrim[2:dim(brownrnatrim)[1],]
rownames(brownrnatrim) <- toupper(rownames(brownrnatrim))

brownrnameta <- data.frame(row.names = colnames(brownrnatrim),
                           "Sex" = rep("Female",length(colnames(brownrnatrim))),
                           "Group" = rep("",length(colnames(brownrnatrim))))
for(i in 1:length(colnames(brownrnatrim))){
  ourid <- colnames(brownrnatrim)[i]
  if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"] == 2){
    brownrnameta[i,"Sex"] <- "Male"
  }
  brownrnameta[i,"Group"] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"]
}
brownrnameta$Cohort <- paste(brownrnameta$Sex,brownrnameta$Group,sep = "_")

dmapTag=tagData(DMAP[, c("MEGA2", "ERY5", "MONO2", "GRAN3",
                         "DENDA1", "TCELLA7",
                         "NKA2", "BCELLA2")], 0.7, max=15,
                ref=brownrnatrim, ref.mean=F)
colnames(dmapTag)[1:2]=c("Megakaryocyte", "Erythrocyte")
colnames(dmapTag)[4] <- "Granulocyte"
colnames(dmapTag)[5] <- "Dendritic.Cell"
irisTag=tagData(IRIS[,c("Neutrophil-Resting","CD4Tcell-N0",
                        "Monocyte-Day0", "Bcell-naïve","NKcell-control" )],2, 
                max=15,ref=brownrnatrim, ref.mean=F)
colnames(irisTag)=c("Neutrophil","T.Cell", "Monocyte",
                    "B.Cell", "NK.Cell" )
tmpTag=combineTags(irisTag, dmapTag[,c(1,2,4,5)])

bloodimmbrownSPVs=getAllSPVs(brownrnatrim, grp=as.factor(brownrnameta$Cohort), tmpTag, method="mixed",plot=T, mix.par=0.3)

adiabdtag <- tagData(adiabdzfinal,max = 15,ref = brownrnatrim,ref.mean = F,cutoff = 0.75)
adiabdannobrownSPVs=getAllSPVs(brownrnatrim, grp=as.factor(brownrnameta$Cohort), adiabdtag, method="mixed",plot=T, mix.par=0.3)

rownames(bloodimmbrownSPVs) <- colnames(brownrnatrim)


brownrnameta[brownrnameta$Group %in% "Eight-week program Control Group","Group"] <- "control"
brownrnameta[brownrnameta$Group %in% "One-week program","Group"] <- "1w"
brownrnameta[brownrnameta$Group %in% "Two-week program","Group"] <- "2w"
brownrnameta[brownrnameta$Group %in% "Four-week program","Group"] <- "4w"
brownrnameta[brownrnameta$Group %in% "Eight-week program Training Group","Group"] <- "8w"

brownrnameta$Group <- factor(brownrnameta$Group,levels = c("control","1w","2w","4w","8w"))
brownrnameta <- brownrnameta[order(brownrnameta$Sex,brownrnameta$Group),]
bloodimmbrownSPVs <- bloodimmbrownSPVs[rownames(brownrnameta),]

rownames(adiabdannobrownSPVs) <- colnames(brownrnatrim)
adiabdannobrownSPVs <- adiabdannobrownSPVs[rownames(brownrnameta),]

# WAT-SC

whiternasymb <- unique(enstosym[rownames(whiternanorm),"Symbol"])
whiternatrim <- matrix(0L,nrow = length(whiternasymb),ncol = dim(whiternanorm)[2])
rownames(whiternatrim) <- whiternasymb
colnames(whiternatrim) <- colnames(whiternanorm)

for(i in 1:length(whiternasymb)){
  if(i%%100 == 0){
    print(i)
  }
  oursymb <- whiternasymb[i]
  ourens <- enstosym[enstosym$Symbol %in% oursymb,"Ensembl"]
  whiternatrim[i,] <- colMeans(whiternanorm[intersect(ourens,rownames(whiternanorm)),])
}

whiternatrim <- whiternatrim[2:dim(whiternatrim)[1],]
rownames(whiternatrim) <- toupper(rownames(whiternatrim))

whiternameta <- data.frame(row.names = colnames(whiternatrim),
                           "Sex" = rep("Female",length(colnames(whiternatrim))),
                           "Group" = rep("",length(colnames(whiternatrim))))
for(i in 1:length(colnames(whiternatrim))){
  ourid <- colnames(whiternatrim)[i]
  if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"] == 2){
    whiternameta[i,"Sex"] <- "Male"
  }
  whiternameta[i,"Group"] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"]
}
whiternameta$Cohort <- paste(whiternameta$Sex,whiternameta$Group,sep = "_")
dmapTag=tagData(DMAP[, c("MEGA2", "ERY5", "MONO2", "GRAN3",
                         "DENDA1", "TCELLA7",
                         "NKA2", "BCELLA2")], 0.7, max=15,
                ref=whiternatrim, ref.mean=F)
colnames(dmapTag)[1:2]=c("Megakaryocyte", "Erythrocyte")
colnames(dmapTag)[4] <- "Granulocyte"
colnames(dmapTag)[5] <- "Dendritic.Cell"
irisTag=tagData(IRIS[,c("Neutrophil-Resting","CD4Tcell-N0",
                        "Monocyte-Day0", "Bcell-naïve","NKcell-control" )],2, 
                max=15,ref=whiternatrim, ref.mean=F)
colnames(irisTag)=c("Neutrophil","T.Cell", "Monocyte",
                    "B.Cell", "NK.Cell" )
tmpTag=combineTags(irisTag, dmapTag[,c(1,2,4,5)])

bloodimmwhiteSPVs=getAllSPVs(whiternatrim, grp=as.factor(whiternameta$Cohort), tmpTag, method="mixed",plot=T, mix.par=0.3)

adiabdtag <- tagData(adiabdzfinal,max = 15,ref = whiternatrim,ref.mean = F,cutoff = 0.75)
adiabdannowhiteSPVs=getAllSPVs(whiternatrim, grp=as.factor(whiternameta$Cohort), adiabdtag, method="mixed",plot=T, mix.par=0.3)

rownames(bloodimmwhiteSPVs) <- colnames(whiternatrim)

whiternameta[whiternameta$Group %in% "Eight-week program Control Group","Group"] <- "control"
whiternameta[whiternameta$Group %in% "One-week program","Group"] <- "1w"
whiternameta[whiternameta$Group %in% "Two-week program","Group"] <- "2w"
whiternameta[whiternameta$Group %in% "Four-week program","Group"] <- "4w"
whiternameta[whiternameta$Group %in% "Eight-week program Training Group","Group"] <- "8w"

whiternameta$Group <- factor(whiternameta$Group,levels = c("control","1w","2w","4w","8w"))
whiternameta <- whiternameta[order(whiternameta$Sex,whiternameta$Group),]
bloodimmwhiteSPVs <- bloodimmwhiteSPVs[rownames(whiternameta),]

rownames(adiabdannowhiteSPVs) <- colnames(whiternatrim)
adiabdannowhiteSPVs <- adiabdannowhiteSPVs[rownames(whiternameta),]

dev.off()

gastrocomboSPVs <- cbind(bloodimmgastroSPVs,adiabdannogastroSPVs[,c("Adipocytes","Endothelial.Cells","PCV.Endothelial.Cells","FAP.Cells","Pericytes")])
heartcomboSPVs <- cbind(bloodimmheartSPVs,adiabdannoheartSPVs[,c("Adipocytes","Endothelial.Cells","PCV.Endothelial.Cells","FAP.Cells","Pericytes")])
hippocomboSPVs <- cbind(bloodimmhippoSPVs,adiabdannohippoSPVs[,c("Adipocytes","Endothelial.Cells","PCV.Endothelial.Cells","FAP.Cells","Pericytes")])
kidneycomboSPVs <- cbind(bloodimmkidneySPVs,adiabdannokidneySPVs[,c("Adipocytes","Endothelial.Cells","PCV.Endothelial.Cells","FAP.Cells","Pericytes")])
livercomboSPVs <- cbind(bloodimmliverSPVs,adiabdannoliverSPVs[,c("Adipocytes","Endothelial.Cells","PCV.Endothelial.Cells","FAP.Cells","Pericytes")])
lungcomboSPVs <- cbind(bloodimmlungSPVs,adiabdannolungSPVs[,c("Adipocytes","Endothelial.Cells","PCV.Endothelial.Cells","FAP.Cells","Pericytes")])
browncomboSPVs <- cbind(bloodimmbrownSPVs,adiabdannobrownSPVs[,c("Adipocytes","Endothelial.Cells","PCV.Endothelial.Cells","FAP.Cells","Pericytes")])
whitecomboSPVs <- cbind(bloodimmwhiteSPVs,adiabdannowhiteSPVs[,c("Adipocytes","Endothelial.Cells","PCV.Endothelial.Cells","FAP.Cells","Pericytes")])

tissuedecontrainsigmat <- matrix(0L,nrow = 8,ncol = length(colnames(whitecomboSPVs)))
rownames(tissuedecontrainsigmat) <- c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC")
colnames(tissuedecontrainsigmat) <- colnames(whitecomboSPVs)
for(i in 1:length(colnames(tissuedecontrainsigmat))){
  ourcelltype <- colnames(tissuedecontrainsigmat)[i]
  ourdf <- data.frame(row.names = rownames(gastrocomboSPVs),
                      "Group" = gastrornameta[,"Group"],
                      "Sex" = gastrornameta[,"Sex"],
                      "Decon" = gastrocomboSPVs[,i])
  ourkwresult <- kruskal.test(Decon ~ Group, data = ourdf)
  tissuedecontrainsigmat["SKM-GN",i] <- ourkwresult$p.value
  ourdf <- data.frame(row.names = rownames(heartcomboSPVs),
                      "Group" = heartrnameta[,"Group"],
                      "Sex" = heartrnameta[,"Sex"],
                      "Decon" = heartcomboSPVs[,i])
  ourkwresult <- kruskal.test(Decon ~ Group, data = ourdf)
  tissuedecontrainsigmat["HEART",i] <- ourkwresult$p.value
  ourdf <- data.frame(row.names = rownames(hippocomboSPVs),
                      "Group" = hippornameta[,"Group"],
                      "Sex" = hippornameta[,"Sex"],
                      "Decon" = hippocomboSPVs[,i])
  ourkwresult <- kruskal.test(Decon ~ Group, data = ourdf)
  tissuedecontrainsigmat["HIPPOC",i] <- ourkwresult$p.value
  ourdf <- data.frame(row.names = rownames(kidneycomboSPVs),
                      "Group" = kidneyrnameta[,"Group"],
                      "Sex" = kidneyrnameta[,"Sex"],
                      "Decon" = kidneycomboSPVs[,i])
  ourkwresult <- kruskal.test(Decon ~ Group, data = ourdf)
  tissuedecontrainsigmat["KIDNEY",i] <- ourkwresult$p.value
  ourdf <- data.frame(row.names = rownames(livercomboSPVs),
                      "Group" = liverrnameta[,"Group"],
                      "Sex" = liverrnameta[,"Sex"],
                      "Decon" = livercomboSPVs[,i])
  ourkwresult <- kruskal.test(Decon ~ Group, data = ourdf)
  tissuedecontrainsigmat["LIVER",i] <- ourkwresult$p.value
  ourdf <- data.frame(row.names = rownames(lungcomboSPVs),
                      "Group" = lungrnameta[,"Group"],
                      "Sex" = lungrnameta[,"Sex"],
                      "Decon" = lungcomboSPVs[,i])
  ourkwresult <- kruskal.test(Decon ~ Group, data = ourdf)
  tissuedecontrainsigmat["LUNG",i] <- ourkwresult$p.value
  ourdf <- data.frame(row.names = rownames(browncomboSPVs),
                      "Group" = brownrnameta[,"Group"],
                      "Sex" = brownrnameta[,"Sex"],
                      "Decon" = browncomboSPVs[,i])
  ourkwresult <- kruskal.test(Decon ~ Group, data = ourdf)
  tissuedecontrainsigmat["BAT",i] <- ourkwresult$p.value
  ourdf <- data.frame(row.names = rownames(whitecomboSPVs),
                      "Group" = whiternameta[,"Group"],
                      "Sex" = whiternameta[,"Sex"],
                      "Decon" = whitecomboSPVs[,i])
  ourkwresult <- kruskal.test(Decon ~ Group, data = ourdf)
  tissuedecontrainsigmat["WAT-SC",i] <- ourkwresult$p.value
}

tissuedeconsexsigmat <- matrix(0L,nrow = 8,ncol = length(colnames(whitecomboSPVs)))
rownames(tissuedeconsexsigmat) <- c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC")
colnames(tissuedeconsexsigmat) <- colnames(whitecomboSPVs)
for(i in 1:length(colnames(tissuedeconsexsigmat))){
  ourcelltype <- colnames(tissuedeconsexsigmat)[i]
  ourdf <- data.frame(row.names = rownames(gastrocomboSPVs),
                      "Group" = gastrornameta[,"Group"],
                      "Sex" = gastrornameta[,"Sex"],
                      "Decon" = gastrocomboSPVs[,i])
  ourkwresult <- kruskal.test(Decon ~ Sex, data = ourdf)
  tissuedeconsexsigmat["SKM-GN",i] <- ourkwresult$p.value
  ourdf <- data.frame(row.names = rownames(heartcomboSPVs),
                      "Group" = heartrnameta[,"Group"],
                      "Sex" = heartrnameta[,"Sex"],
                      "Decon" = heartcomboSPVs[,i])
  ourkwresult <- kruskal.test(Decon ~ Sex, data = ourdf)
  tissuedeconsexsigmat["HEART",i] <- ourkwresult$p.value
  ourdf <- data.frame(row.names = rownames(hippocomboSPVs),
                      "Group" = hippornameta[,"Group"],
                      "Sex" = hippornameta[,"Sex"],
                      "Decon" = hippocomboSPVs[,i])
  ourkwresult <- kruskal.test(Decon ~ Sex, data = ourdf)
  tissuedeconsexsigmat["HIPPOC",i] <- ourkwresult$p.value
  ourdf <- data.frame(row.names = rownames(kidneycomboSPVs),
                      "Group" = kidneyrnameta[,"Group"],
                      "Sex" = kidneyrnameta[,"Sex"],
                      "Decon" = kidneycomboSPVs[,i])
  ourkwresult <- kruskal.test(Decon ~ Sex, data = ourdf)
  tissuedeconsexsigmat["KIDNEY",i] <- ourkwresult$p.value
  ourdf <- data.frame(row.names = rownames(livercomboSPVs),
                      "Group" = liverrnameta[,"Group"],
                      "Sex" = liverrnameta[,"Sex"],
                      "Decon" = livercomboSPVs[,i])
  ourkwresult <- kruskal.test(Decon ~ Sex, data = ourdf)
  tissuedeconsexsigmat["LIVER",i] <- ourkwresult$p.value
  ourdf <- data.frame(row.names = rownames(lungcomboSPVs),
                      "Group" = lungrnameta[,"Group"],
                      "Sex" = lungrnameta[,"Sex"],
                      "Decon" = lungcomboSPVs[,i])
  ourkwresult <- kruskal.test(Decon ~ Sex, data = ourdf)
  tissuedeconsexsigmat["LUNG",i] <- ourkwresult$p.value
  ourdf <- data.frame(row.names = rownames(browncomboSPVs),
                      "Group" = brownrnameta[,"Group"],
                      "Sex" = brownrnameta[,"Sex"],
                      "Decon" = browncomboSPVs[,i])
  ourkwresult <- kruskal.test(Decon ~ Sex, data = ourdf)
  tissuedeconsexsigmat["BAT",i] <- ourkwresult$p.value
  ourdf <- data.frame(row.names = rownames(whitecomboSPVs),
                      "Group" = whiternameta[,"Group"],
                      "Sex" = whiternameta[,"Sex"],
                      "Decon" = whitecomboSPVs[,i])
  ourkwresult <- kruskal.test(Decon ~ Sex, data = ourdf)
  tissuedeconsexsigmat["WAT-SC",i] <- ourkwresult$p.value
}

tissue_cols <- list("Tissue" = c("SKM-GN" = "#088c03",
                                 "HEART" = "#f28b2f",
                                 "HIPPOC" = "#bf7534",
                                 "KIDNEY" = "#7553a7",
                                 "LIVER" = "#da6c75",
                                 "LUNG" = "#04bf8a",
                                 "BAT" = "#8c5220",
                                 "WAT-SC" = "#214da6"))

decontissuemeta <- data.frame(row.names = rownames(tissuedeconsexsigmat),
                              "Tissue" = rownames(tissuedeconsexsigmat))

# Figure 1E
pdf(file = "Figure 1E.pdf",width = 7,height = 6)
pheatmap(t(-log10(tissuedecontrainsigmat[c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC"),])), angle_col = "0",breaks = seq(0,4,length.out = 9),color = brewer.pal(9,"Reds"),annotation_col = tissuemeta,annotation_colors = ann_cols,cluster_cols = F,show_colnames = F,fontsize = 15)
dev.off()

# Figure 1F
pdf(file = "Figure 1F.pdf",width = 7,height = 6)
pheatmap(t(-log10(tissuedeconsexsigmat[c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC"),])), angle_col = "0",breaks = seq(0,4,length.out = 9),color = brewer.pal(9,"Reds"),annotation_col = tissuemeta,annotation_colors = ann_cols,cluster_cols = F,show_colnames = F,fontsize = 15)
dev.off()

####
# Supplemental Figure S5
#####

pdf(file = "Supplemental Figure S5A.pdf",width = 8,height = 3)
pheatmap(t(scale(browncomboSPVs)),cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_colnames = F,annotation_col = brownrnameta[,c("Group","Sex")],annotation_colors = ann_cols,cellwidth = 6.5)
dev.off()

pdf(file = "Supplemental Figure S5B.pdf",width = 8,height = 3)
pheatmap(t(scale(whitecomboSPVs)),cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_colnames = F,annotation_col = whiternameta[,c("Group","Sex")],annotation_colors = ann_cols,cellwidth = 6.5)
dev.off()

pdf(file = "Supplemental Figure S5C.pdf",width = 8,height = 3)
pheatmap(t(scale(gastrocomboSPVs)),cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_colnames = F,annotation_col = gastrornameta[,c("Group","Sex")],annotation_colors = ann_cols,cellwidth = 6.5)
dev.off()

pdf(file = "Supplemental Figure S5D.pdf",width = 8,height = 3)
pheatmap(t(scale(heartcomboSPVs)),cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_colnames = F,annotation_col = heartrnameta[,c("Group","Sex")],annotation_colors = ann_cols,cellwidth = 6.5)
dev.off()

pdf(file = "Supplemental Figure S5E.pdf",width = 8,height = 3)
pheatmap(t(scale(hippocomboSPVs)),cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_colnames = F,annotation_col = hippornameta[,c("Group","Sex")],annotation_colors = ann_cols,cellwidth = 6.5)
dev.off()

pdf(file = "Supplemental Figure S5F.pdf",width = 8,height = 3)
pheatmap(t(scale(kidneycomboSPVs)),cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_colnames = F,annotation_col = kidneyrnameta[,c("Group","Sex")],annotation_colors = ann_cols,cellwidth = 6.5)
dev.off()

pdf(file = "Supplemental Figure S5G.pdf",width = 8,height = 3)
pheatmap(t(scale(livercomboSPVs)),cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_colnames = F,annotation_col = liverrnameta[,c("Group","Sex")],annotation_colors = ann_cols,cellwidth = 6.5)
dev.off()

pdf(file = "Supplemental Figure S5H.pdf",width = 8,height = 3)
pheatmap(t(scale(lungcomboSPVs)),cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_colnames = F,annotation_col = lungrnameta[,c("Group","Sex")],annotation_colors = ann_cols,cellwidth = 6.5)
dev.off()

####
# Figure 1G-L
#####

# Add peak annotations to training data
atactraining <- epigen_atac_seq$training_dea
atactraining$custom_annotation <- ""
atactraining$custom_annotation <- peakanno[atactraining$feature_ID,"custom_annotation"]

tab = table(atactraining$tissue_abbreviation,atactraining$custom_annotation)

#myorder = colnames(tab)
myorder = c("Upstream (<5kb)", "Promoter (<=1kb)","Promoter (1-2kb)", "5' UTR","Exon",
            "Intron", "3' UTR","Downstream (<5kb)", "Distal Intergenic","Overlaps Gene")

mycol = pal_d3()(length(myorder))
names(mycol) = myorder
tab2 = tab[,myorder]
tab2Perc = t(apply(t(tab2), 2, function(x){x*100/sum(x,na.rm=T)}))

# Figure 1G
pdf(file = "Figure 1G.pdf", width=7, height=6)
mosaicplot(tab2Perc,las=1,col=genomic_features_col$Region[colnames(tab2Perc)], main="Genomic Distribution of All Peaks",cex.axis = 1.5)
dev.off()

pcutoff = 0.1
trainingSig = atactraining %>% filter(adj_p_value<pcutoff)

tab = table(trainingSig$tissue_abbreviation,trainingSig$custom_annotation)
tab3 = tab[,myorder[myorder %in% colnames(tab)]]
tab3Perc = t(apply(t(tab3), 2, function(x){x*100/sum(x,na.rm=T)}))

# Figure 1H
pdf(file = "Figure 1H.pdf", width=7, height=6)
mosaicplot(tab3Perc,las=1,col=genomic_features_col$Region[colnames(tab3Perc)], main="Genomic Distribution of Significant Peaks")
dev.off()

tab2 <- tab2[,colnames(tab3)]

tab2rs = rowSums(tab2)
tab3rs = rowSums(tab3)

myres = lapply(1:nrow(tab2), function(i){
  ans2 = lapply(1:ncol(tab2), function(j){
    mydat =c(tab3[i,j],tab2[i,j],tab3rs[i],tab2rs[i])
    ft = fisher.test(matrix(mydat,nrow=2))
    ans = tibble(
      p1=mydat[1],
      p2=mydat[2],
      p3=mydat[3],
      p4=mydat[4],
      tissue=rownames(tab2)[i],feature=colnames(tab2)[j], 
      p_value=ft$p.value, estimate=ft$estimate)
    return(ans)
    
  }) 
  ans3 = do.call(rbind, ans2)
  return(ans3)
})
myres2 = do.call(rbind, myres)

myres2 = myres2 %>% mutate(p_adj=p.adjust(p_value,method="BH"))
myres3 = myres2 %>% mutate(sig=p_adj<0.01, type=ifelse(estimate>1,1,-1))
write.table(myres3,file = "fisher_test_sig_peaks_by_tissue.txt")
#write_tsv(myres3, file.path(outfolder, "fisher_test_sig_peaks_by_tissue.txt"))
myres3sig = myres3 %>% filter(sig)
myres3sig$feature = factor(myres3sig$feature,myorder)



tab4 <- matrix(0L,nrow = dim(tab2)[2],ncol = dim(tab2)[1])
rownames(tab4) <- colnames(tab2)
colnames(tab4) <- rownames(tab2)

for(i in 1:nrow(myres3sig)){
  tab4[myres3sig$feature[i], myres3sig$tissue[i]] = -log2(myres3sig$p_adj[i])
  tab4[myres3sig$feature[i], myres3sig$tissue[i]] = tab4[myres3sig$feature[i], myres3sig$tissue[i]] * myres3sig$type[i]
}

col_fun = colorRamp2(c(-1, 0,1), c("steelblue","white","firebrick"))

# Figure 1K
pdf(file = "Figure 1K.pdf", width=8, height=3)
Heatmap(
  tab4,
  name="Sig",
  col=col_fun,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  rect_gp = gpar(col = "black"),
  show_heatmap_legend=FALSE,
  column_names_rot=TRUE,
  column_names_centered=TRUE,
  row_names_centered=TRUE,
  row_names_side="left"
)
dev.off()

gastroactivegeneactivepeaktab <- table(peakanno[(intersect(rownames(peakanno[peakanno$ensembl_gene %in% rownames(gastrornanorm),]),unique(epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "SKM-GN","feature_ID"]))),"custom_annotation"])
heartactivegeneactivepeaktab <- table(peakanno[(intersect(rownames(peakanno[peakanno$ensembl_gene %in% rownames(heartrnanorm),]),unique(epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "HEART","feature_ID"]))),"custom_annotation"])
hippoactivegeneactivepeaktab <- table(peakanno[(intersect(rownames(peakanno[peakanno$ensembl_gene %in% rownames(hippornanorm),]),unique(epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "HIPPOC","feature_ID"]))),"custom_annotation"])
kidneyactivegeneactivepeaktab <- table(peakanno[(intersect(rownames(peakanno[peakanno$ensembl_gene %in% rownames(kidneyrnanorm),]),unique(epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "KIDNEY","feature_ID"]))),"custom_annotation"])
liveractivegeneactivepeaktab <- table(peakanno[(intersect(rownames(peakanno[peakanno$ensembl_gene %in% rownames(liverrnanorm),]),unique(epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "LIVER","feature_ID"]))),"custom_annotation"])
lungactivegeneactivepeaktab <- table(peakanno[(intersect(rownames(peakanno[peakanno$ensembl_gene %in% rownames(lungrnanorm),]),unique(epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "LUNG","feature_ID"]))),"custom_annotation"])
brownactivegeneactivepeaktab <- table(peakanno[(intersect(rownames(peakanno[peakanno$ensembl_gene %in% rownames(brownrnanorm),]),unique(epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "BAT","feature_ID"]))),"custom_annotation"])
whiteactivegeneactivepeaktab <- table(peakanno[(intersect(rownames(peakanno[peakanno$ensembl_gene %in% rownames(whiternanorm),]),unique(epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "WAT-SC","feature_ID"]))),"custom_annotation"])

tab <- rbind(gastroactivegeneactivepeaktab[c("Distal Intergenic","Promoter (<=1kb)","Upstream (<5kb)","Intron","5' UTR","Exon","Promoter (1-2kb)","3' UTR","Downstream (<5kb)","Overlaps Gene")],
             heartactivegeneactivepeaktab[c("Distal Intergenic","Promoter (<=1kb)","Upstream (<5kb)","Intron","5' UTR","Exon","Promoter (1-2kb)","3' UTR","Downstream (<5kb)","Overlaps Gene")],
             hippoactivegeneactivepeaktab[c("Distal Intergenic","Promoter (<=1kb)","Upstream (<5kb)","Intron","5' UTR","Exon","Promoter (1-2kb)","3' UTR","Downstream (<5kb)","Overlaps Gene")],
             kidneyactivegeneactivepeaktab[c("Distal Intergenic","Promoter (<=1kb)","Upstream (<5kb)","Intron","5' UTR","Exon","Promoter (1-2kb)","3' UTR","Downstream (<5kb)","Overlaps Gene")],
             liveractivegeneactivepeaktab[c("Distal Intergenic","Promoter (<=1kb)","Upstream (<5kb)","Intron","5' UTR","Exon","Promoter (1-2kb)","3' UTR","Downstream (<5kb)","Overlaps Gene")],
             lungactivegeneactivepeaktab[c("Distal Intergenic","Promoter (<=1kb)","Upstream (<5kb)","Intron","5' UTR","Exon","Promoter (1-2kb)","3' UTR","Downstream (<5kb)","Overlaps Gene")],
             brownactivegeneactivepeaktab[c("Distal Intergenic","Promoter (<=1kb)","Upstream (<5kb)","Intron","5' UTR","Exon","Promoter (1-2kb)","3' UTR","Downstream (<5kb)","Overlaps Gene")],
             whiteactivegeneactivepeaktab[c("Distal Intergenic","Promoter (<=1kb)","Upstream (<5kb)","Intron","5' UTR","Exon","Promoter (1-2kb)","3' UTR","Downstream (<5kb)","Overlaps Gene")])
rownames(tab) <- c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC")

tabPerc <- t(apply(t(tab), 2, function(x) {
  x * 100 / sum(x, na.rm = T)
}))

tab = tab[,myorder]
tabPerc = tabPerc[,myorder]

## plots
# Figure 1I
pdf(file = "Figure 1I.pdf", width = 7, height = 6)
mosaicplot(tabPerc, las = 1, col = genomic_features_col$Region[colnames(tabPerc)], main = "Genomic Distribution of Active Genes Percentage")
dev.off()

gastrosiggeneactivepeaktab <- table(peakanno[(intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastrornasig,]),unique(epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "SKM-GN","feature_ID"]))),"custom_annotation"])
heartsiggeneactivepeaktab <- table(peakanno[(intersect(rownames(peakanno[peakanno$ensembl_gene %in% heartrnasig,]),unique(epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "HEART","feature_ID"]))),"custom_annotation"])
hipposiggeneactivepeaktab <- table(peakanno[(intersect(rownames(peakanno[peakanno$ensembl_gene %in% hippornasig,]),unique(epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "HIPPOC","feature_ID"]))),"custom_annotation"])
kidneysiggeneactivepeaktab <- table(peakanno[(intersect(rownames(peakanno[peakanno$ensembl_gene %in% kidneyrnasig,]),unique(epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "KIDNEY","feature_ID"]))),"custom_annotation"])
liversiggeneactivepeaktab <- table(peakanno[(intersect(rownames(peakanno[peakanno$ensembl_gene %in% liverrnasig,]),unique(epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "LIVER","feature_ID"]))),"custom_annotation"])
lungsiggeneactivepeaktab <- table(peakanno[(intersect(rownames(peakanno[peakanno$ensembl_gene %in% lungrnasig,]),unique(epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "LUNG","feature_ID"]))),"custom_annotation"])
brownsiggeneactivepeaktab <- table(peakanno[(intersect(rownames(peakanno[peakanno$ensembl_gene %in% brownrnasig,]),unique(epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "BAT","feature_ID"]))),"custom_annotation"])
whitesiggeneactivepeaktab <- table(peakanno[(intersect(rownames(peakanno[peakanno$ensembl_gene %in% whiternasig,]),unique(epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "WAT-SC","feature_ID"]))),"custom_annotation"])

tab3 <- matrix(0L,nrow = 8,ncol = 10)
colnames(tab3) <- colnames(tab)
rownames(tab3) <- rownames(tab)
for(i in 1:dim(tab3)[2]){
  if(colnames(tab3)[i] %in% names(gastrosiggeneactivepeaktab)){
    tab3["SKM-GN",i] <- gastrosiggeneactivepeaktab[colnames(tab3)[i]]
  }
  if(colnames(tab3)[i] %in% names(heartsiggeneactivepeaktab)){
    tab3["HEART",i] <- heartsiggeneactivepeaktab[colnames(tab3)[i]]
  }
  if(colnames(tab3)[i] %in% names(hipposiggeneactivepeaktab)){
    tab3["HIPPOC",i] <- hipposiggeneactivepeaktab[colnames(tab3)[i]]
  }
  if(colnames(tab3)[i] %in% names(kidneysiggeneactivepeaktab)){
    tab3["KIDNEY",i] <- kidneysiggeneactivepeaktab[colnames(tab3)[i]]
  }
  if(colnames(tab3)[i] %in% names(liversiggeneactivepeaktab)){
    tab3["LIVER",i] <- liversiggeneactivepeaktab[colnames(tab3)[i]]
  }
  if(colnames(tab3)[i] %in% names(lungsiggeneactivepeaktab)){
    tab3["LUNG",i] <- lungsiggeneactivepeaktab[colnames(tab3)[i]]
  }
  if(colnames(tab3)[i] %in% names(brownsiggeneactivepeaktab)){
    tab3["BAT",i] <- brownsiggeneactivepeaktab[colnames(tab3)[i]]
  }
  if(colnames(tab3)[i] %in% names(whitesiggeneactivepeaktab)){
    tab3["WAT-SC",i] <- whitesiggeneactivepeaktab[colnames(tab3)[i]]
  }
}

tab3Perc <- t(apply(t(tab3), 2, function(x) {
  x * 100 / sum(x, na.rm = T)
}))

tab3 = tab3[,myorder]
tab3Perc = tab3Perc[,myorder]

## plots
# Figure 1J
pdf(file = "Figure 1J.pdf", width = 7, height = 6)
mosaicplot(tab3Perc, las = 1, col = genomic_features_col$Region[colnames(tab3Perc)], main = "Genomic Distribution of DEGaP Percentage")
dev.off()

# Figure 1L
tabrs = rowSums(tab)
tab3rs = rowSums(tab3)

myres = lapply(1:nrow(tab), function(i){
  ans2 = lapply(1:ncol(tab), function(j){
    mydat =c(tab3[i,j],tab[i,j],tab3rs[i],tabrs[i])
    ft = fisher.test(matrix(mydat,nrow=2))
    ans = tibble(
      p1=mydat[1],
      p2=mydat[2],
      p3=mydat[3],
      p4=mydat[4],
      tissue=rownames(tab)[i],feature=colnames(tab)[j], 
      p_value=ft$p.value, estimate=ft$estimate)
    return(ans)
    
  }) 
  ans3 = do.call(rbind, ans2)
  return(ans3)
})
myres2 = do.call(rbind, myres)

myres2 = myres2 %>% mutate(p_adj=p.adjust(p_value,method="BH"))
myres3 = myres2 %>% mutate(sig=p_adj<0.01, type=ifelse(estimate>1,1,-1))
write.table(myres3,file = "fisher_test_sig_peaks_by_tissue.txt")
#write_tsv(myres3, file.path(outfolder, "fisher_test_sig_peaks_by_tissue.txt"))
myres3sig = myres3 %>% filter(sig)
myres3sig$feature = factor(myres3sig$feature,myorder)

tab4 <- matrix(0L,nrow = dim(tab)[2],ncol = dim(tab)[1])
rownames(tab4) <- colnames(tab)
colnames(tab4) <- rownames(tab)

for(i in 1:nrow(myres3sig)){
  tab4[myres3sig$feature[i], myres3sig$tissue[i]] = -log2(myres3sig$p_adj[i])
  tab4[myres3sig$feature[i], myres3sig$tissue[i]] = tab4[myres3sig$feature[i], myres3sig$tissue[i]] * myres3sig$type[i]
}

col_fun = colorRamp2(c(-1, 0,1), c("steelblue","white","firebrick"))

# Figure 1L
pdf(file = "Figure 1L.pdf", width=8, height=3)
Heatmap(
  tab4,
  name="Sig",
  col=col_fun,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  rect_gp = gpar(col = "black"),
  show_heatmap_legend=FALSE,
  column_names_rot=TRUE,
  column_names_centered=TRUE,
  row_names_centered=TRUE,
  row_names_side="left"
)
dev.off()


####
# Figure 2
#####

gastroatacsigdf <- data.frame(row.names = gastroatacsig,
                              "ATAC" = gastroatacsig,
                              "RNA" = peakanno[gastroatacsig,"ensembl_gene"])

heartatacsigdf <- data.frame(row.names = heartatacsig,
                             "ATAC" = heartatacsig,
                             "RNA" = peakanno[heartatacsig,"ensembl_gene"])

hippoatacsigdf <- data.frame(row.names = hippoatacsig,
                             "ATAC" = hippoatacsig,
                             "RNA" = peakanno[hippoatacsig,"ensembl_gene"])

kidneyatacsigdf <- data.frame(row.names = kidneyatacsig,
                              "ATAC" = kidneyatacsig,
                              "RNA" = peakanno[kidneyatacsig,"ensembl_gene"])

liveratacsigdf <- data.frame(row.names = liveratacsig,
                             "ATAC" = liveratacsig,
                             "RNA" = peakanno[liveratacsig,"ensembl_gene"])

lungatacsigdf <- data.frame(row.names = lungatacsig,
                            "ATAC" = lungatacsig,
                            "RNA" = peakanno[lungatacsig,"ensembl_gene"])

brownatacsigdf <- data.frame(row.names = brownatacsig,
                             "ATAC" = brownatacsig,
                             "RNA" = peakanno[brownatacsig,"ensembl_gene"])

whiteatacsigdf <- data.frame(row.names = whiteatacsig,
                             "ATAC" = whiteatacsig,
                             "RNA" = peakanno[whiteatacsig,"ensembl_gene"])

atacbarplotdf <- data.frame("Tissue" = c(rep("HIPPOC",2),
                                         rep("SKM-GN",2),
                                         rep("HEART",2),
                                         rep("KIDNEY",2),
                                         rep("LUNG",2),
                                         rep("LIVER",2),
                                         rep("BAT",2),
                                         rep("WAT-SC",2)),
                            "Group" = factor(rep(c("Significant Gene","Not Significant Gene"),8),
                                             levels = c("Significant Gene","Not Significant Gene")),
                            "Count" = rep(0,16))
atacbarplotdf[1,"Count"] <- sum(hippoatacsigdf$RNA %in% hippornasig)
atacbarplotdf[2,"Count"] <- length(hippoatacsigdf$RNA) - sum(hippoatacsigdf$RNA %in% hippornasig)

atacbarplotdf[3,"Count"] <- sum(gastroatacsigdf$RNA %in% gastrornasig)
atacbarplotdf[4,"Count"] <- length(gastroatacsigdf$RNA) - sum(gastroatacsigdf$RNA %in% gastrornasig)

atacbarplotdf[5,"Count"] <- sum(heartatacsigdf$RNA %in% heartrnasig)
atacbarplotdf[6,"Count"] <- length(heartatacsigdf$RNA) - sum(heartatacsigdf$RNA %in% heartrnasig)

atacbarplotdf[7,"Count"] <- sum(kidneyatacsigdf$RNA %in% kidneyrnasig)
atacbarplotdf[8,"Count"] <- length(kidneyatacsigdf$RNA) - sum(kidneyatacsigdf$RNA %in% kidneyrnasig)

atacbarplotdf[9,"Count"] <- sum(lungatacsigdf$RNA %in% lungrnasig)
atacbarplotdf[10,"Count"] <- length(lungatacsigdf$RNA) - sum(lungatacsigdf$RNA %in% lungrnasig)

atacbarplotdf[11,"Count"] <- sum(liveratacsigdf$RNA %in% liverrnasig)
atacbarplotdf[12,"Count"] <- length(liveratacsigdf$RNA) - sum(liveratacsigdf$RNA %in% liverrnasig)

atacbarplotdf[13,"Count"] <- sum(brownatacsigdf$RNA %in% brownrnasig)
atacbarplotdf[14,"Count"] <- length(brownatacsigdf$RNA) - sum(brownatacsigdf$RNA %in% brownrnasig)

atacbarplotdf[15,"Count"] <- sum(whiteatacsigdf$RNA %in% whiternasig)
atacbarplotdf[16,"Count"] <- length(whiteatacsigdf$RNA) - sum(whiteatacsigdf$RNA %in% whiternasig)




rnabarplotdf <- data.frame("Tissue" = c(rep("HIPPOC",2),
                                        rep("SKM-GN",2),
                                        rep("HEART",2),
                                        rep("KIDNEY",2),
                                        rep("LUNG",2),
                                        rep("LIVER",2),
                                        rep("BAT",2),
                                        rep("WAT-SC",2)),
                           "Group" = rep(c("Contains DAR","Does not contain DAR"),8),
                           "Count" = rep(0,16))

rnabarplotdf[1,"Count"] <- sum(hippornasig %in% hippoatacsigdf$RNA)
rnabarplotdf[2,"Count"] <- length(hippornasig) - sum(hippornasig %in% hippoatacsigdf$RNA)
rnabarplotdf[3,"Count"] <- sum(gastrornasig %in% gastroatacsigdf$RNA)
rnabarplotdf[4,"Count"] <- length(gastrornasig) - sum(gastrornasig %in% gastroatacsigdf$RNA)
rnabarplotdf[5,"Count"] <- sum(heartrnasig %in% heartatacsigdf$RNA)
rnabarplotdf[6,"Count"] <- length(heartrnasig) - sum(heartrnasig %in% heartatacsigdf$RNA)
rnabarplotdf[7,"Count"] <- sum(kidneyrnasig %in% kidneyatacsigdf$RNA)
rnabarplotdf[8,"Count"] <- length(kidneyrnasig) - sum(kidneyrnasig %in% kidneyatacsigdf$RNA)
rnabarplotdf[9,"Count"] <- sum(lungrnasig %in% lungatacsigdf$RNA)
rnabarplotdf[10,"Count"] <- length(lungrnasig) - sum(lungrnasig %in% lungatacsigdf$RNA)
rnabarplotdf[11,"Count"] <- sum(liverrnasig %in% liveratacsigdf$RNA)
rnabarplotdf[12,"Count"] <- length(liverrnasig) - sum(liverrnasig %in% liveratacsigdf$RNA)
rnabarplotdf[13,"Count"] <- sum(brownrnasig %in% brownatacsigdf$RNA)
rnabarplotdf[14,"Count"] <- length(brownrnasig) - sum(brownrnasig %in% brownatacsigdf$RNA)
rnabarplotdf[15,"Count"] <- sum(whiternasig %in% whiteatacsigdf$RNA)
rnabarplotdf[16,"Count"] <- length(whiternasig) - sum(whiternasig %in% whiteatacsigdf$RNA)

# Figure 2A
pdf(file = "Figure 2A.pdf",width = 7.5,height = 4)
ggplot(atacbarplotdf, aes(fill=Group, y=Count, x=Tissue)) + 
  geom_bar(position="stack", stat="identity") + scale_fill_manual(values = c("blue3","darkorange2")) + theme_classic() + ggtitle("DARs")
dev.off()

# Figure 2B
pdf(file = "Figure 2B.pdf",width = 7.5,height = 4)
ggplot(rnabarplotdf, aes(fill=Group, y=Count, x=Tissue)) + 
  geom_bar(position="stack", stat="identity") + scale_fill_manual(values = c("blue3","darkorange2")) + theme_classic() + ggtitle("DEGs")
dev.off()

####
# Figure 2C-D
#####

# SKM-GN

fingastrornasig <- intersect(gastrornasig,unique(peakanno$ensembl_gene))

gastropeakdistance <- matrix(0L,nrow = length(gastroatacsig),ncol = length(fingastrornasig))
rownames(gastropeakdistance) <- gastroatacsig
colnames(gastropeakdistance) <- fingastrornasig

gastroatacsiganno <- peakanno[gastroatacsig,c("chrom","start","end")]
gastroatacsiganno$mid <- (gastroatacsiganno$start + gastroatacsiganno$end)/2
gastrornasiganno <- data.frame(row.names = fingastrornasig,"chrom" = rep("1",length(fingastrornasig)),"start" = rep(1,length(fingastrornasig)))
for(i in 1:length(fingastrornasig)){
  if(i%%50 == 0){
    print(i)
  }
  ourgene <- fingastrornasig[i]
  ourgastrornapeaks <- peakanno[peakanno$ensembl_gene %in% ourgene,c("chrom","geneStart")]
  gastrornasiganno[i,"chrom"] <- ourgastrornapeaks$chrom[1]
  gastrornasiganno[i,"start"] <- ourgastrornapeaks$geneStart[1]
}

for(i in 1:length(gastroatacsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(fingastrornasig)){
    gastropeakdistance[i,j] <- (gastrornasiganno$chrom[j] == gastroatacsiganno$chrom[i])*(abs(gastrornasiganno$start[j] - gastroatacsiganno$mid[i]))
  }
}

gastropeakdistancedf <- data.frame(row.names = gastroatacsig,
                                   "Region" = peakanno[gastroatacsig,"custom_annotation"],
                                   "Distance" = rep(0,length(gastroatacsig)))
for(i in 1:length(gastroatacsig)){
  
  ourrow <- gastropeakdistance[i,]
  gastropeakdistancedf[i,"Distance"] <- min(ourrow[ourrow > 0])
  
}

# HEART

finheartrnasig <- intersect(heartrnasig,unique(peakanno$ensembl_gene))

heartpeakdistance <- matrix(0L,nrow = length(heartatacsig),ncol = length(finheartrnasig))
rownames(heartpeakdistance) <- heartatacsig
colnames(heartpeakdistance) <- finheartrnasig

heartatacsiganno <- peakanno[heartatacsig,c("chrom","start","end")]
heartatacsiganno$mid <- (heartatacsiganno$start + heartatacsiganno$end)/2
heartrnasiganno <- data.frame(row.names = finheartrnasig,"chrom" = rep("1",length(finheartrnasig)),"start" = rep(1,length(finheartrnasig)))
for(i in 1:length(finheartrnasig)){
  if(i%%50 == 0){
    print(i)
  }
  ourgene <- finheartrnasig[i]
  ourheartrnapeaks <- peakanno[peakanno$ensembl_gene %in% ourgene,c("chrom","geneStart")]
  heartrnasiganno[i,"chrom"] <- ourheartrnapeaks$chrom[1]
  heartrnasiganno[i,"start"] <- ourheartrnapeaks$geneStart[1]
}

for(i in 1:length(heartatacsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(finheartrnasig)){
    heartpeakdistance[i,j] <- (heartrnasiganno$chrom[j] == heartatacsiganno$chrom[i])*(abs(heartrnasiganno$start[j] - heartatacsiganno$mid[i]))
  }
}

heartpeakdistancedf <- data.frame(row.names = heartatacsig,
                                  "Region" = peakanno[heartatacsig,"custom_annotation"],
                                  "Distance" = rep(0,length(heartatacsig)))
for(i in 1:length(heartatacsig)){
  
  ourrow <- heartpeakdistance[i,]
  heartpeakdistancedf[i,"Distance"] <- min(ourrow[ourrow > 0])
  
}

# HIPPOC

finhippornasig <- intersect(hippornasig,unique(peakanno$ensembl_gene))

hippopeakdistance <- matrix(0L,nrow = length(hippoatacsig),ncol = length(finhippornasig))
rownames(hippopeakdistance) <- hippoatacsig
colnames(hippopeakdistance) <- finhippornasig

hippoatacsiganno <- peakanno[hippoatacsig,c("chrom","start","end")]
hippoatacsiganno$mid <- (hippoatacsiganno$start + hippoatacsiganno$end)/2
hippornasiganno <- data.frame(row.names = finhippornasig,"chrom" = rep("1",length(finhippornasig)),"start" = rep(1,length(finhippornasig)))
for(i in 1:length(finhippornasig)){
  if(i%%50 == 0){
    print(i)
  }
  ourgene <- finhippornasig[i]
  ourhippornapeaks <- peakanno[peakanno$ensembl_gene %in% ourgene,c("chrom","geneStart")]
  hippornasiganno[i,"chrom"] <- ourhippornapeaks$chrom[1]
  hippornasiganno[i,"start"] <- ourhippornapeaks$geneStart[1]
}

for(i in 1:length(hippoatacsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(finhippornasig)){
    hippopeakdistance[i,j] <- (hippornasiganno$chrom[j] == hippoatacsiganno$chrom[i])*(abs(hippornasiganno$start[j] - hippoatacsiganno$mid[i]))
  }
}

hippopeakdistancedf <- data.frame(row.names = hippoatacsig,
                                  "Region" = peakanno[hippoatacsig,"custom_annotation"],
                                  "Distance" = rep(0,length(hippoatacsig)))
for(i in 1:length(hippoatacsig)){
  
  ourrow <- hippopeakdistance[i,]
  hippopeakdistancedf[i,"Distance"] <- min(ourrow[ourrow > 0])
  
}

# KIDNEY

finkidneyrnasig <- intersect(kidneyrnasig,unique(peakanno$ensembl_gene))

kidneypeakdistance <- matrix(0L,nrow = length(kidneyatacsig),ncol = length(finkidneyrnasig))
rownames(kidneypeakdistance) <- kidneyatacsig
colnames(kidneypeakdistance) <- finkidneyrnasig

kidneyatacsiganno <- peakanno[kidneyatacsig,c("chrom","start","end")]
kidneyatacsiganno$mid <- (kidneyatacsiganno$start + kidneyatacsiganno$end)/2
kidneyrnasiganno <- data.frame(row.names = finkidneyrnasig,"chrom" = rep("1",length(finkidneyrnasig)),"start" = rep(1,length(finkidneyrnasig)))
for(i in 1:length(finkidneyrnasig)){
  if(i%%50 == 0){
    print(i)
  }
  ourgene <- finkidneyrnasig[i]
  ourkidneyrnapeaks <- peakanno[peakanno$ensembl_gene %in% ourgene,c("chrom","geneStart")]
  kidneyrnasiganno[i,"chrom"] <- ourkidneyrnapeaks$chrom[1]
  kidneyrnasiganno[i,"start"] <- ourkidneyrnapeaks$geneStart[1]
}

for(i in 1:length(kidneyatacsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(finkidneyrnasig)){
    kidneypeakdistance[i,j] <- (kidneyrnasiganno$chrom[j] == kidneyatacsiganno$chrom[i])*(abs(kidneyrnasiganno$start[j] - kidneyatacsiganno$mid[i]))
  }
}

kidneypeakdistancedf <- data.frame(row.names = kidneyatacsig,
                                   "Region" = peakanno[kidneyatacsig,"custom_annotation"],
                                   "Distance" = rep(0,length(kidneyatacsig)))
for(i in 1:length(kidneyatacsig)){
  
  ourrow <- kidneypeakdistance[i,]
  kidneypeakdistancedf[i,"Distance"] <- min(ourrow[ourrow > 0])
  
}

# LIVER

finliverrnasig <- intersect(liverrnasig,unique(peakanno$ensembl_gene))

liverpeakdistance <- matrix(0L,nrow = length(liveratacsig),ncol = length(finliverrnasig))
rownames(liverpeakdistance) <- liveratacsig
colnames(liverpeakdistance) <- finliverrnasig

liveratacsiganno <- peakanno[liveratacsig,c("chrom","start","end")]
liveratacsiganno$mid <- (liveratacsiganno$start + liveratacsiganno$end)/2
liverrnasiganno <- data.frame(row.names = finliverrnasig,"chrom" = rep("1",length(finliverrnasig)),"start" = rep(1,length(finliverrnasig)))
for(i in 1:length(finliverrnasig)){
  if(i%%50 == 0){
    print(i)
  }
  ourgene <- finliverrnasig[i]
  ourliverrnapeaks <- peakanno[peakanno$ensembl_gene %in% ourgene,c("chrom","geneStart")]
  liverrnasiganno[i,"chrom"] <- ourliverrnapeaks$chrom[1]
  liverrnasiganno[i,"start"] <- ourliverrnapeaks$geneStart[1]
}

for(i in 1:length(liveratacsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(finliverrnasig)){
    liverpeakdistance[i,j] <- (liverrnasiganno$chrom[j] == liveratacsiganno$chrom[i])*(abs(liverrnasiganno$start[j] - liveratacsiganno$mid[i]))
  }
}

liverpeakdistancedf <- data.frame(row.names = liveratacsig,
                                  "Region" = peakanno[liveratacsig,"custom_annotation"],
                                  "Distance" = rep(0,length(liveratacsig)))
for(i in 1:length(liveratacsig)){
  
  ourrow <- liverpeakdistance[i,]
  liverpeakdistancedf[i,"Distance"] <- min(ourrow[ourrow > 0])
  
}

# LUNG

finlungrnasig <- intersect(lungrnasig,unique(peakanno$ensembl_gene))

lungpeakdistance <- matrix(0L,nrow = length(lungatacsig),ncol = length(finlungrnasig))
rownames(lungpeakdistance) <- lungatacsig
colnames(lungpeakdistance) <- finlungrnasig

lungatacsiganno <- peakanno[lungatacsig,c("chrom","start","end")]
lungatacsiganno$mid <- (lungatacsiganno$start + lungatacsiganno$end)/2
lungrnasiganno <- data.frame(row.names = finlungrnasig,"chrom" = rep("1",length(finlungrnasig)),"start" = rep(1,length(finlungrnasig)))
for(i in 1:length(finlungrnasig)){
  if(i%%50 == 0){
    print(i)
  }
  ourgene <- finlungrnasig[i]
  ourlungrnapeaks <- peakanno[peakanno$ensembl_gene %in% ourgene,c("chrom","geneStart")]
  lungrnasiganno[i,"chrom"] <- ourlungrnapeaks$chrom[1]
  lungrnasiganno[i,"start"] <- ourlungrnapeaks$geneStart[1]
}

for(i in 1:length(lungatacsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(finlungrnasig)){
    lungpeakdistance[i,j] <- (lungrnasiganno$chrom[j] == lungatacsiganno$chrom[i])*(abs(lungrnasiganno$start[j] - lungatacsiganno$mid[i]))
  }
}

lungpeakdistancedf <- data.frame(row.names = lungatacsig,
                                 "Region" = peakanno[lungatacsig,"custom_annotation"],
                                 "Distance" = rep(0,length(lungatacsig)))
for(i in 1:length(lungatacsig)){
  
  ourrow <- lungpeakdistance[i,]
  lungpeakdistancedf[i,"Distance"] <- min(ourrow[ourrow > 0])
  
}

# BAT

finbrownrnasig <- intersect(brownrnasig,unique(peakanno$ensembl_gene))

brownpeakdistance <- matrix(0L,nrow = length(brownatacsig),ncol = length(finbrownrnasig))
rownames(brownpeakdistance) <- brownatacsig
colnames(brownpeakdistance) <- finbrownrnasig

brownatacsiganno <- peakanno[brownatacsig,c("chrom","start","end")]
brownatacsiganno$mid <- (brownatacsiganno$start + brownatacsiganno$end)/2
brownrnasiganno <- data.frame(row.names = finbrownrnasig,"chrom" = rep("1",length(finbrownrnasig)),"start" = rep(1,length(finbrownrnasig)))
for(i in 1:length(finbrownrnasig)){
  if(i%%50 == 0){
    print(i)
  }
  ourgene <- finbrownrnasig[i]
  ourbrownrnapeaks <- peakanno[peakanno$ensembl_gene %in% ourgene,c("chrom","geneStart")]
  brownrnasiganno[i,"chrom"] <- ourbrownrnapeaks$chrom[1]
  brownrnasiganno[i,"start"] <- ourbrownrnapeaks$geneStart[1]
}

for(i in 1:length(brownatacsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(finbrownrnasig)){
    brownpeakdistance[i,j] <- (brownrnasiganno$chrom[j] == brownatacsiganno$chrom[i])*(abs(brownrnasiganno$start[j] - brownatacsiganno$mid[i]))
  }
}

brownpeakdistancedf <- data.frame(row.names = brownatacsig,
                                  "Region" = peakanno[brownatacsig,"custom_annotation"],
                                  "Distance" = rep(0,length(brownatacsig)))
for(i in 1:length(brownatacsig)){
  
  ourrow <- brownpeakdistance[i,]
  brownpeakdistancedf[i,"Distance"] <- min(ourrow[ourrow > 0])
  
}

# WAT-SC

finwhiternasig <- intersect(whiternasig,unique(peakanno$ensembl_gene))

whitepeakdistance <- matrix(0L,nrow = length(whiteatacsig),ncol = length(finwhiternasig))
rownames(whitepeakdistance) <- whiteatacsig
colnames(whitepeakdistance) <- finwhiternasig

whiteatacsiganno <- peakanno[whiteatacsig,c("chrom","start","end")]
whiteatacsiganno$mid <- (whiteatacsiganno$start + whiteatacsiganno$end)/2
whiternasiganno <- data.frame(row.names = finwhiternasig,"chrom" = rep("1",length(finwhiternasig)),"start" = rep(1,length(finwhiternasig)))
for(i in 1:length(finwhiternasig)){
  if(i%%50 == 0){
    print(i)
  }
  ourgene <- finwhiternasig[i]
  ourwhiternapeaks <- peakanno[peakanno$ensembl_gene %in% ourgene,c("chrom","geneStart")]
  whiternasiganno[i,"chrom"] <- ourwhiternapeaks$chrom[1]
  whiternasiganno[i,"start"] <- ourwhiternapeaks$geneStart[1]
}

for(i in 1:length(whiteatacsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(finwhiternasig)){
    whitepeakdistance[i,j] <- (whiternasiganno$chrom[j] == whiteatacsiganno$chrom[i])*(abs(whiternasiganno$start[j] - whiteatacsiganno$mid[i]))
  }
}

whitepeakdistancedf <- data.frame(row.names = whiteatacsig,
                                  "Region" = peakanno[whiteatacsig,"custom_annotation"],
                                  "Distance" = rep(0,length(whiteatacsig)))
for(i in 1:length(whiteatacsig)){
  
  ourrow <- whitepeakdistance[i,]
  whitepeakdistancedf[i,"Distance"] <- min(ourrow[ourrow > 0])
  
}

totalpeakdistancedf <- data.frame("Region" = c(gastropeakdistancedf$Region,
                                               heartpeakdistancedf$Region,
                                               hippopeakdistancedf$Region,
                                               kidneypeakdistancedf$Region,
                                               liverpeakdistancedf$Region,
                                               lungpeakdistancedf$Region,
                                               brownpeakdistancedf$Region,
                                               whitepeakdistancedf$Region),
                                  "Distance" = c(gastropeakdistancedf$Distance,
                                                 heartpeakdistancedf$Distance,
                                                 hippopeakdistancedf$Distance,
                                                 kidneypeakdistancedf$Distance,
                                                 liverpeakdistancedf$Distance,
                                                 lungpeakdistancedf$Distance,
                                                 brownpeakdistancedf$Distance,
                                                 whitepeakdistancedf$Distance),
                                  "Tissue" = c(rep("SKM-GN",dim(gastropeakdistancedf)[1]),
                                               rep("HEART",dim(heartpeakdistancedf)[1]),
                                               rep("HIPPOC",dim(hippopeakdistancedf)[1]),
                                               rep("KIDNEY",dim(kidneypeakdistancedf)[1]),
                                               rep("LIVER",dim(liverpeakdistancedf)[1]),
                                               rep("LUNG",dim(lungpeakdistancedf)[1]),
                                               rep("BAT",dim(brownpeakdistancedf)[1]),
                                               rep("WAT-SC",dim(whitepeakdistancedf)[1])))

pdf(file = "Figure 2C.pdf",width = 7,height = 6)
ggplot(totalpeakdistancedf, aes(x = Distance)) +
  geom_histogram(aes(color = Tissue,fill = Tissue),
                 position = "stack", bins = 30) + theme_classic() +
  scale_color_manual(values = c("#8c5220","#f28b2f","#bf7534","#7553a7","#da6c75","#04bf8a","#088c03","#214da6")) + 
  scale_fill_manual(values = c("#8c5220","#f28b2f","#bf7534","#7553a7","#da6c75","#04bf8a","#088c03","#214da6"))  + scale_x_log10(breaks=c(1,100,10000,1000000,100000000)) + xlab("Distance to Nearest DEG TSS") + ylab("Count") + theme(axis.title = element_text(size = 15),axis.text = element_text(size = 12),legend.title = element_text(size = 15),legend.text = element_text(size = 12))
dev.off()

pdf(file = "Figure 2D.pdf",width = 7,height = 6)
ggplot(totalpeakdistancedf, aes(x = Distance)) +
  geom_histogram(aes(color = Region,fill = Region),
                 position = "stack", bins = 30) + theme_classic() +
  scale_color_manual(values = c("#E377C2FF","#D62728FF","#BCBD22FF","#7F7F7FFF","#9467BDFF","#8C564BFF","#2CA02CFF","#FF7F0EFF","#1F77B4FF")) + 
  scale_fill_manual(values = c("#E377C2FF","#D62728FF","#BCBD22FF","#7F7F7FFF","#9467BDFF","#8C564BFF","#2CA02CFF","#FF7F0EFF","#1F77B4FF"))  + scale_x_log10(breaks=c(1,100,10000,1000000,100000000)) + xlab("Distance to Nearest DEG TSS") + ylab("Count")
dev.off()


####
# Figure 2E and Supplemental Figure S6
#####

gastropeakcor <- matrix(0L,nrow = dim(gastropeakdistance)[1],ncol = dim(gastropeakdistance)[2])
rownames(gastropeakcor) <- rownames(gastropeakdistance)
colnames(gastropeakcor) <- colnames(gastropeakdistance)
for(i in 1:dim(gastropeakcor)[1]){
  for(j in 1:dim(gastropeakcor)[2]){
    ourpeak <- rownames(gastropeakdistance)[i]
    ourgene <- colnames(gastropeakdistance)[j]
    gastropeakcor[i,j] <- cor(gastrosigatacl2fc[ourpeak,],gastrol2fcmat[ourgene,])
  }
}


heartpeakcor <- matrix(0L,nrow = dim(heartpeakdistance)[1],ncol = dim(heartpeakdistance)[2])
rownames(heartpeakcor) <- rownames(heartpeakdistance)
colnames(heartpeakcor) <- colnames(heartpeakdistance)
for(i in 1:dim(heartpeakcor)[1]){
  for(j in 1:dim(heartpeakcor)[2]){
    ourpeak <- rownames(heartpeakdistance)[i]
    ourgene <- colnames(heartpeakdistance)[j]
    heartpeakcor[i,j] <- cor(heartsigatacl2fc[ourpeak,],heartl2fcmat[ourgene,])
  }
}


hippopeakcor <- matrix(0L,nrow = dim(hippopeakdistance)[1],ncol = dim(hippopeakdistance)[2])
rownames(hippopeakcor) <- rownames(hippopeakdistance)
colnames(hippopeakcor) <- colnames(hippopeakdistance)
for(i in 1:dim(hippopeakcor)[1]){
  for(j in 1:dim(hippopeakcor)[2]){
    ourpeak <- rownames(hippopeakdistance)[i]
    ourgene <- colnames(hippopeakdistance)[j]
    hippopeakcor[i,j] <- cor(hipposigatacl2fc[ourpeak,],hippol2fcmat[ourgene,])
  }
}


kidneypeakcor <- matrix(0L,nrow = dim(kidneypeakdistance)[1],ncol = dim(kidneypeakdistance)[2])
rownames(kidneypeakcor) <- rownames(kidneypeakdistance)
colnames(kidneypeakcor) <- colnames(kidneypeakdistance)
for(i in 1:dim(kidneypeakcor)[1]){
  for(j in 1:dim(kidneypeakcor)[2]){
    ourpeak <- rownames(kidneypeakdistance)[i]
    ourgene <- colnames(kidneypeakdistance)[j]
    kidneypeakcor[i,j] <- cor(kidneysigatacl2fc[ourpeak,],kidneyl2fcmat[ourgene,])
  }
}


liverpeakcor <- matrix(0L,nrow = dim(liverpeakdistance)[1],ncol = dim(liverpeakdistance)[2])
rownames(liverpeakcor) <- rownames(liverpeakdistance)
colnames(liverpeakcor) <- colnames(liverpeakdistance)
for(i in 1:dim(liverpeakcor)[1]){
  for(j in 1:dim(liverpeakcor)[2]){
    ourpeak <- rownames(liverpeakdistance)[i]
    ourgene <- colnames(liverpeakdistance)[j]
    liverpeakcor[i,j] <- cor(liversigatacl2fc[ourpeak,],liverl2fcmat[ourgene,])
  }
}


lungpeakcor <- matrix(0L,nrow = dim(lungpeakdistance)[1],ncol = dim(lungpeakdistance)[2])
rownames(lungpeakcor) <- rownames(lungpeakdistance)
colnames(lungpeakcor) <- colnames(lungpeakdistance)
for(i in 1:dim(lungpeakcor)[1]){
  for(j in 1:dim(lungpeakcor)[2]){
    ourpeak <- rownames(lungpeakdistance)[i]
    ourgene <- colnames(lungpeakdistance)[j]
    lungpeakcor[i,j] <- cor(lungsigatacl2fc[ourpeak,],lungl2fcmat[ourgene,])
  }
}
lungpeakcor[is.na(lungpeakcor)] <- 0

brownpeakcor <- matrix(0L,nrow = dim(brownpeakdistance)[1],ncol = dim(brownpeakdistance)[2])
rownames(brownpeakcor) <- rownames(brownpeakdistance)
colnames(brownpeakcor) <- colnames(brownpeakdistance)
for(i in 1:dim(brownpeakcor)[1]){
  for(j in 1:dim(brownpeakcor)[2]){
    ourpeak <- rownames(brownpeakdistance)[i]
    ourgene <- colnames(brownpeakdistance)[j]
    brownpeakcor[i,j] <- cor(brownsigatacl2fc[ourpeak,],brownl2fcmat[ourgene,])
  }
}


whitepeakcor <- matrix(0L,nrow = dim(whitepeakdistance)[1],ncol = dim(whitepeakdistance)[2])
rownames(whitepeakcor) <- rownames(whitepeakdistance)
colnames(whitepeakcor) <- colnames(whitepeakdistance)
for(i in 1:dim(whitepeakcor)[1]){
  for(j in 1:dim(whitepeakcor)[2]){
    ourpeak <- rownames(whitepeakdistance)[i]
    ourgene <- colnames(whitepeakdistance)[j]
    whitepeakcor[i,j] <- cor(whitesigatacl2fc[ourpeak,],whitel2fcmat[ourgene,])
  }
}
whitepeakcor[is.na(whitepeakcor)] <- 0

gastropeaknoabsdistance <- matrix(0L,nrow = length(gastroatacsig),ncol = length(fingastrornasig))
rownames(gastropeaknoabsdistance) <- gastroatacsig
colnames(gastropeaknoabsdistance) <- fingastrornasig

for(i in 1:length(gastroatacsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(fingastrornasig)){
    gastropeaknoabsdistance[i,j] <- (gastrornasiganno$chrom[j] == gastroatacsiganno$chrom[i])*((gastrornasiganno$start[j] - gastroatacsiganno$mid[i]))
  }
}

heartpeaknoabsdistance <- matrix(0L,nrow = length(heartatacsig),ncol = length(finheartrnasig))
rownames(heartpeaknoabsdistance) <- heartatacsig
colnames(heartpeaknoabsdistance) <- finheartrnasig

for(i in 1:length(heartatacsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(finheartrnasig)){
    heartpeaknoabsdistance[i,j] <- (heartrnasiganno$chrom[j] == heartatacsiganno$chrom[i])*((heartrnasiganno$start[j] - heartatacsiganno$mid[i]))
  }
}

hippopeaknoabsdistance <- matrix(0L,nrow = length(hippoatacsig),ncol = length(finhippornasig))
rownames(hippopeaknoabsdistance) <- hippoatacsig
colnames(hippopeaknoabsdistance) <- finhippornasig

for(i in 1:length(hippoatacsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(finhippornasig)){
    hippopeaknoabsdistance[i,j] <- (hippornasiganno$chrom[j] == hippoatacsiganno$chrom[i])*((hippornasiganno$start[j] - hippoatacsiganno$mid[i]))
  }
}

kidneypeaknoabsdistance <- matrix(0L,nrow = length(kidneyatacsig),ncol = length(finkidneyrnasig))
rownames(kidneypeaknoabsdistance) <- kidneyatacsig
colnames(kidneypeaknoabsdistance) <- finkidneyrnasig

for(i in 1:length(kidneyatacsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(finkidneyrnasig)){
    kidneypeaknoabsdistance[i,j] <- (kidneyrnasiganno$chrom[j] == kidneyatacsiganno$chrom[i])*((kidneyrnasiganno$start[j] - kidneyatacsiganno$mid[i]))
  }
}

liverpeaknoabsdistance <- matrix(0L,nrow = length(liveratacsig),ncol = length(finliverrnasig))
rownames(liverpeaknoabsdistance) <- liveratacsig
colnames(liverpeaknoabsdistance) <- finliverrnasig

for(i in 1:length(liveratacsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(finliverrnasig)){
    liverpeaknoabsdistance[i,j] <- (liverrnasiganno$chrom[j] == liveratacsiganno$chrom[i])*((liverrnasiganno$start[j] - liveratacsiganno$mid[i]))
  }
}

lungpeaknoabsdistance <- matrix(0L,nrow = length(lungatacsig),ncol = length(finlungrnasig))
rownames(lungpeaknoabsdistance) <- lungatacsig
colnames(lungpeaknoabsdistance) <- finlungrnasig

for(i in 1:length(lungatacsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(finlungrnasig)){
    lungpeaknoabsdistance[i,j] <- (lungrnasiganno$chrom[j] == lungatacsiganno$chrom[i])*((lungrnasiganno$start[j] - lungatacsiganno$mid[i]))
  }
}

brownpeaknoabsdistance <- matrix(0L,nrow = length(brownatacsig),ncol = length(finbrownrnasig))
rownames(brownpeaknoabsdistance) <- brownatacsig
colnames(brownpeaknoabsdistance) <- finbrownrnasig

for(i in 1:length(brownatacsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(finbrownrnasig)){
    brownpeaknoabsdistance[i,j] <- (brownrnasiganno$chrom[j] == brownatacsiganno$chrom[i])*((brownrnasiganno$start[j] - brownatacsiganno$mid[i]))
  }
}

whitepeaknoabsdistance <- matrix(0L,nrow = length(whiteatacsig),ncol = length(finwhiternasig))
rownames(whitepeaknoabsdistance) <- whiteatacsig
colnames(whitepeaknoabsdistance) <- finwhiternasig

for(i in 1:length(whiteatacsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(finwhiternasig)){
    whitepeaknoabsdistance[i,j] <- (whiternasiganno$chrom[j] == whiteatacsiganno$chrom[i])*((whiternasiganno$start[j] - whiteatacsiganno$mid[i]))
  }
}

gastrodistvscordf <- data.frame("Correlation" = as.vector(gastropeakcor),"Distance" = as.vector(gastropeaknoabsdistance))
heartdistvscordf <- data.frame("Correlation" = as.vector(heartpeakcor),"Distance" = as.vector(heartpeaknoabsdistance))
hippodistvscordf <- data.frame("Correlation" = as.vector(hippopeakcor),"Distance" = as.vector(hippopeaknoabsdistance))
kidneydistvscordf <- data.frame("Correlation" = as.vector(kidneypeakcor),"Distance" = as.vector(kidneypeaknoabsdistance))
liverdistvscordf <- data.frame("Correlation" = as.vector(liverpeakcor),"Distance" = as.vector(liverpeaknoabsdistance))
lungdistvscordf <- data.frame("Correlation" = as.vector(lungpeakcor),"Distance" = as.vector(lungpeaknoabsdistance))
browndistvscordf <- data.frame("Correlation" = as.vector(brownpeakcor),"Distance" = as.vector(brownpeaknoabsdistance))
whitedistvscordf <- data.frame("Correlation" = as.vector(whitepeakcor),"Distance" = as.vector(whitepeaknoabsdistance))

totaldistvscordf <- rbind(gastrodistvscordf,heartdistvscordf,hippodistvscordf,kidneydistvscordf,
                          liverdistvscordf,lungdistvscordf,whitedistvscordf,browndistvscordf)
totaldistvscordf$Tissue <- c(rep("SKM-GN",dim(gastrodistvscordf)[1]),
                             rep("HEART",dim(heartdistvscordf)[1]),
                             rep("HIPPOC",dim(hippodistvscordf)[1]),
                             rep("KIDNEY",dim(kidneydistvscordf)[1]),
                             rep("LIVER",dim(liverdistvscordf)[1]),
                             rep("LUNG",dim(lungdistvscordf)[1]),
                             rep("WAT-SC",dim(whitedistvscordf)[1]),
                             rep("BAT",dim(browndistvscordf)[1]))

# Figure 2E
pdf(file = "Figure 2E.pdf",width = 9,height = 6)
ggplot(totaldistvscordf[abs(totaldistvscordf$Distance) > 0 & abs(totaldistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_density2d_filled() + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 18),legend.title = element_text(size = 18),legend.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + geom_point(data = totaldistvscordf[abs(totaldistvscordf$Distance) > 0 & abs(totaldistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation,color=Tissue),size = 2) + scale_color_manual(values = ann_colsheatmaptrim$Tissue)
dev.off()

# Supplemental Figure S6

pdf(file = "Supplemental Figure S6A.pdf",width = 6,height = 6)
ggplot(gastrodistvscordf[abs(gastrodistvscordf$Distance) > 0 & abs(gastrodistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_colsheatmaptrim$Tissue["SKM-GN"]) + geom_density2d(color = ann_colsheatmaptrim$Tissue["SKM-GN"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 18)) + xlim(-500000,500000) + ylim(-1,1)
dev.off()

pdf(file = "Supplemental Figure S6B.pdf",width = 6,height = 6)
ggplot(heartdistvscordf[abs(heartdistvscordf$Distance) > 0 & abs(heartdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_colsheatmaptrim$Tissue["HEART"]) + geom_density2d(color = ann_colsheatmaptrim$Tissue["HEART"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 18)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_colsheatmaptrim$Tissue["HEART"])
dev.off()

pdf(file = "Supplemental Figure S6C.pdf",width = 6,height = 6)
ggplot(hippodistvscordf[abs(hippodistvscordf$Distance) > 0 & abs(hippodistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_colsheatmaptrim$Tissue["HIPPOC"]) + geom_density2d(color = ann_colsheatmaptrim$Tissue["HIPPOC"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 18)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_colsheatmaptrim$Tissue["HIPPOC"])
dev.off()

pdf(file = "Supplemental Figure S6D.pdf",width = 6,height = 6)
ggplot(kidneydistvscordf[abs(kidneydistvscordf$Distance) > 0 & abs(kidneydistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_colsheatmaptrim$Tissue["KIDNEY"]) + geom_density2d(color = ann_colsheatmaptrim$Tissue["KIDNEY"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 18)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_colsheatmaptrim$Tissue["KIDNEY"])
dev.off()

pdf(file = "Supplemental Figure S6E.pdf",width = 6,height = 6)
ggplot(liverdistvscordf[abs(liverdistvscordf$Distance) > 0 & abs(liverdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_colsheatmaptrim$Tissue["LIVER"]) + geom_density2d(color = ann_colsheatmaptrim$Tissue["LIVER"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 18)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_colsheatmaptrim$Tissue["LIVER"])
dev.off()

pdf(file = "Supplemental Figure S6F.pdf",width = 6,height = 6)
ggplot(lungdistvscordf[abs(lungdistvscordf$Distance) > 0 & abs(lungdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_colsheatmaptrim$Tissue["LUNG"]) + geom_density2d(color = ann_colsheatmaptrim$Tissue["LUNG"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 18)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_colsheatmaptrim$Tissue["LUNG"])
dev.off()

pdf(file = "Supplemental Figure S6G.pdf",width = 6,height = 6)
ggplot(browndistvscordf[abs(browndistvscordf$Distance) > 0 & abs(browndistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_colsheatmaptrim$Tissue["BAT"]) + geom_density2d(color = ann_colsheatmaptrim$Tissue["BAT"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 18)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_colsheatmaptrim$Tissue["BAT"])
dev.off()

pdf(file = "Supplemental Figure S6H.pdf",width = 6,height = 6)
ggplot(whitedistvscordf[abs(whitedistvscordf$Distance) > 0 & abs(whitedistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_colsheatmaptrim$Tissue["WAT-SC"]) + geom_density2d(color = ann_colsheatmaptrim$Tissue["WAT-SC"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 18)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_colsheatmaptrim$Tissue["WAT-SC"])
dev.off()

####
# Figure 3
#####

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

####
# Figure 3A - To be added
#####

gastro_gostres <- gost(query = colnames(gastropeakcortrim),organism = "rnorvegicus",sources = c('GO:BP','REAC','WP','KEGG'),custom_bg = rownames(gastrornanorm),significant = F)
heart_gostres <- gost(query = colnames(heartpeakcortrim),organism = "rnorvegicus",sources = c('GO:BP','REAC','WP','KEGG'),custom_bg = rownames(gastrornanorm),significant = F)
kidney_gostres <- gost(query = colnames(kidneypeakcortrim),organism = "rnorvegicus",sources = c('GO:BP','REAC','WP','KEGG'),custom_bg = rownames(kidneyrnanorm),significant = F)
liver_gostres <- gost(query = colnames(liverpeakcortrim),organism = "rnorvegicus",sources = c('GO:BP','REAC','WP','KEGG'),custom_bg = rownames(liverrnanorm),significant = F)
lung_gostres <- gost(query = colnames(lungpeakcortrim),organism = "rnorvegicus",sources = c('GO:BP','REAC','WP','KEGG'),custom_bg = rownames(lungrnanorm),significant = F)
brown_gostres <- gost(query = colnames(brownpeakcortrim),organism = "rnorvegicus",sources = c('GO:BP','REAC','WP','KEGG'),custom_bg = rownames(brownrnanorm),significant = F)

#gastro_gostres$result$padj <- p.adjust(gastro_gostres$result$p_value)
#kidney_gostres$result$padj <- p.adjust(kidney_gostres$result$p_value)
#liver_gostres$result$padj <- p.adjust(liver_gostres$result$p_value)
#lung_gostres$result$padj <- p.adjust(lung_gostres$result$p_value)
#brown_gostres$result$padj <- p.adjust(brown_gostres$result$p_value)

gastrow1sig_gostres <- gost(query = transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "SKM-GN" &
                                                                      transcript_rna_seq$timewise_dea$comparison_group %in% "1w" &
                                                                      transcript_rna_seq$timewise_dea$adj_p_value < 0.05 &
                                                                      abs(transcript_rna_seq$timewise_dea$logFC) > 0.1,"feature_ID"][!is.na(transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "SKM-GN" &
                                                                                                                                                                              transcript_rna_seq$timewise_dea$comparison_group %in% "1w" &
                                                                                                                                                                              transcript_rna_seq$timewise_dea$adj_p_value < 0.05 &
                                                                                                                                                                              abs(transcript_rna_seq$timewise_dea$logFC) > 0.1,"feature_ID"])],
                            organism = "rnorvegicus",sources = c('GO:BP','REAC','WP','KEGG'),custom_bg = rownames(gastrornanorm),significant = F)
gastrow8sig_gostres <- gost(query = transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "SKM-GN" &
                                                                      transcript_rna_seq$timewise_dea$comparison_group %in% "8w" &
                                                                      transcript_rna_seq$timewise_dea$adj_p_value < 0.05 &
                                                                      abs(transcript_rna_seq$timewise_dea$logFC) > 0.1,"feature_ID"][!is.na(transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "SKM-GN" &
                                                                                                                                                                              transcript_rna_seq$timewise_dea$comparison_group %in% "8w" &
                                                                                                                                                                              transcript_rna_seq$timewise_dea$adj_p_value < 0.05 &
                                                                                                                                                                              abs(transcript_rna_seq$timewise_dea$logFC) > 0.1,"feature_ID"])],
                            organism = "rnorvegicus",sources = c('GO:BP','REAC','WP','KEGG'),custom_bg = rownames(gastrornanorm),significant = F)
heartw1sig_gostres <- gost(query = transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "HEART" &
                                                                      transcript_rna_seq$timewise_dea$comparison_group %in% "1w" &
                                                                      transcript_rna_seq$timewise_dea$adj_p_value < 0.05 &
                                                                      abs(transcript_rna_seq$timewise_dea$logFC) > 0.1,"feature_ID"][!is.na(transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "HEART" &
                                                                                                                                                                              transcript_rna_seq$timewise_dea$comparison_group %in% "1w" &
                                                                                                                                                                              transcript_rna_seq$timewise_dea$adj_p_value < 0.05 &
                                                                                                                                                                              abs(transcript_rna_seq$timewise_dea$logFC) > 0.1,"feature_ID"])],
                            organism = "rnorvegicus",sources = c('GO:BP','REAC','WP','KEGG'),custom_bg = rownames(heartrnanorm),significant = F)
heartw8sig_gostres <- gost(query = transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "HEART" &
                                                                      transcript_rna_seq$timewise_dea$comparison_group %in% "8w" &
                                                                      transcript_rna_seq$timewise_dea$adj_p_value < 0.05 &
                                                                      abs(transcript_rna_seq$timewise_dea$logFC) > 0.1,"feature_ID"][!is.na(transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "HEART" &
                                                                                                                                                                              transcript_rna_seq$timewise_dea$comparison_group %in% "8w" &
                                                                                                                                                                              transcript_rna_seq$timewise_dea$adj_p_value < 0.05 &
                                                                                                                                                                              abs(transcript_rna_seq$timewise_dea$logFC) > 0.1,"feature_ID"])],
                            organism = "rnorvegicus",sources = c('GO:BP','REAC','WP','KEGG'),custom_bg = rownames(heartrnanorm),significant = F)

hippow1sig_gostres <- gost(query = transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "HIPPOC" &
                                                                     transcript_rna_seq$timewise_dea$comparison_group %in% "1w" &
                                                                     transcript_rna_seq$timewise_dea$adj_p_value < 0.05 &
                                                                     abs(transcript_rna_seq$timewise_dea$logFC) > 0.1,"feature_ID"][!is.na(transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "HIPPOC" &
                                                                                                                                                                             transcript_rna_seq$timewise_dea$comparison_group %in% "1w" &
                                                                                                                                                                             transcript_rna_seq$timewise_dea$adj_p_value < 0.05 &
                                                                                                                                                                             abs(transcript_rna_seq$timewise_dea$logFC) > 0.1,"feature_ID"])],
                           organism = "rnorvegicus",sources = c('GO:BP','REAC','WP','KEGG'),custom_bg = rownames(hippornanorm),significant = F)
hippow8sig_gostres <- gost(query = transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "HIPPOC" &
                                                                     transcript_rna_seq$timewise_dea$comparison_group %in% "8w" &
                                                                     transcript_rna_seq$timewise_dea$adj_p_value < 0.05 &
                                                                     abs(transcript_rna_seq$timewise_dea$logFC) > 0.1,"feature_ID"][!is.na(transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "HIPPOC" &
                                                                                                                                                                             transcript_rna_seq$timewise_dea$comparison_group %in% "8w" &
                                                                                                                                                                             transcript_rna_seq$timewise_dea$adj_p_value < 0.05 &
                                                                                                                                                                             abs(transcript_rna_seq$timewise_dea$logFC) > 0.1,"feature_ID"])],
                           organism = "rnorvegicus",sources = c('GO:BP','REAC','WP','KEGG'),custom_bg = rownames(hippornanorm),significant = F)

kidneyw1sig_gostres <- gost(query = transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "KIDNEY" &
                                                                     transcript_rna_seq$timewise_dea$comparison_group %in% "1w" &
                                                                     transcript_rna_seq$timewise_dea$adj_p_value < 0.05 &
                                                                     abs(transcript_rna_seq$timewise_dea$logFC) > 0.1,"feature_ID"][!is.na(transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "KIDNEY" &
                                                                                                                                                                             transcript_rna_seq$timewise_dea$comparison_group %in% "1w" &
                                                                                                                                                                             transcript_rna_seq$timewise_dea$adj_p_value < 0.05 &
                                                                                                                                                                             abs(transcript_rna_seq$timewise_dea$logFC) > 0.1,"feature_ID"])],
                           organism = "rnorvegicus",sources = c('GO:BP','REAC','WP','KEGG'),custom_bg = rownames(kidneyrnanorm),significant = F)
kidneyw8sig_gostres <- gost(query = transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "KIDNEY" &
                                                                     transcript_rna_seq$timewise_dea$comparison_group %in% "8w" &
                                                                     transcript_rna_seq$timewise_dea$adj_p_value < 0.05 &
                                                                     abs(transcript_rna_seq$timewise_dea$logFC) > 0.1,"feature_ID"][!is.na(transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "KIDNEY" &
                                                                                                                                                                             transcript_rna_seq$timewise_dea$comparison_group %in% "8w" &
                                                                                                                                                                             transcript_rna_seq$timewise_dea$adj_p_value < 0.05 &
                                                                                                                                                                             abs(transcript_rna_seq$timewise_dea$logFC) > 0.1,"feature_ID"])],
                           organism = "rnorvegicus",sources = c('GO:BP','REAC','WP','KEGG'),custom_bg = rownames(kidneyrnanorm),significant = F)


liverw1sig_gostres <- gost(query = transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "LIVER" &
                                                                      transcript_rna_seq$timewise_dea$comparison_group %in% "1w" &
                                                                      transcript_rna_seq$timewise_dea$adj_p_value < 0.05 &
                                                                      abs(transcript_rna_seq$timewise_dea$logFC) > 0.1,"feature_ID"][!is.na(transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "LIVER" &
                                                                                                                                                                              transcript_rna_seq$timewise_dea$comparison_group %in% "1w" &
                                                                                                                                                                              transcript_rna_seq$timewise_dea$adj_p_value < 0.05 &
                                                                                                                                                                              abs(transcript_rna_seq$timewise_dea$logFC) > 0.1,"feature_ID"])],
                            organism = "rnorvegicus",sources = c('GO:BP','REAC','WP','KEGG'),custom_bg = rownames(liverrnanorm),significant = F)
liverw8sig_gostres <- gost(query = transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "LIVER" &
                                                                      transcript_rna_seq$timewise_dea$comparison_group %in% "8w" &
                                                                      transcript_rna_seq$timewise_dea$adj_p_value < 0.05 &
                                                                      abs(transcript_rna_seq$timewise_dea$logFC) > 0.1,"feature_ID"][!is.na(transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "LIVER" &
                                                                                                                                                                              transcript_rna_seq$timewise_dea$comparison_group %in% "8w" &
                                                                                                                                                                              transcript_rna_seq$timewise_dea$adj_p_value < 0.05 &
                                                                                                                                                                              abs(transcript_rna_seq$timewise_dea$logFC) > 0.1,"feature_ID"])],
                            organism = "rnorvegicus",sources = c('GO:BP','REAC','WP','KEGG'),custom_bg = rownames(liverrnanorm),significant = F)

lungw1sig_gostres <- gost(query = transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "LUNG" &
                                                                     transcript_rna_seq$timewise_dea$comparison_group %in% "1w" &
                                                                     transcript_rna_seq$timewise_dea$adj_p_value < 0.05 &
                                                                     abs(transcript_rna_seq$timewise_dea$logFC) > 0.1,"feature_ID"][!is.na(transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "LUNG" &
                                                                                                                                                                             transcript_rna_seq$timewise_dea$comparison_group %in% "1w" &
                                                                                                                                                                             transcript_rna_seq$timewise_dea$adj_p_value < 0.05 &
                                                                                                                                                                             abs(transcript_rna_seq$timewise_dea$logFC) > 0.1,"feature_ID"])],
                           organism = "rnorvegicus",sources = c('GO:BP','REAC','WP','KEGG'),custom_bg = rownames(lungrnanorm),significant = F)
lungw8sig_gostres <- gost(query = transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "LUNG" &
                                                                     transcript_rna_seq$timewise_dea$comparison_group %in% "8w" &
                                                                     transcript_rna_seq$timewise_dea$adj_p_value < 0.05 &
                                                                     abs(transcript_rna_seq$timewise_dea$logFC) > 0.1,"feature_ID"][!is.na(transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "LUNG" &
                                                                                                                                                                             transcript_rna_seq$timewise_dea$comparison_group %in% "8w" &
                                                                                                                                                                             transcript_rna_seq$timewise_dea$adj_p_value < 0.05 &
                                                                                                                                                                             abs(transcript_rna_seq$timewise_dea$logFC) > 0.1,"feature_ID"])],
                           organism = "rnorvegicus",sources = c('GO:BP','REAC','WP','KEGG'),custom_bg = rownames(lungrnanorm),significant = F)

brownw1sig_gostres <- gost(query = transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "BAT" &
                                                                    transcript_rna_seq$timewise_dea$comparison_group %in% "1w" &
                                                                    transcript_rna_seq$timewise_dea$adj_p_value < 0.05 &
                                                                    abs(transcript_rna_seq$timewise_dea$logFC) > 0.1,"feature_ID"][!is.na(transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "BAT" &
                                                                                                                                                                            transcript_rna_seq$timewise_dea$comparison_group %in% "1w" &
                                                                                                                                                                            transcript_rna_seq$timewise_dea$adj_p_value < 0.05 &
                                                                                                                                                                            abs(transcript_rna_seq$timewise_dea$logFC) > 0.1,"feature_ID"])],
                          organism = "rnorvegicus",sources = c('GO:BP','REAC','WP','KEGG'),custom_bg = rownames(brownrnanorm),significant = F)
brownw8sig_gostres <- gost(query = transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "BAT" &
                                                                    transcript_rna_seq$timewise_dea$comparison_group %in% "8w" &
                                                                    transcript_rna_seq$timewise_dea$adj_p_value < 0.05 &
                                                                    abs(transcript_rna_seq$timewise_dea$logFC) > 0.1,"feature_ID"][!is.na(transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "BAT" &
                                                                                                                                                                            transcript_rna_seq$timewise_dea$comparison_group %in% "8w" &
                                                                                                                                                                            transcript_rna_seq$timewise_dea$adj_p_value < 0.05 &
                                                                                                                                                                            abs(transcript_rna_seq$timewise_dea$logFC) > 0.1,"feature_ID"])],
                          organism = "rnorvegicus",sources = c('GO:BP','REAC','WP','KEGG'),custom_bg = rownames(brownrnanorm),significant = F)

whitew1sig_gostres <- gost(query = transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "WAT-SC" &
                                                                     transcript_rna_seq$timewise_dea$comparison_group %in% "1w" &
                                                                     transcript_rna_seq$timewise_dea$adj_p_value < 0.05 &
                                                                     abs(transcript_rna_seq$timewise_dea$logFC) > 0.1,"feature_ID"][!is.na(transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "WAT-SC" &
                                                                                                                                                                             transcript_rna_seq$timewise_dea$comparison_group %in% "1w" &
                                                                                                                                                                             transcript_rna_seq$timewise_dea$adj_p_value < 0.05 &
                                                                                                                                                                             abs(transcript_rna_seq$timewise_dea$logFC) > 0.1,"feature_ID"])],
                           organism = "rnorvegicus",sources = c('GO:BP','REAC','WP','KEGG'),custom_bg = rownames(whiternanorm),significant = F)
whitew8sig_gostres <- gost(query = transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "WAT-SC" &
                                                                     transcript_rna_seq$timewise_dea$comparison_group %in% "8w" &
                                                                     transcript_rna_seq$timewise_dea$adj_p_value < 0.05 &
                                                                     abs(transcript_rna_seq$timewise_dea$logFC) > 0.1,"feature_ID"][!is.na(transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "WAT-SC" &
                                                                                                                                                                             transcript_rna_seq$timewise_dea$comparison_group %in% "8w" &
                                                                                                                                                                             transcript_rna_seq$timewise_dea$adj_p_value < 0.05 &
                                                                                                                                                                             abs(transcript_rna_seq$timewise_dea$logFC) > 0.1,"feature_ID"])],
                           organism = "rnorvegicus",sources = c('GO:BP','REAC','WP','KEGG'),custom_bg = rownames(whiternanorm),significant = F)








####
# Figure 3B-D
#####

pdf(file = "Figure 3B.pdf",width = 8,height = 5)
pheatmap(gastro50peakcorfin[unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin) & gastro50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),colnames(gastro50peakcorfin[unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin) & gastro50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),colSums(gastro50peakcorfin[unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin) & gastro50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),]) != 0])],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_col = enstosym[colnames(gastro50peakcorfin[unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin) & gastro50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),colSums(gastro50peakcorfin[unique(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% rownames(gastro50peakcorfin) & gastro50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),]) != 0]),"Symbol"],fontsize = 20,cluster_rows = F, cluster_cols = F)
dev.off()

pdf(file = "Figure 3C.pdf",width = 8,height = 5)
pheatmap(lung50peakcorfin[unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin) & lung50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),colnames(lung50peakcorfin[unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin) & lung50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),colSums(lung50peakcorfin[unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin) & lung50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),]) != 0])],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_col = enstosym[colnames(lung50peakcorfin[unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin) & lung50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),colSums(lung50peakcorfin[unique(lung50peakmotifs[lung50peakmotifs$PositionID %in% rownames(lung50peakcorfin) & lung50peakmotifs$Motif.Name %in% "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer","PositionID"]),]) != 0]),"Symbol"],fontsize = 20,cluster_rows = F, cluster_cols = F)
dev.off()

pdf(file = "Figure 3D.pdf",width = 9,height = 11)
pheatmap(liver50peakcorfin[unique(liver50peakmotifs[liver50peakmotifs$PositionID %in% rownames(liver50peakcorfin) & liver50peakmotifs$Motif.Name %in% "Smad3(MAD)/NPC-Smad3-ChIP-Seq(GSE36673)/Homer","PositionID"]),colnames(liver50peakcorfin[unique(liver50peakmotifs[liver50peakmotifs$PositionID %in% rownames(liver50peakcorfin) & liver50peakmotifs$Motif.Name %in% "Smad3(MAD)/NPC-Smad3-ChIP-Seq(GSE36673)/Homer","PositionID"]),colSums(liver50peakcorfin[unique(liver50peakmotifs[liver50peakmotifs$PositionID %in% rownames(liver50peakcorfin) & liver50peakmotifs$Motif.Name %in% "Smad3(MAD)/NPC-Smad3-ChIP-Seq(GSE36673)/Homer","PositionID"]),]) != 0])],angle_col = 315,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_col = enstosym[colnames(liver50peakcorfin[unique(liver50peakmotifs[liver50peakmotifs$PositionID %in% rownames(liver50peakcorfin) & liver50peakmotifs$Motif.Name %in% "Smad3(MAD)/NPC-Smad3-ChIP-Seq(GSE36673)/Homer","PositionID"]),colSums(liver50peakcorfin[unique(liver50peakmotifs[liver50peakmotifs$PositionID %in% rownames(liver50peakcorfin) & liver50peakmotifs$Motif.Name %in% "Smad3(MAD)/NPC-Smad3-ChIP-Seq(GSE36673)/Homer","PositionID"]),]) != 0]),"Symbol"],fontsize = 20,cluster_rows = F, cluster_cols = F)
dev.off()

####
# Supplemental Table 1
#####

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

####
# Figure 3E-H
#####

ourgene <- enstosym[enstosym$Symbol %in% "Igf2","Ensembl"]
igf2_gastro_corpeaks <- rownames(gastropeakcortrim)[abs(gastropeakcortrim[,ourgene]) > 0.5]
ourpeak <- igf2_gastro_corpeaks
ourdf <- data.frame("RNA.L2FC" = gastrol2fcmat[ourgene,],
                    "ATAC.L2FC" = gastrosigatacl2fc[ourpeak,],
                    "Label" = colnames(gastrol2fcmat))
pdf(file = "Figure 3E.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("SKM-GN: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Sall2","Ensembl"]
sall2_gastro_corpeaks <- rownames(gastropeakcortrim)[abs(gastropeakcortrim[,ourgene]) > 0.5]
ourpeak <- sall2_gastro_corpeaks
ourdf <- data.frame("RNA.L2FC" = gastrol2fcmat[ourgene,],
                    "ATAC.L2FC" = gastrosigatacl2fc[ourpeak,],
                    "Label" = colnames(gastrol2fcmat))

pdf(file = "Figure 3F.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("SKM-GN: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Nfkb2","Ensembl"]
nfkb2_lung_corpeaks <- rownames(lungpeakcortrim)[abs(lungpeakcortrim[,ourgene]) > 0.5]
ourpeak <- nfkb2_lung_corpeaks
ourdf <- data.frame("RNA.L2FC" = lungl2fcmat[ourgene,],
                    "ATAC.L2FC" = lungsigatacl2fc[ourpeak,],
                    "Label" = colnames(lungl2fcmat))
pdf(file = "Figure 3G.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LUNG: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Fkbp4","Ensembl"]
fkbp4_liver_corpeaks <- rownames(liverpeakcortrim)[abs(liverpeakcortrim[,ourgene]) > 0.5]
ourpeak <- fkbp4_liver_corpeaks
ourdf <- data.frame("RNA.L2FC" = liverl2fcmat[ourgene,],
                    "ATAC.L2FC" = liversigatacl2fc[ourpeak,],
                    "Label" = colnames(liverl2fcmat))
pdf(file = "Figure 3H.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

####
# Supplemental Figure S7
#####

ourgene <- enstosym[enstosym$Symbol %in% "Glul","Ensembl"]
glul_liver_corpeaks <- rownames(liverpeakcortrim)[abs(liverpeakcortrim[,ourgene]) > 0.5]
ourpeak <- glul_liver_corpeaks[2]
ourdf <- data.frame("RNA.L2FC" = liverl2fcmat[ourgene,],
                    "ATAC.L2FC" = liversigatacl2fc[ourpeak,],
                    "Label" = colnames(liverl2fcmat))
pdf(file = "Supplemental Figure S7A.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Glul","Ensembl"]
glul_liver_corpeaks <- rownames(liverpeakcortrim)[abs(liverpeakcortrim[,ourgene]) > 0.5]
ourpeak <- glul_liver_corpeaks[1]
ourdf <- data.frame("RNA.L2FC" = liverl2fcmat[ourgene,],
                    "ATAC.L2FC" = liversigatacl2fc[ourpeak,],
                    "Label" = colnames(liverl2fcmat))
pdf(file = "Supplemental Figure S7B.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Lpar3","Ensembl"]
lpar3_liver_corpeaks <- rownames(liverpeakcortrim)[abs(liverpeakcortrim[,ourgene]) > 0.5]
ourpeak <- lpar3_liver_corpeaks[5]
ourdf <- data.frame("RNA.L2FC" = liverl2fcmat[ourgene,],
                    "ATAC.L2FC" = liversigatacl2fc[ourpeak,],
                    "Label" = colnames(liverl2fcmat))
pdf(file = "Supplemental Figure S7C.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Lpar3","Ensembl"]
lpar3_liver_corpeaks <- rownames(liverpeakcortrim)[abs(liverpeakcortrim[,ourgene]) > 0.5]
ourpeak <- lpar3_liver_corpeaks[6]
ourdf <- data.frame("RNA.L2FC" = liverl2fcmat[ourgene,],
                    "ATAC.L2FC" = liversigatacl2fc[ourpeak,],
                    "Label" = colnames(liverl2fcmat))
pdf(file = "Supplemental Figure S7D.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Serpina4","Ensembl"]
serpina4_liver_corpeaks <- rownames(liverpeakcortrim)[abs(liverpeakcortrim[,ourgene]) > 0.5]
ourpeak <- serpina4_liver_corpeaks[3]
ourdf <- data.frame("RNA.L2FC" = liverl2fcmat[ourgene,],
                    "ATAC.L2FC" = liversigatacl2fc[ourpeak,],
                    "Label" = colnames(liverl2fcmat))
pdf(file = "Supplemental Figure S7E.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Serpina4","Ensembl"]
serpina4_liver_corpeaks <- rownames(liverpeakcortrim)[abs(liverpeakcortrim[,ourgene]) > 0.5]
ourpeak <- serpina4_liver_corpeaks[4]
ourdf <- data.frame("RNA.L2FC" = liverl2fcmat[ourgene,],
                    "ATAC.L2FC" = liversigatacl2fc[ourpeak,],
                    "Label" = colnames(liverl2fcmat))
pdf(file = "Supplemental Figure S7F.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Abhd2","Ensembl"]
abhd2_liver_corpeaks <- rownames(liverpeakcortrim)[abs(liverpeakcortrim[,ourgene]) > 0.5]
ourpeak <- abhd2_liver_corpeaks[5]
ourdf <- data.frame("RNA.L2FC" = liverl2fcmat[ourgene,],
                    "ATAC.L2FC" = liversigatacl2fc[ourpeak,],
                    "Label" = colnames(liverl2fcmat))
pdf(file = "Supplemental Figure S7G.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Onecut1","Ensembl"]
onecut1_liver_corpeaks <- rownames(liverpeakcortrim)[abs(liverpeakcortrim[,ourgene]) > 0.5]
ourpeak <- onecut1_liver_corpeaks
ourdf <- data.frame("RNA.L2FC" = liverl2fcmat[ourgene,],
                    "ATAC.L2FC" = liversigatacl2fc[ourpeak,],
                    "Label" = colnames(liverl2fcmat))
pdf(file = "Supplemental Figure S7H.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Ccnd1","Ensembl"]
ccnd1_liver_corpeaks <- rownames(liverpeakcortrim)[abs(liverpeakcortrim[,ourgene]) > 0.5]
ourpeak <- ccnd1_liver_corpeaks[2]
ourdf <- data.frame("RNA.L2FC" = liverl2fcmat[ourgene,],
                    "ATAC.L2FC" = liversigatacl2fc[ourpeak,],
                    "Label" = colnames(liverl2fcmat))
pdf(file = "Supplemental Figure S7I.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Xbp1","Ensembl"]
xbp1_liver_corpeaks <- rownames(liverpeakcortrim)[abs(liverpeakcortrim[,ourgene]) > 0.5]
ourpeak <- xbp1_liver_corpeaks[1]
ourdf <- data.frame("RNA.L2FC" = liverl2fcmat[ourgene,],
                    "ATAC.L2FC" = liversigatacl2fc[ourpeak,],
                    "Label" = colnames(liverl2fcmat))
pdf(file = "Supplemental Figure S7J.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Egfr","Ensembl"]
egfr_liver_corpeaks <- rownames(liverpeakcortrim)[abs(liverpeakcortrim[,ourgene]) > 0.5]
ourpeak <- egfr_liver_corpeaks[3]
ourdf <- data.frame("RNA.L2FC" = liverl2fcmat[ourgene,],
                    "ATAC.L2FC" = liversigatacl2fc[ourpeak,],
                    "Label" = colnames(liverl2fcmat))
pdf(file = "Supplemental Figure S7K.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Fkbp5","Ensembl"]
fkbp5_gastro_corpeaks <- rownames(gastropeakcortrim)[abs(gastropeakcortrim[,ourgene]) > 0.5]
ourpeak <- fkbp5_gastro_corpeaks
ourdf <- data.frame("RNA.L2FC" = gastrol2fcmat[ourgene,],
                    "ATAC.L2FC" = gastrosigatacl2fc[ourpeak,],
                    "Label" = colnames(gastrol2fcmat))
pdf(file = "Supplemental Figure S7L.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("SKM-GN: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Ppard","Ensembl"]
ppard_gastro_corpeaks <- rownames(gastropeakcortrim)[abs(gastropeakcortrim[,ourgene]) > 0.5]
ourpeak <- ppard_gastro_corpeaks
ourdf <- data.frame("RNA.L2FC" = gastrol2fcmat[ourgene,],
                    "ATAC.L2FC" = gastrosigatacl2fc[ourpeak,],
                    "Label" = colnames(gastrol2fcmat))
pdf(file = "Supplemental Figure S7M.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("SKM-GN: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Sik1","Ensembl"]
sik1_gastro_corpeaks <- rownames(gastropeakcortrim)[abs(gastropeakcortrim[,ourgene]) > 0.5]
ourpeak <- sik1_gastro_corpeaks
ourdf <- data.frame("RNA.L2FC" = gastrol2fcmat[ourgene,],
                    "ATAC.L2FC" = gastrosigatacl2fc[ourpeak,],
                    "Label" = colnames(gastrol2fcmat))
pdf(file = "Supplemental Figure S7N.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("SKM-GN: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

ourgene <- enstosym[enstosym$Symbol %in% "Igf1","Ensembl"]
igf1_liver_corpeaks <- rownames(liverpeakcortrim)[abs(liverpeakcortrim[,ourgene]) > 0.5]
ourpeak <- igf1_liver_corpeaks[2]
ourdf <- data.frame("RNA.L2FC" = liverl2fcmat[ourgene,],
                    "ATAC.L2FC" = liversigatacl2fc[ourpeak,],
                    "Label" = colnames(liverl2fcmat))
pdf(file = "Supplemental Figure S7O.pdf",width = 6,height = 6)
ggplot(data = ourdf,mapping = aes(x=RNA.L2FC,y=ATAC.L2FC))+geom_text(label = ourdf$Label,size = 5) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle(label = paste("LIVER: ",enstosym[ourgene,"Symbol"]," vs ",ourpeak,"\nPearson Correlation: ",substr(toString(cor(ourdf$RNA.L2FC,ourdf$ATAC.L2FC)),start = 1,stop = 7),sep = "")) + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),plot.title = element_text(size = 18,hjust = 0.5))
dev.off()

# Supplemental Figure S7P - 500x600
pdf(file = "Supplemental Figure S7P.pdf",width = 5,height = 6)
pheatmap(liverl2fcmat[intersect(tfanno[liver50peakmotifs[liver50peakmotifs$PositionID %in% igf1_liver_corpeaks,"Motif.Name"],"Ensembl"],rownames(liverl2fcmat)),],cluster_rows = F,cluster_cols = F,labels_row = enstosym[intersect(tfanno[liver50peakmotifs[liver50peakmotifs$PositionID %in% igf1_liver_corpeaks,"Motif.Name"],"Ensembl"],rownames(liverl2fcmat)),"Symbol"],breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = 0,display_numbers = T)
dev.off()


####
# Figure 4
#####

# Building proteomics dataset

gastropro <- prot_pr$timewise_dea[prot_pr$timewise_dea$tissue %in% "t55-gastrocnemius",]
heartpro <- prot_pr$timewise_dea[prot_pr$timewise_dea$tissue %in% "t58-heart",]
kidneypro <- prot_pr$timewise_dea[prot_pr$timewise_dea$tissue %in% "t59-kidney",]
liverpro <- prot_pr$timewise_dea[prot_pr$timewise_dea$tissue %in% "t68-liver",]
lungpro <- prot_pr$timewise_dea[prot_pr$timewise_dea$tissue %in% "t66-lung",]
whitepro <- prot_pr$timewise_dea[prot_pr$timewise_dea$tissue %in% "t70-white-adipose",]

gastroprol2fc <- matrix(0L,nrow = length(unique(gastropro$feature_ID)),ncol = 8)
rownames(gastroprol2fc) <- unique(gastropro$feature_ID)
colnames(gastroprol2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
for(i in 1:dim(gastroprol2fc)[1]){
  
  if(i%%100 == 0){
    print(paste("Gastro:",toString(i)))
  }
  
  ourpro <- rownames(gastroprol2fc)[i]
  ourpromat <- gastropro[gastropro$feature_ID %in% ourpro,]
  gastroprol2fc[i,"F W1"] <- ourpromat[ourpromat$sex %in% "female" & ourpromat$comparison_group %in% "1w","logFC"]
  gastroprol2fc[i,"F W2"] <- ourpromat[ourpromat$sex %in% "female" & ourpromat$comparison_group %in% "2w","logFC"]
  gastroprol2fc[i,"F W4"] <- ourpromat[ourpromat$sex %in% "female" & ourpromat$comparison_group %in% "4w","logFC"]
  gastroprol2fc[i,"F W8"] <- ourpromat[ourpromat$sex %in% "female" & ourpromat$comparison_group %in% "8w","logFC"]
  
  gastroprol2fc[i,"M W1"] <- ourpromat[ourpromat$sex %in% "male" & ourpromat$comparison_group %in% "1w","logFC"]
  gastroprol2fc[i,"M W2"] <- ourpromat[ourpromat$sex %in% "male" & ourpromat$comparison_group %in% "2w","logFC"]
  gastroprol2fc[i,"M W4"] <- ourpromat[ourpromat$sex %in% "male" & ourpromat$comparison_group %in% "4w","logFC"]
  gastroprol2fc[i,"M W8"] <- ourpromat[ourpromat$sex %in% "male" & ourpromat$comparison_group %in% "8w","logFC"]
}

heartprol2fc <- matrix(0L,nrow = length(unique(heartpro$feature_ID)),ncol = 8)
rownames(heartprol2fc) <- unique(heartpro$feature_ID)
colnames(heartprol2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
for(i in 1:dim(heartprol2fc)[1]){
  
  if(i%%100 == 0){
    print(paste("Heart:",toString(i)))
  }
  
  ourpro <- rownames(heartprol2fc)[i]
  ourpromat <- heartpro[heartpro$feature_ID %in% ourpro,]
  heartprol2fc[i,"F W1"] <- ourpromat[ourpromat$sex %in% "female" & ourpromat$comparison_group %in% "1w","logFC"]
  heartprol2fc[i,"F W2"] <- ourpromat[ourpromat$sex %in% "female" & ourpromat$comparison_group %in% "2w","logFC"]
  heartprol2fc[i,"F W4"] <- ourpromat[ourpromat$sex %in% "female" & ourpromat$comparison_group %in% "4w","logFC"]
  heartprol2fc[i,"F W8"] <- ourpromat[ourpromat$sex %in% "female" & ourpromat$comparison_group %in% "8w","logFC"]
  
  heartprol2fc[i,"M W1"] <- ourpromat[ourpromat$sex %in% "male" & ourpromat$comparison_group %in% "1w","logFC"]
  heartprol2fc[i,"M W2"] <- ourpromat[ourpromat$sex %in% "male" & ourpromat$comparison_group %in% "2w","logFC"]
  heartprol2fc[i,"M W4"] <- ourpromat[ourpromat$sex %in% "male" & ourpromat$comparison_group %in% "4w","logFC"]
  heartprol2fc[i,"M W8"] <- ourpromat[ourpromat$sex %in% "male" & ourpromat$comparison_group %in% "8w","logFC"]
}

kidneyprol2fc <- matrix(0L,nrow = length(unique(kidneypro$feature_ID)),ncol = 8)
rownames(kidneyprol2fc) <- unique(kidneypro$feature_ID)
colnames(kidneyprol2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
for(i in 1:dim(kidneyprol2fc)[1]){
  
  if(i%%100 == 0){
    print(paste("Kidney:",toString(i)))
  }
  
  ourpro <- rownames(kidneyprol2fc)[i]
  ourpromat <- kidneypro[kidneypro$feature_ID %in% ourpro,]
  kidneyprol2fc[i,"F W1"] <- ourpromat[ourpromat$sex %in% "female" & ourpromat$comparison_group %in% "1w","logFC"]
  kidneyprol2fc[i,"F W2"] <- ourpromat[ourpromat$sex %in% "female" & ourpromat$comparison_group %in% "2w","logFC"]
  kidneyprol2fc[i,"F W4"] <- ourpromat[ourpromat$sex %in% "female" & ourpromat$comparison_group %in% "4w","logFC"]
  kidneyprol2fc[i,"F W8"] <- ourpromat[ourpromat$sex %in% "female" & ourpromat$comparison_group %in% "8w","logFC"]
  
  kidneyprol2fc[i,"M W1"] <- ourpromat[ourpromat$sex %in% "male" & ourpromat$comparison_group %in% "1w","logFC"]
  kidneyprol2fc[i,"M W2"] <- ourpromat[ourpromat$sex %in% "male" & ourpromat$comparison_group %in% "2w","logFC"]
  kidneyprol2fc[i,"M W4"] <- ourpromat[ourpromat$sex %in% "male" & ourpromat$comparison_group %in% "4w","logFC"]
  kidneyprol2fc[i,"M W8"] <- ourpromat[ourpromat$sex %in% "male" & ourpromat$comparison_group %in% "8w","logFC"]
}

liverprol2fc <- matrix(0L,nrow = length(unique(liverpro$feature_ID)),ncol = 8)
rownames(liverprol2fc) <- unique(liverpro$feature_ID)
colnames(liverprol2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
for(i in 1:dim(liverprol2fc)[1]){
  
  if(i%%100 == 0){
    print(paste("Liver:",toString(i)))
  }
  
  ourpro <- rownames(liverprol2fc)[i]
  ourpromat <- liverpro[liverpro$feature_ID %in% ourpro,]
  liverprol2fc[i,"F W1"] <- ourpromat[ourpromat$sex %in% "female" & ourpromat$comparison_group %in% "1w","logFC"]
  liverprol2fc[i,"F W2"] <- ourpromat[ourpromat$sex %in% "female" & ourpromat$comparison_group %in% "2w","logFC"]
  liverprol2fc[i,"F W4"] <- ourpromat[ourpromat$sex %in% "female" & ourpromat$comparison_group %in% "4w","logFC"]
  liverprol2fc[i,"F W8"] <- ourpromat[ourpromat$sex %in% "female" & ourpromat$comparison_group %in% "8w","logFC"]
  
  liverprol2fc[i,"M W1"] <- ourpromat[ourpromat$sex %in% "male" & ourpromat$comparison_group %in% "1w","logFC"]
  liverprol2fc[i,"M W2"] <- ourpromat[ourpromat$sex %in% "male" & ourpromat$comparison_group %in% "2w","logFC"]
  liverprol2fc[i,"M W4"] <- ourpromat[ourpromat$sex %in% "male" & ourpromat$comparison_group %in% "4w","logFC"]
  liverprol2fc[i,"M W8"] <- ourpromat[ourpromat$sex %in% "male" & ourpromat$comparison_group %in% "8w","logFC"]
}

lungprol2fc <- matrix(0L,nrow = length(unique(lungpro$feature_ID)),ncol = 8)
rownames(lungprol2fc) <- unique(lungpro$feature_ID)
colnames(lungprol2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
for(i in 1:dim(lungprol2fc)[1]){
  
  if(i%%100 == 0){
    print(paste("Lung:",toString(i)))
  }
  
  ourpro <- rownames(lungprol2fc)[i]
  ourpromat <- lungpro[lungpro$feature_ID %in% ourpro,]
  lungprol2fc[i,"F W1"] <- ourpromat[ourpromat$sex %in% "female" & ourpromat$comparison_group %in% "1w","logFC"]
  lungprol2fc[i,"F W2"] <- ourpromat[ourpromat$sex %in% "female" & ourpromat$comparison_group %in% "2w","logFC"]
  lungprol2fc[i,"F W4"] <- ourpromat[ourpromat$sex %in% "female" & ourpromat$comparison_group %in% "4w","logFC"]
  lungprol2fc[i,"F W8"] <- ourpromat[ourpromat$sex %in% "female" & ourpromat$comparison_group %in% "8w","logFC"]
  
  lungprol2fc[i,"M W1"] <- ourpromat[ourpromat$sex %in% "male" & ourpromat$comparison_group %in% "1w","logFC"]
  lungprol2fc[i,"M W2"] <- ourpromat[ourpromat$sex %in% "male" & ourpromat$comparison_group %in% "2w","logFC"]
  lungprol2fc[i,"M W4"] <- ourpromat[ourpromat$sex %in% "male" & ourpromat$comparison_group %in% "4w","logFC"]
  lungprol2fc[i,"M W8"] <- ourpromat[ourpromat$sex %in% "male" & ourpromat$comparison_group %in% "8w","logFC"]
}

whiteprol2fc <- matrix(0L,nrow = length(unique(whitepro$feature_ID)),ncol = 8)
rownames(whiteprol2fc) <- unique(whitepro$feature_ID)
colnames(whiteprol2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
for(i in 1:dim(whiteprol2fc)[1]){
  
  if(i%%100 == 0){
    print(paste("WhiteAd:",toString(i)))
  }
  
  ourpro <- rownames(whiteprol2fc)[i]
  ourpromat <- whitepro[whitepro$feature_ID %in% ourpro,]
  whiteprol2fc[i,"F W1"] <- ourpromat[ourpromat$sex %in% "female" & ourpromat$comparison_group %in% "1w","logFC"]
  whiteprol2fc[i,"F W2"] <- ourpromat[ourpromat$sex %in% "female" & ourpromat$comparison_group %in% "2w","logFC"]
  whiteprol2fc[i,"F W4"] <- ourpromat[ourpromat$sex %in% "female" & ourpromat$comparison_group %in% "4w","logFC"]
  whiteprol2fc[i,"F W8"] <- ourpromat[ourpromat$sex %in% "female" & ourpromat$comparison_group %in% "8w","logFC"]
  
  whiteprol2fc[i,"M W1"] <- ourpromat[ourpromat$sex %in% "male" & ourpromat$comparison_group %in% "1w","logFC"]
  whiteprol2fc[i,"M W2"] <- ourpromat[ourpromat$sex %in% "male" & ourpromat$comparison_group %in% "2w","logFC"]
  whiteprol2fc[i,"M W4"] <- ourpromat[ourpromat$sex %in% "male" & ourpromat$comparison_group %in% "4w","logFC"]
  whiteprol2fc[i,"M W8"] <- ourpromat[ourpromat$sex %in% "male" & ourpromat$comparison_group %in% "8w","logFC"]
}

gastrotfproids <- intersect(unique(tfproanno$Gastro.Pro.ID),rownames(gastroprol2fc))
hearttfproids <- intersect(unique(tfproanno$Heart.Pro.ID),rownames(heartprol2fc))
kidneytfproids <- intersect(unique(tfproanno$Kidney.Pro.ID),rownames(kidneyprol2fc))
livertfproids <- intersect(unique(tfproanno$Liver.Pro.ID),rownames(liverprol2fc))
lungtfproids <- intersect(unique(tfproanno$Lung.Pro.ID),rownames(lungprol2fc))
whitetfproids <- intersect(unique(tfproanno$WhiteAd.Pro.ID),rownames(whiteprol2fc))


tfsigrna <- list("SKM-GN" = transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue %in% "t55-gastrocnemius" & 
                                                              transcript_rna_seq$training_dea$feature_ID %in% tfanno$Ensemb &
                                                              transcript_rna_seq$training_dea$adj_p_value < 0.05,"feature_ID"],
                 "HEART" = transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue %in% "t58-heart" & 
                                                             transcript_rna_seq$training_dea$feature_ID %in% tfanno$Ensemb &
                                                             transcript_rna_seq$training_dea$adj_p_value < 0.05,"feature_ID"],
                 "HIPPOC" = transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue %in% "t52-hippocampus" & 
                                                              transcript_rna_seq$training_dea$feature_ID %in% tfanno$Ensemb &
                                                              transcript_rna_seq$training_dea$adj_p_value < 0.05,"feature_ID"],
                 "KIDNEY" = transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue %in% "t59-kidney" & 
                                                              transcript_rna_seq$training_dea$feature_ID %in% tfanno$Ensemb &
                                                              transcript_rna_seq$training_dea$adj_p_value < 0.05,"feature_ID"],
                 "LIVER" = transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue %in% "t68-liver" & 
                                                             transcript_rna_seq$training_dea$feature_ID %in% tfanno$Ensemb &
                                                             transcript_rna_seq$training_dea$adj_p_value < 0.05,"feature_ID"],
                 "LUNG" = transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue %in% "t66-lung" & 
                                                            transcript_rna_seq$training_dea$feature_ID %in% tfanno$Ensemb &
                                                            transcript_rna_seq$training_dea$adj_p_value < 0.05,"feature_ID"],
                 "BAT" = transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue %in% "t69-brown-adipose" & 
                                                           transcript_rna_seq$training_dea$feature_ID %in% tfanno$Ensemb &
                                                           transcript_rna_seq$training_dea$adj_p_value < 0.05,"feature_ID"],
                 "WAT-SC" = transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue %in% "t70-white-adipose" & 
                                                              transcript_rna_seq$training_dea$feature_ID %in% tfanno$Ensemb &
                                                              transcript_rna_seq$training_dea$adj_p_value < 0.05,"feature_ID"])



tfsigpro <- list("SKM-GN" = prot_pr$training_dea[prot_pr$training_dea$tissue %in% "t55-gastrocnemius" & 
                                                   prot_pr$training_dea$feature_ID %in% tfproanno$Gastro.Pro.ID &
                                                   prot_pr$training_dea$p_value < 0.06,"feature_ID"],
                 "HEART" = prot_pr$training_dea[prot_pr$training_dea$tissue %in% "t58-heart" & 
                                                  prot_pr$training_dea$feature_ID %in% tfproanno$Heart.Pro.ID &
                                                  prot_pr$training_dea$p_value < 0.05,"feature_ID"],
                 "KIDNEY" = prot_pr$training_dea[prot_pr$training_dea$tissue %in% "t59-kidney" & 
                                                   prot_pr$training_dea$feature_ID %in% tfproanno$Kidney.Pro.ID &
                                                   prot_pr$training_dea$p_value < 0.05,"feature_ID"],
                 "LIVER" = prot_pr$training_dea[prot_pr$training_dea$tissue %in% "t68-liver" & 
                                                  prot_pr$training_dea$feature_ID %in% tfproanno$Liver.Pro.ID &
                                                  prot_pr$training_dea$p_value < 0.05,"feature_ID"],
                 "LUNG" = prot_pr$training_dea[prot_pr$training_dea$tissue %in% "t66-lung" & 
                                                 prot_pr$training_dea$feature_ID %in% tfproanno$Lung.Pro.ID &
                                                 prot_pr$training_dea$p_value < 0.05,"feature_ID"],
                 "WAT-SC" = prot_pr$training_dea[prot_pr$training_dea$tissue %in% "t70-white-adipose" & 
                                                   prot_pr$training_dea$feature_ID %in% tfproanno$WhiteAd.Pro.ID &
                                                   prot_pr$training_dea$p_value < 0.05,"feature_ID"])

tfproanno$tftrim <- gsub("\\(.*","",rownames(tfproanno))

####
# Figure 4A
#####

tfrnapval <- matrix(0L,nrow = length(Reduce(union,tfsigrna)),ncol = 8)
rownames(tfrnapval) <- Reduce(union,tfsigrna)
colnames(tfrnapval) <- c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC")
for(i in 1:dim(tfrnapval)[1]){
  ourrna <- rownames(tfrnapval)[i]
  ourrnatraining <- transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$feature_ID %in% ourrna,]
  if("SKM-GN" %in% ourrnatraining$tissue_abbreviation){
    tfrnapval[i,"SKM-GN"] <- -log10(ourrnatraining[ourrnatraining$tissue %in% "t55-gastrocnemius","p_value"])  
  }
  if("HEART" %in% ourrnatraining$tissue_abbreviation){
    tfrnapval[i,"HEART"] <- -log10(ourrnatraining[ourrnatraining$tissue %in% "t58-heart","p_value"])  
  }
  if("HIPPOC" %in% ourrnatraining$tissue_abbreviation){
    tfrnapval[i,"HIPPOC"] <- -log10(ourrnatraining[ourrnatraining$tissue %in% "t52-hippocampus","p_value"])  
  }
  if("KIDNEY" %in% ourrnatraining$tissue_abbreviation){
    tfrnapval[i,"KIDNEY"] <- -log10(ourrnatraining[ourrnatraining$tissue %in% "t59-kidney","p_value"])  
  }
  if("LIVER" %in% ourrnatraining$tissue_abbreviation){
    tfrnapval[i,"LIVER"] <- -log10(ourrnatraining[ourrnatraining$tissue %in% "t68-liver","p_value"])  
  }
  if("LUNG" %in% ourrnatraining$tissue_abbreviation){
    tfrnapval[i,"LUNG"] <- -log10(ourrnatraining[ourrnatraining$tissue %in% "t66-lung","p_value"])  
  }
  if("BAT" %in% ourrnatraining$tissue_abbreviation){
    tfrnapval[i,"BAT"] <- -log10(ourrnatraining[ourrnatraining$tissue %in% "t69-brown-adipose","p_value"])  
  }
  if("WAT-SC" %in% ourrnatraining$tissue_abbreviation){
    tfrnapval[i,"WAT-SC"] <- -log10(ourrnatraining[ourrnatraining$tissue %in% "t70-white-adipose","p_value"])  
  }
  
}

tfrnamaxl2fc <- matrix(0L,nrow = length(Reduce(union,tfsigrna)),ncol = 8)
rownames(tfrnamaxl2fc) <- Reduce(union,tfsigrna)
colnames(tfrnamaxl2fc) <- c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC")
for(i in 1:dim(tfrnamaxl2fc)[1]){
  ourrna <- rownames(tfrnamaxl2fc)[i]
  ourrnatimewise <- transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$feature_ID %in% ourrna,]
  if("SKM-GN" %in% ourrnatimewise$tissue_abbreviation){
    tfrnamaxl2fc[i,"SKM-GN"] <- ourrnatimewise[ourrnatimewise$tissue %in% "t55-gastrocnemius","logFC"][order(-abs(ourrnatimewise[ourrnatimewise$tissue %in% "t55-gastrocnemius","logFC"]))][1]  
  }
  if("HEART" %in% ourrnatimewise$tissue_abbreviation){
    tfrnamaxl2fc[i,"HEART"] <- ourrnatimewise[ourrnatimewise$tissue %in% "t58-heart","logFC"][order(-abs(ourrnatimewise[ourrnatimewise$tissue %in% "t58-heart","logFC"]))][1] 
  }
  if("HIPPOC" %in% ourrnatimewise$tissue_abbreviation){
    tfrnamaxl2fc[i,"HIPPOC"] <- ourrnatimewise[ourrnatimewise$tissue %in% "t52-hippocampus","logFC"][order(-abs(ourrnatimewise[ourrnatimewise$tissue %in% "t52-hippocampus","logFC"]))][1] 
  }
  if("KIDNEY" %in% ourrnatimewise$tissue_abbreviation){
    tfrnamaxl2fc[i,"KIDNEY"] <- ourrnatimewise[ourrnatimewise$tissue %in% "t59-kidney","logFC"][order(-abs(ourrnatimewise[ourrnatimewise$tissue %in% "t59-kidney","logFC"]))][1]
  }
  if("LIVER" %in% ourrnatimewise$tissue_abbreviation){
    tfrnamaxl2fc[i,"LIVER"] <- ourrnatimewise[ourrnatimewise$tissue %in% "t68-liver","logFC"][order(-abs(ourrnatimewise[ourrnatimewise$tissue %in% "t68-liver","logFC"]))][1] 
  }
  if("LUNG" %in% ourrnatimewise$tissue_abbreviation){
    tfrnamaxl2fc[i,"LUNG"] <- ourrnatimewise[ourrnatimewise$tissue %in% "t66-lung","logFC"][order(-abs(ourrnatimewise[ourrnatimewise$tissue %in% "t66-lung","logFC"]))][1]
  }
  if("BAT" %in% ourrnatimewise$tissue_abbreviation){
    tfrnamaxl2fc[i,"BAT"] <- ourrnatimewise[ourrnatimewise$tissue %in% "t69-brown-adipose","logFC"][order(-abs(ourrnatimewise[ourrnatimewise$tissue %in% "t69-brown-adipose","logFC"]))][1] 
  }
  if("WAT-SC" %in% ourrnatimewise$tissue_abbreviation){
    tfrnamaxl2fc[i,"WAT-SC"] <- ourrnatimewise[ourrnatimewise$tissue %in% "t70-white-adipose","logFC"][order(-abs(ourrnatimewise[ourrnatimewise$tissue %in% "t70-white-adipose","logFC"]))][1]  
  }
  
}

tfrnacluster <- hclust(d = dist(tfrnamaxl2fc))
tfrnaregcluster <- hclust(d = dist(t(tfrnamaxl2fc)))

tfrnapandl2fc <- data.frame("p-val" = as.vector(tfrnapval[tfrnacluster$order,tfrnaregcluster$order]),
                            "L2FC" = as.vector(tfrnamaxl2fc[tfrnacluster$order,tfrnaregcluster$order]),
                            "Tissue" = factor(rep(colnames(tfrnapval[,tfrnaregcluster$order]),each = dim(tfrnapval)[1]),levels = c(colnames(tfrnapval[,tfrnaregcluster$order]))),
                            "TF" = factor(rep(enstosym[rownames(tfrnamaxl2fc[tfrnacluster$order,tfrnaregcluster$order]),"Symbol"],dim(tfrnapval)[2]),levels = enstosym[rownames(tfrnamaxl2fc[tfrnacluster$order,tfrnaregcluster$order]),"Symbol"]))
rownames(tfrnapandl2fc) <- paste(tfrnapandl2fc$TF,"_",tfrnapandl2fc$Tissue,sep = "")

tfrnapandl2fcthresh <- tfrnapandl2fc[tfrnapandl2fc$TF %in% rownames(table(tfrnapandl2fc[tfrnapandl2fc$p.val > 3 & abs(tfrnapandl2fc$L2FC) > 0.25,"TF"]))[table(tfrnapandl2fc[tfrnapandl2fc$p.val > 3 & abs(tfrnapandl2fc$L2FC) > 0.25,"TF"])>0],]

tfrnapandl2fcthresh[tfrnapandl2fcthresh$L2FC > 1,"L2FC"] <- 1
tfrnapandl2fcthresh[tfrnapandl2fcthresh$L2FC < (-1),"L2FC"] <- (-1)

tfrnapandl2fc[tfrnapandl2fc$L2FC > 1,"L2FC"] <- 1
tfrnapandl2fc[tfrnapandl2fc$L2FC < (-1),"L2FC"] <- (-1)

pdf(file = "Figure 4A.pdf",width = 8,height = 12)
ggplot(tfrnapandl2fcthresh,aes(x = Tissue,y = TF, size = p.val,fill = L2FC)) + geom_point(shape = 21,stroke = 0) + theme_classic() +scale_fill_gradientn(colors = c("blue","white","red"),limits = c(-1,1)) + labs(size = "-log10 p-val") + theme(axis.text.x = element_text(size = 15,angle = 45,vjust = 1,hjust = 1),axis.text.y = element_text(size = 15),axis.title = element_blank(),legend.title = element_text(size = 15),legend.text = element_text(size = 10))
dev.off()

####
# Figure 4B
#####

tfproanno$tftrim <- gsub("\\(.*","",rownames(tfproanno))

tfpropval <- matrix(0L,nrow = length(Reduce(union,tfsigpro)),ncol = 6)
rownames(tfpropval) <- Reduce(union,tfsigpro)
colnames(tfpropval) <- c("SKM-GN","HEART","KIDNEY","LIVER","LUNG","WAT-SC")
for(i in 1:dim(tfpropval)[1]){
  ourpro <- rownames(tfpropval)[i]
  ourprotraining <- prot_pr$training_dea[prot_pr$training_dea$feature_ID %in% ourpro,]
  if("SKM-GN" %in% ourprotraining$tissue_abbreviation){
    tfpropval[i,"SKM-GN"] <- -log10(ourprotraining[ourprotraining$tissue %in% "t55-gastrocnemius" & 
                                                     ourprotraining$feature_ID %in% ourpro,"p_value"])  
  }
  if("HEART" %in% ourprotraining$tissue_abbreviation){
    tfpropval[i,"HEART"] <- -log10(ourprotraining[ourprotraining$tissue %in% "t58-heart" & 
                                                    ourprotraining$feature_ID %in% ourpro,"p_value"])  
  }
  if("KIDNEY" %in% ourprotraining$tissue_abbreviation){
    tfpropval[i,"KIDNEY"] <- -log10(ourprotraining[ourprotraining$tissue %in% "t59-kidney" & 
                                                     ourprotraining$feature_ID %in% ourpro,"p_value"])  
  }
  if("LIVER" %in% ourprotraining$tissue_abbreviation){
    tfpropval[i,"LIVER"] <- -log10(ourprotraining[ourprotraining$tissue %in% "t68-liver" & 
                                                    ourprotraining$feature_ID %in% ourpro,"p_value"])  
  }
  if("LUNG" %in% ourprotraining$tissue_abbreviation){
    tfpropval[i,"LUNG"] <- -log10(ourprotraining[ourprotraining$tissue %in% "t66-lung" & 
                                                   ourprotraining$feature_ID %in% ourpro,"p_value"])  
  }
  if("WAT-SC" %in% ourprotraining$tissue_abbreviation){
    tfpropval[i,"WAT-SC"] <- -log10(ourprotraining[ourprotraining$tissue %in% "t70-white-adipose" & 
                                                     ourprotraining$feature_ID %in% ourpro,"p_value"])  
  }
  
}

tfproppvalmeta <- data.frame(row.names = colnames(tfpropval),
                             "Tissue" = colnames(tfpropval))
tfrnappvalmeta <- data.frame(row.names = colnames(tfrnapval),
                             "Tissue" = colnames(tfrnapval))
Reduce(union,list(tfproanno$Gastro.Pro.ID,tfproanno$Heart.Pro.ID,tfproanno$Kidney.Pro.ID,tfproanno$Liver.Pro.ID,tfproanno$Lung.Pro.ID,tfproanno$WhiteAd.Pro.ID))


tfproidanno <- data.frame(row.names = Reduce(union,list(tfproanno$Gastro.Pro.ID,tfproanno$Heart.Pro.ID,tfproanno$Kidney.Pro.ID,tfproanno$Liver.Pro.ID,tfproanno$Lung.Pro.ID,tfproanno$WhiteAd.Pro.ID)),
                          "Protein_ID" = Reduce(union,list(tfproanno$Gastro.Pro.ID,tfproanno$Heart.Pro.ID,tfproanno$Kidney.Pro.ID,tfproanno$Liver.Pro.ID,tfproanno$Lung.Pro.ID,tfproanno$WhiteAd.Pro.ID)),
                          "Ensembl" = "",
                          "Gene.Name" = "")
tfproidanno <- tfproidanno[c(2:dim(tfproidanno)[1]),]
for(i in 1:dim(tfproidanno)[1]){
  ourpro <- rownames(tfproidanno)[i]
  if(ourpro %in% tfproanno$Gastro.Pro.ID){
    tfproidanno[i,"Ensembl"] <- tfproanno[tfproanno$Gastro.Pro.ID %in% ourpro,"Ensembl"][1]
    tfproidanno[i,"Gene.Name"] <- tfproanno[tfproanno$Gastro.Pro.ID %in% ourpro,"Gene.Name"][1]
  } else if(ourpro %in% tfproanno$Heart.Pro.ID){
    tfproidanno[i,"Ensembl"] <- tfproanno[tfproanno$Heart.Pro.ID %in% ourpro,"Ensembl"][1]
    tfproidanno[i,"Gene.Name"] <- tfproanno[tfproanno$Heart.Pro.ID %in% ourpro,"Gene.Name"][1]
  } else if(ourpro %in% tfproanno$Kidney.Pro.ID){
    tfproidanno[i,"Ensembl"] <- tfproanno[tfproanno$Kidney.Pro.ID %in% ourpro,"Ensembl"][1]
    tfproidanno[i,"Gene.Name"] <- tfproanno[tfproanno$Kidney.Pro.ID %in% ourpro,"Gene.Name"][1]
  } else if(ourpro %in% tfproanno$Liver.Pro.ID){
    tfproidanno[i,"Ensembl"] <- tfproanno[tfproanno$Liver.Pro.ID %in% ourpro,"Ensembl"][1]
    tfproidanno[i,"Gene.Name"] <- tfproanno[tfproanno$Liver.Pro.ID %in% ourpro,"Gene.Name"][1]
  } else if(ourpro %in% tfproanno$Lung.Pro.ID){
    tfproidanno[i,"Ensembl"] <- tfproanno[tfproanno$Lung.Pro.ID %in% ourpro,"Ensembl"][1]
    tfproidanno[i,"Gene.Name"] <- tfproanno[tfproanno$Lung.Pro.ID %in% ourpro,"Gene.Name"][1]
  } else if(ourpro %in% tfproanno$WhiteAd.Pro.ID){
    tfproidanno[i,"Ensembl"] <- tfproanno[tfproanno$WhiteAd.Pro.ID %in% ourpro,"Ensembl"][1]
    tfproidanno[i,"Gene.Name"] <- tfproanno[tfproanno$WhiteAd.Pro.ID %in% ourpro,"Gene.Name"][1]
  }
}

tfpropval[11,] <- tfpropval[11,] + tfpropval[64,]
tfpropval[29,] <- colMeans(rbind(tfpropval[29,],tfpropval[62,]))
tfpropval[22,] <- tfpropval[22,] + tfpropval[48,]
tfpropval[27,] <- tfpropval[27,] + tfpropval[55,]
tfpropval[4,] <- tfpropval[4,] + tfpropval[63,]

tfpropval <- tfpropval[c(1:47,49:54,56:61),]

tfpromaxl2fc <- matrix(0L,nrow = length(Reduce(union,tfsigpro)),ncol = 6)
rownames(tfpromaxl2fc) <- Reduce(union,tfsigpro)
colnames(tfpromaxl2fc) <- c("SKM-GN","HEART","KIDNEY","LIVER","LUNG","WAT-SC")
for(i in 1:dim(tfpromaxl2fc)[1]){
  ourpro <- rownames(tfpromaxl2fc)[i]
  ourprotimewise <- prot_pr$timewise_dea[prot_pr$timewise_dea$feature_ID %in% ourpro,]
  if("SKM-GN" %in% ourprotimewise$tissue_abbreviation){
    tfpromaxl2fc[i,"SKM-GN"] <- ourprotimewise[ourprotimewise$tissue %in% "t55-gastrocnemius" & 
                                                 ourprotimewise$feature_ID %in% ourpro,"logFC"][order(-abs(ourprotimewise[ourprotimewise$tissue %in% "t55-gastrocnemius" & 
                                                                                                                            ourprotimewise$feature_ID %in% ourpro,"logFC"]))][1]  
  }
  if("HEART" %in% ourprotimewise$tissue_abbreviation){
    tfpromaxl2fc[i,"HEART"] <- ourprotimewise[ourprotimewise$tissue %in% "t58-heart" & 
                                                ourprotimewise$feature_ID %in% ourpro,"logFC"][order(-abs(ourprotimewise[ourprotimewise$tissue %in% "t58-heart" & 
                                                                                                                           ourprotimewise$feature_ID %in% ourpro,"logFC"]))][1] 
  }
  if("KIDNEY" %in% ourprotimewise$tissue_abbreviation){
    tfpromaxl2fc[i,"KIDNEY"] <- ourprotimewise[ourprotimewise$tissue %in% "t59-kidney" & 
                                                 ourprotimewise$feature_ID %in% ourpro,"logFC"][order(-abs(ourprotimewise[ourprotimewise$tissue %in% "t59-kidney" & 
                                                                                                                            ourprotimewise$feature_ID %in% ourpro,"logFC"]))][1]
  }
  if("LIVER" %in% ourprotimewise$tissue_abbreviation){
    tfpromaxl2fc[i,"LIVER"] <- ourprotimewise[ourprotimewise$tissue %in% "t68-liver" & 
                                                ourprotimewise$feature_ID %in% ourpro,"logFC"][order(-abs(ourprotimewise[ourprotimewise$tissue %in% "t68-liver" & 
                                                                                                                           ourprotimewise$feature_ID %in% ourpro,"logFC"]))][1] 
  }
  if("LUNG" %in% ourprotimewise$tissue_abbreviation){
    tfpromaxl2fc[i,"LUNG"] <- ourprotimewise[ourprotimewise$tissue %in% "t66-lung" & 
                                               ourprotimewise$feature_ID %in% ourpro,"logFC"][order(-abs(ourprotimewise[ourprotimewise$tissue %in% "t66-lung" & 
                                                                                                                          ourprotimewise$feature_ID %in% ourpro,"logFC"]))][1]
  }
  if("WAT-SC" %in% ourprotimewise$tissue_abbreviation){
    tfpromaxl2fc[i,"WAT-SC"] <- ourprotimewise[ourprotimewise$tissue %in% "t70-white-adipose" & 
                                                 ourprotimewise$feature_ID %in% ourpro,"logFC"][order(-abs(ourprotimewise[ourprotimewise$tissue %in% "t70-white-adipose" & 
                                                                                                                            ourprotimewise$feature_ID %in% ourpro,"logFC"]))][1]  
  }
  
}

tfpromaxl2fc[11,] <- tfpromaxl2fc[11,] + tfpromaxl2fc[64,]
tfpromaxl2fc[29,] <- colMeans(rbind(tfpromaxl2fc[29,],tfpromaxl2fc[62,]))
tfpromaxl2fc[22,] <- tfpromaxl2fc[22,] + tfpromaxl2fc[48,]
tfpromaxl2fc[27,] <- tfpromaxl2fc[27,] + tfpromaxl2fc[55,]
tfpromaxl2fc[4,] <- tfpromaxl2fc[4,] + tfpromaxl2fc[63,]

tfpromaxl2fc <- tfpromaxl2fc[c(1:47,49:54,56:61),]

tfprocluster <- hclust(d = dist(tfpromaxl2fc))
tfproregcluster <- hclust(d = dist(t(tfpromaxl2fc)))

tfpropandl2fc <- data.frame("p-val" = as.vector(tfpropval[tfprocluster$order,tfproregcluster$order]),
                            "L2FC" = as.vector(tfpromaxl2fc[tfprocluster$order,tfproregcluster$order]),
                            "Tissue" = factor(rep(colnames(tfpropval[,tfproregcluster$order]),each = dim(tfpropval)[1]),levels = c(colnames(tfpropval[,tfproregcluster$order]))),
                            "TF" = factor(rep(tfproidanno[rownames(tfpromaxl2fc[tfprocluster$order,]),"Gene.Name"],dim(tfpropval)[2]),levels = tfproidanno[rownames(tfpromaxl2fc[tfprocluster$order,]),"Gene.Name"]))
rownames(tfpropandl2fc) <- paste(tfpropandl2fc$TF,"_",tfpropandl2fc$Tissue,sep = "")

pdf(file = "Figure 4B.pdf",width = 7,height = 10)
ggplot(tfpropandl2fc,aes(x = Tissue,y = TF, size = p.val,fill = L2FC)) + geom_point(shape = 21,stroke = 0) + theme_classic() +scale_fill_gradientn(colors = c("blue","white","red"),limits = c(-1,1))
dev.off()

####
# Figure 4C-D
#####

# Correlation plot of TF expression/protein response

tfgeneresponsecor <- matrix(0L,nrow = 8,ncol = 8)
rownames(tfgeneresponsecor) <- colnames(tfrnamaxl2fc)
colnames(tfgeneresponsecor) <- colnames(tfrnamaxl2fc)

tfproresponsecor <- matrix(0L,nrow = 6,ncol = 6)
rownames(tfproresponsecor) <- colnames(tfpromaxl2fc)
colnames(tfproresponsecor) <- colnames(tfpromaxl2fc)

for(i in 1:8){
  for(j in 1:8){
    tis1 <- rownames(tfgeneresponsecor)[i]
    tis2 <- rownames(tfgeneresponsecor)[j]
    ourgenes <- intersect(rownames(tfrnamaxl2fc[tfrnamaxl2fc[,tis1] != 0,]),rownames(tfrnamaxl2fc[tfrnamaxl2fc[,tis2] != 0,]))
    tfgeneresponsecor[i,j] <- cor(tfrnamaxl2fc[ourgenes,tis1],tfrnamaxl2fc[ourgenes,tis2])
  }
}


for(i in 1:6){
  for(j in 1:6){
    tis1 <- rownames(tfproresponsecor)[i]
    tis2 <- rownames(tfproresponsecor)[j]
    ourpros <- intersect(rownames(tfpromaxl2fc[tfpromaxl2fc[,tis1] != 0,]),rownames(tfpromaxl2fc[tfpromaxl2fc[,tis2] != 0,]))
    tfproresponsecor[i,j] <- cor(tfpromaxl2fc[ourpros,tis1],tfpromaxl2fc[ourpros,tis2])
  }
}

# Figure 4C
pdf(file = "Figure 4C.pdf",width = 7, height = 4)
pheatmap(tfgeneresponsecor,angle_col = 0,breaks = seq(-0.5,0.5,length.out = 101),color = colorpanel(101,"blue","white","red"),display_numbers = T,number_color = "black")
dev.off()

# Figure 4D
pdf(file = "Figure 4D.pdf",width = 7, height = 4)
pheatmap(tfproresponsecor,angle_col = 0,breaks = seq(-0.5,0.5,length.out = 101),color = colorpanel(101,"blue","white","red"),display_numbers = T,number_color = "black")
dev.off()

####
# Figure 4E-G and Supplemental Figure S8
#####

tf50sigpro <- list("SKM-GN" = prot_pr$training_dea[prot_pr$training_dea$tissue %in% "t55-gastrocnemius" & 
                                                     prot_pr$training_dea$feature_ID %in% tfproanno$Gastro.Pro.ID &
                                                     prot_pr$training_dea$p_value < 0.1,"feature_ID"],
                   "HEART" = prot_pr$training_dea[prot_pr$training_dea$tissue %in% "t58-heart" & 
                                                    prot_pr$training_dea$feature_ID %in% tfproanno$Heart.Pro.ID &
                                                    prot_pr$training_dea$p_value < 0.1,"feature_ID"],
                   "KIDNEY" = prot_pr$training_dea[prot_pr$training_dea$tissue %in% "t59-kidney" & 
                                                     prot_pr$training_dea$feature_ID %in% tfproanno$Kidney.Pro.ID &
                                                     prot_pr$training_dea$p_value < 0.1,"feature_ID"],
                   "LIVER" = prot_pr$training_dea[prot_pr$training_dea$tissue %in% "t68-liver" & 
                                                    prot_pr$training_dea$feature_ID %in% tfproanno$Liver.Pro.ID &
                                                    prot_pr$training_dea$p_value < 0.1,"feature_ID"],
                   "LUNG" = prot_pr$training_dea[prot_pr$training_dea$tissue %in% "t66-lung" & 
                                                   prot_pr$training_dea$feature_ID %in% tfproanno$Lung.Pro.ID &
                                                   prot_pr$training_dea$p_value < 0.1,"feature_ID"],
                   "WAT-SC" = prot_pr$training_dea[prot_pr$training_dea$tissue %in% "t70-white-adipose" & 
                                                     prot_pr$training_dea$feature_ID %in% tfproanno$WhiteAd.Pro.ID &
                                                     prot_pr$training_dea$p_value < 0.1,"feature_ID"])

tf50propval <- matrix(0L,nrow = length(Reduce(union,tf50sigpro)),ncol = 6)
rownames(tf50propval) <- Reduce(union,tf50sigpro)
colnames(tf50propval) <- c("SKM-GN","HEART","KIDNEY","LIVER","LUNG","WAT-SC")
for(i in 1:dim(tf50propval)[1]){
  ourpro <- rownames(tf50propval)[i]
  ourprotraining <- prot_pr$training_dea[prot_pr$training_dea$feature_ID %in% ourpro,]
  if("SKM-GN" %in% ourprotraining$tissue_abbreviation){
    tf50propval[i,"SKM-GN"] <- -log10(ourprotraining[ourprotraining$tissue %in% "t55-gastrocnemius" & 
                                                       ourprotraining$feature_ID %in% ourpro,"p_value"])  
  }
  if("HEART" %in% ourprotraining$tissue_abbreviation){
    tf50propval[i,"HEART"] <- -log10(ourprotraining[ourprotraining$tissue %in% "t58-heart" & 
                                                      ourprotraining$feature_ID %in% ourpro,"p_value"])  
  }
  if("KIDNEY" %in% ourprotraining$tissue_abbreviation){
    tf50propval[i,"KIDNEY"] <- -log10(ourprotraining[ourprotraining$tissue %in% "t59-kidney" & 
                                                       ourprotraining$feature_ID %in% ourpro,"p_value"])  
  }
  if("LIVER" %in% ourprotraining$tissue_abbreviation){
    tf50propval[i,"LIVER"] <- -log10(ourprotraining[ourprotraining$tissue %in% "t68-liver" & 
                                                      ourprotraining$feature_ID %in% ourpro,"p_value"])  
  }
  if("LUNG" %in% ourprotraining$tissue_abbreviation){
    tf50propval[i,"LUNG"] <- -log10(ourprotraining[ourprotraining$tissue %in% "t66-lung" & 
                                                     ourprotraining$feature_ID %in% ourpro,"p_value"])  
  }
  if("WAT-SC" %in% ourprotraining$tissue_abbreviation){
    tf50propval[i,"WAT-SC"] <- -log10(ourprotraining[ourprotraining$tissue %in% "t70-white-adipose" & 
                                                       ourprotraining$feature_ID %in% ourpro,"p_value"])  
  }
  
}

gastrotfprosig <- rownames(tfproanno[tfproanno$Gastro.Pro.ID %in% rownames(tf50propval)[tf50propval[,"SKM-GN"] > 1],])
hearttfprosig <- rownames(tfproanno[tfproanno$Gastro.Pro.ID %in% rownames(tf50propval)[tf50propval[,"HEART"] > 1],])
kidneytfprosig <- rownames(tfproanno[tfproanno$Gastro.Pro.ID %in% rownames(tf50propval)[tf50propval[,"KIDNEY"] > 1],])
livertfprosig <- rownames(tfproanno[tfproanno$Gastro.Pro.ID %in% rownames(tf50propval)[tf50propval[,"LIVER"] > 1],])
lungtfprosig <- rownames(tfproanno[tfproanno$Gastro.Pro.ID %in% rownames(tf50propval)[tf50propval[,"LUNG"] > 1],])
whitetfprosig <- rownames(tfproanno[tfproanno$Gastro.Pro.ID %in% rownames(tf50propval)[tf50propval[,"WAT-SC"] > 1],])


# SKM-GN 

gastrotftargetpromproxgenesigpct <- matrix(0L,nrow = length(gastrotfprosig),ncol = 1)
rownames(gastrotftargetpromproxgenesigpct) <- gastrotfprosig
colnames(gastrotftargetpromproxgenesigpct) <- "Percent.Significant"

gastrotftargetintrongenesigpct <- matrix(0L,nrow = length(gastrotfprosig),ncol = 1)
rownames(gastrotftargetintrongenesigpct) <- gastrotfprosig
colnames(gastrotftargetintrongenesigpct) <- "Percent.Significant"

gastrotftargetexongenesigpct <- matrix(0L,nrow = length(gastrotfprosig),ncol = 1)
rownames(gastrotftargetexongenesigpct) <- gastrotfprosig
colnames(gastrotftargetexongenesigpct) <- "Percent.Significant"

gastrotftargetintergenicgenesigpct <- matrix(0L,nrow = length(gastrotfprosig),ncol = 1)
rownames(gastrotftargetintergenicgenesigpct) <- gastrotfprosig
colnames(gastrotftargetintergenicgenesigpct) <- "Percent.Significant"


for(i in 1:length(gastrotfprosig)){
  ourtf <- gastrotfprosig[i]
  ourtftargetpeaks <- gastro50peakmotifs[gastro50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
  ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
  ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
  ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
  ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
  ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])
  gastrotftargetpromproxgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetpromproxgenes,gastrornasig))/length(ourtftargetpromproxgenes)
  gastrotftargetintrongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintrongenes,gastrornasig))/length(ourtftargetintrongenes)
  gastrotftargetexongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetexongenes,gastrornasig))/length(ourtftargetexongenes)
  gastrotftargetintergenicgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintergenicgenes,gastrornasig))/length(ourtftargetintergenicgenes)
}

# HEART 

hearttftargetpromproxgenesigpct <- matrix(0L,nrow = length(hearttfprosig),ncol = 1)
rownames(hearttftargetpromproxgenesigpct) <- hearttfprosig
colnames(hearttftargetpromproxgenesigpct) <- "Percent.Significant"

hearttftargetintrongenesigpct <- matrix(0L,nrow = length(hearttfprosig),ncol = 1)
rownames(hearttftargetintrongenesigpct) <- hearttfprosig
colnames(hearttftargetintrongenesigpct) <- "Percent.Significant"

hearttftargetexongenesigpct <- matrix(0L,nrow = length(hearttfprosig),ncol = 1)
rownames(hearttftargetexongenesigpct) <- hearttfprosig
colnames(hearttftargetexongenesigpct) <- "Percent.Significant"

hearttftargetintergenicgenesigpct <- matrix(0L,nrow = length(hearttfprosig),ncol = 1)
rownames(hearttftargetintergenicgenesigpct) <- hearttfprosig
colnames(hearttftargetintergenicgenesigpct) <- "Percent.Significant"


for(i in 1:length(hearttfprosig)){
  ourtf <- hearttfprosig[i]
  ourtftargetpeaks <- heart50peakmotifs[heart50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
  ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
  ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
  ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
  ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
  ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])
  hearttftargetpromproxgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetpromproxgenes,heartrnasig))/length(ourtftargetpromproxgenes)
  hearttftargetintrongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintrongenes,heartrnasig))/length(ourtftargetintrongenes)
  hearttftargetexongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetexongenes,heartrnasig))/length(ourtftargetexongenes)
  hearttftargetintergenicgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintergenicgenes,heartrnasig))/length(ourtftargetintergenicgenes)
}

# KIDNEY

kidneytftargetpromproxgenesigpct <- matrix(0L,nrow = length(kidneytfprosig),ncol = 1)
rownames(kidneytftargetpromproxgenesigpct) <- kidneytfprosig
colnames(kidneytftargetpromproxgenesigpct) <- "Percent.Significant"

kidneytftargetintrongenesigpct <- matrix(0L,nrow = length(kidneytfprosig),ncol = 1)
rownames(kidneytftargetintrongenesigpct) <- kidneytfprosig
colnames(kidneytftargetintrongenesigpct) <- "Percent.Significant"

kidneytftargetexongenesigpct <- matrix(0L,nrow = length(kidneytfprosig),ncol = 1)
rownames(kidneytftargetexongenesigpct) <- kidneytfprosig
colnames(kidneytftargetexongenesigpct) <- "Percent.Significant"

kidneytftargetintergenicgenesigpct <- matrix(0L,nrow = length(kidneytfprosig),ncol = 1)
rownames(kidneytftargetintergenicgenesigpct) <- kidneytfprosig
colnames(kidneytftargetintergenicgenesigpct) <- "Percent.Significant"


for(i in 1:length(kidneytfprosig)){
  ourtf <- kidneytfprosig[i]
  ourtftargetpeaks <- kidney50peakmotifs[kidney50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
  ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
  ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
  ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
  ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
  ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])
  kidneytftargetpromproxgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetpromproxgenes,kidneyrnasig))/length(ourtftargetpromproxgenes)
  kidneytftargetintrongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintrongenes,kidneyrnasig))/length(ourtftargetintrongenes)
  kidneytftargetexongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetexongenes,kidneyrnasig))/length(ourtftargetexongenes)
  kidneytftargetintergenicgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintergenicgenes,kidneyrnasig))/length(ourtftargetintergenicgenes)
}

# LIVER

livertftargetpromproxgenesigpct <- matrix(0L,nrow = length(livertfprosig),ncol = 1)
rownames(livertftargetpromproxgenesigpct) <- livertfprosig
colnames(livertftargetpromproxgenesigpct) <- "Percent.Significant"

livertftargetintrongenesigpct <- matrix(0L,nrow = length(livertfprosig),ncol = 1)
rownames(livertftargetintrongenesigpct) <- livertfprosig
colnames(livertftargetintrongenesigpct) <- "Percent.Significant"

livertftargetexongenesigpct <- matrix(0L,nrow = length(livertfprosig),ncol = 1)
rownames(livertftargetexongenesigpct) <- livertfprosig
colnames(livertftargetexongenesigpct) <- "Percent.Significant"

livertftargetintergenicgenesigpct <- matrix(0L,nrow = length(livertfprosig),ncol = 1)
rownames(livertftargetintergenicgenesigpct) <- livertfprosig
colnames(livertftargetintergenicgenesigpct) <- "Percent.Significant"


for(i in 1:length(livertfprosig)){
  ourtf <- livertfprosig[i]
  ourtftargetpeaks <- liver50peakmotifs[liver50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
  ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
  ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
  ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
  ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
  ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])
  livertftargetpromproxgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetpromproxgenes,liverrnasig))/length(ourtftargetpromproxgenes)
  livertftargetintrongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintrongenes,liverrnasig))/length(ourtftargetintrongenes)
  livertftargetexongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetexongenes,liverrnasig))/length(ourtftargetexongenes)
  livertftargetintergenicgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintergenicgenes,liverrnasig))/length(ourtftargetintergenicgenes)
}

# LUNG

lungtftargetpromproxgenesigpct <- matrix(0L,nrow = length(lungtfprosig),ncol = 1)
rownames(lungtftargetpromproxgenesigpct) <- lungtfprosig
colnames(lungtftargetpromproxgenesigpct) <- "Percent.Significant"

lungtftargetintrongenesigpct <- matrix(0L,nrow = length(lungtfprosig),ncol = 1)
rownames(lungtftargetintrongenesigpct) <- lungtfprosig
colnames(lungtftargetintrongenesigpct) <- "Percent.Significant"

lungtftargetexongenesigpct <- matrix(0L,nrow = length(lungtfprosig),ncol = 1)
rownames(lungtftargetexongenesigpct) <- lungtfprosig
colnames(lungtftargetexongenesigpct) <- "Percent.Significant"

lungtftargetintergenicgenesigpct <- matrix(0L,nrow = length(lungtfprosig),ncol = 1)
rownames(lungtftargetintergenicgenesigpct) <- lungtfprosig
colnames(lungtftargetintergenicgenesigpct) <- "Percent.Significant"


for(i in 1:length(lungtfprosig)){
  ourtf <- lungtfprosig[i]
  ourtftargetpeaks <- lung50peakmotifs[lung50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
  ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
  ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
  ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
  ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
  ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])
  lungtftargetpromproxgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetpromproxgenes,lungrnasig))/length(ourtftargetpromproxgenes)
  lungtftargetintrongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintrongenes,lungrnasig))/length(ourtftargetintrongenes)
  lungtftargetexongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetexongenes,lungrnasig))/length(ourtftargetexongenes)
  lungtftargetintergenicgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintergenicgenes,lungrnasig))/length(ourtftargetintergenicgenes)
}

# WAT-SC

whitetftargetpromproxgenesigpct <- matrix(0L,nrow = length(whitetfprosig),ncol = 1)
rownames(whitetftargetpromproxgenesigpct) <- whitetfprosig
colnames(whitetftargetpromproxgenesigpct) <- "Percent.Significant"

whitetftargetintrongenesigpct <- matrix(0L,nrow = length(whitetfprosig),ncol = 1)
rownames(whitetftargetintrongenesigpct) <- whitetfprosig
colnames(whitetftargetintrongenesigpct) <- "Percent.Significant"

whitetftargetexongenesigpct <- matrix(0L,nrow = length(whitetfprosig),ncol = 1)
rownames(whitetftargetexongenesigpct) <- whitetfprosig
colnames(whitetftargetexongenesigpct) <- "Percent.Significant"

whitetftargetintergenicgenesigpct <- matrix(0L,nrow = length(whitetfprosig),ncol = 1)
rownames(whitetftargetintergenicgenesigpct) <- whitetfprosig
colnames(whitetftargetintergenicgenesigpct) <- "Percent.Significant"


for(i in 1:length(whitetfprosig)){
  ourtf <- whitetfprosig[i]
  ourtftargetpeaks <- white50peakmotifs[white50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
  ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
  ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
  ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
  ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
  ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])
  whitetftargetpromproxgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetpromproxgenes,whiternasig))/length(ourtftargetpromproxgenes)
  whitetftargetintrongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintrongenes,whiternasig))/length(ourtftargetintrongenes)
  whitetftargetexongenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetexongenes,whiternasig))/length(ourtftargetexongenes)
  whitetftargetintergenicgenesigpct[i,"Percent.Significant"] <- length(intersect(ourtftargetintergenicgenes,whiternasig))/length(ourtftargetintergenicgenes)
}


gastroprotargetsigdf <- data.frame("Percent.Significant" = c(gastrotftargetpromproxgenesigpct,
                                                             gastrotftargetintrongenesigpct,
                                                             gastrotftargetexongenesigpct,
                                                             gastrotftargetintergenicgenesigpct),
                                   "Region" = c(rep("Promoter (<=1kb)",length(gastrotftargetpromproxgenesigpct)),
                                                rep("Intron",length(gastrotftargetpromproxgenesigpct)),
                                                rep("Exon",length(gastrotftargetpromproxgenesigpct)),
                                                rep("Distal Intergenic",length(gastrotftargetpromproxgenesigpct))),
                                   "Transcription Factor" = rep(gsub("\\(.*","",rownames(gastrotftargetpromproxgenesigpct)),4))

heartprotargetsigdf <- data.frame("Percent.Significant" = c(hearttftargetpromproxgenesigpct,
                                                            hearttftargetintrongenesigpct,
                                                            hearttftargetexongenesigpct,
                                                            hearttftargetintergenicgenesigpct),
                                  "Region" = c(rep("Promoter (<=1kb)",length(hearttftargetpromproxgenesigpct)),
                                               rep("Intron",length(hearttftargetpromproxgenesigpct)),
                                               rep("Exon",length(hearttftargetpromproxgenesigpct)),
                                               rep("Distal Intergenic",length(hearttftargetpromproxgenesigpct))),
                                  "Transcription Factor" = rep(gsub("\\(.*","",rownames(hearttftargetpromproxgenesigpct)),4))


kidneyprotargetsigdf <- data.frame("Percent.Significant" = c(kidneytftargetpromproxgenesigpct,
                                                             kidneytftargetintrongenesigpct,
                                                             kidneytftargetexongenesigpct,
                                                             kidneytftargetintergenicgenesigpct),
                                   "Region" = c(rep("Promoter (<=1kb)",length(kidneytftargetpromproxgenesigpct)),
                                                rep("Intron",length(kidneytftargetpromproxgenesigpct)),
                                                rep("Exon",length(kidneytftargetpromproxgenesigpct)),
                                                rep("Distal Intergenic",length(kidneytftargetpromproxgenesigpct))),
                                   "Transcription Factor" = rep(gsub("\\(.*","",rownames(kidneytftargetpromproxgenesigpct)),4))

liverprotargetsigdf <- data.frame("Percent.Significant" = c(livertftargetpromproxgenesigpct,
                                                            livertftargetintrongenesigpct,
                                                            livertftargetexongenesigpct,
                                                            livertftargetintergenicgenesigpct),
                                  "Region" = c(rep("Promoter (<=1kb)",length(livertftargetpromproxgenesigpct)),
                                               rep("Intron",length(livertftargetpromproxgenesigpct)),
                                               rep("Exon",length(livertftargetpromproxgenesigpct)),
                                               rep("Distal Intergenic",length(livertftargetpromproxgenesigpct))),
                                  "Transcription Factor" = rep(gsub("\\(.*","",rownames(livertftargetpromproxgenesigpct)),4))

lungprotargetsigdf <- data.frame("Percent.Significant" = c(lungtftargetpromproxgenesigpct,
                                                           lungtftargetintrongenesigpct,
                                                           lungtftargetexongenesigpct,
                                                           lungtftargetintergenicgenesigpct),
                                 "Region" = c(rep("Promoter (<=1kb)",length(lungtftargetpromproxgenesigpct)),
                                              rep("Intron",length(lungtftargetpromproxgenesigpct)),
                                              rep("Exon",length(lungtftargetpromproxgenesigpct)),
                                              rep("Distal Intergenic",length(lungtftargetpromproxgenesigpct))),
                                 "Transcription Factor" = rep(gsub("\\(.*","",rownames(lungtftargetpromproxgenesigpct)),4))

whiteprotargetsigdf <- data.frame("Percent.Significant" = c(whitetftargetpromproxgenesigpct,
                                                            whitetftargetintrongenesigpct,
                                                            whitetftargetexongenesigpct,
                                                            whitetftargetintergenicgenesigpct),
                                  "Region" = c(rep("Promoter (<=1kb)",length(whitetftargetpromproxgenesigpct)),
                                               rep("Intron",length(whitetftargetpromproxgenesigpct)),
                                               rep("Exon",length(whitetftargetpromproxgenesigpct)),
                                               rep("Distal Intergenic",length(whitetftargetpromproxgenesigpct))),
                                  "Transcription Factor" = rep(gsub("\\(.*","",rownames(whitetftargetpromproxgenesigpct)),4))
# Supplemental Figure S8 and Figure 4E
pdf(file = "Supplemental Figure S8A.pdf",width = 7,height = 4.5)
ggplot(gastroprotargetsigdf[!gastroprotargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(gastrornasig)/dim(gastrol2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 315, vjust = 0.5,hjust = 0),axis.title = element_text(size = 13),legend.text = element_text(size = 12)) + ylim(0.00,0.15)
dev.off()

pdf(file = "Figure 4E.pdf",width = 7,height = 4.5)
ggplot(gastroprotargetsigdf[!gastroprotargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(gastrornasig)/dim(gastrol2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 15),axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),axis.title = element_text(size = 15),legend.text = element_text(size = 15)) + ylim(0.00,0.15)
dev.off()

pdf(file = "Supplemental Figure S8B.pdf",width = 7,height = 4.5)
ggplot(heartprotargetsigdf[!heartprotargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(heartrnasig)/dim(heartl2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 315, vjust = 0.5,hjust = 0),axis.title = element_text(size = 13),legend.text = element_text(size = 12)) + ylim(0.00,0.10)
dev.off()

pdf(file = "Supplemental Figure S8C.pdf",width = 7,height = 4.5)
ggplot(kidneyprotargetsigdf[!kidneyprotargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(kidneyrnasig)/dim(kidneyl2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 315, vjust = 0.5,hjust = 0),axis.title = element_text(size = 13),legend.text = element_text(size = 12)) + ylim(0.00,0.05)
dev.off()

pdf(file = "Supplemental Figure S8D.pdf",width = 7,height = 4.5)
ggplot(liverprotargetsigdf[!liverprotargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(liverrnasig)/dim(liverl2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 315, vjust = 0.5,hjust = 0),axis.title = element_text(size = 13),legend.text = element_text(size = 12)) + ylim(0.00,0.10)
dev.off()

pdf(file = "Supplemental Figure S8E.pdf",width = 7,height = 4.5)
ggplot(lungprotargetsigdf[!lungprotargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(lungrnasig)/dim(lungl2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 315, vjust = 0.5,hjust = 0),axis.title = element_text(size = 13),legend.text = element_text(size = 12)) + ylim(0.00,0.15)
dev.off()

pdf(file = "Supplemental Figure S8F.pdf",width = 7,height = 4.5)
ggplot(whiteprotargetsigdf[!whiteprotargetsigdf$Region %in% "Exon",],aes(x=Transcription.Factor,y=Percent.Significant,group=Region,palette="jco")) + geom_hline(yintercept = length(whiternasig)/dim(whitel2fcmat)[1],linetype = "dashed",size = 1) + geom_point(aes(color=Region),size = 3) + geom_line(aes(color=Region),size = 2) + theme_classic() + scale_color_manual(values = c("#BCBD22FF","#8C564BFF","#FF7F0EFF")) + theme(legend.position = "top",axis.text = element_text(size = 12),axis.text.x = element_text(angle = 315, vjust = 0.5,hjust = 0),axis.title = element_text(size = 13),legend.text = element_text(size = 12)) + ylim(0.00,0.20)
dev.off()


# Figure 4F

ourtf <- "Mef2c(MADS)/GM12878-Mef2c-ChIP-Seq(GSE32465)/Homer"
ourtftargetpeaks <- gastro50peakmotifs[gastro50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])

sextimemeta <- data.frame(row.names = colnames(gastrol2fcmat),
                          "Sex" = c(rep("Female",4),rep("Male",4)),
                          "Week" = rep(c("1w","2w","4w","8w"),2))
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

pdf(file = "Figure 4F.pdf",width = 7,height = 4)
pheatmap(gastrol2fcmat[intersect(ourtftargetpromproxgenes,gastrornasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,gastrornasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,annotation_row = ourtftargetpromproxcormeta,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F)
dev.off()

# Figure 4G
ourtf <- "Nur77(NR)/K562-NR4A1-ChIP-Seq(GSE31363)/Homer"
ourtftargetpeaks <- gastro50peakmotifs[gastro50peakmotifs$Motif.Name %in% ourtf,"PositionID"]
ourtftargetgenes <- unique(peakanno[ourtftargetpeaks,"ensembl_gene"])
ourtftargetpromproxgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Promoter (<=1kb)",])),"ensembl_gene"])
ourtftargetintrongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Intron",])),"ensembl_gene"])
ourtftargetexongenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Exon",])),"ensembl_gene"])
ourtftargetintergenicgenes <- unique(peakanno[intersect(ourtftargetpeaks,rownames(peakanno[peakanno$custom_annotation %in% "Distal Intergenic",])),"ensembl_gene"])

sextimemeta <- data.frame(row.names = colnames(gastrol2fcmat),
                          "Sex" = c(rep("Female",4),rep("Male",4)),
                          "Week" = rep(c("1w","2w","4w","8w"),2))
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

pdf(file = "Figure 4G.pdf",width = 7,height = 4)
pheatmap(gastrol2fcmat[intersect(ourtftargetpromproxgenes,gastrornasig),],cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = "0",labels_row = enstosym[intersect(ourtftargetpromproxgenes,gastrornasig),"Symbol"],display_numbers = T,number_color = "black",annotation_col = sextimemeta[,c("Week","Sex")],annotation_colors = ann_cols_procor,cluster_rows = F,annotation_row = ourtftargetpromproxcormeta,show_colnames = F,cellwidth = 30,fontsize = 15,legend = F,annotation_legend = F,fontsize_number = 10,annotation_names_row = F)
dev.off()

####
# Figure 5
#####

####
# We generated a normalized control matrix for the RNAseq data across tissues
#####

bloodrnadata <- read.delim(file = "PASS1B Transcription Factor Paper Data/RNASeq Data/pass1b-06_transcript-rna-seq_raw-data_t30-blood-rna.txt",header = TRUE, row.names = 1,sep = "\t")
colnames(bloodrnadata) <- gsub("X","",colnames(bloodrnadata))
bloodrnadata <- bloodrnadata[,colnames(bloodrnadata) %in% pass1bphenodata$pass1bf0004]
for(i in 1:50){
  rnacol <- colnames(bloodrnadata)[i]
  colnames(bloodrnadata)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 == rnacol,"pass1bf0001"]
}

hippornadata <- read.delim(file = "PASS1B Transcription Factor Paper Data/RNASeq Data/pass1b-06_transcript-rna-seq_raw-data_t52-hippocampus.txt",header = TRUE, row.names = 1,sep = "\t")
colnames(hippornadata) <- gsub("X","",colnames(hippornadata))
hippornadata <- hippornadata[,colnames(hippornadata) %in% pass1bphenodata$pass1bf0004]
for(i in 1:50){
  rnacol <- colnames(hippornadata)[i]
  colnames(hippornadata)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 == rnacol,"pass1bf0001"]
}

cortexrnadata <- read.delim(file = "PASS1B Transcription Factor Paper Data/RNASeq Data/pass1b-06_transcript-rna-seq_raw-data_t53-cortex.txt",header = TRUE, row.names = 1,sep = "\t")
colnames(cortexrnadata) <- gsub("X","",colnames(cortexrnadata))
cortexrnadata <- cortexrnadata[,colnames(cortexrnadata) %in% pass1bphenodata$pass1bf0004]
for(i in 1:50){
  rnacol <- colnames(cortexrnadata)[i]
  colnames(cortexrnadata)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 == rnacol,"pass1bf0001"]
}

hypornadata <- read.delim(file = "PASS1B Transcription Factor Paper Data/RNASeq Data/pass1b-06_transcript-rna-seq_raw-data_t54-hypothalamus.txt",header = TRUE, row.names = 1,sep = "\t")
colnames(hypornadata) <- gsub("X","",colnames(hypornadata))
hypornadata <- hypornadata[,colnames(hypornadata) %in% pass1bphenodata$pass1bf0004]
for(i in 1:50){
  rnacol <- colnames(hypornadata)[i]
  colnames(hypornadata)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 == rnacol,"pass1bf0001"]
}

gastrornadata <- read.delim(file = "PASS1B Transcription Factor Paper Data/RNASeq Data/pass1b-06_transcript-rna-seq_raw-data_t55-gastrocnemius.txt",header = TRUE, row.names = 1,sep = "\t")
colnames(gastrornadata) <- gsub("X","",colnames(gastrornadata))
gastrornadata <- gastrornadata[,colnames(gastrornadata) %in% pass1bphenodata$pass1bf0004]
for(i in 1:50){
  rnacol <- colnames(gastrornadata)[i]
  colnames(gastrornadata)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 == rnacol,"pass1bf0001"]
}

vastusrnadata <- read.delim(file = "PASS1B Transcription Factor Paper Data/RNASeq Data/pass1b-06_transcript-rna-seq_raw-data_t56-vastus-lateralis.txt",header = TRUE, row.names = 1,sep = "\t")
colnames(vastusrnadata) <- gsub("X","",colnames(vastusrnadata))
vastusrnadata <- vastusrnadata[,colnames(vastusrnadata) %in% pass1bphenodata$pass1bf0004]
for(i in 1:50){
  rnacol <- colnames(vastusrnadata)[i]
  colnames(vastusrnadata)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 == rnacol,"pass1bf0001"]
}

heartrnadata <- read.delim(file = "PASS1B Transcription Factor Paper Data/RNASeq Data/pass1b-06_transcript-rna-seq_raw-data_t58-heart.txt",header = TRUE, row.names = 1,sep = "\t")
colnames(heartrnadata) <- gsub("X","",colnames(heartrnadata))
heartrnadata <- heartrnadata[,colnames(heartrnadata) %in% pass1bphenodata$pass1bf0004]
for(i in 1:50){
  rnacol <- colnames(heartrnadata)[i]
  colnames(heartrnadata)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 == rnacol,"pass1bf0001"]
}

kidneyrnadata <- read.delim(file = "PASS1B Transcription Factor Paper Data/RNASeq Data/pass1b-06_transcript-rna-seq_raw-data_t59-kidney.txt",header = TRUE, row.names = 1,sep = "\t")
colnames(kidneyrnadata) <- gsub("X","",colnames(kidneyrnadata))
kidneyrnadata <- kidneyrnadata[,colnames(kidneyrnadata) %in% pass1bphenodata$pass1bf0004]
for(i in 1:50){
  rnacol <- colnames(kidneyrnadata)[i]
  colnames(kidneyrnadata)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 == rnacol,"pass1bf0001"]
}

adrenalrnadata <- read.delim(file = "PASS1B Transcription Factor Paper Data/RNASeq Data/pass1b-06_transcript-rna-seq_raw-data_t60-adrenal.txt",header = TRUE, row.names = 1,sep = "\t")
colnames(adrenalrnadata) <- gsub("X","",colnames(adrenalrnadata))
adrenalrnadata <- adrenalrnadata[,colnames(adrenalrnadata) %in% pass1bphenodata$pass1bf0004]
for(i in 1:50){
  rnacol <- colnames(adrenalrnadata)[i]
  colnames(adrenalrnadata)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 == rnacol,"pass1bf0001"]
}

colonrnadata <- read.delim(file = "PASS1B Transcription Factor Paper Data/RNASeq Data/pass1b-06_transcript-rna-seq_raw-data_t61-colon.txt",header = TRUE, row.names = 1,sep = "\t")
colnames(colonrnadata) <- gsub("X","",colnames(colonrnadata))
colonrnadata <- colonrnadata[,colnames(colonrnadata) %in% pass1bphenodata$pass1bf0004]
for(i in 1:50){
  rnacol <- colnames(colonrnadata)[i]
  colnames(colonrnadata)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 == rnacol,"pass1bf0001"]
}

spleenrnadata <- read.delim(file = "PASS1B Transcription Factor Paper Data/RNASeq Data/pass1b-06_transcript-rna-seq_raw-data_t62-spleen.txt",header = TRUE, row.names = 1,sep = "\t")
colnames(spleenrnadata) <- gsub("X","",colnames(spleenrnadata))
spleenrnadata <- spleenrnadata[,colnames(spleenrnadata) %in% pass1bphenodata$pass1bf0004]
for(i in 1:50){
  rnacol <- colnames(spleenrnadata)[i]
  colnames(spleenrnadata)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 == rnacol,"pass1bf0001"]
}

testesrnadata <- read.delim(file = "PASS1B Transcription Factor Paper Data/RNASeq Data/pass1b-06_transcript-rna-seq_raw-data_t63-testes.txt",header = TRUE, row.names = 1,sep = "\t")
colnames(testesrnadata) <- gsub("X","",colnames(testesrnadata))
testesrnadata <- testesrnadata[,colnames(testesrnadata) %in% pass1bphenodata$pass1bf0004]
for(i in 1:25){
  rnacol <- colnames(testesrnadata)[i]
  colnames(testesrnadata)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 == rnacol,"pass1bf0001"]
}

ovariesrnadata <- read.delim(file = "PASS1B Transcription Factor Paper Data/RNASeq Data/pass1b-06_transcript-rna-seq_raw-data_t64-ovaries.txt",header = TRUE, row.names = 1,sep = "\t")
colnames(ovariesrnadata) <- gsub("X","",colnames(ovariesrnadata))
ovariesrnadata <- ovariesrnadata[,colnames(ovariesrnadata) %in% pass1bphenodata$pass1bf0004]
for(i in 1:24){
  rnacol <- colnames(ovariesrnadata)[i]
  colnames(ovariesrnadata)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 == rnacol,"pass1bf0001"]
}

aortarnadata <- read.delim(file = "PASS1B Transcription Factor Paper Data/RNASeq Data/pass1b-06_transcript-rna-seq_raw-data_t65-aorta.txt",header = TRUE, row.names = 1,sep = "\t")
colnames(aortarnadata) <- gsub("X","",colnames(aortarnadata))
aortarnadata <- aortarnadata[,colnames(aortarnadata) %in% pass1bphenodata$pass1bf0004]
for(i in 1:50){
  rnacol <- colnames(aortarnadata)[i]
  colnames(aortarnadata)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 == rnacol,"pass1bf0001"]
}

lungrnadata <- read.delim(file = "PASS1B Transcription Factor Paper Data/RNASeq Data/pass1b-06_transcript-rna-seq_raw-data_t66-lung.txt",header = TRUE, row.names = 1,sep = "\t")
colnames(lungrnadata) <- gsub("X","",colnames(lungrnadata))
lungrnadata <- lungrnadata[,colnames(lungrnadata) %in% pass1bphenodata$pass1bf0004]
for(i in 1:50){
  rnacol <- colnames(lungrnadata)[i]
  colnames(lungrnadata)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 == rnacol,"pass1bf0001"]
}

smallintrnadata <- read.delim(file = "PASS1B Transcription Factor Paper Data/RNASeq Data/pass1b-06_transcript-rna-seq_raw-data_t67-small-intestine.txt",header = TRUE, row.names = 1,sep = "\t")
colnames(smallintrnadata) <- gsub("X","",colnames(smallintrnadata))
smallintrnadata <- smallintrnadata[,colnames(smallintrnadata) %in% pass1bphenodata$pass1bf0004]
for(i in 1:50){
  rnacol <- colnames(smallintrnadata)[i]
  colnames(smallintrnadata)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 == rnacol,"pass1bf0001"]
}

liverrnadata <- read.delim(file = "PASS1B Transcription Factor Paper Data/RNASeq Data/pass1b-06_transcript-rna-seq_raw-data_t68-liver.txt",header = TRUE, row.names = 1,sep = "\t")
colnames(liverrnadata) <- gsub("X","",colnames(liverrnadata))
liverrnadata <- liverrnadata[,colnames(liverrnadata) %in% pass1bphenodata$pass1bf0004]
for(i in 1:50){
  rnacol <- colnames(liverrnadata)[i]
  colnames(liverrnadata)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 == rnacol,"pass1bf0001"]
}

brownrnadata <- read.delim(file = "PASS1B Transcription Factor Paper Data/RNASeq Data/pass1b-06_transcript-rna-seq_raw-data_t69-brown-adipose.txt",header = TRUE, row.names = 1,sep = "\t")
colnames(brownrnadata) <- gsub("X","",colnames(brownrnadata))
brownrnadata <- brownrnadata[,colnames(brownrnadata) %in% pass1bphenodata$pass1bf0004]
for(i in 1:50){
  rnacol <- colnames(brownrnadata)[i]
  colnames(brownrnadata)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 == rnacol,"pass1bf0001"]
}

whiternadata <- read.delim(file = "PASS1B Transcription Factor Paper Data/RNASeq Data/pass1b-06_transcript-rna-seq_raw-data_t70-white-adipose.txt",header = TRUE, row.names = 1,sep = "\t")
colnames(whiternadata) <- gsub("X","",colnames(whiternadata))
whiternadata <- whiternadata[,colnames(whiternadata) %in% pass1bphenodata$pass1bf0004]
for(i in 1:50){
  rnacol <- colnames(whiternadata)[i]
  colnames(whiternadata)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 == rnacol,"pass1bf0001"]
}

oursamplelist <- colnames(adrenalrnadata)
controlsamples <- as.character(intersect(oursamplelist,unique(pass1bphenodata[pass1bphenodata$pass1bf0011 %in% "Eight-week program Control Group","pass1bf0001"])))

rnacolsummat <- cbind(apply(adrenalrnadata[,controlsamples],2,sum),
                      apply(aortarnadata[,controlsamples],2,sum),
                      apply(bloodrnadata[,controlsamples],2,sum),
                      apply(brownrnadata[,controlsamples],2,sum),
                      apply(colonrnadata[,controlsamples],2,sum),
                      apply(cortexrnadata[,controlsamples],2,sum),
                      apply(gastrornadata[,controlsamples],2,sum),
                      apply(heartrnadata[,controlsamples],2,sum),
                      apply(hippornadata[,controlsamples],2,sum),
                      apply(hypornadata[,controlsamples],2,sum),
                      apply(kidneyrnadata[,controlsamples],2,sum),
                      apply(liverrnadata[,controlsamples],2,sum),
                      apply(lungrnadata[,controlsamples],2,sum),
                      apply(smallintrnadata[,controlsamples],2,sum),
                      apply(spleenrnadata[,controlsamples],2,sum),
                      apply(vastusrnadata[,controlsamples],2,sum),
                      apply(whiternadata[,controlsamples],2,sum))
colnames(rnacolsummat) <- c("Adrenal","Aorta","Blood","Brown Ad","Colon","Cortex","Gastro",
                            "Heart","Hippo","Hypo","Kidney","Liver","Lung","SmallInt",
                            "Spleen","Vastus","White")
controlsamplemeta <- data.frame(row.names = controlsamples,"sex" = rep("Female",10))
for(i in 1:10){
  oursample <- controlsamples[i]
  pass1bphenodata[pass1bphenodata$pass1bf0001 %in% oursample,"pass1bf0027"]
  if(unique(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% oursample,"pass1bf0027"]) == 2){
    controlsamplemeta[oursample,"sex"] <- "Male"
  }
}

sortcontrol <- rownames(controlsamplemeta)[order(controlsamplemeta$sex)]

rnacontroldata <- cbind(adrenalrnadata[,sortcontrol],aortarnadata[,sortcontrol],bloodrnadata[,sortcontrol],brownrnadata[,sortcontrol],
                        colonrnadata[,sortcontrol],cortexrnadata[,sortcontrol],gastrornadata[,sortcontrol],heartrnadata[,sortcontrol],
                        hippornadata[,sortcontrol],hypornadata[,sortcontrol],kidneyrnadata[,sortcontrol],liverrnadata[,sortcontrol],
                        lungrnadata[,sortcontrol],ovariesrnadata[,sortcontrol[1:5]],smallintrnadata[,sortcontrol],
                        spleenrnadata[,sortcontrol],testesrnadata[,sortcontrol[6:10]],vastusrnadata[,sortcontrol],whiternadata[,sortcontrol])
colnames(rnacontroldata) <- c(paste("Adrenal_",sortcontrol,sep = ""),
                              paste("Aorta_",sortcontrol,sep = ""),
                              paste("Blood_",sortcontrol,sep = ""),
                              paste("BrownAd_",sortcontrol,sep = ""),
                              paste("Colon_",sortcontrol,sep = ""),
                              paste("Cortex_",sortcontrol,sep = ""),
                              paste("Gastro_",sortcontrol,sep = ""),
                              paste("Heart_",sortcontrol,sep = ""),
                              paste("Hippo_",sortcontrol,sep = ""),
                              paste("Hypo_",sortcontrol,sep = ""),
                              paste("Kidney_",sortcontrol,sep = ""),
                              paste("Liver_",sortcontrol,sep = ""),
                              paste("Lung_",sortcontrol,sep = ""),
                              paste("Ovaries_",sortcontrol[1:5],sep = ""),
                              paste("SmallInt_",sortcontrol,sep = ""),
                              paste("Spleen_",sortcontrol,sep = ""),
                              paste("Testes_",sortcontrol[6:10],sep = ""),
                              paste("Vastus_",sortcontrol,sep = ""),
                              paste("WhiteAd_",sortcontrol,sep = ""))
rnametadata <- data.frame(row.names = colnames(rnacontroldata),
                          "Tissue" = c(rep("Adrenal",10),
                                       rep("Aorta",10),
                                       rep("Blood",10),
                                       rep("BrownAd",10),
                                       rep("Colon",10),
                                       rep("Cortex",10),
                                       rep("Gastro",10),
                                       rep("Heart",10),
                                       rep("Hippo",10),
                                       rep("Hypo",10),
                                       rep("Kidney",10),
                                       rep("Liver",10),
                                       rep("Lung",10),
                                       rep("Ovaries",5),
                                       rep("SmallInt",10),
                                       rep("Spleen",10),
                                       rep("Testes",5),
                                       rep("Vastus",10),
                                       rep("WhiteAd",10)),
                          "Sex" = c(rep(c(rep("Female",5),rep("Male",5)),13),
                                    rep("Female",5),
                                    rep(c(rep("Female",5),rep("Male",5)),2),
                                    rep("Male",5),
                                    rep(c(rep("Female",5),rep("Male",5)),2)),
                          "Sample" = c(rep(sortcontrol,13),
                                       sortcontrol[1:5],
                                       rep(sortcontrol,2),
                                       sortcontrol[6:10],
                                       rep(sortcontrol,2)))

# raw --> filter (filt)
raw_dge = DGEList(counts=rnacontroldata) # raw counts filtered down to samples in this tissue 
keep = rowSums(cpm(raw_dge) > 0.5) >= 2
filt_dge = raw_dge[keep, , keep.lib.sizes=FALSE]
# filt --> tmm 
dge = calcNormFactors(filt_dge, method='TMM')
rnacontrolnorm = cpm(dge,log=TRUE)

controlcolumns <- c(colnames(rnacontrolnorm)[grep("Gastro",colnames(rnacontrolnorm))],
                    colnames(rnacontrolnorm)[grep("Heart",colnames(rnacontrolnorm))],
                    colnames(rnacontrolnorm)[grep("Hippo",colnames(rnacontrolnorm))],
                    colnames(rnacontrolnorm)[grep("Kidney",colnames(rnacontrolnorm))],
                    colnames(rnacontrolnorm)[grep("Liver",colnames(rnacontrolnorm))],
                    colnames(rnacontrolnorm)[grep("Lung",colnames(rnacontrolnorm))],
                    colnames(rnacontrolnorm)[grep("BrownAd",colnames(rnacontrolnorm))],
                    colnames(rnacontrolnorm)[grep("WhiteAd",colnames(rnacontrolnorm))])

rnacontrolmeta <- data.frame(row.names = colnames(rnacontrolnorm),
                             "Tissue" = c(rep("ADRNL",10),
                                          rep("VENACV",10),
                                          rep("BLOOD",10),
                                          rep("BAT",10),
                                          rep("COLON",10),
                                          rep("CORTEX",10),
                                          rep("SKM-GN",10),
                                          rep("HEART",10),
                                          rep("HIPPOC",10),
                                          rep("HYPOTH",10),
                                          rep("KIDNEY",10),
                                          rep("LIVER",10),
                                          rep("LUNG",10),
                                          rep("OVARY",5),
                                          rep("SMLINT",10),
                                          rep("SPLEEN",10),
                                          rep("TESTES",5),
                                          rep("SKM-VL",10),
                                          rep("WAT-SC",10)),
                             "Sex" = c(rep(c(rep("Female",5),rep("Male",5)),13),
                                       rep("Female",5),
                                       rep(c(rep("Female",5),rep("Male",5)),2),
                                       rep("Male",5),
                                       rep(c(rep("Female",5),rep("Male",5)),2)),
                             "Sample" = gsub(".*_","",colnames(rnacontrolnorm)))

####
# Figure 5A
#####

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

# Figure 5A
pdf(file = "Figure 5A.pdf",width = 6,height = 6)
pheatmap(tfatacsigpmat[finalatacsigtfs,],cluster_rows = F,cluster_cols = F,labels_row = gsub("\\(.*","",finalatacsigtfs),annotation_col = tissuemeta,annotation_colors = ann_cols,show_colnames = F,color = colorpanel(101,"white","firebrick"))
dev.off()


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




# Figure 5B - 592x580
pdf(file = "Figure 5B.pdf",width = 6,height = 6)
pheatmap(tf50logpmat[finalsigtflist,],labels_row = gsub("\\(.*","",finalsigtflist),angle_col = 0,color = colorpanel(101,"white","firebrick"),breaks = seq(0,40,length.out = 101),cluster_cols = F,cluster_rows = F,annotation_col = tissuemeta,annotation_colors = ann_colsheatmaptrim,show_colnames = F)
dev.off()


####
# Figure 5C - TO BE ADDED BY BINGQING
#####

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
motifCountByRegioncurrsigmetatrim <- motifCountByRegioncurrsigmeta[(!motifCountByRegioncurrsigmeta$Region %in% "Overlaps.Gene") & (!motifCountByRegioncurrsigmeta$Tissue %in% c("HIPPOC","WAT-SC")),]

motifCountByRegioncurrsigtrim <- motifCountByRegioncurrsig[,rownames(motifCountByRegioncurrsigmetatrim)]
motifCountByRegioncurrsigtrim <- motifCountByRegioncurrsigtrim[apply(motifCountByRegioncurrsigtrim,1,min) < 0.05,]

pdf("Figure 5C.pdf",width = 7,height = 13)
pheatmap(-log10(motifCountByRegioncurrsigtrim[,rownames(motifCountByRegioncurrsigmetatrim)]),labels_row = gsub("\\/.*","",gsub("\\(.*","",rownames(motifCountByRegioncurrsigtrim))),show_colnames = F,annotation_col = motifCountByRegioncurrsigmetatrim[,c("Region","Tissue")],breaks = seq(0,6,length.out = 101),color = colorpanel(101,"white","firebrick"),cluster_cols = F,annotation_colors = ann_cols)
dev.off()


####
# Figure 5D-E
#####

tf50allregionpctmat <- cbind(as.numeric(gsub("%","",gastrosigpromproxtf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",gastrosigpromfartf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",gastrodisttf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",gastrointtf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",gastroextf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",gastrouptf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",gastrodowntf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",heartsigpromproxtf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",heartsigpromfartf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",heartdisttf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",heartinttf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",heartextf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",heartuptf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",heartdowntf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",hipposigpromproxtf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",hipposigpromfartf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",hippodisttf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",hippointtf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",hippoextf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",hippouptf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",hippodowntf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",kidneysigpromproxtf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",kidneysigpromfartf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",kidneydisttf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",kidneyinttf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",kidneyextf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",kidneyuptf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",kidneydowntf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",liversigpromproxtf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",liversigpromfartf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",liverdisttf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",liverinttf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",liverextf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",liveruptf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",liverdowntf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",lungsigpromproxtf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",lungsigpromfartf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",lungdisttf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",lunginttf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",lungextf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",lunguptf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",lungdowntf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",whitesigpromproxtf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",whitesigpromfartf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",whitedisttf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",whiteinttf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",whiteextf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",whiteuptf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",whitedowntf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",brownsigpromproxtf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",brownsigpromfartf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",browndisttf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",browninttf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",brownextf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",brownuptf[tflist,"X..of.Target.Sequences.with.Motif"])),
                             as.numeric(gsub("%","",browndowntf[tflist,"X..of.Target.Sequences.with.Motif"])))

rownames(tf50allregionpctmat) <- tflist
colnames(tf50allregionpctmat) <- c("SKM-GN_promprox","SKM-GN_promfar","SKM-GN_dist","SKM-GN_int","SKM-GN_ex","SKM-GN_up","SKM-GN_down",
                                   "HEART_promprox","HEART_promfar","HEART_dist","HEART_int","HEART_ex","HEART_up","HEART_down",
                                   "HIPPOC_promprox","HIPPOC_promfar","HIPPOC_dist","HIPPOC_int","HIPPOC_ex","HIPPOC_up","HIPPOC_down",
                                   "KIDNEY_promprox","KIDNEY_promfar","KIDNEY_dist","KIDNEY_int","KIDNEY_ex","KIDNEY_up","KIDNEY_down",
                                   "LIVER_promprox","LIVER_promfar","LIVER_dist","LIVER_int","LIVER_ex","LIVER_up","LIVER_down",
                                   "LUNG_promprox","LUNG_promfar","LUNG_dist","LUNG_int","LUNG_ex","LUNG_up","LUNG_down",
                                   "BAT_promprox","BAT_promfar","BAT_dist","BAT_int","BAT_ex","BAT_up","BAT_down",
                                   "WAT-SC_promprox","WAT-SC_promfar","WAT-SC_dist","WAT-SC_int","WAT-SC_ex","WAT-SC_up","WAT-SC_down")
tf50allregionpctmat[is.na(tf50allregionpctmat)] <- 0

tf50allregionmeta <- data.frame(row.names = colnames(tf50allregionpctmat),
                                "Tissue" = c(rep("SKM-GN",7),rep("HEART",7),rep("HIPPOC",7),rep("KIDNEY",7),rep("LIVER",7),
                                             rep("LUNG",7),rep("BAT",7),rep("WAT-SC",7)),
                                "Region" = rep(c("Promoter.(<=1kb)","Promoter.(1-2kb)","Distal.Intergenic","Intron","Exon",
                                                 "Upstream","Downstream"),8))


# Figure 5D - 698x580
pdf(file = "Figure 5D.pdf",width = 6,height = 6)
pheatmap(cor(t(scale(t(tf50allregionpctmat[apply(tf50allregionpctmat,1,max) > 1,])))),annotation_col = tf50allregionmeta[,c("Region","Tissue")],annotation_row = tf50allregionmeta[,c("Region","Tissue")],breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),annotation_colors = ann_cols,show_rownames = F,show_colnames = F,annotation_legend = F)
dev.off()
pdf(file = "Figure 5DLegend.pdf",width = 6,height = 6)
pheatmap(cor(t(scale(t(tf50allregionpctmat[apply(tf50allregionpctmat,1,max) > 1,])))),annotation_col = tf50allregionmeta,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),annotation_colors = ann_cols,show_rownames = F,show_colnames = F)
dev.off()

tf50subregionpctmat <- tf50allregionpctmat[,rownames(tf50allregionmeta[tf50allregionmeta$Region %in% c("Distal.Intergenic","Exon","Intron","Promoter.(<=1kb)","Promoter.(1-2kb)","Upstream","Downstream"),])]

tf50subregionpctmatzscore <- t(scale(t(tf50subregionpctmat[apply(tf50subregionpctmat,1,max) > 1,])))
tf50subregionavgzscore <- cbind(rowMeans(tf50subregionpctmatzscore[,colnames(tf50subregionpctmatzscore)[grep("promprox",colnames(tf50subregionpctmatzscore))]]),
                                rowMeans(tf50subregionpctmatzscore[,colnames(tf50subregionpctmatzscore)[grep("promfar",colnames(tf50subregionpctmatzscore))]]),
                                rowMeans(tf50subregionpctmatzscore[,colnames(tf50subregionpctmatzscore)[grep("dist",colnames(tf50subregionpctmatzscore))]]),
                                rowMeans(tf50subregionpctmatzscore[,colnames(tf50subregionpctmatzscore)[grep("ex",colnames(tf50subregionpctmatzscore))]]),
                                rowMeans(tf50subregionpctmatzscore[,colnames(tf50subregionpctmatzscore)[grep("int",colnames(tf50subregionpctmatzscore))]]),
                                rowMeans(tf50subregionpctmatzscore[,colnames(tf50subregionpctmatzscore)[grep("up",colnames(tf50subregionpctmatzscore))]]),
                                rowMeans(tf50subregionpctmatzscore[,colnames(tf50subregionpctmatzscore)[grep("down",colnames(tf50subregionpctmatzscore))]]))
colnames(tf50subregionavgzscore) <- c("Promoter.(<=1kb)",
                                      "Promoter.(1-2kb)",
                                      "Distal.Intergenic",
                                      "Exon",
                                      "Intron",
                                      "Upstream",
                                      "Downstream")
tf50subregionavgzscore[is.na(tf50subregionavgzscore)] <- 0
tf50subregionavgmeta <- data.frame(row.names = colnames(tf50subregionavgzscore),
                                   "Region" = colnames(tf50subregionavgzscore))

# Figure 5E - 500x580
pdf(file = "Figure 5E.pdf",width = 5,height = 6)
pheatmap(tf50subregionavgzscore[apply(abs(tf50subregionavgzscore),1,max) > 1.2,],show_rownames = T,labels_row = gsub("\\(.*","",rownames(tf50subregionavgzscore[apply(abs(tf50subregionavgzscore),1,max) > 1.2,])),angle_col = 0,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = T,annotation_col = tf50subregionavgmeta,annotation_colors = ann_cols,show_colnames = F)
dev.off()


####
# Supplemental Figure S9
#####

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

pdf(file = "Supplemental Figure S9A.pdf",width = 6,height = 6)
pheatmap(t(scale(t(tfatacsigpctmat[finalatacsigtfs,]))),cluster_rows = F,cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_row = gsub("\\(.*","",finalatacsigtfs),annotation_col = tissuemeta,annotation_colors = ann_cols,show_colnames = F)
dev.off()

pdf(file = "Supplemental Figure S9B.pdf",width = 6,height = 6)
pheatmap(t(scale(t(rnacontrolnorm[tfanno[finalatacsigtfs,"Ensembl"],controlcolumns[c(1:20,31:70)]]))),cluster_rows = F,cluster_cols = F,breaks = seq(-3,3,length.out = 101),color = colorpanel(101,"blue","white","red"),annotation_col = rnacontrolmeta[,c("Sex","Tissue")],annotation_colors = ann_cols,show_colnames = F,labels_row = enstosym[tfanno[finalatacsigtfs,"Ensembl"],"Symbol"],cellwidth = 4.0)
dev.off()

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

pdf(file = "Supplemental Figure S9C.pdf",width = 6.75,height = 6)
pheatmap(atacsigtfl2fcmat,cluster_rows = F,cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),labels_row = enstosym[rownames(atacsigtfl2fcmat),"Symbol"],annotation_col = atacsigtfl2fcmeta[,c("Group","Sex","Tissue")],annotation_colors = ann_cols,cellwidth = 6,show_colnames = F)
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

pdf(file = "Supplemental Figure S9D.pdf",width = 3.5,height = 5.5)
pheatmap(atactfexprenrichcor["Correlation",],cluster_rows = F,cluster_cols = F,display_numbers = T,number_color = "black",breaks = seq(-0.6,0.6,length.out = 101),color = colorpanel(101,"royalblue4","white","firebrick"),labels_row = colnames(atactfexprenrichcor),angle_col = 0,fontsize = 20,show_rownames = T,cellheight = 60)
dev.off()

heartatacsigexprenrichcordf <- data.frame("TF" = gsub("\\(.*","",finalatacsigtfs),
                                          "Enrichment Z Score" = tfatacsigpctmatz[finalatacsigtfs,"HEART"],
                                          "Expression Z Score" = atactfcontrolmeanmatz[,"HEART"])

lungatacsigexprenrichcordf <- data.frame("TF" = gsub("\\(.*","",finalatacsigtfs),
                                         "Enrichment Z Score" = tfatacsigpctmatz[finalatacsigtfs,"LUNG"],
                                         "Expression Z Score" = atactfcontrolmeanmatz[,"LUNG"])

pdf(file = "Supplemental Figure S9E.pdf",width = 6,height = 6)
ggplot(heartatacsigexprenrichcordf,aes(x=Enrichment.Z.Score,y=Expression.Z.Score,label=TF)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_text_repel(max.overlaps = 15)+geom_smooth(method=lm,  linetype="dashed",color="darkred", fill="blue") + theme_classic() + ggtitle(paste("HEART Correlation: ",toString(cor(tfatacsigpctmatz[finalatacsigtfs,"HEART"],atactfcontrolmeanmatz[,"HEART"])),sep = "")) + theme(text = element_text(size = 15))
dev.off()

pdf(file = "Supplemental Figure S9F.pdf",width = 6,height = 6)
ggplot(lungatacsigexprenrichcordf,aes(x=Enrichment.Z.Score,y=Expression.Z.Score,label=TF)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_text_repel(max.overlaps = 15)+geom_smooth(method=lm,  linetype="dashed",color="darkred", fill="blue") + theme_classic() + ggtitle(paste("LUNG Correlation: ",toString(cor(tfatacsigpctmatz[finalatacsigtfs,"LUNG"],atactfcontrolmeanmatz[,"LUNG"])),sep = "")) + theme(text = element_text(size = 15))
dev.off()

####
# Supplemental Figure S10
#####

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


# Supplemental Figure S10A
pdf(file = "Supplemental Figure S10A.pdf",width = 6,height = 6)
pheatmap(t(scale(t(tf50pctmat[finalsigtflist,]))),labels_row = gsub("\\(.*","",finalsigtflist),angle_col = 0,color = colorpanel(101,"blue","white","red"),breaks = seq(-2,2,length.out = 101),cluster_cols = F,cluster_rows = F,annotation_col = tissuemeta,annotation_colors = ann_colsheatmaptrim,show_colnames = F)
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
ann_control_cols <- list("Tissue" = tissue_cols[unique(rnacontrolmetacolumns$Tissue)],
                         "Sex" = sex_cols[unique(rnacontrolmetacolumns$Sex)])

pdf(file = "Supplemental Figure S10B.pdf",width = 6,height = 6)
pheatmap(t(scale(t(rnacontrolnorm[sigtfensembl,controlcolumns]))),cluster_rows = F,cluster_cols = F,labels_row = tfanno[finalsigtflist,"Gene.Name"],show_colnames = F,annotation_col = rnacontrolmetacolumns[,c("Sex","Tissue")],breaks = seq(-3,3,length.out = 101),color = colorpanel(101,"blue","white","red"),annotation_colors = ann_cols,cellwidth = 3.5)
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

pdf(file = "Supplemental Figure S10C.pdf",width = 7,height = 6)
pheatmap(rnasigtfl2fcmat,cluster_rows = F,angle_col = 315,labels_row = tfanno[finalsigtflist,"Gene.Name"],breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),cluster_cols = F,annotation_col = columnmetadf[,c("Group","Sex","Tissue")],show_colnames = F,annotation_colors = ann_cols,cellwidth = 5)
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

# Supplemental Figure S10D
pdf(file = "Supplemental Figure S10D.pdf",width = 3.5,height = 5.5)
pheatmap(tfexprenrichcor["Correlation",],cluster_rows = F,cluster_cols = F,display_numbers = T,number_color = "black",breaks = seq(0,0.6,length.out = 101),color = colorpanel(101,"white","firebrick"),labels_row = colnames(tfexprenrichcor),angle_col = 0,fontsize = 20,show_rownames = T,cellheight = 43)
dev.off()


liverexprenrichcordf <- data.frame(row.names = gsub("\\(.*","",finalsigtflist),
                                   "TF" = gsub("\\(.*","",finalsigtflist),
                                   "Enrichment Z Score" = tf50pctmatz[finalsigtflist,"LIVER"],
                                   "Expression Z Score" = ourcontrolmeanmatz[,"LIVER"])
lungexprenrichcordf <- data.frame(row.names = gsub("\\(.*","",finalsigtflist),
                                  "TF" = gsub("\\(.*","",finalsigtflist),
                                  "Enrichment Z Score" = tf50pctmatz[finalsigtflist,"LUNG"],
                                  "Expression Z Score" = ourcontrolmeanmatz[,"LUNG"])

# Supplemental Figure S10E
pdf(file = "Supplemental Figure S10E.pdf",width = 6,height = 6)
ggplot(liverexprenrichcordf,aes(x=Enrichment.Z.Score,y=Expression.Z.Score,label=TF)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_text_repel(max.overlaps = 15)+geom_smooth(method=lm,  linetype="dashed",color="darkred", fill="blue") + theme_classic() + ggtitle(paste("LIVER Correlation: ",toString(cor(tf50pctmatz[finalsigtflist,"LIVER"],ourcontrolmeanmatz[,"LIVER"])),sep = "")) + theme(text = element_text(size = 15))
dev.off()

# Supplemental Figure S10F
pdf(file = "Supplemental Figure S10F.pdf",width = 6,height = 6)
ggplot(lungexprenrichcordf,aes(x=Enrichment.Z.Score,y=Expression.Z.Score,label=TF)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_text_repel(max.overlaps = 15)+geom_smooth(method=lm,  linetype="dashed",color="darkred", fill="blue") + theme_classic() + ggtitle(paste("LUNG Correlation: ",toString(cor(tf50pctmatz[finalsigtflist,"LUNG"],ourcontrolmeanmatz[,"LUNG"])),sep = "")) + theme(text = element_text(size = 15))
dev.off()

####
# Supplemental Figure S11
#####

ourtissue <- "SKM-GN"
ourdf <- data.frame(row.names = atacsigtflist,
                    "DEGaP Enrichment" = tf50logpmat[atacsigtflist,ourtissue],
                    "DAR Enrichment" = tfatacsigpmat[atacsigtflist,ourtissue],
                    "TF" = gsub("\\(.*","",atacsigtflist))
pdf(file = "Supplemental Figure S11A.pdf",width = 6.5,height = 5.5)
ggplot(ourdf,aes(x=DEGaP.Enrichment,y=DAR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DEGaP.Enrichment,ourdf$DAR.Enrichment)),sep = ""))
dev.off()


ourtissue <- "HEART"
ourdf <- data.frame(row.names = atacsigtflist,
                    "DEGaP Enrichment" = tf50logpmat[atacsigtflist,ourtissue],
                    "DAR Enrichment" = tfatacsigpmat[atacsigtflist,ourtissue],
                    "TF" = gsub("\\(.*","",atacsigtflist))
pdf(file = "Supplemental Figure S11B.pdf",width = 6.5,height = 5.5)
ggplot(ourdf,aes(x=DEGaP.Enrichment,y=DAR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DEGaP.Enrichment,ourdf$DAR.Enrichment)),sep = ""))
dev.off()


ourtissue <- "KIDNEY"
ourdf <- data.frame(row.names = atacsigtflist,
                    "DEGaP Enrichment" = tf50logpmat[atacsigtflist,ourtissue],
                    "DAR Enrichment" = tfatacsigpmat[atacsigtflist,ourtissue],
                    "TF" = gsub("\\(.*","",atacsigtflist))
pdf(file = "Supplemental Figure S11C.pdf",width = 6.5,height = 5.5)
ggplot(ourdf,aes(x=DEGaP.Enrichment,y=DAR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DEGaP.Enrichment,ourdf$DAR.Enrichment)),sep = ""))
dev.off()


ourtissue <- "LIVER"
ourdf <- data.frame(row.names = atacsigtflist,
                    "DEGaP Enrichment" = tf50logpmat[atacsigtflist,ourtissue],
                    "DAR Enrichment" = tfatacsigpmat[atacsigtflist,ourtissue],
                    "TF" = gsub("\\(.*","",atacsigtflist))
pdf(file = "Supplemental Figure S11D.pdf",width = 6.5,height = 5.5)
ggplot(ourdf,aes(x=DEGaP.Enrichment,y=DAR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DEGaP.Enrichment,ourdf$DAR.Enrichment)),sep = ""))
dev.off()


ourtissue <- "LUNG"
ourdf <- data.frame(row.names = atacsigtflist,
                    "DEGaP Enrichment" = tf50logpmat[atacsigtflist,ourtissue],
                    "DAR Enrichment" = tfatacsigpmat[atacsigtflist,ourtissue],
                    "TF" = gsub("\\(.*","",atacsigtflist))
pdf(file = "Supplemental Figure S11E.pdf",width = 6.5,height = 5.5)
ggplot(ourdf,aes(x=DEGaP.Enrichment,y=DAR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DEGaP.Enrichment,ourdf$DAR.Enrichment)),sep = ""))
dev.off()


ourtissue <- "BAT"
ourdf <- data.frame(row.names = atacsigtflist,
                    "DEGaP Enrichment" = tf50logpmat[atacsigtflist,ourtissue],
                    "DAR Enrichment" = tfatacsigpmat[atacsigtflist,ourtissue],
                    "TF" = gsub("\\(.*","",atacsigtflist))
pdf(file = "Supplemental Figure S11F.pdf",width = 6.5,height = 5.5)
ggplot(ourdf,aes(x=DEGaP.Enrichment,y=DAR.Enrichment,labels=TF))+geom_point()+geom_text_repel(label = ourdf$TF) + theme_classic() + ggtitle(paste(ourtissue," Correlation: ",toString(cor(ourdf$DEGaP.Enrichment,ourdf$DAR.Enrichment)),sep = ""))
dev.off()

####
# Supplemental Figure S12
#####

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

activevsdegcorrmat <- matrix(0L,nrow = 3,ncol = 8)
rownames(activevsdegcorrmat) <- c("Active vs DEG","Active Prom vs DEG Prom","Up-Reg Prom vs Down-Reg Prom")
colnames(activevsdegcorrmat) <- c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC")
activevsdegcorrmat["Active vs DEG","SKM-GN"] <- cor(as.numeric(gsub("%","",gastro_stan_activetf[fintflist,"X..of.Target.Sequences.with.Motif"])),as.numeric(gsub("%","",gastro_stantf[fintflist,"X..of.Target.Sequences.with.Motif"])))
activevsdegcorrmat["Active vs DEG","HEART"] <- cor(as.numeric(gsub("%","",heart_activetf[fintflist,"X..of.Target.Sequences.with.Motif"])),as.numeric(gsub("%","",hearttf[fintflist,"X..of.Target.Sequences.with.Motif"])))
activevsdegcorrmat["Active vs DEG","HIPPOC"] <- cor(as.numeric(gsub("%","",hippo_activetf[fintflist,"X..of.Target.Sequences.with.Motif"])),as.numeric(gsub("%","",hippotf[fintflist,"X..of.Target.Sequences.with.Motif"])))
activevsdegcorrmat["Active vs DEG","KIDNEY"] <- cor(as.numeric(gsub("%","",kidney_activetf[fintflist,"X..of.Target.Sequences.with.Motif"])),as.numeric(gsub("%","",kidneytf[fintflist,"X..of.Target.Sequences.with.Motif"])))
activevsdegcorrmat["Active vs DEG","LIVER"] <- cor(as.numeric(gsub("%","",liver_activetf[fintflist,"X..of.Target.Sequences.with.Motif"])),as.numeric(gsub("%","",livertf[fintflist,"X..of.Target.Sequences.with.Motif"])))
activevsdegcorrmat["Active vs DEG","LUNG"] <- cor(as.numeric(gsub("%","",lung_activetf[fintflist,"X..of.Target.Sequences.with.Motif"])),as.numeric(gsub("%","",lungtf[fintflist,"X..of.Target.Sequences.with.Motif"])))
activevsdegcorrmat["Active vs DEG","BAT"] <- cor(as.numeric(gsub("%","",brown_activetf[fintflist,"X..of.Target.Sequences.with.Motif"])),as.numeric(gsub("%","",browntf[fintflist,"X..of.Target.Sequences.with.Motif"])))
activevsdegcorrmat["Active vs DEG","WAT-SC"] <- cor(as.numeric(gsub("%","",white_stan_activetf[fintflist,"X..of.Target.Sequences.with.Motif"])),as.numeric(gsub("%","",white_stantf[fintflist,"X..of.Target.Sequences.with.Motif"])))

activevsdegcorrmat["Active Prom vs DEG Prom","SKM-GN"] <- cor(tf50activeprompeaks[fintflist,c("SKM.GN.Enrichment.in.Active.Prom.Peaks")],100*gastro50motifdf[fintflist,"Enrichment.in.Sig.Prom.Peaks"])
activevsdegcorrmat["Active Prom vs DEG Prom","HEART"] <- cor(tf50activeprompeaks[fintflist,c("HEART.Enrichment.in.Active.Prom.Peaks")],100*as.numeric(gsub("%","",heart50motifdf[fintflist,"Enrichment.in.Sig.Prom.Peaks"])))
activevsdegcorrmat["Active Prom vs DEG Prom","HIPPOC"] <- cor(tf50activeprompeaks[fintflist,c("HIPPOC.Enrichment.in.Active.Prom.Peaks")],100*as.numeric(gsub("%","",hippo50motifdf[fintflist,"Enrichment.in.Sig.Prom.Peaks"])))
activevsdegcorrmat["Active Prom vs DEG Prom","KIDNEY"] <- cor(tf50activeprompeaks[fintflist,c("KIDNEY.Enrichment.in.Active.Prom.Peaks")],100*as.numeric(gsub("%","",kidney50motifdf[fintflist,"Enrichment.in.Sig.Prom.Peaks"])))
activevsdegcorrmat["Active Prom vs DEG Prom","LIVER"] <- cor(tf50activeprompeaks[fintflist,c("LIVER.Enrichment.in.Active.Prom.Peaks")],100*as.numeric(gsub("%","",liver50motifdf[fintflist,"Enrichment.in.Sig.Prom.Peaks"])))
activevsdegcorrmat["Active Prom vs DEG Prom","LUNG"] <- cor(tf50activeprompeaks[fintflist,c("LUNG.Enrichment.in.Active.Prom.Peaks")],100*as.numeric(gsub("%","",lung50motifdf[fintflist,"Enrichment.in.Sig.Prom.Peaks"])))
activevsdegcorrmat["Active Prom vs DEG Prom","BAT"] <- cor(tf50activeprompeaks[fintflist,c("BAT.Enrichment.in.Active.Prom.Peaks")],100*as.numeric(gsub("%","",brown50motifdf[fintflist,"Enrichment.in.Sig.Prom.Peaks"])))
activevsdegcorrmat["Active Prom vs DEG Prom","WAT-SC"] <- cor(tf50activeprompeaks[fintflist,c("WAT.SC.Enrichment.in.Active.Prom.Peaks")],100*as.numeric(gsub("%","",white50motifdf[fintflist,"Enrichment.in.Sig.Prom.Peaks"])))

activevsdegcorrmat["Up-Reg Prom vs Down-Reg Prom","SKM-GN"] <- cor(gastro50motifdf[,c("Enrichment.in.W8.UpReg.Sig.Prom.Peaks")],gastro50motifdf[,c("Enrichment.in.W8.DownReg.Sig.Prom.Peaks")])
activevsdegcorrmat["Up-Reg Prom vs Down-Reg Prom","HEART"] <- cor(heart50motifdf[,c("Enrichment.in.W8.UpReg.Sig.Prom.Peaks")],heart50motifdf[,c("Enrichment.in.W8.DownReg.Sig.Prom.Peaks")])
activevsdegcorrmat["Up-Reg Prom vs Down-Reg Prom","HIPPOC"] <- cor(hippo50motifdf[,c("Enrichment.in.W8.UpReg.Sig.Prom.Peaks")],hippo50motifdf[,c("Enrichment.in.W8.DownReg.Sig.Prom.Peaks")])
activevsdegcorrmat["Up-Reg Prom vs Down-Reg Prom","KIDNEY"] <- cor(kidney50motifdf[,c("Enrichment.in.W8.UpReg.Sig.Prom.Peaks")],kidney50motifdf[,c("Enrichment.in.W8.DownReg.Sig.Prom.Peaks")])
activevsdegcorrmat["Up-Reg Prom vs Down-Reg Prom","LIVER"] <- cor(liver50motifdf[,c("Enrichment.in.W8.UpReg.Sig.Prom.Peaks")],liver50motifdf[,c("Enrichment.in.W8.DownReg.Sig.Prom.Peaks")])
activevsdegcorrmat["Up-Reg Prom vs Down-Reg Prom","LUNG"] <- cor(lung50motifdf[,c("Enrichment.in.W8.UpReg.Sig.Prom.Peaks")],lung50motifdf[,c("Enrichment.in.W8.DownReg.Sig.Prom.Peaks")])
activevsdegcorrmat["Up-Reg Prom vs Down-Reg Prom","BAT"] <- cor(brown50motifdf[,c("Enrichment.in.W8.UpReg.Sig.Prom.Peaks")],brown50motifdf[,c("Enrichment.in.W8.DownReg.Sig.Prom.Peaks")])
activevsdegcorrmat["Up-Reg Prom vs Down-Reg Prom","WAT-SC"] <- cor(white50motifdf[,c("Enrichment.in.W8.UpReg.Sig.Prom.Peaks")],white50motifdf[,c("Enrichment.in.W8.DownReg.Sig.Prom.Peaks")])

activevsdegcorrdf <- data.frame("Correlation" = as.vector(t(activevsdegcorrmat)),
                                "Comparison" = c(rep("Active vs DEG",8),
                                                 rep("Active Prom vs DEG Prom",8),
                                                 rep("Up-Reg Prom vs Down-Reg Prom",8)),
                                "Tissue" = rep(c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC"),3))
activevsdegcorrdf$Comparison <- factor(activevsdegcorrdf$Comparison,levels = c("Active vs DEG",
                                                                               "Active Prom vs DEG Prom",
                                                                               "Up-Reg Prom vs Down-Reg Prom"))

####
# Supplemental Figure S12
#####

# Supplemental Figure S12A
pdf(file = "Supplemental Figure S12A.pdf",width = 7,height = 5.5)
pheatmap(t(scale(t(gastro50motifdf[apply(gastro50motifdf,1,max) > 0.01,order(motifdfmetadata$Region,motifdfmetadata$Training.Response)]))),labels_row = gsub("\\(.*","",rownames(gastro50motifdf))[apply(gastro50motifdf,1,max) > 0.01],annotation_col = motifdfmetadata[,c("Week","Training.Response","Region")],annotation_colors = ann_colsheatmaptrim,show_colnames = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,cluster_cols = F)
dev.off()
# Supplemental Figure S12B
pdf(file = "Supplemental Figure S12B.pdf",width = 7,height = 5.5)
pheatmap(t(scale(t(heart50motifdf[apply(heart50motifdf,1,max) > 0.01,order(motifdfmetadata$Region,motifdfmetadata$Training.Response)]))),labels_row = gsub("\\(.*","",rownames(heart50motifdf))[apply(heart50motifdf,1,max) > 0.01],annotation_col = motifdfmetadata[,c("Week","Training.Response","Region")],annotation_colors = ann_colsheatmaptrim,show_colnames = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,cluster_cols = F)
dev.off()
# Supplemental Figure S12C
pdf(file = "Supplemental Figure S12C.pdf",width = 7,height = 5.5)
pheatmap(t(scale(t(hippo50motifdf[apply(hippo50motifdf,1,max) > 0.01,order(motifdfmetadata$Region,motifdfmetadata$Training.Response)]))),labels_row = gsub("\\(.*","",rownames(hippo50motifdf))[apply(hippo50motifdf,1,max) > 0.01],annotation_col = motifdfmetadata[,c("Week","Training.Response","Region")],annotation_colors = ann_colsheatmaptrim,show_colnames = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,cluster_cols = F)
dev.off()
# Supplemental Figure S12D
pdf(file = "Supplemental Figure S12D.pdf",width = 7,height = 5.5)
pheatmap(t(scale(t(kidney50motifdf[apply(kidney50motifdf,1,max) > 0.01,order(motifdfmetadata$Region,motifdfmetadata$Training.Response)]))),labels_row = gsub("\\(.*","",rownames(kidney50motifdf))[apply(kidney50motifdf,1,max) > 0.01],annotation_col = motifdfmetadata[,c("Week","Training.Response","Region")],annotation_colors = ann_colsheatmaptrim,show_colnames = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,cluster_cols = F)
dev.off()
# Supplemental Figure S12E
pdf(file = "Supplemental Figure S12E.pdf",width = 7,height = 5.5)
pheatmap(t(scale(t(liver50motifdf[apply(liver50motifdf,1,max) > 0.01,order(motifdfmetadata$Region,motifdfmetadata$Training.Response)]))),labels_row = gsub("\\(.*","",rownames(liver50motifdf))[apply(liver50motifdf,1,max) > 0.01],annotation_col = motifdfmetadata[,c("Week","Training.Response","Region")],annotation_colors = ann_colsheatmaptrim,show_colnames = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,cluster_cols = F)
dev.off()
# Supplemental Figure S12F
pdf(file = "Supplemental Figure S12F.pdf",width = 7,height = 5.5)
pheatmap(t(scale(t(lung50motifdf[apply(lung50motifdf,1,max) > 0.01,order(motifdfmetadata$Region,motifdfmetadata$Training.Response)]))),labels_row = gsub("\\(.*","",rownames(lung50motifdf))[apply(lung50motifdf,1,max) > 0.01],annotation_col = motifdfmetadata[,c("Week","Training.Response","Region")],annotation_colors = ann_colsheatmaptrim,show_colnames = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,cluster_cols = F)
dev.off()
# Supplemental Figure S12G
pdf(file = "Supplemental Figure S12G.pdf",width = 7,height = 5.5)
pheatmap(t(scale(t(brown50motifdf[apply(brown50motifdf,1,max) > 0.01,order(motifdfmetadata$Region,motifdfmetadata$Training.Response)]))),labels_row = gsub("\\(.*","",rownames(brown50motifdf))[apply(brown50motifdf,1,max) > 0.01],annotation_col = motifdfmetadata[,c("Week","Training.Response","Region")],annotation_colors = ann_colsheatmaptrim,show_colnames = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,cluster_cols = F)
dev.off()
# Supplemental Figure S12H
pdf(file = "Supplemental Figure S12H.pdf",width = 7,height = 5.5)
pheatmap(t(scale(t(white50motifdf[apply(white50motifdf,1,max) > 0.01,order(motifdfmetadata$Region,motifdfmetadata$Training.Response)]))),labels_row = gsub("\\(.*","",rownames(white50motifdf))[apply(white50motifdf,1,max) > 0.01],annotation_col = motifdfmetadata[,c("Week","Training.Response","Region")],annotation_colors = ann_colsheatmaptrim,show_colnames = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,cluster_cols = F)
dev.off()


####
# Supplemental Figure S13
#####

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


# Fig S13a - 750x1500
pdf(file = "Supplemental Figure S13A.pdf",width = 7.5,height = 15)
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
         fontsize = 15)
dev.off()

# Fig S13b - 750x1500
pdf(file = "Supplemental Figure S13B.pdf",width = 7.5,height = 15)
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
         fontsize = 15)
dev.off()

# Fig S13c - 750x1500
pdf(file = "Supplemental Figure S13C.pdf",width = 7.5,height = 15)
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
         fontsize = 15)
dev.off()
# Fig S13d - 750x1500
pdf(file = "Supplemental Figure S13D.pdf",width = 7.5,height = 15)
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
         fontsize = 15)
dev.off()
# Fig S13e - 750x1500
pdf(file = "Supplemental Figure S13E.pdf",width = 7.5,height = 15)
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
         fontsize = 15)
dev.off()
# Fig S13f - 750x1500
pdf(file = "Supplemental Figure S13F.pdf",width = 7.5,height = 15)
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
         fontsize = 15)
dev.off()
# Fig S13g - 750x1500
pdf(file = "Supplemental Figure S13G.pdf",width = 7.5,height = 15)
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
         fontsize = 15)
dev.off()
# Fig S13h - 750x1500
pdf(file = "Supplemental Figure S13H.pdf",width = 7.5,height = 15)
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
         fontsize = 15)
dev.off()

####
# Figure 6
#####

####
# Figure 6A-D
#####

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

# Figure 6A
pdf(file = "Figure 6A.pdf",width = 5,height = 6)
pheatmap(tfmotif50prommeanzmatupreg[apply(abs(tfmotif50prommeanzmatupreg[,c(3,7,5,8,2,6,1,4)]),1,max) > 0.88,c(3,7,5,8,2,6,1,4)],show_rownames = T,cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = 0,labels_row = gsub("\\(.*","",rownames(tfmotif50prommeanzmatupreg[apply(abs(tfmotif50prommeanzmatupreg[,c(3,7,5,8,2,6,1,4)]),1,max) > 0.88,])),annotation_col = promzmatmetadf,annotation_colors = tissue_cols,show_colnames = F,cellwidth = 18)
dev.off()

# Figure 6B
pdf(file = "Figure 6B.pdf",width = 5,height = 6)
pheatmap(tfmotif50prommeanzmatdownreg[apply(abs(tfmotif50prommeanzmatdownreg[,c(11,15,10,16,9,12,13,14)]),1,max) > 0.88,c(11,15,10,16,9,12,13,14)],show_rownames = T,cluster_cols = F,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = 0,labels_row = gsub("\\(.*","",rownames(tfmotif50prommeanzmatdownreg[apply(abs(tfmotif50prommeanzmatdownreg[,c(11,15,10,16,9,12,13,14)]),1,max) > 0.88,])),annotation_col = promzmatmetadf,annotation_colors = tissue_cols,show_colnames = F,cellwidth = 18)
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

# Figure 6C
pdf(file = "Figure 6C.pdf",width = 5,height = 4)
pheatmap(prommeanupregtissuecomp50,display_numbers = T,angle_col = 0,color = colorpanel(101,"white","firebrick"),annotation_row = tissuemeta,annotation_col = tissuemeta,annotation_colors = ann_cols,show_rownames = F,show_colnames = F,fontsize = 15,number_color = "black",number_format = "%d",breaks = seq(0,30,length.out = 101),annotation_names_col = F,annotation_names_row = F,annotation_legend = F)
dev.off()

# Figure 6D
pdf(file = "Figure 6D.pdf",width = 5,height = 4)
pheatmap(prommeandownregtissuecomp50,display_numbers = T,angle_col = 0,color = colorpanel(101,"white","firebrick"),annotation_row = tissuemeta,annotation_col = tissuemeta,annotation_colors = ann_cols,show_rownames = F,show_colnames = F,fontsize = 15,number_color = "black",number_format = "%d",breaks = seq(0,30,length.out = 101),annotation_names_col = F,annotation_names_row = F,annotation_legend = F)
dev.off()

####
# Supplemental Figure S14
#####

phenodf <- data.frame(row.names = intersect(phenomeasuredata$pid,pass1bphenodata$pass1bf0001),"pid" = intersect(phenomeasuredata$pid,pass1bphenodata$pass1bf0001))
phenodf$wgt_gain_after_train <- 0
phenodf$pct_body_fat_change <- 0
phenodf$pct_body_lean_change <- 0
phenodf$pct_body_fluid_change <- 0
phenodf$lactate_change_dueto_train <- 0
phenodf$vo2_max_change <- 0
phenodf$sex <- ""
phenodf$group <- ""
for(i in 1:dim(phenodf)[1]){
  ourid <- rownames(phenodf)[i]
  phenodf[ourid,"wgt_gain_after_train"] <- mean(phenomeasuredata[phenomeasuredata$pid %in% ourid,"calculated_variables___wgt_gain_after_train"])
  phenodf[ourid,"pct_body_fat_change"] <- mean(phenomeasuredata[phenomeasuredata$pid %in% ourid,"calculated_variables___pct_body_fat_change"])
  phenodf[ourid,"pct_body_lean_change"] <- mean(phenomeasuredata[phenomeasuredata$pid %in% ourid,"calculated_variables___pct_body_lean_change"])
  phenodf[ourid,"pct_body_fluid_change"] <- mean(phenomeasuredata[phenomeasuredata$pid %in% ourid,"calculated_variables___pct_body_fluid_change"])
  phenodf[ourid,"lactate_change_dueto_train"] <- mean(phenomeasuredata[phenomeasuredata$pid %in% ourid,"calculated_variables___lactate_change_dueto_train"])
  phenodf[ourid,"vo2_max_change"] <- mean(phenomeasuredata[phenomeasuredata$pid %in% ourid,"calculated_variables___vo_2_max_change"])
  phenodf[ourid,"sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0027"][1]]
  phenodf[ourid,"group"] <- pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1]
}

trimphenodf <- phenodf[c(1:36),]
trimphenodf[trimphenodf$group %in% "Eight-week program Control Group","group"] <- "Control"
trimphenodf[trimphenodf$group %in% "Eight-week program Training Group","group"] <- "8w"
trimphenodf[trimphenodf$group %in% "Four-week program","group"] <- "4w"
trimphenodf$group <- factor(trimphenodf$group,levels = c("Control","4w","8w"))
gastrophenomeasuredata <- phenomeasuredata[intersect(gsub("X","",colnames(gastroatacnorm)),rownames(phenomeasuredata)),]
trimmedphenodf <- trimphenodf[trimphenodf$pid %in% gastrophenomeasuredata$pid,]

my_comparisons <- list( c("Control", "4w"), c("4w", "8w"), c("Control", "8w") )

# Supplemental Figure S14A
pdf(file = "Supplemental Figure S14A.pdf",height = 6,width = 7)
ggboxplot(trimmedphenodf,x = "group",y = "wgt_gain_after_train",facet.by = "sex",palette = "jco",add = "jitter",size = 1) + ggtitle("Body Weight",) + theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5)) + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "t.test",size = 8) + ylim(-12,65) + xlab("Experimental Group") + ylab("Change over Experiment")
dev.off()
# Supplemental Figure S14B
pdf(file = "Supplemental Figure S14B.pdf",height = 6,width = 7)
ggboxplot(trimmedphenodf,x = "group",y = "pct_body_fat_change",facet.by = "sex",palette = "jco",add = "jitter",size = 1) + ggtitle("Body Fat") + theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5)) + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "t.test",size = 8) + ylim(-6,9) + xlab("Experimental Group") + ylab("Change over Experiment")
dev.off()
# Supplemental Figure S14C
pdf(file = "Supplemental Figure S14C.pdf",height = 6,width = 7)
ggboxplot(trimmedphenodf,x = "group",y = "vo2_max_change",facet.by = "sex",palette = "jco",add = "jitter",size = 1) + ggtitle("VO2 Max") + theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5)) + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "t.test",size = 8) + ylim(-16,33) + xlab("Experimental Group") + ylab("Change over Experiment")
dev.off()
# Supplemental Figure S14D
pdf(file = "Supplemental Figure S14D.pdf",height = 6,width = 7)
ggboxplot(trimmedphenodf,x = "group",y = "pct_body_lean_change",facet.by = "sex",palette = "jco",add = "jitter",size = 1) + ggtitle("Body Lean") + theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5)) + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "t.test",size = 8) + ylim(-6,12) + xlab("Experimental Group") + ylab("Change over Experiment")
dev.off()
# Supplemental Figure S14E
pdf(file = "Supplemental Figure S14E.pdf",height = 6,width = 7)
ggboxplot(trimmedphenodf,x = "group",y = "lactate_change_dueto_train",facet.by = "sex",palette = "jco",add = "jitter",size = 1) + ggtitle("Lactate") + theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5)) + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "t.test",size = 8) + ylim(2.5,17) + xlab("Experimental Group") + ylab("Change over Experiment")
dev.off()
# Supplemental Figure S14F
pdf(file = "Supplemental Figure S14F.pdf",height = 6,width = 7)
ggboxplot(trimmedphenodf,x = "group",y = "pct_body_fluid_change",facet.by = "sex",palette = "jco",add = "jitter",size = 1) + ggtitle("Body Water") + theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5)) + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "t.test",size = 8) + ylim(-1.6,3) + xlab("Experimental Group") + ylab("Change over Experiment")
dev.off()

####
# Supplemental Figure S15
#####

# SKM-GN

gastrornaphenomeasuredata <- phenomeasuredata[intersect(gsub("X","",colnames(gastrornanorm)),rownames(phenomeasuredata)),]
gastrornaphenomeasuredata <- gastrornaphenomeasuredata[,c("calculated_variables___wgt_gain_after_train","calculated_variables___pct_body_fat_change","calculated_variables___pct_body_lean_change","calculated_variables___pct_body_fluid_change","calculated_variables___lactate_change_dueto_train","calculated_variables___vo_2_max_change")]
gastrornasignorm <- gastrornanorm[gastrornasig,]
gastrornameta <- data.frame(row.names = colnames(gastrornasignorm),
                            "label" = colnames(gastrornasignorm))


gastrophenodf <- data.frame(row.names = colnames(gastrornasignorm),"label" = colnames(gastrornasignorm))
gastrophenodf$wgt_gain_after_train <- 0
gastrophenodf$pct_body_fat_change <- 0
gastrophenodf$pct_body_lean_change <- 0
gastrophenodf$pct_body_fluid_change <- 0
gastrophenodf$lactate_change_dueto_train <- 0
gastrophenodf$vo2_max_change <- 0
gastrophenodf$sex <- ""
gastrophenodf$group <- ""
for(i in 1:dim(gastrophenodf)[1]){
  ourid <- rownames(gastrophenodf)[i]
  gastrophenodf[ourid,"wgt_gain_after_train"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___wgt_gain_after_train"])
  gastrophenodf[ourid,"pct_body_fat_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_fat_change"])
  gastrophenodf[ourid,"pct_body_lean_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_lean_change"])
  gastrophenodf[ourid,"pct_body_fluid_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_fluid_change"])
  gastrophenodf[ourid,"lactate_change_dueto_train"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___lactate_change_dueto_train"])
  gastrophenodf[ourid,"vo2_max_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___vo_2_max_change"])
  gastrophenodf[ourid,"sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"][1]]
  gastrophenodf[ourid,"group"] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1]
}
trimgastrophenodf <- gastrophenodf[gastrophenodf$group %in% c("Four-week program",
                                                              "Eight-week program Training Group",
                                                              "Eight-week program Control Group"),]
trimgastrophenodf[trimgastrophenodf$group %in% "Eight-week program Control Group","group"] <- "Control"
trimgastrophenodf[trimgastrophenodf$group %in% "Eight-week program Training Group","group"] <- "8w"
trimgastrophenodf[trimgastrophenodf$group %in% "Four-week program","group"] <- "4w"
trimgastrophenodf$group <- factor(trimgastrophenodf$group,levels = c("Control","4w","8w"))

gastrornasigvsphenocor <- matrix(0L,nrow = length(gastrornasig),ncol = dim(gastrornaphenomeasuredata)[2])
rownames(gastrornasigvsphenocor) <- gastrornasig
colnames(gastrornasigvsphenocor) <- colnames(gastrornaphenomeasuredata)
for(i in 1:length(gastrornasig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:dim(gastrornaphenomeasuredata)[2]){
    gastrornasigvsphenocor[i,j] <- cor(t(gastrornasignorm[i,rownames(trimgastrophenodf)]),gastrornaphenomeasuredata[rownames(trimgastrophenodf),j])
  }
}

# HEART

heartrnaphenomeasuredata <- phenomeasuredata[intersect(gsub("X","",colnames(heartrnanorm)),rownames(phenomeasuredata)),]
heartrnaphenomeasuredata <- heartrnaphenomeasuredata[,c("calculated_variables___wgt_gain_after_train","calculated_variables___pct_body_fat_change","calculated_variables___pct_body_lean_change","calculated_variables___pct_body_fluid_change","calculated_variables___lactate_change_dueto_train","calculated_variables___vo_2_max_change")]
heartrnasignorm <- heartrnanorm[heartrnasig,]
heartrnameta <- data.frame(row.names = colnames(heartrnasignorm),
                           "label" = colnames(heartrnasignorm))


heartphenodf <- data.frame(row.names = colnames(heartrnasignorm),"label" = colnames(heartrnasignorm))
heartphenodf$wgt_gain_after_train <- 0
heartphenodf$pct_body_fat_change <- 0
heartphenodf$pct_body_lean_change <- 0
heartphenodf$pct_body_fluid_change <- 0
heartphenodf$lactate_change_dueto_train <- 0
heartphenodf$vo2_max_change <- 0
heartphenodf$sex <- ""
heartphenodf$group <- ""
for(i in 1:dim(heartphenodf)[1]){
  ourid <- rownames(heartphenodf)[i]
  heartphenodf[ourid,"wgt_gain_after_train"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___wgt_gain_after_train"])
  heartphenodf[ourid,"pct_body_fat_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_fat_change"])
  heartphenodf[ourid,"pct_body_lean_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_lean_change"])
  heartphenodf[ourid,"pct_body_fluid_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_fluid_change"])
  heartphenodf[ourid,"lactate_change_dueto_train"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___lactate_change_dueto_train"])
  heartphenodf[ourid,"vo2_max_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___vo_2_max_change"])
  heartphenodf[ourid,"sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"][1]]
  heartphenodf[ourid,"group"] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1]
}
trimheartphenodf <- heartphenodf[heartphenodf$group %in% c("Four-week program",
                                                           "Eight-week program Training Group",
                                                           "Eight-week program Control Group"),]
trimheartphenodf[trimheartphenodf$group %in% "Eight-week program Control Group","group"] <- "Control"
trimheartphenodf[trimheartphenodf$group %in% "Eight-week program Training Group","group"] <- "8w"
trimheartphenodf[trimheartphenodf$group %in% "Four-week program","group"] <- "4w"
trimheartphenodf$group <- factor(trimheartphenodf$group,levels = c("Control","4w","8w"))

heartrnasigvsphenocor <- matrix(0L,nrow = length(heartrnasig),ncol = dim(heartrnaphenomeasuredata)[2])
rownames(heartrnasigvsphenocor) <- heartrnasig
colnames(heartrnasigvsphenocor) <- colnames(heartrnaphenomeasuredata)
for(i in 1:length(heartrnasig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:dim(heartrnaphenomeasuredata)[2]){
    heartrnasigvsphenocor[i,j] <- cor(t(heartrnasignorm[i,rownames(trimheartphenodf)]),heartrnaphenomeasuredata[rownames(trimheartphenodf),j])
  }
}

# HIPPOC

hippornaphenomeasuredata <- phenomeasuredata[intersect(gsub("X","",colnames(hippornanorm)),rownames(phenomeasuredata)),]
hippornaphenomeasuredata <- hippornaphenomeasuredata[,c("calculated_variables___wgt_gain_after_train","calculated_variables___pct_body_fat_change","calculated_variables___pct_body_lean_change","calculated_variables___pct_body_fluid_change","calculated_variables___lactate_change_dueto_train","calculated_variables___vo_2_max_change")]
hippornasignorm <- hippornanorm[hippornasig,]
hippornameta <- data.frame(row.names = colnames(hippornasignorm),
                           "label" = colnames(hippornasignorm))


hippophenodf <- data.frame(row.names = colnames(hippornasignorm),"label" = colnames(hippornasignorm))
hippophenodf$wgt_gain_after_train <- 0
hippophenodf$pct_body_fat_change <- 0
hippophenodf$pct_body_lean_change <- 0
hippophenodf$pct_body_fluid_change <- 0
hippophenodf$lactate_change_dueto_train <- 0
hippophenodf$vo2_max_change <- 0
hippophenodf$sex <- ""
hippophenodf$group <- ""
for(i in 1:dim(hippophenodf)[1]){
  ourid <- rownames(hippophenodf)[i]
  hippophenodf[ourid,"wgt_gain_after_train"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___wgt_gain_after_train"])
  hippophenodf[ourid,"pct_body_fat_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_fat_change"])
  hippophenodf[ourid,"pct_body_lean_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_lean_change"])
  hippophenodf[ourid,"pct_body_fluid_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_fluid_change"])
  hippophenodf[ourid,"lactate_change_dueto_train"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___lactate_change_dueto_train"])
  hippophenodf[ourid,"vo2_max_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___vo_2_max_change"])
  hippophenodf[ourid,"sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"][1]]
  hippophenodf[ourid,"group"] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1]
}
trimhippophenodf <- hippophenodf[hippophenodf$group %in% c("Four-week program",
                                                           "Eight-week program Training Group",
                                                           "Eight-week program Control Group"),]
trimhippophenodf[trimhippophenodf$group %in% "Eight-week program Control Group","group"] <- "Control"
trimhippophenodf[trimhippophenodf$group %in% "Eight-week program Training Group","group"] <- "8w"
trimhippophenodf[trimhippophenodf$group %in% "Four-week program","group"] <- "4w"
trimhippophenodf$group <- factor(trimhippophenodf$group,levels = c("Control","4w","8w"))

hippornasigvsphenocor <- matrix(0L,nrow = length(hippornasig),ncol = dim(hippornaphenomeasuredata)[2])
rownames(hippornasigvsphenocor) <- hippornasig
colnames(hippornasigvsphenocor) <- colnames(hippornaphenomeasuredata)
for(i in 1:length(hippornasig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:dim(hippornaphenomeasuredata)[2]){
    hippornasigvsphenocor[i,j] <- cor(t(hippornasignorm[i,rownames(trimhippophenodf)]),hippornaphenomeasuredata[rownames(trimhippophenodf),j])
  }
}

# KIDNEY

kidneyrnaphenomeasuredata <- phenomeasuredata[intersect(gsub("X","",colnames(kidneyrnanorm)),rownames(phenomeasuredata)),]
kidneyrnaphenomeasuredata <- kidneyrnaphenomeasuredata[,c("calculated_variables___wgt_gain_after_train","calculated_variables___pct_body_fat_change","calculated_variables___pct_body_lean_change","calculated_variables___pct_body_fluid_change","calculated_variables___lactate_change_dueto_train","calculated_variables___vo_2_max_change")]
kidneyrnasignorm <- kidneyrnanorm[kidneyrnasig,]
kidneyrnameta <- data.frame(row.names = colnames(kidneyrnasignorm),
                            "label" = colnames(kidneyrnasignorm))


kidneyphenodf <- data.frame(row.names = colnames(kidneyrnasignorm),"label" = colnames(kidneyrnasignorm))
kidneyphenodf$wgt_gain_after_train <- 0
kidneyphenodf$pct_body_fat_change <- 0
kidneyphenodf$pct_body_lean_change <- 0
kidneyphenodf$pct_body_fluid_change <- 0
kidneyphenodf$lactate_change_dueto_train <- 0
kidneyphenodf$vo2_max_change <- 0
kidneyphenodf$sex <- ""
kidneyphenodf$group <- ""
for(i in 1:dim(kidneyphenodf)[1]){
  ourid <- rownames(kidneyphenodf)[i]
  kidneyphenodf[ourid,"wgt_gain_after_train"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___wgt_gain_after_train"])
  kidneyphenodf[ourid,"pct_body_fat_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_fat_change"])
  kidneyphenodf[ourid,"pct_body_lean_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_lean_change"])
  kidneyphenodf[ourid,"pct_body_fluid_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_fluid_change"])
  kidneyphenodf[ourid,"lactate_change_dueto_train"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___lactate_change_dueto_train"])
  kidneyphenodf[ourid,"vo2_max_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___vo_2_max_change"])
  kidneyphenodf[ourid,"sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"][1]]
  kidneyphenodf[ourid,"group"] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1]
}
trimkidneyphenodf <- kidneyphenodf[kidneyphenodf$group %in% c("Four-week program",
                                                              "Eight-week program Training Group",
                                                              "Eight-week program Control Group"),]
trimkidneyphenodf[trimkidneyphenodf$group %in% "Eight-week program Control Group","group"] <- "Control"
trimkidneyphenodf[trimkidneyphenodf$group %in% "Eight-week program Training Group","group"] <- "8w"
trimkidneyphenodf[trimkidneyphenodf$group %in% "Four-week program","group"] <- "4w"
trimkidneyphenodf$group <- factor(trimkidneyphenodf$group,levels = c("Control","4w","8w"))

kidneyrnasigvsphenocor <- matrix(0L,nrow = length(kidneyrnasig),ncol = dim(kidneyrnaphenomeasuredata)[2])
rownames(kidneyrnasigvsphenocor) <- kidneyrnasig
colnames(kidneyrnasigvsphenocor) <- colnames(kidneyrnaphenomeasuredata)
for(i in 1:length(kidneyrnasig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:dim(kidneyrnaphenomeasuredata)[2]){
    kidneyrnasigvsphenocor[i,j] <- cor(t(kidneyrnasignorm[i,rownames(trimkidneyphenodf)]),kidneyrnaphenomeasuredata[rownames(trimkidneyphenodf),j])
  }
}

# LIVER

liverrnaphenomeasuredata <- phenomeasuredata[intersect(gsub("X","",colnames(liverrnanorm)),rownames(phenomeasuredata)),]
liverrnaphenomeasuredata <- liverrnaphenomeasuredata[,c("calculated_variables___wgt_gain_after_train","calculated_variables___pct_body_fat_change","calculated_variables___pct_body_lean_change","calculated_variables___pct_body_fluid_change","calculated_variables___lactate_change_dueto_train","calculated_variables___vo_2_max_change")]
liverrnasignorm <- liverrnanorm[liverrnasig,]
liverrnameta <- data.frame(row.names = colnames(liverrnasignorm),
                           "label" = colnames(liverrnasignorm))


liverphenodf <- data.frame(row.names = colnames(liverrnasignorm),"label" = colnames(liverrnasignorm))
liverphenodf$wgt_gain_after_train <- 0
liverphenodf$pct_body_fat_change <- 0
liverphenodf$pct_body_lean_change <- 0
liverphenodf$pct_body_fluid_change <- 0
liverphenodf$lactate_change_dueto_train <- 0
liverphenodf$vo2_max_change <- 0
liverphenodf$sex <- ""
liverphenodf$group <- ""
for(i in 1:dim(liverphenodf)[1]){
  ourid <- rownames(liverphenodf)[i]
  liverphenodf[ourid,"wgt_gain_after_train"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___wgt_gain_after_train"])
  liverphenodf[ourid,"pct_body_fat_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_fat_change"])
  liverphenodf[ourid,"pct_body_lean_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_lean_change"])
  liverphenodf[ourid,"pct_body_fluid_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_fluid_change"])
  liverphenodf[ourid,"lactate_change_dueto_train"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___lactate_change_dueto_train"])
  liverphenodf[ourid,"vo2_max_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___vo_2_max_change"])
  liverphenodf[ourid,"sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"][1]]
  liverphenodf[ourid,"group"] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1]
}
trimliverphenodf <- liverphenodf[liverphenodf$group %in% c("Four-week program",
                                                           "Eight-week program Training Group",
                                                           "Eight-week program Control Group"),]
trimliverphenodf[trimliverphenodf$group %in% "Eight-week program Control Group","group"] <- "Control"
trimliverphenodf[trimliverphenodf$group %in% "Eight-week program Training Group","group"] <- "8w"
trimliverphenodf[trimliverphenodf$group %in% "Four-week program","group"] <- "4w"
trimliverphenodf$group <- factor(trimliverphenodf$group,levels = c("Control","4w","8w"))

liverrnasigvsphenocor <- matrix(0L,nrow = length(liverrnasig),ncol = dim(liverrnaphenomeasuredata)[2])
rownames(liverrnasigvsphenocor) <- liverrnasig
colnames(liverrnasigvsphenocor) <- colnames(liverrnaphenomeasuredata)
for(i in 1:length(liverrnasig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:dim(liverrnaphenomeasuredata)[2]){
    liverrnasigvsphenocor[i,j] <- cor(t(liverrnasignorm[i,rownames(trimliverphenodf)]),liverrnaphenomeasuredata[rownames(trimliverphenodf),j])
  }
}

# LUNG

lungrnaphenomeasuredata <- phenomeasuredata[intersect(gsub("X","",colnames(lungrnanorm)),rownames(phenomeasuredata)),]
lungrnaphenomeasuredata <- lungrnaphenomeasuredata[,c("calculated_variables___wgt_gain_after_train","calculated_variables___pct_body_fat_change","calculated_variables___pct_body_lean_change","calculated_variables___pct_body_fluid_change","calculated_variables___lactate_change_dueto_train","calculated_variables___vo_2_max_change")]
lungrnasignorm <- lungrnanorm[lungrnasig,]
lungrnameta <- data.frame(row.names = colnames(lungrnasignorm),
                          "label" = colnames(lungrnasignorm))


lungphenodf <- data.frame(row.names = colnames(lungrnasignorm),"label" = colnames(lungrnasignorm))
lungphenodf$wgt_gain_after_train <- 0
lungphenodf$pct_body_fat_change <- 0
lungphenodf$pct_body_lean_change <- 0
lungphenodf$pct_body_fluid_change <- 0
lungphenodf$lactate_change_dueto_train <- 0
lungphenodf$vo2_max_change <- 0
lungphenodf$sex <- ""
lungphenodf$group <- ""
for(i in 1:dim(lungphenodf)[1]){
  ourid <- rownames(lungphenodf)[i]
  lungphenodf[ourid,"wgt_gain_after_train"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___wgt_gain_after_train"])
  lungphenodf[ourid,"pct_body_fat_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_fat_change"])
  lungphenodf[ourid,"pct_body_lean_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_lean_change"])
  lungphenodf[ourid,"pct_body_fluid_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_fluid_change"])
  lungphenodf[ourid,"lactate_change_dueto_train"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___lactate_change_dueto_train"])
  lungphenodf[ourid,"vo2_max_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___vo_2_max_change"])
  lungphenodf[ourid,"sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"][1]]
  lungphenodf[ourid,"group"] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1]
}
trimlungphenodf <- lungphenodf[lungphenodf$group %in% c("Four-week program",
                                                        "Eight-week program Training Group",
                                                        "Eight-week program Control Group"),]
trimlungphenodf[trimlungphenodf$group %in% "Eight-week program Control Group","group"] <- "Control"
trimlungphenodf[trimlungphenodf$group %in% "Eight-week program Training Group","group"] <- "8w"
trimlungphenodf[trimlungphenodf$group %in% "Four-week program","group"] <- "4w"
trimlungphenodf$group <- factor(trimlungphenodf$group,levels = c("Control","4w","8w"))

lungrnasigvsphenocor <- matrix(0L,nrow = length(lungrnasig),ncol = dim(lungrnaphenomeasuredata)[2])
rownames(lungrnasigvsphenocor) <- lungrnasig
colnames(lungrnasigvsphenocor) <- colnames(lungrnaphenomeasuredata)
for(i in 1:length(lungrnasig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:dim(lungrnaphenomeasuredata)[2]){
    lungrnasigvsphenocor[i,j] <- cor(t(lungrnasignorm[i,rownames(trimlungphenodf)]),lungrnaphenomeasuredata[rownames(trimlungphenodf),j])
  }
}


# BAT

brownrnaphenomeasuredata <- phenomeasuredata[intersect(gsub("X","",colnames(brownrnanorm)),rownames(phenomeasuredata)),]
brownrnaphenomeasuredata <- brownrnaphenomeasuredata[,c("calculated_variables___wgt_gain_after_train","calculated_variables___pct_body_fat_change","calculated_variables___pct_body_lean_change","calculated_variables___pct_body_fluid_change","calculated_variables___lactate_change_dueto_train","calculated_variables___vo_2_max_change")]
brownrnasignorm <- brownrnanorm[brownrnasig,]
brownrnameta <- data.frame(row.names = colnames(brownrnasignorm),
                           "label" = colnames(brownrnasignorm))


brownphenodf <- data.frame(row.names = colnames(brownrnasignorm),"label" = colnames(brownrnasignorm))
brownphenodf$wgt_gain_after_train <- 0
brownphenodf$pct_body_fat_change <- 0
brownphenodf$pct_body_lean_change <- 0
brownphenodf$pct_body_fluid_change <- 0
brownphenodf$lactate_change_dueto_train <- 0
brownphenodf$vo2_max_change <- 0
brownphenodf$sex <- ""
brownphenodf$group <- ""
for(i in 1:dim(brownphenodf)[1]){
  ourid <- rownames(brownphenodf)[i]
  brownphenodf[ourid,"wgt_gain_after_train"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___wgt_gain_after_train"])
  brownphenodf[ourid,"pct_body_fat_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_fat_change"])
  brownphenodf[ourid,"pct_body_lean_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_lean_change"])
  brownphenodf[ourid,"pct_body_fluid_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_fluid_change"])
  brownphenodf[ourid,"lactate_change_dueto_train"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___lactate_change_dueto_train"])
  brownphenodf[ourid,"vo2_max_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___vo_2_max_change"])
  brownphenodf[ourid,"sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"][1]]
  brownphenodf[ourid,"group"] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1]
}
trimbrownphenodf <- brownphenodf[brownphenodf$group %in% c("Four-week program",
                                                           "Eight-week program Training Group",
                                                           "Eight-week program Control Group"),]
trimbrownphenodf[trimbrownphenodf$group %in% "Eight-week program Control Group","group"] <- "Control"
trimbrownphenodf[trimbrownphenodf$group %in% "Eight-week program Training Group","group"] <- "8w"
trimbrownphenodf[trimbrownphenodf$group %in% "Four-week program","group"] <- "4w"
trimbrownphenodf$group <- factor(trimbrownphenodf$group,levels = c("Control","4w","8w"))

brownrnasigvsphenocor <- matrix(0L,nrow = length(brownrnasig),ncol = dim(brownrnaphenomeasuredata)[2])
rownames(brownrnasigvsphenocor) <- brownrnasig
colnames(brownrnasigvsphenocor) <- colnames(brownrnaphenomeasuredata)
for(i in 1:length(brownrnasig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:dim(brownrnaphenomeasuredata)[2]){
    brownrnasigvsphenocor[i,j] <- cor(t(brownrnasignorm[i,rownames(trimbrownphenodf)]),brownrnaphenomeasuredata[rownames(trimbrownphenodf),j])
  }
}

# WAT-SC

whiternaphenomeasuredata <- phenomeasuredata[intersect(gsub("X","",colnames(whiternanorm)),rownames(phenomeasuredata)),]
whiternaphenomeasuredata <- whiternaphenomeasuredata[,c("calculated_variables___wgt_gain_after_train","calculated_variables___pct_body_fat_change","calculated_variables___pct_body_lean_change","calculated_variables___pct_body_fluid_change","calculated_variables___lactate_change_dueto_train","calculated_variables___vo_2_max_change")]
whiternasignorm <- whiternanorm[whiternasig,]
whiternameta <- data.frame(row.names = colnames(whiternasignorm),
                           "label" = colnames(whiternasignorm))


whitephenodf <- data.frame(row.names = colnames(whiternasignorm),"label" = colnames(whiternasignorm))
whitephenodf$wgt_gain_after_train <- 0
whitephenodf$pct_body_fat_change <- 0
whitephenodf$pct_body_lean_change <- 0
whitephenodf$pct_body_fluid_change <- 0
whitephenodf$lactate_change_dueto_train <- 0
whitephenodf$vo2_max_change <- 0
whitephenodf$sex <- ""
whitephenodf$group <- ""
for(i in 1:dim(whitephenodf)[1]){
  ourid <- rownames(whitephenodf)[i]
  whitephenodf[ourid,"wgt_gain_after_train"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___wgt_gain_after_train"])
  whitephenodf[ourid,"pct_body_fat_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_fat_change"])
  whitephenodf[ourid,"pct_body_lean_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_lean_change"])
  whitephenodf[ourid,"pct_body_fluid_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_fluid_change"])
  whitephenodf[ourid,"lactate_change_dueto_train"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___lactate_change_dueto_train"])
  whitephenodf[ourid,"vo2_max_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___vo_2_max_change"])
  whitephenodf[ourid,"sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"][1]]
  whitephenodf[ourid,"group"] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1]
}
trimwhitephenodf <- whitephenodf[whitephenodf$group %in% c("Four-week program",
                                                           "Eight-week program Training Group",
                                                           "Eight-week program Control Group"),]
trimwhitephenodf[trimwhitephenodf$group %in% "Eight-week program Control Group","group"] <- "Control"
trimwhitephenodf[trimwhitephenodf$group %in% "Eight-week program Training Group","group"] <- "8w"
trimwhitephenodf[trimwhitephenodf$group %in% "Four-week program","group"] <- "4w"
trimwhitephenodf$group <- factor(trimwhitephenodf$group,levels = c("Control","4w","8w"))

whiternasigvsphenocor <- matrix(0L,nrow = length(whiternasig),ncol = dim(whiternaphenomeasuredata)[2])
rownames(whiternasigvsphenocor) <- whiternasig
colnames(whiternasigvsphenocor) <- colnames(whiternaphenomeasuredata)
for(i in 1:length(whiternasig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:dim(whiternaphenomeasuredata)[2]){
    whiternasigvsphenocor[i,j] <- cor(t(whiternasignorm[i,rownames(trimwhitephenodf)]),whiternaphenomeasuredata[rownames(trimwhitephenodf),j])
  }
}

colnames(gastrornasigvsphenocor) <- gsub("calculated_variables___","",colnames(gastrornasigvsphenocor))
colnames(heartrnasigvsphenocor) <- gsub("calculated_variables___","",colnames(heartrnasigvsphenocor))
colnames(hippornasigvsphenocor) <- gsub("calculated_variables___","",colnames(hippornasigvsphenocor))
colnames(kidneyrnasigvsphenocor) <- gsub("calculated_variables___","",colnames(kidneyrnasigvsphenocor))
colnames(liverrnasigvsphenocor) <- gsub("calculated_variables___","",colnames(liverrnasigvsphenocor))
colnames(lungrnasigvsphenocor) <- gsub("calculated_variables___","",colnames(lungrnasigvsphenocor))
colnames(brownrnasigvsphenocor) <- gsub("calculated_variables___","",colnames(brownrnasigvsphenocor))
colnames(whiternasigvsphenocor) <- gsub("calculated_variables___","",colnames(whiternasigvsphenocor))

sigvsphenocormeta <- data.frame(row.names = colnames(gastrornasigvsphenocor),"Measure" = colnames(gastrornasigvsphenocor))
sigvsphenocormeta$Measure <- c("Body Weight","Body Fat","Body Lean","Body Water","Lactate","VO2 Max")

# Supplemental Figure S15A
pdf(file = "Supplemental Figure S15A.pdf",width = 3,height = 6)
pheatmap(gastrornasigvsphenocor,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,annotation_col = sigvsphenocormeta,show_colnames = F,annotation_colors = ann_cols_sigphenocor,annotation_names_col = F,fontsize = 15,annotation_legend = F)
dev.off()
# Supplemental Figure S15B
pdf(file = "Supplemental Figure S15B.pdf",width = 3,height = 6)
pheatmap(heartrnasigvsphenocor,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,annotation_col = sigvsphenocormeta,show_colnames = F,annotation_colors = ann_cols_sigphenocor,annotation_names_col = F,fontsize = 15,annotation_legend = F)
dev.off()
# Supplemental Figure S15C
pdf(file = "Supplemental Figure S15C.pdf",width = 3,height = 6)
pheatmap(hippornasigvsphenocor,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,annotation_col = sigvsphenocormeta,show_colnames = F,annotation_colors = ann_cols_sigphenocor,annotation_names_col = F,fontsize = 15,annotation_legend = F)
dev.off()
# Supplemental Figure S15D
pdf(file = "Supplemental Figure S15D.pdf",width = 3,height = 6)
pheatmap(kidneyrnasigvsphenocor,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,annotation_col = sigvsphenocormeta,show_colnames = F,annotation_colors = ann_cols_sigphenocor,annotation_names_col = F,fontsize = 15,annotation_legend = F)
dev.off()
# Supplemental Figure S15E
pdf(file = "Supplemental Figure S15E.pdf",width = 3,height = 6)
pheatmap(liverrnasigvsphenocor,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,annotation_col = sigvsphenocormeta,show_colnames = F,annotation_colors = ann_cols_sigphenocor,annotation_names_col = F,fontsize = 15,annotation_legend = F)
dev.off()
# Supplemental Figure S15F
pdf(file = "Supplemental Figure S15F.pdf",width = 3,height = 6)
pheatmap(lungrnasigvsphenocor,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,annotation_col = sigvsphenocormeta,show_colnames = F,annotation_colors = ann_cols_sigphenocor,annotation_names_col = F,fontsize = 15,annotation_legend = F)
dev.off()
# Supplemental Figure S15G
pdf(file = "Supplemental Figure S15G.pdf",width = 3,height = 6)
pheatmap(brownrnasigvsphenocor,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,annotation_col = sigvsphenocormeta,show_colnames = F,annotation_colors = ann_cols_sigphenocor,annotation_names_col = F,fontsize = 15,annotation_legend = F)
dev.off()
# Supplemental Figure S15H
pdf(file = "Supplemental Figure S15H.pdf",width = 3,height = 6)
pheatmap(whiternasigvsphenocor,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,annotation_col = sigvsphenocormeta,show_colnames = F,annotation_colors = ann_cols_sigphenocor,annotation_names_col = F,fontsize = 15,annotation_legend = F)
dev.off()


####
# Figure 7
#####

phenomatrix <- cbind(trimphenodf$wgt_gain_after_train,trimphenodf$pct_body_fat_change,
                     trimphenodf$pct_body_lean_change,trimphenodf$pct_body_fluid_change,
                     trimphenodf$lactate_change_dueto_train,trimphenodf$vo2_max_change)
rownames(phenomatrix) <- rownames(trimphenodf)
colnames(phenomatrix) <- c("Body Weight","Body Fat","Body Lean","Body Water","Lactate","VO2 Max")

# Figure 7A
pdf(file = "Figure 7A.pdf",height = 5,width = 7)
pheatmap(cor(phenomatrix),angle_col = "315",breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),display_numbers = T,number_color = "black",fontsize = 20,show_colnames = T,show_rownames = T)
dev.off()

# SKM-GN

gastro_wgtgainpos <- rownames(gastrornasigvsphenocor[gastrornasigvsphenocor[,"wgt_gain_after_train"] > 0.5,])
gastro_wgtgainneg <- rownames(gastrornasigvsphenocor[gastrornasigvsphenocor[,"wgt_gain_after_train"] < -0.5,])
gastro_bodyfatpos <- rownames(gastrornasigvsphenocor[gastrornasigvsphenocor[,"pct_body_fat_change"] > 0.5,])
gastro_bodyfatneg <- rownames(gastrornasigvsphenocor[gastrornasigvsphenocor[,"pct_body_fat_change"] < -0.5,])
gastro_bodyleanpos <- rownames(gastrornasigvsphenocor[gastrornasigvsphenocor[,"pct_body_lean_change"] > 0.5,])
gastro_bodyleanneg <- rownames(gastrornasigvsphenocor[gastrornasigvsphenocor[,"pct_body_lean_change"] < -0.5,])
gastro_bodyfluidpos <- rownames(gastrornasigvsphenocor[gastrornasigvsphenocor[,"pct_body_fluid_change"] > 0.5,])
gastro_bodyfluidneg <- rownames(gastrornasigvsphenocor[gastrornasigvsphenocor[,"pct_body_fluid_change"] < -0.5,])
gastro_lactatepos <- rownames(gastrornasigvsphenocor[gastrornasigvsphenocor[,"lactate_change_dueto_train"] > 0.5,])
gastro_lactateneg <- rownames(gastrornasigvsphenocor[gastrornasigvsphenocor[,"lactate_change_dueto_train"] < -0.5,])
gastro_vo2maxpos <- rownames(gastrornasigvsphenocor[gastrornasigvsphenocor[,"vo_2_max_change"] > 0.499,])
gastro_vo2maxneg <- rownames(gastrornasigvsphenocor[gastrornasigvsphenocor[,"vo_2_max_change"] < -0.5,])

# HEART

heart_wgtgainpos <- rownames(heartrnasigvsphenocor[heartrnasigvsphenocor[,"wgt_gain_after_train"] > 0.5,])
heart_wgtgainneg <- rownames(heartrnasigvsphenocor[heartrnasigvsphenocor[,"wgt_gain_after_train"] < -0.5,])
heart_bodyfatpos <- rownames(heartrnasigvsphenocor[heartrnasigvsphenocor[,"pct_body_fat_change"] > 0.5,])
heart_bodyfatneg <- rownames(heartrnasigvsphenocor[heartrnasigvsphenocor[,"pct_body_fat_change"] < -0.5,])
heart_bodyleanpos <- rownames(heartrnasigvsphenocor[heartrnasigvsphenocor[,"pct_body_lean_change"] > 0.5,])
heart_bodyleanneg <- rownames(heartrnasigvsphenocor[heartrnasigvsphenocor[,"pct_body_lean_change"] < -0.5,])
heart_bodyfluidpos <- rownames(heartrnasigvsphenocor[heartrnasigvsphenocor[,"pct_body_fluid_change"] > 0.5,])
heart_bodyfluidneg <- rownames(heartrnasigvsphenocor[heartrnasigvsphenocor[,"pct_body_fluid_change"] < -0.5,])
heart_lactatepos <- rownames(heartrnasigvsphenocor[heartrnasigvsphenocor[,"lactate_change_dueto_train"] > 0.5,])
heart_lactateneg <- rownames(heartrnasigvsphenocor[heartrnasigvsphenocor[,"lactate_change_dueto_train"] < -0.5,])
heart_vo2maxpos <- rownames(heartrnasigvsphenocor[heartrnasigvsphenocor[,"vo_2_max_change"] > 0.5,])
heart_vo2maxneg <- rownames(heartrnasigvsphenocor[heartrnasigvsphenocor[,"vo_2_max_change"] < -0.5,])

# HIPPOC

hippo_wgtgainpos <- rownames(hippornasigvsphenocor[hippornasigvsphenocor[,"wgt_gain_after_train"] > 0.5,])
hippo_wgtgainneg <- rownames(hippornasigvsphenocor[hippornasigvsphenocor[,"wgt_gain_after_train"] < -0.5,])
hippo_bodyfatpos <- rownames(hippornasigvsphenocor[hippornasigvsphenocor[,"pct_body_fat_change"] > 0.5,])
hippo_bodyfatneg <- rownames(hippornasigvsphenocor[hippornasigvsphenocor[,"pct_body_fat_change"] < -0.5,])
hippo_bodyleanpos <- rownames(hippornasigvsphenocor[hippornasigvsphenocor[,"pct_body_lean_change"] > 0.5,])
hippo_bodyleanneg <- rownames(hippornasigvsphenocor[hippornasigvsphenocor[,"pct_body_lean_change"] < -0.5,])
hippo_bodyfluidpos <- rownames(hippornasigvsphenocor[hippornasigvsphenocor[,"pct_body_fluid_change"] > 0.5,])
hippo_bodyfluidneg <- rownames(hippornasigvsphenocor[hippornasigvsphenocor[,"pct_body_fluid_change"] < -0.5,])
hippo_lactatepos <- rownames(hippornasigvsphenocor[hippornasigvsphenocor[,"lactate_change_dueto_train"] > 0.5,])
hippo_lactateneg <- rownames(hippornasigvsphenocor[hippornasigvsphenocor[,"lactate_change_dueto_train"] < -0.5,])
hippo_vo2maxpos <- rownames(hippornasigvsphenocor[hippornasigvsphenocor[,"vo_2_max_change"] > 0.5,])
hippo_vo2maxneg <- rownames(hippornasigvsphenocor[hippornasigvsphenocor[,"vo_2_max_change"] < -0.5,])

# KIDNEY

kidney_wgtgainpos <- rownames(kidneyrnasigvsphenocor[kidneyrnasigvsphenocor[,"wgt_gain_after_train"] > 0.5,])
kidney_wgtgainneg <- rownames(kidneyrnasigvsphenocor[kidneyrnasigvsphenocor[,"wgt_gain_after_train"] < -0.5,])
kidney_bodyfatpos <- rownames(kidneyrnasigvsphenocor[kidneyrnasigvsphenocor[,"pct_body_fat_change"] > 0.5,])
kidney_bodyfatneg <- rownames(kidneyrnasigvsphenocor[kidneyrnasigvsphenocor[,"pct_body_fat_change"] < -0.5,])
kidney_bodyleanpos <- rownames(kidneyrnasigvsphenocor[kidneyrnasigvsphenocor[,"pct_body_lean_change"] > 0.5,])
kidney_bodyleanneg <- rownames(kidneyrnasigvsphenocor[kidneyrnasigvsphenocor[,"pct_body_lean_change"] < -0.5,])
kidney_bodyfluidpos <- rownames(kidneyrnasigvsphenocor[kidneyrnasigvsphenocor[,"pct_body_fluid_change"] > 0.5,])
kidney_bodyfluidneg <- rownames(kidneyrnasigvsphenocor[kidneyrnasigvsphenocor[,"pct_body_fluid_change"] < -0.5,])
kidney_lactatepos <- rownames(kidneyrnasigvsphenocor[kidneyrnasigvsphenocor[,"lactate_change_dueto_train"] > 0.5,])
kidney_lactateneg <- rownames(kidneyrnasigvsphenocor[kidneyrnasigvsphenocor[,"lactate_change_dueto_train"] < -0.5,])
kidney_vo2maxpos <- rownames(kidneyrnasigvsphenocor[kidneyrnasigvsphenocor[,"vo_2_max_change"] > 0.5,])
kidney_vo2maxneg <- rownames(kidneyrnasigvsphenocor[kidneyrnasigvsphenocor[,"vo_2_max_change"] < -0.5,])

# LIVER

liver_wgtgainpos <- rownames(liverrnasigvsphenocor[liverrnasigvsphenocor[,"wgt_gain_after_train"] > 0.5,])
liver_wgtgainneg <- rownames(liverrnasigvsphenocor[liverrnasigvsphenocor[,"wgt_gain_after_train"] < -0.5,])
liver_bodyfatpos <- rownames(liverrnasigvsphenocor[liverrnasigvsphenocor[,"pct_body_fat_change"] > 0.5,])
liver_bodyfatneg <- rownames(liverrnasigvsphenocor[liverrnasigvsphenocor[,"pct_body_fat_change"] < -0.5,])
liver_bodyleanpos <- rownames(liverrnasigvsphenocor[liverrnasigvsphenocor[,"pct_body_lean_change"] > 0.5,])
liver_bodyleanneg <- rownames(liverrnasigvsphenocor[liverrnasigvsphenocor[,"pct_body_lean_change"] < -0.5,])
liver_bodyfluidpos <- rownames(liverrnasigvsphenocor[liverrnasigvsphenocor[,"pct_body_fluid_change"] > 0.5,])
liver_bodyfluidneg <- rownames(liverrnasigvsphenocor[liverrnasigvsphenocor[,"pct_body_fluid_change"] < -0.5,])
liver_lactatepos <- rownames(liverrnasigvsphenocor[liverrnasigvsphenocor[,"lactate_change_dueto_train"] > 0.5,])
liver_lactateneg <- rownames(liverrnasigvsphenocor[liverrnasigvsphenocor[,"lactate_change_dueto_train"] < -0.5,])
liver_vo2maxpos <- rownames(liverrnasigvsphenocor[liverrnasigvsphenocor[,"vo_2_max_change"] > 0.5,])
liver_vo2maxneg <- rownames(liverrnasigvsphenocor[liverrnasigvsphenocor[,"vo_2_max_change"] < -0.5,])

# LUNG

lung_wgtgainpos <- rownames(lungrnasigvsphenocor[lungrnasigvsphenocor[,"wgt_gain_after_train"] > 0.5,])
lung_wgtgainneg <- rownames(lungrnasigvsphenocor[lungrnasigvsphenocor[,"wgt_gain_after_train"] < -0.5,])
lung_bodyfatpos <- rownames(lungrnasigvsphenocor[lungrnasigvsphenocor[,"pct_body_fat_change"] > 0.5,])
lung_bodyfatneg <- rownames(lungrnasigvsphenocor[lungrnasigvsphenocor[,"pct_body_fat_change"] < -0.5,])
lung_bodyleanpos <- rownames(lungrnasigvsphenocor[lungrnasigvsphenocor[,"pct_body_lean_change"] > 0.5,])
lung_bodyleanneg <- rownames(lungrnasigvsphenocor[lungrnasigvsphenocor[,"pct_body_lean_change"] < -0.5,])
lung_bodyfluidpos <- rownames(lungrnasigvsphenocor[lungrnasigvsphenocor[,"pct_body_fluid_change"] > 0.5,])
lung_bodyfluidneg <- rownames(lungrnasigvsphenocor[lungrnasigvsphenocor[,"pct_body_fluid_change"] < -0.5,])
lung_lactatepos <- rownames(lungrnasigvsphenocor[lungrnasigvsphenocor[,"lactate_change_dueto_train"] > 0.5,])
lung_lactateneg <- rownames(lungrnasigvsphenocor[lungrnasigvsphenocor[,"lactate_change_dueto_train"] < -0.5,])
lung_vo2maxpos <- rownames(lungrnasigvsphenocor[lungrnasigvsphenocor[,"vo_2_max_change"] > 0.5,])
lung_vo2maxneg <- rownames(lungrnasigvsphenocor[lungrnasigvsphenocor[,"vo_2_max_change"] < -0.5,])

# BAT

brown_wgtgainpos <- rownames(brownrnasigvsphenocor[brownrnasigvsphenocor[,"wgt_gain_after_train"] > 0.5,])
brown_wgtgainneg <- rownames(brownrnasigvsphenocor[brownrnasigvsphenocor[,"wgt_gain_after_train"] < -0.5,])
brown_bodyfatpos <- rownames(brownrnasigvsphenocor[brownrnasigvsphenocor[,"pct_body_fat_change"] > 0.5,])
brown_bodyfatneg <- rownames(brownrnasigvsphenocor[brownrnasigvsphenocor[,"pct_body_fat_change"] < -0.5,])
brown_bodyleanpos <- rownames(brownrnasigvsphenocor[brownrnasigvsphenocor[,"pct_body_lean_change"] > 0.5,])
brown_bodyleanneg <- rownames(brownrnasigvsphenocor[brownrnasigvsphenocor[,"pct_body_lean_change"] < -0.5,])
brown_bodyfluidpos <- rownames(brownrnasigvsphenocor[brownrnasigvsphenocor[,"pct_body_fluid_change"] > 0.5,])
brown_bodyfluidneg <- rownames(brownrnasigvsphenocor[brownrnasigvsphenocor[,"pct_body_fluid_change"] < -0.5,])
brown_lactatepos <- rownames(brownrnasigvsphenocor[brownrnasigvsphenocor[,"lactate_change_dueto_train"] > 0.5,])
brown_lactateneg <- rownames(brownrnasigvsphenocor[brownrnasigvsphenocor[,"lactate_change_dueto_train"] < -0.5,])
brown_vo2maxpos <- rownames(brownrnasigvsphenocor[brownrnasigvsphenocor[,"vo_2_max_change"] > 0.5,])
brown_vo2maxneg <- rownames(brownrnasigvsphenocor[brownrnasigvsphenocor[,"vo_2_max_change"] < -0.5,])

# WAT-SC

white_wgtgainpos <- rownames(whiternasigvsphenocor[whiternasigvsphenocor[,"wgt_gain_after_train"] > 0.5,])
white_wgtgainneg <- rownames(whiternasigvsphenocor[whiternasigvsphenocor[,"wgt_gain_after_train"] < -0.5,])
white_bodyfatpos <- rownames(whiternasigvsphenocor[whiternasigvsphenocor[,"pct_body_fat_change"] > 0.5,])
white_bodyfatneg <- rownames(whiternasigvsphenocor[whiternasigvsphenocor[,"pct_body_fat_change"] < -0.5,])
white_bodyleanpos <- rownames(whiternasigvsphenocor[whiternasigvsphenocor[,"pct_body_lean_change"] > 0.5,])
white_bodyleanneg <- rownames(whiternasigvsphenocor[whiternasigvsphenocor[,"pct_body_lean_change"] < -0.5,])
white_bodyfluidpos <- rownames(whiternasigvsphenocor[whiternasigvsphenocor[,"pct_body_fluid_change"] > 0.5,])
white_bodyfluidneg <- rownames(whiternasigvsphenocor[whiternasigvsphenocor[,"pct_body_fluid_change"] < -0.5,])
white_lactatepos <- rownames(whiternasigvsphenocor[whiternasigvsphenocor[,"lactate_change_dueto_train"] > 0.5,])
white_lactateneg <- rownames(whiternasigvsphenocor[whiternasigvsphenocor[,"lactate_change_dueto_train"] < -0.5,])
white_vo2maxpos <- rownames(whiternasigvsphenocor[whiternasigvsphenocor[,"vo_2_max_change"] > 0.5,])
white_vo2maxneg <- rownames(whiternasigvsphenocor[whiternasigvsphenocor[,"vo_2_max_change"] < -0.5,])

rnaphenocorrcount <- rbind(c(length(gastro_wgtgainpos),length(heart_wgtgainpos),length(hippo_wgtgainpos),length(kidney_wgtgainpos),length(liver_wgtgainpos),length(lung_wgtgainpos),length(brown_wgtgainpos),length(white_wgtgainpos)),
                           c(length(gastro_wgtgainneg),length(heart_wgtgainneg),length(hippo_wgtgainneg),length(kidney_wgtgainneg),length(liver_wgtgainneg),length(lung_wgtgainneg),length(brown_wgtgainneg),length(white_wgtgainneg)),
                           c(length(gastro_bodyfatpos),length(heart_bodyfatpos),length(hippo_bodyfatpos),length(kidney_bodyfatpos),length(liver_bodyfatpos),length(lung_bodyfatpos),length(brown_bodyfatpos),length(white_bodyfatpos)),
                           c(length(gastro_bodyfatneg),length(heart_bodyfatneg),length(hippo_bodyfatneg),length(kidney_bodyfatneg),length(liver_bodyfatneg),length(lung_bodyfatneg),length(brown_bodyfatneg),length(white_bodyfatneg)),
                           c(length(gastro_bodyleanpos),length(heart_bodyleanpos),length(hippo_bodyleanpos),length(kidney_bodyleanpos),length(liver_bodyleanpos),length(lung_bodyleanpos),length(brown_bodyleanpos),length(white_bodyleanpos)),
                           c(length(gastro_bodyleanneg),length(heart_bodyleanneg),length(hippo_bodyleanneg),length(kidney_bodyleanneg),length(liver_bodyleanneg),length(lung_bodyleanneg),length(brown_bodyleanneg),length(white_bodyleanneg)),
                           c(length(gastro_bodyfluidpos),length(heart_bodyfluidpos),length(hippo_bodyfluidpos),length(kidney_bodyfluidpos),length(liver_bodyfluidpos),length(lung_bodyfluidpos),length(brown_bodyfluidpos),length(white_bodyfluidpos)),
                           c(length(gastro_bodyfluidneg),length(heart_bodyfluidneg),length(hippo_bodyfluidneg),length(kidney_bodyfluidneg),length(liver_bodyfluidneg),length(lung_bodyfluidneg),length(brown_bodyfluidneg),length(white_bodyfluidneg)),
                           c(length(gastro_lactatepos),length(heart_lactatepos),length(hippo_lactatepos),length(kidney_lactatepos),length(liver_lactatepos),length(lung_lactatepos),length(brown_lactatepos),length(white_lactatepos)),
                           c(length(gastro_lactateneg),length(heart_lactateneg),length(hippo_lactateneg),length(kidney_lactateneg),length(liver_lactateneg),length(lung_lactateneg),length(brown_lactateneg),length(white_lactateneg)),
                           c(length(gastro_vo2maxpos),length(heart_vo2maxpos),length(hippo_vo2maxpos),length(kidney_vo2maxpos),length(liver_vo2maxpos),length(lung_vo2maxpos),length(brown_vo2maxpos),length(white_vo2maxpos)),
                           c(length(gastro_vo2maxneg),length(heart_vo2maxneg),length(hippo_vo2maxneg),length(kidney_vo2maxneg),length(liver_vo2maxneg),length(lung_vo2maxneg),length(brown_vo2maxneg),length(white_vo2maxneg)))

colnames(rnaphenocorrcount) <- c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC")
rownames(rnaphenocorrcount) <- c("WGTGAIN_POS","WGTGAIN_NEG","BODYFAT_POS","BODYFAT_NEG","BODYLEAN_POS","BODYLEAN_NEG",
                                 "BODYFLUID_POS","BODYFLUID_NEG","LACTATE_POS","LACTATE_NEG","VO2MAX_POS","VO2MAX_NEG")

rnaphenocorrfrac <- rnaphenocorrcount[1:12,]
rnaphenocorrfrac[,1] <- rnaphenocorrcount[,1]/length(gastrornasig)
rnaphenocorrfrac[,2] <- rnaphenocorrcount[,2]/length(heartrnasig)
rnaphenocorrfrac[,3] <- rnaphenocorrcount[,3]/length(hippornasig)
rnaphenocorrfrac[,4] <- rnaphenocorrcount[,4]/length(kidneyrnasig)
rnaphenocorrfrac[,5] <- rnaphenocorrcount[,5]/length(liverrnasig)
rnaphenocorrfrac[,6] <- rnaphenocorrcount[,6]/length(lungrnasig)
rnaphenocorrfrac[,7] <- rnaphenocorrcount[,7]/length(brownrnasig)
rnaphenocorrfrac[,8] <- rnaphenocorrcount[,8]/length(whiternasig)

rownames(rnaphenocorrfrac) <- c("Body Weight Pos Cor",
                                "Body Weight Neg Cor",
                                "Body Fat Pos Cor",
                                "Body Fat Neg Cor",
                                "Body Lean Pos Cor",
                                "Body Lean Neg Cor",
                                "Body Water Pos Cor",
                                "Body Water Neg Cor",
                                "Lactate Pos Cor",
                                "Lactate Neg Cor",
                                "VO2 Max Pos Cor",
                                "VO2 Max Neg Cor")
# Figure 7B
pdf(file = "Figure 7B.pdf",width = 9,height = 4)
pheatmap(rnaphenocorrfrac,angle_col = "0",cluster_rows = F,cluster_cols = F,breaks = seq(0,0.25,length.out = 101),color = colorpanel(101,"white","firebrick"),display_numbers = T,annotation_col = tissuemeta,annotation_colors = ann_colsheatmaptrim,show_colnames = F,number_color = "black",fontsize = 15,cellwidth = 44)
dev.off()

####
# Figure 7C-0
#####

# Figure 7C
gastroactiveprompeaks <- rownames(peakanno[rownames(peakanno) %in% gastroactivepeaks & peakanno$trimanno %in% "Promoter",])

gastrovo2maxpospeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastro_vo2maxpos,]),gastroactivepeaks)
gastrovo2maxposprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastro_vo2maxpos & peakanno$trimanno %in% "Promoter",]),gastroactivepeaks)

gastrovo2maxpospromdf <- data.frame("TF" = rep(rownames(table(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% gastrovo2maxposprompeak,"Motif.Name"]))[(table(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% gastrovo2maxposprompeak,"Motif.Name"]))/length(gastrovo2maxposprompeak) >= 0.0475],2),
                                    "TFabbrev" = rep(gsub("\\(.*","",rownames(table(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% gastrovo2maxposprompeak,"Motif.Name"]))[(table(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% gastrovo2maxposprompeak,"Motif.Name"]))/length(gastrovo2maxposprompeak) >= 0.0475],2)),
                                    "Measure" = c(rep("General.Freq",length(rownames(table(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% gastrovo2maxposprompeak,"Motif.Name"]))[(table(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% gastrovo2maxposprompeak,"Motif.Name"]))/length(gastrovo2maxposprompeak) >= 0.0475])),
                                                  rep("PhenoCor.Freq",length(rownames(table(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% gastrovo2maxposprompeak,"Motif.Name"]))[(table(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% gastrovo2maxposprompeak,"Motif.Name"]))/length(gastrovo2maxposprompeak) >= 0.0475]))))
gastrovo2maxpospromdf$Frequency <- 0
for(i in 1:(length(gastrovo2maxpospromdf$TF)/2)){
  ourtf <- gastrovo2maxpospromdf$TF[i]
  gastrovo2maxpospromdf$Frequency[i] <- length(unique(gastro50peakmotifs[gastro50peakmotifs$Motif.Name %in% ourtf & gastro50peakmotifs$PositionID %in% gastroactiveprompeaks,"PositionID"]))/length(intersect(gastro50peakmotifs$PositionID,gastroactiveprompeaks))
  gastrovo2maxpospromdf$Frequency[(length(gastrovo2maxpospromdf$TF)/2)+i] <- length(unique(gastro50peakmotifs[gastro50peakmotifs$Motif.Name %in% ourtf & gastro50peakmotifs$PositionID %in% gastrovo2maxposprompeak,"PositionID"]))/length(gastrovo2maxposprompeak)
}
gastrovo2maxpospromdf$Measure <- factor(gastrovo2maxpospromdf$Measure,levels = c("PhenoCor.Freq","General.Freq"))

tempanalysis <- data.frame(row.names = gastrovo2maxpospromdf[1:(dim(gastrovo2maxpospromdf)[1]/2),"TF"],
                           "General.Freq" = gastrovo2maxpospromdf[1:(dim(gastrovo2maxpospromdf)[1]/2),"Frequency"],
                           "PhenoCor.Freq" = gastrovo2maxpospromdf[((dim(gastrovo2maxpospromdf)[1]/2)+1):dim(gastrovo2maxpospromdf)[1],"Frequency"])
tempanalysis$Ratio <- tempanalysis$PhenoCor.Freq/tempanalysis$General.Freq
trimtempanalysis <- tempanalysis[tempanalysis$General.Freq > 0.045 | tempanalysis$PhenoCor.Freq > 0.045,]
trimtempanalysis <- trimtempanalysis[!rownames(trimtempanalysis) %in% "Nr5a2(NR)/Pancreas-LRH1-ChIP-Seq(GSE34295)/Homer",]
gastrovo2maxpospromdftrim <- gastrovo2maxpospromdf[gastrovo2maxpospromdf$TF %in% rownames(trimtempanalysis),]
#gastrovo2maxpospromdftrim$Ratio <- c(trimtempanalysis$Ratio,trimtempanalysis$Ratio)
#gastrovo2maxpospromdftrim <- gastrovo2maxpospromdftrim[order(gastrovo2maxpospromdftrim$Ratio),]
gastrovo2maxpospromdftrim$TF <- as.character(gastrovo2maxpospromdftrim$TF)
gastrovo2maxpospromdftrim$TF <- factor(gastrovo2maxpospromdftrim$TF,levels = rownames(trimtempanalysis)[order(trimtempanalysis$Ratio)])
gastrovo2maxpospromdftrim$TFabbrev <- factor(gastrovo2maxpospromdftrim$TFabbrev,levels = gsub("\\(.*","",rownames(trimtempanalysis))[order(trimtempanalysis$Ratio)])

gastrovo2maxpospromdftrim$Significance <- rep("N",dim(gastrovo2maxpospromdftrim)[1])
gastrovo2maxpospromdftrim$MaxFrequency <- 0
for(i in 1:(dim(gastrovo2maxpospromdftrim)[1]/2)){
  gastrovo2maxpospromdftrim$MaxFrequency[i] <- max(gastrovo2maxpospromdftrim$Frequency[i],gastrovo2maxpospromdftrim$Frequency[i+(dim(gastrovo2maxpospromdftrim)[1]/2)])
  gastrovo2maxpospromdftrim$MaxFrequency[i+(dim(gastrovo2maxpospromdftrim)[1]/2)] <- max(gastrovo2maxpospromdftrim$Frequency[i],gastrovo2maxpospromdftrim$Frequency[i+(dim(gastrovo2maxpospromdftrim)[1]/2)])
}

# The two TFs that pass an exact binomial test for % enrichment 
gastrovo2maxpospromdftrim["52","Significance"] <- "Y"
gastrovo2maxpospromdftrim["62","Significance"] <- "Y"

pdf(file = "Figure 7C.pdf",width = 10,height = 5)
ggplot(data = gastrovo2maxpospromdftrim, aes(x = TFabbrev,y = Frequency,group = Measure,color = Measure)) + geom_line(size = 2) + geom_point(size = 3) + theme_classic() + theme(axis.text.x = element_text(size = 15,angle = 315, vjust = 0,hjust = 0.25),axis.text.y = element_text(size = 15),axis.title = element_text(size = 18),legend.text = element_text(size = 15),legend.title = element_text(size = 18),title = element_text(size = 20,hjust = 0.5)) + xlab("Transcription Factor") + ggtitle("TF Enrich in Pos Corr DEGs with VO2 Max in SKM-GN") + geom_point(data = gastrovo2maxpospromdftrim[gastrovo2maxpospromdftrim$Significance == "Y", ], aes(TFabbrev, MaxFrequency + 0.02), shape = "*", size=10, color="black")
dev.off()


# Figure 7E
ourgene <- enstosym[enstosym$Symbol %in% "Me3","Ensembl"]
ourdf <- trimgastrophenodf
ourdf$Gene.Expr <- t(gastrornasignorm[ourgene,rownames(trimgastrophenodf)])
pdf(file = "Figure 7E.pdf",width = 7,height = 4)
ggplot(ourdf,aes(x=Gene.Expr,y=vo2_max_change,shape=sex,color=group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$vo2_max_change,ourdf$Gene.Expr),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("VO2 Max Change") + xlab("Gene Expression")
dev.off()

# Figure 7F
ourgene <- enstosym[enstosym$Symbol %in% "Rora","Ensembl"]
ourdf <- trimgastrophenodf
ourdf$Gene.Expr <- t(gastrornasignorm[ourgene,rownames(trimgastrophenodf)])
pdf(file = "Figure 7F.pdf",width = 7,height = 4)
ggplot(ourdf,aes(x=Gene.Expr,y=vo2_max_change,shape=sex,color=group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$vo2_max_change,ourdf$Gene.Expr),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("VO2 Max Change") + xlab("Gene Expression")
dev.off()

# Figure 7G
ourgene <- enstosym[enstosym$Symbol %in% "Lgi3","Ensembl"]
ourdf <- trimgastrophenodf
ourdf$Gene.Expr <- t(gastrornasignorm[ourgene,rownames(trimgastrophenodf)])
pdf(file = "Figure 7G.pdf",width = 7,height = 4)
ggplot(ourdf,aes(x=Gene.Expr,y=vo2_max_change,shape=sex,color=group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$vo2_max_change,ourdf$Gene.Expr),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("VO2 Max Change") + xlab("Gene Expression")
dev.off()

# Figure 7H
gastrowgtgainpospeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastro_wgtgainpos,]),gastroactivepeaks)
gastrowgtgainposprompeak <- intersect(rownames(peakanno[peakanno$ensembl_gene %in% gastro_wgtgainpos & peakanno$trimanno %in% "Promoter",]),gastroactivepeaks)

gastrowgtgainpospromdf <- data.frame("TF" = rep(rownames(table(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% gastrowgtgainposprompeak,"Motif.Name"]))[(table(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% gastrowgtgainposprompeak,"Motif.Name"]))/length(gastrowgtgainposprompeak) >= 0.05],2),
                                     "TFabbrev" = rep(gsub("\\(.*","",rownames(table(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% gastrowgtgainposprompeak,"Motif.Name"]))[(table(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% gastrowgtgainposprompeak,"Motif.Name"]))/length(gastrowgtgainposprompeak) >= 0.05],2)),
                                     "Measure" = c(rep("General.Freq",length(rownames(table(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% gastrowgtgainposprompeak,"Motif.Name"]))[(table(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% gastrowgtgainposprompeak,"Motif.Name"]))/length(gastrowgtgainposprompeak) >= 0.05])),
                                                   rep("PhenoCor.Freq",length(rownames(table(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% gastrowgtgainposprompeak,"Motif.Name"]))[(table(gastro50peakmotifs[gastro50peakmotifs$PositionID %in% gastrowgtgainposprompeak,"Motif.Name"]))/length(gastrowgtgainposprompeak) >= 0.05]))))
gastrowgtgainpospromdf$Frequency <- 0
for(i in 1:(length(gastrowgtgainpospromdf$TF)/2)){
  ourtf <- gastrowgtgainpospromdf$TF[i]
  gastrowgtgainpospromdf$Frequency[i] <- length(unique(gastro50peakmotifs[gastro50peakmotifs$Motif.Name %in% ourtf & gastro50peakmotifs$PositionID %in% gastroactiveprompeaks,"PositionID"]))/length(intersect(gastro50peakmotifs$PositionID,gastroactiveprompeaks))
  gastrowgtgainpospromdf$Frequency[(length(gastrowgtgainpospromdf$TF)/2)+i] <- length(unique(gastro50peakmotifs[gastro50peakmotifs$Motif.Name %in% ourtf & gastro50peakmotifs$PositionID %in% gastrowgtgainposprompeak,"PositionID"]))/length(gastrowgtgainposprompeak)
}
gastrowgtgainpospromdf$Measure <- factor(gastrowgtgainpospromdf$Measure,levels = c("PhenoCor.Freq","General.Freq"))

tempanalysis <- data.frame(row.names = gastrowgtgainpospromdf[1:(dim(gastrowgtgainpospromdf)[1]/2),"TF"],
                           "General.Freq" = gastrowgtgainpospromdf[1:(dim(gastrowgtgainpospromdf)[1]/2),"Frequency"],
                           "PhenoCor.Freq" = gastrowgtgainpospromdf[((dim(gastrowgtgainpospromdf)[1]/2)+1):dim(gastrowgtgainpospromdf)[1],"Frequency"])
tempanalysis$Ratio <- tempanalysis$PhenoCor.Freq/tempanalysis$General.Freq
trimtempanalysis <- tempanalysis[tempanalysis$General.Freq > 0.05 | tempanalysis$PhenoCor.Freq > 0.05,]
trimtempanalysis <- trimtempanalysis[!rownames(trimtempanalysis) %in% "COUP-TFII(NR)/Artia-Nr2f2-ChIP-Seq(GSE46497)/Homer",]
gastrowgtgainpospromdftrim <- gastrowgtgainpospromdf[gastrowgtgainpospromdf$TF %in% rownames(trimtempanalysis),]
gastrowgtgainpospromdftrim$TF <- as.character(gastrowgtgainpospromdftrim$TF)
gastrowgtgainpospromdftrim$TF <- factor(gastrowgtgainpospromdftrim$TF,levels = rownames(trimtempanalysis)[order(trimtempanalysis$Ratio)])
gastrowgtgainpospromdftrim$TFabbrev <- factor(gastrowgtgainpospromdftrim$TFabbrev,levels = gsub("\\(.*","",rownames(trimtempanalysis))[order(trimtempanalysis$Ratio)])

gastrowgtgainpospromdftrim$Significance <- rep("N",dim(gastrowgtgainpospromdftrim)[1])
gastrowgtgainpospromdftrim$MaxFrequency <- 0
for(i in 1:(dim(gastrowgtgainpospromdftrim)[1]/2)){
  gastrowgtgainpospromdftrim$MaxFrequency[i] <- max(gastrowgtgainpospromdftrim$Frequency[i],gastrowgtgainpospromdftrim$Frequency[i+(dim(gastrowgtgainpospromdftrim)[1]/2)])
  gastrowgtgainpospromdftrim$MaxFrequency[i+(dim(gastrowgtgainpospromdftrim)[1]/2)] <- max(gastrowgtgainpospromdftrim$Frequency[i],gastrowgtgainpospromdftrim$Frequency[i+(dim(gastrowgtgainpospromdftrim)[1]/2)])
}
gastrowgtgainpospromdftrim["2","Significance"] <- "Y"
gastrowgtgainpospromdftrim["52","Significance"] <- "Y"

pdf(file = "Figure 7H.pdf",width = 10,height = 5)
ggplot(data = gastrowgtgainpospromdftrim, aes(x = TFabbrev,y = Frequency,group = Measure,color = Measure)) + geom_line(size = 2) + geom_point(size = 3) + theme_classic() + theme(axis.text.x = element_text(size = 15,angle = 315, vjust = 0,hjust = 0.25),axis.text.y = element_text(size = 15),axis.title = element_text(size = 18),legend.text = element_text(size = 15),legend.title = element_text(size = 18),title = element_text(size = 20,hjust = 0.5)) + xlab("Transcription Factor") + ggtitle("TF Enrich in Pos Corr DEGs with Body Weight in SKM-GN") + geom_point(data = gastrowgtgainpospromdftrim[gastrowgtgainpospromdftrim$Significance == "Y", ], aes(TFabbrev, MaxFrequency + 0.02), shape = "*", size=10, color="black")
dev.off()


# Figure 7J
ourgene <- enstosym[enstosym$Symbol %in% "Chd7","Ensembl"]
ourdf <- trimgastrophenodf
ourdf$Gene.Expr <- t(gastrornasignorm[ourgene,rownames(trimgastrophenodf)])
pdf(file = "Figure 7J.pdf",width = 7,height = 4)
ggplot(ourdf,aes(x=Gene.Expr,y=wgt_gain_after_train,shape=sex,color=group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$wgt_gain_after_train,ourdf$Gene.Expr),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("Body Weight Change") + xlab("Gene Expression")
dev.off()

# Figure 7K
ourgene <- enstosym[enstosym$Symbol %in% "Igf2","Ensembl"]
ourdf <- trimgastrophenodf
ourdf$Gene.Expr <- t(gastrornasignorm[ourgene,rownames(trimgastrophenodf)])
pdf(file = "Figure 7K.pdf",width = 7,height = 4)
ggplot(ourdf,aes(x=Gene.Expr,y=pct_body_fat_change,shape=sex,color=group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$pct_body_fat_change,ourdf$Gene.Expr),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("Body Fat Change") + xlab("Gene Expression")
dev.off()

# Figure 7L
ourgene <- enstosym[enstosym$Symbol %in% "Sall2","Ensembl"]
ourdf <- trimgastrophenodf
ourdf$Gene.Expr <- t(gastrornasignorm[ourgene,rownames(trimgastrophenodf)])
pdf(file = "Figure 7L.pdf",width = 7,height = 4)
ggplot(ourdf,aes(x=Gene.Expr,y=pct_body_fat_change,shape=sex,color=group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$pct_body_fat_change,ourdf$Gene.Expr),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("Body Fat Change") + xlab("Gene Expression")
dev.off()

# Figure 7M
ourgene <- enstosym[enstosym$Symbol %in% "Oas2","Ensembl"]
ourdf <- trimlungphenodf
ourdf$Gene.Expr <- t(lungrnasignorm[ourgene,rownames(trimlungphenodf)])
pdf(file = "Figure 7M.pdf",width = 7,height = 4)
ggplot(ourdf,aes(x=Gene.Expr,y=pct_body_fat_change,shape=sex,color=group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LUNG ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$pct_body_fat_change,ourdf$Gene.Expr),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("Body Fat Change") + xlab("Gene Expression")
dev.off()

# Figure 7N
ourgene <- enstosym[enstosym$Symbol %in% "Nfkb2","Ensembl"]
ourdf <- trimlungphenodf
ourdf$Gene.Expr <- t(lungrnasignorm[ourgene,rownames(trimlungphenodf)])
pdf(file = "Figure 7N.pdf",width = 7,height = 4)
ggplot(ourdf,aes(x=Gene.Expr,y=pct_body_fat_change,shape=sex,color=group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LUNG ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$pct_body_fat_change,ourdf$Gene.Expr),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("Body Fat Change") + xlab("Gene Expression")
dev.off()

# Figure 7O
ourgene <- enstosym[enstosym$Symbol %in% "Fkbp4","Ensembl"]
ourdf <- trimliverphenodf
ourdf$Gene.Expr <- t(liverrnasignorm[ourgene,rownames(trimliverphenodf)])
pdf(file = "Figure 7O.pdf",width = 7,height = 4)
ggplot(ourdf,aes(x=Gene.Expr,y=pct_body_fat_change,shape=sex,color=group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LIVER ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$pct_body_fat_change,ourdf$Gene.Expr),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("Body Fat Change") + xlab("Gene Expression")
dev.off()

####
# Supplemental Figure S16
#####

####
# Supplemental Figure S16A
#####

gastrornatfvsphenocor <- gastrornasigvsphenocor[intersect(rownames(gastrornasigvsphenocor),tfanno$Ensembl),]
heartrnatfvsphenocor <- heartrnasigvsphenocor[intersect(rownames(heartrnasigvsphenocor),tfanno$Ensembl),]
hippornatfvsphenocor <- hippornasigvsphenocor[intersect(rownames(hippornasigvsphenocor),tfanno$Ensembl),]
kidneyrnatfvsphenocor <- kidneyrnasigvsphenocor[intersect(rownames(kidneyrnasigvsphenocor),tfanno$Ensembl),]
liverrnatfvsphenocor <- liverrnasigvsphenocor[intersect(rownames(liverrnasigvsphenocor),tfanno$Ensembl),]
lungrnatfvsphenocor <- lungrnasigvsphenocor[intersect(rownames(lungrnasigvsphenocor),tfanno$Ensembl),]
brownrnatfvsphenocor <- brownrnasigvsphenocor[intersect(rownames(brownrnasigvsphenocor),tfanno$Ensembl),]
whiternatfvsphenocor <- whiternasigvsphenocor[intersect(rownames(whiternasigvsphenocor),tfanno$Ensembl),]

gastrornatfvsphenocorsub <- rownames(gastrornatfvsphenocor)[apply(abs(gastrornatfvsphenocor),1,max) > 0.5]
heartrnatfvsphenocorsub <- rownames(heartrnatfvsphenocor)[apply(abs(heartrnatfvsphenocor),1,max) > 0.5]
hippornatfvsphenocorsub <- rownames(hippornatfvsphenocor)[apply(abs(hippornatfvsphenocor),1,max) > 0.5]
kidneyrnatfvsphenocorsub <- rownames(kidneyrnatfvsphenocor)[apply(abs(kidneyrnatfvsphenocor),1,max) > 0.5]
liverrnatfvsphenocorsub <- rownames(liverrnatfvsphenocor)[apply(abs(liverrnatfvsphenocor),1,max) > 0.5]
lungrnatfvsphenocorsub <- rownames(lungrnatfvsphenocor)[apply(abs(lungrnatfvsphenocor),1,max) > 0.5]
brownrnatfvsphenocorsub <- rownames(brownrnatfvsphenocor)[apply(abs(brownrnatfvsphenocor),1,max) > 0.5]
whiternatfvsphenocorsub <- rownames(whiternatfvsphenocor)[apply(abs(whiternatfvsphenocor),1,max) > 0.5]

combornatfvsphenocor <- rbind(gastrornatfvsphenocor[gastrornatfvsphenocorsub,],
                              kidneyrnatfvsphenocor[kidneyrnatfvsphenocorsub,],
                              liverrnatfvsphenocor[liverrnatfvsphenocorsub,],
                              lungrnatfvsphenocor[lungrnatfvsphenocorsub,],
                              brownrnatfvsphenocor[brownrnatfvsphenocorsub,],
                              whiternatfvsphenocor[whiternatfvsphenocorsub,])
combotfrnaidlist <- c(gastrornatfvsphenocorsub,
                      kidneyrnatfvsphenocorsub,
                      liverrnatfvsphenocorsub,
                      lungrnatfvsphenocorsub,
                      brownrnatfvsphenocorsub,
                      whiternatfvsphenocorsub)
rownames(combornatfvsphenocor)[1:4] <- paste("gastro_",rownames(combornatfvsphenocor)[1:4],sep = "")
rownames(combornatfvsphenocor)[5:6] <- paste("kidney_",rownames(combornatfvsphenocor)[5:6],sep = "")
rownames(combornatfvsphenocor)[7:8] <- paste("liver_",rownames(combornatfvsphenocor)[7:8],sep = "")
rownames(combornatfvsphenocor)[9:11] <- paste("lung_",rownames(combornatfvsphenocor)[9:11],sep = "")
rownames(combornatfvsphenocor)[12:14] <- paste("brown_",rownames(combornatfvsphenocor)[12:14],sep = "")
rownames(combornatfvsphenocor)[15:28] <- paste("white_",rownames(combornatfvsphenocor)[15:28],sep = "")
combornatfvsphenometa <- data.frame(row.names = rownames(combornatfvsphenocor),
                                    "Tissue" = c(rep("SKM-GN",4),
                                                 rep("KIDNEY",2),
                                                 rep("LIVER",2),
                                                 rep("LUNG",3),
                                                 rep("BAT",3),
                                                 rep("WAT-SC",14)))
combotfrnanamelist <- enstosym[combotfrnaidlist,"Symbol"]

# Supplemental Figure S16A
pdf(file = "Supplemental Figure S16A.pdf",height = 7,width = 8)
pheatmap(combornatfvsphenocor,cluster_rows = T,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),annotation_row = combornatfvsphenometa,annotation_col = sigvsphenocormeta,annotation_colors = ann_cols_sigphenocor,show_colnames = F,display_numbers = T,number_color = "black",labels_row = gsub("\\(.*","",combotfrnanamelist),fontsize = 15)
dev.off()

####
# Supplemental Figure S16B
#####

# SKM-GN

gastroprosig <- unique(prot_pr$training_dea[prot_pr$training_dea$adj_p_value < 0.1 & prot_pr$training_dea$tissue_abbreviation %in% "SKM-GN","feature_ID"])

gastropronorm <- read.delim(file = "November 2021 PASS1B Data Freeze/Proteomics Data/Prot_PR/pass1b-06_v1.0_analysis_proteomics-untargeted_prot-pr_normalized-data_motrpac_pass1b-06_t55-gastro.txt",header = T,row.names = 1)
colnames(gastropronorm) <- gsub("X","",colnames(gastropronorm))
gastroprosignorm <- gastropronorm[gastroprosig,]

gastroprosignorm <- gastroprosignorm[rowSums(is.na(gastroprosignorm)) == 0,]
gastroprosig <- rownames(gastroprosignorm)



gastroprophenomeasuredata <- phenomeasuredata[intersect(gsub("X","",colnames(gastropronorm)),rownames(phenomeasuredata)),]
gastroprophenomeasuredata <- gastroprophenomeasuredata[,c("calculated_variables___wgt_gain_after_train","calculated_variables___pct_body_fat_change","calculated_variables___pct_body_lean_change","calculated_variables___pct_body_fluid_change","calculated_variables___lactate_change_dueto_train","calculated_variables___vo_2_max_change")]

gastroprometa <- data.frame(row.names = colnames(gastroprosignorm),
                            "label" = colnames(gastroprosignorm))


gastroprophenodf <- data.frame(row.names = colnames(gastroprosignorm),"label" = colnames(gastroprosignorm))
gastroprophenodf$wgt_gain_after_train <- 0
gastroprophenodf$pct_body_fat_change <- 0
gastroprophenodf$pct_body_lean_change <- 0
gastroprophenodf$pct_body_fluid_change <- 0
gastroprophenodf$lactate_change_dueto_train <- 0
gastroprophenodf$vo2_max_change <- 0
gastroprophenodf$sex <- ""
gastroprophenodf$group <- ""
for(i in 1:dim(gastroprophenodf)[1]){
  ourid <- rownames(gastroprophenodf)[i]
  gastroprophenodf[ourid,"wgt_gain_after_train"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___wgt_gain_after_train"])
  gastroprophenodf[ourid,"pct_body_fat_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_fat_change"])
  gastroprophenodf[ourid,"pct_body_lean_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_lean_change"])
  gastroprophenodf[ourid,"pct_body_fluid_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_fluid_change"])
  gastroprophenodf[ourid,"lactate_change_dueto_train"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___lactate_change_dueto_train"])
  gastroprophenodf[ourid,"vo2_max_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___vo_2_max_change"])
  gastroprophenodf[ourid,"sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"][1]]
  gastroprophenodf[ourid,"group"] <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1]
}
trimgastroprophenodf <- gastroprophenodf[gastroprophenodf$group %in% c("Four-week program",
                                                                       "Eight-week program Training Group",
                                                                       "Eight-week program Control Group"),]
trimgastroprophenodf[trimgastroprophenodf$group %in% "Eight-week program Control Group","group"] <- "Control"
trimgastroprophenodf[trimgastroprophenodf$group %in% "Eight-week program Training Group","group"] <- "8w"
trimgastroprophenodf[trimgastroprophenodf$group %in% "Four-week program","group"] <- "4w"
trimgastroprophenodf$group <- factor(trimgastroprophenodf$group,levels = c("Control","4w","8w"))


gastroprosigvsphenocor <- matrix(0L,nrow = length(gastroprosig),ncol = dim(gastroprophenomeasuredata)[2])
rownames(gastroprosigvsphenocor) <- gastroprosig
colnames(gastroprosigvsphenocor) <- colnames(gastroprophenomeasuredata)
for(i in 1:length(gastroprosig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:dim(gastroprophenomeasuredata)[2]){
    gastroprosigvsphenocor[i,j] <- cor(t(gastroprosignorm[i,rownames(trimgastroprophenodf)]),gastroprophenomeasuredata[rownames(trimgastroprophenodf),j])
  }
}
colnames(gastroprosigvsphenocor) <- gsub("calculated_variables___","",colnames(gastroprosigvsphenocor))


gastroprotfnorm <- gastropronorm[intersect(rownames(gastropronorm),tfproanno$Gastro.Pro.ID),]
gastroprotfnorm <- gastroprotfnorm[rowSums(is.na(gastroprotfnorm)) == 0,]

gastroprotfvsphenocor <- matrix(0L,nrow = dim(gastroprotfnorm)[1],ncol = dim(gastroprophenomeasuredata)[2])
rownames(gastroprotfvsphenocor) <- rownames(gastroprotfnorm)
colnames(gastroprotfvsphenocor) <- colnames(gastroprophenomeasuredata)
for(i in 1:length(rownames(gastroprotfnorm))){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:dim(gastroprophenomeasuredata)[2]){
    gastroprotfvsphenocor[i,j] <- cor(t(gastroprotfnorm[i,rownames(trimgastroprophenodf)]),gastroprophenomeasuredata[rownames(trimgastroprophenodf),j])
  }
}
colnames(gastroprotfvsphenocor) <- gsub("calculated_variables___","",colnames(gastroprotfvsphenocor))


# HEART

heartprosig <- unique(prot_pr$training_dea[prot_pr$training_dea$adj_p_value < 0.1 & prot_pr$training_dea$tissue_abbreviation %in% "HEART","feature_ID"])

heartpronorm <- read.delim(file = "November 2021 PASS1B Data Freeze/Proteomics Data/Prot_PR/pass1b-06_v1.1_analysis_proteomics-untargeted_prot-pr_normalized-data_motrpac_pass1b-06_t58-heart_prot.txt",header = T,row.names = 1)
colnames(heartpronorm) <- gsub("X","",colnames(heartpronorm))
heartprosignorm <- heartpronorm[heartprosig,]

heartprosignorm <- heartprosignorm[rowSums(is.na(heartprosignorm)) == 0,]
heartprosig <- rownames(heartprosignorm)



heartprophenomeasuredata <- phenomeasuredata[intersect(gsub("X","",colnames(heartpronorm)),rownames(phenomeasuredata)),]
heartprophenomeasuredata <- heartprophenomeasuredata[,c("calculated_variables___wgt_gain_after_train","calculated_variables___pct_body_fat_change","calculated_variables___pct_body_lean_change","calculated_variables___pct_body_fluid_change","calculated_variables___lactate_change_dueto_train","calculated_variables___vo_2_max_change")]

heartprometa <- data.frame(row.names = colnames(heartprosignorm),
                           "label" = colnames(heartprosignorm))


heartprophenodf <- data.frame(row.names = colnames(heartprosignorm),"label" = colnames(heartprosignorm))
heartprophenodf$wgt_gain_after_train <- 0
heartprophenodf$pct_body_fat_change <- 0
heartprophenodf$pct_body_lean_change <- 0
heartprophenodf$pct_body_fluid_change <- 0
heartprophenodf$lactate_change_dueto_train <- 0
heartprophenodf$vo2_max_change <- 0
heartprophenodf$sex <- ""
heartprophenodf$group <- ""
for(i in 1:dim(heartprophenodf)[1]){
  ourid <- rownames(heartprophenodf)[i]
  ourpid <- phenomeasuredata[ourid,"pid"]
  heartprophenodf[ourid,"wgt_gain_after_train"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___wgt_gain_after_train"])
  heartprophenodf[ourid,"pct_body_fat_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_fat_change"])
  heartprophenodf[ourid,"pct_body_lean_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_lean_change"])
  heartprophenodf[ourid,"pct_body_fluid_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_fluid_change"])
  heartprophenodf[ourid,"lactate_change_dueto_train"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___lactate_change_dueto_train"])
  heartprophenodf[ourid,"vo2_max_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___vo_2_max_change"])
  heartprophenodf[ourid,"sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourpid,"pass1bf0027"][1]]
  heartprophenodf[ourid,"group"] <- pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourpid,"pass1bf0011"][1]
}
trimheartprophenodf <- heartprophenodf[heartprophenodf$group %in% c("Four-week program",
                                                                    "Eight-week program Training Group",
                                                                    "Eight-week program Control Group"),]
trimheartprophenodf[trimheartprophenodf$group %in% "Eight-week program Control Group","group"] <- "Control"
trimheartprophenodf[trimheartprophenodf$group %in% "Eight-week program Training Group","group"] <- "8w"
trimheartprophenodf[trimheartprophenodf$group %in% "Four-week program","group"] <- "4w"
trimheartprophenodf$group <- factor(trimheartprophenodf$group,levels = c("Control","4w","8w"))


heartprosigvsphenocor <- matrix(0L,nrow = length(heartprosig),ncol = dim(heartprophenomeasuredata)[2])
rownames(heartprosigvsphenocor) <- heartprosig
colnames(heartprosigvsphenocor) <- colnames(heartprophenomeasuredata)
for(i in 1:length(heartprosig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:dim(heartprophenomeasuredata)[2]){
    heartprosigvsphenocor[i,j] <- cor(t(heartprosignorm[i,rownames(trimheartprophenodf)]),heartprophenomeasuredata[rownames(trimheartprophenodf),j])
  }
}
colnames(heartprosigvsphenocor) <- gsub("calculated_variables___","",colnames(heartprosigvsphenocor))


heartprotfnorm <- heartpronorm[intersect(rownames(heartpronorm),tfproanno$Heart.Pro.ID),]
heartprotfnorm <- heartprotfnorm[rowSums(is.na(heartprotfnorm)) == 0,]

heartprotfvsphenocor <- matrix(0L,nrow = dim(heartprotfnorm)[1],ncol = dim(heartprophenomeasuredata)[2])
rownames(heartprotfvsphenocor) <- rownames(heartprotfnorm)
colnames(heartprotfvsphenocor) <- colnames(heartprophenomeasuredata)
for(i in 1:length(rownames(heartprotfnorm))){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:dim(heartprophenomeasuredata)[2]){
    heartprotfvsphenocor[i,j] <- cor(t(heartprotfnorm[i,rownames(trimheartprophenodf)]),heartprophenomeasuredata[rownames(trimheartprophenodf),j])
  }
}
colnames(heartprotfvsphenocor) <- gsub("calculated_variables___","",colnames(heartprotfvsphenocor))


# KIDNEY

kidneyprosig <- unique(prot_pr$training_dea[prot_pr$training_dea$adj_p_value < 0.1 & prot_pr$training_dea$tissue_abbreviation %in% "KIDNEY","feature_ID"])

kidneypronorm <- read.delim(file = "November 2021 PASS1B Data Freeze/Proteomics Data/Prot_PR/pass1b-06_v1.1_analysis_proteomics-untargeted_prot-pr_normalized-data_motrpac_pass1b-06_t59-kidney.txt",header = T,row.names = 1)
colnames(kidneypronorm) <- gsub("X","",colnames(kidneypronorm))
kidneyprosignorm <- kidneypronorm[kidneyprosig,]

kidneyprosignorm <- kidneyprosignorm[rowSums(is.na(kidneyprosignorm)) == 0,]
kidneyprosig <- rownames(kidneyprosignorm)



kidneyprophenomeasuredata <- phenomeasuredata[intersect(gsub("X","",colnames(kidneypronorm)),rownames(phenomeasuredata)),]
kidneyprophenomeasuredata <- kidneyprophenomeasuredata[,c("calculated_variables___wgt_gain_after_train","calculated_variables___pct_body_fat_change","calculated_variables___pct_body_lean_change","calculated_variables___pct_body_fluid_change","calculated_variables___lactate_change_dueto_train","calculated_variables___vo_2_max_change")]

kidneyprometa <- data.frame(row.names = colnames(kidneyprosignorm),
                            "label" = colnames(kidneyprosignorm))


kidneyprophenodf <- data.frame(row.names = colnames(kidneyprosignorm),"label" = colnames(kidneyprosignorm))
kidneyprophenodf$wgt_gain_after_train <- 0
kidneyprophenodf$pct_body_fat_change <- 0
kidneyprophenodf$pct_body_lean_change <- 0
kidneyprophenodf$pct_body_fluid_change <- 0
kidneyprophenodf$lactate_change_dueto_train <- 0
kidneyprophenodf$vo2_max_change <- 0
kidneyprophenodf$sex <- ""
kidneyprophenodf$group <- ""
for(i in 1:dim(kidneyprophenodf)[1]){
  ourid <- rownames(kidneyprophenodf)[i]
  ourpid <- phenomeasuredata[ourid,"pid"]
  kidneyprophenodf[ourid,"wgt_gain_after_train"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___wgt_gain_after_train"])
  kidneyprophenodf[ourid,"pct_body_fat_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_fat_change"])
  kidneyprophenodf[ourid,"pct_body_lean_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_lean_change"])
  kidneyprophenodf[ourid,"pct_body_fluid_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_fluid_change"])
  kidneyprophenodf[ourid,"lactate_change_dueto_train"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___lactate_change_dueto_train"])
  kidneyprophenodf[ourid,"vo2_max_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___vo_2_max_change"])
  kidneyprophenodf[ourid,"sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourpid,"pass1bf0027"][1]]
  kidneyprophenodf[ourid,"group"] <- pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourpid,"pass1bf0011"][1]
}
trimkidneyprophenodf <- kidneyprophenodf[kidneyprophenodf$group %in% c("Four-week program",
                                                                       "Eight-week program Training Group",
                                                                       "Eight-week program Control Group"),]
trimkidneyprophenodf[trimkidneyprophenodf$group %in% "Eight-week program Control Group","group"] <- "Control"
trimkidneyprophenodf[trimkidneyprophenodf$group %in% "Eight-week program Training Group","group"] <- "8w"
trimkidneyprophenodf[trimkidneyprophenodf$group %in% "Four-week program","group"] <- "4w"
trimkidneyprophenodf$group <- factor(trimkidneyprophenodf$group,levels = c("Control","4w","8w"))


kidneyprosigvsphenocor <- matrix(0L,nrow = length(kidneyprosig),ncol = dim(kidneyprophenomeasuredata)[2])
rownames(kidneyprosigvsphenocor) <- kidneyprosig
colnames(kidneyprosigvsphenocor) <- colnames(kidneyprophenomeasuredata)
for(i in 1:length(kidneyprosig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:dim(kidneyprophenomeasuredata)[2]){
    kidneyprosigvsphenocor[i,j] <- cor(t(kidneyprosignorm[i,rownames(trimkidneyprophenodf)]),kidneyprophenomeasuredata[rownames(trimkidneyprophenodf),j])
  }
}
colnames(kidneyprosigvsphenocor) <- gsub("calculated_variables___","",colnames(kidneyprosigvsphenocor))


kidneyprotfnorm <- kidneypronorm[intersect(rownames(kidneypronorm),tfproanno$Kidney.Pro.ID),]
kidneyprotfnorm <- kidneyprotfnorm[rowSums(is.na(kidneyprotfnorm)) == 0,]

kidneyprotfvsphenocor <- matrix(0L,nrow = dim(kidneyprotfnorm)[1],ncol = dim(kidneyprophenomeasuredata)[2])
rownames(kidneyprotfvsphenocor) <- rownames(kidneyprotfnorm)
colnames(kidneyprotfvsphenocor) <- colnames(kidneyprophenomeasuredata)
for(i in 1:length(rownames(kidneyprotfnorm))){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:dim(kidneyprophenomeasuredata)[2]){
    kidneyprotfvsphenocor[i,j] <- cor(t(kidneyprotfnorm[i,rownames(trimkidneyprophenodf)]),kidneyprophenomeasuredata[rownames(trimkidneyprophenodf),j])
  }
}
colnames(kidneyprotfvsphenocor) <- gsub("calculated_variables___","",colnames(kidneyprotfvsphenocor))


# LIVER

liverprosig <- unique(prot_pr$training_dea[prot_pr$training_dea$adj_p_value < 0.1 & prot_pr$training_dea$tissue_abbreviation %in% "LIVER","feature_ID"])

liverpronorm <- read.delim(file = "November 2021 PASS1B Data Freeze/Proteomics Data/Prot_PR/pass1b-06_v1.1_analysis_proteomics-untargeted_prot-pr_normalized-data_motrpac_pass1b-06_t68-liver_prot.txt",header = T,row.names = 1)
colnames(liverpronorm) <- gsub("X","",colnames(liverpronorm))
liverprosignorm <- liverpronorm[liverprosig,]

liverprosignorm <- liverprosignorm[rowSums(is.na(liverprosignorm)) == 0,]
liverprosig <- rownames(liverprosignorm)



liverprophenomeasuredata <- phenomeasuredata[intersect(gsub("X","",colnames(liverpronorm)),rownames(phenomeasuredata)),]
liverprophenomeasuredata <- liverprophenomeasuredata[,c("calculated_variables___wgt_gain_after_train","calculated_variables___pct_body_fat_change","calculated_variables___pct_body_lean_change","calculated_variables___pct_body_fluid_change","calculated_variables___lactate_change_dueto_train","calculated_variables___vo_2_max_change")]

liverprometa <- data.frame(row.names = colnames(liverprosignorm),
                           "label" = colnames(liverprosignorm))


liverprophenodf <- data.frame(row.names = colnames(liverprosignorm),"label" = colnames(liverprosignorm))
liverprophenodf$wgt_gain_after_train <- 0
liverprophenodf$pct_body_fat_change <- 0
liverprophenodf$pct_body_lean_change <- 0
liverprophenodf$pct_body_fluid_change <- 0
liverprophenodf$lactate_change_dueto_train <- 0
liverprophenodf$vo2_max_change <- 0
liverprophenodf$sex <- ""
liverprophenodf$group <- ""
for(i in 1:dim(liverprophenodf)[1]){
  ourid <- rownames(liverprophenodf)[i]
  ourpid <- phenomeasuredata[ourid,"pid"]
  liverprophenodf[ourid,"wgt_gain_after_train"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___wgt_gain_after_train"])
  liverprophenodf[ourid,"pct_body_fat_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_fat_change"])
  liverprophenodf[ourid,"pct_body_lean_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_lean_change"])
  liverprophenodf[ourid,"pct_body_fluid_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_fluid_change"])
  liverprophenodf[ourid,"lactate_change_dueto_train"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___lactate_change_dueto_train"])
  liverprophenodf[ourid,"vo2_max_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___vo_2_max_change"])
  liverprophenodf[ourid,"sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourpid,"pass1bf0027"][1]]
  liverprophenodf[ourid,"group"] <- pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourpid,"pass1bf0011"][1]
}
trimliverprophenodf <- liverprophenodf[liverprophenodf$group %in% c("Four-week program",
                                                                    "Eight-week program Training Group",
                                                                    "Eight-week program Control Group"),]
trimliverprophenodf[trimliverprophenodf$group %in% "Eight-week program Control Group","group"] <- "Control"
trimliverprophenodf[trimliverprophenodf$group %in% "Eight-week program Training Group","group"] <- "8w"
trimliverprophenodf[trimliverprophenodf$group %in% "Four-week program","group"] <- "4w"
trimliverprophenodf$group <- factor(trimliverprophenodf$group,levels = c("Control","4w","8w"))


liverprosigvsphenocor <- matrix(0L,nrow = length(liverprosig),ncol = dim(liverprophenomeasuredata)[2])
rownames(liverprosigvsphenocor) <- liverprosig
colnames(liverprosigvsphenocor) <- colnames(liverprophenomeasuredata)
for(i in 1:length(liverprosig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:dim(liverprophenomeasuredata)[2]){
    liverprosigvsphenocor[i,j] <- cor(t(liverprosignorm[i,rownames(trimliverprophenodf)]),liverprophenomeasuredata[rownames(trimliverprophenodf),j])
  }
}
colnames(liverprosigvsphenocor) <- gsub("calculated_variables___","",colnames(liverprosigvsphenocor))


liverprotfnorm <- liverpronorm[intersect(rownames(liverpronorm),tfproanno$Liver.Pro.ID),]
liverprotfnorm <- liverprotfnorm[rowSums(is.na(liverprotfnorm)) == 0,]

liverprotfvsphenocor <- matrix(0L,nrow = dim(liverprotfnorm)[1],ncol = dim(liverprophenomeasuredata)[2])
rownames(liverprotfvsphenocor) <- rownames(liverprotfnorm)
colnames(liverprotfvsphenocor) <- colnames(liverprophenomeasuredata)
for(i in 1:length(rownames(liverprotfnorm))){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:dim(liverprophenomeasuredata)[2]){
    liverprotfvsphenocor[i,j] <- cor(t(liverprotfnorm[i,rownames(trimliverprophenodf)]),liverprophenomeasuredata[rownames(trimliverprophenodf),j])
  }
}
colnames(liverprotfvsphenocor) <- gsub("calculated_variables___","",colnames(liverprotfvsphenocor))


# LUNG

lungprosig <- unique(prot_pr$training_dea[prot_pr$training_dea$adj_p_value < 0.1 & prot_pr$training_dea$tissue_abbreviation %in% "LUNG","feature_ID"])

lungpronorm <- read.delim(file = "November 2021 PASS1B Data Freeze/Proteomics Data/Prot_PR/pass1b-06_v1.0_analysis_proteomics-untargeted_prot-pr_normalized-data_motrpac_pass1b-06_t66-lung_p.txt",header = T,row.names = 1)
colnames(lungpronorm) <- gsub("X","",colnames(lungpronorm))
lungprosignorm <- lungpronorm[lungprosig,]

lungprosignorm <- lungprosignorm[rowSums(is.na(lungprosignorm)) == 0,]
lungprosig <- rownames(lungprosignorm)



lungprophenomeasuredata <- phenomeasuredata[intersect(gsub("X","",colnames(lungpronorm)),rownames(phenomeasuredata)),]
lungprophenomeasuredata <- lungprophenomeasuredata[,c("calculated_variables___wgt_gain_after_train","calculated_variables___pct_body_fat_change","calculated_variables___pct_body_lean_change","calculated_variables___pct_body_fluid_change","calculated_variables___lactate_change_dueto_train","calculated_variables___vo_2_max_change")]

lungprometa <- data.frame(row.names = colnames(lungprosignorm),
                          "label" = colnames(lungprosignorm))


lungprophenodf <- data.frame(row.names = colnames(lungprosignorm),"label" = colnames(lungprosignorm))
lungprophenodf$wgt_gain_after_train <- 0
lungprophenodf$pct_body_fat_change <- 0
lungprophenodf$pct_body_lean_change <- 0
lungprophenodf$pct_body_fluid_change <- 0
lungprophenodf$lactate_change_dueto_train <- 0
lungprophenodf$vo2_max_change <- 0
lungprophenodf$sex <- ""
lungprophenodf$group <- ""
for(i in 1:dim(lungprophenodf)[1]){
  ourid <- rownames(lungprophenodf)[i]
  ourpid <- phenomeasuredata[ourid,"pid"]
  lungprophenodf[ourid,"wgt_gain_after_train"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___wgt_gain_after_train"])
  lungprophenodf[ourid,"pct_body_fat_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_fat_change"])
  lungprophenodf[ourid,"pct_body_lean_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_lean_change"])
  lungprophenodf[ourid,"pct_body_fluid_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_fluid_change"])
  lungprophenodf[ourid,"lactate_change_dueto_train"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___lactate_change_dueto_train"])
  lungprophenodf[ourid,"vo2_max_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___vo_2_max_change"])
  lungprophenodf[ourid,"sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourpid,"pass1bf0027"][1]]
  lungprophenodf[ourid,"group"] <- pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourpid,"pass1bf0011"][1]
}
trimlungprophenodf <- lungprophenodf[lungprophenodf$group %in% c("Four-week program",
                                                                 "Eight-week program Training Group",
                                                                 "Eight-week program Control Group"),]
trimlungprophenodf[trimlungprophenodf$group %in% "Eight-week program Control Group","group"] <- "Control"
trimlungprophenodf[trimlungprophenodf$group %in% "Eight-week program Training Group","group"] <- "8w"
trimlungprophenodf[trimlungprophenodf$group %in% "Four-week program","group"] <- "4w"
trimlungprophenodf$group <- factor(trimlungprophenodf$group,levels = c("Control","4w","8w"))


lungprosigvsphenocor <- matrix(0L,nrow = length(lungprosig),ncol = dim(lungprophenomeasuredata)[2])
rownames(lungprosigvsphenocor) <- lungprosig
colnames(lungprosigvsphenocor) <- colnames(lungprophenomeasuredata)
for(i in 1:length(lungprosig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:dim(lungprophenomeasuredata)[2]){
    lungprosigvsphenocor[i,j] <- cor(t(lungprosignorm[i,rownames(trimlungprophenodf)]),lungprophenomeasuredata[rownames(trimlungprophenodf),j])
  }
}
colnames(lungprosigvsphenocor) <- gsub("calculated_variables___","",colnames(lungprosigvsphenocor))


lungprotfnorm <- lungpronorm[intersect(rownames(lungpronorm),tfproanno$Lung.Pro.ID),]
lungprotfnorm <- lungprotfnorm[rowSums(is.na(lungprotfnorm)) == 0,]

lungprotfvsphenocor <- matrix(0L,nrow = dim(lungprotfnorm)[1],ncol = dim(lungprophenomeasuredata)[2])
rownames(lungprotfvsphenocor) <- rownames(lungprotfnorm)
colnames(lungprotfvsphenocor) <- colnames(lungprophenomeasuredata)
for(i in 1:length(rownames(lungprotfnorm))){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:dim(lungprophenomeasuredata)[2]){
    lungprotfvsphenocor[i,j] <- cor(t(lungprotfnorm[i,rownames(trimlungprophenodf)]),lungprophenomeasuredata[rownames(trimlungprophenodf),j])
  }
}
colnames(lungprotfvsphenocor) <- gsub("calculated_variables___","",colnames(lungprotfvsphenocor))


# WAT-SC

whiteprosig <- unique(prot_pr$training_dea[prot_pr$training_dea$adj_p_value < 0.1 & prot_pr$training_dea$tissue_abbreviation %in% "WAT-SC","feature_ID"])

whitepronorm <- read.delim(file = "November 2021 PASS1B Data Freeze/Proteomics Data/Prot_PR/pass1b-06_v1.1_analysis_proteomics-untargeted_prot-pr_normalized-data_motrpac_pass1b-06_t70-white-adip.txt",header = T,row.names = 1)
colnames(whitepronorm) <- gsub("X","",colnames(whitepronorm))
whiteprosignorm <- whitepronorm[whiteprosig,]

whiteprosignorm <- whiteprosignorm[rowSums(is.na(whiteprosignorm)) == 0,]
whiteprosig <- rownames(whiteprosignorm)



whiteprophenomeasuredata <- phenomeasuredata[intersect(gsub("X","",colnames(whitepronorm)),rownames(phenomeasuredata)),]
whiteprophenomeasuredata <- whiteprophenomeasuredata[,c("calculated_variables___wgt_gain_after_train","calculated_variables___pct_body_fat_change","calculated_variables___pct_body_lean_change","calculated_variables___pct_body_fluid_change","calculated_variables___lactate_change_dueto_train","calculated_variables___vo_2_max_change")]

whiteprometa <- data.frame(row.names = colnames(whiteprosignorm),
                           "label" = colnames(whiteprosignorm))


whiteprophenodf <- data.frame(row.names = colnames(whiteprosignorm),"label" = colnames(whiteprosignorm))
whiteprophenodf$wgt_gain_after_train <- 0
whiteprophenodf$pct_body_fat_change <- 0
whiteprophenodf$pct_body_lean_change <- 0
whiteprophenodf$pct_body_fluid_change <- 0
whiteprophenodf$lactate_change_dueto_train <- 0
whiteprophenodf$vo2_max_change <- 0
whiteprophenodf$sex <- ""
whiteprophenodf$group <- ""
for(i in 1:dim(whiteprophenodf)[1]){
  ourid <- rownames(whiteprophenodf)[i]
  ourpid <- phenomeasuredata[ourid,"pid"]
  whiteprophenodf[ourid,"wgt_gain_after_train"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___wgt_gain_after_train"])
  whiteprophenodf[ourid,"pct_body_fat_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_fat_change"])
  whiteprophenodf[ourid,"pct_body_lean_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_lean_change"])
  whiteprophenodf[ourid,"pct_body_fluid_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___pct_body_fluid_change"])
  whiteprophenodf[ourid,"lactate_change_dueto_train"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___lactate_change_dueto_train"])
  whiteprophenodf[ourid,"vo2_max_change"] <- mean(phenomeasuredata[rownames(phenomeasuredata) %in% ourid,"calculated_variables___vo_2_max_change"])
  whiteprophenodf[ourid,"sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourpid,"pass1bf0027"][1]]
  whiteprophenodf[ourid,"group"] <- pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourpid,"pass1bf0011"][1]
}
trimwhiteprophenodf <- whiteprophenodf[whiteprophenodf$group %in% c("Four-week program",
                                                                    "Eight-week program Training Group",
                                                                    "Eight-week program Control Group"),]
trimwhiteprophenodf[trimwhiteprophenodf$group %in% "Eight-week program Control Group","group"] <- "Control"
trimwhiteprophenodf[trimwhiteprophenodf$group %in% "Eight-week program Training Group","group"] <- "8w"
trimwhiteprophenodf[trimwhiteprophenodf$group %in% "Four-week program","group"] <- "4w"
trimwhiteprophenodf$group <- factor(trimwhiteprophenodf$group,levels = c("Control","4w","8w"))


whiteprosigvsphenocor <- matrix(0L,nrow = length(whiteprosig),ncol = dim(whiteprophenomeasuredata)[2])
rownames(whiteprosigvsphenocor) <- whiteprosig
colnames(whiteprosigvsphenocor) <- colnames(whiteprophenomeasuredata)
for(i in 1:length(whiteprosig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:dim(whiteprophenomeasuredata)[2]){
    whiteprosigvsphenocor[i,j] <- cor(t(whiteprosignorm[i,rownames(trimwhiteprophenodf)]),whiteprophenomeasuredata[rownames(trimwhiteprophenodf),j])
  }
}
colnames(whiteprosigvsphenocor) <- gsub("calculated_variables___","",colnames(whiteprosigvsphenocor))


whiteprotfnorm <- whitepronorm[intersect(rownames(whitepronorm),tfproanno$WhiteAd.Pro.ID),]
whiteprotfnorm <- whiteprotfnorm[rowSums(is.na(whiteprotfnorm)) == 0,]

whiteprotfvsphenocor <- matrix(0L,nrow = dim(whiteprotfnorm)[1],ncol = dim(whiteprophenomeasuredata)[2])
rownames(whiteprotfvsphenocor) <- rownames(whiteprotfnorm)
colnames(whiteprotfvsphenocor) <- colnames(whiteprophenomeasuredata)
for(i in 1:length(rownames(whiteprotfnorm))){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:dim(whiteprophenomeasuredata)[2]){
    whiteprotfvsphenocor[i,j] <- cor(t(whiteprotfnorm[i,rownames(trimwhiteprophenodf)]),whiteprophenomeasuredata[rownames(trimwhiteprophenodf),j])
  }
}
colnames(whiteprotfvsphenocor) <- gsub("calculated_variables___","",colnames(whiteprotfvsphenocor))


gastroprotfvsphenocorsub <- rownames(gastroprotfvsphenocor)[apply(abs(gastroprotfvsphenocor),1,max) > 0.5]
heartprotfvsphenocorsub <- rownames(heartprotfvsphenocor)[apply(abs(heartprotfvsphenocor),1,max) > 0.5]
kidneyprotfvsphenocorsub <- rownames(kidneyprotfvsphenocor)[apply(abs(kidneyprotfvsphenocor),1,max) > 0.5]
liverprotfvsphenocorsub <- rownames(liverprotfvsphenocor)[apply(abs(liverprotfvsphenocor),1,max) > 0.5]
lungprotfvsphenocorsub <- rownames(lungprotfvsphenocor)[apply(abs(lungprotfvsphenocor),1,max) > 0.5]
whiteprotfvsphenocorsub <- rownames(whiteprotfvsphenocor)[apply(abs(whiteprotfvsphenocor),1,max) > 0.5]

comboprotfvsphenocor <- rbind(gastroprotfvsphenocor[gastroprotfvsphenocorsub,],
                              heartprotfvsphenocor[heartprotfvsphenocorsub,],
                              kidneyprotfvsphenocor[kidneyprotfvsphenocorsub,],
                              liverprotfvsphenocor[liverprotfvsphenocorsub,],
                              lungprotfvsphenocor[lungprotfvsphenocorsub,],
                              whiteprotfvsphenocor[whiteprotfvsphenocorsub,])
rownames(comboprotfvsphenocor)[19] <- rownames(lungprotfvsphenocor)[apply(abs(lungprotfvsphenocor),1,max) > 0.5]
combotfproidlist <- c(gastroprotfvsphenocorsub,
                      heartprotfvsphenocorsub,
                      kidneyprotfvsphenocorsub,
                      liverprotfvsphenocorsub,
                      lungprotfvsphenocorsub,
                      whiteprotfvsphenocorsub)
rownames(comboprotfvsphenocor)[1:2] <- paste("gastro_",rownames(comboprotfvsphenocor)[1:2],sep = "")
rownames(comboprotfvsphenocor)[3:6] <- paste("heart",rownames(comboprotfvsphenocor)[3:6],sep = "")
#rownames(comboprotfvsphenocor)[7] <- paste("kidney_",rownames(comboprotfvsphenocor)[7],sep = "")
rownames(comboprotfvsphenocor)[7:18] <- paste("liver_",rownames(comboprotfvsphenocor)[7:18],sep = "")
rownames(comboprotfvsphenocor)[19] <- paste("lung_",rownames(comboprotfvsphenocor)[19],sep = "")
rownames(comboprotfvsphenocor)[20:31] <- paste("white_",rownames(comboprotfvsphenocor)[20:31],sep = "")
comboprotfvsphenometa <- data.frame(row.names = rownames(comboprotfvsphenocor),
                                    "Tissue" = c(rep("SKM-GN",2),
                                                 rep("HEART",4),
                                                 rep("LIVER",12),
                                                 rep("LUNG",1),
                                                 rep("WAT-SC",12)))

combotfpronamelist <- combotfproidlist
combotfpronamelist[1] <- rownames(tfproanno)[tfproanno$Gastro.Pro.ID %in% combotfproidlist[1]][1]
combotfpronamelist[2] <- rownames(tfproanno)[tfproanno$Gastro.Pro.ID %in% combotfproidlist[2]][1]
combotfpronamelist[3] <- rownames(tfproanno)[tfproanno$Heart.Pro.ID %in% combotfproidlist[3]][1]
combotfpronamelist[4] <- rownames(tfproanno)[tfproanno$Heart.Pro.ID %in% combotfproidlist[4]][1]
combotfpronamelist[5] <- rownames(tfproanno)[tfproanno$Heart.Pro.ID %in% combotfproidlist[5]][1]
combotfpronamelist[6] <- rownames(tfproanno)[tfproanno$Heart.Pro.ID %in% combotfproidlist[6]][1]
combotfpronamelist[7] <- rownames(tfproanno)[tfproanno$Liver.Pro.ID %in% combotfproidlist[7]][1]
combotfpronamelist[8] <- rownames(tfproanno)[tfproanno$Liver.Pro.ID %in% combotfproidlist[8]][1]
combotfpronamelist[9] <- rownames(tfproanno)[tfproanno$Liver.Pro.ID %in% combotfproidlist[9]][1]
combotfpronamelist[10] <- rownames(tfproanno)[tfproanno$Liver.Pro.ID %in% combotfproidlist[10]][1]
combotfpronamelist[11] <- rownames(tfproanno)[tfproanno$Liver.Pro.ID %in% combotfproidlist[11]][1]
combotfpronamelist[12] <- rownames(tfproanno)[tfproanno$Liver.Pro.ID %in% combotfproidlist[12]][1]
combotfpronamelist[13] <- rownames(tfproanno)[tfproanno$Liver.Pro.ID %in% combotfproidlist[13]][1]
combotfpronamelist[14] <- rownames(tfproanno)[tfproanno$Liver.Pro.ID %in% combotfproidlist[14]][1]
combotfpronamelist[15] <- rownames(tfproanno)[tfproanno$Liver.Pro.ID %in% combotfproidlist[15]][1]
combotfpronamelist[16] <- rownames(tfproanno)[tfproanno$Liver.Pro.ID %in% combotfproidlist[16]][1]
combotfpronamelist[17] <- rownames(tfproanno)[tfproanno$Liver.Pro.ID %in% combotfproidlist[17]][1]
combotfpronamelist[18] <- rownames(tfproanno)[tfproanno$Liver.Pro.ID %in% combotfproidlist[18]][1]
combotfpronamelist[19] <- rownames(tfproanno)[tfproanno$Lung.Pro.ID %in% combotfproidlist[19]][1]
combotfpronamelist[20] <- rownames(tfproanno)[tfproanno$WhiteAd.Pro.ID %in% combotfproidlist[20]][1]
combotfpronamelist[21] <- rownames(tfproanno)[tfproanno$WhiteAd.Pro.ID %in% combotfproidlist[21]][1]
combotfpronamelist[22] <- rownames(tfproanno)[tfproanno$WhiteAd.Pro.ID %in% combotfproidlist[22]][1]
combotfpronamelist[23] <- rownames(tfproanno)[tfproanno$WhiteAd.Pro.ID %in% combotfproidlist[23]][1]
combotfpronamelist[24] <- rownames(tfproanno)[tfproanno$WhiteAd.Pro.ID %in% combotfproidlist[24]][1]
combotfpronamelist[25] <- rownames(tfproanno)[tfproanno$WhiteAd.Pro.ID %in% combotfproidlist[25]][1]
combotfpronamelist[26] <- rownames(tfproanno)[tfproanno$WhiteAd.Pro.ID %in% combotfproidlist[26]][1]
combotfpronamelist[27] <- rownames(tfproanno)[tfproanno$WhiteAd.Pro.ID %in% combotfproidlist[27]][1]
combotfpronamelist[28] <- rownames(tfproanno)[tfproanno$WhiteAd.Pro.ID %in% combotfproidlist[28]][1]
combotfpronamelist[29] <- rownames(tfproanno)[tfproanno$WhiteAd.Pro.ID %in% combotfproidlist[29]][1]
combotfpronamelist[30] <- rownames(tfproanno)[tfproanno$WhiteAd.Pro.ID %in% combotfproidlist[30]][1]
combotfpronamelist[31] <- rownames(tfproanno)[tfproanno$WhiteAd.Pro.ID %in% combotfproidlist[31]][1]

# Supplemental Figure S16B
pdf(file = "Supplemental Figure S16B.pdf",height = 7,width = 8)
pheatmap(comboprotfvsphenocor,cluster_rows = T,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),annotation_row = comboprotfvsphenometa,annotation_col = sigvsphenocormeta,annotation_colors = ann_cols_sigphenocor,show_colnames = F,display_numbers = T,number_color = "black",labels_row = gsub("\\(.*","",combotfpronamelist),fontsize = 15)
dev.off()
