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

######
# This file involves loading and processing necessary data into usable form for 
# downstream analysis and figure generation
# - file names used below have been simplified from MoTrPAC data releases
####

# To be edited to your directory containing the folder with all necessary paper data files
setwd("C:/Users/gsmit/Documents/Mount Sinai/Sealfon Laboratory/MoTrPAC/PASS1B Raw Data")

## Collecting phenotypic data and sample metadata
pass1bphenodata <- read.delim(file = "PASS1B Transcription Factor Paper Data/Phenotype Data/motrpac_pass1b-06_pheno_viallabel-data.txt",header = T,sep = "\t")
phenomeasuredata <- read.delim(file = "PASS1B Transcription Factor Paper Data/Phenotype Data/pass1b-06_phenotype_motrpac_pass1b-06_pheno_viallabel_data_merged_v2.0.txt",header = T,sep = "\t",row.names = 4)
# To correct for a missing entry
phenomeasuredata["90412017702",c("calculated_variables___wgt_gain_after_train","calculated_variables___pct_body_fat_change","calculated_variables___pct_body_lean_change","calculated_variables___pct_body_fluid_change","calculated_variables___lactate_change_dueto_train","calculated_variables___vo_2_max_change")] <- phenomeasuredata["90412017014",c("calculated_variables___wgt_gain_after_train","calculated_variables___pct_body_fat_change","calculated_variables___pct_body_lean_change","calculated_variables___pct_body_fluid_change","calculated_variables___lactate_change_dueto_train","calculated_variables___vo_2_max_change")]

saveRDS(pass1bphenodata,file = "pass1bphenodata.rds")
saveRDS(phenomeasuredata,file = "phenomeasuredata.rds")

# Collecting ATACseq peak annotations

peakanno <- read.delim(file = "PASS1B Transcription Factor Paper Data/ATACSeq Data/pass1b-06_epigen-atac-seq_feature-mapping.txt",header = T,row.names = 2)
peakanno$trimanno <- gsub(" \\(.*","",peakanno$custom_annotation)

saveRDS(peakanno,file = "peakanno.RDS")

newmethanno <- read.delim(file = "PASS1B Transcription Factor Paper Data/RRBS Data/pass1b-06_epigen-rrbs-feature-mapping.txt",header = T)
trimmedmethanno <- newmethanno[!duplicated(newmethanno$feature_ID),]
rownames(trimmedmethanno) <- trimmedmethanno$feature_ID
saveRDS(trimmedmethanno,file = "trimmedmethanno.RDS")

#####
# Identifying the significant analytes in each ome
####

load(file = "PASS1B Transcription Factor Paper Data/DEA Analysis/pass1b-06_transcript-rna-seq_dea.RData")
load(file = "PASS1B Transcription Factor Paper Data/DEA Analysis/pass1b-06_epigen-atac-seq_dea.RData")
load(file = "PASS1B Transcription Factor Paper Data/DEA Analysis/pass1b-06_proteomics-prot-pr_dea.RData")
load(file = "PASS1B Transcription Factor Paper Data/DEA Analysis/pass1b-06_proteomics-prot-ph_dea.RData")

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

gastroprosig <- unique(prot_pr$training_dea[prot_pr$training_dea$tissue_abbreviation %in% "SKM-GN" &
                                              prot_pr$training_dea$adj_p_value < 0.1,"feature_ID"])
heartprosig <- unique(prot_pr$training_dea[prot_pr$training_dea$tissue_abbreviation %in% "HEART" &
                                             prot_pr$training_dea$adj_p_value < 0.1,"feature_ID"])
kidneyprosig <- unique(prot_pr$training_dea[prot_pr$training_dea$tissue_abbreviation %in% "KIDNEY" &
                                              prot_pr$training_dea$adj_p_value < 0.1,"feature_ID"])
liverprosig <- unique(prot_pr$training_dea[prot_pr$training_dea$tissue_abbreviation %in% "LIVER" &
                                             prot_pr$training_dea$adj_p_value < 0.1,"feature_ID"])
lungprosig <- unique(prot_pr$training_dea[prot_pr$training_dea$tissue_abbreviation %in% "LUNG" &
                                            prot_pr$training_dea$adj_p_value < 0.1,"feature_ID"])
whiteprosig <- unique(prot_pr$training_dea[prot_pr$training_dea$tissue_abbreviation %in% "WAT-SC" &
                                             prot_pr$training_dea$adj_p_value < 0.1,"feature_ID"])


gastroprophsig <- unique(prot_ph$training_dea[prot_ph$training_dea$tissue_abbreviation %in% "SKM-GN" &
                                             prot_ph$training_dea$adj_p_value < 0.1,"feature_ID"])
heartprophsig <- unique(prot_ph$training_dea[prot_ph$training_dea$tissue_abbreviation %in% "HEART" &
                                            prot_ph$training_dea$adj_p_value < 0.1,"feature_ID"])
kidneyprophsig <- unique(prot_ph$training_dea[prot_ph$training_dea$tissue_abbreviation %in% "KIDNEY" &
                                             prot_ph$training_dea$adj_p_value < 0.1,"feature_ID"])
liverprophsig <- unique(prot_ph$training_dea[prot_ph$training_dea$tissue_abbreviation %in% "LIVER" &
                                            prot_ph$training_dea$adj_p_value < 0.1,"feature_ID"])
lungprophsig <- unique(prot_ph$training_dea[prot_ph$training_dea$tissue_abbreviation %in% "LUNG" &
                                           prot_ph$training_dea$adj_p_value < 0.1,"feature_ID"])
whiteprophsig <- unique(prot_ph$training_dea[prot_ph$training_dea$tissue_abbreviation %in% "WAT-SC" &
                                            prot_ph$training_dea$adj_p_value < 0.1,"feature_ID"])
rm(epigen_atac_seq)
rm(transcript_rna_seq)
rm(prot_ph)
rm(prot_pr)
gc()

load(file = "PASS1B Transcription Factor Paper Data/DEA Analysis/pass1b-06_epigen-rrbs_dea.RData")

gastromethsig <- epigen_rrbs$training_dea[epigen_rrbs$training_dea$adj_p_value < 0.1 & epigen_rrbs$training_dea$tissue_abbreviation %in% "SKM-GN","feature_ID"]
heartmethsig <- epigen_rrbs$training_dea[epigen_rrbs$training_dea$adj_p_value < 0.1 & epigen_rrbs$training_dea$tissue_abbreviation %in% "HEART","feature_ID"]
hippomethsig <- epigen_rrbs$training_dea[epigen_rrbs$training_dea$adj_p_value < 0.1 & epigen_rrbs$training_dea$tissue_abbreviation %in% "HIPPOC","feature_ID"]
kidneymethsig <- epigen_rrbs$training_dea[epigen_rrbs$training_dea$adj_p_value < 0.1 & epigen_rrbs$training_dea$tissue_abbreviation %in% "KIDNEY","feature_ID"]
livermethsig <- epigen_rrbs$training_dea[epigen_rrbs$training_dea$adj_p_value < 0.1 & epigen_rrbs$training_dea$tissue_abbreviation %in% "LIVER","feature_ID"]
lungmethsig <- epigen_rrbs$training_dea[epigen_rrbs$training_dea$adj_p_value < 0.1 & epigen_rrbs$training_dea$tissue_abbreviation %in% "LUNG","feature_ID"]
brownmethsig <- epigen_rrbs$training_dea[epigen_rrbs$training_dea$adj_p_value < 0.1 & epigen_rrbs$training_dea$tissue_abbreviation %in% "BAT","feature_ID"]
whitemethsig <- epigen_rrbs$training_dea[epigen_rrbs$training_dea$adj_p_value < 0.1 & epigen_rrbs$training_dea$tissue_abbreviation %in% "WAT-SC","feature_ID"]

save(list = c("gastrornasig","heartrnasig","hippornasig","kidneyrnasig",
              "liverrnasig","lungrnasig","brownrnasig","whiternasig",
              "gastroatacsig","heartatacsig","hippoatacsig","kidneyatacsig",
              "liveratacsig","lungatacsig","brownatacsig","whiteatacsig",
              "gastromethsig","heartmethsig","hippomethsig","kidneymethsig",
              "livermethsig","lungmethsig","brownmethsig","whitemethsig",
              "gastroprosig","heartprosig","kidneyprosig",
              "liverprosig","lungprosig","whiteprosig",
              "gastroprophsig","heartprophsig","kidneyprophsig",
                "liverprophsig","lungprophsig","whiteprophsig"),file = "omesigdata.RData")

rm(epigen_rrbs)


#####
# Generating L2FC matrices for the significant analytes in each ome
####

load(file = "omesigdata.RData")

load(file = "PASS1B Transcription Factor Paper Data/DEA Analysis/pass1b-06_transcript-rna-seq_dea.RData")
load(file = "PASS1B Transcription Factor Paper Data/DEA Analysis/pass1b-06_epigen-atac-seq_dea.RData")

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

save(list = c("gastrol2fcmat","heartl2fcmat","hippol2fcmat","kidneyl2fcmat",
              "liverl2fcmat","lungl2fcmat","brownl2fcmat","whitel2fcmat"),
     file = "rnal2fcmat.RData")

gastrosigrnal2fc <- gastrol2fcmat[gastrornasig,]
heartsigrnal2fc <- heartl2fcmat[heartrnasig,]
hipposigrnal2fc <- hippol2fcmat[hippornasig,]
kidneysigrnal2fc <- kidneyl2fcmat[kidneyrnasig,]
liversigrnal2fc <- liverl2fcmat[liverrnasig,]
lungsigrnal2fc <- lungl2fcmat[lungrnasig,]
brownsigrnal2fc <- brownl2fcmat[brownrnasig,]
whitesigrnal2fc <- whitel2fcmat[whiternasig,]

save(list = c("gastrosigrnal2fc","heartsigrnal2fc","hipposigrnal2fc","kidneysigrnal2fc",
              "liversigrnal2fc","lungsigrnal2fc","brownsigrnal2fc","whitesigrnal2fc"),
     file = "rnasigl2fcmat.RData")


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


save(list = c("gastrosigatacl2fc","heartsigatacl2fc","hipposigatacl2fc","kidneysigatacl2fc",
              "liversigatacl2fc","lungsigatacl2fc","brownsigatacl2fc","whitesigatacl2fc"),
     file = "atacsigl2fcmat.RData")

rm(epigen_atac_seq)
rm(transcript_rna_seq)
gc()

load(file = "PASS1B Transcription Factor Paper Data/DEA Analysis/pass1b-06_epigen-rrbs_dea.RData")

gastromethsigl2fc <- matrix(0L,nrow = length(gastromethsig),ncol = 8)
rownames(gastromethsigl2fc) <- gastromethsig
colnames(gastromethsigl2fc) <- c("F W1","F W2","F W4","F W8",
                              "M W1","M W2","M W4","M W8")
gastromethtimewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "SKM-GN",]
for(i in 1:length(gastromethsig)){
  ourmeth <- gastromethsig[i]
  ourmethtimewise <- gastromethtimewise[gastromethtimewise$tissue_abbreviation %in% "SKM-GN" & gastromethtimewise$feature_ID %in% ourmeth,]
  gastromethsigl2fc[i,"F W1"] <- ourmethtimewise[ourmethtimewise$sex %in% "female" & ourmethtimewise$comparison_group %in% "1w","logFC"]
  gastromethsigl2fc[i,"F W2"] <- ourmethtimewise[ourmethtimewise$sex %in% "female" & ourmethtimewise$comparison_group %in% "2w","logFC"]
  gastromethsigl2fc[i,"F W4"] <- ourmethtimewise[ourmethtimewise$sex %in% "female" & ourmethtimewise$comparison_group %in% "4w","logFC"]
  gastromethsigl2fc[i,"F W8"] <- ourmethtimewise[ourmethtimewise$sex %in% "female" & ourmethtimewise$comparison_group %in% "8w","logFC"]
  gastromethsigl2fc[i,"M W1"] <- ourmethtimewise[ourmethtimewise$sex %in% "male" & ourmethtimewise$comparison_group %in% "1w","logFC"]
  gastromethsigl2fc[i,"M W2"] <- ourmethtimewise[ourmethtimewise$sex %in% "male" & ourmethtimewise$comparison_group %in% "2w","logFC"]
  gastromethsigl2fc[i,"M W4"] <- ourmethtimewise[ourmethtimewise$sex %in% "male" & ourmethtimewise$comparison_group %in% "4w","logFC"]
  gastromethsigl2fc[i,"M W8"] <- ourmethtimewise[ourmethtimewise$sex %in% "male" & ourmethtimewise$comparison_group %in% "8w","logFC"]
}

heartmethsigl2fc <- matrix(0L,nrow = length(heartmethsig),ncol = 8)
rownames(heartmethsigl2fc) <- heartmethsig
colnames(heartmethsigl2fc) <- c("F W1","F W2","F W4","F W8",
                              "M W1","M W2","M W4","M W8")
heartmethtimewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "HEART",]
for(i in 1:length(heartmethsig)){
  ourmeth <- heartmethsig[i]
  ourmethtimewise <- heartmethtimewise[heartmethtimewise$tissue_abbreviation %in% "HEART" & heartmethtimewise$feature_ID %in% ourmeth,]
  heartmethsigl2fc[i,"F W1"] <- ourmethtimewise[ourmethtimewise$sex %in% "female" & ourmethtimewise$comparison_group %in% "1w","logFC"]
  heartmethsigl2fc[i,"F W2"] <- ourmethtimewise[ourmethtimewise$sex %in% "female" & ourmethtimewise$comparison_group %in% "2w","logFC"]
  heartmethsigl2fc[i,"F W4"] <- ourmethtimewise[ourmethtimewise$sex %in% "female" & ourmethtimewise$comparison_group %in% "4w","logFC"]
  heartmethsigl2fc[i,"F W8"] <- ourmethtimewise[ourmethtimewise$sex %in% "female" & ourmethtimewise$comparison_group %in% "8w","logFC"]
  heartmethsigl2fc[i,"M W1"] <- ourmethtimewise[ourmethtimewise$sex %in% "male" & ourmethtimewise$comparison_group %in% "1w","logFC"]
  heartmethsigl2fc[i,"M W2"] <- ourmethtimewise[ourmethtimewise$sex %in% "male" & ourmethtimewise$comparison_group %in% "2w","logFC"]
  heartmethsigl2fc[i,"M W4"] <- ourmethtimewise[ourmethtimewise$sex %in% "male" & ourmethtimewise$comparison_group %in% "4w","logFC"]
  heartmethsigl2fc[i,"M W8"] <- ourmethtimewise[ourmethtimewise$sex %in% "male" & ourmethtimewise$comparison_group %in% "8w","logFC"]
}

hippomethsigl2fc <- matrix(0L,nrow = length(hippomethsig),ncol = 8)
rownames(hippomethsigl2fc) <- hippomethsig
colnames(hippomethsigl2fc) <- c("F W1","F W2","F W4","F W8",
                              "M W1","M W2","M W4","M W8")
hippomethtimewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "HIPPOC",]
for(i in 1:length(hippomethsig)){
  ourmeth <- hippomethsig[i]
  ourmethtimewise <- hippomethtimewise[hippomethtimewise$tissue_abbreviation %in% "HIPPOC" & hippomethtimewise$feature_ID %in% ourmeth,]
  hippomethsigl2fc[i,"F W1"] <- ourmethtimewise[ourmethtimewise$sex %in% "female" & ourmethtimewise$comparison_group %in% "1w","logFC"]
  hippomethsigl2fc[i,"F W2"] <- ourmethtimewise[ourmethtimewise$sex %in% "female" & ourmethtimewise$comparison_group %in% "2w","logFC"]
  hippomethsigl2fc[i,"F W4"] <- ourmethtimewise[ourmethtimewise$sex %in% "female" & ourmethtimewise$comparison_group %in% "4w","logFC"]
  hippomethsigl2fc[i,"F W8"] <- ourmethtimewise[ourmethtimewise$sex %in% "female" & ourmethtimewise$comparison_group %in% "8w","logFC"]
  hippomethsigl2fc[i,"M W1"] <- ourmethtimewise[ourmethtimewise$sex %in% "male" & ourmethtimewise$comparison_group %in% "1w","logFC"]
  hippomethsigl2fc[i,"M W2"] <- ourmethtimewise[ourmethtimewise$sex %in% "male" & ourmethtimewise$comparison_group %in% "2w","logFC"]
  hippomethsigl2fc[i,"M W4"] <- ourmethtimewise[ourmethtimewise$sex %in% "male" & ourmethtimewise$comparison_group %in% "4w","logFC"]
  hippomethsigl2fc[i,"M W8"] <- ourmethtimewise[ourmethtimewise$sex %in% "male" & ourmethtimewise$comparison_group %in% "8w","logFC"]
}

kidneymethsigl2fc <- matrix(0L,nrow = length(kidneymethsig),ncol = 8)
rownames(kidneymethsigl2fc) <- kidneymethsig
colnames(kidneymethsigl2fc) <- c("F W1","F W2","F W4","F W8",
                              "M W1","M W2","M W4","M W8")
kidneymethtimewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "KIDNEY",]
for(i in 1:length(kidneymethsig)){
  ourmeth <- kidneymethsig[i]
  ourmethtimewise <- kidneymethtimewise[kidneymethtimewise$tissue_abbreviation %in% "KIDNEY" & kidneymethtimewise$feature_ID %in% ourmeth,]
  kidneymethsigl2fc[i,"F W1"] <- ourmethtimewise[ourmethtimewise$sex %in% "female" & ourmethtimewise$comparison_group %in% "1w","logFC"]
  kidneymethsigl2fc[i,"F W2"] <- ourmethtimewise[ourmethtimewise$sex %in% "female" & ourmethtimewise$comparison_group %in% "2w","logFC"]
  kidneymethsigl2fc[i,"F W4"] <- ourmethtimewise[ourmethtimewise$sex %in% "female" & ourmethtimewise$comparison_group %in% "4w","logFC"]
  kidneymethsigl2fc[i,"F W8"] <- ourmethtimewise[ourmethtimewise$sex %in% "female" & ourmethtimewise$comparison_group %in% "8w","logFC"]
  kidneymethsigl2fc[i,"M W1"] <- ourmethtimewise[ourmethtimewise$sex %in% "male" & ourmethtimewise$comparison_group %in% "1w","logFC"]
  kidneymethsigl2fc[i,"M W2"] <- ourmethtimewise[ourmethtimewise$sex %in% "male" & ourmethtimewise$comparison_group %in% "2w","logFC"]
  kidneymethsigl2fc[i,"M W4"] <- ourmethtimewise[ourmethtimewise$sex %in% "male" & ourmethtimewise$comparison_group %in% "4w","logFC"]
  kidneymethsigl2fc[i,"M W8"] <- ourmethtimewise[ourmethtimewise$sex %in% "male" & ourmethtimewise$comparison_group %in% "8w","logFC"]
}

livermethsigl2fc <- matrix(0L,nrow = length(livermethsig),ncol = 8)
rownames(livermethsigl2fc) <- livermethsig
colnames(livermethsigl2fc) <- c("F W1","F W2","F W4","F W8",
                              "M W1","M W2","M W4","M W8")
livermethtimewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "LIVER",]
for(i in 1:length(livermethsig)){
  ourmeth <- livermethsig[i]
  ourmethtimewise <- livermethtimewise[livermethtimewise$tissue_abbreviation %in% "LIVER" & livermethtimewise$feature_ID %in% ourmeth,]
  livermethsigl2fc[i,"F W1"] <- ourmethtimewise[ourmethtimewise$sex %in% "female" & ourmethtimewise$comparison_group %in% "1w","logFC"]
  livermethsigl2fc[i,"F W2"] <- ourmethtimewise[ourmethtimewise$sex %in% "female" & ourmethtimewise$comparison_group %in% "2w","logFC"]
  livermethsigl2fc[i,"F W4"] <- ourmethtimewise[ourmethtimewise$sex %in% "female" & ourmethtimewise$comparison_group %in% "4w","logFC"]
  livermethsigl2fc[i,"F W8"] <- ourmethtimewise[ourmethtimewise$sex %in% "female" & ourmethtimewise$comparison_group %in% "8w","logFC"]
  livermethsigl2fc[i,"M W1"] <- ourmethtimewise[ourmethtimewise$sex %in% "male" & ourmethtimewise$comparison_group %in% "1w","logFC"]
  livermethsigl2fc[i,"M W2"] <- ourmethtimewise[ourmethtimewise$sex %in% "male" & ourmethtimewise$comparison_group %in% "2w","logFC"]
  livermethsigl2fc[i,"M W4"] <- ourmethtimewise[ourmethtimewise$sex %in% "male" & ourmethtimewise$comparison_group %in% "4w","logFC"]
  livermethsigl2fc[i,"M W8"] <- ourmethtimewise[ourmethtimewise$sex %in% "male" & ourmethtimewise$comparison_group %in% "8w","logFC"]
}

lungmethsigl2fc <- matrix(0L,nrow = length(lungmethsig),ncol = 8)
rownames(lungmethsigl2fc) <- lungmethsig
colnames(lungmethsigl2fc) <- c("F W1","F W2","F W4","F W8",
                              "M W1","M W2","M W4","M W8")
lungmethtimewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "LUNG",]
for(i in 1:length(lungmethsig)){
  ourmeth <- lungmethsig[i]
  ourmethtimewise <- lungmethtimewise[lungmethtimewise$tissue_abbreviation %in% "LUNG" & lungmethtimewise$feature_ID %in% ourmeth,]
  lungmethsigl2fc[i,"F W1"] <- ourmethtimewise[ourmethtimewise$sex %in% "female" & ourmethtimewise$comparison_group %in% "1w","logFC"]
  lungmethsigl2fc[i,"F W2"] <- ourmethtimewise[ourmethtimewise$sex %in% "female" & ourmethtimewise$comparison_group %in% "2w","logFC"]
  lungmethsigl2fc[i,"F W4"] <- ourmethtimewise[ourmethtimewise$sex %in% "female" & ourmethtimewise$comparison_group %in% "4w","logFC"]
  lungmethsigl2fc[i,"F W8"] <- ourmethtimewise[ourmethtimewise$sex %in% "female" & ourmethtimewise$comparison_group %in% "8w","logFC"]
  lungmethsigl2fc[i,"M W1"] <- ourmethtimewise[ourmethtimewise$sex %in% "male" & ourmethtimewise$comparison_group %in% "1w","logFC"]
  lungmethsigl2fc[i,"M W2"] <- ourmethtimewise[ourmethtimewise$sex %in% "male" & ourmethtimewise$comparison_group %in% "2w","logFC"]
  lungmethsigl2fc[i,"M W4"] <- ourmethtimewise[ourmethtimewise$sex %in% "male" & ourmethtimewise$comparison_group %in% "4w","logFC"]
  lungmethsigl2fc[i,"M W8"] <- ourmethtimewise[ourmethtimewise$sex %in% "male" & ourmethtimewise$comparison_group %in% "8w","logFC"]
}

brownmethsigl2fc <- matrix(0L,nrow = length(brownmethsig),ncol = 8)
rownames(brownmethsigl2fc) <- brownmethsig
colnames(brownmethsigl2fc) <- c("F W1","F W2","F W4","F W8",
                              "M W1","M W2","M W4","M W8")
brownmethtimewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "BAT",]
for(i in 1:length(brownmethsig)){
  ourmeth <- brownmethsig[i]
  ourmethtimewise <- brownmethtimewise[brownmethtimewise$tissue_abbreviation %in% "BAT" & brownmethtimewise$feature_ID %in% ourmeth,]
  brownmethsigl2fc[i,"F W1"] <- ourmethtimewise[ourmethtimewise$sex %in% "female" & ourmethtimewise$comparison_group %in% "1w","logFC"]
  brownmethsigl2fc[i,"F W2"] <- ourmethtimewise[ourmethtimewise$sex %in% "female" & ourmethtimewise$comparison_group %in% "2w","logFC"]
  brownmethsigl2fc[i,"F W4"] <- ourmethtimewise[ourmethtimewise$sex %in% "female" & ourmethtimewise$comparison_group %in% "4w","logFC"]
  brownmethsigl2fc[i,"F W8"] <- ourmethtimewise[ourmethtimewise$sex %in% "female" & ourmethtimewise$comparison_group %in% "8w","logFC"]
  brownmethsigl2fc[i,"M W1"] <- ourmethtimewise[ourmethtimewise$sex %in% "male" & ourmethtimewise$comparison_group %in% "1w","logFC"]
  brownmethsigl2fc[i,"M W2"] <- ourmethtimewise[ourmethtimewise$sex %in% "male" & ourmethtimewise$comparison_group %in% "2w","logFC"]
  brownmethsigl2fc[i,"M W4"] <- ourmethtimewise[ourmethtimewise$sex %in% "male" & ourmethtimewise$comparison_group %in% "4w","logFC"]
  brownmethsigl2fc[i,"M W8"] <- ourmethtimewise[ourmethtimewise$sex %in% "male" & ourmethtimewise$comparison_group %in% "8w","logFC"]
}

whitemethsigl2fc <- matrix(0L,nrow = length(whitemethsig),ncol = 8)
rownames(whitemethsigl2fc) <- whitemethsig
colnames(whitemethsigl2fc) <- c("F W1","F W2","F W4","F W8",
                              "M W1","M W2","M W4","M W8")
whitemethtimewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "WAT-SC",]
for(i in 1:length(whitemethsig)){
  ourmeth <- whitemethsig[i]
  ourmethtimewise <- whitemethtimewise[whitemethtimewise$tissue_abbreviation %in% "WAT-SC" & whitemethtimewise$feature_ID %in% ourmeth,]
  whitemethsigl2fc[i,"F W1"] <- ourmethtimewise[ourmethtimewise$sex %in% "female" & ourmethtimewise$comparison_group %in% "1w","logFC"]
  whitemethsigl2fc[i,"F W2"] <- ourmethtimewise[ourmethtimewise$sex %in% "female" & ourmethtimewise$comparison_group %in% "2w","logFC"]
  whitemethsigl2fc[i,"F W4"] <- ourmethtimewise[ourmethtimewise$sex %in% "female" & ourmethtimewise$comparison_group %in% "4w","logFC"]
  whitemethsigl2fc[i,"F W8"] <- ourmethtimewise[ourmethtimewise$sex %in% "female" & ourmethtimewise$comparison_group %in% "8w","logFC"]
  whitemethsigl2fc[i,"M W1"] <- ourmethtimewise[ourmethtimewise$sex %in% "male" & ourmethtimewise$comparison_group %in% "1w","logFC"]
  whitemethsigl2fc[i,"M W2"] <- ourmethtimewise[ourmethtimewise$sex %in% "male" & ourmethtimewise$comparison_group %in% "2w","logFC"]
  whitemethsigl2fc[i,"M W4"] <- ourmethtimewise[ourmethtimewise$sex %in% "male" & ourmethtimewise$comparison_group %in% "4w","logFC"]
  whitemethsigl2fc[i,"M W8"] <- ourmethtimewise[ourmethtimewise$sex %in% "male" & ourmethtimewise$comparison_group %in% "8w","logFC"]
}

save(list = c("gastromethsigl2fc","heartmethsigl2fc","hippomethsigl2fc","kidneymethsigl2fc",
              "livermethsigl2fc","lungmethsigl2fc","brownmethsigl2fc","whitemethsigl2fc"),
     file = "new_methsigl2fcmat_51624.RData")

#####
# Setting up normalized data matrices
####

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

gastrornanorm <- gastrornanorm[3:dim(gastrornanorm)[1],]
heartrnanorm <- heartrnanorm[3:dim(heartrnanorm)[1],]
hippornanorm <- hippornanorm[3:dim(hippornanorm)[1],]
kidneyrnanorm <- kidneyrnanorm[3:dim(kidneyrnanorm)[1],]
liverrnanorm <- liverrnanorm[3:dim(liverrnanorm)[1],]
lungrnanorm <- lungrnanorm[3:dim(lungrnanorm)[1],]
whiternanorm <- whiternanorm[3:dim(whiternanorm)[1],]
brownrnanorm <- brownrnanorm[3:dim(brownrnanorm)[1],]

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

save(list = c("gastrornanorm","heartrnanorm","hippornanorm","kidneyrnanorm",
              "liverrnanorm","lungrnanorm","brownrnanorm","whiternanorm"),
     file = "rnanormmatrices.RData")

save(list = c("gastroatacnorm","heartatacnorm","hippoatacnorm","kidneyatacnorm",
              "liveratacnorm","lungatacnorm","brownatacnorm","whiteatacnorm"),
     file = "atacnormmatrices.RData")


gastroMmat <- read.delim(file = "PASS1B Transcription Factor Paper Data/RRBS Data/gastrocnemius_M_matrix.txt",header = T,row.names = 1)
heartMmat <- read.delim(file = "PASS1B Transcription Factor Paper Data/RRBS Data/heart_M_matrix.txt",header = T,row.names = 1)
hippoMmat <- read.delim(file = "PASS1B Transcription Factor Paper Data/RRBS Data/hippocampus_M_matrix.txt",header = T,row.names = 1)
kidneyMmat <- read.delim(file = "PASS1B Transcription Factor Paper Data/RRBS Data/kidney_M_matrix.txt",header = T,row.names = 1)
liverMmat <- read.delim(file = "PASS1B Transcription Factor Paper Data/RRBS Data/liver_M_matrix.txt",header = T,row.names = 1)
lungMmat <- read.delim(file = "PASS1B Transcription Factor Paper Data/RRBS Data/lung_M_matrix.txt",header = T,row.names = 1)
brownMmat <- read.delim(file = "PASS1B Transcription Factor Paper Data/RRBS Data/brown-adipose_M_matrix.txt",header = T,row.names = 1)
whiteMmat <- read.delim(file = "PASS1B Transcription Factor Paper Data/RRBS Data/white-adipose_M_matrix.txt",header = T,row.names = 1)

gastroMmat <- gastroMmat[,grep("FALSE",gastroMmat["is_outlier",])]
heartMmat <- heartMmat[,grep("FALSE",heartMmat["is_outlier",])]
hippoMmat <- hippoMmat[,grep("FALSE",hippoMmat["is_outlier",])]
kidneyMmat <- kidneyMmat[,grep("FALSE",kidneyMmat["is_outlier",])]
liverMmat <- liverMmat[,grep("FALSE",liverMmat["is_outlier",])]
lungMmat <- lungMmat[,grep("FALSE",lungMmat["is_outlier",])]
brownMmat <- brownMmat[,grep("FALSE",brownMmat["is_outlier",])]
whiteMmat <- whiteMmat[,grep("FALSE",whiteMmat["is_outlier",])]

colnames(gastroMmat) <- gastroMmat["pid",]
gastroMmat <- gastroMmat[4:dim(gastroMmat)[1],]

colnames(heartMmat) <- heartMmat["pid",]
heartMmat <- heartMmat[4:dim(heartMmat)[1],]

colnames(hippoMmat) <- hippoMmat["pid",]
hippoMmat <- hippoMmat[4:dim(hippoMmat)[1],]

colnames(kidneyMmat) <- kidneyMmat["pid",]
kidneyMmat <- kidneyMmat[4:dim(kidneyMmat)[1],]

colnames(liverMmat) <- liverMmat["pid",]
liverMmat <- liverMmat[4:dim(liverMmat)[1],]

colnames(lungMmat) <- lungMmat["pid",]
lungMmat <- lungMmat[4:dim(lungMmat)[1],]

colnames(brownMmat) <- brownMmat["pid",]
brownMmat <- brownMmat[4:dim(brownMmat)[1],]

colnames(whiteMmat) <- whiteMmat["pid",]
whiteMmat <- whiteMmat[4:dim(whiteMmat)[1],]

save(list = c("gastroMmat","heartMmat","hippoMmat","kidneyMmat",
              "liverMmat","lungMmat","brownMmat","whiteMmat"),
     file = "methMmatrices.RData")

gastroproprnorm <-read.delim(file = "PASS1B Transcription Factor Paper Data/Prot_PR Data/pass1b-06_proteomics-prot-pr_normalized-data_t55-gastrocnemius.txt",sep = "\t",header = T,row.names = 1)
heartproprnorm <-read.delim(file = "PASS1B Transcription Factor Paper Data/Prot_PR Data/pass1b-06_proteomics-prot-pr_normalized-data_t58-heart.txt",sep = "\t",header = T,row.names = 1)
kidneyproprnorm <-read.delim(file = "PASS1B Transcription Factor Paper Data/Prot_PR Data/pass1b-06_proteomics-prot-pr_normalized-data_t59-kidney.txt",sep = "\t",header = T,row.names = 1)
liverproprnorm <-read.delim(file = "PASS1B Transcription Factor Paper Data/Prot_PR Data/pass1b-06_proteomics-prot-pr_normalized-data_t68-liver.txt",sep = "\t",header = T,row.names = 1)
lungproprnorm <-read.delim(file = "PASS1B Transcription Factor Paper Data/Prot_PR Data/pass1b-06_proteomics-prot-pr_normalized-data_t66-lung.txt",sep = "\t",header = T,row.names = 1)
whiteproprnorm <-read.delim(file = "PASS1B Transcription Factor Paper Data/Prot_PR Data/pass1b-06_proteomics-prot-pr_normalized-data_t70-white-adipose.txt",sep = "\t",header = T,row.names = 1)

colnames(gastroproprnorm) <- gastroproprnorm["pid",]
gastroproprnorm <- gastroproprnorm[3:dim(gastroproprnorm)[1],]
gastroproprnorm[is.na(gastroproprnorm)] <- min(gastroproprnorm[!is.na(gastroproprnorm)])

colnames(heartproprnorm) <- heartproprnorm["pid",]
heartproprnorm <- heartproprnorm[3:dim(heartproprnorm)[1],]
heartproprnorm[is.na(heartproprnorm)] <- min(heartproprnorm[!is.na(heartproprnorm)])

colnames(kidneyproprnorm) <- kidneyproprnorm["pid",]
kidneyproprnorm <- kidneyproprnorm[3:dim(kidneyproprnorm)[1],]
kidneyproprnorm[is.na(kidneyproprnorm)] <- min(kidneyproprnorm[!is.na(kidneyproprnorm)])

colnames(liverproprnorm) <- liverproprnorm["pid",]
liverproprnorm <- liverproprnorm[3:dim(liverproprnorm)[1],]
liverproprnorm[is.na(liverproprnorm)] <- min(liverproprnorm[!is.na(liverproprnorm)])

colnames(lungproprnorm) <- lungproprnorm["pid",]
lungproprnorm <- lungproprnorm[3:dim(lungproprnorm)[1],]
lungproprnorm[is.na(lungproprnorm)] <- min(lungproprnorm[!is.na(lungproprnorm)])

colnames(whiteproprnorm) <- whiteproprnorm["pid",]
whiteproprnorm <- whiteproprnorm[3:dim(whiteproprnorm)[1],]
whiteproprnorm[is.na(whiteproprnorm)] <- min(whiteproprnorm[!is.na(whiteproprnorm)])

save(list = c("gastroproprnorm","heartproprnorm","kidneyproprnorm","liverproprnorm",
              "lungproprnorm","whiteproprnorm"),
     file = "proprnormmatrices.RData")


#### Now phosphoproteome data

gastroprophnormdata <- MotrpacRatTraining6moData::PHOSPHO_SKMGN_NORM_DATA
rownames(gastroprophnormdata) <- gastroprophnormdata$feature_ID
gastroprophnorm <- gastroprophnormdata[,c(5:dim(gastroprophnormdata)[2])]
gastroprophnorm[is.na(gastroprophnorm)] <- min(gastroprophnorm[!is.na(gastroprophnorm)])

heartprophnormdata <- MotrpacRatTraining6moData::PHOSPHO_HEART_NORM_DATA
rownames(heartprophnormdata) <- heartprophnormdata$feature_ID
heartprophnorm <- heartprophnormdata[,c(5:dim(heartprophnormdata)[2])]
heartprophnorm[is.na(heartprophnorm)] <- min(heartprophnorm[!is.na(heartprophnorm)])

kidneyprophnormdata <- MotrpacRatTraining6moData::PHOSPHO_KIDNEY_NORM_DATA
rownames(kidneyprophnormdata) <- kidneyprophnormdata$feature_ID
kidneyprophnorm <- kidneyprophnormdata[,c(5:dim(kidneyprophnormdata)[2])]
kidneyprophnorm[is.na(kidneyprophnorm)] <- min(kidneyprophnorm[!is.na(kidneyprophnorm)])

liverprophnormdata <- MotrpacRatTraining6moData::PHOSPHO_LIVER_NORM_DATA
rownames(liverprophnormdata) <- liverprophnormdata$feature_ID
liverprophnorm <- liverprophnormdata[,c(5:dim(liverprophnormdata)[2])]
liverprophnorm[is.na(liverprophnorm)] <- min(liverprophnorm[!is.na(liverprophnorm)])

lungprophnormdata <- MotrpacRatTraining6moData::PHOSPHO_LUNG_NORM_DATA
rownames(lungprophnormdata) <- lungprophnormdata$feature_ID
lungprophnorm <- lungprophnormdata[,c(5:dim(lungprophnormdata)[2])]
lungprophnorm[is.na(lungprophnorm)] <- min(lungprophnorm[!is.na(lungprophnorm)])

whiteprophnormdata <- MotrpacRatTraining6moData::PHOSPHO_WATSC_NORM_DATA
rownames(whiteprophnormdata) <- whiteprophnormdata$feature_ID
whiteprophnorm <- whiteprophnormdata[,c(5:dim(whiteprophnormdata)[2])]
whiteprophnorm[is.na(whiteprophnorm)] <- min(whiteprophnorm[!is.na(whiteprophnorm)])

save(list = c("gastroprophnorm","heartprophnorm","kidneyprophnorm","liverprophnorm",
              "lungprophnorm","whiteprophnorm"),
     file = "prophnormmatrices.RData")

#####
# Setting up gene ensembl to symbol annotations
####

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
enstosym["ENSRNOG00000046050","Symbol"] <- "Dennd1c"
enstosym["ENSRNOG00000053448","Symbol"] <- "G6pc1"

save(list = c("enstosym"),
     file = "enstosym.RData")

#####
# Generate RNA normalized control matrix
####

pass1bphenodata <- readRDS("pass1bphenodata.rds")

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

save(list = c("rnacontrolnorm","rnacontrolmeta","controlcolumns"),
     file = "rnacontrolnormandmeta.RData")

#####
# Generating ATACseq Normalized Control Matrix
####

#atacseqcontrolnorm <- read.csv(file = "atacseqnorm1922.csv",header = T,row.names = 1)

pass1bphenodata <- readRDS("pass1bphenodata.rds")

gastroatacdata <- read.delim(file = "PASS1B Transcription Factor Paper Data/ATACSeq Data/T55.atac.counts.txt",header = T,sep = "\t")
heartatacdata <- read.delim(file = "PASS1B Transcription Factor Paper Data/ATACSeq Data/T58.atac.counts.txt",header = T,sep = "\t")
hippoatacdata <- read.delim(file = "PASS1B Transcription Factor Paper Data/ATACSeq Data/T52.atac.counts.txt",header = T,sep = "\t")
kidneyatacdata <- read.delim(file = "PASS1B Transcription Factor Paper Data/ATACSeq Data/T59.atac.counts.txt",header = T,sep = "\t")
liveratacdata <- read.delim(file = "PASS1B Transcription Factor Paper Data/ATACSeq Data/T68.atac.counts.txt",header = T,sep = "\t")
lungatacdata <- read.delim(file = "PASS1B Transcription Factor Paper Data/ATACSeq Data/T66.atac.counts.txt",header = T,sep = "\t")
whiteatacdata <- read.delim(file = "PASS1B Transcription Factor Paper Data/ATACSeq Data/T70.atac.counts.txt",header = T,sep = "\t")
brownatacdata <- read.delim(file = "PASS1B Transcription Factor Paper Data/ATACSeq Data/T69.atac.counts.txt",header = T,sep = "\t")

# Gastro - we use location 830 = stanford
gastropeakmeta <- gastroatacdata[,c("chrom","start","end")]
rownames(gastropeakmeta) <- paste(gastropeakmeta$chrom,":",gastropeakmeta$start,"-",gastropeakmeta$end,sep = "")
rownames(gastroatacdata) <- rownames(gastropeakmeta)
colnames(gastroatacdata) <- gsub("X","",colnames(gastroatacdata))

#gastroatacdata <- gastroatacdata[,colnames(gastroatacdata) %in% pass1bphenodata$pass1bf0004]
gastroatacdata <- gastroatacdata[,colnames(gastroatacdata) %in% pass1bphenodata[pass1bphenodata$pass1bf0542 %in% "830","pass1bf0004"]]

for(i in 1:50){
  ataccol <- colnames(gastroatacdata)[i]
  colnames(gastroatacdata)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 == ataccol,"pass1bf0001"]
}

# White Adipose

whitepeakmeta <- whiteatacdata[,c("chrom","start","end")]
rownames(whitepeakmeta) <- paste(whitepeakmeta$chrom,":",whitepeakmeta$start,"-",whitepeakmeta$end)
rownames(whiteatacdata) <- rownames(whitepeakmeta)
colnames(whiteatacdata) <- gsub("X","",colnames(whiteatacdata))

#whiteatacdata <- whiteatacdata[,colnames(whiteatacdata) %in% pass1bphenodata$pass1bf0004]
whiteatacdata <- whiteatacdata[,colnames(whiteatacdata) %in% pass1bphenodata[pass1bphenodata$pass1bf0542 %in% "830","pass1bf0004"]]

for(i in 1:50){
  ataccol <- colnames(whiteatacdata)[i]
  colnames(whiteatacdata)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 == ataccol,"pass1bf0001"]
}

# Heart

heartpeakmeta <- heartatacdata[,c("chrom","start","end")]
rownames(heartpeakmeta) <- paste(heartpeakmeta$chrom,":",heartpeakmeta$start,"-",heartpeakmeta$end)
rownames(heartatacdata) <- rownames(heartpeakmeta)
colnames(heartatacdata) <- gsub("X","",colnames(heartatacdata))

heartatacdata <- heartatacdata[,colnames(heartatacdata) %in% pass1bphenodata$pass1bf0004]
#heartatacdata <- heartatacdata[,colnames(heartatacdata) %in% pass1bphenodata[pass1bphenodata$pass1bf0542 %in% "830","pass1bf0004"]]

for(i in 1:50){
  ataccol <- colnames(heartatacdata)[i]
  colnames(heartatacdata)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 == ataccol,"pass1bf0001"]
}

# Hippocampus

hippopeakmeta <- hippoatacdata[,c("chrom","start","end")]
rownames(hippopeakmeta) <- paste(hippopeakmeta$chrom,":",hippopeakmeta$start,"-",hippopeakmeta$end)
rownames(hippoatacdata) <- rownames(hippopeakmeta)
colnames(hippoatacdata) <- gsub("X","",colnames(hippoatacdata))

hippoatacdata <- hippoatacdata[,colnames(hippoatacdata) %in% pass1bphenodata$pass1bf0004]
#hippoatacdata <- hippoatacdata[,colnames(hippoatacdata) %in% pass1bphenodata[pass1bphenodata$pass1bf0542 %in% "830","pass1bf0004"]]

for(i in 1:50){
  ataccol <- colnames(hippoatacdata)[i]
  colnames(hippoatacdata)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 == ataccol,"pass1bf0001"]
}

# Kidney

kidneypeakmeta <- kidneyatacdata[,c("chrom","start","end")]
rownames(kidneypeakmeta) <- paste(kidneypeakmeta$chrom,":",kidneypeakmeta$start,"-",kidneypeakmeta$end)
rownames(kidneyatacdata) <- rownames(kidneypeakmeta)
colnames(kidneyatacdata) <- gsub("X","",colnames(kidneyatacdata))

kidneyatacdata <- kidneyatacdata[,colnames(kidneyatacdata) %in% pass1bphenodata$pass1bf0004]
#kidneyatacdata <- kidneyatacdata[,colnames(kidneyatacdata) %in% pass1bphenodata[pass1bphenodata$pass1bf0542 %in% "830","pass1bf0004"]]

for(i in 1:50){
  ataccol <- colnames(kidneyatacdata)[i]
  colnames(kidneyatacdata)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 == ataccol,"pass1bf0001"]
}

# Liver

liverpeakmeta <- liveratacdata[,c("chrom","start","end")]
rownames(liverpeakmeta) <- paste(liverpeakmeta$chrom,":",liverpeakmeta$start,"-",liverpeakmeta$end)
rownames(liveratacdata) <- rownames(liverpeakmeta)
colnames(liveratacdata) <- gsub("X","",colnames(liveratacdata))

liveratacdata <- liveratacdata[,colnames(liveratacdata) %in% pass1bphenodata$pass1bf0004]
#liveratacdata <- liveratacdata[,colnames(liveratacdata) %in% pass1bphenodata[pass1bphenodata$pass1bf0542 %in% "830","pass1bf0004"]]

for(i in 1:50){
  ataccol <- colnames(liveratacdata)[i]
  colnames(liveratacdata)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 == ataccol,"pass1bf0001"]
}

# Lung

lungpeakmeta <- lungatacdata[,c("chrom","start","end")]
rownames(lungpeakmeta) <- paste(lungpeakmeta$chrom,":",lungpeakmeta$start,"-",lungpeakmeta$end)
rownames(lungatacdata) <- rownames(lungpeakmeta)
colnames(lungatacdata) <- gsub("X","",colnames(lungatacdata))

lungatacdata <- lungatacdata[,colnames(lungatacdata) %in% pass1bphenodata$pass1bf0004]
#lungatacdata <- lungatacdata[,colnames(lungatacdata) %in% pass1bphenodata[pass1bphenodata$pass1bf0542 %in% "830","pass1bf0004"]]

for(i in 1:50){
  ataccol <- colnames(lungatacdata)[i]
  colnames(lungatacdata)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 == ataccol,"pass1bf0001"]
}

# Brown Adipose

brownpeakmeta <- brownatacdata[,c("chrom","start","end")]
rownames(brownpeakmeta) <- paste(brownpeakmeta$chrom,":",brownpeakmeta$start,"-",brownpeakmeta$end)
rownames(brownatacdata) <- rownames(brownpeakmeta)
colnames(brownatacdata) <- gsub("X","",colnames(brownatacdata))

brownatacdata <- brownatacdata[,colnames(brownatacdata) %in% pass1bphenodata$pass1bf0004]
#brownatacdata <- brownatacdata[,colnames(brownatacdata) %in% pass1bphenodata[pass1bphenodata$pass1bf0542 %in% "830","pass1bf0004"]]

for(i in 1:50){
  ataccol <- colnames(brownatacdata)[i]
  colnames(brownatacdata)[i] <- pass1bphenodata[pass1bphenodata$pass1bf0004 == ataccol,"pass1bf0001"]
}

oursamplelist <- colnames(gastroatacdata)
controlsamples <- as.character(intersect(oursamplelist,unique(pass1bphenodata[pass1bphenodata$pass1bf0011 %in% "Eight-week program Control Group","pass1bf0001"])))
controlsamplemeta <- data.frame(row.names = controlsamples,"sex" = rep("Female",10))
for(i in 1:10){
  oursample <- controlsamples[i]
  pass1bphenodata[pass1bphenodata$pass1bf0001 %in% oursample,"pass1bf0027"]
  if(unique(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% oursample,"pass1bf0027"]) == 2){
    controlsamplemeta[oursample,"sex"] <- "Male"
  }
}
sortcontrol <- rownames(controlsamplemeta)[order(controlsamplemeta$sex)]

atacdata <- cbind(brownatacdata[,sortcontrol],gastroatacdata[,sortcontrol],heartatacdata[,sortcontrol],hippoatacdata[,sortcontrol],
                  kidneyatacdata[,sortcontrol],liveratacdata[,sortcontrol],lungatacdata[,sortcontrol],whiteatacdata[,sortcontrol])
colnames(atacdata) <- c(paste("BrownAd_",sortcontrol,sep = ""),
                        paste("Gastro_",sortcontrol,sep = ""),
                        paste("Heart_",sortcontrol,sep = ""),
                        paste("Hippo_",sortcontrol,sep = ""),
                        paste("Kidney_",sortcontrol,sep = ""),
                        paste("Liver_",sortcontrol,sep = ""),
                        paste("Lung_",sortcontrol,sep = ""),
                        paste("WhiteAd_",sortcontrol,sep = ""))

#site_counts = counts[,site_meta[,viallabel]]

# remove non-auto peaks
atacdata = atacdata[grepl("^chr[0-9]|^chrY|^chrX", rownames(atacdata)),]

# exclude low count peaks in the current dataset
# at least 10 counts in N samples
n_samples = 4
atacfiltdata = atacdata[rowSums(data.frame(lapply(atacdata, function(x) as.numeric(x >= 10)), check.names=F)) >= n_samples,]

# quantile normalize
# this takes a couple of minutes given the size of the peak x sample counts matrix
atacnorm = voom(atacfiltdata,normalize.method = "quantile")
atacsubnorm = round(atacnorm$E,2)

saveRDS(atacsubnorm,file = "atacseqcontrolnorm.RDS")

#####
# Peak Motif Connection
####

allpeakmotifs <- read.table(file = "PASS1B Transcription Factor Paper Data/TF Data/allpeaksmotifs.txt",header = T,sep = "\t")
saveRDS(allpeakmotifs,file = "allpeakmotifs.RDS")

#####
# Specify active peaks in each tissue for downstream analysis
####

load("atacnormmatrices.RData")

activepeakcutoff <- -1

gastroactivepeaks <- rownames(gastroatacnorm[apply(gastroatacnorm,1,median) > activepeakcutoff,])
heartactivepeaks <- rownames(heartatacnorm[apply(heartatacnorm,1,median) > activepeakcutoff,])
hippoactivepeaks <- rownames(hippoatacnorm[apply(hippoatacnorm,1,median) > activepeakcutoff,])
kidneyactivepeaks <- rownames(kidneyatacnorm[apply(kidneyatacnorm,1,median) > activepeakcutoff,])
liveractivepeaks <- rownames(liveratacnorm[apply(liveratacnorm,1,median) > activepeakcutoff,])
lungactivepeaks <- rownames(lungatacnorm[apply(lungatacnorm,1,median) > activepeakcutoff,])
brownactivepeaks <- rownames(brownatacnorm[apply(brownatacnorm,1,median) > activepeakcutoff,])
whiteactivepeaks <- rownames(whiteatacnorm[apply(whiteatacnorm,1,median) > activepeakcutoff,])

save(list = c("gastroactivepeaks","heartactivepeaks","hippoactivepeaks","kidneyactivepeaks",
              "liveractivepeaks","lungactivepeaks","brownactivepeaks","whiteactivepeaks"),
     file = "activepeakfiles.RData")


#####
# TF Annotations
####

tfanno <- read.csv("PASS1B Transcription Factor Paper Data/TF Data/tflist.csv",header = T,row.names = 1)
tfanno$Ensembl <- gsub(" ","",tfanno$Ensembl)
tfanno <- tfanno[!tfanno$Ensembl == "",]
saveRDS(tfanno,file = "tfanno.RDS")

#####
# TF Protein Annotations
####

tfanno <- readRDS("tfanno.RDS")

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

saveRDS(tfproanno,file = "tfproanno.RDS")

#####
# Protein L2FC Matrices
####

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

save(list = c("gastroprol2fc","heartprol2fc","kidneyprol2fc","liverprol2fc","lungprol2fc","whiteprol2fc"),
     file = "prol2fcmat.RData")


#####
# Phospho-prophteome L2FC Matrices
####

# Building prophteomics dataset

gastroproph <- prot_ph$timewise_dea[prot_ph$timewise_dea$tissue %in% "t55-gastrocnemius",]
heartproph <- prot_ph$timewise_dea[prot_ph$timewise_dea$tissue %in% "t58-heart",]
kidneyproph <- prot_ph$timewise_dea[prot_ph$timewise_dea$tissue %in% "t59-kidney",]
liverproph <- prot_ph$timewise_dea[prot_ph$timewise_dea$tissue %in% "t68-liver",]
lungproph <- prot_ph$timewise_dea[prot_ph$timewise_dea$tissue %in% "t66-lung",]
whiteproph <- prot_ph$timewise_dea[prot_ph$timewise_dea$tissue %in% "t70-white-adipose",]

gastroprophl2fc <- matrix(0L,nrow = length(unique(gastroproph$feature_ID)),ncol = 8)
rownames(gastroprophl2fc) <- unique(gastroproph$feature_ID)
colnames(gastroprophl2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
for(i in 1:dim(gastroprophl2fc)[1]){
  
  if(i%%100 == 0){
    print(paste("Gastro:",toString(i)))
  }
  
  ourproph <- rownames(gastroprophl2fc)[i]
  ourprophmat <- gastroproph[gastroproph$feature_ID %in% ourproph,]
  gastroprophl2fc[i,"F W1"] <- ourprophmat[ourprophmat$sex %in% "female" & ourprophmat$comparison_group %in% "1w","logFC"]
  gastroprophl2fc[i,"F W2"] <- ourprophmat[ourprophmat$sex %in% "female" & ourprophmat$comparison_group %in% "2w","logFC"]
  gastroprophl2fc[i,"F W4"] <- ourprophmat[ourprophmat$sex %in% "female" & ourprophmat$comparison_group %in% "4w","logFC"]
  gastroprophl2fc[i,"F W8"] <- ourprophmat[ourprophmat$sex %in% "female" & ourprophmat$comparison_group %in% "8w","logFC"]
  
  gastroprophl2fc[i,"M W1"] <- ourprophmat[ourprophmat$sex %in% "male" & ourprophmat$comparison_group %in% "1w","logFC"]
  gastroprophl2fc[i,"M W2"] <- ourprophmat[ourprophmat$sex %in% "male" & ourprophmat$comparison_group %in% "2w","logFC"]
  gastroprophl2fc[i,"M W4"] <- ourprophmat[ourprophmat$sex %in% "male" & ourprophmat$comparison_group %in% "4w","logFC"]
  gastroprophl2fc[i,"M W8"] <- ourprophmat[ourprophmat$sex %in% "male" & ourprophmat$comparison_group %in% "8w","logFC"]
}

heartprophl2fc <- matrix(0L,nrow = length(unique(heartproph$feature_ID)),ncol = 8)
rownames(heartprophl2fc) <- unique(heartproph$feature_ID)
colnames(heartprophl2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
for(i in 1:dim(heartprophl2fc)[1]){
  
  if(i%%100 == 0){
    print(paste("Heart:",toString(i)))
  }
  
  ourproph <- rownames(heartprophl2fc)[i]
  ourprophmat <- heartproph[heartproph$feature_ID %in% ourproph,]
  heartprophl2fc[i,"F W1"] <- ourprophmat[ourprophmat$sex %in% "female" & ourprophmat$comparison_group %in% "1w","logFC"]
  heartprophl2fc[i,"F W2"] <- ourprophmat[ourprophmat$sex %in% "female" & ourprophmat$comparison_group %in% "2w","logFC"]
  heartprophl2fc[i,"F W4"] <- ourprophmat[ourprophmat$sex %in% "female" & ourprophmat$comparison_group %in% "4w","logFC"]
  heartprophl2fc[i,"F W8"] <- ourprophmat[ourprophmat$sex %in% "female" & ourprophmat$comparison_group %in% "8w","logFC"]
  
  heartprophl2fc[i,"M W1"] <- ourprophmat[ourprophmat$sex %in% "male" & ourprophmat$comparison_group %in% "1w","logFC"]
  heartprophl2fc[i,"M W2"] <- ourprophmat[ourprophmat$sex %in% "male" & ourprophmat$comparison_group %in% "2w","logFC"]
  heartprophl2fc[i,"M W4"] <- ourprophmat[ourprophmat$sex %in% "male" & ourprophmat$comparison_group %in% "4w","logFC"]
  heartprophl2fc[i,"M W8"] <- ourprophmat[ourprophmat$sex %in% "male" & ourprophmat$comparison_group %in% "8w","logFC"]
}

kidneyprophl2fc <- matrix(0L,nrow = length(unique(kidneyproph$feature_ID)),ncol = 8)
rownames(kidneyprophl2fc) <- unique(kidneyproph$feature_ID)
colnames(kidneyprophl2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
for(i in 1:dim(kidneyprophl2fc)[1]){
  
  if(i%%100 == 0){
    print(paste("Kidney:",toString(i)))
  }
  
  ourproph <- rownames(kidneyprophl2fc)[i]
  ourprophmat <- kidneyproph[kidneyproph$feature_ID %in% ourproph,]
  kidneyprophl2fc[i,"F W1"] <- ourprophmat[ourprophmat$sex %in% "female" & ourprophmat$comparison_group %in% "1w","logFC"]
  kidneyprophl2fc[i,"F W2"] <- ourprophmat[ourprophmat$sex %in% "female" & ourprophmat$comparison_group %in% "2w","logFC"]
  kidneyprophl2fc[i,"F W4"] <- ourprophmat[ourprophmat$sex %in% "female" & ourprophmat$comparison_group %in% "4w","logFC"]
  kidneyprophl2fc[i,"F W8"] <- ourprophmat[ourprophmat$sex %in% "female" & ourprophmat$comparison_group %in% "8w","logFC"]
  
  kidneyprophl2fc[i,"M W1"] <- ourprophmat[ourprophmat$sex %in% "male" & ourprophmat$comparison_group %in% "1w","logFC"]
  kidneyprophl2fc[i,"M W2"] <- ourprophmat[ourprophmat$sex %in% "male" & ourprophmat$comparison_group %in% "2w","logFC"]
  kidneyprophl2fc[i,"M W4"] <- ourprophmat[ourprophmat$sex %in% "male" & ourprophmat$comparison_group %in% "4w","logFC"]
  kidneyprophl2fc[i,"M W8"] <- ourprophmat[ourprophmat$sex %in% "male" & ourprophmat$comparison_group %in% "8w","logFC"]
}

liverprophl2fc <- matrix(0L,nrow = length(unique(liverproph$feature_ID)),ncol = 8)
rownames(liverprophl2fc) <- unique(liverproph$feature_ID)
colnames(liverprophl2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
for(i in 1:dim(liverprophl2fc)[1]){
  
  if(i%%100 == 0){
    print(paste("Liver:",toString(i)))
  }
  
  ourproph <- rownames(liverprophl2fc)[i]
  ourprophmat <- liverproph[liverproph$feature_ID %in% ourproph,]
  liverprophl2fc[i,"F W1"] <- ourprophmat[ourprophmat$sex %in% "female" & ourprophmat$comparison_group %in% "1w","logFC"]
  liverprophl2fc[i,"F W2"] <- ourprophmat[ourprophmat$sex %in% "female" & ourprophmat$comparison_group %in% "2w","logFC"]
  liverprophl2fc[i,"F W4"] <- ourprophmat[ourprophmat$sex %in% "female" & ourprophmat$comparison_group %in% "4w","logFC"]
  liverprophl2fc[i,"F W8"] <- ourprophmat[ourprophmat$sex %in% "female" & ourprophmat$comparison_group %in% "8w","logFC"]
  
  liverprophl2fc[i,"M W1"] <- ourprophmat[ourprophmat$sex %in% "male" & ourprophmat$comparison_group %in% "1w","logFC"]
  liverprophl2fc[i,"M W2"] <- ourprophmat[ourprophmat$sex %in% "male" & ourprophmat$comparison_group %in% "2w","logFC"]
  liverprophl2fc[i,"M W4"] <- ourprophmat[ourprophmat$sex %in% "male" & ourprophmat$comparison_group %in% "4w","logFC"]
  liverprophl2fc[i,"M W8"] <- ourprophmat[ourprophmat$sex %in% "male" & ourprophmat$comparison_group %in% "8w","logFC"]
}

lungprophl2fc <- matrix(0L,nrow = length(unique(lungproph$feature_ID)),ncol = 8)
rownames(lungprophl2fc) <- unique(lungproph$feature_ID)
colnames(lungprophl2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
for(i in 1:dim(lungprophl2fc)[1]){
  
  if(i%%100 == 0){
    print(paste("Lung:",toString(i)))
  }
  
  ourproph <- rownames(lungprophl2fc)[i]
  ourprophmat <- lungproph[lungproph$feature_ID %in% ourproph,]
  lungprophl2fc[i,"F W1"] <- ourprophmat[ourprophmat$sex %in% "female" & ourprophmat$comparison_group %in% "1w","logFC"]
  lungprophl2fc[i,"F W2"] <- ourprophmat[ourprophmat$sex %in% "female" & ourprophmat$comparison_group %in% "2w","logFC"]
  lungprophl2fc[i,"F W4"] <- ourprophmat[ourprophmat$sex %in% "female" & ourprophmat$comparison_group %in% "4w","logFC"]
  lungprophl2fc[i,"F W8"] <- ourprophmat[ourprophmat$sex %in% "female" & ourprophmat$comparison_group %in% "8w","logFC"]
  
  lungprophl2fc[i,"M W1"] <- ourprophmat[ourprophmat$sex %in% "male" & ourprophmat$comparison_group %in% "1w","logFC"]
  lungprophl2fc[i,"M W2"] <- ourprophmat[ourprophmat$sex %in% "male" & ourprophmat$comparison_group %in% "2w","logFC"]
  lungprophl2fc[i,"M W4"] <- ourprophmat[ourprophmat$sex %in% "male" & ourprophmat$comparison_group %in% "4w","logFC"]
  lungprophl2fc[i,"M W8"] <- ourprophmat[ourprophmat$sex %in% "male" & ourprophmat$comparison_group %in% "8w","logFC"]
}

whiteprophl2fc <- matrix(0L,nrow = length(unique(whiteproph$feature_ID)),ncol = 8)
rownames(whiteprophl2fc) <- unique(whiteproph$feature_ID)
colnames(whiteprophl2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
for(i in 1:dim(whiteprophl2fc)[1]){
  
  if(i%%100 == 0){
    print(paste("WhiteAd:",toString(i)))
  }
  
  ourproph <- rownames(whiteprophl2fc)[i]
  ourprophmat <- whiteproph[whiteproph$feature_ID %in% ourproph,]
  whiteprophl2fc[i,"F W1"] <- ourprophmat[ourprophmat$sex %in% "female" & ourprophmat$comparison_group %in% "1w","logFC"]
  whiteprophl2fc[i,"F W2"] <- ourprophmat[ourprophmat$sex %in% "female" & ourprophmat$comparison_group %in% "2w","logFC"]
  whiteprophl2fc[i,"F W4"] <- ourprophmat[ourprophmat$sex %in% "female" & ourprophmat$comparison_group %in% "4w","logFC"]
  whiteprophl2fc[i,"F W8"] <- ourprophmat[ourprophmat$sex %in% "female" & ourprophmat$comparison_group %in% "8w","logFC"]
  
  whiteprophl2fc[i,"M W1"] <- ourprophmat[ourprophmat$sex %in% "male" & ourprophmat$comparison_group %in% "1w","logFC"]
  whiteprophl2fc[i,"M W2"] <- ourprophmat[ourprophmat$sex %in% "male" & ourprophmat$comparison_group %in% "2w","logFC"]
  whiteprophl2fc[i,"M W4"] <- ourprophmat[ourprophmat$sex %in% "male" & ourprophmat$comparison_group %in% "4w","logFC"]
  whiteprophl2fc[i,"M W8"] <- ourprophmat[ourprophmat$sex %in% "male" & ourprophmat$comparison_group %in% "8w","logFC"]
}

save(list = c("gastroprophl2fc","heartprophl2fc","kidneyprophl2fc","liverprophl2fc","lungprophl2fc","whiteprophl2fc"),
     file = "prophl2fcmat.RData")


#####
# Generating the Input Files for Homer Analysis
####

# Use tables as input for homer findMotifsGenome.pl command with flags -size 50 -p 10 -preparse

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


#####
# Generating full L2FC matrices for ATACseq and RRBS
####

load("PASS1B Transcription Factor Paper Data/DEA Analysis/pass1b-06_epigen-atac-seq_dea.RData")

gastroatacl2fc <- matrix(0L,nrow = length(unique(epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "SKM-GN","feature_ID"])),ncol = 8)
rownames(gastroatacl2fc) <- unique(epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "SKM-GN","feature_ID"])
colnames(gastroatacl2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
gastroatacFW1timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "SKM-GN" & epigen_atac_seq$timewise_dea$sex %in% "female" & epigen_atac_seq$timewise_dea$comparison_group %in% "1w",]
rownames(gastroatacFW1timewise) <- gastroatacFW1timewise$feature_ID
gastroatacl2fc[,"F W1"] <- gastroatacFW1timewise[rownames(gastroatacl2fc),"logFC"]
rm(gastroatacFW1timewise)
gc()
gastroatacFW2timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "SKM-GN" & epigen_atac_seq$timewise_dea$sex %in% "female" & epigen_atac_seq$timewise_dea$comparison_group %in% "2w",]
rownames(gastroatacFW2timewise) <- gastroatacFW2timewise$feature_ID
gastroatacl2fc[,"F W2"] <- gastroatacFW2timewise[rownames(gastroatacl2fc),"logFC"]
rm(gastroatacFW2timewise)
gc()
gastroatacFW4timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "SKM-GN" & epigen_atac_seq$timewise_dea$sex %in% "female" & epigen_atac_seq$timewise_dea$comparison_group %in% "4w",]
rownames(gastroatacFW4timewise) <- gastroatacFW4timewise$feature_ID
gastroatacl2fc[,"F W4"] <- gastroatacFW4timewise[rownames(gastroatacl2fc),"logFC"]
rm(gastroatacFW4timewise)
gc()
gastroatacFW8timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "SKM-GN" & epigen_atac_seq$timewise_dea$sex %in% "female" & epigen_atac_seq$timewise_dea$comparison_group %in% "8w",]
rownames(gastroatacFW8timewise) <- gastroatacFW8timewise$feature_ID
gastroatacl2fc[,"F W8"] <- gastroatacFW8timewise[rownames(gastroatacl2fc),"logFC"]
rm(gastroatacFW8timewise)
gc()

gastroatacMW1timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "SKM-GN" & epigen_atac_seq$timewise_dea$sex %in% "male" & epigen_atac_seq$timewise_dea$comparison_group %in% "1w",]
rownames(gastroatacMW1timewise) <- gastroatacMW1timewise$feature_ID
gastroatacl2fc[,"M W1"] <- gastroatacMW1timewise[rownames(gastroatacl2fc),"logFC"]
rm(gastroatacMW1timewise)
gc()
gastroatacMW2timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "SKM-GN" & epigen_atac_seq$timewise_dea$sex %in% "male" & epigen_atac_seq$timewise_dea$comparison_group %in% "2w",]
rownames(gastroatacMW2timewise) <- gastroatacMW2timewise$feature_ID
gastroatacl2fc[,"M W2"] <- gastroatacMW2timewise[rownames(gastroatacl2fc),"logFC"]
rm(gastroatacMW2timewise)
gc()
gastroatacMW4timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "SKM-GN" & epigen_atac_seq$timewise_dea$sex %in% "male" & epigen_atac_seq$timewise_dea$comparison_group %in% "4w",]
rownames(gastroatacMW4timewise) <- gastroatacMW4timewise$feature_ID
gastroatacl2fc[,"M W4"] <- gastroatacMW4timewise[rownames(gastroatacl2fc),"logFC"]
rm(gastroatacMW4timewise)
gc()
gastroatacMW8timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "SKM-GN" & epigen_atac_seq$timewise_dea$sex %in% "male" & epigen_atac_seq$timewise_dea$comparison_group %in% "8w",]
rownames(gastroatacMW8timewise) <- gastroatacMW8timewise$feature_ID
gastroatacl2fc[,"M W8"] <- gastroatacMW8timewise[rownames(gastroatacl2fc),"logFC"]
rm(gastroatacMW8timewise)
gc()


heartatacl2fc <- matrix(0L,nrow = length(unique(epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "HEART","feature_ID"])),ncol = 8)
rownames(heartatacl2fc) <- unique(epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "HEART","feature_ID"])
colnames(heartatacl2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
heartatacFW1timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "HEART" & epigen_atac_seq$timewise_dea$sex %in% "female" & epigen_atac_seq$timewise_dea$comparison_group %in% "1w",]
rownames(heartatacFW1timewise) <- heartatacFW1timewise$feature_ID
heartatacl2fc[,"F W1"] <- heartatacFW1timewise[rownames(heartatacl2fc),"logFC"]
rm(heartatacFW1timewise)
gc()
heartatacFW2timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "HEART" & epigen_atac_seq$timewise_dea$sex %in% "female" & epigen_atac_seq$timewise_dea$comparison_group %in% "2w",]
rownames(heartatacFW2timewise) <- heartatacFW2timewise$feature_ID
heartatacl2fc[,"F W2"] <- heartatacFW2timewise[rownames(heartatacl2fc),"logFC"]
rm(heartatacFW2timewise)
gc()
heartatacFW4timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "HEART" & epigen_atac_seq$timewise_dea$sex %in% "female" & epigen_atac_seq$timewise_dea$comparison_group %in% "4w",]
rownames(heartatacFW4timewise) <- heartatacFW4timewise$feature_ID
heartatacl2fc[,"F W4"] <- heartatacFW4timewise[rownames(heartatacl2fc),"logFC"]
rm(heartatacFW4timewise)
gc()
heartatacFW8timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "HEART" & epigen_atac_seq$timewise_dea$sex %in% "female" & epigen_atac_seq$timewise_dea$comparison_group %in% "8w",]
rownames(heartatacFW8timewise) <- heartatacFW8timewise$feature_ID
heartatacl2fc[,"F W8"] <- heartatacFW8timewise[rownames(heartatacl2fc),"logFC"]
rm(heartatacFW8timewise)
gc()

heartatacMW1timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "HEART" & epigen_atac_seq$timewise_dea$sex %in% "male" & epigen_atac_seq$timewise_dea$comparison_group %in% "1w",]
rownames(heartatacMW1timewise) <- heartatacMW1timewise$feature_ID
heartatacl2fc[,"M W1"] <- heartatacMW1timewise[rownames(heartatacl2fc),"logFC"]
rm(heartatacMW1timewise)
gc()
heartatacMW2timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "HEART" & epigen_atac_seq$timewise_dea$sex %in% "male" & epigen_atac_seq$timewise_dea$comparison_group %in% "2w",]
rownames(heartatacMW2timewise) <- heartatacMW2timewise$feature_ID
heartatacl2fc[,"M W2"] <- heartatacMW2timewise[rownames(heartatacl2fc),"logFC"]
rm(heartatacMW2timewise)
gc()
heartatacMW4timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "HEART" & epigen_atac_seq$timewise_dea$sex %in% "male" & epigen_atac_seq$timewise_dea$comparison_group %in% "4w",]
rownames(heartatacMW4timewise) <- heartatacMW4timewise$feature_ID
heartatacl2fc[,"M W4"] <- heartatacMW4timewise[rownames(heartatacl2fc),"logFC"]
rm(heartatacMW4timewise)
gc()
heartatacMW8timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "HEART" & epigen_atac_seq$timewise_dea$sex %in% "male" & epigen_atac_seq$timewise_dea$comparison_group %in% "8w",]
rownames(heartatacMW8timewise) <- heartatacMW8timewise$feature_ID
heartatacl2fc[,"M W8"] <- heartatacMW8timewise[rownames(heartatacl2fc),"logFC"]
rm(heartatacMW8timewise)
gc()


hippoatacl2fc <- matrix(0L,nrow = length(unique(epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "HIPPOC","feature_ID"])),ncol = 8)
rownames(hippoatacl2fc) <- unique(epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "HIPPOC","feature_ID"])
colnames(hippoatacl2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
hippoatacFW1timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "HIPPOC" & epigen_atac_seq$timewise_dea$sex %in% "female" & epigen_atac_seq$timewise_dea$comparison_group %in% "1w",]
rownames(hippoatacFW1timewise) <- hippoatacFW1timewise$feature_ID
hippoatacl2fc[,"F W1"] <- hippoatacFW1timewise[rownames(hippoatacl2fc),"logFC"]
rm(hippoatacFW1timewise)
gc()
hippoatacFW2timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "HIPPOC" & epigen_atac_seq$timewise_dea$sex %in% "female" & epigen_atac_seq$timewise_dea$comparison_group %in% "2w",]
rownames(hippoatacFW2timewise) <- hippoatacFW2timewise$feature_ID
hippoatacl2fc[,"F W2"] <- hippoatacFW2timewise[rownames(hippoatacl2fc),"logFC"]
rm(hippoatacFW2timewise)
gc()
hippoatacFW4timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "HIPPOC" & epigen_atac_seq$timewise_dea$sex %in% "female" & epigen_atac_seq$timewise_dea$comparison_group %in% "4w",]
rownames(hippoatacFW4timewise) <- hippoatacFW4timewise$feature_ID
hippoatacl2fc[,"F W4"] <- hippoatacFW4timewise[rownames(hippoatacl2fc),"logFC"]
rm(hippoatacFW4timewise)
gc()
hippoatacFW8timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "HIPPOC" & epigen_atac_seq$timewise_dea$sex %in% "female" & epigen_atac_seq$timewise_dea$comparison_group %in% "8w",]
rownames(hippoatacFW8timewise) <- hippoatacFW8timewise$feature_ID
hippoatacl2fc[,"F W8"] <- hippoatacFW8timewise[rownames(hippoatacl2fc),"logFC"]
rm(hippoatacFW8timewise)
gc()

hippoatacMW1timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "HIPPOC" & epigen_atac_seq$timewise_dea$sex %in% "male" & epigen_atac_seq$timewise_dea$comparison_group %in% "1w",]
rownames(hippoatacMW1timewise) <- hippoatacMW1timewise$feature_ID
hippoatacl2fc[,"M W1"] <- hippoatacMW1timewise[rownames(hippoatacl2fc),"logFC"]
rm(hippoatacMW1timewise)
gc()
hippoatacMW2timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "HIPPOC" & epigen_atac_seq$timewise_dea$sex %in% "male" & epigen_atac_seq$timewise_dea$comparison_group %in% "2w",]
rownames(hippoatacMW2timewise) <- hippoatacMW2timewise$feature_ID
hippoatacl2fc[,"M W2"] <- hippoatacMW2timewise[rownames(hippoatacl2fc),"logFC"]
rm(hippoatacMW2timewise)
gc()
hippoatacMW4timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "HIPPOC" & epigen_atac_seq$timewise_dea$sex %in% "male" & epigen_atac_seq$timewise_dea$comparison_group %in% "4w",]
rownames(hippoatacMW4timewise) <- hippoatacMW4timewise$feature_ID
hippoatacl2fc[,"M W4"] <- hippoatacMW4timewise[rownames(hippoatacl2fc),"logFC"]
rm(hippoatacMW4timewise)
gc()
hippoatacMW8timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "HIPPOC" & epigen_atac_seq$timewise_dea$sex %in% "male" & epigen_atac_seq$timewise_dea$comparison_group %in% "8w",]
rownames(hippoatacMW8timewise) <- hippoatacMW8timewise$feature_ID
hippoatacl2fc[,"M W8"] <- hippoatacMW8timewise[rownames(hippoatacl2fc),"logFC"]
rm(hippoatacMW8timewise)
gc()


kidneyatacl2fc <- matrix(0L,nrow = length(unique(epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "KIDNEY","feature_ID"])),ncol = 8)
rownames(kidneyatacl2fc) <- unique(epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "KIDNEY","feature_ID"])
colnames(kidneyatacl2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
kidneyatacFW1timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "KIDNEY" & epigen_atac_seq$timewise_dea$sex %in% "female" & epigen_atac_seq$timewise_dea$comparison_group %in% "1w",]
rownames(kidneyatacFW1timewise) <- kidneyatacFW1timewise$feature_ID
kidneyatacl2fc[,"F W1"] <- kidneyatacFW1timewise[rownames(kidneyatacl2fc),"logFC"]
rm(kidneyatacFW1timewise)
gc()
kidneyatacFW2timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "KIDNEY" & epigen_atac_seq$timewise_dea$sex %in% "female" & epigen_atac_seq$timewise_dea$comparison_group %in% "2w",]
rownames(kidneyatacFW2timewise) <- kidneyatacFW2timewise$feature_ID
kidneyatacl2fc[,"F W2"] <- kidneyatacFW2timewise[rownames(kidneyatacl2fc),"logFC"]
rm(kidneyatacFW2timewise)
gc()
kidneyatacFW4timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "KIDNEY" & epigen_atac_seq$timewise_dea$sex %in% "female" & epigen_atac_seq$timewise_dea$comparison_group %in% "4w",]
rownames(kidneyatacFW4timewise) <- kidneyatacFW4timewise$feature_ID
kidneyatacl2fc[,"F W4"] <- kidneyatacFW4timewise[rownames(kidneyatacl2fc),"logFC"]
rm(kidneyatacFW4timewise)
gc()
kidneyatacFW8timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "KIDNEY" & epigen_atac_seq$timewise_dea$sex %in% "female" & epigen_atac_seq$timewise_dea$comparison_group %in% "8w",]
rownames(kidneyatacFW8timewise) <- kidneyatacFW8timewise$feature_ID
kidneyatacl2fc[,"F W8"] <- kidneyatacFW8timewise[rownames(kidneyatacl2fc),"logFC"]
rm(kidneyatacFW8timewise)
gc()

kidneyatacMW1timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "KIDNEY" & epigen_atac_seq$timewise_dea$sex %in% "male" & epigen_atac_seq$timewise_dea$comparison_group %in% "1w",]
rownames(kidneyatacMW1timewise) <- kidneyatacMW1timewise$feature_ID
kidneyatacl2fc[,"M W1"] <- kidneyatacMW1timewise[rownames(kidneyatacl2fc),"logFC"]
rm(kidneyatacMW1timewise)
gc()
kidneyatacMW2timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "KIDNEY" & epigen_atac_seq$timewise_dea$sex %in% "male" & epigen_atac_seq$timewise_dea$comparison_group %in% "2w",]
rownames(kidneyatacMW2timewise) <- kidneyatacMW2timewise$feature_ID
kidneyatacl2fc[,"M W2"] <- kidneyatacMW2timewise[rownames(kidneyatacl2fc),"logFC"]
rm(kidneyatacMW2timewise)
gc()
kidneyatacMW4timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "KIDNEY" & epigen_atac_seq$timewise_dea$sex %in% "male" & epigen_atac_seq$timewise_dea$comparison_group %in% "4w",]
rownames(kidneyatacMW4timewise) <- kidneyatacMW4timewise$feature_ID
kidneyatacl2fc[,"M W4"] <- kidneyatacMW4timewise[rownames(kidneyatacl2fc),"logFC"]
rm(kidneyatacMW4timewise)
gc()
kidneyatacMW8timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "KIDNEY" & epigen_atac_seq$timewise_dea$sex %in% "male" & epigen_atac_seq$timewise_dea$comparison_group %in% "8w",]
rownames(kidneyatacMW8timewise) <- kidneyatacMW8timewise$feature_ID
kidneyatacl2fc[,"M W8"] <- kidneyatacMW8timewise[rownames(kidneyatacl2fc),"logFC"]
rm(kidneyatacMW8timewise)
gc()


liveratacl2fc <- matrix(0L,nrow = length(unique(epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "LIVER","feature_ID"])),ncol = 8)
rownames(liveratacl2fc) <- unique(epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "LIVER","feature_ID"])
colnames(liveratacl2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
liveratacFW1timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "LIVER" & epigen_atac_seq$timewise_dea$sex %in% "female" & epigen_atac_seq$timewise_dea$comparison_group %in% "1w",]
rownames(liveratacFW1timewise) <- liveratacFW1timewise$feature_ID
liveratacl2fc[,"F W1"] <- liveratacFW1timewise[rownames(liveratacl2fc),"logFC"]
rm(liveratacFW1timewise)
gc()
liveratacFW2timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "LIVER" & epigen_atac_seq$timewise_dea$sex %in% "female" & epigen_atac_seq$timewise_dea$comparison_group %in% "2w",]
rownames(liveratacFW2timewise) <- liveratacFW2timewise$feature_ID
liveratacl2fc[,"F W2"] <- liveratacFW2timewise[rownames(liveratacl2fc),"logFC"]
rm(liveratacFW2timewise)
gc()
liveratacFW4timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "LIVER" & epigen_atac_seq$timewise_dea$sex %in% "female" & epigen_atac_seq$timewise_dea$comparison_group %in% "4w",]
rownames(liveratacFW4timewise) <- liveratacFW4timewise$feature_ID
liveratacl2fc[,"F W4"] <- liveratacFW4timewise[rownames(liveratacl2fc),"logFC"]
rm(liveratacFW4timewise)
gc()
liveratacFW8timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "LIVER" & epigen_atac_seq$timewise_dea$sex %in% "female" & epigen_atac_seq$timewise_dea$comparison_group %in% "8w",]
rownames(liveratacFW8timewise) <- liveratacFW8timewise$feature_ID
liveratacl2fc[,"F W8"] <- liveratacFW8timewise[rownames(liveratacl2fc),"logFC"]
rm(liveratacFW8timewise)
gc()

liveratacMW1timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "LIVER" & epigen_atac_seq$timewise_dea$sex %in% "male" & epigen_atac_seq$timewise_dea$comparison_group %in% "1w",]
rownames(liveratacMW1timewise) <- liveratacMW1timewise$feature_ID
liveratacl2fc[,"M W1"] <- liveratacMW1timewise[rownames(liveratacl2fc),"logFC"]
rm(liveratacMW1timewise)
gc()
liveratacMW2timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "LIVER" & epigen_atac_seq$timewise_dea$sex %in% "male" & epigen_atac_seq$timewise_dea$comparison_group %in% "2w",]
rownames(liveratacMW2timewise) <- liveratacMW2timewise$feature_ID
liveratacl2fc[,"M W2"] <- liveratacMW2timewise[rownames(liveratacl2fc),"logFC"]
rm(liveratacMW2timewise)
gc()
liveratacMW4timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "LIVER" & epigen_atac_seq$timewise_dea$sex %in% "male" & epigen_atac_seq$timewise_dea$comparison_group %in% "4w",]
rownames(liveratacMW4timewise) <- liveratacMW4timewise$feature_ID
liveratacl2fc[,"M W4"] <- liveratacMW4timewise[rownames(liveratacl2fc),"logFC"]
rm(liveratacMW4timewise)
gc()
liveratacMW8timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "LIVER" & epigen_atac_seq$timewise_dea$sex %in% "male" & epigen_atac_seq$timewise_dea$comparison_group %in% "8w",]
rownames(liveratacMW8timewise) <- liveratacMW8timewise$feature_ID
liveratacl2fc[,"M W8"] <- liveratacMW8timewise[rownames(liveratacl2fc),"logFC"]
rm(liveratacMW8timewise)
gc()


lungatacl2fc <- matrix(0L,nrow = length(unique(epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "LUNG","feature_ID"])),ncol = 8)
rownames(lungatacl2fc) <- unique(epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "LUNG","feature_ID"])
colnames(lungatacl2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
lungatacFW1timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "LUNG" & epigen_atac_seq$timewise_dea$sex %in% "female" & epigen_atac_seq$timewise_dea$comparison_group %in% "1w",]
rownames(lungatacFW1timewise) <- lungatacFW1timewise$feature_ID
lungatacl2fc[,"F W1"] <- lungatacFW1timewise[rownames(lungatacl2fc),"logFC"]
rm(lungatacFW1timewise)
gc()
lungatacFW2timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "LUNG" & epigen_atac_seq$timewise_dea$sex %in% "female" & epigen_atac_seq$timewise_dea$comparison_group %in% "2w",]
rownames(lungatacFW2timewise) <- lungatacFW2timewise$feature_ID
lungatacl2fc[,"F W2"] <- lungatacFW2timewise[rownames(lungatacl2fc),"logFC"]
rm(lungatacFW2timewise)
gc()
lungatacFW4timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "LUNG" & epigen_atac_seq$timewise_dea$sex %in% "female" & epigen_atac_seq$timewise_dea$comparison_group %in% "4w",]
rownames(lungatacFW4timewise) <- lungatacFW4timewise$feature_ID
lungatacl2fc[,"F W4"] <- lungatacFW4timewise[rownames(lungatacl2fc),"logFC"]
rm(lungatacFW4timewise)
gc()
lungatacFW8timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "LUNG" & epigen_atac_seq$timewise_dea$sex %in% "female" & epigen_atac_seq$timewise_dea$comparison_group %in% "8w",]
rownames(lungatacFW8timewise) <- lungatacFW8timewise$feature_ID
lungatacl2fc[,"F W8"] <- lungatacFW8timewise[rownames(lungatacl2fc),"logFC"]
rm(lungatacFW8timewise)
gc()

lungatacMW1timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "LUNG" & epigen_atac_seq$timewise_dea$sex %in% "male" & epigen_atac_seq$timewise_dea$comparison_group %in% "1w",]
rownames(lungatacMW1timewise) <- lungatacMW1timewise$feature_ID
lungatacl2fc[,"M W1"] <- lungatacMW1timewise[rownames(lungatacl2fc),"logFC"]
rm(lungatacMW1timewise)
gc()
lungatacMW2timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "LUNG" & epigen_atac_seq$timewise_dea$sex %in% "male" & epigen_atac_seq$timewise_dea$comparison_group %in% "2w",]
rownames(lungatacMW2timewise) <- lungatacMW2timewise$feature_ID
lungatacl2fc[,"M W2"] <- lungatacMW2timewise[rownames(lungatacl2fc),"logFC"]
rm(lungatacMW2timewise)
gc()
lungatacMW4timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "LUNG" & epigen_atac_seq$timewise_dea$sex %in% "male" & epigen_atac_seq$timewise_dea$comparison_group %in% "4w",]
rownames(lungatacMW4timewise) <- lungatacMW4timewise$feature_ID
lungatacl2fc[,"M W4"] <- lungatacMW4timewise[rownames(lungatacl2fc),"logFC"]
rm(lungatacMW4timewise)
gc()
lungatacMW8timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "LUNG" & epigen_atac_seq$timewise_dea$sex %in% "male" & epigen_atac_seq$timewise_dea$comparison_group %in% "8w",]
rownames(lungatacMW8timewise) <- lungatacMW8timewise$feature_ID
lungatacl2fc[,"M W8"] <- lungatacMW8timewise[rownames(lungatacl2fc),"logFC"]
rm(lungatacMW8timewise)
gc()


brownatacl2fc <- matrix(0L,nrow = length(unique(epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "BAT","feature_ID"])),ncol = 8)
rownames(brownatacl2fc) <- unique(epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "BAT","feature_ID"])
colnames(brownatacl2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
brownatacFW1timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "BAT" & epigen_atac_seq$timewise_dea$sex %in% "female" & epigen_atac_seq$timewise_dea$comparison_group %in% "1w",]
rownames(brownatacFW1timewise) <- brownatacFW1timewise$feature_ID
brownatacl2fc[,"F W1"] <- brownatacFW1timewise[rownames(brownatacl2fc),"logFC"]
rm(brownatacFW1timewise)
gc()
brownatacFW2timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "BAT" & epigen_atac_seq$timewise_dea$sex %in% "female" & epigen_atac_seq$timewise_dea$comparison_group %in% "2w",]
rownames(brownatacFW2timewise) <- brownatacFW2timewise$feature_ID
brownatacl2fc[,"F W2"] <- brownatacFW2timewise[rownames(brownatacl2fc),"logFC"]
rm(brownatacFW2timewise)
gc()
brownatacFW4timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "BAT" & epigen_atac_seq$timewise_dea$sex %in% "female" & epigen_atac_seq$timewise_dea$comparison_group %in% "4w",]
rownames(brownatacFW4timewise) <- brownatacFW4timewise$feature_ID
brownatacl2fc[,"F W4"] <- brownatacFW4timewise[rownames(brownatacl2fc),"logFC"]
rm(brownatacFW4timewise)
gc()
brownatacFW8timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "BAT" & epigen_atac_seq$timewise_dea$sex %in% "female" & epigen_atac_seq$timewise_dea$comparison_group %in% "8w",]
rownames(brownatacFW8timewise) <- brownatacFW8timewise$feature_ID
brownatacl2fc[,"F W8"] <- brownatacFW8timewise[rownames(brownatacl2fc),"logFC"]
rm(brownatacFW8timewise)
gc()

brownatacMW1timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "BAT" & epigen_atac_seq$timewise_dea$sex %in% "male" & epigen_atac_seq$timewise_dea$comparison_group %in% "1w",]
rownames(brownatacMW1timewise) <- brownatacMW1timewise$feature_ID
brownatacl2fc[,"M W1"] <- brownatacMW1timewise[rownames(brownatacl2fc),"logFC"]
rm(brownatacMW1timewise)
gc()
brownatacMW2timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "BAT" & epigen_atac_seq$timewise_dea$sex %in% "male" & epigen_atac_seq$timewise_dea$comparison_group %in% "2w",]
rownames(brownatacMW2timewise) <- brownatacMW2timewise$feature_ID
brownatacl2fc[,"M W2"] <- brownatacMW2timewise[rownames(brownatacl2fc),"logFC"]
rm(brownatacMW2timewise)
gc()
brownatacMW4timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "BAT" & epigen_atac_seq$timewise_dea$sex %in% "male" & epigen_atac_seq$timewise_dea$comparison_group %in% "4w",]
rownames(brownatacMW4timewise) <- brownatacMW4timewise$feature_ID
brownatacl2fc[,"M W4"] <- brownatacMW4timewise[rownames(brownatacl2fc),"logFC"]
rm(brownatacMW4timewise)
gc()
brownatacMW8timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "BAT" & epigen_atac_seq$timewise_dea$sex %in% "male" & epigen_atac_seq$timewise_dea$comparison_group %in% "8w",]
rownames(brownatacMW8timewise) <- brownatacMW8timewise$feature_ID
brownatacl2fc[,"M W8"] <- brownatacMW8timewise[rownames(brownatacl2fc),"logFC"]
rm(brownatacMW8timewise)
gc()

whiteatacl2fc <- matrix(0L,nrow = length(unique(epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "WAT-SC","feature_ID"])),ncol = 8)
rownames(whiteatacl2fc) <- unique(epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "WAT-SC","feature_ID"])
colnames(whiteatacl2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
whiteatacFW1timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "WAT-SC" & epigen_atac_seq$timewise_dea$sex %in% "female" & epigen_atac_seq$timewise_dea$comparison_group %in% "1w",]
rownames(whiteatacFW1timewise) <- whiteatacFW1timewise$feature_ID
whiteatacl2fc[,"F W1"] <- whiteatacFW1timewise[rownames(whiteatacl2fc),"logFC"]
rm(whiteatacFW1timewise)
gc()
whiteatacFW2timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "WAT-SC" & epigen_atac_seq$timewise_dea$sex %in% "female" & epigen_atac_seq$timewise_dea$comparison_group %in% "2w",]
rownames(whiteatacFW2timewise) <- whiteatacFW2timewise$feature_ID
whiteatacl2fc[,"F W2"] <- whiteatacFW2timewise[rownames(whiteatacl2fc),"logFC"]
rm(whiteatacFW2timewise)
gc()
whiteatacFW4timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "WAT-SC" & epigen_atac_seq$timewise_dea$sex %in% "female" & epigen_atac_seq$timewise_dea$comparison_group %in% "4w",]
rownames(whiteatacFW4timewise) <- whiteatacFW4timewise$feature_ID
whiteatacl2fc[,"F W4"] <- whiteatacFW4timewise[rownames(whiteatacl2fc),"logFC"]
rm(whiteatacFW4timewise)
gc()
whiteatacFW8timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "WAT-SC" & epigen_atac_seq$timewise_dea$sex %in% "female" & epigen_atac_seq$timewise_dea$comparison_group %in% "8w",]
rownames(whiteatacFW8timewise) <- whiteatacFW8timewise$feature_ID
whiteatacl2fc[,"F W8"] <- whiteatacFW8timewise[rownames(whiteatacl2fc),"logFC"]
rm(whiteatacFW8timewise)
gc()

whiteatacMW1timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "WAT-SC" & epigen_atac_seq$timewise_dea$sex %in% "male" & epigen_atac_seq$timewise_dea$comparison_group %in% "1w",]
rownames(whiteatacMW1timewise) <- whiteatacMW1timewise$feature_ID
whiteatacl2fc[,"M W1"] <- whiteatacMW1timewise[rownames(whiteatacl2fc),"logFC"]
rm(whiteatacMW1timewise)
gc()
whiteatacMW2timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "WAT-SC" & epigen_atac_seq$timewise_dea$sex %in% "male" & epigen_atac_seq$timewise_dea$comparison_group %in% "2w",]
rownames(whiteatacMW2timewise) <- whiteatacMW2timewise$feature_ID
whiteatacl2fc[,"M W2"] <- whiteatacMW2timewise[rownames(whiteatacl2fc),"logFC"]
rm(whiteatacMW2timewise)
gc()
whiteatacMW4timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "WAT-SC" & epigen_atac_seq$timewise_dea$sex %in% "male" & epigen_atac_seq$timewise_dea$comparison_group %in% "4w",]
rownames(whiteatacMW4timewise) <- whiteatacMW4timewise$feature_ID
whiteatacl2fc[,"M W4"] <- whiteatacMW4timewise[rownames(whiteatacl2fc),"logFC"]
rm(whiteatacMW4timewise)
gc()
whiteatacMW8timewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "WAT-SC" & epigen_atac_seq$timewise_dea$sex %in% "male" & epigen_atac_seq$timewise_dea$comparison_group %in% "8w",]
rownames(whiteatacMW8timewise) <- whiteatacMW8timewise$feature_ID
whiteatacl2fc[,"M W8"] <- whiteatacMW8timewise[rownames(whiteatacl2fc),"logFC"]
rm(whiteatacMW8timewise)
gc()


save(list = c("gastroatacl2fc","heartatacl2fc","hippoatacl2fc","kidneyatacl2fc",
              "liveratacl2fc","lungatacl2fc","brownatacl2fc","whiteatacl2fc"),file = "atacl2fcmat.RData")

# RRBS L2FC Generation
load("PASS1B Transcription Factor Paper Data/DEA Analysis/pass1b-06_epigen-rrbs_dea.RData")

gastromethl2fc <- matrix(0L,nrow = length(unique(epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "SKM-GN","feature_ID"])),ncol = 8)
rownames(gastromethl2fc) <- unique(epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "SKM-GN","feature_ID"])
colnames(gastromethl2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
gastromethFW1timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "SKM-GN" & epigen_rrbs$timewise_dea$sex %in% "female" & epigen_rrbs$timewise_dea$comparison_group %in% "1w",]
rownames(gastromethFW1timewise) <- gastromethFW1timewise$feature_ID
gastromethl2fc[,"F W1"] <- gastromethFW1timewise[rownames(gastromethl2fc),"logFC"]
rm(gastromethFW1timewise)
gc()
gastromethFW2timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "SKM-GN" & epigen_rrbs$timewise_dea$sex %in% "female" & epigen_rrbs$timewise_dea$comparison_group %in% "2w",]
rownames(gastromethFW2timewise) <- gastromethFW2timewise$feature_ID
gastromethl2fc[,"F W2"] <- gastromethFW2timewise[rownames(gastromethl2fc),"logFC"]
rm(gastromethFW2timewise)
gc()
gastromethFW4timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "SKM-GN" & epigen_rrbs$timewise_dea$sex %in% "female" & epigen_rrbs$timewise_dea$comparison_group %in% "4w",]
rownames(gastromethFW4timewise) <- gastromethFW4timewise$feature_ID
gastromethl2fc[,"F W4"] <- gastromethFW4timewise[rownames(gastromethl2fc),"logFC"]
rm(gastromethFW4timewise)
gc()
gastromethFW8timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "SKM-GN" & epigen_rrbs$timewise_dea$sex %in% "female" & epigen_rrbs$timewise_dea$comparison_group %in% "8w",]
rownames(gastromethFW8timewise) <- gastromethFW8timewise$feature_ID
gastromethl2fc[,"F W8"] <- gastromethFW8timewise[rownames(gastromethl2fc),"logFC"]
rm(gastromethFW8timewise)
gc()

gastromethMW1timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "SKM-GN" & epigen_rrbs$timewise_dea$sex %in% "male" & epigen_rrbs$timewise_dea$comparison_group %in% "1w",]
rownames(gastromethMW1timewise) <- gastromethMW1timewise$feature_ID
gastromethl2fc[,"M W1"] <- gastromethMW1timewise[rownames(gastromethl2fc),"logFC"]
rm(gastromethMW1timewise)
gc()
gastromethMW2timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "SKM-GN" & epigen_rrbs$timewise_dea$sex %in% "male" & epigen_rrbs$timewise_dea$comparison_group %in% "2w",]
rownames(gastromethMW2timewise) <- gastromethMW2timewise$feature_ID
gastromethl2fc[,"M W2"] <- gastromethMW2timewise[rownames(gastromethl2fc),"logFC"]
rm(gastromethMW2timewise)
gc()
gastromethMW4timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "SKM-GN" & epigen_rrbs$timewise_dea$sex %in% "male" & epigen_rrbs$timewise_dea$comparison_group %in% "4w",]
rownames(gastromethMW4timewise) <- gastromethMW4timewise$feature_ID
gastromethl2fc[,"M W4"] <- gastromethMW4timewise[rownames(gastromethl2fc),"logFC"]
rm(gastromethMW4timewise)
gc()
gastromethMW8timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "SKM-GN" & epigen_rrbs$timewise_dea$sex %in% "male" & epigen_rrbs$timewise_dea$comparison_group %in% "8w",]
rownames(gastromethMW8timewise) <- gastromethMW8timewise$feature_ID
gastromethl2fc[,"M W8"] <- gastromethMW8timewise[rownames(gastromethl2fc),"logFC"]
rm(gastromethMW8timewise)
gc()


heartmethl2fc <- matrix(0L,nrow = length(unique(epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "HEART","feature_ID"])),ncol = 8)
rownames(heartmethl2fc) <- unique(epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "HEART","feature_ID"])
colnames(heartmethl2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
heartmethFW1timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "HEART" & epigen_rrbs$timewise_dea$sex %in% "female" & epigen_rrbs$timewise_dea$comparison_group %in% "1w",]
rownames(heartmethFW1timewise) <- heartmethFW1timewise$feature_ID
heartmethl2fc[,"F W1"] <- heartmethFW1timewise[rownames(heartmethl2fc),"logFC"]
rm(heartmethFW1timewise)
gc()
heartmethFW2timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "HEART" & epigen_rrbs$timewise_dea$sex %in% "female" & epigen_rrbs$timewise_dea$comparison_group %in% "2w",]
rownames(heartmethFW2timewise) <- heartmethFW2timewise$feature_ID
heartmethl2fc[,"F W2"] <- heartmethFW2timewise[rownames(heartmethl2fc),"logFC"]
rm(heartmethFW2timewise)
gc()
heartmethFW4timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "HEART" & epigen_rrbs$timewise_dea$sex %in% "female" & epigen_rrbs$timewise_dea$comparison_group %in% "4w",]
rownames(heartmethFW4timewise) <- heartmethFW4timewise$feature_ID
heartmethl2fc[,"F W4"] <- heartmethFW4timewise[rownames(heartmethl2fc),"logFC"]
rm(heartmethFW4timewise)
gc()
heartmethFW8timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "HEART" & epigen_rrbs$timewise_dea$sex %in% "female" & epigen_rrbs$timewise_dea$comparison_group %in% "8w",]
rownames(heartmethFW8timewise) <- heartmethFW8timewise$feature_ID
heartmethl2fc[,"F W8"] <- heartmethFW8timewise[rownames(heartmethl2fc),"logFC"]
rm(heartmethFW8timewise)
gc()

heartmethMW1timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "HEART" & epigen_rrbs$timewise_dea$sex %in% "male" & epigen_rrbs$timewise_dea$comparison_group %in% "1w",]
rownames(heartmethMW1timewise) <- heartmethMW1timewise$feature_ID
heartmethl2fc[,"M W1"] <- heartmethMW1timewise[rownames(heartmethl2fc),"logFC"]
rm(heartmethMW1timewise)
gc()
heartmethMW2timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "HEART" & epigen_rrbs$timewise_dea$sex %in% "male" & epigen_rrbs$timewise_dea$comparison_group %in% "2w",]
rownames(heartmethMW2timewise) <- heartmethMW2timewise$feature_ID
heartmethl2fc[,"M W2"] <- heartmethMW2timewise[rownames(heartmethl2fc),"logFC"]
rm(heartmethMW2timewise)
gc()
heartmethMW4timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "HEART" & epigen_rrbs$timewise_dea$sex %in% "male" & epigen_rrbs$timewise_dea$comparison_group %in% "4w",]
rownames(heartmethMW4timewise) <- heartmethMW4timewise$feature_ID
heartmethl2fc[,"M W4"] <- heartmethMW4timewise[rownames(heartmethl2fc),"logFC"]
rm(heartmethMW4timewise)
gc()
heartmethMW8timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "HEART" & epigen_rrbs$timewise_dea$sex %in% "male" & epigen_rrbs$timewise_dea$comparison_group %in% "8w",]
rownames(heartmethMW8timewise) <- heartmethMW8timewise$feature_ID
heartmethl2fc[,"M W8"] <- heartmethMW8timewise[rownames(heartmethl2fc),"logFC"]
rm(heartmethMW8timewise)
gc()


hippomethl2fc <- matrix(0L,nrow = length(unique(epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "HIPPOC","feature_ID"])),ncol = 8)
rownames(hippomethl2fc) <- unique(epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "HIPPOC","feature_ID"])
colnames(hippomethl2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
hippomethFW1timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "HIPPOC" & epigen_rrbs$timewise_dea$sex %in% "female" & epigen_rrbs$timewise_dea$comparison_group %in% "1w",]
rownames(hippomethFW1timewise) <- hippomethFW1timewise$feature_ID
hippomethl2fc[,"F W1"] <- hippomethFW1timewise[rownames(hippomethl2fc),"logFC"]
rm(hippomethFW1timewise)
gc()
hippomethFW2timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "HIPPOC" & epigen_rrbs$timewise_dea$sex %in% "female" & epigen_rrbs$timewise_dea$comparison_group %in% "2w",]
rownames(hippomethFW2timewise) <- hippomethFW2timewise$feature_ID
hippomethl2fc[,"F W2"] <- hippomethFW2timewise[rownames(hippomethl2fc),"logFC"]
rm(hippomethFW2timewise)
gc()
hippomethFW4timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "HIPPOC" & epigen_rrbs$timewise_dea$sex %in% "female" & epigen_rrbs$timewise_dea$comparison_group %in% "4w",]
rownames(hippomethFW4timewise) <- hippomethFW4timewise$feature_ID
hippomethl2fc[,"F W4"] <- hippomethFW4timewise[rownames(hippomethl2fc),"logFC"]
rm(hippomethFW4timewise)
gc()
hippomethFW8timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "HIPPOC" & epigen_rrbs$timewise_dea$sex %in% "female" & epigen_rrbs$timewise_dea$comparison_group %in% "8w",]
rownames(hippomethFW8timewise) <- hippomethFW8timewise$feature_ID
hippomethl2fc[,"F W8"] <- hippomethFW8timewise[rownames(hippomethl2fc),"logFC"]
rm(hippomethFW8timewise)
gc()

hippomethMW1timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "HIPPOC" & epigen_rrbs$timewise_dea$sex %in% "male" & epigen_rrbs$timewise_dea$comparison_group %in% "1w",]
rownames(hippomethMW1timewise) <- hippomethMW1timewise$feature_ID
hippomethl2fc[,"M W1"] <- hippomethMW1timewise[rownames(hippomethl2fc),"logFC"]
rm(hippomethMW1timewise)
gc()
hippomethMW2timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "HIPPOC" & epigen_rrbs$timewise_dea$sex %in% "male" & epigen_rrbs$timewise_dea$comparison_group %in% "2w",]
rownames(hippomethMW2timewise) <- hippomethMW2timewise$feature_ID
hippomethl2fc[,"M W2"] <- hippomethMW2timewise[rownames(hippomethl2fc),"logFC"]
rm(hippomethMW2timewise)
gc()
hippomethMW4timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "HIPPOC" & epigen_rrbs$timewise_dea$sex %in% "male" & epigen_rrbs$timewise_dea$comparison_group %in% "4w",]
rownames(hippomethMW4timewise) <- hippomethMW4timewise$feature_ID
hippomethl2fc[,"M W4"] <- hippomethMW4timewise[rownames(hippomethl2fc),"logFC"]
rm(hippomethMW4timewise)
gc()
hippomethMW8timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "HIPPOC" & epigen_rrbs$timewise_dea$sex %in% "male" & epigen_rrbs$timewise_dea$comparison_group %in% "8w",]
rownames(hippomethMW8timewise) <- hippomethMW8timewise$feature_ID
hippomethl2fc[,"M W8"] <- hippomethMW8timewise[rownames(hippomethl2fc),"logFC"]
rm(hippomethMW8timewise)
gc()


kidneymethl2fc <- matrix(0L,nrow = length(unique(epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "KIDNEY","feature_ID"])),ncol = 8)
rownames(kidneymethl2fc) <- unique(epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "KIDNEY","feature_ID"])
colnames(kidneymethl2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
kidneymethFW1timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "KIDNEY" & epigen_rrbs$timewise_dea$sex %in% "female" & epigen_rrbs$timewise_dea$comparison_group %in% "1w",]
rownames(kidneymethFW1timewise) <- kidneymethFW1timewise$feature_ID
kidneymethl2fc[,"F W1"] <- kidneymethFW1timewise[rownames(kidneymethl2fc),"logFC"]
rm(kidneymethFW1timewise)
gc()
kidneymethFW2timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "KIDNEY" & epigen_rrbs$timewise_dea$sex %in% "female" & epigen_rrbs$timewise_dea$comparison_group %in% "2w",]
rownames(kidneymethFW2timewise) <- kidneymethFW2timewise$feature_ID
kidneymethl2fc[,"F W2"] <- kidneymethFW2timewise[rownames(kidneymethl2fc),"logFC"]
rm(kidneymethFW2timewise)
gc()
kidneymethFW4timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "KIDNEY" & epigen_rrbs$timewise_dea$sex %in% "female" & epigen_rrbs$timewise_dea$comparison_group %in% "4w",]
rownames(kidneymethFW4timewise) <- kidneymethFW4timewise$feature_ID
kidneymethl2fc[,"F W4"] <- kidneymethFW4timewise[rownames(kidneymethl2fc),"logFC"]
rm(kidneymethFW4timewise)
gc()
kidneymethFW8timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "KIDNEY" & epigen_rrbs$timewise_dea$sex %in% "female" & epigen_rrbs$timewise_dea$comparison_group %in% "8w",]
rownames(kidneymethFW8timewise) <- kidneymethFW8timewise$feature_ID
kidneymethl2fc[,"F W8"] <- kidneymethFW8timewise[rownames(kidneymethl2fc),"logFC"]
rm(kidneymethFW8timewise)
gc()

kidneymethMW1timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "KIDNEY" & epigen_rrbs$timewise_dea$sex %in% "male" & epigen_rrbs$timewise_dea$comparison_group %in% "1w",]
rownames(kidneymethMW1timewise) <- kidneymethMW1timewise$feature_ID
kidneymethl2fc[,"M W1"] <- kidneymethMW1timewise[rownames(kidneymethl2fc),"logFC"]
rm(kidneymethMW1timewise)
gc()
kidneymethMW2timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "KIDNEY" & epigen_rrbs$timewise_dea$sex %in% "male" & epigen_rrbs$timewise_dea$comparison_group %in% "2w",]
rownames(kidneymethMW2timewise) <- kidneymethMW2timewise$feature_ID
kidneymethl2fc[,"M W2"] <- kidneymethMW2timewise[rownames(kidneymethl2fc),"logFC"]
rm(kidneymethMW2timewise)
gc()
kidneymethMW4timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "KIDNEY" & epigen_rrbs$timewise_dea$sex %in% "male" & epigen_rrbs$timewise_dea$comparison_group %in% "4w",]
rownames(kidneymethMW4timewise) <- kidneymethMW4timewise$feature_ID
kidneymethl2fc[,"M W4"] <- kidneymethMW4timewise[rownames(kidneymethl2fc),"logFC"]
rm(kidneymethMW4timewise)
gc()
kidneymethMW8timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "KIDNEY" & epigen_rrbs$timewise_dea$sex %in% "male" & epigen_rrbs$timewise_dea$comparison_group %in% "8w",]
rownames(kidneymethMW8timewise) <- kidneymethMW8timewise$feature_ID
kidneymethl2fc[,"M W8"] <- kidneymethMW8timewise[rownames(kidneymethl2fc),"logFC"]
rm(kidneymethMW8timewise)
gc()


livermethl2fc <- matrix(0L,nrow = length(unique(epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "LIVER","feature_ID"])),ncol = 8)
rownames(livermethl2fc) <- unique(epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "LIVER","feature_ID"])
colnames(livermethl2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
livermethFW1timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "LIVER" & epigen_rrbs$timewise_dea$sex %in% "female" & epigen_rrbs$timewise_dea$comparison_group %in% "1w",]
rownames(livermethFW1timewise) <- livermethFW1timewise$feature_ID
livermethl2fc[,"F W1"] <- livermethFW1timewise[rownames(livermethl2fc),"logFC"]
rm(livermethFW1timewise)
gc()
livermethFW2timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "LIVER" & epigen_rrbs$timewise_dea$sex %in% "female" & epigen_rrbs$timewise_dea$comparison_group %in% "2w",]
rownames(livermethFW2timewise) <- livermethFW2timewise$feature_ID
livermethl2fc[,"F W2"] <- livermethFW2timewise[rownames(livermethl2fc),"logFC"]
rm(livermethFW2timewise)
gc()
livermethFW4timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "LIVER" & epigen_rrbs$timewise_dea$sex %in% "female" & epigen_rrbs$timewise_dea$comparison_group %in% "4w",]
rownames(livermethFW4timewise) <- livermethFW4timewise$feature_ID
livermethl2fc[,"F W4"] <- livermethFW4timewise[rownames(livermethl2fc),"logFC"]
rm(livermethFW4timewise)
gc()
livermethFW8timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "LIVER" & epigen_rrbs$timewise_dea$sex %in% "female" & epigen_rrbs$timewise_dea$comparison_group %in% "8w",]
rownames(livermethFW8timewise) <- livermethFW8timewise$feature_ID
livermethl2fc[,"F W8"] <- livermethFW8timewise[rownames(livermethl2fc),"logFC"]
rm(livermethFW8timewise)
gc()

livermethMW1timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "LIVER" & epigen_rrbs$timewise_dea$sex %in% "male" & epigen_rrbs$timewise_dea$comparison_group %in% "1w",]
rownames(livermethMW1timewise) <- livermethMW1timewise$feature_ID
livermethl2fc[,"M W1"] <- livermethMW1timewise[rownames(livermethl2fc),"logFC"]
rm(livermethMW1timewise)
gc()
livermethMW2timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "LIVER" & epigen_rrbs$timewise_dea$sex %in% "male" & epigen_rrbs$timewise_dea$comparison_group %in% "2w",]
rownames(livermethMW2timewise) <- livermethMW2timewise$feature_ID
livermethl2fc[,"M W2"] <- livermethMW2timewise[rownames(livermethl2fc),"logFC"]
rm(livermethMW2timewise)
gc()
livermethMW4timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "LIVER" & epigen_rrbs$timewise_dea$sex %in% "male" & epigen_rrbs$timewise_dea$comparison_group %in% "4w",]
rownames(livermethMW4timewise) <- livermethMW4timewise$feature_ID
livermethl2fc[,"M W4"] <- livermethMW4timewise[rownames(livermethl2fc),"logFC"]
rm(livermethMW4timewise)
gc()
livermethMW8timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "LIVER" & epigen_rrbs$timewise_dea$sex %in% "male" & epigen_rrbs$timewise_dea$comparison_group %in% "8w",]
rownames(livermethMW8timewise) <- livermethMW8timewise$feature_ID
livermethl2fc[,"M W8"] <- livermethMW8timewise[rownames(livermethl2fc),"logFC"]
rm(livermethMW8timewise)
gc()


lungmethl2fc <- matrix(0L,nrow = length(unique(epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "LUNG","feature_ID"])),ncol = 8)
rownames(lungmethl2fc) <- unique(epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "LUNG","feature_ID"])
colnames(lungmethl2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
lungmethFW1timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "LUNG" & epigen_rrbs$timewise_dea$sex %in% "female" & epigen_rrbs$timewise_dea$comparison_group %in% "1w",]
rownames(lungmethFW1timewise) <- lungmethFW1timewise$feature_ID
lungmethl2fc[,"F W1"] <- lungmethFW1timewise[rownames(lungmethl2fc),"logFC"]
rm(lungmethFW1timewise)
gc()
lungmethFW2timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "LUNG" & epigen_rrbs$timewise_dea$sex %in% "female" & epigen_rrbs$timewise_dea$comparison_group %in% "2w",]
rownames(lungmethFW2timewise) <- lungmethFW2timewise$feature_ID
lungmethl2fc[,"F W2"] <- lungmethFW2timewise[rownames(lungmethl2fc),"logFC"]
rm(lungmethFW2timewise)
gc()
lungmethFW4timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "LUNG" & epigen_rrbs$timewise_dea$sex %in% "female" & epigen_rrbs$timewise_dea$comparison_group %in% "4w",]
rownames(lungmethFW4timewise) <- lungmethFW4timewise$feature_ID
lungmethl2fc[,"F W4"] <- lungmethFW4timewise[rownames(lungmethl2fc),"logFC"]
rm(lungmethFW4timewise)
gc()
lungmethFW8timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "LUNG" & epigen_rrbs$timewise_dea$sex %in% "female" & epigen_rrbs$timewise_dea$comparison_group %in% "8w",]
rownames(lungmethFW8timewise) <- lungmethFW8timewise$feature_ID
lungmethl2fc[,"F W8"] <- lungmethFW8timewise[rownames(lungmethl2fc),"logFC"]
rm(lungmethFW8timewise)
gc()

lungmethMW1timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "LUNG" & epigen_rrbs$timewise_dea$sex %in% "male" & epigen_rrbs$timewise_dea$comparison_group %in% "1w",]
rownames(lungmethMW1timewise) <- lungmethMW1timewise$feature_ID
lungmethl2fc[,"M W1"] <- lungmethMW1timewise[rownames(lungmethl2fc),"logFC"]
rm(lungmethMW1timewise)
gc()
lungmethMW2timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "LUNG" & epigen_rrbs$timewise_dea$sex %in% "male" & epigen_rrbs$timewise_dea$comparison_group %in% "2w",]
rownames(lungmethMW2timewise) <- lungmethMW2timewise$feature_ID
lungmethl2fc[,"M W2"] <- lungmethMW2timewise[rownames(lungmethl2fc),"logFC"]
rm(lungmethMW2timewise)
gc()
lungmethMW4timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "LUNG" & epigen_rrbs$timewise_dea$sex %in% "male" & epigen_rrbs$timewise_dea$comparison_group %in% "4w",]
rownames(lungmethMW4timewise) <- lungmethMW4timewise$feature_ID
lungmethl2fc[,"M W4"] <- lungmethMW4timewise[rownames(lungmethl2fc),"logFC"]
rm(lungmethMW4timewise)
gc()
lungmethMW8timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "LUNG" & epigen_rrbs$timewise_dea$sex %in% "male" & epigen_rrbs$timewise_dea$comparison_group %in% "8w",]
rownames(lungmethMW8timewise) <- lungmethMW8timewise$feature_ID
lungmethl2fc[,"M W8"] <- lungmethMW8timewise[rownames(lungmethl2fc),"logFC"]
rm(lungmethMW8timewise)
gc()


brownmethl2fc <- matrix(0L,nrow = length(unique(epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "BAT","feature_ID"])),ncol = 8)
rownames(brownmethl2fc) <- unique(epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "BAT","feature_ID"])
colnames(brownmethl2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
brownmethFW1timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "BAT" & epigen_rrbs$timewise_dea$sex %in% "female" & epigen_rrbs$timewise_dea$comparison_group %in% "1w",]
rownames(brownmethFW1timewise) <- brownmethFW1timewise$feature_ID
brownmethl2fc[,"F W1"] <- brownmethFW1timewise[rownames(brownmethl2fc),"logFC"]
rm(brownmethFW1timewise)
gc()
brownmethFW2timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "BAT" & epigen_rrbs$timewise_dea$sex %in% "female" & epigen_rrbs$timewise_dea$comparison_group %in% "2w",]
rownames(brownmethFW2timewise) <- brownmethFW2timewise$feature_ID
brownmethl2fc[,"F W2"] <- brownmethFW2timewise[rownames(brownmethl2fc),"logFC"]
rm(brownmethFW2timewise)
gc()
brownmethFW4timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "BAT" & epigen_rrbs$timewise_dea$sex %in% "female" & epigen_rrbs$timewise_dea$comparison_group %in% "4w",]
rownames(brownmethFW4timewise) <- brownmethFW4timewise$feature_ID
brownmethl2fc[,"F W4"] <- brownmethFW4timewise[rownames(brownmethl2fc),"logFC"]
rm(brownmethFW4timewise)
gc()
brownmethFW8timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "BAT" & epigen_rrbs$timewise_dea$sex %in% "female" & epigen_rrbs$timewise_dea$comparison_group %in% "8w",]
rownames(brownmethFW8timewise) <- brownmethFW8timewise$feature_ID
brownmethl2fc[,"F W8"] <- brownmethFW8timewise[rownames(brownmethl2fc),"logFC"]
rm(brownmethFW8timewise)
gc()

brownmethMW1timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "BAT" & epigen_rrbs$timewise_dea$sex %in% "male" & epigen_rrbs$timewise_dea$comparison_group %in% "1w",]
rownames(brownmethMW1timewise) <- brownmethMW1timewise$feature_ID
brownmethl2fc[,"M W1"] <- brownmethMW1timewise[rownames(brownmethl2fc),"logFC"]
rm(brownmethMW1timewise)
gc()
brownmethMW2timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "BAT" & epigen_rrbs$timewise_dea$sex %in% "male" & epigen_rrbs$timewise_dea$comparison_group %in% "2w",]
rownames(brownmethMW2timewise) <- brownmethMW2timewise$feature_ID
brownmethl2fc[,"M W2"] <- brownmethMW2timewise[rownames(brownmethl2fc),"logFC"]
rm(brownmethMW2timewise)
gc()
brownmethMW4timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "BAT" & epigen_rrbs$timewise_dea$sex %in% "male" & epigen_rrbs$timewise_dea$comparison_group %in% "4w",]
rownames(brownmethMW4timewise) <- brownmethMW4timewise$feature_ID
brownmethl2fc[,"M W4"] <- brownmethMW4timewise[rownames(brownmethl2fc),"logFC"]
rm(brownmethMW4timewise)
gc()
brownmethMW8timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "BAT" & epigen_rrbs$timewise_dea$sex %in% "male" & epigen_rrbs$timewise_dea$comparison_group %in% "8w",]
rownames(brownmethMW8timewise) <- brownmethMW8timewise$feature_ID
brownmethl2fc[,"M W8"] <- brownmethMW8timewise[rownames(brownmethl2fc),"logFC"]
rm(brownmethMW8timewise)
gc()


whitemethl2fc <- matrix(0L,nrow = length(unique(epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "WAT-SC","feature_ID"])),ncol = 8)
rownames(whitemethl2fc) <- unique(epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "WAT-SC","feature_ID"])
colnames(whitemethl2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")
whitemethFW1timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "WAT-SC" & epigen_rrbs$timewise_dea$sex %in% "female" & epigen_rrbs$timewise_dea$comparison_group %in% "1w",]
rownames(whitemethFW1timewise) <- whitemethFW1timewise$feature_ID
whitemethl2fc[,"F W1"] <- whitemethFW1timewise[rownames(whitemethl2fc),"logFC"]
rm(whitemethFW1timewise)
gc()
whitemethFW2timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "WAT-SC" & epigen_rrbs$timewise_dea$sex %in% "female" & epigen_rrbs$timewise_dea$comparison_group %in% "2w",]
rownames(whitemethFW2timewise) <- whitemethFW2timewise$feature_ID
whitemethl2fc[,"F W2"] <- whitemethFW2timewise[rownames(whitemethl2fc),"logFC"]
rm(whitemethFW2timewise)
gc()
whitemethFW4timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "WAT-SC" & epigen_rrbs$timewise_dea$sex %in% "female" & epigen_rrbs$timewise_dea$comparison_group %in% "4w",]
rownames(whitemethFW4timewise) <- whitemethFW4timewise$feature_ID
whitemethl2fc[,"F W4"] <- whitemethFW4timewise[rownames(whitemethl2fc),"logFC"]
rm(whitemethFW4timewise)
gc()
whitemethFW8timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "WAT-SC" & epigen_rrbs$timewise_dea$sex %in% "female" & epigen_rrbs$timewise_dea$comparison_group %in% "8w",]
rownames(whitemethFW8timewise) <- whitemethFW8timewise$feature_ID
whitemethl2fc[,"F W8"] <- whitemethFW8timewise[rownames(whitemethl2fc),"logFC"]
rm(whitemethFW8timewise)
gc()

whitemethMW1timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "WAT-SC" & epigen_rrbs$timewise_dea$sex %in% "male" & epigen_rrbs$timewise_dea$comparison_group %in% "1w",]
rownames(whitemethMW1timewise) <- whitemethMW1timewise$feature_ID
whitemethl2fc[,"M W1"] <- whitemethMW1timewise[rownames(whitemethl2fc),"logFC"]
rm(whitemethMW1timewise)
gc()
whitemethMW2timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "WAT-SC" & epigen_rrbs$timewise_dea$sex %in% "male" & epigen_rrbs$timewise_dea$comparison_group %in% "2w",]
rownames(whitemethMW2timewise) <- whitemethMW2timewise$feature_ID
whitemethl2fc[,"M W2"] <- whitemethMW2timewise[rownames(whitemethl2fc),"logFC"]
rm(whitemethMW2timewise)
gc()
whitemethMW4timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "WAT-SC" & epigen_rrbs$timewise_dea$sex %in% "male" & epigen_rrbs$timewise_dea$comparison_group %in% "4w",]
rownames(whitemethMW4timewise) <- whitemethMW4timewise$feature_ID
whitemethl2fc[,"M W4"] <- whitemethMW4timewise[rownames(whitemethl2fc),"logFC"]
rm(whitemethMW4timewise)
gc()
whitemethMW8timewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$tissue_abbreviation %in% "WAT-SC" & epigen_rrbs$timewise_dea$sex %in% "male" & epigen_rrbs$timewise_dea$comparison_group %in% "8w",]
rownames(whitemethMW8timewise) <- whitemethMW8timewise$feature_ID
whitemethl2fc[,"M W8"] <- whitemethMW8timewise[rownames(whitemethl2fc),"logFC"]
rm(whitemethMW8timewise)
gc()

save(list = c("gastromethl2fc","heartmethl2fc","hippomethl2fc","kidneymethl2fc",
              "livermethl2fc","lungmethl2fc","brownmethl2fc","whitemethl2fc"),file = "methl2fcmat.RData")

#####
# Now to divide the genes/peaks/methyl sites into week 8 upreg vs downreg subsets
####

# First ATAC

load("PASS1B Transcription Factor Paper Data/DEA Analysis/pass1b-06_epigen-atac-seq_dea.RData")

load("omesigdata.RData")

gastroatacsigtimewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$feature_ID %in% gastroatacsig & epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "SKM-GN",]
gastroatacsigupordowndf <- data.frame(row.names = gastroatacsig,
                                      "W8 F" = rep(0,length(gastroatacsig)),
                                      "W8 M" = rep(0,length(gastroatacsig)),
                                      "UpReg" = rep(0,length(gastroatacsig)),
                                      "DownReg" = rep(0,length(gastroatacsig)))
for(i in 1:length(gastroatacsig)){
  if(i%%50 == 0){
    print(i)
  }
  ourpeak <- gastroatacsig[i]
  gastroatacsigupordowndf[i,"W8.F"] <- gastroatacsigtimewise[gastroatacsigtimewise$feature_ID %in% ourpeak &
                                                               gastroatacsigtimewise$sex %in% "female" &
                                                               gastroatacsigtimewise$comparison_group %in% "8w","logFC"]
  gastroatacsigupordowndf[i,"W8.M"] <- gastroatacsigtimewise[gastroatacsigtimewise$feature_ID %in% ourpeak &
                                                               gastroatacsigtimewise$sex %in% "male" &
                                                               gastroatacsigtimewise$comparison_group %in% "8w","logFC"]
  if(gastroatacsigupordowndf[i,"W8.F"] > 0.1 & gastroatacsigupordowndf[i,"W8.M"] > 0.1){
    gastroatacsigupordowndf[i,"UpReg"] <- 1  
  }
  if(gastroatacsigupordowndf[i,"W8.F"] < -0.1 & gastroatacsigupordowndf[i,"W8.M"] < -0.1){
    gastroatacsigupordowndf[i,"DownReg"] <- 1  
  }
}

# heart
heartatacsigtimewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$feature_ID %in% heartatacsig & epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "HEART",]
heartatacsigupordowndf <- data.frame(row.names = heartatacsig,
                                     "W8 F" = rep(0,length(heartatacsig)),
                                     "W8 M" = rep(0,length(heartatacsig)),
                                     "UpReg" = rep(0,length(heartatacsig)),
                                     "DownReg" = rep(0,length(heartatacsig)))
for(i in 1:length(heartatacsig)){
  if(i%%50 == 0){
    print(i)
  }
  ourpeak <- heartatacsig[i]
  heartatacsigupordowndf[i,"W8.F"] <- heartatacsigtimewise[heartatacsigtimewise$feature_ID %in% ourpeak &
                                                             heartatacsigtimewise$sex %in% "female" &
                                                             heartatacsigtimewise$comparison_group %in% "8w","logFC"]
  heartatacsigupordowndf[i,"W8.M"] <- heartatacsigtimewise[heartatacsigtimewise$feature_ID %in% ourpeak &
                                                             heartatacsigtimewise$sex %in% "male" &
                                                             heartatacsigtimewise$comparison_group %in% "8w","logFC"]
  if(heartatacsigupordowndf[i,"W8.F"] > 0.1 & heartatacsigupordowndf[i,"W8.M"] > 0.1){
    heartatacsigupordowndf[i,"UpReg"] <- 1  
  }
  if(heartatacsigupordowndf[i,"W8.F"] < -0.1 & heartatacsigupordowndf[i,"W8.M"] < -0.1){
    heartatacsigupordowndf[i,"DownReg"] <- 1  
  }
}

# hippo
hippoatacsigtimewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$feature_ID %in% hippoatacsig & epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "HIPPOC",]
hippoatacsigupordowndf <- data.frame(row.names = hippoatacsig,
                                     "W8 F" = rep(0,length(hippoatacsig)),
                                     "W8 M" = rep(0,length(hippoatacsig)),
                                     "UpReg" = rep(0,length(hippoatacsig)),
                                     "DownReg" = rep(0,length(hippoatacsig)))
for(i in 1:length(hippoatacsig)){
  if(i%%50 == 0){
    print(i)
  }
  ourpeak <- hippoatacsig[i]
  hippoatacsigupordowndf[i,"W8.F"] <- hippoatacsigtimewise[hippoatacsigtimewise$feature_ID %in% ourpeak &
                                                             hippoatacsigtimewise$sex %in% "female" &
                                                             hippoatacsigtimewise$comparison_group %in% "8w","logFC"]
  hippoatacsigupordowndf[i,"W8.M"] <- hippoatacsigtimewise[hippoatacsigtimewise$feature_ID %in% ourpeak &
                                                             hippoatacsigtimewise$sex %in% "male" &
                                                             hippoatacsigtimewise$comparison_group %in% "8w","logFC"]
  if(hippoatacsigupordowndf[i,"W8.F"] > 0.1 & hippoatacsigupordowndf[i,"W8.M"] > 0.1){
    hippoatacsigupordowndf[i,"UpReg"] <- 1  
  }
  if(hippoatacsigupordowndf[i,"W8.F"] < -0.1 & hippoatacsigupordowndf[i,"W8.M"] < -0.1){
    hippoatacsigupordowndf[i,"DownReg"] <- 1  
  }
}

# kidney
kidneyatacsigtimewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$feature_ID %in% kidneyatacsig & epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "KIDNEY",]
kidneyatacsigupordowndf <- data.frame(row.names = kidneyatacsig,
                                      "W8 F" = rep(0,length(kidneyatacsig)),
                                      "W8 M" = rep(0,length(kidneyatacsig)),
                                      "UpReg" = rep(0,length(kidneyatacsig)),
                                      "DownReg" = rep(0,length(kidneyatacsig)))
for(i in 1:length(kidneyatacsig)){
  if(i%%50 == 0){
    print(i)
  }
  ourpeak <- kidneyatacsig[i]
  kidneyatacsigupordowndf[i,"W8.F"] <- kidneyatacsigtimewise[kidneyatacsigtimewise$feature_ID %in% ourpeak &
                                                               kidneyatacsigtimewise$sex %in% "female" &
                                                               kidneyatacsigtimewise$comparison_group %in% "8w","logFC"]
  kidneyatacsigupordowndf[i,"W8.M"] <- kidneyatacsigtimewise[kidneyatacsigtimewise$feature_ID %in% ourpeak &
                                                               kidneyatacsigtimewise$sex %in% "male" &
                                                               kidneyatacsigtimewise$comparison_group %in% "8w","logFC"]
  if(kidneyatacsigupordowndf[i,"W8.F"] > 0.1 & kidneyatacsigupordowndf[i,"W8.M"] > 0.1){
    kidneyatacsigupordowndf[i,"UpReg"] <- 1  
  }
  if(kidneyatacsigupordowndf[i,"W8.F"] < -0.1 & kidneyatacsigupordowndf[i,"W8.M"] < -0.1){
    kidneyatacsigupordowndf[i,"DownReg"] <- 1  
  }
}

# liver
liveratacsigtimewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$feature_ID %in% liveratacsig & epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "LIVER",]
liveratacsigupordowndf <- data.frame(row.names = liveratacsig,
                                     "W8 F" = rep(0,length(liveratacsig)),
                                     "W8 M" = rep(0,length(liveratacsig)),
                                     "UpReg" = rep(0,length(liveratacsig)),
                                     "DownReg" = rep(0,length(liveratacsig)))
for(i in 1:length(liveratacsig)){
  if(i%%50 == 0){
    print(i)
  }
  ourpeak <- liveratacsig[i]
  liveratacsigupordowndf[i,"W8.F"] <- liveratacsigtimewise[liveratacsigtimewise$feature_ID %in% ourpeak &
                                                             liveratacsigtimewise$sex %in% "female" &
                                                             liveratacsigtimewise$comparison_group %in% "8w","logFC"]
  liveratacsigupordowndf[i,"W8.M"] <- liveratacsigtimewise[liveratacsigtimewise$feature_ID %in% ourpeak &
                                                             liveratacsigtimewise$sex %in% "male" &
                                                             liveratacsigtimewise$comparison_group %in% "8w","logFC"]
  if(liveratacsigupordowndf[i,"W8.F"] > 0.1 & liveratacsigupordowndf[i,"W8.M"] > 0.1){
    liveratacsigupordowndf[i,"UpReg"] <- 1  
  }
  if(liveratacsigupordowndf[i,"W8.F"] < -0.1 & liveratacsigupordowndf[i,"W8.M"] < -0.1){
    liveratacsigupordowndf[i,"DownReg"] <- 1  
  }
}

# lung
lungatacsigtimewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$feature_ID %in% lungatacsig & epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "LUNG",]
lungatacsigupordowndf <- data.frame(row.names = lungatacsig,
                                    "W8 F" = rep(0,length(lungatacsig)),
                                    "W8 M" = rep(0,length(lungatacsig)),
                                    "UpReg" = rep(0,length(lungatacsig)),
                                    "DownReg" = rep(0,length(lungatacsig)))
for(i in 1:length(lungatacsig)){
  if(i%%50 == 0){
    print(i)
  }
  ourpeak <- lungatacsig[i]
  lungatacsigupordowndf[i,"W8.F"] <- lungatacsigtimewise[lungatacsigtimewise$feature_ID %in% ourpeak &
                                                           lungatacsigtimewise$sex %in% "female" &
                                                           lungatacsigtimewise$comparison_group %in% "8w","logFC"]
  lungatacsigupordowndf[i,"W8.M"] <- lungatacsigtimewise[lungatacsigtimewise$feature_ID %in% ourpeak &
                                                           lungatacsigtimewise$sex %in% "male" &
                                                           lungatacsigtimewise$comparison_group %in% "8w","logFC"]
  if(lungatacsigupordowndf[i,"W8.F"] > 0.1 & lungatacsigupordowndf[i,"W8.M"] > 0.1){
    lungatacsigupordowndf[i,"UpReg"] <- 1  
  }
  if(lungatacsigupordowndf[i,"W8.F"] < -0.1 & lungatacsigupordowndf[i,"W8.M"] < -0.1){
    lungatacsigupordowndf[i,"DownReg"] <- 1  
  }
}

# brown
brownatacsigtimewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$feature_ID %in% brownatacsig & epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "BAT",]
brownatacsigupordowndf <- data.frame(row.names = brownatacsig,
                                     "W8 F" = rep(0,length(brownatacsig)),
                                     "W8 M" = rep(0,length(brownatacsig)),
                                     "UpReg" = rep(0,length(brownatacsig)),
                                     "DownReg" = rep(0,length(brownatacsig)))
for(i in 1:length(brownatacsig)){
  if(i%%50 == 0){
    print(i)
  }
  ourpeak <- brownatacsig[i]
  brownatacsigupordowndf[i,"W8.F"] <- brownatacsigtimewise[brownatacsigtimewise$feature_ID %in% ourpeak &
                                                             brownatacsigtimewise$sex %in% "female" &
                                                             brownatacsigtimewise$comparison_group %in% "8w","logFC"]
  brownatacsigupordowndf[i,"W8.M"] <- brownatacsigtimewise[brownatacsigtimewise$feature_ID %in% ourpeak &
                                                             brownatacsigtimewise$sex %in% "male" &
                                                             brownatacsigtimewise$comparison_group %in% "8w","logFC"]
  if(brownatacsigupordowndf[i,"W8.F"] > 0.1 & brownatacsigupordowndf[i,"W8.M"] > 0.1){
    brownatacsigupordowndf[i,"UpReg"] <- 1  
  }
  if(brownatacsigupordowndf[i,"W8.F"] < -0.1 & brownatacsigupordowndf[i,"W8.M"] < -0.1){
    brownatacsigupordowndf[i,"DownReg"] <- 1  
  }
}

# white
whiteatacsigtimewise <- epigen_atac_seq$timewise_dea[epigen_atac_seq$timewise_dea$feature_ID %in% whiteatacsig & epigen_atac_seq$timewise_dea$tissue_abbreviation %in% "WAT-SC",]
whiteatacsigupordowndf <- data.frame(row.names = whiteatacsig,
                                     "W8 F" = rep(0,length(whiteatacsig)),
                                     "W8 M" = rep(0,length(whiteatacsig)),
                                     "UpReg" = rep(0,length(whiteatacsig)),
                                     "DownReg" = rep(0,length(whiteatacsig)))
for(i in 1:length(whiteatacsig)){
  if(i%%50 == 0){
    print(i)
  }
  ourpeak <- whiteatacsig[i]
  whiteatacsigupordowndf[i,"W8.F"] <- whiteatacsigtimewise[whiteatacsigtimewise$feature_ID %in% ourpeak &
                                                             whiteatacsigtimewise$sex %in% "female" &
                                                             whiteatacsigtimewise$comparison_group %in% "8w","logFC"]
  whiteatacsigupordowndf[i,"W8.M"] <- whiteatacsigtimewise[whiteatacsigtimewise$feature_ID %in% ourpeak &
                                                             whiteatacsigtimewise$sex %in% "male" &
                                                             whiteatacsigtimewise$comparison_group %in% "8w","logFC"]
  if(whiteatacsigupordowndf[i,"W8.F"] > 0.1 & whiteatacsigupordowndf[i,"W8.M"] > 0.1){
    whiteatacsigupordowndf[i,"UpReg"] <- 1  
  }
  if(whiteatacsigupordowndf[i,"W8.F"] < -0.1 & whiteatacsigupordowndf[i,"W8.M"] < -0.1){
    whiteatacsigupordowndf[i,"DownReg"] <- 1  
  }
}

gastroatacsigupordowndf$W8Reg <- -1*gastroatacsigupordowndf$DownReg + gastroatacsigupordowndf$UpReg
heartatacsigupordowndf$W8Reg <- -1*heartatacsigupordowndf$DownReg + heartatacsigupordowndf$UpReg
hippoatacsigupordowndf$W8Reg <- -1*hippoatacsigupordowndf$DownReg + hippoatacsigupordowndf$UpReg
kidneyatacsigupordowndf$W8Reg <- -1*kidneyatacsigupordowndf$DownReg + kidneyatacsigupordowndf$UpReg
liveratacsigupordowndf$W8Reg <- -1*liveratacsigupordowndf$DownReg + liveratacsigupordowndf$UpReg
lungatacsigupordowndf$W8Reg <- -1*lungatacsigupordowndf$DownReg + lungatacsigupordowndf$UpReg
brownatacsigupordowndf$W8Reg <- -1*brownatacsigupordowndf$DownReg + brownatacsigupordowndf$UpReg
whiteatacsigupordowndf$W8Reg <- -1*whiteatacsigupordowndf$DownReg + whiteatacsigupordowndf$UpReg

save(list = c("gastroatacsigupordowndf","heartatacsigupordowndf","hippoatacsigupordowndf","kidneyatacsigupordowndf",
              "liveratacsigupordowndf","lungatacsigupordowndf","brownatacsigupordowndf","whiteatacsigupordowndf"),file = "atacsigupordowndfs.RData")


# Now RRBS

load("omesigdata.RData")
load("PASS1B Transcription Factor Paper Data/DEA Analysis/pass1b-06_epigen-rrbs_dea.RData")

gastromethsigtimewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$feature_ID %in% gastromethsig & epigen_rrbs$timewise_dea$tissue_abbreviation %in% "SKM-GN",]
gastromethsigupordowndf <- data.frame(row.names = gastromethsig,
                                      "W8 F" = rep(0,length(gastromethsig)),
                                      "W8 M" = rep(0,length(gastromethsig)),
                                      "UpReg" = rep(0,length(gastromethsig)),
                                      "DownReg" = rep(0,length(gastromethsig)))
for(i in 1:length(gastromethsig)){
  if(i%%50 == 0){
    print(i)
  }
  ourpeak <- gastromethsig[i]
  gastromethsigupordowndf[i,"W8.F"] <- gastromethsigtimewise[gastromethsigtimewise$feature_ID %in% ourpeak &
                                                               gastromethsigtimewise$sex %in% "female" &
                                                               gastromethsigtimewise$comparison_group %in% "8w","logFC"]
  gastromethsigupordowndf[i,"W8.M"] <- gastromethsigtimewise[gastromethsigtimewise$feature_ID %in% ourpeak &
                                                               gastromethsigtimewise$sex %in% "male" &
                                                               gastromethsigtimewise$comparison_group %in% "8w","logFC"]
  if(gastromethsigupordowndf[i,"W8.F"] > 0.1 & gastromethsigupordowndf[i,"W8.M"] > 0.1){
    gastromethsigupordowndf[i,"UpReg"] <- 1  
  }
  if(gastromethsigupordowndf[i,"W8.F"] < -0.1 & gastromethsigupordowndf[i,"W8.M"] < -0.1){
    gastromethsigupordowndf[i,"DownReg"] <- 1  
  }
}

# heart
heartmethsigtimewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$feature_ID %in% heartmethsig & epigen_rrbs$timewise_dea$tissue_abbreviation %in% "HEART",]
heartmethsigupordowndf <- data.frame(row.names = heartmethsig,
                                     "W8 F" = rep(0,length(heartmethsig)),
                                     "W8 M" = rep(0,length(heartmethsig)),
                                     "UpReg" = rep(0,length(heartmethsig)),
                                     "DownReg" = rep(0,length(heartmethsig)))
for(i in 1:length(heartmethsig)){
  if(i%%50 == 0){
    print(i)
  }
  ourpeak <- heartmethsig[i]
  heartmethsigupordowndf[i,"W8.F"] <- heartmethsigtimewise[heartmethsigtimewise$feature_ID %in% ourpeak &
                                                             heartmethsigtimewise$sex %in% "female" &
                                                             heartmethsigtimewise$comparison_group %in% "8w","logFC"]
  heartmethsigupordowndf[i,"W8.M"] <- heartmethsigtimewise[heartmethsigtimewise$feature_ID %in% ourpeak &
                                                             heartmethsigtimewise$sex %in% "male" &
                                                             heartmethsigtimewise$comparison_group %in% "8w","logFC"]
  if(heartmethsigupordowndf[i,"W8.F"] > 0.1 & heartmethsigupordowndf[i,"W8.M"] > 0.1){
    heartmethsigupordowndf[i,"UpReg"] <- 1  
  }
  if(heartmethsigupordowndf[i,"W8.F"] < -0.1 & heartmethsigupordowndf[i,"W8.M"] < -0.1){
    heartmethsigupordowndf[i,"DownReg"] <- 1  
  }
}

# hippo
hippomethsigtimewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$feature_ID %in% hippomethsig & epigen_rrbs$timewise_dea$tissue_abbreviation %in% "HIPPOC",]
hippomethsigupordowndf <- data.frame(row.names = hippomethsig,
                                     "W8 F" = rep(0,length(hippomethsig)),
                                     "W8 M" = rep(0,length(hippomethsig)),
                                     "UpReg" = rep(0,length(hippomethsig)),
                                     "DownReg" = rep(0,length(hippomethsig)))
for(i in 1:length(hippomethsig)){
  if(i%%50 == 0){
    print(i)
  }
  ourpeak <- hippomethsig[i]
  hippomethsigupordowndf[i,"W8.F"] <- hippomethsigtimewise[hippomethsigtimewise$feature_ID %in% ourpeak &
                                                             hippomethsigtimewise$sex %in% "female" &
                                                             hippomethsigtimewise$comparison_group %in% "8w","logFC"]
  hippomethsigupordowndf[i,"W8.M"] <- hippomethsigtimewise[hippomethsigtimewise$feature_ID %in% ourpeak &
                                                             hippomethsigtimewise$sex %in% "male" &
                                                             hippomethsigtimewise$comparison_group %in% "8w","logFC"]
  if(hippomethsigupordowndf[i,"W8.F"] > 0.1 & hippomethsigupordowndf[i,"W8.M"] > 0.1){
    hippomethsigupordowndf[i,"UpReg"] <- 1  
  }
  if(hippomethsigupordowndf[i,"W8.F"] < -0.1 & hippomethsigupordowndf[i,"W8.M"] < -0.1){
    hippomethsigupordowndf[i,"DownReg"] <- 1  
  }
}

# kidney
kidneymethsigtimewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$feature_ID %in% kidneymethsig & epigen_rrbs$timewise_dea$tissue_abbreviation %in% "KIDNEY",]
kidneymethsigupordowndf <- data.frame(row.names = kidneymethsig,
                                      "W8 F" = rep(0,length(kidneymethsig)),
                                      "W8 M" = rep(0,length(kidneymethsig)),
                                      "UpReg" = rep(0,length(kidneymethsig)),
                                      "DownReg" = rep(0,length(kidneymethsig)))
for(i in 1:length(kidneymethsig)){
  if(i%%50 == 0){
    print(i)
  }
  ourpeak <- kidneymethsig[i]
  kidneymethsigupordowndf[i,"W8.F"] <- kidneymethsigtimewise[kidneymethsigtimewise$feature_ID %in% ourpeak &
                                                               kidneymethsigtimewise$sex %in% "female" &
                                                               kidneymethsigtimewise$comparison_group %in% "8w","logFC"]
  kidneymethsigupordowndf[i,"W8.M"] <- kidneymethsigtimewise[kidneymethsigtimewise$feature_ID %in% ourpeak &
                                                               kidneymethsigtimewise$sex %in% "male" &
                                                               kidneymethsigtimewise$comparison_group %in% "8w","logFC"]
  if(kidneymethsigupordowndf[i,"W8.F"] > 0.1 & kidneymethsigupordowndf[i,"W8.M"] > 0.1){
    kidneymethsigupordowndf[i,"UpReg"] <- 1  
  }
  if(kidneymethsigupordowndf[i,"W8.F"] < -0.1 & kidneymethsigupordowndf[i,"W8.M"] < -0.1){
    kidneymethsigupordowndf[i,"DownReg"] <- 1  
  }
}

# liver
livermethsigtimewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$feature_ID %in% livermethsig & epigen_rrbs$timewise_dea$tissue_abbreviation %in% "LIVER",]
livermethsigupordowndf <- data.frame(row.names = livermethsig,
                                     "W8 F" = rep(0,length(livermethsig)),
                                     "W8 M" = rep(0,length(livermethsig)),
                                     "UpReg" = rep(0,length(livermethsig)),
                                     "DownReg" = rep(0,length(livermethsig)))
for(i in 1:length(livermethsig)){
  if(i%%50 == 0){
    print(i)
  }
  ourpeak <- livermethsig[i]
  livermethsigupordowndf[i,"W8.F"] <- livermethsigtimewise[livermethsigtimewise$feature_ID %in% ourpeak &
                                                             livermethsigtimewise$sex %in% "female" &
                                                             livermethsigtimewise$comparison_group %in% "8w","logFC"]
  livermethsigupordowndf[i,"W8.M"] <- livermethsigtimewise[livermethsigtimewise$feature_ID %in% ourpeak &
                                                             livermethsigtimewise$sex %in% "male" &
                                                             livermethsigtimewise$comparison_group %in% "8w","logFC"]
  if(livermethsigupordowndf[i,"W8.F"] > 0.1 & livermethsigupordowndf[i,"W8.M"] > 0.1){
    livermethsigupordowndf[i,"UpReg"] <- 1  
  }
  if(livermethsigupordowndf[i,"W8.F"] < -0.1 & livermethsigupordowndf[i,"W8.M"] < -0.1){
    livermethsigupordowndf[i,"DownReg"] <- 1  
  }
}

# lung
lungmethsigtimewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$feature_ID %in% lungmethsig & epigen_rrbs$timewise_dea$tissue_abbreviation %in% "LUNG",]
lungmethsigupordowndf <- data.frame(row.names = lungmethsig,
                                    "W8 F" = rep(0,length(lungmethsig)),
                                    "W8 M" = rep(0,length(lungmethsig)),
                                    "UpReg" = rep(0,length(lungmethsig)),
                                    "DownReg" = rep(0,length(lungmethsig)))
for(i in 1:length(lungmethsig)){
  if(i%%50 == 0){
    print(i)
  }
  ourpeak <- lungmethsig[i]
  lungmethsigupordowndf[i,"W8.F"] <- lungmethsigtimewise[lungmethsigtimewise$feature_ID %in% ourpeak &
                                                           lungmethsigtimewise$sex %in% "female" &
                                                           lungmethsigtimewise$comparison_group %in% "8w","logFC"]
  lungmethsigupordowndf[i,"W8.M"] <- lungmethsigtimewise[lungmethsigtimewise$feature_ID %in% ourpeak &
                                                           lungmethsigtimewise$sex %in% "male" &
                                                           lungmethsigtimewise$comparison_group %in% "8w","logFC"]
  if(lungmethsigupordowndf[i,"W8.F"] > 0.1 & lungmethsigupordowndf[i,"W8.M"] > 0.1){
    lungmethsigupordowndf[i,"UpReg"] <- 1  
  }
  if(lungmethsigupordowndf[i,"W8.F"] < -0.1 & lungmethsigupordowndf[i,"W8.M"] < -0.1){
    lungmethsigupordowndf[i,"DownReg"] <- 1  
  }
}

# brown
brownmethsigtimewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$feature_ID %in% brownmethsig & epigen_rrbs$timewise_dea$tissue_abbreviation %in% "BAT",]
brownmethsigupordowndf <- data.frame(row.names = brownmethsig,
                                     "W8 F" = rep(0,length(brownmethsig)),
                                     "W8 M" = rep(0,length(brownmethsig)),
                                     "UpReg" = rep(0,length(brownmethsig)),
                                     "DownReg" = rep(0,length(brownmethsig)))
for(i in 1:length(brownmethsig)){
  if(i%%50 == 0){
    print(i)
  }
  ourpeak <- brownmethsig[i]
  brownmethsigupordowndf[i,"W8.F"] <- brownmethsigtimewise[brownmethsigtimewise$feature_ID %in% ourpeak &
                                                             brownmethsigtimewise$sex %in% "female" &
                                                             brownmethsigtimewise$comparison_group %in% "8w","logFC"]
  brownmethsigupordowndf[i,"W8.M"] <- brownmethsigtimewise[brownmethsigtimewise$feature_ID %in% ourpeak &
                                                             brownmethsigtimewise$sex %in% "male" &
                                                             brownmethsigtimewise$comparison_group %in% "8w","logFC"]
  if(brownmethsigupordowndf[i,"W8.F"] > 0.1 & brownmethsigupordowndf[i,"W8.M"] > 0.1){
    brownmethsigupordowndf[i,"UpReg"] <- 1  
  }
  if(brownmethsigupordowndf[i,"W8.F"] < -0.1 & brownmethsigupordowndf[i,"W8.M"] < -0.1){
    brownmethsigupordowndf[i,"DownReg"] <- 1  
  }
}

# white
whitemethsigtimewise <- epigen_rrbs$timewise_dea[epigen_rrbs$timewise_dea$feature_ID %in% whitemethsig & epigen_rrbs$timewise_dea$tissue_abbreviation %in% "WAT-SC",]
whitemethsigupordowndf <- data.frame(row.names = whitemethsig,
                                     "W8 F" = rep(0,length(whitemethsig)),
                                     "W8 M" = rep(0,length(whitemethsig)),
                                     "UpReg" = rep(0,length(whitemethsig)),
                                     "DownReg" = rep(0,length(whitemethsig)))
for(i in 1:length(whitemethsig)){
  if(i%%50 == 0){
    print(i)
  }
  ourpeak <- whitemethsig[i]
  whitemethsigupordowndf[i,"W8.F"] <- whitemethsigtimewise[whitemethsigtimewise$feature_ID %in% ourpeak &
                                                             whitemethsigtimewise$sex %in% "female" &
                                                             whitemethsigtimewise$comparison_group %in% "8w","logFC"]
  whitemethsigupordowndf[i,"W8.M"] <- whitemethsigtimewise[whitemethsigtimewise$feature_ID %in% ourpeak &
                                                             whitemethsigtimewise$sex %in% "male" &
                                                             whitemethsigtimewise$comparison_group %in% "8w","logFC"]
  if(whitemethsigupordowndf[i,"W8.F"] > 0.1 & whitemethsigupordowndf[i,"W8.M"] > 0.1){
    whitemethsigupordowndf[i,"UpReg"] <- 1  
  }
  if(whitemethsigupordowndf[i,"W8.F"] < -0.1 & whitemethsigupordowndf[i,"W8.M"] < -0.1){
    whitemethsigupordowndf[i,"DownReg"] <- 1  
  }
}

gastromethsigupordowndf$W8Reg <- -1*gastromethsigupordowndf$DownReg + gastromethsigupordowndf$UpReg
heartmethsigupordowndf$W8Reg <- -1*heartmethsigupordowndf$DownReg + heartmethsigupordowndf$UpReg
hippomethsigupordowndf$W8Reg <- -1*hippomethsigupordowndf$DownReg + hippomethsigupordowndf$UpReg
kidneymethsigupordowndf$W8Reg <- -1*kidneymethsigupordowndf$DownReg + kidneymethsigupordowndf$UpReg
livermethsigupordowndf$W8Reg <- -1*livermethsigupordowndf$DownReg + livermethsigupordowndf$UpReg
lungmethsigupordowndf$W8Reg <- -1*lungmethsigupordowndf$DownReg + lungmethsigupordowndf$UpReg
brownmethsigupordowndf$W8Reg <- -1*brownmethsigupordowndf$DownReg + brownmethsigupordowndf$UpReg
whitemethsigupordowndf$W8Reg <- -1*whitemethsigupordowndf$DownReg + whitemethsigupordowndf$UpReg

save(list = c("gastromethsigupordowndf","heartmethsigupordowndf","hippomethsigupordowndf","kidneymethsigupordowndf",
              "livermethsigupordowndf","lungmethsigupordowndf","brownmethsigupordowndf","whitemethsigupordowndf"),file = "methsigupordowndfs.RData")


# And now RNAseq

load("omesigdata.RData")
load("PASS1B Transcription Factor Paper Data/DEA Analysis/pass1b-06_transcript-rna-seq_dea.RData")

# gastro

gastrornasigtimewise <- transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$feature_ID %in% gastrornasig & transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "SKM-GN",]
gastrornasigupordowndf <- data.frame(row.names = gastrornasig,
                                      "W8 F" = rep(0,length(gastrornasig)),
                                      "W8 M" = rep(0,length(gastrornasig)),
                                      "UpReg" = rep(0,length(gastrornasig)),
                                      "DownReg" = rep(0,length(gastrornasig)))
for(i in 1:length(gastrornasig)){
  if(i%%50 == 0){
    print(i)
  }
  ourpeak <- gastrornasig[i]
  gastrornasigupordowndf[i,"W8.F"] <- gastrornasigtimewise[gastrornasigtimewise$feature_ID %in% ourpeak &
                                                               gastrornasigtimewise$sex %in% "female" &
                                                               gastrornasigtimewise$comparison_group %in% "8w","logFC"]
  gastrornasigupordowndf[i,"W8.M"] <- gastrornasigtimewise[gastrornasigtimewise$feature_ID %in% ourpeak &
                                                               gastrornasigtimewise$sex %in% "male" &
                                                               gastrornasigtimewise$comparison_group %in% "8w","logFC"]
  if(gastrornasigupordowndf[i,"W8.F"] > 0.1 & gastrornasigupordowndf[i,"W8.M"] > 0.1){
    gastrornasigupordowndf[i,"UpReg"] <- 1  
  }
  if(gastrornasigupordowndf[i,"W8.F"] < -0.1 & gastrornasigupordowndf[i,"W8.M"] < -0.1){
    gastrornasigupordowndf[i,"DownReg"] <- 1  
  }
}


# heart

heartrnasigtimewise <- transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$feature_ID %in% heartrnasig & transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "HEART",]
heartrnasigupordowndf <- data.frame(row.names = heartrnasig,
                                     "W8 F" = rep(0,length(heartrnasig)),
                                     "W8 M" = rep(0,length(heartrnasig)),
                                     "UpReg" = rep(0,length(heartrnasig)),
                                     "DownReg" = rep(0,length(heartrnasig)))
for(i in 1:length(heartrnasig)){
  if(i%%50 == 0){
    print(i)
  }
  ourpeak <- heartrnasig[i]
  heartrnasigupordowndf[i,"W8.F"] <- heartrnasigtimewise[heartrnasigtimewise$feature_ID %in% ourpeak &
                                                             heartrnasigtimewise$sex %in% "female" &
                                                             heartrnasigtimewise$comparison_group %in% "8w","logFC"]
  heartrnasigupordowndf[i,"W8.M"] <- heartrnasigtimewise[heartrnasigtimewise$feature_ID %in% ourpeak &
                                                             heartrnasigtimewise$sex %in% "male" &
                                                             heartrnasigtimewise$comparison_group %in% "8w","logFC"]
  if(heartrnasigupordowndf[i,"W8.F"] > 0.1 & heartrnasigupordowndf[i,"W8.M"] > 0.1){
    heartrnasigupordowndf[i,"UpReg"] <- 1  
  }
  if(heartrnasigupordowndf[i,"W8.F"] < -0.1 & heartrnasigupordowndf[i,"W8.M"] < -0.1){
    heartrnasigupordowndf[i,"DownReg"] <- 1  
  }
}


# hippo

hippornasigtimewise <- transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$feature_ID %in% hippornasig & transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "HIPPOC",]
hippornasigupordowndf <- data.frame(row.names = hippornasig,
                                     "W8 F" = rep(0,length(hippornasig)),
                                     "W8 M" = rep(0,length(hippornasig)),
                                     "UpReg" = rep(0,length(hippornasig)),
                                     "DownReg" = rep(0,length(hippornasig)))
for(i in 1:length(hippornasig)){
  if(i%%50 == 0){
    print(i)
  }
  ourpeak <- hippornasig[i]
  hippornasigupordowndf[i,"W8.F"] <- hippornasigtimewise[hippornasigtimewise$feature_ID %in% ourpeak &
                                                             hippornasigtimewise$sex %in% "female" &
                                                             hippornasigtimewise$comparison_group %in% "8w","logFC"]
  hippornasigupordowndf[i,"W8.M"] <- hippornasigtimewise[hippornasigtimewise$feature_ID %in% ourpeak &
                                                             hippornasigtimewise$sex %in% "male" &
                                                             hippornasigtimewise$comparison_group %in% "8w","logFC"]
  if(hippornasigupordowndf[i,"W8.F"] > 0.1 & hippornasigupordowndf[i,"W8.M"] > 0.1){
    hippornasigupordowndf[i,"UpReg"] <- 1  
  }
  if(hippornasigupordowndf[i,"W8.F"] < -0.1 & hippornasigupordowndf[i,"W8.M"] < -0.1){
    hippornasigupordowndf[i,"DownReg"] <- 1  
  }
}


# kidney

kidneyrnasigtimewise <- transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$feature_ID %in% kidneyrnasig & transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "KIDNEY",]
kidneyrnasigupordowndf <- data.frame(row.names = kidneyrnasig,
                                     "W8 F" = rep(0,length(kidneyrnasig)),
                                     "W8 M" = rep(0,length(kidneyrnasig)),
                                     "UpReg" = rep(0,length(kidneyrnasig)),
                                     "DownReg" = rep(0,length(kidneyrnasig)))
for(i in 1:length(kidneyrnasig)){
  if(i%%50 == 0){
    print(i)
  }
  ourpeak <- kidneyrnasig[i]
  kidneyrnasigupordowndf[i,"W8.F"] <- kidneyrnasigtimewise[kidneyrnasigtimewise$feature_ID %in% ourpeak &
                                                             kidneyrnasigtimewise$sex %in% "female" &
                                                             kidneyrnasigtimewise$comparison_group %in% "8w","logFC"]
  kidneyrnasigupordowndf[i,"W8.M"] <- kidneyrnasigtimewise[kidneyrnasigtimewise$feature_ID %in% ourpeak &
                                                             kidneyrnasigtimewise$sex %in% "male" &
                                                             kidneyrnasigtimewise$comparison_group %in% "8w","logFC"]
  if(kidneyrnasigupordowndf[i,"W8.F"] > 0.1 & kidneyrnasigupordowndf[i,"W8.M"] > 0.1){
    kidneyrnasigupordowndf[i,"UpReg"] <- 1  
  }
  if(kidneyrnasigupordowndf[i,"W8.F"] < -0.1 & kidneyrnasigupordowndf[i,"W8.M"] < -0.1){
    kidneyrnasigupordowndf[i,"DownReg"] <- 1  
  }
}


# liver

liverrnasigtimewise <- transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$feature_ID %in% liverrnasig & transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "LIVER",]
liverrnasigupordowndf <- data.frame(row.names = liverrnasig,
                                     "W8 F" = rep(0,length(liverrnasig)),
                                     "W8 M" = rep(0,length(liverrnasig)),
                                     "UpReg" = rep(0,length(liverrnasig)),
                                     "DownReg" = rep(0,length(liverrnasig)))
for(i in 1:length(liverrnasig)){
  if(i%%50 == 0){
    print(i)
  }
  ourpeak <- liverrnasig[i]
  liverrnasigupordowndf[i,"W8.F"] <- liverrnasigtimewise[liverrnasigtimewise$feature_ID %in% ourpeak &
                                                             liverrnasigtimewise$sex %in% "female" &
                                                             liverrnasigtimewise$comparison_group %in% "8w","logFC"]
  liverrnasigupordowndf[i,"W8.M"] <- liverrnasigtimewise[liverrnasigtimewise$feature_ID %in% ourpeak &
                                                             liverrnasigtimewise$sex %in% "male" &
                                                             liverrnasigtimewise$comparison_group %in% "8w","logFC"]
  if(liverrnasigupordowndf[i,"W8.F"] > 0.1 & liverrnasigupordowndf[i,"W8.M"] > 0.1){
    liverrnasigupordowndf[i,"UpReg"] <- 1  
  }
  if(liverrnasigupordowndf[i,"W8.F"] < -0.1 & liverrnasigupordowndf[i,"W8.M"] < -0.1){
    liverrnasigupordowndf[i,"DownReg"] <- 1  
  }
}


# lung

lungrnasigtimewise <- transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$feature_ID %in% lungrnasig & transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "LUNG",]
lungrnasigupordowndf <- data.frame(row.names = lungrnasig,
                                     "W8 F" = rep(0,length(lungrnasig)),
                                     "W8 M" = rep(0,length(lungrnasig)),
                                     "UpReg" = rep(0,length(lungrnasig)),
                                     "DownReg" = rep(0,length(lungrnasig)))
for(i in 1:length(lungrnasig)){
  if(i%%50 == 0){
    print(i)
  }
  ourpeak <- lungrnasig[i]
  lungrnasigupordowndf[i,"W8.F"] <- lungrnasigtimewise[lungrnasigtimewise$feature_ID %in% ourpeak &
                                                             lungrnasigtimewise$sex %in% "female" &
                                                             lungrnasigtimewise$comparison_group %in% "8w","logFC"]
  lungrnasigupordowndf[i,"W8.M"] <- lungrnasigtimewise[lungrnasigtimewise$feature_ID %in% ourpeak &
                                                             lungrnasigtimewise$sex %in% "male" &
                                                             lungrnasigtimewise$comparison_group %in% "8w","logFC"]
  if(is.na(lungrnasigupordowndf[i,"W8.M"])){
    lungrnasigupordowndf[i,"W8.M"] <- 0
  }
  if(is.na(lungrnasigupordowndf[i,"W8.F"])){
    lungrnasigupordowndf[i,"W8.F"] <- 0
  }
  
  if(lungrnasigupordowndf[i,"W8.F"] > 0.1 & lungrnasigupordowndf[i,"W8.M"] > 0.1){
    lungrnasigupordowndf[i,"UpReg"] <- 1  
  }
  if(lungrnasigupordowndf[i,"W8.F"] < -0.1 & lungrnasigupordowndf[i,"W8.M"] < -0.1){
    lungrnasigupordowndf[i,"DownReg"] <- 1  
  }
}


# brown

brownrnasigtimewise <- transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$feature_ID %in% brownrnasig & transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "BAT",]
brownrnasigupordowndf <- data.frame(row.names = brownrnasig,
                                     "W8 F" = rep(0,length(brownrnasig)),
                                     "W8 M" = rep(0,length(brownrnasig)),
                                     "UpReg" = rep(0,length(brownrnasig)),
                                     "DownReg" = rep(0,length(brownrnasig)))
for(i in 1:length(brownrnasig)){
  if(i%%50 == 0){
    print(i)
  }
  ourpeak <- brownrnasig[i]
  brownrnasigupordowndf[i,"W8.F"] <- brownrnasigtimewise[brownrnasigtimewise$feature_ID %in% ourpeak &
                                                             brownrnasigtimewise$sex %in% "female" &
                                                             brownrnasigtimewise$comparison_group %in% "8w","logFC"]
  brownrnasigupordowndf[i,"W8.M"] <- brownrnasigtimewise[brownrnasigtimewise$feature_ID %in% ourpeak &
                                                             brownrnasigtimewise$sex %in% "male" &
                                                             brownrnasigtimewise$comparison_group %in% "8w","logFC"]
  if(brownrnasigupordowndf[i,"W8.F"] > 0.1 & brownrnasigupordowndf[i,"W8.M"] > 0.1){
    brownrnasigupordowndf[i,"UpReg"] <- 1  
  }
  if(brownrnasigupordowndf[i,"W8.F"] < -0.1 & brownrnasigupordowndf[i,"W8.M"] < -0.1){
    brownrnasigupordowndf[i,"DownReg"] <- 1  
  }
}


# white

whiternasigtimewise <- transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$feature_ID %in% whiternasig & transcript_rna_seq$timewise_dea$tissue_abbreviation %in% "WAT-SC",]
whiternasigupordowndf <- data.frame(row.names = whiternasig,
                                     "W8 F" = rep(0,length(whiternasig)),
                                     "W8 M" = rep(0,length(whiternasig)),
                                     "UpReg" = rep(0,length(whiternasig)),
                                     "DownReg" = rep(0,length(whiternasig)))
for(i in 1:length(whiternasig)){
  if(i%%50 == 0){
    print(i)
  }
  ourpeak <- whiternasig[i]
  whiternasigupordowndf[i,"W8.F"] <- whiternasigtimewise[whiternasigtimewise$feature_ID %in% ourpeak &
                                                             whiternasigtimewise$sex %in% "female" &
                                                             whiternasigtimewise$comparison_group %in% "8w","logFC"]
  whiternasigupordowndf[i,"W8.M"] <- whiternasigtimewise[whiternasigtimewise$feature_ID %in% ourpeak &
                                                             whiternasigtimewise$sex %in% "male" &
                                                             whiternasigtimewise$comparison_group %in% "8w","logFC"]
  if(is.na(whiternasigupordowndf[i,"W8.M"])){
    whiternasigupordowndf[i,"W8.M"] <- 0
  }
  if(is.na(whiternasigupordowndf[i,"W8.F"])){
    whiternasigupordowndf[i,"W8.F"] <- 0
  }
  
  if(whiternasigupordowndf[i,"W8.F"] > 0.1 & whiternasigupordowndf[i,"W8.M"] > 0.1){
    whiternasigupordowndf[i,"UpReg"] <- 1  
  }
  if(whiternasigupordowndf[i,"W8.F"] < -0.1 & whiternasigupordowndf[i,"W8.M"] < -0.1){
    whiternasigupordowndf[i,"DownReg"] <- 1  
  }
}


save(list = c("gastrornasigupordowndf","heartrnasigupordowndf","hippornasigupordowndf","kidneyrnasigupordowndf",
              "liverrnasigupordowndf","lungrnasigupordowndf","brownrnasigupordowndf","whiternasigupordowndf"),file = "rnasigupordowndfs.RData")


#####
# We are going to look into training time of day
####

ratsubjecttrainingtime <- data.frame(row.names = unique(MotrpacRatTraining6moData::PHENO$pid),
                                     "TrainingTime_Day1" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day2" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day3" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day4" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day5" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day6" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day7" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day8" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day9" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day10" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day11" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day12" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day13" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day14" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day15" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day16" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day17" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day18" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day19" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day20" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day21" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day22" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day23" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day24" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day25" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day26" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day27" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day28" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day29" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day30" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day31" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day32" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day33" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day34" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day35" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day36" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day37" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day38" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day39" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))),
                                     "TrainingTime_Day40" = rep("",length(unique(MotrpacRatTraining6moData::PHENO$pid))))

for(i in 1:dim(ratsubjecttrainingtime)[1]){
  ourpid <- rownames(ratsubjecttrainingtime)[i]
  ratsubjecttrainingtime[i,"TrainingTime_Day1"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day1time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day2"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day2time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day3"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day3time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day4"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day4time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day5"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day5time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day6"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day6time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day7"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day7time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day8"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day8time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day9"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day9time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day10"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day10time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day11"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day11time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day12"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day12time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day13"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day13time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day14"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day14time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day15"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day15time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day16"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day16time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day17"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day17time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day18"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day18time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day19"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day19time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day20"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day20time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day21"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day21time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day22"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day22time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day23"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day23time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day24"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day24time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day25"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day25time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day26"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day26time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day27"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day27time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day28"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day28time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day29"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day29time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day30"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day30time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day31"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day31time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day32"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day32time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day33"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day33time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day34"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day34time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day35"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day35time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day36"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day36time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day37"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day37time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day38"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day38time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day39"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day39time"][1]
  ratsubjecttrainingtime[i,"TrainingTime_Day40"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid,"training.day40time"][1]
}

ratsubjecttrainingtimenumeric <- ratsubjecttrainingtime
for(i in 1:dim(ratsubjecttrainingtimenumeric)[1]){
  for(j in 1:dim(ratsubjecttrainingtimenumeric)[2]){
    ratsubjecttrainingtimenumeric[i,j] <- as.numeric(gsub(":.*","",ratsubjecttrainingtime[i,j]))*60 + as.numeric(gsub(".*:","",ratsubjecttrainingtime[i,j]))
  }
}

ratsubjecttrainingtimenumeric[is.na(ratsubjecttrainingtimenumeric)] <- 0
ratsubjecttrainingtimenumericmat <- matrix(0L, nrow = dim(ratsubjecttrainingtimenumeric)[1],ncol = dim(ratsubjecttrainingtimenumeric)[2])
rownames(ratsubjecttrainingtimenumericmat) <- rownames(ratsubjecttrainingtimenumeric)
colnames(ratsubjecttrainingtimenumericmat) <- c(1:40)
for(i in 1:dim(ratsubjecttrainingtimenumericmat)[1]){
  for(j in 1:dim(ratsubjecttrainingtimenumericmat)[2]){
    ratsubjecttrainingtimenumericmat[i,j] <- as.numeric(ratsubjecttrainingtimenumeric[i,j])
  }
}

gastrornaviallabels <- colnames(MotrpacRatTraining6moData::TRNSCRPT_SKMGN_NORM_DATA)[5:54]
gastrornapids <- unique(MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$viallabel %in% gastrornaviallabels,"pid"])

trainingtimetable <- ratsubjecttrainingtime[rownames(ratsubjecttrainingtime) %in% gastrornapids,]
colnames(trainingtimetable) <- gsub("TrainingTime_","",colnames(trainingtimetable))
write.csv(trainingtimetable,file = "New_Supplemental_Table_S4_52224.csv")