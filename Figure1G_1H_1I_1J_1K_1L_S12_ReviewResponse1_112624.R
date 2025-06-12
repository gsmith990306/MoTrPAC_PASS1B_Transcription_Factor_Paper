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

setwd("C:/Users/gsmit/Documents/Mount Sinai/Sealfon Laboratory/MoTrPAC/PASS1B Raw Data")

load("PASS1B Transcription Factor Paper Data/DEA Analysis/pass1b-06_epigen-atac-seq_dea.RData")
peakanno <- readRDS("peakanno.RDS")

####
# Figure 1G-L
#####

# Add peak annotations to training data
atactraining <- epigen_atac_seq$training_dea
atactraining$custom_annotation <- ""
atactraining$custom_annotation <- peakanno[atactraining$feature_ID,"custom_annotation"]

tab = table(atactraining$tissue_abbreviation,atactraining$custom_annotation)

myorder = c("Upstream (<5kb)", "Promoter (<=1kb)","Promoter (1-2kb)", "5' UTR","Exon",
            "Intron", "3' UTR","Downstream (<5kb)", "Distal Intergenic","Overlaps Gene")

mycol = pal_d3()(length(myorder))
names(mycol) = myorder
tab2 = tab[,myorder]
tab2Perc = t(apply(t(tab2), 2, function(x){x*100/sum(x,na.rm=T)}))

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

# We want to identify the ATAC peaks that are consistently up or down-regulated at a given time point.
# For simplicity, let's investigate the peaks that are up (>0.25) in both male and female rats at week 8

load("omesigdata.RData")

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

# gastro
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

alltissueatacsigW8upordown <- rbind(c(72,480,111),
                                    c(14,72,5),
                                    c(2,36,0),
                                    c(54,193,85),
                                    c(80,1050,267),
                                    c(7,211,12),
                                    c(12,269,99),
                                    c(0,5,0))
rownames(alltissueatacsigW8upordown) <- c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC")
colnames(alltissueatacsigW8upordown) <- c("W8 Down-Regulated","Neither","W8 Up-Regulated")
tissuemeta <- data.frame(row.names = rownames(alltissueatacsigW8upordown),
                         "Tissue" = rownames(alltissueatacsigW8upordown))

# Figure Aside for Reviewer on Dividing DAR sets by up and down regulation
png(file = "Reviewer Response Fig 1A_112624.png", width=7.5, height=3,units = "in",res = 600)
pheatmap(t(alltissueatacsigW8upordown),cluster_rows = F,cluster_cols = F,color = colorpanel(100,"white","firebrick"),display_numbers = T,number_color = "black",number_format = "%.0f",fontsize = 15,annotation_col = tissuemeta,annotation_colors = ann_cols,show_colnames = F,border_color = NA)
dev.off()

annotationfrequencyacrosstissuesupordown <- matrix(0L,ncol = 10,nrow = 24)
colnames(annotationfrequencyacrosstissuesupordown) <- names(table(peakanno$custom_annotation))
rownames(annotationfrequencyacrosstissuesupordown) <- c("SKM-GN Total",
                                                        "SKM-GN UpReg",
                                                        "SKM-GN DownReg",
                                                        "HEART Total",
                                                        "HEART UpReg",
                                                        "HEART DownReg",
                                                        "HIPPOC Total",
                                                        "HIPPOC UpReg",
                                                        "HIPPOC DownReg",
                                                        "KIDNEY Total",
                                                        "KIDNEY UpReg",
                                                        "KIDNEY DownReg",
                                                        "LIVER Total",
                                                        "LIVER UpReg",
                                                        "LIVER DownReg",
                                                        "LUNG Total",
                                                        "LUNG UpReg",
                                                        "LUNG DownReg",
                                                        "BAT Total",
                                                        "BAT UpReg",
                                                        "BAT DownReg",
                                                        "WAT-SC Total",
                                                        "WAT-SC UpReg",
                                                        "WAT-SC DownReg")

for(i in 1:dim(annotationfrequencyacrosstissuesupordown)[2]){
  ouranno <- colnames(annotationfrequencyacrosstissuesupordown)[i]
  annotationfrequencyacrosstissuesupordown["SKM-GN Total",i] <- sum(peakanno[gastroatacsig,"custom_annotation"] %in% ouranno)/length(gastroatacsig)
  annotationfrequencyacrosstissuesupordown["HEART Total",i] <- sum(peakanno[heartatacsig,"custom_annotation"] %in% ouranno)/length(heartatacsig)
  annotationfrequencyacrosstissuesupordown["HIPPOC Total",i] <- sum(peakanno[hippoatacsig,"custom_annotation"] %in% ouranno)/length(hippoatacsig)
  annotationfrequencyacrosstissuesupordown["KIDNEY Total",i] <- sum(peakanno[kidneyatacsig,"custom_annotation"] %in% ouranno)/length(kidneyatacsig)
  annotationfrequencyacrosstissuesupordown["LIVER Total",i] <- sum(peakanno[liveratacsig,"custom_annotation"] %in% ouranno)/length(liveratacsig)
  annotationfrequencyacrosstissuesupordown["LUNG Total",i] <- sum(peakanno[lungatacsig,"custom_annotation"] %in% ouranno)/length(lungatacsig)
  annotationfrequencyacrosstissuesupordown["BAT Total",i] <- sum(peakanno[brownatacsig,"custom_annotation"] %in% ouranno)/length(brownatacsig)
  annotationfrequencyacrosstissuesupordown["WAT-SC Total",i] <- sum(peakanno[whiteatacsig,"custom_annotation"] %in% ouranno)/length(whiteatacsig)
  
  annotationfrequencyacrosstissuesupordown["SKM-GN UpReg",i] <- sum(peakanno[rownames(gastroatacsigupordowndf)[gastroatacsigupordowndf$UpReg == "1"],"custom_annotation"] %in% ouranno)/length(rownames(gastroatacsigupordowndf)[gastroatacsigupordowndf$UpReg == "1"])
  annotationfrequencyacrosstissuesupordown["SKM-GN DownReg",i] <- sum(peakanno[rownames(gastroatacsigupordowndf)[gastroatacsigupordowndf$DownReg == "1"],"custom_annotation"] %in% ouranno)/length(rownames(gastroatacsigupordowndf)[gastroatacsigupordowndf$DownReg == "1"])
  
  annotationfrequencyacrosstissuesupordown["HEART UpReg",i] <- sum(peakanno[rownames(heartatacsigupordowndf)[heartatacsigupordowndf$UpReg == "1"],"custom_annotation"] %in% ouranno)/length(rownames(heartatacsigupordowndf)[heartatacsigupordowndf$UpReg == "1"])
  annotationfrequencyacrosstissuesupordown["HEART DownReg",i] <- sum(peakanno[rownames(heartatacsigupordowndf)[heartatacsigupordowndf$DownReg == "1"],"custom_annotation"] %in% ouranno)/length(rownames(heartatacsigupordowndf)[heartatacsigupordowndf$DownReg == "1"])
  
  annotationfrequencyacrosstissuesupordown["HIPPOC UpReg",i] <- sum(peakanno[rownames(hippoatacsigupordowndf)[hippoatacsigupordowndf$UpReg == "1"],"custom_annotation"] %in% ouranno)/length(rownames(hippoatacsigupordowndf)[hippoatacsigupordowndf$UpReg == "1"])
  annotationfrequencyacrosstissuesupordown["HIPPOC DownReg",i] <- sum(peakanno[rownames(hippoatacsigupordowndf)[hippoatacsigupordowndf$DownReg == "1"],"custom_annotation"] %in% ouranno)/length(rownames(hippoatacsigupordowndf)[hippoatacsigupordowndf$DownReg == "1"])
  
  annotationfrequencyacrosstissuesupordown["KIDNEY UpReg",i] <- sum(peakanno[rownames(kidneyatacsigupordowndf)[kidneyatacsigupordowndf$UpReg == "1"],"custom_annotation"] %in% ouranno)/length(rownames(kidneyatacsigupordowndf)[kidneyatacsigupordowndf$UpReg == "1"])
  annotationfrequencyacrosstissuesupordown["KIDNEY DownReg",i] <- sum(peakanno[rownames(kidneyatacsigupordowndf)[kidneyatacsigupordowndf$DownReg == "1"],"custom_annotation"] %in% ouranno)/length(rownames(kidneyatacsigupordowndf)[kidneyatacsigupordowndf$DownReg == "1"])
  
  annotationfrequencyacrosstissuesupordown["LIVER UpReg",i] <- sum(peakanno[rownames(liveratacsigupordowndf)[liveratacsigupordowndf$UpReg == "1"],"custom_annotation"] %in% ouranno)/length(rownames(liveratacsigupordowndf)[liveratacsigupordowndf$UpReg == "1"])
  annotationfrequencyacrosstissuesupordown["LIVER DownReg",i] <- sum(peakanno[rownames(liveratacsigupordowndf)[liveratacsigupordowndf$DownReg == "1"],"custom_annotation"] %in% ouranno)/length(rownames(liveratacsigupordowndf)[liveratacsigupordowndf$DownReg == "1"])
  
  annotationfrequencyacrosstissuesupordown["LUNG UpReg",i] <- sum(peakanno[rownames(lungatacsigupordowndf)[lungatacsigupordowndf$UpReg == "1"],"custom_annotation"] %in% ouranno)/length(rownames(lungatacsigupordowndf)[lungatacsigupordowndf$UpReg == "1"])
  annotationfrequencyacrosstissuesupordown["LUNG DownReg",i] <- sum(peakanno[rownames(lungatacsigupordowndf)[lungatacsigupordowndf$DownReg == "1"],"custom_annotation"] %in% ouranno)/length(rownames(lungatacsigupordowndf)[lungatacsigupordowndf$DownReg == "1"])
  
  annotationfrequencyacrosstissuesupordown["BAT UpReg",i] <- sum(peakanno[rownames(brownatacsigupordowndf)[brownatacsigupordowndf$UpReg == "1"],"custom_annotation"] %in% ouranno)/length(rownames(brownatacsigupordowndf)[brownatacsigupordowndf$UpReg == "1"])
  annotationfrequencyacrosstissuesupordown["BAT DownReg",i] <- sum(peakanno[rownames(brownatacsigupordowndf)[brownatacsigupordowndf$DownReg == "1"],"custom_annotation"] %in% ouranno)/length(rownames(brownatacsigupordowndf)[brownatacsigupordowndf$DownReg == "1"])
  
  annotationfrequencyacrosstissuesupordown["WAT-SC UpReg",i] <- sum(peakanno[rownames(whiteatacsigupordowndf)[whiteatacsigupordowndf$UpReg == "1"],"custom_annotation"] %in% ouranno)/length(rownames(whiteatacsigupordowndf)[whiteatacsigupordowndf$UpReg == "1"])
  annotationfrequencyacrosstissuesupordown["WAT-SC DownReg",i] <- sum(peakanno[rownames(whiteatacsigupordowndf)[whiteatacsigupordowndf$DownReg == "1"],"custom_annotation"] %in% ouranno)/length(rownames(whiteatacsigupordowndf)[whiteatacsigupordowndf$DownReg == "1"])
  
}

png(file = "Reviewer Response Fig 1B_112624.png",width = 14,height = 6,units = "in",res = 600)
pheatmap(annotationfrequencyacrosstissuesupordown,cluster_rows = F,cluster_cols = F,color = colorpanel(101,"white","firebrick"),display_numbers = T,number_color = "black",angle_col = 0,border_color = NA)
dev.off()

# Figure 1G
pdf(file = "Figure 1G_112624.pdf", width=8, height=6)
mosaicplot(tab2Perc,las=1,col=genomic_features_col$Region[colnames(tab2Perc)], main="Genomic Distribution of All Peaks",cex.axis = 0.8)
dev.off()

png(file = "Figure 1G_112624.png", width=8, height=6,units = "in",res = 600)
mosaicplot(tab2Perc,las=1,col=genomic_features_col$Region[colnames(tab2Perc)], main="Genomic Distribution of All Peaks")
dev.off()

pcutoff = 0.1
trainingSig = atactraining %>% filter(adj_p_value<pcutoff)

tab = table(trainingSig$tissue_abbreviation,trainingSig$custom_annotation)
tab3 = tab[,myorder[myorder %in% colnames(tab)]]
tab3Perc = t(apply(t(tab3), 2, function(x){x*100/sum(x,na.rm=T)}))

# Figure 1H
pdf(file = "Figure 1H_112624.pdf", width=8, height=6)
mosaicplot(tab3Perc,las=1,col=genomic_features_col$Region[colnames(tab3Perc)], main="Genomic Distribution of Significant Peaks",cex.axis = 0.8)
dev.off()

png(file = "Figure 1H_112624.png", width=8, height=6,units = "in",res = 600)
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

pdf(file = "Figure 1K_112624.pdf", width=8, height=3)
pheatmap(tab4[,c(1:7)],breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"steelblue","white","firebrick"),cluster_cols = F,cluster_rows = F,border_color = NA,legend = F,angle_col = 0)
dev.off()

png(file = "Figure 1K_112624.png", width=8, height=3,units = "in",res = 600)
pheatmap(tab4[,c(1:7)],breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"steelblue","white","firebrick"),cluster_cols = F,cluster_rows = F,border_color = NA,legend = F,angle_col = 0)
dev.off()

load("rnanormmatrices.RData")
load("omesigdata.RData")

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
# Supplemental Figure S12A
pdf(file = "Supplemental Figure S12A_112624.pdf", width = 7, height = 6)
mosaicplot(tabPerc, las = 1, col = genomic_features_col$Region[colnames(tabPerc)], main = "Genomic Distribution of Active Genes Percentage")
dev.off()

png(file = "Supplemental Figure S12A_112624.png", width = 7, height = 6,units = "in",res = 600)
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
# Supplemental Figure S12B
pdf(file = "Supplemental Figure S12B_112624.pdf", width = 7, height = 6)
mosaicplot(tab3Perc, las = 1, col = genomic_features_col$Region[colnames(tab3Perc)], main = "Genomic Distribution of DEGaP Percentage")
dev.off()

png(file = "Supplemental Figure S12B_112624.png", width = 7, height = 6,units = "in",res = 600)
mosaicplot(tab3Perc, las = 1, col = genomic_features_col$Region[colnames(tab3Perc)], main = "Genomic Distribution of DEGaP Percentage")
dev.off()

# Supplemental Figure S7C
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

pdf(file = "Supplemental Figure S12C_112624.pdf", width=8, height=3)
pheatmap(tab4,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"steelblue","white","firebrick"),cluster_cols = F,cluster_rows = F,border_color = NA,legend = F,angle_col = 0)
dev.off()

png(file = "Supplemental Figure S12C_112624.png", width=8, height=3,units = "in",res = 600)
pheatmap(tab4,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"steelblue","white","firebrick"),cluster_cols = F,cluster_rows = F,border_color = NA,legend = F,angle_col = 0)
dev.off()

rm(epigen_atac_seq)
rm(atactraining)

gc()

# Figure 1I

load("PASS1B Transcription Factor Paper Data/DEA Analysis/pass1b-06_epigen-rrbs_dea.RData")

newmethanno <- readRDS("newmethanno.RDS")

# Start using information here

trimmedmethanno <- newmethanno[!duplicated(newmethanno$feature_ID),]
rownames(trimmedmethanno) <- trimmedmethanno$feature_ID

# Add peak annotations to training data
methtraining <- epigen_rrbs$training_dea
methtraining$custom_annotation <- ""
methtraining$custom_annotation <- trimmedmethanno[methtraining$feature_ID,"custom_annotation"]

tab = table(methtraining$tissue_abbreviation,methtraining$custom_annotation)

myorder = c("Upstream (<5kb)", "Promoter (<=1kb)","Promoter (1-2kb)", "5' UTR","Exon",
            "Intron", "3' UTR","Downstream (<5kb)", "Distal Intergenic","Overlaps Gene")

mycol = pal_d3()(length(myorder))
names(mycol) = myorder
tab2 = tab[,myorder]
tab2Perc = t(apply(t(tab2), 2, function(x){x*100/sum(x,na.rm=T)}))

# Figure 1I
pdf(file = "Figure 1I_112624.pdf", width=8, height=6)
mosaicplot(tab2Perc,las=1,col=genomic_features_col$Region[colnames(tab2Perc)], main="Genomic Distribution of All Methyl Sites",cex.axis = 0.8)
dev.off()

png(file = "Figure 1I_112624.png", width=8, height=6,units = "in",res = 600)
mosaicplot(tab2Perc,las=1,col=genomic_features_col$Region[colnames(tab2Perc)], main="Genomic Distribution of All Methyl Sites")
dev.off()

pcutoff = 0.1
trainingSig = methtraining %>% filter(adj_p_value<pcutoff)

tab = table(trainingSig$tissue_abbreviation,trainingSig$custom_annotation)
tab3 = tab[,myorder[myorder %in% colnames(tab)]]
tab3Perc = t(apply(t(tab3), 2, function(x){x*100/sum(x,na.rm=T)}))

# Figure 1J
pdf(file = "Figure 1J_112624.pdf", width=8, height=6)
mosaicplot(tab3Perc,las=1,col=genomic_features_col$Region[colnames(tab3Perc)], main="Genomic Distribution of Significant Methyl Sites",cex.axis = 0.8)
dev.off()

png(file = "Figure 1J_112624.png", width=8, height=6,units = "in",res = 600)
mosaicplot(tab3Perc,las=1,col=genomic_features_col$Region[colnames(tab3Perc)], main="Genomic Distribution of Significant Methyl Sites")
dev.off()

# Figure 1L

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

our_cols = colorRamp2(c(-1, 0,1), c("steelblue","white","firebrick"))

# Figure 1L
pdf(file = "Figure 1L_112624.pdf", width=8, height=3)
pheatmap(tab4,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"steelblue","white","firebrick"),cluster_cols = F,cluster_rows = F,border_color = NA,legend = F,angle_col = 0)
dev.off()

png(file = "Figure 1L_112624.png", width=8, height=3,units = "in",res = 600)
pheatmap(tab4,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"steelblue","white","firebrick"),cluster_cols = F,cluster_rows = F,border_color = NA,legend = F,angle_col = 0)
dev.off()

####
# Investigate up and down-reg in methylation data

# gastro
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

alltissuemethsigW8upordown <- rbind(c(50,83,33),
                                    c(21,70,32),
                                    c(6,43,19),
                                    c(16,53,32),
                                    c(37,70,26),
                                    c(28,85,31),
                                    c(76,272,526),
                                    c(57,242,244))
rownames(alltissuemethsigW8upordown) <- c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC")
colnames(alltissuemethsigW8upordown) <- c("W8 Down-Regulated","Neither","W8 Up-Regulated")
tissuemeta <- data.frame(row.names = rownames(alltissuemethsigW8upordown),
                         "Tissue" = rownames(alltissuemethsigW8upordown))

# Figure Aside for Reviewer on Dividing DAR sets by up and down regulation
png(file = "Reviewer Response Fig 1C_112624.png", width=7.5, height=3,units = "in",res = 600)
pheatmap(t(alltissuemethsigW8upordown),cluster_rows = F,cluster_cols = F,color = colorpanel(100,"white","firebrick"),display_numbers = T,number_color = "black",number_format = "%.0f",fontsize = 15,annotation_col = tissuemeta,annotation_colors = ann_cols,show_colnames = F,border_color = NA)
dev.off()

methannotationfrequencyacrosstissuesupordown <- matrix(0L,ncol = 10,nrow = 24)
colnames(methannotationfrequencyacrosstissuesupordown) <- names(table(trimmedmethanno$custom_annotation))
rownames(methannotationfrequencyacrosstissuesupordown) <- c("SKM-GN Total",
                                                        "SKM-GN UpReg",
                                                        "SKM-GN DownReg",
                                                        "HEART Total",
                                                        "HEART UpReg",
                                                        "HEART DownReg",
                                                        "HIPPOC Total",
                                                        "HIPPOC UpReg",
                                                        "HIPPOC DownReg",
                                                        "KIDNEY Total",
                                                        "KIDNEY UpReg",
                                                        "KIDNEY DownReg",
                                                        "LIVER Total",
                                                        "LIVER UpReg",
                                                        "LIVER DownReg",
                                                        "LUNG Total",
                                                        "LUNG UpReg",
                                                        "LUNG DownReg",
                                                        "BAT Total",
                                                        "BAT UpReg",
                                                        "BAT DownReg",
                                                        "WAT-SC Total",
                                                        "WAT-SC UpReg",
                                                        "WAT-SC DownReg")

for(i in 1:dim(methannotationfrequencyacrosstissuesupordown)[2]){
  ouranno <- colnames(methannotationfrequencyacrosstissuesupordown)[i]
  methannotationfrequencyacrosstissuesupordown["SKM-GN Total",i] <- sum(trimmedmethanno[gastromethsig,"custom_annotation"] %in% ouranno)/length(gastromethsig)
  methannotationfrequencyacrosstissuesupordown["HEART Total",i] <- sum(trimmedmethanno[heartmethsig,"custom_annotation"] %in% ouranno)/length(heartmethsig)
  methannotationfrequencyacrosstissuesupordown["HIPPOC Total",i] <- sum(trimmedmethanno[hippomethsig,"custom_annotation"] %in% ouranno)/length(hippomethsig)
  methannotationfrequencyacrosstissuesupordown["KIDNEY Total",i] <- sum(trimmedmethanno[kidneymethsig,"custom_annotation"] %in% ouranno)/length(kidneymethsig)
  methannotationfrequencyacrosstissuesupordown["LIVER Total",i] <- sum(trimmedmethanno[livermethsig,"custom_annotation"] %in% ouranno)/length(livermethsig)
  methannotationfrequencyacrosstissuesupordown["LUNG Total",i] <- sum(trimmedmethanno[lungmethsig,"custom_annotation"] %in% ouranno)/length(lungmethsig)
  methannotationfrequencyacrosstissuesupordown["BAT Total",i] <- sum(trimmedmethanno[brownmethsig,"custom_annotation"] %in% ouranno)/length(brownmethsig)
  methannotationfrequencyacrosstissuesupordown["WAT-SC Total",i] <- sum(trimmedmethanno[whitemethsig,"custom_annotation"] %in% ouranno)/length(whitemethsig)
  
  methannotationfrequencyacrosstissuesupordown["SKM-GN UpReg",i] <- sum(trimmedmethanno[rownames(gastromethsigupordowndf)[gastromethsigupordowndf$UpReg == "1"],"custom_annotation"] %in% ouranno)/length(rownames(gastromethsigupordowndf)[gastromethsigupordowndf$UpReg == "1"])
  methannotationfrequencyacrosstissuesupordown["SKM-GN DownReg",i] <- sum(trimmedmethanno[rownames(gastromethsigupordowndf)[gastromethsigupordowndf$DownReg == "1"],"custom_annotation"] %in% ouranno)/length(rownames(gastromethsigupordowndf)[gastromethsigupordowndf$DownReg == "1"])
  
  methannotationfrequencyacrosstissuesupordown["HEART UpReg",i] <- sum(trimmedmethanno[rownames(heartmethsigupordowndf)[heartmethsigupordowndf$UpReg == "1"],"custom_annotation"] %in% ouranno)/length(rownames(heartmethsigupordowndf)[heartmethsigupordowndf$UpReg == "1"])
  methannotationfrequencyacrosstissuesupordown["HEART DownReg",i] <- sum(trimmedmethanno[rownames(heartmethsigupordowndf)[heartmethsigupordowndf$DownReg == "1"],"custom_annotation"] %in% ouranno)/length(rownames(heartmethsigupordowndf)[heartmethsigupordowndf$DownReg == "1"])
  
  methannotationfrequencyacrosstissuesupordown["HIPPOC UpReg",i] <- sum(trimmedmethanno[rownames(hippomethsigupordowndf)[hippomethsigupordowndf$UpReg == "1"],"custom_annotation"] %in% ouranno)/length(rownames(hippomethsigupordowndf)[hippomethsigupordowndf$UpReg == "1"])
  methannotationfrequencyacrosstissuesupordown["HIPPOC DownReg",i] <- sum(trimmedmethanno[rownames(hippomethsigupordowndf)[hippomethsigupordowndf$DownReg == "1"],"custom_annotation"] %in% ouranno)/length(rownames(hippomethsigupordowndf)[hippomethsigupordowndf$DownReg == "1"])
  
  methannotationfrequencyacrosstissuesupordown["KIDNEY UpReg",i] <- sum(trimmedmethanno[rownames(kidneymethsigupordowndf)[kidneymethsigupordowndf$UpReg == "1"],"custom_annotation"] %in% ouranno)/length(rownames(kidneymethsigupordowndf)[kidneymethsigupordowndf$UpReg == "1"])
  methannotationfrequencyacrosstissuesupordown["KIDNEY DownReg",i] <- sum(trimmedmethanno[rownames(kidneymethsigupordowndf)[kidneymethsigupordowndf$DownReg == "1"],"custom_annotation"] %in% ouranno)/length(rownames(kidneymethsigupordowndf)[kidneymethsigupordowndf$DownReg == "1"])
  
  methannotationfrequencyacrosstissuesupordown["LIVER UpReg",i] <- sum(trimmedmethanno[rownames(livermethsigupordowndf)[livermethsigupordowndf$UpReg == "1"],"custom_annotation"] %in% ouranno)/length(rownames(livermethsigupordowndf)[livermethsigupordowndf$UpReg == "1"])
  methannotationfrequencyacrosstissuesupordown["LIVER DownReg",i] <- sum(trimmedmethanno[rownames(livermethsigupordowndf)[livermethsigupordowndf$DownReg == "1"],"custom_annotation"] %in% ouranno)/length(rownames(livermethsigupordowndf)[livermethsigupordowndf$DownReg == "1"])
  
  methannotationfrequencyacrosstissuesupordown["LUNG UpReg",i] <- sum(trimmedmethanno[rownames(lungmethsigupordowndf)[lungmethsigupordowndf$UpReg == "1"],"custom_annotation"] %in% ouranno)/length(rownames(lungmethsigupordowndf)[lungmethsigupordowndf$UpReg == "1"])
  methannotationfrequencyacrosstissuesupordown["LUNG DownReg",i] <- sum(trimmedmethanno[rownames(lungmethsigupordowndf)[lungmethsigupordowndf$DownReg == "1"],"custom_annotation"] %in% ouranno)/length(rownames(lungmethsigupordowndf)[lungmethsigupordowndf$DownReg == "1"])
  
  methannotationfrequencyacrosstissuesupordown["BAT UpReg",i] <- sum(trimmedmethanno[rownames(brownmethsigupordowndf)[brownmethsigupordowndf$UpReg == "1"],"custom_annotation"] %in% ouranno)/length(rownames(brownmethsigupordowndf)[brownmethsigupordowndf$UpReg == "1"])
  methannotationfrequencyacrosstissuesupordown["BAT DownReg",i] <- sum(trimmedmethanno[rownames(brownmethsigupordowndf)[brownmethsigupordowndf$DownReg == "1"],"custom_annotation"] %in% ouranno)/length(rownames(brownmethsigupordowndf)[brownmethsigupordowndf$DownReg == "1"])
  
  methannotationfrequencyacrosstissuesupordown["WAT-SC UpReg",i] <- sum(trimmedmethanno[rownames(whitemethsigupordowndf)[whitemethsigupordowndf$UpReg == "1"],"custom_annotation"] %in% ouranno)/length(rownames(whitemethsigupordowndf)[whitemethsigupordowndf$UpReg == "1"])
  methannotationfrequencyacrosstissuesupordown["WAT-SC DownReg",i] <- sum(trimmedmethanno[rownames(whitemethsigupordowndf)[whitemethsigupordowndf$DownReg == "1"],"custom_annotation"] %in% ouranno)/length(rownames(whitemethsigupordowndf)[whitemethsigupordowndf$DownReg == "1"])
  
}

png(file = "Reviewer Response Fig 1D_112624.png",width = 14,height = 6,units = "in",res = 600)
pheatmap(methannotationfrequencyacrosstissuesupordown,cluster_rows = F,cluster_cols = F,color = colorpanel(101,"white","firebrick"),display_numbers = T,number_color = "black",angle_col = 0,border_color = NA)
dev.off()


rm(epigen_rrbs)
gc()

save.image("Figure1G_1H_1I_1J_1K_1L_S12_ReviewResponse1_112624.RData")