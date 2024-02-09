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
peakanno <- readRDS("peakanno.RDS")
trimmedmethanno <- readRDS("trimmedmethanno.RDS")

#####
# Figure 3E
####

peakanno$mid <- round((peakanno$end + peakanno$start)/2)

# gastro

gastromethsiganno <- trimmedmethanno[trimmedmethanno$feature_ID %in% gastromethsig,]
fingastromethsiganno <- data.frame(row.names = gastromethsig,
                                   "assay" = rep("epigen-rrbs",length(gastromethsig)),
                                   "feature_ID" = gastromethsig,
                                   "chrom" = rep(0,length(gastromethsig)),
                                   "start" = rep(0,length(gastromethsig)),
                                   "end" = rep(0,length(gastromethsig)),
                                   "width" = rep(0,length(gastromethsig)),
                                   "chipseeker_annotation" = rep("",length(gastromethsig)),
                                   "custom_annotation" = rep("",length(gastromethsig)),
                                   "distanceToTSS" = rep(0,length(gastromethsig)),
                                   "relationship_to_gene" = rep(0,length(gastromethsig)),
                                   "ensembl_gene" = rep("",length(gastromethsig)),
                                   "geneStart" = rep(0,length(gastromethsig)),
                                   "geneEnd" = rep(0,length(gastromethsig)),
                                   "geneLength" = rep(0,length(gastromethsig)),
                                   "geneStrand" = rep(0,length(gastromethsig)))
for(i in 1:dim(fingastromethsiganno)[1]){
  if(i%%20 == 0){
    print(i)
  }
  ourmeth <- rownames(fingastromethsiganno)[i]
  ourannolist <- gastromethsiganno[gastromethsiganno$feature_ID %in% ourmeth,]
  fingastromethsiganno[i,] <- ourannolist[1,]
}
fingastromethsiganno$sitemid <- (fingastromethsiganno$start + fingastromethsiganno$end)/2

# heart

heartmethsiganno <- trimmedmethanno[trimmedmethanno$feature_ID %in% heartmethsig,]
finheartmethsiganno <- data.frame(row.names = heartmethsig,
                                  "assay" = rep("epigen-rrbs",length(heartmethsig)),
                                  "feature_ID" = heartmethsig,
                                  "chrom" = rep(0,length(heartmethsig)),
                                  "start" = rep(0,length(heartmethsig)),
                                  "end" = rep(0,length(heartmethsig)),
                                  "width" = rep(0,length(heartmethsig)),
                                  "chipseeker_annotation" = rep("",length(heartmethsig)),
                                  "custom_annotation" = rep("",length(heartmethsig)),
                                  "distanceToTSS" = rep(0,length(heartmethsig)),
                                  "relationship_to_gene" = rep(0,length(heartmethsig)),
                                  "ensembl_gene" = rep("",length(heartmethsig)),
                                  "geneStart" = rep(0,length(heartmethsig)),
                                  "geneEnd" = rep(0,length(heartmethsig)),
                                  "geneLength" = rep(0,length(heartmethsig)),
                                  "geneStrand" = rep(0,length(heartmethsig)))
for(i in 1:dim(finheartmethsiganno)[1]){
  if(i%%20 == 0){
    print(i)
  }
  ourmeth <- rownames(finheartmethsiganno)[i]
  ourannolist <- heartmethsiganno[heartmethsiganno$feature_ID %in% ourmeth,]
  finheartmethsiganno[i,] <- ourannolist[1,]
}
finheartmethsiganno$sitemid <- (finheartmethsiganno$start + finheartmethsiganno$end)/2

# hippo

hippomethsiganno <- trimmedmethanno[trimmedmethanno$feature_ID %in% hippomethsig,]
finhippomethsiganno <- data.frame(row.names = hippomethsig,
                                  "assay" = rep("epigen-rrbs",length(hippomethsig)),
                                  "feature_ID" = hippomethsig,
                                  "chrom" = rep(0,length(hippomethsig)),
                                  "start" = rep(0,length(hippomethsig)),
                                  "end" = rep(0,length(hippomethsig)),
                                  "width" = rep(0,length(hippomethsig)),
                                  "chipseeker_annotation" = rep("",length(hippomethsig)),
                                  "custom_annotation" = rep("",length(hippomethsig)),
                                  "distanceToTSS" = rep(0,length(hippomethsig)),
                                  "relationship_to_gene" = rep(0,length(hippomethsig)),
                                  "ensembl_gene" = rep("",length(hippomethsig)),
                                  "geneStart" = rep(0,length(hippomethsig)),
                                  "geneEnd" = rep(0,length(hippomethsig)),
                                  "geneLength" = rep(0,length(hippomethsig)),
                                  "geneStrand" = rep(0,length(hippomethsig)))
for(i in 1:dim(finhippomethsiganno)[1]){
  if(i%%20 == 0){
    print(i)
  }
  ourmeth <- rownames(finhippomethsiganno)[i]
  ourannolist <- hippomethsiganno[hippomethsiganno$feature_ID %in% ourmeth,]
  finhippomethsiganno[i,] <- ourannolist[1,]
}
finhippomethsiganno$sitemid <- (finhippomethsiganno$start + finhippomethsiganno$end)/2

# kidney

kidneymethsiganno <- trimmedmethanno[trimmedmethanno$feature_ID %in% kidneymethsig,]
finkidneymethsiganno <- data.frame(row.names = kidneymethsig,
                                   "assay" = rep("epigen-rrbs",length(kidneymethsig)),
                                   "feature_ID" = kidneymethsig,
                                   "chrom" = rep(0,length(kidneymethsig)),
                                   "start" = rep(0,length(kidneymethsig)),
                                   "end" = rep(0,length(kidneymethsig)),
                                   "width" = rep(0,length(kidneymethsig)),
                                   "chipseeker_annotation" = rep("",length(kidneymethsig)),
                                   "custom_annotation" = rep("",length(kidneymethsig)),
                                   "distanceToTSS" = rep(0,length(kidneymethsig)),
                                   "relationship_to_gene" = rep(0,length(kidneymethsig)),
                                   "ensembl_gene" = rep("",length(kidneymethsig)),
                                   "geneStart" = rep(0,length(kidneymethsig)),
                                   "geneEnd" = rep(0,length(kidneymethsig)),
                                   "geneLength" = rep(0,length(kidneymethsig)),
                                   "geneStrand" = rep(0,length(kidneymethsig)))
for(i in 1:dim(finkidneymethsiganno)[1]){
  if(i%%20 == 0){
    print(i)
  }
  ourmeth <- rownames(finkidneymethsiganno)[i]
  ourannolist <- kidneymethsiganno[kidneymethsiganno$feature_ID %in% ourmeth,]
  finkidneymethsiganno[i,] <- ourannolist[1,]
}
finkidneymethsiganno$sitemid <- (finkidneymethsiganno$start + finkidneymethsiganno$end)/2


# liver

livermethsiganno <- trimmedmethanno[trimmedmethanno$feature_ID %in% livermethsig,]
finlivermethsiganno <- data.frame(row.names = livermethsig,
                                  "assay" = rep("epigen-rrbs",length(livermethsig)),
                                  "feature_ID" = livermethsig,
                                  "chrom" = rep(0,length(livermethsig)),
                                  "start" = rep(0,length(livermethsig)),
                                  "end" = rep(0,length(livermethsig)),
                                  "width" = rep(0,length(livermethsig)),
                                  "chipseeker_annotation" = rep("",length(livermethsig)),
                                  "custom_annotation" = rep("",length(livermethsig)),
                                  "distanceToTSS" = rep(0,length(livermethsig)),
                                  "relationship_to_gene" = rep(0,length(livermethsig)),
                                  "ensembl_gene" = rep("",length(livermethsig)),
                                  "geneStart" = rep(0,length(livermethsig)),
                                  "geneEnd" = rep(0,length(livermethsig)),
                                  "geneLength" = rep(0,length(livermethsig)),
                                  "geneStrand" = rep(0,length(livermethsig)))
for(i in 1:dim(finlivermethsiganno)[1]){
  if(i%%20 == 0){
    print(i)
  }
  ourmeth <- rownames(finlivermethsiganno)[i]
  ourannolist <- livermethsiganno[livermethsiganno$feature_ID %in% ourmeth,]
  finlivermethsiganno[i,] <- ourannolist[1,]
}
finlivermethsiganno$sitemid <- (finlivermethsiganno$start + finlivermethsiganno$end)/2

# lung

lungmethsiganno <- trimmedmethanno[trimmedmethanno$feature_ID %in% lungmethsig,]
finlungmethsiganno <- data.frame(row.names = lungmethsig,
                                 "assay" = rep("epigen-rrbs",length(lungmethsig)),
                                 "feature_ID" = lungmethsig,
                                 "chrom" = rep(0,length(lungmethsig)),
                                 "start" = rep(0,length(lungmethsig)),
                                 "end" = rep(0,length(lungmethsig)),
                                 "width" = rep(0,length(lungmethsig)),
                                 "chipseeker_annotation" = rep("",length(lungmethsig)),
                                 "custom_annotation" = rep("",length(lungmethsig)),
                                 "distanceToTSS" = rep(0,length(lungmethsig)),
                                 "relationship_to_gene" = rep(0,length(lungmethsig)),
                                 "ensembl_gene" = rep("",length(lungmethsig)),
                                 "geneStart" = rep(0,length(lungmethsig)),
                                 "geneEnd" = rep(0,length(lungmethsig)),
                                 "geneLength" = rep(0,length(lungmethsig)),
                                 "geneStrand" = rep(0,length(lungmethsig)))
for(i in 1:dim(finlungmethsiganno)[1]){
  if(i%%20 == 0){
    print(i)
  }
  ourmeth <- rownames(finlungmethsiganno)[i]
  ourannolist <- lungmethsiganno[lungmethsiganno$feature_ID %in% ourmeth,]
  finlungmethsiganno[i,] <- ourannolist[1,]
}
finlungmethsiganno$sitemid <- (finlungmethsiganno$start + finlungmethsiganno$end)/2


# brown

brownmethsiganno <- trimmedmethanno[trimmedmethanno$feature_ID %in% brownmethsig,]
finbrownmethsiganno <- data.frame(row.names = brownmethsig,
                                  "assay" = rep("epigen-rrbs",length(brownmethsig)),
                                  "feature_ID" = brownmethsig,
                                  "chrom" = rep(0,length(brownmethsig)),
                                  "start" = rep(0,length(brownmethsig)),
                                  "end" = rep(0,length(brownmethsig)),
                                  "width" = rep(0,length(brownmethsig)),
                                  "chipseeker_annotation" = rep("",length(brownmethsig)),
                                  "custom_annotation" = rep("",length(brownmethsig)),
                                  "distanceToTSS" = rep(0,length(brownmethsig)),
                                  "relationship_to_gene" = rep(0,length(brownmethsig)),
                                  "ensembl_gene" = rep("",length(brownmethsig)),
                                  "geneStart" = rep(0,length(brownmethsig)),
                                  "geneEnd" = rep(0,length(brownmethsig)),
                                  "geneLength" = rep(0,length(brownmethsig)),
                                  "geneStrand" = rep(0,length(brownmethsig)))
for(i in 1:dim(finbrownmethsiganno)[1]){
  if(i%%20 == 0){
    print(i)
  }
  ourmeth <- rownames(finbrownmethsiganno)[i]
  ourannolist <- brownmethsiganno[brownmethsiganno$feature_ID %in% ourmeth,]
  finbrownmethsiganno[i,] <- ourannolist[1,]
}
finbrownmethsiganno$sitemid <- (finbrownmethsiganno$start + finbrownmethsiganno$end)/2


# white

whitemethsiganno <- trimmedmethanno[trimmedmethanno$feature_ID %in% whitemethsig,]
finwhitemethsiganno <- data.frame(row.names = whitemethsig,
                                  "assay" = rep("epigen-rrbs",length(whitemethsig)),
                                  "feature_ID" = whitemethsig,
                                  "chrom" = rep(0,length(whitemethsig)),
                                  "start" = rep(0,length(whitemethsig)),
                                  "end" = rep(0,length(whitemethsig)),
                                  "width" = rep(0,length(whitemethsig)),
                                  "chipseeker_annotation" = rep("",length(whitemethsig)),
                                  "custom_annotation" = rep("",length(whitemethsig)),
                                  "distanceToTSS" = rep(0,length(whitemethsig)),
                                  "relationship_to_gene" = rep(0,length(whitemethsig)),
                                  "ensembl_gene" = rep("",length(whitemethsig)),
                                  "geneStart" = rep(0,length(whitemethsig)),
                                  "geneEnd" = rep(0,length(whitemethsig)),
                                  "geneLength" = rep(0,length(whitemethsig)),
                                  "geneStrand" = rep(0,length(whitemethsig)))
for(i in 1:dim(finwhitemethsiganno)[1]){
  if(i%%20 == 0){
    print(i)
  }
  ourmeth <- rownames(finwhitemethsiganno)[i]
  ourannolist <- whitemethsiganno[whitemethsiganno$feature_ID %in% ourmeth,]
  finwhitemethsiganno[i,] <- ourannolist[1,]
}
finwhitemethsiganno$sitemid <- (finwhitemethsiganno$start + finwhitemethsiganno$end)/2



#SKM-GN

gastromethatacdistance <- matrix(0L,nrow = length(gastromethsig),ncol = length(gastroatacsig))
rownames(gastromethatacdistance) <- gastromethsig
colnames(gastromethatacdistance) <- gastroatacsig

gastroatacsigpeakanno <- peakanno[gastroatacsig,]

for(i in 1:length(gastromethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(gastroatacsig)){
    gastromethatacdistance[i,j] <- (gastroatacsigpeakanno[gastroatacsig[j],"chrom"] == fingastromethsiganno$chrom[i])*(abs(gastroatacsigpeakanno[gastroatacsig[j],"mid"] - fingastromethsiganno$sitemid[i]))
  }
}

gastromethatacdistancedf <- data.frame(row.names = gastromethsig,
                                       "Region" = fingastromethsiganno[gastromethsig,"custom_annotation"],
                                       "Distance" = rep(0,length(gastromethsig)))
for(i in 1:length(gastromethsig)){
  
  ourrow <- gastromethatacdistance[i,]
  gastromethatacdistancedf[i,"Distance"] <- min(ourrow[ourrow > 0])
  
}


# HEART

heartmethatacdistance <- matrix(0L,nrow = length(heartmethsig),ncol = length(heartatacsig))
rownames(heartmethatacdistance) <- heartmethsig
colnames(heartmethatacdistance) <- heartatacsig

heartatacsigpeakanno <- peakanno[heartatacsig,]

for(i in 1:length(heartmethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(heartatacsig)){
    heartmethatacdistance[i,j] <- (heartatacsigpeakanno[heartatacsig[j],"chrom"] == finheartmethsiganno$chrom[i])*(abs(heartatacsigpeakanno[heartatacsig[j],"mid"] - finheartmethsiganno$sitemid[i]))
  }
}

heartmethatacdistancedf <- data.frame(row.names = heartmethsig,
                                      "Region" = finheartmethsiganno[heartmethsig,"custom_annotation"],
                                      "Distance" = rep(0,length(heartmethsig)))
for(i in 1:length(heartmethsig)){
  
  ourrow <- heartmethatacdistance[i,]
  heartmethatacdistancedf[i,"Distance"] <- min(ourrow[ourrow > 0])
  
}

# HIPPOC

hippomethatacdistance <- matrix(0L,nrow = length(hippomethsig),ncol = length(hippoatacsig))
rownames(hippomethatacdistance) <- hippomethsig
colnames(hippomethatacdistance) <- hippoatacsig

hippoatacsigpeakanno <- peakanno[hippoatacsig,]

for(i in 1:length(hippomethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(hippoatacsig)){
    hippomethatacdistance[i,j] <- (hippoatacsigpeakanno[hippoatacsig[j],"chrom"] == finhippomethsiganno$chrom[i])*(abs(hippoatacsigpeakanno[hippoatacsig[j],"mid"] - finhippomethsiganno$sitemid[i]))
  }
}

hippomethatacdistancedf <- data.frame(row.names = hippomethsig,
                                      "Region" = finhippomethsiganno[hippomethsig,"custom_annotation"],
                                      "Distance" = rep(0,length(hippomethsig)))
for(i in 1:length(hippomethsig)){
  
  ourrow <- hippomethatacdistance[i,]
  hippomethatacdistancedf[i,"Distance"] <- min(ourrow[ourrow > 0])
  
}

# KIDNEY

kidneymethatacdistance <- matrix(0L,nrow = length(kidneymethsig),ncol = length(kidneyatacsig))
rownames(kidneymethatacdistance) <- kidneymethsig
colnames(kidneymethatacdistance) <- kidneyatacsig

kidneyatacsigpeakanno <- peakanno[kidneyatacsig,]

for(i in 1:length(kidneymethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(kidneyatacsig)){
    kidneymethatacdistance[i,j] <- (kidneyatacsigpeakanno[kidneyatacsig[j],"chrom"] == finkidneymethsiganno$chrom[i])*(abs(kidneyatacsigpeakanno[kidneyatacsig[j],"mid"] - finkidneymethsiganno$sitemid[i]))
  }
}

kidneymethatacdistancedf <- data.frame(row.names = kidneymethsig,
                                       "Region" = finkidneymethsiganno[kidneymethsig,"custom_annotation"],
                                       "Distance" = rep(0,length(kidneymethsig)))
for(i in 1:length(kidneymethsig)){
  
  ourrow <- kidneymethatacdistance[i,]
  kidneymethatacdistancedf[i,"Distance"] <- min(ourrow[ourrow > 0])
  
}

# LIVER

livermethatacdistance <- matrix(0L,nrow = length(livermethsig),ncol = length(liveratacsig))
rownames(livermethatacdistance) <- livermethsig
colnames(livermethatacdistance) <- liveratacsig

liveratacsigpeakanno <- peakanno[liveratacsig,]

for(i in 1:length(livermethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(liveratacsig)){
    livermethatacdistance[i,j] <- (liveratacsigpeakanno[liveratacsig[j],"chrom"] == finlivermethsiganno$chrom[i])*(abs(liveratacsigpeakanno[liveratacsig[j],"mid"] - finlivermethsiganno$sitemid[i]))
  }
}

livermethatacdistancedf <- data.frame(row.names = livermethsig,
                                      "Region" = finlivermethsiganno[livermethsig,"custom_annotation"],
                                      "Distance" = rep(0,length(livermethsig)))
for(i in 1:length(livermethsig)){
  
  ourrow <- livermethatacdistance[i,]
  livermethatacdistancedf[i,"Distance"] <- min(ourrow[ourrow > 0])
  
}

# LUNG

lungmethatacdistance <- matrix(0L,nrow = length(lungmethsig),ncol = length(lungatacsig))
rownames(lungmethatacdistance) <- lungmethsig
colnames(lungmethatacdistance) <- lungatacsig

lungatacsigpeakanno <- peakanno[lungatacsig,]

for(i in 1:length(lungmethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(lungatacsig)){
    lungmethatacdistance[i,j] <- (lungatacsigpeakanno[lungatacsig[j],"chrom"] == finlungmethsiganno$chrom[i])*(abs(lungatacsigpeakanno[lungatacsig[j],"mid"] - finlungmethsiganno$sitemid[i]))
  }
}

lungmethatacdistancedf <- data.frame(row.names = lungmethsig,
                                     "Region" = finlungmethsiganno[lungmethsig,"custom_annotation"],
                                     "Distance" = rep(0,length(lungmethsig)))
for(i in 1:length(lungmethsig)){
  
  ourrow <- lungmethatacdistance[i,]
  lungmethatacdistancedf[i,"Distance"] <- min(ourrow[ourrow > 0])
  
}

# BAT

brownmethatacdistance <- matrix(0L,nrow = length(brownmethsig),ncol = length(brownatacsig))
rownames(brownmethatacdistance) <- brownmethsig
colnames(brownmethatacdistance) <- brownatacsig

brownatacsigpeakanno <- peakanno[brownatacsig,]

for(i in 1:length(brownmethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(brownatacsig)){
    brownmethatacdistance[i,j] <- (brownatacsigpeakanno[brownatacsig[j],"chrom"] == finbrownmethsiganno$chrom[i])*(abs(brownatacsigpeakanno[brownatacsig[j],"mid"] - finbrownmethsiganno$sitemid[i]))
  }
}

brownmethatacdistancedf <- data.frame(row.names = brownmethsig,
                                      "Region" = finbrownmethsiganno[brownmethsig,"custom_annotation"],
                                      "Distance" = rep(0,length(brownmethsig)))
for(i in 1:length(brownmethsig)){
  
  ourrow <- brownmethatacdistance[i,]
  brownmethatacdistancedf[i,"Distance"] <- min(ourrow[ourrow > 0])
  
}

# WAT-SC

whitemethatacdistance <- matrix(0L,nrow = length(whitemethsig),ncol = length(whiteatacsig))
rownames(whitemethatacdistance) <- whitemethsig
colnames(whitemethatacdistance) <- whiteatacsig

whiteatacsigpeakanno <- peakanno[whiteatacsig,]

for(i in 1:length(whitemethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(whiteatacsig)){
    whitemethatacdistance[i,j] <- (whiteatacsigpeakanno[whiteatacsig[j],"chrom"] == finwhitemethsiganno$chrom[i])*(abs(whiteatacsigpeakanno[whiteatacsig[j],"mid"] - finwhitemethsiganno$sitemid[i]))
  }
}

whitemethatacdistancedf <- data.frame(row.names = whitemethsig,
                                      "Region" = finwhitemethsiganno[whitemethsig,"custom_annotation"],
                                      "Distance" = rep(0,length(whitemethsig)))
for(i in 1:length(whitemethsig)){
  
  ourrow <- whitemethatacdistance[i,]
  whitemethatacdistancedf[i,"Distance"] <- min(ourrow[ourrow > 0])
  
}

totalmethatacdistancedf <- data.frame("Region" = c(gastromethatacdistancedf$Region,
                                                   heartmethatacdistancedf$Region,
                                                   hippomethatacdistancedf$Region,
                                                   kidneymethatacdistancedf$Region,
                                                   livermethatacdistancedf$Region,
                                                   lungmethatacdistancedf$Region,
                                                   brownmethatacdistancedf$Region,
                                                   whitemethatacdistancedf$Region),
                                      "Distance" = c(gastromethatacdistancedf$Distance,
                                                     heartmethatacdistancedf$Distance,
                                                     hippomethatacdistancedf$Distance,
                                                     kidneymethatacdistancedf$Distance,
                                                     livermethatacdistancedf$Distance,
                                                     lungmethatacdistancedf$Distance,
                                                     brownmethatacdistancedf$Distance,
                                                     whitemethatacdistancedf$Distance),
                                      "Tissue" = c(rep("SKM-GN",dim(gastromethatacdistancedf)[1]),
                                                   rep("HEART",dim(heartmethatacdistancedf)[1]),
                                                   rep("HIPPOC",dim(hippomethatacdistancedf)[1]),
                                                   rep("KIDNEY",dim(kidneymethatacdistancedf)[1]),
                                                   rep("LIVER",dim(livermethatacdistancedf)[1]),
                                                   rep("LUNG",dim(lungmethatacdistancedf)[1]),
                                                   rep("BAT",dim(brownmethatacdistancedf)[1]),
                                                   rep("WAT-SC",dim(whitemethatacdistancedf)[1])))
totalmethatacdistancedf <- totalmethatacdistancedf[!totalmethatacdistancedf$Region %in% "Overlaps Gene",]


peakanno$mid <- round((peakanno$end + peakanno$start)/2)

#SKM-GN

gastromethatacdistance <- matrix(0L,nrow = length(gastromethsig),ncol = length(gastroatacsig))
rownames(gastromethatacdistance) <- gastromethsig
colnames(gastromethatacdistance) <- gastroatacsig

gastroatacsigpeakanno <- peakanno[gastroatacsig,]

for(i in 1:length(gastromethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(gastroatacsig)){
    gastromethatacdistance[i,j] <- (gastroatacsigpeakanno[gastroatacsig[j],"chrom"] == fingastromethsiganno$chrom[i])*(abs(gastroatacsigpeakanno[gastroatacsig[j],"mid"] - fingastromethsiganno$sitemid[i]))
  }
}

gastromethatacdistancedf <- data.frame(row.names = gastromethsig,
                                       "Region" = fingastromethsiganno[gastromethsig,"custom_annotation"],
                                       "Distance" = rep(0,length(gastromethsig)))
for(i in 1:length(gastromethsig)){
  
  ourrow <- gastromethatacdistance[i,]
  gastromethatacdistancedf[i,"Distance"] <- min(ourrow[ourrow > 0])
  
}


# HEART

heartmethatacdistance <- matrix(0L,nrow = length(heartmethsig),ncol = length(heartatacsig))
rownames(heartmethatacdistance) <- heartmethsig
colnames(heartmethatacdistance) <- heartatacsig

heartatacsigpeakanno <- peakanno[heartatacsig,]

for(i in 1:length(heartmethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(heartatacsig)){
    heartmethatacdistance[i,j] <- (heartatacsigpeakanno[heartatacsig[j],"chrom"] == finheartmethsiganno$chrom[i])*(abs(heartatacsigpeakanno[heartatacsig[j],"mid"] - finheartmethsiganno$sitemid[i]))
  }
}

heartmethatacdistancedf <- data.frame(row.names = heartmethsig,
                                      "Region" = finheartmethsiganno[heartmethsig,"custom_annotation"],
                                      "Distance" = rep(0,length(heartmethsig)))
for(i in 1:length(heartmethsig)){
  
  ourrow <- heartmethatacdistance[i,]
  heartmethatacdistancedf[i,"Distance"] <- min(ourrow[ourrow > 0])
  
}

# HIPPOC

hippomethatacdistance <- matrix(0L,nrow = length(hippomethsig),ncol = length(hippoatacsig))
rownames(hippomethatacdistance) <- hippomethsig
colnames(hippomethatacdistance) <- hippoatacsig

hippoatacsigpeakanno <- peakanno[hippoatacsig,]

for(i in 1:length(hippomethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(hippoatacsig)){
    hippomethatacdistance[i,j] <- (hippoatacsigpeakanno[hippoatacsig[j],"chrom"] == finhippomethsiganno$chrom[i])*(abs(hippoatacsigpeakanno[hippoatacsig[j],"mid"] - finhippomethsiganno$sitemid[i]))
  }
}

hippomethatacdistancedf <- data.frame(row.names = hippomethsig,
                                      "Region" = finhippomethsiganno[hippomethsig,"custom_annotation"],
                                      "Distance" = rep(0,length(hippomethsig)))
for(i in 1:length(hippomethsig)){
  
  ourrow <- hippomethatacdistance[i,]
  hippomethatacdistancedf[i,"Distance"] <- min(ourrow[ourrow > 0])
  
}

# KIDNEY

kidneymethatacdistance <- matrix(0L,nrow = length(kidneymethsig),ncol = length(kidneyatacsig))
rownames(kidneymethatacdistance) <- kidneymethsig
colnames(kidneymethatacdistance) <- kidneyatacsig

kidneyatacsigpeakanno <- peakanno[kidneyatacsig,]

for(i in 1:length(kidneymethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(kidneyatacsig)){
    kidneymethatacdistance[i,j] <- (kidneyatacsigpeakanno[kidneyatacsig[j],"chrom"] == finkidneymethsiganno$chrom[i])*(abs(kidneyatacsigpeakanno[kidneyatacsig[j],"mid"] - finkidneymethsiganno$sitemid[i]))
  }
}

kidneymethatacdistancedf <- data.frame(row.names = kidneymethsig,
                                       "Region" = finkidneymethsiganno[kidneymethsig,"custom_annotation"],
                                       "Distance" = rep(0,length(kidneymethsig)))
for(i in 1:length(kidneymethsig)){
  
  ourrow <- kidneymethatacdistance[i,]
  kidneymethatacdistancedf[i,"Distance"] <- min(ourrow[ourrow > 0])
  
}

# LIVER

livermethatacdistance <- matrix(0L,nrow = length(livermethsig),ncol = length(liveratacsig))
rownames(livermethatacdistance) <- livermethsig
colnames(livermethatacdistance) <- liveratacsig

liveratacsigpeakanno <- peakanno[liveratacsig,]

for(i in 1:length(livermethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(liveratacsig)){
    livermethatacdistance[i,j] <- (liveratacsigpeakanno[liveratacsig[j],"chrom"] == finlivermethsiganno$chrom[i])*(abs(liveratacsigpeakanno[liveratacsig[j],"mid"] - finlivermethsiganno$sitemid[i]))
  }
}

livermethatacdistancedf <- data.frame(row.names = livermethsig,
                                      "Region" = finlivermethsiganno[livermethsig,"custom_annotation"],
                                      "Distance" = rep(0,length(livermethsig)))
for(i in 1:length(livermethsig)){
  
  ourrow <- livermethatacdistance[i,]
  livermethatacdistancedf[i,"Distance"] <- min(ourrow[ourrow > 0])
  
}

# LUNG

lungmethatacdistance <- matrix(0L,nrow = length(lungmethsig),ncol = length(lungatacsig))
rownames(lungmethatacdistance) <- lungmethsig
colnames(lungmethatacdistance) <- lungatacsig

lungatacsigpeakanno <- peakanno[lungatacsig,]

for(i in 1:length(lungmethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(lungatacsig)){
    lungmethatacdistance[i,j] <- (lungatacsigpeakanno[lungatacsig[j],"chrom"] == finlungmethsiganno$chrom[i])*(abs(lungatacsigpeakanno[lungatacsig[j],"mid"] - finlungmethsiganno$sitemid[i]))
  }
}

lungmethatacdistancedf <- data.frame(row.names = lungmethsig,
                                     "Region" = finlungmethsiganno[lungmethsig,"custom_annotation"],
                                     "Distance" = rep(0,length(lungmethsig)))
for(i in 1:length(lungmethsig)){
  
  ourrow <- lungmethatacdistance[i,]
  lungmethatacdistancedf[i,"Distance"] <- min(ourrow[ourrow > 0])
  
}

# BAT

brownmethatacdistance <- matrix(0L,nrow = length(brownmethsig),ncol = length(brownatacsig))
rownames(brownmethatacdistance) <- brownmethsig
colnames(brownmethatacdistance) <- brownatacsig

brownatacsigpeakanno <- peakanno[brownatacsig,]

for(i in 1:length(brownmethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(brownatacsig)){
    brownmethatacdistance[i,j] <- (brownatacsigpeakanno[brownatacsig[j],"chrom"] == finbrownmethsiganno$chrom[i])*(abs(brownatacsigpeakanno[brownatacsig[j],"mid"] - finbrownmethsiganno$sitemid[i]))
  }
}

brownmethatacdistancedf <- data.frame(row.names = brownmethsig,
                                      "Region" = finbrownmethsiganno[brownmethsig,"custom_annotation"],
                                      "Distance" = rep(0,length(brownmethsig)))
for(i in 1:length(brownmethsig)){
  
  ourrow <- brownmethatacdistance[i,]
  brownmethatacdistancedf[i,"Distance"] <- min(ourrow[ourrow > 0])
  
}

# WAT-SC

whitemethatacdistance <- matrix(0L,nrow = length(whitemethsig),ncol = length(whiteatacsig))
rownames(whitemethatacdistance) <- whitemethsig
colnames(whitemethatacdistance) <- whiteatacsig

whiteatacsigpeakanno <- peakanno[whiteatacsig,]

for(i in 1:length(whitemethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(whiteatacsig)){
    whitemethatacdistance[i,j] <- (whiteatacsigpeakanno[whiteatacsig[j],"chrom"] == finwhitemethsiganno$chrom[i])*(abs(whiteatacsigpeakanno[whiteatacsig[j],"mid"] - finwhitemethsiganno$sitemid[i]))
  }
}

whitemethatacdistancedf <- data.frame(row.names = whitemethsig,
                                      "Region" = finwhitemethsiganno[whitemethsig,"custom_annotation"],
                                      "Distance" = rep(0,length(whitemethsig)))
for(i in 1:length(whitemethsig)){
  
  ourrow <- whitemethatacdistance[i,]
  whitemethatacdistancedf[i,"Distance"] <- min(ourrow[ourrow > 0])
  
}

totalmethatacdistancedf <- data.frame("Region" = c(gastromethatacdistancedf$Region,
                                                   heartmethatacdistancedf$Region,
                                                   hippomethatacdistancedf$Region,
                                                   kidneymethatacdistancedf$Region,
                                                   livermethatacdistancedf$Region,
                                                   lungmethatacdistancedf$Region,
                                                   brownmethatacdistancedf$Region,
                                                   whitemethatacdistancedf$Region),
                                      "Distance" = c(gastromethatacdistancedf$Distance,
                                                     heartmethatacdistancedf$Distance,
                                                     hippomethatacdistancedf$Distance,
                                                     kidneymethatacdistancedf$Distance,
                                                     livermethatacdistancedf$Distance,
                                                     lungmethatacdistancedf$Distance,
                                                     brownmethatacdistancedf$Distance,
                                                     whitemethatacdistancedf$Distance),
                                      "Tissue" = c(rep("SKM-GN",dim(gastromethatacdistancedf)[1]),
                                                   rep("HEART",dim(heartmethatacdistancedf)[1]),
                                                   rep("HIPPOC",dim(hippomethatacdistancedf)[1]),
                                                   rep("KIDNEY",dim(kidneymethatacdistancedf)[1]),
                                                   rep("LIVER",dim(livermethatacdistancedf)[1]),
                                                   rep("LUNG",dim(lungmethatacdistancedf)[1]),
                                                   rep("BAT",dim(brownmethatacdistancedf)[1]),
                                                   rep("WAT-SC",dim(whitemethatacdistancedf)[1])))
totalmethatacdistancedf <- totalmethatacdistancedf[!totalmethatacdistancedf$Region %in% "Overlaps Gene",]

pdf(file = "Figure 3E.pdf",width = 7,height = 6)
ggplot(totalmethatacdistancedf, aes(x = Distance)) +
  geom_histogram(aes(color = Tissue,fill = Tissue),
                 position = "stack", bins = 30) + theme_classic() +
  scale_color_manual(values = c("#8c5220","#f28b2f","#bf7534","#7553a7","#da6c75","#04bf8a","#088c03","#214da6")) + 
  scale_fill_manual(values = c("#8c5220","#f28b2f","#bf7534","#7553a7","#da6c75","#04bf8a","#088c03","#214da6"))  + scale_x_log10(breaks=c(1,100,10000,1000000,100000000)) + xlab("Distance to Nearest DAR") + ylab("Count") + theme(axis.title = element_text(size = 15),axis.text = element_text(size = 12),legend.title = element_text(size = 15),legend.text = element_text(size = 12))
dev.off()

png(file = "Figure 3E.png",width = 7,height = 6,units = "in",res = 600)
ggplot(totalmethatacdistancedf, aes(x = Distance)) +
  geom_histogram(aes(color = Tissue,fill = Tissue),
                 position = "stack", bins = 30) + theme_classic() +
  scale_color_manual(values = c("#8c5220","#f28b2f","#bf7534","#7553a7","#da6c75","#04bf8a","#088c03","#214da6")) + 
  scale_fill_manual(values = c("#8c5220","#f28b2f","#bf7534","#7553a7","#da6c75","#04bf8a","#088c03","#214da6"))  + scale_x_log10(breaks=c(1,100,10000,1000000,100000000)) + xlab("Distance to Nearest DAR") + ylab("Count") + theme(axis.title = element_text(size = 15),axis.text = element_text(size = 12),legend.title = element_text(size = 15),legend.text = element_text(size = 12))
dev.off()

#####
# Figure 3H
####

load("atacsigl2fcmat.RData")
load("methsigl2fcmat.RData")

# SKM-GN

gastromethataccor <- matrix(0L,nrow = dim(gastromethatacdistance)[1],ncol = dim(gastromethatacdistance)[2])
rownames(gastromethataccor) <- rownames(gastromethatacdistance)
colnames(gastromethataccor) <- colnames(gastromethatacdistance)
for(i in 1:dim(gastromethataccor)[1]){
  if(i%%20 == 0){
    print(i)
  }
  for(j in 1:dim(gastromethataccor)[2]){
    ourmeth <- rownames(gastromethatacdistance)[i]
    ourpeak <- colnames(gastromethatacdistance)[j]
    gastromethataccor[i,j] <- cor(gastromethl2fc[ourmeth,],gastrosigatacl2fc[ourpeak,])
  }
}
gastromethataccor[is.na(gastromethataccor)] <- 0

gastromethatacnoabsdistance <- matrix(0L,nrow = length(gastromethsig),ncol = length(gastroatacsig))
rownames(gastromethatacnoabsdistance) <- gastromethsig
colnames(gastromethatacnoabsdistance) <- gastroatacsig

for(i in 1:length(gastromethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(gastroatacsig)){
    gastromethatacnoabsdistance[i,j] <- (gastroatacsigpeakanno[gastroatacsig[j],"chrom"] == fingastromethsiganno$chrom[i])*(gastroatacsigpeakanno[gastroatacsig[j],"mid"] - fingastromethsiganno$sitemid[i])
  }
}


# HEART

heartmethataccor <- matrix(0L,nrow = dim(heartmethatacdistance)[1],ncol = dim(heartmethatacdistance)[2])
rownames(heartmethataccor) <- rownames(heartmethatacdistance)
colnames(heartmethataccor) <- colnames(heartmethatacdistance)
for(i in 1:dim(heartmethataccor)[1]){
  if(i%%20 == 0){
    print(i)
  }
  for(j in 1:dim(heartmethataccor)[2]){
    ourmeth <- rownames(heartmethatacdistance)[i]
    ourpeak <- colnames(heartmethatacdistance)[j]
    heartmethataccor[i,j] <- cor(heartmethl2fc[ourmeth,],heartsigatacl2fc[ourpeak,])
  }
}
heartmethataccor[is.na(heartmethataccor)] <- 0

heartmethatacnoabsdistance <- matrix(0L,nrow = length(heartmethsig),ncol = length(heartatacsig))
rownames(heartmethatacnoabsdistance) <- heartmethsig
colnames(heartmethatacnoabsdistance) <- heartatacsig

for(i in 1:length(heartmethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(heartatacsig)){
    heartmethatacnoabsdistance[i,j] <- (heartatacsigpeakanno[heartatacsig[j],"chrom"] == finheartmethsiganno$chrom[i])*(heartatacsigpeakanno[heartatacsig[j],"mid"] - finheartmethsiganno$sitemid[i])
  }
}


# HIPPOC

hippomethataccor <- matrix(0L,nrow = dim(hippomethatacdistance)[1],ncol = dim(hippomethatacdistance)[2])
rownames(hippomethataccor) <- rownames(hippomethatacdistance)
colnames(hippomethataccor) <- colnames(hippomethatacdistance)
for(i in 1:dim(hippomethataccor)[1]){
  if(i%%20 == 0){
    print(i)
  }
  for(j in 1:dim(hippomethataccor)[2]){
    ourmeth <- rownames(hippomethatacdistance)[i]
    ourpeak <- colnames(hippomethatacdistance)[j]
    hippomethataccor[i,j] <- cor(hippomethl2fc[ourmeth,],hipposigatacl2fc[ourpeak,])
  }
}
hippomethataccor[is.na(hippomethataccor)] <- 0

hippomethatacnoabsdistance <- matrix(0L,nrow = length(hippomethsig),ncol = length(hippoatacsig))
rownames(hippomethatacnoabsdistance) <- hippomethsig
colnames(hippomethatacnoabsdistance) <- hippoatacsig

for(i in 1:length(hippomethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(hippoatacsig)){
    hippomethatacnoabsdistance[i,j] <- (hippoatacsigpeakanno[hippoatacsig[j],"chrom"] == finhippomethsiganno$chrom[i])*(hippoatacsigpeakanno[hippoatacsig[j],"mid"] - finhippomethsiganno$sitemid[i])
  }
}


# KIDNEY

kidneymethataccor <- matrix(0L,nrow = dim(kidneymethatacdistance)[1],ncol = dim(kidneymethatacdistance)[2])
rownames(kidneymethataccor) <- rownames(kidneymethatacdistance)
colnames(kidneymethataccor) <- colnames(kidneymethatacdistance)
for(i in 1:dim(kidneymethataccor)[1]){
  if(i%%20 == 0){
    print(i)
  }
  for(j in 1:dim(kidneymethataccor)[2]){
    ourmeth <- rownames(kidneymethatacdistance)[i]
    ourpeak <- colnames(kidneymethatacdistance)[j]
    kidneymethataccor[i,j] <- cor(kidneymethl2fc[ourmeth,],kidneysigatacl2fc[ourpeak,])
  }
}
kidneymethataccor[is.na(kidneymethataccor)] <- 0

kidneymethatacnoabsdistance <- matrix(0L,nrow = length(kidneymethsig),ncol = length(kidneyatacsig))
rownames(kidneymethatacnoabsdistance) <- kidneymethsig
colnames(kidneymethatacnoabsdistance) <- kidneyatacsig

for(i in 1:length(kidneymethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(kidneyatacsig)){
    kidneymethatacnoabsdistance[i,j] <- (kidneyatacsigpeakanno[kidneyatacsig[j],"chrom"] == finkidneymethsiganno$chrom[i])*(kidneyatacsigpeakanno[kidneyatacsig[j],"mid"] - finkidneymethsiganno$sitemid[i])
  }
}


# LIVER

livermethataccor <- matrix(0L,nrow = dim(livermethatacdistance)[1],ncol = dim(livermethatacdistance)[2])
rownames(livermethataccor) <- rownames(livermethatacdistance)
colnames(livermethataccor) <- colnames(livermethatacdistance)
for(i in 1:dim(livermethataccor)[1]){
  if(i%%20 == 0){
    print(i)
  }
  for(j in 1:dim(livermethataccor)[2]){
    ourmeth <- rownames(livermethatacdistance)[i]
    ourpeak <- colnames(livermethatacdistance)[j]
    livermethataccor[i,j] <- cor(livermethl2fc[ourmeth,],liversigatacl2fc[ourpeak,])
  }
}
livermethataccor[is.na(livermethataccor)] <- 0

livermethatacnoabsdistance <- matrix(0L,nrow = length(livermethsig),ncol = length(liveratacsig))
rownames(livermethatacnoabsdistance) <- livermethsig
colnames(livermethatacnoabsdistance) <- liveratacsig

for(i in 1:length(livermethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(liveratacsig)){
    livermethatacnoabsdistance[i,j] <- (liveratacsigpeakanno[liveratacsig[j],"chrom"] == finlivermethsiganno$chrom[i])*(liveratacsigpeakanno[liveratacsig[j],"mid"] - finlivermethsiganno$sitemid[i])
  }
}


# LUNG

lungmethataccor <- matrix(0L,nrow = dim(lungmethatacdistance)[1],ncol = dim(lungmethatacdistance)[2])
rownames(lungmethataccor) <- rownames(lungmethatacdistance)
colnames(lungmethataccor) <- colnames(lungmethatacdistance)
for(i in 1:dim(lungmethataccor)[1]){
  if(i%%20 == 0){
    print(i)
  }
  for(j in 1:dim(lungmethataccor)[2]){
    ourmeth <- rownames(lungmethatacdistance)[i]
    ourpeak <- colnames(lungmethatacdistance)[j]
    lungmethataccor[i,j] <- cor(lungmethl2fc[ourmeth,],lungsigatacl2fc[ourpeak,])
  }
}
lungmethataccor[is.na(lungmethataccor)] <- 0

lungmethatacnoabsdistance <- matrix(0L,nrow = length(lungmethsig),ncol = length(lungatacsig))
rownames(lungmethatacnoabsdistance) <- lungmethsig
colnames(lungmethatacnoabsdistance) <- lungatacsig

for(i in 1:length(lungmethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(lungatacsig)){
    lungmethatacnoabsdistance[i,j] <- (lungatacsigpeakanno[lungatacsig[j],"chrom"] == finlungmethsiganno$chrom[i])*(lungatacsigpeakanno[lungatacsig[j],"mid"] - finlungmethsiganno$sitemid[i])
  }
}


# BAT

brownmethataccor <- matrix(0L,nrow = dim(brownmethatacdistance)[1],ncol = dim(brownmethatacdistance)[2])
rownames(brownmethataccor) <- rownames(brownmethatacdistance)
colnames(brownmethataccor) <- colnames(brownmethatacdistance)
for(i in 1:dim(brownmethataccor)[1]){
  if(i%%20 == 0){
    print(i)
  }
  for(j in 1:dim(brownmethataccor)[2]){
    ourmeth <- rownames(brownmethatacdistance)[i]
    ourpeak <- colnames(brownmethatacdistance)[j]
    brownmethataccor[i,j] <- cor(brownmethl2fc[ourmeth,],brownsigatacl2fc[ourpeak,])
  }
}
brownmethataccor[is.na(brownmethataccor)] <- 0

brownmethatacnoabsdistance <- matrix(0L,nrow = length(brownmethsig),ncol = length(brownatacsig))
rownames(brownmethatacnoabsdistance) <- brownmethsig
colnames(brownmethatacnoabsdistance) <- brownatacsig

for(i in 1:length(brownmethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(brownatacsig)){
    brownmethatacnoabsdistance[i,j] <- (brownatacsigpeakanno[brownatacsig[j],"chrom"] == finbrownmethsiganno$chrom[i])*(brownatacsigpeakanno[brownatacsig[j],"mid"] - finbrownmethsiganno$sitemid[i])
  }
}


# WAT-SC

whitemethataccor <- matrix(0L,nrow = dim(whitemethatacdistance)[1],ncol = dim(whitemethatacdistance)[2])
rownames(whitemethataccor) <- rownames(whitemethatacdistance)
colnames(whitemethataccor) <- colnames(whitemethatacdistance)
for(i in 1:dim(whitemethataccor)[1]){
  if(i%%20 == 0){
    print(i)
  }
  for(j in 1:dim(whitemethataccor)[2]){
    ourmeth <- rownames(whitemethatacdistance)[i]
    ourpeak <- colnames(whitemethatacdistance)[j]
    whitemethataccor[i,j] <- cor(whitemethl2fc[ourmeth,],whitesigatacl2fc[ourpeak,])
  }
}
whitemethataccor[is.na(whitemethataccor)] <- 0

whitemethatacnoabsdistance <- matrix(0L,nrow = length(whitemethsig),ncol = length(whiteatacsig))
rownames(whitemethatacnoabsdistance) <- whitemethsig
colnames(whitemethatacnoabsdistance) <- whiteatacsig

for(i in 1:length(whitemethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(whiteatacsig)){
    whitemethatacnoabsdistance[i,j] <- (whiteatacsigpeakanno[whiteatacsig[j],"chrom"] == finwhitemethsiganno$chrom[i])*(whiteatacsigpeakanno[whiteatacsig[j],"mid"] - finwhitemethsiganno$sitemid[i])
  }
}




gastromethatacdistvscordf <- data.frame("Correlation" = as.vector(gastromethataccor),"Distance" = as.vector(gastromethatacnoabsdistance))
heartmethatacdistvscordf <- data.frame("Correlation" = as.vector(heartmethataccor),"Distance" = as.vector(heartmethatacnoabsdistance))
hippomethatacdistvscordf <- data.frame("Correlation" = as.vector(hippomethataccor),"Distance" = as.vector(hippomethatacnoabsdistance))
kidneymethatacdistvscordf <- data.frame("Correlation" = as.vector(kidneymethataccor),"Distance" = as.vector(kidneymethatacnoabsdistance))
livermethatacdistvscordf <- data.frame("Correlation" = as.vector(livermethataccor),"Distance" = as.vector(livermethatacnoabsdistance))
lungmethatacdistvscordf <- data.frame("Correlation" = as.vector(lungmethataccor),"Distance" = as.vector(lungmethatacnoabsdistance))
brownmethatacdistvscordf <- data.frame("Correlation" = as.vector(brownmethataccor),"Distance" = as.vector(brownmethatacnoabsdistance))
whitemethatacdistvscordf <- data.frame("Correlation" = as.vector(whitemethataccor),"Distance" = as.vector(whitemethatacnoabsdistance))

totalmethatacdistvscordf <- rbind(gastromethatacdistvscordf,heartmethatacdistvscordf,hippomethatacdistvscordf,kidneymethatacdistvscordf,
                                  livermethatacdistvscordf,lungmethatacdistvscordf,whitemethatacdistvscordf,brownmethatacdistvscordf)
totalmethatacdistvscordf$Tissue <- c(rep("SKM-GN",dim(gastromethatacdistvscordf)[1]),
                                     rep("HEART",dim(heartmethatacdistvscordf)[1]),
                                     rep("HIPPOC",dim(hippomethatacdistvscordf)[1]),
                                     rep("KIDNEY",dim(kidneymethatacdistvscordf)[1]),
                                     rep("LIVER",dim(livermethatacdistvscordf)[1]),
                                     rep("LUNG",dim(lungmethatacdistvscordf)[1]),
                                     rep("WAT-SC",dim(whitemethatacdistvscordf)[1]),
                                     rep("BAT",dim(brownmethatacdistvscordf)[1]))

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



# New Figure 3H
pdf(file = "New Figure 3H.pdf",width = 9,height = 6)
ggplot(totalmethatacdistvscordf[abs(totalmethatacdistvscordf$Distance) > 0 & abs(totalmethatacdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to Nearest DAR") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_density2d_filled() + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 18),legend.title = element_text(size = 18),legend.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + geom_point(data = totalmethatacdistvscordf[abs(totalmethatacdistvscordf$Distance) > 0 & abs(totalmethatacdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation,color=Tissue),size = 2) + scale_color_manual(values = ann_cols$Tissue)
dev.off()

png(file = "New Figure 3H.png",width = 9,height = 6,units = "in",res = 600)
ggplot(totalmethatacdistvscordf[abs(totalmethatacdistvscordf$Distance) > 0 & abs(totalmethatacdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to Nearest DAR") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_density2d_filled() + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 18),legend.title = element_text(size = 18),legend.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + geom_point(data = totalmethatacdistvscordf[abs(totalmethatacdistvscordf$Distance) > 0 & abs(totalmethatacdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation,color=Tissue),size = 2) + scale_color_manual(values = ann_cols$Tissue)
dev.off()


# New Supplemental Figure S10

pdf(file = "Supplemental Figure S10A.pdf",width = 7,height = 7)
ggplot(gastromethatacdistvscordf[abs(gastromethatacdistvscordf$Distance) > 0 & abs(gastromethatacdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to Nearest DAR") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["SKM-GN"]) + geom_density2d(color = ann_cols$Tissue["SKM-GN"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["SKM-GN"])
dev.off()

pdf(file = "Supplemental Figure S10B.pdf",width = 7,height = 7)
ggplot(heartmethatacdistvscordf[abs(heartmethatacdistvscordf$Distance) > 0 & abs(heartmethatacdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to Nearest DAR") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["HEART"]) + geom_density2d(color = ann_cols$Tissue["HEART"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["HEART"])
dev.off()

pdf(file = "Supplemental Figure S10C.pdf",width = 7,height = 7)
ggplot(kidneymethatacdistvscordf[abs(kidneymethatacdistvscordf$Distance) > 0 & abs(kidneymethatacdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to Nearest DAR") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["KIDNEY"]) + geom_density2d(color = ann_cols$Tissue["KIDNEY"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["KIDNEY"])
dev.off()

pdf(file = "Supplemental Figure S10D.pdf",width = 7,height = 7)
ggplot(livermethatacdistvscordf[abs(livermethatacdistvscordf$Distance) > 0 & abs(livermethatacdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to Nearest DAR") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["LIVER"]) + geom_density2d(color = ann_cols$Tissue["LIVER"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["LIVER"])
dev.off()

pdf(file = "Supplemental Figure S10E.pdf",width = 7,height = 7)
ggplot(lungmethatacdistvscordf[abs(lungmethatacdistvscordf$Distance) > 0 & abs(lungmethatacdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to Nearest DAR") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["LUNG"]) + geom_density2d(color = ann_cols$Tissue["LUNG"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["LUNG"])
dev.off()

pdf(file = "Supplemental Figure S10F.pdf",width = 7,height = 7)
ggplot(brownmethatacdistvscordf[abs(brownmethatacdistvscordf$Distance) > 0 & abs(brownmethatacdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to Nearest DAR") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["BAT"]) + geom_density2d(color = ann_cols$Tissue["BAT"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["BAT"])
dev.off()

png(file = "Supplemental Figure S10A.png",width = 7,height = 7,units = "in",res = 600)
ggplot(gastromethatacdistvscordf[abs(gastromethatacdistvscordf$Distance) > 0 & abs(gastromethatacdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to Nearest DAR") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["SKM-GN"]) + geom_density2d(color = ann_cols$Tissue["SKM-GN"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["SKM-GN"])
dev.off()

png(file = "Supplemental Figure S10B.png",width = 7,height = 7,units = "in",res = 600)
ggplot(heartmethatacdistvscordf[abs(heartmethatacdistvscordf$Distance) > 0 & abs(heartmethatacdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to Nearest DAR") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["HEART"]) + geom_density2d(color = ann_cols$Tissue["HEART"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["HEART"])
dev.off()

png(file = "Supplemental Figure S10C.png",width = 7,height = 7,units = "in",res = 600)
ggplot(kidneymethatacdistvscordf[abs(kidneymethatacdistvscordf$Distance) > 0 & abs(kidneymethatacdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to Nearest DAR") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["KIDNEY"]) + geom_density2d(color = ann_cols$Tissue["KIDNEY"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["KIDNEY"])
dev.off()

png(file = "Supplemental Figure S10D.png",width = 7,height = 7,units = "in",res = 600)
ggplot(livermethatacdistvscordf[abs(livermethatacdistvscordf$Distance) > 0 & abs(livermethatacdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to Nearest DAR") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["LIVER"]) + geom_density2d(color = ann_cols$Tissue["LIVER"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["LIVER"])
dev.off()

png(file = "Supplemental Figure S10E.png",width = 7,height = 7,units = "in",res = 600)
ggplot(lungmethatacdistvscordf[abs(lungmethatacdistvscordf$Distance) > 0 & abs(lungmethatacdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to Nearest DAR") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["LUNG"]) + geom_density2d(color = ann_cols$Tissue["LUNG"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["LUNG"])
dev.off()

png(file = "Supplemental Figure S10F.png",width = 7,height = 7,units = "in",res = 600)
ggplot(brownmethatacdistvscordf[abs(brownmethatacdistvscordf$Distance) > 0 & abs(brownmethatacdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to Nearest DAR") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["BAT"]) + geom_density2d(color = ann_cols$Tissue["BAT"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["BAT"])
dev.off()


save.image("Figure3E_3H_S10.RData")