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

pdf(file = "Figure 3E_112624.pdf",width = 7,height = 6)
ggplot(totalmethatacdistancedf, aes(x = Distance)) +
  geom_histogram(aes(color = Tissue,fill = Tissue),
                 position = "stack", bins = 30) + theme_classic() +
  scale_color_manual(values = c("#8c5220","#f28b2f","#bf7534","#7553a7","#da6c75","#04bf8a","#088c03","#214da6")) + 
  scale_fill_manual(values = c("#8c5220","#f28b2f","#bf7534","#7553a7","#da6c75","#04bf8a","#088c03","#214da6"))  + scale_x_log10(breaks=c(1,100,10000,1000000,100000000)) + xlab("Distance to Nearest DAR") + ylab("Count") + theme(axis.title = element_text(size = 15),axis.text = element_text(size = 12),legend.title = element_text(size = 15),legend.text = element_text(size = 12))
dev.off()

png(file = "Figure 3E_112624.png",width = 7,height = 6,units = "in",res = 600)
ggplot(totalmethatacdistancedf, aes(x = Distance)) +
  geom_histogram(aes(color = Tissue,fill = Tissue),
                 position = "stack", bins = 30) + theme_classic() +
  scale_color_manual(values = c("#8c5220","#f28b2f","#bf7534","#7553a7","#da6c75","#04bf8a","#088c03","#214da6")) + 
  scale_fill_manual(values = c("#8c5220","#f28b2f","#bf7534","#7553a7","#da6c75","#04bf8a","#088c03","#214da6"))  + scale_x_log10(breaks=c(1,100,10000,1000000,100000000)) + xlab("Distance to Nearest DAR") + ylab("Count") + theme(axis.title = element_text(size = 15),axis.text = element_text(size = 12),legend.title = element_text(size = 15),legend.text = element_text(size = 12))
dev.off()

#####
# Figure 3H
####

#load("atacsigl2fcmat.RData")
#load("methsigl2fcmat.RData")

load("new_atacnormmatrices.RData")
load("methMmatrices.RData")

pass1bphenodata <- readRDS("pass1bphenodata.rds")

for(i in 1:dim(gastroatacnorm)[2]){
  ourviallabel <- colnames(gastroatacnorm)[i]
  ourpid <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourviallabel,"pass1bf0001"][1]
  colnames(gastroatacnorm)[i] <- ourpid
}

for(i in 1:dim(heartatacnorm)[2]){
  ourviallabel <- colnames(heartatacnorm)[i]
  ourpid <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourviallabel,"pass1bf0001"][1]
  colnames(heartatacnorm)[i] <- ourpid
}

for(i in 1:dim(hippoatacnorm)[2]){
  ourviallabel <- colnames(hippoatacnorm)[i]
  ourpid <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourviallabel,"pass1bf0001"][1]
  colnames(hippoatacnorm)[i] <- ourpid
}

for(i in 1:dim(kidneyatacnorm)[2]){
  ourviallabel <- colnames(kidneyatacnorm)[i]
  ourpid <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourviallabel,"pass1bf0001"][1]
  colnames(kidneyatacnorm)[i] <- ourpid
}

for(i in 1:dim(liveratacnorm)[2]){
  ourviallabel <- colnames(liveratacnorm)[i]
  ourpid <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourviallabel,"pass1bf0001"][1]
  colnames(liveratacnorm)[i] <- ourpid
}

for(i in 1:dim(lungatacnorm)[2]){
  ourviallabel <- colnames(lungatacnorm)[i]
  ourpid <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourviallabel,"pass1bf0001"][1]
  colnames(lungatacnorm)[i] <- ourpid
}

for(i in 1:dim(brownatacnorm)[2]){
  ourviallabel <- colnames(brownatacnorm)[i]
  ourpid <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourviallabel,"pass1bf0001"][1]
  colnames(brownatacnorm)[i] <- ourpid
}

for(i in 1:dim(whiteatacnorm)[2]){
  ourviallabel <- colnames(whiteatacnorm)[i]
  ourpid <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourviallabel,"pass1bf0001"][1]
  colnames(whiteatacnorm)[i] <- ourpid
}


gastroMsigmat <- gastroMmat[gastromethsig,]
heartMsigmat <- heartMmat[heartmethsig,]
hippoMsigmat <- hippoMmat[hippomethsig,]
kidneyMsigmat <- kidneyMmat[kidneymethsig,]
liverMsigmat <- liverMmat[livermethsig,]
lungMsigmat <- lungMmat[lungmethsig,]
brownMsigmat <- brownMmat[brownmethsig,]
whiteMsigmat <- whiteMmat[whitemethsig,]

pidmeta <- data.frame(row.names = colnames(gastroMsigmat),
                      "pid" = colnames(gastroMsigmat))
pidmeta$sex <- ""
pidmeta$group <- ""

for(i in 1:dim(pidmeta)[1]){
  ourpid <- pidmeta$pid[i]
  pidmeta[i,"sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourpid,"pass1bf0027"][1]]
  pidmeta[i,"group"] <- pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourpid,"pass1bf0011"][1]
}
pidmeta[pidmeta$group %in% "One-week program","group"] <- "1w"
pidmeta[pidmeta$group %in% "Two-week program","group"] <- "2w"
pidmeta[pidmeta$group %in% "Four-week program","group"] <- "4w"
pidmeta[pidmeta$group %in% "Eight-week program Training Group","group"] <- "8w"
pidmeta[pidmeta$group %in% "Eight-week program Control Group","group"] <- "control"



# gastro

gastroMsigmatrix <- matrix(as.numeric(as.matrix(gastroMsigmat)),nrow = dim(gastroMsigmat)[1],ncol = dim(gastroMsigmat)[2])
rownames(gastroMsigmatrix) <- rownames(gastroMsigmat)
colnames(gastroMsigmatrix) <- colnames(gastroMsigmat)

gastrocols <- intersect(colnames(gastroatacnorm),colnames(gastroMsigmatrix))
gastroatacnorm <- gastroatacnorm[,gastrocols]
gastroMsigmatrix <- gastroMsigmatrix[,gastrocols]

gastroatacsignorm <- gastroatacnorm[gastroatacsig,]

gastroatacsigfemalecontrolavg <- apply(gastroatacsignorm[,intersect(colnames(gastroatacsignorm),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% "control",]))],1,mean)
gastroatacsigmalecontrolavg <- apply(gastroatacsignorm[,intersect(colnames(gastroatacsignorm),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% "control",]))],1,mean)

gastroatacfemaletrainingl2fcbysubject <- gastroatacsignorm[,intersect(colnames(gastroatacsignorm),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(gastroatacfemaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(gastroatacfemaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(gastroatacfemaletrainingl2fcbysubject)[j]
    gastroatacfemaletrainingl2fcbysubject[i,j] <- gastroatacsignorm[i,ourcolname]-gastroatacsigfemalecontrolavg[i]
  }
}
gastroatacmaletrainingl2fcbysubject <- gastroatacsignorm[,intersect(colnames(gastroatacsignorm),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(gastroatacmaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(gastroatacmaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(gastroatacmaletrainingl2fcbysubject)[j]
    gastroatacmaletrainingl2fcbysubject[i,j] <- gastroatacsignorm[i,ourcolname]-gastroatacsigmalecontrolavg[i]
  }
}
gastroatactrainingl2fcbysubject <- cbind(gastroatacfemaletrainingl2fcbysubject,gastroatacmaletrainingl2fcbysubject)

gastromethsigfemalecontrolavg <- apply(gastroMsigmatrix[,intersect(colnames(gastroMsigmatrix),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% "control",]))],1,mean)
gastromethsigmalecontrolavg <- apply(gastroMsigmatrix[,intersect(colnames(gastroMsigmatrix),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% "control",]))],1,mean)

gastromethfemaletrainingl2fcbysubject <- gastroMsigmatrix[,intersect(colnames(gastroMsigmatrix),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(gastromethfemaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(gastromethfemaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(gastromethfemaletrainingl2fcbysubject)[j]
    gastromethfemaletrainingl2fcbysubject[i,j] <- gastroMsigmatrix[i,ourcolname]-gastromethsigfemalecontrolavg[i]
  }
}
gastromethmaletrainingl2fcbysubject <- gastroMsigmatrix[,intersect(colnames(gastroMsigmatrix),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(gastromethmaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(gastromethmaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(gastromethmaletrainingl2fcbysubject)[j]
    gastromethmaletrainingl2fcbysubject[i,j] <- gastroMsigmatrix[i,ourcolname]-gastromethsigmalecontrolavg[i]
  }
}
gastromethtrainingl2fcbysubject <- cbind(gastromethfemaletrainingl2fcbysubject,gastromethmaletrainingl2fcbysubject)


gastromethataccor <- matrix(0L,nrow = dim(gastromethatacdistance)[1],ncol = dim(gastromethatacdistance)[2])
rownames(gastromethataccor) <- rownames(gastromethatacdistance)
colnames(gastromethataccor) <- colnames(gastromethatacdistance)
gastromethataccortest <- matrix(0L,nrow = dim(gastromethatacdistance)[1],ncol = dim(gastromethatacdistance)[2])
rownames(gastromethataccortest) <- rownames(gastromethatacdistance)
colnames(gastromethataccortest) <- colnames(gastromethatacdistance)
for(i in 1:dim(gastromethataccor)[1]){
  if(i%%20 == 0){
    print(i)
  }
  for(j in 1:dim(gastromethataccor)[2]){
    ourmeth <- rownames(gastromethatacdistance)[i]
    ourpeak <- colnames(gastromethatacdistance)[j]
    ourout <- cor.test(gastromethtrainingl2fcbysubject[ourmeth,],gastroatactrainingl2fcbysubject[ourpeak,])
    gastromethataccor[i,j] <- ourout$estimate
    gastromethataccortest[i,j] <- ourout$p.value
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


# heart

heartMsigmatrix <- matrix(as.numeric(as.matrix(heartMsigmat)),nrow = dim(heartMsigmat)[1],ncol = dim(heartMsigmat)[2])
rownames(heartMsigmatrix) <- rownames(heartMsigmat)
colnames(heartMsigmatrix) <- colnames(heartMsigmat)

heartcols <- intersect(colnames(heartatacnorm),colnames(heartMsigmatrix))
heartatacnorm <- heartatacnorm[,heartcols]
heartMsigmatrix <- heartMsigmatrix[,heartcols]

heartatacsignorm <- heartatacnorm[heartatacsig,]

heartatacsigfemalecontrolavg <- apply(heartatacsignorm[,intersect(colnames(heartatacsignorm),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% "control",]))],1,mean)
heartatacsigmalecontrolavg <- apply(heartatacsignorm[,intersect(colnames(heartatacsignorm),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% "control",]))],1,mean)

heartatacfemaletrainingl2fcbysubject <- heartatacsignorm[,intersect(colnames(heartatacsignorm),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(heartatacfemaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(heartatacfemaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(heartatacfemaletrainingl2fcbysubject)[j]
    heartatacfemaletrainingl2fcbysubject[i,j] <- heartatacsignorm[i,ourcolname]-heartatacsigfemalecontrolavg[i]
  }
}
heartatacmaletrainingl2fcbysubject <- heartatacsignorm[,intersect(colnames(heartatacsignorm),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(heartatacmaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(heartatacmaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(heartatacmaletrainingl2fcbysubject)[j]
    heartatacmaletrainingl2fcbysubject[i,j] <- heartatacsignorm[i,ourcolname]-heartatacsigmalecontrolavg[i]
  }
}
heartatactrainingl2fcbysubject <- cbind(heartatacfemaletrainingl2fcbysubject,heartatacmaletrainingl2fcbysubject)

heartmethsigfemalecontrolavg <- apply(heartMsigmatrix[,intersect(colnames(heartMsigmatrix),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% "control",]))],1,mean)
heartmethsigmalecontrolavg <- apply(heartMsigmatrix[,intersect(colnames(heartMsigmatrix),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% "control",]))],1,mean)

heartmethfemaletrainingl2fcbysubject <- heartMsigmatrix[,intersect(colnames(heartMsigmatrix),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(heartmethfemaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(heartmethfemaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(heartmethfemaletrainingl2fcbysubject)[j]
    heartmethfemaletrainingl2fcbysubject[i,j] <- heartMsigmatrix[i,ourcolname]-heartmethsigfemalecontrolavg[i]
  }
}
heartmethmaletrainingl2fcbysubject <- heartMsigmatrix[,intersect(colnames(heartMsigmatrix),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(heartmethmaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(heartmethmaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(heartmethmaletrainingl2fcbysubject)[j]
    heartmethmaletrainingl2fcbysubject[i,j] <- heartMsigmatrix[i,ourcolname]-heartmethsigmalecontrolavg[i]
  }
}
heartmethtrainingl2fcbysubject <- cbind(heartmethfemaletrainingl2fcbysubject,heartmethmaletrainingl2fcbysubject)


heartmethataccor <- matrix(0L,nrow = dim(heartmethatacdistance)[1],ncol = dim(heartmethatacdistance)[2])
rownames(heartmethataccor) <- rownames(heartmethatacdistance)
colnames(heartmethataccor) <- colnames(heartmethatacdistance)
heartmethataccortest <- matrix(0L,nrow = dim(heartmethatacdistance)[1],ncol = dim(heartmethatacdistance)[2])
rownames(heartmethataccortest) <- rownames(heartmethatacdistance)
colnames(heartmethataccortest) <- colnames(heartmethatacdistance)
for(i in 1:dim(heartmethataccor)[1]){
  if(i%%20 == 0){
    print(i)
  }
  for(j in 1:dim(heartmethataccor)[2]){
    ourmeth <- rownames(heartmethatacdistance)[i]
    ourpeak <- colnames(heartmethatacdistance)[j]
    ourout <- cor.test(heartmethtrainingl2fcbysubject[ourmeth,],heartatactrainingl2fcbysubject[ourpeak,])
    heartmethataccor[i,j] <- ourout$estimate
    heartmethataccortest[i,j] <- ourout$p.value
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



# hippo

hippoMsigmatrix <- matrix(as.numeric(as.matrix(hippoMsigmat)),nrow = dim(hippoMsigmat)[1],ncol = dim(hippoMsigmat)[2])
rownames(hippoMsigmatrix) <- rownames(hippoMsigmat)
colnames(hippoMsigmatrix) <- colnames(hippoMsigmat)

hippocols <- intersect(colnames(hippoatacnorm),colnames(hippoMsigmatrix))
hippoatacnorm <- hippoatacnorm[,hippocols]
hippoMsigmatrix <- hippoMsigmatrix[,hippocols]

hippoatacsignorm <- hippoatacnorm[hippoatacsig,]

hippoatacsigfemalecontrolavg <- apply(hippoatacsignorm[,intersect(colnames(hippoatacsignorm),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% "control",]))],1,mean)
hippoatacsigmalecontrolavg <- apply(hippoatacsignorm[,intersect(colnames(hippoatacsignorm),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% "control",]))],1,mean)

hippoatacfemaletrainingl2fcbysubject <- hippoatacsignorm[,intersect(colnames(hippoatacsignorm),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(hippoatacfemaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(hippoatacfemaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(hippoatacfemaletrainingl2fcbysubject)[j]
    hippoatacfemaletrainingl2fcbysubject[i,j] <- hippoatacsignorm[i,ourcolname]-hippoatacsigfemalecontrolavg[i]
  }
}
hippoatacmaletrainingl2fcbysubject <- hippoatacsignorm[,intersect(colnames(hippoatacsignorm),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(hippoatacmaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(hippoatacmaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(hippoatacmaletrainingl2fcbysubject)[j]
    hippoatacmaletrainingl2fcbysubject[i,j] <- hippoatacsignorm[i,ourcolname]-hippoatacsigmalecontrolavg[i]
  }
}
hippoatactrainingl2fcbysubject <- cbind(hippoatacfemaletrainingl2fcbysubject,hippoatacmaletrainingl2fcbysubject)

hippomethsigfemalecontrolavg <- apply(hippoMsigmatrix[,intersect(colnames(hippoMsigmatrix),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% "control",]))],1,mean)
hippomethsigmalecontrolavg <- apply(hippoMsigmatrix[,intersect(colnames(hippoMsigmatrix),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% "control",]))],1,mean)

hippomethfemaletrainingl2fcbysubject <- hippoMsigmatrix[,intersect(colnames(hippoMsigmatrix),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(hippomethfemaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(hippomethfemaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(hippomethfemaletrainingl2fcbysubject)[j]
    hippomethfemaletrainingl2fcbysubject[i,j] <- hippoMsigmatrix[i,ourcolname]-hippomethsigfemalecontrolavg[i]
  }
}
hippomethmaletrainingl2fcbysubject <- hippoMsigmatrix[,intersect(colnames(hippoMsigmatrix),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(hippomethmaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(hippomethmaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(hippomethmaletrainingl2fcbysubject)[j]
    hippomethmaletrainingl2fcbysubject[i,j] <- hippoMsigmatrix[i,ourcolname]-hippomethsigmalecontrolavg[i]
  }
}
hippomethtrainingl2fcbysubject <- cbind(hippomethfemaletrainingl2fcbysubject,hippomethmaletrainingl2fcbysubject)


hippomethataccor <- matrix(0L,nrow = dim(hippomethatacdistance)[1],ncol = dim(hippomethatacdistance)[2])
rownames(hippomethataccor) <- rownames(hippomethatacdistance)
colnames(hippomethataccor) <- colnames(hippomethatacdistance)
hippomethataccortest <- matrix(0L,nrow = dim(hippomethatacdistance)[1],ncol = dim(hippomethatacdistance)[2])
rownames(hippomethataccortest) <- rownames(hippomethatacdistance)
colnames(hippomethataccortest) <- colnames(hippomethatacdistance)
for(i in 1:dim(hippomethataccor)[1]){
  if(i%%20 == 0){
    print(i)
  }
  for(j in 1:dim(hippomethataccor)[2]){
    ourmeth <- rownames(hippomethatacdistance)[i]
    ourpeak <- colnames(hippomethatacdistance)[j]
    ourout <- cor.test(hippomethtrainingl2fcbysubject[ourmeth,],hippoatactrainingl2fcbysubject[ourpeak,])
    hippomethataccor[i,j] <- ourout$estimate
    hippomethataccortest[i,j] <- ourout$p.value
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



# kidney

kidneyMsigmatrix <- matrix(as.numeric(as.matrix(kidneyMsigmat)),nrow = dim(kidneyMsigmat)[1],ncol = dim(kidneyMsigmat)[2])
rownames(kidneyMsigmatrix) <- rownames(kidneyMsigmat)
colnames(kidneyMsigmatrix) <- colnames(kidneyMsigmat)

kidneycols <- intersect(colnames(kidneyatacnorm),colnames(kidneyMsigmatrix))
kidneyatacnorm <- kidneyatacnorm[,kidneycols]
kidneyMsigmatrix <- kidneyMsigmatrix[,kidneycols]

kidneyatacsignorm <- kidneyatacnorm[kidneyatacsig,]

kidneyatacsigfemalecontrolavg <- apply(kidneyatacsignorm[,intersect(colnames(kidneyatacsignorm),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% "control",]))],1,mean)
kidneyatacsigmalecontrolavg <- apply(kidneyatacsignorm[,intersect(colnames(kidneyatacsignorm),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% "control",]))],1,mean)

kidneyatacfemaletrainingl2fcbysubject <- kidneyatacsignorm[,intersect(colnames(kidneyatacsignorm),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(kidneyatacfemaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(kidneyatacfemaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(kidneyatacfemaletrainingl2fcbysubject)[j]
    kidneyatacfemaletrainingl2fcbysubject[i,j] <- kidneyatacsignorm[i,ourcolname]-kidneyatacsigfemalecontrolavg[i]
  }
}
kidneyatacmaletrainingl2fcbysubject <- kidneyatacsignorm[,intersect(colnames(kidneyatacsignorm),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(kidneyatacmaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(kidneyatacmaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(kidneyatacmaletrainingl2fcbysubject)[j]
    kidneyatacmaletrainingl2fcbysubject[i,j] <- kidneyatacsignorm[i,ourcolname]-kidneyatacsigmalecontrolavg[i]
  }
}
kidneyatactrainingl2fcbysubject <- cbind(kidneyatacfemaletrainingl2fcbysubject,kidneyatacmaletrainingl2fcbysubject)

kidneymethsigfemalecontrolavg <- apply(kidneyMsigmatrix[,intersect(colnames(kidneyMsigmatrix),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% "control",]))],1,mean)
kidneymethsigmalecontrolavg <- apply(kidneyMsigmatrix[,intersect(colnames(kidneyMsigmatrix),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% "control",]))],1,mean)

kidneymethfemaletrainingl2fcbysubject <- kidneyMsigmatrix[,intersect(colnames(kidneyMsigmatrix),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(kidneymethfemaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(kidneymethfemaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(kidneymethfemaletrainingl2fcbysubject)[j]
    kidneymethfemaletrainingl2fcbysubject[i,j] <- kidneyMsigmatrix[i,ourcolname]-kidneymethsigfemalecontrolavg[i]
  }
}
kidneymethmaletrainingl2fcbysubject <- kidneyMsigmatrix[,intersect(colnames(kidneyMsigmatrix),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(kidneymethmaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(kidneymethmaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(kidneymethmaletrainingl2fcbysubject)[j]
    kidneymethmaletrainingl2fcbysubject[i,j] <- kidneyMsigmatrix[i,ourcolname]-kidneymethsigmalecontrolavg[i]
  }
}
kidneymethtrainingl2fcbysubject <- cbind(kidneymethfemaletrainingl2fcbysubject,kidneymethmaletrainingl2fcbysubject)


kidneymethataccor <- matrix(0L,nrow = dim(kidneymethatacdistance)[1],ncol = dim(kidneymethatacdistance)[2])
rownames(kidneymethataccor) <- rownames(kidneymethatacdistance)
colnames(kidneymethataccor) <- colnames(kidneymethatacdistance)
kidneymethataccortest <- matrix(0L,nrow = dim(kidneymethatacdistance)[1],ncol = dim(kidneymethatacdistance)[2])
rownames(kidneymethataccortest) <- rownames(kidneymethatacdistance)
colnames(kidneymethataccortest) <- colnames(kidneymethatacdistance)
for(i in 1:dim(kidneymethataccor)[1]){
  if(i%%20 == 0){
    print(i)
  }
  for(j in 1:dim(kidneymethataccor)[2]){
    ourmeth <- rownames(kidneymethatacdistance)[i]
    ourpeak <- colnames(kidneymethatacdistance)[j]
    ourout <- cor.test(kidneymethtrainingl2fcbysubject[ourmeth,],kidneyatactrainingl2fcbysubject[ourpeak,])
    kidneymethataccor[i,j] <- ourout$estimate
    kidneymethataccortest[i,j] <- ourout$p.value
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



# liver

liverMsigmatrix <- matrix(as.numeric(as.matrix(liverMsigmat)),nrow = dim(liverMsigmat)[1],ncol = dim(liverMsigmat)[2])
rownames(liverMsigmatrix) <- rownames(liverMsigmat)
colnames(liverMsigmatrix) <- colnames(liverMsigmat)

livercols <- intersect(colnames(liveratacnorm),colnames(liverMsigmatrix))
liveratacnorm <- liveratacnorm[,livercols]
liverMsigmatrix <- liverMsigmatrix[,livercols]

liveratacsignorm <- liveratacnorm[liveratacsig,]

liveratacsigfemalecontrolavg <- apply(liveratacsignorm[,intersect(colnames(liveratacsignorm),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% "control",]))],1,mean)
liveratacsigmalecontrolavg <- apply(liveratacsignorm[,intersect(colnames(liveratacsignorm),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% "control",]))],1,mean)

liveratacfemaletrainingl2fcbysubject <- liveratacsignorm[,intersect(colnames(liveratacsignorm),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(liveratacfemaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(liveratacfemaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(liveratacfemaletrainingl2fcbysubject)[j]
    liveratacfemaletrainingl2fcbysubject[i,j] <- liveratacsignorm[i,ourcolname]-liveratacsigfemalecontrolavg[i]
  }
}
liveratacmaletrainingl2fcbysubject <- liveratacsignorm[,intersect(colnames(liveratacsignorm),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(liveratacmaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(liveratacmaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(liveratacmaletrainingl2fcbysubject)[j]
    liveratacmaletrainingl2fcbysubject[i,j] <- liveratacsignorm[i,ourcolname]-liveratacsigmalecontrolavg[i]
  }
}
liveratactrainingl2fcbysubject <- cbind(liveratacfemaletrainingl2fcbysubject,liveratacmaletrainingl2fcbysubject)

livermethsigfemalecontrolavg <- apply(liverMsigmatrix[,intersect(colnames(liverMsigmatrix),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% "control",]))],1,mean)
livermethsigmalecontrolavg <- apply(liverMsigmatrix[,intersect(colnames(liverMsigmatrix),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% "control",]))],1,mean)

livermethfemaletrainingl2fcbysubject <- liverMsigmatrix[,intersect(colnames(liverMsigmatrix),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(livermethfemaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(livermethfemaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(livermethfemaletrainingl2fcbysubject)[j]
    livermethfemaletrainingl2fcbysubject[i,j] <- liverMsigmatrix[i,ourcolname]-livermethsigfemalecontrolavg[i]
  }
}
livermethmaletrainingl2fcbysubject <- liverMsigmatrix[,intersect(colnames(liverMsigmatrix),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(livermethmaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(livermethmaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(livermethmaletrainingl2fcbysubject)[j]
    livermethmaletrainingl2fcbysubject[i,j] <- liverMsigmatrix[i,ourcolname]-livermethsigmalecontrolavg[i]
  }
}
livermethtrainingl2fcbysubject <- cbind(livermethfemaletrainingl2fcbysubject,livermethmaletrainingl2fcbysubject)


livermethataccor <- matrix(0L,nrow = dim(livermethatacdistance)[1],ncol = dim(livermethatacdistance)[2])
rownames(livermethataccor) <- rownames(livermethatacdistance)
colnames(livermethataccor) <- colnames(livermethatacdistance)
livermethataccortest <- matrix(0L,nrow = dim(livermethatacdistance)[1],ncol = dim(livermethatacdistance)[2])
rownames(livermethataccortest) <- rownames(livermethatacdistance)
colnames(livermethataccortest) <- colnames(livermethatacdistance)
for(i in 1:dim(livermethataccor)[1]){
  if(i%%20 == 0){
    print(i)
  }
  for(j in 1:dim(livermethataccor)[2]){
    ourmeth <- rownames(livermethatacdistance)[i]
    ourpeak <- colnames(livermethatacdistance)[j]
    ourout <- cor.test(livermethtrainingl2fcbysubject[ourmeth,],liveratactrainingl2fcbysubject[ourpeak,])
    livermethataccor[i,j] <- ourout$estimate
    livermethataccortest[i,j] <- ourout$p.value
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



# lung

lungMsigmatrix <- matrix(as.numeric(as.matrix(lungMsigmat)),nrow = dim(lungMsigmat)[1],ncol = dim(lungMsigmat)[2])
rownames(lungMsigmatrix) <- rownames(lungMsigmat)
colnames(lungMsigmatrix) <- colnames(lungMsigmat)

lungcols <- intersect(colnames(lungatacnorm),colnames(lungMsigmatrix))
lungatacnorm <- lungatacnorm[,lungcols]
lungMsigmatrix <- lungMsigmatrix[,lungcols]

lungatacsignorm <- lungatacnorm[lungatacsig,]

lungatacsigfemalecontrolavg <- apply(lungatacsignorm[,intersect(colnames(lungatacsignorm),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% "control",]))],1,mean)
lungatacsigmalecontrolavg <- apply(lungatacsignorm[,intersect(colnames(lungatacsignorm),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% "control",]))],1,mean)

lungatacfemaletrainingl2fcbysubject <- lungatacsignorm[,intersect(colnames(lungatacsignorm),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(lungatacfemaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(lungatacfemaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(lungatacfemaletrainingl2fcbysubject)[j]
    lungatacfemaletrainingl2fcbysubject[i,j] <- lungatacsignorm[i,ourcolname]-lungatacsigfemalecontrolavg[i]
  }
}
lungatacmaletrainingl2fcbysubject <- lungatacsignorm[,intersect(colnames(lungatacsignorm),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(lungatacmaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(lungatacmaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(lungatacmaletrainingl2fcbysubject)[j]
    lungatacmaletrainingl2fcbysubject[i,j] <- lungatacsignorm[i,ourcolname]-lungatacsigmalecontrolavg[i]
  }
}
lungatactrainingl2fcbysubject <- cbind(lungatacfemaletrainingl2fcbysubject,lungatacmaletrainingl2fcbysubject)

lungmethsigfemalecontrolavg <- apply(lungMsigmatrix[,intersect(colnames(lungMsigmatrix),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% "control",]))],1,mean)
lungmethsigmalecontrolavg <- apply(lungMsigmatrix[,intersect(colnames(lungMsigmatrix),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% "control",]))],1,mean)

lungmethfemaletrainingl2fcbysubject <- lungMsigmatrix[,intersect(colnames(lungMsigmatrix),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(lungmethfemaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(lungmethfemaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(lungmethfemaletrainingl2fcbysubject)[j]
    lungmethfemaletrainingl2fcbysubject[i,j] <- lungMsigmatrix[i,ourcolname]-lungmethsigfemalecontrolavg[i]
  }
}
lungmethmaletrainingl2fcbysubject <- lungMsigmatrix[,intersect(colnames(lungMsigmatrix),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(lungmethmaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(lungmethmaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(lungmethmaletrainingl2fcbysubject)[j]
    lungmethmaletrainingl2fcbysubject[i,j] <- lungMsigmatrix[i,ourcolname]-lungmethsigmalecontrolavg[i]
  }
}
lungmethtrainingl2fcbysubject <- cbind(lungmethfemaletrainingl2fcbysubject,lungmethmaletrainingl2fcbysubject)


lungmethataccor <- matrix(0L,nrow = dim(lungmethatacdistance)[1],ncol = dim(lungmethatacdistance)[2])
rownames(lungmethataccor) <- rownames(lungmethatacdistance)
colnames(lungmethataccor) <- colnames(lungmethatacdistance)
lungmethataccortest <- matrix(0L,nrow = dim(lungmethatacdistance)[1],ncol = dim(lungmethatacdistance)[2])
rownames(lungmethataccortest) <- rownames(lungmethatacdistance)
colnames(lungmethataccortest) <- colnames(lungmethatacdistance)
for(i in 1:dim(lungmethataccor)[1]){
  if(i%%20 == 0){
    print(i)
  }
  for(j in 1:dim(lungmethataccor)[2]){
    ourmeth <- rownames(lungmethatacdistance)[i]
    ourpeak <- colnames(lungmethatacdistance)[j]
    ourout <- cor.test(lungmethtrainingl2fcbysubject[ourmeth,],lungatactrainingl2fcbysubject[ourpeak,])
    lungmethataccor[i,j] <- ourout$estimate
    lungmethataccortest[i,j] <- ourout$p.value
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



# brown

brownMsigmatrix <- matrix(as.numeric(as.matrix(brownMsigmat)),nrow = dim(brownMsigmat)[1],ncol = dim(brownMsigmat)[2])
rownames(brownMsigmatrix) <- rownames(brownMsigmat)
colnames(brownMsigmatrix) <- colnames(brownMsigmat)

browncols <- intersect(colnames(brownatacnorm),colnames(brownMsigmatrix))
brownatacnorm <- brownatacnorm[,browncols]
brownMsigmatrix <- brownMsigmatrix[,browncols]

brownatacsignorm <- brownatacnorm[brownatacsig,]

brownatacsigfemalecontrolavg <- apply(brownatacsignorm[,intersect(colnames(brownatacsignorm),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% "control",]))],1,mean)
brownatacsigmalecontrolavg <- apply(brownatacsignorm[,intersect(colnames(brownatacsignorm),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% "control",]))],1,mean)

brownatacfemaletrainingl2fcbysubject <- brownatacsignorm[,intersect(colnames(brownatacsignorm),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(brownatacfemaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(brownatacfemaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(brownatacfemaletrainingl2fcbysubject)[j]
    brownatacfemaletrainingl2fcbysubject[i,j] <- brownatacsignorm[i,ourcolname]-brownatacsigfemalecontrolavg[i]
  }
}
brownatacmaletrainingl2fcbysubject <- brownatacsignorm[,intersect(colnames(brownatacsignorm),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(brownatacmaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(brownatacmaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(brownatacmaletrainingl2fcbysubject)[j]
    brownatacmaletrainingl2fcbysubject[i,j] <- brownatacsignorm[i,ourcolname]-brownatacsigmalecontrolavg[i]
  }
}
brownatactrainingl2fcbysubject <- cbind(brownatacfemaletrainingl2fcbysubject,brownatacmaletrainingl2fcbysubject)

brownmethsigfemalecontrolavg <- apply(brownMsigmatrix[,intersect(colnames(brownMsigmatrix),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% "control",]))],1,mean)
brownmethsigmalecontrolavg <- apply(brownMsigmatrix[,intersect(colnames(brownMsigmatrix),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% "control",]))],1,mean)

brownmethfemaletrainingl2fcbysubject <- brownMsigmatrix[,intersect(colnames(brownMsigmatrix),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(brownmethfemaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(brownmethfemaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(brownmethfemaletrainingl2fcbysubject)[j]
    brownmethfemaletrainingl2fcbysubject[i,j] <- brownMsigmatrix[i,ourcolname]-brownmethsigfemalecontrolavg[i]
  }
}
brownmethmaletrainingl2fcbysubject <- brownMsigmatrix[,intersect(colnames(brownMsigmatrix),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(brownmethmaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(brownmethmaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(brownmethmaletrainingl2fcbysubject)[j]
    brownmethmaletrainingl2fcbysubject[i,j] <- brownMsigmatrix[i,ourcolname]-brownmethsigmalecontrolavg[i]
  }
}
brownmethtrainingl2fcbysubject <- cbind(brownmethfemaletrainingl2fcbysubject,brownmethmaletrainingl2fcbysubject)


brownmethataccor <- matrix(0L,nrow = dim(brownmethatacdistance)[1],ncol = dim(brownmethatacdistance)[2])
rownames(brownmethataccor) <- rownames(brownmethatacdistance)
colnames(brownmethataccor) <- colnames(brownmethatacdistance)
brownmethataccortest <- matrix(0L,nrow = dim(brownmethatacdistance)[1],ncol = dim(brownmethatacdistance)[2])
rownames(brownmethataccortest) <- rownames(brownmethatacdistance)
colnames(brownmethataccortest) <- colnames(brownmethatacdistance)
for(i in 1:dim(brownmethataccor)[1]){
  if(i%%20 == 0){
    print(i)
  }
  for(j in 1:dim(brownmethataccor)[2]){
    ourmeth <- rownames(brownmethatacdistance)[i]
    ourpeak <- colnames(brownmethatacdistance)[j]
    ourout <- cor.test(brownmethtrainingl2fcbysubject[ourmeth,],brownatactrainingl2fcbysubject[ourpeak,])
    brownmethataccor[i,j] <- ourout$estimate
    brownmethataccortest[i,j] <- ourout$p.value
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



# white

whiteMsigmatrix <- matrix(as.numeric(as.matrix(whiteMsigmat)),nrow = dim(whiteMsigmat)[1],ncol = dim(whiteMsigmat)[2])
rownames(whiteMsigmatrix) <- rownames(whiteMsigmat)
colnames(whiteMsigmatrix) <- colnames(whiteMsigmat)

whitecols <- intersect(colnames(whiteatacnorm),colnames(whiteMsigmatrix))
whiteatacnorm <- whiteatacnorm[,whitecols]
whiteMsigmatrix <- whiteMsigmatrix[,whitecols]

whiteatacsignorm <- whiteatacnorm[whiteatacsig,]

whiteatacsigfemalecontrolavg <- apply(whiteatacsignorm[,intersect(colnames(whiteatacsignorm),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% "control",]))],1,mean)
whiteatacsigmalecontrolavg <- apply(whiteatacsignorm[,intersect(colnames(whiteatacsignorm),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% "control",]))],1,mean)

whiteatacfemaletrainingl2fcbysubject <- whiteatacsignorm[,intersect(colnames(whiteatacsignorm),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(whiteatacfemaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(whiteatacfemaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(whiteatacfemaletrainingl2fcbysubject)[j]
    whiteatacfemaletrainingl2fcbysubject[i,j] <- whiteatacsignorm[i,ourcolname]-whiteatacsigfemalecontrolavg[i]
  }
}
whiteatacmaletrainingl2fcbysubject <- whiteatacsignorm[,intersect(colnames(whiteatacsignorm),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(whiteatacmaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(whiteatacmaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(whiteatacmaletrainingl2fcbysubject)[j]
    whiteatacmaletrainingl2fcbysubject[i,j] <- whiteatacsignorm[i,ourcolname]-whiteatacsigmalecontrolavg[i]
  }
}
whiteatactrainingl2fcbysubject <- cbind(whiteatacfemaletrainingl2fcbysubject,whiteatacmaletrainingl2fcbysubject)

whitemethsigfemalecontrolavg <- apply(whiteMsigmatrix[,intersect(colnames(whiteMsigmatrix),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% "control",]))],1,mean)
whitemethsigmalecontrolavg <- apply(whiteMsigmatrix[,intersect(colnames(whiteMsigmatrix),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% "control",]))],1,mean)

whitemethfemaletrainingl2fcbysubject <- whiteMsigmatrix[,intersect(colnames(whiteMsigmatrix),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(whitemethfemaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(whitemethfemaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(whitemethfemaletrainingl2fcbysubject)[j]
    whitemethfemaletrainingl2fcbysubject[i,j] <- whiteMsigmatrix[i,ourcolname]-whitemethsigfemalecontrolavg[i]
  }
}
whitemethmaletrainingl2fcbysubject <- whiteMsigmatrix[,intersect(colnames(whiteMsigmatrix),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(whitemethmaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(whitemethmaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(whitemethmaletrainingl2fcbysubject)[j]
    whitemethmaletrainingl2fcbysubject[i,j] <- whiteMsigmatrix[i,ourcolname]-whitemethsigmalecontrolavg[i]
  }
}
whitemethtrainingl2fcbysubject <- cbind(whitemethfemaletrainingl2fcbysubject,whitemethmaletrainingl2fcbysubject)


whitemethataccor <- matrix(0L,nrow = dim(whitemethatacdistance)[1],ncol = dim(whitemethatacdistance)[2])
rownames(whitemethataccor) <- rownames(whitemethatacdistance)
colnames(whitemethataccor) <- colnames(whitemethatacdistance)
whitemethataccortest <- matrix(0L,nrow = dim(whitemethatacdistance)[1],ncol = dim(whitemethatacdistance)[2])
rownames(whitemethataccortest) <- rownames(whitemethatacdistance)
colnames(whitemethataccortest) <- colnames(whitemethatacdistance)
for(i in 1:dim(whitemethataccor)[1]){
  if(i%%20 == 0){
    print(i)
  }
  for(j in 1:dim(whitemethataccor)[2]){
    ourmeth <- rownames(whitemethatacdistance)[i]
    ourpeak <- colnames(whitemethatacdistance)[j]
    ourout <- cor.test(whitemethtrainingl2fcbysubject[ourmeth,],whiteatactrainingl2fcbysubject[ourpeak,])
    whitemethataccor[i,j] <- ourout$estimate
    whitemethataccortest[i,j] <- ourout$p.value
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



# Figure 3H
pdf(file = "Figure 3H_112624.pdf",width = 9,height = 6)
ggplot(totalmethatacdistvscordf[abs(totalmethatacdistvscordf$Distance) > 0 & abs(totalmethatacdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to Nearest DAR") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_density2d_filled() + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 18),legend.title = element_text(size = 18),legend.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + geom_point(data = totalmethatacdistvscordf[abs(totalmethatacdistvscordf$Distance) > 0 & abs(totalmethatacdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation,color=Tissue),size = 2) + scale_color_manual(values = ann_cols$Tissue)
dev.off()

png(file = "Figure 3H_112624.png",width = 9,height = 6,units = "in",res = 600)
ggplot(totalmethatacdistvscordf[abs(totalmethatacdistvscordf$Distance) > 0 & abs(totalmethatacdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to Nearest DAR") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_density2d_filled() + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 18),legend.title = element_text(size = 18),legend.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + geom_point(data = totalmethatacdistvscordf[abs(totalmethatacdistvscordf$Distance) > 0 & abs(totalmethatacdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation,color=Tissue),size = 2) + scale_color_manual(values = ann_cols$Tissue)
dev.off()


# Supplemental Figure S16

pdf(file = "Supplemental Figure S16A_112624.pdf",width = 7,height = 7)
ggplot(gastromethatacdistvscordf[abs(gastromethatacdistvscordf$Distance) > 0 & abs(gastromethatacdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to Nearest DAR") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["SKM-GN"]) + geom_density2d(color = ann_cols$Tissue["SKM-GN"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["SKM-GN"])
dev.off()

pdf(file = "Supplemental Figure S16B_112624.pdf",width = 7,height = 7)
ggplot(heartmethatacdistvscordf[abs(heartmethatacdistvscordf$Distance) > 0 & abs(heartmethatacdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to Nearest DAR") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["HEART"]) + geom_density2d(color = ann_cols$Tissue["HEART"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["HEART"])
dev.off()

pdf(file = "Supplemental Figure S16C_112624.pdf",width = 7,height = 7)
ggplot(kidneymethatacdistvscordf[abs(kidneymethatacdistvscordf$Distance) > 0 & abs(kidneymethatacdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to Nearest DAR") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["KIDNEY"]) + geom_density2d(color = ann_cols$Tissue["KIDNEY"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["KIDNEY"])
dev.off()

pdf(file = "Supplemental Figure S16D_112624.pdf",width = 7,height = 7)
ggplot(livermethatacdistvscordf[abs(livermethatacdistvscordf$Distance) > 0 & abs(livermethatacdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to Nearest DAR") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["LIVER"]) + geom_density2d(color = ann_cols$Tissue["LIVER"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["LIVER"])
dev.off()

pdf(file = "Supplemental Figure S16E_112624.pdf",width = 7,height = 7)
ggplot(lungmethatacdistvscordf[abs(lungmethatacdistvscordf$Distance) > 0 & abs(lungmethatacdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to Nearest DAR") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["LUNG"]) + geom_density2d(color = ann_cols$Tissue["LUNG"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["LUNG"])
dev.off()

pdf(file = "Supplemental Figure S16F_112624.pdf",width = 7,height = 7)
ggplot(brownmethatacdistvscordf[abs(brownmethatacdistvscordf$Distance) > 0 & abs(brownmethatacdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to Nearest DAR") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["BAT"]) + geom_density2d(color = ann_cols$Tissue["BAT"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["BAT"])
dev.off()

png(file = "Supplemental Figure S16A_112624.png",width = 7,height = 7,units = "in",res = 600)
ggplot(gastromethatacdistvscordf[abs(gastromethatacdistvscordf$Distance) > 0 & abs(gastromethatacdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to Nearest DAR") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["SKM-GN"]) + geom_density2d(color = ann_cols$Tissue["SKM-GN"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["SKM-GN"])
dev.off()

png(file = "Supplemental Figure S16B_112624.png",width = 7,height = 7,units = "in",res = 600)
ggplot(heartmethatacdistvscordf[abs(heartmethatacdistvscordf$Distance) > 0 & abs(heartmethatacdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to Nearest DAR") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["HEART"]) + geom_density2d(color = ann_cols$Tissue["HEART"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["HEART"])
dev.off()

png(file = "Supplemental Figure S16C_112624.png",width = 7,height = 7,units = "in",res = 600)
ggplot(kidneymethatacdistvscordf[abs(kidneymethatacdistvscordf$Distance) > 0 & abs(kidneymethatacdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to Nearest DAR") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["KIDNEY"]) + geom_density2d(color = ann_cols$Tissue["KIDNEY"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["KIDNEY"])
dev.off()

png(file = "Supplemental Figure S16D_112624.png",width = 7,height = 7,units = "in",res = 600)
ggplot(livermethatacdistvscordf[abs(livermethatacdistvscordf$Distance) > 0 & abs(livermethatacdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to Nearest DAR") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["LIVER"]) + geom_density2d(color = ann_cols$Tissue["LIVER"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["LIVER"])
dev.off()

png(file = "Supplemental Figure S16E_112624.png",width = 7,height = 7,units = "in",res = 600)
ggplot(lungmethatacdistvscordf[abs(lungmethatacdistvscordf$Distance) > 0 & abs(lungmethatacdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to Nearest DAR") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["LUNG"]) + geom_density2d(color = ann_cols$Tissue["LUNG"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["LUNG"])
dev.off()

png(file = "Supplemental Figure S16F_112624.png",width = 7,height = 7,units = "in",res = 600)
ggplot(brownmethatacdistvscordf[abs(brownmethatacdistvscordf$Distance) > 0 & abs(brownmethatacdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to Nearest DAR") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["BAT"]) + geom_density2d(color = ann_cols$Tissue["BAT"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["BAT"])
dev.off()


save.image("Figure3E_3H_S16_112624.RData")