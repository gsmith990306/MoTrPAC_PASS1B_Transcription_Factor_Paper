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
# Figure 3D
####

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

fingastrornasig <- intersect(gastrornasig,unique(trimmedmethanno$ensembl_gene))

gastromethdistance <- matrix(0L,nrow = length(gastromethsig),ncol = length(fingastrornasig))
rownames(gastromethdistance) <- gastromethsig
colnames(gastromethdistance) <- fingastrornasig

#gastromethsiganno <- modmethanno[gastromethsig,c("chrom","start","end")]
#gastromethsiganno$mid <- (gastromethsiganno$start + gastromethsiganno$end)/2
gastrornasiganno <- data.frame(row.names = fingastrornasig,"chrom" = rep("1",length(fingastrornasig)),"start" = rep(1,length(fingastrornasig)))
for(i in 1:length(fingastrornasig)){
  if(i%%50 == 0){
    print(i)
  }
  ourgene <- fingastrornasig[i]
  ourgastrornapeaks <- trimmedmethanno[trimmedmethanno$ensembl_gene %in% ourgene,c("chrom","geneStart")]
  gastrornasiganno[i,"chrom"] <- ourgastrornapeaks$chrom[1]
  gastrornasiganno[i,"start"] <- ourgastrornapeaks$geneStart[1]
}

for(i in 1:length(gastromethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(fingastrornasig)){
    gastromethdistance[i,j] <- (gastrornasiganno$chrom[j] == fingastromethsiganno$chrom[i])*(abs(gastrornasiganno$start[j] - fingastromethsiganno$sitemid[i]))
  }
}

gastromethdistancedf <- data.frame(row.names = gastromethsig,
                                   "Region" = fingastromethsiganno[gastromethsig,"custom_annotation"],
                                   "Distance" = rep(0,length(gastromethsig)))
for(i in 1:length(gastromethsig)){
  
  ourrow <- gastromethdistance[i,]
  gastromethdistancedf[i,"Distance"] <- min(ourrow[ourrow > 0])
  
}


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

finheartrnasig <- intersect(heartrnasig,unique(trimmedmethanno$ensembl_gene))

heartmethdistance <- matrix(0L,nrow = length(heartmethsig),ncol = length(finheartrnasig))
rownames(heartmethdistance) <- heartmethsig
colnames(heartmethdistance) <- finheartrnasig

#heartmethsiganno <- modmethanno[heartmethsig,c("chrom","start","end")]
#heartmethsiganno$mid <- (heartmethsiganno$start + heartmethsiganno$end)/2
heartrnasiganno <- data.frame(row.names = finheartrnasig,"chrom" = rep("1",length(finheartrnasig)),"start" = rep(1,length(finheartrnasig)))
for(i in 1:length(finheartrnasig)){
  if(i%%50 == 0){
    print(i)
  }
  ourgene <- finheartrnasig[i]
  ourheartrnapeaks <- trimmedmethanno[trimmedmethanno$ensembl_gene %in% ourgene,c("chrom","geneStart")]
  heartrnasiganno[i,"chrom"] <- ourheartrnapeaks$chrom[1]
  heartrnasiganno[i,"start"] <- ourheartrnapeaks$geneStart[1]
}

for(i in 1:length(heartmethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(finheartrnasig)){
    heartmethdistance[i,j] <- (heartrnasiganno$chrom[j] == finheartmethsiganno$chrom[i])*(abs(heartrnasiganno$start[j] - finheartmethsiganno$sitemid[i]))
  }
}

heartmethdistancedf <- data.frame(row.names = heartmethsig,
                                  "Region" = finheartmethsiganno[heartmethsig,"custom_annotation"],
                                  "Distance" = rep(0,length(heartmethsig)))
for(i in 1:length(heartmethsig)){
  
  ourrow <- heartmethdistance[i,]
  heartmethdistancedf[i,"Distance"] <- min(ourrow[ourrow > 0])
  
}


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

finhippornasig <- intersect(hippornasig,unique(trimmedmethanno$ensembl_gene))

hippomethdistance <- matrix(0L,nrow = length(hippomethsig),ncol = length(finhippornasig))
rownames(hippomethdistance) <- hippomethsig
colnames(hippomethdistance) <- finhippornasig

#hippomethsiganno <- modmethanno[hippomethsig,c("chrom","start","end")]
#hippomethsiganno$mid <- (hippomethsiganno$start + hippomethsiganno$end)/2
hippornasiganno <- data.frame(row.names = finhippornasig,"chrom" = rep("1",length(finhippornasig)),"start" = rep(1,length(finhippornasig)))
for(i in 1:length(finhippornasig)){
  if(i%%50 == 0){
    print(i)
  }
  ourgene <- finhippornasig[i]
  ourhippornapeaks <- trimmedmethanno[trimmedmethanno$ensembl_gene %in% ourgene,c("chrom","geneStart")]
  hippornasiganno[i,"chrom"] <- ourhippornapeaks$chrom[1]
  hippornasiganno[i,"start"] <- ourhippornapeaks$geneStart[1]
}

for(i in 1:length(hippomethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(finhippornasig)){
    hippomethdistance[i,j] <- (hippornasiganno$chrom[j] == finhippomethsiganno$chrom[i])*(abs(hippornasiganno$start[j] - finhippomethsiganno$sitemid[i]))
  }
}

hippomethdistancedf <- data.frame(row.names = hippomethsig,
                                  "Region" = finhippomethsiganno[hippomethsig,"custom_annotation"],
                                  "Distance" = rep(0,length(hippomethsig)))
for(i in 1:length(hippomethsig)){
  
  ourrow <- hippomethdistance[i,]
  hippomethdistancedf[i,"Distance"] <- min(ourrow[ourrow > 0])
  
}


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

finkidneyrnasig <- intersect(kidneyrnasig,unique(trimmedmethanno$ensembl_gene))

kidneymethdistance <- matrix(0L,nrow = length(kidneymethsig),ncol = length(finkidneyrnasig))
rownames(kidneymethdistance) <- kidneymethsig
colnames(kidneymethdistance) <- finkidneyrnasig

#kidneymethsiganno <- modmethanno[kidneymethsig,c("chrom","start","end")]
#kidneymethsiganno$mid <- (kidneymethsiganno$start + kidneymethsiganno$end)/2
kidneyrnasiganno <- data.frame(row.names = finkidneyrnasig,"chrom" = rep("1",length(finkidneyrnasig)),"start" = rep(1,length(finkidneyrnasig)))
for(i in 1:length(finkidneyrnasig)){
  if(i%%50 == 0){
    print(i)
  }
  ourgene <- finkidneyrnasig[i]
  ourkidneyrnapeaks <- trimmedmethanno[trimmedmethanno$ensembl_gene %in% ourgene,c("chrom","geneStart")]
  kidneyrnasiganno[i,"chrom"] <- ourkidneyrnapeaks$chrom[1]
  kidneyrnasiganno[i,"start"] <- ourkidneyrnapeaks$geneStart[1]
}

for(i in 1:length(kidneymethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(finkidneyrnasig)){
    kidneymethdistance[i,j] <- (kidneyrnasiganno$chrom[j] == finkidneymethsiganno$chrom[i])*(abs(kidneyrnasiganno$start[j] - finkidneymethsiganno$sitemid[i]))
  }
}

kidneymethdistancedf <- data.frame(row.names = kidneymethsig,
                                   "Region" = finkidneymethsiganno[kidneymethsig,"custom_annotation"],
                                   "Distance" = rep(0,length(kidneymethsig)))
for(i in 1:length(kidneymethsig)){
  
  ourrow <- kidneymethdistance[i,]
  kidneymethdistancedf[i,"Distance"] <- min(ourrow[ourrow > 0])
  
}


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

finliverrnasig <- intersect(liverrnasig,unique(trimmedmethanno$ensembl_gene))

livermethdistance <- matrix(0L,nrow = length(livermethsig),ncol = length(finliverrnasig))
rownames(livermethdistance) <- livermethsig
colnames(livermethdistance) <- finliverrnasig

#livermethsiganno <- modmethanno[livermethsig,c("chrom","start","end")]
#livermethsiganno$mid <- (livermethsiganno$start + livermethsiganno$end)/2
liverrnasiganno <- data.frame(row.names = finliverrnasig,"chrom" = rep("1",length(finliverrnasig)),"start" = rep(1,length(finliverrnasig)))
for(i in 1:length(finliverrnasig)){
  if(i%%50 == 0){
    print(i)
  }
  ourgene <- finliverrnasig[i]
  ourliverrnapeaks <- trimmedmethanno[trimmedmethanno$ensembl_gene %in% ourgene,c("chrom","geneStart")]
  liverrnasiganno[i,"chrom"] <- ourliverrnapeaks$chrom[1]
  liverrnasiganno[i,"start"] <- ourliverrnapeaks$geneStart[1]
}

for(i in 1:length(livermethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(finliverrnasig)){
    livermethdistance[i,j] <- (liverrnasiganno$chrom[j] == finlivermethsiganno$chrom[i])*(abs(liverrnasiganno$start[j] - finlivermethsiganno$sitemid[i]))
  }
}

livermethdistancedf <- data.frame(row.names = livermethsig,
                                  "Region" = finlivermethsiganno[livermethsig,"custom_annotation"],
                                  "Distance" = rep(0,length(livermethsig)))
for(i in 1:length(livermethsig)){
  
  ourrow <- livermethdistance[i,]
  livermethdistancedf[i,"Distance"] <- min(ourrow[ourrow > 0])
  
}


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

finlungrnasig <- intersect(lungrnasig,unique(trimmedmethanno$ensembl_gene))

lungmethdistance <- matrix(0L,nrow = length(lungmethsig),ncol = length(finlungrnasig))
rownames(lungmethdistance) <- lungmethsig
colnames(lungmethdistance) <- finlungrnasig

#lungmethsiganno <- modmethanno[lungmethsig,c("chrom","start","end")]
#lungmethsiganno$mid <- (lungmethsiganno$start + lungmethsiganno$end)/2
lungrnasiganno <- data.frame(row.names = finlungrnasig,"chrom" = rep("1",length(finlungrnasig)),"start" = rep(1,length(finlungrnasig)))
for(i in 1:length(finlungrnasig)){
  if(i%%50 == 0){
    print(i)
  }
  ourgene <- finlungrnasig[i]
  ourlungrnapeaks <- trimmedmethanno[trimmedmethanno$ensembl_gene %in% ourgene,c("chrom","geneStart")]
  lungrnasiganno[i,"chrom"] <- ourlungrnapeaks$chrom[1]
  lungrnasiganno[i,"start"] <- ourlungrnapeaks$geneStart[1]
}

for(i in 1:length(lungmethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(finlungrnasig)){
    lungmethdistance[i,j] <- (lungrnasiganno$chrom[j] == finlungmethsiganno$chrom[i])*(abs(lungrnasiganno$start[j] - finlungmethsiganno$sitemid[i]))
  }
}

lungmethdistancedf <- data.frame(row.names = lungmethsig,
                                 "Region" = finlungmethsiganno[lungmethsig,"custom_annotation"],
                                 "Distance" = rep(0,length(lungmethsig)))
for(i in 1:length(lungmethsig)){
  
  ourrow <- lungmethdistance[i,]
  lungmethdistancedf[i,"Distance"] <- min(ourrow[ourrow > 0])
  
}


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

finbrownrnasig <- intersect(brownrnasig,unique(trimmedmethanno$ensembl_gene))

brownmethdistance <- matrix(0L,nrow = length(brownmethsig),ncol = length(finbrownrnasig))
rownames(brownmethdistance) <- brownmethsig
colnames(brownmethdistance) <- finbrownrnasig

#brownmethsiganno <- modmethanno[brownmethsig,c("chrom","start","end")]
#brownmethsiganno$mid <- (brownmethsiganno$start + brownmethsiganno$end)/2
brownrnasiganno <- data.frame(row.names = finbrownrnasig,"chrom" = rep("1",length(finbrownrnasig)),"start" = rep(1,length(finbrownrnasig)))
for(i in 1:length(finbrownrnasig)){
  if(i%%50 == 0){
    print(i)
  }
  ourgene <- finbrownrnasig[i]
  ourbrownrnapeaks <- trimmedmethanno[trimmedmethanno$ensembl_gene %in% ourgene,c("chrom","geneStart")]
  brownrnasiganno[i,"chrom"] <- ourbrownrnapeaks$chrom[1]
  brownrnasiganno[i,"start"] <- ourbrownrnapeaks$geneStart[1]
}

for(i in 1:length(brownmethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(finbrownrnasig)){
    brownmethdistance[i,j] <- (brownrnasiganno$chrom[j] == finbrownmethsiganno$chrom[i])*(abs(brownrnasiganno$start[j] - finbrownmethsiganno$sitemid[i]))
  }
}

brownmethdistancedf <- data.frame(row.names = brownmethsig,
                                  "Region" = finbrownmethsiganno[brownmethsig,"custom_annotation"],
                                  "Distance" = rep(0,length(brownmethsig)))
for(i in 1:length(brownmethsig)){
  
  ourrow <- brownmethdistance[i,]
  brownmethdistancedf[i,"Distance"] <- min(ourrow[ourrow > 0])
  
}


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

finwhiternasig <- intersect(whiternasig,unique(trimmedmethanno$ensembl_gene))

whitemethdistance <- matrix(0L,nrow = length(whitemethsig),ncol = length(finwhiternasig))
rownames(whitemethdistance) <- whitemethsig
colnames(whitemethdistance) <- finwhiternasig

#whitemethsiganno <- modmethanno[whitemethsig,c("chrom","start","end")]
#whitemethsiganno$mid <- (whitemethsiganno$start + whitemethsiganno$end)/2
whiternasiganno <- data.frame(row.names = finwhiternasig,"chrom" = rep("1",length(finwhiternasig)),"start" = rep(1,length(finwhiternasig)))
for(i in 1:length(finwhiternasig)){
  if(i%%50 == 0){
    print(i)
  }
  ourgene <- finwhiternasig[i]
  ourwhiternapeaks <- trimmedmethanno[trimmedmethanno$ensembl_gene %in% ourgene,c("chrom","geneStart")]
  whiternasiganno[i,"chrom"] <- ourwhiternapeaks$chrom[1]
  whiternasiganno[i,"start"] <- ourwhiternapeaks$geneStart[1]
}

for(i in 1:length(whitemethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(finwhiternasig)){
    whitemethdistance[i,j] <- (whiternasiganno$chrom[j] == finwhitemethsiganno$chrom[i])*(abs(whiternasiganno$start[j] - finwhitemethsiganno$sitemid[i]))
  }
}

whitemethdistancedf <- data.frame(row.names = whitemethsig,
                                  "Region" = finwhitemethsiganno[whitemethsig,"custom_annotation"],
                                  "Distance" = rep(0,length(whitemethsig)))
for(i in 1:length(whitemethsig)){
  
  ourrow <- whitemethdistance[i,]
  whitemethdistancedf[i,"Distance"] <- min(ourrow[ourrow > 0])
  
}


totalmethdistancedf <- data.frame("Region" = c(gastromethdistancedf$Region,
                                               heartmethdistancedf$Region,
                                               hippomethdistancedf$Region,
                                               kidneymethdistancedf$Region,
                                               livermethdistancedf$Region,
                                               lungmethdistancedf$Region,
                                               brownmethdistancedf$Region,
                                               whitemethdistancedf$Region),
                                  "Distance" = c(gastromethdistancedf$Distance,
                                                 heartmethdistancedf$Distance,
                                                 hippomethdistancedf$Distance,
                                                 kidneymethdistancedf$Distance,
                                                 livermethdistancedf$Distance,
                                                 lungmethdistancedf$Distance,
                                                 brownmethdistancedf$Distance,
                                                 whitemethdistancedf$Distance),
                                  "Tissue" = c(rep("SKM-GN",dim(gastromethdistancedf)[1]),
                                               rep("HEART",dim(heartmethdistancedf)[1]),
                                               rep("HIPPOC",dim(hippomethdistancedf)[1]),
                                               rep("KIDNEY",dim(kidneymethdistancedf)[1]),
                                               rep("LIVER",dim(livermethdistancedf)[1]),
                                               rep("LUNG",dim(lungmethdistancedf)[1]),
                                               rep("BAT",dim(brownmethdistancedf)[1]),
                                               rep("WAT-SC",dim(whitemethdistancedf)[1])))
totalmethdistancedf <- totalmethdistancedf[!totalmethdistancedf$Region %in% "Overlaps Gene",]

pdf(file = "Figure 3D.pdf",width = 7,height = 6)
ggplot(totalmethdistancedf, aes(x = Distance)) +
  geom_histogram(aes(color = Tissue,fill = Tissue),
                 position = "stack", bins = 30) + theme_classic() +
  scale_color_manual(values = c("#8c5220","#f28b2f","#bf7534","#7553a7","#da6c75","#04bf8a","#088c03","#214da6")) + 
  scale_fill_manual(values = c("#8c5220","#f28b2f","#bf7534","#7553a7","#da6c75","#04bf8a","#088c03","#214da6"))  + scale_x_log10(breaks=c(1,100,10000,1000000,100000000)) + xlab("Distance to Nearest DEG TSS") + ylab("Count") + theme(axis.title = element_text(size = 15),axis.text = element_text(size = 12),legend.title = element_text(size = 15),legend.text = element_text(size = 12))
dev.off()

png(file = "Figure 3D.png",width = 7,height = 6,units = "in",res = 600)
ggplot(totalmethdistancedf, aes(x = Distance)) +
  geom_histogram(aes(color = Tissue,fill = Tissue),
                 position = "stack", bins = 30) + theme_classic() +
  scale_color_manual(values = c("#8c5220","#f28b2f","#bf7534","#7553a7","#da6c75","#04bf8a","#088c03","#214da6")) + 
  scale_fill_manual(values = c("#8c5220","#f28b2f","#bf7534","#7553a7","#da6c75","#04bf8a","#088c03","#214da6"))  + scale_x_log10(breaks=c(1,100,10000,1000000,100000000)) + xlab("Distance to Nearest DEG TSS") + ylab("Count") + theme(axis.title = element_text(size = 15),axis.text = element_text(size = 12),legend.title = element_text(size = 15),legend.text = element_text(size = 12))
dev.off()


#####
# Figure 3G
####

load("rnal2fcmat.RData")
load("methsigl2fcmat.RData")

# gastro

gastromethcor <- matrix(0L,nrow = dim(gastromethdistance)[1],ncol = dim(gastromethdistance)[2])
rownames(gastromethcor) <- rownames(gastromethdistance)
colnames(gastromethcor) <- colnames(gastromethdistance)
for(i in 1:dim(gastromethcor)[1]){
  if(i%%20 == 0){
    print(i)
  }
  for(j in 1:dim(gastromethcor)[2]){
    ourmeth <- rownames(gastromethdistance)[i]
    ourgene <- colnames(gastromethdistance)[j]
    gastromethcor[i,j] <- cor(gastromethl2fc[ourmeth,],gastrol2fcmat[ourgene,])
  }
}
gastromethcor[is.na(gastromethcor)] <- 0

gastromethnoabsdistance <- matrix(0L,nrow = length(gastromethsig),ncol = length(fingastrornasig))
rownames(gastromethnoabsdistance) <- gastromethsig
colnames(gastromethnoabsdistance) <- fingastrornasig

for(i in 1:length(gastromethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(fingastrornasig)){
    gastromethnoabsdistance[i,j] <- (gastrornasiganno$chrom[j] == fingastromethsiganno$chrom[i])*((gastrornasiganno$start[j] - fingastromethsiganno$sitemid[i]))
  }
}

# heart

heartmethcor <- matrix(0L,nrow = dim(heartmethdistance)[1],ncol = dim(heartmethdistance)[2])
rownames(heartmethcor) <- rownames(heartmethdistance)
colnames(heartmethcor) <- colnames(heartmethdistance)
for(i in 1:dim(heartmethcor)[1]){
  if(i%%20 == 0){
    print(i)
  }
  for(j in 1:dim(heartmethcor)[2]){
    ourmeth <- rownames(heartmethdistance)[i]
    ourgene <- colnames(heartmethdistance)[j]
    heartmethcor[i,j] <- cor(heartmethl2fc[ourmeth,],heartl2fcmat[ourgene,])
  }
}
heartmethcor[is.na(heartmethcor)] <- 0

heartmethnoabsdistance <- matrix(0L,nrow = length(heartmethsig),ncol = length(finheartrnasig))
rownames(heartmethnoabsdistance) <- heartmethsig
colnames(heartmethnoabsdistance) <- finheartrnasig

for(i in 1:length(heartmethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(finheartrnasig)){
    heartmethnoabsdistance[i,j] <- (heartrnasiganno$chrom[j] == finheartmethsiganno$chrom[i])*((heartrnasiganno$start[j] - finheartmethsiganno$sitemid[i]))
  }
}

# hippo

hippomethcor <- matrix(0L,nrow = dim(hippomethdistance)[1],ncol = dim(hippomethdistance)[2])
rownames(hippomethcor) <- rownames(hippomethdistance)
colnames(hippomethcor) <- colnames(hippomethdistance)
for(i in 1:dim(hippomethcor)[1]){
  if(i%%20 == 0){
    print(i)
  }
  for(j in 1:dim(hippomethcor)[2]){
    ourmeth <- rownames(hippomethdistance)[i]
    ourgene <- colnames(hippomethdistance)[j]
    hippomethcor[i,j] <- cor(hippomethl2fc[ourmeth,],hippol2fcmat[ourgene,])
  }
}
hippomethcor[is.na(hippomethcor)] <- 0

hippomethnoabsdistance <- matrix(0L,nrow = length(hippomethsig),ncol = length(finhippornasig))
rownames(hippomethnoabsdistance) <- hippomethsig
colnames(hippomethnoabsdistance) <- finhippornasig

for(i in 1:length(hippomethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(finhippornasig)){
    hippomethnoabsdistance[i,j] <- (hippornasiganno$chrom[j] == finhippomethsiganno$chrom[i])*((hippornasiganno$start[j] - finhippomethsiganno$sitemid[i]))
  }
}

# kidney

kidneymethcor <- matrix(0L,nrow = dim(kidneymethdistance)[1],ncol = dim(kidneymethdistance)[2])
rownames(kidneymethcor) <- rownames(kidneymethdistance)
colnames(kidneymethcor) <- colnames(kidneymethdistance)
for(i in 1:dim(kidneymethcor)[1]){
  if(i%%20 == 0){
    print(i)
  }
  for(j in 1:dim(kidneymethcor)[2]){
    ourmeth <- rownames(kidneymethdistance)[i]
    ourgene <- colnames(kidneymethdistance)[j]
    kidneymethcor[i,j] <- cor(kidneymethl2fc[ourmeth,],kidneyl2fcmat[ourgene,])
  }
}
kidneymethcor[is.na(kidneymethcor)] <- 0

kidneymethnoabsdistance <- matrix(0L,nrow = length(kidneymethsig),ncol = length(finkidneyrnasig))
rownames(kidneymethnoabsdistance) <- kidneymethsig
colnames(kidneymethnoabsdistance) <- finkidneyrnasig

for(i in 1:length(kidneymethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(finkidneyrnasig)){
    kidneymethnoabsdistance[i,j] <- (kidneyrnasiganno$chrom[j] == finkidneymethsiganno$chrom[i])*((kidneyrnasiganno$start[j] - finkidneymethsiganno$sitemid[i]))
  }
}

# liver

livermethcor <- matrix(0L,nrow = dim(livermethdistance)[1],ncol = dim(livermethdistance)[2])
rownames(livermethcor) <- rownames(livermethdistance)
colnames(livermethcor) <- colnames(livermethdistance)
for(i in 1:dim(livermethcor)[1]){
  if(i%%20 == 0){
    print(i)
  }
  for(j in 1:dim(livermethcor)[2]){
    ourmeth <- rownames(livermethdistance)[i]
    ourgene <- colnames(livermethdistance)[j]
    livermethcor[i,j] <- cor(livermethl2fc[ourmeth,],liverl2fcmat[ourgene,])
  }
}
livermethcor[is.na(livermethcor)] <- 0

livermethnoabsdistance <- matrix(0L,nrow = length(livermethsig),ncol = length(finliverrnasig))
rownames(livermethnoabsdistance) <- livermethsig
colnames(livermethnoabsdistance) <- finliverrnasig

for(i in 1:length(livermethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(finliverrnasig)){
    livermethnoabsdistance[i,j] <- (liverrnasiganno$chrom[j] == finlivermethsiganno$chrom[i])*((liverrnasiganno$start[j] - finlivermethsiganno$sitemid[i]))
  }
}

# lung

lungmethcor <- matrix(0L,nrow = dim(lungmethdistance)[1],ncol = dim(lungmethdistance)[2])
rownames(lungmethcor) <- rownames(lungmethdistance)
colnames(lungmethcor) <- colnames(lungmethdistance)
for(i in 1:dim(lungmethcor)[1]){
  if(i%%20 == 0){
    print(i)
  }
  for(j in 1:dim(lungmethcor)[2]){
    ourmeth <- rownames(lungmethdistance)[i]
    ourgene <- colnames(lungmethdistance)[j]
    lungmethcor[i,j] <- cor(lungmethl2fc[ourmeth,],lungl2fcmat[ourgene,])
  }
}
lungmethcor[is.na(lungmethcor)] <- 0

lungmethnoabsdistance <- matrix(0L,nrow = length(lungmethsig),ncol = length(finlungrnasig))
rownames(lungmethnoabsdistance) <- lungmethsig
colnames(lungmethnoabsdistance) <- finlungrnasig

for(i in 1:length(lungmethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(finlungrnasig)){
    lungmethnoabsdistance[i,j] <- (lungrnasiganno$chrom[j] == finlungmethsiganno$chrom[i])*((lungrnasiganno$start[j] - finlungmethsiganno$sitemid[i]))
  }
}

# brown

brownmethcor <- matrix(0L,nrow = dim(brownmethdistance)[1],ncol = dim(brownmethdistance)[2])
rownames(brownmethcor) <- rownames(brownmethdistance)
colnames(brownmethcor) <- colnames(brownmethdistance)
for(i in 1:dim(brownmethcor)[1]){
  if(i%%20 == 0){
    print(i)
  }
  for(j in 1:dim(brownmethcor)[2]){
    ourmeth <- rownames(brownmethdistance)[i]
    ourgene <- colnames(brownmethdistance)[j]
    brownmethcor[i,j] <- cor(brownmethl2fc[ourmeth,],brownl2fcmat[ourgene,])
  }
}
brownmethcor[is.na(brownmethcor)] <- 0

brownmethnoabsdistance <- matrix(0L,nrow = length(brownmethsig),ncol = length(finbrownrnasig))
rownames(brownmethnoabsdistance) <- brownmethsig
colnames(brownmethnoabsdistance) <- finbrownrnasig

for(i in 1:length(brownmethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(finbrownrnasig)){
    brownmethnoabsdistance[i,j] <- (brownrnasiganno$chrom[j] == finbrownmethsiganno$chrom[i])*((brownrnasiganno$start[j] - finbrownmethsiganno$sitemid[i]))
  }
}

# white

whitemethcor <- matrix(0L,nrow = dim(whitemethdistance)[1],ncol = dim(whitemethdistance)[2])
rownames(whitemethcor) <- rownames(whitemethdistance)
colnames(whitemethcor) <- colnames(whitemethdistance)
for(i in 1:dim(whitemethcor)[1]){
  if(i%%20 == 0){
    print(i)
  }
  for(j in 1:dim(whitemethcor)[2]){
    ourmeth <- rownames(whitemethdistance)[i]
    ourgene <- colnames(whitemethdistance)[j]
    whitemethcor[i,j] <- cor(whitemethl2fc[ourmeth,],whitel2fcmat[ourgene,])
  }
}
whitemethcor[is.na(whitemethcor)] <- 0

whitemethnoabsdistance <- matrix(0L,nrow = length(whitemethsig),ncol = length(finwhiternasig))
rownames(whitemethnoabsdistance) <- whitemethsig
colnames(whitemethnoabsdistance) <- finwhiternasig

for(i in 1:length(whitemethsig)){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:length(finwhiternasig)){
    whitemethnoabsdistance[i,j] <- (whiternasiganno$chrom[j] == finwhitemethsiganno$chrom[i])*((whiternasiganno$start[j] - finwhitemethsiganno$sitemid[i]))
  }
}

gastromethdistvscordf <- data.frame("Correlation" = as.vector(gastromethcor),"Distance" = as.vector(gastromethnoabsdistance))
heartmethdistvscordf <- data.frame("Correlation" = as.vector(heartmethcor),"Distance" = as.vector(heartmethnoabsdistance))
hippomethdistvscordf <- data.frame("Correlation" = as.vector(hippomethcor),"Distance" = as.vector(hippomethnoabsdistance))
kidneymethdistvscordf <- data.frame("Correlation" = as.vector(kidneymethcor),"Distance" = as.vector(kidneymethnoabsdistance))
livermethdistvscordf <- data.frame("Correlation" = as.vector(livermethcor),"Distance" = as.vector(livermethnoabsdistance))
lungmethdistvscordf <- data.frame("Correlation" = as.vector(lungmethcor),"Distance" = as.vector(lungmethnoabsdistance))
brownmethdistvscordf <- data.frame("Correlation" = as.vector(brownmethcor),"Distance" = as.vector(brownmethnoabsdistance))
whitemethdistvscordf <- data.frame("Correlation" = as.vector(whitemethcor),"Distance" = as.vector(whitemethnoabsdistance))

totalmethdistvscordf <- rbind(gastromethdistvscordf,heartmethdistvscordf,hippomethdistvscordf,kidneymethdistvscordf,
                              livermethdistvscordf,lungmethdistvscordf,whitemethdistvscordf,brownmethdistvscordf)
totalmethdistvscordf$Tissue <- c(rep("SKM-GN",dim(gastromethdistvscordf)[1]),
                                 rep("HEART",dim(heartmethdistvscordf)[1]),
                                 rep("HIPPOC",dim(hippomethdistvscordf)[1]),
                                 rep("KIDNEY",dim(kidneymethdistvscordf)[1]),
                                 rep("LIVER",dim(livermethdistvscordf)[1]),
                                 rep("LUNG",dim(lungmethdistvscordf)[1]),
                                 rep("WAT-SC",dim(whitemethdistvscordf)[1]),
                                 rep("BAT",dim(brownmethdistvscordf)[1]))

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


# Figure 3G
pdf(file = "Figure 3G.pdf",width = 9,height = 6)
ggplot(totalmethdistvscordf[abs(totalmethdistvscordf$Distance) > 0 & abs(totalmethdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_density2d_filled() + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 18),legend.title = element_text(size = 18),legend.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + geom_point(data = totalmethdistvscordf[abs(totalmethdistvscordf$Distance) > 0 & abs(totalmethdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation,color=Tissue),size = 2) + scale_color_manual(values = ann_cols$Tissue)
dev.off()

png(file = "Figure 3G.png",width = 9,height = 6,units = "in",res = 600)
ggplot(totalmethdistvscordf[abs(totalmethdistvscordf$Distance) > 0 & abs(totalmethdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_density2d_filled() + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 18),legend.title = element_text(size = 18),legend.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + geom_point(data = totalmethdistvscordf[abs(totalmethdistvscordf$Distance) > 0 & abs(totalmethdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation,color=Tissue),size = 2) + scale_color_manual(values = ann_cols$Tissue)
dev.off()

#####
# Supplemental Figure S9
####

pdf(file = "Supplemental Figure S9A.pdf",width = 7,height = 7)
ggplot(gastromethdistvscordf[abs(gastromethdistvscordf$Distance) > 0 & abs(gastromethdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["SKM-GN"]) + geom_density2d(color = ann_cols$Tissue["SKM-GN"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["SKM-GN"])
dev.off()

pdf(file = "Supplemental Figure S9B.pdf",width = 7,height = 7)
ggplot(heartmethdistvscordf[abs(heartmethdistvscordf$Distance) > 0 & abs(heartmethdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["HEART"]) + geom_density2d(color = ann_cols$Tissue["HEART"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["HEART"])
dev.off()

pdf(file = "Supplemental Figure S9C.pdf",width = 7,height = 7)
ggplot(hippomethdistvscordf[abs(hippomethdistvscordf$Distance) > 0 & abs(hippomethdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["HIPPOC"]) + geom_density2d(color = ann_cols$Tissue["HIPPOC"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["HIPPOC"])
dev.off()

pdf(file = "Supplemental Figure S9D.pdf",width = 7,height = 7)
ggplot(kidneymethdistvscordf[abs(kidneymethdistvscordf$Distance) > 0 & abs(kidneymethdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["KIDNEY"]) + geom_density2d(color = ann_cols$Tissue["KIDNEY"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["KIDNEY"])
dev.off()

pdf(file = "Supplemental Figure S9E.pdf",width = 7,height = 7)
ggplot(livermethdistvscordf[abs(livermethdistvscordf$Distance) > 0 & abs(livermethdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["LIVER"]) + geom_density2d(color = ann_cols$Tissue["LIVER"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["LIVER"])
dev.off()

pdf(file = "Supplemental Figure S9F.pdf",width = 7,height = 7)
ggplot(lungmethdistvscordf[abs(lungmethdistvscordf$Distance) > 0 & abs(lungmethdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["LUNG"]) + geom_density2d(color = ann_cols$Tissue["LUNG"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["LUNG"])
dev.off()

pdf(file = "Supplemental Figure S9G.pdf",width = 7,height = 7)
ggplot(brownmethdistvscordf[abs(brownmethdistvscordf$Distance) > 0 & abs(brownmethdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["BAT"]) + geom_density2d(color = ann_cols$Tissue["BAT"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["BAT"])
dev.off()

pdf(file = "Supplemental Figure S9H.pdf",width = 7,height = 7)
ggplot(whitemethdistvscordf[abs(whitemethdistvscordf$Distance) > 0 & abs(whitemethdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["WAT-SC"]) + geom_density2d(color = ann_cols$Tissue["WAT-SC"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["WAT-SC"])
dev.off()


png(file = "Supplemental Figure S9A.png",width = 7,height = 7,units = "in",res = 600)
ggplot(gastromethdistvscordf[abs(gastromethdistvscordf$Distance) > 0 & abs(gastromethdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["SKM-GN"]) + geom_density2d(color = ann_cols$Tissue["SKM-GN"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["SKM-GN"])
dev.off()

png(file = "Supplemental Figure S9B.png",width = 7,height = 7,units = "in",res = 600)
ggplot(heartmethdistvscordf[abs(heartmethdistvscordf$Distance) > 0 & abs(heartmethdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["HEART"]) + geom_density2d(color = ann_cols$Tissue["HEART"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["HEART"])
dev.off()

png(file = "Supplemental Figure S9C.png",width = 7,height = 7,units = "in",res = 600)
ggplot(hippomethdistvscordf[abs(hippomethdistvscordf$Distance) > 0 & abs(hippomethdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["HIPPOC"]) + geom_density2d(color = ann_cols$Tissue["HIPPOC"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["HIPPOC"])
dev.off()

png(file = "Supplemental Figure S9D.png",width = 7,height = 7,units = "in",res = 600)
ggplot(kidneymethdistvscordf[abs(kidneymethdistvscordf$Distance) > 0 & abs(kidneymethdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["KIDNEY"]) + geom_density2d(color = ann_cols$Tissue["KIDNEY"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["KIDNEY"])
dev.off()

png(file = "Supplemental Figure S9E.png",width = 7,height = 7,units = "in",res = 600)
ggplot(livermethdistvscordf[abs(livermethdistvscordf$Distance) > 0 & abs(livermethdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["LIVER"]) + geom_density2d(color = ann_cols$Tissue["LIVER"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["LIVER"])
dev.off()

png(file = "Supplemental Figure S9F.png",width = 7,height = 7,units = "in",res = 600)
ggplot(lungmethdistvscordf[abs(lungmethdistvscordf$Distance) > 0 & abs(lungmethdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["LUNG"]) + geom_density2d(color = ann_cols$Tissue["LUNG"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["LUNG"])
dev.off()

png(file = "Supplemental Figure S9G.png",width = 7,height = 7,units = "in",res = 600)
ggplot(brownmethdistvscordf[abs(brownmethdistvscordf$Distance) > 0 & abs(brownmethdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["BAT"]) + geom_density2d(color = ann_cols$Tissue["BAT"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["BAT"])
dev.off()

png(file = "Supplemental Figure S9H.png",width = 7,height = 7,units = "in",res = 600)
ggplot(whitemethdistvscordf[abs(whitemethdistvscordf$Distance) > 0 & abs(whitemethdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_cols$Tissue["WAT-SC"]) + geom_density2d(color = ann_cols$Tissue["WAT-SC"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_cols$Tissue["WAT-SC"])
dev.off()

save.image("Figure3D_3G_S9.RData")