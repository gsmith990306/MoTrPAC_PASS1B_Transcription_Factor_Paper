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

load("rnanormmatrices.RData")
load("enstosym.RData")

pass1bphenodata <- readRDS("pass1bphenodata.rds")

adiabdz <- read.csv(file = "PASS1B Transcription Factor Paper Data/Deconvolution Data/adiabdcellz.csv",header = T,row.names = 1)

####
# Cellular Deconvolution Analysis for Figures 1E, 1F
#####

data("GSE20300")
data("GSE20300.grp")
data("cellCounts")
data("Garvan")
data("IRIS")
data("DMAP")

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

adiabdztrim <- adiabdz[intersect(rownames(gastrornatrim),rownames(adiabdz)),]
adiabdzfinal <- adiabdztrim[rowSums(is.na(adiabdztrim)) == 0,]
colnames(adiabdzfinal)[5] <- "Myeloid.Cells"

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

adiabdztrim <- adiabdz[intersect(rownames(heartrnatrim),rownames(adiabdz)),]
adiabdzfinal <- adiabdztrim[rowSums(is.na(adiabdztrim)) == 0,]
colnames(adiabdzfinal)[5] <- "Myeloid.Cells"

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

adiabdztrim <- adiabdz[intersect(rownames(hippornatrim),rownames(adiabdz)),]
adiabdzfinal <- adiabdztrim[rowSums(is.na(adiabdztrim)) == 0,]
colnames(adiabdzfinal)[5] <- "Myeloid.Cells"

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

adiabdztrim <- adiabdz[intersect(rownames(kidneyrnatrim),rownames(adiabdz)),]
adiabdzfinal <- adiabdztrim[rowSums(is.na(adiabdztrim)) == 0,]
colnames(adiabdzfinal)[5] <- "Myeloid.Cells"

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

adiabdztrim <- adiabdz[intersect(rownames(liverrnatrim),rownames(adiabdz)),]
adiabdzfinal <- adiabdztrim[rowSums(is.na(adiabdztrim)) == 0,]
colnames(adiabdzfinal)[5] <- "Myeloid.Cells"

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

adiabdztrim <- adiabdz[intersect(rownames(lungrnatrim),rownames(adiabdz)),]
adiabdzfinal <- adiabdztrim[rowSums(is.na(adiabdztrim)) == 0,]
colnames(adiabdzfinal)[5] <- "Myeloid.Cells"

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

adiabdztrim <- adiabdz[intersect(rownames(brownrnatrim),rownames(adiabdz)),]
adiabdzfinal <- adiabdztrim[rowSums(is.na(adiabdztrim)) == 0,]
colnames(adiabdzfinal)[5] <- "Myeloid.Cells"

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

adiabdztrim <- adiabdz[intersect(rownames(whiternatrim),rownames(adiabdz)),]
adiabdzfinal <- adiabdztrim[rowSums(is.na(adiabdztrim)) == 0,]
colnames(adiabdzfinal)[5] <- "Myeloid.Cells"

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

tissuemeta <- data.frame(row.names = rownames(tissuedeconsexsigmat),
                              "Tissue" = rownames(tissuedeconsexsigmat))

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

# Figure 1E
pdf(file = "Figure 1E.pdf",width = 7,height = 6)
pheatmap(t(-log10(tissuedecontrainsigmat[c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC"),])), angle_col = "0",breaks = seq(0,4,length.out = 9),color = brewer.pal(9,"Reds"),annotation_col = tissuemeta,annotation_colors = ann_cols,cluster_cols = F,show_colnames = F,fontsize = 15)
dev.off()

png(file = "Figure 1E.png",width = 7,height = 6,units = "in",res = 600)
pheatmap(t(-log10(tissuedecontrainsigmat[c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC"),])), angle_col = "0",breaks = seq(0,4,length.out = 9),color = brewer.pal(9,"Reds"),annotation_col = tissuemeta,annotation_colors = ann_cols,cluster_cols = F,show_colnames = F,fontsize = 15)
dev.off()

# Figure 1F
pdf(file = "Figure 1F.pdf",width = 7,height = 6)
pheatmap(t(-log10(tissuedeconsexsigmat[c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC"),])), angle_col = "0",breaks = seq(0,4,length.out = 9),color = brewer.pal(9,"Reds"),annotation_col = tissuemeta,annotation_colors = ann_cols,cluster_cols = F,show_colnames = F,fontsize = 15)
dev.off()

png(file = "Figure 1F.png",width = 7,height = 6,units = "in",res = 600)
pheatmap(t(-log10(tissuedeconsexsigmat[c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC"),])), angle_col = "0",breaks = seq(0,4,length.out = 9),color = brewer.pal(9,"Reds"),annotation_col = tissuemeta,annotation_colors = ann_cols,cluster_cols = F,show_colnames = F,fontsize = 15)
dev.off()

####
# Supplemental Figure S6
#####

pdf(file = "Supplemental Figure S6A.pdf",width = 8,height = 3)
pheatmap(t(scale(browncomboSPVs)),cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_colnames = F,annotation_col = brownrnameta[,c("Group","Sex")],annotation_colors = ann_cols,cellwidth = 6.5)
dev.off()

pdf(file = "Supplemental Figure S6B.pdf",width = 8,height = 3)
pheatmap(t(scale(whitecomboSPVs)),cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_colnames = F,annotation_col = whiternameta[,c("Group","Sex")],annotation_colors = ann_cols,cellwidth = 6.5)
dev.off()

pdf(file = "Supplemental Figure S6C.pdf",width = 8,height = 3)
pheatmap(t(scale(gastrocomboSPVs)),cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_colnames = F,annotation_col = gastrornameta[,c("Group","Sex")],annotation_colors = ann_cols,cellwidth = 6.5)
dev.off()

pdf(file = "Supplemental Figure S6D.pdf",width = 8,height = 3)
pheatmap(t(scale(heartcomboSPVs)),cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_colnames = F,annotation_col = heartrnameta[,c("Group","Sex")],annotation_colors = ann_cols,cellwidth = 6.5)
dev.off()

pdf(file = "Supplemental Figure S6E.pdf",width = 8,height = 3)
pheatmap(t(scale(hippocomboSPVs)),cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_colnames = F,annotation_col = hippornameta[,c("Group","Sex")],annotation_colors = ann_cols,cellwidth = 6.5)
dev.off()

pdf(file = "Supplemental Figure S6F.pdf",width = 8,height = 3)
pheatmap(t(scale(kidneycomboSPVs)),cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_colnames = F,annotation_col = kidneyrnameta[,c("Group","Sex")],annotation_colors = ann_cols,cellwidth = 6.5)
dev.off()

pdf(file = "Supplemental Figure S6G.pdf",width = 8,height = 3)
pheatmap(t(scale(livercomboSPVs)),cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_colnames = F,annotation_col = liverrnameta[,c("Group","Sex")],annotation_colors = ann_cols,cellwidth = 6.5)
dev.off()

pdf(file = "Supplemental Figure S6H.pdf",width = 8,height = 3)
pheatmap(t(scale(lungcomboSPVs)),cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_colnames = F,annotation_col = lungrnameta[,c("Group","Sex")],annotation_colors = ann_cols,cellwidth = 6.5)
dev.off()

png(file = "Supplemental Figure S6A.png",width = 8,height = 3,units = "in", res = 600)
pheatmap(t(scale(browncomboSPVs)),cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_colnames = F,annotation_col = brownrnameta[,c("Group","Sex")],annotation_colors = ann_cols,cellwidth = 6.5)
dev.off()

png(file = "Supplemental Figure S6B.png",width = 8,height = 3,units = "in", res = 600)
pheatmap(t(scale(whitecomboSPVs)),cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_colnames = F,annotation_col = whiternameta[,c("Group","Sex")],annotation_colors = ann_cols,cellwidth = 6.5)
dev.off()

png(file = "Supplemental Figure S6C.png",width = 8,height = 3,units = "in", res = 600)
pheatmap(t(scale(gastrocomboSPVs)),cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_colnames = F,annotation_col = gastrornameta[,c("Group","Sex")],annotation_colors = ann_cols,cellwidth = 6.5)
dev.off()

png(file = "Supplemental Figure S6D.png",width = 8,height = 3,units = "in", res = 600)
pheatmap(t(scale(heartcomboSPVs)),cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_colnames = F,annotation_col = heartrnameta[,c("Group","Sex")],annotation_colors = ann_cols,cellwidth = 6.5)
dev.off()

png(file = "Supplemental Figure S6E.png",width = 8,height = 3,units = "in", res = 600)
pheatmap(t(scale(hippocomboSPVs)),cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_colnames = F,annotation_col = hippornameta[,c("Group","Sex")],annotation_colors = ann_cols,cellwidth = 6.5)
dev.off()

png(file = "Supplemental Figure S6F.png",width = 8,height = 3,units = "in", res = 600)
pheatmap(t(scale(kidneycomboSPVs)),cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_colnames = F,annotation_col = kidneyrnameta[,c("Group","Sex")],annotation_colors = ann_cols,cellwidth = 6.5)
dev.off()

png(file = "Supplemental Figure S6G.png",width = 8,height = 3,units = "in", res = 600)
pheatmap(t(scale(livercomboSPVs)),cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_colnames = F,annotation_col = liverrnameta[,c("Group","Sex")],annotation_colors = ann_cols,cellwidth = 6.5)
dev.off()

png(file = "Supplemental Figure S6H.png",width = 8,height = 3,units = "in", res = 600)
pheatmap(t(scale(lungcomboSPVs)),cluster_cols = F,breaks = seq(-2,2,length.out = 101),color = colorpanel(101,"blue","white","red"),show_colnames = F,annotation_col = lungrnameta[,c("Group","Sex")],annotation_colors = ann_cols,cellwidth = 6.5)
dev.off()

save.image("CellTypeComp_Figure1E_1F_S6.RData")