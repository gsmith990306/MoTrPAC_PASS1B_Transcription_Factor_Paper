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


pdf(file = "Figure 3C_112624.pdf",width = 7,height = 6)
ggplot(totalpeakdistancedf, aes(x = Distance)) +
  geom_histogram(aes(color = Tissue,fill = Tissue),
                 position = "stack", bins = 30) + theme_classic() +
  scale_color_manual(values = c("#8c5220","#f28b2f","#bf7534","#7553a7","#da6c75","#04bf8a","#088c03","#214da6")) + 
  scale_fill_manual(values = c("#8c5220","#f28b2f","#bf7534","#7553a7","#da6c75","#04bf8a","#088c03","#214da6"))  + scale_x_log10(breaks=c(1,100,10000,1000000,100000000)) + xlab("Distance to Nearest DEG TSS") + ylab("Count") + theme(axis.title = element_text(size = 15),axis.text = element_text(size = 12),legend.title = element_text(size = 15),legend.text = element_text(size = 12))
dev.off()

png(file = "Figure 3C_112624.png",width = 7,height = 6,units = "in",res = 600)
ggplot(totalpeakdistancedf, aes(x = Distance)) +
  geom_histogram(aes(color = Tissue,fill = Tissue),
                 position = "stack", bins = 30) + theme_classic() +
  scale_color_manual(values = c("#8c5220","#f28b2f","#bf7534","#7553a7","#da6c75","#04bf8a","#088c03","#214da6")) + 
  scale_fill_manual(values = c("#8c5220","#f28b2f","#bf7534","#7553a7","#da6c75","#04bf8a","#088c03","#214da6"))  + scale_x_log10(breaks=c(1,100,10000,1000000,100000000)) + xlab("Distance to Nearest DEG TSS") + ylab("Count") + theme(axis.title = element_text(size = 15),axis.text = element_text(size = 12),legend.title = element_text(size = 15),legend.text = element_text(size = 12))
dev.off()


load("new_rnanormmatrices.RData")
load("new_atacnormmatrices.RData")

pass1bphenodata <- readRDS("pass1bphenodata.rds")

for(i in 1:dim(gastrornanorm)[2]){
  ourviallabel <- colnames(gastrornanorm)[i]
  ourpid <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourviallabel,"pass1bf0001"][1]
  colnames(gastrornanorm)[i] <- ourpid
}

for(i in 1:dim(gastroatacnorm)[2]){
  ourviallabel <- colnames(gastroatacnorm)[i]
  ourpid <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourviallabel,"pass1bf0001"][1]
  colnames(gastroatacnorm)[i] <- ourpid
}
 
for(i in 1:dim(heartrnanorm)[2]){
  ourviallabel <- colnames(heartrnanorm)[i]
  ourpid <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourviallabel,"pass1bf0001"][1]
  colnames(heartrnanorm)[i] <- ourpid
}  
  
for(i in 1:dim(heartatacnorm)[2]){
  ourviallabel <- colnames(heartatacnorm)[i]
  ourpid <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourviallabel,"pass1bf0001"][1]
  colnames(heartatacnorm)[i] <- ourpid
}
  
for(i in 1:dim(hippornanorm)[2]){
  ourviallabel <- colnames(hippornanorm)[i]
  ourpid <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourviallabel,"pass1bf0001"][1]
  colnames(hippornanorm)[i] <- ourpid
}

for(i in 1:dim(hippoatacnorm)[2]){
  ourviallabel <- colnames(hippoatacnorm)[i]
  ourpid <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourviallabel,"pass1bf0001"][1]
  colnames(hippoatacnorm)[i] <- ourpid
}
  
for(i in 1:dim(kidneyrnanorm)[2]){
  ourviallabel <- colnames(kidneyrnanorm)[i]
  ourpid <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourviallabel,"pass1bf0001"][1]
  colnames(kidneyrnanorm)[i] <- ourpid
}
 
for(i in 1:dim(kidneyatacnorm)[2]){ 
  ourviallabel <- colnames(kidneyatacnorm)[i]
  ourpid <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourviallabel,"pass1bf0001"][1]
  colnames(kidneyatacnorm)[i] <- ourpid
}
 
for(i in 1:dim(liverrnanorm)[2]){ 
  ourviallabel <- colnames(liverrnanorm)[i]
  ourpid <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourviallabel,"pass1bf0001"][1]
  colnames(liverrnanorm)[i] <- ourpid
}

for(i in 1:dim(liveratacnorm)[2]){ 
  ourviallabel <- colnames(liveratacnorm)[i]
  ourpid <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourviallabel,"pass1bf0001"][1]
  colnames(liveratacnorm)[i] <- ourpid
}
  
for(i in 1:dim(lungrnanorm)[2]){
  ourviallabel <- colnames(lungrnanorm)[i]
  ourpid <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourviallabel,"pass1bf0001"][1]
  colnames(lungrnanorm)[i] <- ourpid
}

for(i in 1:dim(lungatacnorm)[2]){
  ourviallabel <- colnames(lungatacnorm)[i]
  ourpid <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourviallabel,"pass1bf0001"][1]
  colnames(lungatacnorm)[i] <- ourpid
}
 
for(i in 1:dim(brownrnanorm)[2]){ 
  ourviallabel <- colnames(brownrnanorm)[i]
  ourpid <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourviallabel,"pass1bf0001"][1]
  colnames(brownrnanorm)[i] <- ourpid
}

for(i in 1:dim(brownatacnorm)[2]){
  ourviallabel <- colnames(brownatacnorm)[i]
  ourpid <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourviallabel,"pass1bf0001"][1]
  colnames(brownatacnorm)[i] <- ourpid
}
 
for(i in 1:dim(whiternanorm)[2]){ 
  ourviallabel <- colnames(whiternanorm)[i]
  ourpid <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourviallabel,"pass1bf0001"][1]
  colnames(whiternanorm)[i] <- ourpid
}

for(i in 1:dim(whiteatacnorm)[2]){ 
  ourviallabel <- colnames(whiteatacnorm)[i]
  ourpid <- pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourviallabel,"pass1bf0001"][1]
  colnames(whiteatacnorm)[i] <- ourpid
}

gastrocols <- intersect(colnames(gastrornanorm),colnames(gastroatacnorm))
gastrornanorm <- gastrornanorm[,gastrocols]
gastroatacnorm <- gastroatacnorm[,gastrocols]

heartcols <- intersect(colnames(heartrnanorm),colnames(heartatacnorm))
heartrnanorm <- heartrnanorm[,heartcols]
heartatacnorm <- heartatacnorm[,heartcols]

hippocols <- intersect(colnames(hippornanorm),colnames(hippoatacnorm))
hippornanorm <- hippornanorm[,hippocols]
hippoatacnorm <- hippoatacnorm[,hippocols]

kidneycols <- intersect(colnames(kidneyrnanorm),colnames(kidneyatacnorm))
kidneyrnanorm <- kidneyrnanorm[,kidneycols]
kidneyatacnorm <- kidneyatacnorm[,kidneycols]

livercols <- intersect(colnames(liverrnanorm),colnames(liveratacnorm))
liverrnanorm <- liverrnanorm[,livercols]
liveratacnorm <- liveratacnorm[,livercols]

lungcols <- intersect(colnames(lungrnanorm),colnames(lungatacnorm))
lungrnanorm <- lungrnanorm[,lungcols]
lungatacnorm <- lungatacnorm[,lungcols]

browncols <- intersect(colnames(brownrnanorm),colnames(brownatacnorm))
brownrnanorm <- brownrnanorm[,browncols]
brownatacnorm <- brownatacnorm[,browncols]

whitecols <- intersect(colnames(whiternanorm),colnames(whiteatacnorm))
whiternanorm <- whiternanorm[,whitecols]
whiteatacnorm <- whiteatacnorm[,whitecols]

gastrornasignorm <- as.matrix(gastrornanorm[fingastrornasig,])
gastroatacsignorm <- gastroatacnorm[gastroatacsig,]
heartrnasignorm <- as.matrix(heartrnanorm[finheartrnasig,])
heartatacsignorm <- heartatacnorm[heartatacsig,]
hippornasignorm <- as.matrix(hippornanorm[finhippornasig,])
hippoatacsignorm <- hippoatacnorm[hippoatacsig,]
kidneyrnasignorm <- as.matrix(kidneyrnanorm[finkidneyrnasig,])
kidneyatacsignorm <- kidneyatacnorm[kidneyatacsig,]
liverrnasignorm <- as.matrix(liverrnanorm[finliverrnasig,])
liveratacsignorm <- liveratacnorm[liveratacsig,]
lungrnasignorm <- as.matrix(lungrnanorm[finlungrnasig,])
lungatacsignorm <- lungatacnorm[lungatacsig,]
brownrnasignorm <- as.matrix(brownrnanorm[finbrownrnasig,])
brownatacsignorm <- brownatacnorm[brownatacsig,]
whiternasignorm <- as.matrix(whiternanorm[finwhiternasig,])
whiteatacsignorm <- whiteatacnorm[whiteatacsig,]

# Now we load metadata so we can calculate average male and female control values
# to use as a comparison for l2fc calculations per sample

pass1bphenodata <- readRDS("pass1bphenodata.rds")
pidmeta <- data.frame(row.names = colnames(gastrornanorm),
                      "pid" = colnames(gastrornanorm))
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
gastrornasigfemalecontrolavg <- apply(gastrornasignorm[,intersect(colnames(gastrornasignorm),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% "control",]))],1,mean)
gastrornasigmalecontrolavg <- apply(gastrornasignorm[,intersect(colnames(gastrornasignorm),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% "control",]))],1,mean)

gastrornafemaletrainingl2fcbysubject <- gastrornasignorm[,intersect(colnames(gastrornasignorm),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(gastrornafemaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(gastrornafemaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(gastrornafemaletrainingl2fcbysubject)[j]
    gastrornafemaletrainingl2fcbysubject[i,j] <- gastrornasignorm[i,ourcolname]-gastrornasigfemalecontrolavg[i]
  }
}
gastrornamaletrainingl2fcbysubject <- gastrornasignorm[,intersect(colnames(gastrornasignorm),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(gastrornamaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(gastrornamaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(gastrornamaletrainingl2fcbysubject)[j]
    gastrornamaletrainingl2fcbysubject[i,j] <- gastrornasignorm[i,ourcolname]-gastrornasigmalecontrolavg[i]
  }
}
gastrornatrainingl2fcbysubject <- cbind(gastrornafemaletrainingl2fcbysubject,gastrornamaletrainingl2fcbysubject)

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


gastropeakcor <- matrix(0L,nrow = dim(gastropeakdistance)[1],ncol = dim(gastropeakdistance)[2])
rownames(gastropeakcor) <- rownames(gastropeakdistance)
colnames(gastropeakcor) <- colnames(gastropeakdistance)
gastropeakcortest <- matrix(0L,nrow = dim(gastropeakdistance)[1],ncol = dim(gastropeakdistance)[2])
rownames(gastropeakcortest) <- rownames(gastropeakdistance)
colnames(gastropeakcortest) <- colnames(gastropeakdistance)
for(i in 1:dim(gastropeakcor)[1]){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:dim(gastropeakcor)[2]){
    ourpeak <- rownames(gastropeakdistance)[i]
    ourgene <- colnames(gastropeakdistance)[j]
    ourout <- cor.test(gastroatactrainingl2fcbysubject[ourpeak,],gastrornatrainingl2fcbysubject[ourgene,])
    gastropeakcor[i,j] <- ourout$estimate
    gastropeakcortest[i,j] <- ourout$p.value
  }
}



# heart
heartrnasigfemalecontrolavg <- apply(heartrnasignorm[,intersect(colnames(heartrnasignorm),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% "control",]))],1,mean)
heartrnasigmalecontrolavg <- apply(heartrnasignorm[,intersect(colnames(heartrnasignorm),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% "control",]))],1,mean)

heartrnafemaletrainingl2fcbysubject <- heartrnasignorm[,intersect(colnames(heartrnasignorm),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(heartrnafemaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(heartrnafemaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(heartrnafemaletrainingl2fcbysubject)[j]
    heartrnafemaletrainingl2fcbysubject[i,j] <- heartrnasignorm[i,ourcolname]-heartrnasigfemalecontrolavg[i]
  }
}
heartrnamaletrainingl2fcbysubject <- heartrnasignorm[,intersect(colnames(heartrnasignorm),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(heartrnamaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(heartrnamaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(heartrnamaletrainingl2fcbysubject)[j]
    heartrnamaletrainingl2fcbysubject[i,j] <- heartrnasignorm[i,ourcolname]-heartrnasigmalecontrolavg[i]
  }
}
heartrnatrainingl2fcbysubject <- cbind(heartrnafemaletrainingl2fcbysubject,heartrnamaletrainingl2fcbysubject)

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


heartpeakcor <- matrix(0L,nrow = dim(heartpeakdistance)[1],ncol = dim(heartpeakdistance)[2])
rownames(heartpeakcor) <- rownames(heartpeakdistance)
colnames(heartpeakcor) <- colnames(heartpeakdistance)
heartpeakcortest <- matrix(0L,nrow = dim(heartpeakdistance)[1],ncol = dim(heartpeakdistance)[2])
rownames(heartpeakcortest) <- rownames(heartpeakdistance)
colnames(heartpeakcortest) <- colnames(heartpeakdistance)
for(i in 1:dim(heartpeakcor)[1]){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:dim(heartpeakcor)[2]){
    ourpeak <- rownames(heartpeakdistance)[i]
    ourgene <- colnames(heartpeakdistance)[j]
    ourout <- cor.test(heartatactrainingl2fcbysubject[ourpeak,],heartrnatrainingl2fcbysubject[ourgene,])
    heartpeakcor[i,j] <- ourout$estimate
    heartpeakcortest[i,j] <- ourout$p.value
  }
}

# hippo
hippornasigfemalecontrolavg <- apply(hippornasignorm[,intersect(colnames(hippornasignorm),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% "control",]))],1,mean)
hippornasigmalecontrolavg <- apply(hippornasignorm[,intersect(colnames(hippornasignorm),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% "control",]))],1,mean)

hippornafemaletrainingl2fcbysubject <- hippornasignorm[,intersect(colnames(hippornasignorm),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(hippornafemaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(hippornafemaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(hippornafemaletrainingl2fcbysubject)[j]
    hippornafemaletrainingl2fcbysubject[i,j] <- hippornasignorm[i,ourcolname]-hippornasigfemalecontrolavg[i]
  }
}
hippornamaletrainingl2fcbysubject <- hippornasignorm[,intersect(colnames(hippornasignorm),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(hippornamaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(hippornamaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(hippornamaletrainingl2fcbysubject)[j]
    hippornamaletrainingl2fcbysubject[i,j] <- hippornasignorm[i,ourcolname]-hippornasigmalecontrolavg[i]
  }
}
hippornatrainingl2fcbysubject <- cbind(hippornafemaletrainingl2fcbysubject,hippornamaletrainingl2fcbysubject)

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


hippopeakcor <- matrix(0L,nrow = dim(hippopeakdistance)[1],ncol = dim(hippopeakdistance)[2])
rownames(hippopeakcor) <- rownames(hippopeakdistance)
colnames(hippopeakcor) <- colnames(hippopeakdistance)
hippopeakcortest <- matrix(0L,nrow = dim(hippopeakdistance)[1],ncol = dim(hippopeakdistance)[2])
rownames(hippopeakcortest) <- rownames(hippopeakdistance)
colnames(hippopeakcortest) <- colnames(hippopeakdistance)
for(i in 1:dim(hippopeakcor)[1]){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:dim(hippopeakcor)[2]){
    ourpeak <- rownames(hippopeakdistance)[i]
    ourgene <- colnames(hippopeakdistance)[j]
    ourout <- cor.test(hippoatactrainingl2fcbysubject[ourpeak,],hippornatrainingl2fcbysubject[ourgene,])
    hippopeakcor[i,j] <- ourout$estimate
    hippopeakcortest[i,j] <- ourout$p.value
  }
}


# kidney
kidneyrnasigfemalecontrolavg <- apply(kidneyrnasignorm[,intersect(colnames(kidneyrnasignorm),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% "control",]))],1,mean)
kidneyrnasigmalecontrolavg <- apply(kidneyrnasignorm[,intersect(colnames(kidneyrnasignorm),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% "control",]))],1,mean)

kidneyrnafemaletrainingl2fcbysubject <- kidneyrnasignorm[,intersect(colnames(kidneyrnasignorm),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(kidneyrnafemaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(kidneyrnafemaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(kidneyrnafemaletrainingl2fcbysubject)[j]
    kidneyrnafemaletrainingl2fcbysubject[i,j] <- kidneyrnasignorm[i,ourcolname]-kidneyrnasigfemalecontrolavg[i]
  }
}
kidneyrnamaletrainingl2fcbysubject <- kidneyrnasignorm[,intersect(colnames(kidneyrnasignorm),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(kidneyrnamaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(kidneyrnamaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(kidneyrnamaletrainingl2fcbysubject)[j]
    kidneyrnamaletrainingl2fcbysubject[i,j] <- kidneyrnasignorm[i,ourcolname]-kidneyrnasigmalecontrolavg[i]
  }
}
kidneyrnatrainingl2fcbysubject <- cbind(kidneyrnafemaletrainingl2fcbysubject,kidneyrnamaletrainingl2fcbysubject)

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


kidneypeakcor <- matrix(0L,nrow = dim(kidneypeakdistance)[1],ncol = dim(kidneypeakdistance)[2])
rownames(kidneypeakcor) <- rownames(kidneypeakdistance)
colnames(kidneypeakcor) <- colnames(kidneypeakdistance)
kidneypeakcortest <- matrix(0L,nrow = dim(kidneypeakdistance)[1],ncol = dim(kidneypeakdistance)[2])
rownames(kidneypeakcortest) <- rownames(kidneypeakdistance)
colnames(kidneypeakcortest) <- colnames(kidneypeakdistance)
for(i in 1:dim(kidneypeakcor)[1]){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:dim(kidneypeakcor)[2]){
    ourpeak <- rownames(kidneypeakdistance)[i]
    ourgene <- colnames(kidneypeakdistance)[j]
    ourout <- cor.test(kidneyatactrainingl2fcbysubject[ourpeak,],kidneyrnatrainingl2fcbysubject[ourgene,])
    kidneypeakcor[i,j] <- ourout$estimate
    kidneypeakcortest[i,j] <- ourout$p.value
  }
}


# liver
liverrnasigfemalecontrolavg <- apply(liverrnasignorm[,intersect(colnames(liverrnasignorm),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% "control",]))],1,mean)
liverrnasigmalecontrolavg <- apply(liverrnasignorm[,intersect(colnames(liverrnasignorm),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% "control",]))],1,mean)

liverrnafemaletrainingl2fcbysubject <- liverrnasignorm[,intersect(colnames(liverrnasignorm),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(liverrnafemaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(liverrnafemaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(liverrnafemaletrainingl2fcbysubject)[j]
    liverrnafemaletrainingl2fcbysubject[i,j] <- liverrnasignorm[i,ourcolname]-liverrnasigfemalecontrolavg[i]
  }
}
liverrnamaletrainingl2fcbysubject <- liverrnasignorm[,intersect(colnames(liverrnasignorm),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(liverrnamaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(liverrnamaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(liverrnamaletrainingl2fcbysubject)[j]
    liverrnamaletrainingl2fcbysubject[i,j] <- liverrnasignorm[i,ourcolname]-liverrnasigmalecontrolavg[i]
  }
}
liverrnatrainingl2fcbysubject <- cbind(liverrnafemaletrainingl2fcbysubject,liverrnamaletrainingl2fcbysubject)

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


liverpeakcor <- matrix(0L,nrow = dim(liverpeakdistance)[1],ncol = dim(liverpeakdistance)[2])
rownames(liverpeakcor) <- rownames(liverpeakdistance)
colnames(liverpeakcor) <- colnames(liverpeakdistance)
liverpeakcortest <- matrix(0L,nrow = dim(liverpeakdistance)[1],ncol = dim(liverpeakdistance)[2])
rownames(liverpeakcortest) <- rownames(liverpeakdistance)
colnames(liverpeakcortest) <- colnames(liverpeakdistance)
for(i in 1:dim(liverpeakcor)[1]){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:dim(liverpeakcor)[2]){
    ourpeak <- rownames(liverpeakdistance)[i]
    ourgene <- colnames(liverpeakdistance)[j]
    ourout <- cor.test(liveratactrainingl2fcbysubject[ourpeak,],liverrnatrainingl2fcbysubject[ourgene,])
    liverpeakcor[i,j] <- ourout$estimate
    liverpeakcortest[i,j] <- ourout$p.value
  }
}


# lung
lungrnasigfemalecontrolavg <- apply(lungrnasignorm[,intersect(colnames(lungrnasignorm),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% "control",]))],1,mean)
lungrnasigmalecontrolavg <- apply(lungrnasignorm[,intersect(colnames(lungrnasignorm),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% "control",]))],1,mean)

lungrnafemaletrainingl2fcbysubject <- lungrnasignorm[,intersect(colnames(lungrnasignorm),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(lungrnafemaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(lungrnafemaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(lungrnafemaletrainingl2fcbysubject)[j]
    lungrnafemaletrainingl2fcbysubject[i,j] <- lungrnasignorm[i,ourcolname]-lungrnasigfemalecontrolavg[i]
  }
}
lungrnamaletrainingl2fcbysubject <- lungrnasignorm[,intersect(colnames(lungrnasignorm),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(lungrnamaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(lungrnamaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(lungrnamaletrainingl2fcbysubject)[j]
    lungrnamaletrainingl2fcbysubject[i,j] <- lungrnasignorm[i,ourcolname]-lungrnasigmalecontrolavg[i]
  }
}
lungrnatrainingl2fcbysubject <- cbind(lungrnafemaletrainingl2fcbysubject,lungrnamaletrainingl2fcbysubject)

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


lungpeakcor <- matrix(0L,nrow = dim(lungpeakdistance)[1],ncol = dim(lungpeakdistance)[2])
rownames(lungpeakcor) <- rownames(lungpeakdistance)
colnames(lungpeakcor) <- colnames(lungpeakdistance)
lungpeakcortest <- matrix(0L,nrow = dim(lungpeakdistance)[1],ncol = dim(lungpeakdistance)[2])
rownames(lungpeakcortest) <- rownames(lungpeakdistance)
colnames(lungpeakcortest) <- colnames(lungpeakdistance)
for(i in 1:dim(lungpeakcor)[1]){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:dim(lungpeakcor)[2]){
    ourpeak <- rownames(lungpeakdistance)[i]
    ourgene <- colnames(lungpeakdistance)[j]
    ourout <- cor.test(lungatactrainingl2fcbysubject[ourpeak,],lungrnatrainingl2fcbysubject[ourgene,])
    lungpeakcor[i,j] <- ourout$estimate
    lungpeakcortest[i,j] <- ourout$p.value
  }
}


# brown
brownrnasigfemalecontrolavg <- apply(brownrnasignorm[,intersect(colnames(brownrnasignorm),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% "control",]))],1,mean)
brownrnasigmalecontrolavg <- apply(brownrnasignorm[,intersect(colnames(brownrnasignorm),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% "control",]))],1,mean)

brownrnafemaletrainingl2fcbysubject <- brownrnasignorm[,intersect(colnames(brownrnasignorm),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(brownrnafemaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(brownrnafemaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(brownrnafemaletrainingl2fcbysubject)[j]
    brownrnafemaletrainingl2fcbysubject[i,j] <- brownrnasignorm[i,ourcolname]-brownrnasigfemalecontrolavg[i]
  }
}
brownrnamaletrainingl2fcbysubject <- brownrnasignorm[,intersect(colnames(brownrnasignorm),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(brownrnamaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(brownrnamaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(brownrnamaletrainingl2fcbysubject)[j]
    brownrnamaletrainingl2fcbysubject[i,j] <- brownrnasignorm[i,ourcolname]-brownrnasigmalecontrolavg[i]
  }
}
brownrnatrainingl2fcbysubject <- cbind(brownrnafemaletrainingl2fcbysubject,brownrnamaletrainingl2fcbysubject)

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


brownpeakcor <- matrix(0L,nrow = dim(brownpeakdistance)[1],ncol = dim(brownpeakdistance)[2])
rownames(brownpeakcor) <- rownames(brownpeakdistance)
colnames(brownpeakcor) <- colnames(brownpeakdistance)
brownpeakcortest <- matrix(0L,nrow = dim(brownpeakdistance)[1],ncol = dim(brownpeakdistance)[2])
rownames(brownpeakcortest) <- rownames(brownpeakdistance)
colnames(brownpeakcortest) <- colnames(brownpeakdistance)
for(i in 1:dim(brownpeakcor)[1]){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:dim(brownpeakcor)[2]){
    ourpeak <- rownames(brownpeakdistance)[i]
    ourgene <- colnames(brownpeakdistance)[j]
    ourout <- cor.test(brownatactrainingl2fcbysubject[ourpeak,],brownrnatrainingl2fcbysubject[ourgene,])
    brownpeakcor[i,j] <- ourout$estimate
    brownpeakcortest[i,j] <- ourout$p.value
  }
}


# white
whiternasigfemalecontrolavg <- apply(whiternasignorm[,intersect(colnames(whiternasignorm),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% "control",]))],1,mean)
whiternasigmalecontrolavg <- apply(whiternasignorm[,intersect(colnames(whiternasignorm),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% "control",]))],1,mean)

whiternafemaletrainingl2fcbysubject <- whiternasignorm[,intersect(colnames(whiternasignorm),rownames(pidmeta[pidmeta$sex %in% "Female" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(whiternafemaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(whiternafemaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(whiternafemaletrainingl2fcbysubject)[j]
    whiternafemaletrainingl2fcbysubject[i,j] <- whiternasignorm[i,ourcolname]-whiternasigfemalecontrolavg[i]
  }
}
whiternamaletrainingl2fcbysubject <- whiternasignorm[,intersect(colnames(whiternasignorm),rownames(pidmeta[pidmeta$sex %in% "Male" & pidmeta$group %in% c("1w","2w","4w","8w"),]))]
for(i in 1:dim(whiternamaletrainingl2fcbysubject)[1]){
  for(j in 1:dim(whiternamaletrainingl2fcbysubject)[2]){
    ourcolname <- colnames(whiternamaletrainingl2fcbysubject)[j]
    whiternamaletrainingl2fcbysubject[i,j] <- whiternasignorm[i,ourcolname]-whiternasigmalecontrolavg[i]
  }
}
whiternatrainingl2fcbysubject <- cbind(whiternafemaletrainingl2fcbysubject,whiternamaletrainingl2fcbysubject)

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


whitepeakcor <- matrix(0L,nrow = dim(whitepeakdistance)[1],ncol = dim(whitepeakdistance)[2])
rownames(whitepeakcor) <- rownames(whitepeakdistance)
colnames(whitepeakcor) <- colnames(whitepeakdistance)
whitepeakcortest <- matrix(0L,nrow = dim(whitepeakdistance)[1],ncol = dim(whitepeakdistance)[2])
rownames(whitepeakcortest) <- rownames(whitepeakdistance)
colnames(whitepeakcortest) <- colnames(whitepeakdistance)
for(i in 1:dim(whitepeakcor)[1]){
  if(i%%10 == 0){
    print(i)
  }
  for(j in 1:dim(whitepeakcor)[2]){
    ourpeak <- rownames(whitepeakdistance)[i]
    ourgene <- colnames(whitepeakdistance)[j]
    ourout <- cor.test(whiteatactrainingl2fcbysubject[ourpeak,],whiternatrainingl2fcbysubject[ourgene,])
    whitepeakcor[i,j] <- ourout$estimate
    whitepeakcortest[i,j] <- ourout$p.value
  }
}





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



pdf(file = "Figure 3F_112624.pdf",width = 9,height = 6)
ggplot(totaldistvscordf[abs(totaldistvscordf$Distance) > 0 & abs(totaldistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_density2d_filled() + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 18),legend.title = element_text(size = 18),legend.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + geom_point(data = totaldistvscordf[abs(totaldistvscordf$Distance) > 0 & abs(totaldistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation,color=Tissue),size = 2) + scale_color_manual(values = ann_colsheatmaptrim$Tissue)
dev.off()

png(file = "Figure 3F_112624.png",width = 9,height = 6,units = "in", res = 600)
ggplot(totaldistvscordf[abs(totaldistvscordf$Distance) > 0 & abs(totaldistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_density2d_filled() + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 18),legend.title = element_text(size = 18),legend.text = element_text(size = 15)) + xlim(-500000,500000) + ylim(-1,1) + geom_point(data = totaldistvscordf[abs(totaldistvscordf$Distance) > 0 & abs(totaldistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation,color=Tissue),size = 2) + scale_color_manual(values = ann_colsheatmaptrim$Tissue)
dev.off()

#####
# Supplemental Figure S14
####

pdf(file = "Supplemental Figure S14A_112624.pdf",width = 6,height = 6)
ggplot(gastrodistvscordf[abs(gastrodistvscordf$Distance) > 0 & abs(gastrodistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_colsheatmaptrim$Tissue["SKM-GN"]) + geom_density2d(color = ann_colsheatmaptrim$Tissue["SKM-GN"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 18)) + xlim(-500000,500000) + ylim(-1,1)
dev.off()

pdf(file = "Supplemental Figure S14B_112624.pdf",width = 6,height = 6)
ggplot(heartdistvscordf[abs(heartdistvscordf$Distance) > 0 & abs(heartdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_colsheatmaptrim$Tissue["HEART"]) + geom_density2d(color = ann_colsheatmaptrim$Tissue["HEART"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 18)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_colsheatmaptrim$Tissue["HEART"])
dev.off()

pdf(file = "Supplemental Figure S14C_112624.pdf",width = 6,height = 6)
ggplot(hippodistvscordf[abs(hippodistvscordf$Distance) > 0 & abs(hippodistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_colsheatmaptrim$Tissue["HIPPOC"]) + geom_density2d(color = ann_colsheatmaptrim$Tissue["HIPPOC"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 18)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_colsheatmaptrim$Tissue["HIPPOC"])
dev.off()

pdf(file = "Supplemental Figure S14D_112624.pdf",width = 6,height = 6)
ggplot(kidneydistvscordf[abs(kidneydistvscordf$Distance) > 0 & abs(kidneydistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_colsheatmaptrim$Tissue["KIDNEY"]) + geom_density2d(color = ann_colsheatmaptrim$Tissue["KIDNEY"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 18)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_colsheatmaptrim$Tissue["KIDNEY"])
dev.off()

pdf(file = "Supplemental Figure S14E_112624.pdf",width = 6,height = 6)
ggplot(liverdistvscordf[abs(liverdistvscordf$Distance) > 0 & abs(liverdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_colsheatmaptrim$Tissue["LIVER"]) + geom_density2d(color = ann_colsheatmaptrim$Tissue["LIVER"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 18)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_colsheatmaptrim$Tissue["LIVER"])
dev.off()

pdf(file = "Supplemental Figure S14F_112624.pdf",width = 6,height = 6)
ggplot(lungdistvscordf[abs(lungdistvscordf$Distance) > 0 & abs(lungdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_colsheatmaptrim$Tissue["LUNG"]) + geom_density2d(color = ann_colsheatmaptrim$Tissue["LUNG"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 18)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_colsheatmaptrim$Tissue["LUNG"])
dev.off()

pdf(file = "Supplemental Figure S14G_112624.pdf",width = 6,height = 6)
ggplot(browndistvscordf[abs(browndistvscordf$Distance) > 0 & abs(browndistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_colsheatmaptrim$Tissue["BAT"]) + geom_density2d(color = ann_colsheatmaptrim$Tissue["BAT"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 18)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_colsheatmaptrim$Tissue["BAT"])
dev.off()

pdf(file = "Supplemental Figure S14H_112624.pdf",width = 6,height = 6)
ggplot(whitedistvscordf[abs(whitedistvscordf$Distance) > 0 & abs(whitedistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_colsheatmaptrim$Tissue["WAT-SC"]) + geom_density2d(color = ann_colsheatmaptrim$Tissue["WAT-SC"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 18)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_colsheatmaptrim$Tissue["WAT-SC"])
dev.off()

png(file = "Supplemental Figure S14A_112624.png",width = 6,height = 6, units = "in", res = 600)
ggplot(gastrodistvscordf[abs(gastrodistvscordf$Distance) > 0 & abs(gastrodistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_colsheatmaptrim$Tissue["SKM-GN"]) + geom_density2d(color = ann_colsheatmaptrim$Tissue["SKM-GN"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 18)) + xlim(-500000,500000) + ylim(-1,1)
dev.off()

png(file = "Supplemental Figure S14B_112624.png",width = 6,height = 6, units = "in", res = 600)
ggplot(heartdistvscordf[abs(heartdistvscordf$Distance) > 0 & abs(heartdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_colsheatmaptrim$Tissue["HEART"]) + geom_density2d(color = ann_colsheatmaptrim$Tissue["HEART"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 18)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_colsheatmaptrim$Tissue["HEART"])
dev.off()

png(file = "Supplemental Figure S14C_112624.png",width = 6,height = 6, units = "in", res = 600)
ggplot(hippodistvscordf[abs(hippodistvscordf$Distance) > 0 & abs(hippodistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_colsheatmaptrim$Tissue["HIPPOC"]) + geom_density2d(color = ann_colsheatmaptrim$Tissue["HIPPOC"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 18)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_colsheatmaptrim$Tissue["HIPPOC"])
dev.off()

png(file = "Supplemental Figure S14D_112624.png",width = 6,height = 6, units = "in", res = 600)
ggplot(kidneydistvscordf[abs(kidneydistvscordf$Distance) > 0 & abs(kidneydistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_colsheatmaptrim$Tissue["KIDNEY"]) + geom_density2d(color = ann_colsheatmaptrim$Tissue["KIDNEY"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 18)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_colsheatmaptrim$Tissue["KIDNEY"])
dev.off()

png(file = "Supplemental Figure S14E_112624.png",width = 6,height = 6, units = "in", res = 600)
ggplot(liverdistvscordf[abs(liverdistvscordf$Distance) > 0 & abs(liverdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_colsheatmaptrim$Tissue["LIVER"]) + geom_density2d(color = ann_colsheatmaptrim$Tissue["LIVER"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 18)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_colsheatmaptrim$Tissue["LIVER"])
dev.off()

png(file = "Supplemental Figure S14F_112624.png",width = 6,height = 6, units = "in", res = 600)
ggplot(lungdistvscordf[abs(lungdistvscordf$Distance) > 0 & abs(lungdistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_colsheatmaptrim$Tissue["LUNG"]) + geom_density2d(color = ann_colsheatmaptrim$Tissue["LUNG"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 18)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_colsheatmaptrim$Tissue["LUNG"])
dev.off()

png(file = "Supplemental Figure S14G_112624.png",width = 6,height = 6, units = "in", res = 600)
ggplot(browndistvscordf[abs(browndistvscordf$Distance) > 0 & abs(browndistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_colsheatmaptrim$Tissue["BAT"]) + geom_density2d(color = ann_colsheatmaptrim$Tissue["BAT"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 18)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_colsheatmaptrim$Tissue["BAT"])
dev.off()

png(file = "Supplemental Figure S14H_112624.png",width = 6,height = 6, units = "in", res = 600)
ggplot(whitedistvscordf[abs(whitedistvscordf$Distance) > 0 & abs(whitedistvscordf$Distance) < 500000,],aes(x=Distance,y=Correlation)) + theme_classic() + xlab("Distance to TSS") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(color = ann_colsheatmaptrim$Tissue["WAT-SC"]) + geom_density2d(color = ann_colsheatmaptrim$Tissue["WAT-SC"]) + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 18)) + xlim(-500000,500000) + ylim(-1,1) + scale_color_manual(values = ann_colsheatmaptrim$Tissue["WAT-SC"])
dev.off()

save.image(file = "Figure3C_3F_S14_112624.RData")
