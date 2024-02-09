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

phenomeasuredata <- readRDS("phenomeasuredata.rds")
pass1bphenodata <- readRDS("pass1bphenodata.rds")
load("atacnormmatrices.RData")

#####
# Supplemental Figure S26
####

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
gastrophenomeasuredata <- phenomeasuredata[intersect(colnames(gastroatacnorm),rownames(phenomeasuredata)),]
trimmedphenodf <- trimphenodf[trimphenodf$pid %in% gastrophenomeasuredata$pid,]

my_comparisons <- list( c("Control", "4w"), c("4w", "8w"), c("Control", "8w") )

# Supplemental Figure S26A
pdf(file = "Supplemental Figure S26A.pdf",height = 6,width = 7)
ggboxplot(trimmedphenodf,x = "group",y = "wgt_gain_after_train",facet.by = "sex",palette = "jco",add = "jitter",size = 1) + ggtitle("Body Weight",) + theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5)) + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "t.test",size = 8) + ylim(-12,65) + xlab("Experimental Group") + ylab("Change over Experiment")
dev.off()
# Supplemental Figure S26B
pdf(file = "Supplemental Figure S26B.pdf",height = 6,width = 7)
ggboxplot(trimmedphenodf,x = "group",y = "pct_body_fat_change",facet.by = "sex",palette = "jco",add = "jitter",size = 1) + ggtitle("Body Fat") + theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5)) + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "t.test",size = 8) + ylim(-6,9) + xlab("Experimental Group") + ylab("Change over Experiment")
dev.off()
# Supplemental Figure S26C
pdf(file = "Supplemental Figure S26C.pdf",height = 6,width = 7)
ggboxplot(trimmedphenodf,x = "group",y = "vo2_max_change",facet.by = "sex",palette = "jco",add = "jitter",size = 1) + ggtitle("VO2 Max") + theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5)) + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "t.test",size = 8) + ylim(-16,33) + xlab("Experimental Group") + ylab("Change over Experiment")
dev.off()
# Supplemental Figure S26D
pdf(file = "Supplemental Figure S26D.pdf",height = 6,width = 7)
ggboxplot(trimmedphenodf,x = "group",y = "pct_body_lean_change",facet.by = "sex",palette = "jco",add = "jitter",size = 1) + ggtitle("Body Lean") + theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5)) + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "t.test",size = 8) + ylim(-6,12) + xlab("Experimental Group") + ylab("Change over Experiment")
dev.off()
# Supplemental Figure S26E
pdf(file = "Supplemental Figure S26E.pdf",height = 6,width = 7)
ggboxplot(trimmedphenodf,x = "group",y = "lactate_change_dueto_train",facet.by = "sex",palette = "jco",add = "jitter",size = 1) + ggtitle("Lactate") + theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5)) + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "t.test",size = 8) + ylim(2.5,17) + xlab("Experimental Group") + ylab("Change over Experiment")
dev.off()
# Supplemental Figure S26F
pdf(file = "Supplemental Figure S26F.pdf",height = 6,width = 7)
ggboxplot(trimmedphenodf,x = "group",y = "pct_body_fluid_change",facet.by = "sex",palette = "jco",add = "jitter",size = 1) + ggtitle("Body Water") + theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5)) + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "t.test",size = 8) + ylim(-1.6,3) + xlab("Experimental Group") + ylab("Change over Experiment")
dev.off()
# Supplemental Figure S26A
png(file = "Supplemental Figure S26A.png",height = 6,width = 7,units = "in",res = 600)
ggboxplot(trimmedphenodf,x = "group",y = "wgt_gain_after_train",facet.by = "sex",palette = "jco",add = "jitter",size = 1) + ggtitle("Body Weight",) + theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5)) + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "t.test",size = 8) + ylim(-12,65) + xlab("Experimental Group") + ylab("Change over Experiment")
dev.off()
# Supplemental Figure S26B
png(file = "Supplemental Figure S26B.png",height = 6,width = 7,units = "in",res = 600)
ggboxplot(trimmedphenodf,x = "group",y = "pct_body_fat_change",facet.by = "sex",palette = "jco",add = "jitter",size = 1) + ggtitle("Body Fat") + theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5)) + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "t.test",size = 8) + ylim(-6,9) + xlab("Experimental Group") + ylab("Change over Experiment")
dev.off()
# Supplemental Figure S26C
png(file = "Supplemental Figure S26C.png",height = 6,width = 7,units = "in",res = 600)
ggboxplot(trimmedphenodf,x = "group",y = "vo2_max_change",facet.by = "sex",palette = "jco",add = "jitter",size = 1) + ggtitle("VO2 Max") + theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5)) + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "t.test",size = 8) + ylim(-16,33) + xlab("Experimental Group") + ylab("Change over Experiment")
dev.off()
# Supplemental Figure S26D
png(file = "Supplemental Figure S26D.png",height = 6,width = 7,units = "in",res = 600)
ggboxplot(trimmedphenodf,x = "group",y = "pct_body_lean_change",facet.by = "sex",palette = "jco",add = "jitter",size = 1) + ggtitle("Body Lean") + theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5)) + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "t.test",size = 8) + ylim(-6,12) + xlab("Experimental Group") + ylab("Change over Experiment")
dev.off()
# Supplemental Figure S26E
png(file = "Supplemental Figure S26E.png",height = 6,width = 7,units = "in",res = 600)
ggboxplot(trimmedphenodf,x = "group",y = "lactate_change_dueto_train",facet.by = "sex",palette = "jco",add = "jitter",size = 1) + ggtitle("Lactate") + theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5)) + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "t.test",size = 8) + ylim(2.5,17) + xlab("Experimental Group") + ylab("Change over Experiment")
dev.off()
# Supplemental Figure S26F
png(file = "Supplemental Figure S26F.png",height = 6,width = 7,units = "in",res = 600)
ggboxplot(trimmedphenodf,x = "group",y = "pct_body_fluid_change",facet.by = "sex",palette = "jco",add = "jitter",size = 1) + ggtitle("Body Water") + theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5)) + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "t.test",size = 8) + ylim(-1.6,3) + xlab("Experimental Group") + ylab("Change over Experiment")
dev.off()

#####
# Supplemental Figure S27
####

load("rnanormmatrices.RData")
load("omesigdata.RData")

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


# Supplemental Figure S27A
pdf(file = "Supplemental Figure S27A.pdf",width = 3,height = 6)
pheatmap(gastrornasigvsphenocor,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,annotation_col = sigvsphenocormeta,show_colnames = F,annotation_colors = ann_cols_sigphenocor,annotation_names_col = F,fontsize = 15,annotation_legend = F)
dev.off()
# Supplemental Figure S27B
pdf(file = "Supplemental Figure S27B.pdf",width = 3,height = 6)
pheatmap(heartrnasigvsphenocor,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,annotation_col = sigvsphenocormeta,show_colnames = F,annotation_colors = ann_cols_sigphenocor,annotation_names_col = F,fontsize = 15,annotation_legend = F)
dev.off()
# Supplemental Figure S27C
pdf(file = "Supplemental Figure S27C.pdf",width = 3,height = 6)
pheatmap(hippornasigvsphenocor,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,annotation_col = sigvsphenocormeta,show_colnames = F,annotation_colors = ann_cols_sigphenocor,annotation_names_col = F,fontsize = 15,annotation_legend = F)
dev.off()
# Supplemental Figure S27D
pdf(file = "Supplemental Figure S27D.pdf",width = 3,height = 6)
pheatmap(kidneyrnasigvsphenocor,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,annotation_col = sigvsphenocormeta,show_colnames = F,annotation_colors = ann_cols_sigphenocor,annotation_names_col = F,fontsize = 15,annotation_legend = F)
dev.off()
# Supplemental Figure S27E
pdf(file = "Supplemental Figure S27E.pdf",width = 3,height = 6)
pheatmap(liverrnasigvsphenocor,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,annotation_col = sigvsphenocormeta,show_colnames = F,annotation_colors = ann_cols_sigphenocor,annotation_names_col = F,fontsize = 15,annotation_legend = F)
dev.off()
# Supplemental Figure S27F
pdf(file = "Supplemental Figure S27F.pdf",width = 3,height = 6)
pheatmap(lungrnasigvsphenocor,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,annotation_col = sigvsphenocormeta,show_colnames = F,annotation_colors = ann_cols_sigphenocor,annotation_names_col = F,fontsize = 15,annotation_legend = F)
dev.off()
# Supplemental Figure S27G
pdf(file = "Supplemental Figure S27G.pdf",width = 3,height = 6)
pheatmap(brownrnasigvsphenocor,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,annotation_col = sigvsphenocormeta,show_colnames = F,annotation_colors = ann_cols_sigphenocor,annotation_names_col = F,fontsize = 15,annotation_legend = F)
dev.off()
# Supplemental Figure S27H
pdf(file = "Supplemental Figure S27H.pdf",width = 3,height = 6)
pheatmap(whiternasigvsphenocor,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,annotation_col = sigvsphenocormeta,show_colnames = F,annotation_colors = ann_cols_sigphenocor,annotation_names_col = F,fontsize = 15,annotation_legend = F)
dev.off()

# Supplemental Figure S27A
png(file = "Supplemental Figure S27A.png",width = 3,height = 6,units = "in",res = 600)
pheatmap(gastrornasigvsphenocor,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,annotation_col = sigvsphenocormeta,show_colnames = F,annotation_colors = ann_cols_sigphenocor,annotation_names_col = F,fontsize = 15,annotation_legend = F)
dev.off()
# Supplemental Figure S27B
png(file = "Supplemental Figure S27B.png",width = 3,height = 6,units = "in",res = 600)
pheatmap(heartrnasigvsphenocor,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,annotation_col = sigvsphenocormeta,show_colnames = F,annotation_colors = ann_cols_sigphenocor,annotation_names_col = F,fontsize = 15,annotation_legend = F)
dev.off()
# Supplemental Figure S27C
png(file = "Supplemental Figure S27C.png",width = 3,height = 6,units = "in",res = 600)
pheatmap(hippornasigvsphenocor,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,annotation_col = sigvsphenocormeta,show_colnames = F,annotation_colors = ann_cols_sigphenocor,annotation_names_col = F,fontsize = 15,annotation_legend = F)
dev.off()
# Supplemental Figure S27D
png(file = "Supplemental Figure S27D.png",width = 3,height = 6,units = "in",res = 600)
pheatmap(kidneyrnasigvsphenocor,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,annotation_col = sigvsphenocormeta,show_colnames = F,annotation_colors = ann_cols_sigphenocor,annotation_names_col = F,fontsize = 15,annotation_legend = F)
dev.off()
# Supplemental Figure S27E
png(file = "Supplemental Figure S27E.png",width = 3,height = 6,units = "in",res = 600)
pheatmap(liverrnasigvsphenocor,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,annotation_col = sigvsphenocormeta,show_colnames = F,annotation_colors = ann_cols_sigphenocor,annotation_names_col = F,fontsize = 15,annotation_legend = F)
dev.off()
# Supplemental Figure S27F
png(file = "Supplemental Figure S27F.png",width = 3,height = 6,units = "in",res = 600)
pheatmap(lungrnasigvsphenocor,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,annotation_col = sigvsphenocormeta,show_colnames = F,annotation_colors = ann_cols_sigphenocor,annotation_names_col = F,fontsize = 15,annotation_legend = F)
dev.off()
# Supplemental Figure S27G
png(file = "Supplemental Figure S27G.png",width = 3,height = 6,units = "in",res = 600)
pheatmap(brownrnasigvsphenocor,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,annotation_col = sigvsphenocormeta,show_colnames = F,annotation_colors = ann_cols_sigphenocor,annotation_names_col = F,fontsize = 15,annotation_legend = F)
dev.off()
# Supplemental Figure S27H
png(file = "Supplemental Figure S27H.png",width = 3,height = 6,units = "in",res = 600)
pheatmap(whiternasigvsphenocor,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),show_rownames = F,annotation_col = sigvsphenocormeta,show_colnames = F,annotation_colors = ann_cols_sigphenocor,annotation_names_col = F,fontsize = 15,annotation_legend = F)
dev.off()

#####
# Figure 8
####

phenomatrix <- cbind(trimphenodf$wgt_gain_after_train,trimphenodf$pct_body_fat_change,
                     trimphenodf$pct_body_lean_change,trimphenodf$pct_body_fluid_change,
                     trimphenodf$lactate_change_dueto_train,trimphenodf$vo2_max_change)
rownames(phenomatrix) <- rownames(trimphenodf)
colnames(phenomatrix) <- c("Body Weight","Body Fat","Body Lean","Body Water","Lactate","VO2 Max")

# Figure 8A
pdf(file = "Figure 8A.pdf",height = 5,width = 7)
pheatmap(cor(phenomatrix),angle_col = "315",breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),display_numbers = T,number_color = "black",fontsize = 20,show_colnames = T,show_rownames = T)
dev.off()

png(file = "Figure 8A.png",height = 5,width = 7,units = "in",res = 600)
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

tissuemeta <- data.frame(row.names = c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC"),
                         "Tissue" = c("SKM-GN","HEART","HIPPOC","KIDNEY","LIVER","LUNG","BAT","WAT-SC"))


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

# Figure 8B
pdf(file = "Figure 8B.pdf",width = 9,height = 4)
pheatmap(rnaphenocorrfrac,angle_col = "0",cluster_rows = F,cluster_cols = F,breaks = seq(0,0.25,length.out = 101),color = colorpanel(101,"white","firebrick"),display_numbers = T,annotation_col = tissuemeta,annotation_colors = ann_colsheatmaptrim,show_colnames = F,number_color = "black",fontsize = 15,cellwidth = 40)
dev.off()
png(file = "Figure 8B.png",width = 9,height = 4,units = "in",res = 600)
pheatmap(rnaphenocorrfrac,angle_col = "0",cluster_rows = F,cluster_cols = F,breaks = seq(0,0.25,length.out = 101),color = colorpanel(101,"white","firebrick"),display_numbers = T,annotation_col = tissuemeta,annotation_colors = ann_colsheatmaptrim,show_colnames = F,number_color = "black",fontsize = 15,cellwidth = 40)
dev.off()

#####
# Figure 7C
####

peakanno <- readRDS("peakanno.RDS")
load("activepeakfiles.RData")
allpeakmotifs <- readRDS("allpeakmotifs.RDS")

gastro50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% gastroactivepeaks,]
heart50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% heartactivepeaks,]
hippo50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% hippoactivepeaks,]
kidney50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% kidneyactivepeaks,]
liver50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% liveractivepeaks,]
lung50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% lungactivepeaks,]
brown50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% brownactivepeaks,]
white50peakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% whiteactivepeaks,]

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

# Figure 8C
pdf(file = "Figure 8C.pdf",width = 10,height = 5)
ggplot(data = gastrovo2maxpospromdftrim, aes(x = TFabbrev,y = Frequency,group = Measure,color = Measure)) + geom_line(size = 2) + geom_point(size = 3) + theme_classic() + theme(axis.text.x = element_text(size = 15,angle = 315, vjust = 0,hjust = 0.25),axis.text.y = element_text(size = 15),axis.title = element_text(size = 18),legend.text = element_text(size = 15),legend.title = element_text(size = 18),title = element_text(size = 20,hjust = 0.5)) + xlab("Transcription Factor") + ggtitle("TF Enrich in Pos Corr DEGs with VO2 Max in SKM-GN") + geom_point(data = gastrovo2maxpospromdftrim[gastrovo2maxpospromdftrim$Significance == "Y", ], aes(TFabbrev, MaxFrequency + 0.02), shape = "*", size=10, color="black")
dev.off()
png(file = "Figure 8C.png",width = 10,height = 5,units = "in",res = 600)
ggplot(data = gastrovo2maxpospromdftrim, aes(x = TFabbrev,y = Frequency,group = Measure,color = Measure)) + geom_line(size = 2) + geom_point(size = 3) + theme_classic() + theme(axis.text.x = element_text(size = 15,angle = 315, vjust = 0,hjust = 0.25),axis.text.y = element_text(size = 15),axis.title = element_text(size = 18),legend.text = element_text(size = 15),legend.title = element_text(size = 18),title = element_text(size = 20,hjust = 0.5)) + xlab("Transcription Factor") + ggtitle("TF Enrich in Pos Corr DEGs with VO2 Max in SKM-GN") + geom_point(data = gastrovo2maxpospromdftrim[gastrovo2maxpospromdftrim$Significance == "Y", ], aes(TFabbrev, MaxFrequency + 0.02), shape = "*", size=10, color="black")
dev.off()

# Figure 8E
load("enstosym.RData")
ourgene <- enstosym[enstosym$Symbol %in% "Me3","Ensembl"]
ourdf <- trimgastrophenodf
ourdf$Gene.Expr <- t(gastrornasignorm[ourgene,rownames(trimgastrophenodf)])
pdf(file = "Figure 8E.pdf",width = 7,height = 4)
ggplot(ourdf,aes(x=Gene.Expr,y=vo2_max_change,shape=sex,color=group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$vo2_max_change,ourdf$Gene.Expr),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("VO2 Max Change") + xlab("Gene Expression")
dev.off()
png(file = "Figure 8E.png",width = 7,height = 4,units = "in",res = 600)
ggplot(ourdf,aes(x=Gene.Expr,y=vo2_max_change,shape=sex,color=group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$vo2_max_change,ourdf$Gene.Expr),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("VO2 Max Change") + xlab("Gene Expression")
dev.off()

# Figure 8F
ourgene <- enstosym[enstosym$Symbol %in% "Rora","Ensembl"]
ourdf <- trimgastrophenodf
ourdf$Gene.Expr <- t(gastrornasignorm[ourgene,rownames(trimgastrophenodf)])
pdf(file = "Figure 8F.pdf",width = 7,height = 4)
ggplot(ourdf,aes(x=Gene.Expr,y=vo2_max_change,shape=sex,color=group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$vo2_max_change,ourdf$Gene.Expr),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("VO2 Max Change") + xlab("Gene Expression")
dev.off()
png(file = "Figure 8F.png",width = 7,height = 4,units = "in",res = 600)
ggplot(ourdf,aes(x=Gene.Expr,y=vo2_max_change,shape=sex,color=group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$vo2_max_change,ourdf$Gene.Expr),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("VO2 Max Change") + xlab("Gene Expression")
dev.off()

# Figure 8G
ourgene <- enstosym[enstosym$Symbol %in% "Lgi3","Ensembl"]
ourdf <- trimgastrophenodf
ourdf$Gene.Expr <- t(gastrornasignorm[ourgene,rownames(trimgastrophenodf)])
pdf(file = "Figure 8G.pdf",width = 7,height = 4)
ggplot(ourdf,aes(x=Gene.Expr,y=vo2_max_change,shape=sex,color=group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$vo2_max_change,ourdf$Gene.Expr),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("VO2 Max Change") + xlab("Gene Expression")
dev.off()
png(file = "Figure 8G.png",width = 7,height = 4,units = "in",res = 600)
ggplot(ourdf,aes(x=Gene.Expr,y=vo2_max_change,shape=sex,color=group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$vo2_max_change,ourdf$Gene.Expr),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("VO2 Max Change") + xlab("Gene Expression")
dev.off()

#####
# Figure 8H
####

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

pdf(file = "Figure 8H.pdf",width = 10,height = 5)
ggplot(data = gastrowgtgainpospromdftrim, aes(x = TFabbrev,y = Frequency,group = Measure,color = Measure)) + geom_line(size = 2) + geom_point(size = 3) + theme_classic() + theme(axis.text.x = element_text(size = 15,angle = 315, vjust = 0,hjust = 0.25),axis.text.y = element_text(size = 15),axis.title = element_text(size = 18),legend.text = element_text(size = 15),legend.title = element_text(size = 18),title = element_text(size = 20,hjust = 0.5)) + xlab("Transcription Factor") + ggtitle("TF Enrich in Pos Corr DEGs with Body Weight in SKM-GN") + geom_point(data = gastrowgtgainpospromdftrim[gastrowgtgainpospromdftrim$Significance == "Y", ], aes(TFabbrev, MaxFrequency + 0.02), shape = "*", size=10, color="black")
dev.off()
png(file = "Figure 8H.png",width = 10,height = 5,units = "in",res = 600)
ggplot(data = gastrowgtgainpospromdftrim, aes(x = TFabbrev,y = Frequency,group = Measure,color = Measure)) + geom_line(size = 2) + geom_point(size = 3) + theme_classic() + theme(axis.text.x = element_text(size = 15,angle = 315, vjust = 0,hjust = 0.25),axis.text.y = element_text(size = 15),axis.title = element_text(size = 18),legend.text = element_text(size = 15),legend.title = element_text(size = 18),title = element_text(size = 20,hjust = 0.5)) + xlab("Transcription Factor") + ggtitle("TF Enrich in Pos Corr DEGs with Body Weight in SKM-GN") + geom_point(data = gastrowgtgainpospromdftrim[gastrowgtgainpospromdftrim$Significance == "Y", ], aes(TFabbrev, MaxFrequency + 0.02), shape = "*", size=10, color="black")
dev.off()

# Figure 8J
ourgene <- enstosym[enstosym$Symbol %in% "Chd7","Ensembl"]
ourdf <- trimgastrophenodf
ourdf$Gene.Expr <- t(gastrornasignorm[ourgene,rownames(trimgastrophenodf)])
pdf(file = "Figure 8J.pdf",width = 7,height = 4)
ggplot(ourdf,aes(x=Gene.Expr,y=wgt_gain_after_train,shape=sex,color=group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$wgt_gain_after_train,ourdf$Gene.Expr),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("Body Weight Change") + xlab("Gene Expression")
dev.off()

# Figure 8K
ourgene <- enstosym[enstosym$Symbol %in% "Igf2","Ensembl"]
ourdf <- trimgastrophenodf
ourdf$Gene.Expr <- t(gastrornasignorm[ourgene,rownames(trimgastrophenodf)])
pdf(file = "Figure 8K.pdf",width = 7,height = 4)
ggplot(ourdf,aes(x=Gene.Expr,y=pct_body_fat_change,shape=sex,color=group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$pct_body_fat_change,ourdf$Gene.Expr),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("Body Fat Change") + xlab("Gene Expression")
dev.off()

# Figure 8L
ourgene <- enstosym[enstosym$Symbol %in% "Sall2","Ensembl"]
ourdf <- trimgastrophenodf
ourdf$Gene.Expr <- t(gastrornasignorm[ourgene,rownames(trimgastrophenodf)])
pdf(file = "Figure 8L.pdf",width = 7,height = 4)
ggplot(ourdf,aes(x=Gene.Expr,y=pct_body_fat_change,shape=sex,color=group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$pct_body_fat_change,ourdf$Gene.Expr),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("Body Fat Change") + xlab("Gene Expression")
dev.off()

# Figure 8M
ourgene <- enstosym[enstosym$Symbol %in% "Oas2","Ensembl"]
ourdf <- trimlungphenodf
ourdf$Gene.Expr <- t(lungrnasignorm[ourgene,rownames(trimlungphenodf)])
pdf(file = "Figure 8M.pdf",width = 7,height = 4)
ggplot(ourdf,aes(x=Gene.Expr,y=pct_body_fat_change,shape=sex,color=group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LUNG ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$pct_body_fat_change,ourdf$Gene.Expr),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("Body Fat Change") + xlab("Gene Expression")
dev.off()

# Figure 8N
ourgene <- enstosym[enstosym$Symbol %in% "Nfkb2","Ensembl"]
ourdf <- trimlungphenodf
ourdf$Gene.Expr <- t(lungrnasignorm[ourgene,rownames(trimlungphenodf)])
pdf(file = "Figure 8N.pdf",width = 7,height = 4)
ggplot(ourdf,aes(x=Gene.Expr,y=pct_body_fat_change,shape=sex,color=group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LUNG ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$pct_body_fat_change,ourdf$Gene.Expr),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("Body Fat Change") + xlab("Gene Expression")
dev.off()

# Figure 8O
ourgene <- enstosym[enstosym$Symbol %in% "Fkbp4","Ensembl"]
ourdf <- trimliverphenodf
ourdf$Gene.Expr <- t(liverrnasignorm[ourgene,rownames(trimliverphenodf)])
pdf(file = "Figure 8O.pdf",width = 7,height = 4)
ggplot(ourdf,aes(x=Gene.Expr,y=pct_body_fat_change,shape=sex,color=group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LIVER ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$pct_body_fat_change,ourdf$Gene.Expr),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("Body Fat Change") + xlab("Gene Expression")
dev.off()


# Figure 8J
ourgene <- enstosym[enstosym$Symbol %in% "Chd7","Ensembl"]
ourdf <- trimgastrophenodf
ourdf$Gene.Expr <- t(gastrornasignorm[ourgene,rownames(trimgastrophenodf)])
png(file = "Figure 8J.png",width = 7,height = 4,units = "in",res = 600)
ggplot(ourdf,aes(x=Gene.Expr,y=wgt_gain_after_train,shape=sex,color=group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$wgt_gain_after_train,ourdf$Gene.Expr),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("Body Weight Change") + xlab("Gene Expression")
dev.off()

# Figure 8K
ourgene <- enstosym[enstosym$Symbol %in% "Igf2","Ensembl"]
ourdf <- trimgastrophenodf
ourdf$Gene.Expr <- t(gastrornasignorm[ourgene,rownames(trimgastrophenodf)])
png(file = "Figure 8K.png",width = 7,height = 4,units = "in",res = 600)
ggplot(ourdf,aes(x=Gene.Expr,y=pct_body_fat_change,shape=sex,color=group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$pct_body_fat_change,ourdf$Gene.Expr),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("Body Fat Change") + xlab("Gene Expression")
dev.off()

# Figure 8L
ourgene <- enstosym[enstosym$Symbol %in% "Sall2","Ensembl"]
ourdf <- trimgastrophenodf
ourdf$Gene.Expr <- t(gastrornasignorm[ourgene,rownames(trimgastrophenodf)])
png(file = "Figure 8L.png",width = 7,height = 4,units = "in",res = 600)
ggplot(ourdf,aes(x=Gene.Expr,y=pct_body_fat_change,shape=sex,color=group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$pct_body_fat_change,ourdf$Gene.Expr),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("Body Fat Change") + xlab("Gene Expression")
dev.off()

# Figure 8M
ourgene <- enstosym[enstosym$Symbol %in% "Oas2","Ensembl"]
ourdf <- trimlungphenodf
ourdf$Gene.Expr <- t(lungrnasignorm[ourgene,rownames(trimlungphenodf)])
png(file = "Figure 8M.png",width = 7,height = 4,units = "in",res = 600)
ggplot(ourdf,aes(x=Gene.Expr,y=pct_body_fat_change,shape=sex,color=group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LUNG ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$pct_body_fat_change,ourdf$Gene.Expr),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("Body Fat Change") + xlab("Gene Expression")
dev.off()

# Figure 8N
ourgene <- enstosym[enstosym$Symbol %in% "Nfkb2","Ensembl"]
ourdf <- trimlungphenodf
ourdf$Gene.Expr <- t(lungrnasignorm[ourgene,rownames(trimlungphenodf)])
png(file = "Figure 8N.png",width = 7,height = 4,units = "in",res = 600)
ggplot(ourdf,aes(x=Gene.Expr,y=pct_body_fat_change,shape=sex,color=group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LUNG ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$pct_body_fat_change,ourdf$Gene.Expr),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("Body Fat Change") + xlab("Gene Expression")
dev.off()

# Figure 8O
ourgene <- enstosym[enstosym$Symbol %in% "Fkbp4","Ensembl"]
ourdf <- trimliverphenodf
ourdf$Gene.Expr <- t(liverrnasignorm[ourgene,rownames(trimliverphenodf)])
png(file = "Figure 8O.png",width = 7,height = 4,units = "in",res = 600)
ggplot(ourdf,aes(x=Gene.Expr,y=pct_body_fat_change,shape=sex,color=group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("LIVER ",enstosym[ourgene,"Symbol"],": Pearson Correlation = ",substr(toString(cor(ourdf$pct_body_fat_change,ourdf$Gene.Expr),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("Body Fat Change") + xlab("Gene Expression")
dev.off()

#####
# Supplemental Figure S28
####

tfanno <- readRDS("tfanno.RDS")

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

# Supplemental Figure S28A
pdf(file = "Supplemental Figure S28A.pdf",height = 7,width = 8)
pheatmap(combornatfvsphenocor,cluster_rows = T,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),annotation_row = combornatfvsphenometa,annotation_col = sigvsphenocormeta,annotation_colors = ann_cols_sigphenocor,show_colnames = F,display_numbers = T,number_color = "black",labels_row = gsub("\\(.*","",combotfrnanamelist),fontsize = 15)
dev.off()
png(file = "Supplemental Figure S28A.png",height = 7,width = 8,units = "in",res = 600)
pheatmap(combornatfvsphenocor,cluster_rows = T,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),annotation_row = combornatfvsphenometa,annotation_col = sigvsphenocormeta,annotation_colors = ann_cols_sigphenocor,show_colnames = F,display_numbers = T,number_color = "black",labels_row = gsub("\\(.*","",combotfrnanamelist),fontsize = 15)
dev.off()

#####
# Supplemental Figure S28B
####

load("PASS1B Transcription Factor Paper Data/DEA Analysis/pass1b-06_proteomics-prot-pr_dea.RData")
tfproanno <- readRDS("tfproanno.RDS")

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

# Supplemental Figure S28B
pdf(file = "Supplemental Figure S28B.pdf",height = 7,width = 8)
pheatmap(comboprotfvsphenocor,cluster_rows = T,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),annotation_row = comboprotfvsphenometa,annotation_col = sigvsphenocormeta,annotation_colors = ann_cols_sigphenocor,show_colnames = F,display_numbers = T,number_color = "black",labels_row = gsub("\\(.*","",combotfpronamelist),fontsize = 15)
dev.off()
png(file = "Supplemental Figure S28B.png",height = 7,width = 8,units = "in",res = 600)
pheatmap(comboprotfvsphenocor,cluster_rows = T,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"blue","white","red"),annotation_row = comboprotfvsphenometa,annotation_col = sigvsphenocormeta,annotation_colors = ann_cols_sigphenocor,show_colnames = F,display_numbers = T,number_color = "black",labels_row = gsub("\\(.*","",combotfpronamelist),fontsize = 15)
dev.off()

save.image("Figure8_S26_S27_S28.RData")