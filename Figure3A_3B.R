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

# Figure 3A
pdf(file = "Figure 3A.pdf",width = 7.5,height = 4)
ggplot(atacbarplotdf, aes(fill=Group, y=Count, x=Tissue)) + 
  geom_bar(position="stack", stat="identity") + scale_fill_manual(values = c("blue3","darkorange2")) + theme_classic() + ggtitle("DARs")
dev.off()

png(file = "Figure 3A.png",width = 7.5,height = 4,units = "in",res = 600)
ggplot(atacbarplotdf, aes(fill=Group, y=Count, x=Tissue)) + 
  geom_bar(position="stack", stat="identity") + scale_fill_manual(values = c("blue3","darkorange2")) + theme_classic() + ggtitle("DARs")
dev.off()

trimmedmethanno <- readRDS("trimmedmethanno.RDS")

gastromethsigdf <- data.frame(row.names = gastromethsig,
                              "meth" = gastromethsig,
                              "RNA" = trimmedmethanno[gastromethsig,"ensembl_gene"])

heartmethsigdf <- data.frame(row.names = heartmethsig,
                             "meth" = heartmethsig,
                             "RNA" = trimmedmethanno[heartmethsig,"ensembl_gene"])

hippomethsigdf <- data.frame(row.names = hippomethsig,
                             "meth" = hippomethsig,
                             "RNA" = trimmedmethanno[hippomethsig,"ensembl_gene"])

kidneymethsigdf <- data.frame(row.names = kidneymethsig,
                              "meth" = kidneymethsig,
                              "RNA" = trimmedmethanno[kidneymethsig,"ensembl_gene"])

livermethsigdf <- data.frame(row.names = livermethsig,
                             "meth" = livermethsig,
                             "RNA" = trimmedmethanno[livermethsig,"ensembl_gene"])

lungmethsigdf <- data.frame(row.names = lungmethsig,
                            "meth" = lungmethsig,
                            "RNA" = trimmedmethanno[lungmethsig,"ensembl_gene"])

brownmethsigdf <- data.frame(row.names = brownmethsig,
                             "meth" = brownmethsig,
                             "RNA" = trimmedmethanno[brownmethsig,"ensembl_gene"])

whitemethsigdf <- data.frame(row.names = whitemethsig,
                             "meth" = whitemethsig,
                             "RNA" = trimmedmethanno[whitemethsig,"ensembl_gene"])

methbarplotdf <- data.frame("Tissue" = c(rep("HIPPOC",2),
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
methbarplotdf[1,"Count"] <- sum(hippomethsigdf$RNA %in% hippornasig)
methbarplotdf[2,"Count"] <- length(hippomethsigdf$RNA) - sum(hippomethsigdf$RNA %in% hippornasig)

methbarplotdf[3,"Count"] <- sum(gastromethsigdf$RNA %in% gastrornasig)
methbarplotdf[4,"Count"] <- length(gastromethsigdf$RNA) - sum(gastromethsigdf$RNA %in% gastrornasig)

methbarplotdf[5,"Count"] <- sum(heartmethsigdf$RNA %in% heartrnasig)
methbarplotdf[6,"Count"] <- length(heartmethsigdf$RNA) - sum(heartmethsigdf$RNA %in% heartrnasig)

methbarplotdf[7,"Count"] <- sum(kidneymethsigdf$RNA %in% kidneyrnasig)
methbarplotdf[8,"Count"] <- length(kidneymethsigdf$RNA) - sum(kidneymethsigdf$RNA %in% kidneyrnasig)

methbarplotdf[9,"Count"] <- sum(lungmethsigdf$RNA %in% lungrnasig)
methbarplotdf[10,"Count"] <- length(lungmethsigdf$RNA) - sum(lungmethsigdf$RNA %in% lungrnasig)

methbarplotdf[11,"Count"] <- sum(livermethsigdf$RNA %in% liverrnasig)
methbarplotdf[12,"Count"] <- length(livermethsigdf$RNA) - sum(livermethsigdf$RNA %in% liverrnasig)

methbarplotdf[13,"Count"] <- sum(brownmethsigdf$RNA %in% brownrnasig)
methbarplotdf[14,"Count"] <- length(brownmethsigdf$RNA) - sum(brownmethsigdf$RNA %in% brownrnasig)

methbarplotdf[15,"Count"] <- sum(whitemethsigdf$RNA %in% whiternasig)
methbarplotdf[16,"Count"] <- length(whitemethsigdf$RNA) - sum(whitemethsigdf$RNA %in% whiternasig)

# Figure 3B
pdf(file = "Figure 3B.pdf",width = 7.5,height = 4)
ggplot(methbarplotdf, aes(fill=Group, y=Count, x=Tissue)) + 
  geom_bar(position="stack", stat="identity") + scale_fill_manual(values = c("blue3","darkorange2")) + theme_classic() + ggtitle("DMRs")
dev.off()

png(file = "Figure 3B.png",width = 7.5,height = 4,units = "in",res = 600)
ggplot(methbarplotdf, aes(fill=Group, y=Count, x=Tissue)) + 
  geom_bar(position="stack", stat="identity") + scale_fill_manual(values = c("blue3","darkorange2")) + theme_classic() + ggtitle("DMRs")
dev.off()

remove(peakanno)
gc()

save.image("Figure3A_3B.RData")