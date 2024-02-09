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

# Figure 1G
pdf(file = "Figure 1G.pdf", width=8, height=6)
mosaicplot(tab2Perc,las=1,col=genomic_features_col$Region[colnames(tab2Perc)], main="Genomic Distribution of All Peaks",cex.axis = 0.8)
dev.off()

png(file = "Figure 1G.png", width=8, height=6,units = "in",res = 600)
mosaicplot(tab2Perc,las=1,col=genomic_features_col$Region[colnames(tab2Perc)], main="Genomic Distribution of All Peaks")
dev.off()

pcutoff = 0.1
trainingSig = atactraining %>% filter(adj_p_value<pcutoff)

tab = table(trainingSig$tissue_abbreviation,trainingSig$custom_annotation)
tab3 = tab[,myorder[myorder %in% colnames(tab)]]
tab3Perc = t(apply(t(tab3), 2, function(x){x*100/sum(x,na.rm=T)}))

# Figure 1H
pdf(file = "Figure 1H.pdf", width=8, height=6)
mosaicplot(tab3Perc,las=1,col=genomic_features_col$Region[colnames(tab3Perc)], main="Genomic Distribution of Significant Peaks",cex.axis = 0.8)
dev.off()

png(file = "Figure 1H.png", width=8, height=6,units = "in",res = 600)
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

pdf(file = "New Figure 1K.pdf", width=8, height=3)
pheatmap(tab4,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"steelblue","white","firebrick"),cluster_cols = F,cluster_rows = F,border_color = "black",legend = F,angle_col = 0)
dev.off()

png(file = "New Figure 1K.png", width=8, height=3,units = "in",res = 600)
pheatmap(tab4,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"steelblue","white","firebrick"),cluster_cols = F,cluster_rows = F,border_color = "black",legend = F,angle_col = 0)
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
# Supplemental Figure S7A
pdf(file = "Supplemental Figure S7A.pdf", width = 7, height = 6)
mosaicplot(tabPerc, las = 1, col = genomic_features_col$Region[colnames(tabPerc)], main = "Genomic Distribution of Active Genes Percentage")
dev.off()

png(file = "Supplemental Figure S7A.png", width = 7, height = 6,units = "in",res = 600)
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
# Supplemental Figure S7B
pdf(file = "Supplemental Figure S7B.pdf", width = 7, height = 6)
mosaicplot(tab3Perc, las = 1, col = genomic_features_col$Region[colnames(tab3Perc)], main = "Genomic Distribution of DEGaP Percentage")
dev.off()

png(file = "Supplemental Figure S7B.png", width = 7, height = 6,units = "in",res = 600)
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

pdf(file = "Supplemental Figure S7C.pdf", width=8, height=3)
pheatmap(tab4,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"steelblue","white","firebrick"),cluster_cols = F,cluster_rows = F,border_color = "black",legend = F,angle_col = 0)
dev.off()

png(file = "Supplemental Figure S7C.png", width=8, height=3,units = "in",res = 600)
pheatmap(tab4,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"steelblue","white","firebrick"),cluster_cols = F,cluster_rows = F,border_color = "black",legend = F,angle_col = 0)
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
pdf(file = "Figure 1I.pdf", width=8, height=6)
mosaicplot(tab2Perc,las=1,col=genomic_features_col$Region[colnames(tab2Perc)], main="Genomic Distribution of All Methyl Sites",cex.axis = 0.8)
dev.off()

png(file = "Figure 1I.png", width=8, height=6,units = "in",res = 600)
mosaicplot(tab2Perc,las=1,col=genomic_features_col$Region[colnames(tab2Perc)], main="Genomic Distribution of All Methyl Sites")
dev.off()

pcutoff = 0.1
trainingSig = methtraining %>% filter(adj_p_value<pcutoff)

tab = table(trainingSig$tissue_abbreviation,trainingSig$custom_annotation)
tab3 = tab[,myorder[myorder %in% colnames(tab)]]
tab3Perc = t(apply(t(tab3), 2, function(x){x*100/sum(x,na.rm=T)}))

# Figure 1J
pdf(file = "Figure 1J.pdf", width=8, height=6)
mosaicplot(tab3Perc,las=1,col=genomic_features_col$Region[colnames(tab3Perc)], main="Genomic Distribution of Significant Methyl Sites",cex.axis = 0.8)
dev.off()

png(file = "Figure 1J.png", width=8, height=6,units = "in",res = 600)
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
pdf(file = "Figure 1L.pdf", width=8, height=3)
pheatmap(tab4,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"steelblue","white","firebrick"),cluster_cols = F,cluster_rows = F,border_color = "black",legend = F,angle_col = 0)
dev.off()

png(file = "Figure 1L.png", width=8, height=3,units = "in",res = 600)
pheatmap(tab4,breaks = seq(-1,1,length.out = 101),color = colorpanel(101,"steelblue","white","firebrick"),cluster_cols = F,cluster_rows = F,border_color = "black",legend = F,angle_col = 0)
dev.off()

rm(epigen_rrbs)
gc()

save.image("Figure1G_1H_1I_1J_1K_1L_S7.RData")