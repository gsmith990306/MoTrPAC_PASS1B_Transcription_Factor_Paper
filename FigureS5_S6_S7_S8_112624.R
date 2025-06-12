# We are going to apply MultiPower to conduct power analysis on the PASS1B Dataset

library(MultiPower)
library(FDRsampsize)
library(lpSolve)
library(ggfortify)

setwd("C:/Users/gsmit/Documents/Mount Sinai/Sealfon Laboratory/MoTrPAC/PASS1B Raw Data")

myfigdir = "C:/Users/gsmit/Documents/Mount Sinai/Sealfon Laboratory/MoTrPAC/PASS1B Raw Data/Power_Analysis"

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
                 "Group" = c("control" = "grey50",
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


miscolores = c("red2", "darkslategray3", "azure4", "sienna2",
               "dodgerblue4")
omicas = c("RNAseq", "ATACseq", "RRBS", "Prot_PR",  "Prot_PH")
names(miscolores) = omicas

load("new_rnanormmatrices.RData")
load("new_atacnormmatrices.RData")
load("methMmatrices.RData")
load("proprnormmatrices.RData")
load("prophnormmatrices.RData")

pass1bphenodata <- readRDS("pass1bphenodata.rds")


# gastro

gastrostatdata = gastrostatdesign = vector("list")

gastrostatdata$RNAseq = gastrornanorm
gastrostatdata$ATACseq = gastroatacnorm
gastrostatdata$RRBS = matrix(as.numeric(as.matrix(gastroMmat)),nrow = dim(gastroMmat)[1],ncol = dim(gastroMmat)[2])
rownames(gastrostatdata$RRBS) <- rownames(gastroMmat)
colnames(gastrostatdata$RRBS) <- colnames(gastroMmat)
gastrostatdata$Prot_PR = gastroproprnorm
gastrostatdata$Prot_PH = gastroprophnorm

rm(gastrornanorm)
rm(gastroatacnorm)
rm(gastroMmat)
rm(gastroproprnorm)
rm(gastroprophnorm)
gc()

gastrostatdata$RNAseq = gastrostatdata$RNAseq - min(gastrostatdata$RNAseq)
gastrostatdata$ATACseq = gastrostatdata$ATACseq - min(gastrostatdata$ATACseq)
gastrostatdata$RRBS = gastrostatdata$RRBS - min(gastrostatdata$RRBS)
gastrostatdata$Prot_PR = gastrostatdata$Prot_PR - min(gastrostatdata$Prot_PR)
gastrostatdata$Prot_PH = gastrostatdata$Prot_PH - min(gastrostatdata$Prot_PH)

#gastroRNAseqdesign <- data.frame(row.names = colnames(gastrostatdata$RNAseq),
#                                 #"Sex" = rep("Female",dim(gastrostatdata$RNAseq)[2]),
#                                 "Group" = rep("exercise",dim(gastrostatdata$RNAseq)[2]))
gastroRNAseqdesign <- data.frame(row.names = colnames(gastrostatdata$RNAseq),
                                  "Sex" = rep("Female",dim(gastrostatdata$RNAseq)[2]),
                                  "Group" = rep("1w",dim(gastrostatdata$RNAseq)[2]))
for(i in 1:dim(gastroRNAseqdesign)[1]){
    ourid <- rownames(gastroRNAseqdesign)[i]
    gastroRNAseqdesign[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"][1]]
    if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Two-week program"){
      gastroRNAseqdesign[i,"Group"] <- "2w"
    } else{
        if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Four-week program"){
          gastroRNAseqdesign[i,"Group"] <- "4w"
        } else{
            if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Training Group"){
              gastroRNAseqdesign[i,"Group"] <- "8w"
            } else {
                if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Control Group"){
                  gastroRNAseqdesign[i,"Group"] <- "control"
                }
            }
        }
    }
}

gastrostatdesign$RNAseq = gastroRNAseqdesign
gastrostatdesign$RNAseq[gastrostatdesign$RNAseq$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
gastrostatdesign$RNAseq <- data.frame(row.names = rownames(gastroRNAseqdesign),
                                      "Group" = gastrostatdesign$RNAseq$Group)


gastroATACseqdesign <- data.frame(row.names = colnames(gastrostatdata$ATACseq),
                                 "Sex" = rep("Female",dim(gastrostatdata$ATACseq)[2]),
                                 "Group" = rep("1w",dim(gastrostatdata$ATACseq)[2]))
for(i in 1:dim(gastroATACseqdesign)[1]){
  ourid <- rownames(gastroATACseqdesign)[i]
  gastroATACseqdesign[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"][1]]
  if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Two-week program"){
    gastroATACseqdesign[i,"Group"] <- "2w"
  } else{
    if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Four-week program"){
      gastroATACseqdesign[i,"Group"] <- "4w"
    } else{
      if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Training Group"){
        gastroATACseqdesign[i,"Group"] <- "8w"
      } else {
        if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Control Group"){
          gastroATACseqdesign[i,"Group"] <- "control"
        }
      }
    }
  }
}

gastrostatdesign$ATACseq = gastroATACseqdesign
gastrostatdesign$ATACseq[gastrostatdesign$ATACseq$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
gastrostatdesign$ATACseq <- data.frame(row.names = rownames(gastroATACseqdesign),
                                      "Group" = gastrostatdesign$ATACseq$Group)

gastroRRBSdesign <- data.frame(row.names = colnames(gastrostatdata$RRBS),
                               "Sex" = rep("Female",dim(gastrostatdata$RRBS)[2]),
                               "Group" = rep("1w",dim(gastrostatdata$RRBS)[2]))
for(i in 1:dim(gastroRRBSdesign)[1]){
  ourid <- rownames(gastroRRBSdesign)[i]
  gastroRRBSdesign[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0027"][1]]
  if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Two-week program"){
    gastroRRBSdesign[i,"Group"] <- "2w"
  } else{
    if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Four-week program"){
      gastroRRBSdesign[i,"Group"] <- "4w"
    } else{
      if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Training Group"){
        gastroRRBSdesign[i,"Group"] <- "8w"
      } else {
        if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Control Group"){
          gastroRRBSdesign[i,"Group"] <- "control"
        }
      }
    }
  }
}

gastrostatdesign$RRBS = gastroRRBSdesign
gastrostatdesign$RRBS[gastrostatdesign$RRBS$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
gastrostatdesign$RRBS <- data.frame(row.names = rownames(gastroRRBSdesign),
                                    "Group" = gastrostatdesign$RRBS$Group)


gastroProt_PRdesign <- data.frame(row.names = colnames(gastrostatdata$Prot_PR),
                                  "Sex" = rep("Female",dim(gastrostatdata$Prot_PR)[2]),
                                  "Group" = rep("1w",dim(gastrostatdata$Prot_PR)[2]))
for(i in 1:dim(gastroProt_PRdesign)[1]){
  ourid <- rownames(gastroProt_PRdesign)[i]
  gastroProt_PRdesign[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0027"][1]]
  if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Two-week program"){
    gastroProt_PRdesign[i,"Group"] <- "2w"
  } else{
    if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Four-week program"){
      gastroProt_PRdesign[i,"Group"] <- "4w"
    } else{
      if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Training Group"){
        gastroProt_PRdesign[i,"Group"] <- "8w"
      } else {
        if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Control Group"){
          gastroProt_PRdesign[i,"Group"] <- "control"
        }
      }
    }
  }
}

gastrostatdesign$Prot_PR = gastroProt_PRdesign
gastrostatdesign$Prot_PR[gastrostatdesign$Prot_PR$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
gastrostatdesign$Prot_PR <- data.frame(row.names = rownames(gastroProt_PRdesign),
                                       "Group" = gastrostatdesign$Prot_PR$Group)

gastroProt_PHdesign <- data.frame(row.names = colnames(gastrostatdata$Prot_PH),
                                  "Sex" = rep("Female",dim(gastrostatdata$Prot_PH)[2]),
                                  "Group" = rep("1w",dim(gastrostatdata$Prot_PH)[2]))
for(i in 1:dim(gastroProt_PHdesign)[1]){
  ourid <- rownames(gastroProt_PHdesign)[i]
  gastroProt_PHdesign[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"][1]]
  if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Two-week program"){
    gastroProt_PHdesign[i,"Group"] <- "2w"
  } else{
    if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Four-week program"){
      gastroProt_PHdesign[i,"Group"] <- "4w"
    } else{
      if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Training Group"){
        gastroProt_PHdesign[i,"Group"] <- "8w"
      } else {
        if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Control Group"){
          gastroProt_PHdesign[i,"Group"] <- "control"
        }
      }
    }
  }
}

gastrostatdesign$Prot_PH = gastroProt_PHdesign
gastrostatdesign$Prot_PH[gastrostatdesign$Prot_PH$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
gastrostatdesign$Prot_PH <- data.frame(row.names = rownames(gastroProt_PHdesign),
                                       "Group" = gastrostatdesign$Prot_PH$Group)



#gastrostatResultsEQ = MultiPower(data = gastrostatdata, groups = gastrostatdesign, type = type1, d0 = 0.8, p1 = p1omics,
#                           omicPower = 0.6, averagePower = 0.8, fdr = 0.05, cost = 1, equalSize = TRUE,
#                           max.size = 200, omicCol = miscolores, dispPerc = 75)

gastrostatResultsRNAEQ = MultiPower(data = list("RNAseq" = gastrostatdata$RNAseq), groups = list("RNAseq" = gastrostatdesign$RNAseq), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)
gastrostatResultsATACEQ = MultiPower(data = list("ATACseq" = gastrostatdata$ATACseq), groups = list("ATACseq" = gastrostatdesign$ATACseq), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)
gastrostatResultsRRBSEQ = MultiPower(data = list("RRBS" = gastrostatdata$RRBS), groups = list("RRBS" = gastrostatdesign$RRBS), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)
gastrostatResultsProt_PREQ = MultiPower(data = list("Prot_PR" = gastrostatdata$Prot_PR), groups = list("Prot_PR" = gastrostatdesign$Prot_PR), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)
gastrostatResultsProt_PHEQ = MultiPower(data = list("Prot_PH" = gastrostatdata$Prot_PH), groups = list("Prot_PH" = gastrostatdesign$Prot_PH), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)

gastrostatResultsaltRNAEQ = MultiPower(data = list("RNAseq" = gastrostatdata$RNAseq), groups = list("RNAseq" = gastrostatdesign$RNAseq), type = 2,omicPower = 0.25,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)

# heart

heartstatdata = heartstatdesign = vector("list")

heartstatdata$RNAseq = heartrnanorm
heartstatdata$ATACseq = heartatacnorm
heartstatdata$RRBS = matrix(as.numeric(as.matrix(heartMmat)),nrow = dim(heartMmat)[1],ncol = dim(heartMmat)[2])
rownames(heartstatdata$RRBS) <- rownames(heartMmat)
colnames(heartstatdata$RRBS) <- colnames(heartMmat)
heartstatdata$Prot_PR = heartproprnorm
heartstatdata$Prot_PH = heartprophnorm

rm(heartrnanorm)
rm(heartatacnorm)
rm(heartMmat)
rm(heartproprnorm)
rm(heartprophnorm)
gc()

heartstatdata$RNAseq = heartstatdata$RNAseq - min(heartstatdata$RNAseq)
heartstatdata$ATACseq = heartstatdata$ATACseq - min(heartstatdata$ATACseq)
heartstatdata$RRBS = heartstatdata$RRBS - min(heartstatdata$RRBS)
heartstatdata$Prot_PR = heartstatdata$Prot_PR - min(heartstatdata$Prot_PR)
heartstatdata$Prot_PH = heartstatdata$Prot_PH - min(heartstatdata$Prot_PH)



heartRNAseqdesign <- data.frame(row.names = colnames(heartstatdata$RNAseq),
                                 "Sex" = rep("Female",dim(heartstatdata$RNAseq)[2]),
                                 "Group" = rep("1w",dim(heartstatdata$RNAseq)[2]))
for(i in 1:dim(heartRNAseqdesign)[1]){
  ourid <- rownames(heartRNAseqdesign)[i]
  heartRNAseqdesign[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"][1]]
  if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Two-week program"){
    heartRNAseqdesign[i,"Group"] <- "2w"
  } else{
    if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Four-week program"){
      heartRNAseqdesign[i,"Group"] <- "4w"
    } else{
      if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Training Group"){
        heartRNAseqdesign[i,"Group"] <- "8w"
      } else {
        if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Control Group"){
          heartRNAseqdesign[i,"Group"] <- "control"
        }
      }
    }
  }
}

heartstatdesign$RNAseq = heartRNAseqdesign
heartstatdesign$RNAseq[heartstatdesign$RNAseq$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
heartstatdesign$RNAseq <- data.frame(row.names = rownames(heartRNAseqdesign),
                                      "Group" = heartstatdesign$RNAseq$Group)


heartATACseqdesign <- data.frame(row.names = colnames(heartstatdata$ATACseq),
                                  "Sex" = rep("Female",dim(heartstatdata$ATACseq)[2]),
                                  "Group" = rep("1w",dim(heartstatdata$ATACseq)[2]))
for(i in 1:dim(heartATACseqdesign)[1]){
  ourid <- rownames(heartATACseqdesign)[i]
  heartATACseqdesign[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"][1]]
  if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Two-week program"){
    heartATACseqdesign[i,"Group"] <- "2w"
  } else{
    if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Four-week program"){
      heartATACseqdesign[i,"Group"] <- "4w"
    } else{
      if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Training Group"){
        heartATACseqdesign[i,"Group"] <- "8w"
      } else {
        if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Control Group"){
          heartATACseqdesign[i,"Group"] <- "control"
        }
      }
    }
  }
}

heartstatdesign$ATACseq = heartATACseqdesign
heartstatdesign$ATACseq[heartstatdesign$ATACseq$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
heartstatdesign$ATACseq <- data.frame(row.names = rownames(heartATACseqdesign),
                                       "Group" = heartstatdesign$ATACseq$Group)


heartRRBSdesign <- data.frame(row.names = colnames(heartstatdata$RRBS),
                              "Sex" = rep("Female",dim(heartstatdata$RRBS)[2]),
                              "Group" = rep("1w",dim(heartstatdata$RRBS)[2]))
for(i in 1:dim(heartRRBSdesign)[1]){
  ourid <- rownames(heartRRBSdesign)[i]
  heartRRBSdesign[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0027"][1]]
  if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Two-week program"){
    heartRRBSdesign[i,"Group"] <- "2w"
  } else{
    if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Four-week program"){
      heartRRBSdesign[i,"Group"] <- "4w"
    } else{
      if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Training Group"){
        heartRRBSdesign[i,"Group"] <- "8w"
      } else {
        if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Control Group"){
          heartRRBSdesign[i,"Group"] <- "control"
        }
      }
    }
  }
}

heartstatdesign$RRBS = heartRRBSdesign
heartstatdesign$RRBS[heartstatdesign$RRBS$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
heartstatdesign$RRBS <- data.frame(row.names = rownames(heartRRBSdesign),
                                   "Group" = heartstatdesign$RRBS$Group)



heartProt_PRdesign <- data.frame(row.names = colnames(heartstatdata$Prot_PR),
                                 "Sex" = rep("Female",dim(heartstatdata$Prot_PR)[2]),
                                 "Group" = rep("1w",dim(heartstatdata$Prot_PR)[2]))
for(i in 1:dim(heartProt_PRdesign)[1]){
  ourid <- rownames(heartProt_PRdesign)[i]
  heartProt_PRdesign[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0027"][1]]
  if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Two-week program"){
    heartProt_PRdesign[i,"Group"] <- "2w"
  } else{
    if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Four-week program"){
      heartProt_PRdesign[i,"Group"] <- "4w"
    } else{
      if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Training Group"){
        heartProt_PRdesign[i,"Group"] <- "8w"
      } else {
        if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Control Group"){
          heartProt_PRdesign[i,"Group"] <- "control"
        }
      }
    }
  }
}

heartstatdesign$Prot_PR = heartProt_PRdesign
heartstatdesign$Prot_PR[heartstatdesign$Prot_PR$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
heartstatdesign$Prot_PR <- data.frame(row.names = rownames(heartProt_PRdesign),
                                      "Group" = heartstatdesign$Prot_PR$Group)


heartProt_PHdesign <- data.frame(row.names = colnames(heartstatdata$Prot_PH),
                                  "Sex" = rep("Female",dim(heartstatdata$Prot_PH)[2]),
                                  "Group" = rep("1w",dim(heartstatdata$Prot_PH)[2]))
for(i in 1:dim(heartProt_PHdesign)[1]){
  ourid <- rownames(heartProt_PHdesign)[i]
  heartProt_PHdesign[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"][1]]
  if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Two-week program"){
    heartProt_PHdesign[i,"Group"] <- "2w"
  } else{
    if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Four-week program"){
      heartProt_PHdesign[i,"Group"] <- "4w"
    } else{
      if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Training Group"){
        heartProt_PHdesign[i,"Group"] <- "8w"
      } else {
        if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Control Group"){
          heartProt_PHdesign[i,"Group"] <- "control"
        }
      }
    }
  }
}

heartstatdesign$Prot_PH = heartProt_PHdesign
heartstatdesign$Prot_PH[heartstatdesign$Prot_PH$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
heartstatdesign$Prot_PH <- data.frame(row.names = rownames(heartProt_PHdesign),
                                       "Group" = heartstatdesign$Prot_PH$Group)




#heartstatResultsEQ = MultiPower(data = heartstatdata, groups = heartstatdesign, type = type1, d0 = 0.8, p1 = p1omics,
#                           omicPower = 0.6, averagePower = 0.8, fdr = 0.05, cost = 1, equalSize = TRUE,
#                           max.size = 200, omicCol = miscolores, dispPerc = 75)

heartstatResultsRNAEQ = MultiPower(data = list("RNAseq" = heartstatdata$RNAseq), groups = list("RNAseq" = heartstatdesign$RNAseq), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)
heartstatResultsATACEQ = MultiPower(data = list("ATACseq" = heartstatdata$ATACseq), groups = list("ATACseq" = heartstatdesign$ATACseq), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)
heartstatResultsRRBSEQ = MultiPower(data = list("RRBS" = heartstatdata$RRBS), groups = list("RRBS" = heartstatdesign$RRBS), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)
heartstatResultsProt_PREQ = MultiPower(data = list("Prot_PR" = heartstatdata$Prot_PR), groups = list("Prot_PR" = heartstatdesign$Prot_PR), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)
heartstatResultsProt_PHEQ = MultiPower(data = list("Prot_PH" = heartstatdata$Prot_PH), groups = list("Prot_PH" = heartstatdesign$Prot_PH), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)



# kidney

kidneystatdata = kidneystatdesign = vector("list")

kidneystatdata$RNAseq = kidneyrnanorm
kidneystatdata$ATACseq = kidneyatacnorm
kidneystatdata$RRBS = matrix(as.numeric(as.matrix(kidneyMmat)),nrow = dim(kidneyMmat)[1],ncol = dim(kidneyMmat)[2])
rownames(kidneystatdata$RRBS) <- rownames(kidneyMmat)
colnames(kidneystatdata$RRBS) <- colnames(kidneyMmat)
kidneystatdata$Prot_PR = kidneyproprnorm
kidneystatdata$Prot_PH = kidneyprophnorm

rm(kidneyrnanorm)
rm(kidneyatacnorm)
rm(kidneyMmat)
rm(kidneyproprnorm)
rm(kidneyprophnorm)
gc()

kidneystatdata$RNAseq = kidneystatdata$RNAseq - min(kidneystatdata$RNAseq)
kidneystatdata$ATACseq = kidneystatdata$ATACseq - min(kidneystatdata$ATACseq)
kidneystatdata$RRBS = kidneystatdata$RRBS - min(kidneystatdata$RRBS)
kidneystatdata$Prot_PR = kidneystatdata$Prot_PR - min(kidneystatdata$Prot_PR)
kidneystatdata$Prot_PH = kidneystatdata$Prot_PH - min(kidneystatdata$Prot_PH)

kidneyRNAseqdesign <- data.frame(row.names = colnames(kidneystatdata$RNAseq),
                                 "Sex" = rep("Female",dim(kidneystatdata$RNAseq)[2]),
                                 "Group" = rep("1w",dim(kidneystatdata$RNAseq)[2]))
for(i in 1:dim(kidneyRNAseqdesign)[1]){
  ourid <- rownames(kidneyRNAseqdesign)[i]
  kidneyRNAseqdesign[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"][1]]
  if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Two-week program"){
    kidneyRNAseqdesign[i,"Group"] <- "2w"
  } else{
    if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Four-week program"){
      kidneyRNAseqdesign[i,"Group"] <- "4w"
    } else{
      if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Training Group"){
        kidneyRNAseqdesign[i,"Group"] <- "8w"
      } else {
        if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Control Group"){
          kidneyRNAseqdesign[i,"Group"] <- "control"
        }
      }
    }
  }
}

kidneystatdesign$RNAseq = kidneyRNAseqdesign
kidneystatdesign$RNAseq[kidneystatdesign$RNAseq$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
kidneystatdesign$RNAseq <- data.frame(row.names = rownames(kidneyRNAseqdesign),
                                      "Group" = kidneystatdesign$RNAseq$Group)


kidneyATACseqdesign <- data.frame(row.names = colnames(kidneystatdata$ATACseq),
                                  "Sex" = rep("Female",dim(kidneystatdata$ATACseq)[2]),
                                  "Group" = rep("1w",dim(kidneystatdata$ATACseq)[2]))
for(i in 1:dim(kidneyATACseqdesign)[1]){
  ourid <- rownames(kidneyATACseqdesign)[i]
  kidneyATACseqdesign[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"][1]]
  if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Two-week program"){
    kidneyATACseqdesign[i,"Group"] <- "2w"
  } else{
    if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Four-week program"){
      kidneyATACseqdesign[i,"Group"] <- "4w"
    } else{
      if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Training Group"){
        kidneyATACseqdesign[i,"Group"] <- "8w"
      } else {
        if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Control Group"){
          kidneyATACseqdesign[i,"Group"] <- "control"
        }
      }
    }
  }
}

kidneystatdesign$ATACseq = kidneyATACseqdesign
kidneystatdesign$ATACseq[kidneystatdesign$ATACseq$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
kidneystatdesign$ATACseq <- data.frame(row.names = rownames(kidneyATACseqdesign),
                                       "Group" = kidneystatdesign$ATACseq$Group)


kidneyRRBSdesign <- data.frame(row.names = colnames(kidneystatdata$RRBS),
                               "Sex" = rep("Female",dim(kidneystatdata$RRBS)[2]),
                               "Group" = rep("1w",dim(kidneystatdata$RRBS)[2]))
for(i in 1:dim(kidneyRRBSdesign)[1]){
  ourid <- rownames(kidneyRRBSdesign)[i]
  kidneyRRBSdesign[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0027"][1]]
  if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Two-week program"){
    kidneyRRBSdesign[i,"Group"] <- "2w"
  } else{
    if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Four-week program"){
      kidneyRRBSdesign[i,"Group"] <- "4w"
    } else{
      if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Training Group"){
        kidneyRRBSdesign[i,"Group"] <- "8w"
      } else {
        if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Control Group"){
          kidneyRRBSdesign[i,"Group"] <- "control"
        }
      }
    }
  }
}

kidneystatdesign$RRBS = kidneyRRBSdesign
kidneystatdesign$RRBS[kidneystatdesign$RRBS$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
kidneystatdesign$RRBS <- data.frame(row.names = rownames(kidneyRRBSdesign),
                                    "Group" = kidneystatdesign$RRBS$Group)



kidneyProt_PRdesign <- data.frame(row.names = colnames(kidneystatdata$Prot_PR),
                                  "Sex" = rep("Female",dim(kidneystatdata$Prot_PR)[2]),
                                  "Group" = rep("1w",dim(kidneystatdata$Prot_PR)[2]))
for(i in 1:dim(kidneyProt_PRdesign)[1]){
  ourid <- rownames(kidneyProt_PRdesign)[i]
  kidneyProt_PRdesign[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0027"][1]]
  if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Two-week program"){
    kidneyProt_PRdesign[i,"Group"] <- "2w"
  } else{
    if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Four-week program"){
      kidneyProt_PRdesign[i,"Group"] <- "4w"
    } else{
      if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Training Group"){
        kidneyProt_PRdesign[i,"Group"] <- "8w"
      } else {
        if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Control Group"){
          kidneyProt_PRdesign[i,"Group"] <- "control"
        }
      }
    }
  }
}

kidneystatdesign$Prot_PR = kidneyProt_PRdesign
kidneystatdesign$Prot_PR[kidneystatdesign$Prot_PR$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
kidneystatdesign$Prot_PR <- data.frame(row.names = rownames(kidneyProt_PRdesign),
                                       "Group" = kidneystatdesign$Prot_PR$Group)


kidneyProt_PHdesign <- data.frame(row.names = colnames(kidneystatdata$Prot_PH),
                                  "Sex" = rep("Female",dim(kidneystatdata$Prot_PH)[2]),
                                  "Group" = rep("1w",dim(kidneystatdata$Prot_PH)[2]))
for(i in 1:dim(kidneyProt_PHdesign)[1]){
  ourid <- rownames(kidneyProt_PHdesign)[i]
  kidneyProt_PHdesign[i,"Sex"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$viallabel %in% ourid,"sex"][1]
  if(MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$viallabel %in% ourid,"group"][1] %in% "2w"){
    kidneyProt_PHdesign[i,"Group"] <- "2w"
  } else{
    if(MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$viallabel %in% ourid,"group"][1] %in% "4w"){
      kidneyProt_PHdesign[i,"Group"] <- "4w"
    } else{
      if(MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$viallabel %in% ourid,"group"][1] %in% "8w"){
        kidneyProt_PHdesign[i,"Group"] <- "8w"
      } else {
        if(MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$viallabel %in% ourid,"group"][1] %in% "control"){
          kidneyProt_PHdesign[i,"Group"] <- "control"
        }
      }
    }
  }
}

kidneystatdesign$Prot_PH = kidneyProt_PHdesign
kidneystatdesign$Prot_PH[kidneystatdesign$Prot_PH$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
kidneystatdesign$Prot_PH <- data.frame(row.names = rownames(kidneyProt_PHdesign),
                                       "Group" = kidneystatdesign$Prot_PH$Group)



#kidneystatResultsEQ = MultiPower(data = kidneystatdata, groups = kidneystatdesign, type = type1, d0 = 0.8, p1 = p1omics,
#                           omicPower = 0.6, averagePower = 0.8, fdr = 0.05, cost = 1, equalSize = TRUE,
#                           max.size = 200, omicCol = miscolores, dispPerc = 75)

kidneystatResultsRNAEQ = MultiPower(data = list("RNAseq" = kidneystatdata$RNAseq), groups = list("RNAseq" = kidneystatdesign$RNAseq), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)
kidneystatResultsATACEQ = MultiPower(data = list("ATACseq" = kidneystatdata$ATACseq), groups = list("ATACseq" = kidneystatdesign$ATACseq), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)
kidneystatResultsRRBSEQ = MultiPower(data = list("RRBS" = kidneystatdata$RRBS), groups = list("RRBS" = kidneystatdesign$RRBS), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)
kidneystatResultsProt_PREQ = MultiPower(data = list("Prot_PR" = kidneystatdata$Prot_PR), groups = list("Prot_PR" = kidneystatdesign$Prot_PR), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)
kidneystatResultsProt_PHEQ = MultiPower(data = list("Prot_PH" = kidneystatdata$Prot_PH), groups = list("Prot_PH" = kidneystatdesign$Prot_PH), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)



# liver

liverstatdata = liverstatdesign = vector("list")

liverstatdata$RNAseq = liverrnanorm
liverstatdata$ATACseq = liveratacnorm
liverstatdata$RRBS = matrix(as.numeric(as.matrix(liverMmat)),nrow = dim(liverMmat)[1],ncol = dim(liverMmat)[2])
rownames(liverstatdata$RRBS) <- rownames(liverMmat)
colnames(liverstatdata$RRBS) <- colnames(liverMmat)
liverstatdata$Prot_PR = liverproprnorm
liverstatdata$Prot_PH = liverprophnorm

rm(liverrnanorm)
rm(liveratacnorm)
rm(liverMmat)
rm(liverproprnorm)
rm(liverprophnorm)
gc()

liverstatdata$RNAseq = liverstatdata$RNAseq - min(liverstatdata$RNAseq)
liverstatdata$ATACseq = liverstatdata$ATACseq - min(liverstatdata$ATACseq)
liverstatdata$RRBS = liverstatdata$RRBS - min(liverstatdata$RRBS)
liverstatdata$Prot_PR = liverstatdata$Prot_PR - min(liverstatdata$Prot_PR)
liverstatdata$Prot_PH = liverstatdata$Prot_PH - min(liverstatdata$Prot_PH)

liverRNAseqdesign <- data.frame(row.names = colnames(liverstatdata$RNAseq),
                                 "Sex" = rep("Female",dim(liverstatdata$RNAseq)[2]),
                                 "Group" = rep("1w",dim(liverstatdata$RNAseq)[2]))
for(i in 1:dim(liverRNAseqdesign)[1]){
  ourid <- rownames(liverRNAseqdesign)[i]
  liverRNAseqdesign[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"][1]]
  if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Two-week program"){
    liverRNAseqdesign[i,"Group"] <- "2w"
  } else{
    if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Four-week program"){
      liverRNAseqdesign[i,"Group"] <- "4w"
    } else{
      if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Training Group"){
        liverRNAseqdesign[i,"Group"] <- "8w"
      } else {
        if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Control Group"){
          liverRNAseqdesign[i,"Group"] <- "control"
        }
      }
    }
  }
}

liverstatdesign$RNAseq = liverRNAseqdesign
liverstatdesign$RNAseq[liverstatdesign$RNAseq$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
liverstatdesign$RNAseq <- data.frame(row.names = rownames(liverRNAseqdesign),
                                      "Group" = liverstatdesign$RNAseq$Group)


liverATACseqdesign <- data.frame(row.names = colnames(liverstatdata$ATACseq),
                                  "Sex" = rep("Female",dim(liverstatdata$ATACseq)[2]),
                                  "Group" = rep("1w",dim(liverstatdata$ATACseq)[2]))
for(i in 1:dim(liverATACseqdesign)[1]){
  ourid <- rownames(liverATACseqdesign)[i]
  liverATACseqdesign[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"][1]]
  if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Two-week program"){
    liverATACseqdesign[i,"Group"] <- "2w"
  } else{
    if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Four-week program"){
      liverATACseqdesign[i,"Group"] <- "4w"
    } else{
      if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Training Group"){
        liverATACseqdesign[i,"Group"] <- "8w"
      } else {
        if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Control Group"){
          liverATACseqdesign[i,"Group"] <- "control"
        }
      }
    }
  }
}

liverstatdesign$ATACseq = liverATACseqdesign
liverstatdesign$ATACseq[liverstatdesign$ATACseq$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
liverstatdesign$ATACseq <- data.frame(row.names = rownames(liverATACseqdesign),
                                       "Group" = liverstatdesign$ATACseq$Group)

liverRRBSdesign <- data.frame(row.names = colnames(liverstatdata$RRBS),
                              "Sex" = rep("Female",dim(liverstatdata$RRBS)[2]),
                              "Group" = rep("1w",dim(liverstatdata$RRBS)[2]))
for(i in 1:dim(liverRRBSdesign)[1]){
  ourid <- rownames(liverRRBSdesign)[i]
  liverRRBSdesign[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0027"][1]]
  if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Two-week program"){
    liverRRBSdesign[i,"Group"] <- "2w"
  } else{
    if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Four-week program"){
      liverRRBSdesign[i,"Group"] <- "4w"
    } else{
      if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Training Group"){
        liverRRBSdesign[i,"Group"] <- "8w"
      } else {
        if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Control Group"){
          liverRRBSdesign[i,"Group"] <- "control"
        }
      }
    }
  }
}

liverstatdesign$RRBS = liverRRBSdesign
liverstatdesign$RRBS[liverstatdesign$RRBS$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
liverstatdesign$RRBS <- data.frame(row.names = rownames(liverRRBSdesign),
                                   "Group" = liverstatdesign$RRBS$Group)


liverProt_PRdesign <- data.frame(row.names = colnames(liverstatdata$Prot_PR),
                                  "Sex" = rep("Female",dim(liverstatdata$Prot_PR)[2]),
                                  "Group" = rep("1w",dim(liverstatdata$Prot_PR)[2]))
for(i in 1:dim(liverProt_PRdesign)[1]){
  ourid <- rownames(liverProt_PRdesign)[i]
  liverProt_PRdesign[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0027"][1]]
  if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Two-week program"){
    liverProt_PRdesign[i,"Group"] <- "2w"
  } else{
    if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Four-week program"){
      liverProt_PRdesign[i,"Group"] <- "4w"
    } else{
      if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Training Group"){
        liverProt_PRdesign[i,"Group"] <- "8w"
      } else {
        if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Control Group"){
          liverProt_PRdesign[i,"Group"] <- "control"
        }
      }
    }
  }
}

liverstatdesign$Prot_PR = liverProt_PRdesign
liverstatdesign$Prot_PR[liverstatdesign$Prot_PR$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
liverstatdesign$Prot_PR <- data.frame(row.names = rownames(liverProt_PRdesign),
                                       "Group" = liverstatdesign$Prot_PR$Group)

liverProt_PHdesign <- data.frame(row.names = colnames(liverstatdata$Prot_PH),
                                 "Sex" = rep("Female",dim(liverstatdata$Prot_PH)[2]),
                                 "Group" = rep("1w",dim(liverstatdata$Prot_PH)[2]))
for(i in 1:dim(liverProt_PHdesign)[1]){
  ourid <- rownames(liverProt_PHdesign)[i]
  liverProt_PHdesign[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"][1]]
  if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Two-week program"){
    liverProt_PHdesign[i,"Group"] <- "2w"
  } else{
    if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Four-week program"){
      liverProt_PHdesign[i,"Group"] <- "4w"
    } else{
      if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Training Group"){
        liverProt_PHdesign[i,"Group"] <- "8w"
      } else {
        if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Control Group"){
          liverProt_PHdesign[i,"Group"] <- "control"
        }
      }
    }
  }
}

liverstatdesign$Prot_PH = liverProt_PHdesign
liverstatdesign$Prot_PH[liverstatdesign$Prot_PH$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
liverstatdesign$Prot_PH <- data.frame(row.names = rownames(liverProt_PHdesign),
                                      "Group" = liverstatdesign$Prot_PH$Group)




#liverstatResultsEQ = MultiPower(data = liverstatdata, groups = liverstatdesign, type = type1, d0 = 0.8, p1 = p1omics,
#                           omicPower = 0.6, averagePower = 0.8, fdr = 0.05, cost = 1, equalSize = TRUE,
#                           max.size = 200, omicCol = miscolores, dispPerc = 75)

liverstatResultsRNAEQ = MultiPower(data = list("RNAseq" = liverstatdata$RNAseq), groups = list("RNAseq" = liverstatdesign$RNAseq), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)
liverstatResultsATACEQ = MultiPower(data = list("ATACseq" = liverstatdata$ATACseq), groups = list("ATACseq" = liverstatdesign$ATACseq), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)
liverstatResultsRRBSEQ = MultiPower(data = list("RRBS" = liverstatdata$RRBS), groups = list("RRBS" = liverstatdesign$RRBS), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)
liverstatResultsProt_PREQ = MultiPower(data = list("Prot_PR" = liverstatdata$Prot_PR), groups = list("Prot_PR" = liverstatdesign$Prot_PR), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)
liverstatResultsProt_PHEQ = MultiPower(data = list("Prot_PH" = liverstatdata$Prot_PH), groups = list("Prot_PH" = liverstatdesign$Prot_PH), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)



# lung

lungstatdata = lungstatdesign = vector("list")

lungstatdata$RNAseq = lungrnanorm
lungstatdata$ATACseq = lungatacnorm
lungstatdata$RRBS = matrix(as.numeric(as.matrix(lungMmat)),nrow = dim(lungMmat)[1],ncol = dim(lungMmat)[2])
rownames(lungstatdata$RRBS) <- rownames(lungMmat)
colnames(lungstatdata$RRBS) <- colnames(lungMmat)
lungstatdata$Prot_PR = lungproprnorm
lungstatdata$Prot_PH = lungprophnorm

rm(lungrnanorm)
rm(lungatacnorm)
rm(lungMmat)
rm(lungproprnorm)
rm(lungprophnorm)
gc()

lungstatdata$RNAseq = lungstatdata$RNAseq - min(lungstatdata$RNAseq)
lungstatdata$ATACseq = lungstatdata$ATACseq - min(lungstatdata$ATACseq)
lungstatdata$RRBS = lungstatdata$RRBS - min(lungstatdata$RRBS)
lungstatdata$Prot_PR = lungstatdata$Prot_PR - min(lungstatdata$Prot_PR)
lungstatdata$Prot_PH = lungstatdata$Prot_PH - min(lungstatdata$Prot_PH)


lungRNAseqdesign <- data.frame(row.names = colnames(lungstatdata$RNAseq),
                                 "Sex" = rep("Female",dim(lungstatdata$RNAseq)[2]),
                                 "Group" = rep("1w",dim(lungstatdata$RNAseq)[2]))
for(i in 1:dim(lungRNAseqdesign)[1]){
  ourid <- rownames(lungRNAseqdesign)[i]
  lungRNAseqdesign[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"][1]]
  if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Two-week program"){
    lungRNAseqdesign[i,"Group"] <- "2w"
  } else{
    if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Four-week program"){
      lungRNAseqdesign[i,"Group"] <- "4w"
    } else{
      if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Training Group"){
        lungRNAseqdesign[i,"Group"] <- "8w"
      } else {
        if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Control Group"){
          lungRNAseqdesign[i,"Group"] <- "control"
        }
      }
    }
  }
}

lungstatdesign$RNAseq = lungRNAseqdesign
lungstatdesign$RNAseq[lungstatdesign$RNAseq$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
lungstatdesign$RNAseq <- data.frame(row.names = rownames(lungRNAseqdesign),
                                      "Group" = lungstatdesign$RNAseq$Group)


lungATACseqdesign <- data.frame(row.names = colnames(lungstatdata$ATACseq),
                                  "Sex" = rep("Female",dim(lungstatdata$ATACseq)[2]),
                                  "Group" = rep("1w",dim(lungstatdata$ATACseq)[2]))
for(i in 1:dim(lungATACseqdesign)[1]){
  ourid <- rownames(lungATACseqdesign)[i]
  lungATACseqdesign[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"][1]]
  if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Two-week program"){
    lungATACseqdesign[i,"Group"] <- "2w"
  } else{
    if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Four-week program"){
      lungATACseqdesign[i,"Group"] <- "4w"
    } else{
      if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Training Group"){
        lungATACseqdesign[i,"Group"] <- "8w"
      } else {
        if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Control Group"){
          lungATACseqdesign[i,"Group"] <- "control"
        }
      }
    }
  }
}

lungstatdesign$ATACseq = lungATACseqdesign
lungstatdesign$ATACseq[lungstatdesign$ATACseq$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
lungstatdesign$ATACseq <- data.frame(row.names = rownames(lungATACseqdesign),
                                       "Group" = lungstatdesign$ATACseq$Group)



lungRRBSdesign <- data.frame(row.names = colnames(lungstatdata$RRBS),
                             "Sex" = rep("Female",dim(lungstatdata$RRBS)[2]),
                             "Group" = rep("1w",dim(lungstatdata$RRBS)[2]))
for(i in 1:dim(lungRRBSdesign)[1]){
  ourid <- rownames(lungRRBSdesign)[i]
  lungRRBSdesign[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0027"][1]]
  if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Two-week program"){
    lungRRBSdesign[i,"Group"] <- "2w"
  } else{
    if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Four-week program"){
      lungRRBSdesign[i,"Group"] <- "4w"
    } else{
      if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Training Group"){
        lungRRBSdesign[i,"Group"] <- "8w"
      } else {
        if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Control Group"){
          lungRRBSdesign[i,"Group"] <- "control"
        }
      }
    }
  }
}

lungstatdesign$RRBS = lungRRBSdesign
lungstatdesign$RRBS[lungstatdesign$RRBS$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
lungstatdesign$RRBS <- data.frame(row.names = rownames(lungRRBSdesign),
                                  "Group" = lungstatdesign$RRBS$Group)



lungProt_PRdesign <- data.frame(row.names = colnames(lungstatdata$Prot_PR),
                                  "Sex" = rep("Female",dim(lungstatdata$Prot_PR)[2]),
                                  "Group" = rep("1w",dim(lungstatdata$Prot_PR)[2]))
for(i in 1:dim(lungProt_PRdesign)[1]){
  ourid <- rownames(lungProt_PRdesign)[i]
  lungProt_PRdesign[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0027"][1]]
  if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Two-week program"){
    lungProt_PRdesign[i,"Group"] <- "2w"
  } else{
    if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Four-week program"){
      lungProt_PRdesign[i,"Group"] <- "4w"
    } else{
      if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Training Group"){
        lungProt_PRdesign[i,"Group"] <- "8w"
      } else {
        if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Control Group"){
          lungProt_PRdesign[i,"Group"] <- "control"
        }
      }
    }
  }
}

lungstatdesign$Prot_PR = lungProt_PRdesign
lungstatdesign$Prot_PR[lungstatdesign$Prot_PR$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
lungstatdesign$Prot_PR <- data.frame(row.names = rownames(lungProt_PRdesign),
                                       "Group" = lungstatdesign$Prot_PR$Group)

lungProt_PHdesign <- data.frame(row.names = colnames(lungstatdata$Prot_PH),
                                "Sex" = rep("Female",dim(lungstatdata$Prot_PH)[2]),
                                "Group" = rep("1w",dim(lungstatdata$Prot_PH)[2]))
for(i in 1:dim(lungProt_PHdesign)[1]){
  ourid <- rownames(lungProt_PHdesign)[i]
  lungProt_PHdesign[i,"Sex"] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$viallabel %in% ourid,"sex"][1]
  if(MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$viallabel %in% ourid,"group"][1] %in% "2w"){
    lungProt_PHdesign[i,"Group"] <- "2w"
  } else{
    if(MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$viallabel %in% ourid,"group"][1] %in% "4w"){
      lungProt_PHdesign[i,"Group"] <- "4w"
    } else{
      if(MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$viallabel %in% ourid,"group"][1] %in% "8w"){
        lungProt_PHdesign[i,"Group"] <- "8w"
      } else {
        if(MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$viallabel %in% ourid,"group"][1] %in% "control"){
          lungProt_PHdesign[i,"Group"] <- "control"
        }
      }
    }
  }
}

lungstatdesign$Prot_PH = lungProt_PHdesign
lungstatdesign$Prot_PH[lungstatdesign$Prot_PH$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
lungstatdesign$Prot_PH <- data.frame(row.names = rownames(lungProt_PHdesign),
                                     "Group" = lungstatdesign$Prot_PH$Group)



#lungstatResultsEQ = MultiPower(data = lungstatdata, groups = lungstatdesign, type = type1, d0 = 0.8, p1 = p1omics,
#                           omicPower = 0.6, averagePower = 0.8, fdr = 0.05, cost = 1, equalSize = TRUE,
#                           max.size = 200, omicCol = miscolores, dispPerc = 75)

lungstatResultsRNAEQ = MultiPower(data = list("RNAseq" = lungstatdata$RNAseq), groups = list("RNAseq" = lungstatdesign$RNAseq), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)
lungstatResultsATACEQ = MultiPower(data = list("ATACseq" = lungstatdata$ATACseq), groups = list("ATACseq" = lungstatdesign$ATACseq), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)
lungstatResultsRRBSEQ = MultiPower(data = list("RRBS" = lungstatdata$RRBS), groups = list("RRBS" = lungstatdesign$RRBS), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)
lungstatResultsProt_PREQ = MultiPower(data = list("Prot_PR" = lungstatdata$Prot_PR), groups = list("Prot_PR" = lungstatdesign$Prot_PR), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)
lungstatResultsProt_PHEQ = MultiPower(data = list("Prot_PH" = lungstatdata$Prot_PH), groups = list("Prot_PH" = lungstatdesign$Prot_PH), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)




# white

whitestatdata = whitestatdesign = vector("list")

whitestatdata$RNAseq = whiternanorm
whitestatdata$ATACseq = whiteatacnorm
whitestatdata$RRBS = matrix(as.numeric(as.matrix(whiteMmat)),nrow = dim(whiteMmat)[1],ncol = dim(whiteMmat)[2])
rownames(whitestatdata$RRBS) <- rownames(whiteMmat)
colnames(whitestatdata$RRBS) <- colnames(whiteMmat)
whitestatdata$Prot_PR = whiteproprnorm
whitestatdata$Prot_PH = whiteprophnorm

rm(whiternanorm)
rm(whiteatacnorm)
rm(whiteMmat)
rm(whiteproprnorm)
rm(whiteprophnorm)
gc()

whitestatdata$RNAseq = whitestatdata$RNAseq - min(whitestatdata$RNAseq)
whitestatdata$ATACseq = whitestatdata$ATACseq - min(whitestatdata$ATACseq)
whitestatdata$RRBS = whitestatdata$RRBS - min(whitestatdata$RRBS)
whitestatdata$Prot_PR = whitestatdata$Prot_PR - min(whitestatdata$Prot_PR)
whitestatdata$Prot_PH = whitestatdata$Prot_PH - min(whitestatdata$Prot_PH)


whiteRNAseqdesign <- data.frame(row.names = colnames(whitestatdata$RNAseq),
                                 "Sex" = rep("Female",dim(whitestatdata$RNAseq)[2]),
                                 "Group" = rep("1w",dim(whitestatdata$RNAseq)[2]))
for(i in 1:dim(whiteRNAseqdesign)[1]){
  ourid <- rownames(whiteRNAseqdesign)[i]
  whiteRNAseqdesign[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"][1]]
  if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Two-week program"){
    whiteRNAseqdesign[i,"Group"] <- "2w"
  } else{
    if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Four-week program"){
      whiteRNAseqdesign[i,"Group"] <- "4w"
    } else{
      if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Training Group"){
        whiteRNAseqdesign[i,"Group"] <- "8w"
      } else {
        if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Control Group"){
          whiteRNAseqdesign[i,"Group"] <- "control"
        }
      }
    }
  }
}

whitestatdesign$RNAseq = whiteRNAseqdesign
whitestatdesign$RNAseq[whitestatdesign$RNAseq$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
whitestatdesign$RNAseq <- data.frame(row.names = rownames(whiteRNAseqdesign),
                                      "Group" = whitestatdesign$RNAseq$Group)


whiteATACseqdesign <- data.frame(row.names = colnames(whitestatdata$ATACseq),
                                  "Sex" = rep("Female",dim(whitestatdata$ATACseq)[2]),
                                  "Group" = rep("1w",dim(whitestatdata$ATACseq)[2]))
for(i in 1:dim(whiteATACseqdesign)[1]){
  ourid <- rownames(whiteATACseqdesign)[i]
  whiteATACseqdesign[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"][1]]
  if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Two-week program"){
    whiteATACseqdesign[i,"Group"] <- "2w"
  } else{
    if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Four-week program"){
      whiteATACseqdesign[i,"Group"] <- "4w"
    } else{
      if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Training Group"){
        whiteATACseqdesign[i,"Group"] <- "8w"
      } else {
        if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Control Group"){
          whiteATACseqdesign[i,"Group"] <- "control"
        }
      }
    }
  }
}

whitestatdesign$ATACseq = whiteATACseqdesign
whitestatdesign$ATACseq[whitestatdesign$ATACseq$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
whitestatdesign$ATACseq <- data.frame(row.names = rownames(whiteATACseqdesign),
                                       "Group" = whitestatdesign$ATACseq$Group)


whiteRRBSdesign <- data.frame(row.names = colnames(whitestatdata$RRBS),
                              "Sex" = rep("Female",dim(whitestatdata$RRBS)[2]),
                              "Group" = rep("1w",dim(whitestatdata$RRBS)[2]))
for(i in 1:dim(whiteRRBSdesign)[1]){
  ourid <- rownames(whiteRRBSdesign)[i]
  whiteRRBSdesign[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0027"][1]]
  if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Two-week program"){
    whiteRRBSdesign[i,"Group"] <- "2w"
  } else{
    if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Four-week program"){
      whiteRRBSdesign[i,"Group"] <- "4w"
    } else{
      if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Training Group"){
        whiteRRBSdesign[i,"Group"] <- "8w"
      } else {
        if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Control Group"){
          whiteRRBSdesign[i,"Group"] <- "control"
        }
      }
    }
  }
}

whitestatdesign$RRBS = whiteRRBSdesign
whitestatdesign$RRBS[whitestatdesign$RRBS$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
whitestatdesign$RRBS <- data.frame(row.names = rownames(whiteRRBSdesign),
                                   "Group" = whitestatdesign$RRBS$Group)



whiteProt_PRdesign <- data.frame(row.names = colnames(whitestatdata$Prot_PR),
                                  "Sex" = rep("Female",dim(whitestatdata$Prot_PR)[2]),
                                  "Group" = rep("1w",dim(whitestatdata$Prot_PR)[2]))
for(i in 1:dim(whiteProt_PRdesign)[1]){
  ourid <- rownames(whiteProt_PRdesign)[i]
  whiteProt_PRdesign[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0027"][1]]
  if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Two-week program"){
    whiteProt_PRdesign[i,"Group"] <- "2w"
  } else{
    if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Four-week program"){
      whiteProt_PRdesign[i,"Group"] <- "4w"
    } else{
      if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Training Group"){
        whiteProt_PRdesign[i,"Group"] <- "8w"
      } else {
        if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Control Group"){
          whiteProt_PRdesign[i,"Group"] <- "control"
        }
      }
    }
  }
}

whitestatdesign$Prot_PR = whiteProt_PRdesign
whitestatdesign$Prot_PR[whitestatdesign$Prot_PR$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
whitestatdesign$Prot_PR <- data.frame(row.names = rownames(whiteProt_PRdesign),
                                       "Group" = whitestatdesign$Prot_PR$Group)

whiteProt_PHdesign <- data.frame(row.names = colnames(whitestatdata$Prot_PH),
                                 "Sex" = rep("Female",dim(whitestatdata$Prot_PH)[2]),
                                 "Group" = rep("1w",dim(whitestatdata$Prot_PH)[2]))
for(i in 1:dim(whiteProt_PHdesign)[1]){
  ourid <- rownames(whiteProt_PHdesign)[i]
  whiteProt_PHdesign[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"][1]]
  if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Two-week program"){
    whiteProt_PHdesign[i,"Group"] <- "2w"
  } else{
    if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Four-week program"){
      whiteProt_PHdesign[i,"Group"] <- "4w"
    } else{
      if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Training Group"){
        whiteProt_PHdesign[i,"Group"] <- "8w"
      } else {
        if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Control Group"){
          whiteProt_PHdesign[i,"Group"] <- "control"
        }
      }
    }
  }
}

whitestatdesign$Prot_PH = whiteProt_PHdesign
whitestatdesign$Prot_PH[whitestatdesign$Prot_PH$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
whitestatdesign$Prot_PH <- data.frame(row.names = rownames(whiteProt_PHdesign),
                                      "Group" = whitestatdesign$Prot_PH$Group)




#whitestatResultsEQ = MultiPower(data = whitestatdata, groups = whitestatdesign, type = type1, d0 = 0.8, p1 = p1omics,
#                           omicPower = 0.6, averagePower = 0.8, fdr = 0.05, cost = 1, equalSize = TRUE,
#                           max.size = 200, omicCol = miscolores, dispPerc = 75)

whitestatResultsRNAEQ = MultiPower(data = list("RNAseq" = whitestatdata$RNAseq), groups = list("RNAseq" = whitestatdesign$RNAseq), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)
whitestatResultsATACEQ = MultiPower(data = list("ATACseq" = whitestatdata$ATACseq), groups = list("ATACseq" = whitestatdesign$ATACseq), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)
whitestatResultsRRBSEQ = MultiPower(data = list("RRBS" = whitestatdata$RRBS), groups = list("RRBS" = whitestatdesign$RRBS), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)
whitestatResultsProt_PREQ = MultiPower(data = list("Prot_PR" = whitestatdata$Prot_PR), groups = list("Prot_PR" = whitestatdesign$Prot_PR), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)
whitestatResultsProt_PHEQ = MultiPower(data = list("Prot_PH" = whitestatdata$Prot_PH), groups = list("Prot_PH" = whitestatdesign$Prot_PH), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)






# hippo

hippostatdata = hippostatdesign = vector("list")

hippostatdata$RNAseq = hippornanorm
hippostatdata$ATACseq = hippoatacnorm
hippostatdata$RRBS = matrix(as.numeric(as.matrix(hippoMmat)),nrow = dim(hippoMmat)[1],ncol = dim(hippoMmat)[2])
rownames(hippostatdata$RRBS) <- rownames(hippoMmat)
colnames(hippostatdata$RRBS) <- colnames(hippoMmat)

rm(hippornanorm)
rm(hippoatacnorm)
rm(hippoMmat)
gc()

hippostatdata$RNAseq = hippostatdata$RNAseq - min(hippostatdata$RNAseq)
hippostatdata$ATACseq = hippostatdata$ATACseq - min(hippostatdata$ATACseq)
hippostatdata$RRBS = hippostatdata$RRBS - min(hippostatdata$RRBS)

hippoRNAseqdesign <- data.frame(row.names = colnames(hippostatdata$RNAseq),
                                 "Sex" = rep("Female",dim(hippostatdata$RNAseq)[2]),
                                 "Group" = rep("1w",dim(hippostatdata$RNAseq)[2]))
for(i in 1:dim(hippoRNAseqdesign)[1]){
  ourid <- rownames(hippoRNAseqdesign)[i]
  hippoRNAseqdesign[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"][1]]
  if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Two-week program"){
    hippoRNAseqdesign[i,"Group"] <- "2w"
  } else{
    if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Four-week program"){
      hippoRNAseqdesign[i,"Group"] <- "4w"
    } else{
      if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Training Group"){
        hippoRNAseqdesign[i,"Group"] <- "8w"
      } else {
        if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Control Group"){
          hippoRNAseqdesign[i,"Group"] <- "control"
        }
      }
    }
  }
}

hippostatdesign$RNAseq = hippoRNAseqdesign
hippostatdesign$RNAseq[hippostatdesign$RNAseq$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
hippostatdesign$RNAseq <- data.frame(row.names = rownames(hippoRNAseqdesign),
                                      "Group" = hippostatdesign$RNAseq$Group)


hippoATACseqdesign <- data.frame(row.names = colnames(hippostatdata$ATACseq),
                                  "Sex" = rep("Female",dim(hippostatdata$ATACseq)[2]),
                                  "Group" = rep("1w",dim(hippostatdata$ATACseq)[2]))
for(i in 1:dim(hippoATACseqdesign)[1]){
  ourid <- rownames(hippoATACseqdesign)[i]
  hippoATACseqdesign[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"][1]]
  if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Two-week program"){
    hippoATACseqdesign[i,"Group"] <- "2w"
  } else{
    if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Four-week program"){
      hippoATACseqdesign[i,"Group"] <- "4w"
    } else{
      if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Training Group"){
        hippoATACseqdesign[i,"Group"] <- "8w"
      } else {
        if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Control Group"){
          hippoATACseqdesign[i,"Group"] <- "control"
        }
      }
    }
  }
}

hippostatdesign$ATACseq = hippoATACseqdesign
hippostatdesign$ATACseq[hippostatdesign$ATACseq$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
hippostatdesign$ATACseq <- data.frame(row.names = rownames(hippoATACseqdesign),
                                       "Group" = hippostatdesign$ATACseq$Group)

hippoRRBSdesign <- data.frame(row.names = colnames(hippostatdata$RRBS),
                               "Sex" = rep("Female",dim(hippostatdata$RRBS)[2]),
                               "Group" = rep("1w",dim(hippostatdata$RRBS)[2]))
for(i in 1:dim(hippoRRBSdesign)[1]){
  ourid <- rownames(hippoRRBSdesign)[i]
  hippoRRBSdesign[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0027"][1]]
  if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Two-week program"){
    hippoRRBSdesign[i,"Group"] <- "2w"
  } else{
    if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Four-week program"){
      hippoRRBSdesign[i,"Group"] <- "4w"
    } else{
      if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Training Group"){
        hippoRRBSdesign[i,"Group"] <- "8w"
      } else {
        if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Control Group"){
          hippoRRBSdesign[i,"Group"] <- "control"
        }
      }
    }
  }
}

hippostatdesign$RRBS = hippoRRBSdesign
hippostatdesign$RRBS[hippostatdesign$RRBS$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
hippostatdesign$RRBS <- data.frame(row.names = rownames(hippoRRBSdesign),
                                    "Group" = hippostatdesign$RRBS$Group)




#hippostatResultsEQ = MultiPower(data = hippostatdata, groups = hippostatdesign, type = type1, d0 = 0.8, p1 = p1omics,
#                           omicPower = 0.6, averagePower = 0.8, fdr = 0.05, cost = 1, equalSize = TRUE,
#                           max.size = 200, omicCol = miscolores, dispPerc = 75)

hippostatResultsRNAEQ = MultiPower(data = list("RNAseq" = hippostatdata$RNAseq), groups = list("RNAseq" = hippostatdesign$RNAseq), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)
hippostatResultsATACEQ = MultiPower(data = list("ATACseq" = hippostatdata$ATACseq), groups = list("ATACseq" = hippostatdesign$ATACseq), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)
hippostatResultsRRBSEQ = MultiPower(data = list("RRBS" = hippostatdata$RRBS), groups = list("RRBS" = hippostatdesign$RRBS), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)




# brown

brownstatdata = brownstatdesign = vector("list")

brownstatdata$RNAseq = brownrnanorm
brownstatdata$ATACseq = brownatacnorm
brownstatdata$RRBS = matrix(as.numeric(as.matrix(brownMmat)),nrow = dim(brownMmat)[1],ncol = dim(brownMmat)[2])
rownames(brownstatdata$RRBS) <- rownames(brownMmat)
colnames(brownstatdata$RRBS) <- colnames(brownMmat)

rm(brownrnanorm)
rm(brownatacnorm)
rm(brownMmat)
gc()

brownstatdata$RNAseq = brownstatdata$RNAseq - min(brownstatdata$RNAseq)
brownstatdata$ATACseq = brownstatdata$ATACseq - min(brownstatdata$ATACseq)
brownstatdata$RRBS = brownstatdata$RRBS - min(brownstatdata$RRBS)

brownRNAseqdesign <- data.frame(row.names = colnames(brownstatdata$RNAseq),
                                 "Sex" = rep("Female",dim(brownstatdata$RNAseq)[2]),
                                 "Group" = rep("1w",dim(brownstatdata$RNAseq)[2]))
for(i in 1:dim(brownRNAseqdesign)[1]){
  ourid <- rownames(brownRNAseqdesign)[i]
  brownRNAseqdesign[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"][1]]
  if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Two-week program"){
    brownRNAseqdesign[i,"Group"] <- "2w"
  } else{
    if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Four-week program"){
      brownRNAseqdesign[i,"Group"] <- "4w"
    } else{
      if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Training Group"){
        brownRNAseqdesign[i,"Group"] <- "8w"
      } else {
        if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Control Group"){
          brownRNAseqdesign[i,"Group"] <- "control"
        }
      }
    }
  }
}

brownstatdesign$RNAseq = brownRNAseqdesign
brownstatdesign$RNAseq[brownstatdesign$RNAseq$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
brownstatdesign$RNAseq <- data.frame(row.names = rownames(brownRNAseqdesign),
                                      "Group" = brownstatdesign$RNAseq$Group)


brownATACseqdesign <- data.frame(row.names = colnames(brownstatdata$ATACseq),
                                  "Sex" = rep("Female",dim(brownstatdata$ATACseq)[2]),
                                  "Group" = rep("1w",dim(brownstatdata$ATACseq)[2]))
for(i in 1:dim(brownATACseqdesign)[1]){
  ourid <- rownames(brownATACseqdesign)[i]
  brownATACseqdesign[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0027"][1]]
  if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Two-week program"){
    brownATACseqdesign[i,"Group"] <- "2w"
  } else{
    if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Four-week program"){
      brownATACseqdesign[i,"Group"] <- "4w"
    } else{
      if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Training Group"){
        brownATACseqdesign[i,"Group"] <- "8w"
      } else {
        if(pass1bphenodata[pass1bphenodata$pass1bf0004 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Control Group"){
          brownATACseqdesign[i,"Group"] <- "control"
        }
      }
    }
  }
}

brownstatdesign$ATACseq = brownATACseqdesign
brownstatdesign$ATACseq[brownstatdesign$ATACseq$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
brownstatdesign$ATACseq <- data.frame(row.names = rownames(brownATACseqdesign),
                                       "Group" = brownstatdesign$ATACseq$Group)

brownRRBSdesign <- data.frame(row.names = colnames(brownstatdata$RRBS),
                               "Sex" = rep("Female",dim(brownstatdata$RRBS)[2]),
                               "Group" = rep("1w",dim(brownstatdata$RRBS)[2]))
for(i in 1:dim(brownRRBSdesign)[1]){
  ourid <- rownames(brownRRBSdesign)[i]
  brownRRBSdesign[i,"Sex"] <- c("Female","Male")[pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0027"][1]]
  if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Two-week program"){
    brownRRBSdesign[i,"Group"] <- "2w"
  } else{
    if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Four-week program"){
      brownRRBSdesign[i,"Group"] <- "4w"
    } else{
      if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Training Group"){
        brownRRBSdesign[i,"Group"] <- "8w"
      } else {
        if(pass1bphenodata[pass1bphenodata$pass1bf0001 %in% ourid,"pass1bf0011"][1] %in% "Eight-week program Control Group"){
          brownRRBSdesign[i,"Group"] <- "control"
        }
      }
    }
  }
}

brownstatdesign$RRBS = brownRRBSdesign
brownstatdesign$RRBS[brownstatdesign$RRBS$Group %in% c("1w","2w","4w","8w"),"Group"] <- "exercise"
brownstatdesign$RRBS <- data.frame(row.names = rownames(brownRRBSdesign),
                                    "Group" = brownstatdesign$RRBS$Group)



#brownstatResultsEQ = MultiPower(data = brownstatdata, groups = brownstatdesign, type = type1, d0 = 0.8, p1 = p1omics,
#                           omicPower = 0.6, averagePower = 0.8, fdr = 0.05, cost = 1, equalSize = TRUE,
#                           max.size = 200, omicCol = miscolores, dispPerc = 75)

brownstatResultsRNAEQ = MultiPower(data = list("RNAseq" = brownstatdata$RNAseq), groups = list("RNAseq" = brownstatdesign$RNAseq), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)
brownstatResultsATACEQ = MultiPower(data = list("ATACseq" = brownstatdata$ATACseq), groups = list("ATACseq" = brownstatdesign$ATACseq), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)
brownstatResultsRRBSEQ = MultiPower(data = list("RRBS" = brownstatdata$RRBS), groups = list("RRBS" = brownstatdesign$RRBS), type = 2,omicPower = 0.6,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = FALSE)

#gastrostatpostRNAEQ = postMultiPower(optResults = gastrostatResultsRNAEQ,max.size = 5, omicCol = miscolores) 
#gastrostatpostATACEQ = postMultiPower(optResults = gastrostatResultsATACEQ,max.size = 5, omicCol = miscolores) 

allstatresults <- rbind(gastrostatResultsRNAEQ$summary,
                        gastrostatResultsATACEQ$summary,
                        gastrostatResultsRRBSEQ$summary,
                        heartstatResultsRNAEQ$summary,
                        heartstatResultsATACEQ$summary,
                        heartstatResultsRRBSEQ$summary,
                        hippostatResultsRNAEQ$summary,
                        hippostatResultsATACEQ$summary,
                        hippostatResultsRRBSEQ$summary,
                        kidneystatResultsRNAEQ$summary,
                        kidneystatResultsATACEQ$summary,
                        kidneystatResultsRRBSEQ$summary,
                        liverstatResultsRNAEQ$summary,
                        liverstatResultsATACEQ$summary,
                        liverstatResultsRRBSEQ$summary,
                        lungstatResultsRNAEQ$summary,
                        lungstatResultsATACEQ$summary,
                        lungstatResultsRRBSEQ$summary,
                        brownstatResultsRNAEQ$summary,
                        brownstatResultsATACEQ$summary,
                        brownstatResultsRRBSEQ$summary,
                        whitestatResultsRNAEQ$summary,
                        whitestatResultsATACEQ$summary,
                        whitestatResultsRRBSEQ$summary)

rownames(allstatresults) <- c("SKM-GN_RNA",
                              "SKM-GN_ATAC",
                              "SKM-GN_RRBS",
                              "HEART-RNA",
                              "HEART-ATAC",
                              "HEART_RRBS",
                              "HIPPOC-RNA",
                              "HIPPOC-ATAC",
                              "HIPPOC_RRBS",
                              "KIDNEY-RNA",
                              "KIDNEY-ATAC",
                              "KIDNEY_RRBS",
                              "LIVER-RNA",
                              "LIVER-ATAC",
                              "LIVER_RRBS",
                              "LUNG-RNA",
                              "LUNG-ATAC",
                              "LUNG_RRBS",
                              "BAT-RNA",
                              "BAT-ATAC",
                              "BAT_RRBS",
                              "WAT-SC-RNA",
                              "WAT-SC-ATAC",
                              "WAT-SC_RRBS")


##########
# Principal Component Analysis Aside
####
gastrolongcode <- "Gastrocnemius"
gastroshortcode <- "SKM-GN"
heartlongcode <- "Heart"
heartshortcode <- "HEART"
hippolongcode <- "Hippocampus"
hipposhortcode <- "HIPPOC"
kidneylongcode <- "Kidney"
kidneyshortcode <- "KIDNEY"
liverlongcode <- "Liver"
livershortcode <- "LIVER"
lunglongcode <- "Lung"
lungshortcode <- "LUNG"
brownlongcode <- "Brown_Adipose"
brownshortcode <- "BAT"
whitelongcode <- "White_Adipose"
whiteshortcode <- "WAT-SC"

####
#gastro

gastrornafin <- gastrostatdata$RNAseq[,!(colnames(gastrostatdata$RNAseq) %in% MotrpacRatTraining6moData::OUTLIERS)]
gastrornapca_res <- prcomp(t(gastrornafin),scale. = F,center = T)
gastroatacfin <- gastrostatdata$ATACseq[,!(colnames(gastrostatdata$ATACseq) %in% MotrpacRatTraining6moData::OUTLIERS)]
gastroatacpca_res <- prcomp(t(gastroatacfin),scale. = F,center = T)
gastroprot_phfin <- gastrostatdata$Prot_PH[,!(colnames(gastrostatdata$Prot_PH) %in% MotrpacRatTraining6moData::OUTLIERS)]
gastroprot_phpca_res <- prcomp(t(as.matrix(gastroprot_phfin)),scale. = F,center = T)

gastrorrbstemp <- gastrostatdata$RRBS
for(i in 1:dim(gastrorrbstemp)[2]){
  ourpid <- colnames(gastrorrbstemp)[i]
  ourviallabel <- MotrpacRatTraining6moData::METHYL_META[MotrpacRatTraining6moData::METHYL_META$PID %in% ourpid & MotrpacRatTraining6moData::METHYL_META$Tissue %in% gastrolongcode,"viallabel"]
  colnames(gastrorrbstemp)[i] <- ourviallabel
}
gastrorrbsfin <- gastrorrbstemp[,!(colnames(gastrorrbstemp) %in% MotrpacRatTraining6moData::OUTLIERS)]
gastrorrbspca_res <- prcomp(t(gastrorrbsfin),scale. = F,center = T)

gastroprot_prtemp <- gastrostatdata$Prot_PR
for(i in 1:dim(gastroprot_prtemp)[2]){
  ourpid <- colnames(gastroprot_prtemp)[i]
  ourviallabel <- intersect(MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid & MotrpacRatTraining6moData::PHENO$tissue %in% gastroshortcode,"viallabel"],MotrpacRatTraining6moData::PROT_META$viallabel)
  colnames(gastroprot_prtemp)[i] <- ourviallabel
}
gastroprot_prfin <- gastroprot_prtemp[,!(colnames(gastroprot_prtemp) %in% MotrpacRatTraining6moData::OUTLIERS)]
gastroprot_prpca_res <- prcomp(t(gastroprot_prfin),scale. = F,center = T)

png(file = "Supplemental Figure S6A_112624.png",width = 4,height = 4,units = "in",res = 600)
autoplot(gastrornapca_res,data = gastroRNAseqdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(gastroshortcode," RNAseq PCA Plot",sep = ""))
dev.off()
pdf(file = "Supplemental Figure S6A_112624.pdf",width = 4,height = 4)
autoplot(gastrornapca_res,data = gastroRNAseqdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(gastroshortcode," RNAseq PCA Plot",sep = ""))
dev.off()

png(file = "Supplemental Figure S7A_112624.png",width = 4,height = 4,units = "in",res = 600)
autoplot(gastroatacpca_res,data = gastroATACseqdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(gastroshortcode," ATACseq PCA Plot",sep = ""))
dev.off()
pdf(file = "Supplemental Figure S7A_112624.pdf",width = 4,height = 4)
autoplot(gastroatacpca_res,data = gastroATACseqdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(gastroshortcode," ATACseq PCA Plot",sep = ""))
dev.off()

png(file = "Supplemental Figure S8A_112624.png",width = 4,height = 4,units = "in",res = 600)
autoplot(gastrorrbspca_res,data = gastroRRBSdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(gastroshortcode," RRBS PCA Plot",sep = ""))
dev.off()
pdf(file = "Supplemental Figure S8A_112624.pdf",width = 4,height = 4)
autoplot(gastrorrbspca_res,data = gastroRRBSdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(gastroshortcode," RRBS PCA Plot",sep = ""))
dev.off()

#png(file = "Supplemental Figure X4A_5724.png",width = 4,height = 4,units = "in",res = 600)
#autoplot(gastroprot_prpca_res,data = gastroProt_PRdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(gastroshortcode," PROT_PR PCA Plot",sep = ""))
#dev.off()
#pdf(file = "Supplemental Figure X4A_5724.pdf",width = 4,height = 4)
#autoplot(gastroprot_prpca_res,data = gastroProt_PRdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(gastroshortcode," PROT_PR PCA Plot",sep = ""))
#dev.off()

#png(file = "Supplemental Figure X5A_5724.png",width = 4,height = 4,units = "in",res = 600)
#autoplot(gastroprot_phpca_res,data = gastroProt_PHdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(gastroshortcode," PROT_PH PCA Plot",sep = ""))
#dev.off()
#pdf(file = "Supplemental Figure X5A_5724.pdf",width = 4,height = 4)
#autoplot(gastroprot_phpca_res,data = gastroProt_PHdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(gastroshortcode," PROT_PH PCA Plot",sep = ""))
#dev.off()


####
#heart

heartrnafin <- heartstatdata$RNAseq[,!(colnames(heartstatdata$RNAseq) %in% MotrpacRatTraining6moData::OUTLIERS)]
heartrnapca_res <- prcomp(t(heartrnafin),scale. = F,center = T)
heartatacfin <- heartstatdata$ATACseq[,!(colnames(heartstatdata$ATACseq) %in% MotrpacRatTraining6moData::OUTLIERS)]
heartatacpca_res <- prcomp(t(heartatacfin),scale. = F,center = T)
heartprot_phfin <- heartstatdata$Prot_PH[,!(colnames(heartstatdata$Prot_PH) %in% MotrpacRatTraining6moData::OUTLIERS)]
heartprot_phpca_res <- prcomp(t(heartprot_phfin),scale. = F,center = T)

heartrrbstemp <- heartstatdata$RRBS
for(i in 1:dim(heartrrbstemp)[2]){
  ourpid <- colnames(heartrrbstemp)[i]
  ourviallabel <- MotrpacRatTraining6moData::METHYL_META[MotrpacRatTraining6moData::METHYL_META$PID %in% ourpid & MotrpacRatTraining6moData::METHYL_META$Tissue %in% heartlongcode,"viallabel"]
  colnames(heartrrbstemp)[i] <- ourviallabel
}
heartrrbsfin <- heartrrbstemp[,!(colnames(heartrrbstemp) %in% MotrpacRatTraining6moData::OUTLIERS)]
heartrrbspca_res <- prcomp(t(heartrrbsfin),scale. = F,center = T)

heartprot_prtemp <- heartstatdata$Prot_PR
for(i in 1:dim(heartprot_prtemp)[2]){
  ourpid <- colnames(heartprot_prtemp)[i]
  ourviallabel <- intersect(MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid & MotrpacRatTraining6moData::PHENO$tissue %in% heartshortcode,"viallabel"],MotrpacRatTraining6moData::PROT_META$viallabel)
  colnames(heartprot_prtemp)[i] <- ourviallabel
}
heartprot_prfin <- heartprot_prtemp[,!(colnames(heartprot_prtemp) %in% MotrpacRatTraining6moData::OUTLIERS)]
heartprot_prpca_res <- prcomp(t(heartprot_prfin),scale. = F,center = T)

png(file = "Supplemental Figure S6B_112624.png",width = 4,height = 4,units = "in",res = 600)
autoplot(heartrnapca_res,data = heartRNAseqdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(heartshortcode," RNAseq PCA Plot",sep = ""))
dev.off()
pdf(file = "Supplemental Figure S6B_112624.pdf",width = 4,height = 4)
autoplot(heartrnapca_res,data = heartRNAseqdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(heartshortcode," RNAseq PCA Plot",sep = ""))
dev.off()

png(file = "Supplemental Figure S7B_112624.png",width = 4,height = 4,units = "in",res = 600)
autoplot(heartatacpca_res,data = heartATACseqdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(heartshortcode," ATACseq PCA Plot",sep = ""))
dev.off()
pdf(file = "Supplemental Figure S7B_112624.pdf",width = 4,height = 4)
autoplot(heartatacpca_res,data = heartATACseqdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(heartshortcode," ATACseq PCA Plot",sep = ""))
dev.off()

png(file = "Supplemental Figure S8B_112624.png",width = 4,height = 4,units = "in",res = 600)
autoplot(heartrrbspca_res,data = heartRRBSdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(heartshortcode," RRBS PCA Plot",sep = ""))
dev.off()
pdf(file = "Supplemental Figure S8B_112624.pdf",width = 4,height = 4)
autoplot(heartrrbspca_res,data = heartRRBSdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(heartshortcode," RRBS PCA Plot",sep = ""))
dev.off()

#png(file = "Supplemental Figure X4B_5724.png",width = 4,height = 4,units = "in",res = 600)
#autoplot(heartprot_prpca_res,data = heartProt_PRdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(heartshortcode," PROT_PR PCA Plot",sep = ""))
#dev.off()
#pdf(file = "Supplemental Figure X4B_5724.pdf",width = 4,height = 4)
#autoplot(heartprot_prpca_res,data = heartProt_PRdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(heartshortcode," PROT_PR PCA Plot",sep = ""))
#dev.off()

#png(file = "Supplemental Figure X5B_5724.png",width = 4,height = 4,units = "in",res = 600)
#autoplot(heartprot_phpca_res,data = heartProt_PHdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(heartshortcode," PROT_PH PCA Plot",sep = ""))
#dev.off()
#pdf(file = "Supplemental Figure X5B_5724.pdf",width = 4,height = 4)
#autoplot(heartprot_phpca_res,data = heartProt_PHdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(heartshortcode," PROT_PH PCA Plot",sep = ""))
#dev.off()


####
#kidney

kidneyrnafin <- kidneystatdata$RNAseq[,!(colnames(kidneystatdata$RNAseq) %in% MotrpacRatTraining6moData::OUTLIERS)]
kidneyrnapca_res <- prcomp(t(kidneyrnafin),scale. = F,center = T)
kidneyatacfin <- kidneystatdata$ATACseq[,!(colnames(kidneystatdata$ATACseq) %in% MotrpacRatTraining6moData::OUTLIERS)]
kidneyatacpca_res <- prcomp(t(kidneyatacfin),scale. = F,center = T)
kidneyprot_phfin <- kidneystatdata$Prot_PH[,!(colnames(kidneystatdata$Prot_PH) %in% MotrpacRatTraining6moData::OUTLIERS)]
kidneyprot_phpca_res <- prcomp(t(kidneyprot_phfin),scale. = F,center = T)

kidneyrrbstemp <- kidneystatdata$RRBS
for(i in 1:dim(kidneyrrbstemp)[2]){
  ourpid <- colnames(kidneyrrbstemp)[i]
  ourviallabel <- MotrpacRatTraining6moData::METHYL_META[MotrpacRatTraining6moData::METHYL_META$PID %in% ourpid & MotrpacRatTraining6moData::METHYL_META$Tissue %in% kidneylongcode,"viallabel"]
  colnames(kidneyrrbstemp)[i] <- ourviallabel
}
kidneyrrbsfin <- kidneyrrbstemp[,!(colnames(kidneyrrbstemp) %in% MotrpacRatTraining6moData::OUTLIERS)]
kidneyrrbspca_res <- prcomp(t(kidneyrrbsfin),scale. = F,center = T)

kidneyprot_prtemp <- kidneystatdata$Prot_PR
for(i in 1:dim(kidneyprot_prtemp)[2]){
  ourpid <- colnames(kidneyprot_prtemp)[i]
  ourviallabel <- intersect(MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid & MotrpacRatTraining6moData::PHENO$tissue %in% kidneyshortcode,"viallabel"],MotrpacRatTraining6moData::PROT_META$viallabel)
  colnames(kidneyprot_prtemp)[i] <- ourviallabel
}
kidneyprot_prfin <- kidneyprot_prtemp[,!(colnames(kidneyprot_prtemp) %in% MotrpacRatTraining6moData::OUTLIERS)]
kidneyprot_prpca_res <- prcomp(t(kidneyprot_prfin),scale. = F,center = T)

png(file = "Supplemental Figure S6D_112624.png",width = 4,height = 4,units = "in",res = 600)
autoplot(kidneyrnapca_res,data = kidneyRNAseqdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(kidneyshortcode," RNAseq PCA Plot",sep = ""))
dev.off()
pdf(file = "Supplemental Figure S6D_112624.pdf",width = 4,height = 4)
autoplot(kidneyrnapca_res,data = kidneyRNAseqdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(kidneyshortcode," RNAseq PCA Plot",sep = ""))
dev.off()

png(file = "Supplemental Figure S7D_112624.png",width = 4,height = 4,units = "in",res = 600)
autoplot(kidneyatacpca_res,data = kidneyATACseqdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(kidneyshortcode," ATACseq PCA Plot",sep = ""))
dev.off()
pdf(file = "Supplemental Figure S7D_112624.pdf",width = 4,height = 4)
autoplot(kidneyatacpca_res,data = kidneyATACseqdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(kidneyshortcode," ATACseq PCA Plot",sep = ""))
dev.off()

png(file = "Supplemental Figure S8D_112624.png",width = 4,height = 4,units = "in",res = 600)
autoplot(kidneyrrbspca_res,data = kidneyRRBSdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(kidneyshortcode," RRBS PCA Plot",sep = ""))
dev.off()
pdf(file = "Supplemental Figure S8D_112624.pdf",width = 4,height = 4)
autoplot(kidneyrrbspca_res,data = kidneyRRBSdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(kidneyshortcode," RRBS PCA Plot",sep = ""))
dev.off()

#png(file = "Supplemental Figure X4C_5724.png",width = 4,height = 4,units = "in",res = 600)
#autoplot(kidneyprot_prpca_res,data = kidneyProt_PRdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(kidneyshortcode," PROT_PR PCA Plot",sep = ""))
#dev.off()
#pdf(file = "Supplemental Figure X4C_5724.pdf",width = 4,height = 4)
#autoplot(kidneyprot_prpca_res,data = kidneyProt_PRdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(kidneyshortcode," PROT_PR PCA Plot",sep = ""))
#dev.off()

#png(file = "Supplemental Figure X5C_5724.png",width = 4,height = 4,units = "in",res = 600)
#autoplot(kidneyprot_phpca_res,data = kidneyProt_PHdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(kidneyshortcode," PROT_PH PCA Plot",sep = ""))
#dev.off()
#pdf(file = "Supplemental Figure X5C_5724.pdf",width = 4,height = 4)
#autoplot(kidneyprot_phpca_res,data = kidneyProt_PHdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(kidneyshortcode," PROT_PH PCA Plot",sep = ""))
#dev.off()


####
#liver

liverrnafin <- liverstatdata$RNAseq[,!(colnames(liverstatdata$RNAseq) %in% MotrpacRatTraining6moData::OUTLIERS)]
liverrnapca_res <- prcomp(t(liverrnafin),scale. = F,center = T)
liveratacfin <- liverstatdata$ATACseq[,!(colnames(liverstatdata$ATACseq) %in% MotrpacRatTraining6moData::OUTLIERS)]
liveratacpca_res <- prcomp(t(liveratacfin),scale. = F,center = T)
liverprot_phfin <- liverstatdata$Prot_PH[,!(colnames(liverstatdata$Prot_PH) %in% MotrpacRatTraining6moData::OUTLIERS)]
liverprot_phpca_res <- prcomp(t(liverprot_phfin),scale. = F,center = T)

liverrrbstemp <- liverstatdata$RRBS
for(i in 1:dim(liverrrbstemp)[2]){
  ourpid <- colnames(liverrrbstemp)[i]
  ourviallabel <- MotrpacRatTraining6moData::METHYL_META[MotrpacRatTraining6moData::METHYL_META$PID %in% ourpid & MotrpacRatTraining6moData::METHYL_META$Tissue %in% liverlongcode,"viallabel"]
  colnames(liverrrbstemp)[i] <- ourviallabel
}
liverrrbsfin <- liverrrbstemp[,!(colnames(liverrrbstemp) %in% MotrpacRatTraining6moData::OUTLIERS)]
liverrrbspca_res <- prcomp(t(liverrrbsfin),scale. = F,center = T)

liverprot_prtemp <- liverstatdata$Prot_PR
for(i in 1:dim(liverprot_prtemp)[2]){
  ourpid <- colnames(liverprot_prtemp)[i]
  ourviallabel <- intersect(MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid & MotrpacRatTraining6moData::PHENO$tissue %in% livershortcode,"viallabel"],MotrpacRatTraining6moData::PROT_META$viallabel)
  colnames(liverprot_prtemp)[i] <- ourviallabel
}
liverprot_prfin <- liverprot_prtemp[,!(colnames(liverprot_prtemp) %in% MotrpacRatTraining6moData::OUTLIERS)]
liverprot_prpca_res <- prcomp(t(liverprot_prfin),scale. = F,center = T)

png(file = "Supplemental Figure S6E_112624.png",width = 4,height = 4,units = "in",res = 600)
autoplot(liverrnapca_res,data = liverRNAseqdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(livershortcode," RNAseq PCA Plot",sep = ""))
dev.off()
pdf(file = "Supplemental Figure S6E_112624.pdf",width = 4,height = 4)
autoplot(liverrnapca_res,data = liverRNAseqdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(livershortcode," RNAseq PCA Plot",sep = ""))
dev.off()

png(file = "Supplemental Figure S7E_112624.png",width = 4,height = 4,units = "in",res = 600)
autoplot(liveratacpca_res,data = liverATACseqdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(livershortcode," ATACseq PCA Plot",sep = ""))
dev.off()
pdf(file = "Supplemental Figure S7E_112624.pdf",width = 4,height = 4)
autoplot(liveratacpca_res,data = liverATACseqdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(livershortcode," ATACseq PCA Plot",sep = ""))
dev.off()

png(file = "Supplemental Figure S8E_112624.png",width = 4,height = 4,units = "in",res = 600)
autoplot(liverrrbspca_res,data = liverRRBSdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(livershortcode," RRBS PCA Plot",sep = ""))
dev.off()
pdf(file = "Supplemental Figure S8E_112624.pdf",width = 4,height = 4)
autoplot(liverrrbspca_res,data = liverRRBSdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(livershortcode," RRBS PCA Plot",sep = ""))
dev.off()

#png(file = "Supplemental Figure X4D_5724.png",width = 4,height = 4,units = "in",res = 600)
#autoplot(liverprot_prpca_res,data = liverProt_PRdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(livershortcode," PROT_PR PCA Plot",sep = ""))
#dev.off()
#pdf(file = "Supplemental Figure X4D_5724.pdf",width = 4,height = 4)
#autoplot(liverprot_prpca_res,data = liverProt_PRdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(livershortcode," PROT_PR PCA Plot",sep = ""))
#dev.off()

#png(file = "Supplemental Figure X5D_5724.png",width = 4,height = 4,units = "in",res = 600)
#autoplot(liverprot_phpca_res,data = liverProt_PHdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(livershortcode," PROT_PH PCA Plot",sep = ""))
#dev.off()
#pdf(file = "Supplemental Figure X5D_5724.pdf",width = 4,height = 4)
#autoplot(liverprot_phpca_res,data = liverProt_PHdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(livershortcode," PROT_PH PCA Plot",sep = ""))
#dev.off()


####
#lung

lungrnafin <- lungstatdata$RNAseq[,!(colnames(lungstatdata$RNAseq) %in% MotrpacRatTraining6moData::OUTLIERS)]
lungrnapca_res <- prcomp(t(lungrnafin),scale. = F,center = T)
lungatacfin <- lungstatdata$ATACseq[,!(colnames(lungstatdata$ATACseq) %in% MotrpacRatTraining6moData::OUTLIERS)]
lungatacpca_res <- prcomp(t(lungatacfin),scale. = F,center = T)
lungprot_phfin <- lungstatdata$Prot_PH[,!(colnames(lungstatdata$Prot_PH) %in% MotrpacRatTraining6moData::OUTLIERS)]
lungprot_phpca_res <- prcomp(t(lungprot_phfin),scale. = F,center = T)

lungrrbstemp <- lungstatdata$RRBS
for(i in 1:dim(lungrrbstemp)[2]){
  ourpid <- colnames(lungrrbstemp)[i]
  ourviallabel <- MotrpacRatTraining6moData::METHYL_META[MotrpacRatTraining6moData::METHYL_META$PID %in% ourpid & MotrpacRatTraining6moData::METHYL_META$Tissue %in% lunglongcode,"viallabel"]
  colnames(lungrrbstemp)[i] <- ourviallabel
}
lungrrbsfin <- lungrrbstemp[,!(colnames(lungrrbstemp) %in% MotrpacRatTraining6moData::OUTLIERS)]
lungrrbspca_res <- prcomp(t(lungrrbsfin),scale. = F,center = T)

lungprot_prtemp <- lungstatdata$Prot_PR
for(i in 1:dim(lungprot_prtemp)[2]){
  ourpid <- colnames(lungprot_prtemp)[i]
  ourviallabel <- intersect(MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid & MotrpacRatTraining6moData::PHENO$tissue %in% lungshortcode,"viallabel"],MotrpacRatTraining6moData::PROT_META$viallabel)
  colnames(lungprot_prtemp)[i] <- ourviallabel
}
lungprot_prfin <- lungprot_prtemp[,!(colnames(lungprot_prtemp) %in% MotrpacRatTraining6moData::OUTLIERS)]
lungprot_prpca_res <- prcomp(t(lungprot_prfin),scale. = F,center = T)

png(file = "Supplemental Figure S6F_112624.png",width = 4,height = 4,units = "in",res = 600)
autoplot(lungrnapca_res,data = lungRNAseqdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(lungshortcode," RNAseq PCA Plot",sep = ""))
dev.off()
pdf(file = "Supplemental Figure S6F_112624.pdf",width = 4,height = 4)
autoplot(lungrnapca_res,data = lungRNAseqdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(lungshortcode," RNAseq PCA Plot",sep = ""))
dev.off()

png(file = "Supplemental Figure S7F_112624.png",width = 4,height = 4,units = "in",res = 600)
autoplot(lungatacpca_res,data = lungATACseqdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(lungshortcode," ATACseq PCA Plot",sep = ""))
dev.off()
pdf(file = "Supplemental Figure S7F_112624.pdf",width = 4,height = 4)
autoplot(lungatacpca_res,data = lungATACseqdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(lungshortcode," ATACseq PCA Plot",sep = ""))
dev.off()

png(file = "Supplemental Figure S8F_112624.png",width = 4,height = 4,units = "in",res = 600)
autoplot(lungrrbspca_res,data = lungRRBSdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(lungshortcode," RRBS PCA Plot",sep = ""))
dev.off()
pdf(file = "Supplemental Figure S8F_112624.pdf",width = 4,height = 4)
autoplot(lungrrbspca_res,data = lungRRBSdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(lungshortcode," RRBS PCA Plot",sep = ""))
dev.off()

#png(file = "Supplemental Figure X4E_5724.png",width = 4,height = 4,units = "in",res = 600)
#autoplot(lungprot_prpca_res,data = lungProt_PRdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(lungshortcode," PROT_PR PCA Plot",sep = ""))
#dev.off()
#pdf(file = "Supplemental Figure X4E_5724.pdf",width = 4,height = 4)
#autoplot(lungprot_prpca_res,data = lungProt_PRdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(lungshortcode," PROT_PR PCA Plot",sep = ""))
#dev.off()

#png(file = "Supplemental Figure X5E_5724.png",width = 4,height = 4,units = "in",res = 600)
#autoplot(lungprot_phpca_res,data = lungProt_PHdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(lungshortcode," PROT_PH PCA Plot",sep = ""))
#dev.off()
#pdf(file = "Supplemental Figure X5E_5724.pdf",width = 4,height = 4)
#autoplot(lungprot_phpca_res,data = lungProt_PHdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(lungshortcode," PROT_PH PCA Plot",sep = ""))
#dev.off()


####
#white

whiternafin <- whitestatdata$RNAseq[,!(colnames(whitestatdata$RNAseq) %in% MotrpacRatTraining6moData::OUTLIERS)]
whiternapca_res <- prcomp(t(whiternafin),scale. = F,center = T)
whiteatacfin <- whitestatdata$ATACseq[,!(colnames(whitestatdata$ATACseq) %in% MotrpacRatTraining6moData::OUTLIERS)]
whiteatacpca_res <- prcomp(t(whiteatacfin),scale. = F,center = T)
whiteprot_phfin <- whitestatdata$Prot_PH[,!(colnames(whitestatdata$Prot_PH) %in% MotrpacRatTraining6moData::OUTLIERS)]
whiteprot_phpca_res <- prcomp(t(whiteprot_phfin),scale. = F,center = T)

whiterrbstemp <- whitestatdata$RRBS
for(i in 1:dim(whiterrbstemp)[2]){
  ourpid <- colnames(whiterrbstemp)[i]
  ourviallabel <- MotrpacRatTraining6moData::METHYL_META[MotrpacRatTraining6moData::METHYL_META$PID %in% ourpid & MotrpacRatTraining6moData::METHYL_META$Tissue %in% whitelongcode,"viallabel"]
  colnames(whiterrbstemp)[i] <- ourviallabel
}
whiterrbsfin <- whiterrbstemp[,!(colnames(whiterrbstemp) %in% MotrpacRatTraining6moData::OUTLIERS)]
whiterrbspca_res <- prcomp(t(whiterrbsfin),scale. = F,center = T)

whiteprot_prtemp <- whitestatdata$Prot_PR
for(i in 1:dim(whiteprot_prtemp)[2]){
  ourpid <- colnames(whiteprot_prtemp)[i]
  ourviallabel <- intersect(MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$pid %in% ourpid & MotrpacRatTraining6moData::PHENO$tissue %in% whiteshortcode,"viallabel"],MotrpacRatTraining6moData::PROT_META$viallabel)
  colnames(whiteprot_prtemp)[i] <- ourviallabel
}
whiteprot_prfin <- whiteprot_prtemp[,!(colnames(whiteprot_prtemp) %in% MotrpacRatTraining6moData::OUTLIERS)]
whiteprot_prpca_res <- prcomp(t(whiteprot_prfin),scale. = F,center = T)

png(file = "Supplemental Figure S6H_112624.png",width = 4,height = 4,units = "in",res = 600)
autoplot(whiternapca_res,data = whiteRNAseqdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(whiteshortcode," RNAseq PCA Plot",sep = ""))
dev.off()
pdf(file = "Supplemental Figure S6H_112624.pdf",width = 4,height = 4)
autoplot(whiternapca_res,data = whiteRNAseqdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(whiteshortcode," RNAseq PCA Plot",sep = ""))
dev.off()

png(file = "Supplemental Figure S7H_112624.png",width = 4,height = 4,units = "in",res = 600)
autoplot(whiteatacpca_res,data = whiteATACseqdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(whiteshortcode," ATACseq PCA Plot",sep = ""))
dev.off()
pdf(file = "Supplemental Figure S7H_112624.pdf",width = 4,height = 4)
autoplot(whiteatacpca_res,data = whiteATACseqdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(whiteshortcode," ATACseq PCA Plot",sep = ""))
dev.off()

png(file = "Supplemental Figure S8H_112624.png",width = 4,height = 4,units = "in",res = 600)
autoplot(whiterrbspca_res,data = whiteRRBSdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(whiteshortcode," RRBS PCA Plot",sep = ""))
dev.off()
pdf(file = "Supplemental Figure S8H_112624.pdf",width = 4,height = 4)
autoplot(whiterrbspca_res,data = whiteRRBSdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(whiteshortcode," RRBS PCA Plot",sep = ""))
dev.off()

#png(file = "Supplemental Figure X4F_5724.png",width = 4,height = 4,units = "in",res = 600)
#autoplot(whiteprot_prpca_res,data = whiteProt_PRdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(whiteshortcode," PROT_PR PCA Plot",sep = ""))
#dev.off()
#pdf(file = "Supplemental Figure X4F_5724.pdf",width = 4,height = 4)
#autoplot(whiteprot_prpca_res,data = whiteProt_PRdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(whiteshortcode," PROT_PR PCA Plot",sep = ""))
#dev.off()

#png(file = "Supplemental Figure X5F_5724.png",width = 4,height = 4,units = "in",res = 600)
#autoplot(whiteprot_phpca_res,data = whiteProt_PHdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(whiteshortcode," PROT_PH PCA Plot",sep = ""))
#dev.off()
#pdf(file = "Supplemental Figure X5F_5724.pdf",width = 4,height = 4)
#autoplot(whiteprot_phpca_res,data = whiteProt_PHdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(whiteshortcode," PROT_PH PCA Plot",sep = ""))
#dev.off()


####
#hippo

hippornafin <- hippostatdata$RNAseq[,!(colnames(hippostatdata$RNAseq) %in% MotrpacRatTraining6moData::OUTLIERS)]
hippornapca_res <- prcomp(t(hippornafin),scale. = F,center = T)
hippoatacfin <- hippostatdata$ATACseq[,!(colnames(hippostatdata$ATACseq) %in% MotrpacRatTraining6moData::OUTLIERS)]
hippoatacpca_res <- prcomp(t(hippoatacfin),scale. = F,center = T)

hipporrbstemp <- hippostatdata$RRBS
for(i in 1:dim(hipporrbstemp)[2]){
  ourpid <- colnames(hipporrbstemp)[i]
  ourviallabel <- MotrpacRatTraining6moData::METHYL_META[MotrpacRatTraining6moData::METHYL_META$PID %in% ourpid & MotrpacRatTraining6moData::METHYL_META$Tissue %in% hippolongcode,"viallabel"]
  colnames(hipporrbstemp)[i] <- ourviallabel
}
hipporrbsfin <- hipporrbstemp[,!(colnames(hipporrbstemp) %in% MotrpacRatTraining6moData::OUTLIERS)]
hipporrbspca_res <- prcomp(t(hipporrbsfin),scale. = F,center = T)

png(file = "Supplemental Figure S6C_112624.png",width = 4,height = 4,units = "in",res = 600)
autoplot(hippornapca_res,data = hippoRNAseqdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(hipposhortcode," RNAseq PCA Plot",sep = ""))
dev.off()
pdf(file = "Supplemental Figure S6C_112624.pdf",width = 4,height = 4)
autoplot(hippornapca_res,data = hippoRNAseqdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(hipposhortcode," RNAseq PCA Plot",sep = ""))
dev.off()

png(file = "Supplemental Figure S7C_112624.png",width = 4,height = 4,units = "in",res = 600)
autoplot(hippoatacpca_res,data = hippoATACseqdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(hipposhortcode," ATACseq PCA Plot",sep = ""))
dev.off()
pdf(file = "Supplemental Figure S7C_112624.pdf",width = 4,height = 4)
autoplot(hippoatacpca_res,data = hippoATACseqdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(hipposhortcode," ATACseq PCA Plot",sep = ""))
dev.off()

png(file = "Supplemental Figure S8C_112624.png",width = 4,height = 4,units = "in",res = 600)
autoplot(hipporrbspca_res,data = hippoRRBSdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(hipposhortcode," RRBS PCA Plot",sep = ""))
dev.off()
pdf(file = "Supplemental Figure S8C_112624.pdf",width = 4,height = 4)
autoplot(hipporrbspca_res,data = hippoRRBSdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(hipposhortcode," RRBS PCA Plot",sep = ""))
dev.off()


####
#brown

brownrnafin <- brownstatdata$RNAseq[,!(colnames(brownstatdata$RNAseq) %in% MotrpacRatTraining6moData::OUTLIERS)]
brownrnapca_res <- prcomp(t(brownrnafin),scale. = F,center = T)
brownatacfin <- brownstatdata$ATACseq[,!(colnames(brownstatdata$ATACseq) %in% MotrpacRatTraining6moData::OUTLIERS)]
brownatacpca_res <- prcomp(t(brownatacfin),scale. = F,center = T)

brownrrbstemp <- brownstatdata$RRBS
for(i in 1:dim(brownrrbstemp)[2]){
  ourpid <- colnames(brownrrbstemp)[i]
  ourviallabel <- MotrpacRatTraining6moData::METHYL_META[MotrpacRatTraining6moData::METHYL_META$PID %in% ourpid & MotrpacRatTraining6moData::METHYL_META$Tissue %in% brownlongcode,"viallabel"]
  colnames(brownrrbstemp)[i] <- ourviallabel
}
brownrrbsfin <- brownrrbstemp[,!(colnames(brownrrbstemp) %in% MotrpacRatTraining6moData::OUTLIERS)]
brownrrbspca_res <- prcomp(t(brownrrbsfin),scale. = F,center = T)

png(file = "Supplemental Figure S6G_112624.png",width = 4,height = 4,units = "in",res = 600)
autoplot(brownrnapca_res,data = brownRNAseqdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(brownshortcode," RNAseq PCA Plot",sep = ""))
dev.off()
pdf(file = "Supplemental Figure S6G_112624.pdf",width = 4,height = 4)
autoplot(brownrnapca_res,data = brownRNAseqdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(brownshortcode," RNAseq PCA Plot",sep = ""))
dev.off()

png(file = "Supplemental Figure S7G_112624.png",width = 4,height = 4,units = "in",res = 600)
autoplot(brownatacpca_res,data = brownATACseqdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(brownshortcode," ATACseq PCA Plot",sep = ""))
dev.off()
pdf(file = "Supplemental Figure S7G_112624.pdf",width = 4,height = 4)
autoplot(brownatacpca_res,data = brownATACseqdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(brownshortcode," ATACseq PCA Plot",sep = ""))
dev.off()

png(file = "Supplemental Figure S8G_112624.png",width = 4,height = 4,units = "in",res = 600)
autoplot(brownrrbspca_res,data = brownRRBSdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(brownshortcode," RRBS PCA Plot",sep = ""))
dev.off()
pdf(file = "Supplemental Figure S8G_112624.pdf",width = 4,height = 4)
autoplot(brownrrbspca_res,data = brownRRBSdesign,colour = "Group",shape = "Sex",size = 3) + theme_classic() + scale_colour_manual(values = ann_cols$Group) + ggtitle(paste(brownshortcode," RRBS PCA Plot",sep = ""))
dev.off()




gastrostatResultsaltRNAEQ = MultiPower(data = list("RNAseq" = gastrostatdata$RNAseq), groups = list("RNAseq" = gastrostatdesign$RNAseq), type = 2,omicPower = 0.25,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = TRUE)
heartstatResultsaltRNAEQ = MultiPower(data = list("RNAseq" = gastrostatdata$RNAseq), groups = list("RNAseq" = gastrostatdesign$RNAseq), type = 2,omicPower = 0.25,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = TRUE)
hippostatResultsaltRNAEQ = MultiPower(data = list("RNAseq" = gastrostatdata$RNAseq), groups = list("RNAseq" = gastrostatdesign$RNAseq), type = 2,omicPower = 0.25,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = TRUE)
kidneystatResultsaltRNAEQ = MultiPower(data = list("RNAseq" = gastrostatdata$RNAseq), groups = list("RNAseq" = gastrostatdesign$RNAseq), type = 2,omicPower = 0.25,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = TRUE)
liverstatResultsaltRNAEQ = MultiPower(data = list("RNAseq" = gastrostatdata$RNAseq), groups = list("RNAseq" = gastrostatdesign$RNAseq), type = 2,omicPower = 0.25,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = TRUE)
lungstatResultsaltRNAEQ = MultiPower(data = list("RNAseq" = gastrostatdata$RNAseq), groups = list("RNAseq" = gastrostatdesign$RNAseq), type = 2,omicPower = 0.25,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = TRUE)
brownstatResultsaltRNAEQ = MultiPower(data = list("RNAseq" = gastrostatdata$RNAseq), groups = list("RNAseq" = gastrostatdesign$RNAseq), type = 2,omicPower = 0.25,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = TRUE)
whitestatResultsaltRNAEQ = MultiPower(data = list("RNAseq" = gastrostatdata$RNAseq), groups = list("RNAseq" = gastrostatdesign$RNAseq), type = 2,omicPower = 0.25,averagePower = 0.8,fdr = 0.05, cost = 1, equalSize = TRUE,max.size = 2000, omicCol = miscolores,powerPlots = TRUE)



gastroRNApower <- data.frame(row.names = names(gastrostatResultsRNAEQ$parameters$RNAseq$d),
                             "Cohens.d" = gastrostatResultsRNAEQ$parameters$RNAseq$d,
                             "Power" = rep(0,length(names(gastrostatResultsRNAEQ$parameters$RNAseq$d))))
for(i in 1:length(names(gastrostatResultsRNAEQ$parameters$RNAseq$d))){
  ouroutput <- pwr.t2n.test(n1 = 10,n2 = 40,d = gastroRNApower[i,"Cohens.d"],alternative = "two.sided")
  gastroRNApower[i,"Power"] <- ouroutput$power
}

gastroATACpower <- data.frame(row.names = names(gastrostatResultsATACEQ$parameters$ATACseq$d),
                             "Cohens.d" = gastrostatResultsATACEQ$parameters$ATACseq$d,
                             "Power" = rep(0,length(names(gastrostatResultsATACEQ$parameters$ATACseq$d))))
for(i in 1:length(names(gastrostatResultsATACEQ$parameters$ATACseq$d))){
  if(i%%10000 == 0){
    print(i)
  }
  ouroutput <- pwr.t2n.test(n1 = 10,n2 = 40,d = gastroATACpower[i,"Cohens.d"],alternative = "two.sided")
  gastroATACpower[i,"Power"] <- ouroutput$power
}

gastroRRBSpower <- data.frame(row.names = names(gastrostatResultsRRBSEQ$parameters$RRBS$d),
                              "Cohens.d" = gastrostatResultsRRBSEQ$parameters$RRBS$d,
                              "Power" = rep(0,length(names(gastrostatResultsRRBSEQ$parameters$RRBS$d))))
for(i in 1:length(names(gastrostatResultsRRBSEQ$parameters$RRBS$d))){
  if(i%%10000 == 0){
    print(i)
  }
  ouroutput <- pwr.t2n.test(n1 = 10,n2 = 40,d = gastroRRBSpower[i,"Cohens.d"],alternative = "two.sided")
  gastroRRBSpower[i,"Power"] <- ouroutput$power
}



heartRNApower <- data.frame(row.names = names(heartstatResultsRNAEQ$parameters$RNAseq$d),
                             "Cohens.d" = heartstatResultsRNAEQ$parameters$RNAseq$d,
                             "Power" = rep(0,length(names(heartstatResultsRNAEQ$parameters$RNAseq$d))))
for(i in 1:length(names(heartstatResultsRNAEQ$parameters$RNAseq$d))){
  ouroutput <- pwr.t2n.test(n1 = 10,n2 = 40,d = heartRNApower[i,"Cohens.d"],alternative = "two.sided")
  heartRNApower[i,"Power"] <- ouroutput$power
}


hippoRNApower <- data.frame(row.names = names(hippostatResultsRNAEQ$parameters$RNAseq$d),
                             "Cohens.d" = hippostatResultsRNAEQ$parameters$RNAseq$d,
                             "Power" = rep(0,length(names(hippostatResultsRNAEQ$parameters$RNAseq$d))))
for(i in 1:length(names(hippostatResultsRNAEQ$parameters$RNAseq$d))){
  ouroutput <- pwr.t2n.test(n1 = 10,n2 = 40,d = hippoRNApower[i,"Cohens.d"],alternative = "two.sided")
  hippoRNApower[i,"Power"] <- ouroutput$power
}


kidneyRNApower <- data.frame(row.names = names(kidneystatResultsRNAEQ$parameters$RNAseq$d),
                             "Cohens.d" = kidneystatResultsRNAEQ$parameters$RNAseq$d,
                             "Power" = rep(0,length(names(kidneystatResultsRNAEQ$parameters$RNAseq$d))))
for(i in 1:length(names(kidneystatResultsRNAEQ$parameters$RNAseq$d))){
  ouroutput <- pwr.t2n.test(n1 = 10,n2 = 40,d = kidneyRNApower[i,"Cohens.d"],alternative = "two.sided")
  kidneyRNApower[i,"Power"] <- ouroutput$power
}


liverRNApower <- data.frame(row.names = names(liverstatResultsRNAEQ$parameters$RNAseq$d),
                             "Cohens.d" = liverstatResultsRNAEQ$parameters$RNAseq$d,
                             "Power" = rep(0,length(names(liverstatResultsRNAEQ$parameters$RNAseq$d))))
for(i in 1:length(names(liverstatResultsRNAEQ$parameters$RNAseq$d))){
  ouroutput <- pwr.t2n.test(n1 = 10,n2 = 40,d = liverRNApower[i,"Cohens.d"],alternative = "two.sided")
  liverRNApower[i,"Power"] <- ouroutput$power
}


lungRNApower <- data.frame(row.names = names(lungstatResultsRNAEQ$parameters$RNAseq$d),
                             "Cohens.d" = lungstatResultsRNAEQ$parameters$RNAseq$d,
                             "Power" = rep(0,length(names(lungstatResultsRNAEQ$parameters$RNAseq$d))))
for(i in 1:length(names(lungstatResultsRNAEQ$parameters$RNAseq$d))){
  ouroutput <- pwr.t2n.test(n1 = 10,n2 = 40,d = lungRNApower[i,"Cohens.d"],alternative = "two.sided")
  lungRNApower[i,"Power"] <- ouroutput$power
}


brownRNApower <- data.frame(row.names = names(brownstatResultsRNAEQ$parameters$RNAseq$d),
                             "Cohens.d" = brownstatResultsRNAEQ$parameters$RNAseq$d,
                             "Power" = rep(0,length(names(brownstatResultsRNAEQ$parameters$RNAseq$d))))
for(i in 1:length(names(brownstatResultsRNAEQ$parameters$RNAseq$d))){
  ouroutput <- pwr.t2n.test(n1 = 10,n2 = 40,d = brownRNApower[i,"Cohens.d"],alternative = "two.sided")
  brownRNApower[i,"Power"] <- ouroutput$power
}


whiteRNApower <- data.frame(row.names = names(whitestatResultsRNAEQ$parameters$RNAseq$d),
                             "Cohens.d" = whitestatResultsRNAEQ$parameters$RNAseq$d,
                             "Power" = rep(0,length(names(whitestatResultsRNAEQ$parameters$RNAseq$d))))
for(i in 1:length(names(whitestatResultsRNAEQ$parameters$RNAseq$d))){
  ouroutput <- pwr.t2n.test(n1 = 10,n2 = 40,d = whiteRNApower[i,"Cohens.d"],alternative = "two.sided")
  whiteRNApower[i,"Power"] <- ouroutput$power
}

save.image("outputpoweranalysis.RData")


heartATACpower <- data.frame(row.names = names(heartstatResultsATACEQ$parameters$ATACseq$d),
                              "Cohens.d" = heartstatResultsATACEQ$parameters$ATACseq$d,
                              "Power" = rep(0,length(names(heartstatResultsATACEQ$parameters$ATACseq$d))))
for(i in 1:length(names(heartstatResultsATACEQ$parameters$ATACseq$d))){
  if(i%%10000 == 0){
    print(i)
  }
  ouroutput <- pwr.t2n.test(n1 = 10,n2 = 40,d = heartATACpower[i,"Cohens.d"],alternative = "two.sided")
  heartATACpower[i,"Power"] <- ouroutput$power
}

heartRRBSpower <- data.frame(row.names = names(heartstatResultsRRBSEQ$parameters$RRBS$d),
                              "Cohens.d" = heartstatResultsRRBSEQ$parameters$RRBS$d,
                              "Power" = rep(0,length(names(heartstatResultsRRBSEQ$parameters$RRBS$d))))
for(i in 1:length(names(heartstatResultsRRBSEQ$parameters$RRBS$d))){
  if(i%%10000 == 0){
    print(i)
  }
  ouroutput <- pwr.t2n.test(n1 = 10,n2 = 40,d = heartRRBSpower[i,"Cohens.d"],alternative = "two.sided")
  heartRRBSpower[i,"Power"] <- ouroutput$power
}

save.image("outputpoweranalysis.RData")

hippoATACpower <- data.frame(row.names = names(hippostatResultsATACEQ$parameters$ATACseq$d),
                              "Cohens.d" = hippostatResultsATACEQ$parameters$ATACseq$d,
                              "Power" = rep(0,length(names(hippostatResultsATACEQ$parameters$ATACseq$d))))
for(i in 1:length(names(hippostatResultsATACEQ$parameters$ATACseq$d))){
  if(i%%10000 == 0){
    print(i)
  }
  ouroutput <- pwr.t2n.test(n1 = 10,n2 = 40,d = hippoATACpower[i,"Cohens.d"],alternative = "two.sided")
  hippoATACpower[i,"Power"] <- ouroutput$power
}

hippoRRBSpower <- data.frame(row.names = names(hippostatResultsRRBSEQ$parameters$RRBS$d),
                              "Cohens.d" = hippostatResultsRRBSEQ$parameters$RRBS$d,
                              "Power" = rep(0,length(names(hippostatResultsRRBSEQ$parameters$RRBS$d))))
for(i in 1:length(names(hippostatResultsRRBSEQ$parameters$RRBS$d))){
  if(i%%10000 == 0){
    print(i)
  }
  ouroutput <- pwr.t2n.test(n1 = 10,n2 = 40,d = hippoRRBSpower[i,"Cohens.d"],alternative = "two.sided")
  hippoRRBSpower[i,"Power"] <- ouroutput$power
}

save.image("outputpoweranalysis.RData")


kidneyATACpower <- data.frame(row.names = names(kidneystatResultsATACEQ$parameters$ATACseq$d),
                              "Cohens.d" = kidneystatResultsATACEQ$parameters$ATACseq$d,
                              "Power" = rep(0,length(names(kidneystatResultsATACEQ$parameters$ATACseq$d))))
for(i in 1:length(names(kidneystatResultsATACEQ$parameters$ATACseq$d))){
  if(i%%10000 == 0){
    print(i)
  }
  ouroutput <- pwr.t2n.test(n1 = 10,n2 = 40,d = kidneyATACpower[i,"Cohens.d"],alternative = "two.sided")
  kidneyATACpower[i,"Power"] <- ouroutput$power
}

kidneyRRBSpower <- data.frame(row.names = names(kidneystatResultsRRBSEQ$parameters$RRBS$d),
                              "Cohens.d" = kidneystatResultsRRBSEQ$parameters$RRBS$d,
                              "Power" = rep(0,length(names(kidneystatResultsRRBSEQ$parameters$RRBS$d))))
for(i in 1:length(names(kidneystatResultsRRBSEQ$parameters$RRBS$d))){
  if(i%%10000 == 0){
    print(i)
  }
  ouroutput <- pwr.t2n.test(n1 = 10,n2 = 40,d = kidneyRRBSpower[i,"Cohens.d"],alternative = "two.sided")
  kidneyRRBSpower[i,"Power"] <- ouroutput$power
}

save.image("outputpoweranalysis.RData")

liverATACpower <- data.frame(row.names = names(liverstatResultsATACEQ$parameters$ATACseq$d),
                              "Cohens.d" = liverstatResultsATACEQ$parameters$ATACseq$d,
                              "Power" = rep(0,length(names(liverstatResultsATACEQ$parameters$ATACseq$d))))
for(i in 1:length(names(liverstatResultsATACEQ$parameters$ATACseq$d))){
  if(i%%10000 == 0){
    print(i)
  }
  ouroutput <- pwr.t2n.test(n1 = 10,n2 = 40,d = liverATACpower[i,"Cohens.d"],alternative = "two.sided")
  liverATACpower[i,"Power"] <- ouroutput$power
}

liverRRBSpower <- data.frame(row.names = names(liverstatResultsRRBSEQ$parameters$RRBS$d),
                              "Cohens.d" = liverstatResultsRRBSEQ$parameters$RRBS$d,
                              "Power" = rep(0,length(names(liverstatResultsRRBSEQ$parameters$RRBS$d))))
for(i in 1:length(names(liverstatResultsRRBSEQ$parameters$RRBS$d))){
  if(i%%10000 == 0){
    print(i)
  }
  ouroutput <- pwr.t2n.test(n1 = 10,n2 = 40,d = liverRRBSpower[i,"Cohens.d"],alternative = "two.sided")
  liverRRBSpower[i,"Power"] <- ouroutput$power
}

save.image("outputpoweranalysis.RData")

lungATACpower <- data.frame(row.names = names(lungstatResultsATACEQ$parameters$ATACseq$d),
                              "Cohens.d" = lungstatResultsATACEQ$parameters$ATACseq$d,
                              "Power" = rep(0,length(names(lungstatResultsATACEQ$parameters$ATACseq$d))))
for(i in 1:length(names(lungstatResultsATACEQ$parameters$ATACseq$d))){
  if(i%%10000 == 0){
    print(i)
  }
  ouroutput <- pwr.t2n.test(n1 = 10,n2 = 40,d = lungATACpower[i,"Cohens.d"],alternative = "two.sided")
  lungATACpower[i,"Power"] <- ouroutput$power
}

lungRRBSpower <- data.frame(row.names = names(lungstatResultsRRBSEQ$parameters$RRBS$d),
                              "Cohens.d" = lungstatResultsRRBSEQ$parameters$RRBS$d,
                              "Power" = rep(0,length(names(lungstatResultsRRBSEQ$parameters$RRBS$d))))
for(i in 1:length(names(lungstatResultsRRBSEQ$parameters$RRBS$d))){
  if(i%%10000 == 0){
    print(i)
  }
  ouroutput <- pwr.t2n.test(n1 = 10,n2 = 40,d = lungRRBSpower[i,"Cohens.d"],alternative = "two.sided")
  lungRRBSpower[i,"Power"] <- ouroutput$power
}

save.image("outputpoweranalysis.RData")

brownATACpower <- data.frame(row.names = names(brownstatResultsATACEQ$parameters$ATACseq$d),
                              "Cohens.d" = brownstatResultsATACEQ$parameters$ATACseq$d,
                              "Power" = rep(0,length(names(brownstatResultsATACEQ$parameters$ATACseq$d))))
for(i in 1:length(names(brownstatResultsATACEQ$parameters$ATACseq$d))){
  if(i%%10000 == 0){
    print(i)
  }
  ouroutput <- pwr.t2n.test(n1 = 10,n2 = 40,d = brownATACpower[i,"Cohens.d"],alternative = "two.sided")
  brownATACpower[i,"Power"] <- ouroutput$power
}

brownRRBSpower <- data.frame(row.names = names(brownstatResultsRRBSEQ$parameters$RRBS$d),
                              "Cohens.d" = brownstatResultsRRBSEQ$parameters$RRBS$d,
                              "Power" = rep(0,length(names(brownstatResultsRRBSEQ$parameters$RRBS$d))))
for(i in 1:length(names(brownstatResultsRRBSEQ$parameters$RRBS$d))){
  if(i%%10000 == 0){
    print(i)
  }
  ouroutput <- pwr.t2n.test(n1 = 10,n2 = 40,d = brownRRBSpower[i,"Cohens.d"],alternative = "two.sided")
  brownRRBSpower[i,"Power"] <- ouroutput$power
}

save.image("outputpoweranalysis.RData")

whiteATACpower <- data.frame(row.names = names(whitestatResultsATACEQ$parameters$ATACseq$d),
                              "Cohens.d" = whitestatResultsATACEQ$parameters$ATACseq$d,
                              "Power" = rep(0,length(names(whitestatResultsATACEQ$parameters$ATACseq$d))))
for(i in 1:length(names(whitestatResultsATACEQ$parameters$ATACseq$d))){
  if(i%%10000 == 0){
    print(i)
  }
  ouroutput <- pwr.t2n.test(n1 = 10,n2 = 40,d = whiteATACpower[i,"Cohens.d"],alternative = "two.sided")
  whiteATACpower[i,"Power"] <- ouroutput$power
}

whiteRRBSpower <- data.frame(row.names = names(whitestatResultsRRBSEQ$parameters$RRBS$d),
                              "Cohens.d" = whitestatResultsRRBSEQ$parameters$RRBS$d,
                              "Power" = rep(0,length(names(whitestatResultsRRBSEQ$parameters$RRBS$d))))
for(i in 1:length(names(whitestatResultsRRBSEQ$parameters$RRBS$d))){
  if(i%%10000 == 0){
    print(i)
  }
  ouroutput <- pwr.t2n.test(n1 = 10,n2 = 40,d = whiteRRBSpower[i,"Cohens.d"],alternative = "two.sided")
  whiteRRBSpower[i,"Power"] <- ouroutput$power
}

save.image("outputpoweranalysis.RData")


load("omesigdata.RData")
library(ggplot2)

gastroRNApower$Signif <- "Non.DEG"
gastroRNApower[gastrornasig,"Signif"] <- "DEG"
heartRNApower$Signif <- "Non.DEG"
heartRNApower[heartrnasig,"Signif"] <- "DEG"
hippoRNApower$Signif <- "Non.DEG"
hippoRNApower[hippornasig,"Signif"] <- "DEG"
kidneyRNApower$Signif <- "Non.DEG"
kidneyRNApower[kidneyrnasig,"Signif"] <- "DEG"
liverRNApower$Signif <- "Non.DEG"
liverRNApower[liverrnasig,"Signif"] <- "DEG"
lungRNApower$Signif <- "Non.DEG"
lungRNApower[lungrnasig,"Signif"] <- "DEG"
brownRNApower$Signif <- "Non.DEG"
brownRNApower[brownrnasig,"Signif"] <- "DEG"
whiteRNApower$Signif <- "Non.DEG"
whiteRNApower[whiternasig,"Signif"] <- "DEG"

gastroATACpower$Signif <- "Non.DAR"
gastroATACpower[gastroatacsig,"Signif"] <- "DAR"
heartATACpower$Signif <- "Non.DAR"
heartATACpower[heartatacsig,"Signif"] <- "DAR"
hippoATACpower$Signif <- "Non.DAR"
hippoATACpower[hippoatacsig,"Signif"] <- "DAR"
kidneyATACpower$Signif <- "Non.DAR"
kidneyATACpower[kidneyatacsig,"Signif"] <- "DAR"
liverATACpower$Signif <- "Non.DAR"
liverATACpower[liveratacsig,"Signif"] <- "DAR"
lungATACpower$Signif <- "Non.DAR"
lungATACpower[lungatacsig,"Signif"] <- "DAR"
brownATACpower$Signif <- "Non.DAR"
brownATACpower[brownatacsig,"Signif"] <- "DAR"
whiteATACpower$Signif <- "Non.DAR"
whiteATACpower[whiteatacsig,"Signif"] <- "DAR"

gastroRRBSpower$Signif <- "Non.DMR"
gastroRRBSpower[gastromethsig,"Signif"] <- "DMR"
heartRRBSpower$Signif <- "Non.DMR"
heartRRBSpower[heartmethsig,"Signif"] <- "DMR"
hippoRRBSpower$Signif <- "Non.DMR"
hippoRRBSpower[hippomethsig,"Signif"] <- "DMR"
kidneyRRBSpower$Signif <- "Non.DMR"
kidneyRRBSpower[kidneymethsig,"Signif"] <- "DMR"
liverRRBSpower$Signif <- "Non.DMR"
liverRRBSpower[livermethsig,"Signif"] <- "DMR"
lungRRBSpower$Signif <- "Non.DMR"
lungRRBSpower[lungmethsig,"Signif"] <- "DMR"
brownRRBSpower$Signif <- "Non.DMR"
brownRRBSpower[brownmethsig,"Signif"] <- "DMR"
whiteRRBSpower$Signif <- "Non.DMR"
whiteRRBSpower[whitemethsig,"Signif"] <- "DMR"


png(file = "gastro_rna_power_hist.png",width = 5,height = 5,units = "in",res = 600)
hist(gastroRNApower$Power,main = "All Gene Power",xlab = "Power")
dev.off()

png(file = "gastro_rna_sig_power_hist.png",width = 5,height = 5,units = "in",res = 600)
hist(gastroRNApower[gastrornasig,"Power"],main = "DEG Power",xlab = "Power")
dev.off()

png(file = "heart_rna_power_hist.png",width = 5,height = 5,units = "in",res = 600)
hist(heartRNApower$Power,main = "All Gene Power",xlab = "Power")
dev.off()

png(file = "heart_rna_sig_power_hist.png",width = 5,height = 5,units = "in",res = 600)
hist(heartRNApower[heartrnasig,"Power"],main = "DEG Power",xlab = "Power")
dev.off()


png(file = "hippo_rna_power_hist.png",width = 5,height = 5,units = "in",res = 600)
hist(hippoRNApower$Power,main = "All Gene Power",xlab = "Power")
dev.off()

png(file = "hippo_rna_sig_power_hist.png",width = 5,height = 5,units = "in",res = 600)
hist(hippoRNApower[hippornasig,"Power"],main = "DEG Power",xlab = "Power")
dev.off()


png(file = "kidney_rna_power_hist.png",width = 5,height = 5,units = "in",res = 600)
hist(kidneyRNApower$Power,main = "All Gene Power",xlab = "Power")
dev.off()

png(file = "kidney_rna_sig_power_hist.png",width = 5,height = 5,units = "in",res = 600)
hist(kidneyRNApower[kidneyrnasig,"Power"],main = "DEG Power",xlab = "Power")
dev.off()


png(file = "liver_rna_power_hist.png",width = 5,height = 5,units = "in",res = 600)
hist(liverRNApower$Power,main = "All Gene Power",xlab = "Power")
dev.off()

png(file = "liver_rna_sig_power_hist.png",width = 5,height = 5,units = "in",res = 600)
hist(liverRNApower[liverrnasig,"Power"],main = "DEG Power",xlab = "Power")
dev.off()


png(file = "lung_rna_power_hist.png",width = 5,height = 5,units = "in",res = 600)
hist(lungRNApower$Power,main = "All Gene Power",xlab = "Power")
dev.off()

png(file = "lung_rna_sig_power_hist.png",width = 5,height = 5,units = "in",res = 600)
hist(lungRNApower[lungrnasig,"Power"],main = "DEG Power",xlab = "Power")
dev.off()


png(file = "brown_rna_power_hist.png",width = 5,height = 5,units = "in",res = 600)
hist(brownRNApower$Power,main = "All Gene Power",xlab = "Power")
dev.off()

png(file = "brown_rna_sig_power_hist.png",width = 5,height = 5,units = "in",res = 600)
hist(brownRNApower[brownrnasig,"Power"],main = "DEG Power",xlab = "Power")
dev.off()


png(file = "white_rna_power_hist.png",width = 5,height = 5,units = "in",res = 600)
hist(whiteRNApower$Power,main = "All Gene Power",xlab = "Power")
dev.off()

png(file = "white_rna_sig_power_hist.png",width = 5,height = 5,units = "in",res = 600)
hist(whiteRNApower[whiternasig,"Power"],main = "DEG Power",xlab = "Power")
dev.off()


gastroRNApowermod <- gastroRNApower
gastroRNApowermod$feature <- rownames(gastroRNApowermod)
rownames(gastroRNApowermod) <- paste("RNA_","gastro_",rownames(gastroRNApowermod),sep = "")

heartRNApowermod <- heartRNApower
heartRNApowermod$feature <- rownames(heartRNApowermod)
rownames(heartRNApowermod) <- paste("RNA_","heart_",rownames(heartRNApowermod),sep = "")

hippoRNApowermod <- hippoRNApower
hippoRNApowermod$feature <- rownames(hippoRNApowermod)
rownames(hippoRNApowermod) <- paste("RNA_","hippo_",rownames(hippoRNApowermod),sep = "")

kidneyRNApowermod <- kidneyRNApower
kidneyRNApowermod$feature <- rownames(kidneyRNApowermod)
rownames(kidneyRNApowermod) <- paste("RNA_","kidney_",rownames(kidneyRNApowermod),sep = "")

liverRNApowermod <- liverRNApower
liverRNApowermod$feature <- rownames(liverRNApowermod)
rownames(liverRNApowermod) <- paste("RNA_","liver_",rownames(liverRNApowermod),sep = "")

lungRNApowermod <- lungRNApower
lungRNApowermod$feature <- rownames(lungRNApowermod)
rownames(lungRNApowermod) <- paste("RNA_","lung_",rownames(lungRNApowermod),sep = "")

brownRNApowermod <- brownRNApower
brownRNApowermod$feature <- rownames(brownRNApowermod)
rownames(brownRNApowermod) <- paste("RNA_","brown_",rownames(brownRNApowermod),sep = "")

whiteRNApowermod <- whiteRNApower
whiteRNApowermod$feature <- rownames(whiteRNApowermod)
rownames(whiteRNApowermod) <- paste("RNA_","white_",rownames(whiteRNApowermod),sep = "")


gastroATACpowermod <- gastroATACpower
gastroATACpowermod$feature <- rownames(gastroATACpowermod)
rownames(gastroATACpowermod) <- paste("ATAC_","gastro_",rownames(gastroATACpowermod),sep = "")

heartATACpowermod <- heartATACpower
heartATACpowermod$feature <- rownames(heartATACpowermod)
rownames(heartATACpowermod) <- paste("ATAC_","heart_",rownames(heartATACpowermod),sep = "")

hippoATACpowermod <- hippoATACpower
hippoATACpowermod$feature <- rownames(hippoATACpowermod)
rownames(hippoATACpowermod) <- paste("ATAC_","hippo_",rownames(hippoATACpowermod),sep = "")

kidneyATACpowermod <- kidneyATACpower
kidneyATACpowermod$feature <- rownames(kidneyATACpowermod)
rownames(kidneyATACpowermod) <- paste("ATAC_","kidney_",rownames(kidneyATACpowermod),sep = "")

liverATACpowermod <- liverATACpower
liverATACpowermod$feature <- rownames(liverATACpowermod)
rownames(liverATACpowermod) <- paste("ATAC_","liver_",rownames(liverATACpowermod),sep = "")

lungATACpowermod <- lungATACpower
lungATACpowermod$feature <- rownames(lungATACpowermod)
rownames(lungATACpowermod) <- paste("ATAC_","lung_",rownames(lungATACpowermod),sep = "")

brownATACpowermod <- brownATACpower
brownATACpowermod$feature <- rownames(brownATACpowermod)
rownames(brownATACpowermod) <- paste("ATAC_","brown_",rownames(brownATACpowermod),sep = "")

whiteATACpowermod <- whiteATACpower
whiteATACpowermod$feature <- rownames(whiteATACpowermod)
rownames(whiteATACpowermod) <- paste("ATAC_","white_",rownames(whiteATACpowermod),sep = "")



gastroRRBSpowermod <- gastroRRBSpower
gastroRRBSpowermod$feature <- rownames(gastroRRBSpowermod)
rownames(gastroRRBSpowermod) <- paste("RRBS_","gastro_",rownames(gastroRRBSpowermod),sep = "")

heartRRBSpowermod <- heartRRBSpower
heartRRBSpowermod$feature <- rownames(heartRRBSpowermod)
rownames(heartRRBSpowermod) <- paste("RRBS_","heart_",rownames(heartRRBSpowermod),sep = "")

hippoRRBSpowermod <- hippoRRBSpower
hippoRRBSpowermod$feature <- rownames(hippoRRBSpowermod)
rownames(hippoRRBSpowermod) <- paste("RRBS_","hippo_",rownames(hippoRRBSpowermod),sep = "")

kidneyRRBSpowermod <- kidneyRRBSpower
kidneyRRBSpowermod$feature <- rownames(kidneyRRBSpowermod)
rownames(kidneyRRBSpowermod) <- paste("RRBS_","kidney_",rownames(kidneyRRBSpowermod),sep = "")

liverRRBSpowermod <- liverRRBSpower
liverRRBSpowermod$feature <- rownames(liverRRBSpowermod)
rownames(liverRRBSpowermod) <- paste("RRBS_","liver_",rownames(liverRRBSpowermod),sep = "")

lungRRBSpowermod <- lungRRBSpower
lungRRBSpowermod$feature <- rownames(lungRRBSpowermod)
rownames(lungRRBSpowermod) <- paste("RRBS_","lung_",rownames(lungRRBSpowermod),sep = "")

brownRRBSpowermod <- brownRRBSpower
brownRRBSpowermod$feature <- rownames(brownRRBSpowermod)
rownames(brownRRBSpowermod) <- paste("RRBS_","brown_",rownames(brownRRBSpowermod),sep = "")

whiteRRBSpowermod <- whiteRRBSpower
whiteRRBSpowermod$feature <- rownames(whiteRRBSpowermod)
rownames(whiteRRBSpowermod) <- paste("RRBS_","white_",rownames(whiteRRBSpowermod),sep = "")

gastroRNApowermod$Tissue <- "SKM-GN"
heartRNApowermod$Tissue <- "HEART"
hippoRNApowermod$Tissue <- "HIPPOC"
kidneyRNApowermod$Tissue <- "KIDNEY"
liverRNApowermod$Tissue <- "LIVER"
lungRNApowermod$Tissue <- "LUNG"
brownRNApowermod$Tissue <- "BAT"
whiteRNApowermod$Tissue <- "WAT-SC"

gastroATACpowermod$Tissue <- "SKM-GN"
heartATACpowermod$Tissue <- "HEART"
hippoATACpowermod$Tissue <- "HIPPOC"
kidneyATACpowermod$Tissue <- "KIDNEY"
liverATACpowermod$Tissue <- "LIVER"
lungATACpowermod$Tissue <- "LUNG"
brownATACpowermod$Tissue <- "BAT"
whiteATACpowermod$Tissue <- "WAT-SC"

gastroRRBSpowermod$Tissue <- "SKM-GN"
heartRRBSpowermod$Tissue <- "HEART"
hippoRRBSpowermod$Tissue <- "HIPPOC"
kidneyRRBSpowermod$Tissue <- "KIDNEY"
liverRRBSpowermod$Tissue <- "LIVER"
lungRRBSpowermod$Tissue <- "LUNG"
brownRRBSpowermod$Tissue <- "BAT"
whiteRRBSpowermod$Tissue <- "WAT-SC"

gastroRNApowermod$Ome <- "RNA"
heartRNApowermod$Ome <- "RNA"
hippoRNApowermod$Ome <- "RNA"
kidneyRNApowermod$Ome <- "RNA"
liverRNApowermod$Ome <- "RNA"
lungRNApowermod$Ome <- "RNA"
brownRNApowermod$Ome <- "RNA"
whiteRNApowermod$Ome <- "RNA"

gastroATACpowermod$Ome <- "ATAC"
heartATACpowermod$Ome <- "ATAC"
hippoATACpowermod$Ome <- "ATAC"
kidneyATACpowermod$Ome <- "ATAC"
liverATACpowermod$Ome <- "ATAC"
lungATACpowermod$Ome <- "ATAC"
brownATACpowermod$Ome <- "ATAC"
whiteATACpowermod$Ome <- "ATAC"

gastroRRBSpowermod$Ome <- "RRBS"
heartRRBSpowermod$Ome <- "RRBS"
hippoRRBSpowermod$Ome <- "RRBS"
kidneyRRBSpowermod$Ome <- "RRBS"
liverRRBSpowermod$Ome <- "RRBS"
lungRRBSpowermod$Ome <- "RRBS"
brownRRBSpowermod$Ome <- "RRBS"
whiteRRBSpowermod$Ome <- "RRBS"

allpowermod <- rbind(gastroRNApowermod,
                     heartRNApowermod,
                     hippoRNApowermod,
                     kidneyRNApowermod,
                     liverRNApowermod,
                     lungRNApowermod,
                     brownRNApowermod,
                     whiteRNApowermod,
                     gastroATACpowermod,
                     heartATACpowermod,
                     hippoATACpowermod,
                     kidneyATACpowermod,
                     liverATACpowermod,
                     lungATACpowermod,
                     brownATACpowermod,
                     whiteATACpowermod,
                     gastroRRBSpowermod,
                     heartRRBSpowermod,
                     hippoRRBSpowermod,
                     kidneyRRBSpowermod,
                     liverRRBSpowermod,
                     lungRRBSpowermod,
                     brownRRBSpowermod,
                     whiteRRBSpowermod)


png(file = "All Analyte Effect Size by Tissue and Ome.png",width = 7,height = 5,units = "in",res = 600)
ggplot(allpowermod,aes(x = Ome,y = Cohens.d,fill = Tissue)) + geom_boxplot() + theme_classic() + scale_fill_manual(values = ann_cols$Tissue)
dev.off()

png(file = "Sig Analyte Effect Size by Tissue and Ome.png",width = 7,height = 5,units = "in",res = 600)
ggplot(allpowermod[allpowermod$Signif %in% c("DAR","DEG","DMR"),],aes(x = Ome,y = Cohens.d,fill = Tissue)) + geom_boxplot() + theme_classic() + scale_fill_manual(values = ann_cols$Tissue)
dev.off()

png(file = "Sig Analyte Power by Tissue and Ome.png",width = 7,height = 5,units = "in",res = 600)
ggplot(allpowermod[allpowermod$Signif %in% c("DAR","DEG","DMR"),],aes(x = Ome,y = Power,fill = Tissue)) + geom_boxplot() + theme_classic() + scale_fill_manual(values = ann_cols$Tissue)
dev.off()

png(file = "All Analyte Power by Tissue and Ome.png",width = 7,height = 5,units = "in",res = 600)
ggplot(allpowermod,aes(x = Ome,y = Power,fill = Tissue)) + geom_boxplot() + theme_classic() + scale_fill_manual(values = ann_cols$Tissue)
dev.off()

save.image("outputpoweranalysis.RData")



######
# Now that we figured out the right design matrix for the power analysis, let's try this again
#####


gastrocrossomepowerout <- MultiPower(data = gastrostatdata,groups = gastrostatdesign,type = c(2,2,2,2,2),omicPower = 0.25,averagePower = 0.5,fdr = 0.05,null.effect = 0,max.size = 1000)
heartcrossomepowerout <- MultiPower(data = heartstatdata,groups = heartstatdesign,type = c(2,2,2,2,2),omicPower = 0.25,averagePower = 0.5,fdr = 0.05,null.effect = 0,max.size = 1000)
kidneycrossomepowerout <- MultiPower(data = kidneystatdata,groups = kidneystatdesign,type = c(2,2,2,2,2),omicPower = 0.25,averagePower = 0.5,fdr = 0.05,null.effect = 0,max.size = 1000)
livercrossomepowerout <- MultiPower(data = liverstatdata,groups = liverstatdesign,type = c(2,2,2,2,2),omicPower = 0.25,averagePower = 0.5,fdr = 0.05,null.effect = 0,max.size = 1000)
lungcrossomepowerout <- MultiPower(data = lungstatdata,groups = lungstatdesign,type = c(2,2,2,2,2),omicPower = 0.25,averagePower = 0.5,fdr = 0.05,null.effect = 0,max.size = 1000)
whitecrossomepowerout <- MultiPower(data = whitestatdata,groups = whitestatdesign,type = c(2,2,2,2,2),omicPower = 0.25,averagePower = 0.5,fdr = 0.05,null.effect = 0,max.size = 1000)
hippocrossomepowerout <- MultiPower(data = hippostatdata,groups = hippostatdesign,type = c(2,2,2),omicPower = 0.25,averagePower = 0.5,fdr = 0.05,null.effect = 0,max.size = 1000)
browncrossomepowerout <- MultiPower(data = brownstatdata,groups = brownstatdesign,type = c(2,2,2),omicPower = 0.25,averagePower = 0.5,fdr = 0.05,null.effect = 0,max.size = 1000)

save.image("outputpoweranalysis.RData")

# Takes too long to run this
#gastrocrossomepowerbyomeout <- MultiPower(data = gastrostatdata,groups = gastrostatdesign,type = c(2,2,2,2,2),omicPower = 0.25,averagePower = 0.5,fdr = 0.05,null.effect = 0,max.size = 1000,equalSize = FALSE)

#save.image("outputpoweranalysis.RData")

# Now postMultiPower computation

gastrocrossomepostMultiout <- postMultiPower(optResults = gastrocrossomepowerout,max.size = 50)
heartcrossomepostMultiout <- postMultiPower(optResults = heartcrossomepowerout,max.size = 50)
hippocrossomepostMultiout <- postMultiPower(optResults = hippocrossomepowerout,max.size = 50)
kidneycrossomepostMultiout <- postMultiPower(optResults = kidneycrossomepowerout,max.size = 50)

save.image("outputpoweranalysis.RData")

livercrossomepostMultiout <- postMultiPower(optResults = livercrossomepowerout,max.size = 50)
lungcrossomepostMultiout <- postMultiPower(optResults = lungcrossomepowerout,max.size = 50)
browncrossomepostMultiout <- postMultiPower(optResults = browncrossomepowerout,max.size = 50)
whitecrossomepostMultiout <- postMultiPower(optResults = whitecrossomepowerout,max.size = 50)

save.image("outputpoweranalysis.RData")

# Let's try making this plot ourselves


####
# gastro

omicCol = colors()[c(554, 89, 111, 512, 17, 586, 132,428, 601, 568, 86, 390)]
omicCol = omicCol[1:length(gastrocrossomepowerout$parameters)]
omicShape = 1:length(gastrocrossomepowerout$parameters)
names(omicCol) = names(omicShape) = names(gastrocrossomepowerout$parameters)

nmax = max(gastrocrossomepowerout$optimalSampleSize$n)
xMin = 2
xMax = round(max(nmax + 20, (3 * nmax - xMin)/2), 0)
xValues = sort(unique(c(round(seq(xMin, xMax, (xMax - xMin)/10)),gastrocrossomepowerout$optimalSampleSize$n)))
yValues = gastrocrossomepowerout$data2plot$PowerVsSsampleSize

png(file = "Supplemental Figure S5A_112624.png",width = 6,height = 6,units = "in",res = 600)
matplot(xValues, yValues, type = "l", lwd = 2, xlab = "Sample size", 
        ylab = "Statistical power", main = "SKM-GN: Power vs Sample Size", 
        col = omicCol, lty = omicShape)
optiSS = gastrocrossomepowerout$optimalSampleSize$n
if (length(optiSS) == 1) 
  optiSS = rep(optiSS, length(gastrocrossomepowerout$parameters))
points(optiSS, diag(yValues[as.character(optiSS), ]), pch = 15, 
       col = omicCol, cex = 1.2)
legend("bottomright", names(gastrocrossomepowerout$parameters), lwd = 2, col = omicCol, 
       lty = omicShape, bty = "n")
dev.off()

pdf(file = "Supplemental Figure S5A_112624.pdf",width = 6,height = 6)
matplot(xValues, yValues, type = "l", lwd = 2, xlab = "Sample size", 
        ylab = "Statistical power", main = "SKM-GN: Power vs Sample Size", 
        col = omicCol, lty = omicShape)
optiSS = gastrocrossomepowerout$optimalSampleSize$n
if (length(optiSS) == 1) 
  optiSS = rep(optiSS, length(gastrocrossomepowerout$parameters))
points(optiSS, diag(yValues[as.character(optiSS), ]), pch = 15, 
       col = omicCol, cex = 1.2)
legend("bottomright", names(gastrocrossomepowerout$parameters), lwd = 2, col = omicCol, 
       lty = omicShape, bty = "n")
dev.off()

xValues2 = seq(0, 0.75, 0.05)
yValues2 = gastrocrossomepowerout$data2plot$PowerVsEffectSize
optiSS = gastrocrossomepowerout$optimalSampleSize$n

png(file = "Supplemental Figure S5B_112624.png",width = 6,height = 6,units = "in",res = 600)
matplot(xValues2 * 100, yValues2, type = "l", lwd = 2, 
        xlab = "Percentiles for effect size cutoff", ylab = "Statistical power", 
        main = "SKM-GN: Power vs Effect size", col = omicCol, lty = omicShape)
points(rep(0, length(gastrocrossomepowerout$parameters)), as.numeric(yValues2[1, 
]), pch = 15, col = omicCol, cex = 1.2)
legend("bottomright", names(gastrocrossomepowerout$parameters), lwd = 2, col = omicCol, 
       lty = omicShape, bty = "n")
dev.off()

pdf(file = "Supplemental Figure S5B_112624.pdf",width = 6,height = 6)
matplot(xValues2 * 100, yValues2, type = "l", lwd = 2, 
        xlab = "Percentiles for effect size cutoff", ylab = "Statistical power", 
        main = "SKM-GN: Power vs Effect size", col = omicCol, lty = omicShape)
points(rep(0, length(gastrocrossomepowerout$parameters)), as.numeric(yValues2[1, 
]), pch = 15, col = omicCol, cex = 1.2)
legend("bottomright", names(gastrocrossomepowerout$parameters), lwd = 2, col = omicCol, 
       lty = omicShape, bty = "n")
dev.off()

####
# heart

omicCol = colors()[c(554, 89, 111, 512, 17, 586, 132,428, 601, 568, 86, 390)]
omicCol = omicCol[1:length(heartcrossomepowerout$parameters)]
omicShape = 1:length(heartcrossomepowerout$parameters)
names(omicCol) = names(omicShape) = names(heartcrossomepowerout$parameters)

nmax = max(heartcrossomepowerout$optimalSampleSize$n)
xMin = 2
xMax = round(max(nmax + 20, (3 * nmax - xMin)/2), 0)
xValues = sort(unique(c(round(seq(xMin, xMax, (xMax - xMin)/10)),heartcrossomepowerout$optimalSampleSize$n)))
yValues = heartcrossomepowerout$data2plot$PowerVsSsampleSize

png(file = "Supplemental Figure S5C_112624.png",width = 6,height = 6,units = "in",res = 600)
matplot(xValues, yValues, type = "l", lwd = 2, xlab = "Sample size", 
        ylab = "Statistical power", main = "HEART: Power vs Sample Size", 
        col = omicCol, lty = omicShape)
optiSS = heartcrossomepowerout$optimalSampleSize$n
if (length(optiSS) == 1) 
  optiSS = rep(optiSS, length(heartcrossomepowerout$parameters))
points(optiSS, diag(yValues[as.character(optiSS), ]), pch = 15, 
       col = omicCol, cex = 1.2)
legend("bottomright", names(heartcrossomepowerout$parameters), lwd = 2, col = omicCol, 
       lty = omicShape, bty = "n")
dev.off()

pdf(file = "Supplemental Figure S5C_112624.pdf",width = 6,height = 6)
matplot(xValues, yValues, type = "l", lwd = 2, xlab = "Sample size", 
        ylab = "Statistical power", main = "HEART: Power vs Sample Size", 
        col = omicCol, lty = omicShape)
optiSS = heartcrossomepowerout$optimalSampleSize$n
if (length(optiSS) == 1) 
  optiSS = rep(optiSS, length(heartcrossomepowerout$parameters))
points(optiSS, diag(yValues[as.character(optiSS), ]), pch = 15, 
       col = omicCol, cex = 1.2)
legend("bottomright", names(heartcrossomepowerout$parameters), lwd = 2, col = omicCol, 
       lty = omicShape, bty = "n")
dev.off()

xValues2 = seq(0, 0.75, 0.05)
yValues2 = heartcrossomepowerout$data2plot$PowerVsEffectSize
optiSS = heartcrossomepowerout$optimalSampleSize$n

png(file = "Supplemental Figure S5D_112624.png",width = 6,height = 6,units = "in",res = 600)
matplot(xValues2 * 100, yValues2, type = "l", lwd = 2, 
        xlab = "Percentiles for effect size cutoff", ylab = "Statistical power", 
        main = "HEART: Power vs Effect size", col = omicCol, lty = omicShape)
points(rep(0, length(heartcrossomepowerout$parameters)), as.numeric(yValues2[1, 
]), pch = 15, col = omicCol, cex = 1.2)
legend("bottomright", names(heartcrossomepowerout$parameters), lwd = 2, col = omicCol, 
       lty = omicShape, bty = "n")
dev.off()

pdf(file = "Supplemental Figure S5D_112624.pdf",width = 6,height = 6)
matplot(xValues2 * 100, yValues2, type = "l", lwd = 2, 
        xlab = "Percentiles for effect size cutoff", ylab = "Statistical power", 
        main = "HEART: Power vs Effect size", col = omicCol, lty = omicShape)
points(rep(0, length(heartcrossomepowerout$parameters)), as.numeric(yValues2[1, 
]), pch = 15, col = omicCol, cex = 1.2)
legend("bottomright", names(heartcrossomepowerout$parameters), lwd = 2, col = omicCol, 
       lty = omicShape, bty = "n")
dev.off()


####
# kidney

omicCol = colors()[c(554, 89, 111, 512, 17, 586, 132,428, 601, 568, 86, 390)]
omicCol = omicCol[1:length(kidneycrossomepowerout$parameters)]
omicShape = 1:length(kidneycrossomepowerout$parameters)
names(omicCol) = names(omicShape) = names(kidneycrossomepowerout$parameters)

nmax = max(kidneycrossomepowerout$optimalSampleSize$n)
xMin = 2
xMax = round(max(nmax + 20, (3 * nmax - xMin)/2), 0)
xValues = sort(unique(c(round(seq(xMin, xMax, (xMax - xMin)/10)),kidneycrossomepowerout$optimalSampleSize$n)))
yValues = kidneycrossomepowerout$data2plot$PowerVsSsampleSize

png(file = "Supplemental Figure S5E_112624.png",width = 6,height = 6,units = "in",res = 600)
matplot(xValues, yValues, type = "l", lwd = 2, xlab = "Sample size", 
        ylab = "Statistical power", main = "KIDNEY: Power vs Sample Size", 
        col = omicCol, lty = omicShape)
optiSS = kidneycrossomepowerout$optimalSampleSize$n
if (length(optiSS) == 1) 
  optiSS = rep(optiSS, length(kidneycrossomepowerout$parameters))
points(optiSS, diag(yValues[as.character(optiSS), ]), pch = 15, 
       col = omicCol, cex = 1.2)
legend("bottomright", names(kidneycrossomepowerout$parameters), lwd = 2, col = omicCol, 
       lty = omicShape, bty = "n")
dev.off()

pdf(file = "Supplemental Figure S5E_112624.pdf",width = 6,height = 6)
matplot(xValues, yValues, type = "l", lwd = 2, xlab = "Sample size", 
        ylab = "Statistical power", main = "KIDNEY: Power vs Sample Size", 
        col = omicCol, lty = omicShape)
optiSS = kidneycrossomepowerout$optimalSampleSize$n
if (length(optiSS) == 1) 
  optiSS = rep(optiSS, length(kidneycrossomepowerout$parameters))
points(optiSS, diag(yValues[as.character(optiSS), ]), pch = 15, 
       col = omicCol, cex = 1.2)
legend("bottomright", names(kidneycrossomepowerout$parameters), lwd = 2, col = omicCol, 
       lty = omicShape, bty = "n")
dev.off()

xValues2 = seq(0, 0.75, 0.05)
yValues2 = kidneycrossomepowerout$data2plot$PowerVsEffectSize
optiSS = kidneycrossomepowerout$optimalSampleSize$n

png(file = "Supplemental Figure S5F_112624.png",width = 6,height = 6,units = "in",res = 600)
matplot(xValues2 * 100, yValues2, type = "l", lwd = 2, 
        xlab = "Percentiles for effect size cutoff", ylab = "Statistical power", 
        main = "KIDNEY: Power vs Effect size", col = omicCol, lty = omicShape)
points(rep(0, length(kidneycrossomepowerout$parameters)), as.numeric(yValues2[1, 
]), pch = 15, col = omicCol, cex = 1.2)
legend("bottomright", names(kidneycrossomepowerout$parameters), lwd = 2, col = omicCol, 
       lty = omicShape, bty = "n")
dev.off()

pdf(file = "Supplemental Figure S5F_112624.pdf",width = 6,height = 6)
matplot(xValues2 * 100, yValues2, type = "l", lwd = 2, 
        xlab = "Percentiles for effect size cutoff", ylab = "Statistical power", 
        main = "KIDNEY: Power vs Effect size", col = omicCol, lty = omicShape)
points(rep(0, length(kidneycrossomepowerout$parameters)), as.numeric(yValues2[1, 
]), pch = 15, col = omicCol, cex = 1.2)
legend("bottomright", names(kidneycrossomepowerout$parameters), lwd = 2, col = omicCol, 
       lty = omicShape, bty = "n")
dev.off()

####
# liver

omicCol = colors()[c(554, 89, 111, 512, 17, 586, 132,428, 601, 568, 86, 390)]
omicCol = omicCol[1:length(livercrossomepowerout$parameters)]
omicShape = 1:length(livercrossomepowerout$parameters)
names(omicCol) = names(omicShape) = names(livercrossomepowerout$parameters)

nmax = max(livercrossomepowerout$optimalSampleSize$n)
xMin = 2
xMax = round(max(nmax + 20, (3 * nmax - xMin)/2), 0)
xValues = sort(unique(c(round(seq(xMin, xMax, (xMax - xMin)/10)),livercrossomepowerout$optimalSampleSize$n)))
yValues = livercrossomepowerout$data2plot$PowerVsSsampleSize

png(file = "Supplemental Figure S5G_112624.png",width = 6,height = 6,units = "in",res = 600)
matplot(xValues, yValues, type = "l", lwd = 2, xlab = "Sample size", 
        ylab = "Statistical power", main = "LIVER: Power vs Sample Size", 
        col = omicCol, lty = omicShape)
optiSS = livercrossomepowerout$optimalSampleSize$n
if (length(optiSS) == 1) 
  optiSS = rep(optiSS, length(livercrossomepowerout$parameters))
points(optiSS, diag(yValues[as.character(optiSS), ]), pch = 15, 
       col = omicCol, cex = 1.2)
legend("bottomright", names(livercrossomepowerout$parameters), lwd = 2, col = omicCol, 
       lty = omicShape, bty = "n")
dev.off()

pdf(file = "Supplemental Figure S5G_112624.pdf",width = 6,height = 6)
matplot(xValues, yValues, type = "l", lwd = 2, xlab = "Sample size", 
        ylab = "Statistical power", main = "LIVER: Power vs Sample Size", 
        col = omicCol, lty = omicShape)
optiSS = livercrossomepowerout$optimalSampleSize$n
if (length(optiSS) == 1) 
  optiSS = rep(optiSS, length(livercrossomepowerout$parameters))
points(optiSS, diag(yValues[as.character(optiSS), ]), pch = 15, 
       col = omicCol, cex = 1.2)
legend("bottomright", names(livercrossomepowerout$parameters), lwd = 2, col = omicCol, 
       lty = omicShape, bty = "n")
dev.off()

xValues2 = seq(0, 0.75, 0.05)
yValues2 = livercrossomepowerout$data2plot$PowerVsEffectSize
optiSS = livercrossomepowerout$optimalSampleSize$n

png(file = "Supplemental Figure S5H_112624.png",width = 6,height = 6,units = "in",res = 600)
matplot(xValues2 * 100, yValues2, type = "l", lwd = 2, 
        xlab = "Percentiles for effect size cutoff", ylab = "Statistical power", 
        main = "LIVER: Power vs Effect size", col = omicCol, lty = omicShape)
points(rep(0, length(livercrossomepowerout$parameters)), as.numeric(yValues2[1, 
]), pch = 15, col = omicCol, cex = 1.2)
legend("bottomright", names(livercrossomepowerout$parameters), lwd = 2, col = omicCol, 
       lty = omicShape, bty = "n")
dev.off()

pdf(file = "Supplemental Figure S5H_112624.pdf",width = 6,height = 6)
matplot(xValues2 * 100, yValues2, type = "l", lwd = 2, 
        xlab = "Percentiles for effect size cutoff", ylab = "Statistical power", 
        main = "LIVER: Power vs Effect size", col = omicCol, lty = omicShape)
points(rep(0, length(livercrossomepowerout$parameters)), as.numeric(yValues2[1, 
]), pch = 15, col = omicCol, cex = 1.2)
legend("bottomright", names(livercrossomepowerout$parameters), lwd = 2, col = omicCol, 
       lty = omicShape, bty = "n")
dev.off()


####
# lung

omicCol = colors()[c(554, 89, 111, 512, 17, 586, 132,428, 601, 568, 86, 390)]
omicCol = omicCol[1:length(lungcrossomepowerout$parameters)]
omicShape = 1:length(lungcrossomepowerout$parameters)
names(omicCol) = names(omicShape) = names(lungcrossomepowerout$parameters)

nmax = max(lungcrossomepowerout$optimalSampleSize$n)
xMin = 2
xMax = round(max(nmax + 20, (3 * nmax - xMin)/2), 0)
xValues = sort(unique(c(round(seq(xMin, xMax, (xMax - xMin)/10)),lungcrossomepowerout$optimalSampleSize$n)))
yValues = lungcrossomepowerout$data2plot$PowerVsSsampleSize

png(file = "Supplemental Figure S5I_112624.png",width = 6,height = 6,units = "in",res = 600)
matplot(xValues, yValues, type = "l", lwd = 2, xlab = "Sample size", 
        ylab = "Statistical power", main = "LUNG: Power vs Sample Size", 
        col = omicCol, lty = omicShape)
optiSS = lungcrossomepowerout$optimalSampleSize$n
if (length(optiSS) == 1) 
  optiSS = rep(optiSS, length(lungcrossomepowerout$parameters))
points(optiSS, diag(yValues[as.character(optiSS), ]), pch = 15, 
       col = omicCol, cex = 1.2)
legend("bottomright", names(lungcrossomepowerout$parameters), lwd = 2, col = omicCol, 
       lty = omicShape, bty = "n")
dev.off()

pdf(file = "Supplemental Figure S5I_112624.pdf",width = 6,height = 6)
matplot(xValues, yValues, type = "l", lwd = 2, xlab = "Sample size", 
        ylab = "Statistical power", main = "LUNG: Power vs Sample Size", 
        col = omicCol, lty = omicShape)
optiSS = lungcrossomepowerout$optimalSampleSize$n
if (length(optiSS) == 1) 
  optiSS = rep(optiSS, length(lungcrossomepowerout$parameters))
points(optiSS, diag(yValues[as.character(optiSS), ]), pch = 15, 
       col = omicCol, cex = 1.2)
legend("bottomright", names(lungcrossomepowerout$parameters), lwd = 2, col = omicCol, 
       lty = omicShape, bty = "n")
dev.off()

xValues2 = seq(0, 0.75, 0.05)
yValues2 = lungcrossomepowerout$data2plot$PowerVsEffectSize
optiSS = lungcrossomepowerout$optimalSampleSize$n

png(file = "Supplemental Figure S5J_112624.png",width = 6,height = 6,units = "in",res = 600)
matplot(xValues2 * 100, yValues2, type = "l", lwd = 2, 
        xlab = "Percentiles for effect size cutoff", ylab = "Statistical power", 
        main = "LUNG: Power vs Effect size", col = omicCol, lty = omicShape)
points(rep(0, length(lungcrossomepowerout$parameters)), as.numeric(yValues2[1, 
]), pch = 15, col = omicCol, cex = 1.2)
legend("bottomright", names(lungcrossomepowerout$parameters), lwd = 2, col = omicCol, 
       lty = omicShape, bty = "n")
dev.off()

pdf(file = "Supplemental Figure S5J_112624.pdf",width = 6,height = 6)
matplot(xValues2 * 100, yValues2, type = "l", lwd = 2, 
        xlab = "Percentiles for effect size cutoff", ylab = "Statistical power", 
        main = "LUNG: Power vs Effect size", col = omicCol, lty = omicShape)
points(rep(0, length(lungcrossomepowerout$parameters)), as.numeric(yValues2[1, 
]), pch = 15, col = omicCol, cex = 1.2)
legend("bottomright", names(lungcrossomepowerout$parameters), lwd = 2, col = omicCol, 
       lty = omicShape, bty = "n")
dev.off()

####
# white

omicCol = colors()[c(554, 89, 111, 512, 17, 586, 132,428, 601, 568, 86, 390)]
omicCol = omicCol[1:length(whitecrossomepowerout$parameters)]
omicShape = 1:length(whitecrossomepowerout$parameters)
names(omicCol) = names(omicShape) = names(whitecrossomepowerout$parameters)

nmax = max(whitecrossomepowerout$optimalSampleSize$n)
xMin = 2
xMax = round(max(nmax + 20, (3 * nmax - xMin)/2), 0)
xValues = sort(unique(c(round(seq(xMin, xMax, (xMax - xMin)/10)),whitecrossomepowerout$optimalSampleSize$n)))
yValues = whitecrossomepowerout$data2plot$PowerVsSsampleSize

png(file = "Supplemental Figure S5K_112624.png",width = 6,height = 6,units = "in",res = 600)
matplot(xValues, yValues, type = "l", lwd = 2, xlab = "Sample size", 
        ylab = "Statistical power", main = "WAT-SC: Power vs Sample Size", 
        col = omicCol, lty = omicShape)
optiSS = whitecrossomepowerout$optimalSampleSize$n
if (length(optiSS) == 1) 
  optiSS = rep(optiSS, length(whitecrossomepowerout$parameters))
points(optiSS, diag(yValues[as.character(optiSS), ]), pch = 15, 
       col = omicCol, cex = 1.2)
legend("bottomright", names(whitecrossomepowerout$parameters), lwd = 2, col = omicCol, 
       lty = omicShape, bty = "n")
dev.off()

pdf(file = "Supplemental Figure S5K_112624.pdf",width = 6,height = 6)
matplot(xValues, yValues, type = "l", lwd = 2, xlab = "Sample size", 
        ylab = "Statistical power", main = "WAT-SC: Power vs Sample Size", 
        col = omicCol, lty = omicShape)
optiSS = whitecrossomepowerout$optimalSampleSize$n
if (length(optiSS) == 1) 
  optiSS = rep(optiSS, length(whitecrossomepowerout$parameters))
points(optiSS, diag(yValues[as.character(optiSS), ]), pch = 15, 
       col = omicCol, cex = 1.2)
legend("bottomright", names(whitecrossomepowerout$parameters), lwd = 2, col = omicCol, 
       lty = omicShape, bty = "n")
dev.off()

xValues2 = seq(0, 0.75, 0.05)
yValues2 = whitecrossomepowerout$data2plot$PowerVsEffectSize
optiSS = whitecrossomepowerout$optimalSampleSize$n

png(file = "Supplemental Figure S5L_112624.png",width = 6,height = 6,units = "in",res = 600)
matplot(xValues2 * 100, yValues2, type = "l", lwd = 2, 
        xlab = "Percentiles for effect size cutoff", ylab = "Statistical power", 
        main = "WAT-SC: Power vs Effect size", col = omicCol, lty = omicShape)
points(rep(0, length(whitecrossomepowerout$parameters)), as.numeric(yValues2[1, 
]), pch = 15, col = omicCol, cex = 1.2)
legend("bottomright", names(whitecrossomepowerout$parameters), lwd = 2, col = omicCol, 
       lty = omicShape, bty = "n")
dev.off()

pdf(file = "Supplemental Figure S5L_112624.pdf",width = 6,height = 6)
matplot(xValues2 * 100, yValues2, type = "l", lwd = 2, 
        xlab = "Percentiles for effect size cutoff", ylab = "Statistical power", 
        main = "WAT-SC: Power vs Effect size", col = omicCol, lty = omicShape)
points(rep(0, length(whitecrossomepowerout$parameters)), as.numeric(yValues2[1, 
]), pch = 15, col = omicCol, cex = 1.2)
legend("bottomright", names(whitecrossomepowerout$parameters), lwd = 2, col = omicCol, 
       lty = omicShape, bty = "n")
dev.off()


####
# hippo

omicCol = colors()[c(554, 89, 111, 512, 17, 586, 132,428, 601, 568, 86, 390)]
omicCol = omicCol[1:length(hippocrossomepowerout$parameters)]
omicShape = 1:length(hippocrossomepowerout$parameters)
names(omicCol) = names(omicShape) = names(hippocrossomepowerout$parameters)

nmax = max(hippocrossomepowerout$optimalSampleSize$n)
xMin = 2
xMax = round(max(nmax + 20, (3 * nmax - xMin)/2), 0)
xValues = sort(unique(c(round(seq(xMin, xMax, (xMax - xMin)/10)),hippocrossomepowerout$optimalSampleSize$n)))
yValues = hippocrossomepowerout$data2plot$PowerVsSsampleSize

png(file = "Supplemental Figure S5M_112624.png",width = 6,height = 6,units = "in",res = 600)
matplot(xValues, yValues, type = "l", lwd = 2, xlab = "Sample size", 
        ylab = "Statistical power", main = "HIPPOC: Power vs Sample Size", 
        col = omicCol, lty = omicShape)
optiSS = hippocrossomepowerout$optimalSampleSize$n
if (length(optiSS) == 1) 
  optiSS = rep(optiSS, length(hippocrossomepowerout$parameters))
points(optiSS, diag(yValues[as.character(optiSS), ]), pch = 15, 
       col = omicCol, cex = 1.2)
legend("bottomright", names(hippocrossomepowerout$parameters), lwd = 2, col = omicCol, 
       lty = omicShape, bty = "n")
dev.off()

pdf(file = "Supplemental Figure S5M_112624.pdf",width = 6,height = 6)
matplot(xValues, yValues, type = "l", lwd = 2, xlab = "Sample size", 
        ylab = "Statistical power", main = "HIPPOC: Power vs Sample Size", 
        col = omicCol, lty = omicShape)
optiSS = hippocrossomepowerout$optimalSampleSize$n
if (length(optiSS) == 1) 
  optiSS = rep(optiSS, length(hippocrossomepowerout$parameters))
points(optiSS, diag(yValues[as.character(optiSS), ]), pch = 15, 
       col = omicCol, cex = 1.2)
legend("bottomright", names(hippocrossomepowerout$parameters), lwd = 2, col = omicCol, 
       lty = omicShape, bty = "n")
dev.off()

xValues2 = seq(0, 0.75, 0.05)
yValues2 = hippocrossomepowerout$data2plot$PowerVsEffectSize
optiSS = hippocrossomepowerout$optimalSampleSize$n

png(file = "Supplemental Figure S5N_112624.png",width = 6,height = 6,units = "in",res = 600)
matplot(xValues2 * 100, yValues2, type = "l", lwd = 2, 
        xlab = "Percentiles for effect size cutoff", ylab = "Statistical power", 
        main = "HIPPOC: Power vs Effect size", col = omicCol, lty = omicShape)
points(rep(0, length(hippocrossomepowerout$parameters)), as.numeric(yValues2[1, 
]), pch = 15, col = omicCol, cex = 1.2)
legend("bottomright", names(hippocrossomepowerout$parameters), lwd = 2, col = omicCol, 
       lty = omicShape, bty = "n")
dev.off()

pdf(file = "Supplemental Figure S5N_112624.pdf",width = 6,height = 6)
matplot(xValues2 * 100, yValues2, type = "l", lwd = 2, 
        xlab = "Percentiles for effect size cutoff", ylab = "Statistical power", 
        main = "HIPPOC: Power vs Effect size", col = omicCol, lty = omicShape)
points(rep(0, length(hippocrossomepowerout$parameters)), as.numeric(yValues2[1, 
]), pch = 15, col = omicCol, cex = 1.2)
legend("bottomright", names(hippocrossomepowerout$parameters), lwd = 2, col = omicCol, 
       lty = omicShape, bty = "n")
dev.off()

####
# brown

omicCol = colors()[c(554, 89, 111, 512, 17, 586, 132,428, 601, 568, 86, 390)]
omicCol = omicCol[1:length(browncrossomepowerout$parameters)]
omicShape = 1:length(browncrossomepowerout$parameters)
names(omicCol) = names(omicShape) = names(browncrossomepowerout$parameters)

nmax = max(browncrossomepowerout$optimalSampleSize$n)
xMin = 2
xMax = round(max(nmax + 20, (3 * nmax - xMin)/2), 0)
xValues = sort(unique(c(round(seq(xMin, xMax, (xMax - xMin)/10)),browncrossomepowerout$optimalSampleSize$n)))
yValues = browncrossomepowerout$data2plot$PowerVsSsampleSize

png(file = "Supplemental Figure S5O_112624.png",width = 6,height = 6,units = "in",res = 600)
matplot(xValues, yValues, type = "l", lwd = 2, xlab = "Sample size", 
        ylab = "Statistical power", main = "BAT: Power vs Sample Size", 
        col = omicCol, lty = omicShape)
optiSS = browncrossomepowerout$optimalSampleSize$n
if (length(optiSS) == 1) 
  optiSS = rep(optiSS, length(browncrossomepowerout$parameters))
points(optiSS, diag(yValues[as.character(optiSS), ]), pch = 15, 
       col = omicCol, cex = 1.2)
legend("bottomright", names(browncrossomepowerout$parameters), lwd = 2, col = omicCol, 
       lty = omicShape, bty = "n")
dev.off()

pdf(file = "Supplemental Figure S5O_112624.pdf",width = 6,height = 6)
matplot(xValues, yValues, type = "l", lwd = 2, xlab = "Sample size", 
        ylab = "Statistical power", main = "BAT: Power vs Sample Size", 
        col = omicCol, lty = omicShape)
optiSS = browncrossomepowerout$optimalSampleSize$n
if (length(optiSS) == 1) 
  optiSS = rep(optiSS, length(browncrossomepowerout$parameters))
points(optiSS, diag(yValues[as.character(optiSS), ]), pch = 15, 
       col = omicCol, cex = 1.2)
legend("bottomright", names(browncrossomepowerout$parameters), lwd = 2, col = omicCol, 
       lty = omicShape, bty = "n")
dev.off()

xValues2 = seq(0, 0.75, 0.05)
yValues2 = browncrossomepowerout$data2plot$PowerVsEffectSize
optiSS = browncrossomepowerout$optimalSampleSize$n

png(file = "Supplemental Figure S5P_112624.png",width = 6,height = 6,units = "in",res = 600)
matplot(xValues2 * 100, yValues2, type = "l", lwd = 2, 
        xlab = "Percentiles for effect size cutoff", ylab = "Statistical power", 
        main = "BAT: Power vs Effect size", col = omicCol, lty = omicShape)
points(rep(0, length(browncrossomepowerout$parameters)), as.numeric(yValues2[1, 
]), pch = 15, col = omicCol, cex = 1.2)
legend("bottomright", names(browncrossomepowerout$parameters), lwd = 2, col = omicCol, 
       lty = omicShape, bty = "n")
dev.off()

pdf(file = "Supplemental Figure S5P_112624.pdf",width = 6,height = 6)
matplot(xValues2 * 100, yValues2, type = "l", lwd = 2, 
        xlab = "Percentiles for effect size cutoff", ylab = "Statistical power", 
        main = "BAT: Power vs Effect size", col = omicCol, lty = omicShape)
points(rep(0, length(browncrossomepowerout$parameters)), as.numeric(yValues2[1, 
]), pch = 15, col = omicCol, cex = 1.2)
legend("bottomright", names(browncrossomepowerout$parameters), lwd = 2, col = omicCol, 
       lty = omicShape, bty = "n")
dev.off()

#####
# Now to plot the postPowerPlot figures
####

max.size = 50

####
# gastro

omicCol = colors()[c(554, 89, 111, 512, 17, 586, 132,428, 601, 568, 86, 390)]
omicCol = omicCol[1:ncol(gastrocrossomepostMultiout$SampleSize)]
names(omicCol) = colnames(gastrocrossomepostMultiout$SampleSize)
if (min(gastrocrossomepostMultiout$SampleSize) > max.size) {
  cat(paste0("The chosen sample size of ", max.size, " is not feasible. \n"))
  max.size = min(gastrocrossomepostMultiout$SampleSize)
  cat(paste0("A sample size of ", max.size, " will be plotted instead. \n"))
}

fff = approxfun(gastrocrossomepostMultiout$SampleSize[, 1], gastrocrossomepostMultiout$d)
myD = round(fff(max.size), 1)
cat(paste0("For having a sample size of ", max.size, 
             " and maintain the desired power, you need to remove features with Cohen's d below ", 
             myD, ". \n"))
cat("The number of remaining features in each omic is: \n")
print(gastrocrossomepostMultiout$NumFeat[paste0("d=", myD), ])

png(file = "gastro_sample_size_vs_effect_size.png",width = 6,height = 6,units = "in",res = 600)
plot(gastrocrossomepostMultiout$d, gastrocrossomepostMultiout$SampleSize[, 1], type = "l", 
       lwd = 2, xlab = "Cohen's d cutoff", ylab = "Number of replicates", 
       main = "SKM-GN: Sample size vs Cohen's d", ylim = c(2, max(gastrocrossomepostMultiout$SampleSize)))
arrows(x0 = 0, y0 = max.size, x1 = myD, y1 = max.size,lty = 2, col = 2)
arrows(x0 = myD, y0 = max.size, x1 = myD, y1 = 2, lty = 2,col = 2)
text(min(gastrocrossomepostMultiout$d) + 1, max(gastrocrossomepostMultiout$SampleSize, 
                                   na.rm = TRUE) - 2, paste0("Cohen's d = ", myD), col = 2)
dev.off()  

png(file = "gastro_power_comparison.png",width = 6,height = 6,units = "in",res = 600)
barplot(gastrocrossomepostMultiout$Power[c(1, min(which(gastrocrossomepostMultiout$d >= myD))), ],
        col = rep(omicCol, each = 2), beside = TRUE, 
        las = 2, ylab = "Statistical power", ylim = c(0, 1), 
        border = rep(omicCol, each = 2), density = rep(c(30,100), ncol(gastrocrossomepostMultiout$Power)),main = "SKM-GN: Power Comparison")
legend(x = 0.5, y = 1.0, c("Optimal SS", "User's SS"), col = 1, 
       density = c(30, 100), ncol = 2, bty = "n")
dev.off()


####
# heart

omicCol = colors()[c(554, 89, 111, 512, 17, 586, 132,428, 601, 568, 86, 390)]
omicCol = omicCol[1:ncol(heartcrossomepostMultiout$SampleSize)]
names(omicCol) = colnames(heartcrossomepostMultiout$SampleSize)
if (min(heartcrossomepostMultiout$SampleSize) > max.size) {
  cat(paste0("The chosen sample size of ", max.size, " is not feasible. \n"))
  max.size = min(heartcrossomepostMultiout$SampleSize)
  cat(paste0("A sample size of ", max.size, " will be plotted instead. \n"))
}

fff = approxfun(heartcrossomepostMultiout$SampleSize[, 1], heartcrossomepostMultiout$d)
myD = round(fff(max.size), 1)
cat(paste0("For having a sample size of ", max.size, 
           " and maintain the desired power, you need to remove features with Cohen's d below ", 
           myD, ". \n"))
cat("The number of remaining features in each omic is: \n")
print(heartcrossomepostMultiout$NumFeat[paste0("d=", myD), ])

png(file = "heart_sample_size_vs_effect_size.png",width = 6,height = 6,units = "in",res = 600)
plot(heartcrossomepostMultiout$d, heartcrossomepostMultiout$SampleSize[, 1], type = "l", 
     lwd = 2, xlab = "Cohen's d cutoff", ylab = "Number of replicates", 
     main = "HEART: Sample size vs Cohen's d", ylim = c(2, max(heartcrossomepostMultiout$SampleSize)))
arrows(x0 = 0, y0 = max.size, x1 = myD, y1 = max.size,lty = 2, col = 2)
arrows(x0 = myD, y0 = max.size, x1 = myD, y1 = 2, lty = 2,col = 2)
text(min(heartcrossomepostMultiout$d) + 1, max(heartcrossomepostMultiout$SampleSize, 
                                                na.rm = TRUE) - 2, paste0("Cohen's d = ", myD), col = 2)
dev.off()  

png(file = "heart_power_comparison.png",width = 6,height = 6,units = "in",res = 600)
barplot(heartcrossomepostMultiout$Power[c(1, min(which(heartcrossomepostMultiout$d >= myD))), ],
        col = rep(omicCol, each = 2), beside = TRUE, 
        las = 2, ylab = "Statistical power", ylim = c(0, 1), 
        border = rep(omicCol, each = 2), density = rep(c(30,100), ncol(heartcrossomepostMultiout$Power)),main = "HEART: Power Comparison")
legend(x = 0.5, y = 1.0, c("Optimal SS", "User's SS"), col = 1, 
       density = c(30, 100), ncol = 2, bty = "n")
dev.off()


####
# hippo

omicCol = colors()[c(554, 89, 111, 512, 17, 586, 132,428, 601, 568, 86, 390)]
omicCol = omicCol[1:ncol(hippocrossomepostMultiout$SampleSize)]
names(omicCol) = colnames(hippocrossomepostMultiout$SampleSize)
if (min(hippocrossomepostMultiout$SampleSize) > max.size) {
  cat(paste0("The chosen sample size of ", max.size, " is not feasible. \n"))
  max.size = min(hippocrossomepostMultiout$SampleSize)
  cat(paste0("A sample size of ", max.size, " will be plotted instead. \n"))
}

fff = approxfun(hippocrossomepostMultiout$SampleSize[, 1], hippocrossomepostMultiout$d)
myD = round(fff(max.size), 1)
cat(paste0("For having a sample size of ", max.size, 
           " and maintain the desired power, you need to remove features with Cohen's d below ", 
           myD, ". \n"))
cat("The number of remaining features in each omic is: \n")
print(hippocrossomepostMultiout$NumFeat[paste0("d=", myD), ])

png(file = "hippo_sample_size_vs_effect_size.png",width = 6,height = 6,units = "in",res = 600)
plot(hippocrossomepostMultiout$d, hippocrossomepostMultiout$SampleSize[, 1], type = "l", 
     lwd = 2, xlab = "Cohen's d cutoff", ylab = "Number of replicates", 
     main = "HIPPOC: Sample size vs Cohen's d", ylim = c(2, max(hippocrossomepostMultiout$SampleSize)))
arrows(x0 = 0, y0 = max.size, x1 = myD, y1 = max.size,lty = 2, col = 2)
arrows(x0 = myD, y0 = max.size, x1 = myD, y1 = 2, lty = 2,col = 2)
text(min(hippocrossomepostMultiout$d) + 1, max(hippocrossomepostMultiout$SampleSize, 
                                                na.rm = TRUE) - 2, paste0("Cohen's d = ", myD), col = 2)
dev.off()  

png(file = "hippo_power_comparison.png",width = 6,height = 6,units = "in",res = 600)
barplot(hippocrossomepostMultiout$Power[c(1, min(which(hippocrossomepostMultiout$d >= myD))), ],
        col = rep(omicCol, each = 2), beside = TRUE, 
        las = 2, ylab = "Statistical power", ylim = c(0, 1), 
        border = rep(omicCol, each = 2), density = rep(c(30,100), ncol(hippocrossomepostMultiout$Power)),main = "HIPPOC: Power Comparison")
legend(x = 0.5, y = 1.0, c("Optimal SS", "User's SS"), col = 1, 
       density = c(30, 100), ncol = 2, bty = "n")
dev.off()


####
# kidney

omicCol = colors()[c(554, 89, 111, 512, 17, 586, 132,428, 601, 568, 86, 390)]
omicCol = omicCol[1:ncol(kidneycrossomepostMultiout$SampleSize)]
names(omicCol) = colnames(kidneycrossomepostMultiout$SampleSize)
if (min(kidneycrossomepostMultiout$SampleSize) > max.size) {
  cat(paste0("The chosen sample size of ", max.size, " is not feasible. \n"))
  max.size = min(kidneycrossomepostMultiout$SampleSize)
  cat(paste0("A sample size of ", max.size, " will be plotted instead. \n"))
}

fff = approxfun(kidneycrossomepostMultiout$SampleSize[, 1], kidneycrossomepostMultiout$d)
myD = round(fff(max.size), 1)
cat(paste0("For having a sample size of ", max.size, 
           " and maintain the desired power, you need to remove features with Cohen's d below ", 
           myD, ". \n"))
cat("The number of remaining features in each omic is: \n")
print(kidneycrossomepostMultiout$NumFeat[paste0("d=", myD), ])

png(file = "kidney_sample_size_vs_effect_size.png",width = 6,height = 6,units = "in",res = 600)
plot(kidneycrossomepostMultiout$d, kidneycrossomepostMultiout$SampleSize[, 1], type = "l", 
     lwd = 2, xlab = "Cohen's d cutoff", ylab = "Number of replicates", 
     main = "KIDNEY: Sample size vs Cohen's d", ylim = c(2, max(kidneycrossomepostMultiout$SampleSize)))
arrows(x0 = 0, y0 = max.size, x1 = myD, y1 = max.size,lty = 2, col = 2)
arrows(x0 = myD, y0 = max.size, x1 = myD, y1 = 2, lty = 2,col = 2)
text(min(kidneycrossomepostMultiout$d) + 1, max(kidneycrossomepostMultiout$SampleSize, 
                                                na.rm = TRUE) - 2, paste0("Cohen's d = ", myD), col = 2)
dev.off()  

png(file = "kidney_power_comparison.png",width = 6,height = 6,units = "in",res = 600)
barplot(kidneycrossomepostMultiout$Power[c(1, min(which(kidneycrossomepostMultiout$d >= myD))), ],
        col = rep(omicCol, each = 2), beside = TRUE, 
        las = 2, ylab = "Statistical power", ylim = c(0, 1), 
        border = rep(omicCol, each = 2), density = rep(c(30,100), ncol(kidneycrossomepostMultiout$Power)),main = "KIDNEY: Power Comparison")
legend(x = 0.5, y = 1.0, c("Optimal SS", "User's SS"), col = 1, 
       density = c(30, 100), ncol = 2, bty = "n")
dev.off()


####
# liver

omicCol = colors()[c(554, 89, 111, 512, 17, 586, 132,428, 601, 568, 86, 390)]
omicCol = omicCol[1:ncol(livercrossomepostMultiout$SampleSize)]
names(omicCol) = colnames(livercrossomepostMultiout$SampleSize)
if (min(livercrossomepostMultiout$SampleSize) > max.size) {
  cat(paste0("The chosen sample size of ", max.size, " is not feasible. \n"))
  max.size = min(livercrossomepostMultiout$SampleSize)
  cat(paste0("A sample size of ", max.size, " will be plotted instead. \n"))
}

fff = approxfun(livercrossomepostMultiout$SampleSize[, 1], livercrossomepostMultiout$d)
myD = round(fff(max.size), 1)
cat(paste0("For having a sample size of ", max.size, 
           " and maintain the desired power, you need to remove features with Cohen's d below ", 
           myD, ". \n"))
cat("The number of remaining features in each omic is: \n")
print(livercrossomepostMultiout$NumFeat[paste0("d=", myD), ])

png(file = "liver_sample_size_vs_effect_size.png",width = 6,height = 6,units = "in",res = 600)
plot(livercrossomepostMultiout$d, livercrossomepostMultiout$SampleSize[, 1], type = "l", 
     lwd = 2, xlab = "Cohen's d cutoff", ylab = "Number of replicates", 
     main = "LIVER: Sample size vs Cohen's d", ylim = c(2, max(livercrossomepostMultiout$SampleSize)))
arrows(x0 = 0, y0 = max.size, x1 = myD, y1 = max.size,lty = 2, col = 2)
arrows(x0 = myD, y0 = max.size, x1 = myD, y1 = 2, lty = 2,col = 2)
text(min(livercrossomepostMultiout$d) + 1, max(livercrossomepostMultiout$SampleSize, 
                                                na.rm = TRUE) - 2, paste0("Cohen's d = ", myD), col = 2)
dev.off()  

png(file = "liver_power_comparison.png",width = 6,height = 6,units = "in",res = 600)
barplot(livercrossomepostMultiout$Power[c(1, min(which(livercrossomepostMultiout$d >= myD))), ],
        col = rep(omicCol, each = 2), beside = TRUE, 
        las = 2, ylab = "Statistical power", ylim = c(0, 1), 
        border = rep(omicCol, each = 2), density = rep(c(30,100), ncol(livercrossomepostMultiout$Power)),main = "LIVER: Power Comparison")
legend(x = 0.5, y = 1.0, c("Optimal SS", "User's SS"), col = 1, 
       density = c(30, 100), ncol = 2, bty = "n")
dev.off()


####
# lung

omicCol = colors()[c(554, 89, 111, 512, 17, 586, 132,428, 601, 568, 86, 390)]
omicCol = omicCol[1:ncol(lungcrossomepostMultiout$SampleSize)]
names(omicCol) = colnames(lungcrossomepostMultiout$SampleSize)
if (min(lungcrossomepostMultiout$SampleSize) > max.size) {
  cat(paste0("The chosen sample size of ", max.size, " is not feasible. \n"))
  max.size = min(lungcrossomepostMultiout$SampleSize)
  cat(paste0("A sample size of ", max.size, " will be plotted instead. \n"))
}

fff = approxfun(lungcrossomepostMultiout$SampleSize[, 1], lungcrossomepostMultiout$d)
myD = round(fff(max.size), 1)
cat(paste0("For having a sample size of ", max.size, 
           " and maintain the desired power, you need to remove features with Cohen's d below ", 
           myD, ". \n"))
cat("The number of remaining features in each omic is: \n")
print(lungcrossomepostMultiout$NumFeat[paste0("d=", myD), ])

png(file = "lung_sample_size_vs_effect_size.png",width = 6,height = 6,units = "in",res = 600)
plot(lungcrossomepostMultiout$d, lungcrossomepostMultiout$SampleSize[, 1], type = "l", 
     lwd = 2, xlab = "Cohen's d cutoff", ylab = "Number of replicates", 
     main = "LUNG: Sample size vs Cohen's d", ylim = c(2, max(lungcrossomepostMultiout$SampleSize)))
arrows(x0 = 0, y0 = max.size, x1 = myD, y1 = max.size,lty = 2, col = 2)
arrows(x0 = myD, y0 = max.size, x1 = myD, y1 = 2, lty = 2,col = 2)
text(min(lungcrossomepostMultiout$d) + 1, max(lungcrossomepostMultiout$SampleSize, 
                                                na.rm = TRUE) - 2, paste0("Cohen's d = ", myD), col = 2)
dev.off()  

png(file = "lung_power_comparison.png",width = 6,height = 6,units = "in",res = 600)
barplot(lungcrossomepostMultiout$Power[c(1, min(which(lungcrossomepostMultiout$d >= myD))), ],
        col = rep(omicCol, each = 2), beside = TRUE, 
        las = 2, ylab = "Statistical power", ylim = c(0, 1), 
        border = rep(omicCol, each = 2), density = rep(c(30,100), ncol(lungcrossomepostMultiout$Power)),main = "LUNG: Power Comparison")
legend(x = 0.5, y = 1.0, c("Optimal SS", "User's SS"), col = 1, 
       density = c(30, 100), ncol = 2, bty = "n")
dev.off()


####
# brown

omicCol = colors()[c(554, 89, 111, 512, 17, 586, 132,428, 601, 568, 86, 390)]
omicCol = omicCol[1:ncol(browncrossomepostMultiout$SampleSize)]
names(omicCol) = colnames(browncrossomepostMultiout$SampleSize)
if (min(browncrossomepostMultiout$SampleSize) > max.size) {
  cat(paste0("The chosen sample size of ", max.size, " is not feasible. \n"))
  max.size = min(browncrossomepostMultiout$SampleSize)
  cat(paste0("A sample size of ", max.size, " will be plotted instead. \n"))
}

fff = approxfun(browncrossomepostMultiout$SampleSize[, 1], browncrossomepostMultiout$d)
myD = round(fff(max.size), 1)
cat(paste0("For having a sample size of ", max.size, 
           " and maintain the desired power, you need to remove features with Cohen's d below ", 
           myD, ". \n"))
cat("The number of remaining features in each omic is: \n")
print(browncrossomepostMultiout$NumFeat[paste0("d=", myD), ])

png(file = "brown_sample_size_vs_effect_size.png",width = 6,height = 6,units = "in",res = 600)
plot(browncrossomepostMultiout$d, browncrossomepostMultiout$SampleSize[, 1], type = "l", 
     lwd = 2, xlab = "Cohen's d cutoff", ylab = "Number of replicates", 
     main = "BAT: Sample size vs Cohen's d", ylim = c(2, max(browncrossomepostMultiout$SampleSize)))
arrows(x0 = 0, y0 = max.size, x1 = myD, y1 = max.size,lty = 2, col = 2)
arrows(x0 = myD, y0 = max.size, x1 = myD, y1 = 2, lty = 2,col = 2)
text(min(browncrossomepostMultiout$d) + 1, max(browncrossomepostMultiout$SampleSize, 
                                                na.rm = TRUE) - 2, paste0("Cohen's d = ", myD), col = 2)
dev.off()  

png(file = "brown_power_comparison.png",width = 6,height = 6,units = "in",res = 600)
barplot(browncrossomepostMultiout$Power[c(1, min(which(browncrossomepostMultiout$d >= myD))), ],
        col = rep(omicCol, each = 2), beside = TRUE, 
        las = 2, ylab = "Statistical power", ylim = c(0, 1), 
        border = rep(omicCol, each = 2), density = rep(c(30,100), ncol(browncrossomepostMultiout$Power)),main = "BAT: Power Comparison")
legend(x = 0.5, y = 1.0, c("Optimal SS", "User's SS"), col = 1, 
       density = c(30, 100), ncol = 2, bty = "n")
dev.off()


####
# white

omicCol = colors()[c(554, 89, 111, 512, 17, 586, 132,428, 601, 568, 86, 390)]
omicCol = omicCol[1:ncol(whitecrossomepostMultiout$SampleSize)]
names(omicCol) = colnames(whitecrossomepostMultiout$SampleSize)
if (min(whitecrossomepostMultiout$SampleSize) > max.size) {
  cat(paste0("The chosen sample size of ", max.size, " is not feasible. \n"))
  max.size = min(whitecrossomepostMultiout$SampleSize)
  cat(paste0("A sample size of ", max.size, " will be plotted instead. \n"))
}

fff = approxfun(whitecrossomepostMultiout$SampleSize[, 1], whitecrossomepostMultiout$d)
myD = round(fff(max.size), 1)
cat(paste0("For having a sample size of ", max.size, 
           " and maintain the desired power, you need to remove features with Cohen's d below ", 
           myD, ". \n"))
cat("The number of remaining features in each omic is: \n")
print(whitecrossomepostMultiout$NumFeat[paste0("d=", myD), ])

png(file = "white_sample_size_vs_effect_size.png",width = 6,height = 6,units = "in",res = 600)
plot(whitecrossomepostMultiout$d, whitecrossomepostMultiout$SampleSize[, 1], type = "l", 
     lwd = 2, xlab = "Cohen's d cutoff", ylab = "Number of replicates", 
     main = "WAT-SC: Sample size vs Cohen's d", ylim = c(2, max(whitecrossomepostMultiout$SampleSize)))
arrows(x0 = 0, y0 = max.size, x1 = myD, y1 = max.size,lty = 2, col = 2)
arrows(x0 = myD, y0 = max.size, x1 = myD, y1 = 2, lty = 2,col = 2)
text(min(whitecrossomepostMultiout$d) + 1, max(whitecrossomepostMultiout$SampleSize, 
                                                na.rm = TRUE) - 2, paste0("Cohen's d = ", myD), col = 2)
dev.off()  

png(file = "white_power_comparison.png",width = 6,height = 6,units = "in",res = 600)
barplot(whitecrossomepostMultiout$Power[c(1, min(which(whitecrossomepostMultiout$d >= myD))), ],
        col = rep(omicCol, each = 2), beside = TRUE, 
        las = 2, ylab = "Statistical power", ylim = c(0, 1), 
        border = rep(omicCol, each = 2), density = rep(c(30,100), ncol(whitecrossomepostMultiout$Power)),main = "WAT-SC: Power Comparison")
legend(x = 0.5, y = 1.0, c("Optimal SS", "User's SS"), col = 1, 
       density = c(30, 100), ncol = 2, bty = "n")
dev.off()

gastropowereffectcomp <- as.data.frame(gastrocrossomepostMultiout$Power)
gastropowereffectcomp$d <- as.numeric(gsub("d=","",rownames(gastropowereffectcomp)))
gastropowereffectcompdf <- data.frame("Power" = c(gastropowereffectcomp$RNAseq,
                                                  gastropowereffectcomp$ATACseq,
                                                  gastropowereffectcomp$RRBS,
                                                  gastropowereffectcomp$Prot_PR,
                                                  gastropowereffectcomp$Prot_PH),
                                      "d" = rep(gastropowereffectcomp$d,length(gastrostatdata)),
                                      "Ome" = c(rep("RNAseq",length(gastropowereffectcomp$d)),
                                                rep("ATACseq",length(gastropowereffectcomp$d)),
                                                rep("RRBS",length(gastropowereffectcomp$d)),
                                                rep("Prot_PR",length(gastropowereffectcomp$d)),
                                                rep("Prot_PH",length(gastropowereffectcomp$d))))


save.image("FigureS5_S6_S7_S8_112624.RData")