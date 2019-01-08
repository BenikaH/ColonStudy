library(biclust)
setwd("/Users/BC/Documents/clinical/")
clinical <- read.table("nationwidechildrens.org_clinical_patient_coad.txt", header=T, sep="\t")
stageI <- read.table("Stage_I_IA.txt", header=T, sep="\t")
stageII <- read.table("Stage_II_ALL.txt", header=T, sep="\t")
stageIII <- read.table("Stage_III_ALL.txt", header=T, sep="\t")
stageIV <- read.table("Stage_IV_ALL.txt", header=T, sep="\t")
miRNA <- read.csv("miRNASeq_list.csv", header=T)
RNASeq <- read.csv("RNASeqV2_list.csv", header=T)

stageIcodes <- stageI$barcode
stageIIcodes <- stageII$barcode
stageIIIcodes <- stageIII$barcode
stageIVcodes <- stageIV$barcode

stageIsamples <- merge(stageI, RNASeq, by.x="barcode", by.y="barcode")
stageIIsamples <- merge(stageII, RNASeq, by.x="barcode", by.y="barcode")
stageIIIsamples <- merge(stageIII, RNASeq, by.x="barcode", by.y="barcode")
stageIVsamples <- merge(stageIV, RNASeq, by.x="barcode", by.y="barcode")

write.table(stageIVsamples, file="stageIV_RNASeq_samples.csv", row.names=F, sep=",")

stageI_miRNA_samples <- merge(stageI, miRNA, by.x="barcode", by.y="barcode")
stageII_miRNA_samples <- merge(stageII, miRNA, by.x="barcode", by.y="barcode")
stageIII_miRNA_samples <- merge(stageIII, miRNA, by.x="barcode", by.y="barcode")
stageIV_miRNA_samples <- merge(stageIV, miRNA, by.x="barcode", by.y="barcode")

write.table(stageIV_miRNA_samples, file="stageIV_miRNA_samples.csv", row.names=F, sep=",")


stageI_miRNA_and_RNASeq_samples <- merge(stageI_miRNA_samples, stageIsamples, by.x="barcode", by.y="barcode")
stageII_miRNA_and_RNASeq_samples <- merge(stageII_miRNA_samples, stageIIsamples, by.x="barcode", by.y="barcode")
stageIII_miRNA_and_RNASeq_samples <- merge(stageIII_miRNA_samples, stageIIIsamples, by.x="barcode", by.y="barcode")
stageIV_miRNA_and_RNASeq_samples <- merge(stageIV_miRNA_samples, stageIVsamples, by.x="barcode", by.y="barcode")

write.table(stageIV_miRNA_and_RNASeq_samples[,1], file="stageIV_miRNA_and_RNASeq_list.csv", row.names=F, sep=",")





