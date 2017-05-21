library(edgeR)
library(limma)

# Repeat ChIP-seq analysis using limma

countdata_filt<-read.table("data.txt.mappPE.mean0.8fl.txt", header=F)
header<-read.table("header.txt",header=F)

metapeaks.hg19.filt.raw<-read.table("comm.hg19.bed.ortho.hg19.mappinput.chrfl.mappPE.mean0.8.bed.mappfl.txt",header=F)
metapeaks.panTro3.filt.raw<-read.table("comm.hg19.bed.ortho.panTro3.mappinput.chrfl.mappPE.mean0.8.bed.mappfl.txt",header=F)

rownames(countdata_filt) <- metapeaks.hg19.filt.raw[,4]

colnames(countdata_filt)<-as.matrix((header[1,]))

colnames(metapeaks.hg19.filt.raw) <- c("Chr", "Start", "End", "Peak")
colnames(metapeaks.panTro3.filt.raw) <- c("Chr", "Start", "End", "Peak")


countdata_filt_noH20961 <- countdata_filt[,c(1,2,3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
                                             21,22,24,25,26,27,28,29,30,31,32,33,34)]

H3K9countdata_noH20961<-countdata_filt_noH20961[,c(2,5,6,7,8,10,13,15,17,18,19,20,23,26,27,32)]

chipcounts100 =  H3K9countdata_noH20961[apply( H3K9countdata_noH20961>0,1,sum) > 10, ]

# Make DGE list
ychip <- DGEList(chipcounts100)

# TMM normalization
ychip <- calcNormFactors(ychip, method="TMM")

# HvC

Ind_chip = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
condition_chip = c("H", "H", "H", "C", "C", "C", "C", "H", "H", "H", "C", "H", "H", "C", "C", "H")

design_chip <- model.matrix(~0 + condition_chip)

colnames(design_chip) <- c("C","H")

# Voom
v_chip <- voom(ychip, design_chip)

# Limma
fit_chip <- lmFit(v_chip, design_chip)

cm <- makeContrasts(
  HvC = H-C,
  levels=design_chip)

fit2_chip <-contrasts.fit(fit_chip, cm)
fit2_chip <-eBayes(fit2_chip)

HvC_chip = topTable(fit2_chip, coef=1, adjust="BH", number=Inf, sort.by="none")

sig.HvC_01_chip = HvC_chip[HvC_chip$adj.P.Val < .01 , ]

write.table(HvC_chip,file='HvC_ChIP_limma.txt'
            , quote=F,row.names=T, col.names=T, sep="\t")
