options(stringsAsFactors=F)

# 0. Prepare file ============
anno_b=read.delim2("EaCon_output/DSFC3_B.TargetGenes.txt")
anno_b$sampleID="DSFC3_B"
anno_br=read.delim2("EaCon_output/DSFC3_BR.TargetGenes.txt")
anno_br$sampleID="DSFC3_BR"
anno_c=read.delim2("EaCon_output/DSFC3_C.TargetGenes.txt")
anno_c$sampleID="DSFC3_C"
anno_cr=read.delim2("EaCon_output/DSFC3_CR.TargetGenes.txt")
anno_cr$sampleID="DSFC3_CR"

all=rbind.data.frame(anno_b, anno_br, anno_c, anno_cr)
all$clonality=abs(as.numeric(all$BAF.Value)*100-50)*2



all$CNV=2*2^as.numeric(all$L2R.Value)

#saveRDS(all, "Rdata/Copy_numbers_all_samples_genes.rds")


# 1. Karyoplots ============

all=readRDS("Rdata/Copy_numbers_all_samples_genes.rds")
library(karyoploteR)



l=read.delim2("EaCon_output/DSFC3_B.SegmentedBAF.txt")
l2=read.delim2("EaCon_output/DSFC3_BR.SegmentedBAF.txt")
l3=read.delim2("EaCon_output/DSFC3_C.SegmentedBAF.txt")
l4=read.delim2("EaCon_output/DSFC3_CR.SegmentedBAF.txt")


p=readRDS("EaCon_output/DSFC3_B_OncoScan_CNV_hg19_processed.RDS")
p=p$data

p2=readRDS("EaCon_output/DSFC3_BR_OncoScan_CNV_hg19_processed.RDS")
p2=p2$data

p3=readRDS("EaCon_output/DSFC3_C_OncoScan_CNV_hg19_processed.RDS")
p3=p3$data

p4=readRDS("EaCon_output/DSFC3_CR_OncoScan_CNV_hg19_processed.RDS")
p4=p4$data


chr=p$SNPpos

baf=p$Tumor_BAF
logr=p$Tumor_LogR
head(logr)


baf2=p2$Tumor_BAF
logr2=p2$Tumor_LogR
head(logr2)

baf3=p3$Tumor_BAF
logr3=p3$Tumor_LogR
head(logr3)

baf4=p4$Tumor_BAF
logr4=p4$Tumor_LogR
head(logr4)


chr=cbind.data.frame(chr, baf$DSFC3_B[match(rownames(chr), rownames(baf))],
                     baf2$DSFC3_BR[match(rownames(chr), rownames(baf2))],
                     baf3$DSFC3_C[match(rownames(chr), rownames(baf3))],
                     baf4$DSFC3_CR[match(rownames(chr), rownames(baf4))],
                     logr$DSFC3_B[match(rownames(chr), rownames(logr))],
                     logr2$DSFC3_BR[match(rownames(chr), rownames(logr2))],
                     logr3$DSFC3_C[match(rownames(chr), rownames(logr3))],
                     logr4$DSFC3_CR[match(rownames(chr), rownames(logr4))])

colnames(chr)=c("chr", "pos", "BAF_DSFC3_B","BAF_DSFC3_BR","BAF_DSFC3_C","BAF_DSFC3_CR", "LR_DSFC3_B","LR_DSFC3_BR","LR_DSFC3_C","LR_DSFC3_CR")
head(chr)



q=readRDS("EaCon_output/DSFC3_B.ASCN.SEQUENZA.RDS")
q=q$data
q$chrom=paste0("chr", q$chrom)

q2=readRDS("EaCon_output/DSFC3_BR.ASCN.SEQUENZA.RDS")
q2=q2$data
q2$chrom=paste0("chr", q2$chrom)

q3=readRDS("EaCon_output/DSFC3_C.ASCN.SEQUENZA.RDS")
q3=q3$data
q3$chrom=paste0("chr", q3$chrom)

q4=readRDS("EaCon_output/DSFC3_CR.ASCN.SEQUENZA.RDS")
q4=q4$data
q4$chrom=paste0("chr", q4$chrom)



chr12=subset(chr, chr=="chr12")


l_12=subset(l, Chrom=="chr12")
l2_12=subset(l2, Chrom=="chr12")
l3_12=subset(l3, Chrom=="chr12")
l4_12=subset(l4, Chrom=="chr12")


q_12=subset(q, chrom=="chr12")

q2_12=subset(q2, chrom=="chr12")

q3_12=subset(q3, chrom=="chr12")

q4_12=subset(q4, chrom=="chr12")


# 1.1 KRAS ===================


pdf("Figures/KRAS_BAF_LR.pdf", useDingbats = F)
kp <- plotKaryotype(chromosomes = c("chr12"), plot.type = 4, genome="hg19")
kpAddBaseNumbers(kp)
kpPlotMarkers(kp, chr="chr12", x=25358180, labels="KRAS",y = -0.22)

kpDataBackground(kp, data.panel = 1, col="white", r0=0,r1=0.20)
kpDataBackground(kp, data.panel = 2, col="white", r0=0.21,r1=0.41)

kpSegments(kp, chr="chr12", x0=0, x1=max(chr12$pos), r0=0,r1=0.20, y0=0.5, y1=0.5, lty=2, lwd=0.5, col="black")
kpSegments(kp, chr="chr12", x0=0, x1=max(chr12$pos), r0=0.21,r1=0.41, y0=0.5, y1=0.5, lty=2, lwd=0.5, col="black")

kpPoints(kp, chr="chr12", x=chr12$pos, y=chr12$BAF_DSFC3_B, data.panel = 1, r0=0,r1=0.20,col = "#fdc500", cex=0.25)
kpPoints(kp, chr="chr12", x=chr12$pos, y=chr12$BAF_DSFC3_BR, data.panel = 1, r0=0,r1=0.20,col = "#7209b7", cex=0.25)
kpPoints(kp, chr="chr12", x=chr12$pos, y=chr12$BAF_DSFC3_C, data.panel = 1,r0=0.21,r1=0.41,col = "#fdc500", cex=0.25)
kpPoints(kp, chr="chr12", x=chr12$pos, y=chr12$BAF_DSFC3_CR, data.panel = 1,r0=0.21,r1=0.41,col = "#7209b7", cex=0.25)


kpRect(kp, chr="chr12", x0=0, x1=max(chr12$pos), r0=0,r1=0.20, y0 = -0.01, y1=1.01)
kpRect(kp, chr="chr12", x0=0, x1=max(chr12$pos), r0=0.21,r1=0.41, y0 = -0.01, y1=1.01)

kpAxis(kp, data.panel=1, r0=0,r1=0.20)
kpAxis(kp, data.panel=2, r0=0.21,r1=0.41)

kpSegments(kp, chr="chr12", x0=l_12$Start, x1=l_12$End,data.panel = 1, r0=0,r1=0.20, y0=as.numeric(l_12$BAF.Value), y1=as.numeric(l_12$BAF.Value),lwd=1.5, col="#fdc500")

kpSegments(kp, chr="chr12", x0=l2_12$Start, x1=l2_12$End,data.panel = 1, r0=0,r1=0.20, y0=as.numeric(l2_12$BAF.Value), y1=as.numeric(l2_12$BAF.Value),lwd=1.5, col="#7209b7")

kpSegments(kp, chr="chr12", x0=l3_12$Start, x1=l3_12$End,data.panel = 1, r0=0.21,r1=0.41, y0=as.numeric(l3_12$BAF.Value), y1=as.numeric(l3_12$BAF.Value),lwd=1.5, col="#fdc500")

kpSegments(kp, chr="chr12", x0=l4_12$Start, x1=l4_12$End,data.panel = 1,r0=0.21,r1=0.41, y0=as.numeric(l4_12$BAF.Value), y1=as.numeric(l4_12$BAF.Value),lwd=1.5, col="#7209b7")


kpDataBackground(kp, data.panel = 1, col="white", r0=0.42,r1=0.63)
kpDataBackground(kp, data.panel = 2, col="white", r0=0.64,r1=0.84)

mi=min(c(min(chr12$LR_DSFC3_B,na.rm = T),min(chr12$LR_DSFC3_BR,na.rm = T),min(chr12$LR_DSFC3_C,na.rm = T),min(chr12$LR_DSFC3_CR,na.rm = T)))
ma=max(c(max(chr12$LR_DSFC3_B,na.rm = T),max(chr12$LR_DSFC3_BR,na.rm = T),max(chr12$LR_DSFC3_C,na.rm = T),max(chr12$LR_DSFC3_CR,na.rm = T)))
mi=-1
ma=2



kpSegments(kp, chr="chr12", x0=0, x1=max(chr12$pos), r0=0.42,r1=0.63, y0=1/3, y1=1/3, lty=2, lwd=0.5, col="black")
kpSegments(kp, chr="chr12", x0=0, x1=max(chr12$pos), r0=0.64,r1=0.84, y0=1/3, y1=1/3, lty=2, lwd=0.5, col="black")


kpPoints(kp, chr="chr12", x=chr12$pos, y=chr12$LR_DSFC3_B, data.panel = 1,  r0=0.42,r1=0.63,col = "#fdc500", ymax=ma, ymin=mi, cex=0.25)
kpPoints(kp, chr="chr12", x=chr12$pos, y=chr12$LR_DSFC3_BR, data.panel = 1, r0=0.42,r1=0.63,col = "#7209b7", ymax=ma, ymin=mi, cex=0.25)


kpPoints(kp, chr="chr12", x=chr12$pos, y=chr12$LR_DSFC3_C, data.panel = 1,r0=0.64,r1=0.84,col = "#fdc500", ymax=ma, ymin=mi, cex=0.25)
kpPoints(kp, chr="chr12", x=chr12$pos, y=chr12$LR_DSFC3_CR, data.panel = 1,r0=0.64,r1=0.84,col = "#7209b7", ymax=ma, ymin=mi, cex=0.25)


kpRect(kp, chr="chr12", data.panel=1, x0=0, x1=max(chr12$pos), r0=0.42,r1=0.63, y0 = 0-0.01, y1=1+0.01)
kpRect(kp, chr="chr12", data.panel=2, x0=0, x1=max(chr12$pos),r0=0.64,r1=0.84, y0 = 0-0.01, y1=1+0.01)


kpAxis(kp, data.panel=1, r0=0.42,r1=0.63, ymax=ma, ymin=mi,tick.pos = c(mi, 0, ma))
kpAxis(kp, data.panel=2, r0=0.64,r1=0.84, ymax=ma, ymin=mi,tick.pos = c(mi, 0, ma))

kpSegments(kp, chr="chr12",  x0=q_12$start.pos, x1=q_12$end.pos, y0=q_12$logR.mean/3+1/3, y1=q_12$logR.mean/3+1/3,lwd=1.5, col="#fdc500",r0=0.42,r1=0.63, data.panel=1)
kpSegments(kp, chr="chr12",  x0=q2_12$start.pos, x1=q2_12$end.pos, y0=q2_12$logR.mean/3+1/3, y1=q2_12$logR.mean/3+1/3,lwd=1.5, col="#7209b7",r0=0.42,r1=0.63, data.panel=1)

kpSegments(kp, chr="chr12",  x0=q3_12$start.pos, x1=q3_12$end.pos, y0=q3_12$logR.mean/3+1/3, y1=q3_12$logR.mean/3+1/3,lwd=1.5, col="#fdc500",r0=0.64,r1=0.84, data.panel=1)
kpSegments(kp, chr="chr12",  x0=q4_12$start.pos, x1=q4_12$end.pos, y0=q4_12$logR.mean/3+1/3, y1=q4_12$logR.mean/3+1/3,lwd=1.5, col= "#7209b7",r0=0.64,r1=0.84, data.panel=1)



kpPoints(kp, chr="chr12", x=25358180,y =0.2, r0=0.86,r1=1 , cex=subset(all, Target.Symbol=="KRAS" & sampleID=="DSFC3_B")$CNV/2, col="#00a6fb")

kpPoints(kp, chr="chr12", x=25358180,y =0.4, r0=0.86,r1=1 , cex=subset(all, Target.Symbol=="KRAS" & sampleID=="DSFC3_BR")$CNV/2, col="#00a6fb")

kpPoints(kp, chr="chr12", x=25358180,y =0.6, r0=0.86,r1=1 , cex=subset(all, Target.Symbol=="KRAS" & sampleID=="DSFC3_C")$CNV/2, col="#00a6fb")

kpPoints(kp, chr="chr12", x=25358180,y =0.8, r0=0.86,r1=1 , cex=subset(all, Target.Symbol=="KRAS" & sampleID=="DSFC3_CR")$CNV/2, col="#00a6fb")



kpPoints(kp, chr="chr12", x=40000000,y =0.2, r0=0.86,r1=1 , cex=1/2, col="black")
kpPoints(kp, chr="chr12", x=40000000,y =0.4, r0=0.86,r1=1 , cex=2/2, col="black")
kpPoints(kp, chr="chr12", x=40000000,y =0.6, r0=0.86,r1=1 , cex=3/2, col="black")
kpPoints(kp, chr="chr12", x=40000000,y =0.8, r0=0.86,r1=1 , cex=4/2, col="black")
dev.off()





