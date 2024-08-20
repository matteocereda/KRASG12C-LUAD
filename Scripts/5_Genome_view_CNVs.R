options(stringsAsFactors=F)
library(readr)
library(ChIPpeakAnno)
library(karyoploteR)
library(dplyr)
library(grid)


cnv = readRDS("Rdata/somatic_CNV.rds")
colnames(cnv)[ grep('patient', colnames(cnv)) ] ='sample'

samples=readRDS("Rdata/metadata_samples.rds")

cnv$sample = factor(cnv$sample, levels=levels(samples$sample))
cnv$barcode  = factor(samples$barcode[match(cnv$sample,samples$sample)], levels=levels(samples$barcode))

dcnv = subset(cnv, selected & ProbCall>.75 & !is.na(seq_cnv))

dcnv_n933 = subset(cnv, selected & gene_name %in% unique(dcnv$gene_name))

genes = read_delim("Supporting_files/GRanges_pc_genes_190218.bed", col_names = F)[,c(1,2,3,7)]
colnames(genes) = c("chromosome", "start", "end", "gene_name")

dcnv = unique(rbind.data.frame(dcnv, dcnv_n933))
dcnv$chromosome=genes$chromosome[match(dcnv$gene_name, genes$gene_name)] 
dcnv$start=genes$start[match(dcnv$gene_name, genes$gene_name)] 
dcnv$end=genes$end[match(dcnv$gene_name, genes$gene_name)] 
 

subset_chrs = unique(dcnv$chromosome)
subset_chrs = subset_chrs[ subset_chrs !='chrY' ]

cnv_genes = unique(subset(dcnv, chromosome%in%subset_chrs & grepl('lung',primary_site))[,c("chromosome", "start", "end", "gene_name")])

oncogenes = unique(subset(dcnv,
                          chromosome%in%subset_chrs
                          & grepl('lung',primary_site)
                          &
                            ( ( grepl('oncogene',cgc_annotation ) & grepl('Amp',exc_cnv))
                              )
                          )[,c("chromosome", "start", "end", "gene_name")])

oncogenes =toGRanges(oncogenes)


tsg = unique(subset(dcnv,
                          chromosome%in%subset_chrs
                          & grepl('lung',primary_site)
                          &
                            ( ( grepl('TSG',cgc_annotation ) & grepl('Loss',exc_cnv))
                            )
)[,c("chromosome", "start", "end", "gene_name")])

tsg =toGRanges(tsg)


x  =  plyr::dlply(dcnv, ~barcode)
max_y = round(max(dcnv$CNF, na.rm = T),0)



get_genome_view_genes <- function(x, oncogenes=NULL, tsg=NULL, title=NULL, subset_chrs=NULL, genome='hg19', EXCAVATOR=T) {

  kp = NULL
  if( is.null(subset_chrs)) {
    kp = plotKaryotype(genome, plot.type = 4, ideogram.plotter = NULL)
    kpAddCytobandsAsLine(kp)
  }else{
    kp=plotKaryotype(genome, plot.type = 4, ideogram.plotter = NULL, chromosomes = subset_chrs)
    kpAddCytobands(kp)
  }

  if(!is.null(oncogenes))  suppressWarnings(kpPlotMarkers(kp, data=oncogenes, labels=oncogenes$gene_name,
                                                          marker.parts = c(0.95,0.025,0.025),
                                                          label.color = 'red', r1 = 1.5))

  if(!is.null(tsg))  suppressWarnings(kpPlotMarkers(kp, data=tsg, labels=tsg$gene_name,
                                                    marker.parts = c(0.95,0.025,0.025), label.color = 'blue', r1 = 1.5))

  for(i in 1:length(x)){
    at <- autotrack(current.track = i, total.tracks = length(x))
    seq   = x[[i]]


    #c("Hom.Loss", "Het.Loss", "Amp", "HL.Ampl")

    if(EXCAVATOR){
      kpAxis(kp, data.panel = 1, ymin = 0, ymax=max_y, numticks = (max_y)+1, side = 1, cex = 0.8, r0=at$r0, r1=at$r1 )
      kpAddLabels(kp, labels = names(x)[i], r0=at$r0, r1=at$r1, srt=90, pos=1, label.margin = 0.08, cex = 0.8)
      kpDataBackground(kp, col="grey98", r0=at$r0, r1=at$r1)
      for(j in 1:max_y) kpAbline(kp, h=j/max_y,r0=at$r0, r1=at$r1, col='grey90')

      A <- with(subset(seq, grepl('Amp',seq$exc_cnv)),
                toGRanges(data.frame(chr=chromosome, start=start, end=end, x=start, y= CNF/max_y  )))
      seqlevelsStyle(A) <- "UCSC"
      kpPoints(kp, data=A,  pch=1, col="red", r0=at$r0, r1=at$r1, cex = 1 )

      A <- with(subset(seq, grepl('Loss',seq$exc_cnv)),
                toGRanges(data.frame(chr=chromosome, start=start, end=end, x=start, y= CNF/max_y )))
      seqlevelsStyle(A) <- "UCSC"
      kpPoints(kp, data=A, pch=1, col="blue", r0=at$r0, r1=at$r1, cex = 1)

    }else{

      seq$max = rowMax(as.matrix(seq[,c("A","B")]))
      seq$min = rowMin(as.matrix(seq[,c("A","B")]))

      max_y = round(max(seq$max),0)

      kpAxis(kp, data.panel = 1, ymin = 0, ymax=max_y, numticks = (max_y)+1, side = 1, cex = 0.8, r0=at$r0, r1=at$r1 )
      kpAddLabels(kp, labels = names(x)[i], r0=at$r0, r1=at$r1, srt=90, pos=1, label.margin = 0.08, cex = 0.8)
      kpDataBackground(kp, col="grey98", r0=at$r0, r1=at$r1)
      for(j in 1:max_y) kpAbline(kp, h=j/max_y,r0=at$r0, r1=at$r1, col='grey90')


      A <- with(subset(seq, grepl('Amp',seq$seq_cnv)),
                toGRanges(data.frame(chr=chromosome, start=start, end=end, x=start, y= max/max_y  )))
      seqlevelsStyle(A) <- "UCSC"
      kpPoints(kp, data=A,  pch=1, col="red", r0=at$r0, r1=at$r1, cex = 1 )


      A <- with(subset(seq, grepl('Loss',seq$seq_cnv)),
                toGRanges(data.frame(chr=chromosome, start=start, end=end, x=start, y= min/max_y )))
      seqlevelsStyle(A) <- "UCSC"
      kpPoints(kp, data=A, pch=1, col="blue", r0=at$r0, r1=at$r1, cex = 1)

    }



  }
}



pdf(file='Figures/genome_view_genes.pdf', height = unit(5, 'cm'), width=unit(10,'cm'), useDingbats = F)
get_genome_view_genes( x, oncogenes = oncogenes, tsg = tsg, subset_chrs = subset_chrs)
dev.off()
