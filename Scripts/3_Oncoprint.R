options(stringsAsFactors=F)
library(plyr)
library(dplyr)
library(viridis)



suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  FS = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "yellow", col = NA))
  },
  NFS = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "orange", col = NA))
  },
  SY = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "pink", col = NA))
  },

  NS = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "blue", col = NA))
  },
  SG = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "red", col = NA))
  },
  SL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "brown", col = NA))
  },

  SP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "green", col = NA))
  }
)

col = c("FS" = "yellow","NFS"='orange',"SY"="pink", "NS" = "blue", "SG" = "red","SL"='brown', "SP" = "green")

mc = c("frameshift substitution"="FS","nonframeshift substitution"="NFS"
       ,"nonsynonymous SNV"='NS',"splicing"='SP',"stopgain"='SG',"stoploss"='SL'
       ,"synonymous SNV"='SY')

get_oncoprint <- function(y, min=NULL) {
  a = y$hgnc_canonical_refseq
  if(sum(a=='')>0) a[which(a=='')] = y$alternative_refseq[which(a=='')]
  a = sapply(strsplit(a, "\\,"), "[[", 1)

  y$key = sapply(strsplit(a, "\\:"), function(x) paste0(x[1],"_",x[length(x)]) )

  message("[*] Counting mutations ...")
  tmp =ddply(y, .(barcode, key), summarise
             , value=paste(unique(ExonicFunc.refGene), collapse=";"))

  suppressWarnings(suppressPackageStartupMessages(library(reshape)))
  mat = cast(tmp, key ~ barcode )
  mat[is.na(mat)] = ""
  rownames(mat) = mat[, 1]
  mat = mat[, -1]
  mat = as.data.frame(mat)

  if(!is.null(min)){
    s=apply(mat,1, function(x) sum(x!=""))
    mat = mat[s>min,]
  }


  vaf = cast(y, key ~ barcode, value = "VAF" )
  vaf[is.na(vaf)] = 0
  rownames(vaf) = vaf[, 1]
  vaf = vaf[, -1]
  vaf = as.data.frame(vaf)


  library(reshape2)

  message("[*] Setting gene categories ...")

  ann = unique(y[,c("key","CancerVar_res","oncoKB",'damaging','cgc_annotation','primary_site')])
  rownames(ann) = ann[,1]; ann[,1] = NULL

  message("[*] Oncoprint ...")

  pal_ann = ggsci::pal_d3()(7)
  oncoPrint(mat[, sort(y$barcode) %>% unique ],
            get_type = function(x) strsplit(x, ";")[[1]],
            alter_fun = alter_fun,
            col = col
            , column_order = sort(y$barcode) %>% unique %>% as.character,
            show_column_names = T,
            row_title_gp= gpar(fontsize=8), column_title_gp= gpar(fontsize=8)
            ,row_names_gp = gpar(fontsize=6), column_names_gp = gpar(fontsize=8)
            ,column_title = ""
            , width = unit(.5*ncol(mat), "cm")
            , pct_gp = gpar(fontsize=8)
            , heatmap_legend_param = list(title = "Alternations",
                                         at = mc,
                                         labels = names(mc)
                                         ,ncol = 1)
            , top_annotation = HeatmapAnnotation(column_bar = anno_oncoprint_barplot(),
                                                 annotation_height = unit(1, "cm"))

  )+HeatmapAnnotation(df=ann[rownames(mat),], width = unit(2, "cm")
                      , which='row'
                      , gp = gpar(col = 'white')
                      , show_annotation_name=T
                      , annotation_name_gp = gpar(fontsize=8)
                      # , col=list(
                      #   actionable=c('1'=pal_ann[1],'0'='grey')
                      #   , cgc=c('1'=pal_ann[2],'0'='grey')
                      #   , sarcoma_CG=c('1'=pal_ann[3],'0'='grey')
                      #   , pediatric_CG=c('1'=pal_ann[4],'0'='grey')
                      #   , pediatric_sarcoma_CG=c('1'=pal_ann[5],'0'='grey')
                      #   , fusion=c('1'=pal_ann[6],'0'='grey')
                      #   , can=c('1'=pal_ann[7],'0'='grey'))
                      , show_legend = T
  )+Heatmap(vaf[rownames(mat),sort(y$barcode) %>% unique], cluster_rows = F, cluster_columns = F,
            width = unit(.5*ncol(mat), "cm")
            ,col=viridis(100),

            show_row_names = F
            ,column_names_gp = gpar(fontsize=8)
            , rect_gp = gpar(color='white')
            )
}


mut = readRDS("Rdata/somatic_mutations.rds")

samples=readRDS("Rdata/metadata_samples.rds")

mut$sample   = factor(mut$sample, levels=levels(samples$sample))
mut$barcode  = factor(mut$barcode, levels=levels(samples$barcode))
dmut = subset(mut, CancerVar_res %in% c('I_strong',"II_potential") | (CancerVar_res=="III_Uncertain" & selected) )

pdf(file="Figures/oncoprint.pdf", height = unit(10, 'cm'), width = unit(8,'cm'), useDingbats = F)
get_oncoprint(dmut)
dev.off()
