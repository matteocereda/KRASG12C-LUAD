options(stringsAsFactors=F)

library(clonevol)
library(grid)


# ClonEvol ==============

x2=read.delim2("PyClone/loci.tsv")
x2$sample_id=factor(x2$sample_id,levels = c("ne933", "nh511", "xeno1", "dsf3c") )
x2$cellular_prevalence=as.numeric(x2$cellular_prevalence)
x2$variant_allele_frequency=as.numeric(x2$variant_allele_frequency)*100
x2$cellular_prevalence_std=as.numeric(x2$cellular_prevalence_std)


x=reshape2::dcast(x2, mutation_id~sample_id, value.var = "variant_allele_frequency")

colnames(x)=c("id", paste0(colnames(x)[2:ncol(x)], ".vaf"))

x$cluster=x2$cluster_id[match(x$id, x2$mutation_id)]

vaf.col.names <- grep('.vaf', colnames(x), value=T)

x <- x[order(x$cluster),]

x$cluster=x$cluster+1


## 01 filter clusters =========================

l=as.data.frame(table(x$cluster))
l=subset(l, Freq>=10)

x3=subset(x, cluster%in%l$Var1)

k=1
for(i in unique(x3$cluster)){
  
  x3$cluster[which(x3$cluster==i)]=k
  
  k=k+1
  
}


df=readRDS("Rdata/somatic_mutations.rds")
x3$gene=df$Gene.refGene[match(x3$id, df$id)]

x3$cancervar=df$CancerVar_res[match(x3$id, df$id)]

figc=read.csv("Tables/Mutations_Figure_1C_annot.csv")

x3$is_driver=x3$id%in%figc$id


clone.colors=c('#30D5C8', "#A5D9FD" ,"#F0027F" ,"#984EA3" ,"#FBE928" ,"#77DD77", "#E78AC3")

x4=reshape2::melt(x3, id.vars = c("id", "cluster","is_driver", "gene", "cancervar"))


x4=ddply(x4, .(cluster, variable), mutate, mean_clus=mean(value))

df$key=paste0(df$id, "_", df$patient)

x4$patient=sapply(strsplit(as.character(x4$variable), "\\."), "[[",1)

x4$key=paste0(x4$id, "_", x4$patient)

x4$selected=df$selected[match(x4$key, df$key)]
x4$selected[is.na(x4$selected)]<-FALSE

x4$patient=factor(x4$patient, levels = c("ne933", "nh511", "xeno1", "dsf3c"))

x4$to_plot=NA
x4$to_plot[which(x4$id%in%figc$id)]=x4$gene[which(x4$id%in%figc$id)]


## 02 infer model =========================

y2 = infer.clonal.models(variants = x3,
                         cluster.col.name = 'cluster',
                         vaf.col.names = vaf.col.names,
                         sample.names = vaf.col.names,
                         cancer.initiation.model='monoclonal',
                         subclonal.test = 'bootstrap',
                         subclonal.test.model = 'non-parametric',
                         num.boots = 10000,
                         founding.cluster = 1,
                         cluster.center = 'median',
                         ignore.clusters = NULL,
                         clone.colors = clone.colors,
                         min.cluster.vaf = 0.01,
                         # min probability that CCF(clone) is non-negative
                         sum.p = 0.05,
                         # alpha level in confidence interval estimate for CCF(clone)
                         alpha = 0.05,
                         random.seed=18494)


w2<- convert.consensus.tree.clone.to.branch(y2, branch.scale = 'sqrt')
w2<- transfer.events.to.consensus.trees(w2,
                                        x3[x3$is_driver,],
                                        cluster.col.name = 'cluster',
                                        event.col.name = 'gene')

## 03 plot  =========================

plot.clonal.models(w2,
                   # box plot parameters
                   box.plot = TRUE,
                   fancy.boxplot = TRUE,
                   fancy.variant.boxplot.highlight = 'is_driver',
                   fancy.variant.boxplot.highlight.shape = 21,
                   fancy.variant.boxplot.highlight.fill.color = 'red',
                   fancy.variant.boxplot.highlight.color = 'black',
                   fancy.variant.boxplot.highlight.note.col.name = 'gene',
                   fancy.variant.boxplot.highlight.note.color = 'blue',
                   fancy.variant.boxplot.highlight.note.size = 2,
                   fancy.variant.boxplot.jitter.alpha = 1,
                   fancy.variant.boxplot.jitter.center.color = 'grey50',
                   fancy.variant.boxplot.base_size = 12,
                   fancy.variant.boxplot.plot.margin = 1,
                   fancy.variant.boxplot.vaf.suffix = '.VAF',
                   # bell plot parameters
                   clone.shape = 'bell',
                   bell.event = TRUE,
                   bell.event.label.color = 'blue',
                   bell.event.label.angle = 60,
                   clone.time.step.scale = 1,
                   bell.curve.step = 2,
                   # node-based consensus tree parameters
                   merged.tree.plot = TRUE,
                   tree.node.label.split.character = NULL,
                   tree.node.shape = 'circle',
                   tree.node.size = 30,
                   tree.node.text.size = 0.5,
                   merged.tree.node.size.scale = 1.25,
                   merged.tree.node.text.size.scale = 2.5,
                   merged.tree.cell.frac.ci = FALSE,
                   # branch-based consensus tree parameters
                   merged.tree.clone.as.branch = TRUE,
                   mtcab.event.sep.char = ',',
                   mtcab.branch.text.size = 1,
                   mtcab.branch.width = 0.75,
                   mtcab.node.size = 3,
                   mtcab.node.label.size = 1,
                   mtcab.node.text.size = 1.5,
                   # cellular population parameters
                   cell.plot = TRUE,
                   num.cells = 100,
                   cell.border.size = 0.25,
                   cell.border.color = 'black',
                   clone.grouping = 'horizontal',
                   #meta-parameters
                   scale.monoclonal.cell.frac = TRUE,
                   show.score = FALSE,
                   cell.frac.ci = TRUE,
                   disable.cell.frac = FALSE,
                   # output figure parameters
                   out.dir = 'Figures',
                   out.format = 'pdf',
                   overwrite.output = TRUE,
                   width = 8,
                   height = 4,
                   # vector of width scales for each panel from left to right
                   panel.widths = c(3,4,2,4,2)
)