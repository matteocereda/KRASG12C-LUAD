options(stringsAsFactors=F)
library(readr)
library(ggplot2)
library(data.table)
library(dplyr)
library(grid)


mut = readRDS("Rdata/somatic_mutations.rds")

samples=readRDS("Rdata/metadata_samples.rds")

mut$sample   = factor(mut$sample, levels=levels(samples$sample))
mut$barcode  = factor(mut$barcode, levels=levels(samples$barcode))

report = read_csv("Tables/Mutations_Figure_1C.csv")
moi = subset(mut, Gene.refGene %in% report$Gene)
moi$hgnc_canonical_refseq[moi$hgnc_canonical_refseq==""] = moi$alternative_refseq[moi$hgnc_canonical_refseq==""]
moi = subset(moi, grepl(paste(unique(report$Mut), collapse="|"), moi$hgnc_canonical_refseq))
moi$key = paste0(moi$Gene.refGene, '_', sapply(strsplit(moi$hgnc_canonical_refseq, 'p.'), last))
report$key = paste0(report$Gene, '_', report$Mut)

moi = subset(moi, key%in%report$key)
ids = unique(moi[,c('id','key')])
report$ids=ids$ids[match(report$key, ids$key)]

#write.csv(ids, "Tables/Mutations_Figure_1C_annot.csv", row.names = F)

# write.csv(
#  moi[,c('key', 'barcode', 'ExonicFunc.refGene', "DP","FREQ","VAF", "CancerVar_res")] %>%
#  arrange(key, barcode)
#  , file='Tables/somatic_variants_Figure1C.csv'
# )

pdf(file="Figures/validations_figure1C.pdf",height = unit(6,'cm'), width = unit(6,'cm'), useDingbats = F)
ggplot(moi, aes(x=barcode, y=VAF, group=key, fill=CancerVar_res))+
  geom_line(aes(y=FREQ), lty='dotted', col='grey50')+
  geom_point(aes(y=FREQ, size=DP) ,shape=21, fill='white', col='grey50')+
  geom_line()+
  geom_point(aes(size=DP),shape=21)+
  facet_wrap(~key, ncol =3)+theme_bw()+theme(axis.text.x = element_text(angle=90))
dev.off()

