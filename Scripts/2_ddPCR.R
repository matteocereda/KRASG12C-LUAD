options(stringsAsFactors=F)
library(grid)
library(ggplot2)


# ddPCR ============

dna=read.csv("ddPCR/gDNA_fraction_G12C_sonda nuova_Analisi_CEREDA.csv")
dna$which="gDNA"

rna=read.csv("ddPCR/cDNA_fraction_G12C_sonda nuova_Analisi_CEREDA.csv")
rna$which="cDNA"
rna=subset(rna, ng==1.5)
rna$ng<-NULL

colnames(dna)=colnames(rna)

all=rbind.data.frame(dna, rna)
all$type="G"
all$type[grep("R", all$Sample.description)]="R"

pdf("Figures/Barplot_qqPCR_G12C.pdf", useDingbats = F, height=unit(5, "cm"), width=unit(4.5, "cm"))
ggplot(all, aes(x=Sample.description, y=Fractional.Abundance, fill=type,color=type))+geom_bar(stat="identity", alpha=0.6)+theme_bw()+
  scale_fill_manual(values=c("G"="#fdc500", "R"="#7209b7"))+
  scale_color_manual(values=c("G"="#fdc500", "R"="#7209b7"))+
  xlab("Sample")+ylab("Fraction G12C")+
  geom_text(aes(label=paste0(Fractional.Abundance, " %")), color="white")+facet_wrap(~which,nrow = 2)
dev.off()

