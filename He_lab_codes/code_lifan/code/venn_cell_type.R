### Plotting venn diagram for dynamic eGenes and 0hr eGenes
library(VennDiagram)
setwd("/project/xinhe/lifan/neuron_stim")

dym_GABA <- read.csv("/path/to/dynamic_testing_GABA")
dym_nmglut <- read.csv("/path/to/dynamic_testing_nmglut")
dym_npglut <- read.csv("/path/to/dynamic_testing_npglut")

topqtl <- readRDS("/project/xinhe/lifan/neuron_stim/mateqtl_100lines_output/final_eqtl_08222023/lead_eQTL_by_context.rds")
pdf("FigS14xx_venn_GABA_egene.pdf",height=3, width=3)
g1 <- venn.diagram(list(`0hr`=topqtl$`0hr_GABA`$gene, `1hr`=topqtl$`1hr_GABA`$gene,
                        `6hr`=topqtl$`6hr_GABA`$gene), filename=NULL,
                   col = "transparent", fill = c("blue", "green", "red"), alpha=0.5)
grid.draw(g1)
dev.off()

pdf("FigS14xx_venn_nmglut_egene.pdf",height=3, width=3)
g1 <- venn.diagram(list(`0hr`=topqtl$`0hr_nmglut`$gene, `1hr`=topqtl$`1hr_nmglut`$gene,
                        `6hr`=topqtl$`6hr_nmglut`$gene), filename=NULL,
                   col = "transparent", fill = c("blue", "green", "red"), alpha=0.5)
grid.draw(g1)
dev.off()

pdf("FigS14xx_venn_npglut_egene.pdf",height=3, width=3)
g1 <- venn.diagram(list(`0hr`=topqtl$`0hr_npglut`$gene, `1hr`=topqtl$`1hr_npglut`$gene,
                        `6hr`=topqtl$`6hr_npglut`$gene), filename=NULL,
                   col = "transparent", fill = c("blue", "green", "red"), alpha=0.5)
grid.draw(g1)
dev.off()

pdf("FigS15xx_venn_dymgene.pdf",height=3, width=3)
g1 <- venn.diagram(list(GABA=dym_GABA$gene[dym_GABA$FDR<0.05], nmglut=dym_nmglut$gene[dym_nmglut$FDR<0.05],
             npglut=dym_npglut$gene[dym_npglut$FDR<0.05]), imagetype = "png", filename=NULL,
             #filename = "venn_dymgene.png", print.mode=c("raw","percent"),
             col = "transparent", fill = c("blue", "green", "red"), alpha=0.3)
grid.draw(g1)
dev.off()
pdf("FigS15xx_venn_0hrgene.pdf",height=3, width=3)
g2 <- venn.diagram(list(GABA=topqtl$`0hr_GABA`$gene, nmglut=topqtl$`0hr_nmglut`$gene,
                  npglut=topqtl$`0hr_npglut`$gene), imagetype = "png", filename=NULL,
             #filename = "venn_0hr_eGene.png", print.mode=c("raw","percent"),
             col = "transparent", fill = c("blue", "green", "red"), alpha=0.3)
grid.newpage()
grid.draw(g2)
dev.off()


### Upset plot with ComplexUpset


