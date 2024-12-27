library(ComplexUpset)
library(patchwork)

### Upset plot for dynamic eGenes
GABA1 <- unique(dym_GABA$gene[dym_GABA$FDR<0.05])
npglut1 <- unique(dym_npglut$gene[dym_npglut$FDR<0.05])
nmglut1 <- unique(dym_nmglut$gene[dym_nmglut$FDR<0.05])
gene.univ <- unique(c(GABA1, npglut1, nmglut1))
dat1 <- as.data.frame(matrix(F,nrow=length(gene.univ),ncol=3),row.names = gene.univ)
colnames(dat1) <- c("GABA","npglut","nmglut")
for(i in 1:nrow(dat1)) {
  dat1$GABA[i] <- ifelse (rownames(dat1)[i] %in% GABA1, TRUE, FALSE) 
  dat1$nmglut[i] <- ifelse (rownames(dat1)[i] %in% nmglut1, TRUE, FALSE) 
  dat1$npglut[i] <- ifelse (rownames(dat1)[i] %in% npglut1, TRUE, FALSE)
}
p1 = upset(dat1, colnames(dat1),height_ratio=0.3,name="cell type sharing",sort_intersections=FALSE,sort_sets=F,
           intersections=list(
             'GABA',
             'npglut',
             'nmglut',
             c('GABA', 'npglut'),
             c('GABA', 'nmglut'),
             c('nmglut', 'npglut'),
             c('GABA', 'nmglut', 'npglut')
           ), set_sizes = F)+ggtitle("dynamic eGenes")


### Upset plot for 0hr eGenes
GABA2 <- unique(topqtl$`0hr_GABA`$gene)
npglut2 <- unique(topqtl$`0hr_npglut`$gene)
nmglut2 <- unique(topqtl$`0hr_nmglut`$gene)
gene.univ <- unique(c(GABA2, npglut2, nmglut2))
dat2 <- as.data.frame(matrix(F,nrow=length(gene.univ),ncol=3),row.names = gene.univ)
colnames(dat2) <- c("GABA","npglut","nmglut")
for(i in 1:nrow(dat2)) {
  dat2$GABA[i] <- ifelse (rownames(dat2)[i] %in% GABA2, TRUE, FALSE) 
  dat2$nmglut[i] <- ifelse (rownames(dat2)[i] %in% nmglut2, TRUE, FALSE) 
  dat2$npglut[i] <- ifelse (rownames(dat2)[i] %in% npglut2, TRUE, FALSE)
}
p2 = upset(dat2, colnames(dat2),height_ratio=0.3,name="cell type sharing",sort_intersections=FALSE,sort_sets=F,
           intersections=list(
             'GABA',
             'npglut',
             'nmglut',
             c('GABA', 'npglut'),
             c('GABA', 'nmglut'),
             c('nmglut', 'npglut'),
             c('GABA', 'nmglut', 'npglut')
           ),set_sizes = F)+ggtitle("0hr eGenes")
### Concatenate the two panels.
p1/p2+plot_layout(nrow=3)

