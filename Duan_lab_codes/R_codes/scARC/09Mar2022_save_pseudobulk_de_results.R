sum(res_GABA_0v6$padj < 0.05, na.rm = T)
sum(is.na(res_GABA_0v6$padj))
res_GABA_0v6$padj[is.na(res_GABA_0v6$padj)] <- Inf
sum(is.na(res_GABA_0v6$log2FoldChange))
res_GABA_0v6$log2FoldChange[is.na(res_GABA_0v6$log2FoldChange)] <- 0

res_GABA_0v6_up <- as.data.frame(res_GABA_0v6[res_GABA_0v6$log2FoldChange > 0 & res_GABA_0v6$padj < 0.05, ])
write.table(res_GABA_0v6_up, file = "res_GABA_0v6_up.csv", quote = F,
            sep = ",", col.names = T)
res_GABA_0v6_do <- as.data.frame(res_GABA_0v6[res_GABA_0v6$log2FoldChange < 0 & res_GABA_0v6$padj < 0.05, ])
write.table(res_GABA_0v6_do, file = "res_GABA_0v6_do.csv", quote = F,
            sep = ",", col.names = T)


sum(res_glut_0v6$padj < 0.05, na.rm = T)
sum(is.na(res_glut_0v6$padj))
res_glut_0v6$padj[is.na(res_glut_0v6$padj)] <- Inf
sum(is.na(res_glut_0v6$log2FoldChange))
res_glut_0v6$log2FoldChange[is.na(res_glut_0v6$log2FoldChange)] <- 0

res_glut_0v6_up <- as.data.frame(res_glut_0v6[res_glut_0v6$log2FoldChange > 0 & res_glut_0v6$padj < 0.05, ])
write.table(res_glut_0v6_up, file = "res_glut_0v6_up.csv", quote = F,
            sep = ",", col.names = T)
res_glut_0v6_do <- as.data.frame(res_glut_0v6[res_glut_0v6$log2FoldChange < 0 & res_glut_0v6$padj < 0.05, ])
write.table(res_glut_0v6_do, file = "res_glut_0v6_do.csv", quote = F,
            sep = ",", col.names = T)


sum(res_NPC_0v6$padj < 0.05, na.rm = T)
sum(is.na(res_NPC_0v6$padj))
res_NPC_0v6$padj[is.na(res_NPC_0v6$padj)] <- Inf
sum(is.na(res_NPC_0v6$log2FoldChange))
res_NPC_0v6$log2FoldChange[is.na(res_NPC_0v6$log2FoldChange)] <- 0

res_NPC_0v6_up <- as.data.frame(res_NPC_0v6[res_NPC_0v6$log2FoldChange > 0 & res_NPC_0v6$padj < 0.05, ])
write.table(res_NPC_0v6_up, file = "res_NPC_0v6_up.csv", quote = F,
            sep = ",", col.names = T)
res_NPC_0v6_do <- as.data.frame(res_NPC_0v6[res_NPC_0v6$log2FoldChange < 0 & res_NPC_0v6$padj < 0.05, ])
write.table(res_NPC_0v6_do, file = "res_NPC_0v6_do.csv", quote = F,
            sep = ",", col.names = T)

