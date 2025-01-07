###Figure 1J Pseudo Time heatmap of selected genes
library(ArchR)
addArchRThreads(8)
ArchR_subset = loadArchRProject("/scratch/midway3/zichengwang/ArchR_subset/")


cell_types = c("GABA","nmglut","npglut")

markerGenes  <- c(
  "FOS", "JUNB", "NR4A1",
  "NR4A3", "BTG2","ATF3",
  "DUSP1","EGR1", "NPAS4",
  "VGF","BDNF","PCSK1", "DUSP4",
  "ATP1B1","SLC7A5", "NPTX1",
  "SCG2","CREM")




traj_heatmap = function(cell_type){
  traj = paste0(cell_type,"U")
  
  trajGSM <- getTrajectory(ArchRProj = ArchR_subset, name = traj, useMatrix = "GeneScoreMatrix")
  rownames(trajGSM) = sub(".*:","",rownames(trajGSM))
  trajGSM = trajGSM[markerGenes,]
  trajGEM <- getTrajectory(ArchRProj = ArchR_subset, name = traj, useMatrix = "GeneExpressionMatrix")
  rownames(trajGEM) = sub(".*:","",rownames(trajGEM))
  trajGEM = trajGEM[markerGenes,]
  
  ht1 <- plotTrajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "solarExtra"),  varCutOff = 0,rowOrder = markerGenes)
  ht2 <- plotTrajectoryHeatmap(trajGEM,  pal = paletteContinuous(set = "solarExtra"), varCutOff = 0,rowOrder = markerGenes)
  
  
  plotPDF(ht1 + ht2, 
          name = paste0("traj_heatmaps_final/",cell_type,"-Heatmap-corGSM_GEM-solarExtra.pdf"), 
          ArchRProj = ArchR_subset, 
          addDOC = FALSE, width = 8, height = 8)
  
  ht3 <- plotTrajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0,rowOrder = markerGenes)
  ht4 <- plotTrajectoryHeatmap(trajGEM,  pal = paletteContinuous(set = "solarExtra"), varCutOff = 0,rowOrder = markerGenes)
  
  
  plotPDF(ht3 + ht4, 
          name = paste0("traj_heatmaps_final/",cell_type,"-Heatmap-corGSM_GEM-blueYellow-solarExtra.pdf"), 
          ArchRProj = ArchR_subset, 
          addDOC = FALSE, width = 10, height = 8)
}


sapply(cell_types,traj_heatmap)


### Figure 3a gene module traj and supplementary
library(tidyverse)
library(gridExtra)

setwd("/project/xinhe/zicheng/neuron_stimulation/script/")

k=15
km.res = readRDS("/project/xinhe/zicheng/neuron_stimulation/output/gene_modules/kmeans_k15_logFC1_unionDEGs.rds")

combined_trajGEX = readRDS("/project/xinhe/zicheng/neuron_stimulation/output/gene_modules/combined_trajGEX.rds")
combined_trajGEX = t(scale(t(combined_trajGEX)))
combined_trajGEX[is.nan(combined_trajGEX)] = 0
rownames(combined_trajGEX) = gsub(".*:","",rownames(combined_trajGEX))

plot(combined_trajGEX["HDAC9",])




plot_cluster_trajectory = function(cluster_num){
  genes_in_cluster = names(km.res$cluster)[km.res$cluster == cluster_num]
  
  subset_gex = combined_trajGEX[genes_in_cluster,]
  
  mean_expr_long_df = data.frame(mean_expr = colMeans(subset_gex),
                                 lower_ci = apply(subset_gex,2,function(x) quantile(x,0.25)),
                                 upper_ci = apply(subset_gex,2,function(x) quantile(x,0.75)),
                                 time = rep(seq(0.5,99.5),3),
                                 cell_type = c(rep('GABA',100),rep('nmglut',100),rep('npglut',100)))
  ggplot(mean_expr_long_df, aes(x = time, y = mean_expr)) + geom_line() +
    facet_grid(~cell_type) + ggtitle(paste0("Cluster ",cluster_num,", nGenes = ",length(genes_in_cluster))) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.3) +
    xlab("Pseudo Time") + ylab("Expression") + theme_bw()
}

cluster_plot_list = lapply(1:k,plot_cluster_trajectory)


pdf("../output/final_figures/supp/gene_module_trajs_all.pdf", width = 10, height = 10)
grid.arrange(grobs = cluster_plot_list,nrow= 5)
dev.off()


pdf("../output/final_figures/fig_3a/gene_module_trajs_early.pdf", width = 6, height = 2)
grid.arrange(grobs = cluster_plot_list[c(6,13)],nrow= 1)
dev.off()

pdf("../output/final_figures/fig_3a/gene_module_trajs_late.pdf", width = 12, height = 2)
grid.arrange(grobs = cluster_plot_list[c(4,15,3,8)],nrow= 1)
dev.off()


pdf("../output/final_figures/fig_3a/gene_module_trajs_repression.pdf", width = 3, height = 2)
grid.arrange(grobs = cluster_plot_list[c(9)],nrow= 1)
dev.off()


############Fig 3B enrichment: See Siwei's code

######## Fig 3C: See Xiaotong's code




###Fig 3D, S13B: plot peak gene correlation
library(ArchR)
library(tidyverse)
addArchRThreads(8)

p2g_res_list = readRDS("../output/scmixedGLM/log_Mito0.01_p2g_res_list.rds")

sig_peak_gene_pairs = lapply(p2g_res_list,function(x){
  x %>% filter(fdr < 0.1,beta > 0) %>% select(Gene,Peak)
}) 


ArchR_subset = loadArchRProject("/scratch/midway3/zichengwang/ArchR_subset/")



lapply(c("GABA","nmglut","npglut"), function(cell_type){
  traj = paste0(cell_type,"U")
  

  trajGEM <- getTrajectory(ArchRProj = ArchR_subset, name = traj, useMatrix = "GeneExpressionMatrix")
  trajPM <- getTrajectory(ArchRProj = ArchR_subset, name = traj, useMatrix = "PeakMatrix")
  
  
  gene_conv = 1:length(rownames(trajGEM))
  names(gene_conv) = sub("chr.*:","",rownames(trajGEM))
  
  peak_conv = 1:length(rownames(trajPM))
  names(peak_conv) = sub("_","-",rownames(trajPM),fixed = T)
  rownames(trajPM) = sub("_","-",rownames(trajPM),fixed = T)
  
  conv_p2g = sig_peak_gene_pairs[[cell_type]] %>% mutate(idx2 = gene_conv[Gene],idx1 = peak_conv[Peak]) %>%
    select(idx1,idx2)
  
  
  trajPM2 <- trajPM[conv_p2g$idx1, ]
  trajGEM2 <- trajGEM[conv_p2g$idx2, ]
  
  trajCombined <- trajPM2
  assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajPM2), 1, scale)) + t(apply(assay(trajGEM2), 1, scale))
  
  combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
  rowOrder <- match(rownames(combinedMat), rownames(trajPM2))
  
  ht1 <- plotTrajectoryHeatmap(trajPM2,  pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0, rowOrder = rowOrder)
  
  ht2 <- plotTrajectoryHeatmap(trajGEM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)
  
  plotPDF(ht1 + ht2, 
          name = paste0("trajectory_heatmaps/",cell_type,"-Heatmap-P2G-scGLM.pdf"), 
          ArchRProj = ArchR_subset, 
          addDOC = FALSE, width = 10, height = 8)
  return("Done!")
})


################Fig 3E: abc and p2g overlaps
library(tidyverse)
library(vroom)
library(GenomicRanges)
library(parallel)


p2g_list = readRDS("../output/scmixedGLM/log_Mito0.01_p2g_res_list.rds") 
lapply(p2g_list,dim)

contexts = list.files("/scratch/midway3/zichengwang/Github_Repo/ABC-Enhancer-Gene-Prediction/results")[1:9]

# Aggregated ABC predictions
slim_agg = lapply(contexts, function(x) vroom(paste0("/scratch/midway3/zichengwang/Github_Repo/ABC-Enhancer-Gene-Prediction/results/",
                                                     x,"/Predictions/EnhancerPredictions_threshold0.021_self_promoter.tsv"))) %>% bind_rows()


p2g_abc_prop_barplot = lapply(c("GABA","nmglut","npglut"),function(x){
  #glm_p2g = p2g_list[[x]] %>% sample_n(5000) %>% separate(Peak, into = c("chr","start","end")) %>% GRanges()
  abc_p2g = slim_agg %>% filter(grepl(x,CellType,fixed = T)) %>% GRanges()
  
  non_sig_overlap_prop = lapply(c(0.1), function(y){
    glm_p2g = p2g_list[[x]] %>% filter(fdr > y,beta>0) %>% separate(Peak, into = c("chr","start","end")) %>% GRanges()
    overlap_peaks = findOverlaps(glm_p2g,abc_p2g)
    
    overlap_pairs = overlap_peaks[glm_p2g[overlap_peaks@from]$Gene == abc_p2g[overlap_peaks@to]$TargetGene]
    
    
    
    return(length(glm_p2g[unique(overlap_pairs@from)])/length(glm_p2g))
  }) %>% unlist()
  
  
  overlap_prop = mclapply(c(0.1,0.05,0.01,0.001,0.0001), function(y){
    glm_p2g = p2g_list[[x]] %>% filter(fdr <= y,beta>0) %>% separate(Peak, into = c("chr","start","end")) %>% GRanges()
    overlap_peaks = findOverlaps(glm_p2g,abc_p2g)
    
    overlap_pairs = overlap_peaks[glm_p2g[overlap_peaks@from]$Gene == abc_p2g[overlap_peaks@to]$TargetGene]    
    
    return(length(glm_p2g[unique(overlap_pairs@from)])/length(glm_p2g))
  },mc.cores = 20) %>% unlist()
  
  fdr_cutoffs = c("FDR > 0.1", "FDR \u2264 0.1", "FDR \u2264 0.05", "FDR \u2264 0.01", "FDR \u2264 0.001", "FDR \u2264 0.0001")
  return(data.frame(fdr_cutoff = factor(fdr_cutoffs,levels = fdr_cutoffs),prop = c(non_sig_overlap_prop,overlap_prop),cell_type=x))
}) %>% bind_rows()

cairo_pdf("../output/final_figures/fig_3e.pdf", width = 8, height = 6)

ggplot(p2g_abc_prop_barplot,aes(x = fdr_cutoff, y = prop*100, fill = cell_type)) +
  geom_bar(stat = "identity",position = "dodge") + 
  labs(x = "FDR", y = "% Overlapped", fill = "Cell Type") +
  theme_bw()
dev.off()

##################fisher's test of p2g-abc overlaps
library(tidyverse)
library(vroom)
library(GenomicRanges)

p2g_list = readRDS("../output/scmixedGLM/log_Mito0.01_p2g_res_list.rds")
lapply(p2g_list,dim)

contexts = list.files("/scratch/midway3/zichengwang/Github_Repo/ABC-Enhancer-Gene-Prediction/results")[1:9]

count_overlaps = function(glm_p2g,abc_p2g){
  overlap_peaks = findOverlaps(glm_p2g,abc_p2g)
  
  overlap_pairs = overlap_peaks[glm_p2g[overlap_peaks@from]$Gene == abc_p2g[overlap_peaks@to]$TargetGene]
  
  # cor(glm_p2g[overlap_pairs@from]$beta,abc_p2g[overlap_pairs@to]$ABC.Score,method = "spearman")
  
  #return(list(glm_p2g[overlap_pairs@from],abc_p2g[overlap_pairs@to]))
  return(length(unique(overlap_pairs@from)))
}



p2g_abc_fisher = function(cell_type){
  cellType_contexts = contexts[grep(cell_type,contexts,fixed = T)]
  
  abc_res = lapply(cellType_contexts, function(x){
    vroom(paste0("/scratch/midway3/zichengwang/Github_Repo/ABC-Enhancer-Gene-Prediction/results/",x,"/Predictions/EnhancerPredictionsAllPutative.tsv.gz"),
          col_select = c("chr", "start", "end","TargetGene","TargetGeneTSS","ABC.Score", "CellType"))
    
  }) %>% bind_rows()
  
  sig_abc = abc_res %>% filter(ABC.Score >= 0.021) %>% GRanges()
  #nonsig_abc = abc_res %>% filter(ABC.Score < 0.021) %>% GRanges()
  
  sig_p2g = p2g_list[[cell_type]] %>% filter(fdr <= 0.05 & beta > 0) %>% separate(Peak, into = c("chr","start","end")) %>% GRanges()
  nonsig_p2g = p2g_list[[cell_type]] %>% filter(!(fdr <= 0.05 & beta > 0)) %>% separate(Peak, into = c("chr","start","end")) %>% GRanges()
  
  
  overlap_count <- count_overlaps(sig_p2g, sig_abc)
  non_overlap_p2g <- length(sig_p2g) - overlap_count
  non_overlap_abc <- count_overlaps(nonsig_p2g,sig_abc)
  neither_count <- length(nonsig_p2g) - non_overlap_abc
  
  # Construct contingency table
  contingency_table <- matrix(
    c(overlap_count, non_overlap_p2g, non_overlap_abc, neither_count),
    nrow = 2,
    dimnames = list(
      "sig_p2g" = c("In sig_abc", "Not in sig_abc"),
      "sig_abc" = c("In sig_p2g", "Not in sig_p2g")
    )
  )
  
  # Perform Fisher's Exact Test
  fisher_result <- fisher.test(contingency_table)
  return(list(table=contingency_table,fisher_result=fisher_result))
}

fisher_res_list = lapply(c("GABA","nmglut","npglut"),p2g_abc_fisher)
names(fisher_res_list) = c("GABA","nmglut","npglut")

fisher_res_list
saveRDS(fisher_res_list,"/project/xinhe/zicheng/neuron_stimulation/output/scmixedGLM/p2g_abc_fisher_res.rds")








#### Figure 3F chromatin with expression traj
### plot gene modules
##### plot traj
library(tidyverse)
library(gridExtra)

k=15
km.res = readRDS("/project/xinhe/zicheng/neuron_stimulation/output/gene_modules/kmeans_k15_logFC1_unionDEGs.rds")


### normalize each row
combined_trajGEX = readRDS("/project/xinhe/zicheng/neuron_stimulation/output/gene_modules/combined_trajGEX.rds")
combined_trajGEX = t(scale(t(combined_trajGEX)))
combined_trajGEX[is.nan(combined_trajGEX)] = 0
rownames(combined_trajGEX) = gsub(".*:","",rownames(combined_trajGEX))

combined_trajPM = readRDS("/project/xinhe/zicheng/neuron_stimulation/output/gene_modules/combined_trajPM.rds")
combined_trajPM = t(scale(t(combined_trajPM)))

p2g_list = readRDS("../output/scmixedGLM/log_Mito0.01_p2g_res_list.rds") %>% lapply(.,function(x) filter(x,fdr<=0.05,beta>0))



plot_cluster_trajectory = function(cluster_num,legend = F){
  genes_in_cluster = names(km.res$cluster)[km.res$cluster == cluster_num]
  
  subset_gex = combined_trajGEX[genes_in_cluster,]
  
  peaks_in_cluster = lapply(p2g_list,function(x) x[x$Gene %in% genes_in_cluster,"Peak"])
  mean_chromatin = lapply(1:3, function(x){
    if (length(peaks_in_cluster[[x]]) > 0) {
      return(colMeans(combined_trajPM[peaks_in_cluster[[x]],((x-1)*100+1):(x*100),drop = FALSE]))
    }else{
      return(rep(NA,100))
    }
  }) %>% unlist()
  
  mean_activity_long_df = data.frame(mean_activity = c(colMeans(subset_gex),mean_chromatin),
                                     activity_type = c(rep('Expression',300),rep('Chromatin Accessibility',300)),
                                     time = rep(seq(0.5,99.5),6),
                                     cell_type = rep(c(rep('GABA',100),rep('nmglut',100),rep('npglut',100)),2)
  )
  if (legend == F) {
    ggplot(mean_activity_long_df, aes(x = time, y = mean_activity, color = activity_type)) + geom_line() +
      facet_grid(~cell_type) + ggtitle(paste0("Cluster ",cluster_num,", nGenes = ",length(genes_in_cluster))) +
      labs(x = "Pseudo Time", y = "Activity", color = "Activity Type") + theme_bw() +  theme(legend.position = "none")
  }else{
    ggplot(mean_activity_long_df, aes(x = time, y = mean_activity, color = activity_type)) + geom_line() +
      facet_grid(~cell_type) + ggtitle(paste0("Cluster ",cluster_num,", nGenes = ",length(genes_in_cluster))) +
      labs(x = "Pseudo Time", y = "Activity", color = "Activity Type") + theme_bw()
  }
  
}
cluster_plot_list = lapply(1:k,plot_cluster_trajectory)


setwd("/project/xinhe/zicheng/neuron_stimulation/script")

### Fig S13A
pdf("../output/final_figures/supp/gene_module_expression_w_chromatin_trajs_all.pdf", width = 10, height = 10)
grid.arrange(grobs = cluster_plot_list,nrow= 5)
dev.off()



pdf("../output/final_figures/fig_3f/gene_module_trajs_early_w_chromatin.pdf", width = 6, height = 2)
grid.arrange(grobs = cluster_plot_list[c(6,13)],nrow= 1)
dev.off()

pdf("../output/final_figures/fig_3f/gene_module_trajs_late_w_chromatin.pdf", width = 12, height = 2)
grid.arrange(grobs = cluster_plot_list[c(3,4,8,15)],nrow= 1)
dev.off()


pdf("../output/final_figures/fig_3f/gene_module_trajs_repression_w_chromatin.pdf", width = 3, height = 2)
grid.arrange(grobs = cluster_plot_list[c(9)],nrow= 1)
dev.off()

plot_cluster_trajectory(1,T)
ggsave("../output/final_figures/fig_3f/figure_legend.pdf",width = 4,height = 2)








### Figure 4A and 4C Upset plot
library(ComplexHeatmap)

load("../output/Diff_Motif/response_list_w_expr_filter.RData")
m1 = make_comb_mat(early_response_tfs)

pdf("../output/final_figures/fig_4a.pdf",width = 6,height = 3)
UpSet(m1,set_order = 1:3,top_annotation = 
        HeatmapAnnotation(
          "TFs for early response" = anno_barplot(
            comb_size(m1), 
            border = FALSE, 
            gp = gpar(fill = "black"), 
            height = unit(2, "cm"),
            axis_param = list(side = "left")
          ), 
          annotation_name_side = "left", 
          annotation_name_rot = 0))
dev.off()

m2 = make_comb_mat(late_response_tfs)

pdf("../output/final_figures/fig_4c.pdf",width = 6,height = 3)
UpSet(m2,set_order = 1:3,top_annotation = 
        HeatmapAnnotation(
          "TFs for late response" = anno_barplot(
            comb_size(m2), 
            border = FALSE, 
            gp = gpar(fill = "black"), 
            height = unit(2, "cm"),
            axis_param = list(side = "left")
          ), 
          annotation_name_side = "left", 
          annotation_name_rot = 0))
dev.off()






### Fig 4B TF expression and motif activity bubble plot
library(tidyverse)
library(ComplexHeatmap)
library(circlize)



library(tidyverse)
avg_expr_contexts = readRDS("/project/xinhe/zicheng/neuron_stimulation/data/diff_motif/Log2CPM_Avg_Expr_9contexts.rds")
avg_motif_contexts = readRDS("/project/xinhe/zicheng/neuron_stimulation/data/diff_motif/Raw_Avg_Motif_9contexts.rds")




tf_regulators = c("FOS", "EGR1", "NPAS4", "BACH2", ### shared early
                  "MEF2C", "ESR1", # shared late
                  "FOSL1","JUND","SMARCC1", # glut late
                  "NR2F2", "ID3", "DLX5", "TCF4", "NR4A2")  # gaba late"



rna_avg_mtx = avg_expr_contexts[tf_regulators,]
mm_avg_mtx = avg_motif_contexts[tf_regulators,]

#####normalize each row
rna_avg_mtx = t(scale(t(rna_avg_mtx)))
mm_avg_mtx = t(scale(t(mm_avg_mtx)))


colnames(rna_avg_mtx) = sub("_"," ",colnames(rna_avg_mtx),fixed = T)
colnames(mm_avg_mtx) = sub("_"," ",colnames(mm_avg_mtx),fixed = T)






rna_long_df = reshape2::melt(rna_avg_mtx,value.name = "Expression")
motif_long_df = reshape2::melt(mm_avg_mtx,value.name = "Motif Activity")
joined_df = full_join(rna_long_df, motif_long_df, by = c("Var1", "Var2"))


ggplot(joined_df, aes(Var2,
                      factor(Var1, levels = rev(rownames(rna_avg_mtx))),
                      fill = Expression, size = `Motif Activity`)) +
  geom_point(shape = 21, stroke = .5) +
  scale_radius(range = c(1, 8)) +
  scale_fill_gradient(low = "white",high = "red",labels = scales::number_format(accuracy = 0.1)) +
  theme_bw() +
  theme(legend.position = "bottom", 
        #panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 60,hjust = 1)) +
  guides(size = guide_legend(override.aes = list(fill = NA, color = "black", stroke = .25), 
                             label.position = "bottom",
                             title.position = "left", 
                             order = 1),
         fill = guide_colorbar(ticks.colour = NA, title.position = "left", label.position = "bottom", order = 2)) +
  labs(size = "Motif Activity", fill = "Expression", x = "Context", y = "TF")


ggsave("../output/final_figures/fig_4b.pdf",width = 4,height = 6)




### Fig 4D: plot TF trajecotries


library(tidyverse)

setwd("/project/xinhe/zicheng/neuron_stimulation/script/")

tf_regulators = c("FOS", "MEF2C", "TCF4", "DLX5") 


combined_trajGEX = readRDS("/project/xinhe/zicheng/neuron_stimulation/output/gene_modules/combined_trajGEX.rds")
combined_trajGEX = combined_trajGEX[rowSums(combined_trajGEX) > 0,]
combined_trajGEX = t(scale(t(combined_trajGEX)))
rownames(combined_trajGEX) = gsub(".*:","",rownames(combined_trajGEX))


combined_trajMM = readRDS("/project/xinhe/zicheng/neuron_stimulation/output/gene_modules/combined_trajMM.rds")
combined_trajMM = t(scale(t(combined_trajMM)))




plot_tf_traj = function(tf){
  expr_traj = combined_trajGEX[tf,]
  motif_traj = combined_trajMM[tf,]
  plot_df = data.frame(activity = c(expr_traj,motif_traj),
                       Type = c(rep("Expression",300),rep("Motif Activity",300)),
                       pseudo_time = rep(rep(0.5:99.5,3),2),
                       cell_type = rep(c(rep("GABA",100),rep("nmglut",100),rep("npglut",100)),2))
  traj_plot = plot_df %>% ggplot(aes(x = pseudo_time, y = activity, color = Type)) + geom_smooth(method = "loess", span = 0.5) +
    facet_grid(~cell_type) + theme_bw() +
    labs(x = "Pseudo Time", y = "Activity", color = "Activity Type")
  ggsave(paste0("../output/final_figures/fig_4d/",tf,"_TF_Traj.pdf"),plot = traj_plot, height = 2,width = 4.5, units = "in")
}


lapply(tf_regulators,plot_tf_traj)



#### ####### ASD GRN
library(tidyverse)

library(openxlsx)

asd_genes = scan("/project/xinhe/zicheng/neuron_stimulation/data/ASD_Fu_et_al_0.05.txt",character())


sig_grn = readRDS("/project/xinhe/zicheng/neuron_stimulation/output/scmixedGLM_GRN/log_Mito0.01_cellType_eGRN_corr_res_tf0.2_corr0.5.rds")

erg_tfs = c("TCF4","MEF2C","MLXIP","RORB")


union_sig_grn = lapply(names(sig_grn), function(x) bind_cols(sig_grn[[x]],cell_type = x)) %>%
  bind_rows() %>% group_by(tf, gene) %>%
  mutate(ObservedIn = case_when(
    all(cell_type == "GABA") ~ "GABA only",
    all(cell_type %in% c("nmglut", "npglut")) ~ "Glut only",
    any(cell_type == "GABA") & any(cell_type %in% c("nmglut", "npglut")) ~ "Both"
  )) %>%
  filter(abs(correlation) == max(abs(correlation))) 


union_sig_grn = union_sig_grn %>% mutate(Color = case_when(
  ObservedIn == "GABA only" ~ "Blue",
  ObservedIn == "Glut only" ~ "Red",
  ObservedIn == "Both" ~ "Purple"
))

union_sig_grn %>% filter(
  tf %in% c("TCF4"),
  gene %in% asd_genes
) 

library(igraph)
library(extrafont)

font_import(prompt = FALSE)
# Check available fonts
fonts()
# Load fonts for the current session
loadfonts(device = "pdf")


### plot network with union GRN
plot_union_network = function(df_cor_2){
  g <- graph_from_data_frame(df_cor_2[,1:3], directed=T)
  
  tf_colors =  rep("orange", length(unique(df_cor_2$tf)))
  names(tf_colors) <- unique(df_cor_2$tf)
  
  gene_colors <- rep("green", length(unique(df_cor_2$gene)))
  names(gene_colors) <- unique(df_cor_2$gene)
  
  vertex.color <- c(tf_colors, gene_colors)
  convert_colors_to_sizes <- function(color) {
    if(color == "orange") {
      return(1)
    } else if(color == "green") {
      return(0.8)
    } else {
      return(NA)
    }
  }
  
  # Apply the conversion function to each element in the vector
  vertex.size <- sapply(vertex.color, convert_colors_to_sizes)
  
  vertex.frame.color <- c(tf_colors, gene_colors)
  
  
  ## compute layout
  set.seed(42)
  co <- layout_with_fr(g, dim = 2, niter = 1000)
  
  ## only show labels for TF
  #V(g)$label <- ifelse(V(g)$name %in% df_cor_2$tf, V(g)$name, NA)
  
  ## for TFs, we increase the size of nodes
  V(g)$size <- ifelse(V(g)$name %in% df_cor_2$tf, 15, 10)
  
  #E(g)$weight <- abs( E(g)$correlation)

  
  
  # pdf("../output/eGRN/eGRN_complete.pdf",width = 25, height = 25)
  plot(g, layout = co, 
       vertex.color = vertex.color,
       vertex.frame.color = vertex.frame.color,
       vertex.label.dist = 0.1,
       edge.width = E(g)$weight*2,
       vertex.label.cex = vertex.size,
       edge.color = df_cor_2$Color,
       vertex.label.family = "Helvetica"# Set font to Helvetica
       
  )
  legend("topleft",
         c("TF", "Gene"),
         col = c("orange", "green"),
         pch = 21,
         pt.bg = c("orange", "green"),
         pt.cex = c(4.5, 3),
         cex = 1.6,
         bty = "n")
  
  legend("topright",
         c("GABA only", "Glut only", "Both"),
         col = c("Blue", "Red", "Purple"),
         lty = 1,
         lwd = 2,
         cex = 1.6,
         bty = "n")
  
  output_plot = recordPlot()
}

pdf("../output/final_figures/fig_4f.pdf",width = 10, height = 8)
union_sig_grn %>% filter(
  tf %in% erg_tfs,
  gene %in% asd_genes
) %>% plot_union_network()
dev.off()




########Fig 4F odds ratio heatmap
library(tidyverse)
library(openxlsx)
library(parallel)

sheet_names = getSheetNames("/project/xinhe/zicheng/neuron_stimulation/output/scmixedGLM_GRN/ASD_Enriched_TFs_by_CellType.xlsx")
sig_asd_enrichment_list = lapply(sheet_names,function(x) read.xlsx("/project/xinhe/zicheng/neuron_stimulation/output/scmixedGLM_GRN/ASD_Enriched_TFs_by_CellType.xlsx",sheet = x))
names(sig_asd_enrichment_list) = sheet_names
top10_tfs = lapply(sig_asd_enrichment_list, function(x) x$TF %>% .[1:10]) %>% unlist() %>% unique()


sig_grn = readRDS("/project/xinhe/zicheng/neuron_stimulation/output/scmixedGLM_GRN/log_Mito0.01_cellType_eGRN_corr_res_tf0.2_corr0.5.rds")
asd_genes = scan("/project/xinhe/zicheng/neuron_stimulation/data/ASD_Fu_et_al_0.05.txt",character())
union_degs = readRDS("/project/xinhe/zicheng/neuron_stimulation/output/gene_modules/DEG_Lexi_union.rds")


asd_enriched_tfs = lapply(c("GABA","nmglut","npglut"), function(cell_type){
  union_targets = sig_grn[[cell_type]]
  union_targets = split(union_targets$gene,union_targets$tf)
  union_targets = lapply(union_targets, function(x) unique(x))
  union_tfs = lapply(sig_grn,function(x) x$tf) %>% unlist() %>% unique()
  
  
  asd_match = union_degs %in% asd_genes
  
  
  asd_enrichment = lapply(top10_tfs, function(x){
    tf_targets_match = union_degs %in% union_targets[[x]]
    
    res_fisher = fisher.test(table(asd_match,tf_targets_match))
    
    res_vector = c(res_fisher$estimate,res_fisher$p.value)
    names(res_vector) = c("odd_ratio","p_value")
    return(res_vector)
  }) %>% bind_rows() %>% cbind(TF = top10_tfs,.) %>% mutate(FDR = p.adjust(p_value, method = "fdr"))
  
  return(asd_enrichment)
})
names(asd_enriched_tfs) = c("GABA","nmglut","npglut")

odds_ratio_mtx = lapply(asd_enriched_tfs, function(x){
  odds_ratio = x$odd_ratio
  names(odds_ratio) = x$TF
  return(odds_ratio)
})  %>% do.call(cbind,.)

library(ComplexHeatmap)
library(circlize)  # For custom color functions


col_fun <- colorRamp2(c(1, 6), c("white", "darkred"))

pdf("../output/final_figures/fig_4g.pdf",width = 6, height = 8)
# Generate the heatmap
Heatmap(
  odds_ratio_mtx,
  name = "Odds Ratio", 
  row_title = "TF", 
  column_title = "Cell Type",
  col = col_fun,
  border = TRUE,  # Add borders to separate cells
  rect_gp = gpar(col = "black", lwd = 0.1)  # Black separator lines with custom width
)
dev.off()



# Fig 7D and S15D: Dynamic caQTL boxplots by genotype
library(tidyverse)
library(vroom)
library(data.table)
library(parallel)

load("/project/xinhe/zicheng/neuron_stimulation/caQTL/output/caQTL_mapping/res_analysis/cPeak_list_25kb.RData")
snp_pos = vroom("/project2/xinhe/lifanl/neuron_stim/geno_info_100lines.txt",col_names = c("chr","pos","a2","a1","rsid")) %>%
  as.data.frame()
rownames(snp_pos) = snp_pos$rsid
snp = data.frame(fread("/project/xinhe/lifan/neuron_stim/mateqtl_input_100lines/snps.txt"),row.names = 1)
logTPM_plot = readRDS("/project/xinhe/zicheng/neuron_stimulation/caQTL/data/count_matrix/TMM_plot/logTPM_List_for_Plot.rds")

boxplot_by_genotype = function(peak_id, snp_id,cluster_to_plot = clusters) {
  
  boxplot_df = mclapply(cluster_to_plot, mc.cores = 3,
                        function(cluster) {
                          subset_snp_pos = unlist(snp_pos[snp_id,])
                          genotype_dict = c(paste0(subset_snp_pos["a2"],"/",subset_snp_pos["a2"]),
                                            paste0(subset_snp_pos["a2"],"/",subset_snp_pos["a1"]),
                                            paste0(subset_snp_pos["a1"],"/",subset_snp_pos["a1"]))
                          names(genotype_dict) = c("0","1","2")
                          subset_tmm_counts = as.matrix(logTPM_plot[[cluster]])[peak_id,]
                          names(subset_tmm_counts) = gsub('\\..*',"",names(subset_tmm_counts))
                          subset_snp = genotype_dict[as.character(snp[snp_id,names(subset_tmm_counts)])]
                          subset_snp = factor(subset_snp,levels = genotype_dict)
                          
                          subset_boxplot_df = data.frame(sample_id = names(subset_tmm_counts),
                                                         Genotype=subset_snp,
                                                         TMM_Counts = subset_tmm_counts,
                                                         from_cluster = gsub("__"," ",cluster))
                          rownames(subset_boxplot_df) = NULL
                          
                          return(subset_boxplot_df)
                        }
  ) %>% bind_rows()
  p = ggplot(boxplot_df,aes(x = Genotype, y = TMM_Counts,fill = Genotype)) +
    geom_boxplot() + geom_jitter() + facet_wrap(~from_cluster) +
    scale_fill_brewer(palette="Paired") + theme_minimal() +
    labs(title = paste0(sub("-",":",peak_id), " ~ ", snp_id)) + ylab("log(CPM) Counts")
  return(p)
}


###GABA 
boxplot_by_genotype("chr8-16810827-16811327","rs35542305",clusters)
ggsave("/project/xinhe/zicheng/neuron_stimulation/output/final_figures/fig_7d.pdf",width = 7,height = 6)

###nmglut
boxplot_by_genotype("chr1-38758687-38759187","rs2483119",clusters)
ggsave("/project/xinhe/zicheng/neuron_stimulation/output/final_figures/supp/fig_s16d/nmglut_dynamic_chr1:38758687-38759187_rs2483119.pdf",width = 7,height = 6)

###npglut
boxplot_by_genotype("chr7-139651308-139651808","rs940543",clusters)
ggsave("/project/xinhe/zicheng/neuron_stimulation/output/final_figures/supp/fig_s16d/npglut_dynamic_chr7:139651308-139651808_rs940543.pdf",width = 7,height = 6)








###Supplementary tables


#########Table S11 table of gene module
library(tidyverse)

km.res = readRDS("/project/xinhe/zicheng/neuron_stimulation/output/gene_modules/kmeans_k15_logFC1_unionDEGs.rds")

data.frame(Gene = names(km.res$cluster), Cluster = paste0("C",km.res$cluster),cluster_num = km.res$cluster) %>%
  arrange(cluster_num, Gene) %>%
  select(Gene, Cluster) %>%
  openxlsx::write.xlsx(.,"/project/xinhe/zicheng/neuron_stimulation/output/supp_tables/table_s13x_genes_in_clusters.xlsx")



######Table S12 GO enrichment table
library(tidyverse)

k = 15
enriched_list = readRDS("/project/xinhe/zicheng/neuron_stimulation/output/gene_modules/kmeans_k15_logFC1_enriched_list.rds")

enrich_datasets = names(enriched_list[[1]])[1:3]
enrich_res = list()
for (dat in enrich_datasets) {
  dat_list = list()
  for (i in 1:k) {
    dat_list[[paste0("Cluster ",i)]] = as.data.frame(enriched_list[[i]][[dat]])  %>% mutate(Cluster = i)
  }
  enrich_res[[dat]] = bind_rows(dat_list) %>% mutate(fdr = p.adjust(P.value, method = "fdr")) %>% filter(Adjusted.P.value <= 0.1) %>%
    arrange(desc(Combined.Score))
  
}


openxlsx::write.xlsx(enrich_res,"/project/xinhe/zicheng/neuron_stimulation/output/gene_modules/kmeans_k15_logFC1_GO_summary.xlsx")

###########Table S13 Peak2Gene FDR 0.1
library(tidyverse)

p2g_list = readRDS("/project/xinhe/zicheng/share/log_Mito0.01_p2g_res_list.rds")

lapply(p2g_list, function(x){
  x %>% filter(fdr <= 0.1,beta>0) %>% arrange(fdr) %>%
    rename(Beta = beta, SE = se, Pvalue = P, FDR = fdr, Qvalue = qval) %>%
    select(Gene, Peak, Beta, SE, Z, Pvalue, FDR)
}) %>% openxlsx::write.xlsx(.,"../output/supp_tables/table_sxx_Peak2Gene_FDR0.1.xlsx")



###Table S16 ABC score tables
library(tidyverse)
library(vroom)

contexts = list.files("/scratch/midway3/zichengwang/Github_Repo/ABC-Enhancer-Gene-Prediction/results")[1:9]

cellType_abc = lapply(c("GABA","nmglut","npglut"), function(cell_type){
  cellType_contexts = contexts[grep(cell_type,contexts, fixed = T)]
  slim_agg = lapply(cellType_contexts, function(x) vroom(paste0("/scratch/midway3/zichengwang/Github_Repo/ABC-Enhancer-Gene-Prediction/results/",
                                                       x,"/Predictions/EnhancerPredictions_threshold0.021_self_promoter.tsv"))) %>%
    bind_rows() %>%
    rename(Context = CellType)
})
names(cellType_abc) = c("GABA","nmglut","npglut")


openxlsx::write.xlsx(cellType_abc,"../output/supp_tables/table_sxx_ABC_EnhancerPredictions_threshold0.021_self_promoter.xlsx")



#########Table S17 motif activity differential test
library(tidyverse)
load("../output/Diff_Motif/motif_diff_pval_fdr.RData")
pval_list = lapply(motif_dff_pval_list,function(x){
  test_res = x$wilcox.p
  colnames(test_res) = paste0("Pval-",colnames(test_res))
  return(test_res)
})

fdr_list = lapply(motif_dff_fdr_list,function(x){
  test_res = x$wilcox.p
  colnames(test_res) = paste0("FDR-",colnames(test_res))
  return(test_res)
})

combined_motif_diff_list = lapply(names(motif_dff_fdr_list), function(x){
  bind_cols(list(TF = rownames(pval_list[[x]]),pval_list[[x]],fdr_list[[x]][rownames(pval_list[[x]]),]))
})

names(combined_motif_diff_list) = names(motif_dff_fdr_list)
openxlsx::write.xlsx(combined_motif_diff_list,"/project/xinhe/zicheng/neuron_stimulation/output/supp_tables/table_s16_motif_diff_test_res.xlsx")






########Table S18 GRN table
sig_grn = readRDS("/project/xinhe/zicheng/neuron_stimulation/output/scmixedGLM_GRN/log_Mito0.01_cellType_eGRN_corr_res_tf0.2_corr0.5.rds")
lapply(sig_grn, function(x){
  colnames(x) = c("TF","Target Gene","Correlation")
  return(x)
}) %>% openxlsx::write.xlsx(.,"/project/xinhe/zicheng/neuron_stimulation/output/supp_tables/table_sxx_all_GRN_cutoffCorr0.5.xlsx")


###table s19 num of targets for 4 ASD TFs
library(tidyverse)
library(ComplexHeatmap)
sig_grn = readRDS("../output/scmixedGLM_GRN/log_Mito0.01_cellType_eGRN_corr_res_tf0.2_corr0.5.rds")

tf_count_list = lapply(sig_grn, function(x) table(x$tf))


key_tfs <- c("MEF2C", "RORB", "MLXIP", "TCF4")

target_counts_mtx = matrix(nrow = length(key_tfs),ncol = 3,dimnames = list(key_tfs,c("GABA","nmglut","npglut")))

for (tf in key_tfs) {
  for (cell_type in c("GABA","nmglut","npglut")) {
    if (tf %in% names(tf_count_list[[cell_type]])) {
      target_counts_mtx[tf,cell_type] = tf_count_list[[cell_type]][tf]
    }else{
      target_counts_mtx[tf,cell_type] = 0
    }
  }
}

rowSums(target_counts_mtx) 

reshape2::melt(target_counts_mtx, id.vars = rownames(target_counts_mtx), variable.name = "Type", value.name = "Value") %>%
  rename(TF = Var1, `Cell Type` = Var2, `Number of Target Genes` = Value) %>%
  openxlsx::write.xlsx(.,"/project/xinhe/zicheng/neuron_stimulation/output/supp_tables/table_s17_ASD_4_TFs_num_targets.xlsx")




######## Table S20 4 ASD ERG TF GO enrichment
library(tidyverse)
library(enrichR)
library(openxlsx)
library(parallel)


sig_grn = readRDS("/project/xinhe/zicheng/neuron_stimulation/output/scmixedGLM_GRN/log_Mito0.01_cellType_eGRN_corr_res_tf0.2_corr0.5.rds")
union_targets = sig_grn %>% bind_rows(.) %>% select(tf, gene) %>% filter(!duplicated(.))
union_targets = split(union_targets$gene,union_targets$tf)

erg_tfs = c("TCF4","MEF2C","MLXIP","RORB")
union_targets = union_targets[erg_tfs]

dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023")

erg_targets_enriched_list <- lapply(union_targets, function(x){
  Sys.sleep(1)
  return(enrichr(unique(x), dbs))
})
saveRDS(erg_targets_enriched_list,"/project/xinhe/zicheng/neuron_stimulation/output/scmixedGLM_GRN/erg_targets_GO_enriched_list.rds")

erg_targets_enriched_list = readRDS("/project/xinhe/zicheng/neuron_stimulation/output/scmixedGLM_GRN/erg_targets_GO_enriched_list.rds")

dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023")

agg_enrichr_res = lapply(dbs, function(x){
  lapply(names(erg_targets_enriched_list), function(y){
    bind_cols(TF= y, erg_targets_enriched_list[[y]][[x]])
  }) %>% bind_rows(.)
})

names(agg_enrichr_res) = dbs


agg_enrichr_res_fdr = lapply(agg_enrichr_res, function(x){
  x[,c(1:4,8:10)] %>% mutate(FDR = p.adjust(P.value, method = "fdr")) %>%
    arrange(FDR) %>%
    relocate(FDR, .after = P.value)
})

openxlsx::write.xlsx(agg_enrichr_res_fdr,"/project/xinhe/zicheng/neuron_stimulation/output/supp_tables/table_s18_ASD_EarlyResponse_risk_TFs_GO_Enrichment.xlsx")


sheet_names = getSheetNames("/project/xinhe/zicheng/neuron_stimulation/output/scmixedGLM_GRN/ASD_Enriched_TFs_by_CellType.xlsx")
sig_asd_enrichment_list = lapply(sheet_names,function(x) read.xlsx("/project/xinhe/zicheng/neuron_stimulation/output/scmixedGLM_GRN/ASD_Enriched_TFs_by_CellType.xlsx",sheet = x))
names(sig_asd_enrichment_list) = sheet_names


lapply(sig_asd_enrichment_list, function(x){
  colnames(x) = c("TF","Odds Ratio","P-value","FDR")
  return(x)
}) %>% write.xlsx(.,"/project/xinhe/zicheng/neuron_stimulation/output/supp_tables/table_s20_asd_enriched_TFs.xlsx")

####Table S21 ASD 3 or more shared target GO enrichment
library(tidyverse)
library(enrichR)
library(openxlsx)

asd_genes = scan("/project/xinhe/zicheng/neuron_stimulation/data/ASD_Fu_et_al_0.05.txt",character())


sig_grn = readRDS("/project/xinhe/zicheng/neuron_stimulation/output/scmixedGLM_GRN/log_Mito0.01_cellType_eGRN_corr_res_tf0.2_corr0.5.rds")
union_targets = sig_grn %>% bind_rows(.)
union_targets = split(union_targets$gene,union_targets$tf)
union_targets = lapply(union_targets, function(x) unique(x))



erg_tfs = c("TCF4","MEF2C","MLXIP","RORB")
overlapped_targets = union_targets[erg_tfs] %>% unlist() %>% table() %>% .[. >=3 ] %>% names(.)
dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023","Reactome_2022","GWAS_Catalog_2023","DisGeNET","MAGMA_Drugs_and_Diseases")
enriched_list <-  enrichr(overlapped_targets, dbs)
write.xlsx(enriched_list,"/project/xinhe/zicheng/neuron_stimulation/output/scmixedGLM_GRN/ASD_EarlyResponse_risk_TFs_Union_Shared_Targets_3_or_More.xlsx")


####TABLE S22  ASD enriched TFs
library(tidyverse)
library(parallel)

sig_grn = readRDS("/project/xinhe/zicheng/neuron_stimulation/output/scmixedGLM_GRN/log_Mito0.01_cellType_eGRN_corr_res_tf0.2_corr0.5.rds")
asd_genes = scan("/project/xinhe/zicheng/neuron_stimulation/data/ASD_Fu_et_al_0.05.txt",character())
union_degs = readRDS("/project/xinhe/zicheng/neuron_stimulation/output/gene_modules/DEG_Lexi_union.rds")


asd_enriched_tfs = mclapply(c("GABA","nmglut","npglut"), function(cell_type){
  union_targets = sig_grn[[cell_type]]
  union_targets = split(union_targets$gene,union_targets$tf)
  union_targets = lapply(union_targets, function(x) unique(x))
  union_tfs = lapply(sig_grn,function(x) x$tf) %>% unlist() %>% unique()
  
  
  asd_match = union_degs %in% asd_genes
  
  
  asd_enrichment = lapply(names(union_targets), function(x){
    tf_targets_match = union_degs %in% union_targets[[x]]
    
    res_fisher = fisher.test(table(asd_match,tf_targets_match))
    
    res_vector = c(res_fisher$estimate,res_fisher$p.value)
    names(res_vector) = c("odd_ratio","p_value")
    return(res_vector)
  }) %>% bind_rows() %>% cbind(TF = names(union_targets),.) %>% mutate(FDR = p.adjust(p_value, method = "fdr"))
  sig_asd_enrichment = asd_enrichment %>% filter(odd_ratio > 1,FDR <= 0.05) %>% arrange(FDR)
  
  return(sig_asd_enrichment)
},mc.cores = 3)
names(asd_enriched_tfs) = c("GABA","nmglut","npglut")

lapply(asd_enriched_tfs, function(x) x$TF) %>% Reduce(union,.) %>% length(.)

library(openxlsx)
write.xlsx(asd_enriched_tfs,"/project/xinhe/zicheng/neuron_stimulation/output/scmixedGLM_GRN/ASD_Enriched_TFs_by_CellType.xlsx")


####table S25 caQTL results
library(tidyverse)
load("/project/xinhe/zicheng/neuron_stimulation/caQTL/output/caQTL_mapping/res_analysis/cPeak_list_25kb.RData")
names(cPeak_list) = sub("__","_",names(cPeak_list))

lapply(cPeak_list, function(x){
  x %>% mutate(phenotype_id = sub("-",":",phenotype_id,fixed = T))
}) %>% openxlsx::write.xlsx("/project/xinhe/zicheng/neuron_stimulation/output/supp_tables/table_s25_cPeak_results.xlsx")



##table S27 dynamic caQTL
library(tidyverse)
library(qqman)
library(qvalue)

load("/project/xinhe/zicheng/neuron_stimulation/caQTL/output/dynamic_QTL/final/LM_Perm_output.RData")


dynamic_qtl_perm = lapply(names(dyn_test_res_list), function(x){
  lapply(dyn_test_res_list[[x]], function(y) y$emp_pvals) %>% bind_rows() %>%
    bind_cols(topQTL_list[[x]],.)
  
}) 
names(dynamic_qtl_perm) = names(dyn_test_res_list)

dynamic_qtl_fdr = lapply(dynamic_qtl_perm, function(x){
  fdr_df = apply(x[,3:8],2,function(y) qvalue(y)$qvalues)
  colnames(fdr_df) = c("emp_1hr_qval","beta_1hr_qval",
                       "emp_6hr_qval","beta_6hr_qval",
                       "emp_anova_qval","beta_anova_qval")
  cbind(x[,1:2],fdr_df)
})




lapply(dynamic_qtl_perm, function(x){
  x %>% mutate(phenotype_id = sub("-",":",phenotype_id,fixed = T),empirical_fdr = p.adjust(emp_anova_pval, method = "fdr")) %>%
    select(phenotype_id,variant_id, empirical_pval = emp_anova_pval,empirical_fdr )
}) %>% openxlsx::write.xlsx("/project/xinhe/zicheng/neuron_stimulation/output/supp_tables/table_sx_dynamic_caQTLs_permutation_anova_results.xlsx")







