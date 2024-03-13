
dir.create ('Plots')

source ('PTENL_study/useful_functions.R')
source ('PTENL_study/load_libraries.R')
source ('PTENL_study/ggplot_aestetics.R')

srt = readRDS ('srt.rds')

# Add pallette
trm_pal = setNames (c('black','grey44', 'darkgreen'), c('CB2_PTENL','CB2_PTENL_C124S','GFP')) 
pal_heatmap = c('#67A9CF','#EF8A62')

# FIGURE 6A-B - Generate UMAPs and celltype proportions plot ####
reductionName = 'sampleID_harmony_umap'
metaGroupName = 'celltype2'

umap1 = DimPlot (srt, group.by = metaGroupName, pt.size=0.01, reduction = reductionName, label=F)+ NoAxes() + NoLegend()
umap2 = DimPlot (srt, group.by = 'treatment', pt.size=0.01, reduction = reductionName)+ NoAxes()

cc_box1 = cellComp (
  seurat_obj = srt, 
  metaGroups = c('sampleID', metaGroupName,'treatment',metaGroupName),
  plot_as = 'box',
  ptable_factor = c(1),
  pal = trm_pal,
  ) + gtheme_no_text

cc_df = cc_box1$data
stat.test <- cc_df %>%
  group_by_at (metaGroupName) %>%
  t_test(Freq ~ treatment) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()
stat.test <- stat.test %>% add_xy_position (x = metaGroupName, step.increase=0.01)

pdf (paste0('Plots/cell_composition_',metaGroupName,'.pdf'), width=13, height=3)
(umap1 | cc_box1) 
plot_layout (widths= c(2,5))
dev.off()

metaGroupName = 'celltype_TIMTAM'
cc_box1 = cellComp (
  seurat_obj = srt, 
  metaGroups = c('sampleID', metaGroupName,'treatment'),
  plot_as = 'box',
  ptable_factor = c(1),
  pal = trm_pal,
  ) + gtheme
order_subtypes = do.call (rbind, lapply (split (cc_box1$data,cc_box1$data$celltype_TIMTAM), function(x) mean (x$Freq)))
order_subtypes = setNames (order_subtypes[,1], rownames(order_subtypes))
cc_box1$data$celltype_TIMTAM = factor (cc_box1$data$celltype_TIMTAM, levels = names(order_subtypes[order(-order_subtypes)]))

library (rstatix)
cc_df = cc_box1$data
stat.test <- cc_df %>%
  group_by_at (metaGroupName) %>%
  t_test(Freq ~ treatment) %>%
  adjust_pvalue(method = "none") %>%
  add_significance()
stat.test = stat.test %>% add_xy_position (x = metaGroupName, step.increase=0.1)

cc_box1 = cc_box1 + stat_pvalue_manual (stat.test, remove.bracket=FALSE,
   bracket.nudge.y = .01, hide.ns = TRUE,
    label = "p.adj.signif")
pdf (paste0('Plots/cell_composition_',metaGroupName,'.pdf'), width=7, height=3)
(umap1 | cc_box1)
plot_layout (widths= c(2,6))
dev.off()


### FIGURE 6C - fgsea analysis of DEG between (PTENL vs controls (PTENL_mut + GFP) identified with muscat ####
fgseaResAll = readRDS ('PTENL_study/data/fgsea_results.rds')
# Plot dotplot of fGSEA annotations per cluster 
top_pathways = Inf
int_pathways = c('interferon','antigen','cytokine')
int_pathways2 = unlist (lapply (fgseaResAll[['c5.bp.v7.1.symbol.gmt']], function(x) x$pathway))
int_pathways3 = unique(unlist(lapply (int_pathways,function(x) int_pathways2[grep (x, int_pathways2)])))
fgseaResAll2 = lapply (fgseaResAll[['c5.bp.v7.1.symbol.gmt']], function(x) x[x$pathway %in% int_pathways3,])
  fgseaResAll_dp = dotGSEA (fgseaResAll2, padj_threshold = 0.05, 
    type = 'fgsea',top_pathways = top_pathways,
    cluster_rows=T,
    cluster_cols=T)
    
png (paste0('Plots/fGSEA_annotation_c5.bp.v7.1.symbol.gmt_dotplots2.png'),2800,1000, res=300)
print(fgseaResAll_dp)
dev.off()



### FIGURE 6D - heatmap of average expression differences between PTENL vs controls of genes in relevant fgsea pathways ####
#pathway = c('GO_antigen_processing_and_presentation','GO_antigen_processing_and_presentation_of_peptide_antigen')
celltype = 'DCs'
pathway = c('GO_cytokine_production','GO_cellular_response_to_interferon_gamma','GO_T_cell_cytokine_production','GO_cytokine_production_involved_in_immune_response','GO_response_to_cytokine','GO_positive_regulation_of_cytokine_production','GO_cellular_response_to_cytokine_stimulus')
celltype = c('MregDCs','TAM2','TIM_NC','DCs','TIM_C')
fgsea1 = as.data.frame (do.call (rbind,fgseaResAll2[celltype]))
fgsea1[fgsea1$pathway %in% pathway,'leadingEdge']

gene2 = unique(unlist(lapply(fgseaResAll2, function(x) x[grep ('antigen', x$pathway), 'leadingEdge'])))
#gene2 = unique (unlist (fgsea1[fgsea1$pathway %in% pathway,'leadingEdge']))

ext_avg = AverageExpression (srt, features = gene2, group.by = c('sampleID','celltype_TIMTAM','treatment2'))
ext_avg = log2(as.data.frame (t(ext_avg[[1]])+1))
ext_avg = split (ext_avg, grepl ('control', rownames(ext_avg)))
celltypes = c('DCs','MregDCs','Neutrophils','TAM1','TAM2','TIM_C','TIM_NC','TNKcells','pDCs')
ext_avg = lapply (ext_avg, function(x) do.call (rbind, lapply (celltypes, function(y) colMeans(x[grepl (y, rownames(x)),]))))
names(ext_avg) = c('PTENL','control')
ext_avg = lapply (names(ext_avg), function(x) {rownames (ext_avg[[x]]) = paste0(x,'_',celltypes); ext_avg[[x]]})
ext_avg = do.call (rbind, ext_avg)
ext_avg = t(ext_avg)
agr3 = ext_avg[,grep ('PTENL',colnames(ext_avg))] - ext_avg[,grep ('control',colnames(ext_avg))] 
rev(RColorBrewer::brewer.pal(3,'RdBu'))
col_fun = colorRamp2(c(-1, 0, 1), c("#67A9CF", "#F7F7F7", "#EF8A62"))

colnames(agr3) = c('DCs','MregDCs','Neutrophils','TAM1','TAM2','TIM_C','TIM_NC','TNKcells','pDCs')
avg_hm = Heatmap (agr3, 
  cluster_rows=T, 
  cluster_columns=T,
  column_split = colnames(agr3),
  width = 1,
  col = col_fun,
  row_names_gp = gpar(fontsize = 7), 
  column_names_gp = gpar (fontsize = 6),
  border=T)
png (paste0('Plots/antigen_genes_heatmap_avg3.png'), width = 1050, height= 1900, res=300)
avg_hm
dev.off() 


### FIGURE 6E - CellphoneDB analysis identifying differentially expressed ligand-receptors between PTENL and CONTROLS (PTENL_MUT + GFP) ####
results.df = readRDS ('PTENL_study/data/cellphonedb_results.rds')

# Filtering low variable interactions
results.df = results.df[results.df$diffprop != 0, ]
results.df = results.df[results.df$log10pval >= .95, ] # Retain interaction with low pvalues (-log10pval)

# Generate dotplot of celltypes x ligand-receptors
results.df$celltype_pair = paste0(results.df$sampleID,'_',results.df$celltype_pair) 
  mat1 = results.df %>% 
    dplyr::select(LR, celltype_pair, diffprop) %>%  
    pivot_wider(names_from = LR, values_from = diffprop) %>% 
    data.frame() %>% replace (.=="NULL", NA)
  row.names(mat1) = mat1$celltype_pair  # put gene in `row`
  mat1 = mat1[,-1] #drop gene column as now in rows
  mat1[is.na (mat1)] = 0
  
  mat2 = results.df %>% 
    dplyr::select(LR, celltype_pair, diffprop) %>%  
    pivot_wider(names_from = celltype_pair, values_from = diffprop) %>% 
    data.frame() %>% replace(.=="NULL", NA) %>% data.frame ()
  row.names(mat2) = mat2$LR  # put gene in `row`
  mat2 = mat2[,-1] #drop gene column as now in rows
  mat2[is.na (mat2)] = 0
  
  clust1 = 1-cor (t(mat1 %>% as.matrix ()))
  clust1[is.na(clust1)] = 0
  clust1 = hclust(as.dist(clust1))# hclust with distance matrix
  
  clust2 = 1-cor (t(mat2 %>% as.matrix ()))
  clust2[is.na(clust2)] = 0
  clust2 = hclust(as.dist(clust2))# hclust with distance matrix
  
  ddgram1 = as.dendrogram (clust1) # create dendrogram
  ddgram2 = as.dendrogram (clust2) # create dendrogram
  ggtree_plot1 = ggtree::ggtree(ddgram1) + ggtree::layout_dendrogram()
  ggtree_plot2 = ggtree::ggtree(ddgram2) 
  
  results.df$LR = factor (results.df$LR, levels = unique(results.df$LR)[clust2$order])
  results.df$celltype_pair = factor (results.df$celltype_pair, levels = unique(results.df$celltype_pair)[clust1$order])
  results.df$log10pval[results.df$log10pval < .95] = 0
  results.df$celltype_pair = gsub ('_','',results.df$celltype_pair)
  dplot = ggplot(data = results.df, mapping = aes(x=celltype_pair, y=LR)) +
    geom_point(shape = 21, aes(fill = diffprop, size = log10pval), color='black') +

    scale_fill_gradient2(
      low = "#67A9CF", 
      mid = "#F7F7F7", 
      high = "#EF8A62", 
      midpoint = 0
    ) +
    theme_minimal() +
    theme(text = element_text(size=12), strip.text = element_text(size=22,face='bold')) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_y_discrete(position = "right") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  ggtree_plot_col = ggtree_plot1 + xlim2(dplot)
  ggtree_plot = ggtree_plot2 + ylim2(dplot)
  
  png (paste0("Plots/dotplot_lig_rec_differential_PTENL_GFP_centered_on_selected_CPI.png"), width=3500, height=2100, res=300, pointsize=30)
  print (plot_spacer() + plot_spacer() + ggtree_plot_col + 
    plot_spacer() + plot_spacer() + plot_spacer() +
    ggtree_plot + plot_spacer() + dplot +
    plot_layout(ncol = 3, widths = c(0.4, 0, 4), heights = c(0.4, 0, 4)))
  dev.off()
    
