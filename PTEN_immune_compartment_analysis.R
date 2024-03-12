conda activate scrnatools 
use UGER
R

source ('PTENL_study/useful_functions.R')
source ('PTENL_study/load_libraries.R')
source ('PTENL_study/ggplot_aestetics.R')

srt = readRDS ('srt.rds')

# Add pallette
trm_pal = setNames (c('black','grey44', 'darkgreen'), c('CB2_PTENL','CB2_PTENL_C124S','GFP')) 
pal_heatmap = c('#67A9CF','#EF8A62')


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

pdf (paste0(projdir,'Plots/cell_composition_',metaGroupName,'.pdf'), width=13, height=3)
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
pdf (paste0(projdir,'Plots/cell_composition_',metaGroupName,'.pdf'), width=7, height=3)
(umap1 | cc_box1)
plot_layout (widths= c(2,6))
dev.off()


results.df = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/Jia_prj/PTENL_demuxEM_seq2_analysis/_cellranger_filtered_Filter_200_500_25/no_harmony/high_quality_subset/sampleID_harmony/immune_subset/sampleID_harmony/CellphoneDBv4_celltypeTIMTAM2/CellphoneDB_results.rds')

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
  ggtree_plot1 = ggtree::ggtree(ddgram1) + layout_dendrogram()
  ggtree_plot2 = ggtree::ggtree(ddgram2) 
  
  results.df3$LR = factor (results.df3$LR, levels = unique(results.df3$LR)[clust2$order])
  results.df3$celltype_pair = factor (results.df3$celltype_pair, levels = unique(results.df3$celltype_pair)[clust1$order])
  results.df3$log10pval[results.df3$log10pval < pValThreshold] = 0

  dplot = ggplot(data = results.df3, mapping = aes(x=celltype_pair, y=LR)) +
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
  if (any (colnames(results.df3) =='sampleID')) dplot = dplot + facet_wrap (~sampleID, ncol=length (unique(results.df3$sampleID)))
  ggtree_plot_col = ggtree_plot1 + xlim2(dplot)
  ggtree_plot = ggtree_plot2 + ylim2(dplot)
  
  if (!is.null(compGroups)) w_plot = 1500+(70*length(unique(results.df3$celltype_pair))) else w_plot =  11500+(70*length(unique(results.df3$celltype_pair)))
  png (paste0(cellphonedb_dir, "Plots/dotplot_lig_rec_differential_PTENL_GFP_centered_on_selected_CPI.png"), width=w_plot, height=500+(60*length(unique(results.df3$LR))), res=300, pointsize=30)
  print (plot_spacer() + plot_spacer() + ggtree_plot_col + 
    plot_spacer() + plot_spacer() + plot_spacer() +
    ggtree_plot + plot_spacer() + dplot +
    plot_layout(ncol = 3, widths = c(0.4, 0, 4), heights = c(0.4, 0, 4)))
  dev.off()
  
  if (!is.null(compGroups)) w_plot = 6+(length(unique(results.df3$celltype_pair))/4) else w_plot =  6 + (length(unique(results.df3$celltype_pair))/4) * length(unique(results.df3$sampleID))
  pdf (paste0("Plots/dotplot_lig_rec_differential_PTENL_GFP_centered_on_selected_CPI.pdf"), width= w_plot, height=4+(length(unique(results.df3$LR))/6))
  print (plot_spacer() + plot_spacer() + ggtree_plot_col + 
    plot_spacer() + plot_spacer() + plot_spacer() +
    ggtree_plot + plot_spacer() + dplot +
    plot_layout(ncol = 3, widths = c(0.4, 0, 4), heights = c(0.4, 0, 4)))
  dev.off()
  

### muscat DS on genotype ###
force = FALSE
do.fgsea = TRUE
logfcThreshold = .5
pvalAdjTrheshold = 0.05
ds_method = "DESeq2" #c("edgeR", "DESeq2", "limma-trend", "limma-voom")
metaGroupName1 = 'sampleID'
metaGroupName2 = 'celltype_TIMTAM'
metaGroupName3 = 'treatment2'
#metaGroupName3 = 'treatment'
#muscatIdents = c('CB2_PTENL','GFP')
#muscatIdents = c('CB2_PTENL','CB2_PTENL_C124S')
muscatIdents = c('CB2_PTENL','control')
pbDS_min_cells = 5
topGenes = 20 
show_genes = c('Cxcr3','Cxcr4','Ccr3', 'Ccr1', 'Cxcl10','Cxcl11','Ccl12','Ccl9','Ccl5','Ccl8','Ccl6') # check genes
#srt2 = srt
#srt = subset (srt, sampleID2_harmony_snn_res.1 %in% c(0,1,2,3,4,5,6,7))
source ('/ahg/regevdata/projects/ICA_Lung/Bruno/scripts/scrna_pipeline/DS_muscat.R')

dhm = diffClustHeat ( 
    deg_genes_df = tbl_df,
    # mat1 = srt_Ident1@assays$RNA@data,
    # mat2 = srt_Ident2@assays$RNA@data,
    # meta1 = srt_Ident1@meta.data[,metaGroupName1],
    # meta2 = srt_Ident2@meta.data[,metaGroupName1],
    sort_by = 'avg_log2FC',
    #sort_by = 'p_val_adj',
    topGenes = topGenes, 
    pvalAdjTrheshold = pvalAdjTrheshold,
    # addGene = addGene,
    #column_title = feat,
    col_limit = 3,
    plotcol = rev(RColorBrewer::brewer.pal(3,'RdBu')),
    #plotcol = rev(viridis::turbo(3)),
    name = paste (muscatIdents, collapse = '-'),
    cluster_columns=T,
    cluster_rows=T,
    #row_title = paste (deg2Ident, collapse = ' vs '),
    border=T)
  max_nrow = nrow (dhm@matrix)
  max_ncol = ncol (dhm@matrix)    
  max_nrow2 = max (unlist(max_nrow)) / 20 + 2
  max_ncol2 = ifelse (max ((unlist(max_ncol)) / 5.5) < 3, 3, max (unlist(max_ncol)) / 10)
  pdf (paste0(projdir_ms,'Plots/',topGenes,'_genes_heatmap2.pdf'), width = max_ncol2 +2, height= max_nrow2)
  draw (dhm, heatmap_legend_side = "right")
  dev.off() 

topGenes=20
vln_p = lapply ("Macs", 
  function(x) plotExpression(sce[, sce$cluster_id == x],
  features = head(top_deg_genes_df$gene[top_deg_genes_df$cluster == x], topGenes),
  x = "sample_id", colour_by = "group_id", ncol = 3) +
  guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle (paste('cluster', x)))
pdf (paste0(projdir_ms,'Plots/top_',topGenes,'_genes_in_Macs_vln_plot2.pdf'), 10,10)
print (vln_p)
dev.off()

reac_terms = unique(unlist(lapply (fgseaResAll[[1]], function(x) x$pathway) ))
mhc_terms = reac_terms[grep ('MHC',reac_terms)]
mhc_genes = unique(unlist(lapply (fgseaResAll [[1]], function(x) x[x$pathway %in% mhc_terms,'leadingEdge'])))
tbl_df2 = tbl_df[!is.na(tbl_df$p_val_adj),]
tbl_df2 = tbl_df2[tbl_df2$p_val_adj < 0.05, ]
mhc_sig_genes = tbl_df2[tbl_df2$gene %in% mhc_genes, ]
write.csv (mhc_sig_genes, paste0(projdir,'mhc_sig_genes.csv'))

# Check genes in pathways
#muscatIdents = c('CB2_PTENL','CB2_PTENL_C124S')
muscatIdents = c('CB2_PTENL','control')
projdir_ms = paste0(projdir,'muscat_',metaGroupName1,'_',metaGroupName2,'_',metaGroupName3,'_',paste(muscatIdents,collapse='_vs_'),'_method_',ds_method,'/')
fgseaResAll = readRDS (paste0(projdir_ms, 'fgsea_annotation_',paste(gmt_annotations,collapse='_'),'.rds'))
fgseaRanks = readRDS (paste0(projdir_ms, 'fgsea_ranks_',paste(gmt_annotations,collapse='_'),'.rds'))

# Plot dotplot of fGSEA annotations per cluster 
top_pathways = Inf
int_pathways = c('interferon','antigen','cytokine')
int_pathways2 = unlist (lapply (fgseaResAll[['c5.bp.v7.1.symbol.gmt']], function(x) x$pathway))
int_pathways3 = unique(unlist(lapply (int_pathways,function(x) int_pathways2[grep (x, int_pathways2)])))
fgseaResAll2 = lapply (fgseaResAll[['c5.bp.v7.1.symbol.gmt']], function(x) x[x$pathway %in% int_pathways3,])
  fgseaResAll_dp = dotGSEA (fgseaResAll2, padj_threshold = pvalAdjTrheshold, 
    type = 'fgsea',top_pathways = top_pathways,
    cluster_rows=T,
    cluster_cols=T)
    
png (paste0(projdir,'Plots/fGSEA_annotation_c5.bp.v7.1.symbol.gmt_dotplots2.png'),2800,1000, res=300)
print(fgseaResAll_dp)
dev.off()

unique (unlist(fgseaResAll_dp$data$pathway))  
pathway = c('GO_antigen_processing_and_presentation','GO_antigen_processing_and_presentation_of_peptide_antigen')
celltype = 'DCs'
pathway = c('GO_cytokine_production','GO_cellular_response_to_interferon_gamma','GO_T_cell_cytokine_production','GO_cytokine_production_involved_in_immune_response','GO_response_to_cytokine','GO_positive_regulation_of_cytokine_production','GO_cellular_response_to_cytokine_stimulus')
celltype = c('MregDCs','TAM2','TIM_NC','DCs','TIM_C')
fgsea1 = as.data.frame (do.call (rbind,fgseaResAll2[celltype]))
fgsea1[fgsea1$pathway %in% pathway,'leadingEdge']

gene2 = unique(unlist(lapply(fgseaResAll2, function(x) x[grep ('antigen', x$pathway), 'leadingEdge'])))
gene2 = unique (unlist (fgsea1[fgsea1$pathway %in% pathway,'leadingEdge']))

p = lapply (gene, function(x) VlnPlot (srt, x, split.by = 'treatment', group.by = 'celltype_TIMTAM'))
pdf (paste0(projdir,'Plots/',pathway[1],'_genes_violin_plot2.pdf'), width=9, height=5)
p
dev.off()

# tbl_df2 = tbl_df[tbl_df$gene %in% gene, ]
# pvalAdjTrheshold=1.1
# dhm = diffClustHeat (
#     deg_genes_df = tbl_df2,
#     # mat1 = srt_Ident1@assays$RNA@data,
#     # mat2 = srt_Ident2@assays$RNA@data,
#     # meta1 = srt_Ident1@meta.data[,metaGroupName1],
#     # meta2 = srt_Ident2@meta.data[,metaGroupName1],
#     sort_by = 'avg_log2FC',
#     #sort_by = 'p_val_adj',
#     topGenes = topGenes, 
#     pvalAdjTrheshold = pvalAdjTrheshold,
#     # addGene = addGene,
#     #column_title = feat,
#     col_limit = 3,
#     plotcol = rev(RColorBrewer::brewer.pal(3,'RdBu')),
#     #plotcol = rev(viridis::turbo(3)),
#     name = paste (muscatIdents, collapse = '-'),
#     cluster_columns=T,
#     cluster_rows=T,
#     #row_title = paste (deg2Ident, collapse = ' vs '),
#     border=T)
#   pdf (paste0(projdir,'Plots/',pathway[1],'_genes_heatmap2.pdf'), width = 6, height= 3)
#   draw (dhm, heatmap_legend_side = "right")
#   dev.off() 

metagroup_df = data.frame (
 barcode = colnames(srt),
 condition = srt@meta.data[,metaGroupName3], 
 cluster = srt@meta.data[,metaGroupName2])
metagroup_df = metagroup_df[metagroup_df$condition %in% c(muscatIdents[1],muscatIdents[2]),]
  
agr = srt@assays$RNA@data[,metagroup_df$barcode]
agr = agr[rownames (agr) %in% gene2, ]
agr = as.data.frame (t(agr))
agr2 = agr
agr2$cluster = paste0(srt$treatment2, '__', srt$celltype_TIMTAM)
agr2 = aggregate (.~cluster, agr2, mean)
rownames (agr2) = agr2$cluster
col_split = sapply (agr2$cluster, function(x) unlist(strsplit (x, '__'))[[2]])
agr2 = agr2[,-1]
#agr = (apply(t(agr), 1, function(x)(x-min(x))/(max(x)-min(x))))

ha2 = HeatmapAnnotation (condition = metagroup_df$condition, col=list(condition=setNames(pal_heatmap,c(muscatIdents[1],muscatIdents[2]))))
allcells_hm = Heatmap (t(agr), 
  cluster_rows=T, 
  cluster_columns=T,
  column_title_rot = 90,
  #column_title = paste0(pathway, ' DEG LFC >',logfcThreshold2, 'pval < ',pvalAdjTrheshold2), 
  #column_title_side = 'bottom',
  column_split = paste0(metagroup_df$condition, ' ',metagroup_df$cluster),
  #column_split = metagroup_df$cluster,
  #top_annotation = ha2,
  top_annotation = ha2,
  width = 1,
  #col=dichromat::colorschemes$DarkRedtoBlue.18,
  #col=RColorBrewer::brewer.pal (9,'Reds'),
  col=viridis::turbo(100),
  row_names_gp = gpar(fontsize = 14), 
  column_names_gp = gpar (fontsize = 0),
  border=T)

pdf (paste0(projdir,'Plots/',pathway[1],'_heatmap_cells2.pdf'), width = 16, height= 10)
  draw (allcells_hm, heatmap_legend_side = "right")
  dev.off()

agr3 = agr2[grep ('PTENL',rownames(agr2)),] - agr2[grep ('control',rownames(agr2)),] 
rev(RColorBrewer::brewer.pal(3,'RdBu'))
col_fun = colorRamp2(c(-1, 0, 1), c("#67A9CF", "#F7F7F7", "#EF8A62"))
#agr3 = max (abs(agr3))
avg_hm = Heatmap (t(agr3), 
  cluster_rows=T, 
  cluster_columns=T,
  column_title_rot = 90,
  #column_title = paste0(pathway, ' DEG LFC >',logfcThreshold2, 'pval < ',pvalAdjTrheshold2), 
  #column_title_side = 'bottom',
  column_split = rownames(agr3),
  #column_split = metagroup_df$cluster,
  #top_annotation = ha2,
  #top_annotation = ha2,
  width = 1,
  #col=dichromat::colorschemes$DarkRedtoBlue.18,
  #col=RColorBrewer::brewer.pal (9,'Reds'),
  #col=viridis::turbo(10),
  col = col_fun,
  row_names_gp = gpar(fontsize = 7), 
  column_names_gp = gpar (fontsize = 6),
  border=T)
png (paste0(projdir,'Plots/',pathway[1],'_genes_heatmap_avg2.png'), width = 1050, height= 15700, res=300)
avg_hm
dev.off() 


gene = c('Cxcr3','Cxcr4','Ccr3', 'Ccr1', 'Cxcl10','Cxcl11','Ccl12','Ccl9','Ccl5','Ccl8')
gene1 = c('Il1b','Cxcl10','Il6','Il10','Ccl22','Ccl4','Ccl5','Ccl3','Ifng')
gene_selection = gene [gene %in% colnames (agr3)]
gene_selection = append (gene_selection,colnames (agr3)[grep ('Ifn', colnames(agr3))])
gene_selection = append (gene_selection,colnames (agr3)[grep ('Tnf', colnames(agr3))])
gene_selection = append (gene_selection,colnames (agr3)[grep ('Il', colnames(agr3))])
gene_selection = append (gene_selection,colnames (agr3)[grep ('Ccl', colnames(agr3))])
gene_selection = append (gene_selection,colnames (agr3)[grep ('Cxcl', colnames(agr3))])
gene_selection = append (gene_selection, gene1)

agr = srt@assays$RNA@data[,metagroup_df$barcode]
agr = agr[rownames (agr) %in% gene_selection, ]
agr = as.data.frame (t(agr))
agr2 = agr
agr2$cluster = paste0(srt$treatment2, '__', srt$celltype_TIMTAM)
agr2 = aggregate (.~cluster, agr2, mean)
rownames (agr2) = agr2$cluster
col_split = sapply (agr2$cluster, function(x) unlist(strsplit (x, '__'))[[2]])
agr2 = agr2[,-1]
agr3 = agr2[grep ('PTENL',rownames(agr2)),] - agr2[grep ('control',rownames(agr2)),] 

#agr = (apply(t(agr), 1, function(x)(x-min(x))/(max(x)-min(x))))

avg_hm = Heatmap (t(agr3[,unique(gene_selection)]), 
  cluster_rows=T, 
  cluster_columns=T,
  column_title_rot = 90,
  #column_title = paste0(pathway, ' DEG LFC >',logfcThreshold2, 'pval < ',pvalAdjTrheshold2), 
  #column_title_side = 'bottom',
  column_split = rownames(agr3),
  #column_split = metagroup_df$cluster,
  #top_annotation = ha2,
  #top_annotation = ha2,
  width = 1,
  #col=dichromat::colorschemes$DarkRedtoBlue.18,
  #col=RColorBrewer::brewer.pal (9,'Reds'),
  #col=viridis::turbo(10),
  col = col_fun,
  row_names_gp = gpar(fontsize = 7), 
  column_names_gp = gpar (fontsize = 6),
  border=T)
png (paste0(projdir,'Plots/Cytokines_genes_heatmap_avg_refined.png'), width = 1050, height= 2700, res=300)
avg_hm
dev.off() 




pathway = 'HALLMARK_INTERFERON_GAMMA_RESPONSE'
fgsea1 = as.data.frame (fgseaResAll[['h.all.v7.1.symbol.gmt']][['DCs']])
fgsea1[fgsea1$pathway == pathway,'leadingEdge']


# Check chemokines
metaGroupName = 'sampleID'
metaGroupName2 = 'treatment'
gene = c('Cxcr3','Cxcr4','Ccr3', 'Ccr1', 'Cxcl10','Cxcl11','Ccl12','Ccl9','Ccl5','Ccl8','Ccl6')
p = geneDot (
  gene = gene,
  y = srt$celltype_TIMTAM,
  x = 'sampleGenotype',
  z = srt$treatment2,
  min_expression = 0)

pdf (paste0(projdir, 'Plots/cytokines_expression_dotplot.pdf'))
p
dev.off()

# Run WGCNA
# Set variables
force=FALSE # force re-running WGCNA 
do.fgsea=TRUE
powerTable=FALSE # plot the powerTable plot. Can take a while to generate
do.plots=TRUE
softPower=5 # Set the softPower
deepSplit=1 # Set this to increase decrease the number of modules idendified. 1-4
mergeCutHeight = 0.25 # Height below which two modules are merged
metacells_k = 30 # number of cells to create pseudobulks
max_shared = 20
metacells_groups = c('sampleID','treatment') # set metagroup in which to find metacells (usually your clustering / celltypes)
metaGroupNames = c('sampleID','celltype_wgcna','treatment') # set of metagroups for generating the boxplots
module_pal = trm_pal
minModuleSize = 30
#genes.keep = VariableFeatures (FindVariableFeatures (srt, nfeat = 10000)) # Genes to use to compute WGCNA
genes.keep = VariableFeatures (FindVariableFeatures (srt, nfeat=5000))
enricher_universe = rownames(srt) # Genes to use as background for pathway enrichments analysis
source (paste0(scrna_pipeline_dir,'scWGCNA.R'))

tbl2 = tbl_df[tbl_df$gene %in% modules[modules$module == 'turquoise','gene_name'],]
tbl2 = tbl2[!is.na(tbl2$p_val_adj),]
tbl2 = tbl2[tbl2$p_val_adj < 0.05,]

tbl2 = tbl_df[tbl_df$gene %in% modules[modules$module == 'pink','gene_name'],]
tbl2 = tbl2[!is.na(tbl2$p_val_adj),]
tbl2 = tbl2[tbl2$p_val_adj < 0.05,]


# Subset for malignant and stromal comp
force=TRUE
subclustername = 'momacs'
#metaGroupName='RNA_snn_res.0.8'
#metaGroupSelection=c(2,3,4,6,12,15,16,18,20,21)
metaGroupName='celltype'
metaGroupSelection=c('Macs','Macs_proliferating','Mono1','Mono2','Mono_Spp1')
exclude=FALSE
source (paste0(scrna_pipeline_dir,'subcluster.R'))

# Subset for malignant and stromal comp
force=FALSE
subclustername = 'TNKcells'
#metaGroupName='RNA_snn_res.0.8'
#metaGroupSelection=c(2,3,4,6,12,15,16,18,20,21)
metaGroupName='celltype'
metaGroupSelection=c('TNKcells')
exclude=FALSE
source (paste0(scrna_pipeline_dir,'subcluster.R'))
