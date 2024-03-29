
### Useful functions ###


# Dotplot showing expression of a gene across two meta groups. Need a Seurat obj
geneDot = function (
  seurat_obj = srt,
	#mat_norm = srt@assays$RNA@data,
	#mat_counts = srt@assays$RNA@counts,
	gene = NULL, 
	x = NULL, # Vector of metagroup1 of equal length to mat columns, If multiple genes are specified this is ignored and genes would make the x axis of dotplot instead
	y = NULL, # vector of metagroup2 of equal length to mat columns
	x_name = 'genes',
  assay='RNA',
	y_name = 'clusters',
	min_expression = 0,
	facet_ncol = 5,
	lim_expression = NULL,
	scale.data = TRUE, # scale data when multiple genes are given
	plotcol = viridis::viridis(100),
	include_NA = TRUE,
  swap_axes = FALSE,
  returnDF = FALSE
	#... # arguments to pass to facet wrap
	)
	{
  require ("scCustomize")
  require ('tidyr')
  if (exists('levels_x')) rm ('levels_x')
  if (exists('levels_y')) rm ('levels_y')
	if (all (grepl ('^\\d', seurat_obj@meta.data[,y]))) seurat_obj@meta.data[,y] = paste0('C',seurat_obj@meta.data[,y])
  if (is.factor (seurat_obj@meta.data[,x])) levels_x = levels (seurat_obj@meta.data[,x]) else
  levels_x = unique (seurat_obj@meta.data[,x])
  if (is.factor (seurat_obj@meta.data[,y])) levels_y = levels (seurat_obj@meta.data[,y]) else 
  levels_y = unique (seurat_obj@meta.data[,y])
  
  seurat_obj@meta.data[,x] = gsub ('[-_ \\+]', '', seurat_obj@meta.data[,x])
  seurat_obj@meta.data[,y] = gsub ('[-_ \\+]', '', seurat_obj@meta.data[,y])
  
  if (exists ('levels_x')) levels_x = gsub ('[-_ \\+]', '', levels_x)
  if (exists ('levels_y')) levels_y = gsub ('[-_ \\+]', '', levels_y)
    
  # Compute percentage expression per cell group
  percent = Percent_Expressing (seurat_object = seurat_obj, assay=assay, threshold = min_expression, features = gene, group_by = x, split_by = y)
  percent$gene = rownames(percent)
  percent = gather (percent, groups, expression, 1:(ncol(percent)-1))
  percent$key = paste0(percent$gene, '_', percent$groups)
  colnames (percent)[colnames(percent) == 'expression'] = 'percent'
  
  # Compute average expression per cell group
  seurat_average = log2 (as.data.frame (AverageExpression (seurat_obj, assay=assay, features= gene, group.by = c(x,y))[[1]]) + 1)
  if (length(gene) == 1) rownames(seurat_average) = gene
  seurat_average$gene = rownames(seurat_average)
  seurat_average = gather (seurat_average, groups, expression, 1:(ncol(seurat_average)-1))
  seurat_average$key = paste0(seurat_average$gene, '_', seurat_average$groups)
  
  # Merge percentage and average expression data.frames
  dot.df = cbind (seurat_average, percent[match (seurat_average$key, percent$key),])
  dot.df = dot.df[,!duplicated (colnames(dot.df))]

  # scale data if scale is TRUE
  if (scale.data) dot.df = transform (dot.df, expression = ave (expression, gene, FUN = function(x) scale(x, scale=T, center=T)))

  dot.df$expression[dot.df$percent == 0] = NA
	dot.df$percent[dot.df$percent == 0] = NA

  dot.df$x_axis = factor (sapply (dot.df$groups, function(x) unlist(strsplit(x, '_'))[1]), levels = levels_x)
  dot.df$y_axis = factor (dot.df$gene, levels = gene) 
  if (!is.null (y)) dot.df$y_axis = factor (sapply (dot.df$groups, function(x) unlist(strsplit(x, '_'))[2]), levels = levels_y)
	if (!is.null (y) & length(gene) > 1) 
    {
    dot.df$y_axis = factor (dot.df$gene, levels = gene) 
    dot.df$z_axis = factor (sapply (dot.df$groups, function(x) unlist(strsplit(x, '_'))[2]), levels = levels_y)
    }

  #dot.df$x_axis = factor (dot.df$x_axis, levels = dot.df$x_axis)
  if (swap_axes)  colnames (dot.df)[match (c('x_axis','y_axis'), colnames (dot.df))] = c('y_axis','x_axis')

	p = ggplot (data = dot.df, aes (x= x_axis, y= y_axis)) +
  	geom_point (shape=21, aes (fill= expression, size = percent), alpha=0.7,colour='black', stroke=0.3) +
    labs (x = x_name, y = y_name, title = ifelse(length(gene) > 1,'',gene), subtitle = paste('Min expression >', min_expression)) +
     scale_shape (solid = FALSE) +
  	theme_classic () +
  	theme (axis.text.x = element_text (angle = 45, vjust = 1, hjust=1))
  if (!is.null (y) & length(gene) > 1) p = p + facet_wrap (~z_axis)  
  if (length (plotcol) > 3) p = p + scale_fill_gradientn (colors = plotcol) else
    p = p + scale_fill_gradient2 (low = plotcol[1], mid = plotcol[2], high=plotcol[3])	

	if (returnDF) return(dot.df) else 
  return (p)
	}

# Generate barplots or boxplots of meta groups proportions specified
# meta_groups: vector or meta group names:
# 1) first meta_group is the group on which proportions are calculated
# 2) second meta_group split first meta_group on x axes  
#	3) third meta_group will group barplots separately
# if splits include only one value runs barplot instead
cellComp = function (
	seurat_obj = NULL, 
	metaGroups = NULL, # vector of at least 3 metaGroups e.g. c('orig.ident','celltypes','celltypes'),
	plot_as = 'box', # box or bar 
	pal = NULL,
	prop = TRUE,
	ptable_factor = 1, # specify which column of the data.frame or seurat object metadata should be used to compute proportions
	facet_ncol = 20,
	facet_scales = 'free',
	subset_prop = NULL, # subset prop table by any group in any column
	removeNA = TRUE,
	returnDF = FALSE
	)
	{
	require (ggplot2)	
	require (ggpubr)	
	if (is.data.frame (seurat_obj))
		{
		meta_groups_df = seurat_obj[,metaGroups]	
		} else {
		meta_groups_df = seurat_obj@meta.data[,metaGroups]
		}
	# Refactor to remove 0 groups
	#meta_groups_df =  as.data.frame(lapply(unclass(meta_groups_df),as.character),stringsAsFactors=T)
	if(is.null(pal)) pal = rainbow (length(unique(meta_groups_df[,2])))
	#if(is.null(pal) & plot_as == 'box') pal = rainbow (length(unique(meta_groups_df[,3])))
 	if (prop)
 		{
 		ccomp_df = as.data.frame (prop.table (table (meta_groups_df),ptable_factor))
 		ccomp_df = na.omit (ccomp_df) # this is to remove NaN somehow produced from the line above 
 		} else {
 		ccomp_df = as.data.frame (table (meta_groups_df))	
 		}

 	if(removeNA) ccomp_df = ccomp_df[ccomp_df$Freq != 0, ] # remove 0s from proportions
 	if (!is.null (subset_prop)) 
 		{
 		subset_col = unlist(sapply (seq(ncol(ccomp_df)), function(x) if(any(ccomp_df[,x] %in% subset_prop)) colnames(ccomp_df)[x]))
 		ccomp_df = ccomp_df[ccomp_df[,subset_col] %in% subset_prop,]
 		}
 	#colnames (ccomp_df) = c(paste0('Var_',seq_along(metaGroups)), 'proportion')  
 	if (plot_as == 'box')
 		{
 		p = ggplot (ccomp_df, aes_string (x= metaGroups[2], y= 'Freq')) +
  			theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    		scale_fill_manual (values= pal) + xlab (metaGroups[2]) + ylab (ifelse (prop, 'proportion','counts'))
    if (length(metaGroups > 2)) p = p + geom_boxplot(aes_string (fill= metaGroups[3]), outlier.size=.2, alpha = 0.7, lwd=.2) 
    else p = p + geom_boxplot(aes_string (fill= metaGroups[2]), outlier.size=.2, alpha = 0.7, lwd=.2)		
  	if (length(metaGroups) > 3) p = p + facet_wrap (as.formula(paste("~", metaGroups[4])), scales=facet_scales, ncol=facet_ncol)
  	}	
  if (plot_as == 'bar')
 		{
 		p = ggplot (ccomp_df, aes_string (x= metaGroups[1], y= 'Freq')) +
  			geom_bar(position="stack", stat="identity", aes_string(fill= metaGroups[2])) +
    		#geom_bar(position="dodge", stat="identity", aes_string(fill= metaGroups[2])) +
    		theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    		scale_fill_manual (values= pal) + xlab (metaGroups[2]) + ylab (ifelse (prop, 'proportion','counts'))
  	if (length(metaGroups) == 3) p = p + facet_wrap (as.formula(paste("~", metaGroups[3])), scales=facet_scales, ncol=facet_ncol)
  	}
  if (returnDF) return(ccomp_df) else 
  return (p)
	}
	



# Build clustered dotplot with GSEA terms on y axes and clusters on x axes
dotGSEA = function (
	enrichmentsTest_list,
	type = c('fgsea','enrich'),
	top_pathways = NULL, 
	padj_threshold = 0.05,
	cluster_rows = T,
	cluster_cols = T,
	remove_ns_modules = TRUE)
	{
	require (tidyr)	
	if (type == 'fgsea')
		{
		sig_terms = unique(unlist(sapply (enrichmentsTest_list, function(x) x[x$padj < padj_threshold,'pathway'])))
		if (length(sig_terms) == 0) return (NULL)
		if (!is.null(top_pathways)) sig_terms = unique(unlist(sapply (enrichmentsTest_list, function(x) 
				{
				x =	na.omit(x)
				x = x[x$padj < padj_threshold,]
				x_top = x[order(-x$NES),]
				x_top = head(x_top[order(x_top$padj),'pathway'], top_pathways)
				x_bot = x[order(x$NES),]
				x_bot = head(x_bot[order(-x_bot$padj),'pathway'], top_pathways)
				c(x_top,x_bot)
				})))
		if (length(sig_terms) == 1)
			{
			mat_pvalue = t(as.matrix (sapply (enrichmentsTest_list, function(x) x$pval[match(sig_terms, x$pathway)])))
			mat_sizeLog = t(as.matrix (sapply (enrichmentsTest_list, function(x) x$NES[match(sig_terms, x$pathway)])))	
			} else {
			mat_pvalue = sapply (enrichmentsTest_list, function(x) x$pval[match(sig_terms, x$pathway)])
			mat_sizeLog = sapply (enrichmentsTest_list, function(x) x$NES[match(sig_terms, x$pathway)])	
			}
		}
	if (type == 'enrich')
		{
		sig_terms = na.omit(unique(unlist(sapply (enrichmentsTest_list, function(x) x[x$p.adjust < padj_threshold,'ID']))))
		if (length(sig_terms) == 0) return (NULL)
		if (!is.null(top_pathways)) sig_terms = unique(as.vector(unlist (sapply (enrichmentsTest_list, function(x) na.omit(head(x[x$p.adjust < padj_threshold,'ID'], top_pathways))))))

		if (length(sig_terms) == 1)
			{
			mat_pvalue = t(as.matrix(sapply (enrichmentsTest_list, function(x) x$p.adjust[match(sig_terms, x$ID)])))
			mat_sizeLog = t(as.matrix(sapply (enrichmentsTest_list, function(x) x$Count[match(sig_terms, x$ID)])))	
			} else {
			mat_pvalue = sapply (enrichmentsTest_list, function(x) x$p.adjust[match(sig_terms, x$ID)])
			mat_sizeLog = sapply (enrichmentsTest_list, function(x) x$Count[match(sig_terms, x$ID)])	
			}
		}
	
	mat_pvalue[is.na(mat_pvalue)] = 1
	mat_sizeLog[is.na(mat_sizeLog)] = 0
	
rownames (mat_pvalue) = sig_terms
colnames (mat_pvalue) = names (enrichmentsTest_list)
rownames (mat_sizeLog) = sig_terms
colnames (mat_sizeLog) = names (enrichmentsTest_list)
mat_sizeLog = as.data.frame (mat_sizeLog)
mat_pvalue = as.data.frame (mat_pvalue)

	if(cluster_rows & min (dim(mat_sizeLog)[1]) >= 2) 
		{
		d_row = dist (mat_sizeLog)
		d_row[is.na(d_row)] = max (d_row)	
		hc1_row = hclust (d_row, method = "ward.D")
		mat_sizeLog = mat_sizeLog[hc1_row$order,]
		mat_pvalue = mat_pvalue[hc1_row$order,]
		}
	if(cluster_cols & min (dim(mat_sizeLog)[2]) >= 2) 
		{
		d_col = dist (t(mat_sizeLog))	
		d_col[is.na(d_col)] = max (d_col)
		hc1_col = hclust (d_col, method = "ward.D")
		mat_sizeLog = mat_sizeLog[,hc1_col$order]
		mat_pvalue = mat_pvalue[,hc1_col$order]
		}
#mat_sizeLog = as.data.frame (mat_sizeLog)
mat_sizeLog[mat_pvalue > padj_threshold] = NA			
mat_pvalue[mat_pvalue > padj_threshold] = NA
mat_pvalue = -log10(as.data.frame (mat_pvalue))	
#colnames (mat_pvalue) = names (enrichmentsTest_list)

mat_pvalue$pathway = rownames (mat_pvalue)
if (remove_ns_modules) # remove modules with no enriched pathways
	{
	remove_cols	= apply (mat_pvalue, 2, function(x) !all (is.na(x)))
	mat_pvalue = mat_pvalue[, remove_cols, drop=F]
	mat_sizeLog = mat_sizeLog[, colnames(mat_pvalue)[colnames(mat_pvalue)!='pathway'], drop=F]
	}
if (type == 'fgsea')
	{
	gsea_df_pvalue = gather (mat_pvalue, cluster, nlogPvalue, colnames(mat_pvalue)[1]:colnames(mat_pvalue)[ncol(mat_pvalue)-1], factor_key=TRUE)
	gsea_df_NES = gather (mat_sizeLog, cluster2, NES, colnames(mat_sizeLog)[1]:colnames(mat_sizeLog)[ncol(mat_sizeLog)], factor_key=TRUE)
	gsea_df = cbind (gsea_df_NES, gsea_df_pvalue)
	#gsea_df$cluster = factor (gsea_df$cluster, levels = unique(gsea_df$cluster))
	gsea_df$pathway = factor (gsea_df$pathway, levels = unique(gsea_df$pathway))
	p = ggplot (data = gsea_df, aes (x=cluster, y=pathway,
	      fill= NES, size=nlogPvalue)) +
	      geom_point (shape=21, color='black') +
	      #scale_fill_gradient2 (low = 'blue',mid='white', high = 'red', midpoint=0) +
	      scale_fill_gradient2(
      	low = "#67A9CF", 
      	mid = "#F7F7F7", 
      	high = "#EF8A62", 
      	midpoint = 0
    		) +
	      labs (x = 'cluster', y = '-log10(pvalue)') +
	      theme_minimal() + 
	      theme(text = element_text(size=11)) +
	      theme(axis.text.y = element_text (angle = 0, vjust = 0.5, hjust=1)) +
	      theme(axis.text.x = element_text (angle = 90, vjust = 0.5, hjust=1))
	return (p)	
	}
if (type == 'enrich')
	{	
	enrich_df_pvalue = gather (mat_pvalue, cluster, nlogPvalue, colnames(mat_pvalue)[1]:colnames(mat_pvalue)[ncol(mat_pvalue)-1], factor_key=TRUE)
	enrich_df_Count = gather (mat_sizeLog, cluster2, Count, colnames(mat_sizeLog)[1]:colnames(mat_sizeLog)[ncol(mat_sizeLog)], factor_key=TRUE)
	enrich_df = cbind (enrich_df_Count, enrich_df_pvalue)
	#gsea_df$cluster = factor (gsea_df$cluster, levels = unique(gsea_df$cluster))
	enrich_df$pathway = factor (enrich_df$pathway, levels = unique(enrich_df$pathway))
	p = ggplot (data = enrich_df, aes (x=cluster, y=pathway,
	      fill= nlogPvalue, size=Count)) +
	      geom_point (shape=21, color='black') +
	      scale_fill_distiller (palette = "Spectral") +
	      labs (x = 'cluster', y = '-log10(pvalue)') +
	      theme_minimal() + 
	      theme(text = element_text(size=11)) +
	      theme(axis.text.y = element_text (angle = 0, vjust = 0.5, hjust=1)) +
	      theme(axis.text.x = element_text (angle = 90, vjust = 0.5, hjust=1))
	return (p)		
	}
	}


# Make heatmap of the expression difference of between two conditions (e.g. KO vs WT) aggregated by cluster and top DEG across the clusters
diffClustHeat = function (
	deg_genes_df = NULL, # it has to contains following columns: "p_val_adj", "avg_log2FC", gene", "cluster"
	topGenes = 20,
	pvalAdjTrheshold = 0.05,
	pvalAdjTrheshold2 = 1e-50,
	pvalAdjTrheshold3 = 1e-200,
	sort_by = 'p_val_adj',
	only_pos = F,
	# mat1 = NULL, # normalised counts mats (in seurat found in seuratobj@assays$RNA@data)
	# mat2 = NULL, # normalised counts mats (in seurat found in seuratobj@assays$RNA@data)
	# meta1 = NULL,
	# meta2 = NULL,
	addGene = NULL, # add a gene of interest in the heatmap
	plotcol = NULL,
	col_limit = NULL,
	return_mat = F,
	...
	)
	{
	require (ComplexHeatmap)
	require (circlize)	
	
	deg_genes_df1 = deg_genes_df
	sig_genes = unique (deg_genes_df1$gene[deg_genes_df1$p_val_adj < pvalAdjTrheshold])
	deg_genes_df = deg_genes_df1[deg_genes_df1$gene %in% sig_genes, ]
	#if (only_pos) deg_genes_df = deg_genes_df1[deg_genes_df1$gene %in% sig_genes & deg_genes_df1$avg_log2FC > 0, ]
	
	if (sort_by == 'p_val_adj')
		{
		geneNames = unique (unlist (lapply (split (deg_genes_df, deg_genes_df$cluster), function(x) 
		{
		x = x[order (x[,sort_by]),]	
		head(x, topGenes)$gene
		})))
		}
if (sort_by == 'avg_log2FC')
		{
		geneNames = unique (unlist (lapply (split (deg_genes_df, deg_genes_df$cluster), function(x) 
		{
		x = x[order (-x[,sort_by]),]
		if (only_pos) 
			{
			head(x, topGenes)$gene
			} else {
			c(head(x, topGenes)$gene, tail(x, 2)$gene)				
			}	
		})))
		}
	# Add gene specified in 'addGene' argument if needed
	if (!is.null(addGene)) geneNames = unique(append (geneNames, addGene))
	# Make mat of pvalues
	pval_mat = split (deg_genes_df1, deg_genes_df1$cluster)
	pval_mat = as.data.frame(sapply (pval_mat, function(x) x$p_val_adj[match (geneNames, x$gene)]))
	if (ncol(pval_mat) == 1) pval_mat = as.data.frame(t(pval_mat))
	rownames(pval_mat) = geneNames
	#pval_mat = pval_mat[!duplicated(pval_mat$gene),]
	#rownames(pval_mat) = pval_mat$gene
	#pval_mat = pval_mat[top_genes, ]
	pval_mat[is.na(pval_mat)] = 1

		lfc_mat = split (deg_genes_df1, deg_genes_df1$cluster)
		lfc_mat = as.data.frame (sapply (lfc_mat, function(x) x$avg_log2FC[match (geneNames, x$gene)]))
		if (ncol(lfc_mat) == 1) lfc_mat = as.data.frame(t(lfc_mat))
		rownames (lfc_mat) = geneNames
		#rownames (mat1) = top_genes
		lfc_mat[is.na(lfc_mat)] = 0

		if (is.null(col_limit)) col_limit = max (abs(lfc_mat)) 
		if (is.null(plotcol)) 
			{	
			plotcol = colorRamp2(c(-max (abs(lfc_mat)), 0, max (abs(lfc_mat))), c("green", "white", "red"))
			} else {
			plotcol = colorRamp2(c(-col_limit,0, col_limit), plotcol)
			}
		#pal_corr2 = colorRamp2(c(-.5, 0, .5), c("blue", "grey33", "red"))
    labels_annotation = HeatmapAnnotation(
        text = anno_text(rownames(lfc_mat), rot = 45, location = unit(1, "npc"), just = "right",gp =gpar(fontsize = 4)),
        annotation_height = max_text_width(rownames(lfc_mat))
    )
		ht = Heatmap (t(lfc_mat), 
		clustering_distance_columns = 'euclidean',
		clustering_distance_rows = 'euclidean',
		column_title = paste('* = < ',pvalAdjTrheshold,'top',topGenes,'genes'),
		col = plotcol,
	
		bottom_annotation = labels_annotation,
		column_names_gp = gpar(fontsize = 0),
		row_names_gp = gpar(fontsize = 5),
		cell_fun = function (j, i, x, y, width, height, fill) 
						{
		       if (t(pval_mat)[i, j] < pvalAdjTrheshold3)
		       		{
		           grid.text("***", x, y, just='center', vjust=.8,
		           	gp = gpar(fontsize = 5, col='black'))
		          } else {
		        	if(t(pval_mat)[i, j] < pvalAdjTrheshold2)
		        			{
    							grid.text("**", x, y, just='center', vjust=.8,
    							gp = gpar(fontsize = 5, col='black'))   
			          	} else {
    							if(t(pval_mat)[i, j] < pvalAdjTrheshold)
    								{
    								grid.text("*", x, y, just='center', vjust=.8,
    								gp = gpar(fontsize = 5, col='black'))   			
			          		}}}
			}, ...)
	if (!return_mat) return (ht) 
	else return (list (lfc = lfc_mat,pval = pval_mat))
	}

