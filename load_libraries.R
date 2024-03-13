# load libraries
message ('Load R packages')
packages = c(
  'Seurat',
  'ggplot2',
  'gplots',
  'ggpubr',
  'patchwork',
  'ComplexHeatmap',
  'ggrepel',
  'fgsea',
  'harmony',
  'clusterProfiler',
  'rstatix',
  'paletteer',
  'dplyr',
  'tidyr',
  'ggtree',
  'aplot'
)
lapply(packages, require, character.only = TRUE)

