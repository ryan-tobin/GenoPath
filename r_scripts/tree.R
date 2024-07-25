args <- commandArgs(trailingOnly = TRUE)

clone_tree <- args[1]
tumor_tree <- args[2]
presence_file <- args[3]
output <- args[4]

if (!requireNamespace("ape", quietly = TRUE)) install.packages("ape")
library(ape)
suppressPackageStartupMessages(library(ape))

row_tree <- read.tree(tumor_tree)
col_tree <- read.tree(clone_tree)

row_tree <- chronos(row_tree)
col_tree <- chronos(col_tree)

dist_matrix <- as.matrix(dist.nodes(row_tree))
min_distances <- apply(dist_matrix, 1 , min, na.rm = TRUE)
isolated_tip_index <- which.max(min_distances)
isolated_tip_label <- row_tree$tip.label[isolated_tip_index]
row_tree <- root(row_tree, outgroup = isolated_tip_label, resolve.root=TRUE)
cat("Rooted Tumor Tree at isolated tip:", isolated_tip_label, '\n')

row_order <- row_tree$tip.label
col_order <- col_tree$tip.label 

presence <- read.table(presence_file,header=TRUE,sep='\t',row.names=1)

presence_matrix <- as.matrix(presence)

ordered_matrix <- presence_matrix[row_order, col_order]

hclust_rows <- as.hclust(row_tree)
hclust_cols <- as.hclust(col_tree)

dendro_rows <- as.dendrogram(hclust_rows)
dendro_cols <- as.dendrogram(hclust_cols)

custom_colors <- colorRampPalette(c('white','blue','red'))(n=256)

png(filename = output, width=800, height=600)
heatmap(x=ordered_matrix, Rowv = dendro_rows, Colv = dendro_cols, na.rm = TRUE, main='Heatmap of Clone/Tumor Phylogenies from Presence Data',col=custom_colors)
dev.off()
