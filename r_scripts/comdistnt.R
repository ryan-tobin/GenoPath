args <- commandArgs(trailingOnly = TRUE)

comm_file <- args[1]
tree_file <- args[2]
abundance_weighted <- as.logical(args[3])
output_file <- args[4]
newick_output_file <- args[5]

if (!requireNamespace("picante", quietly = TRUE)) install.packages("picante")
if (!requireNamespace("ape", quietly = TRUE)) install.packages("ape")
library(picante)
library(ape)
suppressPackageStartupMessages(library(picante))
suppressPackageStartupMessages(library(ape))


comm <- read.table(comm_file, header=TRUE, row.names =1,sep="\t")
tree <- read.tree(tree_file)
comm
phydist <- cophenetic(tree)
phydist

mntdValues.result <- comdistnt(comm, phydist, abundance.weighted=abundance_weighted)

mntdTree <- nj(mntdValues.result)

write.tree(mntdTree, file = newick_output_file)

pdf(output_file)
plot(mntdTree, main="NJ Tree from MNTD (comdistnt) Values")
dev.off()
