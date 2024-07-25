args <- commandArgs(trailingOnly = TRUE)

comm_file <- args[1]
tree_file <- args[2]
abundance_weighted <- as.logical(args[3])
output_file <- args[4]
newick_output_file <- args[5]
newick_output_file

if (!requireNamespace("picante", quietly = TRUE)) install.packages("picante")
if (!requireNamespace("ape", quietly = TRUE)) install.packages("ape")
library(picante)
library(ape)
suppressPackageStartupMessages(library(picante))
suppressPackageStartupMessages(library(ape))

comm <- read.table(comm_file, header=TRUE, row.names =1,sep="\t")
tree <- read.tree(tree_file)
phydist <- cophenetic(tree)

comdistA.result <- comdist(comm, phydist, abundance.weighted=abundance_weighted)

mpdTree <- nj(comdistA.result)
mpdTree

write.tree(mpdTree, file = newick_output_file)

pdf(output_file)
plot(mpdTree, main="NJ Tree from MPD (comdist) Values")
dev.off()
