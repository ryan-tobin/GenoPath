args <- commandArgs(trailingOnly = TRUE)

comm_file <- args[1]
tree_file <- args[2]
output_file <- args[3]
newick_output_file <- args[4]

if (!requireNamespace("picante", quietly = TRUE)) install.packages("picante")
if (!requireNamespace("ape", quietly = TRUE)) install.packages("ape")
library(picante)
library(ape)
suppressPackageStartupMessages(library(picante))
suppressPackageStartupMessages(library(ape))


comm <- read.table(comm_file, header=TRUE, row.names =1,sep="\t")
tree <- read.tree(tree_file)
comm
tree
unifracValues.result <- unifrac(comm, tree)

unifracTree <- nj(unifracValues.result)

write.tree(unifracTree, file = newick_output_file)

pdf(output_file)
plot(unifracTree, main="NJ Tree from Unifrac Values")
dev.off()
