options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!require("BiocManager",quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MutationalPatterns",force=TRUE)