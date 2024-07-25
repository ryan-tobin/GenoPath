options(repos = c(CRAN = "https://cloud.r-project.org"))

install.packages("remotes")
library(remotes)

remotes::install_url("https://www.ncbi.nlm.nih.gov/CBBresearch/Przytycka/software/signatureestimation/SignatureEstimation.tar.gz",force=TRUE)