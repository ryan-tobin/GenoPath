options(repos = c(CRAN = "https://cloud.r-project.org"))

install.packages("remotes")
library(remotes)

remotes::install_github("https://github.com/raerose01/deconstructSigs",force=TRUE)
