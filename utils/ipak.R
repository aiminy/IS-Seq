
ipak <- function(pkg, repository=c('CRAN', 'Bioconductor', 'github')){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  # new.pkg <- pkg
  if (length(new.pkg)) {
    if(repository=='CRAN') {
      if (!requireNamespace(new.pkg,quietly = TRUE))
        install.packages(new.pkg,repos = "http://cran.us.r-project.org",dependencies = TRUE)
    }
    if(repository=='Bioconductor') {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager",repos = "http://cran.us.r-project.org")
      if (!requireNamespace(new.pkg, quietly = TRUE))
        BiocManager::install(new.pkg)
    }
    if(repository=='github') {
      if (!requireNamespace(new.pkg, quietly = TRUE))
        devtools::install_github(new.pkg, build_vignettes = FALSE, dependencies=TRUE)
    }
  }
}