# AB332 2026
# Please make sure that you have the required R packages installed.
# You can get a list of the package already installed on your computer by executing

installed.packages()

# in R alternatively you can just run the installation for each package
# to make sure that you have the latest version.
# Run these installation commands line-by-line in R (or Rstudio)
# and answer "yes" if you are asked to update any previously installed pakages:
# If you are asked if you want to compile from source it is usually easiest to answer "no"
install.packages("ggplot2")   # To do plots
install.packages("tidyverse") # To manipulate and visualize data https://www.tidyverse.org/packages/
install.packages("vegan")     # For ecology analyses
install.packages("ape")            # Phylogenetic tools
install.packages("recluster")      # Clustering tools
install.packages("dendextend")     # To work with dendrograms
install.packages("corrplot")       # makes nice correlation plots
install.packages("RcmdrMisc")      # diverse tools
install.packages("RColorBrewer") 

# Packages from Bioconductor: (https://www.bioconductor.org/)
if (!requireNamespace("BiocManager", quietly = TRUE))  install.packages("BiocManager")
BiocManager::install(c("PCAtools","mixOmics", "phyloseq"))
