### Under construction

We will use R (version 4.0.5 or later) and Rstudio (version 1.4.1 or later) in this course.

Ideally everybody should download and installed R (https://www.r-project.org/) and Rstudio (https://www.rstudio.com/) before this part of the course starts. If you prefer to run R in another IDE (for instance Jupyter Notebook, Visual Studio Code or similar) that is also fine, but Rstudio is probably the easiest to get started with if you have not used R before. 

It is also preferable that you have all packages installed before the course starts (see below or [Install_required_packages.R](Install_required_packages.R). For those of you having trouble with installing the packages, we will try to find some time during the first computer lab to get everybody up and running.

**NB: Some packages require administrator access to be installed. If you are working on a computer administered by somebody else (for instance owned by a University) you might have to ask your IT department to either grant you temporary admin rights or help you install the packages**    

**NB: Setting the working directory to a path inside a One Drive folder can result in errors and I adviced against doing it unless you have done it before and know that it works**

Content of [required packages](Install_required_packages.R):

```R
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
```