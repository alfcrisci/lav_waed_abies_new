##########################################################################
# load libraries

#######################################################################################################


library(tidyverse)
library(ggplot2)
library(ggord)
library(knitr)
library(gridExtra)
library(cowplot)
library(openxlsx)
library(heplots) # Visualizing Hypothesis Tests in Multivariate Linear Models

########################################################################

library(cluster)
library(factoextra)
library(flexclust)
library(hopkins)
library(ClusterR)
library(clusterSim)
library(ordr)
library(ordr.extra)
library(MASS)
library(FactoMineR)         # Multivariate Exploratory Data Analysis and Data Mining
library(CCA)                # Canonical correlation analysis
library(klaR)               # Classification and Visualization
library(caret)
library(matrixTests)
library(multiColl)
library(spls)               # Sparse Partial Least Squares (SPLS) Regression and Classification
library(performance)        # Assessment of Regression Models Performance
library(smacof)             #  Multidimensional Scaling
library(paran)              # Horn's Test of Principal Components/Factors
library(moments)            # skewness, kurtosis and related tests

########################################################################

library(plyr); 
library(dplyr)
library(BBmisc)
library(flextable)          # beautifying tables
library(psych)              # psychological research: descr. stats, FA, PCA etc.
library(skimr)              # summary stats
library(gtsummary)          # publication ready summary tables
library(ggpubr)             # publication ready data visualization in R
library(ggcorrplot)
library(kableExtra)
########################################################################


############################################################################################
# functions

standardize <- function(x) {(x - mean(x))/sd(x)}

stat_groups=function(x,g,xlabv="Groups",ylabv="Compound") {
  myls=list();
  myls[[1]]=aov(x ~ g)
  myls[[2]]=shapiro.test(myls[[1]]$residuals)
  # not normal if 0.05 >
  myls[[3]]=bartlett.test(x ~ g) # not ok if 0.05 >
  # not homogeneous if 0.05 >
  myls[[4]]=kruskal.test(x ~ g)
  # differs 0.05 >
  myls[[5]]=boxplot(x ~ g, xlab=xlabv, ylab=ylabv)
  return(myls)
}

outersect <- function(x, y, ...) {
  big.vec <- c(x, y, ...)
  duplicates <- big.vec[duplicated(big.vec)]
  setdiff(big.vec, unique(duplicates))
}

#######################################################################################################
# useful functions

cal_z_score <- function(x){(x - mean(x)) / sd(x)}



############################################################################################
# references


# [1] https://bioconductor.org/packages/devel/bioc/vignettes/PCAtools/inst/doc/PCAtools.html

