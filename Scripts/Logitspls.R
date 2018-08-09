require(spls)
require(tidyverse)
require(ggfortify)
require(knitr)
require(plsgenomics)
#load data and create response and predictor matrices
load("R_data/allmol_noNAs.rda")
x <- as.matrix(allmol_noNAs[,-1])
y <- allmol_noNAs$VF

### hyper-parameters values to test
lambda.l1.range <- seq(0.05,0.95,by=0.1) # between 0 and 1
ncomp.range <- 1:10
# log-linear range between 0.01 a,d 1000 for lambda.ridge.range
logspace <- function( d1, d2, n) exp(log(10)*seq(d1, d2, length.out=n))
lambda.ridge.range <- signif(logspace(d1 <- -2, d2 <- 3, n=21), digits=3)

cv <- logit.spls.cv(x, y, lambda.ridge.range = lambda.ridge.range, 
              lambda.l1.range = lambda.l1.range, ncomp.range = ncomp.range,
              adapt = TRUE, maxIter = 100, svd.decompose = TRUE,
              return.grid = FALSE, ncores = 1, nfolds = 3, nrun = 1,
              center.X = TRUE, scale.X = TRUE, weighted.center = TRUE, seed = NULL,
              verbose = TRUE)

logit.spls(x, y, lambda.ridge = 1, lambda.l1 = 0.1, 1, Xtest = NULL,
           adapt = TRUE, maxIter = 100, svd.decompose = TRUE, center.X = TRUE,
           scale.X = FALSE, weighted.center = TRUE)

  
