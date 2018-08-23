#load packages
require(spls)
require(tidyverse)
require(ggfortify)
require(knitr)
#load data and create train and test sets
load("R_data/allmol_noNAs.rda")
# ***for predict: train = c(sample(seq(1,9,2), 3), sample(seq(2,10,2), 3))
x <- as.matrix(allmol_noNAs[,-1])
y <- allmol_noNAs$VF
#create logistic model using logistic (alt: LDA) (also used K = 1 in powerpoint)
fit <- splsda(x, y, K = 2, eta = 0.9, kappa=0.5,
              classifier=c('logistic'), scale.x=TRUE)
#print variable selection
print(fit)
print("  plus 148 more...  ")