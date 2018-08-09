#load packages
require(spls)
require(tidyverse)
require(ggfortify)
require(knitr)
#load data and create train and test sets
load("R_data/allmol_noNAs.rda")
lookup <- read_excel("R_data/Plasma_data_lookup_key.xlsx")
# ***for predict: train = c(sample(seq(1,9,2), 3), sample(seq(2,10,2), 3))
x <- as.matrix(allmol_noNAs[,-1])
y <- allmol_noNAs$VF
cvsp <- cv.spls(x,y, K = c(1:3), eta = seq(0.1, 0.9, 0.05), fold = 5)
cv.sgpls(x, y, fold=5, K = c(1:3), eta = seq(0.1, 0.9, 0.1), scale.x=TRUE, 
         plot.it=TRUE, n.core=1)

fit <- sgpls(x, y, K = 1, eta = 0.9)
coef <- coef(fit)
ci <- ci.spls(fit, B = 500, plot.it = TRUE, plot.fix="y")
coef.dat <- as.data.frame(coef)
ci.dat <- as.data.frame(ci[["cibeta"]][[1]])

cis <- data.frame(ID = names(allmol_noNAs)[-1], 
                  coef = ci[["betahat"]][-1], 
                  lower = ci[["lbmat"]], 
                  higher = ci[["ubmat"]])
cis %>% filter(coef != 0)
cis %>% filter(coef != 0) %>% ggplot(aes(cis)) + 
  geom_segment(aes(x = lower, xend=higher,y=1:33, yend = 1:33)) + scale_x_continuous("log10")
predict( object, newx, type = c("fit","coefficient"),
         fit.type = c("class","response"), ... )