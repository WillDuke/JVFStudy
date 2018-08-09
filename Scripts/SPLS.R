#load packages
require(spls)
require(tidyverse)
require(ggfortify)
require(knitr)
#load data and create train and test sets
load("R_data/allmol_noNAs.rda")
# ***for predict: train = c(sample(seq(1,9,2), 3), sample(seq(2,10,2), 3))
x <- as.matrix(allmol_noNAs[,-1])
y <- ifelse(allmol_noNAs$VF == "VF", 1, 0)
#cross validate to determine spls parameters K and eta
cvsp <- cv.spls(x,y, K = c(1:3), eta = seq(0.1, 0.9, 0.05), fold = 5)
#create logistic model using logistic (alt: LDA)
cv.s
fit <- sgpls(x, y, K = 0.8, eta = 3, kappa=0.5, scale.x=TRUE)
fitspls <- spls(x, y, lambda.l1 = )
#print variable selection
print(fitspls)
ci <- ci.spls(fit, B = 400)
ci.f <- ci.spls( fit, plot.it=TRUE, plot.fix='x', plot.var=20 )
print(ci[])
#sanity check
test <- predict.splsda(fit, x, type = c("fit"), fit.type = c("class"))
test
#coefficients
confint.lm(fit)
coefs <- data.frame(ID = rownames(coefs), coef = coefs)
coef.noZ <- coefs %>% filter(coef != 0) %>% arrange(desc(coef))
save(coef.noZ, file = "SPLSDA_SelectVars.rda")
?confint
confint(fit)
ci <- ci.spls(fitsgpls, coverage = 0.95, B = 100)

ecoef.fit <- coef(fit)
coef.noZ <- data.frame(coef = coef.fit[ coef.fit!=0, ], ID = names(coef.fit[ coef.fit!=0, ])) %>%
    arrange(desc(coef))
rownames(coef.noZ) <- NULL
fit
coefs %>% filter(ID == "B004")
summary(fit)

c.coef <- coef(fit)
nonzero <- c.coef[c.coef != 0]
nonzero 
coef <- predict.spls(f, newx = allmol_noNAs[-train,-1], type="coefficient")


citest <- ci.spls(fit, coverage=0.95, B=200, plot.it = TRUE)

correct.spls(fit, plot.it=TRUE )

pred.acc <- ifelse(pred.f >0.5, 1, 0)
table(pred = pred.acc, true = allmol_noNAs$VF[-train])
mean(pred.acc==allmol_noNAs$VF[-train])
#try all train combos
strapped <- replicate(dim(combn(10,2))[[2]], {
  train = c(sample(seq(1,9,2), 3), sample(seq(2,10,2), 3))
  cv <- cv.splsda(allmol_noNAs[train,-1],allmol_noNAs$VF[train], 
                  K = c(1:3), eta = seq(0.1, 0.9, 0.1), fold = 5)
  fit <- splsda(x, y, eta = as.numeric(cv[["eta.opt"]]), K = cv[["K.opt"]], kappa=0.5, 
              classifier = "logistic", scale.x=TRUE, scale.y=FALSE)
  pred.fit <- predict.spls(fit, newx = allmol_noNAs[-train,-1], type="fit")
  pred.acc <- ifelse(pred.f >0.5, 1, 0)
  mean(pred.acc==allmol_noNAs$VF[-train])
})

#try all combos -- raw averages
strapped.raw <- replicate(dim(combn(10,2))[[2]], {
  train = c(sample(seq(1,9,2), 3), sample(seq(2,10,2), 3))
  cv <- cv.spls(allmol_noNAs[train,-1],allmol_noNAs$VF[train], K = c(1:3), eta = seq(0.1, 0.9, 0.1), fold = 5)
  fit <- spls(x, y, eta = as.numeric(cv[["eta.opt"]]), K = cv[["K.opt"]], kappa=0.5, 
              select="pls2", fit="simpls", scale.x=TRUE, scale.y=FALSE)
  predict.spls(fit, newx = allmol_noNAs[-train,-1], type="fit")
})
strapped.raw

fit <- spls(x, y, eta = as.numeric(cv[["eta.opt"]]), K = cv[["K.opt"]], kappa=0.5, 
            select="pls2", fit="simpls", scale.x=TRUE, scale.y=FALSE)
table(strapped)

pred.f <- data.frame(Mouse = rownames(pred.f), Prediction = pred.f, row.names = NULL)
pred.tab <- kable(pred.f)
strapraw <- as.data.frame(strapped.raw)
mean(strapped)
rowMeans(strapraw)
coef <- as.data.frame(coef)
coef.org <- coef %>% filter(V1 != 0) %>% arrange(desc(abs(V1)))
coef.f <- coef(f)
coef.nonzero <- coef.f[coef.f!=0, ]       
coef.nonzero <- as.data.frame(coef.nonzero)

coef.nonzero %>% arrange() %>% head(50)


C <- 20
model_instability <- replicate(C, {
  cv <- cv.splsda(x,y, K = c(1:3), eta = seq(0.1, 0.9, 0.05), fold = 5)
  c(cv[["eta.opt"]], cv[["K.opt"]])
})
sd(model_instability[,1])




