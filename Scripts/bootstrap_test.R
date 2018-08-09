library(MASS)
library(boot)
library(readr)
require(tidyverse)
require(DMwR)
require(ggfortify)
require(caret)
#load data
load("R_data/allmol_noNAs.rda")

#generate more rows with KNN and bootstrapping
allmol.vf <- allmol_noNAs %>% filter(allmol_noNAs$VF == 1)
smote.vf <- SMOTE(VF~., allmol.vf, perc.over = 100, k =3)
smote.all <- SMOTE(VF~., allmol_noNAs, perc.over = 600, k = 3, perc.under = 800)
smote.all <- smote.all %>% filter(!is.na(A001))
smote.all <- unique(smote.all)
smote.all$VF <- as.factor(smote.all$VF)
smote.train <- smote.all %>% anti_join(allmol_noNAs)
allmol_noNAs$VF <- rep(c(1,0), 5)
allmol_noNAs$VF <- as.factor(allmol_noNAs$VF)
duplicated(smote.train)
table(smote.vf$VF)
#create train and test sets
set.seed(3456)
trainIndex <- createDataPartition(smote.all$VF, p = .5, 
                                  list = FALSE, 
                                  times = 2)
#create test and train matrices
x.train <- as.matrix(smote.all[trainIndex,-1])
y.train <- as.numeric(smote.all[trainIndex, 1])
x.test <- as.matrix(smote.all[-trainIndex,-1])
y.test <- as.numeric(smote.all[-trainIndex, 1])
x.train$VF
str(smote.all)
str(unique(smote.all))
#exclude real observations
x.train <- smote.all %>% anti_join(allmol_noNAs)

cv.out <- cv.glmnet(x.train,y.train,alpha= 1,family="binomial",type.measure = "mse")
autoplot(cv.out, label = FALSE)

#fit model with parameters
fit <- glmnet(x.train, y.train, alpha = 0.8, lambda = cv.out[["lambda.1se"]])

#predict test set based on model
lasso_prob <- predict(cv.out,newx = x.test,s=cv.out[["lambda.1se"]], type = "response")

#check accuracy
lasso_predict <- rep("neg",nrow(x.test))
lasso_predict[lasso_prob>.5] <- "pos"
#confusion matrix
table(pred=lasso_predict,true= ifelse(y.test == 2, "pos", "neg"))

#find coefficients of model with lambda = lambda.1se
lasso.coef  <- predict(fit, type = 'coefficients', s =  cv.out[["lambda.1se"]])

#extract coefficients
tmp_coeffs <- coef(fit, s = "lambda.1se")
datcoefs <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], 
                       coefficient = tmp_coeffs@x)

prob <- 1/(1+exp(datcoefs$coefficient))
prob
#ROC curve
pred_y <- as.numeric(lasso_prob > 0.5)
true_y <- (y.test - 1)
true_pos <- (true_y == 1) & (pred_y == 1)
true_neg <- (true_y == 0) & (pred_y == 0)
false_pos <- 
idx <- order(-lasso_prob)
recall <- cumsum()
############
x <- as.matrix(allmol_noNAs[,-1])
y <- allmol_noNAs$VF
x.e <- as.matrix(smote.eico[,-1])
y.e <- as.numeric(smote.eico$VF)
eico <- allmol_noNAs %>% 
  select(VF, grep("A", colnames(allmol_noNAs))) %>% 
  mutate(VF = as.factor(VF))
#subset on eico
smote.ecio <- SMOTE(VF~., eico, perc.over = 400, k = 5, perc.under = 400)
smote.eico <- smote.ecio %>% filter(!is.na(A001))
###########

logistic_model <- glmnet(VF~., alpha = 1)

trainSplit$target <- as.factor(trainSplit$target)
allmol_noNAs[,1:10]
allmol_noNAs$VF <- as.factor(allmol_noNAs$VF)
smote.mol <- SMOTE(VF ~ ., allmol_noNAs[-1,], k = 5, perc.under = 100,  perc.over = 100)

summary(smote.mol)
trainSplit$target <- as.numeric(trainSplit$target)

??smote
my.mod <- glm(VF~., data = allmol_noNAs, family = "binomial")
coef(my.mod)
summary(my.mod)

Coefficients:
  Estimate Std. Error z value Pr(>|z|)    
(Intercept) -3.989979   1.139951  -3.500 0.000465 ***
  gre          0.002264   0.001094   2.070 0.038465 *  
  gpa          0.804038   0.331819   2.423 0.015388 *  
  rank2       -0.675443   0.316490  -2.134 0.032829 *  
  rank3       -1.340204   0.345306  -3.881 0.000104 ***
  rank4       -1.551464   0.417832  -3.713 0.000205 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Set up the non-parametric bootstrap

logit.bootstrap <- function(data, indices) {
  
  d <- data[indices, ]
  fit <- glm(admit ~ gre + gpa + rank, data = d, family = "binomial")
  
  return(coef(fit))
}

set.seed(12345) # seed for the RNG to ensure that you get exactly the same results as here

logit.boot <- boot(data=mydata, statistic=logit.bootstrap, R=10000) # 10'000 samples

logit.boot

#Bootstrap Statistics :
  original        bias    std. error
#t1* -3.989979073 -7.217244e-02 1.165573039
#t2*  0.002264426  4.054579e-05 0.001146039
#t3*  0.804037549  1.440693e-02 0.354361032
#t4* -0.675442928 -8.845389e-03 0.329099277
#t5* -1.340203916 -1.977054e-02 0.359502576
#t6* -1.551463677 -4.720579e-02 0.444998099

# Calculate confidence intervals (Bias corrected ="bca") for each coefficient

boot.ci(logit.boot, type="bca", index=1) # intercept

boot.ci(logit.boot, type="bca", index=2) # gre

boot.ci(logit.boot, type="bca", index=3) # gpa

boot.ci(logit.boot, type="bca", index=4) # rank2

boot.ci(logit.boot, type="bca", index=5) # rank3

boot.ci(logit.boot, type="bca", index=6) # rank4


a <- rnorm(100, mean=0, sd = 1)
b <- rnorm(100, mean=-1, sd =2)
mean(a) - mean(b)
m <- length(a)
n <- length(b)
d <- c(a,b)
e <- sample(d)
a2<- e[1:m]
b2 <-  e[ (m+1) : (m+n) ]
mean(a2) - mean(b2)
 

