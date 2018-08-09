require(spls)
require(tidyverse)
require(ggfortify)
require(knitr)
require(plsgenomics)
require(glmnet)
#load data and create response and predictor matrices
load("R_data/allmol_noNAs.rda")
x <- as.matrix(allmol_noNAs[,-1])
y <- allmol_noNAs$VF


#cross validate for hidden parameters and plot result
cv.out <- cv.glmnet(x,y,alpha= 1,family="binomial",type.measure = "mse")
autoplot(cv.out, label = FALSE)

#fit model with parameters
fit <- glmnet(x, y, alpha = 1, lambda = cv.out[["lambda.1se"]])

#find coefficients of model with lambda = lambda.1se
lasso.coef  <- predict(fit, type = 'coefficients', s =  cv.out[["lambda.1se"]])

#extract coefficients
tmp_coeffs <- coef(fit, s = "lambda.min")
datcoefs <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], 
                       coefficient = tmp_coeffs@x)

write.csv(datcoefs, "Lasso_coefs.csv")

#########


##using train and test sets
train <- c(sample(seq(1,9,2), 3), sample(seq(2,10,2),3))
##I used train: c(1,9,7,8,2,6) for the powerpoint slides
train <- c(1,9,7,8,2,6)

x_train <- x[train,]
y_train <- y[train]
x_test <- x[-train,]
y_test <- y[-train]

#cross validate for hidden parameters and plot result
cv.train <- cv.glmnet(x_train,y_train,alpha= 1,family="binomial",type.measure = "mse")
autoplot(cv.train, label = FALSE)

#fit model with parameters
fit <- glmnet(x_train, y_train, alpha = 1, lambda = cv.train[["lambda.1se"]])

#find coefficients of model with lambda = lambda.1se
lasso.coef  <- predict(fit, type = 'coefficients', s =  cv.train[["lambda.1se"]])

#extract coefficients
tmp_coeffs_train <- coef(fit, s = "lambda.1se")
datcoefs_train <- data.frame(name = tmp_coeffs_train@Dimnames[[1]][tmp_coeffs_train@i + 1], 
                       coefficient = tmp_coeffs_train@x)

#predict test set based on model
lasso_prob <- predict(cv.train,newx = x_test,s=cv.train[["lambda.1se"]], type = "response")

#check accuracy
lasso_predict <- rep("neg",nrow(x_test))
lasso_predict[lasso_prob>.5] <- "pos"
#confusion matrix
table(pred=lasso_predict,true= ifelse(y_test == 1, "pos", "neg"))

######

#run lasso on eicosanoid subset

load("R_data/allmol_noNAs.rda")
eico <- allmol_noNAs %>% dplyr::select(VF, grep("A", colnames(allmol_noNAs))) 

x_eic <- as.matrix(eico[,-1])
y_eic <- ifelse(eico$VF == "VF", 1, 0)
#cross validate for hidden parameters and plot result
cv.out <- cv.glmnet(x_eic,y_eic,alpha= 0.5,family="binomial",type.measure = "mse")
autoplot(cv.out, label = FALSE)

#fit model with parameters
fit <- glmnet(x_eic, y_eic, alpha = 0.5, lambda = cv.out[["lambda.min"]])

#find coefficients of model with lambda = lambda.1se
lasso.coef  <- predict(fit, newx = x_eic,  type = 'response', s =  cv.out[["lambda.min"]])

#extract coefficients
tmp_coeffs <- coef(fit, s = "lambda.min")
datcoefs <- data.frame(Alphanumeric = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], 
                       coefficient = tmp_coeffs@x)

#check sanity
lasso_predict <- rep("neg",nrow(x))
lasso_predict[lasso.coef>.5] <- "pos"
#confusion matrix
table(pred=lasso_predict,true= ifelse(y == 1, "pos", "neg"))

#create list
detail.list <- datcoefs %>% 
  merge(lookup, by = "Alphanumeric") %>%
  mutate(Probability = exp(coefficient)/(1+exp(coefficient)), 
         mz = round(mz, 5), RT = round(RT, 4)) %>% 
  rename(Coefficient = coefficient, ID = IDs) %>%
  unite(mzid, mz, RT, sep = "_") %>%
  dplyr::select(ID, Coefficient, Probability, mzid) %>%
  arrange(-Probability)

write.csv(detail.list, "Figures/Lasso_coefs_eicosanoids.csv")


######old code ^^^^stuff I've actually put in the powerpoint (8/9)

#standard logistic regression with glm

#setting up dataset - make mouse IDs rownames and remove IDs column
allmol_df <- as.data.frame(allmol)
allmol_df$VF <- rep(c("VF", "WT"), 5)
rownames(allmol_df) <- allmol_df$IDs
#remove ID column
allmol_df <- allmol_df[,2:ncol(allmol_df)]
#make all zeros NAs and remove incomplete cases
allmol_df[allmol_df == 0] <- NA
allmol_complete <- allmol_df[ , apply(allmol_df, 2, function(x) !any(is.na(x)))]


#set seed to ensure reproducible results
set.seed(42)
#separate training and test sets
trainset <- allmol_complete[train,]
testset <- allmol_complete[-train,]
head(allmol_complete[train,])
glm_model <- glm(VF~.,data = trainset, family = binomial)
summary(glm_model)
#predict probabilities on testset
glm_prob <- predict.glm(glm_model,testset[,-1],type="response")
glm_predict <- rep(0,nrow(testset))
glm_predict[glm_prob>.5] <- 1
table(pred=glm_predict,true=testset$VF)
mean(glm_predict==testset$VF)

##another option (copied from old plasma markdown file)
library(glmnet)
load(allmol_complete)
#X predictors and Y responses
predic <- model.matrix(VF~., allmol_complete)
#remove columns containing NAs

responses <- allmol$VF
#create lambda sequence to find best
lambda <- 10^seq(10, -2, length = 100)
#create train set
set.seed(489)
train = c(sample(seq(1,9,2), 2), sample(seq(2,10,2),2))
#create test set
test = (-train)
response_test <- responses[test]

#ordinary least squares
mol.lm <- lm(VF~., data = allmol_df)
coef(mol.lm)

#ridge
ridge.mod <- glmnet(predic[train,], responses[train], alpha = 0, lambda = lambda)[1,6,]
predict(ridge.mod, s = 0, exact = T, type = 'coefficients')
ridge.mod
plot(ridge.mod$lambda, ridge.mod$dev.ratio)

lasso.mod <- glmnet(predic[train,], responses[train], alpha = 1, lambda = lambda)
lasso.pred <- predict(lasso.mod, s = bestlam, newx = x[test,])


