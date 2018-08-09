install.packages("sgPLS")
require(sgPLS)

require(spls)
data(yeast)
data(swiss)
install.packages("glmnet")
xsample <- model.matrix(Fertility~., swiss)[,-1]
x <- model.matrix(Fertility~., swiss)[,-1]
train = sample(1:nrow(x), nrow(x)/2)
test = (-train)
ytest = y[test]
swisslm <- lm(Fertility~., data = swiss)
coef(swisslm)
y <- swiss$Fertility
#ridge
ridge.mod <- glmnet(x, y, alpha = 0, lambda = lambda)
predict(ridge.mod, s = 0, exact = T, type = 'coefficients')[1:6,]

swisslm <- lm(Fertility~., data = swiss, subset = train)
ridge.mod <- glmnet(x[train,], y[train], alpha = 0, lambda = lambda)
#find the best lambda from our list via cross-validation
cv.out <- cv.glmnet(x[train,], y[train], alpha = 0)

bestlam <- cv.out$lambda.min

#make predictions
ridge.pred <- predict(ridge.mod, s = bestlam, newx = x[test,])
s.pred <- predict(swisslm, newdata = swiss[test,])
#check MSE
mean((s.pred-ytest)^2)

mean((ridge.pred-ytest)^2)

#a look at the coefficients
out = glmnet(x[train,],y[train],alpha = 0)
predict(ridge.mod, type = "coefficients", s = bestlam)[1:6,]

lasso.mod <- glmnet(x[train,], y[train], alpha = 1, lambda = lambda)
lasso.pred <- predict(lasso.mod, s = bestlam, newx = x[test,])
mean((lasso.pred-ytest)^2)

lasso.coef  <- predict(lasso.mod, type = 'coefficients', s = bestlam)[1:6,]
