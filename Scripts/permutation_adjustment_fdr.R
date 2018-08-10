library(modelr)
load("R_data/allmol_noNAs.rda")
library(qvalue)
require(gtools)

#create list of all other permutations
true <- as.data.frame(matrix(1:10, nrow =2, ncol = 5))
nulls <- as.data.frame(combinations(10,5))
combos <- as.matrix(anti_join(nulls, true))
take.null <- function(x){
    
}
ind2 <- setdiff(1:10, ind1)
y <- x[ind1]
z <- x[ind2]
t <- t.test(y, z, na.action=na.omit)
(t[["statistic"]][["t"]])
combos[2,]
?setdiff
ind1 <- as.vector(combos[1,])



test <- take.null(allmol_noNAs$A001)
allmol_noNAs$Alpha <- rownames(allmol_noNAs)
test2 <- permute(allmol_noNAs, n = 10, .id = "Alpha")

#General Form
for (i in range_of_values) {
  operations_that_use_i_which_is_changing_across_values
}
m = 25
#create an empty vector
s_n <- vector(length = m)
ind1 <- vector(length = 5)
#create for loop to find product for 1:25
for(n in 1:m){
  s_n[n] <- compute_s_n(n)
}

samplev <- 1:10
ind <- vector(length = 10)
test <- for(i in 1:10){
  result <- samplev[i] + 3
  result
}

ts.ps <- mapply(take.t, eico)

take.t <- function(x){
  y <- ind1
  z <- ind2
  t <- t.test(y, z, na.action=na.omit)
  t[[3]]
}

apply.t <- function(x){
  ind1 <- 
  mapply(take.t, allmol_noNAs[,-1])
}
R> apply(M, 1, function(x) 2*x[1]+x[2])


 
















#the true nulls need to be uniformly distributed for an empirical null stat  
pvals <- mapply(take.t, allmol_noNAs[,-1])
#looks good! (ish)
hist(pvals)


test_sub <- allmol_noNAs[, 2:100]
sample.vector <- sample(1:100, 10)

#multiply 5 elements by 5 elements


#bootstrap method
take.t <- function(x){
  ind1 <- sample(1:10, 5)
  ind2 <- setdiff(1:10, ind1)
  t <- t.test(x[ind1], x[ind2], na.action=na.omit)
  abs(t[["statistic"]][["t"]])
}  

take.pvals <- function(x){
  y <- x[c(1,3,5,7,9)]
  z <- x[c(2,4,6,8,10)]
  t <- t.test(y, z, na.action=na.omit)
  abs(t[["p.value"]])
}
take.truet <- function(x){
  y <- x[c(1,3,5,7,9)]
  z <- x[c(2,4,6,8,10)]
  t <- t.test(y, z, na.action=na.omit)
  abs(t[["statistic"]][["t"]])
}

require(qvalue)

B <- 100
null_stats <- replicate(B, mapply(take.t, allmol_noNAs[,-1]))

obs_stats <- mapply(take.truet, allmol_noNAs[,-1])
statnought <- t(null_stats)
pvals <- mapply(take.pvals, allmol_noNAs[,-1])

pvalues <- empPvals(stat = obs_stats, stat0 = null_stats, pool = FALSE)
hist(pvalues)
head(pvalues)
sorted <- sort(pvalues)
head(sorted, 200)
summary(qobj)
max(qvalues[qobj$pvalues <= 0.008])
qvalues <- qobj$qvalues
qobj <- qvalue(p = stats)
hist(qobj[["qvalues"]])
hist(pvalues)
summary(qobj)

hist(qobj)
?empPvals
require(gtable)
require(ggplot2)
