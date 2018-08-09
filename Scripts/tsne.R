require(Rtsne)
load("R_data/allmol_noNAs.rda")
eico <- eico %>% dplyr::select(grep("A", colnames(allmol_noNAs))) 
x <- as.matrix(eico)
y <- allmol_noNAs$VF
vf <- rep(c("VF", "WT"), 5)

train$label<-as.factor(train$label)
## for plotting
colors = rainbow(length(unique(y)))
names(colors) = unique(y)

tsne.dat <- data.frame(Genotype = as.factor(vf), tsne1 = as.numeric(tsne$Y[,1]), tsne2 = as.numeric(tsne$Y[,2]))

tsne.dat %>% ggplot(aes(tsne.dat$tsne1, tsne.dat$tsne2)) + geom_point(aes(color = Genotype))

tsne <- Rtsne(x, dims = 2, perplexity=3, verbose=TRUE, max_iter = 500)
exeTimeTsne<- system.time(Rtsne(x, 
                  dims = 2, perplexity=2, verbose=TRUE, max_iter = 500))
plot(tsne$Y, t='n', main="T-sne Non-Linear Clustering", ylab = "t-sne 2", xlab = "t-sne 1")
text(tsne$Y, labels= vf, col=colors)
exeTimeTsne

