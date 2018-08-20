require(Rtsne)
#load data
load("R_data/allmol_noNAs.rda")
#subset data to eicosanoids
eico <- allmol_noNAs %>% dplyr::select(grep("A", colnames(allmol_noNAs))) 
#create matrix and vectors needed for Rtsne
x <- as.matrix(eico)
y <- allmol_noNAs$VF
vf <- rep(c("VF", "WT"), 5)
#run tsne
tsne <- Rtsne(x, dims = 2, perplexity=3, verbose=TRUE, max_iter = 500)
#create formatted dataframe for tsne
tsne.dat <- data.frame(Genotype = as.factor(vf), 
                       tsne1 = as.numeric(tsne$Y[,1]), 
                       tsne2 = as.numeric(tsne$Y[,2]))
#plot results
tsne.dat %>% ggplot(aes(tsne.dat$tsne1, tsne.dat$tsne2)) + 
    geom_point(aes(color = Genotype)) + xlab("tsne 1") + ylab("tsne 2")

