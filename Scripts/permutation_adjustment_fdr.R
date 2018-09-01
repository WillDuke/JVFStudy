library(ggalt)
library(qvalue)
library(ggrepel)
library(dplyr)
load("R_data/allmol_noNAs.rda")

#functions to take t-tests of each column and return pvals and t statistics
#bootstrap method - computes t-test between random set of 5 and remaining values
take.t <- function(x){
  ind1 <- sample(1:10, 5)
  ind2 <- setdiff(1:10, ind1)
  t <- t.test(x[ind1], x[ind2], na.action=na.omit)
  abs(t[["statistic"]][["t"]])
}  
#takes true pvals
take.pvals <- function(x){
  y <- x[c(1,3,5,7,9)]
  z <- x[c(2,4,6,8,10)]
  t <- t.test(y, z, na.action=na.omit)
  abs(t[["p.value"]])
}
#takes true absolute t statistics
take.truet <- function(x){
  y <- x[c(1,3,5,7,9)]
  z <- x[c(2,4,6,8,10)]
  t <- t.test(y, z, na.action=na.omit)
  abs(t[["statistic"]][["t"]])
}

#the true nulls need to be uniformly distributed for an empirical null stat  
pvals <- mapply(take.pvals, allmol_noNAs[,-1])
#looks good!
hist(pvals)

#create null absolute t stat matrix
B <- 200 #set replication number
null_stats <- replicate(B, mapply(take.t, allmol_noNAs[,-1]))

#create true absolute t stats
obs_stats <- mapply(take.truet, allmol_noNAs[,-1])

#list of pvals from t tests
t.pvals <- mapply(take.pvals, allmol_noNAs[,-1])

#adjusted pvals based on observed and null statistics
adj.pvals <- empPvals(stat = obs_stats, 
                      stat0 = null_stats, 
                      pool = TRUE)

#compute q-values from adjusted p-values, bootstrapping to find pi0
#(smoother method did not find a good pi0)
qobj <- qvalue(p = adj.pvals, pi0.method = "bootstrap")
#show pval cutoffs and their associated qvals
summary(qobj)
#show plot of pi0 v lambda as well as cutoff plots
plot(qobj)
#show adjusted pvals with pi0, qvals, and lfdr
hist(qobj)
#find fdr at p<0.01
max(qobj$qvalues[qobj$pvalues <= 0.01])

###alt method: directly find molecules below fdr with fdr.level
qobj_fdrlevel <- qvalue(p = adj.pvals, pi0.method = "bootstrap", fdr.level = 0.05)

#create dataframe from qobj object with relevant values
siglist <- data.frame(Alphanumeric = colnames(allmol_noNAs[,-1]),
                      `P-value` = qobj$pvalues, 
                      `Q-value` = qobj$qvalues, 
                      LFDR = qobj$lfdr) 
#merge with lookup 
siglist <- siglist %>% merge(lookup, by = "Alphanumeric")
#create top 20 list with mzids and rounded values
orglist <- siglist %>% 
    mutate(`P.value` = round(`P.value`,8), 
           `Q.value` = round(`Q.value`,3), 
           LFDR = round(LFDR, 3), mz = round(mz,5), RT = round(RT,4)) %>% 
    unite(mzid, mz, RT, sep = "_") %>%
    select(IDs, `P.value`, `Q.value`, LFDR, mzid) %>% 
    arrange(`P.value`) %>%
    head(20)

#need to fix lookup so that the Us are included!

#write orglist to a csv for table production
write.csv(orglist, file = "Tables/PA-FDR_Top20_withUs.csv")

##try the same with data subsetted to just eicosanoids
eico <- allmol_noNAs %>% select_if(grepl("A", colnames(allmol_noNAs)))
eico.p <- mapply(take.pvals, eico)
#histogram reveals uneven distribution of nonsignificant p-vals 
#PA-FDR will be unreliable
hist(eico.p, nclass = 20)

siglist %>% filter(Q.value < 0.05, Alphanumeric == "^A") %>% ggplot(aes(P.value, Q.value)) + geom_point()

take.allt <- function(x){
  y <- x[c(1,3,5,7,9)]
  z <- x[c(2,4,6,8,10)]
  t <- t.test(y, z, na.action=na.omit)
  c(t[["estimate"]][["mean of x"]], 
    t[["estimate"]][["mean of y"]],
    t[["conf.int"]][1], 
    t[["conf.int"]][2], t[[3]])
}
means <- mapply(take.allt, eico[,c(-1, -114)])
means <- as.data.frame(t(means))
colnames(means) = c("MeanWT", "MeanVF", "lower", "upper", "t")

voldat <- siglist %>% 
  merge(lookup, by = "Alphanumeric") %>% 
  cbind(means) %>% 
  mutate(logP = -log10(P.value), logFC = log2(MeanVF/MeanWT)) %>% 
  mutate(colors = ifelse(logFC > 2 | logFC < -2, "Black", "Blue"))

voldat %>% ggplot(aes(logFC, logP, label = IDs)) + 
  geom_point(aes(color = colors)) +
  geom_hline(yintercept = -log10(0.0005), linetype = "dotted") +
  theme(legend.position = "none") + xlab("Log2 Fold Change") +
  ylab("log10 of P-value") + xlim(c(-4,4)) +
  geom_text_repel(data = subset(voldat, 
      logP > -log10(0.0005) & (logFC > 1 | logFC < -1)),
      arrow = arrow(length = unit(0.03, "npc"), type = "closed", ends = "first"), 
      segment.size = 0.1) +
  ggtitle("Volcano Plot of Whole Dataset after PA-FDR") +
  theme(plot.title = element_text(hjust = 0.5))

#for eico subset of **already adjusted** pvals
voldat %>% ggplot(aes(logFC, logP, label = IDs)) + 
  geom_point(color = "#00BFC4") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  ylab("log10 of P-value") + xlim(c(-4,4)) +
  annotate("text", label = "P-value < 0.5", x = 2.85, y = 1.35, size = 4, colour = "black") +
  geom_text_repel(data = subset(voldat, 
      logP > -log10(0.009) & (logFC > 0 | logFC < -.8)),
      arrow = arrow(length = unit(0.02, "npc"), type = "closed", ends = "first"), 
      segment.size = 0.2, force = 2, point.padding = 1, ylim = c(2.3, NA)) +
  theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "none") + xlab("Log2 Fold Change") +
  ggtitle("Volcano Plot of Eicosanoids after PA-FDR")
  
ggsave(file = "Figures/Volcano Plot Eicosanoids After PA-FDR.tiff", units = "in", 
       width = 9, height = 7, dpi = 300)
#create dataframe from qobj object with relevant values
siglist <- data.frame(Alphanumeric = colnames(allmol_noNAs[,-1]),
                        `P-value` = qobj$pvalues, 
                        `Q-value` = qobj$qvalues, 
                        LFDR = qobj$lfdr)
siglist <- siglist %>% filter(grepl("A", Alphanumeric))
  