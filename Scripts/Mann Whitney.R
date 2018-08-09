require(tidyverse)
#load dataset
load("~/R_Files/Projects/Jak2VFMetabolomics/R_data/allmol_noNAs.rda")
load("~/R_Files/Projects/Jak2VFMetabolomics/R_data/lookup.rda")
#subset to eicosanoids
eico <- allmol_noNAs %>% select_if(grepl("A", colnames(allmol_noNAs))) 
#construct function to take wilcox tests and calculate means
take.mann <- function(x){
  y <- rep(c(1,0), 5)
  Mean_VF <- mean(x[seq(1,9,2)])
  Mean_WT <- mean(x[c(2,4,6,8,10)])
  mw <- wilcox.test(x~y, conf.int = TRUE)
  c(mw[["estimate"]][["difference in location"]]/10^5, 
    mw[["statistic"]][["W"]], mw[["conf.int"]]/10^5, 
    mw[["p.value"]])
}

#apply to all columns
mw.ps <- mapply(take.mann, eico)

#transpose to create vertical list
list.mw <- as.data.frame(t(mw.ps))
colnames(list.mw) <- c("Difference", "W", 
                       "CI - Lower Bound", "CI - Upper Bound", "P-value")
#add lookup key
list.mw$Alphanumeric <- rownames(list.mw)
list.full <- merge(list.mw, lookup, by = "Alphanumeric")
rownames(list.mw) <- NULL
#correct column order and arrange by p-value
list.full <- list.full %>% 
  unite(mzid, mz, RT, sep = "_") %>%
  dplyr::select(IDs, Difference, W, `CI - Lower Bound`,
                `CI - Upper Bound`, `P-value`, mzid) %>%
  arrange(`P-value`)

list.full %>% ggplot(aes(`P-value`)) + 
  scale_x_continuous(trans = "log10") +
  geom_histogram(bins = 30) +
  geom_vline(xintercept = 0.05)
#write the csv
write.csv(list.full, file = "Figures/MWEicosanoidsNoMeans.csv")



