library(dplyr)
#load dataset
library(tidyverse)
library(tidyr)
library(ggalt)
load("R_data/allmol_noNAs.rda")
load("R_data/lookup.rda")
#subset to eicosanoids
eico <- allmol_noNAs %>% dplyr::select(grep("A", colnames(allmol_noNAs))) 

#construct function to take t tests and calculate means
take.t <- function(x){
    y <- x[c(1,3,5,7,9)]
    z <- x[c(2,4,6,8,10)]
    t <- t.test(y, z, na.action=na.omit)
    c(t[["estimate"]][["mean of x"]], 
      t[["estimate"]][["mean of y"]],
      t[["conf.int"]][1], 
      t[["conf.int"]][2], t[[3]])
}
#apply to all columns
ts.ps <- mapply(take.t, eico)

#adj_pvals (not preferred, not used subsequently)
adj_pvals <- p.adjust(p.lasma, method = "fdr")

#transpose to create vertical list
list.mw <- as.data.frame(t(ts.ps))
colnames(list.mw) <- c("WT Mean", "VF Mean", 
                       "CI - Lower Bound", 
                       "CI - Upper Bound", 
                       "P-value")
#add lookup key
list.mw$Alphanumeric <- rownames(list.mw)
list.full <- merge(list.mw, lookup, by = "Alphanumeric")
rownames(list.mw) <- NULL
#correct column order and arrange by p-value
list.full <- list.full %>% 
  mutate(`WT Mean (x 10^5)` = round(`WT Mean`/10^5,2), 
         `VF Mean (x 10^5)` = round(`VF Mean`/10^5,2),
         `CI - Lower Bound (x 10^5)` = round(`CI - Lower Bound`/10^5, 2),
         `CI - Upper Bound (x 10^5)` = round(`CI - Upper Bound`/10^5, 2),
         `P-value` = round(`P-value`, 4)) %>%
  dplyr::select(IDs,
         `WT Mean (x 10^5)`,
         `VF Mean (x 10^5)`, 
         `CI - Lower Bound (x 10^5)`, 
         `CI - Upper Bound (x 10^5)`, 
         `P-value`) %>% 
  arrange(`P-value`)

#write the csv
write.csv(list.full, file = "Figures/TtestEicosanoids.csv")

#Create graphs showing t-tests with errors:
list.full <- list.full %>% 
        mutate(`WT-VF`= `WT Mean (x 10^5)` - `VF Mean (x 10^5)`) %>% 
        rename(Diff = `WT-VF`, 
               lower = `CI - Lower Bound (x 10^5)`, 
               upper = `CI - Upper Bound (x 10^5)`)
list.full <- list.full %>% 
  mutate(FC = (Diff)/`WT Mean (x 10^5)`, 
         lowerFC = lower/`WT Mean (x 10^5)`, 
         upperFC = upper/`WT Mean (x 10^5)`)

#create factor for dumbell
list.lowp <- list.full %>% 
  filter(`P-value` < 0.1) %>% arrange(`P-value`)
list.lowp$FactorID <- factor(list.lowp$IDs, levels=as.character(list.lowp$IDs))
list.lowp$FactorID <- factor(list.lowp$FactorID, levels=rev(levels(list.lowp$FactorID)))

#create dumbbell plot - save as high res tiff
tiff("Figures/dumbbell_high_res_ttest.tiff", units="in", width=9, height=7, res=300)
list.lowp %>% 
  ggplot(aes(x = lowerFC, xend = upperFC, y = FactorID, group = IDs), fontface = "bold") + 
  geom_dumbbell(aes(color = ifelse(FC < 0, "blue", "red"))) + 
  xlim(-1.5, 2.4) +
  geom_point(aes(FC, color = ifelse(FC < 0, "blue", "red"))) +
  geom_vline(xintercept = 1.9) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Fold Change from WT", y = "", title = "Eicosanoids",
       subtitle = "Molecules with P-value < 0.01 after T-test",
       caption = "Significant hits from metabolomic analysis of murine blood plasma.") + 
  theme(plot.title = element_text(hjust=0.5), 
        plot.subtitle = element_text(hjust=0.5),
        plot.caption = element_text(hjust = 0.5)) +
  coord_cartesian(xlim=c(-1.6, 2.2)) +
  scale_x_continuous(breaks=c(-1,0,1)) +
  geom_rect(aes(xmin=1.9, xmax=2.4, ymin=-Inf, ymax=Inf), fill="white") + 
  geom_text(aes(label=round(`P-value`,3), y=IDs, x=2.15), fontface = "bold",
            color = "black")
dev.off()



############old code####################
#for all p-vals

allpvals <- mapply(take.t, allmol_noNAs[,-1])
allpvals <- as.data.frame(allpvals)
allpvals$Alphanumeric <- rownames(allpvals)
allmol_pvals <- allpvals %>% merge(lookup, by = "Alphanumeric")
top20 <- allmol_pvals %>% select(Alphanumeric, IDs, allpvals) %>% arrange(allpvals) %>% head(20)
write.csv(top20, "R_data/top20.csv")

data.frame(p.lasma) %>% ggplot(aes(p.lasma)) + geom_histogram(binwidth = 0.04) +
  geom_vline(aes(xintercept = 0.01)) + xlab("P-values before correction")

###Table of top raw p-vals

p.lasma <- as.data.frame(p.lasma)
p.lasma$Alphanumeric <- rownames(p.lasma)
raw_pvals <- merge(lookup, p.lasma, by = "Alphanumeric")
raw_pvals <- raw_pvals %>% select(Alphanumeric, IDs, p.lasma) %>% arrange(p.lasma)
candidates_raw <- raw_pvals %>% filter(p.lasma < 0.01)
candidates_raw
write.csv(candidates_raw, "candidates_raw.csv")

###Create histogram of adjusted p-values

#create histogram of adjusted p-values
data.frame(adj_pvals) %>% ggplot(aes(adj_pvals)) + geom_histogram(binwidth = 0.033) + 
  geom_vline(aes(xintercept = 0.15)) + xlab("P-values after FDR correction")


###Create list of candidate molecules with p-values below cutoff

#create list of candidate molecules by merging with lookup file
adj_pvals <- as.data.frame(adj_pvals)

adj_pvals$Alphanumeric <- rownames(adj_pvals)
mol_pvals <- merge(lookup, adj_pvals, by = "Alphanumeric")
mol_pvals <- mol_pvals %>% select(IDs, adj_pvals) %>% arrange(adj_pvals)
candidates <- mol_pvals %>% filter(adj_pvals < 0.15)
candidates
write.csv(cand_df, "candidates.csv")


