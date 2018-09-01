require(ggalt)
require(dplyr)
#load dataset
load("~/R_Files/Projects/Jak2VFMetabolomics/R_data/allmol_noNAs.rda")
load("~/R_Files/Projects/Jak2VFMetabolomics/R_data/lookup.rda")
#subset to eicosanoids
eico <- allmol_noNAs %>% select_if(grepl("A", colnames(allmol_noNAs))) 
#construct function to take wilcox tests and calculate means
take.mann <- function(x){
  y <- rep(c(0,1), 5)
  Median_VF <- median(x[seq(1,9,2)])
  Median_WT <- median(x[c(2,4,6,8,10)])
  mw <- wilcox.test(x~y, conf.int = TRUE)
  c(Median_WT/10^5, Median_VF/10^5, mw[["estimate"]][["difference in location"]]/10^5, 
    mw[["statistic"]][["W"]], mw[["conf.int"]]/10^5, 
    mw[["p.value"]])
}

#apply to all columns
mw.ps <- mapply(take.mann, eico)

#transpose to create vertical list
list.mw <- as.data.frame(t(mw.ps))
colnames(list.mw) <- c("Median_WT", "Median_VF", "Difference", "W", 
                       "CI - Lower Bound", "CI - Upper Bound", "P-value")
#add lookup key
list.mw$Alphanumeric <- rownames(list.mw)
list.full <- merge(list.mw, lookup, by = "Alphanumeric")
rownames(list.mw) <- NULL

#correct column order and arrange by p-value, and plot (skip if making confint)
list.full %>% 
  unite(mzid, mz, RT, sep = "_") %>%
  dplyr::select(IDs, Difference, W, `CI - Lower Bound`,
                `CI - Upper Bound`, `P-value`, mzid) %>%
  arrange(`P-value`)
  ggplot(aes(`P-value`)) + 
  scale_x_continuous(trans = "log10") +
  geom_histogram(bins = 30) +
  geom_vline(xintercept = 0.05)

#write the csv
write.csv(list.full, file = "Figures/MWEicosanoidsNoMeans.csv")

listfordb <- list.full %>% 
  rename(lower = `CI - Lower Bound`,
         upper = `CI - Upper Bound`) %>%
  mutate(FC = (Difference)/Median_WT, 
         lowerFC = lower/Median_WT, 
         upperFC = upper/Median_WT) %>%
  filter(`P-value` < 0.05) %>%
  arrange(`P-value`)
listfordb$FactorID <- factor(listfordb$IDs, levels=as.character(listfordb$IDs))
listfordb$FactorID <- factor(listfordb$FactorID, levels=rev(levels(listfordb$FactorID)))

tiff("Figures/dumbbell_high_res_mw.tiff", units="in", width=9, height=7, res=300)
listfordb %>% 
  ggplot(aes(x = lowerFC, xend = upperFC, y = FactorID, group = IDs), fontface = "bold") + 
  geom_dumbbell(aes(color = ifelse(FC < 0, "blue", "red"))) + 
  geom_point(aes(FC, color = ifelse(FC < 0, "blue", "red"))) +
  geom_vline(xintercept = 11.1) +
  theme_bw() +
  theme(legend.position = "none") + 
  labs(x = "Fold Change of Median from WT", y = "", title = "Eicosanoids",
       subtitle = "Molecules with P-value < 0.05 after Wilcox Test",
       caption = "Significant hits from metabolomic analysis of murine blood plasma.") + 
  theme(plot.title = element_text(hjust=0.5), 
        plot.subtitle = element_text(hjust=0.5),
        plot.caption = element_text(hjust = 0.5)) +
  coord_cartesian(xlim=c(-1, 11.8)) +
  geom_rect(aes(xmin=11.1, xmax=12.5, ymin=-Inf, ymax=Inf), fill="white") +
  geom_text(aes(label=round(`P-value`,3), y=IDs, x=11.75), fontface = "bold",
            color = "black")
dev.off()
listfordb %>% 
  ggplot(aes(x = lower, xend = upper, y = FactorID, group = IDs)) + 
  geom_dumbbell()
