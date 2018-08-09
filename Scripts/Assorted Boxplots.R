require(tidyverse)
load("R_data/allmol_noNAs.rda")
allmol_noNAs$VF <- rep(c("VF", "WT"), 5)
allmol_noNAs %>% ggplot(aes(VF, B004/10^8)) + geom_boxplot() +
  geom_point() + ylim(1,4.8) + xlab("Genotype") + ylab("Intensity") +
  ggtitle("Docosahexaenoic Acid")
allmol_noNAs %>% ggplot(aes(VF, B002/10^8)) + geom_boxplot() + 
  geom_point() + ylim(1,2.5) + xlab("Genotype") + ylab("Intensity") +
  ggtitle("Arachidonic Acid")