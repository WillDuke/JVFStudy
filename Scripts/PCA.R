#require packages
lapply(c("tidyverse", "gdata", 
         "ggfortify", "dplyr", 
         "knitr", "dendextend", 
         "spls", "glmnet"), require, character.only = TRUE)
#load data
load("R_data/allmol_noNAs.rda")
#create subset
eico <- allmol_noNAs %>% 
      dplyr::select(Genotype, grep("A", colnames(allmol_noNAs))) %>% 
      mutate(Sex = c("M", "M", "F", "F", "M", "M", "M", "M", "F", "F")) %>%
      mutate(Sex = as.factor(Sex))

#adjust names for plot output
allmol_noNAs <- allmol_noNAs %>% 
        mutate(VF = rep(c("VF", "WT"), 5)) %>%
        mutate(VF = as.factor(VF)) %>% 
        rename(Genotype = VF) %>%
        mutate(Sex = c("M", "M", "F", "F", "M", "M", "M", "M", "F", "F")) %>%
        mutate(Sex = as.factor(Sex))

#plot PCA on all molecules
autoplot(prcomp(allmol_noNAs[,c(-1, -2620)]), 
         data = allmol_noNAs, 
         size = 2.3, 
         colour = 'Genotype', 
         shape = 'Sex')
ggsave("Figures/PCA_of_allmols.png")

#plot PCA on eicosanoids
autoplot(prcomp(eico[,c(-1,-114)]), 
         data = eico, 
         size = 2.3, 
         colour = 'Genotype',
         shape = 'Sex')
ggsave("Figures/PCA_of_eicosanoids.png")
