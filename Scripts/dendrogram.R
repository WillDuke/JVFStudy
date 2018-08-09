lapply(c("tidyverse", "gdata", 
         "ggfortify", "dplyr", 
         "knitr", "dendextend", 
         "spls", "glmnet"), 
       require, character.only = TRUE)

load("R_data/allmol_noNAs.rda")
plasma_dend <- allmol_noNAs[,-1] %>% dist %>% 
    hclust(method = "complete") %>% 
    as.dendrogram()
plasma_dend %>% set("leaves_pch", c(17)) %>%  # node point type
  set("leaves_cex", 2) %>%  # node point size
  set("branches_k_color", value = c("blue", "red"), k = 2) %>%
  set("leaves_col", c(rep("blue",5), rep("red", 5))) %>% #node point color
  plot()
