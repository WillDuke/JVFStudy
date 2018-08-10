require(tidyverse)
load("R_data/allmol_noNAs.rda")
molwithsex <- allmol_noNAs %>% 
  mutate(Sex = c("M", "M", "F", "F", "M","M","M","M","F","F"))

malesWT <- molwithsex %>% filter(Sex == "M") %>% 
females <- molwithsex %>% filter(Sex == "F")

