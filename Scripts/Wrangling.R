##packages needed for most analyses
lapply(c("tidyverse", "gdata", "ggfortify", "dplyr", "knitr", "dendextend", "spls", "glmnet"), require, character.only = TRUE)

#choose file
kplasma <- as.data.frame(idmol_for_R)
kplasma <- read.xls("~/Documents/R Projects/MullallyLab/Murine Plasma Metabolomics/idmol_for_R.xlsx", header = T)

#setting up dataset - make mouse IDs rownames and remove IDs column
allmol_df <- as.data.frame(allmol_noNAs)
allmol_df$VF <- rep(c("VF", "WT"), 5)
rownames(allmol_df) <- allmol_df$IDs
#remove ID column
allmol_df <- allmol_df[,2:ncol(allmol_df)]
#make all zeros NAs and remove incomplete cases
allmol_df[allmol_df == 0] <- NA
allmol_complete <- allmol_df[ , apply(allmol_df, 2, function(x) !any(is.na(x)))]
allmol_complete$VF <- rep(c(1,0), 5)
save(allmol_complete, file = "R_data/allmol_complete.rda")
save(allmol_df, file = "R_data/allmol_withNAs.rda")


idmol_NAs <- load("~/R_Files/Projects/Jak2VFMetabolomics/R_data/idmol_withNAs.rda")
#add rownames
rownames(kplasma) <- kplasma[,1]
#remove ID column -- kplasma to be used where genotype info is needed
kplasma <- kplasma[,-1]
#remove genotype info 
kplasma_pc <- kplasma[,-1]
#remove columns with NAs for above tests
#remove columns with NAs -- kplasma_NA to be used for PCA and t.tests
kplasma_pc[kplasma_pc == 0] = NA
kplasma_NA <- kplasma_pc %>% select_if(~ !any(is.na(.)))
idmol_withNAs <- kplasma
idmol_noNAforPCA <- kplasma_pc
allmol_noNAs <- allmol_complete
save(allmol_noNAs, file = "R_data/allmol_noNAs.rda")