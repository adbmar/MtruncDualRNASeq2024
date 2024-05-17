#' This script compute jaccard distance

renv::load()
library(tidyverse)
library(janitor)
library(vegan)

gpa <- read_delim(here::here("data/pangenome/gene_presence_absence.Rtab"))
gpa <- clean_names(gpa)

m <- gpa[,c("usda1021", "wsm1022")] %>%
    filter(!(usda1021 == 0 & wsm1022 == 0)) %>%
    as.matrix() %>%
    t()
ncol(m) # 7653 genes in union
sum(apply(m,2,sum) == 2) # 5496 shared genes
apply(m, 1, sum) # 6826 6323
1-vegdist(m, method = "jaccard") # 0.7181497

