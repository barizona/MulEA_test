#xxxxxxxxxxxx
# Initialization ----
#xxxxxxxxxxxx
# MulEA
# devtools::install_github("koralgooll/MulEA")
library(MulEA)

# CPP
library(Rcpp)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp("source/set-based-enrichment-test.cpp")
source("source/set-based-enrichment-test.R")

# GergoR
source("source/Gergo_Enrichment_test_G.r") # required packages: snow, rlecuyer

# Other
library(tidyverse)
#xxxxxxxxxxxx

# AIM: Testing if FDR (empirically corrected p-value) could be smaller than the original p-value

#xxxxxxxx
# Read input files ----
#xxxxxxxx

# GMT tab for MulEA
gmt_tab <- MulEA::readGmtFileAsDataFrame("input/KEGG_filtered.gmt")
# GMT list for CPP and GergoR
# GMT
x <- strsplit(readLines("input/KEGG_filtered.gmt"), "\t")
gmt_list <- lapply(x, tail, n=-2)
names(gmt_list) <- lapply(x, head, n=1)
rm(x)

# gene set
gene_set <- readLines("input/input_select.txt")

# background
gene_bg <- unlist(gmt_list) %>% unique()

#xxxxxxxxxx
# MulEA ORA ----
#xxxxxxxxxx
mulea_ora_results <- MulEA::ORA(
  gmt = gmt_tab, 
  testData = gene_set, 
  pool = gene_bg,
  adjustMethod = "PT",
  numberOfPermutations = 1000) %>% 
  MulEA::runTest() %>% 
  arrange(pValue)

#xxxxxxxxxx
# CPP ORA ----
#xxxxxxxxxx
cpp_results <- set.based.enrichment.test(steps = 1000, 
                                         pool = gene_bg, 
                                         select = gene_set, 
                                         DB = gmt_list) %>% 
  arrange(P)

#xxxxxxxxxx
# GergoR ORA ----
#xxxxxxxxxx
GergoR_results <- HyperGeomFDR(steps = 1000, 
                               pool = gene_bg, 
                               select = gene_set, 
                               DB = gmt_list, 
                               nthreads = 10) %>% 
  arrange(P)

#xxxxxxxxxx
# In which cases pValue or P > adjustedPValueEmpirical or FDR ----
#xxxxxxxxxx
P_greater_result_list <- list()
P_greater_result_list$mulea <- mulea_ora_results %>% 
  filter(pValue > adjustedPValueEmpirical)
P_greater_result_list$cpp <- cpp_results %>% 
  filter(P > FDR)
P_greater_result_list$GergoR <- GergoR_results %>% 
  filter(P > FDR)
  
