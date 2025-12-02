# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html


library(rtracklayer)
library(GenomicRanges)

set.seed(42)

create_test_granges <- function() {
  # [Keep ALL your existing position/pval data EXACTLY the same]
  shared_pos <- list(c(26226,26613), c(27825,28026), c(31245,31711), c(33327,33759), c(91162,91621),
                    c(98247,98626), c(107474,107788), c(116740,116952), c(146317,146646), c(229259,229660),
                    c(234054,234296), c(243963,244066), c(256788,257276), c(288390,288571), c(442418,442875),
                    c(529904,530220), c(571859,572133), c(619177,619419), c(631263,631442), c(670488,670698),
                    c(709571,710061), c(772247,772519), c(776647,776799), c(777573,777720), c(935519,935813))
  
  shared_pvals <- c(2.85e-04,4.25e-03,3.13e-02,3.27e-02,8.11e-07,7.63e-03,7.85e-10,8.01e-08,
                   4.91e-02,4.26e-02,4.79e-02,1.50e-02,1.96e-02,1.32e-02,2.23e-03,2.40e-03,
                   6.19e-06,1.66e-02,1.15e-06,3.64e-02,3.90e-09,1.96e-09,3.93e-02,1.96e-02,9.80e-08)
  shared_pValue <- -log10(shared_pvals)
  
  unique1_pvals <- c(1.57e-03,2.70e-03,5.19e-03,1.99e-02,1.27e-02,3.35e-02,1.24e-02,1.56e-07,
                    1.08e-02,1.26e-02,3.30e-05,1.66e-06,2.36e-05,4.09e-07,2.11e-02,7.12e-06,
                    3.83e-02,9.71e-08,5.05e-05,1.10e-04,3.22e-02,1.13e-02,1.03e-02,3.06e-02,2.32e-02)
  unique1_pValue <- -log10(unique1_pvals)
  
  unique2_pvals <- c(7.64e-07,1.22e-02,3.91e-02,1.24e-02,3.58e-02,3.61e-02,1.94e-07,4.20e-02,
                    4.07e-06,2.85e-02,5.68e-10,1.56e-02,1.87e-05,6.83e-08,6.08e-07,3.15e-02,
                    2.06e-06,1.77e-02,1.15e-03,4.46e-02,1.81e-07,1.73e-02,1.44e-02,6.12e-05,8.79e-09)
  unique2_pValue <- -log10(unique2_pvals)
  
  # BULLETPROOF METHOD: Create with DataFrame, force column names
  shared_mcols <- DataFrame(pValue = shared_pValue)
  shared_gr <- GRanges("chr1", 
                      IRanges(unlist(lapply(shared_pos, `[`, 1)), unlist(lapply(shared_pos, `[`, 2))), 
                      mcols = shared_mcols)
  
  # Verify column name BEFORE combining (debugging)
  stopifnot("pValue" %in% colnames(mcols(shared_gr)))
  
  unique1_mcols <- DataFrame(pValue = unique1_pValue)
  gr1_unique <- GRanges("chr1", 
                       IRanges(unlist(lapply(unique1_pos, `[`, 1)), unlist(lapply(unique1_pos, `[`, 2))), 
                       mcols = unique1_mcols)
  
  unique2_mcols <- DataFrame(pValue = unique2_pValue)
  gr2_unique <- GRanges("chr1", 
                       IRanges(unlist(lapply(unique2_pos, `[`, 1)), unlist(lapply(unique2_pos, `[`, 2))), 
                       mcols = unique2_mcols)
  
  gr1 <- c(shared_gr, gr1_unique)
  gr2 <- c(shared_gr, gr2_unique)
  gr1$pValue <- gr1$pvalue
  gr2$pValue <- gr2$pvalue
  gr1$pvalue <- NULL
  gr2$pvalue <- NULL
  # FINAL VERIFICATION
  stopifnot(all(colnames(mcols(gr1)) == "pValue"))
  stopifnot(all(colnames(mcols(gr2)) == "pValue"))
  
  list(gr1 = gr1, gr2 = gr2)
}








# Test tha filter_by_pvalue works: 

# test_that("Filter by P value works", {
#     #Create test list: 
#     test_list <- list("sample1" = gr1, "sample2" = gr2)

#     #Test 1: Test with known p-value: 
#     filtered <- filter_by_pvalue(test_list, pvalue = 0.01, chromosomes = "chr1")

#     expect_length(filtered, 2)
#     expect_true()

# })
# test_list <- list("sample1" = gr1, "sample2" = gr2)
# filtered <- filter_by_pvalue(test_list, pvalue = 0.01, chromosomes = "chr1")

f <- rtracklayer::import("/Volumes/JM/Analysis_Projects/ChIP.v.CnR/CTCF/NarrowPeakFiles/CHIP/SJHGG059115_C5-LTC42_P21-CTCF-AB2_REP1_peaks.narrowPeak")