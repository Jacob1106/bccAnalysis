# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

library(testthat)
library(rtracklayer)
library(GenomicRanges)

set.seed(123)

n_peaks <- 50
chrom <- "chr1"

# Function to randomly generate peak lengths between 100 and 500 bp
random_peak_lengths <- function(n) sample(100:500, n, replace = TRUE)

# Generate random start positions ensuring no overlaps within one set
generate_starts <- function(n, min_start = 1, max_start = 1e6) {
  starts <- numeric(0)
  lengths <- random_peak_lengths(n)
  i <- 1
  while(length(starts) < n) {
    candidate <- sample(min_start:max_start, 1)
    candidate_end <- candidate + lengths[i] - 1
    if (length(starts) == 0 || 
        !any(starts <= candidate_end & (starts + lengths[1:length(starts)] - 1) >= candidate)) {
      starts <- c(starts, candidate)
    }
    i <- i + 1
  }
  list(starts = starts, lengths = lengths[1:n])
}

# Generate shared peaks (first 25)
shared <- generate_starts(n_peaks / 2)

# Generate unique peaks for each object
unique1 <- generate_starts(n_peaks / 2, min_start = 1e6 + 1, max_start = 2e6)
unique2 <- generate_starts(n_peaks / 2, min_start = 2e6 + 1, max_start = 3e6)

# Combine for each GRanges object
starts1 <- c(shared$starts, unique1$starts)
lengths1 <- c(shared$lengths, unique1$lengths)
ends1 <- starts1 + lengths1 - 1

starts2 <- c(shared$starts, unique2$starts)
lengths2 <- c(shared$lengths, unique2$lengths)
ends2 <- starts2 + lengths2 - 1

# Generate random DNA sequences as character strings
generate_sequences <- function(lengths) {
  sapply(lengths, function(n) paste0(sample(c("A","C","G","T"), n, replace = TRUE), collapse = ""))
}

seqs1 <- generate_sequences(lengths1)
seqs2 <- generate_sequences(lengths2)

# P-values: 25 below 0.01, 25 above 0.01 (shuffled)
make_pvals <- function(n) {
  c(runif(n/2, 1e-6, 0.009), runif(n/2, 0.011, 0.05))
}
pvals1 <- sample(make_pvals(n_peaks))
pvals2 <- sample(make_pvals(n_peaks))

# Create GRanges objects
gr1 <- GRanges(
  seqnames = Rle(chrom),
  ranges = IRanges(start = starts1, end = ends1),
  strand = Rle("*"),
  pvalue = pvals1,
  sequence = seqs1
)

gr2 <- GRanges(
  seqnames = Rle(chrom),
  ranges = IRanges(start = starts2, end = ends2),
  strand = Rle("*"),
  pvalue = pvals2,
  sequence = seqs2
)

# Verify structure
cat("GR1: ", length(gr1), "ranges\n")
cat("Shared peaks (identical ranges):", sum(countOverlaps(gr1, gr2, type = "equal") > 0), "\n")
cat("GR1 pvalues < 0.01:", sum(mcols(gr1)$pvalue < 0.01), "\n")
cat("GR2 pvalues < 0.01:", sum(mcols(gr2)$pvalue < 0.01), "\n")



# Test tha filter_by_pvalue works: 

test_that("Filter by P value works", {
    #Create test list: 
    test_list <- list("sample1" = gr1, "sample2" = gr2)

    #Test 1: Test with known p-value: 
    filtered <- filter_by_pvalue(test_list, pvalue = 0.01, chromosomes = "chr1")

    expect_length(filtered, 2)
})
# test_list <- list("sample1" = gr1, "sample2" = gr2)
# test_list <- GRangesList(test_list)
# filtered <- filter_by_pvalue(test_list, pvalue = 0.01, chromosomes = "chr1")
