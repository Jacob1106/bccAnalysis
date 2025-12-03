# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html
library(bccAnalysis)
library(testthat)
library(rtracklayer)
library(GenomicRanges)

#TEST OBJECT GENERATION 
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

# Generate raw p-values: 25 below 0.01 and 25 above 0.01 (shuffled)
make_pvals <- function(n) {
  c(runif(n/2, 1e-6, 0.009), runif(n/2, 0.011, 0.05))
}

raw_pvals1 <- sample(make_pvals(n_peaks))
raw_pvals2 <- sample(make_pvals(n_peaks))

# Compute negative log10 p-values
neg_log_pvals1 <- -log10(raw_pvals1)
neg_log_pvals2 <- -log10(raw_pvals2)

# Create GRanges objects with pvalue as negative log10 of raw p-value
gr1 <- GRanges(
  seqnames = Rle(chrom),
  ranges = IRanges(start = starts1, end = ends1),
  strand = Rle("*"),
  pvalue = neg_log_pvals1,
  sequence = seqs1
)

gr2 <- GRanges(
  seqnames = Rle(chrom),
  ranges = IRanges(start = starts2, end = ends2),
  strand = Rle("*"),
  pvalue = neg_log_pvals2,
  sequence = seqs2
)

# Check example output
gr1[1:3]
gr2[1:3]
gr1$pValue <- gr1$pvalue
gr2$pValue <- gr2$pvalue
gr1$pvalue <- NULL
gr2$pvalue <- NULL

cat("GR1: ", length(gr1), "ranges\n")
cat("Shared peaks (identical ranges):", sum(countOverlaps(gr1, gr2, type = "equal") > 0), "\n")
cat("GR1 pvalues < 0.01:", sum(mcols(gr1)$pValue < 2), "\n")
cat("GR2 pvalues < 0.01:", sum(mcols(gr2)$pValue < 2), "\n")

#Generation of overlaps and uniques for the example data sets: 

index <- findOverlaps(gr1, gr2)
overlapsGr1 <- gr1[queryHits(index)]
overlapsGr2 <- gr2[subjectHits(index)]
uniqueGr1 <- gr1[-queryHits(index)]
uniqueGr2 <- gr2[-subjectHits(index)]


#--------------------------------------------------------------------
# Test tha filter_by_pvalue works: 

#Test that empty expect error 
#Check it returns correct length with valid p value 
#Test that it returns an error if its a genomic ranges object and not a list
#Test to see if chromosome filter is working

test_that("Filter by P value works", {
  #Create test list: 
  test_list <- GRangesList(list("sample1" = gr1))


  #Test 1: Test with known p-value: 
  filtered <- filter_by_pvalue(test_list, pvalue = 0.01, chromosomes = "chr1")
  expect_equal(length(filtered), 1)
  expect_s4_class(filtered[[1]], "GenomicRanges")
  expect_type(filtered, "list")
  expect_equal(length(filtered[[1]]), 25)

  #Test 2: Test Invalid p value: 
  expect_error(filter_by_pvalue(GRangeList(), pvalue = 0.01, chromosomes = "chr1"))

  #Test 3: Test if chromosome filter works: 
  expect_equal(length(filter_by_pvalue(test_list, pvalue = 0.01, chromosomes = "chr2")[[1]]), 0)
  expect_equal(length(filter_by_pvalue(test_list, pvalue = 0.9, chromosomes = "chr1")[[1]]), 50)
  expect_equal(length(filter_by_pvalue(test_list, pvalue = 0.01, chromosomes = "chr1")[[1]]), 25)

  #Test 4: for genomic range input: 
  expect_error(filter_by_pvalue(gr1, pvalue = 0.01, chromosomes = "chr1"))

})

#OverlapGRangeList TESTS: 

#Test If it ouputs the correct options 
#Test if input is wrong it returns an error 
#Test If the there is an invalid overlap count it returns an error

test_that("OverlapGRangeList Test",{
  #GRANGE List Example: 
  ex <- GRangesList(list("sample1" = gr1, "sample2" = gr2))
  fil <- overlapGRangeList(ex, num_of_overlaps_required = 1)
  fil_noOverlap <- overlapGRangeList(ex, num_of_overlaps_required = 0)
  #Test 1: Outputs the correct class: 
  expect_s4_class(fil, "GenomicRanges")
  expect_s4_class(fil_noOverlap, "GenomicRanges")

  #Test 2: Test if it returns error for wrong input e.g. genomic ranges input: 
  expect_error(overlapGRangeList(gr1, num_of_overlaps_required = 1))

  #Test 3: Test length of outputs for different overlap conditions: 
  expect_equal(length(fil), 25)
  expect_equal(length(fil_noOverlap), 75)

  #Test 4: Invalid overlap number given: 
  expect_error(overlapGRangeList(fil, num_of_overlaps_required = -1)) # Using -1 as invalid input: THERE CANT BE NEGATIVE OVERLAPS 
})


#proportionOverlapTibble Test: 

#Test invalid input
#Test all outputs for correct inputs
#Test that the tibble is generated 

test_that("Tibble Generation Test", {
  
  #Test Examples: 
  ex <- GRangesList(list("sample1" = gr1, "sample2" = gr2))
  f <- proportionOverlapTibble(object1 = gr1, object2 = gr2, object1_unique = uniqueGr1, object2_unique = uniqueGr2, overlapobject1 = overlapsGr1, overlapobject2 = overlapsGr2)

  #Test 1: invalid input
  expect_error(proportionOverlapTibble(object1 = ex, object2 = gr2, object1_unique = uniqueGr1, object2_unique = uniqueGr2, overlapobject1 = overlapsGr1, overlapobject2 = overlapsGr2)) 

  #Test 2: Outputs: 
  expect_s3_class(f, "data.frame")
  expect_all_equal(f$PercentageOverlapPeaks, 50)
  expect_equal(f$sample[1], "gr1")
  expect_equal(f$sample[2], "gr2")
  expect_all_equal(f$TotalPeaks, 50)
  expect_all_equal(f$NumberofOverlapPeaks, 25)
})

#Venn Diagram Test: 
# Test If a plot is generated 
#Test If the input is wrong there is an error: 

test_that("Venn Diagram Test ", {
  example <- custom_venn_diagram("GR1", "GR2", DataA = gr1, DataB = gr2)

  #Test 1: error if wrong input: 
  expect_error(custom_venn_diagram("GR1", "GR2", dataA = GRangesList(list("sample1" = gr1)), dataB = gr2))

})
# ex <- GRangesList(list("sample1" = gr1, "sample2" = gr2))



test_that("WidthPlot", {
  x <- compareGRangeWidth(gr1, gr2)

  #Test if a ggplot is produced: 
  expect_s3_class(x, "ggplot")
})