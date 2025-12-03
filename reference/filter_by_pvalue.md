# Filtering Function that filter by a custom p value and chromosome list

Filtering Function that filter by a custom p value and chromosome list

## Usage

``` r
filter_by_pvalue(
  GRangeList,
  pvalue,
  chromosomes = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
    "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
    "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
)
```

## Arguments

- pvalue:

  float value p value

- chromosomes:

  list of characters e.g. c("chr1", "chr2")

- list:

  List of GRanges Objects

## Value

Returns a list of GRange Objects with a new column for normal p values.
The list is sorted and edited based on the p value and chromosome
inputs.

## Author

Jacob Martin
