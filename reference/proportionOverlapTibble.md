# Filtering Function that filter by a custom p value and chromosome list

Filtering Function that filter by a custom p value and chromosome list

## Usage

``` r
proportionOverlapTibble(
  object1,
  object2,
  overlapobject1,
  overlapobject2,
  object1_unique,
  object2_unique
)
```

## Arguments

- object1:

  GRanges Object

- object2:

  GRanges Object

- overlapobject1:

  The overlap G ranges object of object1

- overlapobject2:

  The overlap G ranges object of object2

- object1_unique:

  The unique peaks as a G ranges object of object1

- object2_unique:

  The unique peaks as a G ranges object of object2

## Value

Returns a tibble that shows: Proportion of Overlap, TotalPeak coverage,
Total number of unique bases, number of overlap peaks and percentage
overlap of peaks for object 1 and 2.

## Author

Jacob Martin
