# Filtering GRange Objects by number of overlaps with other samples. This function is used to merge multiple GRange Objects into one.

Filtering GRange Objects by number of overlaps with other samples. This
function is used to merge multiple GRange Objects into one.

## Usage

``` r
overlapGRangeList(grange_list, num_of_overlaps_required = 1)
```

## Arguments

- grange_list:

  List of GRanges Objects that are similar or part of the same condition

- num_of_overlaps_required:

  Integer value that determines how many overlaps a peak must have with
  the other GRanges objects in order to be outputed by the function

## Value

Returns a single GRange Object that is merged and filtered by overlaps.

## Author

Jacob Martin
