# Importing narrowPeak files and converting into a GRanges Object

Importing narrowPeak files and converting into a GRanges Object

## Usage

``` r
import_samples(named_file_list)
```

## Arguments

- named_file_list:

  A named list with file paths for narrowPeak files. These files will be
  imported as GRanges Objects.

## Value

Returns a list of GRange Objects with a new column for normal p values.
The list is sorted and edited based on the p value and chromosome
inputs.

## Author

Jacob Martin
