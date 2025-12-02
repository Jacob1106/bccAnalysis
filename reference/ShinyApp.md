# Shiny App for Chip-Seq Analysis

Shiny App for Chip-Seq Analysis

## Usage

``` r
ShinyApp(named_file_list, height = 850)
```

## Arguments

- named_file_list:

  A named list with file paths for narrowPeak files. These files will be
  imported as GRanges Objects.

- height:

  Controls Browser window size on first launch

## Value

Returns an Analysis App. Outputs: Overlaps, Unique peak annotations,
width analysis and proportion overlap table

## Author

Jacob Martin
