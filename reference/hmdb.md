# The Human Metabolome DataBase multiplet table

The Human Metabolome DataBase multiplet table

## References

<https://hmdb.ca/>

## Examples

``` r
# Get all the 1-Methylhistidine peaks:
data("hmdb")
hmdb[hmdb$Metabolite == "1-Methylhistidine", ]
#> # A tibble: 4 × 12
#>   Metabolite      pos_in_ppm couple_code J_constant relative_intensity Accession
#>   <chr>                <dbl> <chr>       <chr>                   <dbl> <chr>    
#> 1 1-Methylhistid…       7.67 0           NA                          1 HMDB0000…
#> 2 1-Methylhistid…       7    0           NA                          1 HMDB0000…
#> 3 1-Methylhistid…       3.96 1,1         7.66,4.92                   2 HMDB0000…
#> 4 1-Methylhistid…       3.68 0           NA                          4 HMDB0000…
#> # ℹ 6 more variables: Blood <lgl>, Urine <lgl>, Solvent <chr>,
#> #   ChemShiftRef <chr>, NMRFreq <dbl>, pH <chr>
```
