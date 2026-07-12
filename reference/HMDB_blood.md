# The Human Metabolome DataBase multiplet table: blood metabolites normally found in NMR-based metabolomics

The Human Metabolome DataBase multiplet table: blood metabolites
normally found in NMR-based metabolomics

## References

<https://hmdb.ca/>

## Examples

``` r
data("HMDB_blood")
HMDB_blood[HMDB_blood$Metabolite == "1-Methylhistidine", ]
#>          Metabolite HMDB_code Shift_ppm Type      J_Hz Height
#> 1 1-Methylhistidine HMDB00001     3.695    s           1.0000
#> 2 1-Methylhistidine HMDB00001     3.085   dd           0.1135
#> 3 1-Methylhistidine HMDB00001     7.685    s           0.2230
#> 4 1-Methylhistidine HMDB00001     3.175   dd           0.1115
#> 5 1-Methylhistidine HMDB00001     7.015    s           0.2369
#> 6 1-Methylhistidine HMDB00001     3.975   dd 7.66 4.92 0.1092
#>   Blood_concentration n_reported_in_Blood
#> 1            12.52857                   7
#> 2            12.52857                   7
#> 3            12.52857                   7
#> 4            12.52857                   7
#> 5            12.52857                   7
#> 6            12.52857                   7
```
