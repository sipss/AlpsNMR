# The Human Metabolome DataBase multiplet table: urine metabolites normally found in NMR-based metabolomics

The Human Metabolome DataBase multiplet table: urine metabolites
normally found in NMR-based metabolomics

## References

<https://hmdb.ca/>

## Examples

``` r
data("HMDB_urine")
HMDB_urine[HMDB_urine$Metabolite == "1-Methyladenosine", ]
#>           Metabolite HMDB_code Shift_ppm Type J_Hz Height Urine_concentration
#> 10 1-Methyladenosine HMDB03331     3.855    s      1.0000               5.375
#> 11 1-Methyladenosine HMDB03331     3.945    t 3.33 0.0674               5.375
#> 12 1-Methyladenosine HMDB03331     8.475    s      0.4054               5.375
#> 13 1-Methyladenosine HMDB03331     3.925    t      0.1304               5.375
#> 14 1-Methyladenosine HMDB03331     4.775    t 5.23 0.1088               5.375
#> 15 1-Methyladenosine HMDB03331     3.885    t      0.1264               5.375
#> 16 1-Methyladenosine HMDB03331     6.105    d 5.34 0.1879               5.375
#> 17 1-Methyladenosine HMDB03331     6.025    d 6.05 0.0825               5.375
#> 18 1-Methyladenosine HMDB03331     4.475    t 4.73 0.1294               5.375
#> 19 1-Methyladenosine HMDB03331     4.295    m      0.0988               5.375
#> 20 1-Methyladenosine HMDB03331     8.415    s      0.3359               5.375
#> 21 1-Methyladenosine HMDB03331     8.155    s      0.0508               5.375
#> 22 1-Methyladenosine HMDB03331     8.265    s      0.1741               5.375
#>    n_reported_in_Urine Bouatra_2013
#> 10                   7          yes
#> 11                   7          yes
#> 12                   7          yes
#> 13                   7          yes
#> 14                   7          yes
#> 15                   7          yes
#> 16                   7          yes
#> 17                   7          yes
#> 18                   7          yes
#> 19                   7          yes
#> 20                   7          yes
#> 21                   7          yes
#> 22                   7          yes
```
