# ROIs for cell samples

The template ROI_cell contains the targeted list of metabolites to be
quantified (cell samples)

## References

[github.com/danielcanueto/rDolphin](https://sipss.github.io/AlpsNMR/reference/github.com/danielcanueto/rDolphin)

## Examples

``` r
data("ROI_cell")
ROI_cell[ROI_cell$Metabolite == "Valine", ]
#>   ROI.left.edge..ppm. ROI.right.edge..ppm. Quantification.Mode Metabolite
#> 3                1.02                0.973    Baseline Fitting     Valine
#> 5                1.06                1.015    Baseline Fitting     Valine
#>   Quantification.Signal Chemical.shift.ppm. Chemical.shift.tolerance..ppm.
#> 3                     2              0.9850                           0.02
#> 5                     1              1.0356                           0.02
#>   Half.bandwidth..Hz. Multiplicity J.coupling..Hz. Roof.effect HMDB_code
#> 3                 1.4            2            7.01           0      <NA>
#> 5                 1.4            2            7.05           0      <NA>
```
