# ROIs for urine samples

The template ROI_urine contains the targeted list of metabolites to be
quantified (urine samples)

## References

[github.com/danielcanueto/rDolphin](https://sipss.github.io/AlpsNMR/reference/github.com/danielcanueto/rDolphin)

## Examples

``` r
data("ROI_urine")
ROI_urine[ROI_urine$Metabolite == "Valine", ]
#>   ROI.left.edge ROI.right.edge Quantification.Mode Metabolite
#> 4          1.09           1.03    Baseline Fitting     Valine
#>   Quantification.Signal Chemical.shift Chemical.shift.tolerance Half.bandwidth
#> 4                     1          1.047                    0.005            1.2
#>   Multiplicity J.coupling Roof.effect HMDB.code
#> 4            2          7           0 HMDB00883
```
