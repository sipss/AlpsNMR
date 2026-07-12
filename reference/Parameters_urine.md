# Parameters for urine samples profiling

The template `Parameters_urine` contains the chosen normalization
approach (by default, PQN), the Spectometer Frequency (by default,
600.04MHz), alignment (by default, TSP 0.00 ppm), bucket resolution (by
default, 0.00023)

## References

[github.com/danielcanueto/rDolphin](https://sipss.github.io/AlpsNMR/reference/github.com/danielcanueto/rDolphin)

## Examples

``` r
data("Parameters_urine")
Parameters_urine
#>                                                                   Parameter
#> 1                                                           nmr folder path
#> 2                                                             1D data index
#> 3                                                                   proc_no
#> 4                                                      spectra dataset path
#> 5                                                Metadata path (csv format)
#> 6                                                         ROI patterns file
#> 7  Normalization (0=No;1=Eretic; 2=TSP; 3=Creatinine; 4=Spectra Sum; 5=PQN)
#> 8                              Alignment (0=No;1=Glucose; 2=TSP; 3=Formate)
#> 9                                                      Default Suppressions
#> 10                                              Spectometer Frequency (MHz)
#> 11                                                        Bucket resolution
#> 12                                                                 Biofluid
#> 13                                                                  2D-Path
#> 14                                              Specific dataset parameters
#>                Value
#> 1                   
#> 2                   
#> 3                   
#> 4    NMR_spectra.csv
#> 5  meta_rDolphin.csv
#> 6            ROI.csv
#> 7                  5
#> 8                  2
#> 9                   
#> 10            600.04
#> 11           0.00023
#> 12             Urine
#> 13                  
#> 14                  
```
