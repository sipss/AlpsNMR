# AlpsNMR: Automated spectraL Processing System for NMR

AlpsNMR allows you to import NMR spectra into R and provides automated
and efficient signal processing for untargeted NMR metabolomics.

## Details

The following functions can be combined with the pipe. They create or
modify the
[nmr_dataset](https://sipss.github.io/AlpsNMR/reference/nmr_dataset.md)
object.

- [`nmr_read_samples_dir()`](https://sipss.github.io/AlpsNMR/reference/nmr_read_samples.md)
  or
  [`nmr_read_samples()`](https://sipss.github.io/AlpsNMR/reference/nmr_read_samples.md)

- [`nmr_interpolate_1D()`](https://sipss.github.io/AlpsNMR/reference/nmr_interpolate_1D.md)

- [`nmr_exclude_region()`](https://sipss.github.io/AlpsNMR/reference/nmr_exclude_region.md)

- [`nmr_normalize()`](https://sipss.github.io/AlpsNMR/reference/nmr_normalize.md)

- [plot()](https://sipss.github.io/AlpsNMR/reference/plot.nmr_dataset_1D.md)

There are also functions to extract the metadata and submit the samples
to irods, see the example below.

The
[nmr_dataset](https://sipss.github.io/AlpsNMR/reference/nmr_dataset.md)
object is essentially a list, so it is easy to access its components for
further analysis.

## See also

Useful links:

- <https://sipss.github.io/AlpsNMR/>

- <https://github.com/sipss/AlpsNMR>

- Report bugs at <https://github.com/sipss/AlpsNMR/issues>

## Author

**Maintainer**: Sergio Oller Moreno <sergioller@gmail.com>
([ORCID](https://orcid.org/0000-0002-8994-1549))

Authors:

- Sergio Oller Moreno <sergioller@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-8994-1549))

- Ivan Montoliu Roura <Ivan.MontoliuRoura@rd.nestle.com>

- Francisco Madrid Gambin <fmadrid@ibecbarcelona.eu>
  ([ORCID](https://orcid.org/0000-0001-9333-0014))

- Luis Fernandez <lfernandez@ibecbarcelona.eu>
  ([ORCID](https://orcid.org/0000-0001-9790-6287))

- H\<U+00E9\>ctor Gracia Cabrera <hgracia@ibecbarcelona.eu>

- Santiago Marco Col\<U+00E1\>s <smarco@ibecbarcelona.eu>
  ([ORCID](https://orcid.org/0000-0003-2663-2965))

Other contributors:

- Laura L\<U+00F3\>pez S\<U+00E1\>nchez \[contributor\]

- Nestl\<U+00E9\> Institute of Health Sciences \[copyright holder\]

- Institute for Bioengineering of Catalonia \[copyright holder\]

- Miller Jack <jack.miller@physics.org>
  ([ORCID](https://orcid.org/0000-0002-6258-1299)) (Autophase wrapper,
  ASICS export) \[contributor\]

## Examples

``` r
dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
my_nmr_dataset <- dataset %>%
    nmr_interpolate_1D(axis = c(0.4, 10)) %>%
    nmr_exclude_region(exclude = list(water = c(4.6, 5))) %>%
    nmr_normalize(method = "pqn") %>%
    plot()
#> Warning: There are not enough samples for reliably estimating the median spectra
#> ℹ The Probabalistic Quotient Normalization requires several samples to compute the median spectra. Your number of samples is low
#> ℹ Review your peaks before and after normalization to ensure there are no big distortions
```
