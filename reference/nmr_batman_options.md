# Batman Options helper

Batman Options helper

## Usage

``` r
nmr_batman_options(
  ppmRange = matrix(c(3, 3.1, 3.6, 3.7, 3.9, 4, 4, 4.1, 6.95, 7.05, 7.6, 7.7, 7.8, 7.9),
    ncol = 2, byrow = TRUE),
  specNo = "1",
  paraProc = 4L,
  negThresh = -0.5,
  scaleFac = 1e+06,
  downSamp = 1,
  hiresFlag = 1,
  randSeed = 100025L,
  nItBurnin = 200L,
  nItPostBurnin = 5000L,
  multFile = 2L,
  thinning = 50L,
  cfeFlag = 0,
  nItRerun = 5000L,
  startTemp = 1000,
  specFreq = 600,
  a = 1e-05,
  b = 1e-09,
  muMean = 1.1,
  muVar = 0.2,
  muVar_prop = 0.002,
  nuMVar = 0.0025,
  nuMVarProp = 0.1,
  tauMean = -0.05,
  tauPrec = 2,
  rdelta = 0.02,
  csFlag = 0
)
```

## Arguments

- ppmRange:

  Range of ppm to process

- specNo:

  Index of spectra to process

- paraProc:

  Number of cores to use

- negThresh:

  Truncation threshold for negative intensities

- scaleFac:

  Divide each spectrum by this number

- downSamp:

  Decimate each spectrum by this factor

- hiresFlag:

  Keep High Resolution deconvolved spectra

- randSeed:

  A random seed

- nItBurnin:

  Number of burn-in iterations

- nItPostBurnin:

  Number of iterations after burn-in

- multFile:

  Multiplet file (integer)

- thinning:

  Save MCMC state every thinning iterations

- cfeFlag:

  Same concentration for all spectra (fixed effect)

- nItRerun:

  Number of iterations for a batman rerun

- startTemp:

  Start temperature

- specFreq:

  NMR Spectrometer frequency

- a:

  Shape parameter for the gamma distribution (used for lambda, the
  precision)

- b:

  Rate distribution parameter for the gamma distribution (used for
  lambda, the precision)

- muMean:

  Peak width mean in ln(Hz)

- muVar:

  Peak width variance in ln(Hz)

- muVar_prop:

  Peak width proposed variance in ln(Hz)

- nuMVar:

  Peak width metabolite variance in ln(Hz)

- nuMVarProp:

  Peak width metabolite proposed variance in ln(Hz)

- tauMean:

  mean of the prior on tau

- tauPrec:

  inverse of variance of prior on tau

- rdelta:

  Truncation of the prior on peak shift (ppm)

- csFlag:

  Specify chemical shift for each multiplet in each spectrum?
  (chemShiftperSpectra.csv file)

## Value

A batman_options object with the Batman Options

## See also

Other batman functions:
[`nmr_batman`](https://sipss.github.io/AlpsNMR/reference/nmr_batman.md)

## Examples

``` r
bopts <- nmr_batman_options()
```
