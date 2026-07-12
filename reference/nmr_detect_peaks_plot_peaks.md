# Plot multiple peaks from a peak list

Plot multiple peaks from a peak list

## Usage

``` r
nmr_detect_peaks_plot_peaks(
  nmr_dataset,
  peak_data,
  peak_ids,
  caption = paste("{peak_id}", "(NMRExp.\u00A0{NMRExperiment},",
    "\u03B3(ppb)\u00a0=\u00a0{gamma_ppb},", "\narea\u00a0=\u00a0{area},",
    "nrmse\u00a0=\u00a0{norm_rmse})")
)
```

## Arguments

- nmr_dataset:

  The `nmr_dataset_1D` object with the spectra

- peak_data:

  A data frame, the peak list

- peak_ids:

  The peak ids to plot

- caption:

  The caption for each subplot

## Value

A plot object
