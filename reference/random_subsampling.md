# Random subsampling

Random subsampling

## Usage

``` r
random_subsampling(
  sample_idx,
  iterations = 10L,
  test_size = 0.25,
  keep_together = NULL,
  balance_in_train = NULL
)
```

## Arguments

- sample_idx:

  Typically a numeric vector with sample index to be separated. A
  character vector with sample IDs could also be used

- iterations:

  An integer, the number of iterations in the random subsampling

- test_size:

  A number between 0 and 1. The samples to be included in the test set
  on each interation.

- keep_together:

  Either `NULL` or a factor with the same length as `sample_idx`.
  `keep_together` can be used to ensure that groups of samples are kept
  in together in all iterations (either on training or on test, but
  never split). A typical use case for this is when you have sample
  replicates and you want to keep all replicates together to prevent
  overoptimistic results (having one sample on the train subset and its
  replicate on the test subset would make the prediction easier to
  guess). Another use case for this is when you have a longitudinal
  study and you want to keep some subjects in the same train or test
  group, because you want to use some information in a longitudinal way
  (e.g. a multilevel plsda model).

- balance_in_train:

  Either `NULL` or a factor with the same length as `sample_idx`.
  `balance_in_train` can be used to force that on each iteration, the
  train partition contains the same number of samples of the given
  factor levels. For instance, if we have a dataset with 40 samples of
  class "A" and 20 samples of class "B", using a `test_size = 0.25`, we
  can force to always have 16 samples of class "A" and 16 samples of
  class "B" in the training subset. This is beneficial to those
  algorithms that require that the training groups are balanced.

## Value

A list of length equal to `iterations`. Each element of the list is a
list with two entries (`training` and `test`) containing the
`sample_idx` values that will belong to each subset.

## Examples

``` r
random_subsampling(1:100, iterations = 4, test_size = 0.25)
#> [[1]]
#> [[1]]$training
#>  [1]   1   2   3   5   6   7   8   9  10  12  15  16  18  19  20  23  24  26  27
#> [20]  28  30  32  33  36  37  38  39  40  41  42  43  44  45  47  48  49  51  52
#> [39]  53  54  56  57  58  60  61  63  64  65  66  67  69  70  71  72  73  74  76
#> [58]  77  78  79  80  81  82  83  85  87  88  89  92  94  95  97  98  99 100
#> 
#> [[1]]$test
#>  [1]  4 11 13 14 17 21 22 25 29 31 34 35 46 50 55 59 62 68 75 84 86 90 91 93 96
#> 
#> 
#> [[2]]
#> [[2]]$training
#>  [1]  1  4  5  6  7  9 10 11 12 14 15 16 17 18 19 20 21 22 24 25 26 27 28 29 31
#> [26] 32 34 36 37 38 39 40 41 43 45 47 48 50 52 54 55 56 57 58 61 62 63 66 67 68
#> [51] 69 70 71 72 73 76 77 78 79 80 81 82 84 85 87 88 89 91 92 93 94 95 97 98 99
#> 
#> [[2]]$test
#>  [1]   2   3   8  13  23  30  33  35  42  44  46  49  51  53  59  60  64  65  74
#> [20]  75  83  86  90  96 100
#> 
#> 
#> [[3]]
#> [[3]]$training
#>  [1]   1   2   3   4   5   6   7   9  10  12  13  15  16  18  19  20  21  22  23
#> [20]  25  26  27  28  30  31  32  33  35  36  37  38  41  43  44  45  46  47  48
#> [39]  51  53  54  55  56  58  59  63  64  65  66  68  71  73  74  76  77  78  80
#> [58]  81  82  83  84  85  86  87  88  90  91  92  93  94  96  97  98  99 100
#> 
#> [[3]]$test
#>  [1]  8 11 14 17 24 29 34 39 40 42 49 50 52 57 60 61 62 67 69 70 72 75 79 89 95
#> 
#> 
#> [[4]]
#> [[4]]$training
#>  [1]   1   2   3   4   5   6   8  10  11  13  14  15  16  17  18  19  20  21  22
#> [20]  24  25  26  28  30  31  32  34  35  36  37  39  40  41  42  44  46  47  49
#> [39]  50  52  54  57  58  59  62  63  64  65  66  67  70  71  72  73  74  75  76
#> [58]  77  78  79  80  82  83  85  88  91  92  93  94  95  96  97  98  99 100
#> 
#> [[4]]$test
#>  [1]  7  9 12 23 27 29 33 38 43 45 48 51 53 55 56 60 61 68 69 81 84 86 87 89 90
#> 
#> 

subject_id <- c("Alice", "Bob", "Charlie", "Eve")
random_subsampling(1:4, iterations = 2, test_size = 0.25, keep_together = subject_id)
#> [[1]]
#> [[1]]$training
#> [1] 1 3 4
#> 
#> [[1]]$test
#> [1] 2
#> 
#> 
#> [[2]]
#> [[2]]$training
#> [1] 1 2 3
#> 
#> [[2]]$test
#> [1] 4
#> 
#> 
```
