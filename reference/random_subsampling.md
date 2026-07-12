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
#>  [1]   2   3   5   6   7  11  12  14  15  17  18  19  20  21  22  23  25  26  27
#> [20]  28  29  30  31  35  36  38  40  41  42  43  45  46  47  48  50  52  53  54
#> [39]  55  57  58  59  60  61  62  63  64  65  66  67  68  69  70  73  74  75  77
#> [58]  79  80  81  82  83  86  87  89  90  91  93  94  95  96  97  98  99 100
#> 
#> [[1]]$test
#>  [1]  1  4  8  9 10 13 16 24 32 33 34 37 39 44 49 51 56 71 72 76 78 84 85 88 92
#> 
#> 
#> [[2]]
#> [[2]]$training
#>  [1]   1   2   3   4   5   6   8  11  12  13  14  16  18  19  20  21  22  23  24
#> [20]  25  27  28  29  30  31  32  33  34  35  36  37  39  40  41  42  43  44  45
#> [39]  46  47  48  52  53  54  55  56  57  59  60  61  62  64  67  70  72  73  74
#> [58]  75  78  79  80  81  83  84  85  87  89  91  92  93  95  96  97  98 100
#> 
#> [[2]]$test
#>  [1]  7  9 10 15 17 26 38 49 50 51 58 63 65 66 68 69 71 76 77 82 86 88 90 94 99
#> 
#> 
#> [[3]]
#> [[3]]$training
#>  [1]   1   2   3   4   5   7   9  10  11  12  13  14  15  18  19  20  21  22  23
#> [20]  24  26  27  28  29  32  33  34  35  36  38  41  42  43  45  46  51  52  54
#> [39]  55  56  58  59  60  61  62  63  64  65  66  68  70  71  73  75  76  77  78
#> [58]  79  80  81  82  83  84  85  87  89  90  92  93  94  95  96  97  98 100
#> 
#> [[3]]$test
#>  [1]  6  8 16 17 25 30 31 37 39 40 44 47 48 49 50 53 57 67 69 72 74 86 88 91 99
#> 
#> 
#> [[4]]
#> [[4]]$training
#>  [1]   2   3   7  10  11  12  13  14  15  17  18  19  20  21  22  23  24  25  26
#> [20]  27  28  29  30  31  33  36  37  39  41  43  44  45  47  48  50  51  52  55
#> [39]  56  57  58  62  63  65  66  67  68  69  70  72  74  75  76  77  78  79  80
#> [58]  81  82  83  84  85  86  87  88  89  91  92  93  94  96  97  98  99 100
#> 
#> [[4]]$test
#>  [1]  1  4  5  6  8  9 16 32 34 35 38 40 42 46 49 53 54 59 60 61 64 71 73 90 95
#> 
#> 

subject_id <- c("Alice", "Bob", "Charlie", "Eve")
random_subsampling(1:4, iterations = 2, test_size = 0.25, keep_together = subject_id)
#> [[1]]
#> [[1]]$training
#> [1] 1 2 3
#> 
#> [[1]]$test
#> [1] 4
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
