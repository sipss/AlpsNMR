# Create method for NMR data analysis

Create method for NMR data analysis

## Usage

``` r
new_nmr_data_analysis_method(
  train_evaluate_model,
  train_evaluate_model_params_inner,
  choose_best_inner,
  train_evaluate_model_params_outer,
  train_evaluate_model_digest_outer
)
```

## Arguments

- train_evaluate_model:

  A function. The `train_evaluate_model` must have the following
  signature:

          function(x_train, y_train, identity_train, x_test, y_test, identity_test, ...)

  The `x_train` and `y_train` (and their test counterparts) are
  self-explanatory.

  The `identity_` arguments are expected to be factors. They can be used
  for instance with a callback that uses
  [mixOmics::plsda](https://rdrr.io/pkg/mixOmics/man/plsda.html) in a
  `multilevel` approach for longitudinal studies. In those studies the
  `identity` would be an identifier of the subject.

  The `...` arguments are free to be defined for each
  `train_evaluate_model`.

- train_evaluate_model_params_inner, train_evaluate_model_params_outer:

  A list with additional arguments to pass to `train_evaluate_model`
  either in the inner cv loop or in the outer cv loop.

- choose_best_inner:

  A function with a single argument:

      function(inner_cv_results)

  The argument is a list of `train_evaluate_model` outputs. The return
  value of must be a list with at least an element named
  `train_evaluate_model_args`. `train_evaluate_model_args` must be a
  named list.

  - Each element must be named as one of the `train_evaluate_model`
    arguments.

  - Each element must be a vector as long as the number of outer
    cross-validations.

  - The values of each vector must be the values that the
    `train_evaluate_model` argument must take on each outer
    cross-validation iteration Additional list elements can be returned
    and will be given back to the user

- train_evaluate_model_digest_outer:

  A function with a single argument:

  function(outer_cv_results)

  The argument is a list of `train_evaluate_model` outputs in outer
  cross-validation. The return value is returned by `nmr_data_analysis`

## Value

An object encapsulating the method dependent functions that can be used
with
[nmr_data_analysis](https://sipss.github.io/AlpsNMR/reference/nmr_data_analysis.md)
