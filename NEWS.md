# esc 0.5.1

## General

* Deprecated functions have been removed.
* Update website (https://strengejacke.github.io/esc).

# esc 0.5.0

## General

* Reduce package dependencies.
* Function that convert from one effect size into another are now prefixed with `convert_*()`, instead of `esc_*()` (e.g. `convert_d2r()` instead of `esc_d2r()`).

## Changes to functions

* `esc_mean_sd()`, `esc_mean_se()` and `esc_mean_gain()` can now calculate effect-sizes for within-subject-designs, when argument `r` is given.

# esc 0.4.1

## Bug fixes

* Fix issue with wrong computation of confidence intervals when converting effect sizes from _d_ to _f_ (`esc_d2f()`).

# esc 0.4.0

## General

* Functions now also return effect sizes Cohen's f or eta squared.

## New functions

* More functions to convert effect size into other effect size measures (`cohens_f()`, `odds_ratio()`, `log_odds()`, `pearsons_r()`, `eta_squared()` and `cohens_d()`).

## Bug fixes

* In rare cases, `esc_mean_sd()` tried to calculate the square root of negative values when computing the pooled standard deviance. In such cases, an alternative formular for the pooled SD is used.

# esc 0.3.2

## Bug fixes

* `combine_esc()` did not work for effect sizes computed with `esc_t()`, when the returned effect size was `r`.

# esc 0.3.1

## Bug fixes

* `effect_sizes()` did not always properly extract study names and used numeric indices instead character values.

# esc 0.3.0

## New functions

* `effect_sizes()` to generate a data frame with results of multiple effect size computations, based on data from another data frame.

## Bug fixes
* Fix control statements with condition with greater than one, which currently give a warning, but in future R versions may throw an error.

# esc 0.2.0

## New functions

* `write_esc()` as a convenient wrapper to write results to an Excel csv-file.