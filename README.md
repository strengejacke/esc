esc - Effect Size Computation for Meta Analysis
------------------------------------------------------------------------------
This is an R implementation of the web-based 'Practical Meta-Analysis Effect Size Calculator' from David B. Wilson. The original calculator can be found at http://www.campbellcollaboration.org/escalc/html/EffectSizeCalculator-Home.php.

Based on the input, the effect size can be returned as standardized mean difference (`d`), the correlation coefficient effect size `r`, the odds ratio or log odds effect size.

### Return values

The return value of all functions has the same structure:

* The effect size, whether being `d`, `r`, odds ratios or logits, is always named `es`.
* The standard error of the effect size, `se`.
* The variance of the effect size, `var`.
* The lower and upper confidence limits `ci.lo` and `ci.hi`.
* The weight factor, based on the inverse-variance, `w`.

#### Correlation Effect Size

If the correlation effect size `r` is computed, the transformed Fisher's r and their confidence intervals are also returned. The variance and standard error for the correlation effect size r are always based on Fisher's transformation.

#### Odds Ratio Effect Size

For odds ratios, the variance and standard error are always returned on the log-scale!

**esc** is still under development, i.e. not all effect size computation options are implemented yet. The remaining options will follow in the course of the next days.

## Installation

### Latest development build

To install the latest development snapshot (see latest changes below), type following commands into the R console:

```r
library(devtools)
devtools::install_github("sjPlot/esc")
```

## Citation

In case you want / have to cite my package, please use `citation('esc')` for citation information.
