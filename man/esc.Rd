\encoding{UTF-8}
\name{esc-package}
\alias{esc-package}
\alias{esc}
\docType{package}
\title{Effect Size Computation for Meta Analysis}
\description{
  This is an R implementation of the web-based 'Practical Meta-Analysis Effect Size Calculator' from David B. Wilson. \cr \cr
  Based on the input, the effect size can be returned as standardized mean difference (\code{d}),
  Hedges' \code{g}, correlation coefficient effect size \code{r} or Fisher's transformation \code{z},
  odds ratio or log odds effect size.
  \cr \cr
  \strong{Return values} \cr \cr
  The return value of all functions has the same structure:
  \itemize{
    \item The effect size, whether being \code{d}, \code{g}, \code{r}, (Cox) odds ratios or (Cox) logits, is always named \code{es}.
    \item The standard error of the effect size, \code{se}.
    \item The variance of the effect size, \code{var}.
    \item The lower and upper confidence limits \code{ci.lo} and \code{ci.hi}.
    \item The weight factor, based on the inverse-variance, \code{w}.
    \item The total sample size \code{totaln}.
    \item The effect size measure, \code{measure}, which is typically specified via the \code{es.type}-argument.
    \item Information on the effect-size conversion, \code{info}.
    \item A string with the study name, if the \code{study}-argument was specified in function calls.

  }
  \strong{Correlation Effect Size} \cr \cr
  If the correlation effect size \code{r} is computed, the transformed Fisher's z and their confidence
  intervals are also returned. The variance and standard error for the correlation effect size r are always
  based on Fisher's transformation.
  \cr \cr
  \strong{Odds Ratio Effect Size} \cr \cr
  For odds ratios, the variance and standard error are always returned on the log-scale!
  \cr \cr
  \strong{Preparing an Effect Size Data Frame for Meta-Analysis} \cr \cr
  The results of the effect size calculation functions in this package are returned as list with
  a \code{esc}-class attribute. The \code{\link{combine_esc}}-function takes one or more
  of these \code{esc}-objects and combines them into a \code{data.frame} that can be
  used as argument for further use, for instance with the \code{\link[metafor]{rma}}-function.
  \preformatted{
    e1 <- esc_2x2(grp1yes = 30, grp1no = 50, grp2yes = 40,
                  grp2no = 45, study = "Study 1")
    e2 <- esc_2x2(grp1yes = 30, grp1no = 50, grp2yes = 40, grp2no = 45,
                  es.type = "or", study = "Study 2")
    e3 <- esc_t(p = 0.03, grp1n = 100, grp2n = 150, study = "Study 3")
    e4 <- esc_mean_sd(grp1m = 7, grp1sd = 2, grp1n = 50, grp2m = 9, grp2sd = 3,
                      grp2n = 60, es.type = "logit", study = "Study 4")

    mydat <- combine_esc(e1, e2, e3, e4)

    metafor::rma(yi = es, sei = se, method = "REML", data = mydat)
  }
}
\details{
\tabular{rl}{
Package: \tab esc\cr
Type: \tab Package\cr
Version: \tab 0.2.0\cr
Date: \tab 2016-09-07\cr
License: \tab GPL-3\cr
}
}
\author{
Daniel Lüdecke \email{d.luedecke@uke.de}
}
\references{
  \itemize{
    \item Github: \url{https://github.com/sjPlot/esc}
    \item Bug reports: \url{https://github.com/sjPlot/esc/issues}
  }
}
