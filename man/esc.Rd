\encoding{UTF-8}
\name{esc-package}
\alias{esc-package}
\alias{esc}
\docType{package}
\title{Effect Size Computation for Meta Analysis}
\description{
  This is an R implementation of the web-based 'Practical Meta-Analysis Effect Size Calculator' from David B. Wilson. \cr \cr
  Based on the input, the effect size can be returned as standardized mean difference (\code{d}),
  the correlation coefficient effect size \code{r}, the odds ratio or log odds effect size.
  \cr \cr
  The return value of all functions has the same structure:
  \itemize{
    \item The effect size, whether being \code{d}, \code{r}, odds ratios or logits, is always named \code{es}.
    \item The standard error of the effect size, \code{se}.
    \item The variance of the effect size, \code{var}.
    \item The lower and upper confidence limits \code{ci.lo} and \code{ci.hi}.
    \item The weight factor, based on the inverse-variance, \code{w}.
  }
}
\details{
\tabular{rl}{
Package: \tab esc\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2016-06-30\cr
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
