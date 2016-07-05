#' @title Compute Hedges' g from effect size d
#' @name hedges_g
#'
#' @description Compute Hedges' g from effect size d.
#'
#' @param es The effect size \code{d}.
#' @param totaln The total sample size.
#'
#' @return The Hedges' g coefficient.
#'
#' @references Lipsey MW, Wilson DB. 2001. Practical meta-analysis. Thousand Oaks, Calif: Sage Publications
#'             \cr \cr
#'             Wilson DB. 2016. Formulas Used by the "Practical Meta-Analysis Effect Size Calculator". Unpublished manuscript: George Mason University
#'             \cr \cr
#'             Hedges LV. 1981. Distribution theory for Glass's estimator of effect size and related estimators. Journal of Educational Statistics 6: 107–128. \doi{10.3102/10769986006002107}
#'
#' @examples
#' hedges_g(0.75, 50)
#'
#' @export
hedges_g <- function(es, totaln) es * sssbc(totaln)
