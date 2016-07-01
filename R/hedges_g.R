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
#'
#' @examples
#' hedges_g(0.75, 50)
#'
#' @export
hedges_g <- function(es, totaln) es * sssbc(totaln)
