#' @title Convert correlation coefficient r into Fisher's z
#' @name convert_r2z
#'
#' @description Convert correlation coefficient r into Fisher's z.
#'
#' @param r The correlation coefficient.
#'
#' @return The transformed Fisher's \code{z}.
#'
#' @references Lipsey MW, Wilson DB. 2001. Practical meta-analysis. Thousand Oaks, Calif: Sage Publications
#'             \cr \cr
#'             Wilson DB. 2016. Formulas Used by the "Practical Meta-Analysis Effect Size Calculator". Unpublished manuscript: George Mason University
#'
#' @examples
#' convert_r2z(.03)
#'
#' @export
convert_r2z <- function(r) return(.5 * log((1 + r) / (1 - r)))


#' @title Convert Fisher's z into correlation coefficient r
#' @name convert_z2r
#'
#' @description Convert Fisher's z into correlation coefficient r.
#'
#' @param z Fisher's \code{z}-value.
#'
#' @return The back-transformed correlation coefficient \code{r}.
#'
#' @references Lipsey MW, Wilson DB. 2001. Practical meta-analysis. Thousand Oaks, Calif: Sage Publications
#'             \cr \cr
#'             Wilson DB. 2016. Formulas Used by the "Practical Meta-Analysis Effect Size Calculator". Unpublished manuscript: George Mason University
#'
#' @examples
#' convert_z2r(.03)
#'
#' @export
convert_z2r <- function(z) return((exp(2 * z) - 1) / (exp(2 * z) + 1))



