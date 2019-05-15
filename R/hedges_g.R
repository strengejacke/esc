#' @title Convert effect sizes
#' @name hedges_g
#'
#' @description Convert between different effect sized.
#'
#' @param d,r,f,eta,or,logit A scalar or vector with effect size(s).
#' @param totaln A vector of total sample size(s).
#'
#' @return The requested effect size.
#'
#' @references Lipsey MW, Wilson DB. 2001. Practical meta-analysis. Thousand Oaks, Calif: Sage Publications
#'             \cr \cr
#'             Wilson DB. 2016. Formulas Used by the "Practical Meta-Analysis Effect Size Calculator". Unpublished manuscript: George Mason University
#'             \cr \cr
#'             Hedges LV. 1981. Distribution theory for Glass's estimator of effect size and related estimators. Journal of Educational Statistics 6: 107–128.
#'             \cr \cr
#'             Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. 2009. Introduction to Meta-Analysis. Chichester, West Sussex, UK: Wiley
#'             \cr \cr
#'             Cohen J. 1988. Statistical Power Analysis for the Behavioral Sciences. 2nd ed. Hillsdale, NJ: Erlbaum
#'
#' @examples
#' # convert from d to Hedges' g or odds ratio
#' hedges_g(d = 0.75, totaln = 50)
#' odds_ratio(d = .3)
#'
#' # convert from odds ratio to eta_squared
#' eta_squared(or = 2.3)
#'
#' # convert from f or r to d
#' cohens_d(f = .3)
#' cohens_d(r = .25)
#'
#' # functions are vectorized
#' hedges_g(c(0.75, .3), c(50, 70))
#' cohens_f(r = c(.1, .2, .3))
#'
#' @export
hedges_g <- function(d, totaln) {
  mapply(function(.x, .y) .x * sssbc(.y), d, totaln)
}


#' @rdname hedges_g
#' @export
eta_squared <- function(d, r, f, or, logit) {
  d <- convert_to_d(d = d, r = r, f = f, or = or, logit = logit)
  (d^2 / (1 + d^2)) / 2
}


#' @rdname hedges_g
#' @export
cohens_f <- function(d, r, eta, or, logit) {
  d <- convert_to_d(d = d, r = r, eta = eta, or = or, logit = logit)
  d / 2
}


#' @rdname hedges_g
#' @export
cohens_d <- function(f, r, eta, or, logit) {
  convert_to_d(r = r, eta = eta, f = f, or = or, logit = logit)
}


#' @rdname hedges_g
#' @export
pearsons_r <- function(d, eta, f, or, logit) {
  d <- convert_to_d(d = d, eta = eta, f = f, or = or, logit = logit)
  sqrt(d^2 / (d^2 + 4))
}


#' @rdname hedges_g
#' @export
log_odds <- function(d, eta, f, or, r) {
  d <- convert_to_d(d = d, eta = eta, f = f, or = or, r = r)
  d * pi / sqrt(3)
}


#' @rdname hedges_g
#' @export
odds_ratio <- function(d, eta, f, logit, r) {
  d <- convert_to_d(d = d, eta = eta, f = f, logit = logit, r = r)
  exp(d * pi / sqrt(3) )
}


convert_to_d <- function(d, r, f, eta, or, logit) {
  # convert effect size into d
  if (!missing(d)) d <- d
  if (!missing(r)) d <- ((2 * r) / sqrt(1 - (r^2)))
  if (!missing(f)) d <- 2 * f
  if (!missing(eta)) d <- 2 * (sqrt(eta / (1 - eta)))
  if (!missing(or)) d <- (log(or) / (pi / sqrt(3)))
  if (!missing(logit)) d <- (logit / (pi / sqrt(3)))

  d
}

