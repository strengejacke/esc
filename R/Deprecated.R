#' @rdname convert_d2or
#' @export
esc_d2or <- function(d, se, v, totaln, es.type = c("logit", "cox"), info = NULL, study = NULL) {
  .Deprecated("convert_d2or")
  convert_d2or(d, se, v, totaln, es.type, info, study)
}



#' @rdname convert_d2f
#' @export
esc_d2f <- function(d, se, v, totaln, info = NULL, study = NULL) {
  .Deprecated("convert_d2f")
  convert_d2f(d, se, v, totaln, info, study)
}




#' @rdname convert_d2logit
#' @export
esc_d2logit <- function(d, se, v, totaln, es.type = c("logit", "cox"), info = NULL, study = NULL) {
  .Deprecated("convert_d2logit")
  convert_d2logit(d, se, v, totaln, es.type, info, study)
}



#' @rdname convert_or2d
#' @export
esc_or2d <- function(or, se, v, totaln, es.type = c("d", "cox.d", "g", "f", "eta"), info = NULL, study = NULL) {
  .Deprecated("convert_or2d")
  convert_or2d(or, se, v, totaln, es.type, info, study)
}



#' @rdname convert_d2r
#' @export
esc_d2r <- function(d, se, v, grp1n, grp2n, info = NULL, study = NULL) {
  .Deprecated("convert_d2r")
  convert_d2r(d, se, v, grp1n, grp2n, info, study)
}



#' @rdname convert_d2etasq
#' @export
esc_d2etasq <- function(d, se, v, grp1n, grp2n, info = NULL, study = NULL) {
  .Deprecated("convert_d2etasq")
  convert_d2etasq(d, se, v, grp1n, grp2n, info, study)
}




#' @rdname convert_z2r
#' @export
esc_z2r <- function(z) {
  .Deprecated("convert_z2r")
  convert_z2r(z)
}



#' @rdname convert_r2z
#' @export
esc_r2z <- function(r) {
  .Deprecated("convert_r2z")
  convert_r2z(r)
}
