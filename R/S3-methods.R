#' @export
print.esc <- function(x, ...) {
  if (!is.null(x$study))
    cat(sprintf("\n%s (n=%i)\n\n", x$study, x$totaln))
  else
    cat("\nEffect Size Calculation for Meta Analysis\n\n")
  cat(sprintf("     Conversion: %s\n", x$info))
  cat(sprintf("    Effect Size: %8.4f\n", x$es))
  cat(sprintf(" Standard Error: %8.4f\n", x$se))
  cat(sprintf("       Variance: %8.4f\n", x$var))
  cat(sprintf("       Lower CI: %8.4f\n", x$ci.lo))
  cat(sprintf("       Upper CI: %8.4f\n", x$ci.hi))
  cat(sprintf("         Weight: %8.4f\n", x$w))
  if (!is.null(x$zr)) cat(sprintf("     Fisher's z: %8.4f\n", x$zr))
  if (!is.null(x$ci.lo.zr)) cat(sprintf("      Lower CIz: %8.4f\n", x$ci.lo.zr))
  if (!is.null(x$ci.hi.zr)) cat(sprintf("      Upper CIz: %8.4f\n", x$ci.hi.zr))
}


#' @export
as.data.frame.esc <- function(x, ...) {
  if (any(class(x) == "esc_d2r")) {
    if (is.null(x$study)) x$study <- ""
    return(data.frame(study = x$study, es = x$es, weight = x$w, sample.size = x$totaln, se = x$se,
                      var = x$var, ci.lo = x$ci.lo, ci.hi = x$ci.hi, fishers.z = x$zr,
                      ci.lo.z = x$ci.lo.zr, ci.hi.z = x$ci.hi.zr, measure = x$measure))
  } else {
    if (is.null(x$study)) x$study <- ""
    return(data.frame(study = x$study, es = x$es, weight = x$w, sample.size = x$totaln, se = x$se,
                      var = x$var, ci.lo = x$ci.lo, ci.hi = x$ci.hi, measure = x$measure))
  }
}
