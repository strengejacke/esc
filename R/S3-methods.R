#' @export
print.esc <- function(x, ...) {
  cat("\nEffect Size Calculation for Meta Analysis\n\n")
  cat(sprintf("     Conversion: %s\n", x$info))
  cat(sprintf("    Effect Size: %8.4f\n", x$es))
  cat(sprintf(" Standard Error: %8.4f\n", x$se))
  cat(sprintf("       Variance: %8.4f\n", x$var))
  cat(sprintf("       Lower CI: %8.4f\n", x$ci.lo))
  cat(sprintf("       Upper CI: %8.4f\n", x$ci.hi))
  cat(sprintf("         Weight: %8.4f\n", x$w))
}