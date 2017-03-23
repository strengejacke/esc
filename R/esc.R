#' @title Generate effect size data frame from other data
#' @name effect_sizes
#'
#' @description This method computes any effect size from raw values from a data
#'              frame. Convenient method to compute multiple effect sizes at once,
#'              when the required information to calculate effects sizes are
#'              stored in a table (i.e. data frame).
#'
#' @param data A data frame with columns that contain the values that are passed
#'          to one of the \pkg{esc}-functions.
#' @param ... Named arguments. The name (left-hand side) is the name of one of
#'          \pkg{esc} functions' argument, the argument (right-hand side) is the
#'          name of the column in \code{data} that holds the data values.
#'          See 'Examples'.
#' @param fun Name of one of the \pkg{esc}-functions, as string, where arguments
#'          in \code{...} are passed to. May either be the full function name
#'          (like \code{"esc_t"} or \code{"esc_2x2"}) or the funcion name
#'          \emph{without} the suffix \code{"esc_"} (like \code{"t"} or \code{"2x2"}).
#' @inheritParams esc_beta
#'
#' @details This function rowwise iterates \code{data} and calls the function named
#'         in \code{type} for the values taken from each row of \code{data}.
#'         The column names in \code{data} that contain the necessary values to compute
#'         the effect sizes should be passed as unquoted value for the arguments.
#'         The argument names should match those arguments for the esc-function
#'         that should be called from within \code{effect_sizes()}.
#'         \cr \cr
#'         Example: \cr
#'         If you want to compute effect sizes from chi-squared values, you
#'         would call \code{esc_chisq()}. This function name is used for the
#'         \code{type}-argument: \code{type = "esc_chisq"}. \code{esc_chisq()}
#'         requires one of \code{chisq} or \code{p} as arguments, and \code{totaln}.
#'         Now \code{data} must have columns with values for either \code{chisq}
#'         or \code{p}, and \code{effect_sizes()} automatically selects the
#'         first non-missing value from \code{data} (see 'Examples').
#'
#' @return A data frame with the effect sizes computed for all data from \code{data}.
#'
#' @examples
#' tmp <- data.frame(
#'   tvalue = c(3.3, 2.9, 2.3),
#'   n = c(250, 200, 210),
#'   studyname = c("Study 1", "Study 2", "Study 3")
#' )
#' effect_sizes(tmp, t = tvalue, totaln = n, study = studyname, fun = "esc_t")
#'
#' # missing effect size results are dropped,
#' # shorter function name, calls "esc_t()"
#' tmp <- data.frame(
#'   tvalue = c(3.3, 2.9, NA, 2.3),
#'   n = c(250, 200, 210, 210),
#'   studyname = c("Study 1", "Study 2", NA, "Study 4")
#' )
#' effect_sizes(tmp, t = tvalue, totaln = n, study = studyname, fun = "t")
#'
#'
#' tmp <- data.frame(
#'   coefficient = c(0.4, 0.2, 0.6),
#'   se = c(.15, .1, .2),
#'   treat = c(50, 60, 50),
#'   cntrl = c(45, 70, 40),
#'   author = c("Smith 2000", "Smith 2010 2", "Smith 2012")
#' )
#' effect_sizes(tmp, beta = coefficient, sdy = se, grp1n = treat, grp2n = cntrl,
#'     study = author, fun = "esc_beta", es.type = "or")
#'
#' # the "esc_chisq" function requires *either* the chisq-argument *or*
#' # the pval-argument. If at least one of these values is present,
#' # effect size can be calculated. You can specify both arguments,
#' # and the first non-missing required value from "data" is taken.
#' tmp <- data.frame(
#'   chisqquared = c(NA, NA, 3.3, NA, 2.9),
#'   pval = c(.003, .05, NA, .12, NA),
#'   n = c(250, 200, 210, 150, 180),
#'   studyname = c("Study 1", "Study 2", "Study 3", "Study 4", "Study 5")
#' )
#' effect_sizes(tmp, chisq = chisqquared, p = pval, totaln = n,
#'              study = studyname, fun = "esc_chisq")
#'
#' # if all required information are missing, data will be removed
#' tmp <- data.frame(
#'   chisqquared = c(NA, NA, 3.3, NA, NA),
#'   pval = c(.003, .05, NA, .12, NA),
#'   n = c(250, 200, 210, 150, 180),
#'   studyname = c("Study 1", "Study 2", "Study 3", "Study 4", "Study 5")
#' )
#' effect_sizes(tmp, chisq = chisqquared, p = pval, totaln = n,
#'              study = studyname, fun = "chisq")
#'
#'
#' @importFrom tibble as_tibble
#' @importFrom dplyr slice_
#' @importFrom purrr map_df
#' @export
effect_sizes <- function(data, ..., fun, es.type = c("d", "g", "or", "logit", "r", "cox.or", "cox.log")) {
  es.type <- match.arg(es.type)

  # check function name
  if (substr(fun, 1, 4) != "esc_") fun <- paste0("esc_", fun)

  # get arguments
  params <- match.call(expand.dots = FALSE)$`...`

  # and get columns names
  cols <- unname(unlist(lapply(params, function(x) as.character(x))))

  # make sure all numerics are numeric! and, only
  # return necessary data rows
  data <- suppressWarnings(
    purrr::map_df(data[, cols], function(x) {
      if (!all(is.na(as.numeric(x))))
        as.numeric(x)
      else
        x
    })
  )

  results <- list()

  # iterate data frame rowwise
  for (i in 1:nrow(data)) {
    # get argument names
    argnames <- names(params)
    # get required data row
    datarow <- data[i, ]
    # copy values as list, for do.call
    callparams <- as.list(datarow)
    # give argument proper names
    names(callparams) <- argnames
    # call esc-method
    results[[length(results) + 1]] <-
      suppressWarnings(do.call(eval(fun), args = c(callparams, es.type = es.type)))
  }

  # returned results as combined tibble
  results <- suppressWarnings(purrr::map_df(results, function(x) tibble::as_tibble(x)))

  # do we have any missing effect sizes?
  if (anyNA(results$es) || any(is.infinite(results$es))) {
    # get studies with missing or Inf effect sizes
    mes <- which(is.na(results$es) | is.infinite(results$es))
    # warn user that some entries have been dropped
    warning(sprintf(
      "Row(s) %s in `data` had missing values and were removed.",
      paste0(mes, collapse = ", ")
    ),
    call. = F)
    # remove missings
    results <- dplyr::slice(results, -mes)
  }

  purrr::map_df(results, function(x) if (is.factor(x)) as.character(x) else x)
}