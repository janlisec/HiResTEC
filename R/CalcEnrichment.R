#' @title CalcEnrichment.
#' @description \code{CalcEnrichment} will calculate from a vector of isotope intensities the enrichment of 13C.
#' @param x Vector of intensities.
#' @param sig Round enrichment to this precision.
#' @param robust Provide an intensity threshold that at least one ion has to exceed. Otherwise NA is returned.
#' @return The enrichment, defined as the ratio of 13C/total_C.
#' @keywords internal
#' @noRd
CalcEnrichment <- function(x = c(100, 10, 1), sig = 4, robust = 0) {
  if (any(is.finite(x))) {
    if (!any(x > robust, na.rm = T)) {
      return(NA)
    } else {
      n <- length(x) - 1
      return(100 * round(sum(seq(0, n) * (x / n), na.rm = T) / sum(x, na.rm = T), sig))
    }
  } else {
    return(NA)
  }
}
