#' @title DetermineElementComposition.
#' @description \code{DetermineElementComposition} will compute the Mass
#'    Distribution Vectors of isotopologues based on measured intensities
#'    at [13C]=1\%.
#' @details Not yet.
#' @param int Measured intensities.
#' @param mz Mass ob peak.
#' @param ionization Currently only APCI is supported.
#' @return The most likely combinations of C and Si atoms for this formula.
#' @keywords internal
#' @noRd
DetermineElementComposition <- function(int = NULL, mz = 300, ionization = "APCI") {
  # set nC and nSi range based on mz
  # ToDo
  nC_rng <- 1:ceiling(mz / 12)

  # set reference MID
  ref <- c(100, rep(0, length(int) - 1))

  # find best fitting MID brute force
  if (ionization == "APCI") {
    test <- lapply(nC_rng, function(nCbio) {
      sapply(0:min(c(floor(nCbio / 3), 8)), function(nSi) {
        mid <- CorMID(int = int, fml = paste0("C", nCbio + 3 * nSi, "Si", nSi))
        return(sqrt(sum((mid - ref[1:length(mid)])^2)))
      })
    })
  }
  # return respective nC and nSi
  nC <- which.min(sapply(test, min))
  nSi <- which.min(test[[nC]])
  return(paste0("C", nC, "Si", nSi))
}
