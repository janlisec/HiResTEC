#' @title groupval.
#'
#' @description
#' \code{groupval} is a dirty rewrite of the same function available in xcms.
#'
#' @details Comment on 2019-05-14 by JL: 'xcms' needs to be put to suggest due to mzR
#'     on request by CRAN --> replaced groupval from xcms by own function.
#'
#' @param xg xcmsSet object.
#'
#' @keywords internal
#' @noRd
groupval <- function(xg = NULL) {
  gv <- matrix(NA, nrow = nrow(xg@groups), ncol = nrow(xg@phenoData))
  for (i in 1:nrow(gv)) {
    idx <- xg@peaks[xg@groupidx[[i]], "sample"]
    if (!any(duplicated(idx))) {
      gv[i, idx] <- xg@peaks[xg@groupidx[[i]], "maxo"]
    } else {
      flt <- !duplicated(idx)
      gv[i, idx[flt]] <- xg@peaks[xg@groupidx[[i]][flt], "maxo"]
    }
  }
  return(gv)
}
