#' @title Extract multiple ion chromatograms from mass spectrometry data.
#'
#' @description \code{getMultipleBPC} will extract multiple BPCs from an `xcmsRaw`
#'     or `xcmsRawLike` object for a vector of mz within the limits given by rt,
#'     rt_dev and mz_dev.
#'
#' @details While there are other functions to extract BPC information from raw data,
#'     this one is particularly useful to get all traces belonging to a isotopologue
#'     group. It will attach several derived values to the results object,
#'     i.e. describing the observed mass shift (deviation from expected value) which
#'     is helpful in QC for non-targeted tracer analyses.
#'     While the `mz` and `mz_dev` parameters  can be vectorized, the `rt` and
#'     `rt_dev` values will be consistently used for all ion traces.
#'
#' @param x `xcmsRaw` or `xcmsRawLike` object.
#' @param mz Numeric vector of masses or NULL (default) to return the overall BPC.
#' @param mz_dev Allowed mass deviations (can be a single numeric, a vector, a matrix
#'     with one row (lower bound, upper bound) or a matrix with \code{length(mz)} rows
#'     giving lower and upper bound for each mz).
#' @param rt Target retention time or NULL (default) to use full time range.
#' @param rt_dev Allowed time deviation (if rt is specified).
#' @param zeroVal Set values <=0 to NA or keep as is with NULL.
#' @param smooth Window size for moving average smoother, 0 = no smoothing.
#' @param returnEIC Return EIC instead of BPC?
#'
#' @return A matrix with scan wise (rows) intensities for all requested masses (columns)
#'     as either EIC or BPC.
#'
#' @examples
#' raw <- HiResTEC::raw
#' # search for mz = 556.263 and its isotopic traces
#' mz <- 556.263 + c(0:3) * 1.0034
#' getMultipleBPC(x = raw[[1]], mz = mz, mz_dev = 0.04, rt = 1026)
#'
#' @export
#'
#' @useDynLib HiResTEC, .registration = TRUE
#'
#' @references Uses C code modified from XCMS (see \code{citation("xcms")}).
#'
getMultipleBPC <- function(x, mz = NULL, mz_dev = 0.005, rt = NULL, rt_dev = 2, zeroVal = NA, smooth = 0, returnEIC = FALSE) {
  #

  # use full rt if rt = NULL
  if (is.null(rt)) {
    rt <- median(range(x@scantime))
    rt_dev <- diff(range(x@scantime))/2
  }

  scans <- which(abs(x@scantime - rt) <= rt_dev)
  if (length(scans) < 1) {
    return(NULL)
  }

  scanrange <- range(scans)
  if (!is.double(x@env$mz)) {
    x@env$mz <- as.double(x@env$mz)
  }
  if (!is.double(x@env$intensity)) {
    x@env$intensity <- as.double(x@env$intensity)
  }
  if (!is.integer(x@scanindex)) {
    x@scanindex <- as.integer(x@scanindex)
  }

  # return TIC for mz = NULL
  if (is.null(mz)) {
    mz <- median(x@mzrange)
    mz_dev <- diff(range(x@mzrange))/2
  } else {
    # convert to vector in case user provided a data.frame
    mz <- as.vector(unlist(mz))
    mz <- as.numeric(mz)
  }

  ## prepare mz_dev
  nmz <- length(mz)
  nmzdev <- length(mz_dev)
  isNumeric <- is.numeric(mz_dev)
  if (isNumeric && nmzdev == 1) {
    mzdev_lower <- mzdev_upper <- rep(mz_dev, nmz)
  } else if (isNumeric && nmzdev == nmz) {
    mzdev_lower <- mzdev_upper <- mz_dev
  } else if (isNumeric && is.matrix(mz_dev) && nrow(mz_dev) == 1) {
    mzdev_lower <- rep(mz_dev[, 1], nmz)
    mzdev_upper <- rep(mz_dev[, 2], nmz)
  } else if (isNumeric && is.matrix(mz_dev) && nrow(mz_dev) == nmz) {
    mzdev_lower <- mz_dev[, 1]
    mzdev_upper <- mz_dev[, 2]
  } else {
    stop("mz_dev incorrectly specified")
  }

  # extract for all mz int and mz@int=max for scan range
  tmp <- .Call("getMultipleBPC_C", x@env$mz, x@env$intensity, x@scanindex,
    as.double(mz), as.double(mzdev_lower), as.double(mzdev_upper), as.integer(scanrange),
    as.integer(length(x@scantime)), as.integer(smooth), as.integer(returnEIC),
    PACKAGE = "HiResTEC"
  )
  res <- matrix(tmp$intensity, nrow = length(scans), ncol = length(mz), byrow = TRUE)
  if (!is.null(zeroVal)) res[res == 0] <- zeroVal
  rownames(res) <- round(x@scantime[scans], 2)
  colnames(res) <- round(mz, 4)
  attr(res, "rt") <- x@scantime[scans]
  attr(res, "mz") <- mz
  attr(res, "mz_dev") <- mz_dev
  ln <- which.max(apply(res, 1, sum, na.rm = T)) # scan whose TIC is max
  attr(res, "maxBPC") <- ln
  ln <- max(c(1, ln - 2)):min(c(length(scans), ln + 2)) # expand max scan to left and right (5 scans total)
  # cl <- grep("m",colnames(tmp))
  mzmat <- matrix(tmp$mz, nrow = length(scans), ncol = length(mz), byrow = TRUE)[ln, ] # matrix of accurate mz values [5 x nmz]
  mzmat[mzmat == 0] <- NA
  mzmat <- t(mzmat) - mz # [nmz x 5]
  #mzmat <- round(1000 * Biobase::rowMedians(mzmat, na.rm = TRUE)) # slightly faster, we need Biobase anyway due to xcms
  mzmat <- round(1000 * apply(mzmat, MARGIN=1, FUN=median, na.rm=TRUE))
  # mzmat <- round(1000*apply(mzmat, 1, median, na.rm=TRUE),1) # vector of length nmz
  mzmat[!is.finite(mzmat)] <- NA
  if (length(mzmat) < length(mz)) mzmat <- rep(NA, length(mz))
  attr(res, "mass_defect") <- mzmat
  return(res)
}
