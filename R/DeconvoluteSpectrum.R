#' @title DeconvoluteSpectrum.
#'
#' @description
#' \code{DeconvoluteSpectrum} will evaluate a list of xcmsRaw objects at a given time (rt) and potential mass (mz1).
#' The main purpose is to deconvolute the mass spectrum at rt including mz1.
#'
#' @details
#' Will test all mz at spectrum of base peak within range for co-apex, rt diff and ratio consistency/correlation over a set of samples.
#'
#' @param dat A list of xcmsRaws or an xcmsSet object.
#' @param rt Retention time to search for maxima.
#' @param rt_dev Allowed retention time window.
#' @param mz1 If specified, ensure that this mass is included in the spectrum (assumed base peak). NULL otherwise.
#' @param mz_dev Allowed mz deviation [Da].
#' @param use.mz.adjust Will adjust mz on an experiment wide basis.
#' @param ionization Either APCI or ESI. Choice will modify some internal parameters and checks performed.
#' @param smooth Smoothing parameter passed on to \link{getMultipleBPC}.
#'
#' @return
#' A pseudo spectrum at rt (containing mz1 if specified). Effectively a 2-column matrix (mz, int) with rt as attribute
#'
#' @examples
#' # Please use examples from previous versions as xcms (and xcms objects)
#' # are no longer supported during CRAN checks leading to package rejection
#' # if included (and I do not know a work around). :(
#'
#' @export
#'
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom graphics legend
#' @importFrom stats median
#' @importFrom stats cor
#' @importFrom stats sd
#'

DeconvoluteSpectrum <- function(dat = NULL, rt = NULL, rt_dev = 3, mz1 = NULL, mz_dev = 0.003, use.mz.adjust = FALSE, ionization = c("APCI", "ESI")[1], smooth = 0) {
  # putative parameters
  # rt_dev_max function tries to adjust rt_dev internally; however if this results in extrem values rt_dev_max is applied to limit this
  if (ionization == "APCI") {
    rt_dev_max <- 3 # how far may peaks deviate (between samples)
    rt_dev_min <- 0.5
    allowed_d_rt <- 0.15 # how far can co-apexes be apart from each other (intra sample)
    min_cor_rat <- 0.7 # how well have maxima to be correlated over samples
    nmax <- 1000 # number of peaks to consider in Deconvolution
  }
  if (ionization == "ESI") {
    rt_dev_max <- 6
    # rt_dev_min <- 1
    # allowed_d_rt <- 0.35
    # min_cor_rat <- 0.7
    rt_dev_min <- 2
    allowed_d_rt <- 0.5
    min_cor_rat <- 0.7
    nmax <- 1000
  }

  GetSpectrum <- function(x = NULL, rt = NULL, cutoff = 200, nmax = 150, se = 2, sort_int = FALSE) {
    # this is the scan where target mz1 is at max within all samples
    s <- which.min(abs(x@scantime - rt))
    if (s == length(x@scanindex)) s <- s - se
    if (s < se) s <- s + se + 1
    # these are potential co-eluting masses and their intensities
    m <- x@env$mz[(1 + x@scanindex[s]):x@scanindex[s + 1]]
    i <- x@env$intensity[(1 + x@scanindex[s]):x@scanindex[s + 1]]
    # m, which are considered for the spectrum should be
    # (i) above a cutoff
    flt <- which(i > cutoff)
    i <- i[flt]
    m <- m[flt]
    # # (ii) have there max in close proximity to s, i.e. they should not show higher int s+-se scans up/downstream
    # m_down <- x@env$mz[(1+x@scanindex[s-se]):x@scanindex[s-se+1]]
    # m_up <- x@env$mz[(1+x@scanindex[s+se]):x@scanindex[s+se+1]]
    # i_down <- x@env$intensity[(1+x@scanindex[s-se]):x@scanindex[s-se+1]]
    # i_up <- x@env$intensity[(1+x@scanindex[s+se]):x@scanindex[s+se+1]]
    # flt <- sapply(1:length(m), function(idx) {
    #   i_comp <- c(i_down[abs(m_down-m[idx])<mz_dev],i_up[abs(m_up-m[idx])<mz_dev])
    #   ifelse(length(i_comp)>=1, i_comp<i[idx], FALSE)
    # })
    # i <- i[flt]
    # m <- m[flt]
    # (iii) limited to a managable amount
    flt <- which(rank(-i) <= nmax)
    i <- i[flt]
    m <- m[flt]
    # return mz which fulfill above criteria
    if (sort_int) {
      return(m[order(i, decreasing = TRUE)])
    } else {
      return(m)
    }
  }

  if (is.null(mz1)) {
    # take the median of highest masses found in all provided samples
    mz1 <- median(sapply(dat, function(x) {
      GetSpectrum(x = x, rt = rt, sort_int = TRUE)[1]
    }))
  }

  # determine rt and int of max mz over exp
  if (!is.null(names(dat))) names(dat) <- NULL
  i_mz1 <- sapply(dat, function(x) {
    x <- getMultipleBPC(x = x, mz = mz1, mz_dev = mz_dev, rt = rt, rt_dev = rt_dev, zeroVal = 0, smooth = smooth)
    if (is.null(x)) {
      return(NA)
    } else {
      out <- x[attr(x, "maxBPC"), ]
      names(out) <- names(attr(x, "maxBPC"))
      return(out)
    }
  })

  # remove files which are empty in this region (e.g. no data recorded within time window or for defined mz-window)
  filt_files <- which(i_mz1 == 0 | is.na(i_mz1))
  if (length(filt_files) >= 1) {
    if (length(filt_files) == length(dat)) {
      spec <- matrix(NA, ncol = 2, dimnames = list(NULL, c("mz", "int")))[-1, ]
      attr(spec, "rt") <- rt
      return(spec)
    } else {
      i_mz1 <- i_mz1[-filt_files]
      dat <- dat[-filt_files]
    }
  }

  # readjust rt and rt_dev based on data
  if (!is.finite(median(as.numeric(names(i_mz1))[i_mz1 > 0], na.rm = T))) {
    #browser()
    message("ToDo: The rt should be adjusted.")
  } else {
    rt <- median(as.numeric(names(i_mz1))[i_mz1 > 0], na.rm = T)
  }

  # ==========================================================================================
  # [ToDo] think about readjusting rt and rt_dev based on this result
  # ==========================================================================================
  rt_dev <- ceiling(max(abs(as.numeric(names(i_mz1))[i_mz1 > 0] - rt)))

  # if rt_dev is still high (shouldn't necessary be larger than 2 for well aligned samples) try to fix
  if (rt_dev > rt_dev_max) {
    rt_dev <- rt_dev_max
  }
  if (rt_dev < rt_dev_min) {
    # happens if single samples are processed
    rt_dev <- rt_dev_min
  }

  # reevaluate mz1 with adjusted rt_dev and redo filtering
  i_mz1 <- sapply(dat, function(x) {
    x <- getMultipleBPC(x = x, mz = mz1, mz_dev = mz_dev, rt = rt, rt_dev = rt_dev, zeroVal = 0)
    if (is.null(x)) {
      return(NA)
    } else {
      out <- x[attr(x, "maxBPC"), ]
      names(out) <- names(attr(x, "maxBPC"))
      return(out)
    }
  })

  # remove files which are empty in this region (e.g. no data recorded within time window)
  filt_files <- which(i_mz1 == 0 | is.na(i_mz1))
  if (length(filt_files) >= 1) {
    if (length(filt_files) == length(dat)) {
      spec <- matrix(NA, ncol = 2, dimnames = list(NULL, c("mz", "int")))[-1, ]
      attr(spec, "rt") <- rt
      return(spec)
    } else {
      i_mz1 <- i_mz1[-filt_files]
      dat <- dat[-filt_files]
    }
  }

  # determine masses from spectrum
  bpm <- which.max(i_mz1)

  # ensure that mz1 is still contained
  # browser()
  mz2 <- GetSpectrum(x = dat[[bpm]], rt = as.numeric(names(i_mz1))[bpm], cutoff = min(c(200, 0.1 * max(i_mz1))), nmax = nmax)
  if (!any(abs(mz2 - mz1) < mz_dev)) mz2 <- sort(c(mz1, mz2))
  idx_mz1 <- which.min(abs(mz2 - mz1))

  # determine maxima for these masses from all files
  out <- lapply(dat, function(x) {
    getMultipleBPC(
      x = x, mz = mz2, mz_dev = mz_dev, rt = rt,
      ## rt_dev = if(smooth>0) rt_dev+diff(range(x@scantime))/length(x@scantime)*smooth else rt_dev, # adjust for loss due to smoothing
      rt_dev = rt_dev,
      zeroVal = NULL, smooth = smooth
    )
  })

  # get diffs for rt at max int between candidates (mz2 and base peak mz1)
  d_rt <- round(apply(sapply(1:length(out), function(i) {
    # attr(out[[i]],"rt")[apply(out[[i]],2,which.max)]-as.numeric(names(i_mz1)[i])
    attr(out[[i]], "rt")[apply(out[[i]], 2, which.max)] - attr(out[[i]], "rt")[which.max(out[[i]][, idx_mz1])]
  }), 1, median), 2)

  # get stable ratios and correlation between mz1 and all mz2's if more than 5 samples (as correlation is flawed otherwise)
  if (length(out) >= 5) {
    cor_rat <- round(apply(sapply(1:length(out), function(i) {
      flt <- rank(-out[[i]][, idx_mz1]) <= 10 & out[[i]][, idx_mz1] > 0
      if (sum(flt) >= 6) {
        suppressWarnings(
          cor(out[[i]][flt, , drop = F], out[[i]][flt, idx_mz1], use = "p")
        )
      } else {
        matrix(0, nrow = ncol(out[[i]]), ncol = 1)
      }
    }), 1, median, na.rm = T), 2)
    msg_warn <- NULL
  } else {
    # set cor_rat=1 artificially
    cor_rat <- rep(1, length(d_rt))
    msg_warn <- "Less than 5 samples provided, no correlation testing was performed."
  }

  mz_rat <- round(apply(sapply(1:length(out), function(i) {
    flt <- rank(-out[[i]][, idx_mz1]) <= 10 & out[[i]][, idx_mz1] > 0
    # this fails quite often
    apply(out[[i]][flt, , drop = F] / (out[[i]][flt, idx_mz1]), 2, median)
  }), 1, median), 4)

  # get mz intensity levels by scaling mz2 with median ratios
  mn_int <- round(i_mz1[bpm] * mz_rat)

  # combine these infos into dataframe
  tmp <- cbind(mz2, mn_int, d_rt, cor_rat)
  rownames(tmp) <- 1:nrow(tmp)

  # filter for mz2 belonging to base peak (=mz1) and return pseudo spectrum
  flt <- abs(d_rt) < allowed_d_rt & cor_rat > min_cor_rat
  # browser()
  if (any(flt, na.rm = T)) {
    spec <- matrix(c(mz2[which(flt)], mn_int[which(flt)]), ncol = 2, dimnames = list(NULL, c("mz", "int")))
    spec <- spec[is.finite(spec[, 2]), , drop = FALSE]
    spec <- spec[spec[, 2] > 0, , drop = FALSE]
  } else {
    spec <- matrix(NA, ncol = 2, dimnames = list(NULL, c("mz", "int")))[-1, ]
  }

  attr(spec, "rt") <- rt
  if (!is.null(msg_warn)) attr(spec, "warning") <- msg_warn
  return(spec)
}
