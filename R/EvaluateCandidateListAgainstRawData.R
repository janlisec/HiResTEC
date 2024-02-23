#' @title EvaluateCandidateListAgainstRawData.
#'
#' @description
#' \code{EvaluateCandidateListAgainstRawData} will analyze an xcmsSet result for mass pairs (mz1, mz2) with changes due to 13C incorporation.
#'
#' @details
#' This function will evaluate candidate mz pairs found within an xcmsSet object by \link{EvaluatePairsFromXCMSSet} against the raw measurement data.
#' A special parameter is 'rolp' which can be set to 'non', 'pos', 'neg' or 'all'. It will influence the time performance of the function be determining how many peaks are effectively tested.
#' If 'rolp' is set to 'non', no overlapping peaks will be skipped, every individual mz-pair will be sequentially evaluated (slow but most informative).
#' If it is set to 'pos' or 'neg', overlapping peaks (determined by experiment wide deconvolution) will not be tested aditionally for positive or negative hits ('neg' is standard).
#' If set to 'all' overlapping peaks will always be removed from the list of mz-pairs to be tested (fast).
#'
#' @param x Dataframe of results (output of EvaluatePairsFromXCMSet).
#' @param tp Timepoint.
#' @param gr group, e.g. different genotypes or concentrations.
#' @param dat list of xcmsRaw's for deconvolution and plotting.
#' @param dmz Allowed mass deviation in Da (for BPC extraction).
#' @param drt Allowed rt deviation in seconds (for get extraction).
#' @param dEcut Minimum required change in enrichment before a candidate ID is assigned.
#' @param Pcut Maximum allowed P value before a candidate ID is assigned.
#' @param Icut Minimum required median peak intensity before a candidate ID is assigned.
#' @param method Either APCI or ESI. Choice will modify some internal parameters and checks performed.
#' @param rolp RemoveOverLappingPeaks paramter.
#' @param smooth Smoothing parameter passed to \link{getMultipleBPC}.
#'
#' @return
#' A list of evaluation results.
#'
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#'
#' @examples
#' # Please use examples from previous versions as xcms (and xcms objects)
#' # are no longer supported during CRAN checks leading to package rejection
#' # if included (and I do not know a work around). :(
#'
#' @export
#'
EvaluateCandidateListAgainstRawData <- function(x = NULL, tp = NULL, gr = NULL, dat = NULL, dmz = 0.025, drt = 1, dEcut = 1, Pcut = 0.01, Icut = 1000, method = c("APCI", "ESI")[1], rolp = c("non", "pos", "neg", "all")[2], smooth = 0) {
  # potential paramters
  # ...
  flux_lib <- NULL

  # convert spectra from lib into list of masses (without intensity information)
  if (!is.null(flux_lib)) {
    flux_lib_masses <- lapply(flux_lib[, "Spectrum"], function(x) {
      as.numeric(sapply(strsplit(x, " ")[[1]], function(y) {
        strsplit(y, ":")[[1]][1]
      }))
    })
  } else {
    flux_lib_masses <- NULL
  }

  # check if any data is provided
  if (is.null(x) || nrow(x) == 0) {
    stop("\nNo data provided.")
  }

  # establish mz_shift_corrector
  x[, "n"] <- round(x[, "mz2"] - x[, "mz1"])
  utils::data(mz_shift_corrector, package = "HiResTEC")
  mz_shift_corrector <- mz_shift_corrector[[method]][1:(max(x[, "n"]) + 1)]

  # initialize result list
  res_list <- NULL

  # specify two factors, usually a numeric one (tp) and a categorical (gr)
  if (is.null(tp)) tp <- rep(1, length(dat)) else tp <- as.numeric(tp)
  if (is.null(gr)) gr <- gl(1, length(dat)) else gr <- factor(gr)

  # order timegroup according to max of sum of {mz1,mz2} at any timepoint
  # [ToDo] generalize for different designs
  x <- RankCandidateList(x = x)

  # while there are still candidates fullfilling defined significance criteria do testing
  done_lines <- rep(F, nrow(x))
  pb <- utils::txtProgressBar(max = nrow(x), style = 3)
  for (i in 1:nrow(x)) {
    utils::setTxtProgressBar(pb, i)

    # in case of error store the last line processed in global environment 'last_row_processed'
    on.exit(if (!all(done_lines)) {
      # last_row_processed <- x[i,,drop=FALSE]
      # assign("last_row_processed", envir=.GlobalEnv)
      print("Something went wrong and this is the last line processed where the error most likely occured:")
      print(x[i, , drop = FALSE])
    })

    if (done_lines[i]) {
      # don't do anything if already marked due to Spectra occurence
    } else {
      # evaluate candidate without spectrum deconvolution
      browser()
      dmz_mat <- matrix(c(dmz + mz_shift_corrector[1:(x[i, "n"] + 1)], rep(dmz, x[i, "n"] + 1)), ncol = 2)
      res <- EvaluateCandidate(x = x[i, , drop = FALSE], tp = tp, gr = gr, dmz = dmz_mat, drt = drt, dat = dat, dEcut = dEcut, Pcut = Pcut, Icut = Icut, method = method, flux_lib = flux_lib, flux_lib_masses = flux_lib_masses)
      res <- EvaluateSpectrum(res = res, dat = dat, drt = drt, dmz = dmz, ionization = ifelse(method %in% c("APCI", "test"), "APCI", "ESI"), smooth = smooth)
      res[["row"]] <- rownames(x)[i]

      if (is.null(res[["err_msg"]]) | rolp %in% c("pos", "neg", "all")) {
        # identify overlapping_peaks
        overlapping_peaks <- abs(x[, "rt"] - x[i, "rt"]) < 1
        overlapping_peaks[overlapping_peaks] <- sapply(x[overlapping_peaks, "mz1"], function(mz) {
          any(abs(res[["s"]][, "mz"] - mz) < dmz)
        })
        overlapping_peaks[i] <- TRUE # ensure that current line is TRUE
        res[["OverlappingLine"]] <- paste(unique(x[which(overlapping_peaks), "peak_idx"]), collapse = ", ")
        res[["OverlappingMasses"]] <- paste(unique(x[which(overlapping_peaks), "mz1"]), collapse = ", ")
      }

      if (rolp == "non" | (is.null(res[["err_msg"]]) & rolp == "neg") | (!is.null(res[["err_msg"]]) & rolp == "pos")) {
        overlapping_peaks <- rep(FALSE, nrow(x))
        overlapping_peaks[i] <- TRUE # ensure that current line is TRUE
      }

      # remove overlapping peaks
      done_lines[overlapping_peaks] <- TRUE
      # x <- x[!overlapping_peaks,,drop=FALSE]

      # store evaluation result into list
      res_list[length(res_list) + 1] <- list(res)
    }
  }

  # enumerate positive candidates
  idx <- 1
  for (i in 1:length(res_list)) {
    if (is.null(res_list[[i]][["err_msg"]])) {
      res_list[[i]][["cand_id"]] <- idx
      idx <- idx + 1
    } else {
      res_list[[i]][["cand_id"]] <- ""
    }
  }

  invisible(res_list)
}
