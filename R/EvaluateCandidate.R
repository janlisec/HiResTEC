#' @title EvaluateCandidate.
#'
#' @description
#' \code{EvaluateCandidate} will analyze a mass pair (mz1, mz2) for changes due to 13C incorporation.
#'
#' @details
#' This function will extract relevant information regarding a potential candidate mz-pair (for enrichment) from raw data files.
#'
#' @param x A single row from results dataframe provided by \link{EvaluatePairsFromXCMSSet}.
#' @param tp Timepoint.
#' @param gr Group.
#' @param dat list of xcms raws for deconvolution and plotting.
#' @param mz_iso The isotopic mass shift under investigation (e.g. 1.003355 for 13C experiments)
#' @param dmz Allowed mz deviation in milli Dalton for a BPCs (passed on to getMultipleBPC, can be a matrix).
#' @param drt Allowed rt deviation in seconds for BPCs.
#' @param dEcut Minimum required change in enrichment before a candidate ID is assigned.
#' @param Pcut Maximum allowed P value before a candidate ID is assigned.
#' @param Icut Minimum required median peak intensity (mz1 and mz2) before a candidate ID is assigned.
#' @param flux_lib Flux Library.
#' @param flux_lib_masses Masses from library spectra as a list.
#' @param flux_lib_rt_dev Allow this time deviation to call a metabolite identified.
#' @param method Either APCI or ESI. Choice will modify some internal parameters and checks performed.
#'
#' @return A list with all important information including a deconvoluted
#'     spectrum (from raw data) and linear model analysis.
#'
#' @keywords internal
#' @noRd
EvaluateCandidate <- function(x = NULL, tp = NULL, gr = NULL, dat = NULL, mz_iso = 1.003355, dmz = 0.005, drt = 0.5, dEcut = 0, Pcut = 0.01, Icut = 1000, flux_lib = NULL, flux_lib_masses = NULL, flux_lib_rt_dev = 5, method = c("APCI", "ESI")[1]) {
  # potential parameters
  # assign error message if P_raw exceeds this value

  # measurement precision of the device im mDa (should be a worst case estimate)
  mz_precision <- 4
  if (method == "ESI") mz_precision <- 1000 * unique(dmz)[1]

  res <- vector("list")
  test <- TRUE # initialize test variable
  res[["mz1"]] <- mz1 <- x[, "mz1"]
  res[["mz2"]] <- x[, "mz2"]
  res[["mz_iso"]] <- mz_iso
  res[["rt"]] <- rt <- x[, "rt"]
  res[["ng"]] <- ng <- round(diff(unlist(x[, c("mz1", "mz2")])))
  res[["tp"]] <- tp
  res[["gr"]] <- gr
  if (!is.null(gr)) res[["inter"]] <- interaction(gr, tp, sep = "_", drop = TRUE) else res[["inter"]] <- factor(tp)
  res[["bpc"]] <- NULL
  res[["err_msg"]] <- NULL
  res[["lm"]] <- NA
  res[["P_raw"]] <- 1
  res[["dE"]] <- -100

  # check if cand is present in FluxLib
  if (is.null(flux_lib)) {
    tmp_txt <- ""
  } else {
    rtok <- which(abs(flux_lib[, "RT"] - rt) < flux_lib_rt_dev)
    if (any(rtok)) {
      mzok <- which(sapply(flux_lib_masses[rtok], function(flm) {
        any(abs(flm - mz1) < 0.01)
      }))
      if (any(mzok)) {
        tmp_txt <- paste(flux_lib[rtok[mzok], "Name"], collapse = "; ")
      } else {
        tmp_txt <- "No flux_lib peak in RT range has mz overlap in Spectrum"
      }
    } else {
      tmp_txt <- "No flux_lib peak in RT range"
    }
  }
  res[["FluxLib"]] <- tmp_txt

  # evaluate BPCs
  res[["bpc"]] <- lapply(dat, function(x) {
    getMultipleBPC(x = x, mz_dev = dmz, mz = mz1 + (0:ng) * mz_iso, rt = rt, rt_dev = 1 + 4 * drt, smooth = 0)
  })
  # ion case of empty samples list elements may be NULL and have to be filled with NAs instead
  #browser()
  if (any(sapply(res[["bpc"]], is.null))) {
    emp <- which(sapply(res[["bpc"]], is.null))
    def <- min(which(!(1:length(res[["bpc"]]) %in% emp)))
    def <- res[["bpc"]][[def]]
    def[is.finite(def)] <- NA
    for (k in emp) res[["bpc"]][[k]] <- def
  }
  res[["enr"]] <- t(ldply_base(res[["bpc"]], function(x) {
    x[attr(x, "maxBPC"), , drop = FALSE]
  }))
  attr(res[["enr"]], "Enrichment") <- sapply(res[["bpc"]], function(x) {
    CalcEnrichment(unlist(x[attr(x, "maxBPC"), ]))
  })
  res[["dE"]] <- round(median(attr(res[["enr"]], "Enrichment")[tp == max(tp)], na.rm = T) - median(attr(res[["enr"]], "Enrichment")[tp == min(tp)], na.rm = T), 2)
  ratios <- apply(res[["enr"]][c(1, ng + 1), ], 2, CalcEnrichment)

  # check intensity cutoff on M+0 and M+n
  if (sum(res[["enr"]][1, tp == min(tp)] > Icut, na.rm = T) < 0.5 * sum(tp == min(tp)) & sum(res[["enr"]][nrow(res[["enr"]]), tp == max(tp)] > Icut, na.rm = T) < 0.5 * sum(tp == max(tp))) {
    res[["err_msg"]] <- c(res[["err_msg"]], paste0("mz1 or mz2 < Icut (", Icut, ")"))
  }

  # check dEcut
  if (is.na(res[["dE"]]) || res[["dE"]] < dEcut) {
    res[["err_msg"]] <- c(res[["err_msg"]], paste0("dE < dEcut (", dEcut, ")"))
  }

  # compute the average mass drift from mz1 to mz2
  median_mass_shift <- apply(sapply(res[["bpc"]], attr, "mass_defect"), 1, function(x) {
    sapply(split(x, tp), median, na.rm = T)
  })
  mass_drift <- apply(median_mass_shift, 2, function(y) {
    y[length(y)] - y[1]
  })
  # NA values have to be avoided during below testing
  mass_drift[!is.finite(mass_drift)] <- 0
  median_mass_shift[!is.finite(median_mass_shift)] <- 0
  # for APCI (silylated compounds this is expected to be strongly positive (i.e. >4) in case of labelling for mz=M+2 or higher)
  if (method == "APCI") {
    # masses M+0 and M+1 should not change over time, otherwise a coeluting compound probably caused that
    mass_drift_test <- !all(abs(mass_drift)[1:2] < mz_precision, na.rm = T)
    # all masses should stay within precision (although deviation to lower end are allowed because of Si)
    mass_drift_test <- mass_drift_test | !all(median_mass_shift < mz_precision)
    ratio_change <- diff(range(sapply(split(ratios, tp), median, na.rm = T), na.rm = T))
    if (!is.finite(ratio_change)) ratio_change <- 0
    nm <- length(mass_drift)
    if (nm >= 3) {
      mass_drift_test <- mass_drift_test | !(median_mass_shift[1, nm] <= 0 & median_mass_shift[nrow(median_mass_shift), nm] <= mz_precision)
      # [JL] if sampling starts at t>0 assumption of negative mass shift may be violated in many cases
      if (min(tp) == 0 | ratio_change > 10) {
        # [JL] modified test. The upper line was often too harsh in filtering
        # mass_drift_test <- mass_drift_test | mass_drift[nm]<4 | median_mass_shift[nrow(median_mass_shift),nm]<(-mz_precision)
        mass_drift_test <- mass_drift_test | mass_drift[nm] < ifelse(is.finite(res[["dE"]]) && res[["dE"]] < 10, 0, 3)
      }
    }
  }

  if (method == "test") {
    mass_drift_test <- FALSE
  }

  if (method == "ESI") {
    # check all mz which are >10% base peak intensity in at least 50% of all samples
    flt <- apply(apply(res[["enr"]], 2, function(y) {
      if (all(is.na(y))) {
        return(y)
      } else {
        return(y / max(y, na.rm = T))
      }
    }), 1, function(z) {
      sum(z > 0.1, na.rm = T) >= 0.5 * length(z)
    })
    # check those mz if m/z deviates more than mz_precision allowed
    mass_drift_test <- !all(abs(apply(sapply(res[["bpc"]], attr, "mass_defect")[flt, , drop = FALSE], 1, median, na.rm = T)) <= mz_precision, na.rm = T)
    # check those mz if m/z drift occurs (should be potentially only done for mass_drift[flt])
    mass_drift_test <- mass_drift_test | !all(mass_drift < mz_precision)
  }

  if (mass_drift_test) {
    res[["err_msg"]] <- c(res[["err_msg"]], "unexpected mass drift")
  }

  # compute mean cv ratios (only if >=3 replicates per group are present)
  if (any(table(res[["inter"]]) >= 3)) {
    mean_sd <- mean(sapply(split(ratios, res[["inter"]]), function(x) {
      ifelse(sum(is.finite(x)) >= 3, sd(x, na.rm = T), NA)
    }), na.rm = T)
    mean_sd <- mean_sd / diff(range(ratios, na.rm = T))
    if (!is.finite(mean_sd) || mean_sd > 0.3) {
      res[["err_msg"]] <- c(res[["err_msg"]], "mean_sd > 0.3 of total range")
    }
  }

  # compute a linear model for pseudo enrichments (==ratios)
  tp <- as.numeric(as.character(tp))
  if (length(levels(gr)) > 1) {
    # if groups are specified use:
    ratios_lm <- try(stats::lm(ratios ~ tp * gr), silent = TRUE)
  } else {
    # if only timepoints are specified use:
    ratios_lm <- try(stats::lm(ratios ~ tp), silent = TRUE)
  }
  if (!attr(ratios_lm, "class") == "lm") {
    res[["P_raw"]] <- 1
    c(res[["err_msg"]], "lm on pseudo enrichments failed")
  } else {
    res[["lm"]] <- ratios_lm
    # strip the environment attributes as object blows up in size otherwise
    attr(res[["lm"]]$terms, ".Environment") <- NULL
    attr(attr(res[["lm"]]$model, "terms"), ".Environment") <- NULL
    res[["P_raw"]] <- 1
    alm <- stats::anova(ratios_lm)
    alm_tp <- grep("tp", rownames(alm))
    if (length(alm_tp) >= 1 & any(is.finite(alm[alm_tp, "Pr(>F)"]))) {
      res[["P_raw"]] <- min(alm[alm_tp, "Pr(>F)"], na.rm = T)
    } else {
      # no time dependent P-value could be found using a linear model
      res[["P_raw"]] <- 1
      c(res[["err_msg"]], "lm did not yield a time dependent P-value")
    }
  }

  # check Pcut
  if (res[["P_raw"]] > Pcut) {
    res[["err_msg"]] <- c(res[["err_msg"]], paste0("P > Pcut (", Pcut, ")"))
  }

  return(res)
}
