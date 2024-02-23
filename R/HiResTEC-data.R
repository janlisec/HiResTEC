#' @name sam
#' @title Sample table
#' @description This data frame contains the sample definition of 18 samples from a larger experiment.
#' @docType data
#' @usage data(sam)
#' @source jan.lisec@@bam.de
#' @keywords data
"sam"

#' @name xcms_cand
#' @title Dataframe with putative candidates
#' @description This data frame contains the analysis result of an xcmsSet which can not be provided via CRAN anymore using \link{EvaluatePairsFromXCMSSet} with respect to interesting m/z-pairs.
#' @docType data
#' @usage data(xcms_cand)
#' @source jan.lisec@@bam.de
#' @keywords data
"xcms_cand"

#' @name res_list
#' @title The main results object of a non-targeted search for tracer incorporation.
#' @description This is a list containing the evaluations results established based on processing example data with \link{EvaluateCandidateListAgainstRawData}.
#' @docType data
#' @usage data(res_list)
#' @source jan.lisec@@bam.de
#' @keywords data
"res_list"

#' @name mz_shift_corrector
#' @title Predefined mass search windows to be used internally.
#' @description This is a list defining windows for high res APCI or ESI instrumentation..
#' @docType data
#' @usage data(mz_shift_corrector)
#' @source jan.lisec@@bam.de
#' @keywords data
"mz_shift_corrector"
