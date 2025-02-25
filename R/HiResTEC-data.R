#' @rdname HiResTEC-data
#' @title Different example data files to be used in the help section of `HiResTEC` functions.
#' @source jan.lisec@@bam.de
#' @importClassesFrom CorrectOverloadedPeaks xcmsRawLike
#' @name raw
#' @description `raw` contains a list of `xcmsRawLike` objects, each containing a selected mass range for a small retention time window of 18 samples defined in sam.
#' @docType data
#' @usage data(raw)
#' @keywords data
"raw"

#' @rdname HiResTEC-data
#' @name sam
#' @description `sam` provides a data frame containing the sample definition of 18 samples from a larger experiment.
#' @docType data
#' @usage data(sam)
#' @keywords data
"sam"

#' @rdname HiResTEC-data
#' @name xcms_cand
#' @description `xcms_cand` contains a data frame with the analysis result of an `xcmsSet` obtained using \link{EvaluatePairsFromXCMSSet}.
#' @docType data
#' @usage data(xcms_cand)
#' @keywords data
"xcms_cand"

#' @rdname HiResTEC-data
#' @name res_list
#' @description `res_list` is a list containing the evaluations results established based on processing example data with \link{EvaluateCandidateListAgainstRawData}.
#' @docType data
#' @usage data(res_list)
#' @keywords data
"res_list"

#' @rdname HiResTEC-data
#' @name mz_shift_corrector
#' @description `mz_shift_corrector` contains a list defining windows for high res APCI or ESI instrumentation.
#' @docType data
#' @usage data(mz_shift_corrector)
#' @keywords data
"mz_shift_corrector"
