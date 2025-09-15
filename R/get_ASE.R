#' get_ASE
#'
#' @param reference_path character, path to the ensembl reference folder
#' @param processed_bam_path character, path to the output of the bam processed by spliceWiz
#' @param control character, name of the control sample
#' @param treatment character, name of treatment sample
#' @param replicates numeric, number of replicates
#'
#' @return
#' An object of SummarisedExperiment
#' @export
#' @import magrittr
#' @import SpliceWiz
#' @import SummarizedExperiment
#'
#'
get_ASE<- function(reference_path, processed_bam_path, control, treatment, replicates=3){

  library("magrittr")
  ref_path <- file.path(reference_path)

  pb_path=file.path(processed_bam_path)

  nxtse_path <- file.path(pb_path, "NxtSE_output")

  expr <- SpliceWiz::findSpliceWizOutput(pb_path)

  SpliceWiz::collateData(
        Experiment = expr,
        reference_path = ref_path,
        novelSplicing = TRUE,
        output_path = nxtse_path,

        novelSplicing_requireOneAnnotatedSJ = TRUE,
        # novel junctions must share one annotated splice site

        novelSplicing_minSamples = 3,
        # retain junctions observed in 3+ samples (of any non-zero expression)

        novelSplicing_minSamplesAboveThreshold = 1,
        # only 1 sample required if its junction count exceeds a set threshold
        novelSplicing_countThreshold = 10  ,
        # threshold for previous parameter

        novelSplicing_useTJ = TRUE
        # whether tandem junction reads should be used to identify novel exons
  )

  se <- SpliceWiz::makeSE(nxtse_path, realize = TRUE)

  SummarizedExperiment::colData(se)$condition <- rep(c(control, treatment), each = replicates)
  SummarizedExperiment::colData(se)$batch <- rep(paste0("set",seq(replicates)), 2)

  # se.filtered <- se[SpliceWiz::applyFilters(se),]

  kk <- SummarizedExperiment::colData(se)

  k <- knitr::kable(kk)

  print(k)

  return(se)
}

