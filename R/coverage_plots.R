#' coverage_plots
#'
#' @param res_edgeR output of get_DE_ASE
#' @param se_object output of get_ASE
#' @param n_events number of significant (FDR <0.05, deltaPSI>=0.05) events to be plotted
#'
#' @return
#' a ggplot of coverage of ASE
#' @export
#' @import ggplot2
#' @import purrr
#' @import dplyr
#' @import SpliceWiz
#' @import ggpubr
#' @import magrittr
#'
coverage_plots <- function(res_edgeR, se_object, n_events=3, control,treatment){

  res_edgeR.filtered <- res_edgeR %>%
    dplyr::filter(significance=="S") %>%
    dplyr::arrange(desc(abs_deltaPSI)) %>%
    dplyr::slice(1:n_events) %>%
    dplyr::select(c(EventName, EventType))

  print(res_edgeR.filtered$EventName)

  gg <- res_edgeR.filtered %>%
            dplyr::mutate(cov_plot= purrr::map2(EventName,EventType, function(i,j){

              dataObj <- SpliceWiz::getCoverageData(se_object, Event = i, tracks = colnames(se_object))

              plotObj <- SpliceWiz::getPlotObject(dataObj,Event = i,condition = "condition", tracks = c(control,treatment))

              nn <- SpliceWiz::plotView(plotObj, centerByEvent = TRUE,
                                        trackList = list(c(1,2)), filterByEventTranscripts = TRUE)+ggplot2::ggtitle(j)

              return(nn)

  }))

  x = nrow(gg)

  cols = round(sqrt(x),0)
  rows = ceiling(x/cols)

  g1 <- ggpubr::ggarrange(plotlist = gg$cov_plot, ncol=cols, nrow = rows)

  print(g1)

}


