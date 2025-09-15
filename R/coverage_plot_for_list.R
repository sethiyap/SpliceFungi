
#' coverage_plot_for_list
#'
#' @param event_list vector, of events to be plotted
#' @param res_edgeR output of get_DE_ASE
#' @param se_object output of get_ASE
#' @param control vector, name of control sample
#' @param treatment vector, name of treatment sample
#'
#' @return a coverage plot
#' @export
#' @import tidyr
#' @import dplyr
#' @import ggplot2
#' @import SpliceWiz
#' @import purrr
#' @import ggpubr
#'
#' @examples
#'
#' event_list=c("CNAG_00332/AFR92466_Intron8/clean","CNAG_01940/AFR98135_Intron3/clean", "CNAG_02944/AFR93748_Intron1/clean" )
#' coverage_plot_for_list(event_list = event_list, res_edgeR = res_edgeR_bs181_5h, se_object = se_bs181_5h, control = "WT_5h",treatment = "WT_BS181_5h")
#'
coverage_plot_for_list <- function(event_list, res_edgeR, se_object, control,treatment){

  el <- tidyr::as_tibble(event_list)

res_edgeR.filtered <- res_edgeR %>%
    dplyr::right_join(y=el, by=c("EventName"="value")) %>%
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
