#' splicing_statistics
#'
#' @param se_filtered_object summarised experiment object obtained from get_ASE function
#' @param res_edgeR output of get_DE_ASE
#'
#' @return a barplot
#' @export
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import SummarizedExperiment
#' @import tibble
#' @import magrittr
#'
#'
splicing_statistics <- function(se_filtered_object, res_edgeR){

  # All events Vs signficant events
  AS_tt <- SummarizedExperiment::rowData(se_filtered_object) %>%
    tibble::as_tibble() %>%
    dplyr::group_by(EventType) %>%
    dplyr::tally()

  deg_tt <- res_edgeR %>%
    dplyr::select(EventName, EventType, significance) %>%
    dplyr::filter(significance=="S") %>%
    dplyr::group_by(EventType) %>%
    dplyr::count()

  se_stats <- AS_tt %>%
    dplyr::full_join(deg_tt,by="EventType")

  colnames(se_stats) <- c("EventType", "AllEvents", "SignificantEvents")

  gg_theme <- function(gp){

    gp+
      ggplot2::geom_text(position =ggplot2::position_dodge(width = .9), vjust=-0.5)+
      ggplot2::theme_bw()+
      ggplot2::theme(axis.text = ggplot2::element_text(size=12, color="black"))
  }

  g1 <- se_stats %>%
          tidyr::gather(Class, count, -EventType) %>%
          ggplot2::ggplot(ggplot2::aes(EventType, count, fill=Class, label=count))+
          ggplot2::geom_col(position="dodge")+
          ggplot2::ggtitle("Significant events from all the events")

  print(gg_theme(g1))

  # DEG in signficant events

  DEG_stats <- res_edgeR %>%
    dplyr::select(EventName, EventType,logFC, significance, class) %>%
    dplyr::filter(significance=="S")  %>%
    dplyr::count(class)

  DEG_stats

  up_count <- res_edgeR %>%
              dplyr::filter(class=="up" & significance=="S") %>%
              dplyr::count(EventType)

  down_count <- res_edgeR %>%
                dplyr::filter(class=="down" & significance=="S") %>%
                dplyr::count(EventType)

  deg_count <- up_count %>% dplyr::full_join(y = down_count,by = "EventType")

  colnames(deg_count) <- c("EventType", "up", "down")


  g2 <- deg_count %>%
        tidyr::gather(Class, count, -EventType) %>%
        ggplot2::ggplot(ggplot2::aes(EventType, count, fill=Class, label=count))+
        ggplot2::geom_col(position="dodge")+
        ggplot2::scale_fill_manual(values=c("blue4", "red3")) +
        ggplot2::ggtitle("Differential Expressed Events")

  print(gg_theme(g2))


}
