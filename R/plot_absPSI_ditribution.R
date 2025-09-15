#' plot_absPSI_ditribution
#'
#' @param res_edgeR res_edgeR object
#' @param plot_type character, "densityplot" or "boxplot"
#'
#' @return
#' a ggplot
#' @export
#' @import ggplot2
#' @import dplyr
#' @import magrittr
#'
plot_absPSI_ditribution <- function(res_edgeR, plot_type="boxplot"){

  plot_data <- res_edgeR %>%
    dplyr::filter(significance=="S") %>% dplyr::select(c(EventName, EventType, abs_deltaPSI, class))

  if(plot_type=="boxplot"){

    g1 <- plot_data %>%
      ggplot2::ggplot(ggplot2::aes(class, abs_deltaPSI*100, fill=class))+
      ggplot2::geom_boxplot()+
      ggplot2::labs(title="distribution of abs(∆PSI)", x="", y="% of abs(∆PSI)")+
      ggplot2::scale_y_log10()

  }
  if(plot_type=="densityplot"){

    g1 <- plot_data %>%
      ggplot2::ggplot(ggplot2::aes(abs_deltaPSI*100, fill=class))+
      ggplot2::geom_density()+
      ggplot2::labs(title="distribution of abs(∆PSI)", y="density", x="% of abs(∆PSI)")
  }

  else{
    message("select either boxplot or densityplot")
  }

gg_theme <- function(gg){

  gg+ggplot2::facet_wrap(~EventType)+ggplot2::theme_bw()+
    ggplot2::theme(axis.text = ggplot2::element_text(size=12, color="black"))+
    ggplot2::scale_fill_manual(values = c("#00AFBB","#FC4E07" ))
}

gg_theme(g1)

}
