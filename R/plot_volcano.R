
#' plot_volcano- Volcano plot for Alternate splicing Events
#'
#' @param res_edgeR a res_edgeR object, out of get_DE_ASE
#'
#' @return two ggplots displaying DEG 1. for all the events 2. Classified according to type
#'
#' @export
#'
#' @import ggplot2
#' @examples
#'
#' plot_volcano(res_edgeR_LP)
#'
plot_volcano <- function(res_edgeR){


  # All the events volcano plot

  gg_1 <- ggplot2::ggplot(res_edgeR,
                  ggplot2::aes(x = logFC, y = -log10(FDR), color=significance)) +
    ggplot2::geom_point()  +
    ggplot2::scale_color_manual(values=c("grey", "red"))+
    ggplot2::labs(title = "Differential analysis",
                  x = "Log2-fold change", y = "FDR (-log10)")+ggplot2::theme_bw()+
    ggplot2::theme(axis.text = ggplot2::element_text(size=12, color="black"))


  print(gg_1)

  gg_2 <- ggplot2::ggplot(res_edgeR,
                  ggplot2::aes(x = logFC, y = -log10(FDR), color=significance)) +
    ggplot2::geom_point()  +
    ggplot2::scale_color_manual(values=c("grey", "red"))+
    ggplot2::facet_wrap(vars(EventType), nrow=2)+
    ggplot2::labs(title = "Differential analysis on each type",
                  x = "Log2-fold change", y = "FDR (-log10)")+ggplot2::theme_bw()+
    ggplot2::theme(axis.text = element_text(size=12, color="black"))

  print(gg_2)

}
