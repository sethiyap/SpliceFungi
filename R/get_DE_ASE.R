#' get_DE_ASE
#'
#' @param se_filtered_object summarised experiment object obtained from get_ASE function
#' @param control character, name of the control sample
#' @param treatment character, name of treatment sample
#' @param PSI_cutoff double, cutoff for delta_PSI, default: 0.05
#' @param LogFC_cutoff double, cutoff for Log2 fold-change, default: 1
#' @param FDR_cutoff double, cutoff on False Discovery Rate(FDR), default: 0.05
#' @param xlim double, X-axis limit for volcano plot
#' @param plot_volcano logical, to print volcano plot
#'
#' @return
#' edgeR output
#' @export
#'
#' @import ggplot2
#' @import SpliceWiz
#' @import dplyr
#' @import magrittr
#'
#'
get_DE_ASE <- function(se_filtered_object, control, treatment, PSI_cutoff=0.05,
                       LogFC_cutoff=1, FDR_cutoff=0.05, xlim=6.5, plot_volcano=TRUE){

  res_edgeR <- SpliceWiz::ASE_edgeR(
    se = se_filtered_object,
    test_factor = "condition",
    test_nom = treatment,
    test_denom = control,
    IRmode = "all"
  ) %>%
    dplyr::mutate(class=dplyr::if_else(logFC >=LogFC_cutoff, "up", dplyr::if_else(logFC<= -LogFC_cutoff, "down", "no-change"))) %>%
    dplyr::mutate(significance=dplyr::if_else(FDR< FDR_cutoff & abs_deltaPSI>= PSI_cutoff, "S", "NS")) %>%
    dplyr::arrange(significance)

  gg_theme <- function(gp){

    gp+
      ggplot2::scale_color_manual(values=c("grey", "red")) +
      ggplot2::labs(title = paste("Differential analysis", treatment," Vs ", control),
                    x = "Log2-fold change",
                    y = "FDR (-log10)")+
      ggplot2::theme_bw()+
      ggplot2::theme(axis.text = ggplot2::element_text(size=12, color="black"))
  }

  # volcano plot of log2FC Vs FDR

  if(plot_volcano==TRUE){
    g1 <- ggplot2::ggplot(res_edgeR,ggplot2::aes(x = logFC, y = -log10(FDR), color=significance)) +
      ggplot2::geom_point()+
      ggplot2::xlim(-xlim, xlim)

    print(gg_theme(g1))
    # volcano plot of log2FC Vs FDR for each AS event

    g2 <- ggplot2::ggplot(res_edgeR,ggplot2::aes(x = logFC, y = -log10(FDR), color=significance)) +
      ggplot2::geom_point()  +
      ggplot2::facet_wrap(ggplot2::vars(EventType), nrow=2)

    print(gg_theme(g2))
  }else{
    print("Only PSI plot")
  }


 #PSI plot

 t_PSI <- paste("AvgPSI_", treatment, sep="")
 c_PSI <- paste("AvgPSI_", control, sep="")

 g3 <- ggplot2::ggplot(res_edgeR, ggplot2::aes(x = 100*get(t_PSI),
                                               y = 100*get(c_PSI))) +
                  ggplot2::geom_point(ggplot2::aes(color=significance))+ ggplot2::theme_bw()+
                  ggplot2::xlim(0, 100) + ggplot2::ylim(0, 100) +
                  ggpubr::stat_cor(method = "spearman", label.x = 5, label.y = 75, show.legend = FALSE, na.rm = TRUE)+
                  ggplot2::scale_color_manual(values=c("grey", "red"))+
                  ggplot2::labs(title = "PSI values across conditions",
                               x = paste("PSI of condition",treatment)  , y = paste("PSI of condition",control))
 print(g3)
 return(res_edgeR)
}

