
#' compare_ASE
#'
#' @param res_edgeR_1 res_edgeR tibble of first comparison
#' @param name1 character, name assigned to first comparison
#' @param res_edgeR_2 res_edgeR tibble of second comparison
#' @param name2 character, name assigned to second comparison
#'
#' @return
#' a ggplot
#' @export
#' @import dplyr
#' @import UpSetR
#' @import tidyr
#' @import magrittr
#'
#'
compare_ASE <- function(res_edgeR_1,name1, res_edgeR_2, name2, output_intersection=TRUE){

  # sample 1
  sample1_S_events <- res_edgeR_1 %>%
    dplyr::filter(significance=="S") %>%
    dplyr::filter(class=="up" | class=="down") %>%
    dplyr::select(c(EventName, class)) %>%
    dplyr::group_by(class) %>% dplyr::mutate(class=paste(class,name1, sep="_"))

  sample1_S_events$ID <- 1:nrow(sample1_S_events)

  lst_sample1 <- sample1_S_events %>% tidyr::spread(class, EventName) %>% dplyr::select(-ID) %>% as.list()

  lst_sample1 <- lapply(lst_sample1, function(x) x[!is.na(x)])

  # sample 2
  sample2_S_events <- res_edgeR_2 %>%
    dplyr::filter(significance=="S") %>%
    dplyr::filter(class=="up" | class=="down") %>%
    dplyr::select(c(EventName, class)) %>% dplyr::group_by(class) %>%
    dplyr::mutate(class=paste(class,name2, sep="_"))

  sample2_S_events$ID <- 1:nrow(sample2_S_events)

  lst_sample2 <-sample2_S_events %>% tidyr::spread(class, EventName) %>% dplyr::select(-ID) %>% as.list()

  lst_sample2 <- lapply(lst_sample2, function(x) x[!is.na(x)])

  lst <- append(lst_sample1, lst_sample2)

  U1 <- UpSetR::upset(UpSetR::fromList(lst), order.by = "degree", point.size = 2.5, line.size = 1, sets=names(lst),
                sets.bar.color = rep(c("#FC4E07","#00AFBB"),2),text.scale=2)



  if(output_intersection==TRUE){
    upset_data <- unlist(lst, use.names = FALSE)
    upset_data <- upset_data[ !duplicated(upset_data) ]


    out <- U1$New_data %>%
      tibble::as_tibble() %>%
      dplyr::mutate(us_elem = upset_data) %>%
      tidyr::gather(samples, is_element_present, -us_elem) %>%
      dplyr::filter(is_element_present == 1) %>%
      dplyr::group_by(us_elem) %>%
      tidyr::nest() %>%
      dplyr::mutate(set = purrr::map_chr(data, ~ ..1 %>% dplyr::pull(1) %>%
                                           stringr::str_c(collapse = ","))) %>%
      dplyr::select(-data) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(set) %>%
      tidyr::nest() %>%
      dplyr::mutate(elements = purrr::map(data , ~..1 %>% dplyr::pull(1))) %>%
      dplyr::select(-data) %>%
      dplyr::ungroup() %>%
      tidyr::unnest(cols = "elements")

      return(out)

  }else{
    return(U1)
  }

}
