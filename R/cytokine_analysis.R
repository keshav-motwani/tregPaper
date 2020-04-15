#' Run cytokine analysis and generate figures
#'
#' @param cytokine_data
#' @param features_to_plot
#' @param result_name
#' @param main_result_path
#' @param supplemental_result_path
#' @param data_path
#'
#' @return
#' @export
#'
#' @examples
cytokine_analysis = function(cytokine_data,
                             features_to_plot,
                             result_name,
                             main_result_path,
                             supplemental_result_path,
                             data_path) {
  main_result_path = file.path(main_result_path, result_name)
  dir.create(main_result_path, recursive = TRUE)

  supplemental_result_path = file.path(supplemental_result_path, result_name)
  dir.create(supplemental_result_path, recursive = TRUE)

  data_path = file.path(data_path, result_name)
  dir.create(data_path, recursive = TRUE)

  cytokine_data = flow_cytokine_process_data(cytokine_data) %>%
    dplyr::mutate(group = factor(
      group,
      levels = c("CB_Treg", "APB_Treg", "CB_Tconv", "APB_Tconv")
    ),
    value = log10(value))

  de = flow_cytokine_de(cytokine_data,
                        "group",
                        "donor") %>%
    dplyr::mutate(group1 = level1, group2 = level2) %>%
    rstatix::add_significance("p.value", "p_signif")

  readr::write_csv(de, file.path(data_path, "cytokine_de.csv"))

  cytokine_figure(
    cytokine_data,
    de,
    features_to_plot,
    result_name,
    main_result_path,
    supplemental_result_path,
    data_path
  )
}
