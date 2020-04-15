#' Run flow cytometry analysis and generate figures
#'
#' @param flow_data
#' @param features_to_plot
#' @param result_name
#' @param main_result_path
#' @param supplemental_result_path
#' @param data_result_path
#'
#' @return
#' @export
#'
#' @examples
flow_analysis = function(flow_data,
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

  flow_data = flow_cytokine_process_data(flow_data) %>%
    dplyr::mutate(
      group = factor(
        group,
        levels = c("CB_Treg", "APB_Treg", "CB_Tconv", "APB_Tconv")
      ),
      value = value * 100
    )

  de = flow_cytokine_de(flow_data,
                        "group",
                        "donor") %>%
    dplyr::mutate(group1 = level1, group2 = level2) %>%
    rstatix::add_significance("p.value", "p_signif")

  readr::write_csv(de, file.path(data_path, "flow_de.csv"))

  flow_figure(
    flow_data,
    de,
    features_to_plot,
    result_name,
    main_result_path,
    supplemental_result_path,
    data_path
  )
}
