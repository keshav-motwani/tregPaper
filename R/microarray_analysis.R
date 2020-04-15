#' Run microarray analysis and generate figures
#'
#' @param microarray_data
#' @param features_to_annotate
#' @param annotation_p_value_rank_threshold
#' @param result_name
#' @param main_result_path
#' @param supplemental_result_path
#' @param data_path
#'
#' @return
#' @export
#'
#' @examples
microarray_analysis = function(microarray_data,
                               features_to_annotate,
                               annotation_p_value_rank_threshold = 30,
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

  microarray_data = microarray_process_data(microarray_data)

  de = microarray_de(
    microarray_data,
    attributes(microarray_data)$sample_annotations,
    c("origin", "cell_type"),
    "donor"
  )

  readr::write_csv(de,
                   file.path(data_path, "microarray_de.csv"))

  de %>%
    dplyr::group_by(contrast) %>%
    dplyr::filter(abs(logFC) > 1.5) %>%
    dplyr::top_n(-100, p_val_adj) %>%
    readr::write_csv(file.path(supplemental_result_path, "microarray_top_100_de.csv"))

  microarray_figure(
    microarray_data,
    de,
    features_to_annotate,
    annotation_p_value_rank_threshold,
    main_result_path,
    supplemental_result_path,
    data_path
  )
}
