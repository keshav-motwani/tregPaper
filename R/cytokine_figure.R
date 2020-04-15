#' Generate cytokine figure
#'
#' @param cytokine_data
#' @param de
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
cytokine_figure = function(cytokine_data,
                           de,
                           features_to_plot = unique(cytokine_data$feature),
                           result_name,
                           main_result_path = NULL,
                           supplemental_result_path = NULL,
                           data_path = NULL) {

  pairwise_annotation = ggexp:::.prepare_pairwise_annotation(cytokine_data, de, "group", "value", "p_signif", "ns", groups = "feature")

  plot_data = cytokine_data %>%
    dplyr::filter(feature %in% features_to_plot) %>%
    dplyr::mutate(feature = factor(feature, levels = features_to_plot))

  plot_pairwise_annotation = pairwise_annotation %>%
    dplyr::filter(feature %in% features_to_plot) %>%
    dplyr::mutate(feature = factor(feature, levels = features_to_plot))

  plot = ggexp::plot_distributions(
    data = plot_data,
    pairwise_annotation = plot_pairwise_annotation,
    type = "quasirandom",
    x = "group",
    y = "value",
    color = "group",
    pairwise_annotation_tier_width = 0.1,
    facet_rows = "feature",
    scales = "free",
    facet_type = "wrap",
    nrow = 3,
    pairwise_annotation_label = "p_signif"
  )

  plot = egg::tag_facet(plot,
                        open = "", close = "",
                        tag_pool = LETTERS) +
    ggexp::theme_ggexp()

  plot = plot +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(x = NULL, y = "log10(Concentration)")

  ggplot2::ggsave(file.path(main_result_path, "figure_6.jpeg"), plot, height = 9, width = 14)
}
