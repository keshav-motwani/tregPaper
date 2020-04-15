#' Generate flow cytometry figure
#'
#' @param flow_data
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
flow_figure = function(flow_data,
                       de,
                       features_to_plot = unique(flow_data$feature),
                       result_name,
                       main_result_path = NULL,
                       supplemental_result_path = NULL,
                       data_path = NULL) {
  pairwise_annotation = ggexp:::.prepare_pairwise_annotation(flow_data, de, "group", "value", "p_signif", "ns", groups = "feature")

  plot_data = flow_data %>%
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
                        open = "",
                        close = "",
                        tag_pool = LETTERS) +
    ggexp::theme_ggexp()

  plot = plot + ggplot2::theme(legend.position = "none") + ggplot2::labs(x = NULL, y = "Percent Positive")

  ggplot2::ggsave(
    file.path(main_result_path, "figure_5.jpeg"),
    plot,
    height = 9,
    width = 14
  )

  plot = ggexp::plot_distributions(
    data = flow_data,
    pairwise_annotation = pairwise_annotation,
    type = "quasirandom",
    x = "group",
    y = "value",
    color = "group",
    pairwise_annotation_tier_width = 0.1,
    facet_rows = "feature",
    scales = "free",
    facet_type = "wrap",
    nrow = NULL,
    pairwise_annotation_label = "p_signif"
  )

  plot = egg::tag_facet(
    plot,
    open = "",
    close = "",
    tag_pool = c(LETTERS, paste0(LETTERS, letters))
  ) +
    ggexp::theme_ggexp()

  plot = plot + ggplot2::theme(legend.position = "none") + ggplot2::labs(x = NULL, y = "Percent Positive")

  ggplot2::ggsave(
    file.path(supplemental_result_path, "flow_figure_all.jpeg"),
    plot,
    height = 20,
    width = 20
  )
}
