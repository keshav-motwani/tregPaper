#' Generate cell cycle figure
#'
#' @param sce_list list of SCE objects
#' @param de_data data with output from sc_de
#' @param cluster_resolution name of cluster variable
#' @param contaminant_clusters ids of identified contaminant clusters
#'
#' @import cowplot
#' @import patchwork
#' @import ggplot2
#' @importFrom ggexp plot_barplot
#' @importFrom readr write_csv
#' @importFrom dplyr filter group_by arrange top_n pull
#' @importFrom scanalysis plot_pairwise_features
#'
#' @return
#' @export
#'
#' @examples
sc_cell_cycle_figure = function(sce_list,
                                         main_result_path,
                                         supplemental_result_path,
                                         cache_path,
                                         result_name) {

  ###### DIFFUSION MAP COLORED BY CLUSTERS
  clusters_pca = plot_reduced_dimensions(
    sce_list,
    type = "PCA",
    features = "Phase",
    facet_rows = c("sample"),
    scales = "fixed",
    facet_type = "wrap",
    nrow = 1,
    alpha = 0.5
  ) +
    theme(legend.position = "none") +
    guides(color = guide_legend(ncol = 1, override.aes = aes(alpha = 1, size = 3))) +
    labs(x = "Principle Component 1",
         y = "Principle Component 2") +
    ggsci::scale_color_startrek()


  clusters_diffusion_map = plot_reduced_dimensions(
    sce_list,
    type = "UMAP",
    features = "Phase",
    facet_rows = c("sample"),
    scales = "fixed",
    facet_type = "wrap",
    nrow = 1,
    alpha = 0.5
  ) +
    guides(color = guide_legend(ncol = 1, override.aes = aes(alpha = 1, size = 3))) +
    labs(x = "UMAP Component 1",
         y = "UMAP Component 2") +
    ggsci::scale_color_startrek()

  phase = encode_cell_identity_frequency_long(
    sce_list,
    attributes = "Phase",
    group_by = c("sample"),
    normalize = "none"
  )

  phase = plot_barplot(
    phase,
    x = "sample",
    y = "value",
    color = "Phase",
    label = "value",
    facet_columns = "sample",
    scales = "free"
  ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(y = NULL, x = NULL) +
    theme(
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.text.x = element_blank(),
      legend.position = "none"
    ) +
    ggsci::scale_color_startrek()

  wrap_plots(clusters_pca, clusters_diffusion_map, phase, ncol = 3, widths = c(3, 3, 1.2), byrow = FALSE) + plot_layout(guides = 'collect')
  ggsave(file.path(supplemental_result_path, paste0(result_name, "_cell_cycle_figure.pdf")), height = 4, width = 18)
}
