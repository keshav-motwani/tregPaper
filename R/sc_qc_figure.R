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
sc_qc_figure = function(sce_list,
                                 main_result_path,
                                 supplemental_result_path,
                                 cache_path,
                                 result_name) {
  ###### DIFFUSION MAP COLORED BY CLUSTERS
  pca = plot_reduced_dimensions(
    sce_list,
    type = "PCA",
    features = c("pct_mito", "num_genes_expr", "total_umi"),
    facet_rows = c(".feature"), facet_columns = "sample",
    scales = "fixed",
    # facet_switch = "y",
    facet_type = "grid",
    alpha = 0.5
  ) +
    labs(x = "Principle Component 1",
         y = "Principle Component 2")

  umap = plot_reduced_dimensions(
    sce_list,
    type = "UMAP",
    features = c("pct_mito", "num_genes_expr", "total_umi"),
    facet_rows = c(".feature"), facet_columns = "sample",
    scales = "fixed",
    facet_type = "grid",
    # facet_switch = "y",
    alpha = 0.5
  ) +
    labs(x = "UMAP Component 1",
         y = "UMAP Component 2")

  wrap_plots(
    pca,
    umap,
    ncol = 2,
    byrow = FALSE
  ) + plot_layout(guides = 'collect')

  ggsave(
    file.path(
      supplemental_result_path,
      paste0(result_name, "_qc_figure.pdf")
    ),
    height = 8,
    width = 12
  )
}
