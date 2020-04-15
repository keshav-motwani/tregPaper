#' Generate cell cycle figure
#'
#' @param sce_list list of SCE objects
#' @param de_data data with output from sc_de
#' @param cluster_resolution name of cluster variable
#' @param contaminant_clusters ids of identified contaminant clusters
#'
#' @import patchwork
#' @import ggplot2
#' @importFrom dplyr filter
#' @importFrom scanalysis plot_volcano
#'
#' @return
#' @export
#'
#' @examples
sc_volcano_single_figure = function(de_data,
                                             main_result_path,
                                             supplemental_result_path,
                                             cache_path) {
  apb = de_data %>%
    filter(grepl("APB", .sample)) %>%
    plot_volcano(
      fold_change = "logFC1",
      annotations = c("FOXP3", "IKZF2"),
      p_value = "padj1",
      label = "feature",
      facet_rows = c("group", ".sample"),
      nrow = 5,
      facet_type = "wrap"
    )

  cb = de_data %>%
    filter(grepl("CB", .sample)) %>%
    plot_volcano(
      fold_change = "logFC1",
      annotations = c("FOXP3", "IKZF2"),
      p_value = "padj1",
      label = "feature",
      facet_rows = c("group", ".sample"),
      nrow = 5,
      facet_type = "wrap"
    )

  assembled_plot = cowplot::plot_grid(apb, cb, ncol = 2)

  ggsave(
    file.path(supplemental_result_path, "volcano_single_figure.pdf"),
    assembled_plot,
    height = 20,
    width = 20
  )
}
