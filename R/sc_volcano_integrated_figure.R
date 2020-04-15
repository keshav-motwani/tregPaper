#' Generate cell cycle figure
#'
#' @param de_data data with output from sc_de_cluster
#' @param cluster_resolution name of cluster variable
#' @param contaminant_clusters ids of identified contaminant clusters
#'
#' @import patchwork
#' @importFrom cowplot plot_grid
#' @import ggplot2
#' @importFrom dplyr filter
#' @importFrom scanalysis plot_volcano
#'
#' @return
#' @export
#'
#' @examples
sc_volcano_integrated_figure = function(de_data,
                                        main_result_path,
                                        supplemental_result_path,
                                        cache_path) {
  apbFC = grep("logFC", grep("APB", colnames(de_data), value = TRUE), value = TRUE)
  cbFC = grep("logFC", grep("CB", colnames(de_data), value = TRUE), value = TRUE)

  apbP = grep("padj", grep("APB", colnames(de_data), value = TRUE), value = TRUE)
  cbP = grep("padj", grep("CB", colnames(de_data), value = TRUE), value = TRUE)

  apb = plot_volcano(
    de_data,
    annotations = c("FOXP3", "IKZF2"),
    fold_change = apbFC,
    p_value = apbP,
    label = "feature",
    facet_rows = "group",
    nrow = 5,
    facet_type = "wrap"
  ) +
    labs(subtitle = grep("APB", colnames(de_data), value = TRUE))
  cb = plot_volcano(
    de_data,
    annotations = c("FOXP3", "IKZF2"),
    fold_change = cbFC,
    p_value = cbP,
    label = "feature",
    facet_rows = "group",
    nrow = 5,
    facet_type = "wrap"
  ) +
    labs(subtitle = labs(title = grep("CB", colnames(de_data), value = TRUE)))

  assembled_plot = plot_grid(apb, cb, ncol = 2)

  ggsave(
    file.path(supplemental_result_path, "volcano_integrated_figure.pdf"),
    assembled_plot,
    height = 20,
    width = 20
  )
}
