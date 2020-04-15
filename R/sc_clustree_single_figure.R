#' Generate cell cycle figure
#'
#' @param sce_list list of SCE objects
#' @param de_data data with output from sc_de
#' @param cluster_resolution name of cluster variable
#' @param contaminant_clusters ids of identified contaminant clusters
#'
#' @import patchwork
#' @import ggplot2
#' @import ggraph
#' @importFrom clustree clustree
#' @importFrom dplyr filter
#' @importFrom scanalysis plot_volcano
#'
#' @return
#' @export
#'
#' @examples
sc_clustree_single_figure = function(sce_list,
                                              main_result_path,
                                              supplemental_result_path,
                                              cache_path) {
  samples = names(sce_list)

  plots = list()

  max_foxp3 = c()

  for (sample in samples) {
    plot = clustree(
      sce_list[[sample]],
      "cluster_",
      node_colour = "FOXP3",
      node_colour_aggr = "mean",
      exprs = "logcounts"
    ) + labs(title = sample)
    max_foxp3 = c(max_foxp3, max(plot$data$mean_logcounts_FOXP3))
    plots = c(plots, list(plot))
  }

  for (plot in 1:length(plots)) {
    plots[[plot]] = plots[[plot]] + scale_color_gradient(low = "white",
                                       high = "firebrick",
                                       limits = c(0, max(max_foxp3)))
  }

  assembled_plot = wrap_plots(plots, ncol = 2)

  ggsave(
    file.path(supplemental_result_path, "clustree_single_figure.pdf"),
    assembled_plot,
    height = 10,
    width = 25
  )
}
