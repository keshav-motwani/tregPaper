#' Generate ambient RNA barcode filter figure
#'
#' This figure helps to compare various commonly used filtering methods, including emptyDrops, the algorithm used in cellranger, and the knee and inflection point in the library size vs library size rank curve.
#'
#' @param sce_list list of SingleCellExperiment objects
#' @param barcode_filters result from sc_compute_barcode_filters
#' @param main_result_path path for storring main results
#' @param supplemental_result_path path for storing supplemental results
#' @param cache_path path for caching results
#'
#' @importFrom scanalysis plot_barcode_qc cache
#' @import ggplot2
#'
#' @return
#' @export
#'
#' @examples
sc_barcode_filter_figure = function(sce_list,
                                    barcode_filters,
                                    main_result_path,
                                    supplemental_result_path,
                                    cache_path) {

  plot = plot_barcode_qc(sce_list, barcode_filters, facet_columns = ".sample", facet_type = "wrap", scales = "free",  nrow = 1)

  plot = plot + theme(legend.position = "bottom") + ggsci::scale_color_aaas()
  plot$layers[[1]]$data = filter(plot$layers[[1]]$data, total_umi > 0)
  for (i in 1:4) {
    plot$layers[[i]]$data$.sample = factor(plot$layers[[i]]$data$.sample, levels = c("APB_Treg_Pre", "CB_Treg_Pre", "APB_Treg_Post", "CB_Treg_Post"))
  }

  ggsave(
    file.path(
      supplemental_result_path,
      "barcode_filters.pdf"
    ),
    plot,
    height = 4,
    width = 10
  )

  return(plot)
}