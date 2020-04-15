#' Generate outlier cell filter figure
#'
#' This figure shows the distribution of library sizes, number of genes expressed, percentage of mitochondrial reads, number of VDJ chains, and the relationships between these QC metrics
#'
#' @param sce_list list of SingleCellExperiment objects
#' @param outlier_cell_filters result from sc_compute_outlier_cell_filters
#' @param main_result_path path for storing main results
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
sc_outlier_cell_filter_figure = function(sce_list,
                                         outlier_cell_filters,
                                         barcode_filter_plot,
                                         main_result_path,
                                         supplemental_result_path,
                                         cache_path) {

  total_umi_vs_n_genes_expr = plot_gex_bivariate_qc(
    sce_list,
    x = "total_umi",
    y = "n_genes_expr",
    x_filters = outlier_cell_filters$total_umi_filters,
    y_filters = outlier_cell_filters$n_genes_expr_filters,
    facet_columns = ".sample",
    facet_type = "wrap",
    scales = "free",
    color = "pct_mito",
    nrow = 1
  ) + labs(x = "Total UMI", y = "Number of Genes Expressed")

  total_umi_vs_n_genes_expr$data$.sample = factor(
    total_umi_vs_n_genes_expr$data$.sample,
    levels = c(
      "APB_Treg_Pre",
      "CB_Treg_Pre",
      "APB_Treg_Post",
      "CB_Treg_Post"
    )
  )
  for (i in 3:7) {
    total_umi_vs_n_genes_expr$layers[[i]]$data$.sample = factor(
      total_umi_vs_n_genes_expr$layers[[i]]$data$.sample,
      levels = c(
        "APB_Treg_Pre",
        "CB_Treg_Pre",
        "APB_Treg_Post",
        "CB_Treg_Post"
      )
    )
  }

  total_umi_vs_pct_mito = plot_gex_bivariate_qc(
    sce_list,
    x = "total_umi",
    y = "pct_mito",
    x_filters = outlier_cell_filters$total_umi_filters,
    y_filters = outlier_cell_filters$pct_mito_filters,
    facet_columns = ".sample",
    facet_type = "wrap",
    scales = "free",
    color = "n_genes_expr",
    nrow = 1,
    y_log = FALSE
  ) + labs(x = "Total UMI", y = "Percentage Mitochondrial Reads")

  total_umi_vs_pct_mito$data$.sample = factor(
    total_umi_vs_pct_mito$data$.sample,
    levels = c(
      "APB_Treg_Pre",
      "CB_Treg_Pre",
      "APB_Treg_Post",
      "CB_Treg_Post"
    )
  )
  for (i in 3:7) {
    total_umi_vs_pct_mito$layers[[i]]$data$.sample = factor(
      total_umi_vs_pct_mito$layers[[i]]$data$.sample,
      levels = c(
        "APB_Treg_Pre",
        "CB_Treg_Pre",
        "APB_Treg_Post",
        "CB_Treg_Post"
      )
    )
  }

  total_umi_vs_vdj = plot_vdj_gex_univariate_qc(
    sce_list,
    x_filters = outlier_cell_filters$total_umi_filters,
    vdj_filters = outlier_cell_filters$vdj_filters,
    x = "total_umi",
    facet_columns = ".sample",
    point_size = 0.5,
    alpha = 1,
    scales = "free",
    facet_type = "wrap",
    nrow = 1
  ) + labs(x = "Total UMI", y = "Chain Count: TRA_TRB_IGL_IGK_IGH")

  total_umi_vs_vdj$data$.sample = factor(
    total_umi_vs_vdj$data$.sample,
    levels = c(
      "APB_Treg_Pre",
      "CB_Treg_Pre",
      "APB_Treg_Post",
      "CB_Treg_Post"
    )
  )
  for (i in 2:5) {
    total_umi_vs_vdj$layers[[i]]$data$.sample = factor(
      total_umi_vs_vdj$layers[[i]]$data$.sample,
      levels = c(
        "APB_Treg_Pre",
        "CB_Treg_Pre",
        "APB_Treg_Post",
        "CB_Treg_Post"
      )
    )
  }

  plot = cowplot::plot_grid(
    barcode_filter_plot + theme(legend.position = "right"),
    total_umi_vs_n_genes_expr,
    total_umi_vs_pct_mito,
    total_umi_vs_vdj,
    ncol = 1,
    rel_heights = c(1, 0.8, 0.8, 1.5)
  )

  ggsave(
    file.path(supplemental_result_path, "qc_figure.pdf"),
    plot,
    height = 17,
    width = 13
  )
}