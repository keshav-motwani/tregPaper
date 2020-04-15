#' Compute filters to identify outlier cells
#'
#' Computes filters based on library sizes, number of genes expressed, percentage of mitochondrial reads, and number of VDJ chains
#'
#' @param sce_list list of SingleCellExperiment objects
#' @param cache_path path for caching results
#' @param result_name unique result name used in caching
#'
#' @importFrom scanalysis cache
#'
#' @return
#' @export
#'
#' @examples
#' NULL
sc_compute_outlier_cell_filters = function(sce_list, cache_path, result_name) {
  cache_path = file.path(cache_path, paste0(result_name, ".qs"))

  result = cache(cache_path,
                 function() .sc_compute_outlier_cell_filters(sce_list))

  return(result)
}

#' Helper function for generating outlier cell filters
#'
#' Computes filters based on library sizes, number of genes expressed, percentage of mitochondrial reads, and number of VDJ chains
#'
#' @param sce_list list of SingleCellExperiment objects to compute filters for
#'
#' @importFrom scanalysis filter_total_umi filter_n_genes_expr filter_pct_mito filter_vdj_chain_count
#' @importFrom purrr map
#'
#' @return
#'
#' @examples
#' NULL
.sc_compute_outlier_cell_filters = function(sce_list) {

  total_umi_filters = map(sce_list, ~ filter_total_umi(.x, 3, type = "both"))

  n_genes_expr_filters = map(sce_list, ~ filter_n_genes_expr(.x, 3, type = "both"))

  pct_mito_filters = map(sce_list, ~ filter_pct_mito(.x, 3, type = "higher"))

  vdj_filters = map(sce_list, ~ filter_vdj_chain_count(.x, tra_range = c(0, 2), trb_range = c(0, 1)))

  return(list(total_umi_filters = total_umi_filters,
              n_genes_expr_filters = n_genes_expr_filters,
              pct_mito_filters = pct_mito_filters,
              vdj_filters = vdj_filters))
}
