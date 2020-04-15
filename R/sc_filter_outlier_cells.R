#' Filter out outlier_cells containing only ambient reads
#'
#' Filters outlier_cells and keeps only those that are called cells by emptyDrops and above the inflection point in the library size vs library size rank curve.
#'
#' @param sce_list list of SingleCellExperiment objects
#' @param cache_path path for caching results
#' @param result_name unique result name used in caching
#'
#' @importFrom scanalysis cache
#' @importFrom purrr pmap
#'
#' @return
#' @export
#'
#' @examples
#' NULL
sc_filter_outlier_cells = function(sce_list,
                                   outlier_cell_filters,
                                   cache_path,
                                   result_name) {
  cache_path = file.path(cache_path, paste0(result_name, ".qs"))

  result = cache(cache_path,
                 function()
                   pmap(
                     list(
                       sce_list,
                       outlier_cell_filters$total_umi_filters,
                       outlier_cell_filters$n_genes_expr_filters,
                       outlier_cell_filters$pct_mito_filters,
                       outlier_cell_filters$vdj_filters
                     ),
                     ~ ..1[, ..2 & ..3 & ..4 & ..5]
                   ))

  return(result)
}
