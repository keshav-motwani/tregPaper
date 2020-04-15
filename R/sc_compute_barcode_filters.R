#' Compute barcode filters to identify barcodes containing cells vs ambient RNA
#'
#' Computes filters based on emptyDrops algorithm as well as inflection point in the library size vs library size rank curve.
#'
#' @param sce_list list of SingleCellExperiment objects
#' @param cache_path path for caching results
#' @param result_name unique result name used in caching
#'
#' @importFrom scanalysis cache
#' @importFrom purrr map
#'
#' @return
#' @export
#'
#' @examples
#' NULL
sc_compute_barcode_filters = function(sce_list, cache_path, result_name) {
  cache_path = file.path(cache_path, paste0(result_name, ".qs"))

  result = cache(cache_path,
                 function()
                   map(sce_list, .sc_compute_barcode_filters_single_sample))

  return(result)
}

#' Compute barcode filters for a single sample
#'
#' @param sce SingleCellExperiment object to compute filters for
#'
#' @importFrom scanalysis filter_ambient_barcode
#'
#' @return
#'
#' @examples
#' NULL
.sc_compute_barcode_filters_single_sample = function(sce) {
  set.seed(2387423)
  return(filter_ambient_barcode(sce, 5000, n_iters = 10000))
}