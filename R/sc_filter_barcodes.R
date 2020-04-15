#' Filter out barcodes containing only ambient reads
#'
#' Filters barcodes and keeps only those that are called cells by emptyDrops and above the inflection point in the library size vs library size rank curve.
#'
#' @param sce_list list of SingleCellExperiment objects
#' @param cache_path path for caching results
#' @param result_name unique result name used in caching
#'
#' @importFrom scanalysis cache
#' @importFrom purrr map2
#'
#' @return
#' @export
#'
#' @examples
#' NULL
sc_filter_barcodes = function(sce_list, barcode_filters, cache_path, result_name) {
  cache_path = file.path(cache_path, paste0(result_name, ".qs"))

  result = cache(cache_path,
                 function()
                   map2(sce_list, barcode_filters, .sc_filter_barcodes_single_sample))

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
.sc_filter_barcodes_single_sample = function(sce, barcode_filter) {
  return(sce[, barcode_filter$empty_drops & barcode_filter$inflection])
}