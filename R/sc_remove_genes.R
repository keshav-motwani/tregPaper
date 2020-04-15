#' Remove uninteresting mitochondrial, ribosomal, and TCR genes from single-cell data
#'
#' Any genes matching ^MT-|^MTAT|^MRPL|^MRPS|^RPS|^RPL|^MALAT1$|^TRA|^TRB are removed.
#'
#' @param sce_list list of SingleCellExperiment objects
#' @param cache_path path for caching results
#' @param result_name unique result name used in caching
#'
#' @return
#' @export
#'
#' @examples
#' NULL
sc_remove_genes = function(sce_list, cache_path, result_name) {

  filtered_datasets_path = file.path(cache_path, paste0(result_name, "qs"))
  filtered = cache(filtered_datasets_path, function() map(sce_list, .sc_remove_genes_single_sample))

  return(filtered)
}

#' Helper function to remove genes from single sample
#'
#' @param sce SingleCellExperiment object
#'
#' @importFrom SingleCellExperiment counts
#'
#' @return
#'
#' @examples
#' NULL
.sc_remove_genes_single_sample = function(sce) {
  indices = grep("^MT-|^MTAT|^MRPL|^MRPS|^RPS|^RPL|^MALAT1$|^TRA|^TRB",
                 rownames(sce),
                 invert = TRUE)
  sce = sce[indices, ]
  return(sce)
}
