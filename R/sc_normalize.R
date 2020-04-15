#' Normalize single-cell data
#'
#' Cells are normalized for biases related to sequencing depth using SCTransform. Residuals are scaled to have unit variance.
#'
#' @param sce_list list of SingleCellExperiment objects
#' @param cache_path path for caching results
#' @param result_name unique result name used in caching
#'
#' @importFrom purrr map
#' @importFrom scanalysis cache
#'
#' @return
#' @export
#'
#' @examples
sc_normalize = function(sce_list, cache_path, result_name) {
  result_path = file.path(cache_path, paste0(result_name, ".qs"))
  normalized = cache(result_path, function() map(sce_list, .sc_normalize))
  return(normalized)
}

#' Helper function to normalize data
#'
#' @param sce SingleCellExperiment object
#'
#' @importFrom scanalysis sce_to_seurat seurat_to_sce
#' @importFrom Seurat SCTransform
#'
#' @return
#'
#' @examples
.sc_normalize = function(sce) {
  seurat = sce_to_seurat(sce)
  seurat = SCTransform(seurat, vars.to.regress = c(), return.only.var.genes = FALSE, do.scale = TRUE)
  sce = seurat_to_sce(seurat)
  return(sce)
}
