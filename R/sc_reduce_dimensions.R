#' Run PCA and UMAP
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
sc_reduce_dimensions = function(sce_list, cache_path, result_name) {

  reduced_dimension_datasets_path = file.path(cache_path, paste0(result_name, ".qs"))

  reduced_dimensions = cache(reduced_dimension_datasets_path, function()
    map(sce_list, .sc_reduce_dimensions_single_sample))

  return(reduced_dimensions)
}
#' Run PCA and UMAP on a single sample
#'
#' @param sce SingleCellExperiment object
#'
#' @importFrom Seurat ScaleData RunPCA RunUMAP VariableFeatures
#' @importFrom scater runDiffusionMap
#' @importFrom scanalysis sce_to_seurat seurat_to_sce
#'
#' @return
#'
#' @examples
.sc_reduce_dimensions_single_sample = function(sce) {

  seurat = sce_to_seurat(sce)

  if (all(dim(seurat@assays[[seurat@active.assay]]@scale.data) == c(0, 0))) {
    seurat = ScaleData(seurat, features = c(VariableFeatures(seurat), seurat@misc$cell_cycle_genes))
  }

  reduced = seurat %>%
    RunPCA(npcs = 20, features = VariableFeatures(seurat)) %>%
    RunPCA(npcs = 20, features = seurat@assays[[seurat@active.assay]]@misc$cell_cycle_genes, reduction.name = "cc_pca") %>%
    RunUMAP(dims = 1:20) %>%
    RunUMAP(dims = 1:20, reduction = "cc_pca", reduction.name = "cc_umap") %>%
    seurat_to_sce()

  return(reduced)
}
