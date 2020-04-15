#' Integrate samples using anchor-based method in Seurat v3
#'
#' @param sce_list list of SingleCellExperiment objects
#' @param cache_path path for caching results
#' @param result_name unique result name used in caching
#'
#' @return
#' @export
#'
#' @examples
sc_integrate = function(sce_list, cache_path, result_name) {

  integrated_datasets_path = file.path(cache_path, paste0(result_name, ".qs"))

  integrated = cache(integrated_datasets_path, function() .sc_integrate(sce_list))

  return(list(integrated))
}

#' Helper function to integrate data
#'
#' @param sce_list list of SingleCellExperiment objects
#'
#' @importFrom purrr map
#' @importFrom scanalysis sce_to_seurat seurat_to_sce
#' @importFrom Seurat SelectIntegrationFeatures PrepSCTIntegration FindIntegrationAnchors IntegrateData DefaultAssay<- ScaleData RunPCA RunUMAP
#' @importFrom scater runDiffusionMap
#'
#' @return
#'
#' @examples
.sc_integrate = function(sce_list) {

  seurat_list = map(sce_list, sce_to_seurat)

  features = SelectIntegrationFeatures(object.list = seurat_list)
  seurat_list = PrepSCTIntegration(object.list = seurat_list, anchor.features = features)
  anchors = FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", anchor.features = features)

  combined = IntegrateData(anchorset = anchors, normalization.method = "SCT")

  DefaultAssay(combined) = "integrated"

  combined = combined %>%
    RunPCA(npcs = 20) %>%
    RunUMAP(dims = 1:20) %>%
    seurat_to_sce()

  return(combined)
}
