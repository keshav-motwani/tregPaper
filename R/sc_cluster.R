#' Cluster cells using Louvain community detection on SNN graph
#'
#' FindNeighbors and FindClusters are used for clustering. A range of resolution values are used, but we use 0.4 for the manuscript figures.
#'
#' @param sce_list list of SingleCellExperiment objects
#' @param cache_path path for caching results
#' @param result_name unique result name used in caching
#'
#' @importFrom purrr map
#'
#' @return
#' @export
#'
#' @examples
sc_cluster = function(sce_list, cache_path, result_name) {

  result_path = file.path(cache_path, paste0(result_name, ".rds"))

  clustered = cache(result_path, function() map(sce_list, .sc_cluster_single_sample))

  return(clustered)
}

#' Helper function to cluster single sample
#'
#' @param sce SingleCellExperiment object
#'
#' @importFrom Seurat FindNeighbors FindClusters
#' @importFrom scanalysis seurat_to_sce sce_to_seurat
#' @importFrom stringr str_pad
#'
#' @return
#'
#' @examples
.sc_cluster_single_sample = function(sce) {

  seurat = sce_to_seurat(sce)

  pca = ifelse("PCA" %in% names(seurat@reductions), "PCA", "pca")

  seurat = FindNeighbors(seurat, dims = 1:20, reduction = pca, k.param = 20)

  for (resolution in seq(0.4, 1.2, .1)) {
    seurat = FindClusters(seurat, resolution = resolution, graph.name = "RNA_snn")
    renamed_clusters = paste0("C", str_pad(as.numeric(Seurat::Idents(seurat)), 2, "0", side =
                                                      "left"))
    cluster_identifier = paste0("cluster_", resolution)
    seurat@meta.data[, cluster_identifier] = renamed_clusters
  }

  sce = seurat_to_sce(seurat)

  return(sce)
}
