#' Score cell-cycle phase in single-cell data
#'
#' Score cell-cycle using ScoreCellCycle in Seurat, with genes from Tirosh et al.
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
sc_score_cell_cycle = function(sce_list, cache_path, result_name) {

  cell_cycle_scored_datasets_path = file.path(cache_path, paste0(result_name, ".qs"))

  cell_cycle_scored = cache(cell_cycle_scored_datasets_path, function()
    map(sce_list, .sc_score_cell_cycle_single_sample))

  return(cell_cycle_scored)
}

#' Helper function to score cell-cycle on a single sample
#'
#' @param sce SingleCellExperiment object
#'
#' @importFrom Seurat CellCycleScoring
#' @importFrom scanalysis sce_to_seurat seurat_to_sce
#' @importFrom S4Vectors metadata<-
#'
#' @return
#'
#' @examples
.sc_score_cell_cycle_single_sample = function(sce) {

  seurat = sce_to_seurat(sce)

  s_genes = Seurat::cc.genes$s.genes
  g2m_genes = Seurat::cc.genes$g2m.genes
  seurat = CellCycleScoring(seurat, s.features = s_genes, g2m.features = g2m_genes)

  sce = seurat_to_sce(seurat)

  metadata(sce)$cell_cycle_genes = c(s_genes, g2m_genes)

  return(sce)
}
