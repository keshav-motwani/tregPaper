#' Find differentially expressed genes across clusters
#'
#' presto::wilcoxauc is used for differential expression on each individual cluster, and metap::minimump is used to combine p-values across samples (similar to FindConservedMarkers in Seurat).
#'
#' @param sce_list list of SingleCellExperiment objects
#' @param cluster_resolution cluster resolution to use
#' @param main_result_path path for storing main results
#' @param supplemental_result_path path for storing supplemental results
#' @param cache_path path for caching results
#' @param result_name unique result name used in caching
#'
#' @return
#' @export
#'
#' @examples
sc_de_cluster = function(sce_list,
                         cluster_resolution,
                         main_result_path,
                         supplemental_result_path,
                         cache_path,
                         result_name) {

  de_genes = imap_dfr(
    sce_list,
    ~ .identify_cluster_markers_presto(.x, cluster_resolution, "sample") %>%
      mutate(.sample = .y)
  )

  de_genes %>%
    dplyr::group_by(group) %>%
    dplyr::arrange(combined_p) %>%
    dplyr::top_n(-50, combined_p) %>%
    dplyr::arrange(group, combined_p) %>%
    readr::write_csv(file.path(
      supplemental_result_path,
      paste0(result_name, "_de_cluster_top_50.csv")
    ))

  de_genes %>%
    dplyr::arrange(group, combined_p) %>%
    readr::write_csv(file.path(
      supplemental_result_path,
      paste0(result_name, "_de_cluster.csv")
    ))

  return(de_genes)
}

#' Compute cluster defining genes from a SingleCellExperiment object using presto::wilcoxauc
#'
#' @param sce SingleCellExperiment object
#' @param cluster_resolution resolution to  use
#' @param grouping_vars Column name(s) from metadata to compute DE genes independently on and combine with metap package
#'
#' @return
#'
#' @examples
.identify_cluster_markers_presto = function(sce, cluster_resolution, grouping_vars) {
  seurat = scanalysis::sce_to_seurat(sce)

  seurat@meta.data = seurat@meta.data %>%
    tidyr::unite(group, grouping_vars, remove = FALSE)

  seurat_list = Seurat::SplitObject(seurat, split.by = "group")

  if (length(seurat_list) == 1) {
    names(seurat_list) = c(1)
  }

  if ("SCT" %in% names(seurat_list[[1]]@assays)) {
    seurat_assay = "SCT"
  } else {
    seurat_assay = "RNA"
  }

  de_per_group = purrr::imap(
    seurat_list,
    ~ presto::wilcoxauc(.x, cluster_resolution, seurat_assay = seurat_assay) %>%
      dplyr::rename_at(dplyr::vars(-c(feature, group)), paste0, .y)
  )

  de = purrr::reduce(de_per_group, dplyr::inner_join, by = c("group", "feature"))

  pval_codes = colnames(x = de)[grepl(pattern = "^pval", x = colnames(x = de))]

  if (length(seurat_list) > 1) {
    combined_p = apply(
      X = de[, pval_codes, drop = FALSE],
      MARGIN = 1,
      FUN = function(x) {
        return(metap::minimump(x)$p)
      }
    )
    de$combined_p = combined_p
  } else {
    de$combined_p = de$padj1
  }

  return(de)
}
