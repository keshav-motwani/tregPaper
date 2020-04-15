#' Run merged analysis and generate figures
#'
#' @param sce_list list of SCE objects
#' @param cluster_resolution resolution to use for Louvain community detection
#' @param contaminant_cluster_ids ids of identified contaminant clusters
#' @param result_name
#' @param result_path
#'
#' @import scanalysis
#'
#' @return
#' @export
#'
#' @examples
#' NULL
sc_analysis = function(raw_data_path,
                       result_path,
                       result_name,
                       cluster_resolution,
                       contaminant_clusters,
                       expansion_state = "pre_expansion") {

  main_result_path = file.path(result_path, result_name, expansion_state, "main")
  supplemental_result_path = file.path(result_path, result_name, expansion_state, "supplemental")
  cache_path = file.path(result_path, result_name, expansion_state, "cache")
  processed_data_path = file.path(result_path, result_name, "preprocessing")

  for (path in c(main_result_path, supplemental_result_path, cache_path, processed_data_path)) {
    dir.create(path, recursive = TRUE)
  }

  raw_data = sc_load(raw_data_path, processed_data_path, "raw_data")

  ambient_barcode_filters = sc_compute_barcode_filters(raw_data, processed_data_path, "barcode_filters")

  barcode_figure = sc_barcode_filter_figure(raw_data, ambient_barcode_filters, main_result_path, supplemental_result_path, cache_path)

  filtered_barcodes = sc_filter_barcodes(raw_data, ambient_barcode_filters, processed_data_path, "filtered_barcodes")
  rm(raw_data)
  rm(ambient_barcode_filters)

  outlier_cell_filters = sc_compute_outlier_cell_filters(filtered_barcodes, processed_data_path, "outlier_cell_filters")

  sc_outlier_cell_filter_figure(filtered_barcodes, outlier_cell_filters, barcode_figure, main_result_path, supplemental_result_path, cache_path)

  filtered_cells = sc_filter_outlier_cells(filtered_barcodes, outlier_cell_filters, processed_data_path, "filtered_cells")
  rm(filtered_barcodes)
  rm(outlier_cell_filters)

  subsetted = sc_subset_expansion_state(filtered_cells, expansion_state, cache_path, "subsetted")
  rm(filtered_cells)

  normalized = sc_normalize(subsetted, cache_path, "normalized")
  rm(subsetted)

  cell_cycle_scored = sc_score_cell_cycle(normalized, cache_path, "cell_cycle_scored")
  rm(normalized)

  genes_removed = sc_remove_genes(cell_cycle_scored, cache_path, "genes_removed")
  rm(cell_cycle_scored)

  feature_selected = sc_select_variable_features(genes_removed,
                                                 cache_path,
                                                 "feature_selected")
  rm(genes_removed)

  red_dim = sc_reduce_dimensions(feature_selected, cache_path, "red_dim")

  integrated = sc_integrate(feature_selected, cache_path, "integrated")
  rm(feature_selected)

  clustered_single = sc_cluster(red_dim, cache_path, "clustered_single")
  rm(red_dim)

  clustered_integrated = sc_cluster(integrated, cache_path, "clustered_integrated")
  rm(integrated)

  de_cluster_single = sc_de_cluster(
    clustered_single,
    cluster_resolution,
    main_result_path,
    supplemental_result_path,
    cache_path,
    "single"
  )

  de_cluster_integrated = sc_de_cluster(
    clustered_integrated,
    cluster_resolution,
    main_result_path,
    supplemental_result_path,
    cache_path,
    "integrated"
  )

  sc_volcano_single_figure(de_cluster_single,
                           main_result_path,
                           supplemental_result_path,
                           cache_path)

  sc_clustree_single_figure(clustered_single,
                            main_result_path,
                            supplemental_result_path,
                            cache_path)

  sc_qc_figure(
    clustered_single,
    main_result_path,
    supplemental_result_path,
    cache_path,
    "single"
  )

  sc_cell_cycle_figure(
    clustered_single,
    main_result_path,
    supplemental_result_path,
    cache_path,
    "single"
  )

  sc_cell_cycle_figure(
    clustered_integrated,
    main_result_path,
    supplemental_result_path,
    cache_path,
    "integrated"
  )

  sc_volcano_integrated_figure(de_cluster_integrated,
                               main_result_path,
                               supplemental_result_path,
                               cache_path)

  sc_clustree_integrated_figure(clustered_integrated,
                                main_result_path,
                                supplemental_result_path,
                                cache_path)

  sc_qc_figure(
    clustered_integrated,
    main_result_path,
    supplemental_result_path,
    cache_path,
    "integrated"
  )

  sc_contaminant_integrated_figure(
    clustered_integrated,
    de_cluster_integrated,
    cluster_resolution,
    contaminant_clusters,
    main_result_path,
    supplemental_result_path,
    cache_path
  )

  sc_tcr_figure(
    clustered_integrated,
    cluster_resolution,
    main_result_path,
    supplemental_result_path,
    cache_path,
    "integrated"
  )
}
