#' Subset cells for specific expansion state
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
sc_subset_expansion_state = function(sce_list,
                                     expansion_state,
                                     cache_path,
                                     result_name) {
  sce_list[grep(strsplit(expansion_state, "_")[[1]][[1]],
                names(sce_list),
                ignore.case = TRUE)]
}
