#' Select highly variable features
#'
#' Features are selected for downstream analysis. In pre-expansion data, the features are used from SCTransform, and these same features are also used for post-expansion data.
#'
#' @param raw_data_path path to raw data, must contain folders with cellranger count output
#' @param cache_path path for caching results
#' @param result_name unique result name used in caching
#'
#' @importFrom purrr map
#' @importFrom S4Vectors metadata
#'
#' @return
#' @export
#'
#' @examples
#' NULL
sc_select_variable_features = function(sce_list, cache_path, result_name) {

  result_path = file.path(cache_path, paste0(result_name, ".qs"))
  if (grepl("post_expansion", cache_path)) {
    fn = .sc_select_variable_features_post_expansion
  } else {
    fn = .sc_select_variable_features_pre_expansion
  }

  result = cache(result_path, function() map(sce_list, ~ fn(.x, cache_path)))

  variable_features = unique(unlist(map(result, ~ metadata(.x)$variable_features)))

  fileConn = file(file.path(cache_path, "selected_features.txt"))
  writeLines(variable_features, fileConn)
  close(fileConn)

  return(result)
}

#' Helper function to select features from pre-expansion dataset
#'
#' @param sce SingleCellExperiment object
#' @param cache_path ignored in this function
#'
#' @return
#'
#' @examples
.sc_select_variable_features_pre_expansion = function(sce, cache_path) {
  return(sce)
}

#' Helper function to select features from post-expansion dataset
#'
#' @param sce SingleCellExperiment object
#' @param cache_path ignored in this function
#'
#' @importFrom S4Vectors metadata<-
#'
#' @return
#' @export
#'
#' @examples
.sc_select_variable_features_post_expansion = function(sce, cache_path) {

  fileConn = file(gsub("post_expansion", "pre_expansion", file.path(cache_path, "selected_features.txt")))
  variable_features = readLines(fileConn)
  close(fileConn)

  metadata(sce)$variable_features = variable_features

  return(sce)
}
