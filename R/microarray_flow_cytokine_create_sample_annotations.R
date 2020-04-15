#' Generate sample annotations from data based on column names
#'
#' @param microarray data frame with data
#'
#' @return
#' @export
#'
#' @examples
microarray_flow_cytokine_create_sample_annotations = function(data) {
  sample_annotations = data.frame(origin = factor(substr(colnames(data), start = 1, stop = 2)),
                                  cell_type = factor(substr(colnames(data), start = 3, stop = 4)))
  sample_annotations$cell_type = gsub("TR", "Treg", sample_annotations$cell_type)
  sample_annotations$cell_type = gsub("TC", "Tconv", sample_annotations$cell_type)
  sample_annotations$origin = gsub("PB", "APB", sample_annotations$origin)
  donor_id = stringr::str_sub(colnames(data), 5)
  sample_annotations$donor = paste0(sample_annotations$origin, donor_id)
  return(sample_annotations)
}
