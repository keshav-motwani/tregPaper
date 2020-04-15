#' Process raw microarray data
#'
#' @param microarray data frame with microarray data
#'
#' @return
#' @export
#'
#' @examples
microarray_process_data = function(microarray) {

  cols_keep = grep("Gene Symbol|^CB|^PB", colnames(microarray), value = TRUE)
  cols_keep = grep("UCSF", cols_keep, invert = TRUE, value = TRUE)
  microarray = microarray[, cols_keep]
  colnames(microarray) = gsub("Gene Symbol", "gene", colnames(microarray))
  colnames(microarray) = gsub(" ", "", colnames(microarray))

  genes = microarray$gene
  microarray$gene = NULL
  microarray = as.matrix(microarray)
  rownames(microarray) = make.unique(genes)

  attributes(microarray)$sample_annotations = microarray_flow_cytokine_create_sample_annotations(microarray)

  return(microarray)
}
