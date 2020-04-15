#' Process raw microarray data
#'
#' @param data data frame with microarray data
#'
#' @return
#' @export
#'
#' @examples
flow_cytokine_process_data = function(data) {

  colnames(data) = gsub("X1", "Sample", colnames(data))
  colnames(data) = gsub(" ", "", colnames(data))

  rows_keep = grep("^CB|^PB", data$Sample, value = TRUE)
  rows_keep = grep("GMP|UCSF", rows_keep, invert = TRUE, value = TRUE)

  data = data[data$Sample %in% rows_keep, ]
  data = data[complete.cases(data),]
  data[, "CD38"] = NULL

  sample_names = data$Sample
  data$Sample = NULL

  data = t(as.matrix(data))
  colnames(data) = sample_names

  sample_annotations = microarray_flow_cytokine_create_sample_annotations(data)

  data = cbind(t(data), sample_annotations) %>%
    tidyr::pivot_longer(cols = rownames(data), names_to = "feature") %>%
    tidyr::unite(group, c("origin", "cell_type"))

  return(data)
}
