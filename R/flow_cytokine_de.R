#' Run pairwise comparisons with a mixed effect model based on combined treatment variables and donor as a random effect
#'
#' @param data
#' @param treatment_columns
#' @param donor_column
#'
#' @return
#' @export
#'
#' @examples
flow_cytokine_de = function(data,
                    treatment_columns,
                    donor_column) {
  data = tidyr::unite(data, group, treatment_columns, remove = FALSE) %>%
    dplyr::mutate(donor = data[, donor_column, drop = TRUE])

  data = tidyr::nest(data, -feature)

  result = purrr::map2_dfr(data$data,
                           data$feature,
                           ~ pairwise_comparisons(.x) %>%
                             dplyr::mutate(feature = .y))

  return(result)
}

pairwise_comparisons = function(data, donor_column) {
  model = lme4::lmer(value ~ group + (1 | donor), data = data, REML = FALSE)
  emmeans::emmeans(model, pairwise ~ group, adjust = "none")$contrasts %>%
    broom::tidy()
}
