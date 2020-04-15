#' Run limma on microarray data based on a mixed effects model with fixed effect for treatment and random effect for donor
#'
#' @param data_matrix matrix with an attribute
#' @param subset_col Column of sample annotations to subset on
#' @param subset_value Allowed value for subsetting
#' @param comparison_col Column from sample annotations to compare across
#' @param paired Boolean indicating whether samples are paired or not
#'
#' @return
#' @export
#'
#' @examples
microarray_de = function(data_matrix,
                           sample_annotations,
                           treatment_columns,
                           donor_column) {

  donor = sample_annotations[, donor_column, drop = TRUE]

  group = factor(apply(sample_annotations[, treatment_columns, drop = FALSE], 1, function(x) paste0(x, collapse = "_")))

  design = model.matrix(~ 0 + group)
  colnames(design) = levels(group)

  corfit = limma::duplicateCorrelation(data_matrix, design, block = donor)

  fit = limma::lmFit(
    data_matrix,
    design = design,
    block = donor,
    correlation = corfit$consensus.correlation
  )

  contrasts = c(
    "CB_Treg vs CB_Tconv" = "CB_Treg - CB_Tconv",
    "APB_Treg vs APB_Tconv" = "APB_Treg - APB_Tconv",
    "CB_Treg vs APB_Treg" = "CB_Treg - APB_Treg",
    "CB_Tconv vs APB_Tconv" = "CB_Tconv - APB_Tconv"
  )

  cm = do.call(limma::makeContrasts, c(as.list(contrasts), list(levels = design)))

  fit2 = limma::contrasts.fit(fit, cm)

  fit2 = limma::eBayes(fit2)

  tt = purrr::map_dfr(
    names(contrasts),
    ~ limma::topTable(fit2,
                       coef = .x,
                       number = Inf,
                       sort.by = "none") %>%
    dplyr::mutate(p_val_adj = adj.P.Val,
                  p_val = P.Value,
                  gene = rownames(.),
                  contrast = .x)
  )

  return(tt)
}
