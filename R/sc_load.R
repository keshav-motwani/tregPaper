#' Load single-cell data
#'
#' Loads gene expression and VDJ data, as well as annotate library size, number of genes expressed, and percentage of mitochondrial reads per cell.
#'
#' @param raw_data_path path to raw data, must contain folders with cellranger count output
#' @param cache_path path for caching results
#' @param result_name unique result name used in caching
#'
#' @importFrom purrr map
#' @importFrom scanalysis cache
#'
#' @return
#' @export
#'
#' @examples
#' NULL
sc_load = function(raw_data_path,
                   cache_path,
                   result_name) {

  cache_path = file.path(cache_path, paste0(result_name, ".qs"))

  samples = c("APB_Treg_Pre", "CB_Treg_Pre", "APB_Treg_Post", "CB_Treg_Post")

  result = cache(cache_path,
                    function() map(samples, ~ .sc_load_single_sample(.x, raw_data_path)))
  names(result) = samples

  return(result)
}

#' Load single-cell data for a single sample
#'
#' @param sample sample name
#' @param raw_data_path path to raw data, must contain folders with cellranger count output
#'
#' @importFrom scanalysis read_10x annotate_total_umi_count annotate_n_genes_expr annotate_pct_gene_set annotate_chain_count
#' @importFrom SummarizedExperiment colData colData<-
#'
#' @return
#'
#' @examples
.sc_load_single_sample = function(sample, raw_data_path) {

  sce = read_10x(file.path(raw_data_path, sample, "raw_gene_bc_matrices", "hg19"),
                 file.path(raw_data_path, sample)) %>%
    annotate_total_umi_count() %>%
    annotate_n_genes_expr() %>%
    annotate_pct_gene_set("^MT-", "pct_mito") %>%
    annotate_chain_count() %>%
    assign_clonotypes(tra_range = c(1, 1), trb_range = c(1, 1)) %>%
    annotate_clonotype_count()

  colData(sce)$sample = sample
  colData(sce)$Barcode = paste0(sample, colData(sce)$Barcode)

  return(sce)
}
