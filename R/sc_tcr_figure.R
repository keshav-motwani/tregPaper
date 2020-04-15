#' Generate single-cell TCR figure
#'
#' @param sce_list list of SingleCellExperiment objects
#' @param cluster_resolution name of cluster
#' @param result_path result path
#'
#' @import ggplot2
#' @importFrom ggexp theme_ggexp
#' @importFrom readr write_csv read_csv
#' @importFrom ggedit remove_geom
#' @importFrom dplyr filter mutate group_by
#' @importFrom scanalysis annotate_cdr3
#'
#' @return
#' @export
#'
#' @examples
sc_tcr_figure = function(sce_list,
                                  cluster_resolution,
                                  main_result_path,
                                  supplemental_result_path,
                                  cache_path,
                                  result_name) {

  sce_list[[1]]$sample = map_chr(strsplit(sce_list[[1]]$sample, "_"), 1)

  ###### CLUSTER LEVEL EVENNESS PROFILE
  evenness_profile_cluster = compute_evenness_profile_long(sce_list, group_by = c("sample", cluster_resolution))

  write_csv(evenness_profile_cluster,
                   path = file.path(supplemental_result_path, paste0(result_name, "_evenness_profile_cluster.csv")))

  evenness_profile_cluster_plot = ggplot(
    evenness_profile_cluster,
    aes_string(x = "alpha", y = "evenness", color = cluster_resolution)
  ) +
    geom_jitter(size = 0.2) +
    geom_line(aes(group = group), alpha = 1) +
    theme_ggexp() +
    scale_color_manual(values = get_clusters_palette()) +
    scale_x_continuous(expand = c(0, 0.5)) +
    theme(legend.position = c(1, 0.01),
          legend.justification = c(1, 0.0)) +
    facet_grid(~ sample) +
    ylim(c(0, 1.01)) +
    labs(y = "Evenness", color = "Cluster") +
    guides(color = guide_legend(ncol = 2))

  ###### SAMPLE LEVEL EVENNESS PROFILE
  evenness_profile = compute_evenness_profile_long(sce_list, group_by = c("sample"))

  write_csv(evenness_profile,
                   path = file.path(supplemental_result_path, paste0(result_name, "_evenness_profile_sample.csv")))

  evenness_profile$test = ''
  evenness_profile_plot = ggplot(evenness_profile,
                                 aes_string(x = "alpha", y = "evenness", color = "sample")) +
    geom_jitter(size = 0.2) +
    geom_line(aes(group = group), alpha = 0.5) +
    theme_ggexp() +
    scale_color_manual(values = get_samples_palette()) +
    scale_x_continuous(expand = c(0, 0.5)) +
    theme(
      legend.position = c(1, 0.03),
      legend.justification = c(1, 0.1),
      legend.direction = "vertical"
    ) +
    ylim(c(0, 1.01)) +
    labs(color = "Sample", y = "Evenness") +
    facet_grid(~ test)

  ###### DIFFUSION MAP COLORED BY CLONOTYPE
  clonotype_diffusion = plot_reduced_dimensions(
    sce_list,
    type = "UMAP",
    features = c("clonotype_count_2"),
    facet_rows = c('sample'),
    facet_columns = c(cluster_resolution),
    scales = "fixed",
    alpha = 1
  ) +
    scale_color_manual(values = get_clonotypes_palette()) +
    theme(legend.position = "none") +
    labs(x = "UMAP Component 1",
         y = "UMAP Component 2")

  clonotype_diffusion = remove_geom(clonotype_diffusion, "point") +
    geom_point(data = clonotype_diffusion$data %>% filter(value == "< 2 counts"), alpha = 1, size = 0.3) +
    geom_point(data = clonotype_diffusion$data %>% filter(value != "< 2 counts"), alpha = 1, size = 1.2)

  ###### McPAS SPECIFICITY CATEGORY COMPOSITION
  mcpas = read_csv("data/single_cell/McPAS-TCR.csv") %>%
    filter(Species == "Human")

  for (i in 1:length(sce_list)) {
    sce_list[[i]] = annotate_cdr3(
      sce = sce_list[[i]],
      reference = mcpas,
      reference_cdr3_column = "CDR3.beta.aa",
      reference_annotation_column = "Category",
      chain_match = "TRB",
      max_dist = 0
    )
  }

  specificity_composition_data = encode_vdj_identity_frequency_long(
    sce_list,
    attributes = "Category",
    group_by = "sample",
    normalize = "none"
  )

  specificity_composition_data = specificity_composition_data %>%
    group_by(sample) %>%
    mutate(frequency = 100 * value / sum(value))

  write_csv(specificity_composition_data,
                   path = file.path(supplemental_result_path, paste0(result_name, "_specificity_composition.csv")))

  specificity_composition = ggplot(specificity_composition_data,
                               aes(x = sample,
                                   y = feature,
                                   fill = frequency,
                                   label = paste0(value, " (", round(frequency, 2), "%)"))) +
    geom_tile(color = "white", size = 1) +
    scale_fill_gradient(low = "#E8E8E8", high = "firebrick") +
    theme_ggexp() +
    geom_text(size = 3) +
    theme(
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      legend.position = "none"
    ) +
    labs(y = NULL, x = NULL)


  ###### ASSEMBLE FIGURE
  tcr_figure = plot_grid(
    clonotype_diffusion,
    plot_grid(
      plot_grid(
          plot_grid(
            evenness_profile_plot,
            evenness_profile_cluster_plot,
            rel_widths = c(4, 8)
          ),
        labels = "C",
        label_x = 0.35
      ),
      specificity_composition,
      nrow = 1,
      rel_widths = c(12, 5),
      labels = c("B", "D", "E")
    ),
    labels = c("A", ""),
    rel_heights = c(6, 7),
    ncol = 1
  )

  ggsave(
    file.path(main_result_path, paste0(result_name, "_tcr_figure.jpeg")),
    tcr_figure,
    height = 7,
    width = 9,
    units = "in",
    dpi = 300
  )

}
