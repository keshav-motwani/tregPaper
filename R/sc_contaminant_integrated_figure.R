#' Generate contaminant gene expression figure
#'
#' @param sce_list list of SCE objects
#' @param de_data data with output from sc_de
#' @param cluster_resolution name of cluster variable
#' @param contaminant_clusters ids of identified contaminant clusters
#'
#' @import cowplot
#' @import patchwork
#' @import ggplot2
#' @importFrom ggexp theme_ggexp
#' @importFrom readr write_csv
#' @importFrom dplyr filter group_by arrange top_n pull
#' @importFrom scanalysis plot_pairwise_features plot_features
#' @importFrom ggedit remove_geom
#' @importFrom purrr map
#'
#' @return
#' @export
#'
#' @examples
#' NULL
sc_contaminant_integrated_figure = function(sce_list,
                                            de_data,
                                            cluster_name,
                                            contaminant_clusters,
                                            main_result_path,
                                            supplemental_result_path,
                                            cache_path) {

  sce_list[[1]]$sample = map_chr(strsplit(sce_list[[1]]$sample, "_"), 1)

  ###### DIFFUSION MAP COLORED BY CLUSTERS
  clusters_diffusion_map = plot_reduced_dimensions(
    sce_list,
    type = "UMAP",
    features = cluster_name,
    facet_rows = c("sample"),
    scales = "fixed",
    facet_type = "wrap",
    nrow = 2,
    alpha = 1,
    label = cluster_name
  ) +
    scale_color_manual(values = get_clusters_palette()) +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(ncol = 3, override.aes = aes(alpha = 1, size = 3))) +
    labs(x = "UMAP Component 1",
         y = "UMAP Component 2") +
    theme(plot.margin = unit(c(2, 0.5, 2, 0.5), "lines"))

  ###### CLUSTER FREQUENCY COMPARISON BETWEEN SAMPLES
  cluster_composition_data = encode_cell_identity_frequency_long(
    sce_list,
    attributes = cluster_name,
    group_by = "sample",
    normalize = "none"
  )

  write_csv(
    cluster_composition_data,
    path = file.path(supplemental_result_path, "cluster_composition.csv")
  )

  cluster_composition_data$feature = factor(cluster_composition_data$feature, levels = rev(sort(unique(cluster_composition_data$feature))))

  cluster_composition = ggplot(cluster_composition_data,
                               aes(x = sample,
                                   y = feature,
                                   fill = value,
                                   label = value)) +
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


  genes = de_data %>%
    filter(group %in% contaminant_clusters) %>%
    group_by(group) %>%
    arrange(combined_p) %>%
    top_n(-5, combined_p) %>%
    pull(feature)

  ###### DIFFUSION MAP COLORED BY GENES
  genes_diffusion_map = plot_reduced_dimensions(
    sce_list,
    type = "UMAP",
    features = c(genes, "FOXP3", "IKZF2"),
    facet_columns = c("sample"),
    facet_rows = ".feature",
    switch = "y"
  ) +
    labs(x = "UMAP Component 1",
         y = "UMAP Component 2") +
    theme(legend.position = "none") +
    scale_color_gradient(
      low = "#E8E8E8",
      high = "firebrick",
      breaks = c(0, 1),
      labels = c("min", "max"),
      limits = c(0, 1)
    )

  ###### DIFFUSION COMPONENT 2 VS GENE EXPRESSION
  violin = plot_features(
    sce_list,
    assay = "logcounts",
    x = cluster_name,
    features = c(genes, "FOXP3", "IKZF2"),
    color = cluster_name,
    alpha = 1,
    facet_rows = ".feature",
    facet_columns = "sample",
    alt_exp = NULL,
    switch = "y",
    type = "quasirandom",
    facet_type = "grid",
    annotate_counts = FALSE
  ) +
    scale_color_manual(values = get_clusters_palette()) +
    theme(
      legend.position = "none",
      strip.placement = "outside",
      axis.text.x = element_text(angle = 90, hjust = 1)
    ) +
    labs(y = "Normalized Gene Expression",
         x = "Cluster")

  violin = remove_geom(violin, "point")

  groups = c("sample", ".feature", cluster_name)
  counts = violin$data %>%
    group_by(.dots = groups) %>%
    summarize(n = sum(value > 0))

  violin = violin + geom_jitter(alpha = 1, size = 0.6) + geom_boxplot(
    color = "black",
    alpha = 0,
    outlier.size = 0,
    width = 0.3
  ) +
    geom_text(
      data = counts,
      aes_string(label = "n",
                 x = cluster_name,
                 y = Inf),
      hjust = 0.5,
      vjust = 1,
      size = 2.5,
      color = "black",
      angle = 0
    )

  ###### TREG GENES VS CONTAMINANT GENES
  pairwise = plot_pairwise_features(
    sce_list,
    x = c("FOXP3", "IKZF2"),
    y = genes,
    color = cluster_name,
    alpha = 1,
    facet_rows = "sample",
    alt_exp = NULL,
    switch = "both",
    combination_group = NULL
  ) +
    scale_color_manual(values = get_clusters_palette()) +
    theme(legend.position = "none", strip.placement = "outside")

  counts = pairwise$data %>%
    group_by(.xkey, .ykey, sample) %>%
    summarize(n = sum(.xvalue != 0 & .yvalue != 0))

  pairwise = remove_geom(pairwise, "point") +
    geom_point(
      data = pairwise$data %>%
        filter(.xvalue == 0 | .yvalue == 0),
      alpha = 1,
      size = 0.6
    ) +
    geom_jitter(
      data = pairwise$data %>%
        filter(.xvalue != 0, .yvalue != 0),
      alpha = 1,
      size = 0.6
    ) +
    geom_text(
      data = counts,
      aes(
        label = n,
        x = Inf,
        y = Inf,
        color = NULL
      ),
      hjust = 1,
      vjust = 1
    )

  ###### ASSEMBLED PLOT
  assembled_plot = plot_grid(
    genes_diffusion_map,
    violin,
    plot_grid(pairwise) + geom_hline(yintercept = 0.522),
    rel_widths = c(3, 3, 3),
    nrow = 1,
    labels = c("C", "D", "E")
  )

  assembled_plot = plot_grid(
    plot_grid(
      clusters_diffusion_map + theme(legend.position = "none"),
      cluster_composition,
      ncol = 1,
      labels = c("A", "B")
    ),
    assembled_plot,
    ncol = 2,
    rel_widths = c(2, 9.5)
  )

  ggsave(
    file.path(main_result_path, "contaminant_figure.jpeg"),
    assembled_plot,
    height = 10,
    width = 14,
    units = "in",
    dpi = 300
  )
}
