#' Generate microarray figure
#'
#' @param microarray_data
#' @param de
#' @param features_to_annotate
#' @param annotation_p_value_rank_threshold
#' @param main_result_path
#' @param supplemental_result_path
#' @param data_path
#'
#' @return
#' @export
#'
#' @examples
microarray_figure = function(microarray_data,
                             de,
                             features_to_annotate,
                             annotation_p_value_rank_threshold,
                             main_result_path = NULL,
                             supplemental_result_path = NULL,
                             data_path = NULL) {

  de$contrast = factor(de$contrast, levels = c("APB_Treg vs APB_Tconv", "CB_Treg vs CB_Tconv", "CB_Tconv vs APB_Tconv", "CB_Treg vs APB_Treg"))

  volcano = scanalysis::plot_volcano(de,
                                     "logFC",
                                     "p_val_adj",
                                     annotations = c("FOXP3", "IKZF2"),
                                     annotations_if_threshold = features_to_annotate,
                                     p_value_rank_threshold = annotation_p_value_rank_threshold,
                                     facet_rows = "contrast",
                                     facet_type = "wrap",
                                     nrow = 2,
                                     n_annotate_top = 0)

  canonical_features = intersect(rownames(microarray_data),
                                 c("UTS", "HLA-DRA", "ENTPD1", "RTKN2", "FCRL3", "FOXP3", "IKZF2", "TIGIT", "CTLA4",
                                   "IL2RA", "IL2RB", "CD40LG", "NELL2", "ANK3", "DENDD5A", "THEMIS", "TMEM71"))

  col_anno = ComplexHeatmap::columnAnnotation(
    df = attributes(microarray_data)$sample_annotations[, c("origin", "cell_type")] %>% dplyr::rename("Origin" = "origin", "Cell Type" = "cell_type"),
    col = list("Origin" = c(APB = "#FB6542", CB = "#375E97"), "Cell Type" = get_cell_type_palette()),
    gp = grid::gpar(col = "white"),
    show_annotation_name = TRUE
  )

  heatmap = ComplexHeatmap::Heatmap(t(scale(t(microarray_data)))[canonical_features, ],
                                    name = "Scaled\nExpression",
                                    col = c("royalblue", "cornflowerblue", "white", "gold", "orange"),
                                    cluster_rows = TRUE,
                                    cluster_columns = TRUE,
                                    show_row_dend = FALSE,
                                    show_column_dend = FALSE,
                                    show_row_names = TRUE,
                                    show_column_names = FALSE,
                                    row_names_side = "left",
                                    rect_gp = grid::gpar(col = "white"),
                                    top_annotation = col_anno,
                                    column_names_gp = grid::gpar(fontsize = 7),
                                    row_names_gp = grid::gpar(fontsize = 7))

  heatmap = grid::grid.grabExpr(ComplexHeatmap::draw(heatmap))

  microarray_figure = cowplot::plot_grid(volcano, heatmap, nrow = 1, rel_widths = c(4, 6), labels = "AUTO")

  cowplot::save_plot(file.path(main_result_path, "figure_4.jpeg"), plot = microarray_figure, base_height = 4, base_width = 10)
}
