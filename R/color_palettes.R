#' Get palette for clusters
#'
#' @return
#' @export
#'
#' @examples
get_clusters_palette = function() {
  clusters_palette = ggsci::pal_d3("category20c")(20)
  names(clusters_palette) = paste0("C", stringr::str_pad(0:19, 2, "0", side = "left"))
  clusters_palette
}

#' Get palette for clonotypes
#'
#' @return
#' @export
#'
#' @examples
get_clonotypes_palette = function() {
  clonotypes_palette = c(
    "#D3D3D3",
    ggsci::pal_aaas()(10),
    ggsci::pal_igv()(7)[3:7],
    ggsci::pal_d3(palette = "category20b")(20),
    ggsci::pal_d3(palette = "category20")(20),
    ggsci::pal_uchicago(),
    ggsci::pal_igv()(51)
  )
  c(clonotypes_palette, clonotypes_palette, clonotypes_palette, clonotypes_palette, clonotypes_palette)
}

#' Get palette for samples
#'
#' @return
#' @export
#'
#' @examples
get_samples_palette = function() {
  samples_palette = c(APB_Treg_Pre = "#FB6542",
                      CB_Treg_Pre = "#375E97",
                      APB_Treg_Post = "#FB6542",
                      CB_Treg_Post = "#375E97",
                      APB = "#FB6542",
                      CB = "#375E97")
  samples_palette
}

#' Get palette for cell types
#'
#' @return
#' @export
#'
#' @examples
get_cell_type_palette = function() {
  cell_type_palette = c(Treg = "#EDF5E1",
                        Tconv = "#5CDB95")
  cell_type_palette
}
