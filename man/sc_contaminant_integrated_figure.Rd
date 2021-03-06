% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sc_contaminant_integrated_figure.R
\name{sc_contaminant_integrated_figure}
\alias{sc_contaminant_integrated_figure}
\title{Generate contaminant gene expression figure}
\usage{
sc_contaminant_integrated_figure(
  sce_list,
  de_data,
  cluster_name,
  contaminant_clusters,
  main_result_path,
  supplemental_result_path,
  cache_path
)
}
\arguments{
\item{sce_list}{list of SCE objects}

\item{de_data}{data with output from sc_de}

\item{contaminant_clusters}{ids of identified contaminant clusters}

\item{cluster_resolution}{name of cluster variable}
}
\value{

}
\description{
Generate contaminant gene expression figure
}
\examples{
NULL
}
