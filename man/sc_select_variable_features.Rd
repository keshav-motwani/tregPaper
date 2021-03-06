% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sc_select_variable_features.R
\name{sc_select_variable_features}
\alias{sc_select_variable_features}
\title{Select highly variable features}
\usage{
sc_select_variable_features(sce_list, cache_path, result_name)
}
\arguments{
\item{cache_path}{path for caching results}

\item{result_name}{unique result name used in caching}

\item{raw_data_path}{path to raw data, must contain folders with cellranger count output}
}
\value{

}
\description{
Features are selected for downstream analysis. In pre-expansion data, the features are used from SCTransform, and these same features are also used for post-expansion data.
}
\examples{
NULL
}
