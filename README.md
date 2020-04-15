# tregPaper

All analyses can be reproduced using the following commands:

``` r
library(tregPaper)

sc_analysis(
  raw_data_path = "data/single_cell/",
  result_name = "single_cell",
  result_path = "results",
  cluster_resolution = "cluster_0.4",
  contaminant_clusters = "C06",
  expansion_state = "pre_expansion"
)

sc_analysis(
  raw_data_path = "data/single_cell/",
  result_name = "single_cell",
  result_path = "results",
  cluster_resolution = "cluster_0.4",
  contaminant_clusters = "C08",
  expansion_state = "post_expansion"
)

microarray_data = readr::read_csv("data/microarray/microarray.csv")
microarray_analysis(
  microarray_data,
  features_to_annotate = c("TNFRSF9", "GZMA", "IL7R", "GZMB", "GNLY", "IL18RAP", "LGALS3", "LGMN", "HES1", "SOX4", "GBP1", "STAT1", "GPR55", "DST", "TCF4", "THEMIS", "CTLA4", "TIGIT"),
  annotation_p_value_rank_threshold = 100,
  result_name = "microarray",
  main_result_path = "results/main/",
  supplemental_result_path = "results/supplemental/",
  data_path = "results/data/"
)

cytokine_data = readr::read_csv("data/cytokine/Luminex.csv")
cytokine_analysis(
  cytokine_data,
  features_to_plot = c("IL-20", "IL-22", "IL-26", "IL-27(p28)", "IL-12(p40)", "IL-12(p70)", "IL28A/IFNg2", "IL29/IFNg1", "IL-19", "IL-35", "IL-10", "IL-2"),
  result_name = "cytokine",
  main_result_path = "results/main/",
  supplemental_result_path = "results/supplemental/",
  data_path = "results/data/"
)

library(tregPaper)

flow_data = readr::read_csv("data/flow/PercentPos.csv")
flow_analysis(
  flow_data,
  features_to_plot = c("CD226", "CD73", "CD95L", "CD279", "CD28", "HLA-DR", "TIGIT", "CD95", "CD197", "CD194", "CD183", "CD49b"),
  result_name = "flow",
  main_result_path = "results/main/",
  supplemental_result_path = "results/supplemental/",
  data_path = "results/data/"
)
```

If you are interested in the implementations, please look in the `R` folder. Function names combined with function documentation should make it clear where something is implemented. If you have any questions, please create an issue on this page and I will be more than happy to help.