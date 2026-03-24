library(zellkonverter)
library(SingleCellExperiment)
library(tidyverse)

# raw data available at smb://research.drive.wisc.edu/ksankaran/scdesigner/data/04_differential_expression_benchmark/SRP162188.rds
obj <- readRDS("SRP162188.rds")
sce <- if (methods::is(obj, "SingleCellExperiment")) obj else as(obj, "SingleCellExperiment")

col_data <- as.data.frame(colData(sce)) |>
  select(starts_with("sra")) |>
  mutate(
    treatment = str_extract(sra.sample_attributes, "(?<=treatment;;)[^|]+"),
    row_id = row_number(),
    value = 1L
  ) |>
  pivot_wider(
    id_cols = c(row_id, starts_with("sra")),
    names_from = treatment,
    values_from = value,
    names_prefix = "treatment_",
    values_fill = 0L
  ) |>
  select(-row_id, -any_of("treatment_NA"))

colData(sce) <- DataFrame(col_data)
writeH5AD(sce, "SRP162188.h5ad")
