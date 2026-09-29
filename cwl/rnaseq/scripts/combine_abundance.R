#!/usr/bin/env Rscript

use_data_table <- FALSE

# Check if data.table package is available
if (nzchar(system.file(package = "data.table"))) {
  library(data.table)
  use_data_table <- TRUE
} else {
  message("data.table package is not installed, falling back to base R")
}

if (use_data_table) {
  # data.table approach
  abundance.files <- list.files(pattern = "*.tsv")
  dfs <- lapply(abundance.files, data.table::fread)
  
  accession.ids <- sub(".tsv", "", abundance.files)
  
  rename_cols <- function(dt, suffix) {
    data.table::setnames(dt, old = c("eff_length", "est_counts", "tpm"),
                              new = paste0(suffix, c("_eff_length", "_est_counts", "_tpm")))
  }
  
  for (i in seq_along(dfs)) {
    rename_cols(dfs[[i]], accession.ids[i])
  }
  
  combined_df <- Reduce(function(x, y) merge(x, y, all = TRUE, by = c("target_id", "length")), dfs)
  
  data.table::fwrite(combined_df, "combined_abundance.tsv", sep = "\t", quote = FALSE)
  
} else {
  # Base R fallback
  abundance.files = list.files(pattern=("*.tsv"))
  accession.ids = sub(".tsv", "", abundance.files)
  combined_df = read.csv(abundance.files[1], sep="\t")
  names(combined_df) = c("target_id", "length", 
                        paste0(accession.ids[1], "_eff_length"), 
                        paste0(accession.ids[1], "_est_counts"), 
                        paste0(accession.ids[1], "_tpm"))
  for(i in 2:length(abundance.files)) {
    additional_df = read.csv(abundance.files[i], sep="\t")
    names(additional_df) = c("target_id", "length", 
                            paste0(accession.ids[i], "_eff_length"), 
                            paste0(accession.ids[i], "_est_counts"), 
                            paste0(accession.ids[i], "_tpm"))
    combined_df = merge(combined_df, additional_df, all=T)
  }
  write.table(combined_df, "combined_abundance.tsv", row.names=F, sep="\t", quote=F)
}
