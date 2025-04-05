suppressPackageStartupMessages({
  library(rehh)
  library(tidyverse)
  source("workflow/scripts/util.R")
})

detect_sweeps <- function(df, genes, chromosomes, window_size, window_step, out,
                          val_column, val_name,
                          both_tails = TRUE, highlight = TRUE) {
  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)

  regions <- calc_candidate_regions(df,
    threshold = 2, window_size = window_size,
    overlap = window_step, join_neighbors = FALSE,
    min_n_mrk = 10, min_n_extr_mrk = 0,
    ignore_sign = both_tails
  )
  regions <- add_empirical_pvalue(regions, "MEAN_MRK") %>%
    arrange(CHR, START, END)
  write_csv(regions, str_glue("{out}_regions.csv"))

  regions <- regions %>%
    transmute(
      CHROM = factor(CHR, levels = chromosomes),
      BIN_START = START,
      BIN_END = END,
      N_VARIANTS = N_MRK,
      MEAN = MEAN_MRK,
      EMPIRICAL_PVALUE,
    )

  sweeps <- sweeps_gene_intersection(regions, "MEAN", genes)
  within_sweep_distribution(sweeps$summary, df, chromosomes, val_column, val_name, str_glue("{out}_sweeps_dist/sweeps"), both_tails)
  write_csv(sweeps$summary %>% arrange(desc(MEAN)), str_glue("{out}_sweeps.csv"))

  freq_plot(sweeps$summary, "CHROM", str_glue("{out}_sweeps_freq.png"))
  clusters <- find_clusters(regions)
  write_csv(clusters, str_glue("{out}_top_sweep_clusters.csv"))

  ylabel <- if (both_tails) substitute(group("|", NAME, "|")[w], list(NAME = as.name(val_name))) else substitute(NAME[w], list(NAME = as.name(val_name)))
  if (highlight) {
    res <- manhattan_plot(regions, str_glue("{out}_sweeps.png"), "MEAN", as.expression(ylabel), sweeps$summary)
  } else {
    res <- manhattan_plot(regions, str_glue("{out}_sweeps.png"), "MEAN", as.expression(ylabel), NULL)
  }
}
