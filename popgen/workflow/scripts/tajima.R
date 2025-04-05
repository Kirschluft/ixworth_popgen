# Assessment of nucleotide diversity from VCFtools

suppressPackageStartupMessages({
  library(tidyverse)
  library(argparser)
  library(parallel)
  source("workflow/scripts/util.R")
})



parse_input <- function() {
  p <- arg_parser("Assess and visualise tajima's D")
  p <- add_argument(p, "--files", nargs = Inf, help = "Input Tajima.D files from vcftools")
  p <- add_argument(p, "--cores", default = 16, help = "Maximum number of cores to run")
  p <- add_argument(p, "--out", help = "Output file prefix")
  args <- parse_args(p)
  if (is.nall(args$files)) {
    stop("--files should include a list of file names")
  }
  if (is.nall(args$out)) {
    stop("--out should include an output name")
  }

  dfs <- mclapply(args$files, \(x) read_input_vcftools(x, "Tajima.D", column = "POPULATION", variant_col = "N_SNPS"), mc.cores = args$cores)
  return(list(args = args, dfs = dfs))
}



input <- parse_input()


res <- mclapply(input$dfs, \(df) {
  comp <- df %>%
    pull(POPULATION) %>%
    dplyr::first() %>%
    str_replace_all(pattern = " ", replacement = "_")
  out <- str_glue("{input$args$out}{comp}_hist.png")
  res <- window_hist(df, "TajimaD", "Tajima's D", "Number of Windows", out, x_breaks = seq(0, 1, 0.1))
}, mc.cores = input$args$cores)


dfs <- mclapply(input$dfs, \(x) add_empirical_pvalue(x, "TajimaD"), mc.cores = input$args$cores)
df <- do.call(bind_rows, input$dfs)
res <- windows_summary(df, input$args$out, "TajimaD", "POPULATION", "Tajima's D")
res <- windows_chr_summary(df, input$args$out, "TajimaD", "POPULATION", "Tajima's D", "CHROM")

res <- lapply(dfs, \(window) {
  pop <- window$POPULATION[1] %>%
    str_replace_all(pattern = " ", replacement = "_")
  out_manhattan <- str_glue("{input$args$out}{pop}_manhattan.png")
  res <- manhattan_plot(window, out_manhattan, "TajimaD", "Tajima's D")
})

res <- heatmap_chr(df, "CHROM", "POPULATION", "TajimaD", str_glue("{input$args$out}heatmap_chr.png"), "Tajima's D")
