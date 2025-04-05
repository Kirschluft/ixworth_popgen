# Assessment of nucleotide diversity from VCFtools

suppressPackageStartupMessages({
  library(tidyverse)
  library(argparser)
  library(parallel)
  source("workflow/scripts/util.R")
})



parse_input <- function() {
  p <- arg_parser("Assess and visualise nucleotide diversity (PI)")
  p <- add_argument(p, "--files", nargs = Inf, help = "Input PI files from vcftools")
  p <- add_argument(p, "--cores", default = 16, help = "Maximum number of cores to run")
  p <- add_argument(p, "--out", help = "Output file prefix")
  args <- parse_args(p)
  if (is.nall(args$files)) {
    stop("--files should include a list of file names")
  }
  if (is.nall(args$out)) {
    stop("--out should include an output name")
  }

  dfs <- mclapply(args$files, \(x) read_input_vcftools(x, "windowed.pi", column = "Population"), mc.cores = args$cores)
  return(list(args = args, dfs = dfs))
}



input <- parse_input()


res <- mclapply(input$dfs, \(df) {
  comp <- df %>%
    pull(Population) %>%
    dplyr::first() %>%
    str_replace_all(pattern = " ", replacement = "_")
  out <- str_glue("{input$args$out}{comp}_hist.png")
  res <- window_hist(df, "PI", expression(pi), "Number of Windows", out, x_breaks = seq(0, 1, 0.1))
}, mc.cores = input$args$cores)


dfs <- mclapply(input$dfs, \(x) add_empirical_pvalue(x, "PI"), mc.cores = input$args$cores)
df <- do.call(bind_rows, input$dfs)
res <- windows_summary(df, input$args$out, "PI", "Population", expression(pi))
res <- windows_chr_summary(df, input$args$out, "PI", "Population", expression(pi), "CHROM")

res <- mclapply(dfs, \(window) {
  pop <- window$Population[1] %>%
    str_replace_all(pattern = " ", replacement = "_")
  out_manhattan <- str_glue("{input$args$out}{pop}_manhattan.png")
  res <- manhattan_plot(window, out_manhattan, "PI", expression(pi), round_ylim = FALSE)
}, mc.cores = input$args$cores)

res <- heatmap_chr(df, "CHROM", "Population", "PI", str_glue("{input$args$out}heatmap_chr.png"), expression(pi))
