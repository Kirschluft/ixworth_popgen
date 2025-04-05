## Integrated haplotype score (iHS) post-analysis using rehh

suppressPackageStartupMessages({
  library(tidyverse)
  library(argparser)
  library(rtracklayer)
  source("workflow/scripts/sweeps.R")
})



parse_input <- function() {
  p <- arg_parser("Analyze iHS results based on rehh")
  p <- add_argument(p, "--input", help = "Output of iHS (csv)")
  p <- add_argument(p, "--genes", help = "GFF3 file containing genes")
  p <- add_argument(p, "--chromosomes", nargs = Inf, help = "List of chromosomes to use")
  p <- add_argument(p, "--window_size", default = 4e4, help = "Window size for candidate regions")
  p <- add_argument(p, "--window_step", default = 2e4, help = "Window step size (overlap) for candidate regions")
  p <- add_argument(p, "--out", help = "Output file prefix")
  args <- parse_args(p)
  if (is.nall(args$input)) {
    stop("--input should be a path to a csv file")
  }
  if (is.nall(args$genes)) {
    stop("--genes should be a path to a GFF3 file")
  }
  if (is.nall(args$chromosomes)) {
    stop("--chromosomes should be a list of chromosomes")
  }
  if (is.nall(args$window_size)) {
    stop("--window_size should be an integer")
  }
  if (is.nall(args$window_size)) {
    stop("--window_step should be an integer")
  }
  if (is.nall(args$out)) {
    stop("--out should include an output name")
  }

  genes <- import(args$genes, format = "gff3")
  genes <- genes[genes$type %in% c("gene", "ncRNA_gene", "lnc_RNA", "miRNA", "snoRNA", "snRNA", "scRNA", "rRNA", "tRNA", "Y_RNA", "pseudogene")]
  args$chromosomes <- order_chromosomes(args$chromosomes)
  return(list(args = args, genes = genes))
}

ihs_stats <- function(path, genes, chromosomes, window_size, window_step, out) {
  df <- read.csv(path, header = TRUE) %>%
    mutate(IHS = scale(IHS, center = TRUE, scale = TRUE))
  res <- window_hist(df, "IHS", "iHS", "Number of Windows", str_glue("{out}_hist.png"))

  detect_sweeps(df, genes, chromosomes, window_size, window_step, out,
    val_column = "IHS",
    val_name = "iHS", both_tails = TRUE, highlight = TRUE
  )
  detect_sweeps(df, genes, chromosomes, window_size, window_step, str_glue("{out}_reference/ref"),
    val_column = "IHS", val_name = "iHS", both_tails = FALSE, highlight = FALSE
  )
}


input <- parse_input()
res <- ihs_stats(input$args$input, input$genes, input$args$chromosomes, input$args$window_size, input$args$window_step, input$args$out)
