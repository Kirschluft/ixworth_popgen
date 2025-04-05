# Assessment of FST values from VCFtools
## Determination of top FST windows, overlap related of genes and distribution visualization

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
  library(parallel)
  library(argparser)
  source("workflow/scripts/util.R")
})



parse_input <- function() {
  p <- arg_parser("Determine top FST windows and visualise distribution")
  p <- add_argument(p, "--files", nargs = Inf, help = "Input FST files from vcftools")
  p <- add_argument(p, "--genes", help = "GFF3 file containing genes")
  p <- add_argument(p, "--cores", default = 16, help = "Maximum number of cores to run")
  p <- add_argument(p, "--chromosomes", nargs = Inf, help = "List of chromosomes to use")
  p <- add_argument(p, "--out", help = "Output file prefix")
  args <- parse_args(p)
  if (is.nall(args$files)) {
    stop("--files should include a list of file names")
  }
  if (is.nall(args$genes)) {
    stop("--genes should be a path to a GFF3 file")
  }
  if (is.nall(args$chromosomes)) {
    stop("--chromosomes should be a list of chromosomes")
  }
  if (is.nall(args$out)) {
    stop("--out should include an output name")
  }
  dfs <- mclapply(args$files, \(x) read_input_vcftools(x, "windowed.weir.fst", column = "Comparison"), mc.cores = args$cores)
  genes <- import(args$genes, format = "gff3")
  genes <- genes[genes$type %in% c("gene", "ncRNA_gene", "lnc_RNA", "miRNA", "snoRNA", "snRNA", "scRNA", "rRNA", "tRNA", "Y_RNA", "pseudogene")]
  args$chromosomes <- order_chromosomes(args$chromosomes)
  return(list(args = args, dfs = dfs, genes = genes))
}

fst_analysis <- function(dfs, out, genes, chromosomes, cores) {
  res <- mclapply(dfs, \(df) {
    comp <- df %>%
      pull(Comparison) %>%
      dplyr::first() %>%
      str_replace_all(pattern = " ", replacement = "_")
    out <- str_glue("{out}{comp}_nsnps.png")
    res <- window_hist(df, "N_VARIANTS", "Number of SNPs", "Number of Windows", out)
  }, mc.cores = cores)

  res <- mclapply(dfs, \(df) {
    comp <- df %>%
      pull(Comparison) %>%
      dplyr::first() %>%
      str_replace_all(pattern = " ", replacement = "_")
    out <- str_glue("{out}{comp}_hist.png")
    res <- window_hist(df, "WEIGHTED_FST", expression(F[ST]), "Number of Windows", out, x_breaks = seq(0, 1, 0.1))
  }, mc.cores = cores)


  dfs <- mclapply(dfs, \(x) add_empirical_pvalue(x, "WEIGHTED_FST"), mc.cores = cores)
  df <- do.call(bind_rows, dfs)
  res <- windows_summary(df, out, "WEIGHTED_FST", "Comparison", expression(F[ST]))
  res <- windows_chr_summary(df, out, "WEIGHTED_FST", "Comparison", expression(F[ST]), "CHROM")

  sweeps <- mclapply(dfs, \(x) sweeps_gene_intersection(x, "WEIGHTED_FST", genes), mc.cores = cores)
  res <- lapply(seq_along(sweeps), \(index) {
    sweep <- sweeps %>% pluck(index)
    window <- dfs %>% pluck(index)
    comp <- window$Comparison[1] %>%
      str_replace_all(pattern = " ", replacement = "_")
    write_csv(sweep$summary %>% dplyr::rename(MEAN = WEIGHTED_FST), str_glue("{out}{comp}_sweeps.csv"))

    res <- manhattan_plot(window, str_glue("{out}{comp}_sweeps.png"), "WEIGHTED_FST", expression(F[ST]), sweep$summary)
    clusters <- find_clusters(window, val_column = "WEIGHTED_FST")
    write_csv(clusters, str_glue("{out}{comp}_top_sweep_clusters.csv"))

    summary <- sweep$summary %>%
      mutate(CHROM = factor(CHROM, levels = chromosomes))
    res <- freq_plot(summary, "CHROM", str_glue("{out}{comp}_sweeps_freq.png"))
  })

  res <- heatmap(df, "Comparison", "WEIGHTED_FST", str_glue("{out}heatmap.png"), expression(F[ST]))
  res <- heatmap_chr(df, "CHROM", "Comparison", "WEIGHTED_FST", str_glue("{out}heatmap_chr.png"), expression(F[ST]))
}

input <- parse_input()
fst_analysis(input$dfs, input$args$out, input$genes, input$args$chromosomes, input$args$cores)
