## Integrated haplotype score (iHS) analysis using rehh

suppressPackageStartupMessages({
  library(rehh)
  library(tidyverse)
  library(rtracklayer)
  library(parallel)
  library(argparser)
  source("workflow/scripts/util.R")
})



parse_input <- function() {
  p <- arg_parser("Run iHS analysis based on rehh")
  p <- add_argument(p, "--file", help = "Phased VCF file")
  p <- add_argument(p, "--chromosomes", nargs = Inf, help = "List of chromosomes to use")
  p <- add_argument(p, "--window_size", default = 4e4, help = "Window size for candidate regions")
  p <- add_argument(p, "--window_step", default = 2e4, help = "Window step size (overlap) for candidate regions")
  p <- add_argument(p, "--cores", default = 16, help = "Maximum number of cores to run")
  p <- add_argument(p, "--out", help = "Output file prefix")
  args <- parse_args(p)
  if (is.nall(args$file)) {
    stop("--file should be a path to a VCF or GDS file")
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

  return(args)
}

calculate_ihs_chr <- function(file, chr, threads = 16) {
  haplo <- data2haplohh(hap_file = file, chr.name = chr)
  haplo_scan <- scan_hh(haplo, threads = threads)
  ihs_results <- ihh2ihs(haplo_scan, include_freq = TRUE, min_maf = 0)
  ihs_results$frequency.class <- ihs_results$frequency.class %>%
    mutate(CHR = chr)
  return(ihs_results)
}

calculate_ihs <- function(file, chromosomes, jobs = 4, threads_per_job = 16) {
  ihs_res <- mclapply(chromosomes, \(chr) {
    calculate_ihs_chr(file, chr, threads_per_job)
  }, mc.cores = jobs)
  ihs <- lapply(ihs_res, \(x) x$ihs)
  ihs <- do.call(rbind, ihs)
  frequency_class <- lapply(ihs_res, \(x) x$frequency.class)
  frequency_class <- do.call(rbind, frequency_class)
  return(list(ihs = ihs, frequency_class = frequency_class))
}

run_ihs_analysis <- function(file, chromosomes, cores, out) {
  threads_per_job <- 8
  jobs <- floor(cores / threads_per_job)
  ihs_breed <- calculate_ihs(file, chromosomes, jobs, threads_per_job)
  write_csv(ihs_breed$ihs, str_glue("{out}_ihs.csv.gz"))
  write_csv(ihs_breed$frequency_class, str_glue("{out}_freq.csv.gz"))
  return(ihs_breed)
}

input <- parse_input()
ihs <- run_ihs_analysis(input$file, input$chromosome, input$cores, input$out)
