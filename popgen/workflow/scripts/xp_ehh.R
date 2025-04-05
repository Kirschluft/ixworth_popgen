## Cross-Population Extended Haplotype Homozygosity (XP-EHH) analysis using rehh

suppressPackageStartupMessages({
  library(rehh)
  library(tidyverse)
  library(parallel)
  library(argparser)
  library(rtracklayer)
  source("workflow/scripts/util.R")
})



parse_input <- function() {
  p <- arg_parser("Run XP-EHH analysis based on rehh on two populations")
  p <- add_argument(p, "--file1", help = "Phased VCF file (first population), ending with 'phased.vcf.gz'")
  p <- add_argument(p, "--file2", help = "Phased VCF file (second population), ending with 'phased.vcf.gz'")
  p <- add_argument(p, "--chromosomes", nargs = Inf, help = "List of chromosomes to use")
  p <- add_argument(p, "--cores", default = 16, help = "Maximum number of cores to run")
  p <- add_argument(p, "--out", help = "Output file prefix")
  args <- parse_args(p)
  if (is.nall(args$file1)) {
    stop("--file1 should be a path to a VCF or GDS file")
  }
  if (is.nall(args$file2)) {
    stop("--file2 should be a path to a VCF or GDS file")
  }
  if (is.nall(args$chromosomes)) {
    stop("--chromosomes should be a list of chromosomes")
  }
  if (is.nall(args$out)) {
    stop("--out should include an output name")
  }

  return(args)
}

calculate_xp_ehh_chr <- function(file1, file2, pop1, pop2, chr, threads) {
  hh_pop1 <- data2haplohh(hap_file = file1, chr.name = chr, min_maf = 0)
  hh_pop2 <- data2haplohh(hap_file = file2, chr.name = chr, min_maf = 0)
  scan_pop1 <- scan_hh(hh_pop1, threads = threads)
  scan_pop2 <- scan_hh(hh_pop2, threads = threads)
  xp_ehh_results <- ies2xpehh(scan_pop1, scan_pop2, popname1 = pop1, popname2 = pop2)
  return(xp_ehh_results)
}

calculate_xp_ehh <- function(file1, file2, pop1, pop2, chromosomes, jobs = 4, threads_per_job = 16) {
  xp_ehh_res <- mclapply(chromosomes, \(chr) {
    calculate_xp_ehh_chr(file1, file2, pop1, pop2, chr, threads_per_job)
  }, mc.cores = jobs)
  xp_ehh_res <- do.call(rbind, xp_ehh_res)
  return(xp_ehh_res)
}


run_xp_ehh_analysis <- function(file1, file2, chromosomes, cores, out) {
  threads_per_job <- 8
  jobs <- floor(cores / threads_per_job)
  pop1 <- get_name(file1, "phased.vcf.gz")
  pop2 <- get_name(file2, "phased.vcf.gz")
  xp_ehh <- calculate_xp_ehh(file1, file2, pop1, pop2, chromosomes, jobs, threads_per_job) %>%
    rename_with(.cols = all_of(names(.)[3]), .fn = ~"XPEHH")
  write_csv(xp_ehh, str_glue("{out}_xp_ehh.csv.gz"))
  return(xp_ehh)
}

input <- parse_input()
xpehh <- run_xp_ehh_analysis(input$file1, input$file2, input$chromosome, input$cores, input$out)