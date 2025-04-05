# Calculation of Ne using NeEstimator

suppressPackageStartupMessages({
  library(dartR.base)
  library(dartR.popgen)
  library(argparser)
  library(tidyverse)
  source("workflow/scripts/util.R")
})



parse_input <- function() {
  p <- arg_parser("Calculate effective population size using NeEstimator")
  p <- add_argument(p, "--file", help = "Input vcf file (population)")
  p <- add_argument(p, "--path", help = "Path to NeEstimator binary")
  p <- add_argument(p, "--nchr", help = "Number of chromosomes")
  p <- add_argument(p, "--out", help = "Output file name")
  args <- parse_args(p)
  if (is.nall(args$file)) {
    stop("--file should be a vcf file (single population)")
  }
  if (is.nall(args$path)) {
    stop("--path should be a path to the NeEstimator binary")
  }
  if (is.nall(args$nchr)) {
    stop("--nchr should be the number of chromosomes")
  }
  if (is.nall(args$out)) {
    stop("--out should include an output name")
  }

  return(args)
}

main <- function() {
  args <- parse_input()
  gen <- gl.read.vcf(vcffile = args$file, mode = "genotype")
  res <- gl.LDNe(gen,
    neest.path = dirname(args$path), outfile = basename(args$out), outpath = dirname(args$out), critical = c(0.0, 0.05),
    Waples.correction = "nChromosomes", Waples.correction.value = as.numeric(args$nchr), plot.out = FALSE,
    verbose = 0
  )
  write_csv(res[[1]], args$out)
}

main()
