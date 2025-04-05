## Convert VCF file into GDS

suppressPackageStartupMessages({
  library(SNPRelate)
  library(argparser)
  source("workflow/scripts/util.R")
})


parse_input <- function() {
  p <- arg_parser("Visualise PCA for (pruned) populations")
  p <- add_argument(p, "--vcf", help = "VCF file")
  p <- add_argument(p, "--out", help = "Output file name")
  args <- parse_args(p)
  if (is.nall(args$vcf)) {
    stop("--files should be a path to a VCF file")
  }
  if (is.nall(args$out)) {
    stop("--out should include an output name")
  }

  return(args)
}


input <- parse_input()
snpgdsVCF2GDS(input$vcf, input$out, method = "biallelic.only")