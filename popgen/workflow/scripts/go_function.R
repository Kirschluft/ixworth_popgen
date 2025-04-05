## Gene Enrichment and QTL intersections with XP-EHH and FST Sweeps

suppressPackageStartupMessages({
  library(tidyverse)
  library(argparser)
  source("workflow/scripts/go.R")
  source("workflow/scripts/util.R")
})



parse_input <- function() {
  p <- arg_parser("Functional assessment of selective sweeps (FST, XP-EHH, iHS)")
  p <- add_argument(p, "--file", help = "sweeps file (csv)")
  p <- add_argument(p, "--gprofiler_version", default = "https://biit.cs.ut.ee/gprofiler", help = "gprofiler archive")
  p <- add_argument(p, "--organism", default = "ggallus", help = "gprofiler organism")
  p <- add_argument(p, "--qtls", default = NULL, help = "Animal QTLdb bed file")
  p <- add_argument(p, "--out", help = "Output file prefix")
  args <- parse_args(p)
  if (is.nall(args$file)) {
    stop("--file should be a path to a selective sweeps file (csv)")
  }
  if (is.nall(args$gprofiler_version)) {
    stop("--gprofiler_version indicates whether to use gprofilers archive (e.g. https://biit.cs.ut.ee/gprofiler_archive3/e106_eg53_p16/)")
  }
  if (is.nall(args$organism)) {
    stop("--organism refers to an identifier for an organism used for GO-term retrieval (see https://biit.cs.ut.ee/gprofiler/page/organism-list)")
  }
  if (is.nall(args$out)) {
    stop("--out should include an output name")
  }

  input <- read_input(args$file, args$qtls)
  return(list(args = args, qtls = input$qtls, file = input$file))
}

read_input <- function(file, qtls) {
  f <- read_csv(file, show_col_types = FALSE) %>%
    mutate(CHROM = as.character(CHROM))
  if (!is.null(qtls)) qtls <- read_qtls(qtls)
  return(list(file = f, qtls = qtls))
}



input <- parse_input()
go_terms(input$file, "GENEID", input$args$gprofiler_version, input$args$organism, str_glue("{input$args$out}"))
sweeps <- collapse_sweeps(input$file, 1, str_glue("{input$args$out}"))
intersect_qtls(sweeps, input$qtls, 1, str_glue("{input$args$out}"))
