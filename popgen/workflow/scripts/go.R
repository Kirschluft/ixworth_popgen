suppressPackageStartupMessages({
  library(tidyverse)
  library(gprofiler2)
  library(tidygenomics)
})


read_qtls <- function(file) {
  qtls <- read_tsv(file,
    comment = "#",
    col_names = c(
      "CHROM", "START", "END", "NAME", "SCORE", "STRAND",
      "THICKSTART", "THICKEND", "ITEMRGB", "BLOCKCOUNT",
      "BLOCKSIZES", "BLOCKSTARTS"
    ),
    show_col_types = FALSE
  )
  qtls <- qtls %>%
    mutate(CHROM = str_split(CHROM, pattern = "Chr\\.")) %>%
    mutate(CHROM = unlist(lapply(CHROM, \(x) x[[2]]))) %>%
    mutate(
      START = ifelse(START > END, END, START),
      END = ifelse(START > END, START, END)
    ) %>%
    dplyr::select(CHROM, START, END, NAME)
  return(qtls)
}

get_unique_genes <- function(df, column) {
  genes <- df %>%
    filter(!is.na(!!column)) %>%
    pull(!!column) %>%
    str_split(pattern = ",") %>%
    unlist() %>%
    unique()
  return(genes)
}

go_terms <- function(df, column, gprofiler_version, organism, out) {
  tryCatch(
    {
      genes <- get_unique_genes(df, column)
      res <- retry_with_delay({
        set_base_url(gprofiler_version)
        gost(genes, organism = organism, correction_method = "fdr", evcodes = TRUE)
      })

      p <- gostplot(res, interactive = FALSE)
      res_df <- res$result %>%
        mutate(parents = map(parents, \(x) str_c(x, collapse = ",")) %>%
          unlist()) %>%
        dplyr::select(-query, -significant) %>%
        relocate(term_name, term_id)
      write_csv(res_df, str_glue("{out}_go_terms.csv"))
      ggsave(str_glue("{out}_go_plot.png"), p, units = "px", height = 1080, width = 1920)
    },
    error = function(e) {
      print(str_glue("An error occured: {e}"))
      write_csv(data.frame(), str_glue("{out}_go_terms.csv"))
    }
  )
}

retry_with_delay <- function(expr, max_retries = 20, delay_seconds = 1) {
  for (i in 1:max_retries) {
    result <- tryCatch(expr,
      error = function(e) {
        cat("Attempt", i, "failed:", conditionMessage(e), "\n")
        Sys.sleep(delay_seconds)
        NULL
      }
    )

    if (!is.null(result)) {
      return(result)
    }
  }

  stop("Retry attempts exhausted.")
}

unique_collapse <- function(x) {
  res <- x %>%
    str_split(pattern = ",") %>%
    unlist() %>%
    unique() %>%
    as.character() %>%
    replace_na("") %>%
    str_c(collapse = ",")
  return(res)
}

collapse_sweeps <- function(df, max_distance, out) {
  regions <- df %>%
    genome_cluster(by = c("CHROM", "BIN_START", "BIN_END"), max_distance = max_distance) %>%
    group_by(cluster_id) %>%
    summarise(
      CHROM = dplyr::first(CHROM),
      START = min(BIN_START, na.rm = TRUE),
      END = max(BIN_END, na.rm = TRUE),
      GENEID = unique_collapse(GENEID),
      GENE = unique_collapse(GENE),
      N_VARIANTS = sum(N_VARIANTS),
      MEAN = mean(MEAN)
    ) %>%
    dplyr::select(-cluster_id)
  write_csv(regions, str_glue("{out}_collapsed_sweeps.csv"))
  return(regions)
}

intersect_qtls <- function(df, qtls, max_distance, out) {
  if (is.null(qtls)) {
    print("No qtls file provided ...\nWriting empty data frame")
    write_csv(data.frame(), str_glue("{out}_qtl.csv"))
  }

  qtl_regions <- genome_join_closest(df, qtls,
    by = c("CHROM", "START", "END"),
    max_distance = max_distance, distance_column_name = "DISTANCE",
    mode = "left"
  ) %>%
    transmute(
      CHROM = CHROM.x,
      SWEEP_START = START.x,
      SWEEP_END = END.x,
      GENEID,
      GENE,
      QTL_START = START.y,
      QTL_END = END.y,
      QTL_NAME = NAME,
    )

  qtl_regions_summary <- qtl_regions %>%
    group_by(CHROM, SWEEP_START, SWEEP_END) %>%
    summarise(
      GENEID = dplyr::first(GENEID),
      GENE = dplyr::first(GENE),
      QTL = unique_collapse(QTL_NAME)
    )
  write_csv(qtl_regions_summary, str_glue("{out}_qtl.csv"))
}

combine_split <- function(x, y) {
  split1 <- str_split(x, ",")
  split2 <- str_split(y, ",")
  map2(split1, split2, \(x, y) {
    c(x, y) %>%
      unique() %>%
      unlist() %>%
      str_c(collapse = ",")
  }) %>%
    unlist()
}

common_sweeps <- function(xp_ehh, fst, out) {
  common <- genome_join_closest(xp_ehh, fst,
    by = c("CHROM", "START", "END"),
    max_distance = 1, distance_column_name = "Distance"
  ) %>%
    transmute(
      CHROM = CHROM.x,
      START_XPEHH = START.x,
      END_XPEHH = END.x,
      START_FST = START.y,
      END_FST = END.y,
      GENEID = combine_split(GENEID.x, GENEID.y),
      GENE = combine_split(GENE.x, GENE.y)
    )
  write_csv(common, str_glue("{out}_common.csv"))
  return(common)
}
