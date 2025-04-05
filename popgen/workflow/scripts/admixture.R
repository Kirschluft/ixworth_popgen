## Assess ADMIXTURE

suppressPackageStartupMessages({
  library(tidyverse)
  library(argparser)
  source("workflow/scripts/util.R")
})



parse_input <- function() {
  p <- arg_parser("Evaluate ADMIXTURE runs")
  p <- add_argument(p, "--prefix", help = "prefix to ADMIXTURE runs")
  p <- add_argument(p, "--fam", help = "plink fam file used in ADMIXTURE runs")
  p <- add_argument(p, "--cv", help = "cross-validation error from ADMIXTURE")
  p <- add_argument(p, "--out", help = "Output file prefix")
  args <- parse_args(p)
  if (is.nall(args$prefix)) {
    stop("--prefix should be a file prefix for ADMIXTURE runs")
  }
  if (is.nall(args$fam)) {
    stop("--fam should be a fam file")
  }
  if (is.nall(args$cv)) {
    stop("--cv should be a text file containing the CV error of ADMIXTURE runs")
  }
  if (is.nall(args$out)) {
    stop("--out should include an output name")
  }

  fam <- read_table(args$fam, show_col_types = FALSE, col_names = c("Population", "Sample", "Mother", "Father", "Sex", "Phenotype"))
  cv <- read_lines(args$cv) %>%
    tibble(line = .) %>%
    mutate(
      K = str_extract(line, "K=([0-9]+)") %>% str_remove_all("K=") %>% as.numeric(),
      Error = as.numeric(str_extract(line, "[0-9.]+$"))
    ) %>%
    dplyr::select(K, Error)
  runs <- read_input(args$prefix, fam)
  return(list(args = args, runs = runs, cv = cv))
}


read_input <- function(prefix, fam) {
  files <- list.files(dirname(prefix), pattern = basename(prefix), full.names = TRUE)
  res <- lapply(files, \(x) {
    if (str_ends(x, pattern = ".Q")) {
      run <- str_extract(x, pattern = "^*[0-9]+\\.Q$") %>%
        str_replace(pattern = "\\.Q", replacement = "")
      df <- read_table(x, show_col_types = FALSE, col_names = FALSE) %>%
        mutate(
          Population = pull(fam, Population) %>% str_replace_all(pattern = "_", replacement = " "),
          Sample = pull(fam, Sample),
          K = as.numeric(run)
        ) %>%
        pivot_longer(cols = starts_with("X"), names_to = "Cluster", values_to = "Proportion") %>%
        mutate(Cluster = as.numeric(str_replace(Cluster, pattern = "X", replacement = ""))) %>%
        mutate_at(c("Population", "Sample"), as.factor)
      return(df)
    }
  })

  df <- do.call(bind_rows, res) %>%
    mutate(
      Cluster = factor(Cluster, levels = sort(unique(Cluster))),
      K = factor(K, levels = sort(unique(K)))
    )
  return(df)
}

plot_cross_validation <- function(df, out) {
  p <- ggplot(df, aes(x = K, y = Error)) +
    geom_line(aes(group = 1)) +
    geom_point() +
    white_theme +
    ylab("Cross-Validation Error") +
    xlab("k") +
    scale_x_discrete(limits = as.character(1:10))
  ggsave(str_glue("{out}cross_validation.png"), p, width = 1920, height = 1080, units = "px", dpi = 300)
}

plot_admixture <- function(df, out) {
  p <- ggplot(df, aes(x = Sample, y = Proportion, fill = Cluster)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(K ~ Population, scales = "free_x", space = "free_x", switch = "x") +
    scale_y_continuous(
      name = "",
      sec.axis = dup_axis(name = "k")
    ) +
    white_theme +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      panel.spacing = unit(0.1, "lines"),
      strip.placement = "outside",
      strip.background.x = element_blank(),
      strip.text.x = element_text(size = 10, margin = margin(b = 2), angle = 45),
      strip.text.y = element_text(size = 10),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_blank()
    )
  ggsave(str_glue("{out}admixture.png"), p, width = 1920, height = 1080, units = "px", dpi = 300)
}



input <- parse_input()
plot_cross_validation(input$cv, input$args$out)
plot_admixture(input$runs, input$args$out)
