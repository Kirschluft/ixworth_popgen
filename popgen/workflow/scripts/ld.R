# Visualisation of LD-decay based on PopLdDecay

suppressPackageStartupMessages({
  library(tidyverse)
  library(argparser)
  source("workflow/scripts/util.R")
})


parse_input <- function() {
  p <- arg_parser("Visualize LD decay on a genome and chromosome scale")
  p <- add_argument(p, "--files", nargs = Inf, help = "lD files from PopLdDecay to visualize")
  p <- add_argument(p, "--out", help = "Output file prefix")
  p <- add_argument(p, "--column", default = "Population", help = "Legend name")
  p <- add_argument(p, "--binwidth", default = 1000, help = "Width of binned mean values (after 1000 bp)")
  args <- parse_args(p)
  if (is.nall(args$files)) {
    stop("--files should include a list of file names")
  }
  if (is.nall(args$out)) {
    stop("--out should include an output name")
  }
  if (is.nall(args$column)) {
    stop("--column should be either Population or Chromosome")
  }

  dfs <- lapply(args$files, \(x) {
    df <- read_input(x, args$column) %>%
      transmute(position = `#Dist`, meanr2 = `Mean_r^2`, !!sym(args$column)) %>%
      bin_ld(args$column, args$binwidth)
    return(df)
  })
  df <- do.call(bind_rows, dfs)
  if (args$column == "Chromosome") {
    df <- df %>%
      mutate(Chromosome = order_chromosomes(Chromosome))
  }
  return(list(args = args, df = df))
}

is.nall <- function(x) {
  return(length(x) == 1 && (is.null(x) || is.na(x)))
}

read_input <- function(x, column = "genome") {
  if (column == "Population") {
    population_name <- str_match(x, ".*/([^.]+).stat.gz")[, 2] %>%
      str_replace_all(pattern = "_", " ")
    df <- read_tsv(x, show_col_types = FALSE) %>%
      mutate(Population = population_name)
    return(df)
  } else if (column == "Chromosome") {
    chr_name <- str_match(x, ".*_chr([^.]+).stat.gz")[, 2]
    df <- read_tsv(x, show_col_types = FALSE) %>%
      mutate(Chromosome = chr_name)
  } else {
    stop("column not supported")
  }
}

bin_ld <- function(df, group, binwidth) {
  first_group <- df %>%
    filter(position <= 1000)

  binned_group <- df %>%
    filter(position > 1000) %>%
    mutate(bin = cut(position, breaks = seq(min(position), max(position), by = binwidth))) %>%
    group_by(bin, !!sym(group)) %>%
    summarise(
      position = mean(position),
      meanr2 = mean(meanr2),
    )

  combined <- bind_rows(first_group, binned_group)
  return(combined)
}

ld_decay_populations <- function(df, group, out, xlimits = NULL) {
  ngroups <- df %>%
    pull(!!sym(group)) %>%
    unique()

  p <- df %>%
    ggplot(aes(x = position, y = meanr2)) +
    ylab(expression("Mean LD" ~ (italic(r)^2))) +
    xlab("Distance (bp)") +
    white_theme +
    scale_color_viridis_d()

  line_aes <- NULL
  linewidth <- NULL
  if (length(ngroups) > 1) {
    line_aes <- aes(color = !!sym(group))
    if (length(ngroups) > 10) {
      linewidth <- 0.2
      p <- p + guides(color = guide_legend(ncol = 11)) +
        theme(
          legend.position = "bottom",
          legend.text = element_text(size = 5),
          legend.key.size = unit(0.5, "lines")
        )
    }
  }
  p <- p + geom_line(line_aes, linewidth = linewidth, alpha = 0.7)
  if (!is.null(xlimits)) {
    p <- p + xlim(xlimits)
  } else {
    p <- p +
      scale_x_continuous(
        breaks = seq(0, max(df$position), 1e5),
        labels = seq(0, max(df$position) / 1e3, 100),
        name = "Distance (kb)"
      )
  }
  ggsave(out, p, units = "px", height = 1080, width = 1920)
}


input <- parse_input()
ld_decay_populations(input$df, input$args$column, str_glue("{input$args$out}.png"))
ld_decay_populations(input$df, input$args$column, str_glue("{input$args$out}_close.png"), xlimits = c(0, 1000))
