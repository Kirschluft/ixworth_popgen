# Visualisation of intersections between the populations

suppressPackageStartupMessages({
  library(tidyverse)
  library(UpSetR)
  library(argparser)
})


parse_input <- function() {
  p <- arg_parser("Visualize intersections between population")
  p <- add_argument(p, "--file", help = "Output file from vcf-compare")
  p <- add_argument(p, "--out", help = "Output file name")
  args <- parse_args(p)
  if (is.nall(args$file)) {
    stop("--file should include the output file path by vcf-compare")
  }
  if (is.nall(args$out)) {
    stop("--out should include an output name")
  }

  df <- read_input(args$file)
  return(list(args = args, df = df))
}

is.nall <- function(x) {
  return(length(x) == 1 && (is.null(x) || is.na(x)))
}

get_pop_names <- function(lines) {
  # find line that contains all comparisons
  pops_all <- lines[nchar(lines) == max(nchar(lines))]
  # split and get rid of VN and count
  parts <- str_split(pops_all, "\t")[[1]]
  parts <- parts[-c(1, 2)]
  # extract population names from files
  pop_names <- str_extract(parts, "(.)*(?=\\.vcf)")
  pop_names <- basename(pop_names)
  pop_names_sep <- str_c(pop_names, collapse = "|")
  return(list(sep = pop_names_sep, names = pop_names))
}

read_input <- function(x) {
  lines <- read_lines(x, skip = 7)
  vn <- str_starts(lines, "VN")
  lines <- lines[vn]
  pop_names <- get_pop_names(lines)
  df_list <- list()

  df_list <- lapply(lines, \(line) {
    parts <- str_split(line, "\\s+")[[1]]
    count <- as.numeric(parts[2])
    membership <- as.integer(pop_names$names %in% str_extract(parts, pop_names$sep))
    df_row <- data.frame(t(membership), count = count)
    return(list(df_row))
  })

  df <- bind_rows(df_list)
  colnames(df) <- c(pop_names$names, "count")
  return(df)
}

upset_plot <- function(df, out) {
  pop_columns <- 1:(ncol(df) - 1)
  upset_input <- df %>%
    select(all_of(pop_columns)) %>%
    lapply(\(x) rep(x, round(df$count / 1e3))) %>%
    as.data.frame() %>%
    rename_with(~ str_replace_all(., "_", " "))

  png(out, width = 1920, height = 1080, unit = "px")
  print(
    upset(upset_input,
      order.by = "freq",
      sets = colnames(upset_input)[pop_columns],
      mainbar.y.label = "Intersection Size (10³ SNPs)",
      sets.x.label = "Set Size (10³ SNPs)",
      text.scale = c(4.5, 4, 4.5, 4, 4.5, 3),
      nintersects = 25,
    )
  )
  dev.flush()
  invisible(dev.off())
}


input <- parse_input()
upset_plot(input$df, input$args$out)
