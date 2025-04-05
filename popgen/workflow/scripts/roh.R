## Assess ROH

suppressPackageStartupMessages({
  library(tidyverse)
  library(argparser)
  source("workflow/scripts/util.R")
})



parse_input <- function() {
  p <- arg_parser("Calculate runs of homozygosity")
  p <- add_argument(p, "--input", nargs = Inf, help = "Input roh files from bcftools")
  p <- add_argument(p, "--length", default = 3e9, help = "Genome length (bp)")
  p <- add_argument(p, "--min", default = 3e5, help = "Minimum ROH length")
  p <- add_argument(p, "--out", help = "Output file prefix")
  args <- parse_args(p)
  if (is.nall(args$input)) {
    stop("--input should be a path")
  }
  if (is.nall(args$length)) {
    stop("--length should be an integer describing the genome length in bp")
  }
  if (is.nall(args$out)) {
    stop("--out should be a path prefix")
  }

  files <- lapply(args$input, \(x) {
    df <- read_tsv(x, col_names = c("RG", "Sample", "Chromosome", "Start", "End", "Length", "Markers", "Quality"))
    popname <- x %>%
      basename() %>%
      str_remove_all("\\.[a-zA-Z0-9]*")
    df %>%
      mutate(Population = popname)
  })
  files <- bind_rows(files)

  return(list(files = files, length = args$length, min = args$min, out = args$out))
}

estimate_froh <- function(roh, genome_length, out) {
  froh <- roh %>%
    group_by(Population, Sample) %>%
    summarise(ROH_SUM = sum(Length)) %>%
    mutate(FROH = ROH_SUM / genome_length)

  write_csv(froh, str_glue("{out}/froh_indv.csv"))

  froh_pop <- froh %>%
    group_by(Population) %>%
    summarise(
      FROH_SD = round(sd(FROH), 3),
      FROH = round(mean(FROH), 3)
    )

  write_csv(froh_pop, str_glue("{out}/froh_population.csv"))
}

calculate_stats <- function(roh, out) {
  roh_class <- roh %>%
    mutate(Class = cut(Length,
      breaks = c(0, 3e5, 1e6, 2e6, 4e6, 8e6, 1e7, 1.6e7, 1e9),
      labels = c("0-0.3Mb", "0.3-1Mb", "1-2Mb", "2-4Mb", "4-8Mb", "8-10Mb", "10-16Mb", ">16Mb")
    )) %>%
    group_by(Population, Sample, Class) %>%
    summarise(Number = n(), .groups = "drop") %>%
    complete(nesting(Population, Sample), Class, fill = list(Number = 0)) %>%
    group_by(Population, Class) %>%
    summarise(Number = mean(Number), .groups = "drop")

  write_csv(roh_class, str_glue("{out}/roh_count_mean.csv"))

  p <- ggplot(roh_class, aes(x = Population)) +
    geom_bar(aes(y = Number), stat = "identity") +
    white_theme +
    ylab("Mean number of ROH")
  ggsave(str_glue("{out}/roh_count_mean.png"), p, width = 1920, height = 1080, units = "px")

  p <- ggplot(roh_class, aes(x = Class)) +
    geom_bar(aes(y = Number, fill = Population), stat = "identity", position = "dodge") +
    white_theme +
    ylab("Mean number of ROH")
  ggsave(str_glue("{out}/roh_count_mean_class.png"), p, width = 1920, height = 1080, units = "px")

  roh_class_chr <- roh %>%
    mutate(Class = cut(Length,
      breaks = c(0, 3e5, 1e6, 2e6, 4e6, 8e6, 1e7, 1.6e7, 1e9),
      labels = c("0-0.3Mb", "0.3-1Mb", "1-2Mb", "2-4Mb", "4-8Mb", "8-10Mb", "10-16Mb", ">16Mb")
    )) %>%
    group_by(Chromosome, Population, Sample, Class) %>%
    summarise(Number = n(), .groups = "drop") %>%
    complete(nesting(Chromosome, Population), Sample, Class, fill = list(Number = 0)) %>%
    group_by(Chromosome, Population, Class) %>%
    summarise(Number = mean(Number), .groups = "drop")

  write_csv(roh_class_chr, str_glue("{out}/roh_count_mean_chr.csv"))

  p <- ggplot(roh_class_chr, aes(x = as.factor(Chromosome))) +
    geom_bar(aes(y = Number, fill = Population), stat = "identity", position = "dodge") +
    white_theme +
    ylab("Mean number of ROH") +
    xlab("Chromosome") +
    scale_fill_viridis_d()
  ggsave(str_glue("{out}/roh_count_mean_chr.png"), p, width = 1920, height = 1080, units = "px")

  indiv_number <- roh %>%
    group_by(Population, Sample) %>%
    summarise(Number = n(), Length = sum(Length)) %>%
    group_by(Population) %>%
    summarise(
      MeanNumber = mean(Number),
      MeanIndSum = mean(Length),
      SDIndSum = sd(Length),
      .groups = "drop"
    )

  roh_stats <- roh %>%
    group_by(Population) %>%
    summarise(
      MinimumLength = min(Length),
      MaximumLength = max(Length),
      TotalLength = sum(Length),
      MeanLength = mean(Length),
      .groups = "drop"
    ) %>%
    inner_join(indiv_number, by = "Population")

  write_csv(roh_stats, str_glue("{out}/roh_summary.csv"))
}


main <- function() {
  input <- parse_input()
  roh <- input$files %>%
    filter(Quality >= 10) %>%
    filter(Markers >= 10) %>%
    filter(Length >= input$min)

  estimate_froh(roh, input$length, input$out)
  calculate_stats(roh, input$out)
}


main()
