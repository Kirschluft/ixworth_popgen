suppressPackageStartupMessages({
  library(tidyverse)
})


white_theme <- theme_bw() +
  theme(
    axis.title = element_text(size = 10, face = "bold"),
    axis.text = element_text(size = 7),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key.size = unit(0.5, "lines"),
  )


is.nall <- function(x) {
  return(length(x) == 1 && (is.null(x) || is.na(x)))
}

get_name <- function(x, ext) {
  str_match(x, str_glue(".*/([^.]+).{ext}"))[, 2] %>%
    str_replace_all(pattern = "_", " ")
}

read_input_vcftools <- function(x, ext, nthreshold = 10, column = "Comparison", variant_col = "N_VARIANTS") {
  name <- get_name(x, ext)
  df <- read_tsv(x, show_col_types = FALSE) %>%
    filter(!!sym(variant_col) > nthreshold) %>%
    mutate(!!sym(column) := name)
  return(df)
}

window_hist <- function(df, x, xlab, ylab, out, x_breaks = NULL) {
  p <- ggplot(df, aes(x = !!sym(x))) +
    geom_histogram(bins = 100) +
    xlab(xlab) +
    ylab(ylab) +
    white_theme
  if (!is.null(x_breaks)) {
    p <- p + scale_x_continuous(breaks = x_breaks)
  }
  ggsave(out, p, width = 1920, height = 1080, units = "px")
}

windows_summary <- function(df, out, y, group, ylab) {
  windows_stats <- df %>%
    group_by(!!sym(group)) %>%
    summarise(
      MEAN = mean(!!sym(y), na.rm = TRUE),
      SD = sd(!!sym(y), na.rm = TRUE),
      MIN = min(!!sym(y), na.rm = TRUE),
      MAX = max(!!sym(y), na.rm = TRUE)
    )
  write_csv(windows_stats, str_glue("{out}summary.csv"))

  p <- ggplot(df) +
    geom_boxplot(aes(y = !!sym(y), x = !!sym(group)), fill = "gray35") +
    ylab(ylab) +
    scale_fill_viridis_d() +
    white_theme
  ggsave(str_glue("{out}boxplot.png"), p, width = 1920, height = 1080, units = "px")

  p2 <- ggplot(df) +
    geom_violin(aes(y = !!sym(y), x = !!sym(group)), fill = "gray35", adjust = 0.5, draw_quantiles = 0.5) +
    ylab(ylab) +
    white_theme
  ggsave(str_glue("{out}violinplot.png"), p2, width = 1920, height = 1080, units = "px")
}

order_chromosomes <- function(x) {
  is_num_chrom <- str_remove(x, "chr") %>%
    str_detect(pattern = "^[0-9]+$")
  num_chroms <- unique(x[is_num_chrom])
  char_chroms <- unique(x[!is_num_chrom])
  sorted_chroms <- c(sort(as.numeric(num_chroms)), sort(char_chroms))
  chrom_factor <- factor(x, levels = as.character(sorted_chroms))
  return(chrom_factor)
}

windows_chr_summary <- function(df, out, y, group, ylab, chr) {
  windows_chr_stats <- df %>%
    group_by(!!sym(chr), !!sym(group)) %>%
    summarise(
      SD = sd(!!sym(y), na.rm = TRUE),
      MIN = min(!!sym(y), na.rm = TRUE),
      MAX = max(!!sym(y), na.rm = TRUE),
      MEAN = mean(!!sym(y), na.rm = TRUE)
    ) %>%
    mutate(!!chr := order_chromosomes(!!sym(chr)))
  write_csv(windows_chr_stats, str_glue("{out}chr_summary.csv"))

  num_comps <- df %>%
    dplyr::pull(!!sym(group)) %>%
    unique() %>%
    length()

  chr_bar_plot <- ggplot(windows_chr_stats, aes(x = !!sym(chr), y = MEAN)) +
    geom_bar(stat = "identity", position = "dodge", fill = "gray35") +
    geom_errorbar(aes(ymin = MEAN - SD, ymax = MEAN + SD), width = 0.3, position = position_dodge(0.8)) +
    facet_wrap({{ group }}, nrow = round(num_comps / 2)) +
    xlab("Chromosome") +
    ylab(ylab) +
    white_theme +
    theme(
      axis.text = element_text(angle = 0, size = 4),
      legend.text = element_text(size = 6),
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 6)
    )
  ggsave(str_glue("{out}bar_chr.png"), chr_bar_plot, width = 1920, height = 1080, units = "px")

  chr_dense_plot <- ggplot(df, aes(x = !!sym(y))) +
    geom_density(alpha = 0.5, fill = "gray35") +
    facet_wrap({{ group }}, nrow = round(num_comps / 2)) +
    xlab(ylab) +
    ylab("Density") +
    white_theme +
    theme(
      axis.text = element_text(angle = 0, size = 5),
      legend.text = element_text(size = 6),
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 6)
    )
  ggsave(str_glue("{out}dense_chr.png"), chr_dense_plot, width = 1920, height = 1080, units = "px")
}


add_empirical_pvalue <- function(df, column) {
  df <- df %>%
    mutate(EMPIRICAL_PVALUE = 1 - rank(!!sym(column)) / (nrow(.) + 1))
  return(df)
}

manhattan_plot <- function(windows, out, column, ylab, sweeps = NULL,
                           chr = "CHROM", start = "BIN_START", end = "BIN_END",
                           window_size = 4e4, col = c("#d8b365", "#5ab4ac"), round_ylim = TRUE) {
  require(qqman)

  if (!end %in% colnames(windows)) {
    windows <- windows %>%
      mutate(!!sym(end) := !!sym(start) + window_size)
  }

  df <- windows %>%
    mutate(SNP = paste("SNP", 1:nrow(.))) %>%
    transmute(
      SNP,
      CHR = as.numeric(as.factor(!!sym(chr))),
      BP = (!!sym(end) + !!sym(start) - 1) / 2,
      VAL = !!sym(column)
    )

  if (round_ylim) {
    y_min <- floor(min(df$VAL, na.rm = TRUE))
    y_max <- ceiling(max(df$VAL, na.rm = TRUE))
  } else {
    y_min <- min(df$VAL, na.rm = TRUE)
    y_max <- max(df$VAL, na.rm = TRUE)
  }
  ylim <- c(y_min, y_max)


  if (!is.null(sweeps)) {
    threshold <- min(sweeps %>% dplyr::select(!!sym(column)), na.rm = TRUE)

    png(out, width = 1920, height = 1080, unit = "px")
    print(
      manhattan(df,
        logp = FALSE,
        suggestiveline = FALSE,
        genomewideline = threshold,
        p = "VAL",
        cex.axis = 1.25,
        cex.lab = 1.6,
        cex.main = 1.75,
        cex.sub = 1.4,
        mgp = c(2.5, 1, 0),
        col = col,
        ylim = ylim,
        ylab = ylab
      )
    )
  } else {
    png(out, width = 1920, height = 1080, unit = "px")
    print(
      manhattan(df,
        logp = FALSE,
        suggestiveline = FALSE,
        genomewideline = FALSE,
        p = "VAL",
        cex.axis = 1.25,
        cex.lab = 1.6,
        cex.main = 1.75,
        cex.sub = 1.4,
        mgp = c(2.5, 1, 0),
        col = col,
        ylim = ylim,
        ylab = ylab
      )
    )
  }
  dev.flush()
  invisible(dev.off())
}

sweeps_gene_intersection <- function(df, val_column, genes, max_distance = 0, threshold = 0.001) {
  df <- df %>%
    filter(EMPIRICAL_PVALUE < threshold)

  start(genes) <- start(genes) - max_distance
  end(genes) <- end(genes) + max_distance

  df_ranges <- GRanges(
    seqnames = df$CHROM,
    ranges = IRanges(start = df$BIN_START, end = df$BIN_END),
    N_VARIANTS = df$N_VARIANTS,
    VAL = df %>% pull(!!sym(val_column)),
    EMPIRICAL_PVALUE = df$EMPIRICAL_PVALUE,
  )

  intersection <- findOverlaps(genes, df_ranges)
  query <- genes[queryHits(intersection)] %>%
    as_tibble() %>%
    transmute(
      GENE_ID = str_remove(ID, "gene:|transcript:"),
      GENE = Name,
      GENE_CHROM = seqnames,
      GENE_START = start,
      GENE_END = end,
      GENE_STRAND = strand
    )
  subject <- df_ranges[subjectHits(intersection)] %>%
    as_tibble() %>%
    transmute(
      WINDOW_CHROM = seqnames,
      WINDOW_START = start,
      WINDOW_END = end,
      N_VARIANTS,
      VAL,
      EMPIRICAL_PVALUE
    )
  nonsubject <- df_ranges[-unique(subjectHits(intersection))] %>%
    as_tibble() %>%
    transmute(
      GENE_ID = as.character(NA),
      GENE = as.character(NA),
      GENE_CHROM = as.character(NA),
      GENE_START = as.integer(NA),
      GENE_END = as.integer(NA),
      GENE_STRAND = as.character(NA),
      WINDOW_CHROM = seqnames,
      WINDOW_START = start,
      WINDOW_END = end,
      N_VARIANTS,
      VAL,
      EMPIRICAL_PVALUE,
    )

  hits <- bind_cols(query, subject)

  windows_summary <- bind_rows(hits, nonsubject) %>%
    group_by(WINDOW_CHROM, WINDOW_START, WINDOW_END) %>%
    summarise(
      GENEID = str_c(str_replace_na(GENE_ID, " "), collapse = ","),
      GENE = str_c(str_replace_na(GENE, " "), collapse = ","),
      VAL = dplyr::first(VAL),
      EMPIRICAL_PVALUE = dplyr::first(EMPIRICAL_PVALUE),
      N_VARIANTS = dplyr::first(N_VARIANTS),
      .groups = "drop"
    ) %>%
    ungroup() %>%
    transmute(
      CHROM = WINDOW_CHROM,
      BIN_START = WINDOW_START,
      BIN_END = WINDOW_END,
      N_VARIANTS,
      !!sym(val_column) := VAL,
      EMPIRICAL_PVALUE,
      GENEID,
      GENE
    )

  results <- list(
    hits = hits,
    summary = windows_summary
  )
  return(results)
}

heatmap_chr <- function(df, chr, group, val, out, legend = "") {
  mean_agg <- df %>%
    group_by(!!sym(chr), !!sym(group)) %>%
    summarise(!!sym(val) := mean(!!sym(val), na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(!!sym(chr) := as.factor(!!sym(chr)))

  p <- ggplot(mean_agg, aes(x = !!sym(chr), y = !!sym(group), fill = !!sym(val))) +
    geom_tile() +
    white_theme +
    theme_minimal() +
    scale_fill_distiller(palette = "RdYlBu", direction = -1, name = legend) +
    theme(
      axis.text = element_text(angle = 0, size = 7),
      legend.text = element_text(size = 7),
      axis.text.y = element_text(angle = 0, size = 6),
      axis.text.x = element_text(angle = 45, size = 6),
      axis.title.y = element_blank()
    ) +
    xlab("Chromosome") +
    ylab(group)
  ggsave(out, p, width = 1920, height = 1080, units = "px", dpi = 300)
}

heatmap <- function(df, group, val, out, legend = "") {
  mean_agg <- df %>%
    group_by(!!sym(group)) %>%
    summarise(!!sym(val) := mean(!!sym(val), na.rm = TRUE)) %>%
    separate(!!sym(group), into = c("Population1", "Population2"), sep = " vs ")

  full_tbl <- mean_agg %>%
    mutate(
      temp = Population2,
      Population2 = Population1,
      Population1 = temp
    ) %>%
    dplyr::select(-temp) %>%
    bind_rows(mean_agg) %>%
    bind_rows(
      tibble(
        Population1 = unique(c(mean_agg$Population1, mean_agg$Population2)),
        Population2 = Population1,
        !!sym(val) := 0,
      )
    ) %>%
    mutate(
      pop1 = as.numeric(factor(Population1)),
      pop2 = as.numeric(factor(Population2))
    ) %>%
    mutate(!!sym(val) := if_else(pop1 >= pop2, !!sym(val), NA))

  p <- ggplot(full_tbl, aes(x = Population1, y = Population2, fill = !!sym(val))) +
    geom_tile() +
    geom_text(aes(label = ifelse(is.na(!!sym(val)), "", sprintf("%.2f", !!sym(val)))), na.rm = TRUE, color = "black") +
    white_theme +
    scale_fill_distiller(palette = "RdYlBu", direction = -1, name = legend) +
    theme(
      axis.title = element_blank(),
      axis.text = element_text(angle = 0, size = 12, face = "bold"),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 8),
      axis.text.x = element_text(angle = 30, face = "bold"),
      axis.text.y = element_text(face = "bold")
    )
  ggsave(out, p, width = 1920, height = 1080, units = "px", dpi = 300)
}

freq_plot <- function(df, column, out) {
  counts <- df %>%
    group_by(!!sym(column)) %>%
    summarise(COUNT = n()) %>%
    complete(!!sym(column), fill = list(COUNT = 0))
  p <- ggplot(counts, aes(x = !!sym(column), y = COUNT)) +
    geom_bar(stat = "identity") +
    white_theme +
    xlab("Chromosome") +
    ylab("Number of Sweeps")
  ggsave(out, p, width = 1920, height = 1080, units = "px")
}


within_sweep_distribution <- function(sweeps, raw, chromosomes, val, ylab, out, both_tails = TRUE) {
  require(tidygenomics)
  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)

  raw <- raw %>%
    mutate(
      CHROM = factor(CHR, levels = chromosomes),
      BIN_START = POSITION,
      BIN_END = POSITION
    )

  within_sweeps <- sweeps %>%
    genome_join_closest(raw, by = c("CHROM", "BIN_START", "BIN_END"), max_distance = 1) %>%
    transmute(CHROM = CHROM.x, START = BIN_START.x, END = BIN_END.x, POSITION = BIN_START.y, !!sym(val), MEAN)

  ylims <- range(raw %>% pull(!!sym(val)), na.rm = T)
  res <- within_sweeps %>%
    group_by(CHROM, START, END) %>%
    group_walk(~ {
      mean <- if (both_tails) substitute(group("|", NAME, "|")[w], list(NAME = as.name(ylab))) else substitute(NAME[w], list(NAME = as.name(ylab)))
      p <- ggplot(.x, aes(x = POSITION / 1e6, y = !!sym(val))) +
        geom_point(alpha = 0.5, size = 1) +
        ggtitle(bquote("Selective sweep at chr" ~ .(.y$CHROM) * ":" * .(.y$START / 1e6) * "-" *
          .(.y$END / 1e6) * " Mb with" ~ .(mean) == .(round(.x$MEAN, 2)))) +
        ylab(ylab) +
        xlab("Position (Mb)") +
        ylim(ylims) +
        xlim(.y$START / 1e6, .y$END / 1e6) +
        white_theme +
        theme(plot.title = element_text(size = 8))

      ggsave(str_glue("{out}_chr{.y$CHROM}_{.y$START}-{.y$END}.png"), p, width = 1920, height = 1080, units = "px")
    })
}


find_clusters <- function(df, chr = "CHROM", start = "BIN_START", end = "BIN_END", val_column = "MEAN",
                          threshold = 0.01, max_distance = 4e4 + 1, min_windows = 5) {
  require(tidygenomics)

  res <- df %>%
    slice_max(!!sym(val_column), prop = threshold) %>%
    genome_cluster(by = c(chr, start, end), max_distance = max_distance) %>%
    group_by(cluster_id) %>%
    summarise(CHROM = dplyr::first(!!sym(chr)), START = min(!!sym(start)), END = max(!!sym(end)), N_WINDOWS = n()) %>%
    arrange(desc(N_WINDOWS)) %>%
    filter(N_WINDOWS >= min_windows)
  return(res)
}
