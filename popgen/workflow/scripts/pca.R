## PCA for (pruned) populations

suppressPackageStartupMessages({
  library(SNPRelate)
  library(tidyverse)
  library(GGally)
  library(argparser)
  library(fs)
  source("workflow/scripts/util.R")
})



parse_input <- function() {
  p <- arg_parser("Visualise PCA for (pruned) populations")
  p <- add_argument(p, "--file", help = "VCF or GDS file")
  p <- add_argument(p, "--out", help = "Output file prefix")
  args <- parse_args(p)
  if (is.nall(args$file)) {
    stop("--file should be a path to a VCF or GDS file")
  }
  if (is.nall(args$out)) {
    stop("--out should include an output name")
  }

  geno <- read_input(args$file)
  return(list(args = args, geno = geno))
}

read_input <- function(file) {
  basename <- path_ext_remove(file)
  ext <- str_match(file, "\\.([^/]+)$")[, 2]
  if (ext == "vcf" || ext == "vcf.gz") {
    snpgdsVCF2GDS(file, str_glue("{basename}.gds"), method = "biallelic.only")
    file <- str_glue("{basename}.gds")
  } else if (ext != "gds") {
    stop("file is not in VCF or gds format")
  }
  geno <- snpgdsOpen(file)
  return(geno)
}

visualise_pca <- function(pca_res, out) {
  ids <- str_split(pca_res$sample.id, pattern = "_") %>%
    lapply(\(x) x[[length(x)]]) %>%
    unlist()

  pops <- str_split(pca_res$sample.id, pattern = "_") %>%
    lapply(\(x) x[1:(length(x) - 1)]) %>%
    lapply(\(x) paste(unlist(x), collapse = " ")) %>%
    unlist()

  pca_df <- tibble(ID = ids, Population = pops) %>%
    bind_cols(pca_res$eigenvect[, 1:5])

  pcs <- str_glue("PC{1:5}")
  colnames(pca_df) <- c("ID", "Population", pcs)
  var_exp <- round(pca_res$varprop * 100, 2)
  var_exp <- var_exp[!is.na(var_exp)]
  eigenval <- pca_res$eigenval[!is.na(pca_res$eigenval)]
  eigenvect <- as.data.frame(pca_res$eigenvect)
  write_csv(
    tibble(Component = str_glue("PC{1:length(var_exp)}"), Variance = var_exp, Eigenvalue = eigenval),
    str_glue("{out}variance.csv")
  )
  write_csv(eigenvect, str_glue("{out}eigenvectors.csv"))

  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Population)) +
    geom_point(size = 1.5, alpha = 0.5) +
    labs(x = paste0("PC1 (", var_exp[1], "%)"), y = paste0("PC2 (", var_exp[2], "%)")) +
    white_theme +
    theme(
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 8)
    )
  ggsave(str_glue("{out}pca.png"), p, units = "px", height = 1080, width = 1920, dpi = 300)

  p <- ggduo(
    pca_df,
    aes(color = Population),
    columnsX = pcs,
    columnsY = pcs,
    types = list(
      continuous = wrap("points", size = 0.5, alpha = 0.5)
    ),
    legend = c(5, 3)
  ) +
    white_theme +
    theme(
      legend.position = "bottom",
      axis.title = element_text(size = 8, face = "bold"),
      axis.text = element_text(size = 6),
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 7)
    )
  ggsave(str_glue("{out}pca_top5.png"), p, units = "px", height = 1080, width = 1920, dpi = 300)
}

input <- parse_input()
pca_res <- snpgdsPCA(input$geno, num.thread = 8)
visualise_pca(pca_res, input$args$out)
