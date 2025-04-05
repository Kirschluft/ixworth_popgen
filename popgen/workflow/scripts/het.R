# Assessment of heterozygosity and inbreeding from VCFtools

suppressPackageStartupMessages({
  library(tidyverse)
  library(argparser)
  library(parallel)
  source("workflow/scripts/util.R")
})



parse_input <- function() {
  p <- arg_parser("Assess and visualise heterozygosity and inbreeding (F)")
  p <- add_argument(p, "--files", nargs = Inf, help = "Input het files from vcftools")
  p <- add_argument(p, "--cores", default = 16, help = "Maximum number of cores to run")
  p <- add_argument(p, "--out", help = "Output file prefix")
  args <- parse_args(p)
  if (is.nall(args$files)) {
    stop("--files should include a list of file names")
  }
  if (is.nall(args$out)) {
    stop("--out should include an output name")
  }

  dfs <- mclapply(args$files, \(x) read_input_het(x, "het", column = "POPULATION"), mc.cores = args$cores)
  return(list(args = args, dfs = dfs))
}

read_input_het <- function(x, ext, column) {
  name <- get_name(x, ext)
  df <- read_tsv(x, show_col_types = FALSE) %>%
    mutate(!!sym(column) := name) %>%
    mutate(
      E_HET = 1 - `E(HOM)` / N_SITES,
      O_HET = 1 - `O(HOM)` / N_SITES
    )
  return(df)
}


input <- parse_input()
df <- do.call(bind_rows, input$dfs)

summary <- df %>%
  group_by(POPULATION) %>%
  summarise(
    MEAN_F = round(mean(F), 2),
    SD_F = round(sd(F), 2),
    MEAN_E_HOM = mean(`E(HOM)`),
    MEAN_O_HOM = mean(`O(HOM)`),
    MEAN_E_HOM = mean(E_HET),
    MEAN_O_HET = mean(O_HET)
  )
write_csv(summary, str_glue("{input$args$out}_summary.csv"))


f_dist <- ggplot(df, aes(x = F)) +
  geom_density() +
  white_theme +
  facet_wrap(~POPULATION) +
  ylab("Density") +
  xlab("Inbreeding Coefficient (F)") +
  xlim(-1, 1)
ggsave(str_glue("{input$args$out}_F_density.png"), f_dist, width = 1920, height = 1080, units = "px")

f_dist <- ggplot(df, aes(x = F, col = POPULATION), alpha = 0.5) +
  geom_density() +
  white_theme +
  ylab("Density") +
  xlab("Inbreeding Coefficient (F)") +
  xlim(-1, 1)
ggsave(str_glue("{input$args$out}_F_density_color.png"), f_dist, width = 1920, height = 1080, units = "px")

f_violin <- ggplot(df, aes(y = F, x = POPULATION)) +
  geom_violin(fill = "gray35", adjust = 0.5, draw_quantiles = 0.5) +
  ylab("Inbreeding Coefficient (F)") +
  white_theme
ggsave(str_glue("{input$args$out}_F_violinplot.png"), f_violin, width = 1920, height = 1080, units = "px")
