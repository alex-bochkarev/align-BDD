######################################################################
## Generates stats for the set of optima (simplified problem)
##
## (c) Alexey Bochkarev, Clemson University, 2020

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(optparse)
})

######################################################################
## unpack the command line arguments
option_list = list(
  make_option(c("-d", "--outdir"), type="character", default=null,
              help="output directory to create files", metavar="character"),
  make_option(c("-i", "--input"), type="character", default=null,
              help="input file (run log)", metavar="character"));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
infile <- opt$input

outfile_no_opts <- paste0(opt$outdir, "/figures/no_opts.eps")
outfile_opts_diam <- paste0(opt$outdir, "/figures/opts_diam.eps")
outfile_simscore_vs <- paste0(opt$outdir, "/figures/heuristic_simscore.eps")
outfile_bin2d <- paste0(opt$outdir, "/figures/heuristic_simscore_vs_AB_simscore.eps")
######################################################################
## Number of optima

cat(sprintf("Fig: stats on optima"), "\n")
cat(sprintf("Input file: %s", infile),'\n')
df = read.csv(infile)

cat(sprintf(
  "Instances read: %d, showing %d with no more than 100 optima.",
  nrow(df), nrow(filter(df, no_opts <= 100))
), "\n")

df_less_than_hundred = filter(df, no_opts <= 100)

p_opts =
ggplot(df_less_than_hundred) +
  geom_histogram(aes(x = no_opts), binwidth = 2) +
  xlab("Number of optimal solutions") +
  ylab(paste0("Count of instances (out of ", nrow(df), ")"))+
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 23),
    axis.title.y = element_text(size = 23),
    panel.background = element_blank(),
    panel.grid.major = element_line(
      size = 0.5, linetype = "solid",
      colour = "lightgrey"
    ),
    panel.grid.minor.x = element_line(size = 0.25, linetype = "solid", colour = "lightgrey"),
    panel.grid.minor.y = element_line(size = 0.25, linetype = "solid", colour = "lightgrey")
  )+
  scale_x_continuous(
    breaks = seq(0,100,by = 10),
  )

ggsave(outfile_no_opts, p_opts, width=10, height=10)

cat(sprintf("Optima stats:"),'\n')
cat(sprintf("===================="), "\n")
cat(sprintf("Out of %d instances,", nrow(df)), "\n")
cat(sprintf("More than one optimum: %d", nrow(filter(df, no_opts > 1))), "\n")
cat(sprintf("More than five: %d", nrow(filter(df, no_opts > 5))), "\n")
cat(sprintf("More than fifty: %d", nrow(filter(df, no_opts > 50))), "\n")

######################################################################
## Optima set diameter
cat(sprintf("********************"), "\n")
cat(sprintf("Fig: opts diam"), "\n")

ggplot(df) +
  geom_histogram(aes(x=opts_diam), binwidth=5) +
  xlab("Min. simscore between optima, percent") +
  ylab(paste0("Count of instances (out of ", nrow(df), ")"))+
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 23),
    axis.title.y = element_text(size = 23),
    panel.background = element_blank(),
    panel.grid.major = element_line(
      size = 0.5, linetype = "solid",
      colour = "lightgrey"
    ),
    panel.grid.minor.x = element_line(size = 0.25, linetype = "solid", colour = "lightgrey"),
    panel.grid.minor.y = element_line(size = 0.25, linetype = "solid", colour = "lightgrey")
  ) +
  scale_x_continuous(
      breaks = seq(0, 100, by = 10),
      )

ggsave(outfile_opts_diam, p_opts, width=10, height=10)

cat(sprintf("Opts simscore stats:"), "\n")
cat(sprintf("===================="), "\n")
cat(sprintf("Out of %d instances,", nrow(df)), "\n")
cat(sprintf("Min simscore <= 25: %d (%.1f %%)",
            nrow(filter(df, opts_diam <= 25)),
            nrow(filter(df, opts_diam <= 25)) * 100 / nrow(df)),
    "\n")

######################################################################
## Best / worst simscores for VS
cat(sprintf("********************"), "\n")
cat(sprintf("Fig: varseq vs. optima, simscores."), "\n")

df_m = pivot_longer(select(df,
                           best_VS_simscore,
                           worst_VS_simscore,k),
                    names_to = "score_type",
                    cols=c('best_VS_simscore',
                           'worst_VS_simscore'))
df_m$labels=with(df_m,
                 ifelse(score_type=='best_VS_simscore', 'Best', 'Worst'))

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

p_simscore=
ggplot(df_m) +
  geom_histogram(aes(x = value)) +
  xlab("Similarity score (with an optimum)") +
  ylab(paste0("Count of instances (out of ", nrow(df), ")"))+
  guides(fill=guide_legend("Simscore:"))+
  theme(
    legend.position = c(0.2, 0.8),
    legend.direction = "vertical",
    legend.title = element_text(size=18),
    legend.text = element_text(size=18),
    legend.text.align = 0,
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 23),
    axis.title.y = element_text(size = 23),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 18, face="bold"),
    panel.background = element_blank(),
    panel.grid.major = element_line(
      size = 0.5, linetype = "solid",
      colour = "lightgrey"
    ),
    panel.grid.minor.x = element_line(size = 0.25, linetype = "solid", colour = "lightgrey"),
    panel.grid.minor.y = element_line(size = 0.25, linetype = "solid", colour = "lightgrey")
  ) +
  scale_x_continuous(
    breaks = seq(0, 100, by = 10),
    )+
  facet_wrap(labels~.)+
  coord_flip()

ggsave(outfile_simscore_vs, p_simscore, width = 16, height = 10)

p_bin2d =
ggplot(df) +
  geom_bin2d(aes(AB_simscore, best_VS_simscore), bins = 10) +
  scale_fill_viridis_c() +
  guides(fill = guide_legend("No. instances:")) +
  xlab("Simscore(A,B): initial orders")+
  ylab("Best simscore (heuristic)")+
  theme(
    ## legend.position = c(0.2, 0.8),
    ## legend.direction = "vertical",
    legend.title = element_text(size=18),
    legend.text = element_text(size=18),
    ## legend.text.align = 0,
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 23),
    axis.title.y = element_text(size = 23),
    panel.background = element_blank(),
    panel.grid.major = element_line(
      size = 0.5, linetype = "solid",
      colour = "lightgrey"
    ),
    panel.grid.minor.x = element_line(size = 0.25, linetype = "solid", colour = "lightgrey"),
    panel.grid.minor.y = element_line(size = 0.25, linetype = "solid", colour = "lightgrey")
  )

ggsave(outfile_bin2d, p_bin2d, width = 16, height = 10)

## str(df)
## sum(df$opts_diam < 25)

## str(df)

## sum(df$best_VS_simscore > 75)
