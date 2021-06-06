######################################################################
## Creates ggplot for the original problem (several versions):
## estimates the efficiency of the simplified-problem based heuristic
## as compared to the naive `minAB` heuristic.
##
## (c) Alexey Bochkarev, Clemson University, 2021

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(latex2exp)
  library(optparse)
})

######################################################################
## unpack the command line arguments
option_list = list(
  make_option(c("-o", "--out"), type="character", default="./out.eps",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-i", "--indir"), type="character", default="./out.eps",
              help="directory for input files", metavar="character"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$out)){
    print_help(opt_parser)
    stop("Please specify output file", call.=FALSE)
}

if (is.null(opt$indir)){
  print_help(opt_parser)
  stop("Please specify the input directory", call.=FALSE)
}

######################################################################
## parse the input file
indir = opt$indir
df = read.csv(paste0(indir,"/tUFLP_nat.csv"), stringsAsFactors = FALSE)
df$inst_type = "(a) t-UFLP (natural order)"
cat("Number of t-UFLP (nat) instances:", nrow(df)/2,"\n")

df_rnd <- read.csv(paste0(indir, "/tUFLP_rnd.csv"), stringsAsFactors = FALSE)
df_rnd$inst_type = "(b) t-UFLP (random order)"
cat("Number of t-UFLP (rnd) instances:", nrow(df_rnd)/2,"\n")

dfjUFL <- read.csv(paste0(indir, "/jUFLP.csv"), stringsAsFactors = FALSE)
dfjUFL$inst_type <- "(c) Joint UFLP"
cat("Number of j-UFLP instances:", nrow(dfjUFL)/2,"\n")

df_rnd_dia <- read.csv(paste0(indir, "/rnd_dia.csv"), stringsAsFactors = FALSE)
df_rnd_dia$inst_type <- "(d) Random diagrams"
cat("Number of rnd-dia instances:", nrow(df_rnd_dia)/2,"\n")

df = rbind(df, df_rnd, dfjUFL, df_rnd_dia)

df_wide = pivot_wider(df,
                        id_cols = c("instance", "inst_type"),
                        names_from = "num_type",
                        values_from = "value"
                        )

no_insts = nrow(filter(df_wide, inst_type=="(a) t-UFLP (natural order)"))
cat("Number used for the title: ", no_insts,"\n")

df_wide$simpl_rel_obj = with(df_wide, orig_simpl_obj / orig_minAB_obj)

problems = c("(a) t-UFLP (natural order)", "(b) t-UFLP (random order)",
             "(c) Joint UFLP", "(d) Random diagrams")

df_wide$problem = factor(df_wide$inst_type, levels=problems)

plt_dens =
    ggplot(df_wide, aes(x=simpl_rel_obj))+
    geom_histogram(binwidth = 0.02,position = "identity")+
    geom_vline(xintercept = 1.0, size=0.5, color="red", linetype="dashed")+
    ## styling
    scale_y_continuous(
        paste0("Number of instances, out of ", no_insts),
        labels = scales::number_format(accuracy = 1),
        position = "right"
    )+
    scale_x_continuous(
        "Itersection diagram size (relative to 'Best of A and B')",
        labels = scales::percent,
        limits = c(0.0, 2.5)
    )+
  theme(
    legend.position = c(0.6, 0.8),
    legend.direction = "vertical",
    legend.title = element_text(size=24),
    legend.text = element_text(size=24),
    legend.text.align = 0,
    axis.text.x = element_text(size=15,angle=45,vjust = 0.7),
    axis.text.y = element_text(size=15),
    axis.title.x = element_text(size = 23),
    axis.title.y = element_text(size = 23),
    panel.background = element_blank(),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgrey"),
    panel.grid.minor.x = element_line(size=0.5,linetype = 'solid', colour = "lightgrey"),
    panel.grid.minor.y = element_line(size=0.5,linetype = 'solid', colour = "lightgrey"),
    strip.text.y.left = element_text(size=15, angle=0),
    strip.background = element_blank(),
    strip.text.x.bottom = element_text(size=20, angle=0),
    )+
  facet_grid(. ~ problem, switch="x")+
  coord_flip()
    ## end of styling

cat("************************\n")
cat("Performance summary:\n")
df_wide %>%
  group_by(inst_type) %>%
  summarize(
    total=n(),
    better_than_minAB = sum(simpl_rel_obj < 1.0),
    better_than_minAB_percent = sum(simpl_rel_obj < 1.0)/n()
    )

ggsave(opt$out, plt_dens, width = 16, height = 10)
