######################################################################
## Creates ggplot for the original problem:
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
                help="output file name [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$out)){
    print_help(opt_parser)
    stop("Please specify output file", call.=FALSE)
}

######################################################################
## parse the input file
df = read.csv("./run_logs/jUFL/logfile_tUFL_nat.csv", stringsAsFactors = FALSE)
df$inst_type = "(a) tUFL (natural order)"

df_rnd <- read.csv("./run_logs/jUFL/logfile_tUFL_rnd.csv", stringsAsFactors = FALSE)
df_rnd$inst_type = "(b) tUFL (random order)"

dfjUFL <- read.csv("./run_logs/jUFL/logfile.csv", stringsAsFactors = FALSE)
dfjUFL$inst_type <- "(c) jUFL"


df_rnd_dia <- read.csv("./run_logs/jUFL/logfile_rnd_dia.csv", stringsAsFactors = FALSE)
df_rnd_dia$inst_type <- "(d) Random diagrams"

df = rbind(df, df_rnd, dfjUFL, df_rnd_dia)

df_wide = pivot_wider(filter(df, n==20),
                        id_cols = c("instance", "inst_type", "n", "prob"),
                        names_from = "num_type",
                        values_from = "value"
                        )

df_wide$simpl_rel_obj = with(df_wide, orig_simpl_obj / orig_minAB_obj)

problems = c("(a) tUFL (natural order)", "(b) tUFL (random order)",
             "(c) jUFL", "(d) Random diagrams")

df_wide$problem = factor(df_wide$inst_type, levels=problems)

plt_dens =
    ggplot(df_wide, aes(x=simpl_rel_obj))+
    geom_histogram(binwidth = 0.02,position = "identity")+
    geom_vline(xintercept = 1.0, size=0.5, color="red", linetype="dashed")+
    ## styling
    scale_y_continuous(
        "No. of instances",
        labels = scales::number_format(accuracy = 1),
        position = "right"
    )+
    scale_x_continuous(
        "Itersection diagram size (relative)",
        labels = scales::percent,
        ## breaks = seq(xmin,xmax,length.out=11),
        ## minor_breaks = seq(xmin,xmax,length.out=21),
        limits = c(0.0, 2.5)
    )+
  ggtitle("Intersection diagram size of the proposed heuristic relative to the minAB heuristic.")+
  theme(
    legend.position = c(0.6, 0.8),
    legend.direction = "vertical",
    legend.title = element_text(size=24),
    legend.text = element_text(size=24),
    legend.text.align = 0,
    axis.text.x = element_text(size=15,angle=45,vjust = 0.7),
    axis.text.y = element_text(size=15),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    panel.background = element_blank(),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgrey"),
    panel.grid.minor.x = element_line(size=0.5,linetype = 'solid', colour = "lightgrey"),
    panel.grid.minor.y = element_line(size=0.5,linetype = 'solid', colour = "lightgrey"),
    strip.text.y.left = element_text(size=15, angle=0),
    strip.background = element_blank(),
    strip.text.x.top = element_text(size=15, angle=0),
    )+
  facet_grid(. ~ problem, switch="y")+
  coord_flip()
    ## end of styling

df_wide %>%
  group_by(inst_type) %>%
  summarize(
    total=n(),
    better_than_minAB = sum(simpl_rel_obj < 1.0),
    better_than_minAB_percent = sum(simpl_rel_obj < 1.0)/n()
    )

ggsave(opt$out, plt_dens, width = 16, height = 10)
