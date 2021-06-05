######################################################################
##
## cUFL runtimes overview
##
## (c) Alexey Bochkarev, Clemson University, 2021
## abochka@clemson.edu

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(optparse)
library(viridis)

######################################################################
## unpack the command line arguments
option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL,
                help="input filename (solved instances log)", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="./out.eps",
                help="output file name [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input) | is.null(opt$out)){
    print_help(opt_parser)
    stop("Please specify both input and output files", call.=FALSE)
}

######################################################################
## parse the input file
infile = opt$input

df = read.csv(infile, stringsAsFactors = FALSE)
df$method = NULL
df$run = NULL

df$step = sapply(df$step, function(x) gsub("[+-]",".", x), USE.NAMES=FALSE)

df = pivot_wider(df, id_cols = c("k","n"), names_from = step, values_from = time)

df = df %>%
    mutate(
      naive_MIP = build.solve,
      BDD_MIP = BDD.build + MIP.build + MIP.solve,
      align_vs_SP = BDD.build+VS.build.solve + BDD.align.to.vs + intersection.build + intersection.SP.solve,
      align_gsifts_SP = BDD.build + BDD.align.gsifts + intersection.gsifts.build + intersection.gsifts.SP.solve
    )

## calculate means for each instance size
df_means = df %>%
  group_by(n) %>%
  summarise(
    med_naive_MIP = median(naive_MIP),
    med_BDD_MIP = median(BDD_MIP),
    med_align_vs_SP = median(align_vs_SP),
    med_align_gsifts_SP = median(align_gsifts_SP)
  )

cairo_ps(opt$out, width = 16, height = 10, family="Arial")

ggplot(df)+
  scale_color_viridis(discrete=TRUE) +
  scale_y_log10(labels=scales::comma)+
  annotation_logticks(sides='l') +
  geom_jitter(aes(x=n, y = naive_MIP, color="Naive MIP"),
              size=2, alpha=0.3, width = 0.1, shape=1)+
  geom_jitter(aes(x=n, y = BDD_MIP, color="CPP MIP"),
              size=2, alpha=0.3, width = 0.1, shape=2)+
  geom_jitter(aes(x=n, y = align_vs_SP, color="CPP align-BDD (VS) + SP"),
              size=2, alpha=0.3, width = 0.1, shape=3)+
  geom_jitter(aes(x=n, y = align_gsifts_SP, color="CPP align-BDD (gsifts) + SP"),
              size=2, alpha=0.3, width = 0.1, shape=4)+
  labs(x="No. of variables (nodes in the original graph)",y="Solution time, msec")+
  geom_line(data=df_means, aes(x=n, y=med_naive_MIP,
                               color="Naive MIP", linetype="Naive MIP"), size=2)+
  geom_line(data=df_means, aes(x=n, y=med_BDD_MIP,
                               color="CPP MIP", linetype="CPP MIP"), size=2)+
  geom_line(data=df_means, aes(x=n, y=med_align_vs_SP,
                                color="CPP align-BDD (VS) + SP", linetype="CPP align-BDD (VS) + SP"), size=2)+
  geom_line(data=df_means, aes(x=n, y=med_align_gsifts_SP,
                               color="CPP align-BDD (gsifts) + SP", linetype="CPP align-BDD (gsifts) + SP"), size=2)+
  scale_x_continuous(
    breaks = seq(min(df$n),max(df$n),by = 1),
    minor_breaks = seq(min(df$n),max(df$n), by=1)
  )+
  guides(color=guide_legend("Solution method (median values):"), linetype=guide_legend("Solution method (median values):"))+
  theme(
    legend.position = c(0.2, 0.8),
    legend.direction = "vertical",
    legend.title = element_text(size=16),
    legend.text = element_text(size=16),
    legend.key = element_blank(),
    legend.key.width = unit(2.5,"cm"),
    axis.text.x = element_text(size=18,angle=45,vjust = 0.7),
    axis.text.y = element_text(size=18),
    axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26),
    panel.background = element_blank(),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgrey"),
    panel.grid.minor.x = element_line(size=0.25,linetype = 'solid', colour = "lightgrey"),
    panel.grid.minor.y = element_line(size=0.25,linetype = 'solid', colour = "lightgrey")
  )

dev.off()
