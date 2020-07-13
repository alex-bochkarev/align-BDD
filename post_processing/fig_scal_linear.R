######################################################################
##
## Scalability figures (from the ScalTest - generated instances)
##
## (c) Alexey Bochkarev, Clemson University, 2020
## abochka@clemson.edu

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(optparse)

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
df$comment = NULL
df = filter(df,N>0)
df$value = as.numeric(df$value)

df = pivot_wider(df, id_cols = c("instance","N"), names_from = num_type, values_from = value)

df = df %>%
    mutate(
        rel_obj = orig_simpl_obj / orig_gsifts1p_obj
    )

## calculate means for each instance size
df_means = df %>%
  group_by(N) %>%
  summarise(
    mean_simpl_time = median(orig_simpl_time),
    mean_gsifts_time = median(orig_gsifts1p_time)
  )

## percentage where aux problem outperforms
df_outp = df
df_outp$aux_outp = with(df_outp,
                           orig_simpl_time < orig_gsifts1p_time & rel_obj < 1.0
                           )

df_outp = df_outp %>%
  group_by(N) %>%
  summarize(
    simpl_outp = sum(aux_outp) / length(aux_outp)
  )

p1 =
ggplot(df)+
    geom_point(aes(x=N, y = orig_simpl_time, color="Auxiliary problem / heuristic"),
             shape=3, size=5, alpha=0.3)+
    geom_jitter(aes(x=N, y = orig_gsifts1p_time, color="BDD sifts"),
              size=2, alpha=0.3, width = 0.1)+
    geom_line(data=df_means, aes(x=N, y=mean_simpl_time,
                               color="Auxiliary problem / heuristic", linetype="Auxiliary problem / heuristic"), size=2)+
    geom_line(data=df_means, aes(x=N, y=mean_gsifts_time,
                               color="BDD sifts", linetype="BDD sifts"), size=2)+
    geom_label(data=df_outp, aes(x=N, label=scales::percent(accuracy=1,simpl_outp)),y=1550, size=5)+
    scale_linetype_manual(values = c("solid", "dashed"),
                        labels = c("Auxiliary problem / heuristic", "BDD sifts"))+
    scale_x_continuous(
        breaks = seq(min(df$N),max(df$N),by = 1),
        minor_breaks = seq(min(df$N),max(df$N), by=1)
    )+
    ## scale_y_continuous(limits = c(0, 100))+
    labs(x="No. of variables\n(a)",y="Solution time, sec., linear scale, t")+
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

ggsave("../tmp/fig_scal_linear_all.png",p1,width = 16, height = 10)
