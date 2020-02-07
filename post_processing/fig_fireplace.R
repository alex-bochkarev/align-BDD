######################################################################
## Creates ggplot for the simplified problem:
## Calculation wall-clock time vs. objective for different
## heuristics
##
## (c) Alexey Bochkarev, Clemson University, 2020

library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(stringr)
library(latex2exp)
library(optparse)

## Internal parameters for the figure
ORIG_BASE_COL = "simpl_BB_obj" # col to divide by
Nticks = 20
Y_QUANTILE=0.75

## Heuristics to show (codes)
SHOW_HEU = c("all")

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

df = read.csv(opt$input, stringsAsFactors = FALSE)

df_legend = select(filter(df, num_type=="legend"), value,comment)
df = filter(df, num_type != "legend")
df$value = as.numeric(df$value)

df_wide = pivot_wider(df,
                      id_cols = "instance",
                      names_from = "num_type",
                      values_from = "value"
                      )

simpl_obj_cols = grep("^simpl.*_obj$",colnames(df_wide),value = TRUE)
time_cols = grep(".+time$",colnames(df_wide),value = TRUE)
######################################################################
## draw the figure

for (col in simpl_obj_cols){
    df_wide[[paste(col,"rel",sep="_")]] = df_wide[[col]] / df_wide[[ORIG_BASE_COL]]
}

df_rel = pivot_longer(df_wide,
                      cols = -c("instance"),
                      names_to = "num_type",
                      values_to = "value"
                      )

if (SHOW_HEU[1] == "all"){
    SHOW_HEU = c(grep("^simpl.*_obj_rel$",colnames(df_wide),value = TRUE),
                 grep("^simpl.*_time$",colnames(df_wide),value = TRUE))
}

df_rel = filter(df_rel, num_type %in% SHOW_HEU)

df_rel$entry_type = ifelse(df_rel$num_type %in% time_cols,"time","obj")

df_rel = df_rel %>%
    mutate(
        heuristic = ifelse(entry_type == "obj",
                           substring(num_type, 1, nchar(num_type)-nchar("_obj_rel")),
                           substring(num_type, 1, nchar(num_type)-nchar("_time")))
    )

df_time_o = pivot_wider(select(df_rel,-num_type),names_from = "entry_type",values_from = "value")
df_time_o = merge(x=df_time_o, y=df_legend, by.x = "heuristic", by.y = "value")

ymin = 1.0
ymax = quantile(df_time_o$obj, Y_QUANTILE)
xmin = 0
xmax = max(df_time_o$time)
plt_zoomed =
    ggplot(df_time_o,aes(x=time, y=obj , shape=comment,color=comment))+
    geom_point(alpha=0.4, size=4)+
    scale_y_continuous(
        "Objective value (percent of the exact min)",
        limits = c(ymin,ymax),
        labels = scales::percent
    )+
    scale_x_continuous(
        "Calculation (wall-clock) time per instance, sec.",
        labels = scales::number_format(accuracy = 0.1),
        breaks = seq(xmin,xmax,length.out = Nticks)
    )+
    guides(
        color=guide_legend(title="Heuristic / method:"),
        shape=guide_legend(title="Heuristic / method:")
        )+
    theme(
        legend.position = c(0.8, 0.8),
        legend.direction = "vertical",
        legend.title = element_text(size=24),
        legend.text = element_text(size=24),
        axis.text.x = element_text(size=22,angle=45,vjust = 0.7),
        axis.text.y = element_text(size=22),
        axis.title.x = element_text(size = 26),
        axis.title.y = element_text(size = 26, margin = margin(t=50)),
        panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        color = "darkgrey"),
          )

ggsave(opt$out,plt_zoomed,width = 16, height = 10, device = cairo_ps)
