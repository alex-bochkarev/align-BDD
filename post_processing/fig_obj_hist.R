######################################################################
## Creates ggplot for the original problem:
## benchmarking different heuristics
## (objective values histogram)
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
X_QUANTILE=0.98 # quantile to filter for the histogram(Ox axis)
ORIG_BASE_COL = "orig_gsifts1p_obj" # col to divide by

## Heuristics to show (codes)
SHOW_HEU = c("all")
# SHOW_HEU = c("orig_simpl_obj_rel","orig_bestAB_obj_rel","orig_simpl_rnd_obj_rel")

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

orig_obj_cols = grep("^orig.*_obj$",colnames(df_wide),value = TRUE)
time_cols = grep(".+time$",colnames(df_wide),value = TRUE)

for (col in orig_obj_cols){
    df_wide[[paste(col,"rel",sep="_")]] = df_wide[[col]] / df_wide[[ORIG_BASE_COL]]
}

df_rel = pivot_longer(df_wide,
                      cols = -c("instance"),
                      names_to = "num_type",
                      values_to = "value"
                      )

if (SHOW_HEU[1] == "all"){
    SHOW_HEU = grep("^orig.*_obj_rel$",colnames(df_wide),value = TRUE)
    SHOW_HEU = SHOW_HEU[ SHOW_HEU != paste(ORIG_BASE_COL,"_rel",sep="")]
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

######################################################################
## draw the figure
# X_QUANTILE = 0.98
xmin = min(df_time_o$obj)
xmax = quantile(df_time_o$obj,X_QUANTILE)
latex_label = parse(text = TeX("try $S_A$, $S_B$, choose the best one"))

plt_dens =
    ggplot(df_time_o, aes(x=obj))+
    geom_histogram(alpha=0.4,binwidth = 0.01,position = "identity")+
    #geom_density(alpha=0.1,size=1.5)+
    ## guides(fill=guide_legend(title="Heuristic:"), color = guide_legend(title="Heuristic:"))+
#    ggtitle("Objective values distribution for different heuristics (original problem, 15vars 100k dataset, non-reduced instances)")+
    geom_vline(xintercept = 1.0, size=0.5, color="red", linetype="dashed")+
    annotate("text",x=1.0, y=4.2,label = "100% = greedy BDD sifts", color="red")+
    ## styling
    labs(fill="Heuristic used:", color="Heuristic used:")+
    scale_y_continuous(
        "Density (share of instances)",
        labels = scales::number_format(accuracy = 0.5),
        position = "right"
    )+
    scale_x_continuous(
        "Objective value, relative to greedy BDD-sifts heuristic (taken as 100%)",
        labels = scales::percent,
        breaks = seq(xmin,xmax,length.out=11),
        minor_breaks = seq(xmin,xmax,length.out=21),
        limits = c(xmin,xmax)
    )+
    theme(
        legend.position = c(0.6, 0.8),
        legend.direction = "vertical",
        legend.title = element_text(size=24),
        legend.text = element_text(size=24),
        legend.text.align = 0,
        axis.text.x = element_text(size=22,angle=45,vjust = 0.7),
        axis.text.y = element_text(size=22),
        axis.title.x = element_text(size = 26),
        axis.title.y = element_text(size = 26),
        panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "lightgrey"),
        panel.grid.minor.x = element_line(size=0.5,linetype = 'solid', colour = "lightgrey"),
        panel.grid.minor.y = element_line(size=0.5,linetype = 'solid', colour = "lightgrey"),
        strip.text.y = element_text(size=22, angle=180),
        strip.background = element_blank()
    )+
    facet_grid(comment ~ ., scales="free_y", switch="y")
    ## end of styling

ggsave(opt$out,plt_dens, device = cairo_ps, width = 16, height = 10)

