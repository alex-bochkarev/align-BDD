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

df_legend = select(filter(df, num_type=="legend"), value,comment)
df = filter(df, num_type != "legend")
df$value = as.numeric(df$value)

df_wide = pivot_wider(df,
                        id_cols = "instance",
                        names_from = "num_type",
                        values_from = "value"
                        )

rnd_cols = grep("^orig.*_rnd_obj$", colnames(df_wide), value = TRUE)
nat_cols = grep("^orig.*_nat_obj$", colnames(df_wide), value = TRUE)
ctl_cols = grep("^orig.*_ctl_obj$", colnames(df_wide), value = TRUE)

for (col in rnd_cols) {
    df_wide[[paste(col,"rel",sep="_")]] = df_wide[[col]] / df_wide[["orig_gsifts1p_rnd_obj"]]
}

for (col in nat_cols) {
  df_wide[[paste(col,"rel",sep="_")]] = df_wide[[col]] / df_wide[["orig_gsifts1p_nat_obj"]]
}

for (col in ctl_cols) {
  df_wide[[paste(col,"rel",sep="_")]] = df_wide[[col]] / df_wide[["orig_gsifts1p_ctl_obj"]]
}

df_wide$orig_gsifts1p_nat_obj_rel = NULL
df_wide$orig_gsifts1p_rnd_obj_rel = NULL
df_wide$orig_gsifts1p_ctl_obj_rel = NULL

df_rel = pivot_longer(df_wide,
                      cols = -c("instance"),
                      names_to = "num_type",
                      values_to = "value"
                      )

df_rel = df_rel %>%
    mutate(
      entry_type = ifelse(substring(num_type, 1, nchar(num_type)-nchar("_rel")) %in% nat_cols,
                          'natural_var_order',
                   ifelse(substring(num_type, 1, nchar(num_type)-nchar("_rel")) %in% rnd_cols,
                          'random_var_order',
                          'control_experiment')),
      heuristic = substring(num_type, 1, nchar(num_type)-nchar("_obj_rel")))

df_rel = merge(x=df_rel, y=df_legend, by.x="heuristic", by.y="value")

df_rel = df_rel %>%
  mutate(
    heuristic_label = ifelse(entry_type == 'natural_var_order',
                             substring(comment, 1, nchar(comment)-nchar(" (natural var order)")),
                             ifelse(entry_type == 'random_var_order',
                                    substring(comment, 1, nchar(comment)-nchar(" (random var order)")),
                                    substring(comment, 1, nchar(comment)-nchar(" (control experiment)")))))

######################################################################
## draw the figure
# X_QUANTILE = 0.98
xmin = min(df_rel$value)
xmax = quantile(df_rel$value,X_QUANTILE, na.rm=TRUE)

plt_dens =
    ggplot(df_rel, aes(x=value))+
    geom_histogram(binwidth = 0.02,position = "identity")+
    #geom_density(alpha=0.1,size=1.5)+
    ## guides(fill=guide_legend(title="Heuristic:"), color = guide_legend(title="Heuristic:"))+
#    ggtitle("Objective values distribution for different heuristics (original problem, 15vars 100k dataset, non-reduced instances)")+
    geom_vline(xintercept = 1.0, size=0.5, color="red", linetype="dashed")+
    ## annotate("text",x=1.0, y=4.2,label = "100% = greedy BDD sifts", color="red")+
    ## styling
    labs(fill="Heuristic used:", color="Heuristic used:")+
    scale_y_continuous(
        "No. of instances",
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
        strip.text.y.left = element_text(size=15, angle=0),
        strip.background = element_blank(),
        strip.text.x.top = element_text(size=15, angle=0),
    )+
    facet_grid(heuristic_label ~ entry_type, switch="y")
    ## end of styling

ggsave(opt$out,plt_dens, width = 16, height = 10)


