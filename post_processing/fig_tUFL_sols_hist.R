######################################################################
## Creates ggplot for the original problem:
## benchmarking different heuristics
## (diagram sizes histogram for several approaches.)
##
## (c) Alexey Bochkarev, Clemson University, 2021

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

df_legend = select(filter(df, num_type=="legend"), value, comment)
df = filter(df, num_type != "legend")
df$value = as.numeric(df$value)

df_wide = pivot_wider(df,
                        id_cols = "instance",
                        names_from = "num_type",
                        values_from = "value"
                        )

rnd_cols = grep("^orig.*_rnd_obj$", colnames(df_wide), value = TRUE)
nat_cols = grep("^orig.*_nat_obj$", colnames(df_wide), value = TRUE)

for (col in rnd_cols) {
    df_wide[[paste(col,"rel",sep="_")]] = df_wide[[col]] / df_wide[["orig_gsifts1p_rnd_obj"]]
}

for (col in nat_cols) {
  df_wide[[paste(col,"rel",sep="_")]] = df_wide[[col]] / df_wide[["orig_gsifts1p_nat_obj"]]
}

df_wide$orig_gsifts1p_nat_obj_rel = NULL
df_wide$orig_gsifts1p_rnd_obj_rel = NULL

df_wide$orig_MinSimpl_nat_obj_rel = pmin(df_wide$orig_minAB_nat_obj_rel, df_wide$orig_simpl_nat_obj_rel)
df_wide$orig_MinSimpl_rnd_obj_rel <- pmin(df_wide$orig_minAB_rnd_obj_rel, df_wide$orig_simpl_rnd_obj_rel)

nat_cols = c(nat_cols, 'orig_MinSimpl_nat_obj')
rnd_cols = c(rnd_cols, "orig_MinSimpl_rnd_obj")

df_rel = pivot_longer(df_wide,
                      cols = -c("instance"),
                      names_to = "num_type",
                      values_to = "value"
                      )

df_rel = df_rel %>%
    mutate(
      entry_type = ifelse(substring(num_type, 1, nchar(num_type)-nchar("_rel")) %in% nat_cols,
                          "natural_var_order",
                          ifelse(substring(num_type, 1, nchar(num_type)-nchar("_rel")) %in% rnd_cols,
                                 "random_var_order",
                                 "aux_number")),
      heuristic = substring(num_type, 1, nchar(num_type)-nchar("_dsc_obj_rel")))

######################################################################
## draw the figure
# X_QUANTILE = 0.98
df_rel = filter(df_rel, ! heuristic %in% c('orig_gsifts1p_type2cov', 'orig_gsifts1p_cov2type'))
xmin = min(df_rel$value)
xmax = 2.5 #quantile(df_rel$value,1.0, na.rm=TRUE)

plt_dens =
    ggplot(filter(df_rel, entry_type != "aux_number"), aes(x=value))+
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
        labels = scales::number_format(accuracy = 1),
        position = "right"
    )+
    scale_x_continuous(
        "Obj value, relative to greedy BDD-sifts (taken as 100%)",
        labels = scales::percent,
        ## breaks = seq(xmin,xmax,length.out=11),
        ## minor_breaks = seq(xmin,xmax,length.out=21),
        limits = c(xmin, xmax)
    )+
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
  facet_grid(entry_type ~ heuristic, switch="y")+
  coord_flip()
    ## end of styling

ggsave(opt$out,plt_dens, width = 16, height = 10)
#ggsave("./reports/2021-04-27_Randomized_dia_sizes/rnd_order_10v1k.png",
#       width=16, height=10)  # done with cUFL_randomized_cluster.csv
## ggsave("./reports/2021-04-27_Randomized_dia_sizes/2_rnd_order_w_ctl.png",
##        width=16, height=10)

