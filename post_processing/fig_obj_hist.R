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

# Internal parameters for the figure
X_QUANTILE=0.98 # quantile to filter for the histogram(Ox axis)
ORIG_BASE_COL = "orig_gsifts1p_obj" # col to divide by

# Heuristics to show ()
SHOW_HEU = c("all")
    ## c("orig_from_simpl_red_order_obj_rel_orig",
    ##          "orig_from_simpl_order_obj_rel_orig",
    ##          ## "orig_exact_gsifts_obj_rel_orig",
    ##          "orig_g2sifts_obj_rel_orig",
    ##          "orig_gsifts_obj_rel_orig",
    ##          "orig_gswaps_obj_rel_orig",
    ##          "orig_true_minAB_obj_rel_orig",
    ##          ## "orig_toA_obj_rel_orig",
    ##          ## "orig_toB_obj_rel_orig",
    ##          ## "orig_toRandom_obj_rel_orig",
    ##          "orig_from_simpl_red_order_time",
    ##          "orig_from_simpl_order_time" ,
    ##          ## "orig_exact_gsifts_time",
    ##          "orig_g2sifts_time" ,
    ##          "orig_gsifts_time" ,
    ##          "orig_gswaps_time",
    ##          "orig_true_minAB_time"
    ##          ## "orig_toA_time",
    ##          ## "orig_toB_time",
    ##          ## "orig_toRandom_time"
    ##          )

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

df_legend = filter(df, num_type=="legend")
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

# FIXME: make proper times in the solve-instance.py
## df_wide = df_wide %>%
##     mutate(
##         orig_from_simpl_order_time = simpl_bb_time + orig_from_simpl_checking_time / 8,
##         orig_from_simpl_red_order_time = orig_from_simpl_red_checking_time +
##             simpl_bb_time,
##         orig_g2sifts_time = simpl_heu_g2sifts_time + orig_from_simpl_checking_time / 8,
##         orig_gsifts_time = simpl_heu_gsifts_time + orig_from_simpl_checking_time / 8,
##         orig_gswaps_time = simpl_heu_gswaps_time + orig_from_simpl_checking_time / 8,
##         orig_toA_time = orig_from_simpl_checking_time / 8,
##         orig_toB_time = orig_from_simpl_checking_time / 8,
##         orig_toRandom_time = orig_from_simpl_checking_time / 8
##     )

## df_wide = df_wide %>%
##     mutate(
##         orig_true_minAB_obj_rel_orig = pmin(orig_toA_obj, orig_toB_obj) / orig_exact_gsifts_obj,
##         orig_true_minAB_time = orig_from_simpl_checking_time / 4
##     )

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
        heuristic = ifelse(entry_type == "obj", substring(num_type, 1, nchar(num_type)-nchar("_obj_rel")), substring(num_type, 1, nchar(num_type)-nchar("_time")))
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
    ggplot(df_time_o, aes(x=obj , fill=comment, color=comment))+
    geom_histogram(aes(y=..density..), alpha=0.4,binwidth = 0.01,position = "identity")+
    geom_density(alpha=0.1,size=1.5)+
    guides(fill=guide_legend(title="Heuristic:"), color = guide_legend(title="Heuristic:"))+
#    ggtitle("Objective values distribution for different heuristics (original problem, 15vars 100k dataset, non-reduced instances)")+
    geom_vline(xintercept = 1.0, size=0.5, color="red", linetype="dashed")+
    annotate("text",x=1.0, y=4.2,label = "100% = exact greedy sifts heuristic", color="red")+
    ## styling
    labs(fill="Heuristic used:", color="Heuristic used:")+
    scale_y_continuous(
        "Density (share of instances)",
        ##        limits = c(1,2),
        ## breaks = seq(0.5, 5.0, by=0.5),
        ## minor_breaks = seq(0,5.0,by=0.1),
        labels = scales::number_format(accuracy = 0.5)
    )+
    ## coord_cartesian(xlim = c(min(df_time_o$obj),max(df_time_o$obj)))+
    scale_x_continuous(
        "Objective value, relative to the result of exact-sifts heuristic (taken as 100%)",
        labels = scales::percent,
        breaks = seq(xmin,xmax,length.out=11),
        minor_breaks = seq(xmin,xmax,length.out=21),
        limits = c(xmin,xmax)
    )+
#    coord_cartesian(xlim=c(xmin,xmax))+
    ## scale_fill_manual(values = c(
    ##                       "orig_from_simpl_order" = "red",
    ##                       "orig_from_simpl_red_order" = "blue",
    ##                       "orig_g2sifts" = "yellow",
    ##                       "orig_gsifts" = "darkgreen",
    ##                       "orig_gswaps" = "pink",
    ##                       "orig_true_minAB" = "lightblue"
    ##                   ),labels = c(
    ##                         "orig_from_simpl_order" = "Simplified -> Branch&Bound",
    ##                         "orig_from_simpl_red_order" = "Layer-reduction -> Simplified -> Branch&Bound",
    ##                         "orig_g2sifts" = "Simplified -> Greedy sift pairs",
    ##                         "orig_gsifts" = "Simplified -> Greedy sifts",
    ##                         "orig_gswaps" = "Simplified -> Greedy swaps",
    ##                         "orig_true_minAB" = latex_label
    ##                   ))+
    ## scale_color_manual(values = c(
    ##                        "orig_from_simpl_order" = "red",
    ##                        "orig_from_simpl_red_order" = "blue",
    ##                        "orig_g2sifts" = "yellow",
    ##                        "orig_gsifts" = "darkgreen",
    ##                        "orig_gswaps" = "pink",
    ##                        "orig_true_minAB" = "lightblue"
    ##                    ),labels = c(
    ##                         "orig_from_simpl_order" = "Simplified -> Branch&Bound",
    ##                         "orig_from_simpl_red_order" = "Layer-reduction -> Simplified -> Branch&Bound",
    ##                         "orig_g2sifts" = "Simplified -> Greedy sift pairs",
    ##                         "orig_gsifts" = "Simplified -> Greedy sifts",
    ##                         "orig_gswaps" = "Simplified -> Greedy swaps",
    ##                         "orig_true_minAB" = latex_label
    ##                   ))+
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
        panel.grid.minor.y = element_line(size=0.5,linetype = 'solid', colour = "lightgrey")
          )
    ## end of styling

ggsave(opt$out,plt_dens, device = cairo_ps, width = 16, height = 10)

