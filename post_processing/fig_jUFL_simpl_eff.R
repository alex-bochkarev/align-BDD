######################################################################
## Creates ggplot for the original problem:
## estimates the efficiency of the simplified-problem based heuristic
## as compared to the naive `minAB` heuristic.
##
## (c) Alexey Bochkarev, Clemson University, 2021

library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(stringr)
library(latex2exp)
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

df_wide = pivot_wider(df,
                        id_cols = c("instance", "n", "prob"),
                        names_from = "num_type",
                        values_from = "value"
                        )

df_wide$simpl_rel_obj = with(df_wide, orig_simpl_obj / orig_minAB_obj)

plt_dens =
    ggplot(df_wide, aes(x=simpl_rel_obj))+
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
        "Itersection diagram size: simplified-heuristic relative to minAB",
        labels = scales::percent
        ## breaks = seq(xmin,xmax,length.out=11),
        ## minor_breaks = seq(xmin,xmax,length.out=21),
        ## limits = c(xmin, xmax)
    )+
  ggtitle(paste("Joint UFL problem: Relative efficiency of two heuristics for the original problem. From a total of", nrow(df_wide),"instances (including all combination of parameters)."))+
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
  facet_grid(n ~ prob, switch="y")
  ## coord_flip()
    ## end of styling

ggsave(opt$out,plt_dens, width = 16, height = 10)
#ggsave("./reports/2021-04-27_Randomized_dia_sizes/rnd_order_10v1k.png",
#       width=16, height=10)  # done with cUFL_randomized_cluster.csv
## ggsave("./reports/2021-04-27_Randomized_dia_sizes/2_rnd_order_w_ctl.png",
##        width=16, height=10)

