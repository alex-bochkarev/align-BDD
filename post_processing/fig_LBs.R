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
## xmin = 0.0
## xmax = max(df$gap)
######################################################################
## draw the figure
plt_LBs =
    ggplot(df, aes(x=gap , fill=legend, color=legend))+
    geom_histogram(aes(y=..density..), alpha=0.4,binwidth = 0.01,position = "identity")+
    geom_density(alpha=0.1,size=1.5)+
    guides(fill=guide_legend(title="Lower bound:"), color = guide_legend(title="Lower bound:"))+
    geom_vline(xintercept = 1.0, size=0.5, color="red", linetype="dashed")+
    annotate("text",x=1.0, y=10.0,label = "100% = exact optimum", color="red")+
    geom_vline(xintercept = 0.0, size=0.5, color="red", linetype="dashed")+
    annotate("text",x=0.0, y=10.0,label = "0% = current size, |A|+|B|", color="red")+
    ## styling
    scale_y_continuous(
        "Density (share of instances)",
        labels = scales::number_format(accuracy = 0.5)
    )+
    scale_x_continuous(
        "LB tightness, (current_size - LB) / (current_size - optimum)",
        labels = scales::percent
        ## breaks = seq(xmin,xmax,length.out=11),
        ## minor_breaks = seq(xmin,xmax,length.out=21),
        ## limits = c(xmin,xmax)
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
        panel.grid.minor.y = element_line(size=0.5,linetype = 'solid', colour = "lightgrey")
          )
    ## end of styling

ggsave(opt$out,plt_LBs, device = cairo_ps, width = 16, height = 10)

