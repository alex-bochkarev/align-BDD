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
df_old = df
df = filter(df,legend != "timelog")

## xmin = 0.0
## xmax = max(df$gap)
######################################################################
## draw the figure
plt_LBs =
    ggplot(filter(df, legend != "timelog"), aes(x=gap , fill=legend, color=legend))+
    geom_histogram(aes(y=..density..), alpha=0.4,binwidth = 0.01,position = "identity")+
    geom_density(alpha=0.1,size=1.5)+
    guides(fill=guide_legend(title="Lower bound:"), color = guide_legend(title="Lower bound:"))+
    geom_vline(xintercept = 1.0, size=0.5, color="red", linetype="dashed")+
    annotate("text",x=1.0, y=Inf, vjust=2, label = "optimum", color="red")+
    geom_vline(xintercept = 0.0, size=0.5, color="red", linetype="dashed")+
    annotate("text",x=0.0, y=Inf, vjust=2, label = "|A|+|B|", color="red")+
    ## styling
    scale_y_continuous(
        "Density (share of instances)",
        labels = scales::number_format(accuracy = 0.5)
    )+
    scale_x_continuous(
        "LB tightness, score",
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

df_legends = unique(select(df, LB, legend))

df_times = merge(
    x = select(filter(df_old, legend == "timelog"), LB, gap),
    y = df_legends,
    by.x = "LB", by.y = "LB"
)
xmin = 0
xmax = quantile(df_times$gap,0.997)*1000.0
plt_LBs_time =
    ggplot(df_times, aes(x=gap*1000.0, y=..density.. , fill=legend, color=legend))+
    geom_histogram(, alpha=0.4,binwidth = 0.01,position = "identity")+
    geom_density(alpha=0.1,size=1.5)+
    guides(fill=guide_legend(title="Lower bound:"), color = guide_legend(title="Lower bound:"))+
    ## styling
    scale_y_continuous(
        "Density (share of instances)",
        labels = scales::number_format(accuracy = 0.5)
    )+
    scale_x_continuous(
        "Wall-clock time / instance, msec.",
        labels = scales::number_format(accuracy = 0.5),
        limits = c(0,xmax),
        breaks = seq(xmin,xmax,length.out=11),
        minor_breaks = seq(xmin,xmax,length.out=21),
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

cairo_ps(opt$out, width = 16, height = 10)
grid.arrange(plt_LBs,plt_LBs_time,ncol=2)
dev.off()

