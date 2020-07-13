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
df_legends = unique(select(
  filter(df, !grepl("timelog", legend)), LB, legend))

df = merge(
  x = df,
  y = select(df_legends, LB, legend),
  by.x = "LB", by.y = "LB", suffixes = c(".tech", ".capt")
)


df$caption = str_wrap(df$legend.capt, width=10)
df$caption_f = factor(df$caption, levels = unique(df$caption))

######################################################################
## draw the figure

plt_LBs =
  ggplot(filter(df, !grepl("timelog",legend.tech)), aes(x=gap))+
    geom_histogram(binwidth = 0.01,position = "identity")+
    ## geom_density(alpha=0.1,size=1.5)+
    ## guides(fill=guide_legend(title="Lower bound:"), color = guide_legend(title="Lower bound:"))+
    geom_vline(xintercept = 1.0, size=0.5, color="red", linetype="dashed")+
    annotate("text",x=1.0, y=Inf, vjust=2, label = "optimum", color="red")+
    geom_vline(xintercept = 0.0, size=0.5, color="red", linetype="dashed")+
    annotate("text",x=0.0, y=Inf, vjust=2, label = "|A|+|B|", color="red")+
    ## styling
    scale_y_continuous(
        "No. of instances",
        labels = scales::number_format(accuracy = 0.5)
    )+
    scale_x_continuous(
        "LB tightness, score",
        labels = scales::percent
        ## breaks = seq(xmin,xmax,length.out=11),
        ## minor_breaks = seq(xmin,xmax,length.out=21),
        ## limits = c(xmin,xmax)
    )+
  ## coord_cartesian(ylim=c(0,ymax))+
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
  facet_grid(caption_f ~ ., scales="free_y", switch="y")
  ## end of styling

plt_LBs_time =
  ggplot(filter(df, grepl("timelog",legend.tech)), aes(y=gap*1000,x=0))+
    ## geom_histogram(alpha=0.4,binwidth = 0.01,position = "identity")+
    ## geom_density(alpha=0.1,size=1.5)+
  geom_jitter(color='black',width=0.1, alpha=0.5)+
    geom_boxplot(notch=TRUE)+
    guides(fill=FALSE, color = FALSE)+
    ## styling
    scale_y_continuous(
      "Wall-clock time / instance, msec.",
      labels = scales::number_format(accuracy = 0.1)
      ## breaks = seq(xmin, xmax, length.out = 5),
      ## limits = c(0,xmax)
    )+
    scale_x_discrete(
        " "
        ## labels = scales::number_format(accuracy = 0.5),
        ## breaks = seq(xmin,xmax,length.out=11),
        ## minor_breaks = seq(xmin,xmax,length.out=21)
    )+
  ## coord_cartesian(xlim=c(0,xmax))+
  coord_flip()+
    theme(
        legend.position = c(0.6, 0.8),
        legend.direction = "vertical",
        legend.title = element_text(size=24),
        legend.text = element_text(size=24),
        legend.text.align = 0,
        axis.text.x = element_text(size=18,angle=90,vjust = 0.9),
        axis.text.y = element_text(size=22),
        axis.title.x = element_text(size = 26),
        axis.title.y = element_text(size = 26),
        panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "lightgrey"),
        panel.grid.minor.x = element_line(size=0.5,linetype = 'solid', colour = "lightgrey"),
      panel.grid.minor.y = element_line(size=0.5,linetype = 'solid', colour = "lightgrey"),
      strip.text.y = element_blank(), #element_text(size=22, angle=180),
      strip.background = element_blank()
    )+
  facet_grid(caption_f ~ ., scales="free_y", switch="y")


    ## end of styling

cairo_ps(opt$out, width = 16, height = 10)
grid.arrange(plt_LBs,plt_LBs_time,ncol=2)
dev.off()

