###
### Generates stats for random BDDs (layer sizes)
### produced by the python script
###
### (c) Alexey Bochkarev, Clemson University, 2020

library(ggplot2)
library(reshape)
library(dplyr)
library(optparse)
#library(Hmisc)

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
outfile = opt$out
Nticks = 10
quantile_to_show = 0.9

df = read.csv(file=infile, header=FALSE, sep=",",stringsAsFactors=FALSE)
colnames(df) = c("ID","N","LR",paste('n',seq(1,length(colnames(df))-4),sep="_"),"P")

df_m = melt(df, id.vars=c('ID', 'P','LR',"N"))
max_lw = quantile(df_m$value,quantile_to_show)

plt = ggplot(df_m,aes(x=value, fill=factor(P)))+
    coord_flip(xlim=c(1,max_lw))+
    xlab("Layer widths")+
    ylab("Percent of generated BDDs")+
    scale_y_continuous(labels = scales::percent_format(accuracy = 1L),
                       breaks = c(0.5,1.0))+
    scale_x_continuous(
        breaks = seq(1,max_lw,by = floor((max_lw-1)/Nticks)+1))+
        ## minor_breaks = seq(1,max_lw, length.out = 2*(floor((max_lw-1) / (floor((max_lw-1)/Nticks))))+1))+
    geom_bar(aes(y=..prop..,group=P),alpha=0.9, position='dodge')+
    facet_wrap(~variable, nrow=1)+
    theme(
        panel.background = element_rect(fill = "lightgrey",
                                        colour = "lightgrey",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "darkgrey"),
        legend.position = c(0.1,0.9),
        axis.text.x = element_text(angle=90,hjust=1)
    )+
    guides(fill=guide_legend(title="Generation parameter (p):"))+
    theme(
        legend.position = c(0.2, 0.8),
        legend.direction = "vertical",
        legend.title = element_text(size=24),
        legend.text = element_text(size=24),
        legend.text.align = 0,
        axis.text.x = element_text(size=18,angle=45,vjust = 0.7),
        axis.text.y = element_text(size=22),
        axis.title.x = element_text(size = 26),
        axis.title.y = element_text(size = 26),
        panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "lightgrey"),
        panel.grid.minor.x = element_line(size=0.25,linetype = 'solid', colour = "lightgrey"),
        panel.grid.minor.y = element_line(size=0.25,linetype = 'solid', colour = "lightgrey")
    )

ggsave(plot=plt, device = cairo_ps, outfile,width = 16, height = 9, dpi=300)
