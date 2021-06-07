###
### Generates stats for random BDDs (layer sizes)
### produced by the python script
###
### (c) Alexey Bochkarev, Clemson University, 2021

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape)
  library(dplyr)
  library(optparse)
  library(latex2exp)
  library(stringr)
})

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
Nticks = 15
quantile_to_show = 0.95

df = read.csv(file=infile, header=FALSE, sep=",",stringsAsFactors=FALSE)
## colnames(df) = c("ID","N","LR",paste('n',seq(1,length(colnames(df))-4),sep="_"),"P")
colnames(df) = c("ID","N","LR", as.character(seq(1,length(colnames(df))-4)),"P")

df_m = melt(df, id.vars=c('ID', 'P','LR',"N"))
max_lw = quantile(df_m$value,quantile_to_show)

# df_m$var_label = str_replace(as.character(df_m$variable), "_","[")
# df_m$var_label = paste(df_m$var_label,"]",sep="")
# df_m$layer = factor(df_m$var_label, levels = as.character(unique(df_m$var_label)))

plt =
  ggplot(df_m, aes(x=factor(P, levels=unique(P)), y=value, fill=P))+
    ## coord_flip(xlim=c(1,max_lw))+
    ylab("Layer widths")+
    xlab("Generation parameter value (for each layer)")+
  scale_y_continuous()+
    ## limits = quantile(df_m$value, c(0, quantile_to_show)),
    ## breaks = seq(1,max_lw,by = floor((max_lw-1)/Nticks)+1))+
  scale_x_discrete(breaks = unique(df_m$P))+
  geom_jitter(width=0.07, shape=4, alpha=0.05)+
  geom_boxplot(outlier.shape = NA, alpha=0.5, notch=TRUE)+
  facet_wrap(~variable, nrow=1)+ #, label = "label_parsed"
    theme(
        panel.background = element_rect(fill = "lightgrey",
                                        colour = "lightgrey",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "darkgrey"),
        legend.position = c(0.1,0.9),
        axis.text.x = element_text(angle=90,hjust=1)
    )+
    guides(fill="none")+
    theme(
        legend.position = c(0.2, 0.8),
        legend.direction = "vertical",
        legend.title = element_text(size=24),
        legend.text = element_text(size=24),
        legend.text.align = 0,
        axis.text.x = element_text(size=11,angle=90,vjust = 0.7),
        axis.text.y = element_text(size=22),
        axis.title.x = element_text(size = 26),
        axis.title.y = element_text(size = 26),
        panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "lightgrey"),
        panel.grid.minor.x = element_line(size=0.25,linetype = 'solid', colour = "lightgrey"),
        panel.grid.minor.y = element_line(size=0.25,linetype = 'solid', colour = "lightgrey"),
        strip.text.x = element_text(size = 22)
    )

ggsave(plot=plt, device = cairo_ps(family="Arial"), outfile, width = 16, height = 9)
