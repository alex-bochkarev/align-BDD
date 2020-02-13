######################################################################
##
## Scalability figures (from the ScalTest - generated instances)
##
## (c) Alexey Bochkarev, Clemson University, 2020
## abochka@clemson.edu

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
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
df$comment = NULL
df = filter(df,N>0)
df$value = as.numeric(df$value)

df = pivot_wider(df, id_cols = c("instance","N"), names_from = num_type, values_from = value)

df = df %>%
    mutate(
        rel_obj = orig_simpl_obj / orig_gsifts1p_obj
    )

p1 =
ggplot(df)+
    geom_point(aes(x=N, y = log(orig_simpl_time), color="Auxiliary problem / heuristic"), size=5, alpha=0.7)+
    geom_jitter(aes(x=N, y = log(orig_gsifts1p_time), color="BDD sifts"), size=5, alpha=0.5, width = 0.25)+
    scale_x_continuous(
        breaks = seq(min(df$N),max(df$N),by = 2),
        minor_breaks = seq(min(df$N),max(df$N), by=1)
    )+
    labs(x="No. of variables\n(a)",y="Solution time, sec., logarithmic scale -- ln (t)", color="Solution method:")+
    theme(
        legend.position = c(0.4, 0.8),
        legend.direction = "vertical",
        legend.title = element_text(size=24),
        legend.text = element_text(size=24),
        axis.text.x = element_text(size=22,angle=45,vjust = 0.7),
        axis.text.y = element_text(size=22),
        axis.title.x = element_text(size = 26),
        axis.title.y = element_text(size = 26),
        panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "lightgrey"),
        panel.grid.minor.x = element_line(size=0.25,linetype = 'solid', colour = "lightgrey"),
        panel.grid.minor.y = element_line(size=0.25,linetype = 'solid', colour = "lightgrey")
    )

leg = "# variables (N):"
Nmin = min(df$N)
Nmax = max(df$N)
Ns = unique(df$N)
Nmed = Ns[which.min(abs(Ns - median(Ns)))]

p2 =
ggplot(filter(df, N %in% c(Nmin,Nmed,Nmax)), aes(x=rel_obj, y=..density.., fill=as.factor(N), color=as.factor(N)))+
    geom_histogram(alpha=0.5,position="identity")+
    geom_density(alpha=0.1,position="identity")+
    guides(fill=guide_legend(title=leg), color = guide_legend(title=leg))+
    theme(
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        color = "darkgrey"),
        axis.text.x = element_text(angle=90,hjust=0.8),
        legend.position = c(0.5,0.8),
        )+
    scale_x_continuous(
        "Relative heuristic objective\n(b)",
        label = scales::percent,
        breaks = seq(min(df$rel_obj), max(df$rel_obj), length.out = 11)
    )+
    scale_y_continuous(
        "Density (no. of instances)",
        ## breaks = seq(0, nrow(df), length.out = 11),
        ## minor_breaks = seq(0, nrow(df), length.out = 21),
        label = scales::number_format(accuracy = 1.0)
    )+
    theme(
        legend.title = element_text(size=24),
        legend.text = element_text(size=24),
        axis.text.x = element_text(size=22,angle=45,vjust = 0.7),
        axis.text.y = element_text(size=22),
        axis.title.x = element_text(size = 26),
        axis.title.y = element_text(size = 26),
        panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "lightgrey"),
        panel.grid.minor.x = element_line(size=0.25,linetype = 'solid', colour = "lightgrey"),
        panel.grid.minor.y = element_line(size=0.25,linetype = 'solid', colour = "lightgrey")
    )

cairo_ps(opt$out, width = 16, height = 10)
grid.arrange(p1,p2,ncol=2)
dev.off()
