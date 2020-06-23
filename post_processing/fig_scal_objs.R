######################################################################
##
## Scalability figures (from the ScalTest - generated instances)
## draws objective distributions for different problem sizes
##
## (c) Alexey Bochkarev, Clemson University, 2020
## abochka@clemson.edu

library(ggplot2)
library(tidyr)
library(dplyr)
library(optparse)

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)

    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)

    numPlots = length(plots)

    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    } 

    if (numPlots==1) {
         print(plots[[1]])
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout =
                    grid.layout(nrow(layout),
                    ncol(layout))))

        # Make each plot, in the
        # correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}

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

Ns = sort(unique(df$N))

myplots = lapply(Ns, function(Nval) {
        mydf = filter(df, N==Nval)
        ggplot(mydf, aes(x=rel_obj, y=..density..))+
            geom_histogram(alpha=0.5,position="identity")+
            geom_density(alpha=0.1,position="identity")+
            theme(
            panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        color = "darkgrey"),
            axis.text.x = element_text(angle=90,hjust=0.8)
            )+
            scale_x_continuous(
                "Relative heuristic objective",
                label = scales::percent,
                limits = c(min(df$rel_obj), max(df$rel_obj)),
                breaks = seq(min(df$rel_obj), max(df$rel_obj), length.out = 11)
            )+
            scale_y_continuous(
                "Density (no. of instances)",
                ## breaks = seq(0, nrow(df), length.out = 11),
                ## minor_breaks = seq(0, nrow(df), length.out = 21),
                label = scales::number_format(accuracy = 1.0)
            )+
            ggtitle(paste("N=",Nval,sep=""))
    })

png(opt$out, width=16, height=10, units="in", res=300)
multiplot(plotlist=myplots, cols=5)
dev.off()
