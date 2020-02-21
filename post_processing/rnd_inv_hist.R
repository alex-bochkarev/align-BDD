library(ggplot2)
library(optparse)

#####################################################################
## unpack the command line arguments
option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL,
                    help="input filename (solved instances log)", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="./out.eps",
                    help="output file name [default= %default]",  metavar="character")
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

df = read.csv(infile)

plt = ggplot(df)+
    geom_histogram(aes(x=MIN_I_OPT))

ggsave(opt$out,plt, device = cairo_ps, width = 16, height = 10)
