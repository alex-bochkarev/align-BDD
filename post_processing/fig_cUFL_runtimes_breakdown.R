######################################################################
##
## t-UFLP runtimes breakdown figure.
##
## Relevant experiment: tUFLP_runtimes.py (cUFL_runtimes.py)
##
## (c) Alexey Bochkarev, Clemson University, 2021
## abochka@clemson.edu

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
  library(optparse)
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

df = read.csv(infile, stringsAsFactors = FALSE)
df$ID = paste(df$n,df$k,sep="-")
df$run = NULL
df$method = NULL
df$n = NULL  # all these are n=20

df$step <- sapply(df$step, function(x) gsub("[+-]", ".", x), USE.NAMES = FALSE)

df = pivot_wider(df, id_cols = "ID", names_from = step, values_from = time)

df = df %>%
  mutate(
    build_BDDs = BDD.build,
    CPP_MIP_pipeline = BDD.build + MIP.build + MIP.solve,
    solve_VS = VS.build.solve,
    align_BDDs = BDD.align.to.vs,
    intersect_BDDs = intersection.build,
    solve_SP = intersection.SP.solve,
    VS_pipeline = BDD.build + VS.build.solve + BDD.align.to.vs + intersection.build + intersection.SP.solve,
    naive_MIP = build.solve
  )

df = pivot_longer(select(df, ID,
                         build_BDDs,
                         solve_VS,
                         align_BDDs,
                         intersect_BDDs,
                         solve_SP,
                         VS_pipeline,
                         CPP_MIP_pipeline,
                         naive_MIP),
                         !ID, names_to = "step", values_to = "time")

steps = c("build_BDDs", "solve_VS", "align_BDDs", "intersect_BDDs", "solve_SP", "VS_pipeline", "CPP_MIP_pipeline", "naive_MIP")

captions = list("build_BDDs"="(1) Build BDDs",
                "solve_VS"="(2) Solve VS",
                "align_BDDs"="(3) Align BDDs",
                "intersect_BDDs"="(4) Build\nintersection",
                "solve_SP"="(5) Solve SP",
                "VS_pipeline"="(6) TOTAL\nCPP+SP (VS)",
                "CPP_MIP_pipeline"="(7) CPP MIP",
                "naive_MIP"="(8) Naive MIP")

df = df %>%
  mutate(
    step_factor = factor(step, levels = steps)
  )

levels(df$step_factor) <- sapply(levels(df$step_factor),
                                 function(x) return(captions[[x]]))

breakdown_plot=
ggplot(df, aes(x=time))+
  geom_histogram(position="identity", binwidth = 0.1)+
  scale_x_log10()+
  coord_flip()+
  theme(
    axis.text.x = element_text(angle=90,hjust=0.8),
    axis.text.y = element_text(size=11),
    axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26),
    panel.background = element_blank(),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgrey"),
    panel.grid.minor.x = element_line(size=0.25,linetype = 'solid', colour = "lightgrey"),
    panel.grid.minor.y = element_line(size=0.25,linetype = 'solid', colour = "lightgrey"),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 14, face="bold"))+
  facet_grid(. ~ step_factor, scales="fixed")+
  xlab("Step runtime, msec.")+
  ylab(paste("Instances count (out of ", length(unique(df$ID)), ")", sep=""))

ggsave(opt$out, breakdown_plot, width=16, height=10)
