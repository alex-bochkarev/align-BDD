######################################################################
## Creates ggplot for the original problem:
## benchmarking different heuristics
## (objective values histogram)
##
## (c) Alexey Bochkarev, Clemson University, 2020

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(latex2exp)
  library(optparse)})

## Internal parameters for the figure
BINWIDTH = 0.01

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
df_legends = unique(select(
  filter(df, !grepl("timelog", legend)), LB, legend))

df = merge(
  x = df,
  y = select(df_legends, LB, legend),
  by.x = "LB", by.y = "LB", suffixes = c(".tech", ".capt")
)


df$caption = str_wrap(df$legend.capt, width=10)
df$caption_f = factor(df$caption, levels = unique(df$caption))
df$entry_type = ifelse(grepl("timelog", df$legend.tech), "time","gap")
df = df %>% mutate(
              value = ifelse(entry_type=="gap", gap*100,
                             gap*1000)
                  )
no_instances = nrow(filter(df, entry_type=="time")) / length(unique(df$LB))
# sum(grepl("timelog", df$legend.tech), na.rm=TRUE) / nrow(df_legends)
cat(paste0("Number of instances read: ", no_instances, "\n"))

facets = list("LB_first" = "Min size\nfirst\nelement\naligned",
  "LB_last" = "Min size\nlast\nelement\naligned",
  "LB_levels" = "Inversions-\ndriven LB",
  "gap" = "LB tightness score, percent",
  "time" = "Wall-clock time per instance, msec")

plt_L =
   ggplot(filter(df, entry_type == "gap"),aes(x=value, group=caption))+
    geom_histogram(bins=50)+
    geom_vline(data=filter(df, entry_type=="gap"), aes(xintercept = 100.0), size=0.5, color="red", linetype="dashed")+
    geom_vline(data=filter(df, entry_type=="gap"), aes(xintercept = 0.0), size=0.5, color="red", linetype="dashed")+
    scale_y_continuous(position="right", limits=c(0,no_instances*0.9))+
    guides(fill="none")+
    theme(
      axis.text.x = element_text(size=18),
      axis.text.y = element_text(size=13),
      axis.title.x = element_text(size=26, margin = margin(t=25)),
      axis.title.y = element_blank(),
      panel.background = element_rect(fill = NA, color = "black"),
      panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                      color = "lightgrey"),
      panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                      color = "lightgrey"),
      strip.text.x = element_text(size = 22),
      strip.text.y.left = element_text(size = 22, angle=0),
      strip.background = element_blank()
    )+
  facet_wrap(~LB, labeller=as_labeller(function(x) return(facets[x])),
             strip.position="left", nrow=3)+
  ## ylab(paste0("Number of instances, out of ",
  ##             nrow(filter(df, entry_type=="time")) / length(unique(df$LB)),"\n"))+
  xlab("LB tightness score, percent")

plt_R =
   ggplot(filter(df, entry_type == "time"),aes(x=value, group=caption))+
    geom_histogram(bins=50)+
  scale_y_continuous(position="right", limits=c(0, no_instances*0.9))+
  scale_x_log10()+
  ## annotation_logticks(sides='b') +
    guides(fill="none")+
    theme(
      axis.text.x = element_text(size=18),
      axis.text.y = element_text(size=13),
      axis.title.x = element_text(size=26, margin = margin(t=25)),
      axis.title.y = element_text(size = 26, margin = margin(t=50)),
      panel.background = element_rect(fill = NA, color = "black"),
      panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                      color = "lightgrey"),
      panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                      color = "lightgrey"),
      strip.text.x = element_text(size = 22),
      strip.text.y.left = element_blank(),
      plot.margin = margin(l=90),
      strip.background = element_blank()
    )+
  facet_wrap(~LB, labeller=as_labeller(function(x) return(facets[x])),
             strip.position="left", nrow=3)+
  ylab(paste0("Number of instances, out of ",
              no_instances,"\n"))+
  xlab("Wall-clock time per instance, msec")+


cairo_ps(opt$out, width = 16, height = 10, family="Arial")
grid.arrange(plt_L,plt_R,ncol=2)
dev.off()

