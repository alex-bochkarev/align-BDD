######################################################################
## Creates Figure 7: runtime and objectives for different heuristics (*simplified* problem)
##
## (c) Alexey Bochkarev, Clemson University, 2021

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggridges)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(latex2exp)
  library(optparse)
})

## Internal parameters for the figure
ORIG_BASE_COL = "simpl_BB_obj" # col to divide by
Nticks = 20
Y_QUANTILE=0.98 # show "2-sigmas"
BINWIDTH=0.1
## Heuristics to show (codes)
SHOW_HEU = c("all")

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

df_legend = select(filter(df, num_type=="legend"), value,comment)
df = filter(df, num_type != "legend")
df$value = as.numeric(df$value)

df_wide = pivot_wider(df,
                      id_cols = "instance",
                      names_from = "num_type",
                      values_from = "value"
                      )

simpl_obj_cols = grep("^simpl.*_obj$",colnames(df_wide),value = TRUE)
time_cols = grep(".+time$",colnames(df_wide),value = TRUE)
######################################################################
## draw the figure

for (col in simpl_obj_cols){
    df_wide[[paste(col,"rel",sep="_")]] = df_wide[[col]] / df_wide[[ORIG_BASE_COL]]
}

df_rel = pivot_longer(df_wide,
                      cols = -c("instance"),
                      names_to = "num_type",
                      values_to = "value"
                      )

if (SHOW_HEU[1] == "all"){
    SHOW_HEU = c(grep("^simpl.*_obj_rel$",colnames(df_wide),value = TRUE),
                 grep("^simpl.*_time$",colnames(df_wide),value = TRUE))
}

df_rel = filter(df_rel, num_type %in% SHOW_HEU)

df_rel$entry_type = ifelse(df_rel$num_type %in% time_cols,"time","obj")

df_rel = df_rel %>%
    mutate(
        heuristic = ifelse(entry_type == "obj",
                           substring(num_type, 1, nchar(num_type)-nchar("_obj_rel")),
                           substring(num_type, 1, nchar(num_type)-nchar("_time")))
    )

# df_time_o = pivot_wider(select(df_rel,-num_type),names_from = "entry_type",values_from = "value")
df_time_o = merge(x=df_rel, y=df_legend, by.x = "heuristic", by.y = "value")

facet_names = unique(select(df_time_o, heuristic, comment))
labels = as.list(facet_names$comment)
names(labels) <- facet_names$heuristic

facet_names <- c(labels, list(
  "obj" = "Objective values, percent",
  "time" = "Wall-clock time / instance, msec"
  )
)

mylab = function(x){ return(facet_names[x]) }

my_ticks = annotation_logticks(sides = "b")
my_ticks$data=filter(df_time_o, comment=="Best of A and B")

df_time_o = df_time_o %>%
  mutate(
    value = ifelse(entry_type == "time", value*1000, value*100)
  )

plt_zoomed =
  ggplot(df_time_o,aes(x=value, group=comment))+
  geom_histogram(binwidth=BINWIDTH)+ scale_y_continuous(position="right")+
  guides(fill="none")+
  theme(
    axis.text.x = element_text(size=10,angle = 45,vjust = 0.7),
    axis.text.y = element_text(size=10),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 26, margin = margin(t=50)),
    panel.background = element_rect(fill = NA, color = "black"),
    panel.grid.major = element_line(size = 0.5,
                                    linetype = 'solid',
                                    color = "lightgrey"),
    strip.text.x = element_text(size = 22),
    strip.text.y.left = element_text(size = 22, angle=0),
    strip.background = element_blank()) +
  scale_x_log10()+ my_ticks+
  facet_grid(heuristic ~ entry_type, scales="free_x",
             labeller=as_labeller(mylab), switch="y")+
  ylab(paste0("Number of instances, out of ",
              length(unique(df_rel$instance))))


ggsave(opt$out,plt_zoomed,width = 16, height = 10, device = cairo_ps(family="Arial"))
