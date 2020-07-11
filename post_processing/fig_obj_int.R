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

SHOW_HEU = c("all") # heuristics to show (codes)
ORIG_BASE_COL = "orig_gsifts1p_obj" # col to divide by
xNticks = 15
yNticks = 15
Npts = 50
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

df_legend = select(filter(df, num_type=="legend"), value,comment)
df = filter(df, num_type != "legend")

df$value = as.numeric(df$value)

df_wide = pivot_wider(df,
                        id_cols = "instance",
                        names_from = "num_type",
                        values_from = "value"
                        )

orig_obj_cols = grep("^orig.*_obj$",colnames(df_wide),value = TRUE)
time_cols = grep(".+time$",colnames(df_wide),value = TRUE)

for (col in orig_obj_cols){
    df_wide[[paste(col,"rel",sep="_")]] = df_wide[[col]] / df_wide[[ORIG_BASE_COL]]
}

df_rel = pivot_longer(df_wide,
                      cols = -c("instance"),
                      names_to = "num_type",
                      values_to = "value"
                      )

if (SHOW_HEU[1] == "all"){
    SHOW_HEU = grep("^orig.*_obj_rel$",colnames(df_wide),value = TRUE)
    SHOW_HEU = SHOW_HEU[ SHOW_HEU != paste(ORIG_BASE_COL,"_rel",sep="")]
}

df_rel = filter(df_rel, num_type %in% SHOW_HEU)
df_rel$entry_type = ifelse(df_rel$num_type %in% time_cols,"time","obj")

df_rel = df_rel %>%
    mutate(
        heuristic = ifelse(entry_type == "obj",
                           substring(num_type, 1, nchar(num_type)-nchar("_obj_rel")),
                           substring(num_type, 1, nchar(num_type)-nchar("_time")))
    )


df_time_o = pivot_wider(select(df_rel,-num_type),names_from = "entry_type",values_from = "value")
df_time_o = merge(x=df_time_o, y=df_legend, by.x = "heuristic", by.y = "value")

cat(paste("No. if instances not solved by greedy BDD sifts: ",as.character(nrow(filter(df_time_o, obj <0))), "\n"))

df_time_o = filter(df_time_o, obj > 0)

pfrom = min(df_time_o$obj)
pto = max(df_time_o$obj)

opt_percent = seq(from = pfrom, to = pto,length.out = Npts)

cols = unique(df_rel$num_type)

int_df = lapply(cols,
                     function(col){
                         idf = data.frame(mcol = sapply(opt_percent, function(p){
                             return( sum(df_wide[[col]] <= p) / nrow(df_wide) )
                         }));
                         colnames(idf) = substring(col, 1, nchar(col)-nchar("_obj_rel"));
                         idf;
                     })
int_df = bind_cols(int_df)
int_df$opt_percent = opt_percent
int_df_m = pivot_longer(int_df, names_to = "heuristic", values_to = "shares",cols=-c("opt_percent"))
int_df_m = merge(x=int_df_m, y=df_legend, by.x="heuristic", by.y = "value")

######################################################################
## draw the figure
apprs = int_df_m %>%
  group_by(heuristic) %>%
  group_map(~ approxfun(.x$opt_percent, .x$shares)(1.0))

names(apprs) <- unique(int_df_m$heuristic)
cat("Share of instances outperforming the baseline:")
apprs

apprs = int_df_m %>%
  group_by(heuristic) %>%
  group_map(~ approxfun(.x$opt_percent, .x$shares)(1.3))

names(apprs) <- unique(int_df_m$heuristic)
cat("Share of instances at most 30% of the baseline:")
apprs

plt_integrated =
    ggplot(int_df_m)+
    aes(x=opt_percent,y=shares,color=comment, linetype=comment, shape = comment)+
    geom_line()+
    geom_point()+
    theme_light()+
    theme(
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "lightgrey"),
        axis.text.x = element_text(angle=90,hjust=1),
        legend.position = c(0.8,0.3)
        )+
    scale_x_continuous(
        "Objective value relative to greedy BDD sifts, no worse than",
        breaks = seq(pfrom,pto,length.out = xNticks),
        limits = c(pfrom,pto),
        labels = scales::number_format(accuracy =0.01)
    )+
    scale_y_continuous(
        "Share of instances",
        breaks = seq(0,1.00,length.out = yNticks),
        labels = scales::percent_format(accuracy = 0.1)
    )+
    labs(color = "Heuristic:", linetype = "Heuristic:", shape = "Heuristic:")+
    geom_vline(xintercept = 1.0, size=0.5, color="red", linetype="dashed")+
    annotate("text",x=1.0, y=1.02,label = "100% = greedy BDD sifts", color="red")

ggsave(opt$out,plt_integrated, width = 10, height = 10)
