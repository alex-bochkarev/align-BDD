######################################################################
## Creates ggplot for the simplified problem:
## Calculation wall-clock time vs. objective for different
## heuristics
##
## (c) Alexey Bochkarev, Clemson University, 2020

library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(stringr)
library(latex2exp)
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

df = read.csv(opt$input, stringsAsFactors = FALSE)

df_wide = pivot_wider(df,
                      id_cols = "instance",
                      names_from = "num_type",
                      values_from = "value"
                      )

######################################################################
## draw the figure
simpl_colnames = grep("^simpl.*_obj$",colnames(df_wide),value = TRUE)

for (col in simpl_colnames){
    df_wide[[paste(col,"rel_sim",sep="_")]] = df_wide[[col]] / df_wide[["simpl_obj"]]
}

df_rel_sim = pivot_longer(df_wide,
                      cols = -c("instance"),
                      names_to = "num_type",
                      values_to = "value"
                      )

df_rel_sim = filter(df_rel_sim, num_type %in% c(
                                                  "simpl_heu_g2sifts_obj_rel_sim",
                                                  "simpl_heu_gsifts_obj_rel_sim",
                                                  "simpl_heu_gswaps_obj_rel_sim",
                                                  "simpl_heu_g2sifts_time",
                                                  "simpl_heu_gsifts_time",
                                                  "simpl_heu_gswaps_time",
                                                  "simpl_bb_time"
                                        ))

df_rel_sim$entry_type = ifelse(df_rel_sim$num_type %in% grep(".+time$",colnames(df_wide),value = TRUE),"time","obj")

df_rel_sim = df_rel_sim %>%
    mutate(
        heuristic = ifelse(entry_type == "obj", substring(num_type, 1, nchar(num_type)-12), substring(num_type, 1, nchar(num_type)-5))
    )

df_time = pivot_wider(select(df_rel_sim,-num_type),names_from = "entry_type",values_from = "value")
df_time$obj = ifelse(df_time$heuristic == "simpl_bb",1.0, df_time$obj)

plt_zoomed =
    ggplot(filter(df_time, heuristic != "simpl_heu_gswaps"),
       aes(x=time, y =obj , shape=heuristic,color=heuristic))+
    geom_point(alpha=0.4, size=4)+
    scale_y_continuous(
        "Objective value (percent of the exact min)",
        limits = c(1,1.10),
        labels = scales::percent
    )+
    scale_x_continuous(
        "Calculation (wall-clock) time per instance, sec.",
        labels = scales::number_format(accuracy = 0.1),
        breaks = seq(0,max(df_time$time),length.out = 20)
    )+
    guides(
        color=guide_legend(title="Heuristic / method:"),
        shape=guide_legend(title="Heuristic / method:")
        )+
    scale_color_manual(values = c(
                          "simpl_bb" = "red",
                          "simpl_heu_gsifts" = "blue",
                          "simpl_heu_g2sifts" = "green"
                      ),
                      labels = c(
                          "simpl_bb" = "branch-and-bound",
                          "simpl_heu_gsifts" = "greedy sifts",
                          "simpl_heu_g2sifts" = "greedy sift pairs"
                      ),
                      )+
    scale_shape_manual(values = c(
                           "simpl_bb" = "triangle",
                           "simpl_heu_gsifts" = "cross",
                           "simpl_heu_g2sifts" = "circle"
                       ),
                       labels = c(
                           "simpl_bb" = "branch-and-bound",
                           "simpl_heu_gsifts" = "greedy sifts",
                           "simpl_heu_g2sifts" = "greedy sift pairs"
                       ),
                       )+
    theme(
        legend.position = c(0.8, 0.8),
        legend.direction = "vertical",
        legend.title = element_text(size=24),
        legend.text = element_text(size=24),
        axis.text.x = element_text(size=22,angle=45,vjust = 0.7),
        axis.text.y = element_text(size=22),
        axis.title.x = element_text(size = 26),
        axis.title.y = element_text(size = 26, margin = margin(t=50)),
        panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        color = "darkgrey"),
          )

ggsave(opt$out,plt_zoomed,width = 16, height = 10, device = cairo_ps)
