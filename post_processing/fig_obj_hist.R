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

orig_colnames = grep("^orig.*_obj$",colnames(df_wide),value = TRUE)

for (col in orig_colnames){
    df_wide[[paste(col,"rel_orig",sep="_")]] = df_wide[[col]] / df_wide[["orig_exact_gsifts_obj"]]
}

# FIXME: make proper times in the solve-instance.py
df_wide = df_wide %>%
    mutate(
        orig_from_simpl_order_time = simpl_bb_time + orig_from_simpl_checking_time / 8,
        orig_from_simpl_red_order_time = orig_from_simpl_red_checking_time +
            simpl_bb_time,
        orig_g2sifts_time = simpl_heu_g2sifts_time + orig_from_simpl_checking_time / 8,
        orig_gsifts_time = simpl_heu_gsifts_time + orig_from_simpl_checking_time / 8,
        orig_gswaps_time = simpl_heu_gswaps_time + orig_from_simpl_checking_time / 8,
        orig_toA_time = orig_from_simpl_checking_time / 8,
        orig_toB_time = orig_from_simpl_checking_time / 8,
        orig_toRandom_time = orig_from_simpl_checking_time / 8
    )

df_wide = df_wide %>%
    mutate(
        orig_true_minAB_obj_rel_orig = pmin(orig_toA_obj, orig_toB_obj) / orig_exact_gsifts_obj,
        orig_true_minAB_time = orig_from_simpl_checking_time / 4
    )

df_rel_orig = pivot_longer(df_wide,
                      cols = -c("instance"),
                      names_to = "num_type",
                      values_to = "value"
                      )

df_rel_orig = filter(df_rel_orig, num_type %in% c("orig_from_simpl_red_order_obj_rel_orig",
                                                      "orig_from_simpl_order_obj_rel_orig",
                                                      ## "orig_exact_gsifts_obj_rel_orig",
                                                      "orig_g2sifts_obj_rel_orig",
                                                      "orig_gsifts_obj_rel_orig",
                                                      "orig_gswaps_obj_rel_orig",
                                                      "orig_true_minAB_obj_rel_orig",
                                                      ## "orig_toA_obj_rel_orig",
                                                      ## "orig_toB_obj_rel_orig",
                                                      ## "orig_toRandom_obj_rel_orig",

                                                      "orig_from_simpl_red_order_time",
                                                      "orig_from_simpl_order_time" ,
                                                      ## "orig_exact_gsifts_time",
                                                      "orig_g2sifts_time" ,
                                                      "orig_gsifts_time" ,
                                                      "orig_gswaps_time",
                                                      "orig_true_minAB_time"
                                                      ## "orig_toA_time",
                                                      ## "orig_toB_time",
                                                      ## "orig_toRandom_time"
                                        ))
df_rel_orig$entry_type = ifelse(df_rel_orig$num_type %in% grep(".+time$",colnames(df_wide),value = TRUE),"time","obj")

df_rel_orig = df_rel_orig %>%
    mutate(
        heuristic = ifelse(entry_type == "obj", substring(num_type, 1, nchar(num_type)-13), substring(num_type, 1, nchar(num_type)-5))
    )


df_time_o = pivot_wider(select(df_rel_orig,-num_type),names_from = "entry_type",values_from = "value")

######################################################################
## draw the figure
plt_dens =
    ggplot(df_time_o, aes(x=obj , fill=heuristic, color=heuristic))+
    geom_histogram(aes(y=..density..), alpha=0.4,binwidth = 0.01,position = "identity")+
    geom_density(alpha=0.1,size=1.5)+
    ## guides(fill=guide_legend(title="Heuristic code:"), color = FALSE)+
#    ggtitle("Objective values distribution for different heuristics (original problem, 15vars 100k dataset, non-reduced instances)")+
    geom_vline(xintercept = 1.0, size=0.5, color="red", linetype="dashed")+
    annotate("text",x=1.0, y=4.2,label = "100% = exact greedy sifts heuristic", color="red")+
    ## styling
    labs(fill="Heuristic used:", color="Heuristic used:")+
    scale_y_continuous(
        "Density (share of instances) for the corresponding objective values",
        ##        limits = c(1,2),
        ## breaks = seq(0.5, 5.0, by=0.5),
        ## minor_breaks = seq(0,5.0,by=0.1),
        labels = scales::number_format(accuracy = 0.5)
    )+
    ## coord_cartesian(xlim = c(min(df_time_o$obj),max(df_time_o$obj)))+
    scale_x_continuous(
        "Objective value, relative to the result of exact-sifts heuristic (taken as 100%)",
        labels = scales::percent,
        breaks = seq(min(df_time_o$obj),max(df_time_o$obj),length.out=11),
        minor_breaks = seq(min(df_time_o$obj),max(df_time_o$obj),length.out=21)
    )+
    scale_fill_manual(values = c(
                          "orig_from_simpl_order" = "red",
                          "orig_from_simpl_red_order" = "blue",
                          "orig_g2sifts" = "yellow",
                          "orig_gsifts" = "darkgreen",
                          "orig_gswaps" = "pink",
                          "orig_true_minAB" = "lightblue"
                      ),labels = c(
                            "orig_from_simpl_order" = "Simplified -> Branch&Bound ",
                            "orig_from_simpl_red_order" = "Layer-reduction -> Simplified -> Branch&Bound",
                            "orig_g2sifts" = "Simplified -> Greedy sift pairs ",
                            "orig_gsifts" = "Simplified -> Greedy sifts ",
                            "orig_gswaps" = "Simplified -> Greedy swaps ",
                            "orig_true_minAB" = parse(text = TeX("try $S_A$, $S_B$, choose the best one"))
                      ))+
    scale_color_manual(values = c(
                           "orig_from_simpl_order" = "red",
                           "orig_from_simpl_red_order" = "blue",
                           "orig_g2sifts" = "yellow",
                           "orig_gsifts" = "darkgreen",
                           "orig_gswaps" = "pink",
                           "orig_true_minAB" = "lightblue"
                       ),labels = c(
                            "orig_from_simpl_order" = "Simplified -> Branch&Bound ",
                            "orig_from_simpl_red_order" = "Layer-reduction -> Simplified -> Branch&Bound",
                            "orig_g2sifts" = "Simplified -> Greedy sift pairs ",
                            "orig_gsifts" = "Simplified -> Greedy sifts ",
                            "orig_gswaps" = "Simplified -> Greedy swaps ",
                            "orig_true_minAB" = parse(text = TeX("try $S_A$, $S_B$, choose the best one"))
                      ),
                      )+
    theme(
        legend.position = c(0.6, 0.8),
        legend.direction = "vertical",
        legend.title = element_text(size=24),
        legend.text = element_text(size=24),
        legend.text.align = 0,
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
    ## end of styling

ggsave(opt$out,plt_dens, device = cairo_ps, width = 16, height = 10)

## ## integrated figure
## pfrom = 0.50
## Nticks = 70
## opt_percent = seq(from = pfrom, to = opt_xl,length.out = Nticks)
## int_15var = data.frame(opt_perc = opt_percent,
##                        orig_share = sapply(opt_percent, function(p){
##                            return(sum(df_wide$orig_from_simpl_order_obj_rel_orig <= p) / nrow(df_wide))
##                        }),
##                        orig_red_share = sapply(opt_percent, function(p){
##                            return(sum(df_wide$orig_from_simpl_red_order_obj_rel_orig <= p) / nrow(df_wide))
##                        }),

##                        orig_g2sifts_share = sapply(opt_percent, function(p){
##                            return(sum(df_wide$orig_g2sifts_obj_rel_orig <= p) / nrow(df_wide))
##                        }),
##                        orig_gsifts_share = sapply(opt_percent, function(p){
##                            return(sum(df_wide$orig_gsifts_obj_rel_orig <= p) / nrow(df_wide))
##                        }),
##                        orig_gswaps_share = sapply(opt_percent, function(p){
##                            return(sum(df_wide$orig_gswaps_obj_rel_orig <= p) / nrow(df_wide))
##                        }),
##                        toMinAB_share =sapply(opt_percent, function(p){
##                            return(sum(df_wide$orig_true_minAB_obj_rel_orig <= p) / nrow(df_wide))
##                        })
##                        )

## int15_m = pivot_longer(int_15var, names_to = "heuristic", values_to = "shares",cols = -c("opt_perc"))
## str(int15_m)

## labels_df = data.frame(
##     heuristic = c("orig_share",
##                   "orig_red_share",
##                   "orig_g2sifts_share",
##                   "orig_gsifts_share",
##                   "orig_gswaps_share",
##                   "toMinAB_share"),
##     labels = c("Simplified problem (BB, no reduction)",
##                "Simplified problem (BB, prior reduction)",
##                "Simplified problem (greedy 2-sifts)",
##                "Simplified problem (greedy sifts)",
##                "Simplified problem (greedy swaps)",
##                "Best of {to A, to B} (original problem)"),
##     stringsAsFactors = FALSE
## )

## int15_m = merge(x = int15_m, y = labels_df, by.x = "heuristic", by.y = "heuristic")

## plt_integrated =
##     ggplot(int15_m)+
##     aes(x=opt_perc,y=shares,color=labels, linetype=labels, shape = labels)+
##     geom_line()+
##     geom_point()+
##     theme_light()+
##     theme(
##         panel.grid.major = element_line(size = 0.5, linetype = 'solid',
##                                         colour = "darkgrey"),
##         axis.text.x = element_text(angle=90,hjust=1),
##         legend.position = c(0.8,0.3)
##         )+
##     scale_x_continuous(
##         "Objective value, relative to the greedy sifts heuristic for the original problem (taken as 100%)",
##         breaks = seq(pfrom,opt_xl,length.out = Nticks),
##         limits = c(pfrom,opt_xl),
##         labels = scales::percent
##     )+
##     scale_y_continuous(
##         "Share of instances with heuristic performance no worse than the objective along the Ox axis",
##         breaks = seq(0,1.00,length.out = 50),
##         labels = scales::number_format(accuracy = 0.01)
##     )+
##     ggtitle("Heuristics quality for the original problem (integral): 100k 15-vars non-reduced instances")+
##     labs(color = "Target variable order", linetype = "Target variable order", shape = "Target variable order")

## ggsave(paste(out_dir, "fig_orig_perf_15v_integr_obj.png",sep=""),plt_integrated, width = 10, height = 10)
