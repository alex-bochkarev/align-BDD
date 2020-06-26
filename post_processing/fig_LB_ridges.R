######################################################################
## Creates ggplot for the simplified problem:
## B&B convergence figure
##
## (c) Alexey Bochkarev, Clemson University, 2020
## abochka@clemson.edu

library(ggplot2)
library(dplyr)
library(tidyr)
library(optparse)
library(ggridges)

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
df = read.csv(infile,
              stringsAsFactors = FALSE)

steplog = filter(df, num_type == "steplog")

## extract optima
opts = filter(df, comment == "status:optimal")
rownames(opts) = opts$instance

## adjust UB and LB (to make them comparable)
steplog = steplog %>%
    mutate(
        LB_rel = LB / opts[as.character(instance),"LB"],
        UB_rel = UB / opts[as.character(instance),"UB"]
    )

## s_min = 10
## s_max = ceiling(max(steplog$step / 100))*100
## step_vals = seq(s_min,s_max, length.out = 21)
#    c(40,80,100,140,200,240,300,400)
step_vals = unique(filter(steplog,grepl("--none--",comment))$step)

steps_df = steplog %>%
    filter(grepl("--none--",comment)) %>%
    select(instance,step,LB_rel,UB_rel)

    ## filter(step %in% step_vals) ## omit 'optimal' timestamps

full_steps_df = steps_df %>%
    pivot_wider(id_cols = "instance", values_from = c("LB_rel","UB_rel"),
                names_from = "step",
                names_sep = ":",
                values_fill = list(LB_rel=1.0, UB_rel=1.0)
                )

full_steps_df = full_steps_df %>%
    pivot_longer(
        cols = -c("instance"),
        names_to = c("bound_type", "step"),
        names_sep = ":",
        values_to = "rel_value",
        )

full_steps_df$step = as.numeric(full_steps_df$step)

######################################################################
## draw the figure

## df = pivot_wider(full_steps_df, names_from = c("bound_type"), values_from = "rel_value")

## N_inst = length(unique(df$instance))

## gaps_df = data.frame(step = step_vals,
##                      instshares_0.5 = sapply(step_vals, function(s){
##                          with(filter(df, step == s),
##                               sum(UB_rel - LB_rel <= 0.5) / N_inst
##                      )}
##                      ),
##                      instshares_0.10 = sapply(step_vals, function(s){
##                          with(filter(df, step == s),
##                               sum(UB_rel - LB_rel <= 0.10) / N_inst
##                               )}
##                          ),
##                      instshares_0.05 = sapply(step_vals, function(s){
##                          with(filter(df, step == s),
##                               sum(UB_rel - LB_rel <= 0.05) / N_inst
##                               )}
##                          ),
##                      instshares_0.0 = sapply(step_vals, function(s){
##                          with(filter(df, step == s),
##                               sum(UB_rel - LB_rel == 0.0) / N_inst
##                               )}
##                          )
##                      )

## gaps_df = pivot_longer(gaps_df,
##                        cols = -step,
##                        names_to = c("parameter","gap_size"),
##                        names_sep = "_"
##                        )
## gaps_df$parameter = NULL

df = filter(full_steps_df, (bound_type %in% c("UB_rel")))

plt_convergence =
  ggplot(df, aes(x=as.factor(step), y=rel_value))+
  geom_violin()
  ## geom_density_ridges(scale=2)+
  ## coord_cartesian(xlim=c(0.9,1.2))
    ## theme(
    ##     legend.position = c(0.8, 0.5),
    ##     legend.direction = "vertical",
    ##     legend.title = element_text(size=24),
    ##     legend.text = element_text(size=24),
    ##     axis.text.x = element_text(size=22,angle=45,vjust = 0.7),
    ##     axis.text.y = element_text(size=22),
    ##     axis.title.x = element_text(size = 26),
    ##     axis.title.y = element_text(size = 26),
    ##     panel.background = element_blank(),
    ##     panel.grid.major = element_line(size = 0.5, linetype = 'solid',
    ##                                     colour = "lightgrey"),
    ##     panel.grid.minor.x = element_line(size=0.25,linetype = 'solid', colour = "lightgrey"),
    ##     panel.grid.minor.y = element_line(size=0.25,linetype = 'solid', colour = "lightgrey")
    ##     )
##     scale_x_continuous(
##         "Branch-and-bound step number (search tree node expansion)",
##         breaks = seq(0,max(step_vals),length.out = 21),
##         minor_breaks = seq(0,max(step_vals),length.out = 41),
##         labels = scales::number_format(accuracy = 1)
##     )+
##     scale_y_continuous(
##         "Share of instances",
##         breaks = seq(0,1,length.out = 11),
##         minor_breaks = seq(0,1,length.out = 51),
##         labels = scales::percent
##     )+
##     scale_linetype_manual(values = c(
##                            "0.0" = "solid",
##                            "0.05" = "longdash",
##                            "0.10" = "dotted",
##                            "0.5" = "twodash"
##                        ),
##                        labels = c(
##                            "0.0" = "0% (exact match)",
##                            "0.05" = "5% of optimum",
##                            "0.10" = "10% of optimum",
##                            "0.5" = "50% of optimum"
##                        ),
##                        )+
##     scale_color_manual(values = c(
##                               "0.0" = "blue",
##                               "0.05" = "red",
##                               "0.10" = "darkgreen",
##                               "0.5" = "black"
##                           ),
##                           labels = c(
##                               "0.0" = "0% (exact match)",
##                               "0.05" = "5% of optimum",
##                               "0.10" = "10% of optimum",
##                               "0.5" = "50% of optimum"
##                           ), guide="UB/LB gap size:"
##                           )+
## #    ggtitle("Convergence of the Branch and Bound procedure for the simplified problem")+
##     guides(color=guide_legend("UB/LB gap size:"), linetype=guide_legend("UB/LB gap size:"))+
##     geom_hline(yintercept = 1.0, color="red",linetype="dotted",size=1)


ggsave(opt$out,plt_convergence, device=cairo_ps,width = 16, height = 10)
