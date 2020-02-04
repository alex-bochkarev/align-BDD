######################################################################
## B&B performance figures
##
## (c) Alexey Bochkarev, 2019
## abochka@clemson.edu

library(ggplot2)
library(dplyr)
library(tidyr)

out_dir = "../../results/100k_sets/BB_performance/"
input_file = "../../results/100k_sets/BB_performance/15vars100kNR_BB.full_log"

df = read.csv(input_file,
              stringsAsFactors = FALSE)

str(df)
unique(df$num_type)
unique(df$comment)

steplog = filter(df, num_type == "steplog")

## extract optima
opts = filter(df, comment == "status:optimal")
rownames(opts) = opts$instance
str(opts)

## adjust UB and LB (to make them comparable)
steplog = steplog %>%
    mutate(
        LB_rel = LB / opts[as.character(instance),"LB"],
        UB_rel = UB / opts[as.character(instance),"UB"]
    )

str(steplog)
summary(steplog$UB_rel)
summary(steplog$LB_rel)

steps = c(40,80,100,140,200,240,300,400)

steps_df = steplog %>%
    filter(grepl("--none--",comment)) %>%
    select(instance,step,LB_rel,UB_rel) %>%
    filter(as.character(step) %in% steps) ## omit 'optimal' timestamps

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

str(full_steps_df)

full_steps_df$step = as.numeric(full_steps_df$step)

full_steps_df = full_steps_df[
    order(full_steps_df$step, full_steps_df$instance),
]

full_steps_df$bound_label = ifelse(full_steps_df$bound_type == "LB_rel", "Lower bound", "Upper bound")

plt_BB =
    ggplot(full_steps_df)+
    geom_histogram(
        aes(x=rel_value,y=..count..,fill=bound_label),
        position="identity",
        binwidth = 0.01, alpha=0.5)+
    facet_grid(. ~ step)+
    geom_vline(xintercept = 1.0,color="red",size=0.5,linetype="dashed")+
    ggtitle("Upper and lower bounds distributions at different BB-search steps; one facet = a snapshot for one step (search tree node expansion) number, 100k instances, 15 vars, non-reduced")+
    theme(
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "darkgrey"),
        axis.text.x = element_text(angle=90,hjust=1),
        legend.position = c(0.9,0.9),
        )+
    scale_x_continuous(
        "Bound values, relative to the respective optimal value",
        breaks = seq(0.75,1.15,length.out = 25),
        labels = scales::percent
    )+
    scale_y_continuous(
        "Instances count for resp. bounds value (out of 100,000 instances)",
        labels = scales::number_format(accuracy = 100),
        breaks = seq(0,100000,length.out = 9)
    )+
    coord_flip(xlim=c(0.75, 1.15))

plt_BB
ggsave(paste(out_dir,"BB_hists.png",sep=""),plt_BB, width = 16, height = 10)

steptime_df = filter(df, grepl("timelog",num_type))
str(steptime_df)

sps = data.frame(
    instance = steptime_df$instance,
    time_per_step = (steptime_df$LB / steptime_df$step) * 1000 # in milliseconds
)

ggplot(data=sps)+
    geom_density(aes(x=time_per_step))+
    theme(
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "darkgrey"),
        axis.text.x = element_text(angle=90,hjust=1),
        )+
    scale_x_continuous(
        "Mean time per one search tree node expansion step (one entry = mean time for an instance), ms (0.001 sec)",
        breaks = seq(0,max(sps$time_per_step),length.out = 100),
        labels = scales::number_format(accuracy = 0.01)
    )+
    ylab("Density for the corresponding value")+
    ggtitle("Mean time per search tree node expansion step, 100k VarSeq instances")

step_vals = c(40,80,100,140,200,240,300,400)

df = pivot_wider(select(full_steps_df, -bound_label), names_from = c("bound_type"), values_from = "rel_value")
str(df)

N_inst = length(unique(df$instance))

gaps_df = data.frame(step = step_vals,
                     instshares_0.5 = sapply(step_vals, function(s){
                         with(filter(df, step == s),
                              sum(UB_rel - LB_rel <= 0.5) / N_inst
                     )}
                     ),
                     ## instshares_0.25 = sapply(step_vals, function(s){
                     ##     with(filter(df, step == s),
                     ##          sum(UB_rel - LB_rel <= 0.25) / N_inst
                     ##          )}
                     ##     ),
                     ## instshares_0.15 = sapply(step_vals, function(s){
                     ##     with(filter(df, step == s),
                     ##          sum(UB_rel - LB_rel <= 0.15) / N_inst
                     ##          )}
                     ##     ),
                     instshares_0.10 = sapply(step_vals, function(s){
                         with(filter(df, step == s),
                              sum(UB_rel - LB_rel <= 0.10) / N_inst
                              )}
                         ),
                     instshares_0.05 = sapply(step_vals, function(s){
                         with(filter(df, step == s),
                              sum(UB_rel - LB_rel <= 0.05) / N_inst
                              )}
                         ),
                     instshares_0.0 = sapply(step_vals, function(s){
                         with(filter(df, step == s),
                              sum(UB_rel - LB_rel == 0.0) / N_inst
                              )}
                         )
                     )

gaps_df = pivot_longer(gaps_df,
                       cols = -step,
                       names_to = c("parameter","gap_size"),
                       names_sep = "_"
                       )
gaps_df$parameter = NULL
str(gaps_df)


plt_convergence =
    ggplot(gaps_df)+
    geom_line(aes(x = step, y = value, color=gap_size, linetype = gap_size), size=0.7)+
    theme(
        legend.position = c(0.8, 0.5),
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
        )+
    scale_x_continuous(
        "Branch-and-bound step number (search tree node expansion)",
        breaks = seq(0,400,length.out = 21),
        minor_breaks = seq(0,400,length.out = 41),
        labels = scales::number_format(accuracy = 1)
    )+
    scale_y_continuous(
        "Share of instances",
        breaks = seq(0,1,length.out = 11),
        minor_breaks = seq(0,1,length.out = 51),
        labels = scales::percent
    )+
    scale_linetype_manual(values = c(
                           "0.0" = "solid",
                           "0.05" = "longdash",
                           "0.10" = "dotted",
                           "0.5" = "twodash"
                       ),
                       labels = c(
                           "0.0" = "0% (exact match)",
                           "0.05" = "5% of optimum",
                           "0.10" = "10% of optimum",
                           "0.5" = "50% of optimum"
                       ),
                       )+
    scale_color_manual(values = c(
                              "0.0" = "blue",
                              "0.05" = "red",
                              "0.10" = "darkgreen",
                              "0.5" = "black"
                          ),
                          labels = c(
                              "0.0" = "0% (exact match)",
                              "0.05" = "5% of optimum",
                              "0.10" = "10% of optimum",
                              "0.5" = "50% of optimum"
                          ), guide="UB/LB gap size:"
                          )+
#    ggtitle("Convergence of the Branch and Bound procedure for the simplified problem")+
    guides(color=guide_legend("UB/LB gap size:"), linetype=guide_legend("UB/LB gap size:"))+
    geom_hline(yintercept = 1.0, color="red",linetype="dotted",size=1)


ggsave(paste(out_dir,"BB_conv.eps",sep=""),plt_convergence, device=cairo_ps,width = 16, height = 10)
