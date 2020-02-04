######################################################################
## Creates ggplots reg runtimes for different instance sizes
##
## (c) Alexey Bochkarev, Clemson University, 2019

library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(stringr)

######################################################################
## PROCESS 6-var non-reduced instances
out_dir = "../../results/100k_sets/non-reduced-dataset/"

df = read.csv("../../results/100k_sets/6vars100k_RED.full_log",
              stringsAsFactors = FALSE)

str(df)

str(df15)

df_wide = pivot_wider(df,
                      id_cols = "instance",
                      names_from = "num_type",
                      values_from = "value"
                      )

str(df_wide)

######################################################################
## runtimes per instance
## >>> 6-var instances
plt_full = ggplot(filter(df, num_type %in% grep(".+_time",colnames(df_wide),value = TRUE)))+
    geom_boxplot(aes(x=num_type,y=value))+
    geom_jitter(aes(x=num_type, y=value), alpha=0.2)+
    xlab("Operation code (instance generation and analysis)")+
    ylab("Time per operation, sec.")+
    ggtitle("Runtime of running different heuristics, per instance (wall-clock time): 100k instances, 6 variables, non-reduced")


plt_wo_bf = ggplot(filter(df, num_type %in% grep(".+_time",colnames(df_wide),value = TRUE), num_type != "orig_bf_time"))+
    geom_boxplot(aes(x=num_type,y=value))+
    geom_jitter(aes(x=num_type, y=value), alpha=0.2)+
    xlab("Operation code (instance generation and analysis)")+
    ylab("Time per operation, sec.")+
    ggtitle("Runtime of running different heuristics, per instance (wall-clock time): 100k instances, 6 variables, non-reduced (brute-force enumeration excluded)")


ggsave(paste(out_dir,"time_summary_6v100kNR_full.png",sep=""),plt_full,width = 16, height = 10)
ggsave(paste(out_dir,"time_summary_6v100kNR_wo_orig_bf.png",sep=""),plt_wo_bf,width = 16, height = 10)

## >>> 15-var instances
plt_full = ggplot(filter(df15, num_type %in% grep(".+_time",colnames(df15_wide),value = TRUE)))+
    geom_boxplot(aes(x=num_type,y=value))+
    geom_jitter(aes(x=num_type, y=value), alpha=0.2)+
    xlab("Operation code (instance generation and analysis)")+
    ylab("Time per operation, sec.")+
    ggtitle("Runtime of running different heuristics, per instance (wall-clock time): 100k instances, 15 variables, non-reduced")

ggsave(paste(out_dir,"time_summary_15v100kNR_full.png",sep=""),plt_full,width = 16, height = 10)

######################################################################
## Solution quality analysis

######################################################################
## 0) FIG 9: check if the original quality figure is the same
## >>> 6 variables case
orig_colnames = grep("^orig.*_obj$",colnames(df_wide),value = TRUE)

for (col in orig_colnames){
    df_wide[[paste(col,"rel_ex",sep="_")]] = df_wide[[col]] / df_wide[["orig_bf_obj"]]
}


df_rel_ex = pivot_longer(df_wide,
                      cols = -c("instance"),
                      names_to = "num_type",
                      values_to = "value"
                      )

df_rel_ex_f = filter(df_rel_ex, num_type %in% grep(".+_rel_ex",unique(df_rel_ex$num_type),value = TRUE))

df_rel_ex_f = filter(df_rel_ex, num_type %in% c("orig_from_simpl_order_obj_rel_ex",
                                                "orig_toA_obj_rel_ex",
                                                "orig_toB_obj_rel_ex"))

var_meds = df_rel_ex_f %>%
    group_by(num_type) %>%
    summarise(var_median = median(value)) %>%
    ungroup()

var_meds[["labels"]] = ifelse(var_meds$num_type == "orig_from_simpl_order_obj_rel_ex","Take the order from the simplified problem",
                          ifelse(var_meds$num_type == "orig_toA_obj_rel_ex","Align to A (BDD)","Align to B (BDD)"))

## str(var_means)

xl = max(df_rel_ex_f$value)*1.00

df_rel_ex_f[["labels"]] = ifelse(df_rel_ex_f$num_type == "orig_from_simpl_order_obj_rel_ex","Take the order from the simplified problem",
                          ifelse(df_rel_ex_f$num_type == "orig_toA_obj_rel_ex","Align to A (BDD)","Align to B (BDD)"))

ggplot(df_rel_ex_f, aes(x=value, fill=labels))+
    geom_density(alpha=0.5)+
    theme(
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "darkgrey"),
        axis.text.x = element_text(angle=90,hjust=1),
        legend.position = c(0.7,0.8),
    )+
    scale_x_continuous(
        "Objective value, relative to the exact minimum of the original problem (taken as 100%)",
        breaks = seq(1,xl,length.out = 50),
        limits = c(1,xl),
        labels = scales::percent
    )+
    scale_y_continuous(
        "Density of corresponding values",
        breaks = seq(0,15,length.out = 40),
        labels = scales::number_format(accuracy = 0.1)
    )+
    geom_vline(data=var_meds, aes(xintercept=var_median,color=labels),linetype="dashed",size=1)+
    guides(fill=guide_legend(title="Objectives for various heuristics:"), color=FALSE)+
    ggtitle("Heuristics quality for the original problem (density): 100k 6-vars, non-reduced instances")

ggsave(paste(out_dir,"fig9_6var_orig_performance.png",sep=""),plt_fig9_old,width = 16, height = 10)

opt_xl = 2.00
opt_percent = seq(from = 1.00, to = opt_xl,length.out = 50)

integrated = data.frame(opt_perc = opt_percent,
                        simpl_share = sapply(opt_percent, function(p){
                            return(sum(df_wide$orig_from_simpl_order_obj_rel_ex <= p) / nrow(df_wide))
                        }),
                        toA_share = sapply(opt_percent, function(p){
                            return(sum(df_wide$orig_toA_obj_rel_ex <= p) / nrow(df_wide))
                        }),
                        toB_share = sapply(opt_percent, function(p){
                            return(sum(df_wide$orig_toB_obj_rel_ex <= p) / nrow(df_wide))
                        }),
                        toMinAB_share =sapply(opt_percent, function(p){
                            return(sum(df_wide$orig_true_minAB_obj_rel_ex <= p) / nrow(df_wide))
                        })
                        )

plt = ggplot(integrated)+
    geom_line(aes(x=opt_percent,y=simpl_share,color="Simplified problem"), size=1)+
    geom_line(aes(x=opt_percent,y=toA_share,color="A order"), size=1,linetype="dashed")+
    geom_line(aes(x=opt_percent,y=toB_share,color="B order"), size=1,linetype="dotted")+
    ## geom_line(aes(x=opt_percent,y=toMinAB_share,color="best of A and B orders"), size=1,linetype="dotted")+
    theme(
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "darkgrey"),
        axis.text.x = element_text(angle=90,hjust=1),
        legend.position = c(0.7,0.8),
        )+
    scale_x_continuous(
        "Objective value, relative to the exact minimum (taken as 100%)",
        breaks = seq(1,opt_xl,length.out = 50),
        limits = c(1,opt_xl),
        labels = scales::percent
    )+
    scale_y_continuous(
        "Share of corresponding instances (probability estimate)",
        breaks = seq(0,1.00,length.out = 50),
        labels = scales::number_format(accuracy = 0.01)
    )+
    guides(color=guide_legend(title="Target order taken from:"))+
    ggtitle("Heuristics quality for the original problem (integral): 100k 6-vars, non-reduced instances")

ggsave(paste(out_dir,"fig9p_6var_orig_performance_integrated.png",sep=""),plt,width = 10, height = 10)

######################################################################
## A note on this -- based on the exact min{A,B} heuristic

df_wide[["orig_true_minAB_obj"]] = pmin(df_wide[["orig_toA_obj"]], df_wide[["orig_toB_obj"]])
df_wide[["orig_true_minAB_obj_rel_ex"]] = df_wide[["orig_true_minAB_obj"]] / df_wide[["orig_bf_obj"]]

df_rel_ex = pivot_longer(df_wide,
                      cols = -c("instance"),
                      names_to = "num_type",
                      values_to = "value"
                      )

df_rel_ex_f = filter(df_rel_ex, num_type %in% c("orig_from_simpl_order_obj_rel_ex",
                                                "orig_true_minAB_obj_rel_ex",
                                                "orig_exact_gsifts_obj_rel_ex"
                                                ))

plt_context_dens = ggplot(df_rel_ex_f, aes(x=value, fill=num_type))+
    geom_density(data=filter(df_rel_ex_f, num_type == "orig_from_simpl_order_obj_rel_ex"),
                 aes(fill = "Simplified problem"),
                 alpha=0.5)+
    geom_density(data=filter(df_rel_ex_f, num_type == "orig_true_minAB_obj_rel_ex"),
                 aes(fill = "Best of {A,B} orders"),
                 alpha=0.5)+
    geom_density(data=filter(df_rel_ex_f, num_type == "orig_exact_gsifts_obj_rel_ex"),
                 aes(fill = "Exact (BDD) sifts"),
                 alpha=0.5)+
    theme(
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "darkgrey"),
        axis.text.x = element_text(angle=90,hjust=1),
        legend.position = c(0.7,0.8),
    )+
    scale_x_continuous(
        "Objective value, relative to the exact minimum of the original problem (taken as 100%)",
        breaks = seq(1,xl,length.out = 50),
        limits = c(1,xl),
        labels = scales::percent
    )+
    scale_y_continuous(
        "Density of corresponding values",
        breaks = seq(0,100,length.out = 40),
        labels = scales::number_format(accuracy = 0.1)
    )+
    guides(fill=guide_legend(title="Objectives for various heuristics:"), color=FALSE)+
    ggtitle("Heuristics quality for the original problem (density): 100k 6-vars, non-reduced instances")

integrated2 = data.frame(opt_perc = opt_percent,
                        simpl_share = sapply(opt_percent, function(p){
                            return(sum(df_wide$orig_from_simpl_order_obj_rel_ex <= p) / nrow(df_wide))
                        }),
                        toMinAB_share =sapply(opt_percent, function(p){
                            return(sum(df_wide$orig_true_minAB_obj_rel_ex <= p) / nrow(df_wide))
                        }),
                        exact_sifts_share = sapply(opt_percent, function(p){
                            return(sum(df_wide$orig_exact_gsifts_obj_rel_ex <= p) / nrow(df_wide))
                        })
                        )

plt_integrated = ggplot(integrated2)+
    geom_line(aes(x=opt_percent,y=simpl_share,color="Simplified problem"), size=1)+
    geom_line(aes(x=opt_percent,y=toMinAB_share,color="Best of {A,B} orders"), size=1,linetype="dashed")+
    geom_line(aes(x=opt_percent,y=exact_sifts_share,color="Exact (BDD) sifts"),size=2,linetype="dotted")+
    theme_light()+
    theme(
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "darkgrey"),
        axis.text.x = element_text(angle=90,hjust=1),
        legend.position = c(0.7,0.8)
        )+
    scale_x_continuous(
        "Objective value, relative to the exact minimum (taken as 100%)",
        breaks = seq(1,opt_xl,length.out = 50),
        limits = c(1,opt_xl),
        labels = scales::percent
    )+
    scale_y_continuous(
        "Share of instances no worse than Ox value (probability estimate)",
        breaks = seq(0,1.00,length.out = 50),
        labels = scales::number_format(accuracy = 0.01)
    )+
    guides(color=guide_legend(title="Target order taken from:"))+
    ggtitle("Heuristics quality for the original problem (integral): 100k 6-vars non-reduced instances")

pdf(paste(out_dir,"NOTE_context.pdf",sep=""), width = 16, height = 10)
grid.arrange(plt_context_dens, plt_integrated,ncol=2)
dev.off()

######################################################################
## 1) simplified solution quality: performance of the SIMPLIFIED
##    heuristics' results against the SIMPLIFIED problem

## 1.1) guessing rates
## >>> 6-vars case
simpl_colnames = grep("^simpl.*_obj$",colnames(df_wide),value = TRUE)

guessing_power = data.frame(heuristic = grep(".*obj$",simpl_colnames,value = TRUE),stringsAsFactors = FALSE)
guessing_power$exact_guesses = sapply(guessing_power$heuristic, function(h){
    sum(df_wide[[h]] == df_wide[[ "simpl_obj" ]])
})

guessing_power$good_guesses = sapply(guessing_power$heuristic, function(h){
    sum(df_wide[[h]] <= df_wide[["simpl_obj"]]*1.10)
})

guessing_power = pivot_longer(guessing_power, cols = c("good_guesses","exact_guesses"), values_to = "no_insts", names_to = "guess_type")

heu_labels = data.frame(heuristic = c(
                            "simpl_BF_obj",
                            "simpl_heu_gswaps_obj",
                            "simpl_heu_gsifts_obj",
                            "simpl_heu_g2sifts_obj",
                            "simpl_toAB_obj",
                            "simpl_toR_obj"
                        ),
                        labels = c(
                            "Total instances",
                            "Greedy swaps",
                            "Greedy sifts",
                            "Greedy 2-sifts",
                            "Best of {A,B}",
                            "Random order"
                        )
                        )

guessing_power = merge(guessing_power,heu_labels, by="heuristic")

plt = ggplot(filter(guessing_power, heuristic %in% c(
                                                 "simpl_BF_obj",
                                                 "simpl_heu_gswaps_obj",
                                                 "simpl_heu_gsifts_obj",
                                                 "simpl_heu_g2sifts_obj",
                                                 "simpl_toAB_obj",
                                                 "simpl_toR_obj"
                                             )),
       aes(x=labels, y =no_insts, fill=guess_type))+
    geom_bar(stat = "identity",position="dodge")+
    geom_text(aes(label=no_insts),vjust=-0.25, position = position_dodge(width = 0.9))+
    geom_hline(yintercept = nrow(df_wide), color="red",linetype="dotted", size=0.3)+
    ggtitle("Heuristics quality for the *simplified* problem: 100k 6-vars, non-reduced instances")+
    xlab("Heuristic")+
    ylab("No. of instances when the true simplified objective was guessed pretty good (within 5%) or exactly")

ggsave(paste(out_dir,"fig12_6vars_100k.png",sep=""),plt, width = 16, height = 10)

## >>> 15 vars case
simpl_colnames15 = grep("^simpl.*_obj$",colnames(df15_wide),value = TRUE)

guessing_power15 = data.frame(heuristic = grep(".*obj$",simpl_colnames15,value = TRUE),stringsAsFactors = FALSE)
guessing_power15$exact_guesses = sapply(guessing_power15$heuristic, function(h){
    sum(df15_wide[[h]] == df15_wide[[ "simpl_obj" ]])
})

guessing_power15$good_guesses = sapply(guessing_power15$heuristic, function(h){
    sum(df15_wide[[h]] <= df15_wide[["simpl_obj"]]*1.10)
})

guessing_power15 = pivot_longer(guessing_power15, cols = c("good_guesses","exact_guesses"), values_to = "no_insts", names_to = "guess_type")

heu_labels = data.frame(heuristic = c(
                            "simpl_obj",
                            "simpl_heu_gswaps_obj",
                            "simpl_heu_gsifts_obj",
                            "simpl_heu_g2sifts_obj",
                            "simpl_toAB_obj",
                            "simpl_toR_obj"
                        ),
                        labels = c(
                            "Total instances",
                            "Greedy swaps",
                            "Greedy sifts",
                            "Greedy 2-sifts",
                            "Best of {A,B}",
                            "Random order"
                        )
                        )

guessing_power15 = merge(guessing_power15,heu_labels, by="heuristic")

plt = ggplot(filter(guessing_power15, heuristic %in% c(
                                                 "simpl_obj",
                                                 "simpl_heu_gswaps_obj",
                                                 "simpl_heu_gsifts_obj",
                                                 "simpl_heu_g2sifts_obj",
                                                 "simpl_toAB_obj",
                                                 "simpl_toR_obj"
                                             )),
       aes(x=labels, y =no_insts, fill=guess_type))+
    geom_bar(stat = "identity",position="dodge")+
    geom_text(aes(label=no_insts),vjust=-0.25, position = position_dodge(width = 0.9))+
    geom_hline(yintercept = nrow(df15_wide), color="red",linetype="dotted", size=0.3)+
    ggtitle("Heuristics quality for the *simplified* problem: 100k 15-vars, non-reduced instances")+
    xlab("Heuristic")+
    ylab("No. of instances when the true simplified objective was guessed pretty good (within 5%) or exactly")

ggsave(paste(out_dir,"fig12_15vars_100k.png",sep=""),plt, width = 16, height = 10)

######################################################################
## 1.2) time vs. objective trade-off for the SIMPLIFIED problem
## a ``2D'' figure

## >>> 6-var case
for (col in simpl_colnames){
    df_wide[[paste(col,"rel_sim",sep="_")]] = df_wide[[col]] / df_wide[["simpl_obj"]]
}

str(df_wide)

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
                                        ## "simpl_toAB_obj_rel_sim",
                                        ## "simpl_toR_obj_rel_sim"
                                        ))
unique(df_rel_sim$num_type)
## making a 2D graph

df_rel_sim$entry_type = ifelse(df_rel_sim$num_type %in% grep(".+time$",colnames(df_wide),value = TRUE),"time","obj")

df_rel_sim = df_rel_sim %>%
    mutate(
        heuristic = ifelse(entry_type == "obj", substring(num_type, 1, nchar(num_type)-12), substring(num_type, 1, nchar(num_type)-5))
    )
str(df_rel_sim)
unique(df_rel_sim$heuristic)

df_time = pivot_wider(select(df_rel_sim,-num_type),names_from = "entry_type",values_from = "value")
df_time$obj = ifelse(df_time$heuristic == "simpl_bb",1.0, df_time$obj)
str(df_time)

## full fig13-like figure
plt_full = ggplot(df_time, aes(x=time, y =obj , color=heuristic))+
    geom_point(alpha=0.4)+
    theme(
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "darkgrey"),
        axis.text.x = element_text(angle=90,hjust=1),
        legend.position = c(0.7,0.8),
        )+
    scale_y_continuous(
        "Objective value, relative to the exact minimum of the simplified problem (taken as 100%)",
        limits = c(1,2),
        labels = scales::percent
    )+
    scale_x_continuous(
        "Wall-time of calculation per instance, sec.",
        labels = scales::number_format(accuracy = 0.01)
    )+
    guides(color=guide_legend(title="Heuristics:"), color=FALSE)+
    ggtitle("Heuristics speed vs. quality, simplified problem: 100k 6-var instances, non-reduced")

ggsave(paste(out_dir, "fig13_simpl_perf_6v.png",sep=""),plt_full,width = 16, height = 10)

## "zoomed" version of fig 13
plt_zoomed = ggplot(filter(df_time, heuristic != "simpl_heu_gswaps"),
       aes(x=time, y =obj , color=heuristic))+
    geom_point(alpha=0.4)+
    theme(
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "darkgrey"),
        axis.text.x = element_text(angle=90,hjust=1),
        legend.position = c(0.7,0.8),
        )+
    scale_y_continuous(
        "Objective value, relative to the exact minimum of the simplified problem (taken as 100%)",
        limits = c(1,1.2),
        labels = scales::percent
    )+
    scale_x_continuous(
        "Wall-time of calculation per instance, sec.",
        labels = scales::number_format(accuracy = 0.01)
    )+
    guides(color=guide_legend(title="Heuristics:"), color=FALSE)+
    ggtitle("Heuristics speed vs. quality, simplified problem: 100k 6-var instances, non-reduced")

ggsave(paste(out_dir, "fig13_simpl_perf_6v_zoomed.png",sep=""),plt_full,width = 16, height = 10)

## >>> 15-var case
for (col in simpl_colnames15){
    df15_wide[[paste(col,"rel_sim",sep="_")]] = df15_wide[[col]] / df15_wide[["simpl_obj"]]
}

df15_rel_sim = pivot_longer(df15_wide,
                      cols = -c("instance"),
                      names_to = "num_type",
                      values_to = "value"
                      )

df15_rel_sim = filter(df15_rel_sim, num_type %in% c(
                                                  "simpl_heu_g2sifts_obj_rel_sim",
                                                  "simpl_heu_gsifts_obj_rel_sim",
                                                  "simpl_heu_gswaps_obj_rel_sim",
                                                  "simpl_heu_g2sifts_time",
                                                  "simpl_heu_gsifts_time",
                                                  "simpl_heu_gswaps_time",
                                                  "simpl_bb_time"
                                        ))
## making a 2D graph

df15_rel_sim$entry_type = ifelse(df15_rel_sim$num_type %in% grep(".+time$",colnames(df15_wide),value = TRUE),"time","obj")

df15_rel_sim = df15_rel_sim %>%
    mutate(
        heuristic = ifelse(entry_type == "obj", substring(num_type, 1, nchar(num_type)-12), substring(num_type, 1, nchar(num_type)-5))
    )

df_time15 = pivot_wider(select(df15_rel_sim,-num_type),names_from = "entry_type",values_from = "value")
df_time15$obj = ifelse(df_time15$heuristic == "simpl_bb",1.0, df_time15$obj)

## full fig13-like figure
plt_full = ggplot(df_time15, aes(x=time, y =obj , color=heuristic))+
    geom_point(alpha=0.4)+
    theme(
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "darkgrey"),
        axis.text.x = element_text(angle=90,hjust=1),
        legend.position = c(0.7,0.8),
        )+
    scale_y_continuous(
        "Objective value, relative to the exact minimum of the simplified problem (taken as 100%)",
        limits = c(1,2),
        labels = scales::percent
    )+
    scale_x_continuous(
        "Wall-time of calculation per instance, sec.",
        labels = scales::number_format(accuracy = 0.01)
    )+
    guides(color=guide_legend(title="Heuristics:"), color=FALSE)+
    ggtitle("Heuristics speed vs. quality, simplified problem: 100k 15-var instances, non-reduced")

ggsave(paste(out_dir, "fig13_simpl_perf_15v.png",sep=""),plt_full,width = 16, height = 10)

## "zoomed" version of fig 13
plt_zoomed = ggplot(filter(df_time15, heuristic != "simpl_heu_gswaps"),
       aes(x=time, y =obj , color=heuristic))+
    geom_point(alpha=0.4)+
    theme(
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "darkgrey"),
        axis.text.x = element_text(angle=90,hjust=1),
        legend.position = c(0.7,0.8),
        )+
    scale_y_continuous(
        "Objective value, relative to the exact minimum of the simplified problem (taken as 100%)",
        limits = c(1,1.12),
        labels = scales::percent
    )+
    scale_x_continuous(
        "Wall-time of calculation per instance, sec.",
        labels = scales::number_format(accuracy = 0.01)
    )+
    guides(color=guide_legend(title="Heuristics:"), color=FALSE)+
    ggtitle("Heuristics speed vs. quality, simplified problem: 100k 15-var instances, non-reduced (greedy swaps excluded)")

ggsave(paste(out_dir, "fig13_simpl_perf_15v_zoomed_fireplace.png",sep=""),plt_zoomed,width = 16, height = 10)

######################################################################
## simplified 2-D figure
## >>> 6-var case
heu_dict = data.frame(heu_time = c(
                          "simpl_bb",
                          "simpl_heu_g2sifts",
                          "simpl_heu_gsifts",
                          "simpl_heu_gswaps"
                      ),
                      heu_guessing = c(
                          "simpl_BF_obj",
                          "simpl_heu_g2sifts_obj",
                          "simpl_heu_gsifts_obj",
                          "simpl_heu_gswaps_obj"
                      ), stringsAsFactors = FALSE)

guess_df = merge(x=filter(guessing_power, guess_type=="exact_guesses"), y=heu_dict, by.x = "heuristic", by.y = "heu_guessing")

df_time_ex = merge(x=df_time, y = guess_df, by.x = "heuristic", by.y = "heu_time")
df_time_ex$labels = ifelse(df_time_ex$heuristic == "simpl_bb", "Branch and bound", as.character(df_time_ex$labels))

plt_2d_simpl = ggplot(filter(df_time_ex, heuristic %in% heu_dict$heu_time))+
    geom_jitter(aes(x = time, y = no_insts, color = as.character(labels)), alpha=0.3, width = 0, height = 500)+
    theme(
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "darkgrey"),
        legend.position = c(0.7,0.8),
        )+
    scale_y_continuous(
        "No. of (simplified problem) instances where an optimum was guessed exactly (noise added for readability)",
        labels = scales::number_format(accuracy = 100)
    )+
    scale_x_continuous(
        "Wall-time of calculation per instance, sec.",
        labels = scales::number_format(accuracy = 0.005)
    )+
    guides(color=guide_legend(title="Heuristics:"), color=FALSE)+
    ggtitle("Heuristics speed vs. quality measure, simplified problem: 100k 6-var instances, non-reduced")

ggsave(paste(out_dir, "fig_simpl_perf_6v_simpl.png",sep=""),plt_2d_simpl,width = 16, height = 10)

######################################################################
## >>> 15-var case

heu_dict = data.frame(heu_time = c(
                          "simpl_bb",
                          "simpl_heu_g2sifts",
                          "simpl_heu_gsifts",
                          "simpl_heu_gswaps"
                      ),
                      heu_guessing = c(
                          "simpl_obj",
                          "simpl_heu_g2sifts_obj",
                          "simpl_heu_gsifts_obj",
                          "simpl_heu_gswaps_obj"
                      ), stringsAsFactors = FALSE)

guess_df15 = merge(x=filter(guessing_power15, guess_type=="exact_guesses"), y=heu_dict, by.x = "heuristic", by.y = "heu_guessing")

df_time_ex15 = merge(x=df_time15, y = guess_df15, by.x = "heuristic", by.y = "heu_time")
df_time_ex15$labels = ifelse(df_time_ex15$heuristic == "simpl_bb", "Branch and bound", as.character(df_time_ex15$labels))

plt_2d_simpl = ggplot(filter(df_time_ex15, heuristic %in% heu_dict$heu_time))+
    geom_jitter(aes(x = time, y = no_insts, color = as.character(labels)), alpha=0.3, width = 0, height = 500)+
    theme(
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "darkgrey"),
        legend.position = c(0.7,0.8),
        )+
    scale_y_continuous(
        "No. of (simplified problem) instances where an optimum was guessed exactly (noise added for readability)",
        labels = scales::number_format(accuracy = 100)
    )+
    scale_x_continuous(
        "Wall-time of calculation per instance, sec.",
        labels = scales::number_format(accuracy = 0.005)
    )+
    guides(color=guide_legend(title="Heuristics:"), color=FALSE)+
    ggtitle("Heuristics speed vs. quality measure, simplified problem: 100k 15-var instances, non-reduced")

ggsave(paste(out_dir, "fig_simpl_perf_15v_simpl.png",sep=""),plt_2d_simpl,width = 16, height = 10)
######################################################################

######################################################################
######################################################################
## Analysis of the original problem solution quality:
## benchmarking heuristics for the original objective

## 1) fireplace-like figure

## >>> 6-var case
for (col in orig_colnames){
    df_wide[[paste(col,"rel_orig",sep="_")]] = df_wide[[col]] / df_wide[["orig_exact_gsifts_obj"]]
}

df_wide = df_wide %>%
    mutate(
        orig_from_simpl_order_time = simpl_bb_time + orig_from_simpl_checking_time / 8,
        orig_g2sifts_time = simpl_heu_g2sifts_time + orig_from_simpl_checking_time / 8,
        orig_gsifts_time = simpl_heu_gsifts_time + orig_from_simpl_checking_time / 8,
        orig_gswaps_time = simpl_heu_gswaps_time + orig_from_simpl_checking_time / 8,
        orig_toA_time = orig_from_simpl_checking_time / 8,
        orig_toB_time = orig_from_simpl_checking_time / 8,
        orig_toRandom_time = orig_from_simpl_checking_time / 8
    )

str(df_wide)

df_rel_orig = pivot_longer(df_wide,
                      cols = -c("instance"),
                      names_to = "num_type",
                      values_to = "value"
                      )

df_rel_orig = filter(df_rel_orig, num_type %in% c(
                                                    "orig_from_simpl_order_obj_rel_orig",
                                                    "orig_bf_obj_rel_orig",
                                                    "orig_exact_gsifts_obj_rel_orig",
                                                    "orig_g2sifts_obj_rel_orig",
                                                    "orig_gsifts_obj_rel_orig",
                                                    "orig_gswaps_rel_orig",
                                                    "orig_toA_obj_rel_orig",
                                                    "orig_toB_obj_rel_orig",
                                                    "orig_toRandom_obj_rel_orig",

                                                    "orig_from_simpl_order_time" ,
                                                    "orig_bf_time",
                                                    "orig_exact_gsifts_time",
                                                    "orig_g2sifts_time" ,
                                                    "orig_gsifts_time" ,
                                                    "orig_gswaps_time" ,
                                                    "orig_toA_time",
                                                    "orig_toB_time",
                                                    "orig_toRandom_time"
                                        ))
unique(df_rel_orig$num_type)
## making a 2D graph

df_rel_orig$entry_type = ifelse(df_rel_orig$num_type %in% grep(".+time$",colnames(df_wide),value = TRUE),"time","obj")

df_rel_orig = df_rel_orig %>%
    mutate(
        heuristic = ifelse(entry_type == "obj", substring(num_type, 1, nchar(num_type)-13), substring(num_type, 1, nchar(num_type)-5))
    )
str(df_rel_orig)

unique(df_rel_orig$heuristic)

df_time_o = pivot_wider(select(df_rel_orig,-num_type),names_from = "entry_type",values_from = "value")
str(df_time_o)

## full figure
plt_full =ggplot(df_time_o, aes(x=time, y=obj , color=heuristic))+
    geom_point(alpha=0.4)+
    theme(
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "darkgrey"),
        axis.text.x = element_text(angle=90,hjust=1),
        legend.position = c(0.7,0.8),
        )+
    scale_y_continuous(
        "Objective value, relative to the result of 'exact greedy sifts' heuristic (of the original problem), taken as 100%",
#        limits = c(1,2),
        labels = scales::percent
    )+
    scale_x_continuous(
        "Wall-time of calculation per instance, sec.",
        labels = scales::number_format(accuracy = 0.01)
    )+
    guides(color=guide_legend(title="Heuristics:"), color=FALSE)+
    ggtitle("Heuristics speed vs. quality, original problem: 100k 6-var instances, non-reduced")

ggsave(paste(out_dir, "fig_orig_perf_6v_full.png",sep=""),plt_full,width = 16, height = 10)

ggplot(df_time_o),
    ## filter(df_time_o,
    ##           !grepl("orig_bf", heuristic)),
    ##    aes(x=time, y=obj , color=heuristic))+
    geom_point(alpha=0.4)+
    theme(
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        color = "darkgrey"),
        axis.text.x = element_text(angle=90,hjust=1),
        legend.position = c(0.7,0.8),
        )+
    scale_y_continuous(
        "Objective value, relative to the result of 'exact greedy sifts' heuristic (of the original problem), taken as 100%",
        ## limits = c(0.8,2.5),
        labels = scales::percent
    )+
    scale_x_continuous(
        "Wall-time of calculation per instance, sec.",
        ## limits = c(0,0.05),
        labels = scales::number_format(accuracy = 0.01)
    )+
    guides(color=guide_legend(title="Heuristics:"), color=FALSE)+
    ggtitle("Heuristics speed vs. quality, original problem: 100k 6-var instances, non-reduced (brute-force enumeration excluded)")+
    facet_wrap(. ~ heuristic, ncol = 3)

