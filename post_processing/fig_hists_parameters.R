## analyzing the heuristic solutions structure
## (inspired by reviewers' notes)
## (c) A. Bochkarev, Clemson University, 2021

library(ggplot2)
library(dplyr)
library(tidyr)

df = read.csv("./run_logs/tUFL_rnd_DDs_30vars.csv", stringsAsFactors = FALSE)

# data preparation
legend_df = filter(df, num_type == "legend")
df = filter(df, num_type != "legend")

df$value = as.numeric(df$value)

df_wide = pivot_wider(df, id_cols = "instance", names_from = num_type, values_from = value)
str(df_wide)

df_wide = df_wide %>%
  mutate(
    simpl_nat_rel = orig_simpl_nat_obj / orig_minAB_nat_obj,
    simpl_rnd_rel = orig_simpl_rnd_obj / orig_minAB_rnd_obj
  )

ggplot(df_wide)+
  geom_point(aes(x=log(orig_simpl_rnd_obj), y= simpl_rnd_rel))

dfl = pivot_longer(select(df_wide, simpl_nat_rel, simpl_rnd_rel, instance), cols=c("simpl_nat_rel", "simpl_rnd_rel"), names_to="heuristic", values_to = "value")

ggplot(data=dfl)+
  geom_histogram(aes(x=value), binwidth=0.01)+
  facet_grid(heuristic ~ .)

str(dfl)
str(df_wide)

sum(df_wide$simpl_nat_rel < 1.0) / nrow(df_wide)
sum(df_wide$simpl_rnd_rel < 1.0) / nrow(df_wide)

## ggsave("./reports//2021-04-27_Randomized_dia_sizes/3_ind_dia_sizes.png",
       ## width = 16, height = 10)


## second experiment: varying p
dfp = read.csv("./run_logs/tUFL_varying_p.csv", stringsAsFactors = FALSE)

legend_dfp = filter(dfp, num_type == "legend")
dfp = filter(dfp, num_type != "legend")

dfp$value = as.numeric(dfp$value)

df_wide = pivot_wider(dfp, id_cols = "instance", names_from = num_type, values_from = value)
df_wide = df_wide %>%
  mutate(
    simpl_nat_rel = orig_simpl_nat_obj / orig_minAB_nat_obj,
    simpl_rnd_rel = orig_simpl_rnd_obj / orig_minAB_rnd_obj
  )

dfl = pivot_longer(select(df_wide, simpl_nat_rel, simpl_rnd_rel, instance), cols=c("simpl_nat_rel", "simpl_rnd_rel"), names_to="heuristic", values_to = "value")

dfl = dfl %>%
  mutate(
    p = sapply(strsplit(instance,'-'), `[`, 1)
  )

ggplot(data=filter(dfl, heuristic == "simpl_nat_rel"))+
  geom_histogram(aes(x=value), binwidth=0.01)+
  facet_grid(. ~ p)+
  coord_flip()+
  xlim(0.8,1.5)

ggplot(data=filter(dfl, heuristic == "simpl_rnd_rel"))+
  geom_histogram(aes(x=value), binwidth=0.01)+
  facet_grid(. ~ p)+
  coord_flip()+
  xlim(0.6,1.8)

### old code (to be deleted)
str(dfl)
df_wide = df_wide %>%
  mutate(
    simpl2AB_nat = orig_simpl_nat_obj / orig_minAB_nat_obj
  )

## simpl_to_AB = filter(df_wide, n %in% c(5,10, 12, 14, 18, 20))

#ggplot(data = simpl_to_AB)+

ggplot(data = filter(df_wide, n != 5))+
# ggplot(data=df_wide)+
  geom_histogram(aes(x=simpl2AB_nat), binwidth=0.02)+
  coord_flip()+
  facet_grid(. ~ n)+
  xlab("Heuristic objective / minAB objective (both wrt original problem)")+
  ylab("Instances count")

## ggsave("./reports//2021-04-27_Randomized_dia_sizes/4_gen_parameter_p.png",
# ggsave("./reports//2021-04-27_Randomized_dia_sizes/5_gen_parameter_colors.png",
ggsave("./reports//2021-04-27_Randomized_dia_sizes/6_dia_sizes_vs_problem_sizes.png",
       width = 16, height=10)
