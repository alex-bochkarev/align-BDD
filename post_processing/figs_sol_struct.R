## analyzing the heuristic solutions structure
## (inspired by reviewers' notes)
## (c) A. Bochkarev, Clemson University, 2021

library(ggplot2)
library(dplyr)
library(tidyr)

df = read.csv("./run_logs/rnd_orders_nosifts.csv", stringsAsFactors = FALSE)
## df = read.csv("./run_logs/gen_test_p.csv", stringsAsFactors = FALSE)
# df = read.csv("./run_logs/gen_test_colors.csv", stringsAsFactors = FALSE)

# data preparation
legend_df = filter(df, num_type == "legend")
df = filter(df, num_type != "legend")

df$value = as.numeric(df$value)

df = df %>%
  mutate(
    n = sapply(strsplit(instance,'-'), `[`, 1)
  )

sizes = filter(df, (n == 20) & num_type %in% c('color_size_nat',
                                               'cover_size_nat',
                                               'color_size_rnd',
                                               'cover_size_rnd',
                                                 'color_size_ctl',
                                                 'cover_size_ctl'))
ggplot(data=sizes)+
  geom_histogram(aes(x=value, fill=num_type), binwidth=20)+
  facet_grid(. ~ num_type)+
  coord_flip()+
  xlab("Diagram size (no. of nodes)")+
  ylab("Instances count")+
  theme(legend.position="none")

ggsave("./reports//2021-04-27_Randomized_dia_sizes/3_ind_dia_sizes.png",
       width = 16, height = 10)
df_wide = pivot_wider(df, id_cols = c("instance","n"), names_from = num_type, values_from = value)
str(df_wide)

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
