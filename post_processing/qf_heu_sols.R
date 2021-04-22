library(ggplot2)
library(dplyr)

df = read.csv("./run_logs/heu_sol_struct_6vars.csv")
nrow(df)

ggplot(df) +
  geom_histogram(aes(x = no_opts), binwidth = 2) +
  xlab("Number of optima") +
  ylab("Count of instances") +
  ggtitle(paste("Number of different optimal solutions (out of ", nrow(df), ")", sep=""))

ggsave("./figures/experimental/no_opts_all.png")

ggplot(filter(df, no_opts <= 100)) +
  geom_histogram(aes(x = no_opts), binwidth = 2) +
  xlab("Number of optima, no more than a hundred") +
  ylab("Count of instances") +
  ggtitle(paste("Number of different optimal solutions (out of ",
                nrow(filter(df, no_opts <= 100)), ")", sep = ""))

ggsave("./figures/experimental/no_opts_reasonable.png")

ggplot(df) +
  geom_histogram(aes(x=opts_diam), binwidth=5) +
  xlab("Min similarity score between optima (i.e., 'reciprocal diameter' of the solutions set)") +
  ylab("Count of instances") +
  ggtitle(paste("Histogram of the solutions set diameter values (out of ",
          nrow(df), " instances)",
          sep = ""
      ))

ggsave("./figures/experimental/opts_diam.png")

ggplot(df) +
  geom_histogram(aes(x = best_VS_simscore), fill='grey', color='black') +
  geom_histogram(aes(x = worst_VS_simscore), fill="orange", alpha=0.2) +
  xlab("Similarity score (best and worst), 0 = fully inverted, 100 = same order") +
  ggtitle("Best and worst (across true optima) simscore for the heuristic solution.")

ggsave("./figures/experimental/heuristic_simscore.png")

ggplot(df) +
  geom_jitter(aes(x=AB_simscore, y=best_VS_simscore), width=2, height=2) +
  xlab("Simscore between A and B (original BDD orders)") +
  ylab("Simscore between the heuristic solution and a true optimum (max across optima)") +
  ggtitle("Simscore between the heuristic solution and a closest true optimum.") +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgrey"),
    panel.grid.minor.x = element_line(size=0.5,linetype = 'solid', colour = "lightgrey"),
    panel.grid.minor.y = element_line(size=0.5,linetype = 'solid', colour = "lightgrey")
  )

ggsave("./figures/experimental/heuristic_simscore_vs_AB_simscore.png")

str(df)
sum(df$opts_diam < 25)

str(df)

sum(df$best_VS_simscore > 75)
