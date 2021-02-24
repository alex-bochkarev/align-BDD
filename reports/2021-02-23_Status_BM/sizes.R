library(ggplot2)

df = read.csv("./run_logs/dia_sizes.csv")

ggplot(df) +
    geom_jitter(aes(x = n, y = A_size, color = "|A|"), width = 0.1) +
    geom_jitter(aes(x = n, y = C_size, color = "|C|"), width = 0.1) +
    geom_jitter(aes(x = n, y = plain_MIP_vars, color = "no. of variables in plain MIP"), width = 0.1) +
    theme(
        legend.position = c(0.2, 0.8),
        axis.text.x = element_text(vjust = 0.5, hjust = 1, size = 18),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 26),
        axis.title.y = element_text(size = 26, margin = margin(t = 50)),
        panel.background = element_rect(fill = NA, color = "black"),
        panel.grid.major = element_line(
            size = 0.5, linetype = "solid",
            color = "lightgrey"
        )
    )+
  ylab("|A|, |C|, and # of vars")

ggsave("./reports/2021-02-23_Status_BM/AandCsizes.eps")

ggplot(df) +
  geom_jitter(aes(x = n, y = log(int_size), color = "|A∩C| (aux problem)"), width = 0.1) +
  geom_jitter(aes(x = n, y = log(int_size_A2C), color = "|A∩C| (A-to-C)"), width = 0.1)+
  theme(
    legend.position = c(0.2, 0.8),
    axis.text.x = element_text(vjust = 0.5, hjust = 1, size = 18),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26, margin = margin(t = 50)),
    panel.background = element_rect(fill = NA, color = "black"),
    panel.grid.major = element_line(
      size = 0.5, linetype = "solid",
      color = "lightgrey"
    )
  ) +
  ylab("*LOG* # of nodes (!)")

ggsave("./reports/2021-02-23_Status_BM/inter_size.eps")

ggplot(df) +
  geom_jitter(aes(x = n, y = exp_time), width = 0.1)+
  theme(
      legend.position = c(0.2, 0.8),
      axis.text.x = element_text(vjust = 0.5, hjust = 1, size = 18),
      axis.text.y = element_text(size = 13),
      axis.title.x = element_text(size = 26),
      axis.title.y = element_text(size = 26, margin = margin(t = 50)),
      panel.background = element_rect(fill = NA, color = "black"),
      panel.grid.major = element_line(
          size = 0.5, linetype = "solid",
          color = "lightgrey"
      )
  ) +
  ylab("Experiment time, sec")

ggsave("./reports/2021-02-23_Status_BM/runtimes.eps")
