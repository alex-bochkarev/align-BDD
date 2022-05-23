library(ggplot2)

df = read.csv("./run_logs/dcloud_jUFLP_cavemen.csv", stringsAsFactors = FALSE)

str(df)


df$t_VS_rel = df$t_VS / df$t_toA
df$size_VS_rel = df$intsize_VS / df$intsize_toA


ggplot(data=df)+
  geom_histogram(aes(x=t_VS_rel, fill="runtime, sec"), fill='black', alpha=0.3)+
  geom_histogram(aes(x=size_VS_rel, fill="int. size, sec"), fill='red', alpha=0.3)+
  geom_vline(xintercept=1.0, color='red')+
  theme(
    axis.text.x = element_text(vjust = 0.5, hjust = 1, size=18),
    axis.text.y = element_text(size=13),
    axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26, margin = margin(t=50)),
    panel.background = element_rect(fill = NA, color = "black"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    color = "lightgrey"),
    strip.text.x = element_text(size = 22),
    strip.text.y = element_text(size=22),
    strip.background = element_blank()
  )+
  xlim(0,5.0)+
  xlab("Runtime and intersection dia size, relative to 'to-A' heuristic.")+
  ylab(paste("count ( out of", nrow(df), ")"))

ggsave("./jUFLP_hist.png", width = 16, height = 10)
sum(df$t_VS_rel < 1.0) / nrow(df)

ggplot(data=df)+
  geom_point(aes(y=t_VS, x=t_toA))+
  geom_abline(slope = 1.0, color='red')+
  theme(
    axis.text.x = element_text(vjust = 0.5, hjust = 1, size=18),
    axis.text.y = element_text(size=13),
    axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26, margin = margin(t=50)),
    panel.background = element_rect(fill = NA, color = "black"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    color = "lightgrey"),
    strip.text.x = element_text(size = 22),
    strip.text.y = element_text(size=22),
    strip.background = element_blank()
  )+
  xlab("Runtime: 'to A' heuristic")+
  ylab("Runtime: 'VS' heuristic")

ggsave("./jUFLP_runtimes.png", width = 16, height = 10)

ggplot(data=df)+
  geom_point(aes(y=intsize_VS, x=intsize_toA))+
  geom_abline(slope = 1.0, color='red')+
  theme(
    axis.text.x = element_text(vjust = 0.5, hjust = 1, size=18),
    axis.text.y = element_text(size=13),
    axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26, margin = margin(t=50)),
    panel.background = element_rect(fill = NA, color = "black"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    color = "lightgrey"),
    strip.text.x = element_text(size = 22),
    strip.text.y = element_text(size=22),
    strip.background = element_blank()
  )+
  xlab("Int. size: 'to A' heuristic")+
  ylab("Int. size: VS heuristic")

ggsave("./jUFLP_sizes.png", width = 16, height = 10)
