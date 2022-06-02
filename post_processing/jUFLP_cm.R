library(ggplot2)
library(dplyr)
library(viridis)

df = read.csv("./run_logs/jUFLP_cm.csv", stringsAsFactors = FALSE)
str(df)


df$isize_VS_rel = df$int_VS / pmin(df$int_toA, df$int_toB)
df$tVS_rel = df$tDD_VS / (df$tDD_toA + df$tDD_toB)
df$tVS_relA = df$tDD_VS / df$tDD_toA
df$tVS_relB = df$tDD_VS / df$tDD_toB

ggplot(filter(df, n %in% c(5, 10, 14)))+
  geom_histogram(aes(x=tDD_VS / tMIP))+
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
    strip.background = element_blank())+
  geom_vline(xintercept = 1.0, color="red")+
  xlab("tVS / tMIP")+
  facet_wrap(n ~ .)

ggplot(df)+
  geom_histogram(aes(x=tVS_rel, color="rel to (toA or toB)"))+
  geom_histogram(aes(x=tVS_relA, color="rel to toA"))+
  geom_histogram(aes(x=tVS_relB, color="rel to toB"))+
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
    strip.background = element_blank())+
  xlab("Intersection DD sizes: VS / best of 'to-A' and 'to-B' alignment.")+
  xlim(c(0, 2.0))

with(filter(df, n > 13),
     sum(isize_VS_rel < 1.0) / sum(isize_VS_rel > 0.0))

nrow(df)


df_means = df %>%
  group_by(n) %>%
  summarise(
    median_VS_time = median(tDD_VS / tMIP),
    median_toA_time = median(tDD_toA / tMIP)
    ## median_toB_time = median(tDD_toB),
    ## median_MIP_time = median(tMIP)
  )

ggplot(df)+
  geom_jitter(aes(x=n-0.15, y=tDD_VS / tMIP, color="VS heuristic"), shape=3, size=5, alpha=0.3, width=0.1)+
  geom_jitter(aes(x=n+0.15, y=tDD_toA / tMIP, color="to-A heuristic"), shape=1, size=5, alpha=0.3, width=0.1)+
  ## geom_jitter(aes(x=n, y=tMIP, color="MIP"), shape=2, size=5, alpha=0.3, width=0.1)+
  scale_x_continuous(
    breaks = seq(min(df$n),max(df$n),by = 1),
    minor_breaks = seq(min(df$n),max(df$n), by=1)
  )+
  scale_color_viridis(discrete=TRUE) +
  annotation_logticks(sides = "l")+
  scale_y_log10(labels=scales::comma)+
  labs(x="No. of variables",y="Solution time, relative to MIP.")+
  geom_line(data=df_means, aes(x=n, y=median_VS_time,
                               color="VS heuristic",
                               linetype="VS heuristic"), size=2)+
  ## geom_line(data=df_means, aes(x=n, y=median_MIP_time,
  ##                              color="MIP",
  ##                              linetype="MIP"), size=2)+
  geom_line(data=df_means, aes(x=n, y=median_toA_time,
                               color="to-A heuristic",
                               linetype="to-A heuristic"), size=2)+
  guides(color=guide_legend("Solution method (median values):"),
         linetype=guide_legend("Solution method (median values):"))+
  geom_hline(yintercept = 1.0, color="red")+
  theme(
    legend.position = c(0.4, 0.8),
    legend.direction = "vertical",
    legend.title = element_text(size=16),
    legend.text = element_text(size=16),
    legend.key = element_blank(),
    legend.key.width = unit(2.5,"cm"),
    axis.text.x = element_text(size=18,angle=45,vjust = 0.7),
    axis.text.y = element_text(size=10),
    axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26),
    panel.background = element_blank(),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgrey"),
    panel.grid.minor.x = element_line(size=0.25,linetype = 'solid', colour = "lightgrey"),
    panel.grid.minor.y = element_line(size=0.25,linetype = 'solid', colour = "lightgrey"))
