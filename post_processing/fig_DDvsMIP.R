library(ggplot2)
library(viridis)
library(dplyr)

df = read.csv("./run_logs/2022-07-19_jUFLP.csv", stringsAsFactors = FALSE)

df = filter(df, experiment <= 30)  # just to have 10 points per size

# Figure: relative runtimes
ggplot(df) +
  geom_jitter(aes(y=tMIP / tDD_VS, x = tMIP_CPP / tDD_VS, shape=as.factor(M)),
              width=0.2, size=7)+
  scale_shape_manual(values = c(1,4,16)) +
  scale_color_viridis(discrete = TRUE)+
  geom_abline(intercept = 1.0, slope = 0.0, color='red')+
  geom_vline(xintercept = 1.0, color='red')+
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
    strip.background = element_blank(),
    legend.position = c(0.8,0.9),
    legend.text = element_text(size=25),
    legend.title = element_text(size=25)
  )+
  ylab("t(naive MiP) / t(DD with VS)")+
  xlab("t(CPP MIP) / t(DD with VS)")+
  guides(shape=guide_legend("M, points per cluster"))

ggsave("figures/tMIP_tMIPCPP_tDD.eps", width = 16, height = 10)

## Figure: intersection diagram sizes
ggplot(df) +
  geom_jitter(aes(y=int_VS, x = int_VS_toA, shape=as.factor(M)),
              width=0.2, size=7)+
  scale_shape_manual(values = c(1,4,16)) +
  geom_abline(intercept = 0.0, slope = 1.0, color='red')+
  theme(
    axis.text.x = element_text(vjust = 0.5, hjust = 1, size=20),
    axis.text.y = element_text(size=20),
    axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26, margin = margin(t=50)),
    panel.background = element_rect(fill = NA, color = "black"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    color = "lightgrey"),
    strip.text.x = element_text(size = 22),
    strip.text.y = element_text(size=22),
    strip.background = element_blank(),
    legend.position = c(0.25,0.9),
    legend.text = element_text(size=30),
    legend.title = element_text(size=30)
  )+
  ylab("Intersection DD size: DD with VS (VarSeq).")+
  xlab("Intersection DD size: for 'align-to-A.'")+
  annotation_logticks()+
  scale_x_log10(labels=scales::comma)+
  scale_y_log10(labels=scales::comma)+
  guides(shape=guide_legend("M, points / cluster"))

ggsave("figures/intDD_VS_vs_toA.eps", width = 10, height = 10)

# Figure: toA vs VS runtimes
ggplot(df) +
  geom_point(aes(x=tDD_VS, y = tDD_toA, shape=as.factor(M)),
              width=0.2, size=7)+
  scale_shape_manual(values = c(1,4,16)) +
  geom_abline(intercept = 0.0, slope = 1.0, color='red')+
  ## xlim(c(0, 400))+
  ## ylim(c(0, 450))+
  theme(
    axis.text.x = element_text(vjust = 0.5, hjust = 1, size=20),
    axis.text.y = element_text(size=20),
    axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26, margin = margin(t=50)),
    panel.background = element_rect(fill = NA, color = "black"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    color = "lightgrey"),
    strip.text.x = element_text(size = 22),
    strip.text.y = element_text(size=22),
    strip.background = element_blank(),
    legend.position = c(0.25,0.9),
    legend.text = element_text(size=30),
    legend.title = element_text(size=30)
  )+
  ylab("t(DD with 'align-to-A'), sec")+
  xlab("t(DD with VS), sec")+
  guides(shape=guide_legend("M, points / cluster"))

ggsave("figures/t_VS_vs_toA.eps", width = 10, height = 10)
