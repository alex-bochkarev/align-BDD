library(ggplot2)

df = read.csv("./run_logs/2022-07-07_DDvsMIP.csv", stringsAsFactors = FALSE)

ggplot(df) +
  geom_jitter(aes(y=tMIP / tDD_VS, x = tMIP_CPP / tDD_VS, color=as.factor(M), shape=as.factor(M)),
              width=0.2, size=3)+
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
    strip.background = element_blank()
  )+
  ylab("t(naive MiP) / t(DD with VS)")+
  xlab("t(MiP CPP) / t(DD with VS)")+
  guides(color=guide_legend("M"), shape=guide_legend("M"))

ggsave("figures/tMIP_tMIPCPP_tDD.png", width = 16, height = 10)

ggplot(df) +
  geom_jitter(aes(y=int_VS, x = int_VS_toA, color=as.factor(M), shape=as.factor(M)),
              width=0.2, size=3)+
  geom_abline(intercept = 0.0, slope = 1.0, color='red')+
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
  ylab("Intersection DD size: for VarSeq method.")+
  xlab("Intersection DD size: for 'align-to-A.'")+
  scale_y_continuous(labels = scales::comma)+
  scale_x_continuous(labels = scales::comma)+
  guides(color=guide_legend("M"), shape=guide_legend("M"))

ggsave("figures/intDD_VS_vs_toA.png", width = 16, height = 10)


ggplot(df) +
  geom_point(aes(x=tDD_VS, y = tDD_toA, color=as.factor(M), shape=as.factor(M)),
              width=0.2, size=3)+
  geom_abline(intercept = 0.0, slope = 1.0, color='red')+
  xlim(c(0, 400))+
  ylim(c(0, 450))+
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
  ylab("t(DD with 'toA'), sec")+
  xlab("t(DD with VS), sec")+
  guides(color=guide_legend("M"), shape=guide_legend("M"))

ggsave("figures/t_VS_vs_toA.png", width = 10, height = 10)
