library(ggplot2)
library(dplyr)

df = read.csv("./run_logs/darkcloud_BDD_vs_MIP_longMIP.csv", stringsAsFactors = FALSE)

str(df)
df$tMIPvstDD = df$t_MIP / df$t_BDD_

ggplot(data=df)+
  geom_point(aes(x=t_BDD_, y=t_MIP))+
  geom_abline(slope=1.0, intercept = 0.0, color='red')+
  xlab("Solving with BDDs, sec.")+
  ylab("Solving with MIP, sec.")+
  ggtitle("untyped UFLP: runtimes CPP/BDD+SP vs naive MIP")+
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
  )

ggsave("BDDvsMIP_untyped.png", width = 16, height=10)


dft = read.csv("./run_logs/dclouds_typed.csv", stringsAsFactors = FALSE)

str(dft)
df$tMIPvstDD = df$t_MIP / df$t_BDD_

ggplot(data=dft)+
  geom_point(aes(x=tTDD, y=tTMIP, color='typed version'))+
  geom_point(aes(x=tUDD, y=tUMIP, color='untyped version'))+
  geom_abline(slope=1.0, intercept = 0.0, color='red')+
  xlab("Solving with BDDs, sec.")+
  ylab("Solving with MIP, sec.")+
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
  )

ggsave("BDDvsMIP_wtyped.png", width = 16, height=10)
