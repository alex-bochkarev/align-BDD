library(ggplot2)

## df = read.csv("./run_logs/2023-01-03_jUFLP_simscores.csv",
##               stringsAsFactors = FALSE)

df = read.csv("./run_logs/2022-12-06_jUFLP_simscores.csv",
              stringsAsFactors = FALSE)

ggplot(df)+
  geom_jitter(aes(x = param, y = tDD_VS / tDD_toA), width=0.05, size=3)+
  theme(
    axis.text.x = element_text(size=22),
    axis.text.y = element_text(size=22),
    axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26),
    panel.background = element_blank(),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgrey"),
    panel.grid.minor.x = element_line(size=0.5,linetype = 'solid', colour = "lightgrey"),
    panel.grid.minor.y = element_line(size=0.5,linetype = 'solid', colour = "lightgrey"))+
  xlab("Linking order parameter")+
  ylab("t (DD with VS) / t (DD with 'align-to-A')")+
  facet_grid(M ~ .)
    ## end of styling

ggsave("figures/fig_jUFLP_simscores.eps", width = 16, height = 10)

ggplot(df)+
  geom_jitter(aes(x = param, y = CPP_simscore), width=0.05, size=3)+
  theme(
    axis.text.x = element_text(size=22),
    axis.text.y = element_text(size=22),
    axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26),
    panel.background = element_blank(),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgrey"),
    panel.grid.minor.x = element_line(size=0.5,linetype = 'solid', colour = "lightgrey"),
    panel.grid.minor.y = element_line(size=0.5,linetype = 'solid', colour = "lightgrey"))
## end of styling

ggplot(df)+
  geom_point(aes(x=CPP_simscore, y=tDD_VS / tDD_toA, color=param), size=3)+
  theme(
    axis.text.x = element_text(size=22),
    axis.text.y = element_text(size=22),
    axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26),
    panel.background = element_blank(),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgrey"),
    panel.grid.minor.x = element_line(size=0.5,linetype = 'solid', colour = "lightgrey"),
    panel.grid.minor.y = element_line(size=0.5,linetype = 'solid', colour = "lightgrey"))
## end of styling
