library(ggplot2)

df = read.csv("./run_logs//darkcloud_rnd_cover.csv", stringsAsFactors = FALSE)

str(df)

df$int_size_VS_rel = df$size_int_VS / df$size_int_toC

ggplot(df)+
  geom_histogram(aes(x=int_size_VS_rel))+
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
  xlab("Intersection DD sizes: VS / 'to-A' alignment.")

sum(df$int_size_VS_rel < 1.0) / sum(df$int_size_VS_rel > 0.0)
sum(df$int_size_VS_rel == 1.0) / sum(df$int_size_VS_rel > 0.0)

nrow(df)
ggsave("./reports/2022-05-20_special_classes/int_sizes_hist.png",
       width = 16, height = 10)

ggplot(df)+
  geom_point(aes(x=log(size_int_toC), y=log(size_int_VS)))+
  geom_abline(intercept = 0.0, slope=1.0, color='red')+
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
    strip.background = element_blank())

ggsave("./reports/2022-05-20_special_classes/int_sizes_points.png",
       width = 16, height = 10)
