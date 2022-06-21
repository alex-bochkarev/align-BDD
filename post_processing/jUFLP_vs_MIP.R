library(ggplot2)

df = read.csv("./run_logs/2022-06-21_jUFLP_vs_MIPs.csv", stringsAsFactors = FALSE)
str(df)

ggplot(df)+
  geom_point(aes(x=tMIP, y=tDD_VS))+
  geom_abline(intercept = 0, slope=1.0, color='red')

sum(df$)
