library(ggplot2)
library(dplyr)

df = read.csv("./run_logs/darkcloud_BDD_vs_MIP.csv", stringsAsFactors = FALSE)

str(df)
df$tMIPvstDD = df$t_MIP / df$t_BDD_

ggplot(data=df)+
  geom_jitter(aes(x=total_nodes, y=log(t_MIP / t_BDD_), color=L), width=5)


with(filter(df, M > 10), sum(tMIPvstDD > 1.0) / length(tMIPvstDD))
