library(ggplot2)

df5 = read.csv("./run_logs/softcover_goldfish_large_M5.csv", stringsAsFactors = FALSE)
df15 = read.csv("./run_logs/softcover_goldfish_large2_M15.csv", stringsAsFactors = FALSE)
str(df15)

ggplot(df5)+
  geom_point(aes(x=N, y=log(t_BDD/t_MIP)))+
  ggtitle("Cluster size M=5")

ggsave("./figures/adhoc_df5.png", width=16,height = 10)

ggplot(df15)+
  geom_point(aes(x=N, y=log(t_BDD/t_MIP)))+
  ggtitle("Cluster size M=15")

ggsave("./figures/adhoc_df15.png", width=16,height = 10)


df = read.csv("./run_logs/softcover_goldfish_large.csv", stringsAsFactors = FALSE)
str(df)

ggplot(df)+
  geom_point(aes(x=N, y=log(t_BDD/t_MIP)))+
  ggtitle("Misc")
