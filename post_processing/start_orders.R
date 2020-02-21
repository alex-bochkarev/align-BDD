library(ggplot2)

df = read.csv("../run_logs/start_orders.log", skip = 7)
str(df)

ggplot(df)+
    geom_jitter(aes(x=objective, y=MIN_I_OPT), width = 0.3, alpha=0.5,color="blue")
