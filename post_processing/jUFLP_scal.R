######################################################################
##
## Scalability figures (from the ScalTest - jUFLP)
##
## (c) Alexey Bochkarev, Clemson University, 2022
## abochka@g.clemson.edu

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
  library(viridis)})

######################################################################
## parse the input file
infile = "./run_logs/jUFLP_cm.csv"

df = read.csv(infile, stringsAsFactors = FALSE)
df = filter(df, n<15)
## calculate means for each instance size
df_medians = df %>%
  group_by(n) %>%
  summarise(
    mean_MIP_time = median(tMIP),
    mean_DD_toA_time = median(tDD_toA),
    mean_DD_toB_time = median(tDD_toB),
    mean_DD_VS_time = median(tDD_VS))

## percentage where aux problem outperforms
df_outp = df
df_outp$aux_outp= with(df_outp, tDD_VS < tMIP*0.9)

df_outp = df_outp %>%
  group_by(n) %>%
  summarize(
    simpl_outp = sum(aux_outp) / length(aux_outp))

p1 =
ggplot(df)+
  geom_jitter(aes(x=n-0.15, y = tDD_VS,
                 color="VS-heuristic"),
             shape=3, size=5, alpha=0.3, width=0.1)+
  geom_jitter(aes(x=n+0.15, y = tMIP,
                  color="MIP"),
              size=2, alpha=0.3, width = 0.1)+
  geom_line(data=df_medians, aes(x=n, y=mean_DD_VS_time,
                               color="VS-heuristic",
                               linetype="VS-heuristic"), size=2)+
  geom_line(data=df_medians, aes(x=n, y=mean_MIP_time,
                               color="MIP", linetype="MIP"), size=2)+
  geom_label(data=df_outp, aes(x=n, label=scales::percent(accuracy=1,simpl_outp)),
             y=1.5, size=5)+
  scale_linetype_manual(values = c("solid", "dashed"),
                        labels = c("MIP", "VS-heuristic"))+
  scale_x_continuous(
    breaks = seq(min(df$n),max(df$n),by = 1),
    minor_breaks = seq(min(df$n),max(df$n), by=1)
  )+
  scale_color_viridis(discrete=TRUE) +
  annotation_logticks(sides = "l")+
  scale_y_log10(labels=scales::comma)+
  labs(x="No. of variables\n(a)",y="Solution time, sec")+
  guides(color=guide_legend("Solution method (median values):"),
         linetype=guide_legend("Solution method (median values):"))+
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
    panel.grid.minor.x = element_line(size=0.25,linetype = 'solid',
                                      colour = "lightgrey"),
    panel.grid.minor.y = element_line(size=0.25,linetype = 'solid',
                                      colour = "lightgrey"))

######################################################################
## The right panel

Ns = c(7,10, 14)
dfn = filter(df, n %in% Ns)

xmin = min(dfn$tDD_VS)
xmax = quantile(dfn$tMIP, 0.95)

dfn$nlabel = factor(paste0("n = ", dfn$n),
                    levels = c("n = 7", "n = 10", "n = 14"))

p2 =
ggplot(dfn)+
  geom_histogram(aes(x=tDD_toA / tDD_VS, fill="to A"), position="identity",
                 alpha=0.5, binwidth = 0.05)+
  geom_histogram(aes(x=tMIP / tDD_VS, fill="MIP"), position="identity",
                 alpha=0.5, binwidth = 0.05)+
  scale_x_continuous(
    "Runtime, relative to VS-heuristic\n(b)",
    label = scales::percent,
    breaks = seq(xmin, xmax, length.out = 5),
    limits = c(xmin,xmax))+
  scale_y_continuous(
    paste0("Number of instances, out of ",
           nrow(filter(df, n == Ns[1])),"\n"),
    position="right")+
  scale_fill_viridis(discrete=TRUE) +
  theme(
    legend.title = element_text(size=24),
    legend.text = element_text(size=24),
    legend.position = c(0.5,0.9),
    axis.text.x = element_text(size=22,angle=45,vjust = 0.7),
    axis.text.y = element_text(size=11),
    axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26),
    panel.background = element_blank(),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgrey"),
    panel.grid.minor.x = element_line(size=0.25,linetype = 'solid', colour = "lightgrey"),
    panel.grid.minor.y = element_line(size=0.25,linetype = 'solid', colour = "lightgrey"),
    strip.text.y = element_text(size=25, angle=180))+
  guides(fill=guide_legend("Methods based on:"))+
  geom_vline(xintercept = 1, color='red')+
  facet_grid(nlabel ~ ., scales="fixed", switch="y")

if(length(unique(table(dfn$N)))>1){
  cat("WARNING: unbalanced dataset: different number of elements for different Ns in the figure! Please check the data\n")
}

outfile = "figures/jUFLP.eps"
cairo_ps(outfile, width = 16, height = 10, family="Arial")
grid.arrange(p1,p2,ncol=2)
dev.off()
