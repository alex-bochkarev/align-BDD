######################################################################
##
## Analyses the varsplit effects
## (based on the log file from the varsplit_test.py file)
##
##
## (c) Alexey Bochkarev, Clemson University, 2019

library(ggplot2)
library(dplyr)

sps = read.csv("../results/VStest.log",stringsAsFactors = FALSE)
str(sps)

pl1 =
    ggplot(filter(sps, num_type == "no_splits"))+
    geom_histogram(aes(x=value),fill="blue")+
    geom_hline(yintercept = nrow(filter(sps, num_type == "no_splits")), color="red",linetype="dotted", size=0.3)+
    stat_bin(binwidth = 1, geom='text', color='black',size=8, aes(x=value, label=..count..),
             position=position_stack(vjust = 1.05))+
    ggtitle("Varsplit property test: distribution of the number of sub-instances generated from an initial instance (100k, 15-var non-reduced instances)")+
    scale_x_continuous(
        "No. of sub-instances generated",
        labels = scales::number_format(accuracy = 1.0),
        breaks = seq(1,10,by=1)
    )+
    ylab("Count")

pl2 =
    ggplot(filter(sps, num_type == "split_len"))+
    geom_histogram(aes(x=value),fill="blue")+
    stat_bin(binwidth = 1, geom='text', color='black', size=8, aes(x=value, label=..count..),
             position=position_stack(vjust = 1.05))+
    ggtitle("Varsplit property test: lengths of sub-instances (generated from 100k 15-var non-reduced instances)")+
    scale_x_continuous(
        "Sub-instance length, variables",
        labels = scales::number_format(accuracy = 1.0),
        breaks = seq(1,15,by=1)
    )+
    ylab("Count of such sub-instances")
