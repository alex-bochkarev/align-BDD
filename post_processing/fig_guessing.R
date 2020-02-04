######################################################################
## Creates ggplot for the simplified problem:
## solution guessing quality with different
## heuristics
##
## (c) Alexey Bochkarev, Clemson University, 2020

library(ggplot2)
library(dplyr)
library(tidyr)
library(latex2exp)
library(optparse)

######################################################################
## unpack the command line arguments
option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL,
                help="input filename (solved instances log)", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="./out.eps",
                help="output file name [default= %default]", metavar="character"),
    make_option(c("-a","--allowance"),type="numeric",
                help="a multiple for 'good enough' solution [default = %default]", default=1.10)
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input) | is.null(opt$out)){
    print_help(opt_parser)
    stop("Please specify both input and output files", call.=FALSE)
}

######################################################################
## parse the input file

df = read.csv(opt$input, stringsAsFactors = FALSE)

df_wide = pivot_wider(df,
                      id_cols = "instance",
                      names_from = "num_type",
                      values_from = "value"
                      )

simpl_colnames = grep("^simpl.*_obj$",colnames(df_wide),value = TRUE)

guessing_power = data.frame(heuristic = grep(".*obj$",simpl_colnames,value = TRUE),stringsAsFactors = FALSE)
guessing_power$exact_guesses = sapply(guessing_power$heuristic, function(h){
    sum(df_wide[[h]] == df_wide[[ "simpl_obj" ]])
})

guessing_power$good_guesses = sapply(guessing_power$heuristic, function(h){
    sum(df_wide[[h]] <= df_wide[["simpl_obj"]]*opt$allowance)
})

guessing_power = pivot_longer(guessing_power, cols = c("good_guesses","exact_guesses"), values_to = "no_insts", names_to = "guess_type")

######################################################################
## draw the figure
plt =
    ggplot(filter(guessing_power, heuristic %in% c(
                                                 "simpl_obj",
                                                 "simpl_heu_gswaps_obj",
                                                 "simpl_heu_gsifts_obj",
                                                 "simpl_heu_g2sifts_obj",
                                                 "simpl_toAB_obj",
                                                 "simpl_toR_obj"
                                                 )),
           aes(x=heuristic, y =no_insts, fill=guess_type))+
    geom_bar(stat = "identity",position="dodge")+
    geom_text(aes(label=no_insts),vjust=-0.25, position = position_dodge(width = 0.9), size=7)+
    geom_hline(yintercept = nrow(df_wide), color="red",linetype="dotted", size=2)+
    labs(y="Number of instances", x="Heuristic",fill="Solution quality:")+
    scale_x_discrete(labels=c(
                         "simpl_obj" = "Branch & bound",
                         "simpl_heu_gswaps_obj" = "G. swaps",
                         "simpl_heu_gsifts_obj" = "G. sifts",
                         "simpl_heu_g2sifts_obj" = "G. sift pairs",
                         "simpl_toAB_obj" = parse(text = TeX("Best of $S_A$, $S_B$")),
                         "simpl_toR_obj" = "Random order"
                     ))+
    scale_fill_manual(values = c(
                          "exact_guesses" = "blue",
                          "good_guesses" = "lightblue"
                      ),
                      labels = c(
                          "exact_guesses" = "an optimum",
                          "good_guesses" = "within 10% of optimum"
                      ),
                      )+
    theme(
        legend.position = c(0.8, 0.8),
        legend.direction = "vertical",
        legend.title = element_text(size=24),
        legend.text = element_text(size=24),
        axis.text.x = element_text(size=22,angle=45,vjust = 0.7),
        axis.text.y = element_text(size=22),
        axis.title.x = element_text(size = 26),
        axis.title.y = element_text(size = 26),
        panel.background = element_blank()
          )


ggsave(opt$out, plt, device=cairo_ps, width = 16, height = 10)
