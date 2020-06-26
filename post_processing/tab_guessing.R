######################################################################
## Creates ggplot for the simplified problem:
## solution guessing quality with different
## heuristics
##
## (c) Alexey Bochkarev, Clemson University, 2020

library(xtable)
library(dplyr)
library(tidyr)
library(latex2exp)
library(optparse)

## Internal parameters for the figure
SIMPL_BASE_COL = "simpl_BB_obj" # exact solution
## Heuristics to show (codes)
SHOW_HEU = c("all")

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

df_legend = select(filter(df, num_type=="legend"), value,comment)
df = filter(df, num_type != "legend")
df$value = as.numeric(df$value)

df_wide = pivot_wider(df,
                      id_cols = "instance",
                      names_from = "num_type",
                      values_from = "value"
                      )

simpl_obj_cols = grep("^simpl.*_obj$",colnames(df_wide),value = TRUE)
time_cols = grep(".+time$",colnames(df_wide),value = TRUE)

if (SHOW_HEU[1] == "all"){
    SHOW_HEU = c(grep("^simpl.*_obj$",simpl_obj_cols,value = TRUE))
}

guessing_power = data.frame(heuristic = grep(".*obj$",SHOW_HEU,value = TRUE),stringsAsFactors = FALSE)

guessing_power$exact_guesses = sapply(guessing_power$heuristic, function(h){
    sum(df_wide[[h]] == df_wide[[ SIMPL_BASE_COL ]])
})

guessing_power$good_guesses = sapply(guessing_power$heuristic, function(h){
    sum(df_wide[[h]] <= df_wide[[ SIMPL_BASE_COL ]]*opt$allowance)
})

## guessing_power = pivot_longer(guessing_power, cols = c("good_guesses","exact_guesses"), values_to = "no_insts", names_to = "guess_type")

guessing_power = guessing_power %>%
    mutate(
        heuristic = substring(heuristic, 1, nchar(heuristic)-nchar("_obj"))
    )

guessing_power = merge(x=guessing_power, y=df_legend, by.x = "heuristic", by.y = "value")

######################################################################
## draw the figure
print(xtable(select(guessing_power, c("comment","exact_guesses", "good_guesses")), type = "latex"), file=opt$out)
## plt =
##     ggplot(guessing_power,
##            aes(x=comment, y =no_insts, fill=guess_type))+
##     geom_bar(stat = "identity",position="dodge")+
##     geom_text(aes(label=no_insts),vjust=-0.25, position = position_dodge(width = 0.9), size=7)+
##     geom_hline(yintercept = nrow(df_wide), color="red",linetype="dotted", size=2)+
##     labs(y="Number of instances", x="Heuristic",fill="Objective value:")+
##     scale_fill_manual(values = c(
##                           "exact_guesses" = "blue",
##                           "good_guesses" = "lightblue"
##                       ),
##                       labels = c(
##                           "exact_guesses" = "optimal",
##                           "good_guesses" = paste("up to",opt$allowance,"of the optimum",sep=" ")
##                       ),
##                       )+
##     theme(
##         legend.position = c(0.2, 0.8),
##         legend.direction = "vertical",
##         legend.title = element_text(size=24),
##         legend.text = element_text(size=24),
##         axis.text.x = element_text(size=22,angle=45,vjust = 0.7),
##         axis.text.y = element_text(size=22),
##         axis.title.x = element_text(size = 26),
##         axis.title.y = element_text(size = 26),
##         panel.background = element_blank()
##           )


## ggsave(opt$out, plt, device=cairo_ps, width = 16, height = 10)
