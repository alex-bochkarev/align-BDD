######################################################################
# Logfiles testing script
# tests BB_bounds and solved_R logfiles (for consistency)
# (c) A. Bochkarev, Clemson University, 2020
######################################################################

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(optparse))

######################################################################
## unpack the command line arguments
option_list = list(
  make_option(c("-b", "--BBlog"), type="character", default="run_logs/BB_bounds_R.log",
              help="branch-and-bound log [default= %default]", metavar="character"),
  make_option(c("-s", "--sollog"), type="character", default="run_logs/solved_R.log",
              help="general solution blog [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$BBlog) | is.null(opt$sollog)){
  print_help(opt_parser)
  stop("Please specify both files", call.=FALSE)
}

######################################################################
## parse the input file
BB_log = opt$BBlog
sol_log = opt$sollog
BBs = read.csv(BB_log, stringsAsFactors = FALSE)
sols= read.csv(sol_log, stringsAsFactors = FALSE)
sols = filter(sols, instance >= 0 & num_type == "simpl_BB_obj")
sols$num_type <- NULL
sols$comment <- NULL
sols$value = as.numeric(sols$value)

BBopts = filter(BBs, num_type == "steplog" & comment == "status:optimal")
cat("Checking correctness of 'BB_bounds' logfile vs 'solved' logfile...\n")

## assertion: LB = UB
neqs = filter(BBopts, LB != UB)
if (nrow(neqs) > 0){
  cat(paste("LB != UB found:", nrow(neqs),"instances.\n"))
  neqs
} else {
  cat("BB_bounds log: [ OK ] All LBs equal UBs for optimal points\n")
}

cross_check = merge(x=BBopts, y=sols, by="instance")

neqs = filter(cross_check, UB != value)
if (nrow(neqs) > 0){
  cat(paste("'UB != optimum' found:", nrow(neqs),"instances.\n"))
  neqs
} else {
  cat("BB_bounds vs solved log: [ OK ] optima are consistent between the two logfiles\n")
}
