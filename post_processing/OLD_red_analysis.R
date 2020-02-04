######################################################################
##
## (c) Bochkarev, 2019
## Clemson University, abochka@clemson.edu
library(tidyr)

df = read.csv("../../results/100k_sets/6vars100k_RED.full_log",
              stringsAsFactors = FALSE)
str(df)
df_wide = pivot_wider(df, id_cols = "instance",names_from = "num_type", values_from = "value")

str(df_wide)
sum(df_wide$orig_exact_gsifts_obj < df_wide$orig_bf_obj, na.rm = TRUE)
