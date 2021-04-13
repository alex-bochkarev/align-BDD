#!/bin/bash

python -m experiments.cUFL_dia_sizes -H > run_logs/cUFL_int_sizes.csv
# for n in 10 15 20 30 35 40; do
python -m experiments.cUFL_dia_sizes -K 10 -n 20 >> run_logs/cUFL_int_sizes.csv
# done
