#!/bin/bash

oname="run_logs/benchmark/bm"
n=20
k=50
PARFLAG=$1
mkdir -p run_logs/benchmark
parallel python -m experiments.tUFL_benchmark_one_size -n $n -K $k -p "{}" ">" "$oname"_{}.csv ::: $(seq 1 $PARFLAG)

logfile="run_logs/benchmark/logfile.csv"
python -m experiments.tUFL_benchmark_one_size -H > $logfile
cat "$oname"_*.csv >> $logfile
