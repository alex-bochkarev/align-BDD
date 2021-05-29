#!/bin/bash

odir="run_logs/jUFL"
k=50

PARFLAG=$1

mkdir -p $odir
parallel -j $PARFLAG python -m experiments.jUFL_hist_sizes -n {3} -K $k -p {2} --prefix "{1}_{2}_{3}" "|" tee $odir/tmp.{3}_{2}_{1}.csv ::: $(seq 1 20) ::: 0.2 0.5 0.8 ::: 15 20 25

python -m experiments.jUFL_hist_sizes -H > $odir/logfile.csv

cat $odir/tmp.* >> $odir/logfile.csv
rm $odir/tmp.*
