#!/bin/bash

python -m experiments.cUFL_dia_sizes -H
for n in 5 7 10 15 20 25; do
    python -m experiments.cUFL_dia_sizes -K 5 -n $n
done
