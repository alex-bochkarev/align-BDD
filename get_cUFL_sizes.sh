#!/bin/bash

python -m experiments.cUFL_dia_sizes -H
for n in 20 25 30; do
    python -m experiments.cUFL_dia_sizes -K 10 -n $n
done
