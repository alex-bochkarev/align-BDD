#!/bin/bash

python -m experiments.cUFL_dia_sizes -H
for n in 10 15 20; do
    python -m experiments.cUFL_dia_sizes -K 3 -n $n
done
