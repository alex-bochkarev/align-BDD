#!/bin/bash

python -m experiments.UFL_dia_sizes -H
for n in 3 4 5 6; do
    python -m experiments.UFL_dia_sizes -K 15 -n $n
done
