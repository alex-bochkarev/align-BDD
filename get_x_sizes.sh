#!/bin/bash
experiment() {
    python -m experiments.UFL_x_sizes "$@"
}

experiment -H
for n in 5 6 7 8 10 11 12 13 14 15 17 20 25; do
    experiment -K 15 -n $n
done
