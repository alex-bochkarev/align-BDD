#!/bin/bash

python -m experiments.tUFL_generation_tests -H
for p in 0.1 0.2 0.4 0.5 0.6 0.7 0.8 0.9; do
    python -m experiments.tUFL_generation_tests -n 15 -K 100 -P $p -p $p
done
