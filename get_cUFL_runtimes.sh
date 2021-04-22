#!/bin/bash

python -m experiments.cUFL_runtimes -H
for n in 20 25 30; do
    python -m experiments.cUFL_runtimes -K 5 -n $n
done
