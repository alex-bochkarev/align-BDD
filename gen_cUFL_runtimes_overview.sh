#!/bin/bash

python -m experiments.cUFL_runtimes -H
for n in $(seq 5 25); do
    python -m experiments.cUFL_runtimes -K 100 -n $n --with-gsifts
done
