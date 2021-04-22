#!/bin/bash

python -m experiments.cUFL_runtimes -H
for n in $(seq 5 21); do
    python -m experiments.cUFL_runtimes -K 100 -n $n --with-gsifts
done
