#!/bin/bash

python -m experiments.tUFL_generation_tests -H
for npc in 3 5 7 1; do
    python -m experiments.tUFL_generation_tests -n 15 -K 100 --nodes-per-color $npc --prefix $npc
done
