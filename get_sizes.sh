#!/bin/bash

python ./UFL.py -H
for n in 3 4 5 6; do
    python ./UFL.py -K 15 -n $n
done
