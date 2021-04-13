#!/bin/bash

python -m experiments.cUFL_runtimes -H  > run_logs/cUFL_runtimes.csv
python -m experiments.cUFL_runtimes -K 10 -n 20 >> run_logs/cUFL_runtimes.csv
