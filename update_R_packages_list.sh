#!/bin/bash

grep -e "library" post_processing/*.R | sed 's/^[^:]*://' | sort | uniq | sed 's/library(/"/' | sed 's/)/"/' > aux/packages.list
