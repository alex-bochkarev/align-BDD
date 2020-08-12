"""

Generates different BDDs for an independent set instance
(illustration of the fact that the order of variables does
affect the BDD size)

NOTE: Supposed to be run interactively

(c) A. Bochkarev, Clemson University, 2020
a@bochkarev.io
"""

import sys
import numpy as np
import pandas as pd
from importlib import reload
from itertools import permutations

sys.path.append('../..')
import BDD

reload(BDD)

B = BDD.BDD()
B.load("./sample_5var_inst.bdd")

B.show()

print("Order -- size")
min_o = None; max_o = None
min_s = 0; max_s = 0


for o in permutations(B.vars):
    Bp = deepcopy(B)
    Bp.align_to(o, inplace = True)
    print("{} -- {}".format(Bp.vars, Bp.size()))

    if min_o is None or Bp.size() < min_s:
        min_o = o
        min_s = Bp.size()

    if max_o is None or Bp.size() > max_s:
        max_o = o
        max_s = Bp.size()

print("Min order: {} -- |B|={}".format(min_o, min_s))
print("Max order: {} -- |B|={}".format(max_o, max_s))

Bp = B.align_to(max_o)
Bp.dump_gv().view("max_BDD")
B.dump_gv().view("orig_BDD")

print("B is reduced: {}".format(B.is_reduced()))
print("Bp is reduced: {}".format(Bp.is_reduced()))
print("B and Bp are equivalent: {}".format(are_equivalent(B,Bp)))
