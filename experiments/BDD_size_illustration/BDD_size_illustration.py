"""

Generates different BDDs for an independent set instance
(illustration of the fact that the order of variables does
affect the BDD size)

(c) A. Bochkarev, Clemson University, 2020
a@bochkarev.io
"""

import sys
import numpy as np
import pandas as pd
from importlib import reload

from copy import deepcopy


sys.path.append('../..')
import BDD

def get_value(B, x):
    """Finds the terminal node (T or F) corresponding to the variable choices in x

    Args:
        B: a BDD
        x: a dict of {var_name: value}, where value is in {0,1}.

    Returns:
        True/False/-1: terminal node corresponding to the path implied by x, encoded as Boolean (or -1 if error)
    """

    (node,) = B.layers[0]

    for i in range(len(x)):
        if x[B.vars[i]] == 0:
            node = node.lo
        elif x[B.vars[i]] == 1:
            node = node.hi
        else:
            print("Wrong value ({}) for variabe {}: 0 or 1 expected".format(x[B.vars[i]], B.vars[i]))
            return -1

    if node.id == B.T.id:
        return True
    elif node.id == B.F.id:
        return False
    else:
        print("Error not a True or False node reached!")
        print("node is {}, T is {}, F is {}".format(node.id, B.T.id, B.F.id))
        return -1

# x = dict(zip([x for x in "12345"], [0,0,0,0,0]))
# print("Value for {} is {}".format(x, get_value(B,x)))

def get_truth_table(B):
    """Prints a truth table for the Boolean function defined by BDD B"""
    tt = []
    ind = []
    for x_num in range(2**len(B)):
        x = [int(j) for j in np.binary_repr(x_num, width = len(B.vars))]
        tt.append(x + [ get_value(B,dict( zip(B.vars, x) )) ])
        ind.append(np.binary_repr(x_num, width = len(B.vars)))

    tt = pd.DataFrame(tt, columns = B.vars + ["Value"], index = ind)
    return tt

# print("")
# tt = get_truth_table(B)

def are_equivalent(A,B):
    """Checks if two BDDs (A and B) are equivalent in the sense that they define the same Boolean function,
    by checking if the corresponding truth tables coincide. Returns Bool
    """

    tt_A = get_truth_table(A); tt_B = get_truth_table(B)
    tt_B = tt_B[A.vars + ['Value']]
    tt_B['new_index'] = [c for c in tt_B.apply(lambda row: "".join([str(x) for x in row[:-1]]), axis=1) ]
    tt_B.set_index('new_index', inplace=True)

    for idx, row in tt_A.iterrows():
        if row['Value'] != tt_B.loc[idx]['Value']:
            print("Discrepancy found: for x={}, value(A)={}, value(B)={}".format(idx, row['Value'], tt_B.loc[idx]['Value']))
            return False

    return True

reload(BDD)

B = BDD.BDD()
B.load("./sample_5var_inst.bdd")

def test_BDD_transformations(B,no_iters=250):
    # let's test the basic things
    print("First: for Bp = B,")
    Bp = deepcopy(B)
    print("Variables:\n B:{}\nBp:{}".format(B.vars,Bp.vars))
    eq = are_equivalent(B,Bp))
    print("B ~ Bp: {}".format(eq))
    if not eq:
          print("Non-equivalent BDDs: STOP")
          exit(-1)

    N = len(B.vars)
    ##
    print("Randomly sifting variables...")
    for i in range(no_iters):
        Bp.sift(B.vars[np.random.randint(N)], np.random.randint(N))
        print(".",end="")

    print("Variables:\n B:{}\nBp:{}".format(B.vars,Bp.vars))
    eq = are_equivalent(B,Bp))
    print("B ~ Bp: {}".format(eq))
    if not eq:
          print("Non-equivalent BDDs: STOP")
          exit(-1)

    ##
    print("Randomly aligning variables (inplace)...")
    for i in range(no_iters):
        Bp.align_to(np.random.permutation(Bp.vars),inplace=True)
        print(".",end="")

    print("Variables:\n B:{}\nBp:{}".format(B.vars,Bp.vars))
    eq = are_equivalent(B,Bp))
    print("B ~ Bp: {}".format(eq))
    if not eq:
          print("Non-equivalent BDDs: STOP")
          exit(-1)

    ##

    print("Randomly aligning variables (non inplace)...")
    for i in range(no_iters):
        Bp = Bp.align_to(np.random.permutation(Bp.vars))
        print(".",end="")

    print("Variables:\n B:{}\nBp:{}".format(B.vars,Bp.vars))
    eq = are_equivalent(B,Bp))
    print("B ~ Bp: {}".format(eq))
    if not eq:
          print("Non-equivalent BDDs: STOP")
          exit(-1)

for i in range(100):
    A = BDD.BDD.random(10)
    test_BDD_transformations(A)
