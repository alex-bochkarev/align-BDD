"""Test unit: equivalence of BDDs after operations.

Tests if BDDs encode the same Boolean functions after the manipulations (swap, sift, align_to)
against randomly generated instances (not necessarily unique)

(c) A. Bochkarev, Clemson University, 2020
a@bochkarev.io

"""

import sys
import numpy as np
import pandas as pd
import argparse as ap
from copy import deepcopy
from time import time

sys.path.append('..')
import BDD as BDD

N = 8
NO_SWAPS = 20
NO_SIFTS = 20
NO_ALIGNS = 20
NO_GSIFTS = 20
TIMEOUT = 600 # in seconds
PROGRESS_STEPS = 100
P_STEP = TIMEOUT / PROGRESS_STEPS

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

    msg = None
    tt_A = get_truth_table(A); tt_B = get_truth_table(B)
    tt_B = tt_B[A.vars + ['Value']]
    tt_B['new_index'] = [c for c in tt_B.apply(lambda row: "".join([str(x) for x in row[:-1]]), axis=1) ]
    tt_B.set_index('new_index', inplace=True)

    for idx, row in tt_A.iterrows():
        if row['Value'] != tt_B.loc[idx]['Value']:
            msg = "\nINFO: Discrepancy found. For x={}, value(A)={}, value(B)={}".format(idx, row['Value'], tt_B.loc[idx]['Value'])
            return [False, msg]

    return [True, msg]

def test_BDD_transformations(B,no_iters=250):
    """Implements a simple ad-hoc test"""
    # let's test the basic things
    print("First: for Bp = B,")
    Bp = deepcopy(B)
    print("Variables:\n B:{}\nBp:{}".format(B.vars,Bp.vars))
    eq, _ = are_equivalent(B,Bp)
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
    eq, _ = are_equivalent(B,Bp)
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
    eq, _ = are_equivalent(B,Bp)
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
    eq, _ = are_equivalent(B,Bp)
    print("B ~ Bp: {}".format(eq))
    if not eq:
          print("Non-equivalent BDDs: STOP")
          exit(-1)

def adhoc_test():
    for i in range(100):
        A = BDD.BDD.random(10)
        test_BDD_transformations(A)

def print_pair(B,Bp):
    """Prints details on a pair of BDDs"""

    print("Var(B)={} with weights {}".format(B.vars, [B.n(i) for i in range(len(B))]))
    print("Var(Bp)={} with weights {}".format(Bp.vars, [Bp.n(i) for i in range(len(Bp))]))
    print("Profiles: B:\n{}\nBp:\n{}".format(B.profile(), Bp.profile()))

if __name__ == '__main__':
    print(f"Running BDD equivalence test for {N} variables ({PROGRESS_STEPS} progress steps)")
    print(f"Approximate wallclock: {TIMEOUT} sec.")
    print(f"Per instance: {NO_SWAPS} swaps, {NO_SIFTS} sifts, {NO_ALIGNS} aligns inplace and not inplace, {NO_GSIFTS} greedy sifts with random BDDs")
    print("[Ctrl-C] to interrupt")

    time_ctr = 0
    test_ctr = 0
    t0 = time()
    t1 = time()

    print("["+"".join(["⋅" for _ in range(PROGRESS_STEPS)])+"]")
    print("[", end="")
    while t1 - t0 < TIMEOUT:
        B = BDD.BDD.random(N)

        # trivial tests first
        # self-equivalence
        Bp = deepcopy(B)
        eq, msg = are_equivalent(B,Bp)
        if not eq:
            print("ERROR: B not equivalent to Bp=B!")
            print(msg)
            print_pair(B,Bp)

        # non-eq to another random BDD
        keep_trying = True
        Bp = BDD.BDD.random(N)

        while keep_trying:
            eq, msg = are_equivalent(B,Bp)
            if eq:
                Bp.make_reduced()
                B.make_reduced()
                if B.profile() != Bp.profile():
                    print("ERROR: B is equivalent to another random BDD with a different profile!")
                    print_pair(B,Bp)
                else:
                    Bp = BDD.BDD.random(N)
            else:
                keep_trying = False

        ## Random swaps
        Bp = deepcopy(B)
        for i in range(NO_SWAPS):
            Bp.swap_up(np.random.randint(1,N))

        eq, msg = are_equivalent(B,Bp)
        if not eq:
            print("ERROR: B is not equivalent to Bp after random swaps")
            print(msg)
            print_pair(B,Bp)

        ## Random sifts
        Bp = deepcopy(B)
        for i in range(NO_SIFTS):
            Bp.sift(B.vars[np.random.randint(N)], np.random.randint(N))

        eq, msg = are_equivalent(B,Bp)
        if not eq:
            print("ERROR: B is not equivalent to Bp after random sifts")
            print(msg)
            print_pair(B,Bp)

        ## Random aligns (not inplace)
        Bp = deepcopy(B)
        for i in range(NO_ALIGNS):
            Bp = Bp.align_to(np.random.permutation(B.vars), inplace=False)

        eq, msg = are_equivalent(B,Bp)
        if not eq:
            print("ERROR: B is not equivalent to Bp after random aligns")
            print(msg)
            print_pair(B,Bp)

        ## Random aligns (inplace)
        Bp = deepcopy(B)
        for i in range(NO_ALIGNS):
            Bp.align_to(np.random.permutation(B.vars), inplace=True)

        eq, msg = are_equivalent(B,Bp)
        if not eq:
            print("ERROR: B is not equivalent to Bp after random aligns (inplace)")
            print(msg)
            print_pair(B,Bp)

        ## Gsifts
        Bp = deepcopy(B)
        for i in range(NO_ALIGNS):
            A = BDD.BDD.random(N)
            A.rename_vars(dict(zip(Bp.vars,np.random.permutation(Bp.vars))))
            Bp.gsifts(A)

        eq, msg = are_equivalent(B,Bp)
        if not eq:
            print("ERROR: B is not equivalent to Bp after gsifts with a random BDD")
            print(msg)
            print_pair(B,Bp)

        t1 = time()
        if (t1-t0) // P_STEP > time_ctr:
            time_ctr += 1
            print("▮", end="")
            sys.stdout.flush()

        test_ctr += 1
    print(f"]\nFinished {test_ctr} full test passes in {t1-t0} seconds ({PROGRESS_STEPS} steps)")
