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

def test_BDD_transformations(B,no_iters=250):
    """Implements a simple ad-hoc test"""
    # let's test the basic things
    print("First: for Bp = B,")
    Bp = deepcopy(B)
    print("Variables:\n B:{}\nBp:{}".format(B.vars,Bp.vars))
    eq, _ = B.is_equivalent(Bp)
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
    eq, _ = B.is_equivalent(Bp)
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
    eq, _ = B.is_equivalent(Bp)
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
    eq, _ = B.is_equivalent(Bp)
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
        eq, msg = B.is_equivalent(Bp)
        if not eq:
            print("ERROR: B not equivalent to Bp=B!")
            print(msg)
            print_pair(B,Bp)

        # non-eq to another random BDD
        keep_trying = True
        Bp = BDD.BDD.random(N)

        while keep_trying:
            eq, msg = B.is_equivalent(Bp)
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

        eq, msg = B.is_equivalent(Bp)
        if not eq:
            print("ERROR: B is not equivalent to Bp after random swaps")
            print(msg)
            print_pair(B,Bp)

        ## Random sifts
        Bp = deepcopy(B)
        for i in range(NO_SIFTS):
            Bp.sift(B.vars[np.random.randint(N)], np.random.randint(N))

        eq, msg = B.is_equivalent(Bp)
        if not eq:
            print("ERROR: B is not equivalent to Bp after random sifts")
            print(msg)
            print_pair(B,Bp)

        ## Random aligns (not inplace)
        Bp = deepcopy(B)
        for i in range(NO_ALIGNS):
            Bp = Bp.align_to(np.random.permutation(B.vars), inplace=False)

        eq, msg = B.is_equivalent(Bp)
        if not eq:
            print("ERROR: B is not equivalent to Bp after random aligns")
            print(msg)
            print_pair(B,Bp)

        ## Random aligns (inplace)
        Bp = deepcopy(B)
        for i in range(NO_ALIGNS):
            Bp.align_to(np.random.permutation(B.vars), inplace=True)

        eq, msg = B.is_equivalent(Bp)
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

        eq, msg = B.is_equivalent(Bp)
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
