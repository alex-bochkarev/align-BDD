"""
Test unit: BB search correctness.

Compares the objective obtained by the BB search
vs. the brute-force enumerated, true-optimal objective

An infinite loop, prints unsuccessful instances only.
Randomly generated instances (not necessarily unique)

(c) A. Bochkarev, Clemson University, 2020
abochka@clemson.edu
"""

import sys
sys.path.append('..')

import numpy as np
import varseq as vs
import BB_search as bb

print("Running a simple BB search correctness test (1 period = 1 instance)")
print("[Ctrl-C] to stop")

while True:
    A = vs.VarSeq.random(N = 20)
    B = vs.VarSeq.random(N = 20)

    Ap1 = A.greedy_sort(B.layer_var)
    Ap2 = A.q_align_to(B.layer_var)

    if not (np.array_equal(Ap1.layer_var,Ap2.layer_var) and \
            np.array_equal(Ap1.n, Ap2.n)):
        print("MISMATCH DETECTED:")
        print("Instance:\n{}\nvs.\n{}".format(A,B), flush = True)
    # else:
    #     print(".",end="",flush=True)
