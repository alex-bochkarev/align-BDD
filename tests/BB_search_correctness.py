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

import varseq as vs
import BB_search as bb

STEPS_CHUNK = 1
print("Running a simple BB search correctness test (1 period = {} instances)".format(STEPS_CHUNK))
print("[Ctrl-C] to stop")

ctr = 0
while True:
    A = vs.VarSeq.random(N = 8)
    B = vs.VarSeq.random(N = 8)

    b = bb.BBSearch(A,B)

    alts, Ap, Bp = A.OA_bruteforce(B)

    opt_size = Ap[0].size() + Bp[0].size()

    status = b.search()

    BB_size = b.Ap_cand.size() + b.Bp_cand.size()

    ctr = (ctr + 1) % STEPS_CHUNK
    if BB_size != opt_size:
        print("MISMATCH DETECTED:")
        print("Instance:\n{}\nvs.\n{}\nBB size={}, opt size ={}".format(A,B,BB_size,opt_size), flush = True)
    elif ctr == 0:
        print(".",end="",flush=True)
