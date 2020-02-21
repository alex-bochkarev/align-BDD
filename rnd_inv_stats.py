"""
Generates inversions statistics.

(c) A. Bochkarev, Clemson University, 2020
abochka@clemson.edu
"""
from BDD import BDD
from varseq import VarSeq
import matplotlib.pyplot as plt
import numpy as np
from time import time

N = 8 # no. of vars
n = 50 # no. of instances for the histogram

MIN_I_OPT = []
i = 0
print("INST_NO,MIN_I_OPT,PROCESSING_TIME")
for i in range(n):
    t0 = time()
    A = BDD.random(N=N); B = BDD.random(N=N)
    A.rename_vars(dict(zip([i for i in range(1,N+1)],np.random.permutation([i for i in range(1,N+1)]))))
    B.rename_vars(dict(zip([i for i in range(1,N+1)],np.random.permutation([i for i in range(1,N+1)]))))
    A.make_reduced()
    B.make_reduced()
    no_opts, As, Bs = A.OA_bruteforce(B)
    vsA = VarSeq(A.vars,[A.n(i) for i in range(N)])
    vsB = VarSeq(B.vars,[B.n(i) for i in range(N)])
    MIN_I_OPT = min([vsA.count_inversions_to(D.vars) for D in As]+[vsB.count_inversions_to(D.vars) for D in As])
    i += 1
    t1 = time()
    print("{},{},{}".format(i,MIN_I_OPT,t1-t0),flush=True)
