"""
Performs an experiment: summarizes objectives
derived from the simplified problem wrt
different starting orderings of diagrams.
"""

from varseq import VarSeq
from BDD import BDD
import numpy as np
from math import factorial
from itertools import permutations
from copy import deepcopy
from BB_search import BBSearch
import matplotlib.pyplot as plt

N = 8
# if __name__ == "__main__":
print("# Generating a problem...")
A = BDD.random(N=N); B = BDD.random(N=N)
A.rename_vars(dict(zip([i for i in range(1,N+1)],np.random.permutation([i for i in range(1,N+1)]))))
B.rename_vars(dict(zip([i for i in range(1,N+1)],np.random.permutation([i for i in range(1,N+1)]))))
A.make_reduced()
B.make_reduced()
print("# Generated. Initial size {}".format(A.size()+B.size()))
print("# A: {}|{}; B:{}|{}".format(A.vars, [A.n(i) for i in range(N)],
    B.vars, [B.n(i) for i in range(N)]))

print("# Solving...")
no_opts, As, Bs = A.OA_bruteforce(B)
print("# done. Optimal objective: {}".format(As[0].size()+Bs[0].size()))
print("# Building all the possible {} orderings...".format(factorial(N)))
perms = permutations(A.vars)
orders = [o for o in perms]
print("# Solving for all the possible {} problem configurations...".format(len(orders)*(len(orders)+1)/2))
print("MIN_I_OPT,I_AB,objective")
for i in range(len(orders)):
    for j in range(i,len(orders)):
        c = [orders[i],orders[j]]

        Ap = deepcopy(A); Bp = deepcopy(B)
        Ap.align_to(c[0],inplace=True); B.align_to(c[1],inplace=True)
        vsA = VarSeq(Ap.vars,[Ap.n(i) for i in range(N)])
        vsB = VarSeq(Bp.vars,[Bp.n(i) for i in range(N)])
        b = BBSearch(vsA,vsB)
        b.search()
        Ap.align_to(b.Ap_cand.layer_var, inplace=True); Bp.align_to(b.Ap_cand.layer_var, inplace=True)
        obj = Ap.size()+Bp.size()
        print("{},{},{}".format(min([vsA.count_inversions_to(D.vars) for D in As]+[vsB.count_inversions_to(D.vars) for D in As]),vsA.count_inversions_to(vsB),obj))

