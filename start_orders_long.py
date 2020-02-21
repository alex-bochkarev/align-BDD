"""
Performs an experiment: summarizes objectives
derived from the simplified problem wrt
different starting orderings of diagrams.

(c) A. Bochkarev, Clemson University, 2020
abochka@clemson.edu
"""

from varseq import VarSeq
from BDD import BDD
import numpy as np
from math import factorial
from itertools import permutations
from copy import deepcopy,copy
from BB_search import BBSearch
from argparse import ArgumentsParser, ArgumentDefaultsHelpFormatter

N = 3
if __name__ == "__main__":
    parser = ArgumentParser(description="Iterates through all permutations for a random align-BDD instance (c) A. Bochkarev, Clemson University, 2020",
                               formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-N", "--no_vars", action="store", dest="N", help="help",
                        type=int,default=3)
    parser.add_argument("-d","--directory","out_dir", help="output directory (to save bdds)")

    args = parser.parse_args()
    N = args.N
    out_dir = args.out_dir
    A = BDD.random(N=N); B = BDD.random(N=N)
    A.rename_vars(dict(zip([i for i in range(1,N+1)],np.random.permutation([i for i in range(1,N+1)]))))
    B.rename_vars(dict(zip([i for i in range(1,N+1)],np.random.permutation([i for i in range(1,N+1)]))))
    A.make_reduced()
    B.make_reduced()
    print("# Generated. Initial size {}".format(A.size()+B.size()))
    A.save(f"{out_dir}/random_A.bdd")
    B.save(f"{out_dir}/random_B.bdd")

    print("varsA,varsB,aligned,size,orig_simpl_objective")
    A_vars = copy(A.vars)
    B_vars = copy(B.vars)

    iters_A = permutations(A_vars)
    iters_B = deepcopy(iters_A)
    for pA in iters_A:
        for pB in iters_B:
            A.align_to(pA,inplace=True); B.align_to(pB,inplace=True)
            vsA = VarSeq(A.vars,[A.n(i) for i in range(N)])
            vsB = VarSeq(B.vars,[B.n(i) for i in range(N)])
            b = BBSearch(vsA,vsB)
            b.search()
            aligned = A.is_aligned(B)
            size = A.size() + B.size()
            A.align_to(b.Ap_cand.layer_var, inplace=True); B.align_to(b.Ap_cand.layer_var, inplace=True)
            orig_simpl_objective = A.size() + B.size()
            vA = "-".join([str(x) for x in pA])
            vB = "-".join([str(x) for x in pB])
            print(f"{vA},{vB},{aligned},{size},{orig_simpl_objective}",flush=True)
        iters_B = deepcopy(iters_A)
