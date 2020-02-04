"""
"Scalability" test:
Aux script: time-tracks some instances of different size:
- generate (unique) original and simplified problems;
- solve simplified problem: heuristics, BB, brute-force;
- solve exact problem: heuristics from the simplified prob, sifting;

(c) A. Bochkarev, Clemson University, 2020
abochka@clemson.edu
"""

import varseq as simpl
import BDD as exact
import numpy as np
import BB_search as bbs
from experiments.misc import log
import sys
from time import time
from copy import deepcopy
import gc as garbage
import argparse as ap

start_id = 0
Ns_short = [5,6,7,8] # 9,10,11,12,13,14,15]

if __name__ == "__main__":

    # create an exact instance
    parser = ap.ArgumentParser(description="Performs scalability test. (c) A. Bochkarev, Clemson University, 2020",
                               formatter_class=ap.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-v","--no_variables",action="store", dest="V",
                        help="number of variables per instance (0 for the 'short' set, see the source for details)",
                        type=int, default=0)
    parser.add_argument("-k","--no_insts",action="store", dest="K",
                        help="number of instances per instance size",
                        type=int, default=10)
    parser.add_argument("-p", "--prob_exp", dest="p",help="tree expansion probability parameter",
                        action="store", type=float, default=0.6)
    parser.add_argument("-H","--header", help="print column headers only and exit",action="store_true")
    parser.add_argument("-t","--tmpdir", dest="out_dir", help="temporary output directory")
    args = parser.parse_args()

    if args.header:
        log("instance","N","num_type","value",comment="comment")
        exit(0)

    inst_profiles = set()
    inst_id = 0
    if args.V <= 0:
        Ns = Ns_short
    else:
        Ns = [args.V]

    K = args.K
    p = args.p

for N in Ns:
    for k in range(K):
        inst_accepted = False
        trials = 0
        t0 = time()
        while not inst_accepted:
            # this loop is needed to ensure instances are unique
            trials += 1
            A = exact.BDD.random(N=N,p=p)
            B = exact.BDD.random(N=N,p=p)
            B.rename_vars(dict(zip([i for i in range(1,N+1)],np.random.permutation([i for i in range(1,N+1)]))))
            # A.make_reduced()
            # B.make_reduced()

            # check if the instance is unique / never seen before
            prof_A = A.profile()
            prof_B = B.profile()
            if not ( ( (prof_A,prof_B) in inst_profiles ) or ( (prof_B, prof_A) in inst_profiles )):
                inst_profiles.add((prof_A,prof_B))
                inst_accepted = True

        A.save(args.out_dir+"A{}.bdd".format(inst_id))
        B.save(args.out_dir+"B{}.bdd".format(inst_id))

        t1 = time()
        log(inst_id, N, "trials",trials)
        log(inst_id, N,"orig_gen_time",t1 - t0)

        # generate corresponding varseq-instances
        t0 = time()
        vsA = simpl.VarSeq(A.vars, [len(l) for l in A.layers[:-1]])
        vsB = simpl.VarSeq(B.vars, [len(l) for l in B.layers[:-1]])
        t1 = time()
        log(inst_id, N, "vs_gen_time",t1 - t0)

        t0 = time()
        bb = bbs.BBSearch(vsA,vsB)
        bb.search()
        t1 = time()
        log(inst_id, N, "BBsearch_time", t1-t0)

        t0 = time()
        order = bb.Ap_cand.layer_var
        simpl_obj = A.align_to(order).size() + B.align_to(order).size()
        t1 = time()
        log(inst_id, N, "alig_to_order_time", t1-t0)
        log(inst_id, N, "orig_aux_obj", simpl_obj)

        ## solve with exact sifts (with and without preliminary reduction)
        t0 = time()
        if A.align_to(B.vars).size() + B.size() <= A.size()+B.align_to(A.vars).size():
            A.gsifts(B)
        else:
            B.gsifts(A)

        t1 = time()

        log(inst_id, N, "orig_sifts_red_time", t1-t0)
        log(inst_id, N, "orig_sifts_red_obj", A.size()+B.size())
        inst_id += 1
