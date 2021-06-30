"""
"Scalability" test:
Aux script: time-tracks some instances of different size:
- open an instance;
- solve the simplified problem: heuristics, BB;
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
    parser = ap.ArgumentParser(description="Performs scalability test. (c) A. Bochkarev, Clemson University, 2020",
                               formatter_class=ap.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-H","--header", help="print column headers only and exit",action="store_true")
    parser.add_argument("-d","--dir", dest="inst_dir", help="temporary output directory")
    parser.add_argument("-l","--list_inst", dest="inst_list", help="list of instance IDs (filename)")

    args = parser.parse_args()

    if args.header:
        log("instance","N","num_type","value",comment="comment")
        exit(0)

    inst_profiles = set()
    inst_id = 0

with open(args.inst_list,"r") as inst_list:
    for inst_id in inst_list:
        inst_id = int(inst_id.rstrip())
        A = exact.BDD(); B = exact.BDD()
        A.load(args.inst_dir+"/A{}.bdd".format(inst_id))
        B.load(args.inst_dir+"/B{}.bdd".format(inst_id))

        N = len(A.vars)

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
