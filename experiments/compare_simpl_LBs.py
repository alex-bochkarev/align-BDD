"""Compares different lower bounds for the simplified problem.

Given a list of (align-BDD) problems, builds and solves a simplified
problem corresponding to each one to calculate the 'LB tightness'
parameter and keep track of the runtime.
"""
import varseq as simpl
import BDD as exact
import BB_search as bbs
import argparse as ap
from experiments.misc import log
from time import time

if __name__ == "__main__":
    parser = ap.ArgumentParser(description="Compares various LB performance for the simplified problem. (c) A. Bochkarev, Clemson University, 2020",
                               formatter_class=ap.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-d","--dir", dest="inst_dir", help="temporary output directory")
    parser.add_argument("-l","--list_inst", dest="inst_list", help="list of instance IDs (filename)")
    parser.add_argument("-H", "--header", action="store_true", dest="header", help="show header only and exit")
    args = parser.parse_args()

    if args.header:
        log("LB,gap",comment="legend")
        exit(0)

    with open(args.inst_list,"r") as inst_list:
        for inst_id in inst_list:
            inst_id = int(inst_id.rstrip())
            A = exact.BDD(); B = exact.BDD()
            A.load(args.inst_dir+"/A{}.bdd".format(inst_id))
            B.load(args.inst_dir+"/B{}.bdd".format(inst_id))

            N = len(A.vars)

            vsA = simpl.VarSeq(A.vars, [len(l) for l in A.layers[:-1]])
            vsB = simpl.VarSeq(B.vars, [len(l) for l in B.layers[:-1]])
            t0 = time()
            bb = bbs.BBSearch(vsA,vsB)
            bb.search()
            t1 = time()
            log("simpl_BB",t1-t0,comment="timelog")

            curr_size = vsA.size()+vsB.size()
            opt_size = bb.Ap_cand.size()+bb.Bp_cand.size()
            if opt_size == curr_size:
                continue # not an interesting instance

            for LB in bbs.LOWER_BOUNDS:
                t0 = time()
                log(LB[0],(LB[1](vsA,vsB) - curr_size)/(opt_size - curr_size), comment=LB[2])
                t1 = time()
                log(LB[0],t1-t0,comment="timelog")
