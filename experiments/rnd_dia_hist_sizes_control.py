"""Benchmarks DD sizes for different solution methods (typed UFL)


---
(c) A. Bochkarev, Clemson University, 2021
abochka@clemson.edu
"""
from time import time
import argparse as ap
import sys
from gurobipy import GRB

import cUFL
import varseq as vs
import BDD as DD
import BB_search as bb
import numpy as np


def main():
    """Main code (to be run from the command line)."""
    parser = ap.ArgumentParser(
        description= # noqa
        "Random diagrams -- a control for t-UFLP (c) A. Bochkarev, 2021",
        formatter_class=ap.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-n",
                        "--facilities",
                        action="store",
                        dest="n",
                        default="3",
                        help="no. of facilities/customers")
    parser.add_argument("-p",
                        "--probs",
                        action="store",
                        dest="prob",
                        default="0.25",
                        help="prob parameter for instance generation")
    parser.add_argument("-P",
                        "--prefix",
                        action="store",
                        dest="prefix",
                        default="test",
                        help="prefix to be printed")
    parser.add_argument("-K",
                        "--instances",
                        action="store",
                        dest="K",
                        default="50",
                        help="no. of instances to generate")
    parser.add_argument("-H",
                        "--header",
                        action="store_true",
                        dest="header",
                        help="show header only and exit")
    parser.add_argument("-l",
                        "--instance-log",
                        action="store",
                        dest="logdir",
                        default="none",
                        help="dir name to save the instances (bdd format)")

    args = parser.parse_args()

    if args.header:
        print("instance,n,prob,num_type,value")
        exit(0)

    for k in range(int(args.K)):
        A = DD.BDD.random(N=int(args.n), p=float(args.prob), weighted=True)
        B = DD.BDD.random(N=int(args.n), p=float(args.prob), weighted=True)
        B.shuffle_vars()

        A.make_reduced()
        B.make_reduced()

        if args.logdir != "none":
            A.save(f"{args.logdir}/A{args.prefix}-{k}.bdd")
            B.save(f"{args.logdir}/B{args.prefix}-{k}.bdd")

        A2B = A.align_to(B.vars, inplace=False)
        B2A = B.align_to(A.vars, inplace=False)

        int_A2B = DD.intersect(A2B, B)
        int_A2B.make_reduced()

        int_B2A = DD.intersect(A, B2A)
        int_B2A.make_reduced()

        print(f"{args.prefix}-{k},{args.n},{args.prob},orig_minAB_obj,{min(int_A2B.size(), int_B2A.size())}")

        vsA = vs.VarSeq(A.vars, [len(L) for L in A.layers[:-1]])
        vsB = vs.VarSeq(B.vars, [len(L) for L in B.layers[:-1]])

        assert set(vsA.layer_var) == set(vsB.layer_var), f"A:{vsA.layer_var}, B:{vsB.layer_var}"
        b = bb.BBSearch(vsA, vsB)

        bb.TIMEOUT_ITERATIONS=15000
        status = b.search()
        assert status == "optimal"

        Ap = A.align_to(b.Ap_cand.layer_var, inplace=False)
        Bp = B.align_to(b.Ap_cand.layer_var, inplace=False)

        int_DD_VS = DD.intersect(Ap, Bp)
        int_DD_VS.make_reduced()

        print(f"{args.prefix}-{k},{args.n},{args.prob},orig_simpl_obj,{int_DD_VS.size()}")

        sys.stdout.flush()

if __name__ == '__main__':
    main()
