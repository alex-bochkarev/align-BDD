"""Benchmarks DD sizes for different solution methods (colored UFL)

An experiment concerning the Uncapacitated Facility Location with types (typed UFL)

---
(c) A. Bochkarev, Clemson University, 2021
abochka@clemson.edu
"""
from time import time
import argparse as ap
import sys
from gurobipy import GRB

import jUFL
import cUFL
import varseq as vs
import BDD as DD
import BB_search as bb
import numpy as np


def main():
    """Main code (to be run from the command line)."""
    parser = ap.ArgumentParser(
        description= # noqa
        "Joint Uncapacitated Facility Location benchmarking (c) A. Bochkarev, 2021",
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
    parser.add_argument("-R",
                        "--random-order",
                        action="store_true",
                        dest="randomize",
                        default=False,
                        help="randomize the order of the two diagrams.")

    args = parser.parse_args()

    if args.header:
        print("instance,n,prob,num_type,value")
        exit(0)

    for k in range(int(args.K)):
        # S_1, S_2, f_1, f_2 = jUFL.generate_instance(int(args.n), p=float(args.prob))
        S, f, fc, kb = cUFL.generate_test_instance(int(args.n),
                                                   p=float(args.prob))

        # C_1, _ = cUFL.build_randomized_cover_DD(S_1, f_1)
        # C_2, _ = cUFL.build_randomized_cover_DD(S_2, f_2)
        cover_DD, _ = cUFL.build_cover_DD(S, f)
        pref_order = [int(x[1:]) for x in cover_DD.vars]
        type_DD, _ = cUFL.build_color_DD(f, fc, kb,
                                         preferred_order=pref_order)

        cover_DD.make_reduced()
        type_DD.make_reduced()

        if args.randomize:
            cover_DD.align_to(np.random.permutation(cover_DD.vars), inplace=True)
            type_DD.align_to(np.random.permutation(type_DD.vars), inplace=True)

        cover_DD_to_type = cover_DD.align_to(type_DD.vars, inplace=False)
        type_to_cover_DD = type_DD.align_to(cover_DD.vars, inplace=False)

        int_DD_cover_DD_to_type = DD.intersect(cover_DD_to_type, type_DD)
        int_DD_cover_DD_to_type.make_reduced()

        int_DD_type_to_cover_DD = DD.intersect(cover_DD, type_to_cover_DD)
        int_DD_type_to_cover_DD.make_reduced()

        print(f"{args.prefix}-{k},{args.n},{args.prob},orig_minAB_obj,{min(int_DD_cover_DD_to_type.size(), int_DD_type_to_cover_DD.size())}")

        vs_1 = vs.VarSeq(cover_DD.vars, [len(L) for L in cover_DD.layers[:-1]])
        vs_2 = vs.VarSeq(type_DD.vars, [len(L) for L in type_DD.layers[:-1]])

        assert set(vs_1.layer_var) == set(vs_2.layer_var), f"A:{vs_1.layer_var}, C:{vs_2.layer_var}"
        b = bb.BBSearch(vs_1, vs_2)

        bb.TIMEOUT_ITERATIONS=15000
        status = b.search()
        assert status == "optimal"

        coverp = cover_DD.align_to(b.Ap_cand.layer_var, inplace=False)
        typep = type_DD.align_to(b.Ap_cand.layer_var, inplace=False)

        int_DD_VS = DD.intersect(coverp, typep)
        int_DD_VS.make_reduced()

        print(f"{args.prefix}-{k},{args.n},{args.prob},orig_simpl_obj,{int_DD_VS.size()}")

        sys.stdout.flush()

if __name__ == '__main__':
    main()
