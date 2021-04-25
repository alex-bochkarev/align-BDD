"""Benchmarks DD sizes for different solution methods (colored UFL)

An experiment concerning the Uncapacitated Facility Location with colors (colored UFL):
- compares diagram sizes for color, covering, and intersection DDs vs the number of
  variables in 'naive' MIPs.
- experiment runner file: `get_cUFL_sizes.sh`

---
(c) A. Bochkarev, Clemson University, 2021
abochka@clemson.edu
"""
from copy import deepcopy
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
        "Uncapacitated Facility Location benchmarking (c) A. Bochkarev, 2021",
        formatter_class=ap.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-n",
                        "--facilities",
                        action="store",
                        dest="n",
                        default="3",
                        help="no. of facilities/customers")
    parser.add_argument("-K",
                        "--instances",
                        action="store",
                        dest="K",
                        default="50",
                        help="no. of instances to generate")
    parser.add_argument("-p",
                        "--prefix",
                        action="store",
                        dest="prefix",
                        default="0",
                        help="prefix for the instance number")
    parser.add_argument("-H",
                        "--header",
                        action="store_true",
                        dest="header",
                        help="show header only and exit")

    args = parser.parse_args()
    prefix = args.prefix

    if args.header:
        print("instance,num_type,value,comment")
        print("-1,legend,orig_simpl_nat,Simplified problem (natural var order)")
        print("-1,legend,orig_simpl_rnd,Simplified problem (random var order)")
        print("-1,legend,orig_minAB_nat,Best of A and B (natural var order)")
        print("-1,legend,orig_minAB_rnd,Best of A and B (random var order)")
        print("-1,legend,orig_5random_nat,Best of 5 random (natural var order)")
        print("-1,legend,orig_5random_rnd,Best of 5 random (random var order)")
        print("-1,legend,orig_5random_ctl,Best of 5 random (control experiment)")
        print("-1,legend,orig_minAB_ctl,Best of A and B (control experiment)")
        print("-1,legend,orig_simpl_ctl,Simplified problem (control experiment)")
        exit(0)

    for k in range(int(args.K)):
        S, f, fc, kb = cUFL.generate_test_instance(int(args.n))

        cover_nat, _ = cUFL.build_cover_DD(S, f)
        pref_order = [int(x[1:]) for x in cover_nat.vars]
        color_nat, _ = cUFL.build_color_DD(f, fc, kb, pref_order)

        color_nat.make_reduced()
        cover_nat.make_reduced()

        cover_rnd, _ = cUFL.build_randomized_cover_DD(S, f)
        color_rnd, _ = cUFL.build_randomized_color_DD(f, fc, kb)

        color_rnd.make_reduced()
        cover_rnd.make_reduced()

        A = DD.BDD.random(N = int(args.n), p=0.6, weighted=True)
        B = DD.BDD.random(N = int(args.n), p=0.6, weighted=True)
        B.rename_vars(dict(zip([i for i in range(1,int(args.n)+1)],
                               np.random.permutation([i for i in range(1,int(args.n)+1)]))))

        A.make_reduced()
        B.make_reduced()

        experiment_versions = [(color_nat, cover_nat, "nat"),
                               (color_rnd, cover_rnd, "rnd"),
                               (A, B, "ctl")]

        for experiment in experiment_versions:
            color, cover, dsc = experiment

            print(f"{prefix}-{k},color_size_{dsc},{color.size()},--none--")
            print(f"{prefix}-{k},cover_size_{dsc},{cover.size()},--none--")
            vs_color = vs.VarSeq(color.vars, [len(L) for L in color.layers[:-1]])
            vs_cover = vs.VarSeq(cover.vars, [len(L) for L in cover.layers[:-1]])

            bb.TIMEOUT_ITERATIONS = 10000
            b = bb.BBSearch(vs_color, vs_cover)
            status = b.search()
            assert status in ["optimal", "timeout"]

            color_p = color.align_to(b.Ap_cand.layer_var, inplace=False)
            cover_p = cover.align_to(b.Ap_cand.layer_var, inplace=False)

            int_DD_VS = DD.intersect(color_p, cover_p)
            int_DD_VS.make_reduced()

            print(f"{prefix}-{k},orig_simpl_{dsc}_obj,{int_DD_VS.size()},--none--")

            cover_to_color = cover.align_to(color.vars, inplace=False)
            color_to_cover = color.align_to(cover.vars, inplace=False)

            int_DD_cov2col = DD.intersect(color, cover_to_color)
            int_DD_cov2col.make_reduced()

            int_DD_col2cov = DD.intersect(color_to_cover, cover)
            int_DD_col2cov.make_reduced()

            print(f"{prefix}-{k},orig_minAB_{dsc}_obj,{min(int_DD_cov2col.size(), int_DD_col2cov.size())},--none--")

            # int_DD_rnd = None
            # for _ in range(5):
            #     perm = np.random.permutation(color.vars)
            #     cov_rnd = cover.align_to(perm, inplace=False)
            #     col_rnd = color.align_to(perm, inplace=False)

            #     int_rnd = DD.intersect(cov_rnd, col_rnd)
            #     int_rnd.make_reduced()
            #     if int_DD_rnd is None or int_DD_rnd > int_rnd.size():
            #         int_DD_rnd = int_rnd.size()

            # print(f"{prefix}-{k},orig_5random_{dsc}_obj,{int_DD_rnd},--none--")

            sys.stdout.flush()

if __name__ == '__main__':
    main()
