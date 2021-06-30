"""Benchmarks DD sizes for different solution methods (typed UFL)

An experiment concerning the Uncapacitated Facility Location with types (types UFL):
- compares intersection diagram sizes in several settings
- experiment runner file: `tUFL_bm_static.sh`

---
(c) A. Bochkarev, Clemson University, 2021
abochka@clemson.edu
"""
from copy import deepcopy
from time import time
import argparse as ap
import sys
from gurobipy import GRB

import tUFLP
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

        print("-1,legend,orig_gsifts1p_nat,Greedy BDD sifts (natural var order)")
        print("-1,legend,orig_gsifts1p_rnd,Greedy BDD sifts (random var order)")
        exit(0)

    for k in range(int(args.K)):
        S, f, fc, kb = tUFLP.generate_test_instance(int(args.n))

        cover_nat, _ = tUFLP.build_cover_DD(S, f)
        pref_order = [int(x[1:]) for x in cover_nat.vars]
        type_nat, _ = tUFLP.build_type_DD(f, fc, kb, pref_order)

        type_nat.make_reduced()
        cover_nat.make_reduced()

        cover_rnd = cover_nat.align_to(np.random.permutation(cover_nat.vars),
                                          inplace=False)
        type_rnd = type_nat.align_to(np.random.permutation(type_nat.vars),
                                        inplace=False)

        experiment_versions = [(type_nat, cover_nat, "nat"),
                               (type_rnd, cover_rnd, "rnd")]

        for experiment in experiment_versions:
            type_DD, cover_DD, dsc = experiment

            # create and solve the simplified problem
            vs_type = vs.VarSeq(type_DD.vars, [len(L) for L in type_DD.layers[:-1]])
            vs_cover = vs.VarSeq(cover_DD.vars, [len(L) for L in cover_DD.layers[:-1]])

            bb.TIMEOUT_ITERATIONS = 10000
            b = bb.BBSearch(vs_type, vs_cover)
            status = b.search()
            assert status in ["optimal", "timeout"]

            type_p = type_DD.align_to(b.Ap_cand.layer_var, inplace=False)
            cover_p = cover_DD.align_to(b.Ap_cand.layer_var, inplace=False)

            int_DD_VS = DD.intersect(type_p, cover_p)
            int_DD_VS.make_reduced()

            print(f"{prefix}-{k},orig_simpl_{dsc}_obj,{int_DD_VS.size()},--none--")

            # baseline 1: A-to-B and B-to-A
            cover_to_type = cover_DD.align_to(type_DD.vars, inplace=False)
            type_to_cover = type_DD.align_to(cover_DD.vars, inplace=False)

            int_DD_cov2type = DD.intersect(type_DD, cover_to_type)
            int_DD_cov2type.make_reduced()

            int_DD_type2cov = DD.intersect(type_to_cover, cover_DD)
            int_DD_type2cov.make_reduced()

            print(
                f"{prefix}-{k},orig_cov2type_{dsc}_obj,{int_DD_cov2type.size()},--none--"
            )
            print(
                f"{prefix}-{k},orig_type2cov_{dsc}_obj,{int_DD_type2cov.size()},--none--"
            )
            print(
                f"{prefix}-{k},orig_minAB_{dsc}_obj,{min(int_DD_cov2type.size(), int_DD_type2cov.size())},--none--"
            )

            int_size_rnd = None
            for _ in range(5):
                perm = np.random.permutation(type_DD.vars)
                cov_rnd = cover_DD.align_to(perm, inplace=False)
                type_rnd = type_DD.align_to(perm, inplace=False)

                int_rnd = DD.intersect(cov_rnd, type_rnd)
                int_rnd.make_reduced()

                if int_size_rnd is None or int_size_rnd > int_rnd.size():
                    int_size_rnd = int_rnd.size()

            print(f"{prefix}-{k},orig_5random_{dsc}_obj,{int_size_rnd},--none--")

            cover2 = deepcopy(cover_DD)
            type2 = deepcopy(type_DD)

            cover_DD.gsifts(type_DD, start_order=type_DD.vars)
            cover2.gsifts(type2, start_order=cover2.vars)

            int_DD_gsifts = DD.intersect(cover_DD, type_DD)
            int_DD_gsifts.make_reduced()

            int_DD_gsifts2 = DD.intersect(cover2, type2)
            int_DD_gsifts2.make_reduced()

            print(f"{prefix}-{k},orig_gsifts1p_cov2type_{dsc}_obj,{int_DD_gsifts.size()},--none--")
            print(f"{prefix}-{k},orig_gsifts1p_type2cov_{dsc}_obj,{int_DD_gsifts2.size()},--none--")
            print(f"{prefix}-{k},orig_gsifts1p_{dsc}_obj,{min(int_DD_gsifts.size(), int_DD_gsifts2.size())},--none--")

            sys.stdout.flush()

if __name__ == '__main__':
    main()
