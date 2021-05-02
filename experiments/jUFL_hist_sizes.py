"""Benchmarks DD sizes for different solution methods (colored UFL)

An experiment concerning the Uncapacitated Facility Location with colors (colored UFL):
- compares diagram sizes for color, covering, and intersection DDs vs the number of
  variables in 'naive' MIPs.
- experiment runner file: `get_cUFL_sizes.sh`

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

    args = parser.parse_args()

    if args.header:
        print("instance,num_type,value,comment")
        print("-1,legend,orig_simpl,Simplified problem")
        print("-1,legend,orig_minAB,Best of A and B")
        print("-1,legend,orig_5random,Best of 5 random")
        print("-1,legend,orig_gsifts1p,Greedy BDD sifts")
        exit(0)

    for k in range(int(args.K)):
        S_1, S_2, f_1, f_2 = jUFL.generate_instance(int(args.n), p=float(args.prob))

        C_1, _ = cUFL.build_randomized_cover_DD(S_1, f_1)
        C_2, _ = cUFL.build_randomized_cover_DD(S_2, f_2)
        # C_1, _ = cUFL.build_cover_DD(S_1, f_1)
        # C_2, _ = cUFL.build_cover_DD(S_2, f_2)

        C_1.make_reduced()
        C_2.make_reduced()

        C1_to_C2 = C_1.align_to(C_2.vars, inplace=False)
        C2_to_C1 = C_2.align_to(C_1.vars, inplace=False)

        int_DD_C1_to_C2 = DD.intersect(C1_to_C2, C_2)
        int_DD_C1_to_C2.make_reduced()

        int_DD_C2_to_C1 = DD.intersect(C_1, C2_to_C1)
        int_DD_C2_to_C1.make_reduced()

        print(f"{k},orig_minAB_obj,{min(int_DD_C1_to_C2.size(), int_DD_C2_to_C1.size())},--none--")

        vs_1 = vs.VarSeq(C_1.vars, [len(L) for L in C_1.layers[:-1]])
        vs_2 = vs.VarSeq(C_2.vars, [len(L) for L in C_2.layers[:-1]])

        assert set(vs_1.layer_var) == set(vs_2.layer_var), f"A:{vs_1.layer_var}, C:{vs_2.layer_var}"
        b = bb.BBSearch(vs_1, vs_2)

        bb.TIMEOUT_ITERATIONS=15000
        status = b.search()
        assert status == "optimal"

        C1p = C_1.align_to(b.Ap_cand.layer_var, inplace=False)
        C2p = C_2.align_to(b.Ap_cand.layer_var, inplace=False)

        int_DD_VS = DD.intersect(C1p, C2p)
        int_DD_VS.make_reduced()

        print(f"{k},orig_simpl_obj,{int_DD_VS.size()},--none--")

        # int_DD_rnd = None
        # for _ in range(5):
        #     perm = np.random.permutation(C_1.vars)
        #     C1_rnd = C_1.align_to(perm, inplace=False)
        #     C2_rnd = C_2.align_to(perm, inplace=False)

        #     assert C1_rnd.is_reduced()
        #     assert C2_rnd.is_reduced()

        #     int_rnd = DD.intersect(C1_rnd, C2_rnd)
        #     int_rnd.make_reduced()

        #     if int_DD_rnd is None or int_DD_rnd > int_rnd.size():
        #         int_DD_rnd = int_rnd.size()

        # print(f"{k},orig_5random_obj,{int_DD_rnd},--none--")

        # C_1.gsifts(C_2)
        # int_DD_gsifts = DD.intersect(C_1, C_2)
        # int_DD_gsifts.make_reduced()
        # print(f"{k},orig_gsifts1p_obj,{int_DD_gsifts.size()},--none--")

        sys.stdout.flush()

if __name__ == '__main__':
    main()
