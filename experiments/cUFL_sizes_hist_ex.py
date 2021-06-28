"""Benchmarks DD sizes for different solution methods (colored UFL)

An experiment concerning the Uncapacitated Facility Location with colors (colored UFL):

    - compares diagram sizes for color, covering, and intersection DDs
        vs the number of variables in 'naive' MIPs.
    - experiment runner file: `get_cUFL_sizes.sh`

(c) A. Bochkarev, Clemson University, 2021, abochka@clemson.edu
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
    parser.add_argument("-H",
                        "--header",
                        action="store_true",
                        dest="header",
                        help="show header only and exit")

    args = parser.parse_args()

    if args.header:
        print("instance,num_type,value,comment")
        print("-1,legend,orig_simpl,Simplified problem: |A∩B|")
        print("-1,legend,orig_minAB,Best of A and B: |A∩B|")
        print("-1,legend,orig_5random,Best of 5 random: |A∩B|")
        print("-1,legend,orig_gsifts1p,Greedy BDD sifts: |A∩B|")
        print("-1,legend,orig_simpl_AB,Simplified problem: |A|+|B|")
        print("-1,legend,orig_minAB_AB,Best of A and B: |A|+|B|")
        print("-1,legend,orig_5random_AB,Best of 5 random: |A|+|B|")
        print("-1,legend,orig_gsifts1p_AB,Greedy BDD sifts: |A|+|B|")
        exit(0)

    for k in range(int(args.K)):
        S, f, fc, kb = cUFL.generate_test_instance(int(args.n))

        cover, _ = cUFL.build_cover_DD(S, f)
        pref_order = [int(x[1:]) for x in cover.vars]
        color, _ = cUFL.build_color_DD(f, fc, kb, pref_order)

        color.make_reduced()
        cover.make_reduced()

        vs_color = vs.VarSeq(color.vars, [len(L) for L in color.layers[:-1]])
        vs_cover = vs.VarSeq(cover.vars, [len(L) for L in cover.layers[:-1]])

        bb.TIMEOUT_ITERATIONS = 10000
        b = bb.BBSearch(vs_color, vs_cover)
        status = b.search()
        assert status == "optimal" or status == "timeout"

        color_p = color.align_to(b.Ap_cand.layer_var, inplace=False)
        cover_p = cover.align_to(b.Ap_cand.layer_var, inplace=False)

        # assert color_p.is_reduced()
        # assert cover_p.is_reduced()

        print(f"{k},orig_simpl_AB_obj,{color_p.size() + cover_p.size()},--none--")

        int_DD_VS = DD.intersect(color_p, cover_p)
        int_DD_VS.make_reduced()

        print(f"{k},orig_simpl_obj,{int_DD_VS.size()},--none--")

        cover_to_color = cover.align_to(color.vars, inplace=False)
        color_to_cover = color.align_to(cover.vars, inplace=False)

        # assert cover_to_color.is_reduced()
        # assert color_to_cover.is_reduced()

        minAB = min(color.size()+cover_to_color.size(),
                    cover.size()+color_to_cover.size())

        print(f"{k},orig_minAB_AB_obj,{minAB},--none--")

        int_DD_cov2col = DD.intersect(color, cover_to_color)
        int_DD_cov2col.make_reduced()

        int_DD_col2cov = DD.intersect(color_to_cover, cover)
        int_DD_col2cov.make_reduced()

        print(f"{k},orig_minAB_obj,{min(int_DD_cov2col.size(), int_DD_col2cov.size())},--none--")
        int_DD_rnd = None
        rnd_5_size = None

        for _ in range(5):
            perm = np.random.permutation(color.vars)
            cov_rnd = cover.align_to(perm, inplace=False)
            col_rnd = color.align_to(perm, inplace=False)

            # assert cov_rnd.is_reduced()
            # assert col_rnd.is_reduced()

            cur_rnd_size = cov_rnd.size() + col_rnd.size()

            if rnd_5_size is None or rnd_5_size > cur_rnd_size:
                rnd_5_size = cur_rnd_size

            int_rnd = DD.intersect(cov_rnd, col_rnd)
            int_rnd.make_reduced()
            if int_DD_rnd is None or int_DD_rnd > int_rnd.size():
                int_DD_rnd = int_rnd.size()

            if rnd_5_size is None or rnd_5_size > cur_rnd_size:
                rnd_5_size = cur_rnd_size

        print(f"{k},orig_5random_obj,{int_DD_rnd},--none--")
        print(f"{k},orig_5random_AB_obj,{rnd_5_size},--none--")

        cover.gsifts(color)
        # assert cover.is_reduced()
        # assert color.is_reduced()
        print(f"{k},orig_gsifts1p_AB_obj,{cover.size()+color.size()},--none--")

        int_DD_gsifts = DD.intersect(cover, color)
        int_DD_gsifts.make_reduced()
        print(f"{k},orig_gsifts1p_obj,{int_DD_gsifts.size()},--none--")

        sys.stdout.flush()

if __name__ == '__main__':
    main()
