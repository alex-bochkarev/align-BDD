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

import cUFL
import varseq as vs
import BDD as DD
import BB_search as bb
from copy import deepcopy

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
    parser.add_argument("-p",
                        "--prefix",
                        action="store",
                        dest="prefix",
                        default="0",
                        help="prefix string (ID of the run)")
    parser.add_argument("-t",
                        "--type",
                        action="store",
                        dest="inst_type",
                        default="ER",
                        help="instance type (ER or STRING)")

    args = parser.parse_args()

    if args.header:
        # print("run,k,n,m,method,step,duration")
        print("n, instance, PREORDER, color_size, cover_size, int_size, int_cov2col_size, int_col2cov_size, int_cov2col_gsifts, int_col2cov_gsifts, plain_MIP_vars, plain_MIP_constrs, exp_time")
        exit(0)

    for instance in range(int(args.K)):
        t0 = time()
        if args.inst_type == "ER":
            S, f, fc, kb = cUFL.generate_test_instance(int(args.n))
        elif args.inst_type == "STRING":
            S, f, fc, kb = cUFL.generate_string_instance(int(args.n))
        else:
            print(f"FATAL ERROR: {args.inst_type} -- wrong instance type. Expected `ER` or `STRING`.")
            exit(-1)

        for preorder in [True, False]:
            cover, _ = cUFL.build_cover_DD(S, f)
            if preorder:
                pref_order = [int(x[1:]) for x in cover.vars]
            else:
                pref_order = None

            color, _ = cUFL.build_color_DD(f, fc, kb, pref_order)

            vs_color = vs.VarSeq(color.vars, [len(L) for L in color.layers[:-1]])
            vs_cover = vs.VarSeq(cover.vars, [len(L) for L in cover.layers[:-1]])

            bb.TIMEOUT_ITERATIONS = 5000
            b = bb.BBSearch(vs_color, vs_cover)
            status = b.search()
            assert status == "optimal" or status == "timeout"

            color_p = color.align_to(b.Ap_cand.layer_var, inplace=False)
            cover_p = cover.align_to(b.Ap_cand.layer_var, inplace=False)

            cover_to_color = cover.align_to(color.vars, inplace=False)
            color_to_cover = color.align_to(cover.vars, inplace=False)

            int_DD = DD.intersect(color_p, cover_p)

            int_DDp = DD.intersect(color, cover_to_color)
            int_DDpp = DD.intersect(color_to_cover, cover)

            cov_c = deepcopy(cover)
            col_c = deepcopy(color)

            cov_c.gsifts(col_c)
            int_cov2col_gsifts = DD.intersect(cov_c, col_c)

            cover.gsifts(color)
            int_col2cov_gsifts = DD.intersect(cover, color)

            model = cUFL.build_cUFL_MIP(S, f, fc, kb)
            t1 = time()
            print(f"{args.n}, {args.prefix}-{instance}, {preorder}, {color.size()}, {cover.size()}, {int_DD.size()}, {int_DDp.size()}, {int_DDpp.size()}, {int_cov2col_gsifts.size()}, {int_col2cov_gsifts.size()}, {len(model.getVars())}, {len(model.getConstrs())}, {(t1-t0):.1f}")
            sys.stdout.flush()

if __name__ == '__main__':
    main()
