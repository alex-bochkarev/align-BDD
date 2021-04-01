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

    args = parser.parse_args()

    if args.header:
        # print("run,k,n,m,method,step,duration")
        print("n, color_size, cover_size, int_size, int_cov2col_size, int_col2cov_size, plain_MIP_vars, exp_time")
        exit(0)

    for _ in range(int(args.K)):
        t0 = time()
        S, f, fc, kb = cUFL.generate_test_instance(int(args.n))
        color, _ = cUFL.build_color_DD(f, fc, kb)
        cover, _ = cUFL.build_cover_DD(S, f)

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

        model = cUFL.build_cUFL_MIP(S, f, fc, kb)
        t1 = time()
        print(f"{args.n}, {color.size()}, {cover.size()}, {int_DD.size()}, {int_DDp.size()}, {int_DDpp.size()}, {len(model.getVars())}, {(t1-t0):.1f}")
        sys.stdout.flush()

if __name__ == '__main__':
    main()
