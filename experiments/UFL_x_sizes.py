"""Benchmarks DD sizes for different solution methods (UFL): test x-only DDs.

An experiment concerning the Uncapacitated Facility Location (UFL):
- compares diagram sizes for availability, covering, and DP-inspired DDs vs the number of
  variables in 'naive' MIPs.
- experiment runner file: `get_sizes.sh`

---
(c) A. Bochkarev, Clemson University, 2021
abochka@clemson.edu
"""
from time import time
import argparse as ap
import sys
from gurobipy import GRB

import UFL
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
                        help="no. of facilities")
    parser.add_argument("-m",
                        "--customers",
                        action="store",
                        dest="m",
                        default="4",
                        help="no. of customers")
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
        print("n, A_size, C_size, DP_size, plain_MIP_vars, DP_NF_vars, exp_time")
        exit(0)

    for _ in range(int(args.K)):
        t0 = time()
        # S, f, g = UFL.generate_test_instance(n=int(args.n), m=2*int(args.n))
        S, f, g = UFL.generate_dense_instance(n=int(args.n), m=2*int(args.n), covering=0.95)
        A = UFL.create_availability_BDD(S, f)
        C = UFL.create_covering_BDD_wg(S, g)

        DP_DD, _ = UFL.build_DP_DD(S, f, g)
        DP_model, _, _ = UFL.create_NF(DP_DD)

        model = UFL.build_MIP(S, f, g)
        t1 = time()
        print(f"{args.n}, {A.size()}, {C.size()}, {DP_DD.size()}, {len(model.getVars())}, {len(DP_model.getVars())}, {(t1-t0):.1f}")

if __name__ == '__main__':
    main()
