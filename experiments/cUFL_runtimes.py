"""Benchmarks runtimes for different solution methods (colored-UFL)

An experiment concerning the colored Uncapacitated Facility Location (UFL).

---
(c) A. Bochkarev, Clemson University, 2021
abochka@clemson.edu
"""
from time import time
import argparse as ap
import sys
from gurobipy import GRB
import UFL
import cUFL
import varseq as vs
import BDD as DD
import BB_search as bb


def show_header():
    """Shows the table header."""
    print("run,k,n,method,step,time")


def benchmark(K=10, TOL=1e-3, n=5, prefix=0):
    """Runs the solution for three different methods.

    Compares plain MIP, CPP MIP, and network-flow from DD (LP)
    in terms of objectives and time.

    Args:
        K (int): number of instances to generate (default 10)
        TOL (float): objective tolerance (for comparison) (default 1e-3)
        n (int): number of facilities (default 5)

    Returns:
        Nothing (outputs the data to `stdout`)

    Output format (stdout, comma-separated):
        prefix(str): run ID
        k (int): instance number
        n (int): no. of customers
        method (str): method code
        step (str): step code
        time (float): duration in msec.
    """
    for k in range(K):
        t0 = time()
        S, f, fc, kb = cUFL.generate_test_instance(n=n)
        t1 = time()
        print(f"{prefix},{k},{n},gen_instance,all, {(t1-t0)*1000:.3f}")

        # Generate and solve plain MIP
        t0 = time()
        model = cUFL.build_cUFL_MIP(S, f, fc, kb)
        model.setParam("OutputFlag", 0)
        model.optimize()
        t1 = time()
        print(f"{prefix},{k},{n},plain_MIP,all,{(t1-t0)*1000:.3f}")

        if model.status != GRB.OPTIMAL:
            print(f"Plain MIP status is: {model.status}")
            print(f"\nS={S}; f={f}; fc={fc}; kb={kb}")
            continue

        plain_MIP_obj = model.objVal

        # Generate and solve CPP MIP
        t0 = time()
        cover_DD, cover_nl = cUFL.build_cover_DD(S, f)
        pref_order = [int(x[1:]) for x in cover_DD.vars]
        color_DD, color_nl = cUFL.build_color_DD(f, fc, kb, pref_order)
        t1 = time()
        DD_build_time = t1 - t0

        t0 = time()
        m, c, v, x = UFL.add_BDD_to_MIP(cover_DD, prefix="cover")
        m, c, v, x = UFL.add_BDD_to_MIP(color_DD, m, x, "color")
        m.update()
        m.setParam("OutputFlag", 0)
        m.optimize()
        t1 = time()
        print(f"{prefix},{k},{n},CPP_MIP,build_and_solve,{(DD_build_time+t1-t0)*1000:.1f}")

        if m.status != GRB.OPTIMAL:
            print(f"CPP MIP status is: {m.status}")
            print(f"\nS={S}; f={f}; fc={fc}; kb={kb}")
            continue

        CPP_MIP_obj = m.objVal

        t0 = time()
        vs_cover = vs.VarSeq(cover_DD.vars, [len(L) for L in cover_DD.layers[:-1]])
        vs_color = vs.VarSeq(color_DD.vars, [len(L) for L in color_DD.layers[:-1]])

        assert set(vs_cover.layer_var) == set(vs_color.layer_var), f"cover:{vs_cover.layer_var}, color:{vs_color.layer_var}"
        b = bb.BBSearch(vs_cover, vs_color)

        # bb.TIMEOUT_ITERATIONS=10000
        status = b.search()
        assert status == "optimal" or status == "timeout"

        cover_p = cover_DD.align_to(b.Ap_cand.layer_var, inplace=False)
        color_p = color_DD.align_to(b.Ap_cand.layer_var, inplace=False)

        int_DD = DD.intersect(cover_p, color_p)
        m, c, v = UFL.create_NF(int_DD)
        m.setParam("OutputFlag", 0)
        m.optimize()
        t1 = time()

        print(f"{prefix},{k},{n},int_DD_based,build_and_solve,{(DD_build_time+t1-t0)*1000:.1f}")
        NF_obj = m.objVal

        if abs(CPP_MIP_obj - plain_MIP_obj) >= TOL:
            print("! (plain vs CPP)")
            print(f"\nS={S}; f={f}; fc={fc}; kb={kb}")
        elif (abs(plain_MIP_obj - NF_obj) >= TOL):
            print("! (plain vs aDD)")
            print(f"plain obj={plain_MIP_obj}, while aBDD obj={NF_obj}")
            print(f"\nS={S}; f={f}; fc={fc}; kb={kb}")

        sys.stdout.flush()

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
        show_header()
        exit(0)

    benchmark(int(args.K), n=int(args.n), prefix="test")

if __name__ == '__main__':
    main()
