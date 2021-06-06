"""Benchmarks runtimes for different solution methods (colored-UFL)

An experiment concerning the colored Uncapacitated Facility Location (UFL).

---
(c) A. Bochkarev, Clemson University, 2021
abochka@clemson.edu
"""
from time import time
import argparse as ap
import sys
import cUFL
from gurobipy import GRB
import UFL
import varseq as vs
import BDD as DD
import BB_search as bb
from copy import deepcopy
import json


def show_header():
    """Shows the table header."""
    print("run,k,n,method,step,time")


def benchmark(K=10, TOL=1e-3, n=5, prefix=0, do_sifts=False, logging=None):
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
        print(f"{prefix},{k},{n},gen_instance,all, {(t1-t0)*1000:.1f}")
        if logging != None:
            logging.write(json.dumps({
                'S':S, 'f':f, 'fc':fc, 'kb':kb})+"\n")

        # Generate and solve plain MIP
        t0 = time()
        model = cUFL.build_cUFL_MIP(S, f, fc, kb)
        model.setParam("OutputFlag", 0)
        model.optimize()
        t1 = time()
        print(f"{prefix},{k},{n},plain_MIP,build+solve,{(t1-t0)*1000:.1f}")

        if model.status != GRB.OPTIMAL:
            print(f"Plain MIP status is: {model.status}")
            print(f"\nS={S}; f={f}; fc={fc}; kb={kb}")
            continue

        plain_MIP_obj = model.objVal

        # Generate and solve CPP MIP
        t0 = time()
        cover_DD, _ = cUFL.build_cover_DD(S, f)
        cover_DD.make_reduced()
        pref_order = [int(x[1:]) for x in cover_DD.vars]
        color_DD, _ = cUFL.build_color_DD(f, fc, kb, pref_order)
        color_DD.make_reduced()
        t1 = time()
        DD_build_time = t1 - t0
        print(f"{prefix},{k},{n},CPP,BDD-build,{DD_build_time*1000:.1f}")

        t0 = time()
        m, c, v, x = UFL.add_BDD_to_MIP(cover_DD, prefix="cover")
        m, c, v, x = UFL.add_BDD_to_MIP(color_DD, m, x, "color")
        m.setParam("OutputFlag", 0)
        t1 = time()
        print(f"{prefix},{k},{n},CPP,MIP-build,{(t1-t0)*1000:.1f}")

        t0 = time()
        m.optimize()
        t1 = time()

        print(f"{prefix},{k},{n},CPP,MIP-solve,{(t1-t0)*1000:.1f}")

        if m.status != GRB.OPTIMAL:
            print(f"CPP MIP status is: {m.status}")
            print(f"\nS={S}; f={f}; fc={fc}; kb={kb}")
            continue

        CPP_MIP_obj = m.objVal

        t0 = time()
        vs_cover = vs.VarSeq(cover_DD.vars, [len(L) for L in cover_DD.layers[:-1]])
        vs_color = vs.VarSeq(color_DD.vars, [len(L) for L in color_DD.layers[:-1]])
        t1 = time()

        assert set(vs_cover.layer_var) == set(vs_color.layer_var), f"cover:{vs_cover.layer_var}, color:{vs_color.layer_var}"

        b = bb.BBSearch(vs_cover, vs_color)
        bb.TIMEOUT_ITERATIONS=10000
        status = b.search()
        assert status == "optimal" or status == "timeout"
        t1 = time()
        print(f"{prefix},{k},{n},CPP,VS-build+solve,{(t1-t0)*1000:.1f}")

        t0 = time()
        cover_p = cover_DD.align_to(b.Ap_cand.layer_var, inplace=False)
        color_p = color_DD.align_to(b.Ap_cand.layer_var, inplace=False)
        t1 = time()
        print(f"{prefix},{k},{n},CPP,BDD-align-to-vs,{(t1-t0)*1000:.1f}")

        t0 = time()
        int_DD = DD.intersect(cover_p, color_p)
        t1 = time()
        print(f"{prefix},{k},{n},CPP,intersection-build,{(t1-t0)*1000:.1f}")

        t0 = time()
        nl = int_DD.shortest_path()
        t1 = time()

        print(f"{prefix},{k},{n},CPP,intersection-SP-solve,{(t1-t0)*1000:.1f}")
        NF_obj = nl[DD.NROOT]

        if abs(CPP_MIP_obj - plain_MIP_obj) >= TOL:
            print("! (plain vs CPP)")
            print(f"\nS={S}; f={f}; fc={fc}; kb={kb}")
        elif (abs(plain_MIP_obj - NF_obj) >= TOL):
            print("! (plain vs aDD)")
            print(f"plain obj={plain_MIP_obj}, while aBDD obj={NF_obj}")
            print(f"\nS={S}; f={f}; fc={fc}; kb={kb}")

        if do_sifts:
            t0 = time()
            cover_pp = deepcopy(cover_DD)
            color_pp = deepcopy(color_DD)
            color_pp.gsifts(cover_pp)
            t1 = time()
            print(f"{prefix},{k},{n},CPP,BDD-align-gsifts,{(t1-t0)*1000:.1f}")

            t0 = time()
            int_DD_gsifts = DD.intersect(cover_pp, color_pp)
            t1 = time()
            print(f"{prefix},{k},{n},CPP,intersection-gsifts-build,{(t1-t0)*1000:.1f}")

            t0 = time()
            nl = int_DD_gsifts.shortest_path()
            t1 = time()

            print(f"{prefix},{k},{n},CPP,intersection-gsifts-SP-solve,{(t1-t0)*1000:.1f}")
            NF_obj = nl[DD.NROOT]

            if (abs(plain_MIP_obj - NF_obj) >= TOL):
                print("! (plain vs gsifts)")
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
    parser.add_argument("-s",
                        "--with-gsifts",
                        action="store_true",
                        dest="do_sifts",
                        help="perform greedy sifts (slow)")
    parser.add_argument("-p",
                        "--prefix",
                        action="store",
                        dest="prefix",
                        default="0",
                        help="prefix string (ID of the run)")
    parser.add_argument("-l",
                        "--instance-log",
                        action="store",
                        dest="inst_log",
                        default="none",
                        help="file to save the instances (json)")

    args = parser.parse_args()

    if args.header:
        show_header()
        exit(0)

    if args.inst_log == "none":
        inst_log = None
    else:
        inst_log = open(args.inst_log, "w")

    benchmark(int(args.K), n=int(args.n), prefix="test", do_sifts=args.do_sifts, logging=inst_log)

    if not inst_log is None:
        inst_log.close()

if __name__ == '__main__':
    main()
