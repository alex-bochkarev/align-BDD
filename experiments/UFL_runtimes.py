"""Benchmarks runtimes for different solution methods (UFL)

An experiment concerning the Uncapacitated Facility Location (UFL).

---
(c) A. Bochkarev, Clemson University, 2021
abochka@clemson.edu
"""
from time import time
import argparse as ap
import sys
from gurobipy import GRB
import UFL


def benchmark(K=100, TOL=1e-3, n=5, m=7, prefix=0):
    """Runs the solution for three different methods.

    Compares plain MIP, CPP MIP, and DP DD (LP) in terms of
    objectives and time.

    Args:
        K (int): number of instances to generate (default 500)
        TOL (float): objective tolerance (for comparison) (default 1e-3)
        n (int): number of facilities (default 3)
        m (type): number of customers (default 4)

    Returns:
        Nothing (outputs the data to `stdout`)

    Output format (stdout, comma-separated):
        prefix(str): run ID
        k (int): instance number
        n (int): no. of customers
        m (int): no. of facilities
        method (str): method code
        step (str): step code
        duration (float): duration in msec.
    """
    for k in range(K):
        t0 = time()
        S, f, g = UFL.generate_test_instance(n=n, m=m)
        t1 = time()
        print(f"{prefix},{k},{n},{m},gen_instance,--,{(t1-t0)*1000:.3f}")

        # Generate and solve plain MIP
        t0 = time()
        model = UFL.build_MIP(S, f, g)
        model.setParam("OutputFlag", 0)
        t1 = time()
        print(f"{prefix},{k},{n},{m},plain_MIP,model_setup,{(t1-t0)*1000:.3f}")
        t0 = time()
        model.optimize()
        t1 = time()
        print(f"{prefix},{k},{n},{m},plain_MIP,solve,{(t1-t0)*1000:.3f}")

        if model.status != GRB.OPTIMAL:
            print(f"Plain MIP status is: {model.status}")
            print(f"\nS={S}; f={f}; g={g}")
            continue

        plain_MIP_obj = model.objVal

        # Generate and solve CPP MIP
        t0 = time()
        C = UFL.create_covering_BDD_wg(S, g)
        A = UFL.create_availability_BDD(S, f)
        t1 = time()
        print(f"{prefix},{k},{n},{m},CPP,BDD_setup,{(t1-t0)*1000:.3f}")

        t0 = time()
        model, _, _, x = UFL.add_BDD_to_MIP(A, prefix="A_")
        model, _, _, x = UFL.add_BDD_to_MIP(C, model=model, x=x, prefix="C_")
        model.update()
        model.setParam("OutputFlag", 0)
        t1 = time()
        print(f"{prefix},{k},{n},{m},CPP_MIP,model_setup,{(t1-t0)*1000:.3f}")

        t0 = time()
        model.optimize()
        t1 = time()
        print(f"{prefix},{k},{n},{m},CPP_MIP,model_solve,{(t1-t0)*1000:.3f}")

        if model.status != GRB.OPTIMAL:
            print(f"CPP MIP status is: {model.status}")
            print(f"\nS={S}; f={f}; g={g}")
            continue

        CPP_MIP_obj = model.objVal

        # t0 = time()
        # vs_A = vs.VarSeq(A.vars, [len(L) for L in A.layers[:-1]])
        # vs_C = vs.VarSeq(C.vars, [len(L) for L in C.layers[:-1]])
        # assert set(vs_A.layer_var) == set(vs_C.layer_var)
        # t1 = time()
        # print(f"{prefix},{k},{n},{m},CPP_intBDD,varseq_setup,{(t1-t0)*1000:.3f}")

        # t0 = time()
        # b = bb.BBSearch(vs_A, vs_C)
        # status = b.search()
        # assert status == "optimal" or status == "timeout"
        # t1 = time()
        # print(f"{prefix},{k},{n},{m},CPP_intBDD,varseq_solve,{(t1-t0)*1000:.3f}")

        # t0 = time()
        # Ap = A.align_to(b.Ap_cand.layer_var, inplace=False)
        # Cp = C.align_to(b.Ap_cand.layer_var, inplace=False)
        # t1 = time()
        # print(f"{prefix},{k},{n},{m},CPP_intBDD,BDD_revise,{(t1-t0)*1000:.3f}")

        # t0 = time()
        # int_DD = DD.intersect(Ap, Cp)
        # t1 = time()
        # print(f"{prefix},{k},{n},{m},CPP_intBDD,intBDD_setup,{(t1-t0)*1000:.3f}")

        t0 = time()
        D, _ = UFL.build_DP_DD(S, f, g)
        t1 = time()
        print(f"{prefix},{k},{n},{m},DP_BDD,DP_DD_setup,{(t1-t0)*1000:.3f}")

        t0 = time()
        model, c, v = UFL.create_NF(D)
        model.setParam("OutputFlag", 0)
        t1 = time()
        print(f"{prefix},{k},{n},{m},DP_BDD,NF_model_setup,{(t1-t0)*1000:.3f}")

        t0 = time()
        model.optimize()
        t1 = time()
        print(f"{prefix},{k},{n},{m},DP_BDD,NF_solve,{(t1-t0)*1000:.3f}")

        DP_obj = model.objVal

        if abs(CPP_MIP_obj - plain_MIP_obj) >= TOL:
            print("! (plain vs CPP)")
            print(f"\nS={S}; f={f}; g={g}")
        elif (abs(plain_MIP_obj - DP_obj) >= TOL):
            print("! (plain vs aDD)")
            print(f"plain obj={plain_MIP_obj}, while aBDD obj={DP_obj}")
            print(f"\nS={S}; f={f}; g={g}")

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
        print("run,k,n,m,method,stage,time")
        exit(0)

    benchmark(int(args.K), n=int(args.n), prefix="test")

if __name__ == '__main__':
    main()
