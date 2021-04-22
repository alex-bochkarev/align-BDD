"""
Examines heuristic solutions to a random align-BDD instance.
(inspired by the reviewer's comment)

---
(c) A. Bochkarev, Clemson University, 2021
abochka@clemson.edu
"""
import BDD as DD
import varseq as vs
import BB_search as bb
import numpy as np
from copy import deepcopy
import argparse as ap
import sys
from time import time
import cUFL


def run_experiment(N=5, k=0, inst_type="rnd"):
    """Runs the experiment with diagrams of `N` layers."""
    # prepare an align-BDD instance
    t0 = time()
    if inst_type == "rnd":
        A = DD.BDD.random(N=N)
        B = DD.BDD.random(N=N)
        B.rename_vars(
            dict(
                zip([i for i in range(1, N + 1)],
                    np.random.permutation([i for i in range(1, N + 1)]))))
    elif inst_type == "tUFL":
        S, f, fc, kb = cUFL.generate_test_instance(n=N)
        A, _ = cUFL.build_cover_DD(S, f)
        pref_order = [int(x[1:]) for x in A.vars]
        B, _ = cUFL.build_color_DD(f, fc, kb, pref_order)
    else:
        print(f"Wrong instance inst_type: {inst_type} (expected: 'rnd' or 'tUFL')")

    A.make_reduced()
    B.make_reduced()

    AB_simscore = A.simscore(B)

    # create and solve the simplified problem
    vsA = vs.VarSeq(A.vars, [A.n(i) for i in range(N)])
    vsB = vs.VarSeq(B.vars, [B.n(i) for i in range(N)])
    b = bb.BBSearch(vsA, vsB)
    b.search()

    Ap = deepcopy(A)
    Bp = deepcopy(B)

    Ap.align_to(b.Ap_cand.layer_var, inplace=True); Bp.align_to(b.Ap_cand.layer_var, inplace=True)
    obj = Ap.size()+Bp.size()

    # calculate (brute-force) the real optima
    no_opts, A_aligned, B_aligned = A.OA_bruteforce(B)
    opts_diam = None
    best_VS_simscore = None
    worst_VS_simscore = None

    seq_VS = b.Ap_cand.layer_var

    if no_opts == 1:
        opts_diam = 0.0
        VS_simscore = DD.simscore(seq_VS, A_aligned[0].vars)
        best_VS_simscore = VS_simscore
        worst_VS_simscore = best_VS_simscore

    for i in range(no_opts-1):
        for j in range(i, no_opts):
            dist = 1-DD.simscore(A_aligned[i].vars, A_aligned[j].vars)
            if opts_diam is None or dist > opts_diam:
                opts_diam = dist

            VS_simscore = DD.simscore(seq_VS, A_aligned[i].vars)
            if best_VS_simscore is None or best_VS_simscore < VS_simscore:
                best_VS_simscore = VS_simscore

            if worst_VS_simscore is None or worst_VS_simscore > VS_simscore:
                worst_VS_simscore = VS_simscore

    t1 = time()

    print(f"{k}, {b.status}, {AB_simscore*100:.1f}, {no_opts}, {opts_diam*100:.1f}, {best_VS_simscore*100:.1f}, {worst_VS_simscore*100:.1f}, {t1-t0:.1f}")

if __name__ == '__main__':
    parser = ap.ArgumentParser(description=" Experiment: analyzing heuristic solutions. (c) A. Bochkarev, 2021",
                               formatter_class=ap.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-n", action="store", dest="n", help="number of layers")
    parser.add_argument("-k", action="store", dest="k", help="number of instances to generate")
    parser.add_argument("-t", action="store", dest="inst_type", help="instance type 'rnd' (random) or 'tUFL", default="rnd")

    args = parser.parse_args()
    assert int(args.k) > 0, "Please specify a positive no of instances to generate."
    assert int(args.n) > 0, "Please specify number of layers per instance (e.g., 5))"
    assert args.inst_type in ['rnd', 'tUFL'], "wrong instance type ('rnd' or 'tUFL' expected)"

    print("k, VS_status, AB_simscore, no_opts, opts_diam, best_VS_simscore, worst_VS_simscore, exp_time")
    for k in range(int(args.k)):
        run_experiment(N=int(args.n), k=k)
        sys.stdout.flush()
