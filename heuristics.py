"""Implements several heuristics
to be used with a separate script to generate ``.csv`` logs.

Heuristics for the original problem generally have ``orig_`` prefix
in the name, for the simplified problem -- ``simpl_`` prefix.

The module also contains some code related to possible further
research and testing.
"""

import sys
import numpy as np

import opt_parser as als
import varseq as vs
import BB_search as BB

from glob import glob
from copy import deepcopy, copy

# for debugging // profiling
import cProfile
import pstats


def toA(A, B):
    """Helper: align to ``A`` heuristic."""
    Bp = B.align_to(A)
    return [A.size() + Bp.size(), A.layer_var]


def toB(A, B):
    """Helper: align to ``B`` heuristic."""
    return toA(B, A)


# min { toA, toB}
def minAB(A, B):
    """Best of A and B heuristic (simplified problem)."""
    sA, vA = toA(A, B)
    sB, vB = toA(B, A)

    if sA <= sB:
        return sA, vA
    else:
        return sB, vB


def simpl_greedy_swaps(A, B):
    """finds a good alignment target with greedy swaps

    Makes swaps while it can improve the objective
    (picking the best improvement at each step)
    """
    cur_size, T = minAB(A, B)
    T = vs.VarSeq(T, [1 for _ in range(len(T))])
    improvement_possible = True

    while improvement_possible:
        improvement_possible = False
        Tp = deepcopy(T)

        for i in range(len(T) - 1):
            a = T.layer_var[i]
            Ap = A.align_to(T.slide(a, i - 1))
            Bp = B.align_to(T.slide(a, i - 1))

            if Ap.size() + Bp.size() < cur_size:
                cur_size = Ap.size() + Bp.size()
                Tp = T.slide(a, i - 1)
                improvement_possible = True
        T = Tp

    return [A.align_to(T).size() + B.align_to(T).size(), T]


# greedy sifts
def simpl_greedy_sifts(A, B, passes=-1):
    """Finds a good alignment target with greedy sifts

    Makes sifts while it can improve the objective
    or for a given number of passes,
    picking the best improvement at each step.
    (one pass = examining all the variables)
    """
    cur_size, T = minAB(A, B)
    T = vs.VarSeq(T, [1 for _ in range(len(T))])
    improvement_possible = True
    pass_no = 0

    while improvement_possible:
        improvement_possible = False
        pass_no += 1
        Tp = deepcopy(T)

        for i in range(len(T)):
            a = T.layer_var[i]
            for j in range(len(T)):
                Ap = A.align_to(T.slide(a, j))
                Bp = B.align_to(T.slide(a, j))

                if Ap.size() + Bp.size() < cur_size:
                    cur_size = Ap.size() + Bp.size()
                    Tp = T.slide(a, j)
                    improvement_possible = True
        T = Tp
        if pass_no == passes:
            break

    return [A.align_to(T).size() + B.align_to(T).size(), T]


def simpl_gsifts_1p(A, B):
    """Greedy sifts: one pass."""
    return simpl_greedy_sifts(A, B, passes=1)


def simpl_gsifts_2p(A, B):
    """Greedy sifts: two passes."""
    return simpl_greedy_sifts(A, B, passes=2)


def simpl_gsifts_inf(A, B):
    """Greedy sifts: all passes."""
    return simpl_greedy_sifts(A, B, passes=-1)


def simpl_gsifts_3p(A, B):
    """Greedy sifts: three passes."""
    return simpl_greedy_sifts(A, B, passes=3)


def simpl_5random(A, B):
    """Implements 'Best of 5 random orders' heuristic."""
    best_size = -1
    for i in range(5):
        o = np.random.permutation(A.layer_var)
        Ap = A.align_to(o)
        Bp = B.align_to(o)
        if best_size < 0 or Ap.size() + Bp.size() < best_size:
            best_size = Ap.size() + Bp.size()
            best_o = copy(Ap.layer_var)
    return [best_size, best_o]


def simpl_greedy_2sifts(A, B, passes=-1):
    """Finds a good alignment target with greedy pairs of sifts.

    Makes sifts while it can improve the objective,
    or for a given number of passes,
    picking the best improvement at each step.
    (one pass = examining all the variables.)
    """
    cur_size, T = minAB(A, B)
    T = vs.VarSeq(T, [1 for i in range(len(T))])
    improvement_possible = True
    pass_no = 0

    while improvement_possible:
        improvement_possible = False
        pass_no += 1
        Tp = deepcopy(T)

        for i1 in range(len(T)):
            a1 = T.layer_var[i1]
            for j1 in range(len(T)):
                for i2 in range(len(T)):
                    a2 = T.layer_var[i2]
                    for j2 in range(len(T)):
                        Ap = A.align_to(T.slide(a1, j1).slide(a2, j2))
                        Bp = B.align_to(T.slide(a1, j1).slide(a2, j2))

                        if Ap.size() + Bp.size() < cur_size:
                            cur_size = Ap.size() + Bp.size()
                            Tp = T.slide(a1, j1).slide(a2, j2)
                            improvement_possible = True
        T = Tp
        if pass_no == passes:
            break

    return [A.align_to(T).size() + B.align_to(T).size(), T]


def simpl_g2sifts_1p(A, B):
    """Experimental: greedy pairs of sifts (simplified problem), one pass."""
    return simpl_greedy_2sifts(A, B, 1)


def simpl_g2sifts_2p(A, B):
    """Experimental: greedy pairs of sifts (simplified problem), two passes."""
    return simpl_greedy_2sifts(A, B, 2)


######################################################################
## an optimized version for ``sift`` operation.
def fastslide(X, a, p):
    """An experimental function implementing faster 'sift' operation."""
    if X[p] == a:
        return X

    X = [x for x in X]
    # so we definitely change something

    for i in range(len(X)):
        if X[i] != a and p > i:
            continue  # nothing to do, yet

        elif X[i] == a:
            # so, p>i -- slide down
            j = i
            while j < p:
                X[j] = X[j + 1]
                j += 1

            X[p] = a
            break

        elif p == i:
            # so, p(a)>i -- slide up
            buf1 = X[i]
            X[i] = a
            j = i + 1
            while X[j] != a:
                buf2 = X[j]
                X[j] = buf1
                buf1 = buf2
                j += 1

            X[j] = buf1
            break

    return X


def fast_greedy_sifts(A, B):
    """Finds a good alignment target with greedy sifts.

    Makes sifts while it can improve the objective
    (picking the best improvement at each step)
    """
    cur_size, T = minAB(A, B)
    improvement_possible = True
    no_of_ops = 0

    while improvement_possible:
        improvement_possible = False
        Tp = T

        for i in range(len(T)):
            a = T[i]
            for j in range(len(T)):
                Tpp = fastslide(T, a, j)
                Ap = A.align_to(Tpp)
                Bp = B.align_to(Tpp)

                no_of_ops += 2 * abs(i - j) + 2 * len(A)

                cand_size = Ap.size() + Bp.size()
                if cand_size < cur_size:
                    cur_size = cand_size
                    Tp = Tpp
                    no_of_ops += abs(i - j)
                    improvement_possible = True
        T = Tp

    return [T, A.align_to(T).size() + B.align_to(T).size(), no_of_ops]


def fast_greedy_2sifts(A, B):
    """finds a good alignment target with greedy pairs of sifts

    Makes sifts while it can improve the objective
    (picking the best improvement at each step)
    """
    cur_size, T = minAB(A, B)
    improvement_possible = True
    no_of_ops = 0

    while improvement_possible:
        # print("obj:{}".format(cur_size))
        improvement_possible = False
        Tp = T
        for i1 in range(len(T) - 1):
            a1 = T[i1]
            for j1 in range(i1 + 1, len(T)):
                Tpi = fastslide(T, a1, j1)

                for i2 in range(i1 + 1, len(T) - 1):
                    a2 = T[i2]
                    for j2 in range(i2 + 1, len(T)):
                        Tpp = fastslide(Tpi, a2, j2)
                        Ap = A.align_to(Tpp)
                        Bp = B.align_to(Tpp)
                        no_of_ops += 2 * len(A)
                        cand_size = Ap.size() + Bp.size()
                        if cand_size < cur_size:
                            cur_size = cand_size
                            Tp = Tpp
                            no_of_ops += abs(i1 - j1) + abs(i2 - j2)
                            improvement_possible = True
        T = Tp

    return [T, A.align_to(T).size() + B.align_to(T).size(), no_of_ops]


def orig_simpl(A, B, simpl):
    """A heuristic for the **original** problem involving simplified problem."""
    Ap = deepcopy(A)
    Bp = deepcopy(B)
    Ap.align_to(simpl["simpl_BB"][2], inplace=True)
    Bp.align_to(simpl["simpl_BB"][2], inplace=True)
    assert Ap.is_aligned(Bp)
    return [Ap.size() + Bp.size(), simpl["simpl_BB"][1], simpl["simpl_BB"][2]]


def orig_gsifts1p(A, B, simpl):
    """Greedy sifts heuristic for the **original** problem."""
    ApA = deepcopy(A)
    BpA = deepcopy(B)
    ApB = deepcopy(A)
    BpB = deepcopy(B)
    ApB.gsifts(BpB)
    BpA.gsifts(ApA)
    sA = ApA.size() + BpA.size()
    sB = ApB.size() + BpB.size()

    # make sure sequences are aligned
    if not ApA.is_aligned(BpA):
        sA = -1
    if not ApB.is_aligned(BpB):
        sB = -1

    if sA <= sB and sA > 0:
        return [sA, 0, ApA.vars]
    elif sB < sA and sB > 0:
        return [sB, 0, ApB.vars]
    else:
        return [-1, 0, ApA.vars]


def orig_bestAB(A, B, simpl):
    """Best of ``A`` and ``B`` (original problem)."""
    Ap = deepcopy(A)
    Bp = deepcopy(B)
    Ap.align_to(B.vars, inplace=True)
    Bp.align_to(A.vars, inplace=True)
    if Ap.size() + B.size() <= A.size() + Bp.size():
        return [Ap.size() + B.size(), 0, Ap.vars]
    else:
        return [A.size() + Bp.size(), 0, Bp.vars]


def orig_5random(A, B, simpl):
    """Five random orders (original problem)."""
    best_size = -1
    for i in range(5):
        Ap = deepcopy(A)
        Bp = deepcopy(B)
        o = np.random.permutation(A.vars)
        Ap.align_to(o, inplace=True)
        Bp.align_to(o, inplace=True)
        if best_size < 0 or Ap.size() + Bp.size() < best_size:
            best_size = Ap.size() + Bp.size()
            best_o = copy(Ap.vars)
    return [best_size, 0, best_o]


def orig_interleaved(A, B, simpl):
    """Experimental heuristic (further research)."""
    N = len(A)
    Nlast = N
    Ac = deepcopy(A)
    Bc = deepcopy(B)

    for i in range(Nlast):
        if Ac.vars[i] == Bc.vars[i]:
            continue

        vsA = vs.VarSeq(Ac.vars[i:], [Ac.n(j) for j in range(i, N)])
        vsB = vs.VarSeq(Bc.vars[i:], [Bc.n(j) for j in range(i, N)])
        b = BB.BBSearch(vsA, vsB)
        b.search()
        cand_var = b.Ap_cand.layer_var[0]
        Ac.sift(cand_var, i)
        Bc.sift(cand_var, i)

    i = Nlast
    vsAF = vs.VarSeq(Ac.vars[i:], [Ac.n(j) for j in range(i, N)])
    vsBF = vs.VarSeq(Bc.vars[i:], [Bc.n(j) for j in range(i, N)])
    b = BB.BBSearch(vsAF, vsBF)
    b.search()
    for i in range(Nlast, N):
        cand_var = b.Ap_cand.layer_var[i - Nlast]
        Ac.sift(cand_var, i)
        Bc.sift(cand_var, i)

    assert Ac.is_aligned(Bc)

    return [Ac.size() + Bc.size(), 0, Ac.vars]


def orig_meta(A, B, simpl):
    """Another experimental heuristic (further research)."""
    N = len(A)
    Ac = deepcopy(A)
    Bc = deepcopy(B)
    vsA = vs.VarSeq(Ac.vars, [Ac.n(j) for j in range(N)])
    vsB = vs.VarSeq(Bc.vars, [Bc.n(j) for j in range(N)])
    b = BB.BBSearch(vsA, vsB)
    b.search()

    Ac.gsifts(Bc, start_order=b.Ap_cand.layer_var)

    assert Ac.is_aligned(Bc)
    return [Ac.size() + Bc.size(), 0, Ac.vars]


def orig_interleave_when_diverge(A, B, simpl):
    """Another experimental heuristic (further research)."""
    ALLOWANCE = 0.99
    Ap = deepcopy(A)
    Bp = deepcopy(B)
    N = len(Ap)
    vsA = vs.VarSeq(A.vars, [A.n(i) for i in range(N)])
    vsB = vs.VarSeq(B.vars, [B.n(i) for i in range(N)])

    simpl_order = simpl["simpl_BB"][2]
    i = 0
    while i < N:
        if Ap.vars[i] == Bp.vars[i]:
            i += 1
            continue
        Ap.sift(simpl_order[i], i)
        Bp.sift(simpl_order[i], i)
        vsA.slide(simpl_order[i], i, inplace=True)
        vsB.slide(simpl_order[i], i, inplace=True)
        if (Ap.size() + Bp.size() < ALLOWANCE *
            (vsA.size() + vsB.size())) and i < 8:
            # we diverged too much -- need to re-solve
            # rebuild varseqs
            vsA = vs.VarSeq(Ap.vars, [Ap.n(i) for i in range(N)])
            vsB = vs.VarSeq(Bp.vars, [Bp.n(i) for i in range(N)])
            # solve the simplified problem
            b = BB.BBSearch(vsA, vsB)
            b.search()
            simpl_order = b.Ap_cand.layer_var

        i += 1

    assert Ap.is_aligned(Bp)
    return [Ap.size() + Bp.size(), simpl["simpl_BB"][1], simpl_order]


def orig_rnd_starts(A, B, simpl):
    """Another experimental heuristic (further research)."""
    best_o = simpl["simpl_BB"][2]  # retrieve opt order
    best_size = deepcopy(A).align_to(best_o).size() + deepcopy(B).align_to(
        best_o).size()

    N = len(A)
    for i in range(5):
        Ap = deepcopy(A)
        Bp = deepcopy(B)

        Ap = Ap.align_to(np.random.permutation(A.vars))
        Bp = Bp.align_to(np.random.permutation(B.vars))

        vsA = vs.VarSeq(Ap.vars, [Ap.n(i) for i in range(N)])
        vsB = vs.VarSeq(Bp.vars, [Bp.n(i) for i in range(N)])
        # solve the simplified problem
        b = BB.BBSearch(vsA, vsB)
        b.search()
        simpl_order = b.Ap_cand.layer_var
        Ap = Ap.align_to(simpl_order)
        Bp = Bp.align_to(simpl_order)

        if (Ap.size() + Bp.size()) < best_size:
            best_size = Ap.size() + Bp.size()
            best_o = copy(Ap.vars)

    return [best_size, 0, best_o]


# heuristic structure
# CODE - FUNC - LEGEND
SIMPL_HEU = [
    # ["simpl_toA",toA, "align to A"],
    # ["simpl_toB",toB, "align to B"],
    ["simpl_minAB", minAB, "Best of A and B"],
    ["simpl_gswaps", simpl_greedy_swaps, "Greedy swaps"],
    ["simpl_5random", simpl_5random, "Best/five random orders"],
    ["simpl_gsifts_1p", simpl_gsifts_1p, "Greedy sifts (one pass)"],
    ["simpl_gsifts_2p", simpl_gsifts_2p, "Greedy sifts (two passes)"],
    ["simpl_gsifts_inf", simpl_gsifts_inf, "Greedy sifts (all)"]
    # ["simpl_g2sifts_1p",simpl_g2sifts_1p,"greedy 2-sifts (1 pass)"]
]
"""List of heuristics for the **simplified** problem."""

ORIG_HEU = [
    ["orig_simpl", orig_simpl, "Simplified problem"],
    ["orig_gsifts1p", orig_gsifts1p, "Greedy BDD sifts (one pass)"],
    ["orig_bestAB", orig_bestAB, "Best of A and B"],
    ["orig_5random", orig_5random, "Best of five random orders"]
    # ["orig_ild",orig_interleave_when_diverge,"interleave-when-diverge"],
    #["orig_simpl_rnd",orig_rnd_starts,"Simplified problem w/randomized starts"]
    # ["orig_meta",orig_meta,"greedy sifts + simplified problem"]
]
"""List of heuristics for the **original** problem."""

# NOTE: Legacy code below may no longer work
#       (i.e., the file is to be used with experiment scripts only)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("USAGE: {} <inst_directory>/".format(sys.argv[0]))
        exit()

    inst_dir = sys.argv[1]

    print("inst_no,opt_obj," + "_obj,".join(heu_description) + "_obj" + "," +
          "_steps,".join(heu_description) + "_steps")

    pr = cProfile.Profile()
    tg_sendtext(
        "LOCAL JOB: heuristics testing started, instance dir: {}".format(
            inst_dir))
    n = 0
    pr.enable()

    for file in glob(inst_dir + "vs*.md"):
        inst_no = file.split('.')[-2]
        inst_no = inst_no.split('/vs')[-1]

        inst = als.alignments(file)
        inst.load(verbose=False)

        A = inst.As[0]
        B = inst.Bs[0]

        objs = []
        steps = []

        for h in heuristics:
            T, obj, st = h(A, B)
            objs.append(str(obj))
            steps.append(str(st))

        print(inst_no + "," + str(inst.As_aligned[0][0].size() +
                                  inst.Bs_aligned[0][0].size()) + "," +
              ",".join(objs) + "," + ",".join(steps))

        sys.stdout.flush()
        n += 1

    pr.disable()
    ps = pstats.Stats(pr).sort_stats('time')

    ps.dump_stats("./run_profile.dmp")
    tg_sendtext("LOCAL (heuristics): done. {} files processed".format(n))
