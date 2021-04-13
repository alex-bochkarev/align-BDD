"""
Preorder effect on the intersection DD siez: a quick experiment.

---
(c) A. Bochkarev, Clemson University, 2021
abochka@clemson.edu
"""
from BDD import intersect
import cUFL
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import varseq as vs
import BB_search as bb
import numpy as np
import sys


def get_sizes(n=10, K=1):
    """Gets intersection DD sizes with and without color DD preordering."""
    int_sizes = pd.DataFrame(columns=['preorder', 'no_preorder', 'pre_factor'])

    for k in range(K):
        S, f, fc, kb = cUFL.generate_test_instance(n)
        cov, _ = cUFL.build_cover_DD(S, f)
        cov.make_reduced()

        # ==================================================
        # (1) no preordering
        col_nopre, _ = cUFL.build_color_DD(f, fc, kb)
        col_nopre.make_reduced()

        vs_col = vs.VarSeq(col_nopre.vars, [len(L) for L in col_nopre.layers[:-1]])
        vs_cov = vs.VarSeq(cov.vars, [len(L) for L in cov.layers[:-1]])

        bb.TIMEOUT_ITERATIONS = 10000
        b = bb.BBSearch(vs_col, vs_cov)
        status = b.search()
        assert status == "optimal" or status == "timeout"

        color_p = col_nopre.align_to(b.Ap_cand.layer_var, inplace=False)
        cover_p = cov.align_to(b.Ap_cand.layer_var, inplace=False)

        int_nopre = intersect(color_p, cover_p)
        int_nopre.make_reduced()

        # ==================================================
        # (2) with preordering
        pref_order = [int(x[1:]) for x in cov.vars]
        col_pre, _ = cUFL.build_color_DD(f, fc, kb, pref_order)
        col_pre.make_reduced()

        vs_col = vs.VarSeq(col_pre.vars, [len(L) for L in col_pre.layers[:-1]])
        vs_cov = vs.VarSeq(cov.vars, [len(L) for L in cov.layers[:-1]])

        bb.TIMEOUT_ITERATIONS = 10000
        b = bb.BBSearch(vs_col, vs_cov)
        status = b.search()
        assert status == "optimal" or status == "timeout"

        color_p = col_pre.align_to(b.Ap_cand.layer_var, inplace=False)
        cover_p = cov.align_to(b.Ap_cand.layer_var, inplace=False)

        int_pre = intersect(color_p, cover_p)
        int_pre.make_reduced()

        int_sizes = int_sizes.append({'preorder': int_pre.size(), 'no_preorder': int_nopre.size(), 'pre_factor': int_pre.size() / int_nopre.size()}, ignore_index=True)
        if k % 10 == 0:
            print(".", end="")
            sys.stdout.flush()

    return int_sizes


if __name__ == '__main__':
    int_sizes = get_sizes(K=500)
    print(f"No of cases where preorder made sense: {np.sum(int_sizes['pre_factor'] < 1.0)} out of {len(int_sizes['pre_factor'])}")
    sns.histplot(int_sizes, x="pre_factor", binwidth=0.1)
    plt.show()
