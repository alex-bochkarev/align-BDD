"""Tests BB search correctness (:py:mod:`BB_search`)

Compares the objective obtained by the BB search vs. the brute-force
enumerated, true-optimal objective (to the **simplified** problem).
Tests against randomly generated instances (not necessarily unique)
"""
import pytest
import sys
sys.path.append('..')

import varseq as vs
import BB_search as bb

@pytest.mark.parametrize("N", [6 for _ in range(100)])
def test_BB_search(N):
    """Implements the test."""
    A = vs.VarSeq.random(N = N)
    B = vs.VarSeq.random(N = N)

    b = bb.BBSearch(A,B)

    _, Ap, Bp = A.OA_bruteforce(B)

    opt_size = Ap[0].size() + Bp[0].size()

    status = b.search()
    assert status in ['optimal', 'timeout']

    BB_size = b.Ap_cand.size() + b.Bp_cand.size()

    assert BB_size == opt_size, \
        f"Instance:\n{A}\nvs.\n{B}\nBB size={BB_size}, opt size ={opt_size}"
