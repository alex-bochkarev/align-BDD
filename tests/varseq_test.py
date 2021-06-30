"""Tests the :py:class:`varseq.VarSeq` procedures."""
import pytest
import sys
sys.path.append('..')

import numpy as np
import varseq as vs

@pytest.mark.parametrize("N", [np.random.randint(5,20) for i in range(500)])
def test_align_to(N):
    """Tests the correctness of :py:func:`varseq.VarSeq.align_to`.

    Generates random variable sequences (of length ``N``) and tests
    :py:func:`VarSeq.align_to` against a naive implementation of
    :py:func:`VarSeq.greedy_sort`.
    """
    A = vs.VarSeq.random(N=N)
    B = vs.VarSeq.random(N=N)

    Ap1 = A.greedy_sort(B.layer_var)
    Ap2 = A.align_to(B.layer_var)

    assert (np.array_equal(Ap1.layer_var,Ap2.layer_var) and \
            np.array_equal(Ap1.n, Ap2.n)), \
            f"Instance:\n{A}\nvs.\n{B}"
