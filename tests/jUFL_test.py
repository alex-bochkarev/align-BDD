"""Tests joint-UFLP code (:py:mod:`jUFL`).

Essentially, tests if the construction of diagrams is implemented
correctly.
"""
import pytest
from jUFL import *


@pytest.mark.parametrize("test_inst", [generate_instance(np.random.randint(5, 15))
                                       for _ in range(100)])
def test_jUFL(test_inst):
    """Implements the testing code.

    Generates random jUFLP instances and compares
    "Naive MIP" and "CPP MIP".
    """
    m, _ = solve_with_naive_MIP(test_inst)
    o_MIP = m.objVal
    assert m.status == GRB.OPTIMAL

    m, _ = solve_with_DD_MIP(test_inst)
    o_DD_MIP = m.objVal
    assert m.status == GRB.OPTIMAL

    assert abs(o_MIP - o_DD_MIP) < 0.00001, f"Naive MIP ({o_MIP}) vs DD MIP ({o_DD_MIP}): a discrepancy"
