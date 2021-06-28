"""Tests jUFL.py module (joint-UFLP code).

---
(c) A. Bochkarev, Clemson University, 2021
abochka@g.clemson.edu
"""
import pytest
from jUFL import *


@pytest.mark.parametrize("test_inst", [generate_instance(np.random.randint(5, 15))
                                       for _ in range(100)])
def test_jUFL(test_inst):
    m, _ = solve_with_naive_MIP(test_inst)
    o_MIP = m.objVal
    assert m.status == GRB.OPTIMAL

    m, _ = solve_with_DD_MIP(test_inst)
    o_DD_MIP = m.objVal
    assert m.status == GRB.OPTIMAL

    assert abs(o_MIP - o_DD_MIP) < 0.00001, f"Naive MIP ({o_MIP}) vs DD MIP ({o_DD_MIP}): a discrepancy"
