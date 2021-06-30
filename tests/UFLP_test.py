"""Tests the UFLP-related machinery (:py:mod:`tUFLP`).

Note:
    The module requires Gurobi solver.
"""
import pytest
import numpy as np
from gurobipy import GRB

import UFL


def gen_UFL_instance(n, m):
    """Generates a test facility location instance.

    Args:
       n (int): number of facilities
       m (int): number of customers

    Returns:
        A tuple with the following values.

            - S (list): neighborhood list
            - f (dict): costs of facility location
              (generated uniformly random int from `f_min` to `f_max`)
            - g (dict): overlap costs, keys: (customer, overlap)
    """
    N = [i for i in range(1, n+1)]
    M = [j for j in range(1, m+1)]

    good_instance = False
    while not good_instance:
        S = [np.random.choice(N, np.random.randint(1, n+1),
                              replace=False).tolist()
             for j in M]
        f = {i: np.random.randint(5, 10) for i in N}

        g = dict()
        for j in M:
            g[(j, 0)] = np.random.randint(2*n, 4*n)
            g[(j, 1)] = 0
            for k in range(2, len(S[j-1])+1):
                g[(j, k)] = g[(j, k-1)] + \
                    np.random.randint(low=1, high=2+int(5*len(S[j-1])/j))

        Sf = UFL.build_Sf(S)
        if np.sum([1 for i in N if i not in Sf]) == 0:
            good_instance = True

    return S, f, g


@pytest.mark.parametrize("test_inst", [gen_UFL_instance(np.random.randint(1, 15),
                                                        np.random.randint(1, 15))
                                       for _ in range(100)])
def test_MIPs(test_inst):
    """Tests that MIPs return the same objectives."""
    S, f, g = test_inst
    print(f"S={S}; f={f}; g={g}")

    model = UFL.build_MIP(S, f, g)
    model.setParam("OutputFlag", 0)
    model.optimize()

    assert model.status == GRB.OPTIMAL
    plain_MIP_obj = model.objVal

    # Generate and solve CPP MIP
    C = UFL.create_covering_BDD_wg(S, g)
    A = UFL.create_availability_BDD(S, f)

    model, c, v, x = UFL.add_BDD_to_MIP(A, prefix="A_")
    model, c, v, x = UFL.add_BDD_to_MIP(C, model=model, x=x, prefix="C_")
    model.update()
    model.setParam("OutputFlag", 0)
    model.optimize()

    assert model.status == GRB.OPTIMAL

    CPP_MIP_obj = model.objVal
    assert abs(CPP_MIP_obj - plain_MIP_obj) < 1e-3


@pytest.mark.parametrize("test_inst", [gen_UFL_instance(np.random.randint(1, 15),
                                                        np.random.randint(1, 15))
                                       for _ in range(100)])
def test_DPs(test_inst):
    """Tests the DP-inspired diagrams."""
    S, f, g = test_inst
    print(f"S={S}; f={f}; g={g}")

    model = UFL.build_MIP(S, f, g)
    model.setParam("OutputFlag", 0)
    model.optimize()

    assert model.status == GRB.OPTIMAL
    plain_MIP_obj = model.objVal

    # Generate and solve CPP MIP
    D, _ = UFL.build_DP_DD(S, f, g)

    model, c, v = UFL.create_NF(D)
    model.update()
    model.setParam("OutputFlag", 0)
    model.optimize()

    assert model.status == GRB.OPTIMAL
    assert abs(model.objVal - plain_MIP_obj) < 1e-3
