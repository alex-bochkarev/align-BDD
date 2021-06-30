"""Joint UFL example: another proof-of-concept
for aligning BDDs.

Tests coverage: :py:mod:`jUFL_test`.
"""

import pytest
import tUFLP
import BDD as DD
import numpy as np
import varseq as vs
import BB_search as bb
import gurobipy as gp
from gurobipy import GRB

def solve_with_DD_MIP(inst):
    """Solves a jUFL instance with Gurobi.

    Returns:
        m (gurobipy Model): the model after solution,
        x (gurobipy var): location-decision varaibles.
    """
    S_1, S_2, f_1, f_2 = inst
    cov_1, _ = tUFLP.build_cover_DD(S_1, f_1)
    cov_2, _ = tUFLP.build_cover_DD(S_2, f_2)

    m, c, v, x = tUFLP.add_BDD_to_MIP(cov_1)
    m, c, v, x = tUFLP.add_BDD_to_MIP(cov_2, m, x)

    m.optimize()
    return m, x


def solve_with_naive_MIP(instance):
    """Solves `instance` with a naive MIP."""
    S_1, S_2, f_1, f_2 = instance
    m = gp.Model()
    m.modelSense = GRB.MINIMIZE

    x = dict()

    V_1 = range(1, len(S_1)+1)
    V_2 = range(1, len(S_2)+1)

    for i in V_1:
        x[i] = m.addVar(vtype=GRB.BINARY, name=f"x_{i}", obj=f_1[i] + f_2[i])

    # covering constraints
    for j in V_1:
        m.addConstr(gp.quicksum(x[i] for i in S_1[j-1]) >= 1)

    for j in V_2:
        m.addConstr(gp.quicksum(x[i] for i in S_2[j-1]) >= 1)

    m.update()
    m.optimize()
    return m, x


def generate_instance(n, p=0.25):
    """Generates a joint-UFL instance of size `n`."""
    status = GRB.INFEASIBLE

    while status != GRB.OPTIMAL:
        S_1, f_1, _, _ = tUFLP.generate_test_instance(n, p=p)
        S_2, f_2, _, _ = tUFLP.generate_test_instance(n, p=p)
        m, _ = solve_with_naive_MIP((S_1, S_2, f_1, f_2))
        status = m.status

    return S_1, S_2, f_1, f_2


def solve_with_align_BDD(instance):
    """Solves the problem by generating two BDDs, aligning them, and solving a NF.
    """
    S_1, S_2, f_1, f_2 = instance
    C_1, _ = tUFLP.build_randomized_cover_DD(S_1, f_1)
    C_2, _ = tUFLP.build_randomized_cover_DD(S_2, f_2)

    C_1.make_reduced()
    C_2.make_reduced()

    vs_1 = vs.VarSeq(C_1.vars, [len(L) for L in C_1.layers[:-1]])
    vs_2 = vs.VarSeq(C_2.vars, [len(L) for L in C_2.layers[:-1]])

    assert set(vs_1.layer_var) == set(vs_2.layer_var), f"A:{vs_1.layer_var}, C:{vs_2.layer_var}"
    b = bb.BBSearch(vs_1, vs_2)

    # bb.TIMEOUT_ITERATIONS=10000
    status = b.search()
    assert status == "optimal" or status == "timeout"

    C_1.align_to(b.Ap_cand.layer_var, inplace=True)
    C_2.align_to(b.Ap_cand.layer_var, inplace=True)

    int_DD = DD.intersect(C_1, C_2)
    int_DD.make_reduced()

    nl = int_DD.shortest_path()
    return nl[DD.NROOT]
