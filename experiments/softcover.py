"""Generates and solves facility location instances
with ``soft cover'' constraints.
"""
import numpy as np
import gurobipy as gp
from gurobipy import GRB
from time import time
import argparse as ap

def make_MIP(S, c, f):
    """Creates a MIP model (for gurobi) from an instance specs.

    Args:
        S (list): list of adjacency lists,
        c (list): location costs per point,
        f (list): overlap cost function, f[j][a]

    Returns:
        (m, x, y): gurobi model and vars
    """
    m = gp.Model()
    m.modelSense = GRB.MINIMIZE
    m.setParam("OutputFlag", 0)
    x = dict()
    y = dict()

    # create variables
    for j in range(len(S)):
        x[j] = m.addVar(vtype=GRB.BINARY, name=f"x_{j}",
                        obj=c[j])
        for a in range(1, len(S[j])+1):
            y[(j, a)] = m.addVar(vtype=GRB.BINARY, name=f"y_{j}_{a}",
                                 obj=f[j][a]-f[j][a-1])

    # create constraints
    for j in range(len(S)):
        m.addConstr(gp.quicksum(x[k] for k in S[j]) ==
                    gp.quicksum(y[(j, a)] for a in range(1, len(S[j])+1)))
        for a in range(1, len(S[j])):
            m.addConstr(y[(j,a)] >= y[(j, a+1)])

    m.update()
    return m, x, y

def generate_S(n, p=0.25):
    """Generates adjacency lists `S`.

    Args:
        n (int): number of points,
        p (float): probability of an edge.

    Returns:
        S: list of lists.
    """
    S = [[j] for j in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if np.random.uniform() <= p:
                S[i].append(j)
                S[j].append(i)

    return S


def generate_overlaps(S):
    """Generates overlaps values `f`.

    Args:
        S (list): list of lists (adjacencies)

    Returns:
        f (list): list of lists, f[j][a], a=0,..,|S_j|
    """
    Fmax = 10.0
    f = [[Fmax, 0.0] for _ in range(len(S))]
    for j in range(len(S)):
        for a in range(1, len(S[j])):
            f[j].append(Fmax * (np.random.uniform()-0.5))

    return f

def make_instance(n, p=0.25):
    """Generates an instance of `n` points and connectivity `p`."""
    S = generate_S(n,p)
    f = generate_overlaps(S)
    Cmax = 5.0
    c = [Cmax*np.random.uniform() for _ in range(len(S))]
    print(f"S={S}, f={f}, c={c}")
    return S, f, c

# testing the code ############################################################
def test_make_MIP():
    """Tries a trivial instance."""
    S = [[0, 1], [0, 1, 2], [1, 2]]
    c = [2, 5, 2]
    f = [[10, 0, 0], [10, 0, 0, 0], [10, 0, 0]]

    m, x, y = make_MIP(S, c, f)
    m.optimize()

    assert m.status == GRB.OPTIMAL
    aobj = m.objVal + sum(fs[0] for fs in f)
    print(f"(Adjusted) Objective is: {aobj}")
    print(f"Decisions are: {[x[j].x for j in range(len(S))]}")
    assert abs(aobj - 4.0)<1e-5

    f = [[10, 0, 0], [10, 0, 10, 0], [10, 0, 0]]
    m, x, y = make_MIP(S, c, f)
    m.optimize()

    assert m.status == GRB.OPTIMAL
    aobj = m.objVal + sum(fs[0] for fs in f)
    assert abs(aobj - 5.0) < 1e-5
    print(f"(Adjusted) Objective is: {m.objVal + sum(fs[0] for fs in f)}")
    print(f"Decisions are: {[x[j].x for j in range(len(S))]}")

def test_MIP_example():
    """Tries a specific, simple instance."""
    S = [[0, 3, 4, 5],
         [1, 0, 2],
         [2, 1, 6, 7, 8],
         [3, 0], [4, 0], [5, 0],
         [6, 2], [7, 2], [8, 2]]

    c = [4, 1, 4, 1, 1, 1, 1, 1, 1]

    f = [[10, 0, 0, 0, 0, 0],
         [10, 0, 0, 0],
         [10, 0, 0, 0, 0, 0],
         [10, 0, 0],
         [10, 0, 0],
         [10, 0, 0],
         [10, 0, 0],
         [10, 0, 0],
         [10, 0, 0]]

    m, x, y = make_MIP(S, c, f)
    m.optimize()

    assert m.status == GRB.OPTIMAL
    aobj = m.objVal + sum(fs[0] for fs in f)
    print(f"(Adjusted) Objective is: {aobj}")
    print(f"Decisions are: {[x[j].x for j in range(len(S))]}")
    assert abs(aobj - 7.0) < 1e-5

    print("Updating f to screen overlaps at 0 and 2")
    f[0] = [10, 0, 10, 10, 10, 10, 10]
    f[2] = [10, 0, 10, 10, 10, 10, 10]
    m, x, y = make_MIP(S, c, f)
    m.optimize()

    assert m.status == GRB.OPTIMAL
    aobj = m.objVal + sum(fs[0] for fs in f)
    print(f"(Adjusted) Objective is: {aobj}")
    print(f"Decisions are: {[x[j].x for j in range(len(S))]}")
    assert abs(aobj - 8.0) < 1e-5

    print("Updating f to screen overlap at 1")
    f[1][2] = 10
    m, x, y = make_MIP(S, c, f)
    m.optimize()

    assert m.status == GRB.OPTIMAL
    aobj = m.objVal + sum(fs[0] for fs in f)
    print(f"(Adjusted) Objective is: {aobj}")
    print(f"Decisions are: {[x[j].x for j in range(len(S))]}")

# main run (if module is started by itself) ###################################
if __name__ == '__main__':
    parser = ap.ArgumentParser(description="''Soft overlap'' instance generator. (c) A. Bochkarev, 2022",
                               formatter_class=ap.ArgumentDefaultsHelpFormatter)

    parser.add_argument("N", action="store",
                        help="number of points in the network")
    args = parser.parse_args()

    S, f, c = make_instance(int(args.N))
    m, x, y = make_MIP(S, c, f)
    t0 = time()
    m.optimize()
    t1 = time()
    assert m.status == GRB.OPTIMAL
    aobj = m.objVal + sum(fs[0] for fs in f)
    print(f"(Adjusted) Objective is: {aobj}")
    print(f"Decisions are: {[x[j].x for j in range(len(S))]}")
    print(f"Finished in {t1 - t0:.2f} sec.")
