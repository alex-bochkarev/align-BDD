"""Performs the sum-of-products linearize/solve exercise."""

import numpy as np
from math import comb
from time import time
import gurobipy as gp
from gurobipy import GRB
import argparse as ap

def gen_instance(n, M=3, jmax=3, Cmax=10.0):
    """Generates an instance.

    Args:
        n (int): number of variables
        M (int): number of terms in the sum
        jmax (int): number of vars in a term, at most
        Cmax (float): cost coefficient will be -Cmax..Cmax

    Returns:
        I (list): list of lists of indices
        c (list): list of costs (per sum term)
    """
    I = []
    V = [i for i in range(n)]
    terms = {()}
    assert comb(n, jmax) > M

    for m in range(M):
        cand = ()
        while tuple(cand) in terms:
            cand = np.sort(np.random.choice(V, size=np.random.randint(1, jmax+1),
                                            replace=False))
        I.append(list(cand))
        terms.add(tuple(cand))


    c = [(np.random.uniform()-0.5)*Cmax*2 for _ in range(len(I))]
    print("Instance generated:")
    print("min ", end="")
    for j in range(len(I)):
        print(f" {c[j]:+.2f}·"+"·".join([f"x{k}" for k in I[j]]), end="")
    print()
    return I, c


def make_MIP(I, c):
    """Makes a MIP model for the problem, with Gurobi.

    Args:
        I (list): list of indices per sum term.
        c (list): corresponding cost coefficients.

    Returns:
        m, x, y: the resulting model and variables.
    """
    assert len(I) == len(c)
    m = gp.Model()
    m.modelSense = GRB.MINIMIZE

    x = {j: m.addVar(vtype=GRB.BINARY, name=f"x_{j}")
         for j in np.unique(sum(I, []))}
    y = {j: m.addVar(vtype=GRB.BINARY, name=f"y_{j}")
         for j in range(len(I)) if len(I[j])>0}

    obj = gp.LinExpr(0.0)

    for j, term in enumerate(I):
        if len(term) > 1:
            m.addConstrs((y[j] <= x[k] for k in term), name=f"term_{j}")
            m.addConstr(y[j] >= gp.quicksum(x[k] for k in term) - len(term)+1)
            obj += c[j] * y[j]
        else:
            obj += c[j] * x[term[0]]

    m.setObjective(obj)
    m.setParam("OutputFlag", 0)
    m.update()
    return m, x, y

# main run (if module is started by itself) ###################################
if __name__ == '__main__':
    parser = ap.ArgumentParser(description="''Sum-of-products'' instance generator. (c) A. Bochkarev, 2022",
                               formatter_class=ap.ArgumentDefaultsHelpFormatter)

    parser.add_argument("N", action="store",
                        help="number of vars")
    parser.add_argument("M", action="store",
                        help="number of terms")
    parser.add_argument("J", action="store",
                        help="number of vars per term (max)")
    args = parser.parse_args()

    I, c = gen_instance(int(args.N), int(args.M), int(args.J))
    m, x, y = make_MIP(I, c)
    t = time()
    m.optimize()

    assert m.status == GRB.OPTIMAL
    print(f"Objective is: {m.objVal}")
    print(f"Decisions are: x={[x[j].x for j in x.keys()]}, y={[y[j].x for j in y.keys()]}")
    print(f"Finished in {time() - t} sec.")
