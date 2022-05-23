"""A specialized code for typed UFLP over cavemen instances."""
import pytest
from copy import deepcopy
import numpy as np
import gurobipy as gp

from darkcloud import ptscloud, generate_overlaps, DDSolver
from BDD import intersect
from varseq import VarSeq
from BB_search import BBSearch

from time import time


def gen_caveman_inst(n=10, M=7, L=0.25, verbose=False):
    """Generates an instance with the related metadata (info on caves).

    Args:
      n (int): number of caves,
      M (int): number of points in a cave,
      L (float): edge sparsity parameter (share of missing edges)
      verbose (Bool): print debug info

    Returns:
      S, S2, f, f2, c, caves: instance and caves description.

    Note:
      The implementation is based on :py:func:`darkcloud.gen_caveman_inst`,
      but with the 'types' condition replaced with another cover-like
      condition. Note that caves are the same for the two networks:
      sets of nodes coincide, but the connection within each cave
      may be different.
    """
    # creating the first graph topology
    S = [[j+1] for j in range(n * M)]  # creating the nodes first
    entry_point = (None, None)
    caves = []
    for k in range(n):
        lastcave = [k*M + 1]
        if len(caves) > 0:
            caves[-1].e2 = entry_point
            S[entry_point[0]-1].append(entry_point[1])
            S[entry_point[1]-1].append(entry_point[0])

        # create the necessary number of connected nodes
        while len(lastcave) < M:
            connect_to = np.random.choice(lastcave)
            lastcave.append(lastcave[-1]+1)
            S[lastcave[-1]-1].append(connect_to)
            S[connect_to-1].append(lastcave[-1])

        caves.append(ptscloud(entry_point,
                              (None, None),
                              [j for j in range(k*M+1, (k+1)*M+1)]))

        entry_point = (np.random.choice(lastcave), (k+1)*M+1)

    S2 = deepcopy(S)

    # add edges to the first graph
    for cave in caves:
        n_edges = 0
        for i in cave.S:
            for j in S[i - 1]:
                if (i != j) and ((i, j) not in [cave.e1, cave.e2]) and (
                        (j, i) not in [cave.e1, cave.e2]):
                    n_edges += 1

        n_edges /= 2

        while (1 - 2*n_edges / (M*(M-1))) > L:
            n1 = np.random.choice(cave.S)
            n2 = np.random.choice(cave.S)

            if n1 not in S[n2-1]:
                S[n2-1].append(n1)
                S[n1-1].append(n2)
                n_edges += 1

    # add edges to the second graph
    caves2 = deepcopy(caves)

    for cave in caves2:
        n_edges = 0
        for i in cave.S:
            for j in S2[i - 1]:
                if (i != j) and ((i, j) not in [cave.e1, cave.e2]) and (
                        (j, i) not in [cave.e1, cave.e2]):
                    n_edges += 1

        n_edges /= 2

        while (1 - 2*n_edges / (M*(M-1))) > L:
            n1 = np.random.choice(cave.S)
            n2 = np.random.choice(cave.S)

            if n1 not in S2[n2-1]:
                S2[n2-1].append(n1)
                S2[n1-1].append(n2)
                n_edges += 1

    # creating costs info (location and overlap costs)
    f = generate_overlaps(S)
    f2 = generate_overlaps(S2)

    Cmax = 5.0
    c = [Cmax*np.random.uniform() for _ in range(len(S))]
    if verbose:
        print(f"S={S};\nf={f}\n;c={c}")
    return S, S2, f, f2, c, caves


def dump_instance(S, caves, filename="tmp/S.dot"):
    """Dumps a graph implied by S into a `.dot` file. """
    added = set([])
    with open(filename, "w") as fout:
        fout.write("graph G {\n")
        for i in range(len(S)):
            for j in S[i]:
                if ((i+1) != j) and not (((j, (i+1)) in added)
                                         or ((i+1, j) in added)):
                    fout.write(f"    n{i+1} -- n{j};\n")
                    added.add(((i+1), j))

        fout.write("}")


def solve_with_MIP(S, S2, f, f2, c, caves):
    """Generates and solves a MIP for the supplied jUFL cavemen instance.

    Note:
      The instance is parameterized as per :py:func:`gen_caveman_inst`.
      Most of the code is adapted form :py:func:`darkcloud.solve_with_MIP`.
    """
    m = gp.Model()
    m.modelSense = gp.GRB.MINIMIZE
    m.setParam("OutputFlag", 0)
    x = dict()
    y = dict()
    y2 = dict()

    # create variables
    for j in range(1, len(S)+1):
        x[j] = m.addVar(vtype=gp.GRB.BINARY, name=f"x_{j}",
                        obj=c[j-1])

        for a in range(1, len(S[j-1])+1):
            y[(j, a)] = m.addVar(vtype=gp.GRB.BINARY, name=f"y_{j}_{a}",
                                 obj=f[j-1][a]-f[j-1][a-1])

        for a in range(1, len(S2[j-1])+1):
            y2[(j, a)] = m.addVar(vtype=gp.GRB.BINARY, name=f"y2_{j}_{a}",
                                  obj=f2[j-1][a]-f2[j-1][a-1])

    # create constraints
    for j in range(1, len(S)+1):
        # first-graph constraints
        m.addConstr(gp.quicksum(x[k] for k in S[j-1]) ==
                    gp.quicksum(y[(j, a)] for a in range(1, len(S[j-1])+1)))

        for a in range(1, len(S[j-1])):
            m.addConstr(y[(j, a)] >= y[(j, a+1)])

        # second-graph constraints
        m.addConstr(gp.quicksum(x[k] for k in S2[j-1]) ==
                    gp.quicksum(y2[(j, a)] for a in range(1, len(S2[j-1])+1)))

        for a in range(1, len(S2[j-1])):
            m.addConstr(y2[(j, a)] >= y2[(j, a+1)])

    m.update()
    print("Optimizing the model...", end="", flush=True)
    m.optimize()
    print("done.")
    assert m.status == gp.GRB.OPTIMAL
    return m, (m.objVal + sum(fs[0] for fs in f) +
               sum(fs2[0] for fs2 in f2)), x, y, y2


def solve_with_DDs(S, S2, f, f2, c, caves,
                   intmode="toA"):
    """Solves the jUFLP cavemen instance with DDs.

    Args:
      intmode (str): intersection mode, 'toA' or 'VS'

    Notes:
      Instance is parameterized as per :py:func:`gen_caveman_inst`,
      The diagrams are built with :py:class:`darkcloud.DDSolver`.
    """
    print("Building the first DD...", end="", flush=True)
    sol = DDSolver(S, f, c, caves)
    B1 = sol.build_cover_DD()
    print("done.")

    print("Building the second DD...", end="", flush=True)
    sol = DDSolver(S2, f2, [0.0 for _ in c], caves)
    B2 = sol.build_cover_DD()
    print("done.")

    if intmode == "toA":
        target = B1.vars
    elif intmode == 'VS':
        vs1 = VarSeq(B1.vars, [len(L) for L in B1.layers[:-1]])
        vs2 = VarSeq(B2.vars, [len(L) for L in B2.layers[:-1]])

        assert set(vs1.layer_var) == set(
            vs2.layer_var), f"1:{vs1.layer_var}, 2:{vs2.layer_var}"

        b = BBSearch(vs1, vs2)
        status = b.search()
        assert status == "optimal" or status == "timeout"
        target = b.Ap_cand.layer_var
    else:
        print(f"Wrong intersection mode '{intmode}'. Expected: 'toA' or 'VS'.")

    B1.align_to(target, inplace=True)
    B2.align_to(target, inplace=True)

    int_DD = intersect(B1, B2)
    sp = int_DD.shortest_path()
    return sp[0]

# Testing code ######################################################
@pytest.mark.parametrize("test_inst", [gen_caveman_inst()
                                       for _ in range(5)])
def test_jUFL_DDs(test_inst):
    S, S2, f, f2, c, caves = test_inst
    obj1 = solve_with_DDs(S,S2,f,f2,c,caves, intmode='toA')
    obj2 = solve_with_DDs(S,S2,f,f2,c,caves, intmode='VS')
    assert abs(obj1 - obj2) < 0.001

@pytest.mark.parametrize("test_inst", [gen_caveman_inst()
                                       for _ in range(5)])
def test_jUFL_DDvsMIP(test_inst):
    S, S2, f, f2, c, caves = test_inst
    obj1 = solve_with_DDs(S,S2,f,f2,c,caves, intmode='VS')
    m, obj2, _, _, _ = solve_with_MIP(S, S2, f, f2, c, caves)
    assert m.status == gp.GRB.OPTIMAL
    assert abs(obj1 - obj2) < 0.001

def compare_runtimes():
    """Performs a quick runtimes comparison (toA vs VS)."""
    s,s2,f,f2,c,caves = gen_caveman_inst()
    sol = DDSolver(s, f, c, caves)
    B1 = sol.build_cover_DD()

    sol = DDSolver(s2, f2, [0.0 for _ in c], caves)
    B2 = sol.build_cover_DD()
    B2.shuffle_vars()  # renames the variables, without reordering

    # B1.make_reduced()
    # B2.make_reduced()

    t0 = time()
    B1p = B1.align_to(B2.vars, inplace=False)
    int_toA = intersect(B1p, B2)
    obj_toA = int_toA.shortest_path()[0]
    t_toA = time() - t0

    print(f"{t_toA:.2f}, {int_toA.size()},", end="", flush=True)

    t0 = time()
    vs1 = VarSeq(B1.vars, [len(L) for L in B1.layers[:-1]])
    vs2 = VarSeq(B2.vars, [len(L) for L in B2.layers[:-1]])

    assert set(vs1.layer_var) == set(
        vs2.layer_var), f"1:{vs1.layer_var}, 2:{vs2.layer_var}"

    b = BBSearch(vs1, vs2)
    status = b.search()
    assert status == "optimal" or status == "timeout"
    target = b.Ap_cand.layer_var
    B1p = B1.align_to(target, inplace=False)
    B2p = B2.align_to(target, inplace=False)

    int_VS = intersect(B1p, B2p)
    obj_VS = int_VS.shortest_path()[0]
    t_VS = time() - t0
    print(f"{t_VS:.2f}, {int_VS.size()}")
    assert (abs(obj_VS - obj_toA) < 0.01), f"{obj_VS:.2f} vs {obj_toA:.2f}"

def main():
    print("t_toA, intsize_toA, t_VS, intsize_VS")
    for _ in range(250):
        compare_runtimes()

if __name__ == '__main__':
   main()
