"""Generates and solves FLP with ``soft cover'' constraints."""
import numpy as np
import gurobipy as gp
from gurobipy import GRB
from time import time
import argparse as ap
from tUFLP import DegreeKeeper
import BDD as DD
from copy import copy
from tUFLP import add_BDD_to_MIP
import pytest

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
        x[j+1] = m.addVar(vtype=GRB.BINARY, name=f"x_{j+1}",
                          obj=c[j])

        for a in range(1, len(S[j])+1):
            y[(j+1, a)] = m.addVar(vtype=GRB.BINARY, name=f"y_{j+1}_{a}",
                                   obj=f[j][a]-f[j][a-1])

    # create constraints
    for j in range(1, len(S)+1):
        m.addConstr(gp.quicksum(x[k] for k in S[j-1]) ==
                    gp.quicksum(y[(j, a)] for a in range(1, len(S[j-1])+1)))

        for a in range(1, len(S[j-1])):
            m.addConstr(y[(j, a)] >= y[(j, a+1)])

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
    S = [[j+1] for j in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if np.random.uniform() <= p:
                S[i].append(j+1)
                S[j].append(i+1)

    return S


def generate_overlaps(S):
    """Generates overlaps values (costs) `f` given adjacencies `S`.

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


def dump_instance(S, filename="tmp/S.dot"):
    """Dumps a graph implied by S into a `.dot` file. """
    added = set([])
    with open(filename, "w") as fout:
        fout.write("graph G {\n")
        for i in range(len(S)):
            for j in S[i]:
                if ((i+1) != j) and not (((j, (i+1)) in added)
                                         or ((i+1, j) in added)):
                    fout.write(f"n{i+1} -- n{j};\n")
                    added.add(((i+1), j))

        fout.write("}")

def make_instance(n, p=0.25, verbose=True):
    """Generates an instance of `n` points and connectivity `p`."""
    S = generate_S(n, p)
    f = generate_overlaps(S)
    Cmax = 5.0
    c = [Cmax*np.random.uniform() for _ in range(len(S))]
    if verbose:
        print(f"S={S}, f={f}, c={c}")
    return S, f, c


def make_string_inst(n, verbose=True):
    """Generates a special-type instance ('string')."""
    S = [[j+1] for j in range(n)]
    for j in range(1, n-1):
        S[j].append(j)
        S[j].append(j+2)

    S[-1].append(len(S)-1)
    S[0].append(2)

    f = generate_overlaps(S)
    Cmax = 5.0
    c = [Cmax*np.random.uniform() for _ in range(len(S))]
    if verbose:
        print(f"S={S}, f={f}, c={c}")
    return S, f, c


def make_organic_inst(n, verbose=True):
    """Generates a special-type instance ('organic molecule').

    Note: `n` here is length of the "main" line.
    """
    good_inst = False
    while not good_inst:
        groups = [np.random.randint(5) for _ in range(n-2)]
        N = 1
        while (N+1) + sum(groups[:(N+1)]) < n:
            N += 1

        if N>2:
            good_inst = True

    S = [[j+1] for j in range(N)]
    for j in range(1, N-1):
        S[j].append(j)
        S[j].append(j+2)

    S[-1].append(len(S)-1)
    S[0].append(2)

    newnode = N+1
    for j in range(1, N-1):
        for _ in range(groups[j-1]):
            S.append([newnode])
            S[j].append(newnode)
            S[-1].append(j+1)
            newnode += 1

    f = generate_overlaps(S)
    Cmax = 5.0
    c = [Cmax*np.random.uniform() for _ in range(len(S))]
    if verbose:
        print(f"S={S}, f={f}, c={c}")
    return S, f, c


def make_caveman_inst(n=10, M=5, pcave=0.8, verbose=True):
    """Generates a special type instance ('cavemen')."""
    S = [[j+1] for j in range(n)]
    for j in range(1, len(S)-1):
        S[j].append(j)
        S[j].append(j+2)

    S[-1].append(len(S)-1)
    S[0].append(2)

    ncaves = n // M
    for k in range(1, ncaves+1):
        for i in range((k-1)*M, min(k*M, n)):
            for j in range(i+2, min(k*M, n)):
                if np.random.uniform() <= pcave:
                    S[i].append(j+1)
                    S[j].append(i+1)

    f = generate_overlaps(S)
    Cmax = 5.0
    c = [Cmax*np.random.uniform() for _ in range(len(S))]
    if verbose:
        print(f"S={S};\nf={f}\n;c={c}")
    return S, f, c


def assert_instance(S, f, c):
    """Check that the instance is technically correct."""
    assert min(sum(S, [])) == 1
    assert max(sum(S, [])) == len(S)

    assert len(S) == len(f)
    assert len(S) == len(c)
    for Sj in S:
        assert len(Sj) == len(set(Sj))  # no repeats within each adj list
    for j in range(len(S)):
        assert (j+1) in S[j]

        for k in S[j]:
            assert (j+1) in S[k-1]

        assert len(f[j]) == len(S[j])+1

# generating a Cover BDD for the problem ######################################


def make_label(state, pts):
    """Helper: formats a node label out of the state."""
    return "\n".join([f"{pts[j]}:{state[j]}" for j in range(len(state))])


def mknode(state, freedoms):
    """Makes a node key.

    Args:
        state (list): no. of overlaps per point,
        freedoms (DegreeKeeper): no. of point (residual) degrees.

    Returns:
        String label.
    """
    return ";".join([f"{j+1}:{freedoms[j+1]}" for j in range(len(state))
                     if state[j] > 0])


def build_soft_cover_DD(S, f, c, next_node_type='min'):
    """Builds a BDD for the soft-UFL problem with overlaps.

    Args:
        S (list): neighborhood list,
        f (list): overlap costs, `f[j][a]`, a=0,..,|S_j|
        c (list): location costs.
            next_node_type (str): `min`, `max`, or `rnd` -- see `DegreeKeeper`
            docstring for details.

    Returns:
        The resulting BDD.

    Notes:
        - employs a DP approach with state being the Boolean
          covering at each customer.

    """
    N = len(S)
    B = DD.BDD(N=N, vars=[f"stub{i}" for i in range(1, N+1)],
               weighted=True)

    freedoms = DegreeKeeper(S, next_node_type)  # node degrees

    root_state = tuple(0 for _ in range(len(S)))
    pts = [j for j in range(1, len(S)+1) if freedoms.has_freedom(j)]
    assert len(pts) == len(root_state)
    node_labels = dict({DD.NROOT: make_label(root_state, pts)})

    next_layer = {tuple(root_state): B.addnode(None)}

    k = 0  # created layers counter

    while k < N-1:
        # we take node `i` and try to add it and all its neighbors
        # to the diagram (unless they are already added)
        i = freedoms.get_next()  # current 'central' node to process
        for j in S[i-1]:
            if f"x{j}" in B.vars or k == N-1:
                continue

            # adding j-th point (of the original graph)
            current_layer = copy(next_layer)
            next_layer = dict()

            for q in S[j-1]:
                if freedoms.has_freedom(q):
                    freedoms.decrement(q)

            for state in current_layer:
                node = current_layer[state]

                next_state = tuple([state[k] for k in range(len(state))
                                    if freedoms.has_freedom(pts[k])])

                if next_state in next_layer:
                    B.llink(node, next_layer[next_state], "lo")
                else:
                    newnode = B.addnode(node, "lo")
                    next_layer.update({next_state: newnode})
                    node_labels.update({newnode.id: make_label(next_state,
                                                               pts)})

                next_state = tuple(state[k] + (pts[k] in S[j-1])
                                   for k in range(len(state))
                                   if freedoms.has_freedom(pts[k]))

                a = {pts[j]: state[j] for j in range(len(pts))}
                overlap_cost = 0

                for q in S[j-1]:
                    overlap_cost += f[q - 1][a[q]+1] - f[q - 1][a[q]]

                if next_state in next_layer:
                    B.llink(node, next_layer[next_state],
                            "hi", edge_weight=c[j-1] + overlap_cost)
                else:
                    newnode = B.addnode(node, "hi",
                                        edge_weight=c[j-1]+overlap_cost)
                    next_layer.update({next_state: newnode})
                    node_labels.update({newnode.id: make_label(next_state,
                                                               pts)})

            ptsp = [j for j in pts if freedoms.has_freedom(j)]
            pts = ptsp

            B.rename_vars({f"stub{k+1}": f"x{j}"})
            k += 1

    # process the last node in S[i-1] separately
    i = -1

    while ((f"x{i}" in B.vars) or (i == -1)):
        i, _ = freedoms.pop()


    current_layer = copy(next_layer)

    const_cost = sum(fj[0] for fj in f)

    for state in current_layer:
        node = current_layer[state]
        B.link(node.id, DD.NTRUE, "lo", const_cost)

        a = {pts[j]: state[j] for j in range(len(pts))}
        overlap_cost = 0

        for q in S[i-1]:
            overlap_cost += f[q - 1][a[q]+1] - f[q - 1][a[q]]

        B.link(node.id, DD.NTRUE, "hi", c[i-1] +
               overlap_cost + const_cost)

    B.rename_vars({f"stub{N}": f"x{i}"})

    return B, node_labels


# experiments code ############################################################
def dia_sizes(n1=5, n2=10, K=5, igen=make_instance):
    """Prints random diagram sizes for different instance size."""
    print("experiment, N, Bsize_nonred, Bsize_reduced, naive_MIP_vars, t_MIP, t_BDD")
    expn = 1
    for n in range (n1, n2+1):
        for _ in range(K):
            S, f, c = igen(n, verbose=False)
            assert_instance(S, f, c)
            # naive MIP approach
            t0 = time()
            m, x, y = make_MIP(S, c, f)
            m.optimize()
            t1 = time()
            t_MIP = t1-t0
            assert m.status == GRB.OPTIMAL
            aobj = m.objVal + sum(fs[0] for fs in f)
            nvars = len(x.keys())+len(y.keys())

            # BDD approach
            t0 = time()
            B, _ = build_soft_cover_DD(S, f, c, next_node_type='min')
            B_as_is = B.size()
            B.make_reduced()
            B_red = B.size()
            m, _, _, x = add_BDD_to_MIP(B)
            m.update()
            m.optimize()
            t1 = time()
            t_BDD = t1-t0
            assert m.status == GRB.OPTIMAL
            assert (abs(aobj - m.objVal) < 1e-2)
            print(
                f"{expn}, {n}, {B_as_is}, {B_red}, {nvars}, {t_MIP:.2f}, {t_BDD:.2f}",
                flush=True)
            expn += 1

# testing the code ############################################################

def try_softcover_inst(S, c, f):
    """Solves a softcover inst with naive MIP and BDD MIP (both w/Gurobi)."""
    assert_instance(S, f, c)
    m, x, y = make_MIP(S, c, f)
    m.optimize()

    assert m.status == GRB.OPTIMAL
    aobj = m.objVal + sum(fs[0] for fs in f)
    print(f"(Adjusted) Objective is: {aobj}")
    print(f"Decisions are: {[x[j].x for j in range(1, len(S)+1)]}")

    x1 = [x[j].x for j in range(1, len(S)+1)]
    print("Alternative approach: building a BDD...")
    B, _ = build_soft_cover_DD(S, f, c)
    B.make_reduced()
    print(f"BDD built, {B.size()} nodes.")
    m, _, _, x = add_BDD_to_MIP(B)
    print("Added to the model.")
    m.update()
    m.optimize()
    print("Model solved.")
    assert m.status == GRB.OPTIMAL
    print(f"LP objective is: {m.objVal}")
    print(f"Decisions are: {[x[xvar].x for xvar in x]}")
    x2 = [x[xvar].x for xvar in x]

    return (abs(aobj - m.objVal) < 1e-2)


def test_build_soft_cover_DD_simple1():
    """Tests the softcover DD construction procedure."""
    S = [[1, 2], [1, 2, 3], [2, 3]]
    c = [2, 5, 2]
    f = [[10, 0, 0], [10, 0, 0, 0], [10, 0, 0]]

    assert try_softcover_inst(S, c, f)


def test_build_soft_cover_DD_simple2():
    """Tests the softcover DD construction procedure."""
    S = [[1, 2, 4, 5, 6],
         [2, 1, 3],
         [3, 2, 7, 8, 9],
         [4, 1], [5, 1], [6, 1],
         [7, 3], [8, 3], [9, 3]]

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

    assert try_softcover_inst(S, c, f)


@pytest.mark.parametrize("test_inst", [make_instance(np.random.randint(5, 10),
                                                     verbose=False)
                                       for _ in range(1000)])
def test_build_soft_cover_DD(test_inst):
    """Asserts naive MIP ~ BDD MIP."""
    S, f, c = test_inst
    print(f"S={S}; c={c}; f={f}")
    assert try_softcover_inst(S, c, f)


def test_make_MIP():
    """Tries a trivial instance."""
    S = [[1, 2], [1, 2, 3], [2, 3]]
    c = [2, 5, 2]
    f = [[10, 0, 0], [10, 0, 0, 0], [10, 0, 0]]

    m, x, y = make_MIP(S, c, f)
    m.optimize()

    assert m.status == GRB.OPTIMAL
    aobj = m.objVal + sum(fs[0] for fs in f)
    print(f"(Adjusted) Objective is: {aobj}")
    print(f"Decisions are: {[x[j].x for j in range(1, len(S)+1)]}")
    assert abs(aobj - 4.0)<1e-5

    f = [[10, 0, 0], [10, 0, 10, 0], [10, 0, 0]]
    m, x, y = make_MIP(S, c, f)
    m.optimize()

    assert m.status == GRB.OPTIMAL
    aobj = m.objVal + sum(fs[0] for fs in f)
    assert abs(aobj - 5.0) < 1e-5
    print(f"(Adjusted) Objective is: {m.objVal + sum(fs[0] for fs in f)}")
    print(f"Decisions are: {[x[j].x for j in range(1, len(S)+1)]}")


def test_MIP_example():
    """Tries a specific, simple instance."""
    S = [[1, 4, 5, 6],
         [2, 1, 3],
         [3, 2, 7, 8, 9],
         [4, 1], [5, 1], [6, 1],
         [7, 3], [8, 3], [9, 3]]

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
    print(f"Decisions are: {[x[j].x for j in range(1, len(S)+1)]}")
    assert abs(aobj - 7.0) < 1e-5

    print("Updating f to screen overlaps at 0 and 2")
    f[0] = [10, 0, 10, 10, 10, 10, 10]
    f[2] = [10, 0, 10, 10, 10, 10, 10]
    m, x, y = make_MIP(S, c, f)
    m.optimize()

    assert m.status == GRB.OPTIMAL
    aobj = m.objVal + sum(fs[0] for fs in f)
    print(f"(Adjusted) Objective is: {aobj}")
    print(f"Decisions are: {[x[j].x for j in range(1, len(S)+1)]}")
    assert abs(aobj - 8.0) < 1e-5

    print("Updating f to screen overlap at 1")
    f[1][2] = 10
    m, x, y = make_MIP(S, c, f)
    m.optimize()

    assert m.status == GRB.OPTIMAL
    aobj = m.objVal + sum(fs[0] for fs in f)
    print(f"(Adjusted) Objective is: {aobj}")
    print(f"Decisions are: {[x[j].x for j in range(1, len(S)+1)]}")


# main run (if module is started by itself) ###################################
if __name__ == '__main__':
    parser = ap.ArgumentParser(description="''Soft overlap'' experiments. (c) A. Bochkarev, 2022",
                               formatter_class=ap.ArgumentDefaultsHelpFormatter)

    parser.add_argument("EXPERIMENT", action="store",
                        help="experiment name (rnd_inst, dia_sizes)")
    parser.add_argument('-N', dest="N", action="store", default="-1",
                        help="number of points in the network")
    args = parser.parse_args()

    if args.EXPERIMENT == "rnd_inst":
        if args.N == "-1":
            print("N is required!")
            exit(1)

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
    elif args.EXPERIMENT == "dia_sizes":
        dia_sizes(5, 20, K=5)
    elif args.EXPERIMENT == "dia_sizes_string":
        dia_sizes(5, 15, K=1, igen=make_string_inst)
    elif args.EXPERIMENT == "dia_sizes_organic":
        dia_sizes(5, 25, K=1, igen=make_organic_inst)
    elif args.EXPERIMENT == "dia_sizes_caves":
        dia_sizes(25, 175, K=5, igen=make_caveman_inst)
    elif args.EXPERIMENT == "MIP_time_cavemen":
        if args.N == "-1":
            print("N is required!")
            exit(1)

        S, f, c = make_caveman_inst(int(args.N))
        m, x, y = make_MIP(S, c, f)
        t0 = time()
        m.optimize()
        t1 = time()
        assert m.status == GRB.OPTIMAL
        aobj = m.objVal + sum(fs[0] for fs in f)
        print(f"(Adjusted) Objective is: {aobj}")
        print(f"Decisions are: {[x[j].x for j in range(1, len(S)+1)]}")
        print(f"Finished in {t1 - t0:.2f} sec.")
