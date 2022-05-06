"""
Implements the 'cloud' algorithm that builds a BDD for a UFL over a relaxed
cavemen graph.

---
(c) A. Bochkarev, Clemson University, 2022
a@bochkarev.io
"""
import numpy as np
from dataclasses import dataclass
import gurobipy as gp
from experiments.softcover import generate_overlaps, assert_instance
from BDD import BDD, NTRUE
import pytest


@dataclass
class ptscloud:
    """Keeps current 'cloud' of points.

    Attributes:
        e1, e2 (tuple): the two external entry points of the cloud (edges),
        S (list[int]): a list of points within the cloud.
    """
    e1: tuple[int]
    e2: tuple[int]
    S: list[int]


def gen_caveman_inst(n=10, M=5, L=0.5, verbose=False):
    """Generates an instance with the related metadata (info on caves).

    Args:
      n (int): number of caves,
      M (int): number of points in a cave,
      L (float): edge sparsity parameter (share of missing edges)
      verbose (Bool): print debug info

    Returns:
      S, f, c, caves: instance (S,f,c) and caves description.

    Note:
      The parameter ``L`` is calculated as follows.

      .. math::
          \\Lambda = 1 - 2\\frac{#\\textrm{existing arcs}}{N(N-1)}

      (See, e.g., Sefair and Smith 2016 on DSPI for a similar approach.)
    """
    # creating the graph topology first
    S = [[j+1] for j in range(n * M)]  # creating the nodes first
    entry_point = (None, None)
    caves = []
    for k in range(n):
        lastcave = [k*M + 1]
        if len(caves) > 0:
            caves[-1].e2 = entry_point
            S[entry_point[0]-1].append(entry_point[1])
            S[entry_point[1]-1].append(entry_point[0])

        n_edges = 0
        # create the necessary number of connected nodes
        while len(lastcave) < M:
            connect_to = np.random.choice(lastcave)
            lastcave.append(lastcave[-1]+1)
            S[lastcave[-1]-1].append(connect_to)
            S[connect_to-1].append(lastcave[-1])
            n_edges += 1

        # add edges to ensure the required sparsity
        while (1- 2*n_edges / (M*(M-1))) > L:
            n1 = np.random.choice(lastcave)
            n2 = np.random.choice(lastcave)

            if n1 not in S[n2-1]:
                S[n2-1].append(n1)
                S[n1-1].append(n2)
                n_edges += 1

        caves.append(ptscloud(entry_point,
                              (None, None),
                              [j for j in range(k*M+1, (k+1)*M+1)]))

        entry_point = (np.random.choice(lastcave), (k+1)*M+1)

    # creating costs info (location and overlap costs)
    f = generate_overlaps(S)
    Cmax = 5.0
    c = [Cmax*np.random.uniform() for _ in range(len(S))]
    if verbose:
        print(f"S={S};\nf={f}\n;c={c}")
    return S, f, c, caves


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


def gen_simple_cavemen_inst0():
    """Generates a simple instance with metadata.

    Returns:
      S, f, c, caves
    """
    S = [[1, 2, 4, 5],
         [2, 1, 3, 8, 10],
         [3, 2, 12, 13, 15],
         [4, 1, 5, 6],
         [5, 1, 4, 6, 7],
         [6, 4, 5, 7],
         [7, 6, 5],
         [8, 2, 9, 10, 11],
         [9, 8, 11],
         [10, 2, 8],
         [11, 8, 9],
         [12, 3, 13, 14],
         [13, 3, 12, 14],
         [14, 12, 13, 15, 16],
         [15, 3, 14],
         [16, 14]]

    f = generate_overlaps(S)

    c = [5.0*np.random.uniform() for _ in range(len(S))]
    caves = [ptscloud(2, 0, [1,4,5,6,7]),
             ptscloud(1, 3, [2,8,9,10,11]),
             ptscloud(2, 0, [3,12,13,14,15,16])]

    return S, f, c, caves


def gen_simple_cavemen_inst1():
    """Generates (another) simple instance with metadata.

    Returns:
      S, f, c, caves
    """
    S = [[1,2,3],
         [2,1,3,4],
         [3,1,2,4],
         [4,2,3,5],
         [5,4,6,7],
         [6,5,8],
         [7,5,8],
         [8, 6,7,9],
         [9, 8,11],
         [10, 11],
         [11, 9,10,12],
         [12, 11,13],
         [13, 12,14,15],
         [14, 13,15],
         [15, 13,14,18],
         [16, 18],
         [17, 18,20],
         [18, 15,16,17],
         [19, 20,22],
         [20, 17, 19,21],
         [21, 20,22],
         [22, 19,21,23],
         [23, 22,25,26],
         [24, 26],
         [25, 23, 26],
         [26, 23, 24,25]]

    caves = [ptscloud((None, None), (4,5), [1,2,3,4]),
             ptscloud((4,5), (8,9), [5,6,7,8]),
             ptscloud((8,9), (12,13), [9,10,11,12]),
             ptscloud((12,13), (15,18), [13,14,15]),
             ptscloud((15,18), (17,20), [16,17,18]),
             ptscloud((17,20), (22,23), [19,20,21,22]),
             ptscloud((22,23), (None, None), [23,24,25,26])]

    f = generate_overlaps(S)

    c = [5.0*np.random.uniform() for _ in range(len(S))]

    return S, f, c, caves


class DDSolver:

    def __init__(self, S, f, c, caves):
        self.S = S
        self.f = f
        self.c = c
        self.caves = caves
        self.B = None

    def _add_interim_point(self, x, current_state, new_state,  current_layer,
                           add_costs=None):
        """Adds variable after the current layer."""
        assert self.B is not None
        assert new_state[-1] == x, f"Last var, {new_state[-1]}, is not {x}"

        next_layer = dict()
        state_pos = {var: idx for idx, var in enumerate(current_state)}
        state_pos.update({x: len(state_pos)})

        for state in current_layer:
            for (xval, arc) in [(True, "hi"), (False, "lo")]:
                astate = state + (xval,)
                newstate = tuple(astate[state_pos[j]] for j in new_state)
                if add_costs is None:
                    cost = 0.0
                else:
                    fnodes = {current_state[i]: astate[i]
                              for i in range(len(current_state))}
                    fnodes.update({x: astate[-1]})
                    cost = self._calc_cave(add_costs, fixed_nodes=fnodes)

                if newstate in next_layer:
                    self.B.llink(current_layer[state], next_layer[newstate],
                                 arc, cost)
                else:
                    newnode = self.B.addnode(current_layer[state], arc,
                                             edge_weight=cost)
                    next_layer.update({newstate: newnode})

        return next_layer

    def _calc_cave(self, cave, fixed_nodes):
        """Calculates part of the objective due to the given cave."""
        m = gp.Model()
        m.modelSense = gp.GRB.MINIMIZE

        x = dict()
        y = dict()

        endpoints = [x for x in cave.e1] + [x for x in cave.e2]
        endpoints = [x for x in endpoints if x is not None]

        V = list(np.unique([j for j in cave.S] + endpoints))

        # variables
        for i in V:
            # decision (location) variables
            if i in cave.S:
                x[i] = m.addVar(vtype=gp.GRB.BINARY, name=f"x{i}",
                                obj=self.c[i-1])
                # overlap (aux) variables
                for a in range(1, len(self.S[i-1])+1):
                    y[(i, a)] = m.addVar(vtype=gp.GRB.BINARY,
                                         name=f"y_{i}_{a}",
                                         obj=(self.f[i-1][a] -
                                              self.f[i-1][a-1]))
            else:
                x[i] = m.addVar(vtype=gp.GRB.BINARY, name=f"x{i}",
                                obj=0.0)

        # constraints
        for j in V:
            if j not in cave.S:
                continue

            m.addConstr(gp.quicksum(x[k] for k in self.S[j-1]) ==
                        gp.quicksum(y[(j, a)]
                                    for a in range(1, len(self.S[j-1])+1)))

            for a in range(1, len(self.S[j-1])):
                m.addConstr(y[(j, a)] >= y[(j, a+1)])

        # fixing variables
        print(f"Calc: S={cave.S}, fixed={fixed_nodes}, e1={cave.e1}, e2={cave.e2}")
        for var in fixed_nodes:
            m.addConstr(x[var] == fixed_nodes[var]*1)

        m.display()
        m.update()
        m.optimize()
        assert m.status == gp.GRB.OPTIMAL

        fterm = sum(self.f[j-1][0] for j in cave.S)
        return m.objVal + fterm

    def build_cover_DD(self):
        """Builds a cover/overlap DD for the instance."""
        assert len(self.caves) > 1

        vars_to_add = sum([[c.e1[0], c.e1[1], c.e2[0], c.e2[1]]
                          for c in self.caves], [])
        vars_to_add = np.unique([v for v in vars_to_add if v is not None])
        self.B = BDD(N=len(vars_to_add), weighted=True)
        varnames = []

        current_layer = {tuple(): self.B.addnode(parent_node=None)}

        self.curr_state = []

        for cavenum, cave in enumerate(self.caves[:-1]):
            new_points = []
            drop_points = []

            new_points = list(np.unique([j for j in (cave.e1 + cave.e2)
                                         if (j not in self.curr_state) and (
                                                 j is not None)]))
            drop_points = [j for j in cave.e1 if (j not in cave.e2) and (
                j is not None)]

            # for x in (cave.e1 + cave.e2):
            #     if x is None:
            #         continue

            #     if x not in self.curr_state:
            #         if x not in new_points:  # corner case: e1=(a,b), e2=(b,c)
            #             new_points.append(x)

            #     if (x in cave.e1) and (
            #             x not in drop_points) and (x not in cave.e2):
            #         drop_points.append(x)

            print(f"new={new_points}, drop={drop_points}")
            for x in new_points[:-1]:
                current_layer = self._add_interim_point(x, self.curr_state,
                                                        self.curr_state + [x],
                                                        current_layer)
                varnames.append(x)
                self.curr_state += [x]

            x = new_points[-1]
            if cavenum < len(self.caves)-2:
                # processing the last point
                newstate = [y for y in (self.curr_state + [x])
                            if y not in drop_points]
                current_layer = self._add_interim_point(
                    x, self.curr_state, newstate,
                    current_layer,
                    add_costs=cave)
            else:
                # that's the last point of the pre-last cloud
                # connecting everything to T and F nodes
                # (no new nodes needed, just calculate costs)
                for state in current_layer:
                    for (xval, arc) in [(True, "hi"), (False, "lo")]:
                        fixed_vals = {
                            pt: state[i]
                            for i, pt in enumerate(self.curr_state)}

                        fixed_vals.update({new_points[-1]: xval})

                        fvals2 = {var: fixed_vals[var] for var in fixed_vals
                                  if var in (self.caves[-1].S +
                                             list(self.caves[-1].e1) +
                                             list(self.caves[-1].e2))}
                        cost = self._calc_cave(cave,
                                               fixed_vals) + self._calc_cave(
                                                   self.caves[-1], fvals2)

                        self.B.link(current_layer[state].id, NTRUE, arc, cost)

            varnames.append(x)
            self.curr_state = newstate

        self.B.rename_vars({i: varnames[i-1] for i in self.B.vars})
        return self.B


## Testing code ######################################################

# First: test that instances are correct

@pytest.mark.parametrize("_", [None for _ in range(30)])
def test_inst_gen(_):
    """Tests that instance generator works correctly."""
    S, f, c, caves = gen_simple_cavemen_inst0()
    assert_instance(S, f, c)
    S, f, c, caves = gen_simple_cavemen_inst1()
    assert_instance(S, f, c)
    n = np.random.randint(5, 10)
    M = np.random.randint(3, 10)
    L = np.random.uniform(0.1, 0.9)
    S, f, c, caves = gen_caveman_inst(n, M, L)
    assert_instance(S, f, c)


def solve_with_MIP(S, f, c):
    """Creates a MIP model (for gurobi) from an instance specs.

    Args:
        S (list): list of adjacency lists,
        f (list): overlap cost function, f[j][a],
        c (list): location costs per point.

    Returns:
        objective value, x, and y
    """
    m = gp.Model()
    m.modelSense = gp.GRB.MINIMIZE
    m.setParam("OutputFlag", 0)
    x = dict()
    y = dict()

    # create variables
    for j in range(1, len(S)+1):
        x[j] = m.addVar(vtype=gp.GRB.BINARY, name=f"x_{j}",
                        obj=c[j-1])

        for a in range(1, len(S[j-1])+1):
            y[(j, a)] = m.addVar(vtype=gp.GRB.BINARY, name=f"y_{j}_{a}",
                                 obj=f[j-1][a]-f[j-1][a-1])

    # create constraints
    for j in range(1, len(S)+1):
        m.addConstr(gp.quicksum(x[k] for k in S[j-1]) ==
                    gp.quicksum(y[(j, a)] for a in range(1, len(S[j-1])+1)))

        for a in range(1, len(S[j-1])):
            m.addConstr(y[(j, a)] >= y[(j, a+1)])

    m.update()
    m.optimize()
    assert m.status == gp.GRB.OPTIMAL
    return m, (m.objVal + sum(fs[0] for fs in f)), x, y


@pytest.mark.parametrize("_", [None for _ in range(10)])
def test_BDD_vs_MIP_simple(_):
    """Tests BDD vs MIP over a single topology (random costs)."""
    S, f, c, caves = gen_simple_cavemen_inst1()
    assert compare_BDD_vs_MIP(S, f, c, caves)


@pytest.mark.parametrize("_", [None for _ in range(100)])
def test_BDD_vs_MIP_random(_):
    n = np.random.randint(5, 10)
    M = np.random.randint(3, 10)
    L = np.random.uniform(0.1, 0.9)
    S, f, c, caves = gen_caveman_inst(n, M, L)
    assert compare_BDD_vs_MIP(S, f, c, caves)

def compare_BDD_vs_MIP(S,f,c,caves):
    sol = DDSolver(S, f, c, caves)
    B = sol.build_cover_DD()
    sp = B.shortest_path()

    _, obj, x, y = solve_with_MIP(S, f, c)
    print(f"obj={obj}, sp={sp[0]}")
    return abs(obj - sp[0]) < 1e-3
