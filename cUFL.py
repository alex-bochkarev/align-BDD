"""Uncapacitated Facility Location with color-constraints.

Problem-solving machinery.
---
(c) A. Bochkarev, Clemson University, 2021
abochka@clemson.edu
"""
from UFL import add_BDD_to_MIP
from copy import copy as cpy
import heapq
import numpy as np
from graphviz import Graph
import gurobipy as gp
import BDD as DD
import UFL
from gurobipy import GRB  # noqa

#######################################################################
# 1. Building the BDD representation
def make_label(state):
    """Helper: formats a node label out of the state."""
    return "["+",".join([str(j+1) for j in range(len(state)) if state[j]])+"]"


class DegreeKeeper:
    """Keeps a heap of node degrees."""
    def __init__(self, S=None):
        """Initializes the heap and index."""
        self.dh = []
        self.index = dict()
        if S is None:
            heapq.heapify(self.dh)
        else:
            for j in range(len(S)):
                heapq.heappush(self.dh, (len(S[j]), j+1))
                self.index.update({j+1: len(S[j])})

    def __len__(self):
        """Returns no. of elements in the heap."""
        return len(self.dh)

    def __getitem__(self, key):
        """Returns a node's degree."""
        return self.index[key]

    def get_next(self):
        """Returns the next node ID for processing.

        Currently: a one with the minimal degree.
        """
        return self.dh[0][1]

    def decrement(self, j):
        """Decrements a degree of node `j`."""
        entry = (self.index[j], j)
        self.dh.remove(entry)
        heapq.heapify(self.dh)
        if entry[0] > 1:
            heapq.heappush(self.dh, (entry[0]-1, j))
            self.index.update({j: entry[0]-1})
            return None
        else:
            self.index.pop(j)
            return j

    def push(self, i, d):
        """Pushes a node `i` with degree `d`."""
        heapq.heappush(self.dh, (d, i))
        self.index.update({i: d})

    def pop(self):
        """Pops a node from the heap, returns (`j`, `degree`)."""
        entry = heapq.heappop(self.dh)
        self.index.pop(entry[1])
        return entry[1], entry[0]

    def has_freedom(self, key):
        """Checks whether a node can be further covered."""
        return key in self.index.keys()


def process_node(i, S, f, res_deg=None, B=None, fixed_nodes=None, trash_pipe=None, next_layer=None, node_labels=None):
    """Processes the given node.

    Args:
        i (int): node ID (one-based),
        S (list): customer neighborhood lists,
        f (dict): costs of location,
        fixed_nodes (set): the ones mentioned (fixed) earlier in `B`,
        res_deg (class DegreeKeeper): residual degrees data.

    Returns:
        B (class BDD): updated diagram.
    """
    assert B is not None

    if B is None:
        root_state = np.array([False for _ in range(len(S))], dtype=bool)
        B = DD.BDD(N=len(S[i-1]), vars=[f"stub{j+1}" for j in range(len(S))], weighted=True)
        node_labels = dict({DD.NROOT: make_label(root_state)})

    if res_deg is None:
        res_deg = DegreeKeeper(S)

    if next_layer is None:
        next_layer = {tuple(root_state): B.addnode(None)}

    if fixed_nodes is None:
        fixed_nodes = set()

    for j in S[i-1]:
        print(f"at {j}, fixed nodes are {fixed_nodes}")
        if j in fixed_nodes:
            continue

        fixed_nodes.add(j)

        current_layer = cpy(next_layer)
        next_layer = dict()

        critical_nodes = set()
        for k in S[j-1]:
            res_deg.decrement(k)
            if not res_deg.has_freedom(k):
                critical_nodes.add(k)

        print(f"At {j}, critical nodes = {critical_nodes}")
        if trash_pipe is not None:
            # create a new "false" running node.
            new_tp = B.addnode(trash_pipe, "lo")
            B.llink(trash_pipe, new_tp, "hi")
            trash_pipe = new_tp
            node_labels.update({trash_pipe.id: "💀"})

        for state in current_layer:
            # processing nodes of the BDD, introducing another node
            # of the *original* graph
            node = current_layer[tuple(state)]
            if state in next_layer:
                B.llink(node, next_layer[state], "lo")
            else:
                print(f"state for {critical_nodes} are: {[q for idx, q in enumerate(state) if (idx+1) in critical_nodes]}")
                if False in [q for idx, q in enumerate(state)
                             if (idx+1) in critical_nodes]:
                    if trash_pipe is None:
                        trash_pipe = B.addnode(node, "lo")
                        node_labels.update({trash_pipe.id: "💀"})
                    else:
                        B.llink(node, trash_pipe, "lo")
                else:
                    newnode = B.addnode(node, "lo")
                    next_layer.update({state: newnode})
                    node_labels.update({newnode.id: make_label(state)})

            next_state = list(state)
            for k in S[j-1]:
                next_state[k-1] = True
            next_state = tuple(next_state)

            if next_state in next_layer:
                B.llink(node, next_layer[next_state], "hi",
                        edge_weight=f[j])
            else:
                newnode = B.addnode(node, "hi", edge_weight=f[j])
                next_layer.update({next_state: newnode})
                node_labels.update({newnode.id: make_label(next_state)})

    # process the last node in S[i-1] separately
    # current_layer = cpy(next_layer)

    # for state in current_layer:
    #     node = current_layer[tuple(state)]
    #     if not state[i-1]:
    #         B.link(node.id, DD.NFALSE, "lo")
    #     else:
    #         B.link(node.id, DD.NTRUE, "lo")

    #     next_state = list(state)
    #     for k in S[j-1]:
    #         next_state[k-1] = True

    #     if not next_state[i-1]:
    #         B.link(node.id, DD.NFALSE, "hi", f[i])
    #     else:
    #         B.link(node.id, DD.NTRUE, "hi", f[i])

    # if trash_pipe is not None:
    #     B.link(trash_pipe.id, DD.NFALSE, "lo")
    #     B.link(trash_pipe.id, DD.NFALSE, "hi")

    return B, node_labels, res_deg, fixed_nodes, trash_pipe, next_layer


def build_cover_DD(S, f):  # pylint: disable=too-many-locals,invalid-name
    """Builds a BDD for the UFL problem.

    Introduces `x`-variables only; encodes the condition of covering
    each customer at least once (includes location costs).

    Args:
       S (list): neighborhood list,
       f (dict): location costs.

    Returns:
        The resulting BDD.

    Notes:
        - employs a DP approach with state being the Boolean
            covering at each customer.
    """
    trash_pipe = None
    # prepare the diagram
    N = len(S)
    B = DD.BDD(N=N, vars=[f"stub{i}" for i in range(N)],
               weighted=True)
    freedoms = DegreeKeeper(S)  # node degrees

    root_state = np.array([False for _ in range(len(S))], dtype=bool)
    node_labels = dict({DD.NROOT: make_label(root_state)})

    next_layer = {tuple(root_state): B.addnode(None)}

    k = 1  # layers counter

    while k < N:
        i = freedoms.get_next()  # current 'central' node to process
        print(f"running at {i}, degrees are: {freedoms.index}:")
        print(f"Si to go: {S[i-1]}")
        for j in S[i-1]:
            if f"x{j}" in B.vars:
                continue

            print(f"Introducing node {j}")
            current_layer = cpy(next_layer)
            next_layer = dict()

            # define 'critical' nodes -- that got the last 'freedom'
            # eliminated during the current iteration,
            # so they can affect the grand outcome (`true`/`false` terminal)
            critical_nodes = set()
            for q in S[j-1]:
                if freedoms.has_freedom(q):
                    new_critical = freedoms.decrement(q)
                    if new_critical is not None:
                        critical_nodes.add(new_critical)

            if trash_pipe is not None:
                # create a new "false" running node.
                new_tp = B.addnode(trash_pipe, "lo")
                B.llink(trash_pipe, new_tp, "hi")
                trash_pipe = new_tp
                node_labels.update({trash_pipe.id: "💀"})

            for state in current_layer:
                # processing nodes of the BDD, introducing another node
                # of the *original* graph
                node = current_layer[tuple(state)]

                if state in next_layer:
                    B.llink(node, next_layer[state], "lo")
                else:
                    print(f"state for {critical_nodes} are: {[q for idx, q in enumerate(state) if (idx+1) in critical_nodes]}")
                    if False in [q for idx, q in enumerate(state)
                                 if (idx+1) in critical_nodes]:
                        if trash_pipe is None:
                            trash_pipe = B.addnode(node, "lo")
                            node_labels.update({trash_pipe.id: "💀"})
                        else:
                            B.llink(node, trash_pipe, "lo")
                    else:
                        newnode = B.addnode(node, "lo")
                        next_layer.update({state: newnode})
                        node_labels.update({newnode.id: make_label(state)})

                next_state = list(state)
                for q in S[j-1]:
                    next_state[q-1] = True
                next_state = tuple(next_state)

                if next_state in next_layer:
                    B.llink(node, next_layer[next_state], "hi",
                            edge_weight=f[j])
                else:
                    newnode = B.addnode(node, "hi", edge_weight=f[j])
                    next_layer.update({next_state: newnode})
                    node_labels.update({newnode.id: make_label(next_state)})
            B.rename_vars({f"stub{k-1}": f"x{j}"})
            print(f"renamed k={k-1} to x{j}")
            k += 1

    # process the last node in S[i-1] separately
    print(f"freedoms index={freedoms.index}")
    i = -1
    while f"x{i}" in B.vars or i == -1:
        i, _ = freedoms.pop()

    print(f"Processing the last layer {k}, introducing node {i}")
    current_layer = cpy(next_layer)

    for state in current_layer:
        node = current_layer[tuple(state)]
        if not state[i-1]:
            B.link(node.id, DD.NFALSE, "lo")
        else:
            B.link(node.id, DD.NTRUE, "lo")

        next_state = list(state)
        for q in S[i-1]:
            next_state[q-1] = True

        if not next_state[i-1]:
            B.link(node.id, DD.NFALSE, "hi", f[i])
        else:
            B.link(node.id, DD.NTRUE, "hi", f[i])

    if trash_pipe is not None:
        B.link(trash_pipe.id, DD.NFALSE, "lo")
        B.link(trash_pipe.id, DD.NFALSE, "hi")

    B.rename_vars({f"stub{k-1}": f"x{i}"})
    print(f"renamed k={k-1} to x{i}")

    return B, node_labels


def build_color_DD(f, f_color, k_bar):  # pylint: disable=invalid-name
    """Builds a BDD encoding location-level constraints.

    Deals with no. locations per color and location costs.

    Args:
        f (dict): location costs,
        colors (dict): color codes per facility
        k_bar (np.array): max. number of locations per color.

    Returns:
        D (class BDD): resulting diagram.
        nl (dict): string node labels ("states"), `id` -> `label` (string).
    """
    D = DD.BDD(N=len(f_color), vars=[f"stub{i}" for i in range(1, len(f_color)+1)],  # noqa pylint: disable=invalid-name, unnecessary-comprehension
               weighted=True)
    # a *state* is the number of located facilities
    # for the *current* color (a single number)
    n = 1  # (original) nodes counter
    N = len(f_color)

    node_labels = dict({DD.NROOT: 0})
    next_layer = {0: D.addnode(None)}
    trash_pipe = None

    # TODO: here I'd need to decide the order of colors
    #       and order of nodes within a color
    colors = [c for c in range(len(k_bar))]
    customers = [[C+1 for (C, f_c) in enumerate(f_color) if f_c == c]
                 for c in colors]

    for c in colors:
        for customer in customers[c][:-1]:
            if n == N:
                break

            current_layer = cpy(next_layer)
            next_layer = dict()
            if trash_pipe is not None:
                # create a new "false" running node.
                new_tp = D.addnode(trash_pipe, "lo")
                D.llink(trash_pipe, new_tp, "hi")
                trash_pipe = new_tp
                node_labels.update({trash_pipe.id: "💀"})

            for state in current_layer:
                node = current_layer[state]
                if state in next_layer:
                    D.llink(node, next_layer[state], "lo")
                else:
                    newnode = D.addnode(node, "lo")
                    next_layer.update({state: newnode})
                    node_labels.update({newnode.id: str(state)})

                if (state+1) in next_layer:
                    D.llink(node, next_layer[state+1], "hi",
                            edge_weight=f[customer])
                else:
                    if (state+1) > k_bar[c]:
                        if trash_pipe is None:
                            trash_pipe = D.addnode(node, "hi")
                            node_labels.update({trash_pipe.id: "💀"})
                        else:
                            D.llink(node, trash_pipe, "hi")
                    else:
                        newnode = D.addnode(node, "hi")
                        next_layer.update({state+1: newnode})
                        node_labels.update({newnode.id: str(state+1)})
            D.rename_vars({f"stub{n}": f"x{customer}"})
            n += 1

        # Processing the last customer separately
        # (within a color)

        if n < N:
            current_layer = cpy(next_layer)
            next_layer = dict()

            if trash_pipe is not None:
                # create a new "false" running node.
                new_tp = D.addnode(trash_pipe, "lo")
                D.llink(trash_pipe, new_tp, "hi")
                trash_pipe = new_tp
                node_labels.update({trash_pipe.id: "💀"})

            for state in current_layer:
                node = current_layer[state]
                new_state = 0  # we 'reset' the color counter
                if new_state in next_layer:
                    D.llink(node, next_layer[new_state], "lo")
                else:
                    newnode = D.addnode(node, "lo")
                    next_layer.update({new_state: newnode})
                    node_labels.update({newnode.id: str(new_state)})

                if (state+1) > k_bar[c]:
                    if trash_pipe is None:
                        trash_pipe = D.addnode(node, "hi")
                        node_labels.update({trash_pipe.id: "💀"})
                    else:
                        D.llink(node, trash_pipe, "hi")
                else:
                    D.llink(node, next_layer[new_state], "hi")

            D.rename_vars({f"stub{n}": f"x{customers[c][-1]}"})
            n += 1
        else:
            # the last layer of the DD
            if trash_pipe is not None:
                D.link(trash_pipe.id, DD.NFALSE, "hi")
                D.link(trash_pipe.id, DD.NFALSE, "lo")

            for state in next_layer:
                node_id = next_layer[state].id
                if state+1 > k_bar[c]:
                    y_target = DD.NFALSE
                else:
                    y_target = DD.NTRUE

                D.link(node_id, y_target, "hi")
                D.link(node_id, DD.NTRUE, "lo")

            D.rename_vars({f"stub{n}": f"{customers[c][-1]}"})
    return D, node_labels


######################################################################
# 2. Building a MIP


def build_cUFL_MIP(S, f, f_color, k_bar):
    """Builds 'plain-vanilla' MIP instance for 'colorful' UFL.

    Args:
        S (list): list of customer neighborhoods,
        f (dict): facility costs,
        f_color (list): facility colors,
        k_bar (list): max no. of locations per color,

    Returns:
        `gurobipy` model (`gp.Model`)
    """
    m = gp.Model()
    m.modelSense = GRB.MINIMIZE

    x = dict()

    V = range(1, len(S)+1)

    for i in V:
        x[i] = m.addVar(vtype=GRB.BINARY, name=f"x{i}", obj=f[i])

    # covering constraints
    for j in V:
        m.addConstr(gp.quicksum(x[i] for i in S[j-1]) >= 1)

    # color-budget constraints
    for c in range(len(k_bar)):
        m.addConstr(gp.quicksum(x[i] for i in V if f_color[i-1] == c) <= k_bar[c])

    m.update()
    return m


######################################################################
# 3. Quick-testing code
def solve_with_MIP(S, f, fc, kb):
    """Solves the problem with a naive MIP approach."""
    m_naive = build_cUFL_MIP(S, f, fc, kb)
    m_naive.update()
    m_naive.optimize()
    return m_naive


def generate_test_instance(n):
    """Generates a (colored) facility location instance.

    Args:
        n (int): number of nodes (customers / facilities)

    Returns:
        S (list): neighborhood list,
        f (dict): location costs
        f_colors (list): location colors,
        k_bar (list): budget per color.
    """
    F_MIN = 5
    F_MAX = 10
    C_MIN = 2
    C_MAX = 10
    MAX_DEGREE = min(5, n)
    MAX_BUDGET = 5

    status = -1
    while status != GRB.OPTIMAL:
        N = [i for i in range(1,n+1)]  # a set of nodes

        S = [[i] + np.random.choice([v for v in N if v!=i],
                                    np.random.randint(1, MAX_DEGREE),
                                    replace=False).tolist()
            for i in N]

        f = {i: np.random.randint(F_MIN, F_MAX) for i in N}

        no_colors = np.random.randint(C_MIN, C_MAX)
        f_colors = [np.random.randint(0, no_colors) for _ in N]

        k_bar = [np.random.randint(1, 1+MAX_BUDGET) for _ in range(no_colors)]
        status = solve_with_MIP(S, f, f_colors, k_bar).status

    return S, f, f_colors, k_bar


def make_simple_problem():
    """Generates a simple problem instance.

    Returns:
    S, f, f_colors, k_bar.
    """
    S = [[1, 2], [1, 2, 3, 5], [2, 3, 4], [3, 4, 5], [2, 4, 5, 6, 7], [5, 6], [5, 7]]
    f = {1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7}
    f_colors = [0, 1, 2, 2, 1, 0, 0]
    k_bar = [5, 1, 3]
    return S, f, f_colors, k_bar

def draw_problem_dia(S, f, f_colors, k_bar,
                     filename="run_logs/c_problem_dia.gv"):
    """Draws the bipartite graph describing the problem."""
    assert len(np.unique(f_colors)) <= 6
    cols = {0: 'red', 1: 'blue', 2: 'green', 3: 'yellow',
            4: 'black', 5: 'purple', 6: 'orange', 7: 'white', 8: 'gray'}

    dia = Graph('G', comment="Uncapacitated Facility Location")

    edges = set()
    for i in range(len(f_colors)):
        dia.node(f"{i+1}", label=f"{i+1} ({f[i+1]})", style='filled', color=cols[f_colors[i]])
        for j in S[i]:
            if i+1 != j and (i+1, j) not in edges and (j, i+1) not in edges:
                dia.edge(f"{i+1}", f"{j}")
                edges.add((i+1, j))

    print("Color limits: " + ", ".join([f"{cols[c]}: {k_bar[c]}" for c in range(len(k_bar))]))
    dia.view(filename=filename)


def check_simple_example():
    """Shows a simple example (for visual inspection)."""
    S, f, fc, kb = make_simple_problem()
    draw_problem_dia(S, f, fc, kb, "run_logs/test")
    build_cUFL_MIP(S, f, fc, kb).display()

    oD, nl_oD = build_cover_DD(S, f)
    lD, nl_lD = build_color_DD(f, fc, kb)

    oD.dump_gv(node_labels=nl_oD).view(filename="overlap")
    lD.dump_gv(node_labels=nl_lD).view(filename="location")

def solve_with_BDD_MIP(S, f, fc, kb):
    """Sovles the problem by building a MIP from two BDDs."""
    cover_DD, cover_nl = build_cover_DD(S, f)
    color_DD, color_nl = build_color_DD(f, fc, kb)

    m, c, v, x = UFL.add_BDD_to_MIP(cover_DD, prefix="cover")
    m, c, v, x = UFL.add_BDD_to_MIP(color_DD, m, x, "color")
    m.update()
    m.optimize()
    return m

def test_color_UFL():
    """Tests the formulation for color-UFL (overlap DD)."""
    S, f, fc, kb = make_simple_problem()

    m_naive = solve_with_MIP(S, f, fc, kb)
    m = solve_with_BDD_MIP(S, f, fc, kb)
    print(f"Naive model: status={m_naive.status}, obj={m_naive.objVal}")
    print(f"BDD model: status={m.status}, obj={m.objVal}")
    assert m_naive.objVal == m.objVal, f"Naive: {m_naive.objVal} (status {m_naive.status}), while BDD: {m.objVal} (status {m.status})"

def test_random_UFL():
    """Tests the formulation for color-UFL (overlap DD) -- random instance."""
    S, f, fc, kb = generate_test_instance(5)

    m_naive = solve_with_MIP(S, f, fc, kb)
    m = solve_with_BDD_MIP(S, f, fc, kb)
    assert m_naive.objVal == m.objVal, f"Naive: {m_naive.objVal} (status {m_naive.status}), while BDD: {m.objVal} (status {m.status})"
