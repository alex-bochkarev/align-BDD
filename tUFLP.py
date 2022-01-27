"""Implements the basic machinery for Uncapacitated Facility Location
problem with type-constraints.

Building the diagrams (cover DD and type DD) and solving the
corresponding CPP instance.

Tests coverage: :py:mod:`tUFLP_test` and :py:mod:`UFLP_test`.
"""
from UFL import add_BDD_to_MIP, generate_test_instance
from copy import copy as cpy
import heapq
import numpy as np
from graphviz import Graph
from itertools import permutations
import networkx as nx
import gurobipy as gp
import BDD as DD
import BB_search as bb
import varseq as vs
import UFL
from gurobipy import GRB  # noqa
import pytest


#######################################################################
# 1. Building the BDD representation
def make_label(state):
    """Helper: formats a node label out of the state."""
    return "["+",".join([str(j+1) for j in range(len(state)) if state[j]])+"]"


class DegreeKeeper:
    """Keeps a heap of node degrees.

    Note: it works correctly iff nodes are
    numbered starting with one.
    """
    def __init__(self, S=None, next_node_type="min"):
        """Initializes the heap and index.

        EXPERIMENT BRACH version:
        Allows for different approaches to `get_next` node,
        parameterized by `next_node_type`:

        - `min`: minimum residual degree,
        - `max`: resp., maximum,
        - `rnd`: getting random node.
        """
        self.dh = []
        self.index = dict()

        assert next_node_type in ['min', 'max', 'rnd']
        self.next_type = next_node_type

        if S is None:
            heapq.heapify(self.dh)
        else:
            for j in range(len(S)):
                if next_node_type == 'min' or next_node_type == 'rnd':
                    heapq.heappush(self.dh, (len(S[j]), j+1))
                    self.index.update({j+1: len(S[j])})
                elif next_node_type == 'max':
                    heapq.heappush(self.dh, (-len(S[j]), j+1))
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
        if self.next_type in ['min', 'max']:
            return self.dh[0][1]
        else:
            return np.random.choice(list(self.index.keys()))

    def decrement(self, j):
        """Decrements a degree of node `j`."""
        if self.next_type in ['min', 'rnd']:
            entry = (self.index[j], j)
        else:
            entry = (-self.index[j], j)

        self.dh.remove(entry)
        heapq.heapify(self.dh)
        if (self.next_type == "min" or self.next_type == "rnd") and entry[0] > 1:
            heapq.heappush(self.dh, (entry[0]-1, j))
            self.index.update({j: entry[0]-1})
            return None
        elif self.next_type == "max" and entry[0] < -1:
            heapq.heappush(self.dh, (entry[0]+1, j))
            self.index.update({j: -(entry[0]+1)})
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
        if self.next_type in ['min', 'rnd']:
            return entry[1], entry[0]
        else:
            return entry[1], -entry[0]  # this is 'max' type

    def has_freedom(self, key):
        """Checks whether a node can be further covered."""
        return key in self.index.keys()


class TypeSorter:
    """Sorts type blocks (to derive a good var order for the type BDD)."""
    def __init__(self, f_types, node_order):
        self.f_types = f_types
        self.node_order = node_order
        self.no_types = len(np.unique(f_types))
        self.pos = {node_order[i]: i for i in range(len(node_order))}

    def rel_weight(self, elements_in_A, elements_in_B):
        """Calculates a 'relative weight' of the type.

        Args:
            elements_in_A (list): elements with type A,
            elements_in_B (list): elements with type B.
        Returns:
            RW (int): 'relative weight' of type A (relative to B)

        Notes:
            In fact, if RW > 0, then configuration (A,B) will have fewer
            inverstions with `node_order` (and the other way around if
            `RW` < 0).
        """
        RW = 0
        for a in elements_in_A:
            for b in elements_in_B:
                if self.pos[a] > self.pos[b]:
                    RW += 1
                elif self.pos[a] < self.pos[b]:
                    RW -= 1
        return RW

    def sort_types(self):
        """Sorts type blocks and nodes within each type.

        The aim is to have the minimum number of inversions between the resulting
        typing diagram and the `node_order` parameter.

        Returns:
            types (list): ordered type numbers.
            customers (dict): (ordered) list of customer per type.
        """
        types = [set([idx+1 for idx, ptype in enumerate(self.f_types)
                       if ptype == c]) for c in range(self.no_types)]
        col_idx = list(range(len(types)))

        for _ in range(self.no_types**2 +1):
            for i in range(self.no_types-1):
                RW = self.rel_weight(types[i], types[i+1])
                if RW >= 0:
                    c = cpy(types[i])
                    types[i] = cpy(types[i+1])
                    types[i+1] = c

                    c = col_idx[i]
                    col_idx[i] = col_idx[i+1]
                    col_idx[i+1] = c

        customers = {c: [] for c in range(self.no_types)}

        for j in self.node_order:
            customers[self.f_types[j-1]].append(j)

        return col_idx, customers


# pylint: disable=too-many-locals,invalid-name
def build_cover_DD(S, f, next_node_type='min'):
    """Builds a BDD for the UFL problem.

    Introduces `x`-variables only; encodes the condition of covering
    each customer at least once (includes location costs).

    Args:
        S (list): neighborhood list,
        f (dict): location costs.
        next_node_type (str): `min`, `max`, or `rnd` -- see `DegreeKeeper`
                                docstring for details.
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
    freedoms = DegreeKeeper(S, next_node_type)  # node degrees

    root_state = np.array([False for _ in range(len(S))], dtype=bool)
    node_labels = dict({DD.NROOT: make_label(root_state)})

    next_layer = {tuple(root_state): B.addnode(None)}

    k = 1  # layers counter

    while k < N:
        i = freedoms.get_next()  # current 'central' node to process
        for j in S[i-1]:
            if f"x{j}" in B.vars or k == N:
                continue

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
                B.llink(trash_pipe, new_tp, "hi", edge_weight=f[j])
                trash_pipe = new_tp
                node_labels.update({trash_pipe.id: "💀"})

            for state in current_layer:
                # processing nodes of the BDD, introducing another node
                # of the *original* graph
                node = current_layer[tuple(state)]

                if state in next_layer:
                    B.llink(node, next_layer[state], "lo")
                else:
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
            k += 1

    # process the last node in S[i-1] separately
    i = -1
    while f"x{i}" in B.vars or i == -1:
        i, _ = freedoms.pop()

    current_layer = cpy(next_layer)

    for state in current_layer:
        node = current_layer[tuple(state)]
        if not np.all([state[j-1] for j in S[i-1]]):
            B.link(node.id, DD.NFALSE, "lo")
        else:
            B.link(node.id, DD.NTRUE, "lo")

        next_state = list(state)
        for q in S[i-1]:
            next_state[q-1] = True

        if not np.all([next_state[j-1] for j in S[i-1]]):
            B.link(node.id, DD.NFALSE, "hi", f[i])
        else:
            B.link(node.id, DD.NTRUE, "hi", f[i])

    if trash_pipe is not None:
        B.link(trash_pipe.id, DD.NFALSE, "lo")
        B.link(trash_pipe.id, DD.NFALSE, "hi", f[i])

    B.rename_vars({f"stub{k-1}": f"x{i}"})

    return B, node_labels


def build_randomized_cover_DD(S, f):
    """Builds a BDD for the typed-UFL problem.

    Introduces `x`-variables only; encodes the condition of covering
    each customer at least once (includes location costs).

    *Processes the nodes in random order.*

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
    freedoms = {i+1: len(S[i]) for i in range(N)}

    root_state = np.array([False for _ in range(len(S))], dtype=bool)
    node_labels = dict({DD.NROOT: make_label(root_state)})

    next_layer = {tuple(root_state): B.addnode(None)}

    residual_nodes = set(i for i in range(1, N+1))

    for k in range(1, N):
        # pick a random node of the original graph
        i = np.random.choice([x for x in list(freedoms.keys())
                              if freedoms[x] > 0 and x in residual_nodes])

        current_layer = cpy(next_layer)
        next_layer = dict()

        if trash_pipe is not None:
            # create a new "false" running node.
            new_tp = B.addnode(trash_pipe, "lo")
            B.llink(trash_pipe, new_tp, "hi", edge_weight=f[i])
            trash_pipe = new_tp
            node_labels.update({trash_pipe.id: "💀"})

        for q in S[i-1]:
            freedoms[q] = max(freedoms[q]-1, 0)

        for state in current_layer:
            # processing nodes of the BDD, introducing another node
            # of the *original* graph
            node = current_layer[tuple(state)]

            if state in next_layer:
                B.llink(node, next_layer[state], "lo")
            else:
                if False in [state[r] for r in range(N) if freedoms[r+1] == 0]:
                    # we won't be able to cover the node further
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
            for q in S[i-1]:
                next_state[q-1] = True

            next_state = tuple(next_state)

            if next_state in next_layer:
                B.llink(node, next_layer[next_state], "hi",
                        edge_weight=f[i])
            else:
                newnode = B.addnode(node, "hi", edge_weight=f[i])
                next_layer.update({next_state: newnode})
                node_labels.update({newnode.id: make_label(next_state)})

        B.rename_vars({f"stub{k-1}": f"x{i}"})
        residual_nodes.remove(i)

    # process the last layer separately
    i = residual_nodes.pop()

    current_layer = cpy(next_layer)

    for state in current_layer:
        node = current_layer[tuple(state)]
        if not np.all([state[j-1] for j in S[i-1]]):
            B.link(node.id, DD.NFALSE, "lo")
        else:
            B.link(node.id, DD.NTRUE, "lo")

        next_state = list(state)
        for q in S[i-1]:
            next_state[q-1] = True

        if not np.all([next_state[j-1] for j in S[i-1]]):
            B.link(node.id, DD.NFALSE, "hi", f[i])
        else:
            B.link(node.id, DD.NTRUE, "hi", f[i])

    if trash_pipe is not None:
        B.link(trash_pipe.id, DD.NFALSE, "lo")
        B.link(trash_pipe.id, DD.NFALSE, "hi", f[i])

    B.rename_vars({f"stub{N-1}": f"x{i}"})

    return B, node_labels


def build_type_DD(f, f_types, k_bar, preferred_order=None):  # pylint: disable=invalid-name
    """Builds a BDD encoding location-level constraints.

    Deals with no. locations per type and location costs.

    Args:
        f (dict): location costs,
        types (dict): type codes per facility
        k_bar (np.array): max. number of locations per type.
        preferred_order (list): order of nodes to align to, as possible
                                (e.g., the one of the already built cover DD)

    Returns:
        D (class BDD): resulting diagram.
        nl (dict): string node labels ("states"), `id` -> `label` (string).

    Notes:
        If `preferred_order` is given, randomly shuffles types, but makes sensible
        order of customers within each type (corresponding to `preferred_order`).
    """
    D = DD.BDD(N=len(f_types), vars=[f"stub{i}" for i in range(1, len(f_types)+1)],  # noqa pylint: disable=invalid-name, unnecessary-comprehension
               weighted=True)
    # a *state* is the number of located facilities
    # for the *current* type (a single number)
    n = 1  # (original) nodes counter
    N = len(f_types)

    node_labels = dict({DD.NROOT: 0})
    next_layer = {0: D.addnode(None)}
    trash_pipe = None

    if preferred_order is None:
        types = [c for c in range(len(k_bar))]
        customers = [[C+1 for (C, f_c) in enumerate(f_types) if f_c == c]
                    for c in types]
    else:
        types = [c for c in range(len(k_bar))]
        customers = {c: [] for c in range(len(k_bar))}

        for j in preferred_order:
            customers[f_types[j-1]].append(j)

        types = list(np.random.permutation(types))

    for c in types:
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
        # (within a type)

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
                new_state = 0  # we 'reset' the type counter
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

            D.rename_vars({f"stub{n}": f"x{customers[c][-1]}"})
    return D, node_labels


def build_randomized_type_DD(f, f_types, k_bar):  # pylint: disable=invalid-name
    """Builds a BDD encoding location-level constraints.

    Deals with no. locations per type and location costs.
    *Uses random order of types, and random order
    of nodes within each type*

    Args:
        f (dict): location costs,
        type (dict): type codes per facility
        k_bar (np.array): max. number of locations per type.

    Returns:
        D (class BDD): resulting diagram.
        nl (dict): string node labels ("states"), `id` -> `label` (string).
    """
    D, _ = build_type_DD(f, f_types, k_bar)
    D.align_to(np.random.permutation(D.vars), inplace=True)
    return D, None


######################################################################
# 2. Building a MIP
def build_tUFLP_MIP(S, f, f_types, k_bar):
    """Builds 'plain-vanilla' MIP instance for typed UFLP.

    Args:
        S (list): list of customer neighborhoods,
        f (dict): facility costs,
        f_types (list): facility types,
        k_bar (list): max no. of locations per type,

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

    # type-budget constraints
    for c in range(len(k_bar)):
        m.addConstr(gp.quicksum(x[i] for i in V if f_types[i-1] == c) <= k_bar[c])

    m.update()
    return m


######################################################################
# 3. Quick-testing code
def solve_with_MIP(S, f, fc, kb):
    """Solves the problem with a naive MIP approach."""
    m_naive = build_tUFLP_MIP(S, f, fc, kb)
    m_naive.update()
    m_naive.optimize()
    return m_naive


def generate_test_instance(n, p=0.3):
    """Generates a (typed) facility location instance.

    Args:
        n (int): number of nodes (customers / facilities)
        p (float): graph generation parameter
            (prob of an edge between a pair of nodes)

    Returns:
        S (list): neighborhood list,
        f (dict): location costs
        f_types (list): location types,
        k_bar (list): budget per type.
    """
    F_MIN = 5
    F_MAX = 10
    C_MIN = 2
    C_MAX = min(10, n)
    MAX_DEGREE = min(5, n)
    MAX_BUDGET = 5

    status = -1
    N = [i for i in range(1,n+1)]  # a set of nodes

    while status != GRB.OPTIMAL:

        G = nx.gnp_random_graph(n, p, directed=False)
        AM = nx.adjacency_matrix(G)
        S = [[i+1] for i in range(n)]

        for i in range(n):
            for j in range(n):
                if AM[(i,j)] == 1:
                    if j+1 not in S[i]: S[i].append(j+1)
                    if i+1 not in S[j]: S[j].append(i+1)

        f = {i: np.random.randint(F_MIN, F_MAX) for i in N}

        no_types = np.random.randint(C_MIN, C_MAX)
        N_res = N
        f_types = [0] * n

        for c in range(1, no_types):
            facilities_c = np.random.choice(N_res,
                                            np.random.randint(1, 1+ len(N_res) - (no_types-c)),
                                            replace=False)
            for facility in facilities_c:
                f_types[facility-1] = c

            N_res = [v for v in N_res if v not in facilities_c]

        k_bar = [np.random.randint(1, 1+MAX_BUDGET) for _ in range(no_types)]
        status = solve_with_MIP(S, f, f_types, k_bar).status

    return (S, f, f_types, k_bar)


def generate_test_instance_alt(n, p=0.3, nodes_per_type=5):
    """Generates a (typed) facility location instance.

    Args:
        n (int): number of nodes (customers / facilities)
        p (float): graph generation parameter
            (prob of an edge between a pair of nodes)
        nodes_per_type(int): type generation parameter

    Returns:
        S (list): neighborhood list,
        f (dict): location costs
        f_types (list): location types,
        k_bar (list): budget per type.
    """
    F_MIN = 5
    F_MAX = 10
    MAX_BUDGET = 5

    status = -1
    N = [i for i in range(1,n+1)]  # a set of nodes

    while status != GRB.OPTIMAL:

        G = nx.gnp_random_graph(n, p, directed=False)
        AM = nx.adjacency_matrix(G)
        S = [[i+1] for i in range(n)]

        for i in range(n):
            for j in range(n):
                if AM[(i,j)] == 1:
                    if j+1 not in S[i]: S[i].append(j+1)
                    if i+1 not in S[j]: S[j].append(i+1)

        f = {i: np.random.randint(F_MIN, F_MAX) for i in N}

        no_types = round(n / nodes_per_type)
        N_res = N
        np.random.shuffle(N_res)
        N_res = list(N_res)

        f_types = [0] * n

        for c in range(no_types):
            cur_node = N_res.pop()
            f_types[cur_node - 1] = c

        while len(N_res) > 0:
            f_types[N_res.pop()-1] = np.random.randint(0, no_types)

        k_bar = [np.random.randint(1, 1+MAX_BUDGET) for _ in range(no_types)]
        status = solve_with_MIP(S, f, f_types, k_bar).status

    return (S, f, f_types, k_bar)


def generate_string_instance(n):
    """Generates a (typed) facility location instance: string.

    Args:
        n (int): number of nodes (customers / facilities)

    Returns:
        S (list): neighborhood list,
        f (dict): location costs
        f_types (list): location types,
        k_bar (list): budget per type.
    """
    F_MIN = 5
    F_MAX = 10
    C_MIN = 2
    C_MAX = min(10, n)
    MAX_BUDGET = 5

    f = {i: np.random.randint(F_MIN, F_MAX + 1)
         for i in range(1, n+1)}

    status = -1
    S = [[] for _ in range(n)]
    S[0] = [1, 2]
    S[-1] = [n-1, n]

    for k in range(1, n-1):
        S[k] = [k, k+1, k+2]  # so, e.g. "3" is accessible from "2","3","4"

    N = [i for i in range(1,n+1)]  # a set of nodes
    while status != GRB.OPTIMAL:
        no_types = np.random.randint(C_MIN, C_MAX)
        N_res = N
        f_types = [0] * n

        for c in range(1, no_types):
            facilities_c = np.random.choice(N_res,
                                            np.random.randint(1, 1+ len(N_res) - (no_types-c)),
                                            replace=False)
            for facility in facilities_c:
                f_types[facility-1] = c

            N_res = [v for v in N_res if v not in facilities_c]

        k_bar = [np.random.randint(1, 1+MAX_BUDGET) for _ in range(no_types)]
        status = solve_with_MIP(S, f, f_types, k_bar).status

    return (S, f, f_types, k_bar)


def generate_organic_instance(n):
    """Generates a (typed) facility location instance: 'organic' thing.

    Args:
        n (int): number of nodes (customers / facilities)

    Returns:
        S (list): neighborhood list,
        f (dict): location costs
        f_types (list): location types,
        k_bar (list): budget per type.
    """
    F_MIN = 5
    F_MAX = 10
    C_MIN = 2
    C_MAX = min(10, n)
    MAX_BUDGET = 5

    f = {i: np.random.randint(F_MIN, F_MAX + 1)
         for i in range(1, n+1)}

    status = -1
    S = [[] for _ in range(n)]

    k = 0
    S[0] = [1, 2]

    prev_node = 1
    k = 2
    while k < n-1:
        if np.random.uniform() <= 0.5:
            # degree-2 node
            S[k - 1] = [prev_node, k, k+1]
            prev_node = k
        else:
            # deg-3 node + deg-1 node
            S[k - 1] = [prev_node, k, k+1, k+2]
            S[(k+1) - 1] = [k, k+1]
            prev_node = k
            k += 1

        k += 1

    if k == n-1:
        S[k-1] = [prev_node, n-1, n]
        S[k] = [n-1, n]
    elif k == n:
        S[k-1] = [n-2, n]

    print(f"k={k}, n={n}")
    print(f"S={S}")
    # 'typing' the nodes
    N = [i for i in range(1,n+1)]  # a set of nodes
    while status != GRB.OPTIMAL:
        no_types = np.random.randint(C_MIN, C_MAX)
        N_res = N
        f_types = [0] * n

        for c in range(1, no_types):
            facilities_c = np.random.choice(N_res,
                                            np.random.randint(1, 1+ len(N_res) - (no_types-c)),
                                            replace=False)
            for facility in facilities_c:
                f_types[facility-1] = c

            N_res = [v for v in N_res if v not in facilities_c]

        k_bar = [np.random.randint(1, 1+MAX_BUDGET) for _ in range(no_types)]
        status = solve_with_MIP(S, f, f_types, k_bar).status

    return (S, f, f_types, k_bar)


def generate_d4_instance(n):
    """Generates a (typed) facility location instance: degree-4 nodes.

    Args:
        n (int): number of nodes (customers / facilities)

    Returns:
        S (list): neighborhood list,
        f (dict): location costs
        f_types (list): location types,
        k_bar (list): budget per type.
    """
    F_MIN = 5
    F_MAX = 10
    C_MIN = 2
    C_MAX = min(10, n)
    MAX_BUDGET = 5

    f = {i: np.random.randint(F_MIN, F_MAX + 1)
         for i in range(1, n+1)}

    status = -1
    S = [[] for _ in range(n)]

    k = 0
    S[0] = [1, 2]

    prev_node = 1
    k = 2
    while k < n-3:
        if np.random.uniform() <= 0.5:
            # degree-2 node
            S[k - 1] = [prev_node, k, k+1]
            prev_node = k
        else:
            # deg-4 node
            S[k - 1] = [prev_node, k, k+1, k+2, k+3]
            S[(k+1) - 1] = [k, k+1]
            S[(k+2) - 1] = [k, k+2]
            prev_node = k
            k += 2

        k += 1

    while k < n:
        S[k-1] = [prev_node, k, k+1]
        prev_node = k
        k += 1

    S[-1]=[n-1, n]

    print(f"k={k}, n={n}")
    print(f"S={S}")
    # typing the nodes
    N = [i for i in range(1,n+1)]  # a set of nodes
    while status != GRB.OPTIMAL:
        no_types = np.random.randint(C_MIN, C_MAX)
        N_res = N
        f_types = [0] * n

        for c in range(1, no_types):
            facilities_c = np.random.choice(N_res,
                                            np.random.randint(1, 1+ len(N_res) - (no_types-c)),
                                            replace=False)
            for facility in facilities_c:
                f_types[facility-1] = c

            N_res = [v for v in N_res if v not in facilities_c]

        k_bar = [np.random.randint(1, 1+MAX_BUDGET) for _ in range(no_types)]
        status = solve_with_MIP(S, f, f_types, k_bar).status

    return (S, f, f_types, k_bar)


def generate_simple_problem():
    """Generates a simple problem instance.

    Returns:
    S, f, f_types, k_bar.
    """
    S = [[1, 2], [1, 2, 3, 5], [2, 3, 4], [3, 4, 5], [2, 4, 5, 6, 7], [5, 6], [5, 7]]
    f = {1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7}
    f_types = [0, 1, 2, 2, 1, 0, 0]
    k_bar = [5, 1, 3]
    return S, f, f_types, k_bar

def generate_5d_problem(n=None):
    """Generates a simple problem instance, with two nodes of degree 5.

    Returns:
    S, f, f_types, k_bar.
    """
    S = [[1, 2, 3, 4, 5,6], [2, 1], [3, 1], [4, 1], [5, 1], [6, 1, 7, 8, 9, 10],
         [7, 6], [8, 6], [9, 6], [10, 6]]
    f = {1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, 10: 10}
    f_types = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    k_bar = [11]
    return S, f, f_types, k_bar


def draw_problem_dia(S, f, f_types, k_bar,
                     filename="run_logs/c_problem_dia.gv"):
    """Draws the bipartite graph describing the problem."""
    assert len(np.unique(f_types)) <= 10
    cols = {0: 'red', 1: 'blue', 2: 'green', 3: 'yellow',
            4: 'black', 5: 'purple', 6: 'orange', 7: 'white', 8: 'gray',
            9: 'cyan', 10: 'magenta'}

    dia = Graph('G', comment="Uncapacitated Facility Location")

    edges = set()
    for i in range(len(f_types)):
        dia.node(f"{i+1}", label=f"{i+1} ({f[i+1]})", style='filled', color=cols[f_types[i]])
        for j in S[i]:
            if i+1 != j and (i+1, j) not in edges and (j, i+1) not in edges:
                dia.edge(f"{i+1}", f"{j}")
                edges.add((i+1, j))

    print("Type limits: " + ", ".join([f"{cols[c]}: {k_bar[c]}" for c in range(len(k_bar))]))
    dia.view(filename=filename)


def check_simple_example():
    """Shows a simple example (for visual inspection)."""
    S, f, fc, kb = generate_simple_problem()
    draw_problem_dia(S, f, fc, kb, "run_logs/test")
    build_tUFLP_MIP(S, f, fc, kb).display()

    oD, nl_oD = build_cover_DD(S, f)
    lD, nl_lD = build_type_DD(f, fc, kb)

    oD.dump_gv(node_labels=nl_oD).view(filename="overlap")
    lD.dump_gv(node_labels=nl_lD).view(filename="location")

def solve_with_BDD_MIP(S, f, fc, kb):
    """Sovles the problem by building a MIP from two BDDs."""
    cover_DD, cover_nl = build_cover_DD(S, f)
    pref_order = [int(x[1:]) for x in cover_DD.vars]
    type_DD, type_nl = build_type_DD(f, fc, kb, pref_order)

    m, c, v, x = UFL.add_BDD_to_MIP(cover_DD, prefix="cover")
    m, c, v, x = UFL.add_BDD_to_MIP(type_DD, m, x, "type")
    m.update()
    m.optimize()
    return m

def test_type_UFL():
    """Tests the formulation for typed-UFLP (overlap DD)."""
    S, f, fc, kb = generate_simple_problem()

    m_naive = solve_with_MIP(S, f, fc, kb)
    m = solve_with_BDD_MIP(S, f, fc, kb)
    print(f"Naive model: status={m_naive.status}, obj={m_naive.objVal}")
    print(f"BDD model: status={m.status}, obj={m.objVal}")
    assert m_naive.objVal == m.objVal, f"Naive: {m_naive.objVal} (status {m_naive.status}), while BDD: {m.objVal} (status {m.status})"


def solve_with_align_BDD(S, f, fc, kb):
    """Solves the problem by generating two BDDs, aligning them, and solving a NF.

    Args:
        S (list): neighborhood,
        f (dict): location costs,
        fc (list): facility types,
        kb (list): type budgets.

    Returns:
        m (class gurobipy.Model): the resulting model.
    """
    C, cover_nl = build_cover_DD(S, f)
    pref_order = [int(x[1:]) for x in C.vars]
    A, type_nl = build_type_DD(f, fc, kb, pref_order)


    vs_A = vs.VarSeq(A.vars, [len(L) for L in A.layers[:-1]])
    vs_C = vs.VarSeq(C.vars, [len(L) for L in C.layers[:-1]])

    assert set(vs_A.layer_var) == set(vs_C.layer_var), f"A:{vs_A.layer_var}, C:{vs_C.layer_var}"
    b = bb.BBSearch(vs_A, vs_C)

    # bb.TIMEOUT_ITERATIONS=10000
    status = b.search()
    assert status == "optimal" or status == "timeout"

    Ap = A.align_to(b.Ap_cand.layer_var, inplace=False)
    Cp = C.align_to(b.Ap_cand.layer_var, inplace=False)

    int_DD = DD.intersect(Ap, Cp)
    m, c, v = UFL.create_NF(int_DD)
    m.setParam("OutputFlag", 0)
    m.optimize()
    return m

# Separate tesing code for TypeSorter

def no_invs(order, blocks, target):
    """Returns no. inversions b/w `blocks` in `order`, with the `target`.

    Examples:
        1) If I have a target of [1,2,3,4,5,6,7],
        but [1,2,3], [4,5,6], and [7] have to be "glued"
        in such blocks (numbered `0`, `1`, and `2`), then the order
        (2,0,1) will imply the element order [7,1,2,3,4,5,6] and 6 inversions:

        >>> no_invs((2,0,1), [[1,2,3], [4,5,6], [7]], [1,2,3,4,5,6,7])
        6

        Similarly,

        >>> no_invs((0,2,1), [[1,2,3],[4,5,6],[7]], [1,2,3,4,5,6,7])
        3

    2) The order of variables within a block, of course, does not matter:
    >>> no_invs((0,2,1), [[1,2,3],[6,5,4],[7]], [1,2,3,4,5,6,7])
    3

    """
    no_invs = 0
    pos = {target[i]: i for i in list(range(len(target)))}

    for first in range(len(order)-1):
        for second in range(first+1, len(order)):
            for a in blocks[order[first]]:
                for b in blocks[order[second]]:
                    if pos[a] > pos[b]:
                        no_invs += 1

    return no_invs


def bruteforce_correct_order(f_types, target_order):
    """Calculates the correct order with bruteforce."""
    min_inv = len(target_order)**2 + 1  # cannot be more than this
    best_order = None
    no_types = len(np.unique(f_types))

    types = [[] for c in range(no_types)]
    for e in target_order:
        types[f_types[e-1]].append(e)

    for order in permutations(list(range(len(types)))):
        invs = no_invs(order, types, target_order)
        if invs < min_inv:
            min_inv = invs
            best_order = order

    customers = {c: [] for c in best_order}

    for e in target_order:
        customers[f_types[e-1]].append(e)

    return sum([customers[c] for c in best_order], [])

def find_correct_order(f_types, target_order):
    """Finds the correct order with the code above."""

    cs = TypeSorter(f_types, target_order);
    types, customers = cs.sort_types()

    return sum([customers[c] for c in types], [])

def make_instance(n):
    """Makes a quick test instance for TypeSorter"""
    S, f, f_types, k_bar = generate_test_instance(n)
    target_order = np.random.permutation(list(range(1,len(f_types)+1)))
    return (f_types, target_order)

def get_score(A, target):
    """Calculates no. of inversions between `A` and `target`"""
    pos = {target[i]: i for i in list(range(len(target)))}
    invs = 0
    for i in range(len(A)-1):
        for j in range(i, len(A)):
            if pos[A[i]] > pos[A[j]]:
                invs += 1
    return invs

def show_cov(n, gen_func):
    """Generates an instance, generates and shows a *cover* diagram.

    Available generators:
    - generate_simple_problem
    - generate_test_instance
    - generate_string_instance
    - generate_organic_instance
    - generate_d4_instance
    - generate_5d_problem
    """
    S,f, fc, kb = gen_func(n)
    draw_problem_dia(S,f,fc,kb)
    c, nl = build_cover_DD(S,f)
    c.dump_gv(node_labels=nl).view()
