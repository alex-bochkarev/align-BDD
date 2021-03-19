"""Uncapacitated Facility Location with color-constraints.

Problem-solving machinery.
---
(c) A. Bochkarev, Clemson University, 2021
abochka@clemson.edu
"""
from copy import copy
from graphviz import Digraph
import numpy as np
import gurobipy as gp
from gurobipy import GRB  # noqa
import BDD as DD
import UFL

#######################################################################
# 1. Building the BDD representation


def build_overlap_DD(S, f, g):  # pylint: disable=too-many-locals,invalid-name
    """Builds a BDD for the UFL problem (with x-variables only).

    Args:
       S (list): neighborhood list,
       f (dict): location costs,
       g (dict): overlap costs.

    Returns:
        The resulting BDD.

    Notes:
        - employs a DP approach with state being the number of overlaps
            for each customer.
    """
    Sf = UFL.build_Sf(S)  # pylint: disable=invalid-name
    D = DD.BDD(N=len(Sf), vars=[i for i in range(1, len(Sf)+1)],  # noqa pylint: disable=invalid-name, unnecessary-comprehension
               weighted=True)

    # a *state* is the number of overlaps for each customer
    root_state = [0 for _ in range(len(S))]
    i = 1
    node_labels = dict({DD.NROOT: root_state})
    next_layer = {tuple(root_state): D.addnode(None)}

    while i < len(Sf):
        current_layer = copy(next_layer)
        next_layer = dict()
        for state in current_layer:
            node = current_layer[tuple(state)]
            if state in next_layer:
                D.llink(node, next_layer[state], "lo")
            else:
                newnode = D.addnode(node, "lo")
                next_layer.update({state: newnode})
                node_labels.update({newnode.id: str(state)})

            next_state = list(state)
            for j in Sf[i]:
                next_state[j-1] += 1
            next_state = tuple(next_state)

            if next_state in next_layer:
                D.llink(node, next_layer[next_state], "hi",
                        edge_weight=f[i])
            else:
                newnode = D.addnode(node, "hi", edge_weight=f[i])
                next_layer.update({next_state: newnode})
                node_labels.update({newnode.id: str(next_state)})
        i += 1

    for state in next_layer:
        next_state = list(state)
        for j in Sf[i]:
            next_state[j-1] += 1
        next_state = tuple(next_state)

        w_hi = sum(g[(j + 1, state[j] + 1)]
                          if (j + 1) in Sf[i] else g[(j + 1, state[j])]  # noqa pylint: disable=superfluous-parens
                          for j in range(len(state)))

        w_lo = sum(g[(j+1, state[j])]
                   for j in range(len(state)))

        node_id = next_layer[state].id
        if np.any(np.array(state) == 0):
            n_target = DD.NFALSE
            w_lo = 0.0
            if np.any(np.array(next_state) == 0):
                y_target = DD.NFALSE
                w_hi = 0.0
            else:
                y_target = DD.NTRUE
        else:
            y_target = DD.NTRUE
            n_target = DD.NTRUE

        D.link(node_id, y_target, "hi", edge_weight=w_hi)
        D.link(node_id, n_target, "lo", edge_weight=w_lo)

    return D, node_labels


def build_location_DD(f, f_color, k_bar):  # pylint: disable=invalid-name
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
    D = DD.BDD(N=len(f_color), vars=[i for i in range(1, len(f_color)+1)],  # noqa pylint: disable=invalid-name, unnecessary-comprehension
               weighted=True)
    # a *state* is the number of located facilities per color
    root_state = [0 for _ in range(len(k_bar))]
    i = 1
    node_labels = dict({DD.NROOT: root_state})
    next_layer = {tuple(root_state): D.addnode(None)}

    while i < len(f_color):
        current_layer = copy(next_layer)
        next_layer = dict()
        for state in current_layer:
            node = current_layer[tuple(state)]
            if state in next_layer:
                D.llink(node, next_layer[state], "lo")
            else:
                newnode = D.addnode(node, "lo")
                next_layer.update({state: newnode})
                node_labels.update({newnode.id: str(state)})

            next_state = list(state)
            next_state[f_color[i-1]] += 1
            next_state = tuple(next_state)

            if next_state in next_layer:
                D.llink(node, next_layer[next_state], "hi",
                        edge_weight=f[i])
            else:
                newnode = D.addnode(node, "hi", edge_weight=f[i])
                next_layer.update({next_state: newnode})
                node_labels.update({newnode.id: str(next_state)})
        i += 1

    for state in next_layer:
        node_id = next_layer[state].id
        if np.any(np.array(state) > k_bar):
            y_target = DD.NFALSE
            n_target = DD.NFALSE
        else:
            n_target = DD.NTRUE
            yes_state = list(state)
            yes_state[f_color[i-1]] += 1
            if np.any(np.array(yes_state) > k_bar):
                y_target = DD.NFALSE
            else:
                y_target = DD.NTRUE

        D.link(node_id, y_target, "hi", edge_weight=f[i])
        D.link(node_id, n_target, "lo")

    return D, node_labels


######################################################################
# 2. Building a MIP


def build_cUFL_MIP(S, f, f_color, k_bar, g):
    """Builds 'plain-vanilla' MIP instance for 'colorful' UFL.

    Encodes arbitrary penalty function g(·).

    Args:
        S (list): list of customer neighborhoods,
        f (dict): facility costs,
        f_color (list): facility colors,
        k_bar (list): max no. of locations per color,
        g (dict): values for overlap penalties, `(j, n)`
                    for `n` overlaps at consumer `j`

    Returns:
        `gurobipy` model (`gp.Model`)
    """
    m = gp.Model()
    m.modelSense = GRB.MINIMIZE

    Sf = UFL.build_Sf(S)

    x = dict()
    z = dict()
    b = dict()
    v = dict()

    F = range(1, len(Sf.keys())+1)
    C = range(1, len(S)+1)

    for i in F:
        x[i] = m.addVar(vtype=GRB.BINARY, name=f"x{i}")
        for j in Sf[i]:
            z[(i, j)] = m.addVar(vtype=GRB.BINARY, name=f"z{i}_{j}")
            m.addConstr(z[(i, j)] == x[i])

    for j in C:
        b[j] = m.addVar(name=f"b{j}")
        m.addConstr(b[j] == gp.quicksum(z[(i, j)] for i in S[j-1]))

    # set the Objective
    obj = gp.quicksum(f[i]*x[i] for i in F)

    # overlap penalty
    for j in C:
        for k in range(1, len(S[j-1])+1):
            v[(j, k)] = m.addVar(vtype=GRB.BINARY, name=f"v{j}_{k}")
            if k > 1:
                m.addConstr(v[(j, k)] <= v[(j, k-1)])

        m.addConstr(b[j] == gp.quicksum(v[(j, k)]
                    for k in range(1, len(S[j-1])+1)))

        obj += g[(j, 0)] + gp.quicksum((g[(j, k)] - g[(j, k-1)])*v[(j, k)]
                                       for k in range(1, len(S[j-1])+1))

    m.setObjective(obj)
    m.update()
    return m
######################################################################
# 3. Quick-testing code


def make_simple_problem():
    """Generates a simple problem instance."""
    S = [[1], [1, 2], [1, 2, 3], [2, 3]]
    f = {1: 1, 2: 2, 3: 3}
    g = {
        (1, 0): 11,  (2, 0): 12,  (3, 0): 13,  (4, 0): 14,
        (1, 1): 0,  (2, 1): 0,  (3, 1): 0,  (4, 1): 0,
        (1, 2): 2.1,  (2, 2): 2.2,  (3, 2): 2.3,  (4, 2): 2.4,
        (3, 3): 3.1}
    f_colors = [0, 1, 0]
    k_bar = [1, 1, 1]
    return S, g, f, f_colors, k_bar

def draw_problem_dia(S, g, f, f_colors, k_bar, filename="run_logs/c_problem_dia.gv"):
    """Draws the bipartite graph describing the problem."""
    n = len(UFL.build_Sf(S).keys())
    cols = {0: 'red', 1: 'blue', 2: 'green', 3: 'yellow', 4: 'black', 5: 'purple'}

    dia = Digraph('G', comment="Uncapacitated Facility Location")
    for i in range(n):
        dia.node(f"F{i+1}", shape='doublecircle',
                 style='filled', color=cols[f_colors[i]],
                 label=f"F{i+1} ({f[i+1]}) /c{f_colors[i]}")

    for j in range(len(S)):
        dia.node(f"C{j+1}", shape='circle',
                 label=f"C{j+1}\ng=({','.join([str(g[(k, l)]) for (k, l) in g.keys() if k==(j+1)])})")  # noqa: E501

    for j, S_j in enumerate(S):
        for i in S_j:
            dia.edge(f"F{i}", f"C{j+1}")

    print("Color limits: " + ", ".join([f"{cols[c]}: {k_bar[c]}" for c in range(len(k_bar))]))
    dia.view(filename=filename)


def show_dias():
    S, g, f, f_color, K = make_simple_problem()
    oD, nl_oD = build_overlap_DD(S, f, g)
    lD, nl_lD = build_location_DD(f, f_color, K)

    oD.dump_gv(node_labels=nl_oD).view(filename="overlap")
    lD.dump_gv(node_labels=nl_lD).view(filename="location")
    draw_problem_dia(S, g, f, f_color, K)
