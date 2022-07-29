"""Implements basic functions for UFL -- Uncapacitated Facility Location.

Used as a basis for :py:mod:`tUFLP` (some foundational UFLP-related code).

Tests coverage: :py:mod:`UFLP_test`
"""
from graphviz import Digraph
import numpy as np
import gurobipy as gp
from gurobipy import GRB
import sys
from time import time
from copy import copy

import varseq as vs
import BB_search as bb
import BDD as DD


def draw_problem_dia(S, f, g, filename="run_logs/problem_dia.gv"):
    """Draws the bipartite graph describing the problem."""
    n = len(build_Sf(S).keys())

    dia = Digraph('G', comment="Uncapacitated Facility Location")
    for i in range(n):
        dia.node(f"F{i+1}", shape='doublecircle',
                 style='filled', color='lightgrey',
                 label=f"F{i+1} ({f[i+1]})")

    for j in range(len(S)):
        dia.node(f"C{j+1}", shape='circle',
                 label=f"C{j+1}\ng=({','.join([str(g[(k, l)]) for (k, l) in g.keys() if k==(j+1)])})")  # noqa: E501

    for j, S_j in enumerate(S):
        for i in S_j:
            dia.edge(f"F{i}", f"C{j+1}")

    dia.view(filename=filename)


def create_covering_BDD_wg(S, g):  # pylint: disable=all
    """Creates and returns "covering" BDD given an UFL instance.

    Args:
        S (list): of adjacency lists
        g (dict): overlapping costs, g(k,n) is a cost of `n`
            overlaps for consumer `k`)

    Returns:
        C (class BDD): covering BDD
    """
    m = len(S)

    var_names = [f"z{i}-{j+1}" for j in range(m) for i in S[j]]
    C = DD.BDD(
        N=len(var_names),  # number of z-variables
        vars=var_names,
        weighted=True)

    root = C.addnode(None)
    current_layer = [root]

    for j in range(1, m+1):
        assert len(S[j-1]) > 0

        for i in range(1, len(S[j-1])):  # all *but* the last nearby facility
            next_layer = [C.addnode(n, "lo") for n in current_layer] + \
                [C.addnode(current_layer[-1], "hi")]

            for q in range(len(current_layer)-1):
                C.link(current_layer[q].id, next_layer[q+1].id, "hi")

            current_layer = next_layer

        # process the last layer of the subnetwork corresponding
        # to consumer `j`: this is always a special case
        if j < m:
            target = C.addnode(current_layer[0], "lo")
        else:
            target = C.T

        for idx, n in enumerate(current_layer):
            C.link(n.id, target.id, "hi",
                   edge_weight=g[(j, idx+1)])
            C.link(n.id, target.id, "lo",
                   edge_weight=g[(j, idx)])

        current_layer = [target]

    return C


def create_covering_BDD(S, c):
    """Create and return the diagram given the instance.

    Args:
        S (list): of adjacency lists. len(S) = no. of consumers
        c (dict): costs of covering, keyed by (facility_id, consumer_id)

    Notes:
        - encodes the fact that all customers need to be covered
            by at least one facility.
        - assumes no edge costs, does *not* allow
            for arbitrary g(·) function.
    """
    m = len(S)

    C = DD.BDD(
        N=len(np.sum(S)),  # number of z-variables
        vars=[f"z{i}-{j+1}" for j in range(m) for i in S[j]],
        weighted=True
    )

    root = C.addnode(None)

    # current nodes
    s_maybe = root
    s_no = None

    for j in range(m):
        assert len(S[j]) > 0

        s_yes = None
        if j == m-1:
            S_j = S[j][:-1]
        else:
            S_j = S[j]

        for idx, i in enumerate(S_j):
            new_maybe = C.addnode(s_maybe, "lo")
            new_yes = C.addnode(s_maybe, "hi", edge_weight=c[(i, j+1)])

            if s_yes is not None:
                C.link(s_yes.id, new_yes.id, "hi", edge_weight=c[(i, j+1)])
                C.link(s_yes.id, new_yes.id, "lo")

            if s_no is not None:
                if idx < len(S[j])-1:
                    new_no = C.addnode(s_no, edge_weight=c[(i, j+1)])
                else:
                    new_no = new_maybe
                    C.link(s_no.id, new_no.id, "hi", edge_weight=c[(i, j+1)])

                C.link(s_no.id, new_no.id, "lo")
                s_no = new_no

            s_yes = new_yes

            if idx == len(S[j])-1:
                # the last layer of the subnetwork
                # (related to customer j)
                s_no = new_maybe
                s_maybe = new_yes
                C.link(s_yes.id, new_yes.id, "hi", edge_weight=c[(i, j+1)])
                C.link(s_yes.id, new_yes.id, "lo")
            else:
                s_maybe = new_maybe

    C.link(s_maybe.id, DD.NTRUE, "hi", edge_weight=c[(i, j+1)])
    C.link(s_maybe.id, DD.NFALSE, "lo")

    if not (s_yes is None) and (s_yes.id != s_maybe.id):
        C.link(s_yes.id, DD.NTRUE, "hi", edge_weight=c[(i, j+1)])
        C.link(s_yes.id, DD.NTRUE, "lo")

    if s_no is not None:
        C.link(s_no.id, DD.NFALSE, "hi", edge_weight=c[(i, j+1)])
        C.link(s_no.id, DD.NFALSE, "lo")

    return C


def build_Sf(S):
    """Returns the list of facility neighborhoods.

    Args:
        S (list): of consumer neighborhoods (facility IDs)
    """
    Sf = dict()
    for j in range(len(S)):
        for facility in S[j]:
            if facility in Sf.keys():
                Sf[facility].add(j+1)
            else:
                Sf.update({facility: set([j+1])})

    return Sf


def create_availability_BDD(S, f):
    """Encodes consistency of "cover"(z) decisions.

    Either all "yes", or all "no" for a specific facility.

    Args:
        S (list): neighborhood list
        f (dict): facility costs
    """
    Sf = build_Sf(S)
    n = len(Sf.keys())
    F = range(1, n+1)

    var_names = [f"z{facility}-{j}"
                 for facility in F for j in Sf[facility]]

    A = DD.BDD(N=len(var_names), vars=var_names, weighted=True)

    s_main = A.addnode(None)
    s_false = None

    for i in F:
        s_yes = None
        s_no = None
        y_weight = f[i]
        consumers_i = list(Sf[i])

        for j in consumers_i:
            if i == n and j == consumers_i[-1]:
                # this is the very last layer
                if s_yes is not None:
                    A.link(s_yes.id, DD.NTRUE, "hi", edge_weight=y_weight)
                    A.link(s_yes.id, DD.NFALSE, "lo")

                if s_no is not None:
                    A.link(s_no.id, DD.NTRUE, "lo")
                    A.link(s_no.id, DD.NFALSE, "hi", edge_weight=y_weight)

                if s_main is not None:
                    A.link(s_main.id, DD.NTRUE, "hi", edge_weight=y_weight)
                    A.link(s_main.id, DD.NTRUE, "lo")

                if s_false is not None:
                    A.link(s_false.id, DD.NFALSE, "hi")
                    A.link(s_false.id, DD.NFALSE, "lo")
            else:
                # this is not the very last layer
                # ( so we need to create a layer )
                if s_main is not None:
                    # this is the first decision for a facility
                    new_yes = A.addnode(s_main, "hi", edge_weight=y_weight)
                else:
                    new_yes = A.addnode(s_yes, "hi", edge_weight=y_weight)

                if j == consumers_i[-1]:
                    # this is the last decision for a facility
                    new_no = new_yes
                    if s_main is not None:
                        A.link(s_main.id, new_no.id, "lo")
                    s_main = new_yes
                else:
                    if s_main is not None:
                        new_no = A.addnode(s_main, "lo")
                        s_main = None
                    else:
                        new_no = A.addnode(s_no, "lo")

                if s_false is None and s_yes is not None:
                    new_f = A.addnode(s_yes, "lo")
                elif s_false is not None:
                    new_f = A.addnode(s_false, "hi")
                    A.link(s_false.id, new_f.id, "lo")
                else:
                    new_f = None

                if s_yes is not None:
                    A.link(s_yes.id, new_yes.id, "hi", edge_weight=y_weight)
                    A.link(s_yes.id, new_f.id, "lo")

                if s_no is not None:
                    A.link(s_no.id, new_no.id, "lo")
                    A.link(s_no.id, new_f.id, "hi")

                s_yes = new_yes
                s_no = new_no
                s_false = new_f

            y_weight = 0

    return A


# testing
def test_DD_creation(S, f, c, name):
    """Testing unit: BDD generation."""
    # unpack the problem from S
    draw_problem_dia(S, f, c, filename="run_logs/"+name+"_problem.gv")

    create_covering_BDD(S, c).dump_gv().view(filename="run_logs/" +
                                             name+"_covering.gv")
    create_availability_BDD(S, f).dump_gv().view(filename="run_logs/" +
                                                 name+"_avail.gv")


def build_MIP(S, f, g):
    """Builds 'plain-vanilla' MIP instance for UFL.

    Encodes arbitrary penalty function g(·).

    Args:
        S (list): list of customer neighborhoods
        f (dict): facility costs
        g (dict): values for overlap penalties, `(j, n)`
                    for `n` overlaps at consumer `j`

    Returns:
        `gurobipy` model (`gp.Model`)
    """
    m = gp.Model()
    m.modelSense = GRB.MINIMIZE

    Sf = build_Sf(S)

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


def add_BDD_to_MIP(D, model=None, x=None, prefix=""):
    """Builds a MIP (network flow) from a given BDD.

    Args:
        D (class BDD): the diagram object
        model (gurobipy Model):
            a model to amend with vars and constraints(default None)
        x (dict): linking binary variables, if available (default None)
        prefix (str): a string prefix for variable names

    Returns:
        m (gurobipy Model): resulting network flow
        c (list): list of constraints added
        v (dict): flow variables added, keyed as (from_id, to_id, arc_type)
        x (dict): binary variables added (one per layer), keyed by _layer_no_.
    """
    if model is None:
        model = gp.Model()
        model.modelSense = GRB.MINIMIZE

    constraints = []
    v = dict()
    if x is None:
        x = dict()
    incoming_flows = dict()

    l_no = 0

    for L in D.layers:
        l_no += 1

        if l_no <= len(D):
            if not D.vars[l_no-1] in x.keys():
                x[D.vars[l_no-1]] = model.addVar(vtype=GRB.BINARY,
                                                 name=f"link_{D.vars[l_no-1]}")
            yes_sum = gp.LinExpr(0.0)

        inc_nodes_new = dict()
        for n in L:
            # determine node inflow
            inflow = gp.LinExpr(0.0)
            if n.id == DD.NROOT:
                inflow = gp.LinExpr(1.0)
            else:
                if n.id in incoming_flows.keys():
                    inflow = gp.quicksum(v_i for v_i in incoming_flows[n.id])

            # determine node outflow and add a linking constraint
            if n.id == DD.NTRUE:
                outflow = gp.LinExpr(1.0)

            elif n.id == DD.NFALSE:
                outflow = gp.LinExpr(0.0)
            else:
                v[(n.id, n.hi.id, "hi")] = model.addVar(
                    name=f"{prefix}vh{n.id}_{n.hi.id}",
                    obj=D.weights[(n.id, n.hi.id, "hi")])
                v[(n.id, n.lo.id, "lo")] = model.addVar(
                    name=f"{prefix}vl{n.id}_{n.lo.id}",
                    obj=D.weights[(n.id, n.lo.id, "lo")])
                outflow = gp.LinExpr(v[(n.id, n.hi.id, "hi")] +
                                     v[(n.id, n.lo.id, "lo")])
                yes_sum += v[(n.id, n.hi.id, "hi")]

            constraints.append(model.addConstr(inflow == outflow,
                                               name=f"{prefix}cont-at-{n.id}"))

            if n.id != DD.NTRUE and n.id != DD.NFALSE:
                if n.hi.id in inc_nodes_new.keys():
                    inc_nodes_new[n.hi.id].append(v[(n.id, n.hi.id, "hi")])
                else:
                    inc_nodes_new[n.hi.id] = [v[(n.id, n.hi.id, "hi")]]

            if n.id != DD.NTRUE and n.id != DD.NFALSE:
                if n.lo.id in inc_nodes_new.keys():
                    inc_nodes_new[n.lo.id].append(v[(n.id, n.lo.id, "lo")])
                else:
                    inc_nodes_new[n.lo.id] = [v[(n.id, n.lo.id, "lo")]]

        incoming_flows = inc_nodes_new

        if l_no <= len(D):
            constraints.append(model.addConstr(
                x[D.vars[l_no-1]] == yes_sum, name=f"{prefix}bin-link-{l_no}"))

    model.update()
    return model, constraints, v, x


def create_NF(D):
    """Generates a network flow instance for diagram `D`.

    Args:
       D (class BDD): the underlying diagram.

    Returns:
        model (`gurobipy` model): resulting model,
        constraints (list),
        v(dict): variables.
    """
    model = gp.Model()
    model.modelSense = GRB.MINIMIZE

    constraints = []
    v = dict()
    incoming_flows = dict()

    l_no = 0

    for L in D.layers:
        l_no += 1

        inc_nodes_new = dict()
        for n in L:
            # determine node inflow
            inflow = gp.LinExpr(0.0)
            if n.id == DD.NROOT:
                inflow = gp.LinExpr(1.0)
            else:
                if n.id in incoming_flows.keys():
                    inflow = gp.quicksum(v_i for v_i in incoming_flows[n.id])

            # determine node outflow and add a linking constraint
            if n.id == DD.NTRUE:
                outflow = gp.LinExpr(1.0)

            elif n.id == DD.NFALSE:
                outflow = gp.LinExpr(0.0)
            else:
                v[(n.id, n.hi.id, "hi")] = model.addVar(
                    name=f"vh{n.id}_{n.hi.id}",
                    obj=D.weights[(n.id, n.hi.id, "hi")])
                v[(n.id, n.lo.id, "lo")] = model.addVar(
                    name=f"vl{n.id}_{n.lo.id}",
                    obj=D.weights[(n.id, n.lo.id, "lo")])
                outflow = gp.LinExpr(v[(n.id, n.hi.id, "hi")] +
                                     v[(n.id, n.lo.id, "lo")])

            constraints.append(model.addConstr(inflow == outflow,
                                               name=f"cont-at-{n.id}"))

            if n.id != DD.NTRUE and n.id != DD.NFALSE:
                if n.hi.id in inc_nodes_new.keys():
                    inc_nodes_new[n.hi.id].append(v[(n.id, n.hi.id, "hi")])
                else:
                    inc_nodes_new[n.hi.id] = [v[(n.id, n.hi.id, "hi")]]

            if n.id != DD.NTRUE and n.id != DD.NFALSE:
                if n.lo.id in inc_nodes_new.keys():
                    inc_nodes_new[n.lo.id].append(v[(n.id, n.lo.id, "lo")])
                else:
                    inc_nodes_new[n.lo.id] = [v[(n.id, n.lo.id, "lo")]]

        incoming_flows = inc_nodes_new

    model.update()
    return model, constraints, v


def solve_with_intBDD(S, f, g):
    """Solves the UFL instance by aligning two BDDs.

    Args:
       S (list): neighborhood list
       f (dict): facility location costs
       g (dict): overlap penalty dict.

    Returns:
        Resulting model, constraints, and variables.
    """
    A = create_availability_BDD(S, f)
    C = create_covering_BDD_wg(S, g)
    vs_A = vs.VarSeq(A.vars, [len(L) for L in A.layers[:-1]])
    vs_C = vs.VarSeq(C.vars, [len(L) for L in C.layers[:-1]])

    assert set(vs_A.layer_var) == set(vs_C.layer_var)
    b = bb.BBSearch(vs_A, vs_C)

    # bb.TIMEOUT_ITERATIONS=10000
    status = b.search()
    assert status == "optimal" or status == "timeout"

    Ap = A.align_to(b.Ap_cand.layer_var, inplace=False)
    Cp = C.align_to(b.Ap_cand.layer_var, inplace=False)

    int_DD = DD.intersect(Ap, Cp)
    m, c, v = create_NF(int_DD)
    m.setParam("OutputFlag", 0)
    m.optimize()
    return m, c, v


def build_DP_DD(S, f, g):
    """Builds a BDD for the UFL problem (with x-variables only)

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
    Sf = build_Sf(S)
    D = DD.BDD(N=len(Sf), vars=[i for i in range(1, len(Sf)+1)],
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
        w_hi = f[i] + sum(g[(j + 1, state[j] + 1)]
                          if (j + 1) in Sf[i] else g[(j + 1, state[j])]
                          for j in range(len(state)))

        w_lo = sum(g[(j+1, state[j])]
                   for j in range(len(state)))

        node_id = next_layer[state].id
        D.link(node_id, DD.NTRUE, "hi", edge_weight=w_hi)
        D.link(node_id, DD.NTRUE, "lo", edge_weight=w_lo)

    return D, node_labels
# =====================================================================
# Testing code
# =====================================================================

def make_simple_problem():
    """Creates a 'toy' UFLP instance (for debugging)."""
    S = [[1], [1, 2], [1, 2, 3], [2, 3]]
    f = {1: 1, 2: 2, 3: 3}
    g = {
        (1, 0): 11,  (2, 0): 12,  (3, 0): 13,  (4, 0): 14,
        (1, 1): 0,  (2, 1): 0,  (3, 1): 0,  (4, 1): 0,
        (1, 2): 2.1,  (2, 2): 2.2,  (3, 2): 2.3,  (4, 2): 2.4,
        (3, 3): 3.1}
    return S, f, g


def show_BDD_to_MIP():
    """Tests making a MIP from a single BDD."""
    # adhoc experiments
    A = DD.BDD()
    A.load("tests/simple_DD.bdd")
    A.show(filename="A.dot")

    B = DD.BDD()
    B.load("tests/simple_BDD_different_vars.bdd")
    B.show(filename="B.dot")

    m, c, v, x = add_BDD_to_MIP(A, prefix="A_")
    m, c, v, x = add_BDD_to_MIP(D=B, model=m, x=x, prefix="B_")
    m.display()


def show_build_MIP():
    """Tests a simple MIP building function."""
    S = [[1], [1, 2], [1, 2, 3], [2, 3]]
    f = {1: 0.1, 2: 0.2, 3: 0.3}
    g = {
        (1, 0): 1.0,  (2, 0): 2.0,  (3, 0): 3.0,  (4, 0): 4.0,
        (1, 1): 1.1,  (2, 1): 2.1,  (3, 1): 3.1,  (4, 1): 4.1,
        (1, 2): 1.21,  (2, 2): 2.22,  (3, 2): 3.23,  (4, 2): 4.24,
        (3, 3): 3.33}

    print(f"S={S},\n f={f},\n g={g}")
    print("===")
    print("The resulting MIP is:")
    m = build_MIP(S, f, g)
    m.display()


def show_BDD_build():
    """Runs a simple test against a toy problem."""
    S = [[1], [1, 2], [1, 2], [2]]
    f = {1: 0.5, 2: 0.7}
    c = {(1, 1): 0.11,  (1, 2): 0.12,  (2, 2): 0.22,  (1, 3): 0.13,
         (2, 3): 0.23,  (2, 4): 0.24}
    g = {
        (1, 0): 0,  (2, 0): 0,  (3, 0): 0,  (4, 0): 0,
        (1, 1): 1,  (2, 1): 1,  (3, 1): 1,  (4, 1): 1,
        (1, 2): 2,  (2, 2): 2,  (3, 2): 2,  (4, 2): 2}

    draw_problem_dia(S, f, c, g)
    C = create_covering_BDD_wg(S, g)
    C.show(x_prefix='', filename="covering.dot", dir="run_logs")
    A = create_availability_BDD(S, f)
    A.show(x_prefix='', filename="availability.dot", dir="run_logs")


def show_BDD_to_MIP_wg(S, f, g):
    """Runs a simple test against a toy problem."""
    draw_problem_dia(S, f, g)
    C = create_covering_BDD_wg(S, g)
    C.show(x_prefix='', filename="covering.dot", dir="run_logs")
    A = create_availability_BDD(S, f)
    A.show(x_prefix='', filename="availability.dot", dir="run_logs")

    m, c, v, x = add_BDD_to_MIP(A, prefix="A_")
    m, c, v, x = add_BDD_to_MIP(D=C, model=m, x=x, prefix="C_")
    m.display()


def generate_test_figures():
    """Generates a simple test instance and creates PDFs."""
    S = [[1, 2], [1, 2, 3], [2, 3, 4], [2, 4]]
    f = {1: 0.1, 2: 0.2, 3: 0.3, 4: 0.4}
    c = {(1, 1): 0.11, (1, 2): 0.12,
         (2, 1): 0.21, (2, 2): 0.22, (2, 3): 0.23, (2, 4): 0.24,
         (3, 2): 0.32, (3, 3): 0.33, (3, 4): 0.34,
         (4, 3): 0.43, (4, 4): 0.44}
    g = {
        (1, 0): 0,  (2, 0): 0,  (3, 0): 0,  (4, 0): 0,
        (1, 1): 1,  (2, 1): 1,  (3, 1): 1,  (4, 1): 1,
        (1, 3): 1,  (2, 3): 1,  (3, 3): 1,  (4, 3): 1,
        (1, 2): 2,  (2, 2): 2,  (3, 2): 2,  (4, 2): 2}

    draw_problem_dia(S, f, c, g)

    C = create_covering_BDD_wg(S, g)
    C.show(x_prefix='', filename="tst_covering.dot", dir="run_logs")
    A = create_availability_BDD(S, f)
    A.show(x_prefix='', filename="tst_availability.dot", dir="run_logs")

    m, c, v, x = add_BDD_to_MIP(A, prefix="A_")
    m, c, v, x = add_BDD_to_MIP(D=C, model=m, x=x, prefix="C_")
    m.display()


def generate_test_instance(n, m):
    """Generates a test facility location instance.

    Args:
       n (int): number of facilities
       m (int): number of customers

    Returns:
        A tuple of the following values.

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

        Sf = build_Sf(S)
        if np.sum([1 for i in N if i not in Sf]) == 0:
            good_instance = True

    return S, f, g


def generate_dense_instance(n, m, covering=0.95):
    """Generates a dense test facility location instance.

    In the sense that every customer can be covered by (almost)
    every facility.

    Args:
        n (int): number of facilities
        m (int): number of customers
        covering (float): share of facilities covering each customer

    Returns:
        A tuple of the following values.

            - S (list): neighborhood list
            - f (dict): costs of facility location
              (generated uniformly random int from `f_min` to `f_max`)
            - g (dict): overlap costs, keys: (customer, overlap)
    """
    N = [i for i in range(1, n+1)]
    M = [j for j in range(1, m+1)]

    good_instance = False
    while not good_instance:
        S = [np.random.choice(N, round(n*covering),
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

        Sf = build_Sf(S)
        if np.sum([1 for i in N if i not in Sf]) == 0:
            good_instance = True

    return S, f, g


def test_MIPs_protocol():
    """Runs a series of cross-checks.

    Uses :py:func:`generate_test_instance`.

    Covers functions:
        - build_MIP
        - create_covering_BDD_wg
        - create_availability_BDD
        - add_BDD_to_MIP
    """
    test_protocol = [(1000, 10, 20), (1000, 20, 10),
                     (500, 20, 20), (500, 50, 100),
                     (100, 100, 50), (100, 100, 100)]

    for test in test_protocol:
        print(f"Testing {test[0]} instances with n={test[1]}, m={test[2]}")
        test_BDD_and_plain_MIPs(K=test[0], n=test[1], m=test[2])

    print("Test finished.")


def test_BDD_and_plain_MIPs(K=500, TOL=1e-3, n=3, m=4):
    """Tests that objectives for the two MIPs coincide.

    Args:
        K (int): number of instances to generate (default 500)
        TOL (float): objective tolerance (for comparison) (default 1e-3)
        n (int): number of facilities (default 3)
        m (type): number of customers (default 4)

    Returns:
        Nothing (outputs the success rate to the screen)
    """
    # define a problem
    no_success = 0

    plain_MIP_time = 0
    CPP_MIP_time = 0

    for k in range(K):
        S, f, g = generate_test_instance(n=n, m=m)

        # Generate and solve plain MIP
        t0 = time()
        model = build_MIP(S, f, g)
        model.setParam("OutputFlag", 0)
        model.optimize()
        t1 = time()
        plain_MIP_time += (t1 - t0)

        if model.status != GRB.OPTIMAL:
            print(f"Plain MIP status is: {model.status}")
            print(f"\nS={S}; f={f}; g={g}")
            continue
        plain_MIP_obj = model.objVal

        # Generate and solve CPP MIP
        t0 = time()
        C = create_covering_BDD_wg(S, g)
        A = create_availability_BDD(S, f)

        model, c, v, x = add_BDD_to_MIP(A, prefix="A_")
        model, c, v, x = add_BDD_to_MIP(C, model=model, x=x, prefix="C_")
        model.update()
        model.setParam("OutputFlag", 0)
        model.optimize()
        t1 = time()
        CPP_MIP_time += (t1 - t0)
        if model.status != GRB.OPTIMAL:
            print(f"CPP MIP status is: {model.status}")
            print(f"\nS={S}; f={f}; g={g}")
            continue

        CPP_MIP_obj = model.objVal
        if abs(CPP_MIP_obj - plain_MIP_obj) < TOL:
            no_success += 1
            print(".", end="")
        else:
            print("!")
            print(f"\nS={S}; f={f}; g={g}")

        sys.stdout.flush()

    print(f"Testing finished. Success rate: {no_success*100/K:.1f}%")
    print(f"Mean times (per instance), construction + solving:")
    print(f"Plain MIP: {plain_MIP_time / K:.1f} sec.")
    print(f"CPP MIP: {CPP_MIP_time / K:.1f} sec.")

    # benchmark(K=int(args.K), n=int(args.n), m=int(args.m), prefix=args.prefix)
    # Procedure that is run if executed from the command line

    # S = [[1], [1, 2], [1, 2, 3], [2, 3]]
    # f = {1: 0.1, 2: 0.2, 3: 0.3}
    # g = {
    #     (1, 0): 0,  (2, 0): 0,  (3, 0): 0,  (4, 0): 0,
    #     (1, 1): 1,  (2, 1): 1,  (3, 1): 1,  (4, 1): 1,
    #     (1, 2): 2,  (2, 2): 2,  (3, 2): 2,  (4, 2): 2,
    #     (3, 3): 3}

    # m = build_MIP(S, f, c, g)
    # print("*==============================================================*")
    # print("* A MIP for the problem:                                       *")
    # print("*==============================================================*")
    # m.display()

    # test_MIPs_protocol()
    # test_build_MIP()
    # test_BDD_to_MIP_wg(S, f, g)
    # generate_test_figures()
    # test_BDD_to_MIP()
