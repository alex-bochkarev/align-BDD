"""
UFL -- Uncapacitated Facility Location
(testing align-BDD machinery for specific applications)

General notes:
--------------
- on indices: `j` is used for consumers, `i` -- for facilities
---
(c) A. Bochkarev, Clemson University, 2021
abochka@clemson.edu
"""
from graphviz import Digraph
import numpy as np
import gurobipy as gp
from gurobipy import GRB

import BDD as DD


def draw_problem_dia(S, f, filename="run_logs/problem_dia.gv"):
    """Draws the bipartite graph describing the probelm"""

    n = len(build_Sf(S).keys())

    dia = Digraph('G', comment="Uncapacitated Facility Location")
    for i in range(n):
        dia.node(f"F{i+1}", shape='doublecircle',
                 style='filled', color='lightgrey',
                 label=f"F{i+1} ({f[i+1]})")

    for j, S_j in enumerate(S):
        dia.node(f"C{j+1}", shape='circle')
        for i in S_j:
            dia.edge(f"F{i}", f"C{j+1}")

    dia.view(filename=filename)



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
        N=len(np.sum(S)), # number of z-variables
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

    if not s_no is None:
        C.link(s_no.id, DD.NFALSE, "hi", edge_weight=c[(i, j+1)])
        C.link(s_no.id, DD.NFALSE, "lo")

    return C


def build_Sf(S):
    """returns the list of facility neighborhoods

    Arguments
    ---------
    S: list -- of consumer neighborhoods (facility IDs)
    """

    Sf = dict()
    for j in range(len(S)):
        for facility in S[j]:
            if facility in Sf.keys():
                Sf[facility].add(j+1)
            else:
                Sf.update({facility : set([j+1])})

    return Sf



def create_availability_BDD(S, f):
    """Encodes the fact that "turn-on" (x) and "cover"(z) decisions
    need to be consistent.
    Args:
        f (dict): facility costs
    """

    Sf = build_Sf(S)
    n = len(Sf.keys())

    var_names = sum([[f"x{facility}"] + [f"z{facility}-{j}"
                                            for j in Sf[facility]]
                        for facility in Sf],[])

    A = DD.BDD(N=len(var_names), vars=var_names, weighted=True)

    s_main = A.addnode(None);

    s_false = None
    for i in range(1,n+1):
        # create a block per facility
        s_yes = A.addnode(s_main, "hi", edge_weight=f[i])
        s_no = A.addnode(s_main, "lo")

        if s_false is not None:
            new_f = A.addnode(s_false, "hi", edge_weight=f[i])
            A.link(s_false.id, new_f.id, "lo")
            s_false = new_f

        consumers_i = list(Sf[i])
        for j in consumers_i[:-1]:
            new_yes = A.addnode(s_yes, "hi")
            new_no = A.addnode(s_no, "lo")

            new_f = A.addnode(s_yes, "lo")
            A.link(s_no.id, new_f.id, "hi")

            if s_false is not None:
               A.link(s_false.id, new_f.id, "hi")
               A.link(s_false.id, new_f.id, "lo")

            s_yes = new_yes
            s_no = new_no
            s_false = new_f

        if i < n:
            s_main = A.addnode(s_yes, "hi")
            A.link(s_no.id, s_main.id, "lo")

            new_f = A.addnode(s_false, "hi")
            A.link(s_false.id, new_f.id, "lo")
            A.link(s_yes.id, new_f.id, "lo")
            A.link(s_no.id, new_f.id, "hi")
            s_false = new_f
        else:
            # this is a special case:
            # we'd need to link to the terminal nodes
            A.link(s_yes.id,DD.NTRUE, "hi");
            A.link(s_yes.id,DD.NFALSE, "lo")
            A.link(s_no.id,DD.NTRUE, "lo")
            A.link(s_no.id,DD.NFALSE, "hi")
            A.link(s_false.id,DD.NFALSE,"hi")
            A.link(s_false.id,DD.NFALSE,"lo")

    return A


# testing
def test_DD_creation(S, f, c, name):
    """Testing unit: BDD generation"""

    # unpack the problem from S
    draw_problem_dia(S, f, filename="run_logs/"+name+"_problem.gv")

    create_covering_BDD(S, c).dump_gv().view(filename="run_logs/"+name+"_covering.gv")
    create_availability_BDD(S, f).dump_gv().view(filename="run_logs/"+name+"_avail.gv")


def build_MIP(S, f, c, λ):
    """Builds 'plain-vanilla' MIP instance for UFL
    (linear overlap penalty)

    Args:
    ----------
    S (list): -- list of customer neighborhoods
    f (dict): -- facility costs
    c (dict): -- covering costs ((i,j) for covering customer j with facility i)
    λ (num): -- coefficient for overlapping costs.

    Returns:
    ---------
    Gurobipy model (`gp.Model`)
    """

    m = gp.Model()
    m.modelSense = GRB.MINIMIZE

    Sf = build_Sf(S)

    x = dict()
    z = dict()
    b = dict()

    F = range(1, len(Sf.keys())+1)
    C = range(1, len(S)+1)
    for i in F:
        x[i] = m.addVar(vtype=GRB.BINARY, name=f"x{i}")
        for j in Sf[i]:
            z[(i, j)] = m.addVar(vtype=GRB.BINARY, name=f"z{i}_{j}")
            m.addConstr(z[(i, j)] <= x[i])

    for j in C:
        b[j] = m.addVar(name=f"b{j}")
        m.addConstr(b[j] == gp.quicksum(z[(i, j)] for i in S[j-1]))
        m.addConstr(b[j] >= 1)

    # set the Objective
    obj = gp.quicksum(f[i]*x[i] for i in F)
    obj += gp.quicksum(c[(i, j)]*z[(i, j)] for j in C for i in S[j-1])

    # linear overlap penalty:
    obj += λ*gp.quicksum(b[j] for j in C)
    obj += - λ * len(S)

    m.setObjective(obj)
    m.update()
    return m


def add_BDD_to_MIP(D, model=None, x=None, prefix=""):
    """Builds a MIP (network flow) from a given BDD

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
                x[D.vars[l_no-1]] = model.addVar(vtype=GRB.BINARY, name=f"x{D.vars[l_no-1]}")
            yes_sum = gp.LinExpr(0.0)

        inc_nodes_new = dict()
        for n in L:

            if n.id == DD.NROOT:
                inflow = gp.LinExpr(1.0)
            else:
                inflow = gp.quicksum(v_i for v_i in incoming_flows[n.id])

            if n.id == DD.NTRUE:
                outflow = gp.LinExpr(1.0)
            elif n.id == DD.NFALSE:
                outflow = gp.LinExpr(0.0)
            else:
                v[(n.id, n.hi.id, "hi")] = model.addVar(name=f"{prefix}vh{n.id}_{n.hi.id}")
                v[(n.id, n.lo.id, "lo")] = model.addVar(name=f"{prefix}vl{n.id}_{n.lo.id}")
                outflow = gp.LinExpr(v[(n.id, n.hi.id, "hi")] + v[(n.id, n.lo.id, "lo")])
                yes_sum += v[(n.id, n.hi.id, "hi")]

            constraints.append(model.addConstr(inflow == outflow, name=f"{prefix}cont-at-{n.id}"))

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

def test_BDD_to_MIP():
    """Tests making a MIP from a single BDD"""

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

def test_build_simple_MIP():
    """Tests a simple MIP building function"""
    S = [[1], [1,2], [1,2], [2]]
    f = { 1: 0.5, 2:0.7 }
    c = {(1,1):1.1, (1,2):1.2, (2,2):2.2, (1,3):1.3, (2,3):2.3, (2,4):2.4 }
    λ = 0.5

    test_DD_creation(S, f, c, "toy")

    m = build_MIP(S, f, c, λ)
    m.display()


if __name__ == '__main__':
    """Runs when called from a command line"""
    # Procedure that is run if executed from the command line

    # define a toy problem
    # S = [[1], [1,2], [1,2], [2]]
    # test_DD_creation(S, "toy")

    # # define a slightly more complicated one
    # S = [[1, 4], [1, 4], [1, 4],
    #      [1, 2], [1, 2],
    #      [2, 3],
    #      [2, 3, 5],
    #      [3, 5]]
    # test_DD_creation(S, "simple")

    test_build_simple_MIP()
    ##test_BDD_to_MIP()
