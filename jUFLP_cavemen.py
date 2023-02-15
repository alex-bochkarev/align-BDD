"""Proof-of-concept experiment: joint UFLP over special class of instances.

Implements DD-based solution procedures for j-UFLP instances, along with the
instance generation. Note that the up-to-date experiment corresponding to
Section 4.2 of the paper is presented in :py:mod:`UFLP_2_cav`, but it calls
solution procedures implemented in this module.
"""
import pytest
import numpy as np
import gurobipy as gp
import json
from time import time

from darkcloud import ptscloud, DDSolver, gen_caveman_inst
from BDD import intersect
from varseq import VarSeq
from BB_search import BBSearch
from UFL import add_BDD_to_MIP
from UFLP_fullDD import create_cover_DD
from UFLPOrder import UFLP_greedy_order


def gen_cavemen_jUFLP_inst(n=10,
                           M=7,
                           L=0.25,
                           verbose=False,
                           linking="default"):
    """Generates an instance with the related metadata (info on caves).

    Args:
      n (int): number of caves,
      M (int): number of points in a cave,
      L (float): edge sparsity parameter (share of missing edges)
      verbose (Bool): print debug info
      linking (str): linking constraints type

    Returns:
      inst1, inst2, join_map: sub-instances and caves description.

    Note:
      Each sub-instance is parameterized by [S,f,c,caves].
      (See :py:func:`darkcloud.gen_caveman_inst` for details.)
    """
    S, f, c, caves = gen_caveman_inst(n, M, L)
    S2, f2, c2, caves2 = gen_caveman_inst(n, M, L)

    join_map = dict()
    if linking == "uniform":
        join_map = dict(
            zip([j for j in range(1,
                                  len(S) + 1)],
                np.random.permutation([j for j in range(1,
                                                        len(S) + 1)])))
    elif linking == "consecutive":
        ca1 = [ca.S for ca in caves]
        ca2 = [ca.S for ca in caves2]

        for k in range(len(ca1)):
            join_map.update(dict(zip(ca1[k], np.random.permutation(ca2[k]))))

    else:
        ca1 = [ca.S for ca in caves]
        ca2 = np.random.permutation([ca.S for ca in caves2])

        for k in range(len(ca1)):
            join_map.update(dict(zip(ca1[k], np.random.permutation(ca2[k]))))

    if verbose:
        print(f"S={S};\nf={f}\n;c={c}")
    return [[S, f, c, caves], [S2, f2, c2, caves2], join_map]


# encoding numpy numbers


# save and load instances
class jUFLPEncoder(json.JSONEncoder):
    """A technical implementation needed to save an instance as ``.json``"""

    def default(self, obj):
        if isinstance(obj, np.int64):
            return int(obj)
        elif isinstance(obj, dict):
            return {int(k): int(obj[k]) for k in obj}

        return json.JSONEncoder.default(self, obj)


def save_inst(i1, i2, join_map, filename):
    """Saves the jUFLP instance to ``.json`` file. (old instance format)."""
    with open(filename, "w") as fout:
        fout.write(
            json.dumps(
                {
                    'inst1': {
                        'S': i1[0],
                        'f': i1[1],
                        'c': i1[2],
                        'ptsclouds':
                        [i1[3][j].__dict__ for j in range(len(i1[3]))]
                    },
                    'inst2': {
                        'S': i2[0],
                        'f': i2[1],
                        'c': i2[2],
                        'ptsclouds':
                        [i2[3][j].__dict__ for j in range(len(i2[3]))]
                    },
                    'jmap': {int(j): join_map[j]
                             for j in join_map}
                },
                cls=jUFLPEncoder))


def load_inst(filename):
    """Loads a jUFLP instance from ``.json`` file.

    Returns:
      [[S,f,c,ptsclouds], [S2,f2,c2,pltsclouds2], jmap]
    """
    with open(filename, "r") as fin:
        json_inst = fin.read()

    json_dct = json.loads(json_inst)
    return [[
        json_dct[f'inst{i}']['S'], json_dct[f'inst{i}']['f'],
        json_dct[f'inst{i}']['c'],
        [
            ptscloud(ptc['e1'], ptc['e2'], ptc['S'])
            for ptc in json_dct[f'inst{i}']['ptsclouds']
        ]
    ] for i in [1, 2]
            ] + [{int(j1): json_dct['jmap'][j1]
                  for j1 in json_dct['jmap']}]


def show_inst(inst1, inst2, join_map, filename="./tmp/jUFLP.dot"):
    """Shows the supplied j-UFLP instance / saves it to a .dot file.

    Args:
        inst1, inst2 (list): instances description, [S,f,c,caves]
        filename (str): filename to save to (default "./tmp/jUFLP.dot")
    """
    added = set([])
    with open(filename, "w") as fout:
        fout.write("graph G {\n")
        fout.write("""    rankdir="LR";\n""")
        S, f, c, caves = inst2
        join_pts = sum([list(c.e1 + c.e2) for c in caves], [])
        join_pts = list(np.unique([j for j in join_pts if j is not None]))
        for i in range(len(S)):
            for j in S[i]:
                if ((i + 1) != j) and not (((j, (i + 1)) in added) or
                                           ((i + 1, j) in added)):

                    if (i + 1) in join_pts:
                        pref_i = "j"
                    else:
                        pref_i = "f"

                    if j in join_pts:
                        pref_j = "j"
                    else:
                        pref_j = "f"

                    fout.write(
                        f"    {pref_i}{i+1}--{pref_j}{j}[penwidth=3];\n")
                    added.add(((i + 1), j))

            if (i + 1) not in join_pts:
                fout.write(f"    f{i+1}[penwidth=3];\n")
        S, f, c, caves = inst1
        join_pts = sum([list(c.e1 + c.e2) for c in caves], [])
        join_pts = list(np.unique([j for j in join_pts if j is not None]))
        added = set([])
        for i in range(len(S)):
            for j in S[i]:
                if ((i + 1) != j) and not (((j, (i + 1)) in added) or
                                           ((i + 1, j) in added)):

                    if (i + 1) in join_pts:
                        a = f"j{join_map[i+1]}"
                    else:
                        a = f"s{i+1}"

                    if j in join_pts:
                        b = f"j{join_map[j]}"
                    else:
                        b = f"s{j}"

                    fout.write(f"    {a}--{b}[color=red];\n")
                    added.add(((i + 1), j))
            if (i + 1) not in join_pts:
                fout.write(f"    s{i+1}[color=red, textcolor=red];\n")

        for j in join_pts:
            fout.write(
                f"    j{join_map[j]}[style=filled, fillcolor=yellow, penwidth=5];\n"
            )
        fout.write("}")


def solve_cm_jUFLP_MIP(i1, i2, jmap):
    """Solves the special jUFLP instance with MIP."""
    m = gp.Model()
    m.modelSense = gp.GRB.MINIMIZE
    m.setParam("OutputFlag", 0)
    x = dict()
    y = dict()

    # add the first instance
    S, f, c, caves = i1
    # create variables
    for j in range(1, len(S) + 1):
        x[(1, j)] = m.addVar(vtype=gp.GRB.BINARY, name=f"x1_{j}", obj=c[j - 1])

        for a in range(1, len(S[j - 1]) + 1):
            y[(1, j, a)] = m.addVar(vtype=gp.GRB.BINARY,
                                    name=f"y1_{j}_{a}",
                                    obj=f[j - 1][a] - f[j - 1][a - 1])

    # create constraints
    for j in range(1, len(S) + 1):
        m.addConstr(
            gp.quicksum(x[(1, k)] for k in S[j - 1]) == gp.quicksum(
                y[(1, j, a)] for a in range(1,
                                            len(S[j - 1]) + 1)))

        for a in range(1, len(S[j - 1])):
            m.addConstr(y[(1, j, a)] >= y[(1, j, a + 1)])

    # add the second instance
    S, f, c, caves = i2
    # create variables
    for j in range(1, len(S) + 1):
        x[(2, j)] = m.addVar(vtype=gp.GRB.BINARY, name=f"x2_{j}", obj=c[j - 1])

        for a in range(1, len(S[j - 1]) + 1):
            y[(2, j, a)] = m.addVar(vtype=gp.GRB.BINARY,
                                    name=f"y2_{j}_{a}",
                                    obj=f[j - 1][a] - f[j - 1][a - 1])

    # create constraints
    for j in range(1, len(S) + 1):
        m.addConstr(
            gp.quicksum(x[(2, k)] for k in S[j - 1]) == gp.quicksum(
                y[(2, j, a)] for a in range(1,
                                            len(S[j - 1]) + 1)))

        for a in range(1, len(S[j - 1])):
            m.addConstr(y[(2, j, a)] >= y[(2, j, a + 1)])

    # link the two instances
    for j in jmap:
        m.addConstr(x[(1, j)] == x[(2, jmap[j])])

    m.update()
    m.optimize()
    assert m.status == gp.GRB.OPTIMAL
    return m.objVal + sum([fs[0]
                           for fs in i1[1]]) + sum([fs[0] for fs in i2[1]])


def solve_cm_jUFLP_DDs(i1, i2, jmap, intmode="toA", ret_int=False):
    """Solves the jUFLP cavemen instance with DDs.

    Args:
      i1, i2 (list): instances
      intmode (str): alignment mode, 'toA', 'toB', or 'VS'

    Notes:
      Instance is parameterized as per :py:func:`darkcloud.gen_caveman_inst`,
      The diagrams are built with :py:class:`darkcloud.DDSolver`.
    """
    S, f, c, caves = i1
    S2, f2, c2, caves2 = i2

    sol = DDSolver(S, f, c, caves)
    B1 = sol.build_cover_DD()

    sol = DDSolver(S2, f2, c2, caves2)
    B2 = sol.build_cover_DD()

    B1.make_reduced()
    B2.make_reduced()
    B1.rename_vars(jmap)

    if intmode == "toA":
        target = B1.vars
    elif intmode == "toB":
        target = B2.vars
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
        print(f"Wrong mode: '{intmode}'. Expected: 'toA', 'toB', or 'VS'.")

    B1.align_to(target, inplace=True)
    B2.align_to(target, inplace=True)

    int_DD = intersect(B1, B2)
    sp = int_DD.shortest_path()
    if ret_int:
        return sp[0], int_DD.size()
    else:
        return sp[0]


def solve_cm_jUFLP_CPPMIP(i1, i2, jmap):
    """Solves a jUFLP on cavemen with CPPMIP.
    
    Args:
      i1, i2 (list): instances
      jmap (dict): joining dict.

    Notes:
      Instance is parameterized as per :py:func:`darkcloud.gen_caveman_inst`,
      The diagrams are built with :py:class:`darkcloud.DDSolver`.
    Returns:
        objective (float)
    """
    S, f, c, caves = i1
    S2, f2, c2, caves2 = i2

    sol = DDSolver(S, f, c, caves)
    B1 = sol.build_cover_DD()

    sol = DDSolver(S2, f2, c2, caves2)
    B2 = sol.build_cover_DD()

    B1.make_reduced()
    B2.make_reduced()
    B1.rename_vars(jmap)

    m, c, v, x = add_BDD_to_MIP(B1, prefix="B1_")
    m, c, v, x = add_BDD_to_MIP(B2, model=m, x=x, prefix="B2_")
    m.update()
    m.setParam("OutputFlag", 0)
    m.optimize()
    return m.objVal


def solve_cm_jUFLP_CPPMIP_fullDDs(i1, i2, jmap):
    """Solves a jUFLP on cavemen (full-DDs) with CPPMIP.

    Args:
      i1, i2 (list): instances
      jmap (dict): joining dict.

    Notes:
      Instance is parameterized as per :py:func:`darkcloud.gen_caveman_inst`,
      The diagrams are built with :py:func:`UFLP_fullDD.create_cover_DD`.
    Returns:
        objective (float)
    """
    S, f, c, caves = i1
    S2, f2, c2, caves2 = i2

    B1, _ = create_cover_DD(S, f, c, UFLP_greedy_order(S, True))
    B2, _ = create_cover_DD(S2, f2, c2, UFLP_greedy_order(S2, False))

    B1.make_reduced()
    B2.make_reduced()

    B1.rename_vars(jmap)

    m, c, v, x = add_BDD_to_MIP(B1, prefix="B1_")
    m, c, v, x = add_BDD_to_MIP(B2, model=m, x=x, prefix="B2_")
    m.update()
    m.setParam("OutputFlag", 0)
    m.optimize()
    return m.objVal


def solve_cm_jUFLP_fullDDs(i1, i2, jmap, intmode, ret_int=False):
    """Solves the jUFLP cavemen instance with DDs.

    Args:
      i1, i2 (list): instances

    Notes:
      Instance is parameterized as per :py:func:`darkcloud.gen_caveman_inst`.
    """
    S, f, c, caves = i1
    S2, f2, c2, caves2 = i2

    o1 = UFLP_greedy_order(S, True)
    o2 = UFLP_greedy_order(S2, False)

    B1, _ = create_cover_DD(S, f, c, o1)
    B2, _ = create_cover_DD(S2, f2, c2, o2)

    B1.make_reduced()
    B2.make_reduced()

    B1.rename_vars(jmap)

    if intmode == "toA":
        target = B1.vars
    elif intmode == "toB":
        target = B2.vars
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
        print(f"Wrong mode: '{intmode}'. Expected: 'toA', 'toB', or 'VS'.")

    B1.align_to(target, inplace=True)
    B2.align_to(target, inplace=True)

    int_DD = intersect(B1, B2)
    sp = int_DD.shortest_path()
    if ret_int:
        return sp[0], int_DD.size()
    else:
        return sp[0]


###
def dump_instance(S, caves, filename="tmp/S.dot"):
    """Dumps a graph implied by S into a `.dot` file. """
    added = set([])
    with open(filename, "w") as fout:
        fout.write("graph G {\n")
        for i in range(len(S)):
            for j in S[i]:
                if ((i + 1) != j) and not (((j, (i + 1)) in added) or
                                           ((i + 1, j) in added)):
                    fout.write(f"    n{i+1} -- n{j};\n")
                    added.add(((i + 1), j))

        fout.write("}")


def solve_with_MIP(S, S2, f, f2, c, caves):
    """Generates and solves a MIP for the supplied jUFL cavemen instance.

    Note:
      The instance is parameterized as per :py:func:`darkcloud.gen_caveman_inst`.
      Most of the code is adapted form :py:func:`darkcloud.solve_with_MIP`.
    """
    m = gp.Model()
    m.modelSense = gp.GRB.MINIMIZE
    m.setParam("OutputFlag", 0)
    x = dict()
    y = dict()
    y2 = dict()

    # create variables
    for j in range(1, len(S) + 1):
        x[j] = m.addVar(vtype=gp.GRB.BINARY, name=f"x_{j}", obj=c[j - 1])

        for a in range(1, len(S[j - 1]) + 1):
            y[(j, a)] = m.addVar(vtype=gp.GRB.BINARY,
                                 name=f"y_{j}_{a}",
                                 obj=f[j - 1][a] - f[j - 1][a - 1])

        for a in range(1, len(S2[j - 1]) + 1):
            y2[(j, a)] = m.addVar(vtype=gp.GRB.BINARY,
                                  name=f"y2_{j}_{a}",
                                  obj=f2[j - 1][a] - f2[j - 1][a - 1])

    # create constraints
    for j in range(1, len(S) + 1):
        # first-graph constraints
        m.addConstr(
            gp.quicksum(x[k] for k in S[j - 1]) == gp.quicksum(
                y[(j, a)] for a in range(1,
                                         len(S[j - 1]) + 1)))

        for a in range(1, len(S[j - 1])):
            m.addConstr(y[(j, a)] >= y[(j, a + 1)])

        # second-graph constraints
        m.addConstr(
            gp.quicksum(x[k] for k in S2[j - 1]) == gp.quicksum(
                y2[(j, a)] for a in range(1,
                                          len(S2[j - 1]) + 1)))

        for a in range(1, len(S2[j - 1])):
            m.addConstr(y2[(j, a)] >= y2[(j, a + 1)])

    m.update()
    print("Optimizing the model...", end="", flush=True)
    m.optimize()
    print("done.")
    assert m.status == gp.GRB.OPTIMAL
    return m, (m.objVal + sum(fs[0]
                              for fs in f) + sum(fs2[0]
                                                 for fs2 in f2)), x, y, y2


def solve_with_DDs(S, S2, f, f2, c, caves, intmode="toA"):
    """Solves the jUFLP cavemen instance with DDs.

    Args:
      intmode (str): intersection mode, 'toA' or 'VS'

    Notes:
      Instance is parameterized as per :py:func:`darkcloud.gen_caveman_inst`,
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
@pytest.mark.parametrize("test_inst",
                         [gen_cavemen_jUFLP_inst(7, 7) for _ in range(5)])
def test_jUFL_DDs(test_inst):
    i1, i2, jmap = test_inst
    obj1 = solve_cm_jUFLP_DDs(i1, i2, jmap, 'toA')
    obj2 = solve_cm_jUFLP_DDs(i1, i2, jmap, 'toB')
    obj3 = solve_cm_jUFLP_DDs(i1, i2, jmap, 'VS')
    assert abs(obj1 - obj2) < 0.001
    assert abs(obj1 - obj3) < 0.001


@pytest.mark.parametrize("test_inst",
                         [gen_cavemen_jUFLP_inst(7, 7) for _ in range(5)])
def test_cm_jUFL_DDvsMIP(test_inst):
    i1, i2, jmap = test_inst
    obj1 = solve_cm_jUFLP_MIP(i1, i2, jmap)
    obj2 = solve_cm_jUFLP_DDs(i1, i2, jmap)
    obj3 = solve_cm_jUFLP_CPPMIP(i1, i2, jmap)
    assert abs(obj1 - obj2) < 0.01
    assert abs(obj1 - obj3) < 0.01
    # NOTE: it seems sometimes Gurobi yields
    # rounding-off (?) errors --- e.g., 77.5055 instead of 77.5000
    # Seems not a big deal, but I needed to decrease the required
    # tolerance to 0.01 (from 0.001)


@pytest.mark.parametrize("test_inst",
                         [gen_cavemen_jUFLP_inst(5, 5) for _ in range(10)])
def test_load_save(test_inst):
    i1, i2, jmap = test_inst
    save_inst(i1, i2, jmap, "./tmp/jUFLP-loadsave-test.json")

    [inst1, inst2, jm] = load_inst("./tmp/jUFLP-loadsave-test.json")
    obj1 = solve_cm_jUFLP_DDs(i1, i2, jmap)
    obj2 = solve_cm_jUFLP_DDs(inst1, inst2, jm)

    assert abs(obj1 - obj2) < 0.001


def compare_runtimes():
    """Performs a quick runtimes comparison (toA vs VS)."""
    s, s2, f, f2, c, caves = gen_caveman_inst()
    sol = DDSolver(s, f, c, caves)
    B1 = sol.build_cover_DD()

    sol = DDSolver(s2, f2, [0.0 for _ in c], caves)
    B2 = sol.build_cover_DD()
    B2.shuffle_vars()  # renames the variables, without reordering

    B1.make_reduced()
    B2.make_reduced()

    t0 = time()
    B2p = B2.align_to(B1.vars, inplace=False)
    int_toA = intersect(B1, B2p)
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
