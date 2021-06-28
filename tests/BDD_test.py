"""Tests some key procedures for :py:mod:`BDD` .

"""
import BDD
import pytest
import numpy as np
from random import randint
import copy

@pytest.mark.parametrize("D", [BDD.BDD.random(N = 1 + np.random.randint(1, 10),
                                              p = np.random.uniform(),
                                              weighted=np.random.choice([True, False]))
                               for _ in range(500)])
def test_load_save(D):
    """Tests that load/save functionality works as intended."""
    D.save("tests/BDD_load_save.bdd")
    print(f"D is: {D}")
    D2 = BDD.BDD()
    D2.load("tests/BDD_load_save.bdd")
    assert D.profile() == D2.profile()


@pytest.mark.parametrize("test_inst", [(1+np.random.randint(5, 15),
                                        np.random.uniform())
                                       for _ in range(500)])
def test_shortest_path(test_inst):
    """Tests the shortest path procedure (with Gurobi model)"""
    N, p = test_inst
    A = BDD.BDD.random(N, p, weighted=True)
    from UFL import add_BDD_to_MIP

    m, c, v, x = add_BDD_to_MIP(A)
    m.update()
    m.optimize()

    nl = A.shortest_path()
    if m.status == 2:
        # optimal
        assert nl[BDD.NROOT] == m.objVal, "SP={nl[NROOT], while Gurobi gave {m.objVal}}"
    else:
        assert nl[BDD.NROOT] == "∞", "with Gurobi status {m.status}, root node label is {nl[NROOT]}"

@pytest.mark.parametrize("test_inst", [(1+np.random.randint(3, 7),
                                        np.random.randint(1,20),
                                        np.random.uniform(),
                                        'swaps')
                                       for _ in range(50)]+ \
                         [(1+np.random.randint(3, 7),
                           np.random.randint(1,20),
                           np.random.uniform(),
                           'sifts')
                          for _ in range(50)]
                         )
def test_random_swapping(test_inst):
    """Tests if :py:func:`BDD.BDD.swap` work correctly (random problems).

    Generates a random BDD, makes some random swaps (or sifts),
    and then checks that the function encoded by the BDD remained the same
    (that is, every set of var values result in the same terminal node for
    the original BDD and for the randomly changed BDD).

    Args:
        N: no. of variables in the diagram
        m: no. of consecutive test swaps per BDD
        p: BDD expansion prob parameter
        mode: ``swaps`` or ``sifts``, determines what kind of events
                are to be generated

    Notes:
        - Arguments are packed into a single tuple (for compatibility)
            with ``pytest`` module.
        - Involves brute-force enumaration of all the 2^n possible decisions
            (concerning all vars), so choose ``n`` wisely (<=8 seems OK)

    Example:
        test_random_swapping((8,20,0.8,sys.argv[1]))
    """
    N, m, p, mode = test_inst

    bdd = BDD.BDD.random(N, p)
    print(f"BDD proflie is: {bdd}")  # so we could reconstruct the instance
    bdd_s = copy.deepcopy(bdd)

    # make some random layer swaps
    assert mode in ['swaps', 'sifts'], f"Only 'swap' and 'sift' modes are allowed (got {mode})"

    if mode == 'swaps':
        for _ in range(m):
            bdd_s.swap_up(randint(1, len(bdd.vars) - 1))
    else:
        for _ in range(m):
            bdd_s.sift(bdd_s.vars[randint(0,
                                            len(bdd.vars) - 1)],
                        randint(0,
                                len(bdd.vars) - 1))

    # now check that everything remained the same
    no_of_bits = len(bdd.vars)

    for bit_path in range(0, 2**no_of_bits):

        # unpack the var values in the original order (1=hi, 0=lo)
        var_vals = [(bit_path >> bit) & 1
                    for bit in range(no_of_bits - 1, -1, -1)]

        # get the correct answer (T or F)
        cur_node_src = list(bdd.layers[0])[0]
        cur_node_swd = list(bdd_s.layers[0])[0]

        i = 0
        while i < len(bdd.vars):
            if var_vals[i] == 0:
                cur_node_src = cur_node_src.lo
            else:
                cur_node_src = cur_node_src.hi

            if var_vals[bdd.p(bdd_s.vars[i])] == 0:
                cur_node_swd = cur_node_swd.lo
            else:
                cur_node_swd = cur_node_swd.hi

            i += 1

        corr_answer = cur_node_src.id
        assert corr_answer == BDD.NTRUE or corr_answer == BDD.NFALSE

        assert cur_node_swd.id == corr_answer, \
            "ERROR: {}, not {}".format(cur_node_swd.id,
                                         cur_node_src) + \
                                         "path {}".format(var_vals)


def test_swaps_w():
    """Tests :py:func:`BDD.BDD.swap` for a weighted BDD (a simple case).

    Uses a simple problem instance (``tests/simple_DD.wbdd``)
    to run a test: all paths costs and terminals must coincide
    after a series of swaps.
    """
    A = BDD.BDD()
    A.load("./tests/simple_DD.wbdd", weighted=True)

    B = BDD.BDD()
    B.load("./tests/simple_DD.wbdd", weighted=True)
    eq, msg = A.is_equivalent(B)
    assert eq, msg

    B.swap_up(1)
    eq, msg = A.is_equivalent(B)
    assert eq, msg

    B.swap_up(1)
    eq, msg = A.is_equivalent(B)
    assert eq, msg

    B.swap_up(1)
    B.swap_up(2)
    eq, msg = A.is_equivalent(B)
    assert eq, msg


@pytest.mark.parametrize("i", [i for i in range(1000)])
def test_swaps_weighted(i):
    """Tests swap operation for weighted diagrams (random instances)."""
    A = BDD.BDD.random(N=5, weighted=True)
    B = copy.deepcopy(A)

    B.swap_up(1)
    eq, msg = A.is_equivalent(B)
    assert eq, msg

    B.swap_up(len(A)-1)
    eq, msg = A.is_equivalent(B)
    assert eq, msg

    for k in range(5):
        C = copy.deepcopy(B)
        what = np.random.randint(2, len(A))
        B.swap_up(what)
        eq, msg = A.is_equivalent(B)
        assert eq, f"{i}: for swap at {pos}: {msg}"


@pytest.mark.parametrize("i", [i for i in range(1000)])
def test_swaps_uweighted(i):
    A = BDD.BDD.random(N=7, weighted=False)
    B = copy.deepcopy(A)

    B.swap_up(1)
    eq, msg = A.is_equivalent(B)
    assert eq, msg

    B.swap_up(len(A)-1)
    eq, msg = A.is_equivalent(B)
    assert eq, msg

    for _ in range(5):
        what = np.random.randint(2, len(A))
        print(f"{i}: swapping up {what}")
        B.swap_up(what)
        eq, msg = A.is_equivalent(B)
        assert eq, msg


@pytest.mark.parametrize("i", [i for i in range(500)])
def test_intersect(i):
    """Tests intersection function, :py:func:`BDD.intersect`."""
    A = BDD.BDD.random(N=7, weighted=True)
    B = BDD.BDD.random(N=7, weighted=True)

    C = BDD.intersect(A, B)

    tt_A = A.truth_table()
    tt_B = B.truth_table()
    tt_C = C.truth_table()

    for idx in tt_C.index:
        assert (tt_C["Terminal"][idx] == (tt_A["Terminal"][idx] and tt_B["Terminal"][idx])), f"{i} failed"
        assert abs(tt_C["Cost"][idx] - (tt_A["Cost"][idx] + tt_B["Cost"][idx]))<1e-5, f"{i} failed"

#################################################################
## Ad-hoc testing / visual inspection

def save_load(bdd):
    """Quickly tests the load-save functionality (for visual inspection).

    must result in two equivalent graphs (BDDs) on the screen
    """
    bdd.save("./test_save.bdd")
    bdd2 = BDD()
    bdd2.load("./test_save.bdd")
    bdd.dump_gv().view("./before_save.dot",cleanup=True)
    bdd2.dump_gv().view("./after_save.dot",cleanup=True)
    eq, msg = bdd.is_equivalent(bdd2)
    assert eq, msg


def gen_4BDD():
    """creates a simple 4-var BDD"""

    bdd = BDD(4)
    root = bdd.addnode(None)
    nh = bdd.addnode(root,"hi")
    nl = bdd.addnode(root,"lo")

    nh2 = bdd.addnode(nh,"hi")
    bdd.nodes[nh.id].lo = nh2

    nl2 = bdd.addnode(nl, "hi")
    nl.link(nl2,"lo")

    nl = bdd.addnode(nl2,"hi")
    nl2.link(nl,"lo")

    nh = bdd.addnode(nh2,"lo")
    nh2.link(nh,"hi")

    nl.link(bdd.T,"lo")
    nl.link(bdd.F,"hi")

    nh.link(bdd.T,"lo")
    nh.link(bdd.F, "hi")
    return bdd


def test_save_load_noargs():
    """Tests load/save for a specific situation."""
    bdd = gen_4BDD()
    bdd.sift(3,0)
    save_load(bdd)


# BDD creation and rendering (with graphviz)
def show_create_render():
    """Tests the BDD creation and rendering code."""
    bdd = BDD.BDD(4)
    root = bdd.addnode(None)
    nh = bdd.addnode(root, "hi")
    nl = bdd.addnode(root, "lo")

    nh2 = bdd.addnode(nh, "hi")
    bdd.nodes[nh.id].lo = nh2

    nl2 = bdd.addnode(nl, "hi")
    nl.link(nl2,"lo")

    nl = bdd.addnode(nl2, "hi")
    nl2.link(nl, "lo")

    nh = bdd.addnode(nh2,"lo")
    nh2.link(nh,"hi")

    nl.link(bdd.T,"lo")
    nl.link(bdd.F,"hi")

    nh.link(bdd.T,"lo")
    nh.link(bdd.F, "hi")

    g = bdd.dump_gv()

    g.render("./initial.dot",view=True)

def show_swap_sift(bdd):
    """quick test of swap and sift operations"""
    ## swap operation
    bdd.swap_up(2)

    g = bdd.dump_gv()
    g.render("./swapped.dot", view=True)

    bdd.swap_up(1)

    g = bdd.dump_gv()
    g.render("./swapped_top.dot", view=True)

    ## sift operatio
    bdd.sift(3,3)
    g = bdd.dump_gv()
    g.render("./sifted.dot", view=True)

    ### more complex sifting example

    bdd = BDD(4)
    root = bdd.addnode(None)
    n1 = bdd.addnode(root, "hi")
    n2 = bdd.addnode(root, "lo")

    n3 = bdd.addnode(n1,"hi")
    n4 = bdd.addnode(n1,"lo")
    n5 = bdd.addnode(n2,"hi")
    n2.link(n4,"lo")

    n6 = bdd.addnode(n3,"hi")
    n3.link(n6,"lo")

    n7 = bdd.addnode(n4,"hi")
    n4.link(n7,"lo")

    n8 = bdd.addnode(n5, "lo")
    n5.link(n7,"hi")

    n6.link(bdd.F,"hi")
    n6.link(bdd.F,"lo")

    n7.link(bdd.F,"hi")
    n7.link(bdd.F,"lo")

    n8.link(bdd.F,"lo")
    n8.link(bdd.T,"hi")

    # g = bdd.dump_gv()
    # g.view()

    bdd.sift(4,1)

    g = bdd.dump_gv()
    g.render("./after_sifting.dot",view=True)


def show_intersections():
    """Tests intersection function."""
    A = BDD()
    A.load("./tests/int_A.wbdd", weighted=True)
    B = BDD()
    B.load("./tests/int_B.wbdd", weighted=True)

    C = intersect(A, B)
    A.show(dir="testing", filename="A.dot")
    B.show(dir="testing", filename="B.dot")
    C.show(dir="testing", filename="C.dot")


def show_rnd():
    bdd = BDD.BDD.random(N=5,p=0.8)
    bdd.dump_gv().view("./random.dot",directory="./testing",cleanup = True)

def show_align():
    bdd_A = BDD.BDD.random(N=4,p=0.8)
    bdd_B = BDD.BDD.random(N=4,p=0.8)

    alts, Aa, Ba = bdd_A.OA_bruteforce(bdd_B)
    print("Opt order: {}, opt size:{}".format(Aa[0].vars, Aa[0].size()+Ba[0].size()))

    bdd_A.dump_gv().render("./testing/A.dot",view=True)
    bdd_B.dump_gv().render("./testing/B.dot",view=True)
    Aa[0].dump_gv().render("./testing/Aa.dot",view=True)
    Ba[0].dump_gv().render("./testing/Ba.dot", view=True)


def show_bruteforcing():
    bdd = BDD.BDD.random(N=4, p=0.8) # it works for gen_4BDD(), though

    perms = iters.permutations(bdd.vars)

    for p in perms:
        print("Aligning to: {}".format(p))
        bdd = bdd.align_to(list(p))
        bdd.show()
        input("Press Enter to continue...")

def show_rnd_naming():
    mybdd = BDD.BDD.random(N=4,p=0.8)
    mybdd.show()

    mybdd_al = copy.deepcopy(mybdd)
    mybdd_al.show()

    mybdd_al.swap_up(2)
    # mybdd_al = mybdd.align_to([2,1,3,4])
    mybdd_al.show()
    mybdd_al.swap_up(mybdd_al.p(3))
    mybdd_al.show()

def show_swapping_2():
    mybdd = BDD.BDD.random(N=4,p=0.8)
    mybdd.show()

    bdds = copy.deepcopy(mybdd)
    bdds.swap_up(3)
    bdds.dump_gv().render("./testing/after_swap.dot",view=True, cleanup=True)


    mybdd = BDD.BDD.random(N=4,p=0.8)
    mybdd.show()

    bdds = copy.deepcopy(mybdd)
    bdds.swap_up(3)
    bdds.dump_gv().render("./testing/after_swap.dot", view=True, cleanup=True)
