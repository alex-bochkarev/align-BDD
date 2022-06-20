"""Solves UFLP over cavemen instances using a full-DD approach."""
from BDD import BDD, NROOT, NTRUE
from darkcloud import solve_with_MIP, gen_caveman_inst
from copy import copy
import pytest
import numpy as np
from UFLPOrder import UFLP_greedy_order

def make_label(state):
    """Helper: formats a node label according to the ``state``."""
    return "\n"+"\n".join([f"{j+1}:{state[j]}" for j in range(len(state))
                      if (state[j])])


def calc_cost(pts_list, state, prev_state, last_decision, S, f, c):
    """Helper: calculates the added cost during the BDD creation.

    Args:
        S, f, c: UFLP instance parameters,
        pts_list (list): points to generate cost for,
        state: the new state
        last_decision(tuple): a var (int) - value(int) pair.

    Returns:
        Added cost according to the ``state``.
    """
    cost = 0.0
    for j in pts_list:
        a = sum([1 * ((state[i-1] == 1)
                      or (prev_state[i-1] == 1)
                      or ((i == last_decision[0]) and (last_decision[1] == 1)))
                 for i in S[j-1]])

        cost += f[j - 1][a]

        if (state[j - 1] == 1) or (prev_state[j - 1]
                                   == 1) or ((j == last_decision[0])
                                             and (last_decision[1] == 1)):
            cost += c[j - 1]

    return cost


def create_cover_DD(S, f, c, node_order):
    """Generates a BDD for an UFLP with cavemen structure.

    Args:
        S, f, c: instance params (see :py:func:`darkcloud.gen_caveman_inst`),
        node_order (list[int]): order of nodes to add.

    Returns:
        B (class BDD): resulting diagram

    Note:
        Implements 'forgetting' nodes along the way.
    """
    N = len(S)
    assert N == len(node_order), f"node_order must contain exactly {N} nodes."

    B = BDD(N=N, vars=[n for n in node_order],
            weighted=True)

    # how many node costs do I need this point for?
    node_value = [len(S[j-1]) for j in range(1, len(S)+1)]

    # how many more points do I need to add to calculate costs?
    # ("cost uncertainty")
    cost_unc = [len(S[j-1]) for j in range(1, len(S)+1)]

    root_state = tuple(0 for _ in range(len(S)))
    # STATE CODES:
    # 1 = true (set), -1 = false (unset), 0 = doesn't matter

    node_labels = dict({NROOT: make_label(root_state)})

    next_layer = {root_state: B.addnode(None)}

    k = 0  # created layers counter
    nodes = copy(node_order)

    while k < N:
        i = nodes.pop(0)
        add_costs = []  # list of points to calc costs for during the step

        for j in S[i-1]:
            cost_unc[j-1] -= 1
            if cost_unc[j-1] == 0:
                add_costs.append(j)
                for s in S[j-1]:
                    node_value[s-1] -= 1

        # adding i-th point (of the original graph)
        current_layer = copy(next_layer)
        if k == N-1:
            # last layer
            next_layer = {tuple([0 for _ in root_state]): B.nodes[NTRUE]}
        else:
            next_layer = dict()

        for state in current_layer:
            # a no-arc
            node = current_layer[state]
            next_state = tuple([
                (state[k] - 1 * ((k + 1) == i)) * (node_value[k] > 0)
                for k in range(len(state))
            ])

            arc_cost = calc_cost(add_costs, next_state, state, (i, -1),
                                 S, f, c)

            if next_state in next_layer:
                B.llink(node, next_layer[next_state], "lo",
                        edge_weight=arc_cost)
            else:
                newnode = B.addnode(node, "lo", edge_weight=arc_cost)
                next_layer.update({next_state: newnode})
                node_labels.update({newnode.id: make_label(next_state)})

            # a yes-arc
            next_state = tuple([(state[q] + ((q+1) == i)) * (node_value[q] > 0)
                                for q in range(len(state))])

            arc_cost = calc_cost(add_costs, next_state, state, (i, 1),
                                 S, f, c)

            if next_state in next_layer:
                B.llink(node, next_layer[next_state], "hi",
                        edge_weight=arc_cost)
            else:
                newnode = B.addnode(node, "hi", edge_weight=arc_cost)
                next_layer.update({next_state: newnode})
                node_labels.update({newnode.id: make_label(next_state)})

        k += 1

    return B, node_labels


# Testing code ######################################################
@pytest.mark.parametrize("inst", [
    gen_caveman_inst(n=np.random.randint(5, 8),
                     M=np.random.randint(5, 8),
                     L=max(0.1 + np.random.rand(), 0.9)) for _ in range(100)])
def test_full_DD(inst):
    """Tests the full-DD-based approach against a MiP."""
    S, f, c, _ = inst
    _, objMIP, _, _ = solve_with_MIP(S, f, c)
    order = UFLP_greedy_order(S)
    print(f"order = {order}")
    B, nl = create_cover_DD(S, f, c, order)
    objDD = B.shortest_path()[0]
    assert abs(objMIP - objDD) < 0.01


def test_cost_calc():
    """Tests :py:func:`UFLP_fullDD.calc_cost` function."""
    S = [[1, 2,3,4,5], [2, 1,3], [3, 1,2,6], [4, 1], [5, 1,6], [6, 5,3]]
    f = []
    for k in range(len(S)):
        f.append([10, 0] + [1+j+0.1*(k+1) for j in range(len(S[k]))])

    c = [1, 2, 3, 4, 5, 6]

    # testing for a specific instance
    state = (0,0,0,0,0,0)
    prev_state = (0,0,0,0,0,0)
    ld = (1, -1)
    forget = [1]
    assert calc_cost(forget, state, prev_state, ld, S, f, c) == 10.0

    forget = [1, 3, 5]
    assert calc_cost(forget, state, prev_state, ld, S, f, c) == 30.0

    forget = [4, 6]
    assert calc_cost(forget, state, prev_state, ld, S, f, c) == 20.0

    forget = [1]
    ld = (1, 1)
    assert calc_cost(forget, state, prev_state, ld, S, f, c) == 1.0

    forget = [1]
    ld = (3, -1)
    prev_state = (1,0,0,0,0,0)
    assert calc_cost(forget, state, prev_state, ld, S, f, c) == 1.0

    state = (1, -1, -1, 1, 0, 0)
    assert calc_cost(forget, state, prev_state, ld, S, f, c) == 2.1

    ld = (1, 1)
    assert calc_cost(forget, state, prev_state, ld, S, f, c) == 2.1


def test_fullDD_simple():
    S = [[1, 2,3,4,5], [2, 1,3], [3, 1,2,6], [4, 1], [5, 1,6], [6, 5,3]]
    f = []
    for k in range(len(S)):
        f.append([10, 0] + [1+j+0.1*(k+1) for j in range(len(S[k]))])

    c = [1, 2, 3, 4, 5, 6]

    return create_cover_DD(S, f, c, UFLP_greedy_order(S))
