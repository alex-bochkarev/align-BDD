import numpy as np
from copy import copy
from itertools import permutations
import pytest
from cUFL import generate_test_instance

class ColorSorter:
    def __init__(self, f_colors, node_order):
        self.f_colors = f_colors
        self.node_order = node_order
        self.no_colors = len(np.unique(f_colors))
        self.pos = {node_order[i]: i for i in range(len(node_order))}

    def rel_weight(self, elements_in_A, elements_in_B):
        """Calculates a 'relative weight' of the color.

        Args:
            elements_in_A (list): elements with color A,
            elements_in_B (list): elements with color B.
        Returns:
            RW (int): 'relative weight' of color A (relative to B)

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

    def sort_colors(self):
        """Sorts color blocks and nodes within each color.

        The aim is to have the minimum number of inversions between the resulting
        coloring diagram and the `node_order` parameter.

        Returns:
            colors (list): ordered color numbers.
            customers (dict): (ordered) list of customer per color.
        """
        colors = [set([idx+1 for idx, color in enumerate(self.f_colors)
                       if color == c]) for c in range(self.no_colors)]
        col_idx = list(range(len(colors)))
        changed = True

        while changed:
            changed = False
            for i in range(self.no_colors-1):
                for j in range(i, self.no_colors):
                    if self.rel_weight(colors[i], colors[j]) > 0:
                        c = copy(colors[i])
                        colors[i] = copy(colors[j])
                        colors[j] = c

                        c = col_idx[i]
                        col_idx[i] = col_idx[j]
                        col_idx[j] = c
                        changed = True
                        break

        print(f"final order={col_idx}")
        customers = {c: [] for c in range(self.no_colors)}

        for j in self.node_order:
            customers[self.f_colors[j-1]].append(j)

        return col_idx, customers

# ============================================================================
# Testing code
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


def bruteforce_correct_order(f_color, target_order):
    """Calculates the correct order with bruteforce."""
    min_inv = len(target_order)**2 + 1  # cannot be more than this
    best_order = None
    no_colors = len(np.unique(f_color))

    colors = [[] for c in range(no_colors)]
    for e in target_order:
        colors[f_color[e-1]].append(e)

    for order in permutations(list(range(len(colors)))):
        invs = no_invs(order, colors, target_order)
        if invs < min_inv:
            min_inv = invs
            best_order = order

    customers = {c: [] for c in best_order}

    for e in target_order:
        customers[f_color[e-1]].append(e)

    return sum([customers[c] for c in best_order], [])

def find_correct_order(f_color, target_order):
    """Finds the correct order with the code above."""

    cs = ColorSorter(f_color, target_order);
    colors, customers = cs.sort_colors()

    print(f"colors={colors}, cust={customers}")

    return sum([customers[c] for c in colors], [])

def make_instance(n):
    """Makes a quick test instance for ColorSorter"""
    S, f, f_colors, k_bar = generate_test_instance(n)
    target_order = np.random.permutation(list(range(1,len(f_colors)+1)))
    return (f_colors, target_order)

def get_score(A, target):
    """Calculates no. of inversions between `A` and `target`"""
    pos = {target[i]: i for i in list(range(len(target)))}
    invs = 0
    for i in range(len(A)-1):
        for j in range(i, len(A)):
            if pos[A[i]] > pos[A[j]]:
                invs += 1
    return invs

@pytest.mark.parametrize("test_inst", [make_instance(7)
                                       for _ in range(1000)])
def test_ColorSorter(test_inst):
    """Tests that ColorSorter finds an optimal order (compares to bruteforce)."""
    f_colors, target_order = test_inst
    print(f"The instance is: f_colors={f_colors}, target={target_order}")
    assert get_score(find_correct_order(f_colors, target_order), target_order) == get_score(
        bruteforce_correct_order(f_colors, target_order), target_order)
