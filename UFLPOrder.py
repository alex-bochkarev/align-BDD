"""Implements an algorithm to find  a var order for BDD representing an UFLP.

Assuming overlap costs and instance parameterization with S, f, c.
"""
from heapq import heappop, heappush
from itertools import count
import numpy as np

REMOVED = '<DEL>'


class N2RList:
    """Implements a heap of N2-neighborhoods keyed by size.

    Note:
        Based on implementation suggested in ``heapq``,
        see `heapq docs <https://docs.python.org/3/library/heapq.html>`_
        for more info.
    """
    def __init__(self, S):
        """Initializes the prio que structure."""
        self.S = S
        self.N2 = []
        self.counter = count()
        self.pq = []
        self.set_index = {}
        self.removed = []

        for j in range(len(S)):
            N2 = []
            for i in S[j]:
                N2 += [s for s in S[i-1] if s not in N2]

            self.N2.append(N2)

            curcount = next(self.counter)
            entry = [len(N2), curcount, (j+1)]
            self.set_index[j+1] = entry
            heappush(self.pq, entry)

    def _remove(self, N2):
        """Removes points ``N2`` from the heap (adjusting it as necessary)."""
        tbu = {}  # entries to be updated

        for k in N2:
            for i in self.N2[k-1]:
                # all such elements are affected (the key to be decremented.)
                if i in tbu:
                    tbu[i] += 1
                else:
                    tbu[i] = 1
            self.removed.append(k)

        for j in tbu:
            entry = self.set_index[j]
            newsize = entry[0] - tbu[j]
            if newsize > 0:
                newentry = [newsize, next(self.counter), j]
                entry[-1] = REMOVED
                heappush(self.pq, newentry)
                self.set_index[j] = newentry
            else:
                del self.set_index[j]

    def pop(self):
        """Pops a smallest-2-neighborhood point (adjusting the heap)."""
        while self.pq:
            _, _, j = heappop(self.pq)
            if j is not REMOVED:
                N2 = [s for s in self.N2[j-1] if s not in self.removed]
                self._remove(N2)
                return N2

        return []

    def __len__(self):
        return len(self.pq)


def UFLP_greedy_order(S):
    """Finds a good order for BDD representing UFLP.

    Args:
        S (list): adjacencies,
        f (list): overlap costs,
        c (list): location costs.

    Note:
        As usual, point  numbers are 1-based.

    Returns:
        A list of nodes to add to a BDD (according to a greedy algo).
    """
    N2 = N2RList(S)
    order = []
    while (len(N2) > 0) and len(order) < len(S):
        order += N2.pop()

    return order


# Testing code #######################################################
def test_greedy_order_simple():
    """Testing over a toy instance."""
    S = [[1, 2, 3],
         [1, 2, 3, 4],
         [1, 2, 3, 6],
         [2, 4, 5],
         [4, 5, 6],
         [3, 5, 6, 7],
         [6, 7]]

    t = N2RList(S)
    assert t.S == S
    assert [list(np.sort(row)) for row in t.N2] == [[1, 2, 3, 4, 6],
                                                    [1, 2, 3, 4, 5, 6],
                                                    [1, 2, 3, 4, 5, 6, 7],
                                                    [1, 2, 3, 4, 5, 6],
                                                    [2, 3, 4, 5, 6, 7],
                                                    [1, 2, 3, 4, 5, 6, 7],
                                                    [3, 5, 6, 7]]

    assert UFLP_greedy_order(S) == [3, 5, 6, 7, 2, 4, 1]


def test_greedy_order_toy2():
    """Another toy-instance test."""
    S = [[1, 2], [1,2,3,4], [2,3], [2,4,5,6],
         [4,5,10], [4,6,7], [6,7,8], [7,8,9,10],
         [8,9], [5,8,10], [11]]

    t = N2RList(S)
    assert t.S == S
    assert [list(np.sort(row)) for row in t.N2] == [[1,2,3,4],
                                                    [1,2,3,4,5,6],
                                                    [1,2,3,4],
                                                    [1,2,3,4,5,6,7,10],
                                                    [2,4,5,6,8,10],
                                                    [2,4,5,6,7,8],
                                                    [4,6,7,8,9,10],
                                                    [5,6,7,8,9,10],
                                                    [7,8,9,10],
                                                    [4,5,7,8,9,10],
                                                    [11]]
    order = UFLP_greedy_order(S)
    assert [list(np.sort(order[a[0]:a[1]])) for a in [(0,1), (1, 5),
                                                       (5, 7), (7,9),
                                                       (9,10), (10,11)]] == [
                                                           [11], [1,2,3,4],
                                                           [5,6], [7,10], [8],
                                                           [9]]
