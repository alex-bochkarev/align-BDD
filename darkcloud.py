"""
Implements the 'cloud' algorithm that builds a BDD for a UFL over a relaxed
cavemen graph.

---
(c) A. Bochkarev, Clemson University, 2022
a@bochkarev.io
"""
import numpy as np
from dataclasses import dataclass
from experiments.softcover import generate_overlaps
from BDD import BDD, NTRUE


@dataclass
class ptscloud:
    """Keeps current 'cloud' of points.

    Attributes:
        e1, e2 (tuple): the two external entry points of the cloud (edges),
        S (list[int]): a list of points within the cloud.
    """
    e1: tuple[int]
    e2: tuple[int]
    S: list[int]


# def gen_caveman_inst(n=10, M=5, pcave=0.8, verbose=True):
#     """Generates an instance with the related metadata (info on caves).

#     Args:
#       n (int): number of caves,
#       M (int): number of points in a cave0,
#       pcave (float): probability of connection within the cave0,
#       verbose (Bool): print debug info

#     Returns:
#       S, f, c, caves: instance (S,f,c) and cave0 description.
#     """
#     # Create a connected 'string' first
#     S = [[j+1] for j in range(n)]
#     for j in range(1, len(S)-1):
#         S[j].append(j)
#         S[j].append(j+2)

#     S[-1].append(len(S)-1)
#     S[0].append(2)

#     ncaves = n // M
#     for k in range(1, ncaves+1):
#         for i in range((k-1)*M, min(k*M, n)):
#             for j in range(i+2, min(k*M, n)):
#                 if np.random.uniform() <= pcave:
#                     S[i].append(j+1)
#                     S[j].append(i+1)

#     f = generate_overlaps(S)
#     Cmax = 5.0
#     c = [Cmax*np.random.uniform() for _ in range(len(S))]
#     if verbose:
#         print(f"S={S};\nf={f}\n;c={c}")
#     return S, f, c

def gen_simple_cavemen_inst1():
    """Generates a simple instance with metadata.

    Returns:
      S, f, c, caves
    """
    S = [[1, 2, 4, 5],
         [2, 1, 3, 8, 10],
         [3, 2, 12, 13, 15],
         [4, 1, 5, 6],
         [5, 1, 4, 6, 7],
         [6, 4, 5, 7],
         [7, 6, 5],
         [8, 2, 9, 10, 11],
         [9, 8, 11],
         [10, 2, 8],
         [11, 8, 9],
         [12, 3, 13, 14],
         [13, 3, 12, 14],
         [14, 12, 13, 15, 16],
         [15, 3, 14],
         [16, 14]]

    f = generate_overlaps(S)

    c = [5.0*np.random.uniform() for _ in range(len(S))]
    caves = [ptscloud(2, 0, [1,4,5,6,7]),
             ptscloud(1, 3, [2,8,9,10,11]),
             ptscloud(2, 0, [3,12,13,14,15,16])]

    return S, f, c, caves


def gen_simple_cavemen_inst1():
    """Generates (another) simple instance with metadata.

    Returns:
      S, f, c, caves
    """
    S = [[1,2,3],
         [2,1,3,4],
         [3,1,2,4],
         [4,2,3,5],
         [5,4,6,7],
         [6,5,8],
         [7,5,8],
         [8, 6,7,9],
         [9, 8,11],
         [10, 11],
         [11, 9,10,12],
         [12, 11,13],
         [13, 12,14,15],
         [14, 13,15],
         [15, 13,14,18],
         [16, 18],
         [17, 18,20],
         [18, 15,16,17],
         [19, 20,22],
         [20, 17, 19,21],
         [21, 20,22],
         [22, 19,21,23],
         [23, 22,25,26],
         [24, 26],
         [25, 23, 26],
         [26, 23, 24,25]]

    caves = [ptscloud((None, None), (4,5), [1,2,3,4]),
             ptscloud((4,5), (8,9), [5,6,7,8]),
             ptscloud((8,9), (12,13), [9,10,11,12]),
             ptscloud((12,13), (15,18), [13,14,15]),
             ptscloud((15,18), (17,20), [16,17,18]),
             ptscloud((17,20), (22,23), [19,20,21,22]),
             ptscloud((22,23), (None, None), [23,24,25,26])]

    f = generate_overlaps(S)

    c = [5.0*np.random.uniform() for _ in range(len(S))]

    return S, f, c, caves


class DDSolver:

    def __init__(self, S, f, c, caves):
        self.S = S
        self.f = f
        self.c = c
        self.caves = caves
        self.B = None

    def _add_interim_point(self, x, current_state, new_state,  current_layer,
                           add_costs=None):
        """Adds variable after the current layer."""
        assert self.B is not None
        assert new_state[-1] == x, f"Last var, {new_state[-1]}, is not {x}"

        next_layer = dict()
        state_pos = {var: idx for idx, var in enumerate(current_state)}
        state_pos.update({x: len(state_pos)})

        for state in current_layer:
            for (xval, arc) in [(True, "hi"), (False, "lo")]:
                astate = state + (xval,)
                newstate = tuple(astate[state_pos[j]] for j in new_state)
                if add_costs is None:
                    cost = 0.0
                else:
                    cost = self._calc_cave(add_costs,
                                           fixed_nodes={
                                               new_state[i]: newstate[i]
                                               for i in range(len(newstate))
                                           })
                if newstate in next_layer:
                    self.B.llink(current_layer[state], next_layer[newstate],
                                 arc, cost)
                else:
                    newnode = self.B.addnode(current_layer[state], arc,
                                             edge_weight=cost)
                    next_layer.update({newstate: newnode})

        print(f"Added: point {x}.")
        return next_layer

    def _calc_cave(self, cave, fixed_nodes):
        """Calculates part of the objective due to the given cave."""
        print(f"Calc: {cave} with the decisions fixed:\n{fixed_nodes}")
        return -1.0

    def build_cover_DD(self):
        """Builds a cover/overlap DD for the instance."""
        assert len(caves) > 1

        vars_to_add = sum([[c.e1[0], c.e1[1], c.e2[0], c.e2[1]]
                          for c in self.caves], [])
        vars_to_add = np.unique([v for v in vars_to_add if v is not None])
        print(f"vars = {vars_to_add}")
        self.B = BDD(N=len(vars_to_add), weighted=True)
        varnames = []

        current_layer = {tuple(): self.B.addnode(parent_node=None)}

        self.curr_state = []

        for cavenum, cave in enumerate(self.caves[:-1]):
            print(f"=== Adding cave {cave}===")
            new_points = []
            drop_points = []
            for x in (cave.e1 + cave.e2):
                if x is None:
                    continue

                if x not in self.curr_state:
                    if x not in new_points:  # corner case: e1=(a,b), e2=(b,c)
                        new_points.append(x)
                else:
                    drop_points.append(x)

            print(f"new: {new_points}, drop: {drop_points}")
            for x in new_points[:-1]:
                current_layer = self._add_interim_point(x, self.curr_state,
                                                        self.curr_state + [x],
                                                        current_layer)
                varnames.append(x)
                self.curr_state += [x]

            x = new_points[-1]
            if cavenum < len(caves)-2:
                # processing the last point
                newstate = [y for y in (self.curr_state + [x])
                            if y not in drop_points]
                current_layer = self._add_interim_point(
                    x, self.curr_state, newstate,
                    current_layer,
                    add_costs=cave)
            else:
                # that's the last point of the pre-last cloud
                # connecting everything to T and F nodes
                # (no new nodes needed, just calculate costs)
                print("Linking the last layer")
                for state in current_layer:
                    for (xval, arc) in [(True, "hi"), (False, "lo")]:
                        fixed_vals = {
                            pt: state[i]
                            for i, pt in enumerate(self.curr_state)}

                        fixed_vals.update({new_points[-1]: xval})
                        cost = self._calc_cave(cave,
                                               fixed_vals) + self._calc_cave(
                                                   caves[-1], fixed_vals)

                        self.B.link(current_layer[state].id, NTRUE, arc, cost)

            varnames.append(x)
            self.curr_state = newstate

        print(f"Renaming varnames from {self.B.vars} to {varnames}")
        self.B.rename_vars({i: varnames[i-1] for i in self.B.vars})
        return self.B
