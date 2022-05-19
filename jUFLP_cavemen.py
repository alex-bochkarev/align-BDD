"""A specialized code for typed UFLP over cavemen instances."""
from copy import copy, deepcopy
import numpy as np

from darkcloud import ptscloud, generate_overlaps


def gen_caveman_inst(n=10, M=5, L=0.5, verbose=False):
    """Generates an instance with the related metadata (info on caves).

    Args:
      n (int): number of caves,
      M (int): number of points in a cave,
      L (float): edge sparsity parameter (share of missing edges)
      verbose (Bool): print debug info

    Returns:
      S, S2, f, f2, c, caves: instance and caves description.

    Note:
      The implementation is based on :py:func:`darkcloud.gen_caveman_inst`,
      but with the 'types' condition replaced with another cover-like
      condition. Note that caves are the same for the two networks.
    """
    # creating the first graph topology
    S = [[j+1] for j in range(n * M)]  # creating the nodes first
    entry_point = (None, None)
    caves = []
    for k in range(n):
        lastcave = [k*M + 1]
        if len(caves) > 0:
            caves[-1].e2 = entry_point
            S[entry_point[0]-1].append(entry_point[1])
            S[entry_point[1]-1].append(entry_point[0])

        n_edges = 0
        # create the necessary number of connected nodes
        while len(lastcave) < M:
            connect_to = np.random.choice(lastcave)
            lastcave.append(lastcave[-1]+1)
            S[lastcave[-1]-1].append(connect_to)
            S[connect_to-1].append(lastcave[-1])
            n_edges += 1

        # add edges to ensure the required sparsity
        while (1 - 2*n_edges / (M*(M-1))) > L:
            n1 = np.random.choice(lastcave)
            n2 = np.random.choice(lastcave)

            if n1 not in S[n2-1]:
                S[n2-1].append(n1)
                S[n1-1].append(n2)
                n_edges += 1

        caves.append(ptscloud(entry_point,
                              (None, None),
                              [j for j in range(k*M+1, (k+1)*M+1)]))

        entry_point = (np.random.choice(lastcave), (k+1)*M+1)

    # creating the second graph
    caves2 = deepcopy(caves)
    S2 = copy(S)

    n_edges = 0
    for cave in caves2:
        if cave.e1 != (None, None) and cave.e1[1] not in S2[cave.e1[0]-1]:
            S2[cave.e1[0]-1].append(cave.e1[1])
            S2[cave.e1[1]-1].append(cave.e1[0])

        if cave.e2 != (None, None) and cave.e2[1] not in S2[cave.e2[0]-1]:
            S2[cave.e2[0]-1].append(cave.e2[1])
            S2[cave.e2[1]-1].append(cave.e2[0])

        n_edges = 0
        while (1 - 2 * n_edges / (M * (M - 1))) > L:
            n1 = np.random.choice(cave.S)
            n2 = np.random.choice(cave.S)
            if n1 not in S2[n2-1]:
                S2[n2-1].append(n1)
                S2[n1-1].append(n2)
                n_edges += 1

    # creating costs info (location and overlap costs)
    f = generate_overlaps(S)
    f2 = generate_overlaps(S2)

    Cmax = 5.0
    c = [Cmax*np.random.uniform() for _ in range(len(S))]
    if verbose:
        print(f"S={S};\nf={f}\n;c={c}")
    return S, S2, f, f2, c, caves
