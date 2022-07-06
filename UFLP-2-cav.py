"""Generates cavemen-like instances with node-links.

... -- [ ...-⊙-] -- ⊙ -- [-⊙-... ] -- ...

Therefore, after adding the first point of the next cluster, all the points of
the previous clusters can be 'forgotten'. Theoretically, this should imply
really compact BDD representations, somewhat of :math:`O(n 2^M)`.
"""
import numpy as np
from copy import copy
from time import time

from experiments.softcover import generate_overlaps
from jUFLP_cavemen import solve_cm_jUFLP_MIP, solve_cm_jUFLP_CPPMIP_fullDDs
from jUFLP_cavemen import solve_cm_jUFLP_fullDDs
from darkcloud import gen_caveman_inst

# UFLP instance description structure
# UFLP_inst = (S, f, c, caves)
COL_S = 0
COL_f = 1
COL_c = 2
COL_caves = 3


def gen_nlinks_cavemen_inst(n=10, M=5, L=0.5):
    """Generates an instance with the related metadata (info on caves).

    Args:
      n (int): number of clusters,
      M (int): number of points in a cluster,
      L (float): edge sparsity parameter (share of missing edges)

    Returns:
      S, f, c, caves: instance (S,f,c) and caves description.

    Note:
      The parameter ``L`` is calculated as follows.

      .. math::
          \\Lambda = 1 - 2\\frac{#\\textrm{existing arcs}}{N(N-1)}

      (See, e.g., Sefair and Smith 2016 on DSPI for a similar approach.)

      Caves description in this edition is just a list of lists of points
      (which are within each cluster).
    """
    S = [[j+1] for j in range((n-1) * (M+1) + M)]  # creating the nodes first
    caves = []
    for k in range(n):
        lastcave = [k*(M+1) + 1]

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

        if len(caves) > 0:
            i = np.random.choice(caves[-1])
            j = np.random.choice(lastcave)
            m = k*(M+1)
            S[i-1].append(m)
            S[m-1].append(i)
            S[m-1].append(j)
            S[j-1].append(m)

        caves.append(copy(lastcave))

    # creating costs info (location and overlap costs)
    f = generate_overlaps(S)
    C0 = 200.0
    Cw = 10.0
    c = [C0 + Cw*(np.random.uniform() - 0.5) for _ in range(len(S))]
    return S, f, c, caves


def gen_special_jUFLP(n, M, L, linking="consecutive", inst_type="cavemen"):
    """Generate a special instance of jUFLP.

    Args:
        n, M, L: instance parameters (see :py:func:`gen_nlinks_cavemen_inst`.)
        linking (str): linking constraints type, one of:

            - "uniform": link randomly, regardless of clusters,
            - "consecutive": link randomly within consecutive clusters
              ('ladder linking')
            - "by-cluster": link randomly, but cluster-to-cluster.
            - "cluster-reverse": consecutive clusters, in reverse order
              within each pair of clusters.
            - "literal": trivial linking, 1-to-1.

        inst_type (str): instance topology, one of:
            - "1-link": ... (cluster) - (node) -(cluster) ...
            - "cavemen": ... (cluster) - (cluster) ...

    Returns:
      i1, i2, link: two instances (S,f,c, caves) and linking dict.
    """
    if inst_type == "1-link":
        i1 = gen_nlinks_cavemen_inst(n, M, L)
        i2 = gen_nlinks_cavemen_inst(n, M, L)
    elif inst_type == "cavemen":
        i1 = gen_caveman_inst(n, M, L)
        i2 = gen_caveman_inst(n, M, L)
        caves1 = [ca.S for ca in i1[3]]
        caves2 = [ca.S for ca in i2[3]]
        i1 = [i1[0], i1[1], i1[2], caves1]
        i2 = [i2[0], i2[1], i2[2], caves1]
    else:
        raise ValueError(f"Wrong value of `inst_type`({inst_type})" +
                         "expected: '1-link' or 'cavemen'")

    assert len(i1[COL_S]) == len(i2[COL_S])  # must be the same no. of nodes.

    if linking == "uniform":
        link = dict(zip([j for j in range(1, len(i1[COL_S])+1)],
                        np.random.permutation(
                            [j for j in range(1, len(i2[COL_S])+1)])))
    elif linking == "consecutive":
        ca1 = [S for S in i1[COL_caves]]
        ca2 = [S for S in i2[COL_caves]]

        link = dict()
        for k in range(len(ca1)):
            link.update(dict(zip(ca1[k], np.random.permutation(ca2[k]))))

        in_clusters1 = np.unique(sum([], ca1))
        in_clusters2 = np.unique(sum([], ca2))

        origins = [j for j in range(1, len(i1[0])+1)
                   if j not in in_clusters1]

        targets = [j for j in range(1, len(i2[0])+1)
                   if j not in in_clusters2]

        link.update(dict(zip(origins, targets)))
    elif linking == "by-cluster":
        ca1 = [S for S in i1[COL_caves]]
        ca2 = [S for S in i2[COL_caves]]

        link = dict()
        clusters = np.random.permutation([k for k in range(len(ca2))])
        for k in range(len(ca1)):
            link.update(dict(zip(ca1[k],
                                 np.random.permutation(ca2[clusters[k]]))))

        in_clusters1 = np.unique(sum([], ca1))
        in_clusters2 = np.unique(sum([], ca2))

        origins = [j for j in range(1, len(i1[0])+1)
                   if j not in in_clusters1]

        targets = [j for j in range(1, len(i2[0])+1)
                   if j not in in_clusters2]

        link.update(dict(zip(origins, targets)))
    elif linking == "cluster-reverse":
        ca1 = [S for S in i1[COL_caves]]
        ca2 = [S for S in i2[COL_caves]]

        link = dict()
        # clusters = [k for k in np.random.permutation(range(len(ca2)))]
        # clusters = [k for k in reversed(range(len(ca2)))]
        clusters = [k for k in range(len(ca2))]
        for k in range(len(ca1)):
            # link.update(dict(zip(ca1[k],
            #                      np.random.permutation(ca2[clusters[k]]))))
            link.update(dict(zip(ca1[k],
                                 reversed(ca2[clusters[k]]))))

        in_clusters1 = np.unique(sum([], ca1))
        in_clusters2 = np.unique(sum([], ca2))

        origins = [j for j in range(1, len(i1[0])+1)
                   if j not in in_clusters1]

        targets = [j for j in range(1, len(i2[0])+1)
                   if j not in in_clusters2]

        link.update(dict(zip(origins, targets)))
    elif linking == "literal":
        link = dict(zip([j for j in range(1, len(i1[0])+1)],
                        [j for j in range(1, len(i2[0])+1)]))
    else:
        raise ValueError(f"Wrong value of `linking`({linking})" +
                         "expected 'uniform,' 'by-cluster,' or 'consecutive")

    return i1, i2, link


def draw_jUFLP_inst(i1, i2, link, filename="tmp/jUFLP.dot"):
    """Saves an instance to a ``.dot`` file."""
    with open(filename, "w") as fout:
        fout.write("graph G {\n")

        for (inst, no, pref) in [(i1, 1, 'f'), (i2, 2, 's')]:
            S, f, c, caves = inst
            fout.write(f"    subgraph cluster_{no-1}" +
                       " {\n")  # encoding a sub-instance
            fout.write(f'        color=blue; label="sub-UFLP-{no}";')
            added = set([])
            for i in range(len(S)):
                for j in S[i]:
                    if ((i+1) != j) and not (((j, (i+1)) in added)
                                             or ((i+1, j) in added)):

                        fout.write(f"        {pref}{i+1}--{pref}{j};\n")
                        added.add(((i+1), j))

            fout.write("    };\n")  # end of sub-instance

        for j in link:
            fout.write(f"    f{j} -- s{link[j]}[color=red, style=dashed, penwidth=1];\n")
        fout.write("}")


def shuffle_inst(instance):
    """Shuffles the graph nodes randomly."""
    S, f, c, caves = instance
    inew = [j for j in range(len(S))]
    np.random.shuffle(inew)
    Sp = [[inew[j-1]+1 for j in S[inew[i]]] for i in range(len(S))]
    fp = [f[inew[i]] for i in range(len(S))]
    cp = [c[inew[i]] for i in range(len(S))]
    cavesp = [[inew[j-1] for j in cave] for cave in caves]

    return (Sp, fp, cp, cavesp)


if __name__ == '__main__':
    print("experiment, n, M, L, N, A, inst_type, linking, tMIP, tMIP_CPP, tDD_VS, tDD_toA, int_VS, int_VS_toA")
    M = 14
    L = 0.35
    n = 3
    linking = "cluster-reverse"
    inst_type = "cavemen"

    for i in range(1, 100+1):
        i1, i2, jm = gen_special_jUFLP(n, M, L, linking, inst_type)
        # save_inst(i1, i2, jm, f"instances/jUFLP_cm/inst_wMIP_{i}.json")
        # i2 = shuffle_inst(i2)

        print("---")
        t0 = time()
        objMIP = solve_cm_jUFLP_MIP(i1, i2, jm)
        tMIP = time() - t0
        print(f"✅ MIP in {tMIP:.2f} sec", flush=True)

        # solve with CPP MIP
        t0 = time()
        objMIP_CPP = solve_cm_jUFLP_CPPMIP_fullDDs(i1, i2, jm)
        tMIP_CPP = time() - t0
        print(f"✅ CPP MIP in {tMIP_CPP:.2f} sec", flush=True)

        t0 = time()
        objDD_VS, int_VS = solve_cm_jUFLP_fullDDs(i1, i2, jm, "VS", True)
        tDD_VS = time() - t0
        print(f"✅ Full DDs VS in {tDD_VS:.2f} sec", flush=True)

        t0 = time()
        objDD_toA, int_VS_toA = solve_cm_jUFLP_fullDDs(i1, i2, jm, "toA", True)
        tDD_toA = time() - t0
        print(f"✅ Full DDs toA in {tDD_toA:.2f} sec", flush=True)


        # assert abs(objMIP - objMIP_CPP) < 0.01, f"objMIP = {objMIP:.2f}, objMIP_CPP={objMIP_CPP:.2f}"
        # assert abs(objMIP - objDD_VS) < 0.01, f"objMIP = {objMIP:.2f}, objDD_VS={objDD_VS:.2f}"
        # assert abs(objMIP - objDD_toA) < 0.01, f"objMIP = {objMIP:.2f}, objDD_toA={objDD_toA:.2f}"


        A = sum([len(s)-1 for s in i1[0]])/2+sum([len(s)-1 for s in i2[0]])/2
        print(f"{i}, {n}, {M}, {L}, {len(i1[0])+len(i2[0])}, {A}, " +
              f"{inst_type}, {linking}, " +
              f"{tMIP:.2f}, {tMIP_CPP:.2f}, {tDD_VS:.2f}, {tDD_toA:.2f}, {int_VS}, {int_VS_toA}",
              flush=True)
