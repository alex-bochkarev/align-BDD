"""Generates and solves the special class j-UFLP instances.

Instances and the experiment as discussed in the paper, in Section 4.2 and
Appendix F.2. Instances are generated using
:py:func:`UFLP_2_cav.gen_special_jUFLP`, which is a wrapper for function
:py:func:`darkcloud.gen_caveman_inst` for this class of instances. The
experiment is implemented in :py:func:`UFLP_2_cav.main`.

The rest of the code implemented alternative experiments (left out from the
second revision of the paper).

"""
import numpy as np
from copy import copy
from time import time

from experiments.softcover import generate_overlaps
from jUFLP_cavemen import solve_cm_jUFLP_MIP, solve_cm_jUFLP_CPPMIP_fullDDs
from jUFLP_cavemen import solve_cm_jUFLP_fullDDs
from darkcloud import gen_caveman_inst
from jUFLP_utils import save_inst

from BDD import simscore

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
      The parameter ``L`` assumed to be:

      .. math::
          L = 1 - 2\\frac{\\textrm{Number of existing arcs}}{N(N-1)}

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


def make_cluster_reverse_custom_matching(ca1, ca2, ss):
    """Makes clusters matching given simscore parameter `ss`."""
    link = dict()
    clusters = [k for k in range(len(ca2))]  # cluster-to-cluster matching

    for k in range(len(ca1)):
        # within-cluster matching
        assert len(ca1[k]) == len(ca2[clusters[k]])

        src_order = [j for j in range(len(ca1[k]))]
        dest_order = copy(src_order)

        while simscore(src_order, dest_order) > ss:
            # Make a random swap: introduce a single inversion
            swapped = False
            while not swapped:
                swap_pos = np.random.randint(1, len(dest_order))
                if dest_order[swap_pos] > dest_order[swap_pos-1]:
                    s1, s2 = dest_order[swap_pos], dest_order[swap_pos-1]
                    dest_order[swap_pos] = s2
                    dest_order[swap_pos-1] = s1
                    swapped = True

        link.update(dict(zip(ca1[k],
                                [ca2[clusters[k]][j] for j in dest_order])))

    return link


def gen_special_jUFLP(n, M, L, linking="consecutive", inst_type="cavemen",
                      param=None):
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
            - "cluster-reverse-custom": consecutive clusters, in the orders
              with a given similarity score (up to an error caused by the
              integer number of possible inversions).
            - "literal": trivial linking, 1-to-1.

        inst_type (str): instance topology, one of:
            - "1-link": ... (cluster) - (node) -(cluster) ...
            - "cavemen": ... (cluster) - (cluster) ...

        param (float): instance generation parameter for
        "cluster-reverse-custom" linking, simscore tier number.
        (1.0 corresponds to 100% simscore, 0.0 -- to 0%)

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
        i2 = [i2[0], i2[1], i2[2], caves2]
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
        clusters = [k for k in range(len(ca2))]
        for k in range(len(ca1)):
            link.update(dict(zip(ca1[k],
                                 reversed(ca2[clusters[k]]))))

        in_clusters1 = np.unique(sum([], ca1))
        in_clusters2 = np.unique(sum([], ca2))

        origins = [j for j in range(1, len(i1[0])+1)
                   if j not in in_clusters1]

        targets = [j for j in range(1, len(i2[0])+1)
                   if j not in in_clusters2]

        link.update(dict(zip(origins, targets)))
    elif linking == "cluster-reverse-custom":
        ca1 = [S for S in i1[COL_caves]]
        ca2 = [S for S in i2[COL_caves]]

        link = make_cluster_reverse_custom_matching(ca1, ca2, param)
    elif linking == "literal":
        link = dict(zip([j for j in range(1, len(i1[0])+1)],
                        [j for j in range(1, len(i2[0])+1)]))
    else:
        raise ValueError(f"Wrong value of `linking`({linking})" +
                         "expected 'uniform,' 'by-cluster,' or 'consecutive")

    return i1, i2, link


def main():
    """Implements the j-UFLP experiment discussed in Section 4.2."""
    print("experiment, n, M, L, N, A, inst_type, linking, tMIP, tMIP_CPP, tDD_VS, tDD_toA, int_VS, int_VS_toA")
    M = 15
    L = 0.35
    n = 2
    linking = "cluster-reverse"
    inst_type = "cavemen"

    k = 0
    for i in range(1, 100+1):
        for M in [10, 12, 15]:
            k += 1
            i1, i2, jm = gen_special_jUFLP(n, M, L, linking, inst_type)
            save_inst(i1, i2, jm, f"instances/jUFLP_cm/inst_{k}.json")

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

            assert abs(objMIP - objMIP_CPP) < 0.01, f"objMIP = {objMIP:.2f}, objMIP_CPP={objMIP_CPP:.2f}"
            assert abs(objMIP - objDD_VS) < 0.01, f"objMIP = {objMIP:.2f}, objDD_VS={objDD_VS:.2f}"
            assert abs(objMIP - objDD_toA) < 0.01, f"objMIP = {objMIP:.2f}, objDD_toA={objDD_toA:.2f}"

            A = sum([len(s)-1 for s in i1[0]])/2+sum([len(s)-1 for s in i2[0]])/2
            print(f"{k}, {n}, {M}, {L}, {len(i1[0])+len(i2[0])}, {A}, " +
                  f"{inst_type}, {linking}, " +
                  f"{tMIP:.2f}, {tMIP_CPP:.2f}, {tDD_VS:.2f}, {tDD_toA:.2f}, {int_VS}, {int_VS_toA}",
                  flush=True)


if __name__ == '__main__':
    main()
