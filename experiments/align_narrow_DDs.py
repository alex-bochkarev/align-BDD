"""Trying to figure if the heuristic performs worse for narrow DDs."""
import numpy as np
from BDD import BDD, intersect
from varseq import VarSeq
from BB_search import BBSearch

def main():
    print("p=pA/pB: N, A_width, B_width, sim_AB, intVS_size, intA_size")
    N = 20
    pA = 0.80
    pB = 0.15

    for _ in range(15):
        A = BDD.random(N=N, p=pA, weighted=True)
        B = BDD.random(N=N, p=pB, weighted=True)
        B.rename_vars(dict(zip([i for i in range(1, N+1)],
                            np.random.permutation([i for i in range(1, N+1)]))))

        A.make_reduced()
        B.make_reduced()

        Bpp = B.align_to(A.vars, inplace=False)

        intA = intersect(A, Bpp)
        intA.make_reduced()

        vs_A = VarSeq(A.vars, [len(L) for L in A.layers[:-1]])
        vs_B = VarSeq(B.vars, [len(L) for L in B.layers[:-1]])

        assert set(vs_A.layer_var) == set(
            vs_B.layer_var), f"A:{vs_A.layer_var}, B:{vs_B.layer_var}"

        b = BBSearch(vs_A, vs_B)

        # bb.TIMEOUT_ITERATIONS=5000
        status = b.search()
        assert status == "optimal" or status == "timeout", f"status: {status}"

        Ap = A.align_to(b.Ap_cand.layer_var, inplace=False)
        Bp = B.align_to(b.Ap_cand.layer_var, inplace=False)

        intVS = intersect(Ap, Bp)
        intVS.make_reduced()

        print(f"p={pA:.2f}/{pB:.2f}: {N}, {A.width()}, {B.width()}, {A.simscore(B):.2f}, {intVS.size()}, {intA.size()}")

if __name__ == '__main__':
    main()
