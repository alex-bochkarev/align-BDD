"""Investigates a special class of instances: randomized cover DD."""

from time import time
import numpy as np

from darkcloud import gen_typed_cavemen_inst
from darkcloud import DDSolver, DDTypedSolver

from BDD import intersect
from varseq import VarSeq
from BB_search import BBSearch


def main():
    """Generates a dataset of BDD vs MIP runtimes and objectives."""
    print("exp_num, n, M, L, K_types, kmax, gen_iters, objU, objT, t_novsA, " +
          "tTDD, size_int_VS, sim_C_VS, size_int_toC, sim_C_toC")
    for s in range(1, 250):
        generated = False
        gen_iters = 0
        while not generated:
            gen_iters += 1
            n = 10
            M = 10
            L = 0.45
            K = 3
            kbmax = 3
            S, f, c, caves, k, kbar = gen_typed_cavemen_inst(n, M, L, K, kbmax)

            # first: solve the untyped version
            t0 = time()
            sol = DDSolver(S, f, c, caves)
            B = sol.build_cover_DD()
            sp = B.shortest_path()
            objU = sp[0]
            tUDD = time() - t0

            # second: solve the typed version
            t0 = time()
            sol = DDTypedSolver(S, f, c, caves, k, kbar)
            ###

            C = sol.build_cover_DD()
            T, _ = sol.build_type_DD()

            C.make_reduced()
            T.make_reduced()

            # randomize the order of the cover DD -- to remove special
            # structure, so it would not be necessarily optimal
            # to align type DD to cover DD.
            C.align_to(np.random.permutation(C.vars), inplace=True)

            vs_C = VarSeq(C.vars, [len(L) for L in C.layers[:-1]])
            vs_T = VarSeq(T.vars, [len(L) for L in T.layers[:-1]])

            assert set(vs_C.layer_var) == set(
                vs_T.layer_var), f"T:{vs_T.layer_var}, C:{vs_C.layer_var}"

            b = BBSearch(vs_C, vs_T)

            # bb.TIMEOUT_ITERATIONS=5000
            status = b.search()
            assert status == "optimal" or status == "timeout", f"status: {status}"

            Cp = C.align_to(b.Ap_cand.layer_var, inplace=False)
            Tp = T.align_to(b.Ap_cand.layer_var, inplace=False)

            int_DD = intersect(Cp, Tp)
            # assert int_DD.is_reduced(), "not reduced int DD!"
            objT = int_DD.shortest_path()[0]
            ###
            tTDD = time() - t0

            size_int_VS = int_DD.size()
            sim_C_VS = Cp.simscore(C)

            if abs(objU - objT) > 1:
                generated = True  # type constraints are binding

        t0 = time()
        sol = DDTypedSolver(S, f, c, caves, k, kbar)
        ###
        # C = sol.build_cover_DD()
        # T, _ = sol.build_type_DD()

        # T.make_reduced()
        # C.make_reduced()

        Tpp = T.align_to(C.vars, inplace=True)

        int_DD = intersect(C, Tpp)
        objT_novsA = int_DD.shortest_path()[0]
        ###
        tT_novsA = time() - t0
        assert abs(objT - objT_novsA) < 0.01

        size_int_toC = int_DD.size()
        sim_C_toC = T.simscore(C)

        # t0 = time()
        # sol = DDTypedSolver(S, f, c, caves, k, kbar)
        # objT_novsB = sol.solve_with_DDs_noVS(TtoC=False)
        # tT_novsB = time() - t0
        # assert abs(objT - objT_novsB) < 0.01

        print(f"{s}, {n}, {M}, {L}, {K}, {kbmax}, {gen_iters}," +
              f" {objU:.2f}, {objT:.2f}, {tT_novsA:.2f}, {tTDD:.2f}," +
              f" {size_int_VS}, {sim_C_VS:.2f}, {size_int_toC}, {sim_C_toC:.2f}",
              flush=True)


if __name__ == '__main__':
    main()
