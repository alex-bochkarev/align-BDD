from time import time
from darkcloud import gen_typed_cavemen_inst, solve_typed_with_MIP
from darkcloud import solve_with_MIP, DDSolver, DDTypedSolver


def main():
    """Generates a dataset of BDD vs MIP runtimes and objectives."""
    print("exp_num, n, M, L, K_types, kmax, gen_iters, objU, objT, t_novsA, " +
          "tTDD")
    for s in range(1, 250):
        generated = False
        gen_iters = 0
        while not generated:
            gen_iters += 1
            n = 30
            M = 3
            L = 0.45
            K = 10
            kbmax = 15
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
            objT = sol.solve_with_DDs()
            tTDD = time() - t0

            if abs(objU - objT) > 1:
                generated = True  # type constraints are binding

        t0 = time()
        sol = DDTypedSolver(S, f, c, caves, k, kbar)
        objT_novsA = sol.solve_with_DDs_noVS(TtoC=True)
        tT_novsA = time() - t0
        assert abs(objT - objT_novsA) < 0.01

        # t0 = time()
        # sol = DDTypedSolver(S, f, c, caves, k, kbar)
        # objT_novsB = sol.solve_with_DDs_noVS(TtoC=False)
        # tT_novsB = time() - t0
        # assert abs(objT - objT_novsB) < 0.01

        print(f"{s}, {n}, {M}, {L}, {K}, {kbmax}, {gen_iters}," +
              f" {objU:.2f}, {objT:.2f}, {tT_novsA:.2f}, {tTDD:.2f}",
              flush=True)


if __name__ == '__main__':
    main()
