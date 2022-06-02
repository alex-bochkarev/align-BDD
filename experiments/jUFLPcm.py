"""An experiment for joint UFLP + special instance type (cavemen)."""
from jUFLP_cavemen import gen_cavemen_jUFLP_inst, solve_cm_jUFLP_MIP
from jUFLP_cavemen import solve_cm_jUFLP_DDs, save_inst
from time import time


def main():
    print("experiment, n, M, L, tMIP, tDD_toA, tDD_toB, tDD_VS, int_toA, int_toB, int_VS")
    M = 5
    L = 0.25

    i = 1

    while True:
        for n in range(15, 18):
            i1, i2, jm = gen_cavemen_jUFLP_inst(n, M, L)
            t0 = time()
            objMIP = solve_cm_jUFLP_MIP(i1, i2, jm)
            tMIP = time() - t0

            save_inst(i1, i2, jm, f"instances/jUFLP_cm/inst_{i}.json")

            t0 = time()
            objDD, int_toA = solve_cm_jUFLP_DDs(i1, i2, jm, "toA", True)
            tDD_toA = time() - t0

            t0 = time()
            objDD2, int_toB = solve_cm_jUFLP_DDs(i1, i2, jm, "toB", True)
            tDD_toB = time() - t0

            t0 = time()
            objDD3, int_VS = solve_cm_jUFLP_DDs(i1, i2, jm, "VS", True)
            tDD_VS = time() - t0

            assert abs(objMIP - objDD)< 0.01
            assert abs(objMIP - objDD2)< 0.01
            assert abs(objMIP - objDD3)< 0.01

            print(f"{i}, {n}, {M}, {L}, {tMIP:.2f}, {tDD_toA:.2f}, {tDD_toB:.2f}, {tDD_VS:.2f}, {int_toA}, {int_toB}, {int_VS}", flush=True)

            i += 1


if __name__ == '__main__':
    main()
