"""An experiment for joint UFLP + special instance type (cavemen).

A version with CPP MIP.
"""
from jUFLP_cavemen import gen_cavemen_jUFLP_inst, solve_cm_jUFLP_MIP
from jUFLP_cavemen import solve_cm_jUFLP_CPPMIP_fullDDs
from jUFLP_cavemen import save_inst, solve_cm_jUFLP_fullDDs
from time import time


def main():
    print("experiment, n, M, L, tMIP, tMIP_CPP, tDD_VS, int_VS")
    M = 5
    L = 0.50

    i = 1

    while True:
        for n in range(3, 8):
            i1, i2, jm = gen_cavemen_jUFLP_inst(n, M, L)
            # save_inst(i1, i2, jm, f"instances/jUFLP_cm/inst_wCPPMIP_{i}.json")

            print("---")
            t0 = time()
            objMIP = solve_cm_jUFLP_MIP(i1, i2, jm)
            tMIP = time() - t0
            print(f"✅ MIP in {tMIP:.2f} sec")

            # solve with CPP MIP
            t0 = time()
            objMIP_CPP = solve_cm_jUFLP_CPPMIP_fullDDs(i1, i2, jm)
            tMIP_CPP = time() - t0
            print(f"✅ CPP MIP in {tMIP_CPP:.2f} sec")

            t0 = time()
            objDD3, int_VS = solve_cm_jUFLP_fullDDs(i1, i2, jm, "VS", True)
            tDD_VS = time() - t0
            print(f"✅ Full DDs VS in {tDD_VS:.2f} sec")

            assert abs(objMIP - objMIP_CPP) < 0.01, f"objMIP = {objMIP:.2f}, objMIP_CPP={objMIP_CPP:.2f}"
            assert abs(objMIP - objDD3)< 0.01, f"objMIP = {objMIP:.2f}, objDD3={objDD3:.2f}"


            print(f"{i}, {n}, {M}, {L}, {tMIP:.2f}, {tMIP_CPP:.2f}, {tDD_VS:.2f}, {int_VS}", flush=True)

            i += 1


if __name__ == '__main__':
    main()
