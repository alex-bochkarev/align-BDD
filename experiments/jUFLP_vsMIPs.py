"""An experiment for joint UFLP + special instance type (cavemen).

A version with CPP MIP.
"""
from jUFLP_cavemen import gen_cavemen_jUFLP_inst, solve_cm_jUFLP_MIP
from jUFLP_cavemen import solve_cm_jUFLP_CPPMIP_fullDDs
from jUFLP_cavemen import save_inst, solve_cm_jUFLP_fullDDs
from time import time

def main():
    print("experiment, n, M, L, N, A, tMIP, tMIP_CPP, tDD_VS, int_VS")
    M = 11
    L = 0.35
    n = 2

    i = 1

    for _ in range(100):
        i1, i2, jm = gen_cavemen_jUFLP_inst(n, M, L)
        save_inst(i1, i2, jm, f"instances/jUFLP_cm/inst_wMIP_{i}.json")

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
        objDD3, int_VS = solve_cm_jUFLP_fullDDs(i1, i2, jm, "VS", True)
        tDD_VS = time() - t0
        print(f"✅ Full DDs VS in {tDD_VS:.2f} sec", flush=True)

        assert abs(objMIP - objMIP_CPP) < 0.01, f"objMIP = {objMIP:.2f}, objMIP_CPP={objMIP_CPP:.2f}"
        assert abs(objMIP - objDD3)< 0.01, f"objMIP = {objMIP:.2f}, objDD3={objDD3:.2f}"


        A = sum([len(s)-1 for s in i1[0]])/2+sum([len(s)-1 for s in i2[0]])/2
        print(f"{i}, {n}, {M}, {L}, {len(i1[0])+len(i2[0])}, {A}, {tMIP:.2f}, {tMIP_CPP:.2f}, {tDD_VS:.2f}, {int_VS}", flush=True)

        i += 1


if __name__ == '__main__':
    main()
