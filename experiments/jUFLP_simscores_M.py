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
from time import time

from jUFLP_cavemen import solve_cm_jUFLP_MIP, solve_cm_jUFLP_CPPMIP_fullDDs
from jUFLP_cavemen import solve_cm_jUFLP_fullDDs
from jUFLP_utils import save_inst

from BDD import simscore
from UFLP_2_cav import gen_special_jUFLP, make_cluster_reverse_custom_matching
from UFLP_2_cav import COL_caves

from UFLP_fullDD import create_cover_DD
from UFLPOrder import UFLP_greedy_order


def main():
    """The j-UFLP experiment involving instances w/different simscores."""
    print("experiment, n, M, L, N, A, inst_type, linking, param, CPP_simscore, tMIP, tMIP_CPP, tDD_VS, tDD_toA, int_VS, int_VS_toA")
    L = 0.35
    n = 2
    linking = "cluster-reverse-custom"
    inst_type = "cavemen"

    k = 0
    for M in [5, 8, 10, 13, 15]:
        for i in range(1, 10+1):
            i1, i2, jm = gen_special_jUFLP(n, M, L, linking, inst_type, 0.0)
            for param in [0.00, 0.25, 0.50, 0.75, 1.00]:
                k += 1
                ca1 = [S for S in i1[COL_caves]]
                ca2 = [S for S in i2[COL_caves]]
                jm = make_cluster_reverse_custom_matching(ca1, ca2, param)
                save_inst(i1, i2, jm, f"instances/jUFLP_ss_M/inst_{M}_{param}_{k}.json")

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

                # check the resulting CPP BDDs simscore
                S, f, c, caves = i1
                S2, f2, c2, caves2 = i2

                B1, _ = create_cover_DD(S, f, c, UFLP_greedy_order(S, True))
                B2, _ = create_cover_DD(S2, f2, c2, UFLP_greedy_order(S2, False))
                B1.make_reduced()
                B2.make_reduced()
                B1.rename_vars(jm)

                print(f"{k}, {n}, {M}, {L}, {len(i1[0])+len(i2[0])}, {A}, " +
                    f"{inst_type}, {linking}, {param}, " +
                    f"{simscore(B1.vars, B2.vars)}, " +
                    f"{tMIP:.2f}, {tMIP_CPP:.2f}, {tDD_VS:.2f}, {tDD_toA:.2f}, {int_VS}, {int_VS_toA}",
                    flush=True)


if __name__ == '__main__':
    main()
