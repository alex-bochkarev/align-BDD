"""An experiment for joint UFLP + special instance type (cavemen).

A version with CPP MIP.
"""
from darkcloud import gen_caveman_inst
from UFLP_fullDD import create_cover_DD
from UFLPOrder import UFLP_greedy_order

from jUFLP_cavemen import gen_cavemen_jUFLP_inst, solve_cm_jUFLP_MIP
from jUFLP_cavemen import solve_cm_jUFLP_CPPMIP_fullDDs
from jUFLP_cavemen import save_inst, solve_cm_jUFLP_fullDDs
from time import time
import numpy as np

def main():
    print("GreedyOrderSize, RandomSize")
    M = 14
    L = 0.35
    n = 1

    i = 1

    for _ in range(100):
        S, f, c, caves = gen_caveman_inst(n, M, L)

        BG, _ = create_cover_DD(S, f, c, UFLP_greedy_order(S))
        BR, _ = create_cover_DD(S, f, c, list(np.random.permutation([j for j in range(1, len(S)+1)])))

        print(f"{BG.size()}, {BR.size()}", flush=True)

        i += 1


if __name__ == '__main__':
    main()
