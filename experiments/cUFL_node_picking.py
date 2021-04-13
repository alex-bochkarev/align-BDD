"""
An experiment to compare different ways of picking nodes
in (colored) UFL (random, min, and max degrees).

---
(c) A. Bochkarev, Clemson University, 2021
abochka@clemson.edu
"""

import cUFL
import seaborn as sns
import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
from copy import deepcopy


def get_sizes(n=10, K=1):
    """Gets cover DD sizes (list of lists, of size (3, `K`)) for given `n` variables."""
    DD_sizes = pd.DataFrame(columns = ['min', 'max', 'rnd'])

    for k in range(K):
        S, f, fc, kb = cUFL.generate_test_instance(n)
        cov_min, _ = cUFL.build_cover_DD(S, f)
        cm = deepcopy(cov_min)
        cov_min.make_reduced()
        assert cm.is_equivalent(cov_min)[0], f"S={S}; f={f} # cov_min problem"

        cov_max, _ = cUFL.build_cover_DD(S, f, next_node_type="max")
        cov_max.make_reduced()
        assert cov_max.is_equivalent(cov_min)[0], f"S={S}; f={f} # cov_max problem"

        cov_rnd, _ = cUFL.build_cover_DD(S, f, next_node_type="rnd")
        cov_rnd.make_reduced()
        assert cov_rnd.is_equivalent(cov_min)[0], f"S={S}; f={f} # cov_rnd problem"

        DD_sizes = DD_sizes.append({'min':cov_min.size(), 'max': cov_max.size(), 'rnd': cov_rnd.size()}, ignore_index=True)
        if k % 10 == 0:
            print(".", end="")
            sys.stdout.flush()

    return DD_sizes

if __name__ == '__main__':
    DDs = get_sizes(n=7, K=1000)
    DDs["max_rel"] = DDs["max"] / DDs["min"]
    DDs["rnd_rel"] = DDs["rnd"] / DDs["min"]

    sns.histplot(DDs, x="max_rel", binwidth=0.1)
    print(f"No of cases where max is better than min: {np.sum(DDs['max_rel'] < 1.0)} out of {len(DDs['max_rel'])}")
    plt.show()


    sns.histplot(DDs, x="rnd_rel", binwidth=0.1)
    plt.show()
    print(f"No of cases where rnd is better than min: {np.sum(DDs['rnd_rel'] < 1.0)} out of {len(DDs['rnd_rel'])}")
