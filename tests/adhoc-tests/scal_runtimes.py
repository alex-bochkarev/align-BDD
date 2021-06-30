"""
Aux script: tests for a runtime of a single
instance in a `scal_test`

(c) A. Bochkarev, Clemson University, 2020
abochka@clemson.edu
"""
import timeit
import sys
import numpy as np
import cProfile
import pstats
import BDD as exact
import BB_search as BB
import varseq as vs

sys.path.append('..')

def summarize_runtime(N):
    SETUP_CODE = '''import varseq as vs
import BDD as exact
import BB_search as BB
import numpy as np'''

    TEST_CODE='''N = {}
A = exact.BDD.random(N=N)
B = exact.BDD.random(N=N)
B.rename_vars(dict(zip([i for i in range(1,N+1)],np.random.permutation([i for i in range(1,N+1)]))))
vsA = vs.VarSeq(A.vars, [len(l) for l in A.layers[:-1]])
vsB = vs.VarSeq(B.vars, [len(l) for l in B.layers[:-1]])
b = BB.BBSearch(vsA,vsB)
b.verbose = True
status = b.search()
order = b.Ap_cand.layer_var
target_obj = A.align_to(order).size()+B.align_to(order).size()
A.gsifts(B)
'''.format(N)
    times = timeit.repeat(setup = SETUP_CODE,
                        stmt = TEST_CODE,
                        repeat = 10,
                                number = 1)
    mean = np.mean(times)
    print("for N={}, mean time is {} s. per instance".format(N, mean))
    print("which makes for {} minutes for 250 instances".format(mean*250/60))
    print("All times for {} runs: {}".format(len(times),times))

## profiling
N=25
pr = cProfile.Profile()
pr.enable()
A = exact.BDD.random(N=N)
B = exact.BDD.random(N=N)
B.rename_vars(dict(zip([i for i in range(1,N+1)],np.random.permutation([i for i in range(1,N+1)]))))
vsA = vs.VarSeq(A.vars, [len(l) for l in A.layers[:-1]])
vsB = vs.VarSeq(B.vars, [len(l) for l in B.layers[:-1]])
b = BB.BBSearch(vsA,vsB)
status = b.search()
order = b.Ap_cand.layer_var
target_obj = A.align_to(order).size()+B.align_to(order).size()
A.gsifts(B)
pr.disable()
ps = pstats.Stats(pr).strip_dirs().sort_stats('time')
ps.print_stats()
