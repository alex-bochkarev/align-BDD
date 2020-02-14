"""
Aux script: draws a sample B&B search tree
for the auxiliary problem (with varseqs)

(c) A. Bochkarev, Clemson University, 2019
abochka@clemson.edu
"""
import timeit
import sys
sys.path.append('..')

SETUP_CODE = '''import varseq as vs'''

TEST_CODE_R='''A = vs.VarSeq.random(N = 15)
Ap = A.align_to([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])'''

TEST_CODE_Q='''A = vs.VarSeq.random(N = 15)
Ap = A.q_align_to([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])'''

print("Old version:")
print(timeit.repeat(setup = SETUP_CODE,
                    stmt = TEST_CODE_R,
                    repeat = 3,
                    number = 1000))

print("Quick version:")

print(timeit.repeat(setup = SETUP_CODE,
                    stmt = TEST_CODE_Q,
                    repeat = 3,
                    number = 1000))
