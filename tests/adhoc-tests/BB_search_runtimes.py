"""
Aux script: draws a sample B&B search tree
for the auxiliary problem (with varseqs)

(c) A. Bochkarev, Clemson University, 2019
abochka@clemson.edu
"""
import timeit
import sys
sys.path.append('..')

SETUP_CODE = '''import varseq as vs
import BB_search as BB'''

TEST_CODE='''A = vs.VarSeq.random(N = 10)
B = vs.VarSeq.random(N = 10)
b = BB.BBSearch(A,B)
status = b.search()'''

print(timeit.repeat(setup = SETUP_CODE,
                    stmt = TEST_CODE,
                    repeat = 5,
                    number = 20))
