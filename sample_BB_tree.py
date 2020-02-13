"""
Aux script: draws a sample B&B search tree
for the auxiliary problem (with varseqs)

(c) A. Bochkarev, Clemson University, 2019
abochka@clemson.edu
"""
import timeit

SETUP_CODE = '''import varseq as vs
import BB_search as BB
from importlib import reload'''

TEST_CODE='''
A = vs.VarSeq.random(N = 8)
B = vs.VarSeq.random(N = 8)
b = BB.BBSearch(A,B)
status = b.search()
'''
#b.verbose = True

print(timeit.repeat(setup = SETUP_CODE,
                          stmt = TEST_CODE,
                          repeat = 5,
                          number = 10))

# b.dump("./run_logs/sample_BB_tree.dot")

# with open("./run_logs/sample_BB_instance.txt","w") as f:
#     f.write("A:\n")
#     f.write(str(A))
#     f.write("\nB:\n")
#     f.write(str(B))
