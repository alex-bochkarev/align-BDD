"""
Aux script: draws a sample B&B search tree
for the auxiliary problem (with varseqs)

(c) A. Bochkarev, Clemson University, 2019
abochka@clemson.edu
"""

import varseq as vs
import BB_search as BB
from importlib import reload

reload(BB)
A = vs.VarSeq.random(N = 7)
B = vs.VarSeq.random(N = 7)

print(B)

b = BB.BBSearch(A,B)
b.verbose = True

status = b.search()

b.dump("./run_logs/sample_BB_tree.dot")

with open("./run_logs/sample_BB_instance.txt","w") as f:
    f.write("A:\n")
    f.write(str(A))
    f.write("\nB:\n")
    f.write(str(B))
