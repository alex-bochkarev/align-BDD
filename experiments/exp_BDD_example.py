import BDD
import itertools
from copy import deepcopy

if __name__ == '__main__':
    A = BDD.BDD()
    A.load("experiments/exp_BDD.bdd")
    A.dump_gv().save("./run_logs/exp_BDD_example/exp_BDD.gv")

    perms = itertools.permutations(A.vars[0:3])
    B = deepcopy(A)
    for v in perms:
        B.align_to(list(v)+A.vars[3:], inplace=True)
        B.show(dir="./run_logs/exp_BDD_example", filename="_".join(list(v)+A.vars[3:])+".gv")
        print(f"v={list(v)+A.vars[3:]}, |B|={B.size()} (reduced={B.is_reduced()}), cross-check: {A.is_equivalent(B)[0]}")

    B.show()
