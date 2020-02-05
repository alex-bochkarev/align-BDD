"""
Aux script: generates random align-BDD instances
(``original problems'') given dataset parameters
and prints key instance parameters to the stdout (in a .csv format):
    instance     -- instance ID
    inversions   -- number of inversions between var(A) and var(B)
    reduced[A/B] -- whether A (resp., B) is layer-reduced
    gen_trial    -- number of trials until a unique entry is generated

(c) A. Bochkarev, Clemson University, 2020
abochka@clemson.edu
"""

import BDD as exact
import numpy as np
import sys
import os
import argparse as ap

def count_inversions(X, Y):
    invs = 0
    p = dict()

    # O(n): index all the vars of Y
    for i in range(len(Y)):
        p.update({Y[i] : i})

    # O(n^2): check all pairs
    # NOTE: O(n log n) is possible with MergeSort

    for i in range(len(X)-1):
        for j in range(i+1,len(X)):
            if p[X[i]] >= p[X[j]]:
                invs += 1

    return invs


if __name__ == "__main__":
    parser = ap.ArgumentParser(description="Generates align-BDD instances. (c) A. Bochkarev, Clemson University, 2020",
                               formatter_class=ap.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-n", "--no_instances", action="store", dest="n", help="number of instances to generate",
                        type=int,default=1)
    parser.add_argument("-v","--no_variables",action="store", dest="V",
                        help="number of variables per instance",
                        type=int, default=15)
    parser.add_argument("-s","--start_id", action="store",dest="start_id",
                        help="start instance ID number to use",
                        type=int,default=1)
    parser.add_argument("-p", "--prob_exp", dest="p",help="tree expansion probability parameter",
                        action="store", type=float, default=0.6)
    parser.add_argument("-U","--unique", help="generate unique instances",action="store_true")
    parser.add_argument("out_dir", help="output directory")

    group = parser.add_mutually_exclusive_group()
    group.add_argument("-R","--reduced", help="generate reduced instances",action="store_true")
    group.add_argument("-N","--nonreduced", help="generate reduced instances",action="store_true")

    args = parser.parse_args()

    start_id = args.start_id
    n = args.n
    N = args.V
    p = args.p
    create_reduced = args.reduced
    create_unique = args.unique
    out_dir = args.out_dir

    if not os.path.isdir(out_dir):
        print("Error: '{}' -- not a directory".format(out_dir))
        exit(1)

    print("instance,inversions,reducedA,reducedB,gen_trial")

    inst_profiles = set()

    for inst_id in range(start_id, start_id+n):
        # create an exact instance
        inst_accepted = False
        trials = 0
        while not inst_accepted:
            trials += 1
            bdd_A = exact.BDD.random(N=N,p=p)
            bdd_B = exact.BDD.random(N=N,p=p)
            bdd_B.rename_vars(dict(zip([i for i in range(1,N+1)],np.random.permutation([i for i in range(1,N+1)]))))

            if create_reduced:
                bdd_A.make_reduced()
                bdd_B.make_reduced()

            if create_unique:
                # check if the instance is unique / never seen before
                prof_A = bdd_A.profile()
                prof_B = bdd_B.profile()
                if not ( ( (prof_A,prof_B) in inst_profiles ) or ( (prof_B, prof_A) in inst_profiles )):
                    inst_profiles.add((prof_A,prof_B))
                    inst_accepted = True
            else:
                inst_accepted = True

        bdd_A.save(out_dir+"/A{}.bdd".format(inst_id))
        bdd_B.save(out_dir+"/B{}.bdd".format(inst_id))

        print("{},{},{},{},{}".format(inst_id,count_inversions(bdd_A.vars, bdd_B.vars), bdd_A.is_reduced(), bdd_B.is_reduced(), trials))
        sys.stdout.flush() # otherwise we are risking not to have anyting on job kill...
