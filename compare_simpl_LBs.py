"""
Compares different lower bounds for the simplified problem.
(given a list of problems and a solution log)

(c) A. Bochkarev, Clemson University, 2020
abochka@clemson.edu
"""

import varseq as simpl
import BDD as exact
import BB_search as bbs
import argparse as ap
from experiments.misc import log

if __name__ == "__main__":
    parser = ap.ArgumentParser(description="Compares various LB performance for the simplified problem. (c) A. Bochkarev, Clemson University, 2020",
                               formatter_class=ap.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-d","--dir", dest="inst_dir", help="temporary output directory")
    parser.add_argument("-l","--list_inst", dest="inst_list", help="list of instance IDs (filename)")
    args = parser.parse_args()

    log("LB,gap",comment="legend")
    with open(args.inst_list,"r") as inst_list:
        for inst_id in inst_list:
            inst_id = int(inst_id.rstrip())
            A = exact.BDD(); B = exact.BDD()
            A.load(args.inst_dir+"/A{}.bdd".format(inst_id))
            B.load(args.inst_dir+"/B{}.bdd".format(inst_id))

            N = len(A.vars)

            vsA = simpl.VarSeq(A.vars, [len(l) for l in A.layers[:-1]])
            vsB = simpl.VarSeq(B.vars, [len(l) for l in B.layers[:-1]])
            bb = bbs.BBSearch(vsA,vsB)
            bb.search()

            curr_size = vsA.size()+vsB.size()
            opt_size = bb.Ap_cand.size()+bb.Bp_cand.size()
            if opt_size == curr_size:
                continue # not an interesting instance

            for LB in bbs.LOWER_BOUNDS:
                log(LB[0],(LB[1](vsA,vsB) - curr_size)/(opt_size - curr_size), comment=LB[2])

