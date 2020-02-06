"""
Experiment: Branch-and-bound algorithm ``convergence'':
logs upper and lower bounds at certain step numbers

(c) Alexey Bochkarev, Clemson University, 2020
abochka@clemson.edu
"""

import varseq as vs
import BDD as exact
import BB_search as bb
import sys
import os
from experiments.misc import log
from time import time
import argparse as ap

snapshot_steps = [i*10 for i in range(1,100)]

if __name__ == "__main__":
    parser = ap.ArgumentParser(description="Creates B&B logs for the simplified problem. (c) A. Bochkarev, Clemson University, 2020",
                               formatter_class=ap.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-d","--dir", dest="inst_dir", help="directory to take instances from",default="./instances/raw/reduced/")
    parser.add_argument("-i","--input_list", dest="inst_list", help="list of instance IDs to process",default="./instances/raw/reduced/instances.list")
    parser.add_argument("-o","--out_file", dest="logfile", help="output log file name",default="./BB_bounds.log")
    parser.add_argument("-H", "--header", action="store_true", dest="header", help="show header only and exit")

    args = parser.parse_args()
    if args.header:
        log("instance","num_type","step","LB","UB",comment="comment")
        exit(0)

    inst_dir = args.inst_dir

    if not os.path.isdir(inst_dir):
        print("Error: '{}' -- not a directory".format(out_dir))
        exit(1)

    with open(args.inst_list,"r") as inst_list:
        with open(args.logfile,"w") as logf:
            log("instance","num_type","step","LB","UB",comment="comment",outfile = logf)
            for inst_id in inst_list:
                inst_id = inst_id.rstrip()
                if inst_id == "":
                    continue

                # construct a corresponding varseq instance
                fnameA = "".join([inst_dir, "A",inst_id,".bdd"])
                fnameB = "".join([inst_dir, "B",inst_id,".bdd"])
                bdd_A = exact.BDD(); bdd_A.load(fnameA)
                bdd_B = exact.BDD(); bdd_B.load(fnameB)

                ## run some heuristics for the simplified problem
                vs_A = vs.VarSeq(bdd_A.vars, [len(l) for l in bdd_A.layers[:-1]])
                vs_B = vs.VarSeq(bdd_B.vars, [len(l) for l in bdd_B.layers[:-1]])

                ## run a logged bb-search
                t0 = time()
                b = bb.BBSearch(vs_A,vs_B)
                b.set_logging(logf, "{},steplog".format(inst_id), snapshot_steps)

                o = b.search()
                t1 = time()
                log(inst_id,"timelog",b.step,t1-t0,0, outfile = logf)
