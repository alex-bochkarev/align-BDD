"""
Creates a log for scaling test
(comparing runtimes for different number of variables)

(c) A. Bochkarev, Clemson University, 2020
abochka@clemson.edu
"""

import BDD as exact
import varseq as simpl
from time import time
import BB_search as bb
import heuristics as heu
import numpy as np
import sys
from copy import deepcopy
import argparse as ap
from experiments.misc import log

if __name__ == "__main__":
    parser = ap.ArgumentParser(description="Processes (solves) a list of align-BDD instances (c) A. Bochkarev, Clemson University, 2020",
                               formatter_class=ap.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-i", "--in", action="store", dest="inst_list", help="filename for instances list to process")
    parser.add_argument("-o", "--out", action="store", dest="logfile", help="filename for the log")
    parser.add_argument("-d", "--directory", action="store", dest="inst_dir", help="directory with instances")
    parser.add_argument("-H", "--header", action="store_true", dest="header", help="show header only and exit")
    args = parser.parse_args()

    if args.header:
        log("N","instance","num_type","value",comment="comment")
        # needed for the legend
        log("-1,-1,legend,orig_simpl", comment="Simplified problem-based")
        log("-1,-1,legend,orig_gsifts1p", comment="Greedy sifts")
        exit(0)

    inst_dir = args.inst_dir

    with open(args.inst_list,"r") as inst_list:
        with open(args.logfile,"w") as logf:
            # process instances
            for inst_id in inst_list:
                inst_id = inst_id.rstrip()
                if inst_id == "":
                    continue

                fnameA = "".join([inst_dir, "A",inst_id,".bdd"])
                fnameB = "".join([inst_dir, "B",inst_id,".bdd"])
                bdd_A = exact.BDD(); bdd_A.load(fnameA)
                bdd_B = exact.BDD(); bdd_B.load(fnameB)

                ######################################################
                ## The simplified problem-based heuristic

                t0=time()
                vs_A = simpl.VarSeq(bdd_A.vars, [len(l) for l in bdd_A.layers[:-1]])
                vs_B = simpl.VarSeq(bdd_B.vars, [len(l) for l in bdd_B.layers[:-1]])

                N = len(vs_A)

                # find opts for varseq instance with BB-search
                b = bb.BBSearch(vs_A,vs_B)
                status = b.search()
                o = b.Ap_cand.layer_var # optimal order found by the BBSearch
                orig_simpl_obj = bdd_A.align_to(o).size() + bdd_B.align_to(o).size()
                t1 = time()

                log(N,inst_id,"simpl_opt_status",1-int(status=="optimal"),outfile=logf, comment=status)
                log(N,inst_id, "orig_simpl_time",t1-t0, outfile=logf)
                log(N,inst_id, "orig_simpl_obj",orig_simpl_obj, outfile=logf)

                ######################################################
                ## Greedy sifts heuristic
                t0 = time()
                A = deepcopy(bdd_A)
                B = deepcopy(bdd_B)

                A.gsifts(B)
                if not A.is_aligned(B):
                    A = bdd_A; B = bdd_B
                    B.gsifts(A)

                t1 = time()
                if A.is_aligned(B):
                    gs_status="aligned"
                else:
                    gs_status="non-aligned"

                log(N,inst_id, "gsifts_status",gs_status, outfile=logf)
                log(N,inst_id, "orig_gsifts1p_time",t1-t0, outfile=logf)
                log(N,inst_id, "orig_gsifts1p_obj",A.size()+B.size(), outfile=logf)
