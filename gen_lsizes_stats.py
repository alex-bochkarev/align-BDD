"""
Auxiliary script: generates a dataset for BDDs
(given a list of BDD filenames)

(c) A. Bochkarev, Clemson University, 2020
abochka@clemson.edu
"""

import BDD as exact
from experiments.misc import log
import argparse as ap
import sys

if __name__ == "__main__":
    parser = ap.ArgumentParser(description="Produces a layer widths statistics (dataset). (c) A. Bochkarev, Clemson University, 2020",
                               formatter_class=ap.ArgumentDefaultsHelpFormatter)
    parser.add_argument(dest="filename", help="list of instance names to process (filenames with paths)")
    parser.add_argument("-s","--suffix",help="a suffix to append to each row (such as p value)",
                        action="store")
    parser.add_argument("-H", "--header", action="store_true", dest="header", help="show header only and exit")

    args = parser.parse_args()

    if args.header:
        log("# instance","N", "is_reduced", comment="<comma-separated layer-widths>,<suffix>")
        exit(0)

    bdd_no = 0
    with open(args.filename,"r") as inst_list:
        for inst_name in inst_list:
            inst_name = inst_name.rstrip()
            bdd = exact.BDD()
            bdd.load(inst_name)
            N = len(bdd.vars)
            print(",".join([str(bdd_no), str(N), str(bdd.is_reduced())] + [str(bdd.n(i)) for i in range(N)]+[str(args.suffix)]))
            bdd_no += 1
