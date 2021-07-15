"""Aux script: draws a sample B&B search tree
for the auxiliary problem (with :py:class:`varseq.VarSeq`)
"""
# import timeit

import varseq as vs
import BB_search as BB
import argparse as ap

def main():
    """Generates a random align-VS instance and saves the BB search tree."""
    parser = ap.ArgumentParser(description=" Produces a sample BB search tree (a .dot file) (c) A. Bochkarev, Clemson University, 2020",
                               formatter_class=ap.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-n", "--nodes", action="store", dest="n", help="max. number of search tree nodes ",
                        type=int,default=15)
    parser.add_argument("-V", "--vars", action="store", dest="V", help="no. of variables",
                        type=int,default=15)
    parser.add_argument("-o","--output", dest="filename", help=".dot filename to save")
    parser.add_argument("-v","--verbose",dest="verbose", action="store_true")

    args = parser.parse_args()

    tree_size = args.n*2
    while tree_size > args.n:
        A = vs.VarSeq.random(N = args.V)
        B = vs.VarSeq.random(N = args.V)
        b = BB.BBSearch(A,B)
        b.verbose = args.verbose
        status = b.search()
        tree_size = b.tree_size

    b.dump(args.filename)
    print("Instance generated:")
    print(A)
    print(B)

if __name__ == "__main__":
    main()
