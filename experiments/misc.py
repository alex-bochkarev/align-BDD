"""
Contains common misc utilities for numerical experiments
(mostly run on the cluster)

(c) A. Bochkarev, Clemson University, 2019
abochka@clemson.edu
"""
import sys

def log(instance, *what_to_print, outfile = sys.stdout, comment="--none--"):
    outfile.write(str(instance))
    for val in list(what_to_print):
        outfile.write(",")
        outfile.write(str(val))

    outfile.write(",{}\n".format(comment))
    if outfile == sys.stdout:
        outfile.flush()


######################################################################
## module testing code
#
# log(5,"hey",0.3,0.4,0.5,0.7)
# with open("./test.txt","w") as f:
#     log(5,"hey",0.3,0.4,0.5,0.7,comment="mycomm",outfile = f)
