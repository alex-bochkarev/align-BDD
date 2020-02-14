# Branch and bound type algorithm implementation
# A. Bochkarev, 2019

import varseq as vs
#import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import copy as copy
import math as math
import anytree as at
from anytree.exporter import DotExporter
import myheap as mh
import heapq as heap
from experiments.misc import log
import heuristics as heu

# Selected algorithm parameters
TIMEOUT_ITERATIONS = 100000
UB_UPDATE_FREQ = 5000 # in iterations
ALWAYS_FAST = True

######################################################################
## HEURISTICS for UB calculation

def fix_effect(A,B,a,pos):
    return A.align_to(B.layer_var).size() - A.align_to(B.slide(a,pos).layer_var).size() - B.S(a,pos)

def stable_algt(A,B,algt):

    Ap = A.align_to(algt)
    Bp = B.align_to(algt)

    for i in range(len(A)):
        for j in range(len(A)):
            if fix_effect(A,Bp,A.layer_var[i],j) > 0 or fix_effect(B,Ap,B.layer_var[i],j) > 0:
                return False

    return True

def greedy_fix(A,B,track=False):
    X = vs.VarSeq(layer_vars = B.layer_var, layer_sizes = [1 for i in range(len(A))])
    tr = [copy.copy(X.layer_var)]
    while not stable_algt(A,B,X.layer_var):
        for i in range(len(A)):
            for j in range(len(A)):
                if fix_effect(B,A.align_to(X.layer_var),B.layer_var[i],j) > 0:
                    X = X.slide(B.layer_var[i],j)
                    tr.append(X.layer_var)

        for i in range(len(A)):
            for j in range(len(A)):
                if fix_effect(A,B.align_to(X.layer_var),A.layer_var[i],j) > 0:
                    X = X.slide(A.layer_var[i],j)
                    tr.append(X.layer_var)

    if track:
        return tr
    else:
        return [A.align_to(X.layer_var), B.align_to(X.layer_var)]

## a simpler approach -- neighborhood search
def nsearch(A,B,track=False):
    X = vs.VarSeq(layer_vars = B.layer_var, layer_sizes = [1 for i in range(len(A))])
    tr = [copy.copy(X.layer_var)]
    done = False
    while not done:
        cur_size = A.align_to(X).size() + B.align_to(X).size()
        for i in range(len(X)):
            for j in range(len(X)):
                Xc = X.slide(X.layer_var[i])

## END of HEURISTICS section
######################################################################

######################################################################
## Lower bound calculations
def LB_current(A,B):
    """current size (before the alignment)"""
    return A.size()+B.size()

def LB_first_aligned(A,B):
    """min current size after aligning the first element only"""
    return min(A.size()+B.slide(A.layer_var[0],0).size(), A.slide(B.layer_var[0],0).size()+B.size())

def LB_last_aligned(A,B):
    """min current size after aligning the last element (enumeration)"""
    N = len(A)-1
    return min([A.slide(A.layer_var[i],N).size() + B.slide(A.layer_var[i],N).size() for i in range(N+1)])

def LB_by_level(A,B):
    N = len(A)
    LB = LB_current(A,B)

    for i in range(N-1):
        for j in range(i+1,N):
            Ai = A.layer_var[i]; Aj = A.layer_var[j]
            if B.p[Ai] > B.p[Aj]:
                # it is an inversion
                LB += min(
                    2*A.n[i] - A.n[i+1],
                    2*B.var_size(Aj) - B.n[ B.p[Aj]+1 ]
                )

    return LB

def LB_by_level_complicated(A,B):
    N = len(A)
    LB = LB_current(A,B)

    for i in range(N-1):
        Bs = []
        for j in range(i+1, N):
            Ai = A.layer_var[i]; Aj = A.layer_var[j]
            if B.p[Ai] > B.p[Aj]:
                # an inversion
                Bs.append(2*B.var_size(Aj) - B.n[ B.p[Aj]+1 ])
        B_srt = -1*np.sort(-1 * np.array(Bs))
        A_srt = np.array([A.n[i] * (2**k) - A.n[i+k] for k in range(1,len(Bs)+1)])
        LB += np.sum(np.minimum(A_srt,B_srt))
    return LB

def LB_lvl_compl_symm(A,B):
    return max(LB_by_level_complicated(A,B), LB_by_level_complicated(B,A))

def LB_lvl_symm(A,B):
    return max(LB_by_level(A,B), LB_by_level(B,A))

# CODE - FUNC - LEGEND
LOWER_BOUNDS = [
    ["LB_first",LB_first_aligned,"min size first element aligned"],
    ["LB_last",LB_last_aligned,"min size last element aligned"],
    ["LB_levels",LB_by_level,"inversions-driven LB"]
#    ["LB_levels_amd",LB_by_level_complicated,"inversions-driven LB (amd)"],
#    ["LB_lvl_symm",LB_lvl_symm,"inversion-driven (symmetric)"]
#    ["LB_lvl_symm_amd",LB_lvl_compl_symm,"inversion-driven (amd, symm)"]
]
######################################################################
## auxiliary functions

###
# DOT graph (tree) export-specific
# adding labels (for DOT export)
def nodenamefunc(node):
    fixed_A_start = node.A_tail_start #len(node.A)-node.depth
    fixed_B_start = node.B_tail_start #len(node.B)-node.depth
    return "{}[{}]\nA:{}{}({: >3d})\nB:{}{}({: >3d})\n|A|+|B|={}, LB:{}, UB:{}".format(node.name, node.status, node.A.layer_var[:fixed_A_start],node.A.layer_var[fixed_A_start:],node.A.size(), node.B.layer_var[:fixed_B_start],node.B.layer_var[fixed_B_start:], node.B.size(),node.size(), node.LB, node.UB)

def nodeattrfunc(node):
    nlabel = "xlabel=\"Tree: LB={}, UB={}, obj={}\"".format(node.tree_LB,node.tree_UB, node.best_obj)
    ncolor = ""

    if node.status == "P":
        ncolor = "fillcolor=grey, style=filled"
    if node.status == "?":
        ncolor = "fillcolor=lightgrey, style=dashed"
    if node.status == "O":
        ncolor = "fillcolor=orange, style=filled"
    if node.status == "T":
        ncolor = "penwidth=5"

    return "{},{}".format(nlabel,ncolor)
#
###

######################################################################
## class SearchNode
## node-level search logic
## node status chars: I - init-d, ? - open, E - expanded, P - pruned,
##                    T - terminal
class SearchNode(at.NodeMixin):
    def __init__(self, name, parent, Aresid, Bresid, A_tail_start, B_tail_start,A,B):
        self.name = name
        self.parent = parent
        self.Ar = copy.deepcopy(Aresid)
        self.Br = copy.deepcopy(Bresid)
        self.A_tail_start = A_tail_start
        self.B_tail_start = B_tail_start
        self.A = copy.deepcopy(A)
        self.B = copy.deepcopy(B)
        self.status = "I"

        # tree-specific notes
        self.tree_UB = None
        self.tree_LB = None
        self.best_obj = None
        self.LB = None
        self.UB = None

    def size(self):
        return self.A.size() + self.B.size()

    # upper bound at the current node
    def calculate_UB(self, t='fast'):
        if self.status == "T":
            self.UB = self.A.size()+self.B.size()
            order = self.A.layer_var
        else:
            if t=='fast':
                self.UB, order = heu.minAB(self.A,self.B)
            else:
                self.UB, order = heu.simpl_greedy_swaps(self.A,self.B)

        return [self.UB, order]

    # lower bound at the current node
    def calculate_LB(self):
        """returns a lower bound for the current search tree node"""
        if self.status=="T":
            self.LB = self.size()
        else:
            self.LB = LB_by_level(self.A,self.B) # inversions-driven bound

        return self.LB


    # a special method that tells which node to choose
    # if lower bounds are equal
    # approach: break ties arbitrarily
    def __lt__(self, other):
        return np.random.ranf() > 0.5

######################################################################
## class BBSearch
## data structure to keep the search tree and the tree-level logic
class BBSearch:
    # technical: tree initialization procedure
    def __init__(self, A,B):
        self.A = A
        self.B = B
        self.step = 0
        self.tree_size = 1
        self.verbose = False


        self.logs_UB = []
        self.logs_LB = []
        self.logs_step = []

        self.logging = False

        self.status = "initialized"

        ## create the root node
        self.root = SearchNode("root", None, A,B,len(A),len(B),A,B)

        self.LB = self.root.calculate_LB()
        self.UB, order = self.root.calculate_UB('fast')
        self.Ap_cand = self.A.align_to(order)
        self.Bp_cand = self.B.align_to(order)
        self.node_cand = None
        self.cand_parent = self.root

        self.open_nodes = []

        if self.LB == self.UB:
            # a hypothetical case of good bounds
            opt_node = SearchNode("step {}, node {}: UB=LB".format(self.step, self.tree_size),self.root,[],[],0,0,self.A.align_to(order), self.B.align_to(order))
            opt_node.status = "O"
            self.status = "optimal" # no need to proceed
        else:
            self.open_nodes = [(self.LB, self.root)]
            heap.heapify(self.open_nodes) # create a heap, key = lower bound

    def current_best(self):
        """returns current upper bound (best objective seen)"""

        return self.Ap_cand.size()+self.Bp_cand.size()

    # technical: dump the tree into a DOT-file (graphivz)
    def dump(self, filename=None):
        if filename and self.root:
            DotExporter(self.root, nodenamefunc=nodenamefunc,nodeattrfunc = nodeattrfunc).to_dotfile(filename)

    def make_graph_gap(self, ax=None, trueOpt=None):
        """figure: produces an LB/UB vs step no. figure"""
        ax = ax or plt.gca()

        graph_df = pd.DataFrame(list(zip(self.logs_step, self.logs_UB, self.logs_LB)))
        graph_df.columns = ['step','UB','LB']

        sns.lineplot(x='step', y ='UB', data=graph_df, ax = ax, color="red",markers=True, ci=None)
        sns.lineplot(x='step', y ='LB', data=graph_df, ax = ax, color="blue", markers=True,ci=None )

        if trueOpt:
            ax.axhline(trueOpt,ls="--")
        # Add titles (main and on axis)
        ax.set_xlabel("step number (node expansions)")
        ax.set_ylabel("Upper (red) and lower (blue) bounds")
        ax.set_title("LB/UB gap")
        ax.text(ax.get_xlim()[1]*0.05,ax.get_ylim()[1]*0.8,"A:\n{}".format(self.A))
        ax.text(ax.get_xlim()[1]*0.55,ax.get_ylim()[1]*0.8,"B:\n{}".format(self.B))

    def set_logging(self, logfile, prefix, steps_list):
        """sets up the logging process for BB search"""
        self.logging = True
        self.log_file = logfile
        self.log_prefix = prefix
        self.snapshot_steps = steps_list

    ## MAIN LOGIC
    ## search stop criterion: timeout or gap closed
    def stop_criterion(self):
        """Stop criterion: timeout, no nodes left to process or LB=UB"""

        if self.step > TIMEOUT_ITERATIONS:
            self.status="timeout"
            return True

        if len(self.open_nodes)==0 or (self.current_best() == self.LB):
            if self.verbose:
                print("|open nodes|={}, curr_best({}) = LB({}): {}".format(len(self.open_nodes),self.current_best(), self.LB, (self.current_best() == self.LB)))
            self.status = "optimal"
            return True

        return False

    ## !!!
    ## the main search procedure
    def search(self):
        # calculate the initial UB and LB
        if self.status=="optimal" or len(self.open_nodes)==0:
            if self.verbose:
                print("No nodes to explore / optimal state reached")

            return self.status

        if self.verbose: print("Initial LB={}".format(self.LB))

        self.step = 0

        while not self.stop_criterion():
            node_LB, next_node = self.open_nodes.pop() # pick a node with the smallest LB
            self.expand(next_node, node_LB, self.step)     # expand a node
            self.step += 1
            if self.verbose:
                print("node expanded. Step {}, LB={}, UB={}".format(self.step, self.LB, self.UB))

            if self.logging:
                if self.step in self.snapshot_steps:
                    log(self.log_prefix, str(self.step), str(self.LB), str(self.UB), outfile = self.log_file)
        if (self.node_cand is None) and not (self.cand_parent is None):
            # create an optimal node
            order = self.Ap_cand.layer_var
            opt_node = SearchNode("step {}, node {}: UB=LB".format(self.step, self.tree_size),self.cand_parent,[],[],0,0,self.A.align_to(order), self.B.align_to(order))
            opt_node.status = "O"
            if self.verbose:
                print("optimal node created")

        if self.verbose:
            print("\nSearch done in {} iterations. UB={}, LB={}, status: {}, tree size: {}".format(self.step, self.UB, self.LB, self.status, self.tree_size))
            if self.node_cand is None:
                print("No candidate node found!")

        if self.status == "optimal":
            # self.LB = self.current_best()
            if not self.node_cand is None:
                self.node_cand.status="O"

            while len(self.open_nodes) > 0:
                node_LB, next_node = self.open_nodes.pop()
                if next_node.status == "?":
                    next_node.status = "P"

        if self.logging:
            self.logs_LB.append(self.LB)
            self.logs_UB.append(self.UB)
            self.logs_step.append(self.step)
            log(self.log_prefix, str(self.step), str(self.LB), str(self.UB), outfile = self.log_file, comment = "status:{}".format(self.status))

        return self.status

    ## "expanding" a search tree node
    def expand(self, node, node_LB, step=0):
        """Process the node =node="""
        if self.verbose:
            print("Node expansion: {}".format(node.name))

        # set tree bounds to be shown on the graph
        # (sample BB search tree figure)
        node.tree_LB = self.LB
        node.tree_UB = self.UB
        node.best_obj = self.current_best()

        if node_LB > self.UB:
            # we can prune the node
            node.status ="P"
            return

        node.status = "E" # actually processing the node

        for i in reversed( range(len(node.Ar)) ):
            # enumerating all the possible last elements
            # starting from the last one
            if vs.non_dominated(node.Ar.layer_var[i],node.Ar,node.Br):
                ## add a node
                ## sift the element under consideration to the last position
                A_new = copy.deepcopy(node.A)
                B_new = copy.deepcopy(node.B)
                a = node.Ar.layer_var[i]
                pos = len(node.Ar)-1
                A_new = A_new.slide(a,pos)
                B_new = B_new.slide(a,pos)
                ## reshuffle the rest ``covered'' elements
                xA = node.Ar.layer_var
                xB = node.Br.layer_var
                A_covered_els = xA[i:]
                i_A = i
                if A_covered_els != []:
                    for j in range(len(xB)):
                        if xB[j] in A_covered_els and xB[j]!=xA[i]:
                            A_new.layer_var[i_A] = xB[j]
                            A_new.p[xB[j]] = i_A
                            i_A += 1

                B_covered_els = xB[node.B.p[xA[i]]:]
                if B_covered_els != []:
                    i_B = node.B.p[xA[i]]
                    for j in range(len(xA)):
                        if xA[j] in B_covered_els and j!=i:
                            B_new.layer_var[i_B] = node.A.layer_var[j]
                            B_new.p[node.A.layer_var[j]] = i_B
                            i_B += 1

                ## save the results / create a node
                Ar_new = vs.VarSeq(A_new.layer_var[:len(node.Ar)-1], A_new.n[:len(node.Ar)-1])
                Br_new = vs.VarSeq(B_new.layer_var[:len(node.Br)-1], B_new.n[:len(node.Br)-1])

                newnode = SearchNode("step {}, node {}: {}to{}".format(step, self.tree_size, node.Ar.layer_var[i],len(node.Ar)-1), node, Ar_new, Br_new, pos, pos,A_new, B_new)

                self.tree_size += 1

                ## check if the new one is a terminal node
                if A_new.is_aligned(B_new):
                    newnode.status = "T"
                else:
                    newnode.status = "?" # node to be expanded

                ## calculate the node lower bound
                newnode.LB = newnode.calculate_LB()

                ## update UPPER BOUND
                t = 'fast'
                if step % UB_UPDATE_FREQ == 0 and not ALWAYS_FAST:
                    t = 'slow'

                newnode.UB, order = newnode.calculate_UB(t)

                if self.UB > newnode.UB:
                    self.UB = newnode.UB

                if self.current_best() > newnode.UB:
                    self.Ap_cand = self.A.align_to(order)
                    self.Bp_cand = self.B.align_to(order)
                    if newnode.status == "T":
                        self.node_cand = newnode
                    else:
                        self.cand_parent = newnode

                if self.LB < self.UB and newnode.status != "T":
                    heap.heappush(self.open_nodes,(newnode.calculate_LB(), newnode))

        ## update LOWER BOUND
        if self.open_nodes:
            self.LB = min(self.current_best(), self.open_nodes[0][0])
        else:
            self.LB = self.current_best()

        if self.logging:
            self.logs_step.append(step)
            self.logs_UB.append(self.UB)
            self.logs_LB.append(self.LB)

    ## auxiliary functions
