""" A branch-and-bound search scheme
implementation for the align-sequences problem:
alignts two :py:class:`varseq.VarSeq` objects.

Tests coverage: :py:mod:`BB_search_test`.
"""
import varseq as vs
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import copy as copy
import anytree as at
from anytree.exporter import DotExporter
import heapq as heap
from experiments.misc import log
import heuristics as heu

# Selected algorithm parameters
TIMEOUT_ITERATIONS = 1000
UB_UPDATE_FREQ = 5000 # in iterations
ALWAYS_FAST = True

######################################################################
## Lower bound calculations
def LB_current(A,B):
    """LB: current size (before the alignment)."""
    return A.size()+B.size()

def LB_first_aligned(A,B):
    """LB: min current size after aligning the first element only."""
    return min(A.size()+B.slide(A.layer_var[0],0).size(), A.slide(B.layer_var[0],0).size()+B.size())

def LB_last_aligned(A,B):
    """LB: min current size after aligning the last element (enumeration)."""
    N = len(A)-1
    return min([A.slide(A.layer_var[i],N).size() + B.slide(A.layer_var[i],N).size() for i in range(N+1)])

def LB_by_level(A,B):
    """LB: Inversions-driven heuristic.

    The one based on the lemma presented in the paper.
    """
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
    """LB: Another variant of the inversions-driven heuristic (improved)"""
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
    """LB: Symmetric version of the previous one."""
    return max(LB_by_level_complicated(A,B), LB_by_level_complicated(B,A))

def LB_lvl_symm(A,B):
    """Symmetric version of the :py:func:`LB_by_level`."""

    return max(LB_by_level(A,B), LB_by_level(B,A))

# List of lower bounds to examine
# CODE - FUNC - LEGEND
LOWER_BOUNDS = [
    ["LB_first",LB_first_aligned,"min size first element aligned"],
    ["LB_last",LB_last_aligned,"min size last element aligned"],
    ["LB_levels",LB_by_level,"inversions-driven LB"]
#    ["LB_levels_amd",LB_by_level_complicated,"inversions-driven LB (amd)"],
#    ["LB_lvl_symm",LB_lvl_symm,"inversion-driven (symmetric)"]
#    ["LB_lvl_symm_amd",LB_lvl_compl_symm,"inversion-driven (amd, symm)"]
]
"""List of lower bounds to examine (used in a separate experiment)."""

######################################################################
## auxiliary functions

## DOT graph (tree) export-specific
## adding labels (for DOT export)
def nodenamefunc(node):
    """Helper: Generates node names for graphviz export."""

    fixed_A_start = node.A_tail_start #len(node.A)-node.depth
    fixed_B_start = node.B_tail_start #len(node.B)-node.depth
    return "{}[{}]\nA:{}{}({: >3d})\nB:{}{}({: >3d})\n|A|+|B|={}, LB:{}, UB:{}".format(node.name, node.status, node.A.layer_var[:fixed_A_start],node.A.layer_var[fixed_A_start:],node.A.size(), node.B.layer_var[:fixed_B_start],node.B.layer_var[fixed_B_start:], node.B.size(),node.size(), node.LB, node.UB)

def nodeattrfunc(node):
    """Helper: Generates node attributes for graphviz export."""

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

######################################################################
class SearchNode(at.NodeMixin):
    """Implements search tree node.

    Attributes:
        name (str): node name,
        parent (:py:class:`SearchNode`): parent node,
        A, B (:py:class:`varseq.VarSeq`): sequences to align,
        Ar, Br (:py:class:`varseq.VarSeq`): not-yet-aligned parts of ``A`` and ``B``.
        A_tail_start, B_tail_start (int) : position for the already-aligned tail start.
        status (str): node type (see the note below),
        tree_UB, tree_LB (int): current (tree) upper / lower bound,
        best_obj (int): current best objective seen,

    Note:

        - possible node types are: ``T`` = Terminal, ``O`` = Optimal,
          ``I`` = initialized, ``P`` = pruned, ``?`` = unknown
          ('open', not processed), ``E`` = processed ('expanded')

    """
    def __init__(self, name, parent, Aresid, Bresid, A_tail_start, B_tail_start,A,B):
        """Class constructor"""
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
        self.UB = None

    def size(self):
        """Returns node size (total no. of nodes in both diagrams)."""
        return self.A.size() + self.B.size()

    # upper bound at the current node
    def calculate_UB(self, t='fast'):
        """Returns an upper bound for the specific node."""
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
        """Returns a lower bound for the current search tree node."""
        if self.status=="T":
            LB = self.size()
        else:
            LB = LB_by_level(self.A,self.B) # inversions-driven bound

        return LB


    # a special method that tells which node to choose
    # if lower bounds are equal
    # approach: break ties arbitrarily
    def __lt__(self, other):
        """A special methods used to break ties (if LBs are equal) -- randomly."""
        return np.random.ranf() > 0.5

######################################################################
class BBSearch:
    """Implements the search tree.

    Attributes:
        A, B (:py:class:`varseq.VarSeq`): first sequence,
        step (int): current step number,
        tree_size (int): current number of nodes in the tree,
        verbose (Boolean): debug information printing flag,
        logging (Boolean): log bounds flag,
        logs_UB, logs_LB (list): upper/lower bounds by step (if ``verbose``),
        logs_step (list): step numbers (if ``verbose``),
        status (str): search status,
        root (class SearchNode): root node,
        Ap_cand, Bp_cand (:py:class:`varseq.VarSeq`): ``A`` / ``B`` aligned to the current candidate optimum,
        node_cand, cand_parent (:py:class:`SearchNode`): node corresponding to the current candidate and its parent,
        open_nodes (list): list of nodes to process during the search,
            each entry is a tuple (``<lower bound>``,``SearchNode``).
    """

    def __init__(self, A,B):
        """Class constructor"""
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

        self.UB, order = self.root.calculate_UB('fast')
        self.Ap_cand = self.A.align_to(order)
        self.Bp_cand = self.B.align_to(order)
        self.node_cand = None
        self.cand_parent = self.root

        LB = self.root.calculate_LB()

        self.open_nodes = [(LB, self.root)]
        heap.heapify(self.open_nodes) # create a heap, key = lower bound

    def current_best(self):
        """Returns current upper bound (best objective seen)."""

        return self.Ap_cand.size()+self.Bp_cand.size()

    def dump(self, filename=None):
        """Technical: dump the tree into a DOT-file (graphivz)."""
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
        """Sets up the logging process for BB search."""
        self.logging = True
        self.log_file = logfile
        self.log_prefix = prefix
        self.snapshot_steps = steps_list

    ## MAIN LOGIC
    def search(self):
        """Performs BB search (the main procedure).

        Note:
            After the call, the following attributes allow to recover the
            solution:

                - ``status`` (str): ``optimal``, ``timeout`` (``initialized``
                  before),
                - ``Ap_cand``, ``Bp_cand``
                  (:py:class:`varseq.VarSeq`): ``A`` and ``B`` after aligning
                  to the candidate (optimum if ``status`` = ``optimal``).
        """

        # calculate the initial UB and LB
        if self.status=="optimal" or len(self.open_nodes)==0:
            if self.verbose:
                print("No nodes to explore / optimal state reached")

            return self.status

        if self.verbose: print("Initial LB={}".format(LB))

        self.step = 0

        while (len(self.open_nodes) > 0 and self.step <= TIMEOUT_ITERATIONS):
            LB, node = heap.heappop(self.open_nodes) # pick a node with the smallest LB
            ######################################################################
            # expanding the node

            if self.verbose:
                print("Node expansion: {}, node LB={}, UB={}".format(node.name, LB, node.UB))

            # set tree bounds to be shown on the graph
            # (sample BB search tree figure)
            node.tree_LB = LB
            node.tree_UB = self.UB
            node.best_obj = self.current_best()

            if LB >= self.UB:
                # we can prune the node
                node.status ="P"
                self.status="optimal"
                LB = self.UB
                break

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
                    A_new.slide(a,pos,inplace=True)
                    B_new.slide(a,pos,inplace=True)
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

                    newnode = SearchNode("step {}, node {}: {}to{}".format(self.step, self.tree_size, node.Ar.layer_var[i],len(node.Ar)-1), node, Ar_new, Br_new, pos, pos,A_new, B_new)

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
                    if self.step % UB_UPDATE_FREQ == 0 and not ALWAYS_FAST:
                        t = 'slow'

                    newnode.UB, order = newnode.calculate_UB(t)

                    if self.UB > newnode.UB:
                        self.UB = newnode.UB
                        self.Ap_cand = self.A.align_to(order)
                        self.Bp_cand = self.B.align_to(order)
                        if newnode.status == "T":
                            self.node_cand = newnode
                        else:
                            self.cand_parent = newnode

                    if newnode.LB < self.UB and newnode.status != "T":
                        heap.heappush(self.open_nodes,(newnode.LB, newnode))

            ## update LOWER BOUND

            if self.logging:
                self.logs_step.append(self.step)
                self.logs_UB.append(self.UB)
                self.logs_LB.append(LB)
            # end of node expansion
            ######################################################################
            self.step += 1
            if self.verbose:
                print("node expanded. Step {}, LB={}, UB={}".format(self.step, LB, self.UB))

            if self.logging:
                if self.step in self.snapshot_steps:
                    log(self.log_prefix, str(self.step), str(LB), str(self.UB), outfile = self.log_file)

        if self.step > TIMEOUT_ITERATIONS:
            self.status="timeout"
        elif (self.node_cand is None) and not (self.cand_parent is None):
            # create an optimal node
            order = self.Ap_cand.layer_var
            opt_node = SearchNode("step {}, node {}: UB=LB".format(self.step, self.tree_size),self.cand_parent,[],[],0,0,self.A.align_to(order), self.B.align_to(order))
            opt_node.status = "O"
            if self.verbose:
                print("optimal node created")
            self.status="optimal"
            LB = self.UB
        elif (self.node_cand is None and self.cand_parent is None):
            print("ERROR: optimum found, but both node candidate and candidate parent are None (please report a bug)")
            print("Optimal order: {}, objective={}".format(self.Ap_cand.layer_var, self.current_best()))
            print("...while aligning: \nA:\n{}\nB:\n{}".format(self.A,self.B))
            exit(1)
        else:
            self.status="optimal"
            LB = self.UB

        if self.verbose:
            print("\nSearch done in {} iterations. UB={}, LB={}, status: {}, tree size: {}".format(self.step, self.UB, LB, self.status, self.tree_size))
            if self.node_cand is None:
                print("No candidate node found!")

        if self.status == "optimal":
            if not self.node_cand is None:
                self.node_cand.status="O"

            while len(self.open_nodes) > 0:
                node_LB, next_node = self.open_nodes.pop()
                if next_node.status == "?":
                    next_node.status = "P"

        if self.logging:
            self.logs_LB.append(LB)
            self.logs_UB.append(self.UB)
            self.logs_step.append(self.step)
            log(self.log_prefix, str(self.step), str(LB), str(self.UB), outfile = self.log_file, comment = "status:{}".format(self.status))

        return self.status
