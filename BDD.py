"""
BDD.py -- implementation of the BDD-related machinery

Involves actual variable ordering and BDDs alignment classes
and functions (as opposed to the simplification defined
in the `varseq.py`)

(c) A. Bochkarev, Clemson University, 2020
"""

import numpy as np
import pandas as pd
from numpy.random import random as rnd
from random import choice as rnd_choose
from random import randint
import sys

import itertools as iters
from graphviz import Digraph
from collections import deque

import copy
import math

# Module-level constants
GSIFTS_MAX_INCREASE_MUL = 100000  # max increase of the BDD size (e.g., 100 times) during greedy sifts
WEIGHT_TOL = 3  # number of decimal points for the arc weights tolerance (swap_up)
C_MAX = 50  # max edge cost (for random generation, non-inclusive)

# special node IDs (see the node class docstring)
NROOT = 0
NTRUE = -1
NFALSE = -2
######################################################################


class node(object):
    """Encodes a node of the BDD.

    Attributes:
        id: an int node ID
        hi: a pointer to the hi-node
        lo: a pointer to the lo-node
        layer: layer number (id)

    Note:
    There are special node IDs:
        0   for the root
        -1  for True sink
        -2  for the False sink
    """

    def __init__(self, id, hi=None, lo=None, layer=-1):
        self.id = id
        self.hi = hi
        self.lo = lo
        self.layer = layer

    def link(self, to_whom, arc="hi"):
        """Helper function: links the node to another one."""
        if arc == "hi":
            self.hi = to_whom
        elif arc == "lo":
            self.lo = to_whom
        else:
            print("ERROR: wrong arc specifier while creating an edge -- {} (expected 'hi' or 'lo')".format(arc))

######################################################################


class BDD(object):
    """Encodes a BDD and implements layer-ordering (revision) methods.

    Attributes:
        layers (list): a list of layers (sets/hashes), including {T, F}-nodes
        nodes (dict): a set/hash of nodes, key by node ID
        vars (list): variables associated with layers
        T,F (node): pointers to `True` and `False` sink nodes
        weighted (bool): a flag whether arcs have weights
        weights (dict): arc weights, keys: `(id_from, id_to, arc_type)`
                        (with the latter either 'hi', or 'lo')
    """

    def __init__(self, N=2, vars=None, weighted=False):
        """Defines an empty BDD with N variables."""
        if not vars:
            vars = [i for i in range(1, N+1)]

        self.vars = vars
        self.weighted = weighted
        if weighted:
            self.weights = dict()
        self.layers = [set() for i in range(N)]
        self.T = node(NTRUE)
        self.F = node(NFALSE)
        self.layers.append(set([self.T, self.F]))
        self.max_id = -1
        self.av_node_ids = deque()
        self.nodes = dict()
        self.var_pos = dict(zip(vars, [i for i in range(len(self.layers))]))
        self.nodes.update({NTRUE: self.T, NFALSE: self.F})

    # some helper functions

    def link(self, parent, child, etype="hi", edge_weight=0.0):
        """Creates a link between two nodes.

        (Saving auxiliary info as appropriate)

        Args:
           self (type): description
           parent (int): parent node id
           child (int): child node id
           etype (str): edge type, "hi" or "lo" (default "hi")
           edge_weight (float): edge weight to e imposed (default 0.0)
        """
        assert parent in self.nodes.keys()
        assert child in self.nodes.keys()
        self.nodes[parent].link(self.nodes[child], etype)
        if self.weighted:
            self.weights[(parent, child, etype)] = edge_weight

    def llink(self, parent, child, etype="hi", edge_weight=0.0):
        """Creates a link between two nodes.

        (Saving auxiliary info as appropriate)

        Args:
           self (type): description
           parent (class node): parent node
           child (class node): child node
           etype (str): edge type, "hi" or "lo" (default "hi")
           edge_weight (float): edge weight to e imposed (default 0.0)
        """
        parent.link(child, etype)
        if self.weighted:
            self.weights[(parent.id, child.id, etype)] = edge_weight

    def __len__(self):
        """Returns no. of variables (layers, except the terminal)."""
        return len(self.vars)

    def n(self, i):
        """Returns size of the i-th layer."""
        return len(self.layers[i])

    def p(self, a):
        """Returns the position (layer index) of the variable `a`."""
        return self.var_pos[a]

    def size(self):
        """Returns the size of the BDD (total no. of nodes)."""
        return len(self.nodes)-2

    def new_node_name(self):
        """Returns a unique name for the next node to be created.

        Picks the one to re-cycle after some node deletions, when available;
        when it is not -- just takes the `max_id+1` (adjusting the `max_id` itself)
        """
        if len(self.av_node_ids) > 0:
            return self.av_node_ids.popleft()
        else:
            self.max_id += 1
            return self.max_id

    def dump_gv(self, layerCapt=True, x_prefix="x", node_labels=None):
        """Exports the BDD to the Graphviz format (`.dot`).
        Args:
            layerCapt (bool): whether to generate layer captions.
            x_prefix (str): a prefix to be shown for layer name (default 'x')
        Returns: a Digraph object (see `graphviz` module).
        """
        g = Digraph()

        for i, layer in enumerate(self.layers):
            with g.subgraph(name="cluster_{}".format(i)) as s:
                for n in layer:
                    if n.id == NTRUE:
                        s.node(str(NTRUE), label="T",fillcolor = "orange",style="filled")
                    elif n.id == NFALSE:
                        s.node(str(NFALSE),label="F", fillcolor = "gray",style="filled")
                    else:
                        if node_labels is None:
                            s.node( str(n.id) )
                        else:
                            s.node(str(n.id), label=f"{n.id}: {node_labels[n.id]}")

                if i!=len(self.layers)-1:
                    if layerCapt:
                        s.attr(label=f"{x_prefix}{self.vars[i]}, sz={len(layer)}", color="lightgrey")
                    else:
                        s.attr(color="lightgrey")
                else:
                    s.attr(color="white")

        for layer in self.layers:
            for n in layer:
                if not n.hi is None:
                    if self.weighted:
                        elabel = f"{self.weights[(n.id, n.hi.id, 'hi')]}"
                    else:
                        elabel = ""
                    g.edge(str(n.id), str(n.hi.id), label=elabel)

                if not n.lo is None:
                    if self.weighted:
                        elabel = f"{self.weights[(n.id, n.lo.id, 'lo')]}"
                    else:
                        elabel = ""
                    g.edge(str(n.id), str(n.lo.id), style="dashed", label=elabel)

        return g

    def addnode(self, parent_node, arc="hi", node_id=None, edge_weight=0.0):
        """Adds a node and updates aux info as necessary."""
        if node_id is None:
            node_id = self.new_node_name()

        newnode = node(node_id)

        if parent_node is None:
            newnode.layer = 0
            self.layers[0].add(newnode)

        elif arc == "hi":
            parent_node.hi = newnode
            self.layers[parent_node.layer+1].add(newnode)
            newnode.layer = parent_node.layer+1
            if self.weighted:
                self.weights[(parent_node.id,
                              parent_node.hi.id, 'hi')] = edge_weight

        elif arc == "lo":
            parent_node.lo = newnode
            self.layers[parent_node.layer+1].add(newnode)
            newnode.layer = parent_node.layer+1
            if self.weighted:
                self.weights[(parent_node.id,
                              parent_node.lo.id, 'lo')] = edge_weight

        else:
            print("ERROR: Wrong arc type: {}\n(allowed values are 'hi' and 'lo')".format(arc))

        self.nodes.update({newnode.id: newnode})

        return newnode

    def swap_up(self, layer_idx):
        """Swaps the layer with the one immediately above it (in-place).

        (by _index_! in-place operation)

        Args:
            layer_idx (int): number (index) of the layer to be ''bubbled up''

        Note:
            Operation takes O (no-of-nodes in the upper layer) time.
            Drops the layer being swapped. Then iterates through the nodes
            of the layer immediately above it and creates nodes as necessary
            (re-cycling new nodes to avoid redundancy in terms of logical
            functions)
        """
        assert layer_idx >= 1  # otherwise there's no layer to swap to
        assert layer_idx <= len(self.layers)-2  # can't swap-up the last layer

        new_nodes = dict()

        for n in self.layers[layer_idx]:
            del self.nodes[n.id]
            self.av_node_ids.append(n.id)

        self.layers[layer_idx] = set()

        prev_layer = self.layers[layer_idx-1]

        if not self.weighted:
            for F in prev_layer:
                # iterating through the nodes of the upper layer

                # save the pointers
                F_hi_hi = F.hi.hi
                F_hi_lo = F.hi.lo
                F_lo_hi = F.lo.hi
                F_lo_lo = F.lo.lo

                # creating the "hi"-node
                if not (F_hi_hi.id, F_lo_hi.id) in new_nodes.keys():
                    F_hi = self.addnode(F, "hi")
                    F_hi.link(F_hi_hi, "hi")
                    F_hi.link(F_lo_hi, "lo")
                    new_nodes.update({(F_hi.hi.id, F_hi.lo.id): F_hi})
                else:
                    F_hi = new_nodes[(F_hi_hi.id,
                                    F_lo_hi.id)]  # re-cycle the old node
                    F.link(F_hi, "hi")

                # creating the "lo"-node
                if not (F_hi_lo.id, F_lo_lo.id) in new_nodes.keys():
                    F_lo = self.addnode(F, "lo")
                    F_lo.link(F_hi_lo, "hi")
                    F_lo.link(F_lo_lo, "lo")
                    new_nodes.update({(F_lo.hi.id, F_lo.lo.id): F_lo})
                else:
                    F_lo = new_nodes[(F_hi_lo.id,
                                    F_lo_lo.id)]  # re-cycle the old node
                    F.link(F_lo, "lo")

            self.layers[layer_idx] = set(new_nodes.values())

        else:
            # weighted version (with arc costs)
            old_weights = copy.copy(self.weights)
            for F in prev_layer:
                # save the pointers
                F_hi_hi = F.hi.hi
                F_hi_lo = F.hi.lo
                F_lo_hi = F.lo.hi
                F_lo_lo = F.lo.lo

                # calculating the necessary costs
                cF_hi = old_weights[(F.id, F.hi.id, "hi")]
                cF_lo = old_weights[(F.id, F.lo.id, "lo")]
                cF_hi_hi = old_weights[(F.hi.id, F.hi.hi.id, "hi")]
                cF_lo_hi = old_weights[(F.lo.id, F.lo.hi.id, "hi")]
                cF_hi_lo = old_weights[(F.hi.id, F.hi.lo.id, "lo")]
                cF_lo_lo = old_weights[(F.lo.id, F.lo.lo.id, "lo")]

                cH_min = np.round(min(cF_hi + cF_hi_hi, cF_lo + cF_lo_hi),
                                  decimals=WEIGHT_TOL)
                cH_hi = np.round((cF_hi + cF_hi_hi) - cH_min,
                                 decimals=WEIGHT_TOL)
                cH_lo = np.round((cF_lo + cF_lo_hi) - cH_min,
                                 decimals=WEIGHT_TOL)

                cL_min = np.round(min(cF_hi + cF_hi_lo, cF_lo + cF_lo_lo),
                                  decimals=WEIGHT_TOL)
                cL_hi = np.round((cF_hi + cF_hi_lo) - cL_min,
                                 decimals=WEIGHT_TOL)
                cL_lo = np.round((cF_lo + cF_lo_lo) - cL_min,
                                 decimals=WEIGHT_TOL)

                # ----------------------------------------------------
                # creating the "hi"-node

                if not (F_hi_hi.id, F_lo_hi.id,
                        cH_hi, cH_lo) in new_nodes.keys():
                    F_hi = self.addnode(F, "hi", edge_weight=cH_min)
                    self.llink(F_hi, F_hi_hi, "hi", edge_weight=cH_hi)
                    self.llink(F_hi, F_lo_hi, "lo", edge_weight=cH_lo)
                    new_nodes.update({(F_hi.hi.id, F_hi.lo.id,
                                       cH_hi, cH_lo): F_hi})
                else:
                    F_hi = new_nodes[(F_hi_hi.id, F_lo_hi.id,
                                      cH_hi, cH_lo)]  # re-cycle the old node
                    self.llink(F, F_hi, "hi", edge_weight=cH_min)

                # ----------------------------------------------------
                # creating the "lo"-node

                if (F_hi_lo.id, F_lo_lo.id,
                        cL_hi, cL_lo) not in new_nodes.keys():
                    F_lo = self.addnode(F, "lo", edge_weight=cL_min)
                    self.llink(F_lo, F_hi_lo, "hi", edge_weight=cL_hi)
                    self.llink(F_lo, F_lo_lo, "lo", edge_weight=cL_lo)

                    new_nodes.update({(F_hi_lo.id, F_lo_lo.id,
                                       cL_hi, cL_lo): F_lo})
                else:
                    # re-cycle the old node
                    F_lo = new_nodes[(F_hi_lo.id, F_lo_lo.id, cL_hi, cL_lo)]
                    self.llink(F, F_lo, "lo", edge_weight=cL_min)

            self.layers[layer_idx] = set(new_nodes.values())
        # NOTE: old nodes in the source layer must be (physically) deleted
        # by the Python garbage collection

        # rename layers and update aux structures as necessary
        self.var_pos.update({
            self.vars[layer_idx]: layer_idx - 1,
            self.vars[layer_idx - 1]: layer_idx
        })

        v_1 = self.vars[layer_idx-1]
        self.vars[layer_idx-1] = self.vars[layer_idx]
        self.vars[layer_idx] = v_1

    def sift(self, var, pos):
        """Sifts variable `var` (name) to position `pos`, in-place."""
        assert pos >= 0 and pos < len(self.layers)-1
        assert var in self.vars

        if pos < self.p(var):
            # this is a sift-up
            while self.p(var) > pos:
                self.swap_up(self.p(var))

        elif pos == self.p(var):
            # a trivial case -- nothing to do
            return

        else:
            # this is a sift-down
            while pos > self.p(var):
                self.swap_up(self.p(var)+1)

    def gsifts(self, with_whom, start_order=None):
        """Perform greedy sifts to align with `with_whom`.

        Runs a simplistic implementation of Rudell'93
        sifting alorithm extension to minimize |A|+|B|

        starts with aligning to self or start_order (if given)
        """
        N = len(self.layers)-1
        if start_order is None:
            start_order = self.vars
        else:
            for i in range(N):
                self.sift(start_order[i], i)

        for i in range(N):
            with_whom.sift(start_order[i], i)

        # now the pair is aligned, but possibly huge
        # let us see if we can compress it:

        processed_vars = set()

        while len(processed_vars) < N:
            i = 0
            while self.vars[i] in processed_vars:
                i += 1  # skip processed variables

            best_pos = i
            active_var = self.vars[i]
            best_size = self.size() + with_whom.size()

            cur_size = best_size

            # try moving the var down
            for j in range(i+1, N):
                #if self.size()+with_whom.size() > cur_size*GSIFTS_MAX_INCREASE_MUL:
                #    break

                self.swap_up(j)
                with_whom.swap_up(j)
                cur_size = self.size()+with_whom.size()

                if cur_size < best_size:
                    best_size = cur_size
                    best_pos = j

            # try moving the var up
            self.sift(active_var, i)
            with_whom.sift(active_var, i)
            cur_size = self.size()+with_whom.size()

            for j in reversed(range(1,i+1)):
                #if self.size()+with_whom.size() > cur_size*GSIFTS_MAX_INCREASE_MUL:
                #    break

                self.swap_up(j)
                with_whom.swap_up(j)
                cur_size = self.size()+with_whom.size()

                if cur_size < best_size:
                    best_size = cur_size
                    best_pos = j-1

            # now choose the best position for the variable
            self.sift(active_var, best_pos)
            with_whom.sift(active_var, best_pos)

            processed_vars.add(active_var)

    # file manipulation procedures
    def save(self, filename):
        """Saves a BDD to an ASCII(text) file.

        Note:
            File format is as follows.
            File header:
                N=<no-of-layers>
                vars=<var1, var2, var3, ...> (variables' names)
            Then, one node description = one line of the format:
                id:hi,lo

            where `id` is node's ID, `hi` and `lo` are IDs(!) of the
            nodes corresponding to hi and lo pointers of the node "ID"
            The procedure performs breadth-first search and just saves
            all the nodes.
        """
        S = deque()  # a deque of node IDs to explore (FIFO)
        V = set()  # set of visited (saved) nodes (by ID)

        for n in self.layers[0]:
            S.append(n.id)

        with open(filename, "w") as f:
            f.write("N={}\n".format(len(self.layers)))
            f.write("vars={}\n".format(', '.join([str(v) for v in self.vars])))

            while len(S) > 0:
                n_id = S.popleft()
                # this is a new node
                if not (self.nodes[n_id].hi is None
                        or self.nodes[n_id].hi is None):
                    n_hi = self.nodes[n_id].hi.id
                    n_lo = self.nodes[n_id].lo.id
                    f.write("{}:{},{}\n".format(n_id, n_hi, n_lo))
                    if not (n_hi in V or n_hi in S or n_hi is None):
                        S.append(n_hi)

                    if not (n_lo in V or n_lo in S or n_lo is None):
                        S.append(n_lo)

                else:
                    f.write("{}:None,None\n".format(n_id))

                V.add(n_id)

    def load(self, filename, weighted=False):
        """Loads a BDD from an ASCII(text) file.

        Note:
            The format corresponds to the one described in the `save` function
        """
        with open(filename, "r") as f:
            line = f.readline()
            while line[0] == "#" or line == "":
                line = f.readline()

            assert line[:2] == 'N='
            N = int(line[2:])
            assert N > 1

            line = f.readline()
            assert line[:5] == 'vars='
            line = line[5:].split(',')
            assert len(line) == N-1

            # initialize the attributes
            self.layers = [set() for i in range(N)]
            self.nodes = dict()
            self.vars = [v.strip() for v in line]
            self.var_pos = dict(zip(self.vars, [i for i in range(N)]))

            self.weighted = weighted
            if weighted:
                self.weights = dict()

            line = f.readline()

            current_layer = 0
            next_layer = set()
            while line:
                if weighted:
                    id, line, weights = line.split(':')
                    weights = weights.split(",")
                    w_hi = float(weights[0])
                    w_lo = float(weights[1])
                else:
                    id, line = line.split(':')

                hi_id, lo_id = line.split(',')
                id = int(id)
                if id != NTRUE and id != NFALSE:
                    hi_id = int(hi_id)
                    lo_id = int(lo_id)
                    if weighted:
                        self.weights[(id, hi_id, "hi")] = w_hi
                        self.weights[(id, lo_id, "lo")] = w_lo
                else:
                    line = f.readline()
                    continue

                if id in next_layer:
                    current_layer += 1
                    next_layer = set()
                elif id not in self.nodes.keys():
                    self.addnode(None, node_id=id)
                    if id != 0:
                        print(f"WARNING: during load, just added an orphan node with id: {id}")

                if hi_id in next_layer:
                    F_hi = self.nodes[hi_id]
                    self.nodes[id].link(F_hi, "hi")
                else:
                    if self.weighted:
                        e_w = w_hi
                    else:
                        e_w = 0.0
                    F_hi = self.addnode(self.nodes[id],
                                        "hi",
                                        node_id=hi_id,
                                        edge_weight=e_w)
                    next_layer.add(hi_id)

                if lo_id in next_layer:
                    F_lo = self.nodes[lo_id]
                    self.nodes[id].link(F_lo, "lo")
                else:
                    if self.weighted:
                        e_w = w_lo
                    else:
                        e_w = 0.0
                    F_lo = self.addnode(self.nodes[id],
                                        "lo",
                                        node_id=lo_id,
                                        edge_weight=e_w)
                    next_layer.add(lo_id)

                if id not in self.nodes.keys():
                    if current_layer > 0:
                        print("WARNING: a node with no source at layer {}".
                              format(current_layer))

                    F = node(id, F_hi, F_lo)
                    self.layers[current_layer].add(F)
                    self.nodes.update({id: F})

                if id == NTRUE:
                    self.T = self.nodes[id]

                if id == NFALSE:
                    self.F = self.nodes[id]

                line = f.readline()

                self.max_id = max([n.id for n in self.nodes.values()])
                self.av_node_ids = deque([
                    i for i in range(1, self.max_id + 1)
                    if i not in self.nodes.keys()
                ])

    @classmethod
    def random(cls, N=5, p=0.5, weighted=False):
        """Generates a random BDD with `N` variables (N+1 layers).

        Args:
            N (int):   number of variables (results in N+1 layers, incl. T,F)
            p (float):   tree size parameter
                    0 will generate a non-random exponential-sized diagram,
                    1 will result in a single node per layer.
            weighted (bool): whether to generate a weighted diagram.
        Returns:
            A generated BDD
        """
        bdd = BDD(N, weighted=weighted)

        assert N > 1
        bdd.addnode(parent_node=None)  # create a root node

        for layer in range(1, N):
            # generate nodes (and edges) for the layer

            for n in bdd.layers[layer-1]:
                ch = np.random.randint(0, 1+C_MAX)
                cl = np.random.randint(0, 1+C_MAX)
                if rnd() <= p or len(bdd.layers[layer]) == 0:
                    newnode = bdd.addnode(n, arc="hi", edge_weight=ch)
                else:
                    bdd.link(n.id, rnd_choose(tuple(bdd.layers[layer])).id,
                             "hi", edge_weight=ch)

                if rnd() <= p:
                    newnode = bdd.addnode(n, arc="lo", edge_weight=cl)
                else:
                    bdd.link(n.id, rnd_choose(tuple(bdd.layers[layer])).id,
                           "lo", edge_weight=cl)

        # separately process the last layer
        for n in bdd.layers[-2]:
            bdd.link(n.id, rnd_choose(tuple(bdd.layers[-1])).id,
                     "hi", edge_weight=np.random.randint(0, 1+C_MAX))
            bdd.link(n.id, rnd_choose(tuple(bdd.layers[-1])).id,
                     "lo", edge_weight=np.random.randint(0, 1+C_MAX))

        return bdd

    def profile(self):
        """Returns a BDD ``profile'' -- an (almost) BFS-ordered list of nodes.

        A string of the format <vars>:<BFS-nodes>, where:
        <vars>      --  comma-separated variable names by layer;
        <BFS-nodes> --  comma-separated list of node-numbers (not IDs!)
                        in a BFS-traverse order
        """
        Q = deque()
        V = set()
        node_nums = dict()
        cur_node = 0

        p = []

        for n in self.layers[0]:
            node_nums.update({n.id: cur_node})
            Q.append(n)
            p.append(str(cur_node))
            cur_node += 1

        p = [",".join(p) + ";"]

        while(len(Q) > 0):
            n = Q.popleft()

            for successor in [n.hi, n.lo]:
                if successor is None:
                    continue

                if successor.id not in node_nums.keys():
                    node_nums.update({successor.id: cur_node})
                    cur_node += 1

                p.append(str(node_nums[successor.id]))
                if successor.id not in V:
                    Q.append(successor)
                    V.add(successor.id)

        return p[0] + ",".join(p[1:])

    def show(self,
             dir="testing",
             filename="showfunc.dot",
             layerCapt=True,
             x_prefix="x"):
        """Shows the diagram.

        Generates a `.dot` file and compiles it to `.pdf`.

        Args:
            dir (str): directory to place files (default: "testing")
            filename (str): `.dot` filename (default: "showfunc.dot")
            layerCapt (bool): whether to show layer captions (default: True)
            x_prefix(str): a prefix for variable captions (default:'x')
        """
        self.dump_gv(layerCapt, x_prefix).view(filename,
                                               directory=dir,
                                               cleanup=True)

    def align_to(self, vars_order, inplace=False):
        """Revises the BDD to a given order.

        (with a series of consecutive sifts)
        """
        if inplace:
            aligned = self
        else:
            aligned = copy.deepcopy(self)

        for i in range(len(vars_order)):
            aligned.sift(vars_order[i], i)

        return aligned

    def OA_bruteforce(self, with_what, LR=False):
        """Finds the best (smallest) alignment with BDD `with_what`.

        Uses brute-force enumeration.

        Args:
            with_what (BDD): target diagram.
            LR (bool): if set, layer-reduces each candidate BDD.
        """
        # generate all possible permutations
        perms = iters.permutations(self.vars)
        min_size = math.inf
        A_aligned = []
        B_aligned = []
        alternatives = 0

        for perm in perms:
            A_c = self.align_to(list(perm))
            B_c = with_what.align_to(list(perm))

            if LR:
                A_c.make_reduced()
                B_c.make_reduced()

            if min_size > A_c.size() + B_c.size():
                A_aligned = [A_c]
                B_aligned = [B_c]
                min_size = A_c.size()+B_c.size()
                alternatives = 1
            elif min_size == A_c.size() + B_c.size():
                A_aligned.append(A_c)
                B_aligned.append(B_c)
                alternatives += 1

        return [alternatives, A_aligned, B_aligned]

    def is_reduced(self):
        """Checks if the BDD is reduced (no ``redundant'' nodes)."""
        checked_nodes = set()

        for layer in range(len(self.layers)-2, -1, -1):
            for n in self.layers[layer]:
                if (n.hi.id, n.lo.id) in checked_nodes:
                    return False
                else:
                    checked_nodes.add((n.hi.id, n.lo.id))

        return True

    def make_reduced(self):
        """Makes the BDD reduced.

        Swaps each layer back and forth, starting from the last one.
        """
        for i in range(len(self.vars)-1, 0, -1):
            self.swap_up(i)
            self.swap_up(i)

    def rename_vars(self, ren_dict):
        """Helper function: renames variables.

        Note: expects the *full* dict of variables (of length `N`)

        Args:
          ren_dict (dict): a dict of labels in the form {before: after}
        """
        assert len(ren_dict) == len(self)
        new_vars = [ren_dict[v] for v in self.vars]
        self.vars = new_vars
        self.var_pos = dict(zip(self.vars, [i for i in range(len(self.vars))]))

    def is_aligned(self, to_what):
        """Helper function: checks if the BDD is aligned w/ to_what."""
        return np.array_equal(self.vars, to_what.vars)

    # functions related to testing equivalence and finding truth tables
    def get_value(self, x):
        """Finds the terminal node (T or F) -- an endpoint for `x`.

        Finds a terminal node corresponding to the variable choices in `x`.

        Args:
            x (dict): of {var_name: value}, where value is in {0,1}.

        Returns:
            terminal (bool or `-1`): terminal node corresponding to the path
                            implied by x, encoded as Boolean (or -1 if error)
            cost (float): if the BDD is weighted, cost of the path.
        """
        (node,) = self.layers[0]
        cost = 0.0

        for i in range(len(x)):
            if x[self.vars[i]] == 0:
                if self.weighted:
                    cost += self.weights[(node.id, node.lo.id, "lo")]
                node = node.lo
            elif x[self.vars[i]] == 1:
                if self.weighted:
                    cost += self.weights[(node.id, node.hi.id, "hi")]
                node = node.hi
            else:
                print(
                    "Wrong value ({}) for variabe {}: 0 or 1 expected".format(
                        x[self.vars[i]], self.vars[i]))
                return -1

        if node.id == self.T.id:
            if self.weighted:
                return [True, cost]
            else:
                return True
        elif node.id == self.F.id:
            if self.weighted:
                return [False, cost]
            else:
                return False
        else:
            print("Error not a True or False node reached!")
            print("node is {}, T is {}, F is {}".format(
                node.id, self.T.id, self.F.id))
            return -1

    def truth_table(self):
        """Returns a truth table for the BDD (as a Boolean function)."""
        tt = []
        ind = []
        for x_num in range(2**len(self)):
            x = [int(j) for j in np.binary_repr(x_num, width=len(self.vars))]
            if not self.weighted:
                tt.append(x + [self.get_value(dict(zip(self.vars, x)))])
            else:
                tt.append(x + self.get_value(dict(zip(self.vars, x))))
            ind.append(np.binary_repr(x_num, width=len(self.vars)))

        if not self.weighted:
            return pd.DataFrame(tt, columns=[str(x) for x in self.vars] + ["Terminal"], index=ind)
        else:
            return pd.DataFrame(tt,
                                columns=[str(x) for x in self.vars] + ["Terminal", "Cost"],
                                index=ind)

    def is_equivalent(self, B):
        """Checks if two BDDs (A and B) are equivalent.

        Equivalent in the sense that they define the same Boolean function,
        by checking:
        - if the corresponding truth tables coincide.
        - (if the BDD is weighted) if all paths have the same costs.

        Returns Bool.
        """
        assert self.weighted == B.weighted

        msg = None
        tt_self = self.truth_table()
        tt_B = B.truth_table()
        if self.weighted:
            tt_B = tt_B[[str(x) for x in self.vars] + ['Terminal', 'Cost']]
        else:
            tt_B = tt_B[[str(x) for x in self.vars] + ['Terminal']]

        if self.weighted:
            tt_B['new_index'] = [
                c for c in tt_B.apply(
                    lambda row: "".join([str(x) for x in row[:-2]]), axis=1)
            ]
        else:
            tt_B['new_index'] = [
                c for c in tt_B.apply(
                    lambda row: "".join([str(x) for x in row[:-1]]), axis=1)
            ]
        tt_B.set_index('new_index', inplace=True)

        for idx, row in tt_self.iterrows():
            if row['Terminal'] != tt_B.loc[idx]['Terminal']:
                msg = "\nINFO: Discrepancy found. For x={}, value(self)={}, value(B)={}".format(
                    idx, row['Terminal'], tt_B.loc[idx]['Terminal'])
                return [False, msg]
            if self.weighted:
                if abs(row['Cost'] - tt_B.loc[idx]['Cost']) >= WEIGHT_TOL:
                    msg = "\nINFO: Discrepancy found. For x={}, cost(self)={}, cost(B)={}".format(
                        idx, row['Cost'], tt_B.loc[idx]['Cost'])
                    return [False, msg]

        return [True, msg]


def intersect(A, B):
    """Produces an intersection BDD of `A` and `B`.

    Notes:
        BDDs need to be order-associated (same vars in the same order).
    """
    assert A.vars == B.vars
    assert A.weighted == B.weighted

    weighted = A.weighted
    C = BDD(vars=A.vars, N=len(A.vars), weighted=weighted)
    root = C.addnode(parent_node=None)
    current_layer = {(NROOT, NROOT): root}

    # for all layers but the last one
    for i in range(len(A.vars)-1):
        new_layer = dict()
        for n in current_layer:
            niA, niB = n
            nA = A.nodes[niA]
            nB = B.nodes[niB]
            nC = current_layer[n]
            # yes-arc
            new_yes = (nA.hi.id, nB.hi.id)
            wA = A.weights[(nA.id, nA.hi.id, "hi")]
            wB = B.weights[(nB.id, nB.hi.id, "hi")]
            if new_yes in new_layer:
                C.link(nC.id, new_layer[new_yes].id,
                       "hi", edge_weight=wA+wB)
            else:
                new_layer[new_yes] = C.addnode(nC,
                                               "hi", edge_weight=wA+wB)
            # no-arc
            new_no = (nA.lo.id, nB.lo.id)
            wA = A.weights[(nA.id, nA.lo.id, "lo")]
            wB = B.weights[(nB.id, nB.lo.id, "lo")]
            if new_no in new_layer:
                C.link(nC.id, new_layer[new_no].id,
                       "lo", edge_weight=wA+wB)
            else:
                new_layer[new_no] = C.addnode(nC,
                                              "lo", edge_weight=wA+wB)

        current_layer = new_layer

    # a special case of the last layer
    for n in current_layer:
        niA, niB = n
        nA = A.nodes[niA]
        nB = B.nodes[niB]
        nC = current_layer[n]
        # yes-arc
        wA = A.weights[(nA.id, nA.hi.id, "hi")]
        wB = B.weights[(nB.id, nB.hi.id, "hi")]
        if nA.hi.id == NTRUE and nB.hi.id == NTRUE:
            C.link(nC.id, NTRUE, "hi", edge_weight=wA+wB)
        else:
            C.link(nC.id, NFALSE, "hi", edge_weight=wA+wB)
        # no-arc
        wA = A.weights[(nA.id, nA.lo.id, "lo")]
        wB = B.weights[(nB.id, nB.lo.id, "lo")]
        if nA.lo.id == NTRUE and nB.lo.id == NTRUE:
            C.link(nC.id, NTRUE, "lo", edge_weight=wA+wB)
        else:
            C.link(nC.id, NFALSE, "lo", edge_weight=wA+wB)

    return C


# =====================================================================
# testing code
# =====================================================================
def test_intersections():
    """Tests intersection function."""
    A = BDD()
    A.load("./tests/int_A.wbdd", weighted=True)
    B = BDD()
    B.load("./tests/int_B.wbdd", weighted=True)

    C = intersect(A, B)
    A.show(dir="testing", filename="A.dot")
    B.show(dir="testing", filename="B.dot")
    C.show(dir="testing", filename="C.dot")


# BDD creation and rendering (with graphviz)
def test_create_render():
    """Tests the BDD creation and rendering code."""
    bdd = BDD(4)
    root = bdd.addnode(None)
    nh = bdd.addnode(root, "hi")
    nl = bdd.addnode(root, "lo")

    nh2 = bdd.addnode(nh, "hi")
    bdd.nodes[nh.id].lo = nh2

    nl2 = bdd.addnode(nl, "hi")
    nl.link(nl2,"lo")

    nl = bdd.addnode(nl2, "hi")
    nl2.link(nl, "lo")

    nh = bdd.addnode(nh2,"lo")
    nh2.link(nh,"hi")

    nl.link(bdd.T,"lo")
    nl.link(bdd.F,"hi")

    nh.link(bdd.T,"lo")
    nh.link(bdd.F, "hi")

    g = bdd.dump_gv()

    g.render("./initial.dot",view=True)

def gen_4BDD():
    """creates a simple 4-var BDD"""

    bdd = BDD(4)
    root = bdd.addnode(None)
    nh = bdd.addnode(root,"hi")
    nl = bdd.addnode(root,"lo")

    nh2 = bdd.addnode(nh,"hi")
    bdd.nodes[nh.id].lo = nh2

    nl2 = bdd.addnode(nl, "hi")
    nl.link(nl2,"lo")

    nl = bdd.addnode(nl2,"hi")
    nl2.link(nl,"lo")

    nh = bdd.addnode(nh2,"lo")
    nh2.link(nh,"hi")

    nl.link(bdd.T,"lo")
    nl.link(bdd.F,"hi")

    nh.link(bdd.T,"lo")
    nh.link(bdd.F, "hi")
    return bdd

def test_swap_sift():
    """quick test of swap and sift operations"""
    ## swap operation
    bdd.swap_up(2)

    g = bdd.dump_gv()
    g.render("./swapped.dot", view=True)

    bdd.swap_up(1)

    g = bdd.dump_gv()
    g.render("./swapped_top.dot", view=True)

    ## sift operatio
    bdd.sift(3,3)
    g = bdd.dump_gv()
    g.render("./sifted.dot", view=True)

    ### more complex sifting example

    bdd = BDD(4)
    root = bdd.addnode(None)
    n1 = bdd.addnode(root, "hi")
    n2 = bdd.addnode(root, "lo")

    n3 = bdd.addnode(n1,"hi")
    n4 = bdd.addnode(n1,"lo")
    n5 = bdd.addnode(n2,"hi")
    n2.link(n4,"lo")

    n6 = bdd.addnode(n3,"hi")
    n3.link(n6,"lo")

    n7 = bdd.addnode(n4,"hi")
    n4.link(n7,"lo")

    n8 = bdd.addnode(n5, "lo")
    n5.link(n7,"hi")

    n6.link(bdd.F,"hi")
    n6.link(bdd.F,"lo")

    n7.link(bdd.F,"hi")
    n7.link(bdd.F,"lo")

    n8.link(bdd.F,"lo")
    n8.link(bdd.T,"hi")

    # g = bdd.dump_gv()
    # g.view()

    bdd.sift(4,1)

    g = bdd.dump_gv()
    g.render("./after_sifting.dot",view=True)


def test_save_load(bdd):
    """quickly tests the load-save functionality

    must result in two equivalent graphs (BDDs) on the screen
    """
    bdd.save("./test_save.bdd")
    bdd2 = BDD()
    bdd2.load("./test_save.bdd")
    bdd.dump_gv().view("./before_save.dot",cleanup=True)
    bdd2.dump_gv().view("./after_save.dot",cleanup=True)

def test_save_load_noargs():
    bdd = gen_4BDD()
    bdd.sift(3,0)
    test_save_load(bdd)


def test_rnd():
    bdd = BDD.random(N=5,p=0.8)
    bdd.dump_gv().view("./random.dot",directory="./testing",cleanup = True)

def test_align():
    bdd_A = BDD.random(N=4,p=0.8)
    bdd_B = BDD.random(N=4,p=0.8)

    alts, Aa, Ba = bdd_A.OA_bruteforce(bdd_B)
    print("Opt order: {}, opt size:{}".format(Aa[0].vars, Aa[0].size()+Ba[0].size()))

    bdd_A.dump_gv().render("./testing/A.dot",view=True)
    bdd_B.dump_gv().render("./testing/B.dot",view=True)
    Aa[0].dump_gv().render("./testing/Aa.dot",view=True)
    Ba[0].dump_gv().render("./testing/Ba.dot", view=True)


def test_bruteforcing():
    bdd = BDD.random(N=4, p=0.8) # it works for gen_4BDD(), though

    perms = iters.permutations(bdd.vars)

    for p in perms:
        print("Aligning to: {}".format(p))
        bdd = bdd.align_to(list(p))
        bdd.show()
        input("Press Enter to continue...")

def test_rnd_naming():
    mybdd = BDD.random(N=4,p=0.8)
    mybdd.show()

    mybdd_al = copy.deepcopy(mybdd)
    mybdd_al.show()

    mybdd_al.swap_up(2)
    # mybdd_al = mybdd.align_to([2,1,3,4])
    mybdd_al.show()
    mybdd_al.swap_up(mybdd_al.p(3))
    mybdd_al.show()

def test_swapping_2():
    mybdd = BDD.random(N=4,p=0.8)
    mybdd.show()

    bdds = copy.deepcopy(mybdd)
    bdds.swap_up(3)
    bdds.dump_gv().render("./testing/after_swap.dot",view=True, cleanup=True)


    mybdd = BDD.random(N=4,p=0.8)
    mybdd.show()

    bdds = copy.deepcopy(mybdd)
    bdds.swap_up(3)
    bdds.dump_gv().render("./testing/after_swap.dot", view=True, cleanup=True)


def test_random_swapping(N, k, m, p=0.8, mode="swaps"):
    """TEST: checks if swaps work correctly (random problems).

    Generates a random BDD, makes some random swaps (or sifts),
    and then checks that the function encoded by the BDD remained the same
    (that is, every set of var values result in the same terminal node for
    the original BDD and for the randomly changed BDD).

    Involves brute-force enumaration of all the 2^n possible decisions
    (concerning all vars).

    Args:
        N: no. of variables in the diagram
        k: no. of random BDDs generated)
        m: no. of consecutive test swaps per BDD
        p: BDD expansion prob parameter
        mode: 'swaps' or 'sifts', determines what kind of events
                are to be generated

    Example:
        test_random_swapping(8,500,20,0.8,sys.argv[1])
    """
    status = "OK"

    for n in range(k):
        bdd = BDD.random(N, p)
        bdd_s = copy.deepcopy(bdd)

        # make some random layer swaps
        if mode == 'swaps':
            for t in range(m):
                bdd_s.swap_up(randint(1, len(bdd.vars) - 1))
        elif mode == 'sifts':
            for t in range(m):
                bdd_s.sift(bdd_s.vars[randint(0,
                                              len(bdd.vars) - 1)],
                           randint(0,
                                   len(bdd.vars) - 1))
        else:
            print(
                "ERROR: wrong mode -- {}. Only 'swaps' and 'sifts' are allowed"
                .format(mode))

        # now check that everything remained the same
        no_of_bits = len(bdd.vars)

        for bit_path in range(0, 2**no_of_bits):

            # unpack the var values in the original order (1=hi, 0=lo)
            var_vals = [(bit_path >> bit) & 1
                        for bit in range(no_of_bits - 1, -1, -1)]

            # get the correct answer (T or F)
            cur_node_src = list(bdd.layers[0])[0]
            cur_node_swd = list(bdd_s.layers[0])[0]

            i = 0
            while i < len(bdd.vars):
                if var_vals[i] == 0:
                    cur_node_src = cur_node_src.lo
                else:
                    cur_node_src = cur_node_src.hi

                if var_vals[bdd.p(bdd_s.vars[i])] == 0:
                    cur_node_swd = cur_node_swd.lo
                else:
                    cur_node_swd = cur_node_swd.hi

                i += 1

            corr_answer = cur_node_src.id
            assert corr_answer == NTRUE or corr_answer == NFALSE
            if cur_node_swd.id != corr_answer:
                status = "\nERROR: {}, not {}".format(cur_node_swd.id,
                                                      cur_node_src)
                print("path {}: {}".format(var_vals, status))

        if status != "OK":
            print("ERROR encountered!")
        else:
            print(".", end="")
    print("\n{} instances processed, status: {}".format(k, status))


def test_swaps_w():
    """Tests if swaps are implemented correctly for a weighted BDD.

    (All paths costs and terminals must coincide after a series of swaps.)

    Returns:
        Nothing (prints on the screen)
    """
    A = BDD()
    A.load("./tests/simple_DD.wbdd", weighted=True)

    B = BDD()
    B.load("./tests/simple_DD.wbdd", weighted=True)
    eq, msg = A.is_equivalent(B)
    assert eq, msg

    B.swap_up(1)
    eq, msg = A.is_equivalent(B)
    assert eq, msg

    B.swap_up(1)
    eq, msg = A.is_equivalent(B)
    assert eq, msg

    B.swap_up(1)
    B.swap_up(2)
    eq, msg = A.is_equivalent(B)
    assert eq, msg

def test_swaps_weighted():
    for i in range(1000):
        A = BDD.random(N=5, weighted=True)
        B = copy.deepcopy(A)

        B.swap_up(1)
        eq, msg = A.is_equivalent(B)
        assert eq, msg

        B.swap_up(len(A)-1)
        eq, msg = A.is_equivalent(B)
        assert eq, msg

        for k in range(5):
            C = copy.deepcopy(B)
            what = np.random.randint(2, len(A))
            B.swap_up(what)
            eq, msg = A.is_equivalent(B)
            assert eq, f"for swap at {pos}: {msg}"

        print(".", end="")
        sys.stdout.flush()

def test_swaps_uweighted():
    for i in range(100):
        A = BDD.random(N=10, weighted=False)
        B = copy.deepcopy(A)

        B.swap_up(1)
        eq, msg = A.is_equivalent(B)
        assert eq, msg

        B.swap_up(len(A)-1)
        eq, msg = A.is_equivalent(B)
        assert eq, msg

        for k in range(5):
            C = copy.deepcopy(B)
            what = np.random.randint(2, len(A))
            print(f"swapping up {what}")
            B.swap_up(what)
            eq, msg = A.is_equivalent(B)
            assert eq, msg

        print(".", end="")
        sys.stdout.flush()

def main():
    test_swaps_weighted()

if __name__ == '__main__':
    main()
