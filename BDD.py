"""Implements the Binary Decision Diagram data type.

Defines the "original" BDD-related machinery (as opposed to the
simplification defined in :py:class:`varseq.VarSeq`), including
actual simultaneous variable reordering (inspired by Rudell'93) and
other BDD alignment-related functions, along with technical utilities
(such as load/save to file, etc.). A simple auxiliary class
:py:class:`node` keeps the necessary data related to BDD nodes.

Tests coverage: :py:mod:`BDD_test`

"""

import numpy as np
import pandas as pd
from numpy.random import random as rnd
from random import choice as rnd_choose
import sys
import pytest

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

def simscore(seq_A, seq_B, p_B = None):
    """Calculates a 'similarity score' between `A` and `B`.

    The function has nothing to do with other data structures,
    technically: it just calculates the normalized number of
    inversions between the two ordered lists of variables
    (which we call 'similarity score').

    Args:
        seq_A (list): sequence A
        seq_B (list): sequence B to compare
        p_B (dict): dict of element positions in B (default: None, will be built)

    Returns:
        number (between 0 and 1): share of possible inversions.

    Examples:
        >>> simscore([1,2,3], [1,2,3])
        1.0

        >>> simscore([1,2,3], [3,2,1])
        0.0

        >>> simscore([1,2,3], [1,3,2])
        0.6666666666666667

        >>> simscore([1,2,3,4,5,6,7], [5,6,3,4,2,7,1])
        0.33333333333333337
    """
    N = len(seq_A)
    no_inversions = 0
    if p_B is None:
        p_B = {seq_B[i]: i for i in range(N)}

    assert set(seq_A) == set(seq_B)
    for i in range(N-1):
        for j in range(i, N):
            if p_B[seq_A[i]] > p_B[seq_A[j]]:
                no_inversions += 1

    return 1 - no_inversions / (N*(N-1)/2)


class node(object):
    """Represents a BDD node, with the following attributes.

    Attributes:
        id: an integer node ID
        hi: a pointer to the ``hi``-node
        lo: a pointer to the ``lo``-node
        layer: layer number (id)

    Note:
        There are special node IDs:
        ``0``   for the root,
        ``-1``  for True sink,
        ``-2``  for the False sink.
    """

    def __init__(self, id, hi=None, lo=None, layer=-1):
        """A simple class constructor."""
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
        N (int): number of variables (equals number of layers - 1)
        layers (list): a list of layers (sets/hashes), including {T, F}-nodes
        nodes (dict): a set/hash of nodes, key by node ID
        vars (list): variables associated with layers
        T,F (node): pointers to `True` and `False` sink nodes
        weighted (bool): a flag whether arcs have weights
        weights (dict): arc weights, keys: `(id_from, id_to, arc_type)`
                        (with the latter either 'hi', or 'lo')
    """

    def __init__(self, N=2, vars=None, weighted=False):
        """Defines an empty BDD with N variables.

        Args:
            N (int): number of variables ( = number of layers except the terminal),
            vars (list): variable names (strings or ints),
            weighted (bool): a flag for weighted BDD (assigning arc weights).
        """
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


    def width(self):
        """Returns BDD width (max no. nodes in a layer)."""
        return max([len(L) for L in self.layers])
    # some helper functions

    def link(self, parent, child, etype="hi", edge_weight=0.0):
        """Creates a link between two nodes in the BDD.

        (Saving auxiliary info as appropriate)

        Args:
           parent (int): parent node id
           child (int): child node id
           etype (str): edge type, ``hi`` or ``lo`` (default ``hi``)
           edge_weight (float): edge weight (default 0.0)
        """
        assert parent in self.nodes.keys()
        assert child in self.nodes.keys()
        self.nodes[parent].link(self.nodes[child], etype)
        if self.weighted:
            self.weights[(parent, child, etype)] = edge_weight

    def llink(self, parent, child, etype="hi", edge_weight=0.0):
        """Creates a link between two nodes.

        A version of :py:func:`BDD.BDD.link` that takes
        :py:class:`BDD.node` objects as arguments instead of node
        IDs.

        """
        parent.link(child, etype)
        if self.weighted:
            self.weights[(parent.id, child.id, etype)] = edge_weight

    def __len__(self):
        """Returns no of variables (layers, except the terminal)."""
        return len(self.vars)

    def n(self, i):
        """Returns size of the i-th layer."""
        return len(self.layers[i])

    def p(self, a):
        """Returns position (0-based layer index) of variable ``a``."""
        return self.var_pos[a]

    def size(self):
        """Returns the size of the BDD (total no. of nodes)."""
        return len(self.nodes)-2

    def new_node_name(self):
        """Returns a unique name for the next node to be created.

        Picks the one to re-cycle after some node deletions, when available;
        when it is not -- just takes the ``max_id+1`` (adjusting the ``max_id``)
        """
        if len(self.av_node_ids) > 0:
            return self.av_node_ids.popleft()
        else:
            self.max_id += 1
            return self.max_id

    def dump_gv(self, layerCapt=True, x_prefix="x", node_labels=None):
        """Exports the BDD to the Graphviz format (``.dot``).

        Args:
            layerCapt (bool): whether to generate layer captions.
            x_prefix (str): a prefix to be shown for layer name (default 'x')

        Returns:
            a Digraph object (see ``graphviz`` package).
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
                        elabel = f"{self.weights[(n.id, n.hi.id, 'hi')]:.2f}"
                    else:
                        elabel = ""
                    g.edge(str(n.id), str(n.hi.id), label=elabel)

                if not n.lo is None:
                    if self.weighted:
                        elabel = f"{self.weights[(n.id, n.lo.id, 'lo')]:.2f}"
                    else:
                        elabel = ""
                    g.edge(str(n.id), str(n.lo.id), style="dashed", label=elabel)

        return g

    def addnode(self, parent_node, arc="hi", node_id=None, edge_weight=0.0):
        """Adds a node and updates aux info as necessary.

        Args:
            parent_node (node): parent node pointer,
            arc (str): arc type, `hi` or `lo`,
            node_id: node name. If `None` is provided, the value is derived
                from :py:func:`BDD.BDD.new_node_name`,
            edge_weight(float): edge weight

        Returns:
            The created node as a :py:class:`BDD.node` object.
        """
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

        (by *index* ! in-place operation)

        Args:
            layer_idx (int): number (0-based index) of the layer to be ''bubbled up''

        Note:
            Operation takes O (`no-of-nodes in the upper layer`)
            time. Drops the layer being swapped. Then iterates
            through the nodes of the layer immediately above it and
            creates nodes as necessary (re-cycling new nodes to avoid
            redundancy: a quasi-reduced BDD keeps this property)

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
        """Sifts variable ``var`` (by name) to (0-based) position ``pos``.

        Uses :py:func:`BDD.sift` under the hood, in-place operation."""
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
        """Performs greedy sifts to align with ``with_whom``.

        Runs a simplistic implementation of Rudell'93
        sifting alorithm extension to minimize \|A\|+\|B\|

        starts with aligning to ``self`` or ``start_order`` (if given)
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
        """Saves a BDD to a text file.

        Args:
            filename (string): destination file name.

        Note:
            File format is as follows.

                - File header, the first two lines
                    | N= `<no-of-layers>`
                    | vars= `<comma-separated layer labels>`

                - Then, one node description = one line of the format:
                    `id` : `hi` , `lo`

            where ``id`` is node's ID, ``hi`` and ``lo`` are IDs of
            the nodes corresponding to one- and zero- pointers of the
            node ``ID``. The procedure performs breadth-first search
            and just saves all the nodes.

            Three nodes have reserved ID values:

                - **Root** node: 0
                - **True** terminal: -1
                - **False** terminal: -2

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
        """Loads a BDD from a (text) file.

        Note:
            The format is described in :py:func:`save`.
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

        Layers are numbered with consecutive integers by default,
        starting with ``1``.

        Args:
            N (int):   number of variables (results in N+1 layers, incl. T,F)
            p (float):   tree size parameter
                    0 will generate a non-random exponential-sized diagram,
                    1 will result in a single node per layer.
            weighted (bool): whether to generate a weighted diagram.
        Returns:
            A generated :py:class:`BDD`.

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
        """Returns a BDD "profile" -- an (almost) BFS-ordered list of nodes.

        Note:
            This is used to determine, e.g., if the generated
            instance is unique. A 'profile' is a string of the format
            ``<vars>``:``<BFS-nodes>``, where:

            - ``<vars>``      :  comma-separated variable names by layer;
            - ``<BFS-nodes>`` :  comma-separated list of node-numbers (not IDs!)
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
             dir="tmp",
             filename="showfunc.dot",
             layerCapt=True,
             x_prefix="x",
             node_labels=None):
        """Shows the diagram on the screen (saving the auxiliary pdf).

        Generates a ``.dot`` file and compiles it to ``.pdf``.
        Requires ``graphviz`` (``dot`` program)

        Args:
            dir (str): directory to place files (default: ``tmp``)
            filename (str): `.dot` filename (default: ``showfunc.dot``)
            layerCapt (bool): whether to show layer caption (default: ``True``)
            x_prefix(str): a prefix for variable captions (default: ``x``)
            node_labels(dict): node labels to display (optional).

        """
        self.dump_gv(layerCapt, x_prefix, node_labels).view(filename,
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
        """Finds the best (smallest) alignment with BDD ``with_what``.

        Uses brute-force enumeration of all possible orders.

        Args:
            with_what (BDD): target diagram.
            LR (bool): if set, layer-reduces each candidate BDD.

        Returns:
            A tuple of values

                - alternatives (int): number of alternative optima
                - A_aligned (list): of A-aligned (in each optimum)
                - B_aligned (list): of B-aligned (in respective optimum)
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
        """Checks if the BDD is quasi-reduced (no "redundant" nodes)."""
        checked_nodes = set()

        # quasi-reduced BDD for weighted diagrams must also check the
        # uniqueness of weights, not just uniqueness of (hi, lo) tuples.
        # Not a big issue, just a separate procedure (need to take care
        # of the precision, though).
        assert (not self.weighted), "ERROR: not implemented for weighted BDDs."

        for layer in range(len(self.layers)-2, -1, -1):
            for n in self.layers[layer]:
                if (n.hi.id, n.lo.id) in checked_nodes:
                    return False
                else:
                    checked_nodes.add((n.hi.id, n.lo.id))

        return True

    def make_reduced(self):
        """Makes the BDD quasi-reduced.

        Swaps each layer back and forth, starting from the last one,
        (which makes sure the results is quasi-reduced due to
        :py:func:`BDD.BDD.swap` implementation.)

        """
        for i in range(len(self.vars)-1, 0, -1):
            self.swap_up(i)
            self.swap_up(i)

    def rename_vars(self, ren_dict):
        """Helper function: renames variables.

        Args:
            ren_dict (dict): a dict of <= :py:attr:`BDD.N` labels
                in the form {``before``: ``after``}
        """
        for v in ren_dict.keys():
            assert v in self.vars, f"{v} is not in {self.vars}"

        new_vars = [(lambda v: ren_dict[v] if v in ren_dict.keys() else v)(v)
                    for v in self.vars]
        self.vars = new_vars
        self.var_pos = dict(zip(self.vars, [i for i in range(len(self.vars))]))

    def shuffle_vars(self):
        """Randomly renames the variables (shuffles their order)."""
        self.rename_vars(
            dict(zip(self.vars, np.random.permutation(self.vars))))

    def is_aligned(self, to_what):
        """Checks if the BDD is aligned with ``to_what``."""
        return np.array_equal(self.vars, to_what.vars)

    # functions related to testing equivalence and finding truth tables
    def get_value(self, x):
        """Finds the terminal node (T or F) -- an endpoint for ``x``.

        Returns a terminal node ID corresponding to the variable choices in ``x``.

        Args:
            x (dict): of {var_name: value}, where value is in {0,1}.

        Returns:
            A tuple of values

                - terminal (bool or `-1`):
                    terminal node corresponding to the path
                    implied by x, encoded as Boolean (or -1 if error)

                - cost (float):
                    if the BDD is weighted, cost of the path.
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
        """Returns a "truth table" for the encoded Boolean function.

        Returns:

        ``pandas`` dataframe: column numbers correspond to
        variable names with two additional columns: ``Terminal`` for
        the function value (of Boolean type) and ``Cost`` if the
        diagram is :py:attr:`BDD.weighted`.

        """
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

        Equivalent in the sense that they define the same Boolean
        function and the corresponding path costs coincide. Checks:

            - if the corresponding truth tables coincide.
            - (if the BDD is weighted) if all paths have the same costs.

        Returns:
            A tuple of

                - result (bool): whether the diagrams are equivalent,
                - message (str): information regarding the difference, if found.
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

    def shortest_path(self):
        """Finds a shortest-path in O(`no-of-nodes`) time.

        Returns:
            dict: node labels (using node IDs as keys)
            corresponding to a shortest path from the current node to
            **True** terminal.

        Note:
            Expects the diagram to be :py:attr:`BDD.weighted`. (An
            unweighted version, in terms of 'number of hops', does
            not make any sense for a quasi-reduced BDD: the answer
            would be ``len(self)``.)

        """
        assert self.weighted

        node_labels = {NTRUE: 0}

        for L in reversed(range(len(self.vars))):
            for node in self.layers[L]:
                node_label = None
                if node.hi.id in node_labels.keys():
                    node_label = node_labels[node.hi.id] + self.weights[
                        (node.id, node.hi.id, "hi")]

                if node.lo.id in node_labels.keys():
                    if node_label is not None:
                        node_label = min(
                            node_label, node_labels[node.lo.id] +
                            self.weights[(node.id, node.lo.id, "lo")])
                    else:
                        node_label = node_labels[node.lo.id] + self.weights[
                            (node.id, node.lo.id, "lo")]

                if node_label is not None:
                    node_labels.update({node.id: node_label})

        for node in self.nodes:
            if node not in node_labels.keys():
                node_labels[node] = "∞"

        return node_labels

    def simscore(self, B):
        """Computes a simscore with another diagram, `B` (see :py:func:`simscore`)."""
        return simscore(self.vars, B.vars, B.var_pos)

def intersect(A, B):
    """Produces an intersection BDD of ``A`` and ``B``.

    Notes:
        BDDs need to be order-associated (same vars in the same order).

    Returns:
        :py:class:`BDD` : an intersection BDD.
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
