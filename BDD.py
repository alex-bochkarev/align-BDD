"""
Implementation of the (exact) BDD-related machinery

Involves actual variable ordering and BDDs alignment classes
and functions (as opposed to the approximation of `varseq`)

(c) A. Bochkarev, Clemson University, 2019
"""

import numpy as np
from numpy.random import random as rnd
from random import choice as rnd_choose
from random import randint
import sys

import itertools as iters
from graphviz import Digraph
from collections import deque

import copy
import math

## special node IDs (see the node class docstring)
NROOT = 0
NTRUE = -1
NFALSE = -2

######################################################################
class node(object):
    """
    Encodes a node of the BDD

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

    def __init__(self, id, hi=None,lo=None, layer=-1):
        self.id = id
        self.hi = hi
        self.lo = lo
        self.layer = layer

    def link(self, to_whom, arc="hi"):
        """Helper function: links the node to another one"""
        if arc=="hi":
            self.hi = to_whom
        elif arc=="lo":
            self.lo = to_whom
        else:
            print("ERROR: wrong arc specifier while creating an edge -- {}".format(arc))

######################################################################
class BDD(object):
    """
    Encodes a BDD and implements layer-ordering methods

    Attributes:
        layers: a list of layers (sets/hashes)
        nodes: a set/hash of nodes
        vars: variables associated with layers
        T,F: pointers to `True` and `False` sink nodes
    """

    def __init__(self,N=2, vars = None):
        """defines an empty BDD with N variables"""

        if not vars:
            vars = [i for i in range(1,N+1)]

        self.vars = vars
        self.layers = [set() for i in range(N)]
        self.T = node(NTRUE)
        self.F = node(NFALSE)
        self.layers.append(set([self.T,self.F]))
        self.max_id = -1
        self.av_node_ids = deque()
        self.nodes = dict()
        self.var_pos = dict(zip(vars,[i for i in range(len(self.layers))]))
        self.nodes.update({NTRUE:self.T, NFALSE:self.F})

    ## some helper functions

    def __len__(self):
        """returns the no. of the layer-encoding variables (i.e., excluding the terminal one)"""
        return len(self.vars)

    def n(self, i):
        """returns size of the i-th layer"""

        return len(self.layers[i])

    def p(self, a):
        """returns the position (layer index) of the variable `a`"""

        return self.var_pos[a]

    def size(self):
        """returns the size of the BDD (total no. of nodes)"""
        return len(self.nodes)-2

    def new_node_name(self):
        """returns a unique name for the next node to be created

        Picks the one to re-cycle after some node deletions, when available;
        when it is not -- just takes the `max_id+1` (adjusting the `max_id` itself)
        """

        if len(self.av_node_ids)>0:
            return self.av_node_ids.popleft()
        else:
            self.max_id += 1
            return self.max_id

    def dump_gv(self,layerCapt = True):
        """Exports the BDD to the Graphviz format"""
        g = Digraph()

        for i, layer in enumerate(self.layers):
            with g.subgraph(name="cluster_{}".format(i)) as s:
                for n in layer:
                    if n.id == NTRUE:
                        s.node(str(NTRUE), label="T",fillcolor = "orange",style="filled")
                    elif n.id == NFALSE:
                        s.node(str(NFALSE),label="F", fillcolor = "gray",style="filled")
                    else:
                        s.node( str(n.id) )

                if i!=len(self.layers)-1:
                    if layerCapt:
                        s.attr(label="x{}, sz={}".format(self.vars[i], len(layer)), color="lightgrey")
                    else:
                        s.attr(color="lightgrey")
                else:
                    s.attr(color="white")

        for layer in self.layers:
            for n in layer:
                if not n.hi is None:
                    g.edge(str(n.id), str(n.hi.id))

                if not n.lo is None:
                    g.edge(str(n.id), str(n.lo.id), style="dashed")

        return g

    def addnode(self, parent_node, arc="hi", node_id=None):
        """Adds a node and updates aux info as necessary"""
        if node_id is None:
            node_id = self.new_node_name()

        newnode = node(node_id)

        if parent_node is None:
            newnode.layer = 0
            self.layers[0].add(newnode)

        elif arc=="hi":
            parent_node.hi = newnode
            self.layers[parent_node.layer+1].add(newnode)
            newnode.layer = parent_node.layer+1

        elif arc=="lo":
            parent_node.lo = newnode
            self.layers[parent_node.layer+1].add(newnode)
            newnode.layer = parent_node.layer+1

        else:
            print("ERROR: Wrong arc type: {}".format(arc))


        self.nodes.update({newnode.id: newnode})

        return newnode

    def swap_up(self, layer):
        """Swaps the layer with the one immediately above it (by _index_!)

        Arguments:
            layer: number (index) of the layer to be ''bubbled up''

        Note:
            Operation takes O (no-of-nodes in the upper layer) time.
            Drops the layer being swapped. Then iterates through the nodes of the layer
            immediately above it and creates nodes as necessary (re-cycling new nodes to avoid
            redundancy in terms of logical functions)
        """

        assert layer >= 1 # otherwise there's no layer to swap to
        assert layer <= len(self.layers)-2 # we can't swap-up the terminal layer

        new_nodes = dict()

        for n in self.layers[layer]:
            if not n.id in self.nodes.keys():
                print("Node {} not found in layer {}, among:".format(n.id,layer))
                for k in self.layers[layer]:
                    print(k.id)
                print("The set is: {}".format(self.layers[layer]))

            del self.nodes[n.id]
            self.av_node_ids.append(n.id)

        self.layers[layer] = set()

        for F in self.layers[layer-1]:
            # iterating throuhg the nodes of the upper layer

            # creating the "hi"-node
            if not ( F.hi.hi.id, F.lo.hi.id ) in new_nodes.keys():
                F_hi = node(self.new_node_name(), F.hi.hi, F.lo.hi)
                new_nodes.update({( F_hi.hi.id, F_hi.lo.id ):F_hi})
            else:
                F_hi = new_nodes[( F.hi.hi.id, F.lo.hi.id )] # re-cycle the old node

            # creating the "lo"-node
            if not ( F.hi.lo.id, F.lo.lo.id ) in new_nodes.keys():
                F_lo = node(self.new_node_name(), F.hi.lo, F.lo.lo)
                new_nodes.update({( F_lo.hi.id, F_lo.lo.id ):F_lo})
            else:
                F_lo = new_nodes[( F.hi.lo.id, F.lo.lo.id )] # re-cycle the old node

            # now, implement the change to the BDD for this source F
            F.link(F_hi,"hi")
            F.link(F_lo,"lo")

            if not F_hi.id in self.nodes.keys():
                self.nodes.update({F_hi.id: F_hi})

            if not F_lo.id in self.nodes.keys():
                self.nodes.update({F_lo.id: F_lo})

        self.layers[layer] = set(new_nodes.values())

        # NOTE: old nodes in the source layer must be (physically) deleted
        # by the Python garbage collection

        # rename layers and update aux structures as necessary
        self.var_pos.update({self.vars[layer]:layer-1, self.vars[layer-1]:layer})

        v_1 = self.vars[layer-1]
        self.vars[layer-1] = self.vars[layer]
        self.vars[layer] = v_1


    def sift(self, var, pos):
        """Sifts a variable `var`(name, not index) to position `pos`"""

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

    def gsifts(self, with_whom):
        """
        Runs a simplistic implementation of Rudell'93
        sifting alorithm extension to minimize |A|+|B|

        Note: it is in-place! no BDD copy is created
        """

        ## align with_whom to self
        N = len(self.layers)-1
        for i in range(N):
            with_whom.sift(self.vars[i],i)

        # now the pair is aligned, but possibly huge
        # let us see if we can compress it:

        processed_vars = set()

        while len(processed_vars) < N:
            i = 0
            while self.vars[i] in processed_vars:
                i += 1 # skip processed variables

            best_pos = i
            active_var = self.vars[i]
            best_size = self.size() + with_whom.size()

            for j in range(i+1,N):
                self.swap_up(j)
                with_whom.swap_up(j)
                cur_size = self.size()+with_whom.size()

                if cur_size < best_size:
                    best_size = cur_size
                    best_pos = j

            self.sift(active_var,best_pos)
            with_whom.sift(active_var,best_pos)

            processed_vars.add(active_var)

    # file manipulation procedures
    def save(self, filename):
        """saves a BDD to an ASCII(text) file

        Note:
            File format is as follows.
            File header:
                N=<no-of-layers>
                vars=<var1, var2, var3, ...> (variables' names)
            Then, one node description = one line of the format:
                id:hi,lo

            where `id` is node's ID, `hi` and `lo` are IDs(!) of the nodes corresponding
            to hi and lo pointers of the node "ID"
            The procedure performs breadth-first search and just saves all the nodes
        """

        S = deque() # a deque of node IDs to explore (FIFO)
        V = set() # set of visited (saved) nodes (by ID)

        for n in self.layers[0]:
            S.append(n.id)

        with open(filename,"w") as f:
            f.write("N={}\n".format(len(self.layers)))
            f.write("vars={}\n".format(', '.join([str(v) for v in self.vars])))

            while len(S)>0:
                n_id = S.popleft()
                # this is a new node
                if not (self.nodes[n_id].hi is None or self.nodes[n_id].hi is None):
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

    def load(self, filename):
       """loads a BDD from an ASCII(text) file

       Note:
            The format corresponds to the one described in the `save` function
       """

       with open(filename, "r") as f:
            line = f.readline()
            while line[0]=="#" or line=="":
                line = f.readline()

            assert line[:2]=='N='
            N = int(line[2:])
            assert N>1

            line = f.readline()
            assert line[:5]=='vars='
            line = line[5:].split(',')
            assert len(line) == N-1

            # initialize the attributes
            self.layers = [set() for i in range(N)]
            self.nodes = dict()
            self.vars = [v.strip() for v in line]
            self.var_pos = dict(zip(self.vars, [i for i in range(N)]))

            line = f.readline()

            current_layer = 0
            next_layer = set()
            while line:
                id, line = line.split(':')
                hi_id,lo_id = line.split(',')
                id = int(id);
                if id!=NTRUE and id!=NFALSE:
                    hi_id = int(hi_id); lo_id = int(lo_id)
                else:
                    line = f.readline()
                    continue


                if id in next_layer:
                    current_layer += 1
                    next_layer = set()

                if hi_id in next_layer:
                    F_hi = self.nodes[hi_id]
                else:
                    F_hi = node(hi_id)
                    self.nodes.update({hi_id:F_hi})
                    self.layers[current_layer+1].add(F_hi)
                    F_hi.layer = current_layer+1
                    next_layer.add(hi_id)

                if lo_id in next_layer:
                    F_lo = self.nodes[lo_id]
                else:
                    F_lo = node(lo_id)
                    self.nodes.update({lo_id:F_lo})
                    self.layers[current_layer+1].add(F_lo)
                    F_lo.layer = current_layer+1
                    next_layer.add(lo_id)

                if id in self.nodes.keys():
                    self.nodes[id].link(F_hi,"hi")
                    self.nodes[id].link(F_lo,"lo")
                else:
                    if current_layer>0:
                        print("WARNING: a node with no source at layer {}".format(current_layer))

                    F = node(id,F_hi,F_lo)
                    self.layers[current_layer].add(F)
                    self.nodes.update({id:F})

                if id == NTRUE:
                    self.T = self.nodes[id]

                if id == NFALSE:
                    self.F = self.nodes[id]

                line = f.readline()

                self.max_id = max([n.id for n in self.nodes.values()])
                self.av_node_ids = deque([i for i in range(1,self.max_id+1) if not i in self.nodes.keys()])

    @classmethod
    def random(cls, N=5,p=0.5):
        """generates a random BDD with `N` variables (N+1 layers)

        Arguments:
            N:   number of variables (to result in N+1 layers, including the T,F-layer)
            p:   tree size parameter
                    0 will generate a non-random exponential-sized diagram,
                    1 will result in a single node per layer
        Returns:
            A generated BDD
        """

        bdd = BDD(N)

        assert N>1
        bdd.addnode(parent_node=None) # create a root node

        profile = []
        for layer in range(1,N):
            # generate nodes (and edges) for the layer

            for n in bdd.layers[layer-1]:
                if rnd() <= p or len(bdd.layers[layer])==0:
                    newnode = bdd.addnode(n,arc="hi")
                else:
                    n.link(rnd_choose(tuple(bdd.layers[layer])),"hi")

                if rnd() <= p:
                    newnode = bdd.addnode(n,arc="lo")
                else:
                    n.link(rnd_choose(tuple(bdd.layers[layer])),"lo")

        # separately process the last layer
        for n in bdd.layers[-2]:
            n.link(rnd_choose(tuple(bdd.layers[-1])),"hi")
            n.link(rnd_choose(tuple(bdd.layers[-1])),"lo")

        return bdd


    def profile(self):
        """returns a BDD ``profile'' -- an (almost) BFS-ordered list of nodes

        A string of the format <vars>:<BFS-nodes>, where:
        <vars>      --  comma-separated variable names by layer;
        <BFS-nodes> --  comma-separated list of node-numbers (not IDs!)
                        in a BFS-traverse order
        """

        Q = deque()
        V = set()
        node_nums = dict();
        cur_node = 0

        p = []

        for n in self.layers[0]:
            node_nums.update({n.id:cur_node})
            Q.append(n)
            p.append(str(cur_node))
            cur_node += 1

        p = [",".join(p) + ";"]

        while(len(Q) > 0):
            n = Q.popleft()

            for successor in [n.hi, n.lo]:
                if successor is None:
                    continue

                if not successor.id in node_nums.keys():
                    node_nums.update({successor.id: cur_node})
                    cur_node += 1

                p.append(str(node_nums[successor.id]))
                if not successor.id in V:
                    Q.append(successor)
                    V.add(successor.id)

        return p[0] + ",".join(p[1:])

    def show(self, dir="testing",layerCapt=True):
        """shows the diagram (a .dot, compiled to PDF)"""
        self.dump_gv(layerCapt).view("showfunc.dot", directory=dir, cleanup=True)

    def align_to(self, vars_order):
        aligned = copy.deepcopy(self);

        for i in range(len(vars_order)):
            aligned.sift(vars_order[i],i)

        return aligned;

    def OA_bruteforce(self, with_what):
        """finds the best (smallest) alignment with `with_what` BDD by bruteforce enumeration"""

        # generate all possible permutations
        perms = iters.permutations(self.vars)
        min_size = math.inf
        A_aligned = []
        B_aligned = []
        alternatives = 0

        for perm in perms:
            A_c = self.align_to(list(perm))
            B_c = with_what.align_to(list(perm))

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
        """quickly checks if the BDD is reduced (no ``redundant'' nodes)"""

        checked_nodes = set()

        for layer in range(len(self.layers)-2,-1,-1):
            for n in self.layers[layer]:
                if (n.hi.id, n.lo.id) in checked_nodes:
                    return False
                else:
                    checked_nodes.add((n.hi.id, n.lo.id))

        return True

    def make_reduced(self):
        """quickly makes the BDD reduced

        Swaps each layer back and forth, starting from the last ones
        """

        for i in range(len(self.vars)-1,0,-1):
            self.swap_up(i)
            self.swap_up(i)

    def rename_vars(self, ren_dict):
        new_vars = [ren_dict[v] for v in self.vars]
        self.vars = new_vars
        self.var_pos = dict(zip(self.vars,[i for i in range(len(self.vars))]))

######################################################################
## testing code

## BDD creation and rendering (with graphviz)
def test_create_render():
    """tests the BDD creation and rendering code"""
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
    bdds.dump_gv().render("./testing/after_swap.dot",view=True, cleanup=True)

def test_random_swapping(N, k, m, p=0.8, mode="swaps"):
    """TEST: checks if swaps work correctly (random problems)

    Generates a random BDD, makes some random swaps (or sifts),
    and then checks that the function encoded by the BDD remained the same
    (that is, every set of var values result in the same terminal node for
    the original BDD and for the randomly changed BDD).

    Involves brute-force enumaration of all the 2^n possible decisions (concerning all vars).

    Arguments:
        N: no. of variables in the diagram
        k: no. of random BDDs generated)
        m: no. of consecutive test swaps per BDD
        p: BDD expansion prob parameter
        mode: 'swaps' or 'sifts', determines what kind of events are to be generated

    Example:
        test_random_swapping(8,500,20,0.8,sys.argv[1])
    """

    status = "OK"

    for n in range(k):
        bdd = BDD.random(N,p)
        bdd_s = copy.deepcopy(bdd)

        # make some random layer swaps
        if mode=='swaps':
            for t in range(m):
                bdd_s.swap_up(randint(1,len(bdd.vars)-1))
        elif mode=='sifts':
            for t in range(m):
                bdd_s.sift(bdd_s.vars[randint(0,len(bdd.vars)-1)],randint(0,len(bdd.vars)-1))
        else:
            print("ERROR: wrong mode -- {}. Only 'swaps' and 'sifts' are allowed".format(mode))

        # now check that everything remained the same
        no_of_bits = len(bdd.vars)

        for bit_path in range(0, 2**no_of_bits):

            # unpack the var values in the original order (1=hi, 0=lo)
            var_vals = [(bit_path >> bit) & 1 for bit in range(no_of_bits-1,-1,-1)]

            # get the correct answer (T or F)
            cur_node_src = list(bdd.layers[0])[0]
            cur_node_swd = list(bdd_s.layers[0])[0]

            i = 0
            while i<len(bdd.vars):
                if var_vals[i]==0:
                    cur_node_src = cur_node_src.lo
                else:
                    cur_node_src = cur_node_src.hi

                if var_vals[ bdd.p(bdd_s.vars[i]) ]==0:
                    cur_node_swd = cur_node_swd.lo
                else:
                    cur_node_swd = cur_node_swd.hi

                i += 1

            corr_answer = cur_node_src.id
            assert corr_answer==NTRUE or corr_answer==NFALSE
            if cur_node_swd.id != corr_answer:
                status = "\nERROR: {}, not {}".format(cur_node_swd.id, cur_node_src)
                print("path {}: {}".format(var_vals,status))

        if status!="OK":
            print("ERROR encountered!")
        else:
            print(".",end="")
    print("\n{} instances processed, status: {}".format(k,status))
