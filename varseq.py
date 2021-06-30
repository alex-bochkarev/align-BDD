"""Implements weighted variable sequences data structure.

Mostly mirrors the corresponding key methods of :py:class:`BDD.BDD`
class.

Tests coverage: :py:mod:`varseq_test`
"""
import numpy as np
import copy
import itertools as iters
import math

class VarSeq:
    """Implements weighted variable sequence data structure."""

    def __init__(self, layer_vars, layer_sizes):
        """Class constructor.

        Args:
            layer_vars (list): variable labels (strings or ints),
            layer_sizes (list): variable weights (integers)
        """
        assert len(layer_vars)==len(layer_sizes) # obviously, each layer should have a size
        # NOTE: we do not enforce n_1 = 1

        self.layer_var = copy.copy(layer_vars)
        self.n = copy.copy(layer_sizes)
        self.p = dict(zip(layer_vars, range(len(layer_vars)))) # dictionary of var positions

    # generating a random varseq
    @classmethod
    def generate_weights(cls, N):
        """Generates a random valid sequence of ``N`` weights.

        That is, respects the requirement: :math:`n_{i+1} \leq 2n_i`.
        """
        ns = [1]
        for i in range(N-1):
            ns.append(np.random.randint(1,2*ns[-1]+1))
        return np.array(ns)

    @classmethod
    def random(cls, vars = None, N = 7):
        """Generates a random variable sequence.

        Returns:
            :py:class:`VarSeq`: an object representing a weighted
            variable sequence with `N` variables and random element
            weights. If no variable labels are given (in ``vars``),
            uses a random permutation of integer labels 1,...,N.

        """
        if vars is None: vars = np.random.permutation([(i+1) for i in range(N)])
        return cls(layer_vars=vars,layer_sizes = cls.generate_weights(N))

    def size(self):
        """Returns the total sequence size (sum of element weights)."""
        return np.sum(self.n)

    def var_size(self, var):
        """Returns the element weight (given its label ``var``)."""

        return self.n[ self.p[var] ]

    def S(self, a,j):
        """Returns cost of sliding element ``a`` (label) to position ``j``."""
        if self.p[a] == j:
            return 0 # nothing to do

        if self.p[a] > j:  # "slide forward"
             return np.sum([self.n[i] for i in range(j,self.p[a])]) + self.n[j] - self.n[self.p[a]]

        # otherwise, it's p(a) < j -- "slide back"
        return (2**(j-self.p[a]+1)-1) * self.n[self.p[a]] - np.sum([self.n[i] for i in range(self.p[a],j+1)])

    # performs an actual "slide"
    # of the element a(int) to position j
    def slide(self, a, j, inplace=False):
        """Performs a sift of element ``a`` (label) to position ``j``.

        Returns:
            :py:class:`VarSeq`: a revised copy of a variable sequence
            or a self-pointer (if ``inplace`` = *True*)
        """

        if inplace:
            s = self
        else:
            s = copy.deepcopy(self)

        if self.p[a] == j:
            return s

        if self.p[a] > j: # slide forward
            for i in reversed( range(j,self.p[a]) ):
                s.n[i+1] = s.n[i]*2 # update layer sizes
                s.layer_var[i+1] = s.layer_var[i] # update layer vars
                s.p[s.layer_var[i]] += 1 # positions dictionary

            s.layer_var[j] = a;
            s.p[a] = j;

        if self.p[a] < j: # slide backward
            for i in range(self.p[a]+1, j+1):
                s.n[i] = 2*s.n[i-1] # update layer sizes
                s.layer_var[i-1] = s.layer_var[i] # update layer vars
                s.p[s.layer_var[i]] -= 1 # update positions dict

            s.layer_var[j] = a;
            s.p[a] = j;

        return s

    ## auxiliary functions
    def __str__(self):
        """Converts the object to a string (used to print a sequence)."""
        return 'Vars: {}\nn   : {} (sz={})'.format(np.array(self.layer_var),np.array(self.n),np.sum(self.n))

    def __len__(self):
        """Returns a sequence length (number of elements)"""
        return len(self.layer_var)

    def count_inversions_to(self, to_what):
        """Counts inversions with another list (or :py:class:`VarSeq`).
        """
        if type(to_what)==list:
            to_what = VarSeq(to_what,[1 for i in range(len(to_what))])

        invs = 0
        for i in range(len(self)-1):
            for j in range(i+1,len(self)):
                if to_what.p[self.layer_var[i]] >= to_what.p[self.layer_var[j]]:
                    invs += 1

        return invs

    # weighted variable sequence transformations (revisions)
    def greedy_sort(self, to_what = None):
        """Revises the :py:class:`VarSeq` to a given order of labels ``to_what``.

        Args:
            to_what (list / :py:class:`VarSeq`): target order of variables
                (or a varseq where it can be extracted from.)

        Note:
            A naive implementation / cross-check, "sift-align" procedure, :math:`O(N^2)`.
        """
        if to_what is None:
            to_what = np.sort(self.layer_var)

        galigned = copy.deepcopy(self);
        for i in range(len(self)):
            galigned = galigned.slide(to_what[i],i)
        return galigned;

    def align_to(self, to_what = None):
        """Revises the varseq to given variable order in :math:`O(N)` time.

        Args:
            to_what: a list/tuple/np.darray of target var names,
            or a variable sequence (to extract variable labels from).

        Returns:
            :py:class:`VarSeq`: a revised variable sequence

        Note:
            The original sequence in unchanged (the operation is NOT in-place!)
        """
        if to_what is None:
            to_what = np.sort(self.layer_var)

        if not (isinstance(to_what,list) or isinstance(to_what,tuple) \
                or isinstance(to_what, np.ndarray)):
            to_what = to_what.layer_var

        # to-be var indices for each variable
        # in the current sequence
        var_ind = [0 for i in range(len(self))]

        for i in range(len(to_what)):
            var_ind[self.p[to_what[i]]] = i

        # var_ind[i] now is the index of the i-th variable  in the *target* alignment

        # initalize an empty sequence (-1 to indicate un-init value)
        lsorted = VarSeq(to_what, [-1]*len(self))

        sliding_up = set()
        iB = 0
        for iA in range(len(self)):
            if not self.layer_var[iA] in sliding_up:
                ## a "sinking" element found

                lsorted.n[iB] = self.n[iA]*(2**len(sliding_up))

                while iB < var_ind[iA]:
                    sliding_up.add(lsorted.layer_var[iB])
                    iB += 1
                    lsorted.n[iB] = lsorted.n[iB-1]*2
                iB += 1

            else:
                sliding_up.remove(self.layer_var[iA])

        return lsorted

    def OA_bruteforce(self, with_what):
        """Enumerates all possible alignments with another sequence.

        Args:
            with_what (:py:class:`VarSeq`): target ordering.

        Returns:
            list: [``<no. of optima>``, ``<list of A-aligned>``, ``<list of B-aligned>``]

        Example:

            Having the returned object ``res``, ``res[0]`` gives the
            number of optima. Then, assuming it is at least 3,
            ``res[1][2]`` gives A-aligned in the third optimum,
            ``res[2][2]`` gives B-aligned in the third optimum.
        """

        # generate all possible permutations
        perms = iters.permutations(self.layer_var)
        min_size = math.inf
        A_aligned = []
        B_aligned = []
        alternatives = 0

        for perm in perms:
            A_c = self.align_to(perm)
            B_c = with_what.align_to(perm)
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

    def is_aligned(self, to_what):
        """Helper function: checks if varseq is aligned with ``to_what``."""
        return np.array_equal(self.layer_var, to_what.layer_var)

######################################################################
## Auxiliary functions (outside of the varseq class)

def non_dominated(e, A, B):
    """Checks if the element "e" is non-dominated.

    Considering ``e`` as a candidate for the last position
    in the target alignment of :py:class:`VarSeq` -s ``A`` and ``B``,
    check if it is non-dominated.

    Returns:
        *True* = non-dominated, *False* = dominated.
    """
    if A.p[e]==len(A)-1 or B.p[e]==len(B)-1:
        return True
    else:
        return set(A.layer_var[A.p[e]+1:]).isdisjoint(set(B.layer_var[B.p[e]+1:]))
